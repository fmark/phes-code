#!/usr/bin/env python

try:
    from osgeo import gdal
    from osgeo.gdalconst import *
except ImportError:
    import gdal
    from gdalconst import *

try:
    import numpy
except ImportError:
    import Numeric as numpy


import sys
import locale
import time
import csv
from itertools import product
from utils import MultiBandBlockIO
from blist import sortedlist

# =============================================================================
# Config
CACHE_BLOCKS = 48000
DEFAULT_MAX_HEIGHT_DIFF = 6.4
DEFAULT_MIN_POND_PIXELS = 4
#
POND_NO_DATA = -2**31 + 1
# =============================================================================
def Usage():
    print('Usage: build_ponds.py indem inpixelpairs outraster outcsv\n')
    sys.exit( 1 )

# =============================================================================

def is_no_data(nd, val):
    if isinstance(nd, float):
        return nd == float(val)
    if isinstance(nd, int):
        return nd == int(val)
    if isinstance(nd, long):
        return nd == long(val)
    print("Unsupported no data type, %s of %s" % (
            str(nd), str(type(nd))))
    sys.exit(1)

# Find all 8 neighbours 

def fringe8((px, py), extent=None):
    nonzero = (pair for pair in product((-1,0,1),(-1,0,1)) if pair != (0, 0))
    f = ((px+dx, py+dy) for (dx,dy) in nonzero)
    if extent is None:
        return f
    x1, y1, x2, y2 = extent
    return [(fx, fy) for (fx, fy) in f if x1 <= fx < x2 and y1 <= fy < y2]

def find_closest_list_head(lists, val_func, compare_val):
    closest = None
    closest_diff = max_height_diff + 1
    for list in lists:
        try:
            diff = abs(val_func(list) - compare_val)
#            print diff
            if diff < closest_diff:
                closest = list
                closest_diff = diff
        except IndexError:
            continue
    return closest

def rebalance_lists(lower, upper, pivot):
    while len(lower) > 0 and lower[0][1] > pivot: 
        upper.add(lower.pop())
    while (len(upper) > 0) and (upper[0][1] < pivot):
        lower.add(upper.pop())


# Parse command line arguments.
if len(sys.argv) != 5 and len(sys.argv) != 6:
    Usage()

indem = sys.argv[1]
inpxpairs = sys.argv[2]
outfile = sys.argv[3]
outcsv = sys.argv[4]
if len(sys.argv) == 6:
    max_height_diff = float(sys.argv[5])
else:
    max_height_diff = DEFAULT_MAX_HEIGHT_DIFF

print("Maximum height difference %.2f" % max_height_diff)

indemds = gdal.Open(indem, GA_ReadOnly)
if indemds is None:
    print("Input file \"%s\" could not be found." % indem)
    sys.exit(1)
inpxpairsds = gdal.Open(inpxpairs, GA_ReadOnly)
if inpxpairsds is None:
    print("Input file \"%s\" could not be found." % inpxpairs)
    sys.exit(1)

out_driver = indemds.GetDriver()
outdataset = out_driver.Create(outfile, indemds.RasterXSize, indemds.RasterYSize, indemds.RasterCount, GDT_Int32)
if outdataset is None:
    print("Could not create file \"%s\"." % outfile)
    sys.exit(1)

outdataset.SetProjection(indemds.GetProjection())
outdataset.SetGeoTransform(indemds.GetGeoTransform())

indemband = indemds.GetRasterBand(1)
inpxpairsband = inpxpairsds.GetRasterBand(1)
outband = outdataset.GetRasterBand(1)

# Get some stats
adf = indemds.GetGeoTransform()
xcellsize, ycellsize = (abs(adf[1]), 
                        abs(adf[5]))
dem_no_data = indemband.GetNoDataValue()
pxpairs_no_data = inpxpairsband.GetNoDataValue()

outband.SetNoDataValue(POND_NO_DATA)

print "Initialising empty new raster."
# Cast to float, as band. From the docs, "The fill value is 
# passed in as a double but this will be converted to the 
# underlying type before writing to the file. "
outband.Fill(float(POND_NO_DATA))
io = MultiBandBlockIO((indemband, inpxpairsband, outband), CACHE_BLOCKS, True)

total_pixel_count = indemband.XSize * indemband.YSize
locale.setlocale(locale.LC_ALL, "")
start_time = time.time()
world_extent = (0, 0, indemband.XSize, indemband.YSize)
pond_num = -1
for (demblock, pxpairsblock, outblock), block_x, block_y, world_x, world_y in ( 
    io.extent_iterator(*world_extent)):
    # numpy arrays are accessed array[y][x]
    pair_val = pxpairsblock.data[block_y][block_x]
    # skip empty pixel pairs
    if is_no_data(pxpairs_no_data, pair_val):
        continue

    out_val = outblock.data[block_y][block_x]
    # if the pixel has been assigned to a pond already, skip it
    if not is_no_data(POND_NO_DATA, out_val):
        continue
    elevation = demblock.data[block_y][block_x]    
    new_pixel = (world_x, world_y)
    # Make a pond
    pond_num += 1
    pond_pixels = set()
    lower_pond = sortedlist(key=(lambda x: -1 * x[1]))
    upper_pond = sortedlist(key=(lambda x: x[1]))
    equal_pond = [(new_pixel, elevation)]
    pond_pixels.add(new_pixel)
    # Keep stats on contents of pond
    pond_size = 1
    pond_mean_elevation = elevation
    pond_sum_elevation = elevation
    print "Added pixel 1, elevation %f." % pond_sum_elevation
    pond_total_diff = 0.0
    lower_fringe = sortedlist(key=(lambda x: -1 * x[1]))
    upper_fringe = sortedlist(key=(lambda x: x[1]))
    fringe_pixels = set()
    bad_fringe = set()
    # there is no do while, so we'll manually break out of the loop
#    print "New pond, number %4d" % pond_num
    while True: 
        # add surrounding 8 pixels to the fringe if valid
        for f in fringe8(new_pixel, world_extent):
            # screen out any pixels we don't want in the fringe
            if f in fringe_pixels or f in bad_fringe or f in pond_pixels:
                continue
            f_demval, f_pairval, f_outval = io.get_pixel(f)
            # if its not in a pair, skip it
            if is_no_data(pxpairs_no_data, f_pairval):
                bad_fringe.add(f)
                continue
            # if pixel is already in a pond, skip it
            if not is_no_data(POND_NO_DATA, f_outval):
                bad_fringe.add(f)
                continue
            # add valid pixel to fringe
            fringe_pixels.add(f)
            if f_demval > pond_mean_elevation:
                upper_fringe.add((f, f_demval))
            else:
                lower_fringe.add((f, f_demval))
        closest_fringe = find_closest_list_head(
            [upper_fringe, lower_fringe],
            (lambda x: x[0][1]), 
            pond_mean_elevation)

        # if the fringes are empty, then we have finished this pond
        if closest_fringe is None:
            break
        new_pixel, new_elevation = closest_fringe.pop()
        fringe_pixels.remove(new_pixel)
        # if the closest pixel is too far away, then we are done
        diff = new_elevation - pond_mean_elevation
        if pond_total_diff + abs(diff) > (
            (pond_size + 1) * max_height_diff):
            break
        # now recalculate pond statistics
        pond_size += 1
        pond_sum_elevation += new_elevation
        old_mean_elevation = pond_mean_elevation
        pond_mean_elevation = float(pond_sum_elevation) / pond_size
        print "Added pixel %d, elevation %f.  Total diff %.2f, old mean elevation %.2f, new mean elevation %.2f." % (
            pond_size, new_elevation, pond_total_diff, old_mean_elevation, pond_mean_elevation)
        #add the new pixel
        pond_pixels.add(new_pixel)
        pond_total_diff += abs(pond_mean_elevation - new_elevation)
        print "Added new pixel's elevation diff %.2f, total diff %.2f" % (
            abs(pond_mean_elevation - new_elevation), pond_total_diff)
        # do we need to move any pixels between ponds?
        # set aside equal pond to deal with later
        old_equal_pond = equal_pond
        equal_pond = []
        # do any required pond rebalancing
        if new_elevation > pond_mean_elevation:
            # print "Elevation rising"
            print upper_pond, lower_pond
#            pond_total_diff += len(lower_pond) * abs(pond_mean_elevation - old_mean_elevation)
#            print "Added %d lower pond member's elevation to diff, total diff %.2f" % (
#                len(lower_pond), pond_total_diff)
            # transfer any from upper pond to lower or equal ponds if the rising mean has trapped them
            while len(upper_pond) > 0:
                if upper_pond[0][1] < pond_mean_elevation:
                    print "Need to move"
                break
        # which pond to add the new pixel to?
        if new_elevation == pond_mean_elevation:
            equal_pond.append(((new_pixel), new_elevation))
            # restore the old_equal pond as we haven't changed mean
            equal_pond.extend(old_equal_pond)
        else:
            # Add the contents of the old equal pond to the diff
            pond_total_diff += len(old_equal_pond) * abs(old_mean_elevation - pond_mean_elevation)
            print "Added %d old equal pond member's elevation to diff, total diff %.2f" % (
                len(old_equal_pond), pond_total_diff)
            if new_elevation > pond_mean_elevation:
                upper_pond.add(((new_pixel), new_elevation))
                # add the old equal to lower
                lower_pond |= old_equal_pond
            else:
                lower_pond.add(((new_pixel), new_elevation))
                upper_pond |= old_equal_pond
        # validate
        d_tot = 0.0
        d_i = 0
        for _, e in upper_pond:
            d_tot += e - pond_mean_elevation
            d_i += 1
        for _, e in lower_pond:
            d_tot += pond_mean_elevation - e
            d_i += 1
        for _, e in equal_pond:
            d_tot += abs(pond_mean_elevation - e)
            d_i += 1
        d_mean = d_tot / d_i
        print "Counted %d pixels to validate" % d_i
        if abs(pond_total_diff / pond_size - d_mean) > 0.001:
            print "Fail: %f != %f   (total diff %.2f, target %.2f, diff diff %.6f)" % (
                pond_total_diff / pond_size, d_mean, pond_total_diff, d_tot, abs(d_tot - pond_total_diff))
            sys.exit(1)
        else:
            print "Pass   -   total diff %.2f" % pond_total_diff
    print "Fin pond\n"
    # end of for loop.  Following code is kept for reference
    continue
    outblock.data[block_j][block_i] |= numpy.uint8(BOTTOM_RESERVOIR)
    outblock.dirty = True
    origin_outblock.data[origin_block_iy][origin_block_ix] |= (
        numpy.uint8(TOP_RESERVOIR))
    origin_outblock.dirty = True
    
# write cached dirty blocks to disk
io.write_flush()
now = time.time()
print "Search complete, in %.2f minutes." % (
    (now - start_time) / 60)

