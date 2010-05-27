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
from utils import *
from blist import sortedlist

# =============================================================================
# Config
CACHE_BLOCKS = 48000
DEFAULT_MAX_HEIGHT_DIFF   = 6.4
DEFAULT_MIN_POND_SIZE_SQM = 32400

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

def main():

    indem = sys.argv[1]
    inpxpairs = sys.argv[2]
    outfile = sys.argv[3]
    outcsv = sys.argv[4]
    if len(sys.argv) > 5:
        max_height_diff = float(sys.argv[5])
    else:
        max_height_diff = DEFAULT_MAX_HEIGHT_DIFF

    if len(sys.argv) > 6:
        min_pond_size_sqm = float(sys.argv[6])
    else:
        min_pond_size_sqm = DEFAULT_MIN_POND_SIZE_SQM


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
    xcellsize, ycellsize = get_cellsize(indemds)
    min_pond_pixels = int(ceil_div(min_pond_size_sqm, (xcellsize * ycellsize)))
    print "Using minimum pond size of %.2f sq.m, which requires %d pixels of %.2f x %.2f m." % (
        min_pond_size_sqm, min_pond_pixels, xcellsize, ycellsize)

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
    pond_num = 0
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
        pond_pixels = {new_pixel: elevation}
        # Keep stats on contents of pond
        pond_size = 1
        mean_diff = 0.0
        pond_mean_elevation = elevation
        pond_sum_elevation = elevation
    #    print "Added pixel 1, elevation %f." % pond_sum_elevation
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
                if (f in fringe_pixels or 
                    f in bad_fringe or 
                    pond_pixels.has_key(f)):
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
            # now recalculate pond statistics
            pond_size += 1
            pond_sum_elevation += new_elevation
            pond_mean_elevation = float(pond_sum_elevation) / pond_size
            pond_pixels[new_pixel] = new_elevation
            # calculate total pond earth to move
            old_mean_diff = mean_diff
            total_diff = sum(map(lambda x: abs(x - pond_mean_elevation), pond_pixels.values()))
            mean_diff = total_diff / mean_diff
            # If 
            if mean_diff > max_height_diff:
                del pond_pixels[new_pixel]
                mean_diff = old_mean_diff
                break

        if pond_size >= min_pond_pixels:
            print "Found pond %d, %d pixels\n" % (pond_num, pond_size)        
            #save the found pond to buffer
            for save_pixel, save_elevation in pond_pixels.items():
                io.set_pixel(outband, save_pixel, pond_num)
            pond_num += 1
        # end of for loop.  Following code is kept for reference
        continue

    # write cached dirty blocks to disk
    io.write_flush()
    now = time.time()
    print "Search complete, in %.2f minutes." % (
        (now - start_time) / 60)



# Parse command line arguments.
if not (4 < len(sys.argv) < 8):
    Usage()
else:
    main()
