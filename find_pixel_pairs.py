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
from utils import MultiBandBlockIO

# =============================================================================
# Config
CACHE_BLOCKS = 48000
MAX_SLOPE = 1.5
# Classes
# TOP and BOTTOM must be perfect powers of two, e.g. 1, 2, 4, 8, 16, 32, etc.
TOP_RESERVOIR = 1
BOTTOM_RESERVOIR = 2
BOTH_RESERVOIRS = TOP_RESERVOIR | BOTTOM_RESERVOIR
NEITHER_RESERVOIR = 0
# =============================================================================
def Usage():
    print('Usage: find_pixel_pairs.py indem inslopeslope outfile\n')
    sys.exit( 1 )

# =============================================================================

def is_no_data(nd, val):
    if isinstance(nd, float):
        return nd == float(val)
    if isinstance(nd, int):
        return nd == int(val)
    print("Unsupported no data type, %s of %s" % (
            str(nd), str(type(nd))))
    sys.exit(1)
       

# Parse command line arguments.
if len(sys.argv) != 4:
    Usage()

indem = sys.argv[1]
inslope = sys.argv[2]
outfile = sys.argv[3]

indemds = gdal.Open(indem, GA_ReadOnly)
if indemds is None:
    print("Input file \"%s\" could not be found." % indem)
    sys.exit(1)
inslopeds = gdal.Open(inslope, GA_ReadOnly)
if inslopeds is None:
    print("Input file \"%s\" could not be found." % inslope)
    sys.exit(1)

out_driver = indemds.GetDriver()
outdataset = out_driver.Create(outfile, indemds.RasterXSize, indemds.RasterYSize, indemds.RasterCount, GDT_Byte)
if outdataset is None:
    print("Could not create file \"%s\"." % outfile)
    sys.exit(1)

outdataset.SetProjection(indemds.GetProjection())
outdataset.SetGeoTransform(indemds.GetGeoTransform())

indemband = indemds.GetRasterBand(1)
inslopeband = inslopeds.GetRasterBand(1)
outband = outdataset.GetRasterBand(1)

# Get some stats
adf = indemds.GetGeoTransform()
xcellsize, ycellsize = (abs(adf[1]), 
                        abs(adf[5]))
dem_no_data = indemband.GetNoDataValue()
slope_no_data = inslopeband.GetNoDataValue()

outband.SetNoDataValue(NEITHER_RESERVOIR)

print "Initialising empty new raster."
# Cast to float, as band. From the docs, "The fill value is 
# passed in as a double but this will be converted to the 
# underlying type before writing to the file. "
outband.Fill(float(NEITHER_RESERVOIR))
io = MultiBandBlockIO((indemband, inslopeband, outband), CACHE_BLOCKS, True)

x_search_dist_constant = 10.0 / xcellsize
y_search_dist_constant = 10.0 / ycellsize

match_count = 0
pixel_count = 0
total_pixel_count = indemband.XSize * indemband.YSize
locale.setlocale(locale.LC_ALL, "")
start_time = time.time()
for (origin_demblock, origin_slopeblock, origin_outblock), origin_block_ix, origin_block_iy, x, y in ( 
    io.extent_iterator(0, 0, indemband.XSize, indemband.YSize)):

                # progress bar related stuff
                pixel_count += 1
                if pixel_count % 1000 == 0:
                    s = ("Scanned %s of %s pixels (%.2f%%), total %d matches. " %
                         (locale.format('%d', pixel_count, True), 
                          locale.format('%d', total_pixel_count, True), 
                          float(pixel_count)/total_pixel_count * 100,
                          match_count))
                    if pixel_count >= 10000:
                        now = time.time()
                        seconds_remaining = ((total_pixel_count - pixel_count) * 
                                             ((now - start_time) / pixel_count))
                        fin = time.localtime(now + seconds_remaining)
                        s += ("ETA %2d:%2d, %2d/%02d/%d." % 
                              (fin.tm_hour, fin.tm_min, fin.tm_mday,
                               fin.tm_mon, fin.tm_year))
                    print s

                # numpy arrays are accessed array[y][x]
                elevation = origin_demblock.data[origin_block_iy][origin_block_ix]
                origin_slope = origin_slopeblock.data[origin_block_iy][origin_block_ix]
                if is_no_data(dem_no_data, elevation):
                    continue
                if is_no_data(slope_no_data, origin_slope) or (
                    float(origin_slope) > MAX_SLOPE):
                    continue
                # elevation is in meters.  We want a 1 in 10 slope, so there is no point
                # searching farther than 10 times the elevation (unless we have
                # land below sea level).
                
                x_search_dist = elevation * x_search_dist_constant + 1
                y_search_dist = elevation * y_search_dist_constant + 1
                extent = io.search_extent_rect(x, y, x_search_dist, y_search_dist)
                
                # we want to iterate over the search area one block at a time in order
                # to minimise the amount of IO.  First we find which blocks to search:
                block_extent = io.find_block_extent(extent)
                for block_x in xrange(block_extent[0], block_extent[2] + 1):
                    for block_y in xrange(block_extent[1], block_extent[3] + 1):
                        # Get this block in the registered bands
                        demblock, slopeblock, outblock = (
                            io.get_block_in_bands(block_x, block_y) )
                        # Find the subset of pixels in this block that area within
                        # the search radius
                        cur_block_extent = io.block_to_pixel_extent(
                            block_x, block_y)
                        subset = io.intersection((extent, cur_block_extent))
                        # translate the subset coords back to the coords of the 
                        # array within the block
                        xoff = subset[0] - cur_block_extent[0]
                        yoff = subset[1] - cur_block_extent[1]
                        xcount = subset[2] - subset[0]
                        ycount = subset[3] - subset[1]
                        # Blocks are reverse indexed, i.e. instead of being
                        # block[x][y] they are block[y][x]
                        for block_j in xrange(yoff, yoff + ycount):
                            for block_i in xrange(xoff, xoff + xcount):
                                # For each test_pixel in the search radius:
                                #   if output is already 1 (matched), then continue
                                if (outblock.data[block_j][block_i] == 
                                    numpy.uint8(NEITHER_RESERVOIR) or 
                                    origin_outblock.data[origin_block_iy][origin_block_ix] == (
                                        numpy.uint8(NEITHER_RESERVOIR))):
                                    # if slope is too great, then set 0
                                    slope = slopeblock.data[block_j][block_i]
                                    if slope <= MAX_SLOPE:
                                        # Look for downhill gradient of >= 0.1 
                                        # for a 1 in 10 slope
                                        new_e = demblock.data[block_j][block_i]
                                        elevation_delta = elevation - new_e
                                        if elevation_delta > 100:
                                            new_x, new_y = io.block_to_pixel_coord(
                                                block_x, block_y)
                                            new_x += block_i
                                            new_y += block_j
                                            dist = numpy.sqrt(
                                                (((x - new_x) * xcellsize)**2) + 
                                                (((y - new_y) * ycellsize)**2))
                                            gradient = (float(elevation_delta) / 
                                                        float(dist))
                                            if gradient >= 0.1:
                                                outblock.data[block_j][block_i] |= numpy.uint8(BOTTOM_RESERVOIR)
                                                outblock.dirty = True
                                                origin_outblock.data[origin_block_iy][origin_block_ix] |= (
                                                    numpy.uint8(TOP_RESERVOIR))
                                                origin_outblock.dirty = True
                                                match_count += 1
                        # Try to ensure memory is released
                        demblock, outblock, slopeblock = None, None, None
                        del demblock, outblock, slopeblock

# write cached dirty blocks to disk
io.write_flush()
now = time.time()
print "Search complete, found %d candidate pixel pairs in %.2f minutes." % (
    match_count, (now - start_time) / 60)

