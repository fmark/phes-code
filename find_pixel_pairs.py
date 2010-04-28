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
PIXEL_SIZE = 250
DEM_NO_DATA = -3.4028235e+38
CACHE_BLOCKS = 32000
MAX_SLOPE = 15.0
# =============================================================================
def Usage():
    print('Usage: find_pixel_pairs.py indem inslopeslope outfile')
    print('')
    sys.exit( 1 )

# =============================================================================


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
print "Initialising empty new raster to 0."
outband.Fill(0.0, 0.0)
io = MultiBandBlockIO((indemband, inslopeband, outband), CACHE_BLOCKS, True)

search_radius_constant = 10.0 / PIXEL_SIZE
count = 0
total_count = indemband.XSize * indemband.YSize
locale.setlocale(locale.LC_ALL, "")
start_time = time.time()
for y in xrange(indemband.YSize - 1, -1, -1):
    scanline_dem = indemband.ReadAsArray(0, y, indemband.XSize, 1, indemband.XSize, 1)[0]
    scanline_slope = inslopeband.ReadAsArray(0, y, indemband.XSize, 1, indemband.XSize, 1)[0]
    for x in xrange(0, len(scanline_dem)):
        count += 1
        if count % 1000 == 0:
            s = ("Scanned neighbourhood in %s of %s pixels (%.2f%%). " %
                 (locale.format('%d', count, True), 
                  locale.format('%d', total_count, True), 
                  float(count)/total_count * 100))
            if count >= 10000:
                now = time.time()
                seconds_remaining = ((total_count - count) * 
                                     ((now - start_time) / count))
                fin = time.localtime(now + seconds_remaining)
                s += ("ETA %2d:%2d, %2d/%02d/%d." % 
                      (fin.tm_hour, fin.tm_min, fin.tm_mday,
                      fin.tm_mon, fin.tm_year))
            print s
        elevation = scanline_dem[x]
        origin_slope = scanline_slope[x]
        if elevation == numpy.float32(DEM_NO_DATA):
            continue
        if (origin_slope == numpy.float32(DEM_NO_DATA) or
            origin_slope > numpy.float32(MAX_SLOPE)):
            continue
        # elevation is in meters.  We want a 1 in 10 slope, so there is no point
        # searching farther than 10 times the elevation (unless we have
        # land below sea level).
        search_radius = elevation * search_radius_constant + 1
        extent = io.search_extent(x, y, search_radius)
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
                        if outblock.data[block_j][block_i] != (
                            numpy.uint8(1) ):
                            # if slope is too great, then set 0
                            slope = slopeblock.data[block_j][block_i]
                            if slope <= MAX_SLOPE:
                                # Look for downhill gradient of >= 0.1 
                                # for a 1 in 10 slope
                                new_e = demblock.data[block_j][block_i]
                                elevation_delta = elevation - new_e
                                if elevation_delta > 0:
                                    new_x, new_y = io.block_to_pixel_coord(
                                        block_x, block_y)
                                    new_x += block_i
                                    new_y += block_j
                                    dist = (((x - new_x) * PIXEL_SIZE)**2 + 
                                            ((y - new_y) * PIXEL_SIZE)**2)
                                    gradient = (float(elevation_delta) / 
                                                numpy.sqrt(dist))
                                    
                                    if gradient >= 0.1:
                                        outblock.data[block_j][block_i] = 1
                                        outblock.dirty = True
                # Try to ensure memory is released
                demblock, outblock, slopeblock = None, None, None
                del demblock, outblock, slopeblock
        #      
print (x, y, elevation, search_radius, extent, cur_block_extent, subset)
io.write_flush()
