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
from utils import MultiBandBlockIO
from orderedset import OrderedSet

# =============================================================================
# Config
CACHE_BLOCKS = 48000
DEFAULT_MAX_HEIGHT_DIFF = 2.0
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
    print("Unsupported no data type, %s of %s" % (
            str(nd), str(type(nd))))
    sys.exit(1)

# Parse command line arguments.
if len(sys.argv) != 5 and len(sys.argv) != 6:
    Usage()

indem = sys.argv[1]
inpxpairs = sys.argv[2]
outfile = sys.argv[3]
outcsv = sys.argv[4]
if len(sys.argv != 6):
    max_height_diff = float(sys.argv[5])
else:
    max_height_diff = DEFAULT_MAX_HEIGHT_DIFF

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
for (demblock, pxpairsblock, outblock), block_x, block_y, world_x, world_y in ( 
    io.extent_iterator(0, 0, indemband.XSize, indemband.YSize)):
    # numpy arrays are accessed array[y][x]
    pair_val = pxpairsblock.data[block_y][block_x]
    # skip empty pixel pairs
    if is_no_data(pxpairs_no_data, pair_val):
        continue
    out_val = outblock.data[block_y][block_x]
    # if the pixel has been assigned to a pond already, skip it
    if not is_no_data(POND_NO_DATA, out_val):
        continue
    
    # Make an ordered set
    pond = [(world_x, world_y)]
    pond_size = len(pond)
    while True:
        fringe = get_fringe(pond)
        # add from fringe to pond if okay
        # for f in fringe:
        #   if f is okay:
        #     pond += f
        # if we haven't added to our pond this iteration, break out
        if len(pond) = pond_size:
            break;
        else:
            pond_size = len(pond)
    
    x_search_dist = elevation * x_search_dist_constant + 1
    y_search_dist = elevation * y_search_dist_constant + 1
    extent = io.search_extent_rect(x, y, x_search_dist, y_search_dist)

    # we want to iterate over the search area one block at a time in order
    # to minimise the amount of IO.  First we find which blocks to search:
    for (demblock, pxpairsblock, outblock), block_i, block_j, new_x, new_y in ( 
        io.extent_iterator(*extent)):
        # For each test pixel in the search radius:
        #   if output is already 1 (matched), then continue
        if (outblock.data[block_j][block_i] == 
            numpy.uint8(NEITHER_RESERVOIR) or 
            origin_outblock.data[origin_block_iy][origin_block_ix] == (
                numpy.uint8(NEITHER_RESERVOIR))):
            # if slope is too great, then set 0
            pxpairs = pxpairsblock.data[block_j][block_i]
            if pxpairs <= MAX_SLOPE:
                # Look for downhill gradient of >= 0.1 
                # for a 1 in 10 pxpairs
                new_e = demblock.data[block_j][block_i]
                elevation_delta = elevation - new_e
                # Minimum 100m head required for PHES
                if elevation_delta > 100:
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
    demblock, outblock, pxpairsblock = None, None, None
    del demblock, outblock, pxpairsblock

# write cached dirty blocks to disk
io.write_flush()
now = time.time()
print "Search complete, found %d candidate pixel pairs in %.2f minutes." % (
    match_count, (now - start_time) / 60)

