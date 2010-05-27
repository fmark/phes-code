#!/usr/bin/env python

from osgeo import gdal, ogr, osr
from osgeo.gdalconst import *
import numpy
import sys
import locale
import time
import os.path
from itertools import product
from utils import *
from blist import sortedlist


# =============================================================================
# Config
CACHE_BLOCKS = 48000
DEFAULT_MAX_HEIGHT_DIFF   = 6.4
DEFAULT_MIN_POND_SIZE_SQM = 32400

POND_NO_DATA = 0
# =============================================================================
def Usage():
    print('Usage: build_ponds.py indem inpixelpairs outraster outshp\n')
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

#repartition already partitioned, sorted lists when the pivot has changed
def repartition_lists(lower, upper, pivot):
    while len(lower) > 0 and lower[0][1] > pivot: 
        upper.add(lower.pop())
    while (len(upper) > 0) and (upper[0][1] < pivot):
        lower.add(upper.pop())


def find_closest_list_head(lists, val_func, compare_val):
    closest = None
    closest_diff = float("inf")
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

class PondBuilder():
    def parse_args(self):
        self.indem = sys.argv[1]
        self.inpxpairs = sys.argv[2]
        self.outfile = sys.argv[3]
        self.outshp = sys.argv[4]
        if len(sys.argv) > 5:
            self.max_height_diff = float(sys.argv[5])
        else:
            self.max_height_diff = DEFAULT_MAX_HEIGHT_DIFF

        if len(sys.argv) > 6:
            self.min_pond_size_sqm = float(sys.argv[6])
        else:
            self.min_pond_size_sqm = DEFAULT_MIN_POND_SIZE_SQM


        print("Maximum height difference %.2f" % self.max_height_diff)

        self.indemds = gdal.Open(self.indem, GA_ReadOnly)
        if self.indemds is None:
            print("Input file \"%s\" could not be found." % self.indem)
            sys.exit(1)
        self.inpxpairsds = gdal.Open(self.inpxpairs, GA_ReadOnly)
        if self.inpxpairsds is None:
            print("Input file \"%s\" could not be found." % self.inpxpairs)
            sys.exit(1)

        out_driver = self.indemds.GetDriver()
        self.outdataset = out_driver.Create(self.outfile, self.indemds.RasterXSize, self.indemds.RasterYSize, self.indemds.RasterCount, GDT_Int32)
        if self.outdataset is None:
            print("Could not create file \"%s\"." % outfile)
            sys.exit(1)

        self.outdataset.SetProjection(self.indemds.GetProjection())
        self.outdataset.SetGeoTransform(self.indemds.GetGeoTransform())

        self.indemband = self.indemds.GetRasterBand(1)
        self.inpxpairsband = self.inpxpairsds.GetRasterBand(1)
        self.outband = self.outdataset.GetRasterBand(1)



    def __init__(self):
        self.parse_args()

    def prepare_no_data(self):
        self.dem_no_data = self.indemband.GetNoDataValue()
        self.pxpairs_no_data = self.inpxpairsband.GetNoDataValue()

        self.outband.SetNoDataValue(POND_NO_DATA)

        print "Initialising empty new raster."
        # Cast to float, as band. From the docs, "The fill value is 
        # passed in as a double but this will be converted to the 
        # underlying type before writing to the file. "
        self.outband.Fill(float(POND_NO_DATA))

    # we can skip a pixel if it's not in a pair (not useful)
    # or it has already been assigned to a pond
    def can_skip_pixel(self, pairblock, outblock, block_x, block_y):
        # numpy arrays are accessed array[y][x]
        pair_val = pairblock.data[block_y][block_x]
        # skip empty pixel pairs
        if is_no_data(self.pxpairs_no_data, pair_val):
            return True
        out_val = outblock.data[block_y][block_x]
        # if the pixel has been assigned to a pond already, skip it
        return (not is_no_data(POND_NO_DATA, out_val))

    def delete_shapefile(self, path):
        basename = os.path.splitext(path)[0]
        for ext in ['.shp', '.dbf', '.prj', '.shx']:
            try:
                os.remove(basename + ext)
            except:
                pass

    def polygonize(self):
        print ("Converting raster to %d polygons." % self.pond_num);
        shapefile_driver = ogr.GetDriverByName('ESRI Shapefile')
        self.delete_shapefile(self.outshp)
        self.shapefile_ds = shapefile_driver.CreateDataSource(self.outshp)
        out_srs = None
        if self.outdataset.GetProjectionRef() != '':
            out_srs = osr.SpatialReference()
            out_srs.ImportFromWkt(self.outdataset.GetProjectionRef())
        self.shapefile_polygons_layer = self.shapefile_ds.CreateLayer('ponds', srs = out_srs)

        fields = [('id', ogr.OFTInteger),
                  ('size',   ogr.OFTReal),
                  ('mean_h', ogr.OFTReal),
                  ('earth_m3', ogr.OFTReal)]
        for name, type in fields:
            fd = ogr.FieldDefn(name, type)
            self.shapefile_polygons_layer.CreateField( fd )
        gdal.Polygonize(self.outband, self.outband, self.shapefile_polygons_layer, 0)

    def add_polygon_attributes(self):
        # update attributes
        self.shapefile_polygons_layer.ResetReading()
        feat = self.shapefile_polygons_layer.GetNextFeature()
        while feat is not None:
            id = feat.GetFieldAsInteger(0)
            feat_size, feat_mean_h, feat_mean_diff = self.pond_attributes[id]
            feat.SetField('size', feat_size)
            feat.SetField('mean_h', feat_mean_h)
            feat.SetField('earth_m3', feat_mean_diff)
            self.shapefile_polygons_layer.SetFeature(feat)
            feat = self.shapefile_polygons_layer.GetNextFeature()



    def main(self):
        # Get some stats
        self.xcellsize, self.ycellsize = get_cellsize(self.indemds)
        self.min_pond_pixels = int(ceil_div(self.min_pond_size_sqm, (self.xcellsize * self.ycellsize)))
        print "Using minimum pond size of %.2f sq.m, which requires %d pixels of %.2f x %.2f m." % (
            self.min_pond_size_sqm, self.min_pond_pixels, self.xcellsize, self.ycellsize)
        self.prepare_no_data()
        
        self.io = MultiBandBlockIO((self.indemband, self.inpxpairsband, self.outband), CACHE_BLOCKS, True)

        self.total_pixel_count = self.indemband.XSize * self.indemband.YSize

        locale.setlocale(locale.LC_ALL, "")
        start_time = time.time()
        self.world_extent = (0, 0, self.indemband.XSize, self.indemband.YSize)
        self.pond_num = POND_NO_DATA + 1
        self.pond_attributes = {}
        for (demblock, pxpairsblock, outblock), block_x, block_y, world_x, world_y in ( 
            self.io.extent_iterator(*self.world_extent)):
            if self.can_skip_pixel(pxpairsblock, outblock, block_x, block_y):
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
                for f in fringe8(new_pixel, self.world_extent):
                    # screen out any pixels we don't want in the fringe
                    if (f in fringe_pixels or 
                        f in bad_fringe or 
                        pond_pixels.has_key(f)):
                        continue
                    f_demval, f_pairval, f_outval = self.io.get_pixel(f)
                    # if its not in a pair, skip it
                    if is_no_data(self.pxpairs_no_data, f_pairval):
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
                mean_diff = total_diff / pond_size
                # If 
                if mean_diff > self.max_height_diff:
                    del pond_pixels[new_pixel]
                    mean_diff = old_mean_diff
                    break
                repartition_lists(lower_fringe, upper_fringe, pond_mean_elevation)

            if pond_size >= self.min_pond_pixels:
#                print "Found pond %d, %d pixels" % (self.pond_num, pond_size)        
                #save the found pond to buffer
                for save_pixel, save_elevation in pond_pixels.items():
                    self.io.set_pixel(self.outband, save_pixel, self.pond_num)
                pond_m2 = pond_size * self.xcellsize * self.ycellsize
                self.pond_attributes[self.pond_num] = (pond_m2,
                                                  pond_mean_elevation,
                                                  total_diff * pond_m2)
                self.pond_num += 1
            # end of for loop.  Following code is kept for reference
            continue

        # write cached dirty blocks to disk
        self.io.write_flush()
        now = time.time()
        print "Search complete, in %.2f minutes." % (
            (now - start_time) / 60)
        polygonize_start_time = now
        self.polygonize()
        now = time.time()
        print "Polygonization complete, in %.2f minutes." % (
            (now - polygonize_start_time) / 60)
        attributes_start_time = now
        self.add_polygon_attributes()
        now = time.time()
        print "Adding attributes complete, in %.2f minutes." % (
            (now - attributes_start_time) / 60)
        print "Total build time %.2f mintutes." % (
            (now - start_time) / 60)

# Parse command line arguments.
if not (4 < len(sys.argv) < 8):
    Usage()
else:
    bp = PondBuilder()
    bp.main()
