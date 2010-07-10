#!/usr/bin/env python

from osgeo import gdal, ogr, osr
from osgeo.gdalconst import *
import numpy
import sys
import locale
import time
import os.path
from itertools import product
from collections import namedtuple
from blist import sortedlist
from utils import *



# =============================================================================
# Config
CACHE_BLOCKS = 48000
DEFAULT_MAX_HEIGHT_DIFF   = 6.4
DEFAULT_MIN_POND_SIZE_SQM = 32400
MAX_POND_SIZE = 5500

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
def neighbouring_pixels((px, py), n=8, extent=None):
    if n == 8:
        nonzero = (pair for pair in product((-1,0,1),(-1,0,1)) if pair != (0, 0))
    elif n == 4:
        nonzero = [(1, 0), (0, 1), (-1, 0), (0, -1)]
    else:
        raise ValueError("Can only get 8 or 4 neighbours.")
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

Pixel = namedtuple('Pixel', 'x y elevation n_neighbours is_site')

class Fringe(object):
    def __init__(self, initial_elevation, pond, favour_neighbours=True):
        self.mean_elevation = initial_elevation
        self._lower_fringe = sortedlist(key=(lambda x: -1 * x.elevation))
        self._upper_fringe = sortedlist(key=(lambda x: x.elevation))
        self._fringe_pixels = set()
        self._bad_fringe = set()
        self._pond = pond
        self._favour_neighbour = favour_neighbours
        
    def __contains__(self, pixel):
        return pixel in self._fringe_pixels or pixel in self._bad_fringe

    def add_bad(self, pixel):
        self._bad_fringe.add(pixel)

    def add(self, pixel, elevation):
        self._fringe_pixels.add(pixel)
        if elevation > self.mean_elevation:
            self._upper_fringe.add(Pixel(pixel[0], pixel[1], elevation, 0, False))
        else:
            self._lower_fringe.add(Pixel(pixel[0], pixel[1], elevation, 0, False))

    def __iter__(self):
        return self

    def next(self):
        
        closest = self._find_closest_list_head(
                    [self._upper_fringe, self._lower_fringe],
                    (lambda x: x[0].elevation), 
                    self.mean_elevation)
        if closest is None:
            raise StopIteration
        else:
            self._fringe_pixels.remove((closest[0].x, closest[0].y))
            c = closest.pop()
            return (c.x, c.y), c.elevation

    def update_mean(self, new_mean):
        self.mean_elevation = new_mean
        repartition_lists(self._lower_fringe, self._upper_fringe, self.mean_elevation)

    def _find_closest_list_head(self, lists, val_func, compare_val):
        closest = None
        closest_diff = float("inf")
        for list in lists:
            try:
                diff = abs(val_func(list) - compare_val)
                if diff < closest_diff:
                    closest = list
                    closest_diff = diff
            except IndexError:
                continue
            return closest

class PondBuilder(object):
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
        try:
            self.parse_args()
        except:
            print "Could not open all files.  Aborting."
            sys.exit(1)

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
    # or it has already been assigned to a pond,
    # or it has no elevation data
    def can_skip_pixel(self, demblock, pairblock, outblock, block_x, block_y):
        # numpy arrays are accessed array[y][x]
        pair_val = pairblock.data[block_y][block_x]
        # skip empty pixel pairs
        if is_no_data(self.pxpairs_no_data, pair_val):
            return True
        dem_val = demblock.data[block_y][block_x]
        # skip empty pixel pairs
        if is_no_data(self.dem_no_data, dem_val):
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
        pixel_check_count = 0
        for (demblock, pxpairsblock, outblock), block_x, block_y, world_x, world_y in ( 
            self.io.extent_iterator(*self.world_extent)):
            if self.can_skip_pixel(demblock, pxpairsblock, outblock, block_x, block_y):
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
            fringe = Fringe(elevation, pond_pixels)
            # there is no do while, so we'll manually break out of the loop
        #    print "New pond, number %4d" % pond_num
            while True: 
                # add surrounding 8 pixels to the fringe if valid
                for f in neighbouring_pixels(new_pixel, 4, self.world_extent):
                    # screen out any pixels we don't want in the fringe
                    if (f in fringe or 
                        f in pond_pixels):
                        continue
                    f_demval, f_pairval, f_outval = self.io.get_pixel(f)
                    # if it doesn't have an elevation value, skip it
                    if is_no_data(self.dem_no_data, f_demval):
                        fringe.add_bad(f)
                        continue
                    # if pixel is already in a pond, skip it
                    if not is_no_data(POND_NO_DATA, f_outval):
                        fringe.add_bad(f)
                        continue
                    # add valid pixel to fringe
                    fringe.add(f, f_demval)
                    
                # if the fringes are empty, then we have finished this pond
                try:
                    new_pixel, new_elevation = fringe.next()
                except StopIteration:
                    break
                # if the closest pixel is too far away, then we are done
                # now recalculate pond statistics
                pond_size += 1
                pond_sum_elevation += new_elevation
                pond_mean_elevation = float(pond_sum_elevation) / pond_size
                pond_pixels[new_pixel] = new_elevation
                # calculate total pond earth to move
                old_mean_diff = mean_diff
                # do we want gross mean diff or net mean diff
                # net mean diff means that any cutting and filling can be combined
                # gross:
                total_diff = sum([abs(x - pond_mean_elevation) for x in pond_pixels.values()])
                #net
                # total_diff = abs(sum(map(lambda x: x - pond_mean_elevation, pond_pixels.values())))
                mean_diff = total_diff / pond_size
                # If the closest pond isn't close enough
                if mean_diff > self.max_height_diff:
                    del pond_pixels[new_pixel]
                    mean_diff = old_mean_diff
                    pond_size -= 1
                    pond_sum_elevation -= new_elevation
                    pond_mean_elevation = float(pond_sum_elevation) / pond_size
                    break
                fringe.update_mean(pond_mean_elevation)

            if pond_size >= self.min_pond_pixels:
#                print "Found pond %d, %d pixels" % (self.pond_num, pond_size)        
                #save the found pond to buffer
                for save_pixel, save_elevation in pond_pixels.items():
                    self.io.set_pixel(self.outband, save_pixel, self.pond_num)
                pond_m2 = pond_size * self.xcellsize * self.ycellsize
                self.pond_attributes[self.pond_num] = (pond_m2,
                                                  pond_mean_elevation,
                                                  total_diff * pond_m2)
                print "Found pond %d, %d pixels (%.2f m sq.), mean elevation %.2f, total_earth_to_move %.2f" % (self.pond_num, pond_size, pond_m2,
                                                  pond_mean_elevation,
                                                  total_diff * pond_m2)
                self.io.write_flush() 
                self.pond_num += 1
            pixel_check_count += 1
            period = int(self.total_pixel_count / 100)
            if self.total_pixel_count % pixel_check_count == period:
                print "Completed %.2f%%" % (float(self.total_pixel_count) / pixel_check_count * 100)
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
