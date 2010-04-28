
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

from lrucache import LRUCache

class UnequalBandException(Exception):
    pass

class MultiBandBlockIO:
        
    _block_cache = None
    bands = None
    xsize, ysize, xblockcount, yblockcount, xblocksize, yblocksize = (None, None, None, None, None, None)
    verbose = False
    
    class BlockWrapper:
        def __init__(self, band, xoff, yoff, width, height):
            self.band = band
            self._xoff = xoff
            self._yoff = yoff
            self._width = width
            self._height = height
            self.dirty = False
#            if 
            self.data = band.ReadAsArray(xoff,
                                      yoff,
                                      width,
                                      height,
                                      width,
                                      height)
        def finalise(self):
            if self.dirty:
                self.band.WriteArray(self.data, self._xoff, self._yoff)
                self.dirty = False

    @staticmethod
    def _ceil_div(numerator, denominator):
        return (numerator + denominator - 1) / denominator

    def __init__(self, bands, cache_size=24000, verbose=False):
        self._block_cache = LRUCache(cache_size)
        self.bands = bands
        self.verbose = verbose
        bsize = None
        i = 0
        b = bands[i]
        self.xsize = b.XSize
        self.ysize = b.YSize
        (self.xblocksize, self.yblocksize) = b.GetBlockSize()
        self.xblockcount = self._ceil_div(self.xsize, self.xblocksize)
        self.yblockcount = self._ceil_div(self.ysize, self.yblocksize)
        # Ensure all bands are of equal size and block size
        while True:
            i += 1
            if i >= len(bands):
                break
            b = bands[i]
            if ((self.xsize != b.XSize) or
                (self.ysize != b.YSize) or
                ((self.xblocksize, self.yblocksize) != 
                 tuple(b.GetBlockSize()))):
               
               raise UnequalBandException(((self.xsize, self.ysize, 
                                            self.xblocksize, self.yblocksize), 
                                           (b.XSize, b.YSize) + 
                                           tuple(b.GetBlockSize())))

        if self.verbose:
            print ("Files are %d (%d x %d) pixels, block size %d x %d, "
                   "meaning %d x %d blocks.  Caching %d blocks." %
                   (self.xsize * self.ysize, 
                    self.xsize, self.ysize, 
                    self.xblocksize, self.yblocksize,
                    self._ceil_div(self.xsize, self.xblocksize), 
                    self._ceil_div(self.ysize, self.yblocksize),
                    cache_size))

    def search_extent(self, x, y, r):
        if isinstance(r, float):
            r = int(r + 0.5)
        x1 = max(x - r, 0)
        x2 = min(x + r, self.xsize)
        y1 = max(y - r, 0)
        y2 = min(y + r, self.ysize)
        return (x1, y1, x2, y2)

    def pixel_coord_to_block_coord(self, x, y):
        return (x / self.xblocksize, y / self.yblocksize)

    def block_to_pixel_coord(self, block_x, block_y):
        return (block_x * self.xblocksize, 
                block_y * self.yblocksize)

    def block_to_pixel_extent(self, block_x, block_y):
        (x1, y1) = self.block_to_pixel_coord(block_x, block_y)
        x2 = x1 + min(self.xblocksize, self.xsize - x1)
        y2 = y1 + min(self.yblocksize, self.ysize - y1)
        return (x1, y1, x2, y2)

    def find_block_extent(self, (x1, y1, x2, y2)):
        return (self.pixel_coord_to_block_coord(x1, y1) + 
                self.pixel_coord_to_block_coord(x2, y2))

    

    # adapted from http://www.python.org/download/releases/src/lib1.4.tar.gz
    # /Lib/stdwin/rect.py
    @staticmethod
    def intersection(list):
        def is_empty(r):
            (left, top, right, bottom) = r
            return left >= right or top >= bottom
        
        empty = (0, 0, 0, 0)
        if not list: raise error, 'intersect called with empty list'
        if is_empty(list[0]): return empty
        (left, top, right, bottom) = list[0]
        for rect in list[1:]:
                if is_empty(rect):
                        return empty
                (l, t, r, b) = rect
                if left < l: left = l
                if top < t: top = t
                if right > r: right = r
                if bottom > b: bottom = b
                if is_empty(((left, top, right, bottom))):
                        return empty
        return (left, top, right, bottom)

    def get_block_in_bands(self, block_x, block_y):
        res = []
        # Array size is not always self.xblocksize * self.yblocksize,
        # because we may be on edge of extent
        (xoff, yoff, x2, y2) = self.block_to_pixel_extent(block_x, block_y)
        width = x2 - xoff
        height = y2 - yoff
        for b in self.bands:
            key = (b, xoff, yoff, width, height)
            wrapped_array = self._block_cache.get(key)
            if wrapped_array is None:
#                print "Miss"
                wrapped_array = self.BlockWrapper(b, xoff, yoff, width, height)
                self._block_cache.put(key, wrapped_array)
            res.append(wrapped_array)
        return tuple(res)

    def write_flush(self):
        self._block_cache.flush()

    