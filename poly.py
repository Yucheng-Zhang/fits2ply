'''
Polygon, mangle, fits and Healpy map.
'''
import healpy as hp
import numpy as np
import collections
from astropy.io import fits


class poly:
    '''Polygon, mangle, fits and Healpy map.'''

    def __init__(self):
        # polygons
        self.plys = collections.OrderedDict()
        self.plys['nplys'] = None

    def read_fits(self, fn, wattr):
        '''Read polygons from file.'''
        hdu = fits.open(fn)
        data = hdu[1].data

        self.plys['nplys'] = len(data)
        print('>> Number of polygons: {0:d}'.format(self.plys['nplys']))

    def read_ply(self, fn):
        '''Read polygons from .ply file.'''
        pass

    def write_ply(self, fo):
        '''Write polygons to .ply file.'''
        pass

    def make_map(self, nside):
        '''Make Healpy map for mask.'''
        pass
