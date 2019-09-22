'''
Polygon, mangle, fits and Healpy map.
'''
import healpy as hp
import numpy as np
import collections
from astropy.io import fits
import time
from . import utils
import re


class poly:
    '''Polygon, mangle, fits and Healpy map.'''

    def __init__(self):
        # polygons
        self.plys = collections.OrderedDict()
        self.plys['nplys'] = None
        # self.plys[i] info for polygon i

    def read_fits(self, fn, wattr):
        '''Read polygons from file.'''
        hdu = fits.open(fn)
        data = hdu[1].data

        self.plys['nplys'] = len(data)
        print('>> Number of polygons: {0:d}'.format(self.plys['nplys']))

        for i in range(self.plys['nplys']):
            self.plys[i] = collections.OrderedDict()
            self.plys[i]['ncaps'] = ncaps_ = data['NCAPS'][i]
            # get directions and sizes of the caps
            self.plys[i]['caps'] = np.column_stack((data['XCAPS'][i][:ncaps_],
                                                    data['CMCAPS'][i][:ncaps_]))
            self.plys[i]['weight'] = data[wattr][i]
            self.plys[i]['pixel'] = data['PIXEL'][i]
            self.plys[i]['str'] = data['STR'][i]

    def read_ply(self, fn):
        '''Read polygons from .ply file.'''
        f = open(fn, 'r')
        line = f.readline()
        s = re.match(r'(\d+)\s+polygons', line)
        self.plys['nplys'] = int(s.group(1))

        # identify polygon
        ex1 = re.compile(r'polygon\s+(\d+)\s+\(\s*(\d+)\s+caps')
        # check weight
        ex2 = re.compile(r'(\d*\.?\d+)\s+weight')
        # check area
        ex3 = re.compile(r"(\d*\.?\d+)\s+str")

        s = ex1.match(line)
        # working
        f.close()
        pass

    def write_ply(self, fo):
        '''Write polygons to .ply file.'''
        fw = open(fo, 'w')
        fw.write('{0:d} polygons\n'.format(self.plys['nplys']))

        for i in range(self.plys['nplys']):
            fw.write('polygon {0:10d} ( {1:d} caps, {2:.7f} weight, {3:d} pixel, {4:25.20f} str):\n'.format(
                i, self.plys[i]['ncaps'], self.plys[i]['weight'], self.plys[i]['pixel'], self.plys[i]['str']))
            caps = self.plys[i]['caps']
            for j in range(self.plys[i]['ncaps']):
                fw.write('{0:25.20f} {1:25.20f} {2:25.20f} {3:25.20f}\n'.format(
                    caps[j][0], caps[j][1], caps[j][2], caps[j][3]))

        fw.close()
        print(':: Written: {0:s}'.format(fo))

    def _in_cap(self, cap, vec):
        '''Check if vectors in vec array are in cap.'''
        cd = 1. - np.sum(cap[:3] * vec, axis=1)
        if cap[3] < 0.:
            return cd > np.abs(cap[3])
        else:
            return cd < cap[3]

    def _in_ply(self, caps, vec):
        '''Check if vectors in vec array are in polygon.'''
        idx = np.arange(vec.shape[0], dtype=np.int32)
        for cap in caps:
            idx = idx[self._in_cap(cap, vec[idx])]
            if idx.shape[0] == 0:
                break

        res = np.full(vec.shape[0], False, dtype=np.bool)
        res[idx] = True

        return res

    def make_map(self, nside, fo=None, unseen=hp.UNSEEN, ralim=None, declim=None):
        '''Make Healpy map for mask.'''
        npix = hp.nside2npix(nside)
        ipix = np.arange(npix)
        mask = np.full(npix, unseen)
        theta_g, phi_g = hp.pix2ang(nside, ipix)  # theta, phi in G
        ra, dec = utils.get_ra_dec(theta_g, phi_g)

        # set RA (0,360) and DEC (-90,90) range to speed up
        if ralim != None or declim != None:
            # cut pixels
            if ralim != None:
                idx1 = (ra >= ralim[0]) & (ra <= ralim[1])
            if declim != None:
                idx2 = (dec >= declim[0]) & (dec <= declim[1])
            idx = idx1 & idx2
            ipix, ra, dec = ipix[idx], ra[idx], dec[idx]

        # theta, phi in C
        theta, phi = np.deg2rad(90.-dec), np.deg2rad(ra)

        # x, y, z coords of pixels on unit sphere
        vec = np.column_stack((np.sin(theta) * np.cos(phi),
                               np.sin(theta) * np.sin(phi),
                               np.cos(theta)))

        # loop over polygons
        print('>> looping over polygons...')
        t0 = time.time()
        for i in range(self.plys['nplys']):
            print(i)
            idx = self._in_ply(self.plys[i]['caps'], vec)
            # set mask value to polygon weight
            mask[ipix[idx]] = self.plys[i]['weight']
            # cut off the pixels in this polygon to speed up
            ipix = ipix[np.invert(idx)]
            vec = vec[np.invert(idx)]
        print('<< time elapsed: {0:.2f}'.format(time.time()-t0))

        if fo is not None:
            hp.write_map(fo, mask, overwrite=True)
            print(':: Written: {0:s}'.format(fo))

        return mask
