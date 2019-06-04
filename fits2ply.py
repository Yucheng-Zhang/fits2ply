'''
Generate .ply polygon file from .fits.
'''
from astropy.io import fits
import numpy as np
import argparse
import collections

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Generate .ply polygon file from .fits.')
    parser.add_argument('-fits', type=str, default='',
                        help='Input fits file.')
    parser.add_argument('-wattr', type=str, default='WEIGHT',
                        help='Attribute that gives the sector completeness.')
    parser.add_argument('-fo', type=str, default='output.ply',
                        help='Output polygon file.')

    args = parser.parse_args()


def write_ply(fo, plys):
    '''Write polygons to file.'''
    fw = open(fo, 'w')
    nplys = plys['nplys']
    fw.write('{0:d} polygons\n'.format(nplys))

    for i in range(nplys):
        ncaps = plys[i]['ncaps']
        weight = plys[i]['weight']
        pixel = plys[i]['pixel']
        str_ = plys[i]['str']
        caps = plys[i]['caps']

        fw.write('polygon {0:10d} ( {1:d} caps, {2:.7f} weight, {3:d} pixel, {4:25.20f} str):\n'.format(
            i, ncaps, weight, pixel, str_))
        for j in range(ncaps):
            fw.write(
                '{0:25.20f} {1:25.20f} {2:25.20f} {3:25.20f}\n'.format(caps[j]))

    fw.close()
    print(':: Written: {}'.format(fo))


if __name__ == "__main__":

    hdu = fits.open(args.fits)
    data = hdu[1].data

    plys = collections.OrderedDict()
    plys['nplys'] = len(data)
    print('>> Number of polygons: {0:d}'.format(plys['nplys']))

    for i in range(plys['nplys']):
        ncaps = data['NCAPS'][i]
        caps = data['XCAPS'][i][:ncaps]
        caps = np.column_stack((caps, data['CMCAPS'][i][:ncaps]))
        weight = data[args.wattr]
        pixel = data['PIXEL']
        str_ = data['STR']

        plys[i] = collections.OrderedDict()
        plys[i]['ncaps'] = ncaps
        plys[i]['caps'] = caps
        plys[i]['weight'] = weight
        plys[i]['pixel'] = pixel
        plys[i]['str'] = str_

    write_ply(args.fo, plys)
