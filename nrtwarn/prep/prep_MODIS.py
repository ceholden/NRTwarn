#!/usr/bin/env python
""" Preprocess MODIS MOD09GA/MYD09GA and MOD09GQ/MYD09GQ

Usage:
    prep_MODIS.py [options] <input_location> <output_location>

Options:
    --pattern=<pattern>     Pattern for files to preprocess [default: M*09*hdf]
    -n --ncpu=<n>           Number of CPUs to use [default: 1]
    --nodata=<ndv>          Output NODATA value [default: -9999]
    -v --verbose            Show verbose debugging options
    -q --quiet              Do not show extra information
    -h --help               Show help

"""
from __future__ import print_function, division

import fnmatch
import logging
import multiprocessing
import os
import sys

from docopt import docopt
import numpy as np
from osgeo import gdal
import scipy.ndimage

_version = '0.1.0'

gdal.AllRegister()
gdal.UseExceptions()

logging.basicConfig(format='%(asctime)s %(levelname)s: %(message)s',
                    level=logging.INFO,
                    datefmt='%H:%M:%S')
logger = logging.getLogger(__name__)


def find_MODIS_pairs(location, pattern='M*09*hdf'):
    """ Finds matching sets of M[OY]D09GA and M[OY]D09GQ within location

    Args:
      location (str): directory of stored data
      pattern (str, optional): glob pattern to limit search

    Returns:
      pairs (list): list of tuples containing M[OY]D09GQ and M[OY]D09GA

    """
    files = [os.path.join(location, f) for f in
             fnmatch.filter(os.listdir(location), pattern)]

    if len(files) < 2:
        raise IOError('Could not find any MODIS image pairs')

    # Parse out product and acquisition date
    products = []
    dates = []
    for f in files:
        s = os.path.basename(f).split('.')
        products.append(s[0])
        dates.append(s[1])

    products = np.array(products)
    dates = np.array(dates)

    # Retain dates if there are matching MOD09GA/MOD09GQ or MYD09GA/MYD09GQ
    pairs = []
    for d in np.unique(dates):
        i = np.where(dates == d)[0]
        prods = products[i]

        i_mod09ga = np.core.defchararray.startswith(prods, 'MOD09GA')
        i_mod09gq = np.core.defchararray.startswith(prods, 'MOD09GQ')
        i_myd09ga = np.core.defchararray.startswith(prods, 'MYD09GA')
        i_myd09gq = np.core.defchararray.startswith(prods, 'MYD09GQ')

        if i_mod09gq.sum() == 1 and i_mod09ga.sum() == 1:
            pairs.append((files[i[i_mod09gq]], files[i[i_mod09ga]]))
        if i_myd09gq.sum() == 1 and i_myd09ga.sum() == 1:
            pairs.append((files[i[i_myd09gq]], files[i[i_myd09ga]]))

    logger.debug('Found {n} pairs of M[OY]D09GQ and M[OY]D09GA'.format(
        n=len(pairs)))

    return pairs


def get_mask(modQA, dilate=7):
    """ Return a mask image from an input QA band from M[OY]D09G[AQ]

    Args:
      modQA (ndarray): input QA/QC band image
      dilate (int, optional): pixels around aerosols and clouds to buffer

    Returns:
      mask (ndarray): output mask image with only good observations masked

    Porting note:
        MATLAB's 'bitshift' shifts to the left

    """
    # Identify land from water
    land = (np.mod(np.right_shift(modQA, 3) + 6, 8) / 7).astype(np.uint8)
    # Identify cloud
    cloud = (np.mod(modQA, 8) |  # unsure!
             np.mod(np.right_shift(modQA, 8), 4) |  # cirrus == '00' (none)
             np.mod(np.right_shift(modQA, 10), 2) |  # cloud mask == '0'
             np.mod(np.right_shift(modQA, 13), 2)) > 0  # adjacent to cloud

    cloud_buffer = scipy.ndimage.morphology.binary_dilation(
        cloud, structure=np.ones((dilate, dilate)))

    return ((cloud_buffer == 0) * land)


def enlarge(array, scaling):
    """ Enlarge an array by a scaling factor

    Args:
      array (ndarray): array to be scaled
      scaling (int): amount of scaling

    Returns:
      scaled (ndarray): scaled array

    """
    return np.kron(array, np.ones((scaling, scaling)))


if __name__ == '__main__':
    args = docopt(__doc__, version=_version)

    logger.setLevel(logging.INFO)
    if args['--quiet']:
        logger.setLevel(logging.WARNING)
    if args['--verbose']:
        logger.setLevel(logging.DEBUG)

    location = args['<input_location>']

    pattern = args['--pattern']
    ncpu = int(args['--ncpu'])
    ndv = int(args['--nodata'])

    # TEST
    logger.warning('JUST TESTING ON EXAMPLE')

    here = os.path.dirname(os.path.abspath(__file__))
    hdf = '../../tests/MOD09G/MOD09GA.A2000055.h10v08.005.2006268014216.hdf'
    ds = gdal.Open(os.path.join(here, hdf), gdal.GA_ReadOnly)

    sds_qa = ds.GetSubDatasets()[1]
    ds = None

    ds = gdal.Open(sds_qa[0], gdal.GA_ReadOnly)
    qa_band = ds.GetRasterBand(1).ReadAsArray()

    get_mask(qa_band)
