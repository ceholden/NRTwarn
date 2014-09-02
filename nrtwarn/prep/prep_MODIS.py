#!/usr/bin/env python
""" Preprocess MODIS MOD09GA/MYD09GA and MOD09GQ/MYD09GQ

Usage:
    prep_MODIS.py [options] <location>

Options:
    --pattern=<pattern>     Glob pattern for files to preprocess [default: *]
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

    location = args['<location>']

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

    enlarge(get_mask(qa_band), 4).shape
