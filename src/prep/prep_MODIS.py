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

_version = '0.1.0'

gdal.AllRegister()
gdal.UseExceptions()

logging.basicConfig(format='%(asctime)s %(levelname)s: %(message)s',
                    level=logging.INFO,
                    datefmt='%H:%M:%S')
logger = logging.getLogger(__name__)


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
