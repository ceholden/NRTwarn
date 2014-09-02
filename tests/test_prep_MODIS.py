""" Test pre-processing workflow
"""
import os
import sys
import unittest

import numpy as np
from osgeo import gdal

gdal.UseExceptions()
gdal.AllRegister()

here = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(here, '../'))

from nrtwarn.prep import prep_MODIS


class PreprocessTest(unittest.TestCase):
    """ Tests for `prep_MODIS.py`

    Tests are run against pre-existing MATLAB code results
    """

    def setUp(self):
        """ Load in test data """
        mod09ga = 'MOD09G/MOD09GA.A2000055.h10v08.005.2006268014216.hdf'
        ds = gdal.Open(os.path.join(here, mod09ga), gdal.GA_ReadOnly)
        state_sds = ds.GetSubDatasets()[1]
        ds = None
        self.state_ds = gdal.Open(state_sds[0], gdal.GA_ReadOnly)
        self.state = self.state_ds.GetRasterBand(1).ReadAsArray()
        self.state_ds = None

    def test_setUp(self):
        """ Test full sum of QA band """
        total = 14906462501
        self.assertEqual(self.state.sum(), total)

    def test_mask_sum(self):
        """ Test total number of good pixels in mask result """
        total = 1048607
        mask = prep_MODIS.get_mask(self.state, dilate=7)
        self.assertEqual(mask.sum(), total)

    def test_mask_locations(self):
        """ Test assortment of pixel values in mask """
        px = [49, 349, 999, 998, 749, 1199]
        py = [0, 98, 199, 399, 998, 0]

        values = [0, 0, 1, 1, 1, 1]

        mask = prep_MODIS.get_mask(self.state, dilate=7)

        for _px, _py, _v in zip(px, py, values):
            self.assertEqual(
                mask[_py, _px], _v,
                'Pixel x/y {px}/{py} is wrong ({m} vs. {t})'.format(
                    px=_px, py=_py, m=mask[_py, _px], t=_v))

    def test_enlarge(self):
        """ Test image rescaling """
        output = np.array([[1, 1, 0, 0],
                           [1, 1, 0, 0],
                           [0, 0, 1, 1],
                           [0, 0, 1, 1]])

        scaled = prep_MODIS.enlarge(np.eye(2), 2)

        np.testing.assert_array_equal(output, scaled)

    def test_find_MODIS_pairs(self):
        """ Test finding of MODIS image pairs from example dataset """
        test_location = os.path.join(here, 'MOD09G')

        # Test against included dataset
        images = [('MOD09GQ.A2000055.h10v08.005.2006268014216.hdf',
                  'MOD09GA.A2000055.h10v08.005.2006268014216.hdf')]

        found_images = prep_MODIS.find_MODIS_pairs(test_location, 'MOD09G*')
        found_images = [(os.path.basename(ga), os.path.basename(gq))
                        for ga, gq in found_images]

        self.assertListEqual(images, found_images)

    def test_find_MODIS_pairs_bad_pattern(self):
        """ Test finding of MODIS image pairs with bad pattern """
        test_location = os.path.join(here, 'MOD09G')

        # Test against included dataset with bad pattern
        with self.assertRaises(IOError):
            prep_MODIS.find_MODIS_pairs(test_location, 'MYD09G*')




if __name__ == '__main__':
    unittest.main()
