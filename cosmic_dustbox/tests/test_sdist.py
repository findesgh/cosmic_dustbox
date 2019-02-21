from unittest import TestCase
from .. import sdist
import numpy as np
import astropy.units as u


class TestSdist(TestCase):

    @classmethod
    def setUpClass(cls):
        cls._sizes = np.logspace(-11, -5)*u.m
        cls._unitLessSizes = cls._sizes.value
        return

    @classmethod
    def tearDownClass(cls):
        return

    def test_init(self):
        """
        Make sure we can create an instance.
        """
        sdist.SizeDist(3.5*u.angstrom, 1*u.micron, lambda x: 1/x.unit)

    def test_call(self):
        low = 4.0*u.angstrom
        high = 1*u.micron
        a = sdist.SizeDist(low, high, lambda x: 1/x.unit)
        r = a(self.__class__._sizes)
        self.assertTrue(
            np.all(r[self.__class__._sizes > high] == 0.0))
        self.assertTrue(
            np.all(r[self.__class__._sizes < low] == 0.0))
        self.assertTrue(
            np.all(r[np.where(r > 0)] == 1.0/u.m))
        return

    def test_unitless_sizes(self):
        a = sdist.SizeDist(3.5*u.angstrom, 1*u.micron, lambda x: 1/x.unit)
        with self.assertRaises(AttributeError) as context:
            a(self.__class__._unitLessSizes)
        return

    def test_unitless_func(self):
        a = sdist.SizeDist(3.5*u.angstrom, 1*u.micron, lambda x: 1)
        with self.assertRaises(u.core.UnitConversionError) as context:
            a(self.__class__._sizes)
        return

    def test_scalar_sizes(self):
        """
        This is something to be fixed in the future.
        """
        a = sdist.SizeDist(3.5*u.angstrom, 1*u.micron, lambda x: 1/x.unit)
        with self.assertRaises(TypeError) as context:
            a(0.5*u.micron)
        return

    def test_add_scalar_quantity_int(self):
        scalar = 1/u.micron
        a = sdist.SizeDist(3.5*u.angstrom, 1*u.micron, lambda x: 1/x.unit)
        b = a + scalar
        r = a(self.__class__._sizes)
        r[np.where(r != 0)] = r[np.where(r != 0)] + scalar
        self.assertTrue(
            np.all(
                np.isclose(b(self.__class__._sizes), r, rtol=1e-10, atol=0)
            ))
        return

    def test_add_scalar_quantity_float(self):
        scalar = 2.5/u.micron
        a = sdist.SizeDist(3.5*u.angstrom, 1*u.micron, lambda x: 1/x.unit)
        b = a + scalar
        r = a(self.__class__._sizes)
        r[np.where(r != 0)] = r[np.where(r != 0)] + scalar
        self.assertTrue(
            np.all(
                np.isclose(b(self.__class__._sizes), r, rtol=1e-10, atol=0)
            ))
        return


class TestWD01(TestCase):

    def test_smoke_test(self):
        sdist.WD01(3.1, 0.0, 'A')
        return

    def test_wrong_params(self):
        with self.assertRaises(ValueError) as context:
            sdist.WD01(3.1, 8.0, 'A')
        with self.assertRaises(ValueError) as context:
            sdist.WD01(3.1, 8.0, 'C')
        return

    def test_out_of_bounds(self):
        for sd in sdist.WD01(3.1, 0.0, 'A'):
            self.assertEqual(
                sd(np.array([1])*u.angstrom)[0],
                0.0/u.angstrom
            )
            self.assertEqual(
                sd(np.array([11])*u.micron)[0],
                0.0/u.angstrom
            )
        for sd in sdist.WD01(3.1, 6.0, 'A'):
            self.assertEqual(
                sd(np.array([1])*u.angstrom)[0],
                0.0/u.angstrom
            )
            self.assertEqual(
                sd(np.array([11])*u.micron)[0],
                0.0/u.angstrom
            )
