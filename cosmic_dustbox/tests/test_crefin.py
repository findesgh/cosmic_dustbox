from unittest import TestCase
from .. import crefin
import numpy as np
import astropy.units as u


class TestFromData(TestCase):

    @classmethod
    def setUpClass(cls):
        cls.data = {
            'a': [1e-9],
            'lam': np.linspace(1, 10, num=10),
            'n': np.random.random(10) + 1j * np.random.random(10),
        }
        return

    @classmethod
    def tearDownClass(cls):
        return

    def test_smoke_test(self):
        a = crefin.Crefin.fromData(self.__class__.data)
        return

    def test_interp(self):
        a = crefin.Crefin.fromData(self.__class__.data)
        a(np.array([1.0])*u.micron, np.array([1.0])*u.micron)
        return
