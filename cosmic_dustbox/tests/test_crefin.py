from unittest import TestCase
from .. import crefin
import numpy as np
import astropy.units as u
import os


class TestFromData(TestCase):

    @classmethod
    def setUpClass(cls):
        cls.a1 = None
        cls.lam1 = np.linspace(1, 10, num=10)
        cls.n1= np.random.random(10) + 1j * np.random.random(10)

        cls.a2 = np.array([1e-9, 1e-8])
        cls.lam2 = np.linspace(1, 10, num=10)
        cls.n2 = (
            np.random.random(20) + 1j * np.random.random(20)).reshape(10, 2)
        return

    @classmethod
    def tearDownClass(cls):
        return

    def test_no_a(self):
        c = self.__class__
        a = crefin.Crefin.fromData(c.a1, c.lam1, c.n1)
        return

    def test_normal(self):
        c = self.__class__
        a = crefin.Crefin.fromData(c.a2, c.lam2, c.n2)
        return

    def test_normal_no_bounds(self):
        c = self.__class__
        a = crefin.Crefin.fromData(c.a2, c.lam2, c.n2, bounds_error=False)
        return


class TestSGPAHCrefin(TestCase):

    @classmethod
    def setUpClass(cls):
        cls.basepath = os.path.join(
            os.path.dirname(os.path.abspath(__file__)), '../data/crefin')
        return

    @classmethod
    def tearDownClass(cls):
        return

    def test_parseCrefinFile_asil(self):
        path = os.path.join(self.__class__.basepath, 'callindex.out_silD03')
        size, data = crefin.SGPAHCrefin.parseCrefinFile(path)
        self.assertAlmostEqual(0.1, size)
        self.assertEqual(data.shape[0], 837)
        self.assertEqual(data.shape[1], 5)
        return

    def test_parseCrefinFile_cpar(self):
        path = os.path.join(self.__class__.basepath,
                            'callindex.out_CpaD03_0.01')
        size, data = crefin.SGPAHCrefin.parseCrefinFile(path)
        self.assertAlmostEqual(0.01, size)
        self.assertEqual(data.shape[0], 386)
        self.assertEqual(data.shape[1], 5)
        return

    def test_parseCrefinFile_cperp(self):
        path = os.path.join(self.__class__.basepath,
                            'callindex.out_CpeD03_0.01')
        size, data = crefin.SGPAHCrefin.parseCrefinFile(path)
        self.assertAlmostEqual(0.01, size)
        self.assertEqual(data.shape[0], 383)
        self.assertEqual(data.shape[1], 5)
        return

    def test_fromFiles_asil(self):
        path = os.path.join(self.__class__.basepath, 'callindex.out_silD03')
        c = crefin.SGPAHCrefin.fromFiles([path])
        a = np.array([0.1])*u.micron
        wave = np.array([1.239840e5, 5.481167e2, 6.199200e-5])*u.micron
        nre_expect = [2.4350+1, 2.4040+1, -1.9450e-6+1]
        nim_expect = [1.1190e-3, 8.2550e-2, 1.7830e-8]
        n = c(a, wave)

        self.assertEqual(n.shape[0], len(a))
        self.assertEqual(n.shape[1], len(wave))
        for j, r in enumerate(n[0, :]):
            self.assertAlmostEqual(r.real, nre_expect[j])
            self.assertAlmostEqual(r.imag, nim_expect[j])
        return

    def test_fromFiles_asil_multiple_a(self):
        path = os.path.join(self.__class__.basepath, 'callindex.out_silD03')
        c = crefin.SGPAHCrefin.fromFiles([path])
        a = np.array([1.0, 0.1, 0.01])*u.micron
        wave = np.array([1.239840e5, 5.481167e2, 6.199200e-5])*u.micron
        nre_expect = [2.4350+1, 2.4040+1, -1.9450e-6+1]
        nim_expect = [1.1190e-3, 8.2550e-2, 1.7830e-8]
        n = c(a, wave)

        self.assertEqual(n.shape[0], len(a))
        self.assertEqual(n.shape[1], len(wave))
        for i in range(len(a)):
            for j, r in enumerate(n[i, :]):
                self.assertAlmostEqual(r.real, nre_expect[j])
                self.assertAlmostEqual(r.imag, nim_expect[j])
        return

    def test_fromFiles_cperp(self):
        paths = [
            os.path.join(self.__class__.basepath, 'callindex.out_CpeD03_0.01'),
            os.path.join(self.__class__.basepath, 'callindex.out_CpeD03_0.10')
        ]
        c = crefin.SGPAHCrefin.fromFiles(paths)
        a = np.array([1.0, 0.1, 0.01, 0.001])*u.micron
        wave = np.array([1.239840E+05, 6.322489E+00, 6.199200E-05])*u.micron
        nre_expect = np.array([
            [1.0596E+03, 4.5430E+00, 1.5405E-10],
            [1.0596E+03, 4.5430E+00, 1.5405E-10],
            [3.4137E+02, 4.5568E+00, 1.5405E-10],
            [3.4137E+02, 4.5568E+00, 1.5405E-10],
        ])+1
        nim_expect = [
            [1.0636E+03, 4.2997E+00, 1.4471E-17],
            [1.0636E+03, 4.2997E+00, 1.4471E-17],
            [3.4149E+02, 4.3122E+00, 1.3963E-16],
            [3.4149E+02, 4.3122E+00, 1.3963E-16],
        ]
        n = c(a, wave)

        self.assertEqual(n.shape[0], len(a))
        self.assertEqual(n.shape[1], len(wave))
        for i in [0, 1]:
            for j, r in enumerate(n[i, :]):
                self.assertAlmostEqual(r.real, nre_expect[i][j],
                                       msg='Failed at i='+str(i)+' j='+str(j))
                self.assertAlmostEqual(r.imag, nim_expect[i][j])
        return

    def test_fromFiles_cpar(self):
        paths = [
            os.path.join(self.__class__.basepath, 'callindex.out_CpaD03_0.01'),
            os.path.join(self.__class__.basepath, 'callindex.out_CpaD03_0.10')
        ]
        c = crefin.SGPAHCrefin.fromFiles(paths)
        a = np.array([1.0, 0.1, 0.01, 0.001])*u.micron
        wave = np.array([1.239840E+05, 6.199200E+00, 6.199200E-05])*u.micron
        nre_expect = np.array([
            [1.5051E+02, 1.1266E+00, 2.3165E-10],
            [1.5051E+02, 1.1266E+00, 2.3165E-10],
            [1.4329E+02, 1.0261E+00, 2.3165E-10],
            [1.4329E+02, 1.0261E+00, 2.3165E-10],
        ])+1
        nim_expect = [
            [1.5156E+02, 2.3038E-02, 1.4161E-17],
            [1.5156E+02, 2.3038E-02, 1.4161E-17],
            [1.4434E+02, 3.0850E+00, 1.5614E-17],
            [1.4434E+02, 3.0850E+00, 1.5614E-17],
        ]
        n = c(a, wave)

        self.assertEqual(n.shape[0], len(a))
        self.assertEqual(n.shape[1], len(wave))
        for i in [0, 1]:
            for j, r in enumerate(n[i, :]):
                self.assertAlmostEqual(r.real, nre_expect[i][j],
                                       msg='Failed at i='+str(i)+' j='+str(j))
                self.assertAlmostEqual(r.imag, nim_expect[i][j])
        return

    def test_fromFiles_cpar(self):
        paths = [
            os.path.join(self.__class__.basepath, 'callindex.out_CpaD03_0.01'),
            os.path.join(self.__class__.basepath, 'callindex.out_CpaD03_0.10')
        ]
        a = crefin.SGPAHCrefin.fromFiles(paths)
        return

    def test_fromFiles_diff_wavelengths(self):
        paths = [
            os.path.join(self.__class__.basepath, 'callindex.out_CpeD03_0.01'),
            os.path.join(self.__class__.basepath, 'callindex.out_CpaD03_0.10')
        ]
        with self.assertRaises(ValueError) as context:
            a = crefin.SGPAHCrefin.fromFiles(paths)
        return

    def test_fromFiles_same_size(self):
        paths = [
            os.path.join(self.__class__.basepath, 'callindex.out_CpeD03_0.01'),
            os.path.join(self.__class__.basepath, 'callindex.out_CpeD03_0.01')
        ]
        with self.assertRaises(ValueError) as context:
            a = crefin.SGPAHCrefin.fromFiles(paths)
        return
