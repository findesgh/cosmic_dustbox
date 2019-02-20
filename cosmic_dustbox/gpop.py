import astropy.units as _u
import numpy as _numpy
import warnings as _warnings
import os as _os
from scipy.interpolate import interp1d
from scipy import integrate
import sys
Path = _os.path.dirname(_os.path.dirname(_os.path.realpath(__file__))) + \
                 '//Solver//'
sys.path.insert(0, Path)
import sdist as sd
import GrainSpecies as gs
###############################################################################
###############################################################################


class Dp(object):

    def __init__(self, en):
        self.en = en
        return

    def Sigma_abs_E(en):
        sizes = _np.logspace(self.sizeMin, self.sizeMax,
                             num=10000)

    def func1(sizes):
        if issubclass(self.__class__, sd.SizeDist):
            dist_gr = self.__class__(sizes).sdist_gr
            dist_sil = self.__class__(sizes).sdist_sil
        return dist_gr, dist_sil

    def func2(sizes):
        sigma_cr = gs.SpGraphite.Sigma_abs(en, sizes)[0]
        sigma_sl = gs.SpAsil.Sigma_abs(en, sizes)[0]
        return sigma_cr, sigma_sl

    f1_gr = interp1d(sizes, func1(sizes)[0], kind='linear')
    f2_sl = interp1d(sizes, func1(sizes)[1], kind='linear')

    f3_gr = interp1d(sizes, func2(sizes)[0], kind='linear')
    f4_sl = interp1d(sizes, func2(sizes)[1], kind='linear')

    Integ1 =integrate.quad(f1_gr(Sizes) * f3_gr(Sizes),self.sizeMin,
                           self.sizeMax)
    Integ2 = integrate.quad(f1_gr(Sizes), self.sizeMin,self.sizeMax)

    Integ3 = integrate.quad(f2_sl(Sizes) * f4_sl(Sizes),self.sizeMin,
                           self.sizeMax)
    Integ4 = integrate.quad(f2_sl(Sizes), self.sizeMin,self.sizeMax)

    sigma_abs_1 = Integ1/Integ2
    sigma_abs_2 = Integ3/Integ4
    sigma_abs_E = sigma_abs_1 + sigma_abs_2

    return sigma_abs_E
