from scipy import interpolate as _interp
import astropy.units as _u
import numpy as _np


class Crefin(object):

    def __init__(self, f):
        self.f = f
        return

    def __call__(self, a, wave):
        return self.f(
            a.to(_u.micron, equivalencies=_u.spectral()),
            wave.to(_u.micron, equivalencies=_u.spectral()))

    @classmethod
    def fromData(cls, data, a_bounds_error=False, lam_bounds_error=True):
        if len(data['a']) == 1:
            interpolation = _interp.interp1d(
                data['lam'], data['n'], bounds_error=lam_bounds_error)

            def f(a, lam):
                n = interpolation(lam)
                return _np.broadcast_to(n, (len(a), len(lam)))
        elif len(data['a']) == 2:
            interpolation = _interp.interp2d(
                data['a'], data['lam'], data['n'], bounds_error=True)

            def f(a, lam):
                r = _np.clip(a, data['a'][0], data['a'][1])
                return interpolation(a, lam)
        else:
            raise ValueError("Currently no method implemented!")
        return cls(f)
