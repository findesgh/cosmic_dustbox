from scipy import interpolate as _interp
import astropy.units as _u
import numpy as _np
import os as _os


class Crefin(object):

    def __init__(self, f):
        self.f = f
        return

    def __call__(self, a, wave):
        return self.f(
            a.to(_u.micron, equivalencies=_u.spectral()).value,
            wave.to(_u.micron, equivalencies=_u.spectral()).value)

    @classmethod
    def fromData(cls, a, lam, n, bounds_error=True, a_bounds_error=False):
        """
        Interpolate complex refractive index from given data.

        Parameters
        ----------
        a, lam : array
            1-D arrays of grain sizes and wavelengths in micron. Must have
            length >= 2. `a` can also be None if there is no size dependence or
            it is not known

        n : array
            Array of the complex refractive index. Note: if `n` is 2-D, its
            first index is wavelength and its second index is grain size. For
            more details see the documentation of scipy.interpolate.interp2d.

        bounds_error : bool, optional
            Whether to raise out of bounds errors when interpolating.

        a_bounds_error : bool, optional
            If `bounds_error` is True, this can be used to specify if bounds
            errors on grain size should be ignored by setting it False. This is
            useful because often `n` is not known at many grain sizes. If
            `bounds_error` is False or if `a` is None, this is ignored.

        Note
        ----
        scipy.interpolate.interp2d implicitly sorts its `x` and `y` arguments
        and returns an array that's sorted accordingly. Here we shuffle the
        rows and columns of the returned array around so that their
        row-(column-)order corresponds to the order of the input arrays. This
        is highly inefficient but convenient. Like this the interpolation
        behaves exactly as a 2D function evaluated on a meshgrid of `x` and `y`
        would.
        """
        if a is None:
            interpolation = _interp.interp1d(lam, n, bounds_error=bounds_error)

            def f(a, lam):
                return _np.broadcast_to(interpolation(lam), (len(a), len(lam)))
        else:
            rinterpolation = _interp.interp2d(a, lam, n.real,
                bounds_error=bounds_error)
            iinterpolation = _interp.interp2d(a, lam, n.imag,
                bounds_error=bounds_error)

            if bounds_error and not a_bounds_error:
                def f(aa, lam):
                    r = _np.clip(aa, _np.amin(a), _np.amax(a))
                    xind = _np.argsort(r)
                    yind = _np.argsort(lam)
                    return (
                        rinterpolation(r, lam)
                        + 1j*iinterpolation(r, lam))[yind][:, xind].T
            else:
                def f(a, lam):
                    xind = _np.argsort(r)
                    yind = _np.argsort(lam)
                    return (
                        rinterpolation(a, lam)
                        + 1j*iinterpolation(a, lam))[yind][:, xind].T
        return cls(f)


class SGPAHCrefin(Crefin):

    @classmethod
    def parseCrefinFile(cls, path):
        with open(path) as f:
            data = _np.loadtxt(f, skiprows=5)
            f.seek(0)
            next(f)
            size = float(f.readline().split('=')[0])
        return size, data

    @classmethod
    def fromFiles(cls, paths):
        a = []
        n = []
        for j, path in enumerate(paths):
            size, data = cls.parseCrefinFile(path)
            if size in a:
                raise ValueError("Grain size "
                                 + str(size) + "provided twice.")
            else:
                a.append(size)
            if j == 0:
                lam = data[:, 0]
            else:
                if not _np.array_equal(lam, data[:, 0]):
                    raise ValueError(
                        "Wavelengths in file " + path
                        + "differ from those in the other files visited "
                        "so far.")
            nt = data[:, 3] + 1 + 1j*data[:, 4]
            n.append(nt)
        if len(a) > 1:
            a = _np.array(a)
            n = _np.array(n).T
        elif len(a) == 1:
            a = None
            n = _np.array(n)
        return cls.fromData(a, lam, n, bounds_error=True, a_bounds_error=False)


cpaD03 = SGPAHCrefin.fromFiles([
    _os.path.join(_os.path.dirname(__file__),
                  'data/crefin/callindex.out_CpaD03_0.01'),
    _os.path.join(_os.path.dirname(__file__),
                  'data/crefin/callindex.out_CpaD03_0.10')
])
cpeD03 = SGPAHCrefin.fromFiles([
    _os.path.join(_os.path.dirname(__file__),
                  'data/crefin/callindex.out_CpeD03_0.01'),
    _os.path.join(_os.path.dirname(__file__),
                  'data/crefin/callindex.out_CpeD03_0.10')
])
silD03 = SGPAHCrefin.fromFiles([
    _os.path.join(_os.path.dirname(__file__),
                  'data/crefin/callindex.out_silD03')
])
