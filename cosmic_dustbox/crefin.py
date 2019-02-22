from scipy import interpolate as _interp
import astropy.units as _u
import numpy as _np
import os as _os


class Crefin(object):
    """
    Base class used to represent the complex refractive index of a certain
    grain material possibly along a certain axis.

    This class implements a `__call__` method which computes the crefin on a
    grid of grain size and photon energy (or equivalent).

    Parameters
    ----------
    f : callable
        A callable object that takes as first argument a 1-D array of grain
        sizes in micron and as second argument a 1-D array of wavelengths in
        micron. It shall return a 2-D array of the complex refractive index
        evaluated at the given sizes and wavelengths. The first axis of this
        array shall correspond to grain size and the second axis to
        wavelength. Note that this axis ordering is reversed wrt to the output
        of `np.meshgrid`.

    Attributes
    ----------
    _f : callable
        Directly taken from parameters.
    """

    def __init__(self, f):
        """
        See class docstring.
        """
        self._f = f
        return

    def __call__(self, a, wave):
        """
        Compute the complex refractive index as a function of grain size and
        wavelength or equivalent.

        Parameters
        ----------
        a : array-valued Quantity [L]
            1-D array of grain sizes as which to evaluate crefin.

        wave : array-valued Quantity
            1-D array of photon wavelengths or equivalent at which to evaluate
            crefin.

        Returns
        -------
         : array
            2-D array of crefin values. The first axis corresponds to grain
            size and the second to photon wavelength.
        """
        return self._f(
            a.to(_u.micron, equivalencies=_u.spectral()).value,
            wave.to(_u.micron, equivalencies=_u.spectral()).value)

    @classmethod
    def fromData(cls, a, lam, n, bounds_error=True, a_bounds_error=False):
        """
        Linearly interpolate complex refractive index from given data.

        Create an object which interpolates the complex refractive index in
        grain size and wavelength from the given data in its `__call__` method.

        Parameters
        ----------
        a, lam : array
            1-D arrays of grain sizes and wavelengths in micron. Must have
            length >= 2. `a` can also be None if there is no size dependence or
            it is not known.

        n : array
            2-D array of the complex refractive index. Note: its first index is
            grain size and its second index is wavelength, so if len(a) == l
            and len(lam) == k, n.shape == (l, k). If a is None, `n` shall be
            1-D with the same length as `lam`.

        bounds_error : bool, optional
            Whether to raise out of bounds errors when interpolating. If False,
            extrapolation is used. See also `a_bounds_error`.

        a_bounds_error : bool, optional
            If `bounds_error` is True, this can be used to specify if bounds
            errors on grain size should be ignored by setting it False. This is
            useful because often `n` is not known at many grain sizes. Note
            that if the bounds in grain size are to be ignored, the values of
            `n` for the largest/smallest grain size will be used outside the
            bounds instead of extrapolating. If `bounds_error` is False,
            `a_bounds_error` will be ignored, i.e. in this case extrapolation
            is used as for the wavelengths. If `a` is None, it is assumed that
            there is no size dependence of `n`.

        Returns
        -------
         : cls
            An instance of the calling class the `__call__` method of which
            interpolates the crefin on the given data.
        """
        if a is None:
            interpolation = _interp.interp1d(lam, n, bounds_error=bounds_error)

            def f(a, lam):
                return _np.broadcast_to(interpolation(lam), (len(a), len(lam)))
        else:
            rinterpolation = _interp.interp2d(lam, a, n.real,
                bounds_error=bounds_error)
            iinterpolation = _interp.interp2d(lam, a, n.imag,
                bounds_error=bounds_error)

            # scipy.interpolate.interp2d implicitly sorts its `x` and `y`
            # arguments and returns an array that's sorted accordingly. Here we
            # shuffle the rows and columns of the returned array around so that
            # their row-(column-)order corresponds to the order of the input
            # arrays. This is done for convenience, however it is inefficient
            # and might be improved in the future.
            if bounds_error and not a_bounds_error:
                def f(aa, lam):
                    r = _np.clip(aa, _np.amin(a), _np.amax(a))
                    n = rinterpolation(lam, r) + 1j*iinterpolation(lam, r)
                    if _np.prod(r.shape)*_np.prod(lam.shape) > 1:
                        rowind_unsorted = _np.argsort(r)
                        colind_unsorted = _np.argsort(lam)
                        return n[rowind_unsorted][:, colind_unsorted]
                    else:
                        return _np.atleast_2d(n)
            else:
                def f(a, lam):
                    n = rinterpolation(lam, a) + 1j*iinterpolation(lam, a)
                    if _np.prod(r.shape)*_np.prod(lam.shape) > 1:
                        rowind_unsorted = _np.argsort(lam)
                        colind_unsorted = _np.argsort(r)
                        return n[rowind_unsorted][:, colind_unsorted]
                    else:
                        return _np.atleast_2d(n)
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
