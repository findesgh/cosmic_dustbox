###############################################################################
import astropy.units as _u
import numpy as _np
import os as _os
from scipy.special import erf as _erf
###############################################################################


class SizeDist(object):
    """
    Grain size distribution.

    Computes 1/n_H * dn_gr/da [1/L], where n_H is the hydrogen nuclei number
    density, n_gr(a) is the number density of grains <= a and a is a grain
    size. Note that this is not a distribution in the mathematical sense. Since
    it is commonly referred to as such in the community, though, and for lack
    of a better term, this convention is maintained here. It might, however, be
    changed in the future. See [1]_ and specifically Sec. 1 for an example.

    Parameters
    ----------
    sizeMin : scalar Quantity [L]
        Low end size cutoff of the distribution.
    sizeMax : scalar Quantity [L]
        High end size cutoff of the distribution.
    func : callable
        Must take a single Quantity with array value of sizes and return
        1/n_H * dn_gr/da [1/L]. It is called in the limits set by `sizeMin`
        and `sizeMax` when an instance of this class is called.

    Attributes
    ----------
    sizeMin : scalar Quantity [L]
        Directly taken from parameters.
    sizeMax : scalar Quantity [L]
        Directly taken from parameters.
    func : callable
        Directly taken from parameters.

    References
    ----------
    .. [1] Weingartner & Draine 2001ApJ...548..296W

    Examples
    --------
    >>> import astropy.units as u
    >>> def f(s): \
            return 1.0/s.unit
    >>> a = SizeDist(3.5*u.angstrom, 1*u.micron, f)
    >>> a(_np.logspace(-11, -5, 10)*u.m)
    <Quantity [0., 0., 0., 1., 1., 1., 1., 1., 0., 0.] 1 / m>
    """
    def __init__(self, sizeMin, sizeMax, func):
        """
        See class docstring.
        """
        self.sizeMin = sizeMin
        self.sizeMax = sizeMax
        self.func = func
        return

    def __call__(self, sizes):
        """
        Compute 1/n_H * dn_gr/da [1/L].

        Returns func(sizes) in limits set by sizeMin and sizeMax and zero
        elsewhere.

        Parameters
        ----------
        sizes : array-valued Quantity [L]
            Grain sizes at which to evaluate func.

        Returns
        -------
        r : array-valued Quantity [1/L]
        """
        r = _np.zeros(len(sizes))/sizes.unit
        ind = _np.where((sizes >= self.sizeMin) & (sizes <= self.sizeMax))
        r[ind] = self.func(sizes[ind])
        return r

    def __add__(self, other):
        """
        Commutative addition of size distributions.

        Overloading of ``+`` operator. Can add an instance of ``SizeDist`` (or
        subclass thereof), any callable or a scalar to ``self.func``. No check
        on units is performed. Returns an instance of ``SizeDist``, the
        ``func`` attribute of which is defined to return the corresponding sum.

        If ``other`` is an instance of ``SizeDist`` (or subclass), take maximum
        of the two ``sizeMin`` attributes and minimum of the two ``sizeMax``
        attributes.

        Parameters
        ----------
        other : SizeDist, callable or scalar-valued Quantity [1/L]
            Is added to ``self.func``.

        Returns
        -------
         : ``self.__class__``
            Instance of own class with corresponding ``func`` attribute.

        Examples
        --------
        >>> import astropy.units as u

        >>> def f(s): \
                return 1.0/s.unit
        >>> a = SizeDist(1*u.angstrom, 1*u.micron, f)
        >>> b = SizeDist(10*u.angstrom, 10*u.micron, f)
        >>> c = a + b
        >>> c(_np.logspace(-11, -4, 10)*u.m)
        <Quantity [0., 0., 0., 2., 2., 2., 2., 0., 0., 0.] 1 / m>
        """
        if issubclass(other.__class__, SizeDist):
            # find new size limits
            sizeMin = max(self.sizeMin, other.sizeMin)
            sizeMax = min(self.sizeMax, other.sizeMax)

            # new differential number density is sum
            def func(sizes):
                return self.func(sizes) + other.func(sizes)

            return SizeDist(sizeMin, sizeMax, func)
        elif callable(other):
            def func(sizes):
                return self.func(sizes) + other(sizes)
            return SizeDist(self.sizeMin, self.sizeMax, func)
        else:
            def func(sizes):
                return other + self.func(sizes)
            return SizeDist(self.sizeMin, self.sizeMax, func)

    # make addition commutative
    __radd__ = __add__

    def __mul__(self, other):
        if callable(other):
            def func(sizes):
                return self.function(sizes) * other(sizes)
            return self.__class__(self.sizeMin, self.sizeMax, func)
        else:
            def func(sizes):
                return other * self.function(sizes)
            return self.__class__(self.sizeMin, self.sizeMax, func)

    # make multiplication commutative
    __rmul__ = __mul__


class PowerLaw(SizeDist):
    """
    Parameters
    ----------
    sizeMin : scalar Quantity [L]
        Low end size cutoff of the distribution.
    sizeMax : scalar Quantity [L]
        High end size cutoff of the distribution.
    power : float
        Log-slope of the size distribution.
    C : scalar Quantity [L**(-1-power)]
        Normalization of the size distribution.
    """
    def __init__(self, sizeMin, sizeMax, power, C):
        def f(a):
            return C*a**power
        super().__init__(sizeMin, sizeMax, f)
        return


class LogNormal(SizeDist):
    """
    Log normal distribution as in Eq. (2) of [1]_.

    References
    ----------
    .. [1] Weingartner & Draine 2001ApJ...548..296W
    """
    def __init__(self, sizeMin, sizeMax, rho, sgma, bc, a0):

        # mass of carbon atom
        m_C = 12.0107*_u.u

        nominator = (3 * _np.exp(-4.5 * sgma**2) * bc * m_C)
        denominator = (
            (2*_np.pi**2)**1.5 * rho * a0**3 * sgma *
            (1 +
             _erf((3 * sgma/_np.sqrt(2)) +
                  (_np.log(a0/3.5/_u.angstrom) /
                   (sgma * _np.sqrt(2)))
             )
            )
        )
        # FIXME: what do we do if the denominator is zero
        if denominator != 0:
            B = (nominator/denominator).decompose()

        def f(a):
            return B/a * _np.exp(-0.5*(
                _np.log(a/a0)/sgma)**2)

        super().__init__(sizeMin, sizeMax, f)
        return


class WD01ExpCutoff(SizeDist):
    """
    Power law distribution with exponential cutoff as in Eqs. (4) and (5) of
    [1]_.

    References
    ----------
    .. [1] Weingartner & Draine 2001ApJ...548..296W
    """
    def __init__(self, sizeMin, sizeMax, alpha, beta, a_t, a_c, C):

        if beta >= 0:
            def F(a):
                return 1.0 + (beta * a) / a_t
        else:
            def F(a):
                return 1.0 / (1 - (beta * a) / a_t)

        def exp_func(a):
            r = _np.ones_like(a.value)
            ind = _np.where(a > a_t)
            r[ind] = _np.exp(-(((a[ind] - a_t)/a_c))**3)
            return r

        def f(a):
            return C/a*(a / a_t)**alpha \
                * F(a) * exp_func(a)

        super().__init__(sizeMin, sizeMax, f)
        return


def WD01(RV, bc, case):
    # pandas dataframe would be far better suited for this but don't know if
    # pulling in another dependency makes sense just for this
    data = _np.genfromtxt(
        _os.path.join(
            _os.path.dirname(__file__),
            'data/sdist/WD01-2001ApJ...548..296W-Table1.txt'
        ),
        dtype=None,
        encoding=None
    )
    # doing == float comparison here, might come back to bite us at some point
    params = data[
        [
            j for j, x in enumerate(data)
            if (x[0] == RV and x[1] == bc and x[3] == case)
        ]
    ]

    if len(params) > 1:
        # this should never happen
        raise ValueError("Could not uniquely identify parameter set.")
    elif len(params) == 0:
        raise ValueError("Could not identify parameter set.")

    params = params[0]

    sizeMin = 3.5*_u.angstrom
    sizeMax = 10*_u.micron
    rho = 2.24*_u.g/_u.cm**3
    s_car = LogNormal(
        sizeMin,
        sizeMax,
        rho,
        0.4,
        0.75*params[2],
        3.5*_u.angstrom
    ) + \
    LogNormal(
        sizeMin,
        sizeMax,
        rho,
        0.4,
        0.25*params[2],
        30*_u.angstrom
    )
    l_car = WD01ExpCutoff(
        sizeMin,
        sizeMax,
        params[4],
        params[5],
        params[6]*_u.m,
        params[7]*_u.m,
        params[8]
    )
    sil = WD01ExpCutoff(
        sizeMin,
        sizeMax,
        params[9],
        params[10],
        params[11]*_u.m,
        params[12]*_u.m,
        params[13]
    )
    return s_car + l_car, sil


###############################################################################
if __name__ == "__main__":
    import doctest
    doctest.testmod()
###############################################################################
