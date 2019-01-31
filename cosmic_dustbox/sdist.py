###############################################################################
import astropy.units as _u
import numpy as _np
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
        on units is performed. Returns an instance of ``self.__class__``, the
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

            return self.__class__(sizeMin, sizeMax, func)
        elif callable(other):
            def func(sizes):
                return self.func(sizes) + other(sizes)
            return self.__class__(self.sizeMin, self.sizeMax, func)
        else:
            def func(sizes):
                return other + self.func(sizes)
            return self.__class__(self.sizeMin, self.sizeMax, func)

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
    C : scalar Quantity [L]**(-1-power)
        Normalization of the size distribution.
    """
    def __init__(self, sizeMin, sizeMax, power, C):
        def f(a):
            return C*a**power
        super(self).__init__(sizeMin, sizeMax, f)
        return


class Log_Normal(SizeDist):

    def __init__(self, sizeMin, sizeMax, rho, sgma, bc, b, a0):

        # mass of carbon atom
        m_C = 12.0107*_u.u

        def B_val(i):
            nominator = (3 * _np.exp(-4.5 * sgma**2) * bc[i] * m_C)
            denominator = (
                (2*_np.pi**2)**1.5 * rho * a0[i]**3 * sgma *
                (1 + _np.erf((3 * sgma/_np.sqrt(2)) +
                             (_np.log((a0[i]/3.5/_u.angstrom).decompose().value) / (sgma * _np.sqrt(2))))))
            # FIXME: what do we do if the denominator is zero
            if denominator != 0:
                B = (nominator/denominator).decompose().value
            return B

        _B = [B_val(i) for i in range(2)]

        def f(a):
            return sum([_B(i)/a * _np.exp(-0.5*(
                _np.log((a/a0[i]).decompose().value)/sgma)**2)
                        for i in range(2)])

        super(self).__init__(sizeMin, sizeMax, f)
        return


class WD01_dst(SizeDist):

        def __init__(self, sizeMin, sizeMax, a_t, beta, a_c, alpha, C):

            if beta >= 0:
                def F(a):
                    return 1 + (beta * a) / a_t
            else:
                def F(a):
                    return 1 / (1 - (beta * a) / a_t)

            def exp_func(a):
                r = _np.ones_like(a.value)
                ind = _np.where(a > a_t)
                r[ind] = _np.exp(-(((a[ind] - a_t)/a_c).decompose().value)**3)
                return r

            def f(a):
                return C*(a / a_t)**alpha/a \
                        * F(a) * exp_func(a)

            super(self).__init__(sizeMin, sizeMax, f)
            return


###############################################################################
if __name__ == "__main__":
    import doctest
    doctest.testmod()
###############################################################################
