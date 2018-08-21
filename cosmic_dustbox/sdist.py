###############################################################################
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
    <Quantity [ 0., 0., 0., 1., 1., 1., 1., 1., 0., 0.] 1 / m>
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


###############################################################################
if __name__ == "__main__":
    import doctest
    doctest.testmod()
###############################################################################
