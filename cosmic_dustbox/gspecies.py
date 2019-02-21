###############################################################################
import numpy as _np
import astropy.units as u
###############################################################################


class GrainSpecies(object):
    """
    Grain species base class.
    """
    def __init__(self):
        return

    def scatt(self, sizes, wave, **kwargs):
        return self._scatt(
            sizes.to(_u.micron),
            wave.to(_u.micron, equivalencies=_u.spectral()),
            **kwargs
        )


class SphericalGS(GrainSpecies):
    """
    Base class for spherical grain species.
    """
    def __init__(self):
        return

    def volume(self, sizes):
        return 4*_np.pi/3 * size**3

    def geomCS(self, sizes):
        return _np.pi * sizes**2

    def surfaceArea(self, sizes):
        return 4*_np.pi * sizes**2


class HomSphericalGS(SphericalGS):
    """
    Base class for homogeneous spherical grain species.

    Represents grain species for which Mie theory and/or anomalous diffraction
    theory can be used to compute optical properties.
    """
    def __init__(self):
        return

    def mass(self, sizes):
        return self.volume(sizes)*self.massDensity


class IsoHomSphericalGS(HomSphericalGS):
    """
    Base class for isotropic and homogeneous spherical grain species.

    For these grain species, Mie theory and anomalous diffraction theory can be
    directly employed.
    """
    def __init__(self, crefin, mieSolver):
        self.crefin = crefin
        return

    def _scatt(self, sizes, wave):
        return


class cAxisSphericalGS(HomSphericalGS):
    """
    Base class for homogeneous spherical grain species composed of anisotropic
    material.

    Represents grain species that are spherical but composed of a material that
    has crystal axes so that its properties depend on orientation. It is
    assumed that they can be treated in an "effective" way using Mie
    theory/anomalous diffraction theory by averaging properties along different
    axes appropriately. This is then supposed to represent the average
    properties of an ensemble of randomly oriented grains.

    The classic example of this are graphite grains where the crystal axis is
    perpendicular to the so-called basal plane and properties are computed
    using the "1/3 - 2/3 approximation", which has been shown to work well
    [1]_.

    References
    ----------
    .. [1] Draine & Malhotra 1993ApJ...414..632D
    """
    def __init__(self, crefins, axesWeights, mieSolver):
        self.crefins = crefins
        self.axesWeights = axesWeights
        self.solver = mieSolver
        return

    def _scatt(self, sizes, wave):
        return
