from cosmic_dustbox import sdist
import astropy.units as _u

sdists = {
    'gra': sdist.PowerLaw(
        50*_u.angstrom, 0.25*_u.micron, 10**-25.13*_u.cm**2.5, -3.5),
    'sil': sdist.PowerLaw(
        50*_u.angstrom, 0.25*_u.micron, 10**-25.11*_u.cm**2.5, -3.5)
}
