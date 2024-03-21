'''
Datasets related to the asteroid (433) Eros.

Shape
-----
NLR_shape  :  Wieczorek (2024)
SPC_shape  :  Wieczorek (2024)

Gravity
-------
JGE15A01   :  Miller et al. (2002)
'''
from pooch import os_cache as _os_cache
from pooch import retrieve as _retrieve
from pooch import HTTPDownloader as _HTTPDownloader
from pooch import DOIDownloader as _DOIDownloader
from ..shclasses import SHGravCoeffs as _SHGravCoeffs
from ..shclasses import SHCoeffs as _SHCoeffs
from ..constants.Eros import angular_velocity as _omega


def NLR_shape(lmax=719):
    '''
    NLR_shape is a spherical harmonic model of the shape of asteroid (433)
    Eros based on laser altimer data obtained by the Near Earth Asteroid
    Rendezvous mission. The maximum spherical harmonic degree of the model is
    719, which has an effective spatial resolution of 8 pixels per degree. The
    coefficients are in units of meters.

    Parameters
    ----------
    lmax : int, optional, default = 719
        The maximum spherical harmonic degree to return.

    References
    ----------
    Wieczorek, M. (2024). Spherical harmonic models of the shape of (433) Eros
        (1.0.0) [Data set]. Zenodo. https://doi.org/10.5281/zenodo.10820812
    Zuber, M. T., Smith, D. E., Cheng, A. F., Garvin, J. B., Aharonson, O.,
        Cole, T. D., Dunn, P. J., Guo, Y., Lemoine, F. G., Neumann, G. A.,
        Rowlands, D. D., & Torrence, M. H. (2000). The Shape of 433 Eros from
        the NEAR-Shoemaker Laser Rangefinder. Science, 289(5487), 2097–2101.
        https://doi.org/10.1126/science.289.5487.2097
    '''
    fname = _retrieve(
        url="doi:10.5281/zenodo.10820812/Eros_NLR_shape_719.bshc.gz",  # noqa: E501
        known_hash="sha256:42563d6083319b0a96b0f931a2bbefcf8ebf9a98b835c78a1e258fc164965963",  # noqa: E501
        downloader=_DOIDownloader(progressbar=True),
        path=_os_cache('pyshtools'),
    )
    return _SHCoeffs.from_file(fname, lmax=lmax, name='NLR_shape (Eros)',
                               units='m', format='bshc')


def SPC_shape(lmax=511):
    '''
    SPC_shape is a spherical harmonic model of the shape of asteroid (433)
    Eros based on stere-photoclinometry shape model derived from imagery
    obtained by the Near Earth Asteroid Rendezvous mission. The maximum
    spherical harmonic degree of the model is 511, which has an effective
    spatial resolution of ~5.7 pixels per degree. The coefficients are in
    units of meters.

    Parameters
    ----------
    lmax : int, optional, default = 511
        The maximum spherical harmonic degree to return.

    References
    ----------
    Wieczorek, M. (2024). Spherical harmonic models of the shape of (433) Eros
        (1.0.0) [Data set]. Zenodo. https://doi.org/10.5281/zenodo.10820812
    Gaskell, R.W., Gaskell Eros Shape Model V1.1.
        urn:nasa:pds:gaskell.ast-eros.shape-model::1.1. NASA Planetary Data
        System, 2021; doi: 10.26033/d0gq-9427.
    '''
    fname = _retrieve(
        url="doi:10.5281/zenodo.10820812/Eros_SPC_shape_511.bshc.gz",  # noqa: E501
        known_hash="sha256:13c467661c96da9927d04e5c026afd1d6381786a295b23fb70c526daf995c05f",  # noqa: E501
        downloader=_DOIDownloader(progressbar=True),
        path=_os_cache('pyshtools'),
    )
    return _SHCoeffs.from_file(fname, lmax=lmax, name='SPC_shape (Eros)',
                               units='m', format='bshc')


def JGE15A01(lmax=15):
    '''
    JGE15A01 is a JPL 15 degree and order spherical harmonic model of the
    gravitational potential of asteroid (433) Eros.

    Parameters
    ----------
    lmax : int, optional, default = 15
        The maximum spherical harmonic degree to return.

    Reference
    ---------
    Miller, J. K., Konopliv, A. S., Antreasian, P. G., Bordi, J. J., Chesley,
        S., Helfrich, C. E., Owen, W. M., Wang, T. C., Williams, B. G.,
        Yeomans, D. K., & Scheeres, D. J. (2002). Determination of Shape,
        Gravity, and Rotational State of Asteroid 433 Eros. Icarus, 155(1),
        3–17. https://doi.org/10.1006/icar.2001.6753

    '''
    fname = _retrieve(
        url="https://sbnarchive.psi.edu/pds3/near/NEAR_A_5_COLLECTED_MODELS_V1_0/data/rss/n15acoeff.tab",  # noqa: E501
        known_hash="sha256:e08068e2ea5167bee685ae00a8596144964e1da71ab16c51f2328f642d0be90e",  # noqa: E501
        downloader=_HTTPDownloader(progressbar=True),
        path=_os_cache('pyshtools'),
    )
    return _SHGravCoeffs.from_file(fname, lmax=lmax, gm=4.46275472004e5,
                                   r0=16.e3, header=None, errors=True,
                                   name='JGE15A01 (Eros)', encoding='utf-8',
                                   omega=_omega.value)


__all__ = ['NLR_shape', 'SPC_shape', 'JGE15A01']
