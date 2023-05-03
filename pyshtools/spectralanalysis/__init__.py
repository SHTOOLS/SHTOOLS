"""
pyshtools Global and Localized Spectral Analysis Routines.

This subpackage of pyshtools defines the following functions:

Global spectral analysis
------------------------
spectrum               Calculate the spectrum of a real or complex function.
cross_spectrum         Calculate the cross-spectrum of two real or complex
                       functions.
SHAdmitCorr            Calculate the admittance and correlation spectra
                       of two functions.
SHConfidence           Compute the probability that two sets of
                       spherical harmonic coefficients are correlated at
                       a given degree and for a given correlation coefficient.

Multitaper spectral estimation (spherical cap domain)
-----------------------------------------------------
SHMultiTaperSE         Perform a localized multitaper spectral analysis.
SHMultiTaperCSE        Perform a localized multitaper cross-spectral analysis.
SHLocalizedAdmitCorr   Calculate the localized admittance and correlation
                       spectra of two functions at a given location.
SHReturnTapers         Calculate the eigenfunctions of the spherical-cap
                       concentration problem.
SHReturnTapersM        Calculate the eigenfunctions of the spherical-cap
                       concentration problem for a single angular order.
ComputeDm              Compute the space-concentration kernel of a spherical
                       cap.
ComputeDG82            Compute the tridiagonal matrix of Grunbaum et al. (1982)
                       that commutes with the space-concentration kernel of a
                       spherical cap.
SHFindLWin             Determine the spherical-harmonic bandwidth that is
                       necessary to achieve a certain concentration factor.
SHBiasK                Calculate the multitaper (cross-)power spectrum
                       expectation of a windowed function.
SHMTCouplingMatrix     Compute the multitaper coupling matrix for a given set
                       of localization windows.
SHBiasAdmitCorr        Calculate the expected multitaper admittance and
                       correlation spectra associated with the input global
                       cross-power spectra of two functions.
SHMTDebias             Invert for the global power spectrum given a localized
                       multitaper spectrum estimate.
SHMTVarOpt             Calculate the minimum variance and corresponding optimal
                       weights of a localized multitaper spectral estimate.
SHMTVar                Calculate the theoretical variance of a multitaper
                       spectral estimate for a given input power spectrum.
SHSjkPG                Calculate the expectation of the product of two
                       functions, each multiplied by a different data taper,
                       for a given spherical harmonic degree and two different
                       angular orders.
SHRotateTapers         Rotate the spherical cap tapers by three Euler angles.

Localization windows (arbitrary domain)
---------------------------------------
SHReturnTapersMap      Calculate the eigenfunctions of the concentration
                       problem for an arbitrary concentration region.
SHMultiTaperMaskSE     Perform a localized multitaper spectral analysis using
                       arbitrary windows.
SHMultiTaperMaskCSE    Perform a localized multitaper cross-spectral analysis
                       using arbitrary windows.
SHBiasKMask            Calculate the multitaper (cross-)power spectrum
                       expectation of a function localized by arbitrary
                       windows derived from a mask.
ComputeDMap            Compute the space-concentration kernel of a mask
                       defined on the sphere.
Curve2Mask             Given a set of latitude and longitude coordinates
                       representing a closed curve, output a gridded mask.

Localization Bias (General)
---------------------------
SHBias                 Calculate the (cross-)power spectrum expectation of a
                       windowed function.

Slepian function expansions
---------------------------
SlepianCoeffs          Determine the expansion coefficients of a function for
                       a given set of input Slepian functions.
SlepianCoeffsToSH      Convert a function expressed in Slepian coefficients to
                       spherical harmonic coefficients.
SHSCouplingMatrix      Compute the spherical harmonic coupling matrix for a
                       given set of Slepian functions.
SHSCouplingMatrixCap   Compute the spherical harmonic coupling matrix for a
                       given set of spherical-cap Slepian functions.
SHSlepianVar           Calculate the theoretical variance of the power of a
                       function expanded in spherical-cap Slepian functions.

Other
-----
SphericalCapCoef       Calculate the spherical harmonic coefficients of a
                       spherical cap.
"""
from ..backends.shtools import SHAdmitCorr
from ..backends.shtools import SHConfidence
from ..backends.shtools import SHMultiTaperSE
from ..backends.shtools import SHMultiTaperCSE
from ..backends.shtools import SHLocalizedAdmitCorr
from ..backends.shtools import SHReturnTapers
from ..backends.shtools import SHReturnTapersM
from ..backends.shtools import ComputeDm
from ..backends.shtools import ComputeDG82
from ..backends.shtools import SHFindLWin
from ..backends.shtools import SHBiasK
from ..backends.shtools import SHMTCouplingMatrix
from ..backends.shtools import SHBiasAdmitCorr
from ..backends.shtools import SHMTDebias
from ..backends.shtools import SHMTVarOpt
from ..backends.shtools import SHMTVar
from ..backends.shtools import SHSjkPG
from ..backends.shtools import SHMultiTaperMaskSE
from ..backends.shtools import SHMultiTaperMaskCSE
from ..backends.shtools import SHReturnTapersMap
from ..backends.shtools import SHBiasKMask
from ..backends.shtools import ComputeDMap
from ..backends.shtools import Curve2Mask
from ..backends.shtools import SHBias
from ..backends.shtools import SphericalCapCoef
from ..backends.shtools import SHRotateTapers
from ..backends.shtools import SlepianCoeffs
from ..backends.shtools import SlepianCoeffsToSH
from ..backends.shtools import SHSCouplingMatrix
from ..backends.shtools import SHSlepianVar
from ..backends.shtools import SHSCouplingMatrixCap

from .spectrum import spectrum
from .cross_spectrum import cross_spectrum


# ---- Define __all__ for use with: from pyshtools import * ----
__all__ = ['SHAdmitCorr', 'SHConfidence', 'spectrum', 'cross_spectrum',
           'SHMultiTaperSE', 'SHMultiTaperCSE', 'SHLocalizedAdmitCorr',
           'SHReturnTapers', 'SHReturnTapersM', 'ComputeDm', 'ComputeDG82',
           'SHFindLWin', 'SHBiasK', 'SHMTCouplingMatrix', 'SHBiasAdmitCorr',
           'SHMTDebias', 'SHMTVarOpt', 'SHMTVar', 'SHSjkPG',
           'SHMultiTaperMaskSE', 'SHMultiTaperMaskCSE', 'SHReturnTapersMap',
           'SHBiasKMask', 'ComputeDMap', 'Curve2Mask', 'SHBias',
           'SphericalCapCoef', 'SHRotateTapers', 'SlepianCoeffs',
           'SlepianCoeffsToSH', 'SHSCouplingMatrix', 'SHSlepianVar',
           'SHSCouplingMatrixCap']
