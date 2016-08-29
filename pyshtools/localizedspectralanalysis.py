"""
pyshtools Localized Spectral Analysis Routines.

This submodule of pyshtools defines the following functions:

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
SHMTCouplingMatrix     Calculate the multitaper coupling matrix for a given
                       set of localization windows.
SHBiasAdmitCorr        Calculate the expected multitaper admittance and
                       correlation spectra associated with the input global
                       cross-power spectra of two functions.
SHMTDebias             Invert for the global power spectrum given a localized
                       multitaper spectrum estimate.
SHMTVarOpt             Calculate the minimum variance and corresponding optimal
                       weights of a localized multitaper spectral estimate.
SHSjkPG                Calculate the expectation of the product of two
                       functions, each multiplied by a different data taper,
                       for a given spherical harmonic degree and two different
                       angular orders.

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

Other
-----
SphericalCapCoef       Calculate the spherical harmonic coefficients of a
                       spherical cap.
"""

from __future__ import absolute_import as _absolute_import
from __future__ import division as _division
from __future__ import print_function as _print_function

from ._SHTOOLS import SHMultiTaperSE
from ._SHTOOLS import SHMultiTaperCSE
from ._SHTOOLS import SHLocalizedAdmitCorr
from ._SHTOOLS import SHReturnTapers
from ._SHTOOLS import SHReturnTapersM
from ._SHTOOLS import ComputeDm
from ._SHTOOLS import ComputeDG82
from ._SHTOOLS import SHFindLWin
from ._SHTOOLS import SHBiasK
from ._SHTOOLS import SHMTCouplingMatrix
from ._SHTOOLS import SHBiasAdmitCorr
from ._SHTOOLS import SHMTDebias
from ._SHTOOLS import SHMTVarOpt
from ._SHTOOLS import SHSjkPG
from ._SHTOOLS import SHMultiTaperMaskSE
from ._SHTOOLS import SHMultiTaperMaskCSE
from ._SHTOOLS import SHReturnTapersMap
from ._SHTOOLS import SHBiasKMask
from ._SHTOOLS import ComputeDMap
from ._SHTOOLS import Curve2Mask
from ._SHTOOLS import SHBias
from ._SHTOOLS import SphericalCapCoef
