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

from ._SHTOOLS import SHMultiTaperSE, SHMultiTaperCSE, SHLocalizedAdmitCorr
from ._SHTOOLS import SHReturnTapers, SHReturnTapersM, ComputeDm, ComputeDG82
from ._SHTOOLS import SHFindLWin, SHBiasK, SHMTCouplingMatrix, SHBiasAdmitCorr
from ._SHTOOLS import SHMTDebias, SHMTVarOpt, SHSjkPG
from ._SHTOOLS import SHReturnTapersMap, ComputeDMap, Curve2Mask, SHBias
from ._SHTOOLS import SphericalCapCoef
