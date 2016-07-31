"""
pyshtools Global Spectral Analysis Routines.

This submodule of pyshtools defines the following functions:

Real spectral analysis
----------------------
SHPowerL                      Compute the power of a real function for a single
                              spherical harmonic degree.
SHPowerDensityL               Compute the power spectral density of a real
                              function for a single spherical harmonic degree.
SHCrossPowerL                 Compute the cross-power of two functions for a
                              single spherical harmonic degree.
SHCrossPowerDensityL          Compute the cross-power spectral density of two
                              functions for a single spherical harmonic degree.
SHPowerSpectrum               Compute the power spectrum of a function.
SHPowerSpectrumDensity        Compute the power spectral density of a function.
SHCrossPowerSpectrum          Compute the cross-power spectrum of two
                              functions.
SHCrossPowerSpectrumDensity   Compute the cross-power spectral density of two
                              functions.
SHAdmitCorr                   Calculate the admittance and correlation spectra
                              of two functions.
SHConfidence                  Compute the probability that two sets of
                              spherical harmonic coefficients are correlated at
                              a given degree and for a given correlation
                              coefficient.

Complex spectral analysis
-------------------------
SHPowerLC                     Compute the power of a complex function for a
                              single spherical harmonic degree.
SHPowerDensityLC              Compute the power spectral density of a complex
                              function for a single spherical harmonic degree.
SHCrossPowerLC                Compute the cross-power of two complex functions
                              for a single spherical harmonic degree.
SHCrossPowerDensityLC         Compute the cross-power spectral density of two
                              complex functions for a single spherical harmonic
                              degree.
SHPowerSpectrumC              Compute the power spectrum of a complex function.
SHPowerSpectrumDensityC       Compute the power spectral density of a complex
                              function.
SHCrossPowerSpectrumC         Compute the cross-power spectrum of two complex
                              functions.
SHCrossPowerSpectrumDensityC  Compute the cross-power spectral density of two
                              complex functions.
"""

from __future__ import absolute_import as _absolute_import
from __future__ import division as _division
from __future__ import print_function as _print_function

from ._SHTOOLS import SHPowerL
from ._SHTOOLS import SHPowerDensityL
from ._SHTOOLS import SHCrossPowerL
from ._SHTOOLS import SHCrossPowerDensityL
from ._SHTOOLS import SHPowerSpectrum
from ._SHTOOLS import SHPowerSpectrumDensity
from ._SHTOOLS import SHCrossPowerSpectrum
from ._SHTOOLS import SHCrossPowerSpectrumDensity
from ._SHTOOLS import SHAdmitCorr
from ._SHTOOLS import SHConfidence
from ._SHTOOLS import SHPowerLC
from ._SHTOOLS import SHPowerDensityLC
from ._SHTOOLS import SHCrossPowerLC
from ._SHTOOLS import SHCrossPowerDensityLC
from ._SHTOOLS import SHPowerSpectrumC
from ._SHTOOLS import SHPowerSpectrumDensityC
from ._SHTOOLS import SHCrossPowerSpectrumC
from ._SHTOOLS import SHCrossPowerSpectrumDensityC
