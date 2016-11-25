"""
pyshtools Global Spectral Analysis Routines.

This supackage of pyshtools defines the following functions:

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

from ..shtools import SHPowerL
from ..shtools import SHPowerDensityL
from ..shtools import SHCrossPowerL
from ..shtools import SHCrossPowerDensityL
from ..shtools import SHPowerSpectrum
from ..shtools import SHPowerSpectrumDensity
from ..shtools import SHCrossPowerSpectrum
from ..shtools import SHCrossPowerSpectrumDensity
from ..shtools import SHAdmitCorr
from ..shtools import SHConfidence
from ..shtools import SHPowerLC
from ..shtools import SHPowerDensityLC
from ..shtools import SHCrossPowerLC
from ..shtools import SHCrossPowerDensityLC
from ..shtools import SHPowerSpectrumC
from ..shtools import SHPowerSpectrumDensityC
from ..shtools import SHCrossPowerSpectrumC
from ..shtools import SHCrossPowerSpectrumDensityC
