# SHCrossPowerDensityL

Compute the cross-power spectral density of two real functions for a single spherical harmonic degree.

# Usage

`cpsd` = pyshtools.SHCrossPowerDensityL (`cilm1`, `cilm2`, `l`)

# Returns

`cpsd` : float
:   Cross-power spectral density of the two real functions for spherical harmonic degree `l`.

# Parameters

`cilm1` : float, dimension (2, `lmaxin1`+1, `lmaxin1`+1)
:   The spherical harmonic coefficients of the first real function.

`cilm2` : float, dimension (2, `lmaxin2`+1, `lmaxin2`+1)
:   The spherical harmonic coefficients of the second real function.

`l` : integer
:   The spherical harmonic degree. This must be less than or equal to the minimum of `lmaxin1` and `lmaxin2`.

# Description

`SHCrossPowerDensityL` will calculate the cross-power spectral density of two real functions expressed in real spherical harmonics for a single spherical harmonic degree `l`. This is explicitly calculated as:

`cpsd = Sum_{i=0}^1 Sum_{m=0}^l cilm1[i, l, m] * cilm2[i, l, m] / (2l + 1)`.

# See also

[shpowerl](pyshpowerl.html), [shpowerdensityl](pyshpowerdensityl.html), [shcrosspowerl](pyshcrosspowerl.html), [shpowerspectrum](pyshpowerspectrum.html), [shpowerspectrumdensity](pyshpowerspectrumdensity.html), [shcrosspowerspectrum](pyshcrosspowerspectrum.html), [shcrosspowerspectrumdensity](pyshcrosspowerspectrumdensity.html), [shadmitcorr](pyshadmitcorr.html)
