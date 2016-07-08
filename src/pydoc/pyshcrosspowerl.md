# SHCrossPowerL

Compute the cross-power of two real functions for a single spherical harmonic degree.

# Usage

`cpower` = pyshtools.SHCrossPowerL (`cilm1`, `cilm2`, `l`)

# Returns

`cpower` : float
:   The cross power of the two functions for spherical harmonic degree `l`.

# Parameters

`cilm1` : float, dimension (2, `lmaxin1`+1, `lmaxin1`+1)
:   The spherical harmonic coefficients of the first function.

`cilm2` : float, dimension (2, `lmaxin2`+1, `lmaxin2`+1)
:   The spherical harmonic coefficients of the second function.

`l` : integer
:   The spherical harmonic degree. This must be less than or equal to the minimum of `lmaxin1` and `lmaxin2`.

# Description

`SHCrossPowerL` will calculate the cross-power of two functions for a single spherical harmonic degree `l`. This is explicitly calculated as:

`cpower = Sum_{i=0}^1 Sum_{m=0}^l cilm1[i, l, m] * cilm2[i, l, m]`.

# See also

[shpowerl](pyshpowerl.html), [shpowerdensityl](pyshpowerdensityl.html), [shcrosspowerdensityl](pyshcrosspowerdensityl.html), [shpowerspectrum](pyshpowerspectrum.html), [shpowerspectrumdensity](pyshpowerspectrumdensity.html), [shcrosspowerspectrum](pyshcrosspowerspectrum.html), [shcrosspowerspectrumdensity](pyshcrosspowerspectrumdensity.html), [shadmitcorr](pyshadmitcorr.html)
