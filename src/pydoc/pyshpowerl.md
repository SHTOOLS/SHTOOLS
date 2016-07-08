# SHPowerL

Compute the power of a real function for a single spherical harmonic degree.

# Usage

`power` = pyshtools.SHPowerL (`cilm`, `l`)

# Returns

`power` : float
:   Power of the function for a single spherical harmonic degree `l`.

# Parameters

`cilm` : float, dimension (2, `lmaxin`+1, `lmaxin`+1)
:   The spherical harmonic coefficients of the function.

`l` : integer
:   The spherical harmonic degree.

# Description

`SHPowerL` will calculate the power of a function expressed in real spherical harmonics for a single degree `l`. This is explicitly calculated as:

`power = Sum_{i=0}^1 Sum_{m=0}^l cilm[i, l, m]**2` .

# See also

[shpowerdensityl](pyshpowerdensityl.html), [shcrosspowerl](pyshcrosspowerl.html), [shcrosspowerdensityl](pyshcrosspowerdensityl.html), [shpowerspectrum](pyshpowerspectrum.html), [shpowerspectrumdensity](pyshpowerspectrumdensity.html), [shcrosspowerspectrum](pyshcrosspowerspectrum.html), [shcrosspowerspectrumdensity](pyshcrosspowerspectrumdensity.html), [shadmitcorr](pyshadmitcorr.html)
