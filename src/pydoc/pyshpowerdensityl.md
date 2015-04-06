# SHPowerDensityL

Compute the power spectral density of a real function for a single spherical harmonic degree.

# Usage

`psd` = pyshtools.SHPowerDensityL (`cilm`, `l`)

# Returns

`psd` : float
:   Power spectral density of the function for the spherical harmonic degree `l`.

# Parameters

`cilm` : float, dimension (2, `lmaxin`+1, `lmaxin`+1)
:   The function expressed in real spherical harmonics.
	
`l` : integer
:   The spherical harmonic degree. 

# Description

`SHPowerDensityL` will calculate the power spectral density of a function expressed in real spherical harmonics for a single degree `l`. This is explicitly calculated as:

`psd = Sum_{i=0}^1 Sum_{m=0}^l cilm[i, l, m]**2 / (2l + 1)`.

# See also

[shpowerl](pyshpowerl.html), [shcrosspowerl](pyshcrosspowerl.html), [shcrosspowerdensityl](pyshcrosspowerdensityl.html), [shpowerspectrum](pyshpowerspectrum.html), [shpowerspectrumdensity](pyshpowerspectrumdensity.html), [shcrosspowerspectrum](pyshcrosspowerspectrum.html), [shcrosspowerspectrumdensity](pyshcrosspowerspectrumdensity.html), [shadmitcorr](pyshadmitcorr.html)
