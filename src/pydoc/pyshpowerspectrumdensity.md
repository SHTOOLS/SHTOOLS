# SHPowerSpectrumDensity

Compute the power spectral density of a real function.

# Usage

`pspectrum` = SHPowerSpectrumDensity (`cilm`, [`lmax`])

# Returns

`pspectrum` : float, dimension (`lmax`+1)
:   The power spectral density of the function.

# Parameters

`cilm` : float, dimension (2, `lmaxin`+1, `lmaxin`+1)
:   The real function expressed in real spherical harmonics.

`lmax` : optional, integer, default = `lmaxin`
:   The maximum spherical harmonic degree used in calculating the power spectrum. This must be less than or equal to `lmaxin`.

# Description

`SHPowerSpectrumDensity` will calculate the power spectral density of a function expressed in real spherical harmonics. For a given spherical harmonic degree `l`, this is explicitly calculated as:

`pspectrum(l) = Sum_{i=0}^1 Sum_{m=0}^l cilm[i, l, m]**2 / (2l + 1)`.

# See also

[shpowerl](pyshpowerl.html), [shpowerdensityl](pyshpowerdensityl.html), [shcrosspowerl](pyshcrosspowerl.html), [shcrosspowerdensityl](pyshcrosspowerdensityl.html), [shpowerspectrum](pyshpowerspectrum.html), [shcrosspowerspectrum](pyshcrosspowerspectrum.html), [shcrosspowerspectrumdensity](pyshcrosspowerspectrumdensity.html), [shadmitcorr](pyshadmitcorr.html)
