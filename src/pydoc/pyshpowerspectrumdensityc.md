# SHPowerSpectrumDensityC

Compute the power spectral density of a complex function.

# Usage

`pspectrum` = pyshtools.SHPowerSpectrumDensityC (`cilm`, [`lmax`])

# Returns

`pspectrum` : float, dimension (`lmax`+1)
:   The power spectral density of the function.

# Parameters

`cilm` : complex, dimension (2, `lmaxin`+1, `lmaxin`+1)
:   The complex function expressed in complex spherical harmonics.

`lmax` : optional, integer, default = `lmaxin`
:   The maximum spherical harmonic degree of the power spectral density. This must be less than or equal to `lmaxin`.

# Description

`SHPowerSpectrumDensityC` will calculate the power spectral density of a complex function expressed in complex spherical harmonics. For a given spherical harmonic degree `l`, this is calculated as:

`pspectrum(l) = Sum_{i=0}^1 Sum_{m=0}^l | cilm[i, l, m] |**2 / (2l + 1)`.

# See also

[shpowerlc](pyshpowerlc.html), [shpowerdensitylc](pyshpowerdensitylc.html), [shcrosspowerlc](pyshcrosspowerlc.html), [shcrosspowerdensitylc](pyshcrosspowerdensitylc.html), [shpowerspectrumc](pyshpowerspectrumc.html), [shcrosspowerspectrumc](pyshcrosspowerspectrumc.html), [shcrosspowerspectrumdensityc](pyshcrosspowerspectrumdensityc.html)
