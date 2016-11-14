# SHPowerSpectrumC

Compute the power spectrum of a complex function.

# Usage

`pspectrum` = SHPowerSpectrumC (`cilm`, [`lmax`])

# Returns

`pspectrum` : float, dimension (`lmax`+1)
:   The power spectrum of the complex function.

# Parameters

`cilm` : complex, dimension (2, `lmaxin`+1, `lmaxin`+1)
:   The complex function expressed in complex spherical harmonics.

`lmax` : optional, integer, default = `lmaxin`
:   The maximum spherical harmonic degree of the power spectrum. This must be less than or equal to `lmaxin`.

# Description

`SHPowerSpectrumC` will calculate the power spectrum of a complex function expressed in complex spherical harmonics. For a given spherical harmonic degree `l`, this is  calculated as:

`pspectrum(l) = Sum_{i=0}^1 Sum_{m=0}^l | cilm[i, l, m] |**2`.

# See also

[shpowerlc](pyshpowerlc.html), [shpowerdensitylc](pyshpowerdensitylc.html), [shcrosspowerlc](pyshcrosspowerlc.html), [shcrosspowerdensitylc](pyshcrosspowerdensitylc.html), [shpowerspectrumdensityc](pyshpowerspectrumdensityc.html), [shcrosspowerspectrumc](pyshcrosspowerspectrumc.html), [shcrosspowerspectrumdensityc](pyshcrosspowerspectrumdensityc.html)
