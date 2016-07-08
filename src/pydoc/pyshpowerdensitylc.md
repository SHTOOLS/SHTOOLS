# SHPowerDensityLC

Compute the power spectral density of a complex function for a single spherical harmonic degree.

# Usage

`psd` = pyshtools.SHPowerDensityLC (`cilm`, `l`)

# Returns

`psd` : float
:   Power spectral density of the complex function for spherical harmonic degree `l`.

# Parameters

`cilm` : complex, dimension (2, `lmaxin`+1, `lmaxin`+1)
:   The complex function expressed in complex spherical harmonics.

`l` : integer
:   The spherical harmonic degree. This must be less than or equal to `lmaxin`

# Description

`SHPowerDensityLC` will calculate the power spectral density of a complex function expressed in complex spherical harmonics for a single spherical harmonic degree `l`. This is calculated as:

`psd = Sum_{i=0}^1 Sum_{m=0}^l | cilm[i, l, m] |**2 / (2l + 1)`.

# See also

[shpowerlc](pyshpowerlc.html), [shcrosspowerlc](pyshcrosspowerlc.html), [shcrosspowerdensitylc](pyshcrosspowerdensitylc.html), [shpowerspectrumc](pyshpowerspectrumc.html), [shpowerspectrumdensityc](pyshpowerspectrumdensityc.html), [shcrosspowerspectrumc](pyshcrosspowerspectrumc.html), [shcrosspowerspectrumdensityc](pyshcrosspowerspectrumdensityc.html)
