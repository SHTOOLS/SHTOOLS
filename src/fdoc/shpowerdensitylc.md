# SHPowerDensityLC

Compute the power spectral density of a complex function for a single spherical harmonic degree.

# Usage

`psd` = SHPowerDensityLC (`cilm`, `l`)

# Parameters

`psd` : output, real\*8
:   Power spectral density of the complex function for spherical harmonic degree `l`.

`cilm` : input, complex\*16, dimension (2, `lmaxin`+1, `lmaxin`+1)
:   The complex function expressed in complex spherical harmonics.

`l` : input, integer
:   The spherical harmonic degree. This must be less than or equal to `lmaxin`.

# Description

`SHPowerDensityLC` will calculate the power spectral density of a complex function expressed in complex 4-pi normalized spherical harmonics for a single spherical harmonic degree `l`. This is calculated as:

`psd = Sum_{i=1}^2 Sum_{m=0}^l | cilm(i, l+1, m+1) |**2 / (2l + 1)`.

# See also

[shpowerlc](shpowerlc.html), [shcrosspowerlc](shcrosspowerlc.html), [shcrosspowerdensitylc](shcrosspowerdensitylc.html), [shpowerspectrumc](shpowerspectrumc.html), [shpowerspectrumdensityc](shpowerspectrumdensityc.html), [shcrosspowerspectrumc](shcrosspowerspectrumc.html), [shcrosspowerspectrumdensityc](shcrosspowerspectrumdensityc.html)
