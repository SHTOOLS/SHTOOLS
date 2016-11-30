# SHPowerLC

Compute the power of a complex function for a single spherical harmonic degree.

# Usage

`power` = SHPowerLC (`cilm`, `l`)

# Parameters

`power` : output, real\*8
:   Power of the complex function for spherical harmonic degree `l`.

`cilm` : input, complex\*16, dimension (2, `lmaxin`+1, `lmaxin`+1)
:   The complex spherical harmonics of the complex function.

`l` : input, integer
:   The spherical harmonic degree. This must be less than or equal to `lmaxin`.

# Description

`SHPowerLC` will calculate the power of a complex function expressed in complex 4-pi normalized spherical harmonics for a single spherical harmonic degree `l`. This is calculated as:

`power = Sum_{i=1}^2 Sum_{m=0}^l | cilm(i, l+1, m+1) |**2`.

# See also

[shpowerdensitylc](shpowerdensitylc.html), [shcrosspowerlc](shcrosspowerlc.html), [shcrosspowerdensitylc](shcrosspowerdensitylc.html), [shpowerspectrumc](shpowerspectrumc.html), [shpowerspectrumdensityc](shpowerspectrumdensityc.html), [shcrosspowerspectrumc](shcrosspowerspectrumc.html), [shcrosspowerspectrumdensityc](shcrosspowerspectrumdensityc.html)
