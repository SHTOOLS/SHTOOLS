# SHPowerLC

Compute the power of a complex function for a single spherical harmonic degree.

# Usage

`power` = SHPowerLC (`cilm`, `l`)

# Returns

`power` : float
:   Power of the complex function for spherical harmonic degree `l`.

# Parameters

`cilm` : complex, dimension (2, `lmaxin`+1, `lmaxin`+1)
:   The complex spherical harmonics of the complex function.

`l` : integer
:   The spherical harmonic degree. This must be less than or equal to `lmaxin`.

# Description

`SHPowerLC` will calculate the power of a complex function expressed in complex 4-pi normalized spherical harmonics for a single spherical harmonic degree `l`. This is calculated as:

`power = Sum_{i=0}^1 Sum_{m=0}^l | cilm[i, l, m] |**2`.

# See also

[shpowerdensitylc](pyshpowerdensitylc.html), [shcrosspowerlc](pyshcrosspowerlc.html), [shcrosspowerdensitylc](pyshcrosspowerdensitylc.html), [shpowerspectrumc](pyshpowerspectrumc.html), [shpowerspectrumdensityc](pyshpowerspectrumdensityc.html), [shcrosspowerspectrumc](pyshcrosspowerspectrumc.html), [shcrosspowerspectrumdensityc](pyshcrosspowerspectrumdensityc.html)
