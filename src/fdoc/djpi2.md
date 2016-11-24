# djpi2

Compute the rotation matrix d(pi/2) used in rotating data expressed in spherical harmonics.

# Usage

call djpi2 (`dj`, `lmax`, `exitstatus`)

# Parameters

`dj` : output, real\*8, dimension (`lmax`+1, `lmax`+1, `lmax`+1)
:   The rotation matrix dj(pi/2).

`lmax` : input, integer
:   The maximum spherical harmonic degree of the spherical harmonic rotation.

`exitstatus` : output, optional, integer
:   If present, instead of executing a STOP when an error is encountered, the variable exitstatus will be returned describing the error. 0 = No errors; 1 = Improper dimensions of input array; 2 = Improper bounds for input variable; 3 = Error allocating memory; 4 = File IO error.

# Description

`djpi2` will calculate the rotation matrix `d_{mM}^j (pi/2)` that is used in rotating spherical harmonics in the routines `SHRotateRealCoef` and `SHRotateCoef`.

This routine is based on code originally written by Guy Masters.

# See also

[shrotatecoef](shrotatecoef.html), [shrotaterealcoef](shrotaterealcoef.html)
