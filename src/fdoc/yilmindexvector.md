# YilmIndexVector

Compute the index of an 1D array of spherical harmonic coefficients corresponding to `i`, `l`, and `m`.

# Usage

`index` = YilmIndex (`i`, `l`, `m`)

# Parameters

`index` : output, integer 
:   Index of an 1D array of spherical harmonic coefficients corresponding to `i`, `l`, and `m`.

`i` : input, integer
:   1 corresponds to the cosine coefficient `cilm(1,:,:)`, and 2 corresponds to the sine coefficient `cilm(2,:,:)`.

`l` : input, integer
:   The spherical harmonic degree.

`m` : input, integer
:   The angular order.

# Description

`YilmIndexVector` will calculate the index of a 1D vector of spherical harmonic coefficients corresponding to degree `l`, angular order `m` and `i` (1 = cosine, 2 = sine). The index is given by `l**2+(i-1)*l+m+1`.

# See also

[shcilmtovector](shcilmtovector.html), [shvectortocilm](shvectortocilm.html)
