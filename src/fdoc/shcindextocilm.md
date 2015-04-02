# SHCindexToCilm

Convert a two-dimensional indexed array of spherical harmonic coefficients to a three-dimensional array.

# Usage

call subroutine SHCindexToCilm (`cindex`, `cilm`, `degmax`)

# Parameters

`cindex` : input, real\*8, dimension (2, :)
:   The indexed spherical harmonic coefficients. The second dimension of this array must be greater than `(lmax+1)*(lmax+2)/2` where `lmax` is the maximum degree of `cilm`.

`cilm` : output, real\*8, dimension (2, :, :)
:   The input spherical harmonic coefficients. `cilm(1,:,:)` and `cilm(2,:,:)` correspond to either the real and imaginary components, or cosine and sine coefficients, respectively.

`degmax` : input, integer, optional
:   The maximum degree of the output coefficients. By default, the entire array `cilm` will be converted.

# Description

`SHCindexToCilm` will convert a two-dimensional indexed array of spherical harmonic coefficients to a three-dimensional array of complex spherical harmonic coefficients.  The degree `l` and order `m` corresponds to the index `l*(l+1)/2+m+1`. The default is to convert the entire array `cindex`, but a subset of this array can be converted by specifying the optional argument `degmax`.

# See also

[shcindextocilm](shcindextocilm.html), [shcindextovector](pyshcindextovector.html), [shvectortocilm](pyshvectortocilm.html)

