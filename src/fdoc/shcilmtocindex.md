# SHCilmToCindex

Convert a three-dimensional array of spherical harmonic coefficients to a two-dimensional indexed array.

# Usage

call SHCilmToCindex (`cilm`, `cindex`, `degmax`)

# Parameters

`cilm` : input, real\*8, dimension (2, `lmaxin`+1, `lmaxin`+1)
:   The input spherical harmonic coefficients. `cilm(1,:,:)` and `cilm(2,:,:)` correspond to either the real and imaginary components, or cosine and sine coefficients, respectively.
	
`cindex` : output, real\*8, dimension (2, (`degmax`+1)\*(`degmax`+2)/2)
:   The indexed output spherical harmonic coefficients.

`degmax` : input, integer, optional, default = `lmaxin`
:   The maximum degree of the output coefficients. 

# Description

`SHCilmToCindex` will convert a three-dimensional array of spherical harmonic coefficients to a two-dimensional indexed array.  The degree `l` and order `m` corresponds to the index `l*(l+1)/2+m+1`. The default is to convert the entire array `cilm`, but a subset of this array can be converted by specifying the optional argument `degmax`.

# See also

[shcindextocilm](shcindextocilm.html), [shcindextovector](pyshcindextovector.html), [shvectortocilm](pyshvectortocilm.html)
