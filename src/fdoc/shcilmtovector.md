# SHCilmToVector

Convert a three-dimensional array of real spherical harmonic coefficients to a 1-dimensional indexed vector.

# Usage

call SHCilmToVector (`cilm`, `vector`, [`lmax`])

# Parameters

`cilm` : input, real\*8, dimension (2, `lmax`+1, `lmax`+1)
:   The input real spherical harmonic coefficients.
	
`vector` : output, real\*8, dimension ( (`lmax`+1)\*\*2 )
:   The indexed output real spherical harmonic coefficients.

`lmax` : input, integer
:   The maximum degree of the output coefficients to convert.

# Description

`SHCilmToVector` will convert a three-dimensional array of real spherical harmonic coefficients to a 1-dimensional indexed array.  The degree `l`, order `m`, and `i` (1 for cosine, 2 for sine) corresponds to the index `l**2+(i-1)*l+M+1`.

# See also

[shvectortocilm](shvectortocilm.html), [yilmindex](yilmindex.html), [shcindextocilm](pyshcindextocilm.html), [shcilmtocindex](pyshcilmtocindex.html)
