# SHCilmToVector

Convert a three-dimensional array of real spherical harmonic coefficients to a 1-dimensional indexed vector.

# Usage

call SHCilmToVector (`cilm`, `vector`, `lmax`, `exitstatus`)

# Parameters

`cilm` : input, real\*8, dimension (2, `lmax`+1, `lmax`+1)
:   The input real spherical harmonic coefficients.

`vector` : output, real\*8, dimension ( (`lmax`+1)\*\*2 )
:   The indexed output real spherical harmonic coefficients.

`lmax` : input, integer
:   The maximum degree of the output coefficients to convert.

`exitstatus` : output, optional, integer
:   If present, instead of executing a STOP when an error is encountered, the variable exitstatus will be returned describing the error. 0 = No errors; 1 = Improper dimensions of input array; 2 = Improper bounds for input variable; 3 = Error allocating memory; 4 = File IO error.

# Description

`SHCilmToVector` will convert a three-dimensional array of real spherical harmonic coefficients to a 1-dimensional indexed array.  The degree `l`, order `m`, and `i` (1 for cosine, 2 for sine) corresponds to the index `l**2+(i-1)*l+m+1`.

# See also

[shvectortocilm](shvectortocilm.html), [yilmindexvector](yilmindexvector.html), [shcindextocilm](shcindextocilm.html), [shcilmtocindex](shcilmtocindex.html)
