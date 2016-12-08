# SHVectorToCilm

Convert a 1-dimensional indexed vector of real spherical harmonic coefficients to a three-dimensional array.

# Usage

call SHVectorToCilm (`vector`, `cilm`, `lmax`, `exitstatus`)

# Parameters

`vector` : input, real\*8, dimension ( (`lmax`+1)\*\*2 )
:   The input 1-D indexed array of real spherical harmonic coefficients.

`cilm` : output, real\*8, dimension (2, `lmax`+1, `lmax`+1)
:   The 3-D arrary of output real spherical harmonic coefficients.

`lmax` : input, integer
:   The maximum degree of the output coefficients.

`exitstatus` : output, optional, integer
:   If present, instead of executing a STOP when an error is encountered, the variable exitstatus will be returned describing the error. 0 = No errors; 1 = Improper dimensions of input array; 2 = Improper bounds for input variable; 3 = Error allocating memory; 4 = File IO error.

# Description

`SHVectorToCilm` will convert a 1-dimensional indexed vector of real spherical harmonic coefficients to a three-dimensional array.  The degree `l`, order `m`, and `i` (1 = cosine, 2 = sine) corresponds to the index `l**2+(i-1)*l+m+1.

# See also

[shvectortocilm](shvectortocilm.html), [yilmindexvector](yilmindexvector.html), [shcilmtocindex](shcilmtocindex.html), [shcindextocilm](shcindextocilm.html)
