# SHCilmToCindex

Convert a three-dimensional array of spherical harmonic coefficients to a two-dimensional indexed array.

# Usage

call SHCilmToCindex (`cilm`, `cindex`, `degmax`, `exitstatus`)

# Parameters

`cilm` : input, real\*8, dimension (2, `lmaxin`+1, `lmaxin`+1)
:   The input spherical harmonic coefficients. `cilm(1,:,:)` and `cilm(2,:,:)` correspond to either the real and imaginary components, or cosine and sine coefficients, respectively.

`cindex` : output, real\*8, dimension (2, (`degmax`+1)\*(`degmax`+2)/2)
:   The indexed output spherical harmonic coefficients.

`degmax` : input, optional, integer, default = `lmaxin`
:   The maximum degree of the output coefficients.

`exitstatus` : output, optional, integer
:   If present, instead of executing a STOP when an error is encountered, the variable exitstatus will be returned describing the error. 0 = No errors; 1 = Improper dimensions of input array; 2 = Improper bounds for input variable; 3 = Error allocating memory; 4 = File IO error.

# Description

`SHCilmToCindex` will convert a three-dimensional array of spherical harmonic coefficients to a two-dimensional indexed array.  The degree `l` and order `m` corresponds to the index `l*(l+1)/2+m+1`. The default is to convert the entire array `cilm`, but a subset of this array can be converted by specifying the optional argument `degmax`.

# See also

[shcindextocilm](shcindextocilm.html), [shcilmtovector](shcilmtovector.html), [shvectortocilm](shvectortocilm.html)
