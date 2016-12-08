# SHCindexToCilm

Convert a two-dimensional indexed array of spherical harmonic coefficients to a three-dimensional array.

# Usage

call SHCindexToCilm (`cindex`, `cilm`, `degmax`, `exitstatus`)

# Parameters

`cindex` : input, real\*8, dimension (2, (`lmaxin`+1)\*(`lmaxin`+2)/2)
:   The indexed spherical harmonic coefficients.

`cilm` : output, real\*8, dimension (2, `degmax`+1, `degmax`+1)
:   The input spherical harmonic coefficients. `cilm(1,:,:)` and `cilm(2,:,:)` correspond to either the real and imaginary components, or cosine and sine coefficients, respectively.

`degmax` : input, optional, integer, default = `lmaxin`
:   The maximum degree of the output coefficients.

`exitstatus` : output, optional, integer
:   If present, instead of executing a STOP when an error is encountered, the variable exitstatus will be returned describing the error. 0 = No errors; 1 = Improper dimensions of input array; 2 = Improper bounds for input variable; 3 = Error allocating memory; 4 = File IO error.

# Description

`SHCindexToCilm` will convert a two-dimensional indexed array of spherical harmonic coefficients to a three-dimensional array of complex spherical harmonic coefficients.  The degree `l` and order `m` corresponds to the index `l*(l+1)/2+m+1`. The default is to convert the entire array `cindex`, but a subset of this array can be converted by specifying the optional argument `degmax`.

# See also

[shcilmtocindex](shcilmtocindex.html), [shcilmtovector](shcilmtovector.html), [shvectortocilm](shvectortocilm.html)

