# SlepianCoeffs

Determine the expansion coefficients of a function for a given set of input Slepian functions.

# Usage

call SlepianCoeffs(`falpha`, `galpha`, `film`, `lmax`, `nmax`, `exitstatus`)

# Parameters

`falpha` : output, real(dp), dimension (`nmax`)
:   A vector containing the Slepian coefficients of the input function `film`.

`galpha` : input, real(dp), dimension ((`lmax`+1)**2, `nmax`)
:   An array containing the spherical harmonic coefficients of the Slepian functions. Each column corresponds to a single function of which the spherical harmonic coefficients can be unpacked with `SHVectorToCilm`.

`film` : input, real(dp), dimension (2, `lmax`+1, `lmax`+1)
:   The spherical harmonic coefficients of the global function to be expanded in Slepian functions.

`lmax` : input, integer
:   The spherical harmonic bandwidth of the Slepian functions.

`nmax` : input, integer
:   The number of expansion coefficients to compute. This must be less than or equal to (`lmax`+1)\*\*2.

`exitstatus` : output, optional, integer
:   If present, instead of executing a STOP when an error is encountered, the variable exitstatus will be returned describing the error. 0 = No errors; 1 = Improper dimensions of input array; 2 = Improper bounds for input variable; 3 = Error allocating memory; 4 = File IO error.

# Description

`SlepianCoeffs` will compute the Slepian coefficients of a global input function `film` given the Slepian functions `galpha`. The Slepian functions are determined by a call to either (1) `SHReturnTapers` and then `SHRotateTapers`, or (2) `SHReturnTapersMap`. Each row of `galpha` contains the (`lmax`+1)**2 spherical harmonic coefficients of a Slepian function that can be unpacked using `SHVectorToCilm`. The Slepian functions must be normalized to have unit power (that is the sum of the coefficients squared is 1), and the Slepian coefficients are calculated as

`f_alpha = sum_{ilm}^{lmax} f_ilm g(alpha)_lm`  

# See also

[slepiancoeffstosh](slepiancoeffstosh.html), [shreturntapers](shreturntapers.html), [shslepianvar](shslepianvar.html), [shreturntapersmap](shreturntapersmap.html), [shrotatetapers](shrotatetapers.html), [shvectortocilm](shvectortocilm.html), [shscouplingmatrix](shscouplingmatrix.html), [shscouplingmatrixcap](shscouplingmatrixcap.html), [shmtcouplingmatrixcap](shmtcouplingmatrixcap.html)
