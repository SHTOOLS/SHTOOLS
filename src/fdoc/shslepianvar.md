# SHSlepianVar

Calculate the theoretical variance of the power of a function expanded in spherical-cap Slepian functions for a given spherical harmonic degree.

# Usage

call SHSlepianVar (`l`, `galpha`, `galpha_order`, `lmax`, `kmax`, `sff`, `variance`, `exitstatus`)

# Parameters

`l` : input, integer
:   The spherical harmonic degree used to calculate the theoretical variance.

`galpha` : input, real\*8, dimension (`lmax`+1, `kmax`)
:   A matrix of Slepian functions obtained from `SHReturnTapers` or `SHReturnTapersM`.

`galpha_order` : input, integer, dimension (`kmax`)
:   The angular order of the Slepian functions in `galpha`.

`lmax` : input, integer
:   The spherical harmonic bandwidth of the Slepian functions.

`kmax` : input, integer
:   The maximum number of Slepian functions to use when calculating the variance.

`sff` : input, real\*8, dimension (`lmax`+1)
:   The global power spectrum of the function.

`variance` : output, real\*8
:   The theoretical variance of the spectral estimate for degree `l`.

`exitstatus` : output, optional, integer
:   If present, instead of executing a STOP when an error is encountered, the variable exitstatus will be returned describing the error. 0 = No errors; 1 = Improper dimensions of input array; 2 = Improper bounds for input variable; 3 = Error allocating memory; 4 = File IO error.

# Description

`SHSlepianVar` will compute the theoretical variance of the power of a function expanded in spherical-cap Slepian functions for a given spherical harmonic degree. This routine takes as input the spherical harmonic coefficients of spherical-cap Slepian functions as obtained by a call to `SHReturnTapers`, and only the first `KMAX` Slepian functions in the matrix `GALPHA` are used to compute the variance.

# See also

[shreturntapers](shreturntapers.html), [shreturntapersm](shreturntapersm.html), [slepiancoeffs](slepiancoeffs.html), [slepiancoeffstosh](slepiancoeffstosh.html)
