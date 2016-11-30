# ComputeDMap

Compute the space-concentration kernel of an arbitrary mask on the sphere.

# Usage

call ComputeDMap (`dij`, `dh_mask`, `n`, `lmax`, `sampling`, `exitstatus`)

# Parameters

`dij` : output, real\*8, dimension ( (`lmax`+1)\*\*2, (`lmax`+1)\*\*2 )
:   The space-concentration kernel corresponding to the mask dh_mask.

`dh_mask` : input, integer, dimension (`n`, `sampling`\*`n`)
:   A Driscoll and Healy (1994) sampled grid describing the concentration region R. All elements should either be 1 (for inside the concentration region) or 0 (for outside R).

`n` : input, integer
:   The number of latitudinal samples in `dh_mask`. The effective spherical harmonic bandwidth of this grid is `L=n/2-1`.

`lmax` : input, integer
:   The maximum spherical harmonic degree of the matrix `dij`.

`sampling` : input, optional, integer, default = 1
:   For 1 (default), `dh_mask` has `n` x `n` samples. For 2, `dh_mask` has `n` x `2n` samples.

`exitstatus` : output, optional, integer
:   If present, instead of executing a STOP when an error is encountered, the variable exitstatus will be returned describing the error. 0 = No errors; 1 = Improper dimensions of input array; 2 = Improper bounds for input variable; 3 = Error allocating memory; 4 = File IO error.

# Description

`ComputeDMap` will calculate the space-concentration kernel for a generic mask defined on the sphere. The input mask `dh_mask` must be sampled according to the Driscoll and Healy (1994) sampling theorem with `n` samples in latitude, and possess a value of 1 inside the concentration region, and 0 elsewhere. `dh_mask` can either possess `n` samples in longitude (`sampling=1`) or `2n` samples in longitude (`sampling=2`). Given the approximate way in which the elements of `dij` are calculated (see below), `sampling=2` should be preferred. `dij` is symmetric, and the elements are ordered according to the scheme described in `YilmIndexVector`. See Simons et al. (2006) for further details.

The elements of DIJ are explicitly given by

`Dlm,l'm' = 1/(4pi) Integral_R Ylm Yl'm' dOmega`,

where `R` is the concentration region. In this routine, all values of `l'm'` are calculated in a single spherical harmonic transform for a given value of `lm` according to

`Dl'm' = 1/(4pi) Integral_Omega F Yl'm' dOmega`.

where

`F = Ylm dh_mask`.

The function `F` is in general not a polynomial, and thus the coefficients `Dl'm'` should not be expected to be exact. For this reason, the effective spherical harmonic degree of the input mask (`L=n/2-1`) should be greater than `lmax`. The exact value of `n` should be chosen such that further increases in `n` do not alter the returned eigenvalues. The routine prints out the fractional area of the mask computed in the pixel domain divided by `D(1,1)` (the fractional area computed by the spherical harmonic transforms), and the ratio of the two should be close to 1. Experience suggests that `l` should be about 5 times `lmax`.

# References

Driscoll, J.R. and D.M. Healy, Computing Fourier transforms and convolutions on the 2-sphere, Adv. Appl. Math., 15, 202-250, 1994.

Simons, F.J., F.A. Dahlen, and M.A. Wieczorek, Spatiospectral concentration on a sphere, SIAM Review, 48, 504-536, 2006.

# See also

[shreturntapersmap](shreturntapersmap.html), [yilmindexvector](yilmindexvector.html)
