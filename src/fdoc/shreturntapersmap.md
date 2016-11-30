# SHReturnTapersMap

Calculate the eigenfunctions and eigenvalues of the space-concentration problem for an arbitrary region.

# Usage

call SHReturnTapersMap (`tapers`, `eigenvalues`, `dh_mask`, `n`, `lmax`, `sampling`, `ntapers`, `exitstatus`)

# Parameters

`tapers` : input, real\*8, dimension ((`lmax`+1)\*\*2, `ntapers`)
:   The spherical harmonic coefficients of the tapers, arranged in columns, from best to worst concentrated. The spherical harmonic coefficients in each column are indexed according to the scheme described in `YilmIndexVector`.

`eigenvalues` : input, real\*8, dimension (`ntapers`)
:   The concentration factor for each localization window specified in the columns of `tapers`.

`dh_mask` : input, integer, dimension (`n`, `n`\*`sampling`)
:   A Driscoll and Healy (1994) sampled grid describing the concentration region R. All elements should either be 1 (for inside the concentration region) or 0 (for outside R).

`n` : input, integer
:   The number of latitudinal samples in `dh_mask`. The effective spherical harmonic bandwidth of this grid is `L=n/2-1`.

`lmax` : input, integer
:   The spherical harmonic bandwidth of the localization windows.

`sampling` : input, integer, default = 1
:   For 1 (default), `dh_mask` has `n x n` samples. For 2, `dh_mask` has `n x 2n` samples.

`ntapers` : input, optional, integer, default = (`lmax`+1)\*\*2
:   The number of best concentrated eigenvalues and corresponding eigenfunctions to return in `tapers` and `eigenvalues`. The default value is to return all tapers.

`exitstatus` : output, optional, integer
:   If present, instead of executing a STOP when an error is encountered, the variable exitstatus will be returned describing the error. 0 = No errors; 1 = Improper dimensions of input array; 2 = Improper bounds for input variable; 3 = Error allocating memory; 4 = File IO error.

# Description

`SHReturnTapersMap` will calculate the eigenfunctions (i.e., localization windows) of the space-concentration problem for an arbitrary concentration region specified in `dh_mask` (see Simons et al. (2006) for further details). The input mask `dh_mask` must be sampled according to the Driscoll and Healy (1994) sampling theorem with `n` samples in latitude, and possess a value of 1 inside the concentration region, and 0 elsewhere. `dh_mask` can either possess `n` samples in longitude (`sampling=1`) or `2n` samples in longitude (`sampling=2`). Given the approximate way in which the elements of the space-concentration kernel are calculated (see `ComputeDMap` for details), `sampling=2` should be preferred. The effective spherical harmonic bandwidth (L=N/2-1) of the grid `dh_mask` determines the accuracy of the results, and experience shows that this should be about 4 times larger than `lmax`.

The spherical harmonic coefficients of each window are given in the columns of `tapers`, and the corresponding concentration factors are given in `eigenvaules`. The spherical harmonic coefficients are ordered according to the scheme described in `YilmIndexVector`, which can be converted to matrix form using `SHVectorToCilm`, and the columns of `tapers` are ordered from best to worst concentrated. The localization windows are normalized such that they have unit power. If the optional parameter `ntapers` is specified, then only the `ntapers` largest eigenvalues and corresponding eigenfunctions will be calculated and returned.

# References

Driscoll, J.R. and D.M. Healy, Computing Fourier transforms and convolutions on the 2-sphere, Adv. Appl. Math., 15, 202-250, 1994.

Simons, F.J., F.A. Dahlen, and M.A. Wieczorek, Spatiospectral concentration on a sphere, SIAM Review, 48, 504-536, 2006.

# See also

[computedmap](computedmap.html), [yilmindexvector](yilmindexvector.html), [shvectortocilm](shvectortocilm.html)
