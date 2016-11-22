# SHReturnTapersM

Calculate the eigenfunctions of the spherical-cap concentration problem for a single angular order.

# Usage

`tapers`, `eigenvalues` = SHReturnTapersM (`theta0`, `lmax`, `m`)

# Returns

`tapers` : float, dimension (`lmax`+1, `lmax`+1)
:   The spherical harmonic coefficients of the `lmax+1` localization windows, arranged in columns. The first and last rows of each column correspond to spherical harmonic degrees 0 and `lmax`, respectively, and the columns are arranged from best to worst concentrated. 

`eigenvalues` : float, dimension (`lmax`+1)
:   The concentration factors of the localization windows.

# Parameters

`theta0` : float
:   The angular radius of the spherical cap in radians.

`lmax` : integer
:   The spherical harmonic bandwidth of the localization windows.

`m` : integer
:   The angular order of the localization windows.

# Description

`SHReturnTapersM` will calculate the eigenfunctions (i.e., localization windows) of the spherical-cap concentration problem for a singular angular order. The spherical harmonic coefficients of each window are given in the columns of `tapers`, and the corresponding concentration factors are given in `eigenvaules`. The columns of `tapers` are ordered from best to worst concentrated, and the first and last rows of each column correspond to spherical harmonic degrees 0 and `lmax`, respectively. The localization windows are normalized such that they have unit power.

# References

Wieczorek, M. A. and F. J. Simons, Localized spectral analysis on the sphere, 
`Geophys. J. Int.`, 162, 655-675, 2005.

Simons, F.J., F.A. Dahlen, and M.A. Wieczorek, Spatiospectral concentration on a sphere, `SIAM Review`, 48, 504-536, 2006.

# See also

[shreturntapers](pyshreturntapers.html), [computedg82](pycomputedg82.html), [computedm](pycomputedm.html)
