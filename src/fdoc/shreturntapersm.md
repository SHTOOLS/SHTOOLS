# SHReturnTapersM

Calculate the eigenfunctions of the spherical-cap concentration problem for a single angular order.

# Usage

call SHReturnTapersM (`theta0`, `lmax`, `m`, `tapers`, `eigenvalues`, `shannon`, `degrees`, `ntapers`, `exitstatus`)

# Parameters

`theta0` : input, real\*8
:   The angular radius of the spherical cap in radians.

`lmax` : input, integer
:   The spherical harmonic bandwidth of the localization windows.

`m` : input, integer
:   The angular order of the localization windows.

`tapers` : output, real\*8, dimension (`lmax`+1, `lmax`+1)
:   The spherical harmonic coefficients of the `lmax+1` localization windows, arranged in columns. The first and last rows of each column correspond to spherical harmonic degrees 0 and `lmax`, respectively, and the columns are arranged from best to worst concentrated. Only the first `ntapers` columns are non-zero.

`eigenvalues` : output, real\*8, dimension (`lmax`+1)
:   The concentration factors of the localization windows.

`shannon` : output, optional, real\*8
:   The Shannon number, which is the trace of the concentration kernel.

`degrees` : input, integer, optional, dimension (`lmax`+1)
:   List of degrees to use when computing the eigenfunctions. Only those degrees where `degrees(l+1)` is non-zero will be employed.

`ntapers` : output, integer, optional
:   The number of non-zero tapers.

`exitstatus` : output, optional, integer
:   If present, instead of executing a STOP when an error is encountered, the variable exitstatus will be returned describing the error. 0 = No errors; 1 = Improper dimensions of input array; 2 = Improper bounds for input variable; 3 = Error allocating memory; 4 = File IO error.

# Description

`SHReturnTapersM` will calculate the eigenfunctions of the spherical-cap concentration problem for a singule angular order. The 4pi normalized spherical harmonic coefficients of each window are given in the columns of `tapers`, and the corresponding concentration factors are given in `eigenvaules`. The columns of `tapers` are ordered from best to worst concentrated, and the first and last rows of each column correspond to spherical harmonic degrees 0 and `lmax`, respectively. The localization windows are normalized such that they have unit power. If the optional vector `degrees` is specified, then the eigenfunctions will be computed using only those degrees `l` where `degrees(l+1)` is not zero.

When possible, the eigenfunctions are calculated using the kernel of Grunbaum et al. 1982 and the eigenvalues are then calculated by integration using the definition of the space-concentration problem. Use of the Grunbaum et al. kernel is prefered over the space-concentration kernel as the eigenfunctions of the later are unreliable when there are several eigenvalues identical (within machine precision) to either 1 or zero. If, the optional parameter `degrees` is specified, and at least one element is zero for degrees greater or equal to abs(m), then the eigenfunctions and eigenvalues will instead be computed directly using the space-concentration kernel.

# References

Grunbaum, F. A., L. Longhi, and M. Perlstadt, Differential operators commuting with finite convolution integral operators: Some non-abelian examples, SIAM, J. Appl. Math. 42, 941–955, 1982.

Simons, F. J., F. A. Dahlen, and M. A. Wieczorek, Spatiospectral concentration on a sphere, `SIAM Review`, 48, 504-536, 2006.

Wieczorek, M. A. and F. J. Simons, Localized spectral analysis on the sphere, 
`Geophys. J. Int.`, 162, 655-675, 2005.

# See also

[shreturntapers](shreturntapers.html), [computedg82](computedg82.html), [computedm](computedm.html)
