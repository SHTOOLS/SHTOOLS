# ComputeDM

Compute the space-concentration kernel of a spherical cap.

# Usage

call ComputeDM (`dm`, `lmax`, `m`, `theta0`, `exitstatus`)

# Parameters

`dm` : output, real\*8, dimension (`lmax`+1, `lmax`+1)
:   The space-concentration kernel or angular order `m`.

`lmax` : input, integer
:   The spherical harmonic bandwidth of the windows.

`m` : input, integer
:   The angular order of the concentration problem.

`theta0` : input, real\*8
:   The angular radius of the spherical cap in radians.

`exitstatus` : output, optional, integer
:   If present, instead of executing a STOP when an error is encountered, the variable exitstatus will be returned describing the error. 0 = No errors; 1 = Improper dimensions of input array; 2 = Improper bounds for input variable; 3 = Error allocating memory; 4 = File IO error.

# Description

`ComputeDM` will calculate the space-concentration kernel of angular order `m` for the spherical-cap concentration problem. The eigenfunctions of this matrix correspond to a family of orthogonal windowing functions, and the eigenvalues correspond to the window's concentration factor (i.e., the power of the window within `theta0` divided by the total power of the function). It is assumed that the employed spherical harmonic functions are normalized to the same value for all degrees and angular orders, which is the case for both the geodesy 4-pi and orthonormalized harmonics. This kernel is symmetric and is computed exactly by Gauss-Legendre quadrature.

# References

Simons, F.J., F.A. Dahlen, and M.A. Wieczorek, Spatiospectral concentration on a sphere, SIAM Review, 48, 504-536, 2006.

# See also

[computedg82](computedg82.html), [shreturntapers](shreturntapers.html), [shreturntapersm](shreturntapersm.html)
