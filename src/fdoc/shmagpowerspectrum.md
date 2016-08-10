# SHMagPowerSpectrum

Compute the power spectrum of the magnetic field given the Schmidt seminormalized magnetic potential spherical harmonic coefficients.

# Usage

call SHMagPowerSpectrum (`c`, `a`, `r`, `lmax`, `spectrum`)

# Parameters

`c` : input, real\*8, dimension (2, `lmax`+1, `lmax`+1)
:   The Schmidt seminormalized spherical harmonic coefficients of the magnetic potential.

`a` : input, real\*8
:   The reference radius of the magnetic potential spherical harmonic coefficients.

`r` : input, real\*8
:   The radius to evaluate the magnetic field.

`lmax` : input, integer
:   The maximum spherical harmonic degree to calculate the power spectrum.

`spectrum` : output, real\*8, dimension (`lmax`+1)
:   The power spectrum of the magnetic field.

# Description

`SHMagPowerSpectrum` will calculate the power spectrum of the magnetic field at radius `r` given the magnetic potential Schmidt seminormalized spherical harmonic coefficients `c` evaluated at radius `a`. For a given degree `l`, this is explicitly calculated as (Lowes 1966):

`S(l) = (l+1) (a/r)**(2l+4) Sum_{m=0}^l [ c(1, l+1, m+1)**2 + c(2, l+1, m+1)**2 ].`

# Reference

Lowes, F. J., Mean-square values on sphere of spherical harmonic fields, J. Geophys. Res., 71(8), 2179, 1966.

# See also

[shmagpowerl](shmagpowerl.html)
