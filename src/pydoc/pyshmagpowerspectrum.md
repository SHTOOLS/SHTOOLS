# SHMagPowerSpectrum

Compute the power spectrum of the magnetic field given the Schmidt seminormalized magnetic potential spherical harmonic coefficients.

# Usage

`spectrum` = pyshtools.SHMagPowerSpectrum (`c`, `a`, [`r`, `lmax`])

# Returns

`spectrum` : float, dimension (`lmax`+1)
:   The power spectrum of the magnetic field.

# Parameters

`c` : float, dimension (2, `lmaxin`+1, `lmaxin`+1)
:   The Schmidt seminormalized spherical harmonic coefficients of the magnetic potential.

`a` : float
:   The reference radius of the magnetic potential spherical harmonic coefficients.

`r` : float, default = `a`
:   The radius to evaluate the magnetic field.

`lmax` : optional, integer, default = `lmaxin`
:   The maximum spherical harmonic degree to calculate the power spectrum.

# Description

`SHMagPowerSpectrum` will calculate the power spectrum of the magnetic field at radius `r` given the magnetic potential Schmidt seminormalized spherical harmonic coefficients `c` evaluated at radius `a`. For a given degree `l`, this is explicitly calculated as:

`S(l) = (l+1) (a/r)**(2l+4) Sum_{m=0}^l [ c[0,l,m]**2 + c[1,l,m]**2 ].`  

# See also

[shmagpowerl](pyshmagpowerl.html)
