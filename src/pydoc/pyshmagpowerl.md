# SHMagPowerL

Compute the power of the magnetic field for a single degree `l` given the Schmidt seminormalized magnetic potential spherical harmonic coefficients.

# Usage

`power` = SHMagPowerL (`c`, `a`, `l`, [`r`])

# Returns

`power` : float
:   The power at degree `l`

# Parameters

`c` : float, dimension (2, l+1, l+1)
:   The Schmidt seminormalized spherical harmonic coefficients of the magnetic potential.

`a` : float
:   The reference radius of the magnetic potential spherical harmonic coefficients.

`l` : integer
:   The spherical harmonic degree for which the power will be calculated.

`r` : optional, float
:   The radius to evaluate the magnetic field.

# Description

`SHMagPowerL` will calculate the power of the magnetic field at radius `r` for a single degree `l` given the magnetic potential Schmidt seminormalized spherical harmonic coefficients `c` evaluated at radius `a`. This is explicitly calculated as:

`S(l) = (l+1) (a/r)**(2l+4) Sum_{m=0}^l [ c[0,l,m]**2 + c[1,l,m)**2 ].`  

# See also

[shmagpowerspectrum](pyshmagpowerspectrum.html)
