# PlanetsConstants

SHTOOLS Module containing planetary constants.

# Description

`PlanetsConstants` defines several constants that are useful when working with gravity and topography data of the Earth, Mars, Venus, and Moon. All units are SI. Confer with the source file or the Python documentation for exact values and references.

# Constants

## Fundamental Constants

`Grav_constant`
:   Gravitational Constant.

`pi_constant`
:   Pi.

`mu0_constant` 
:   Magnetic constant.


## The Moon

`R_Moon`
:   Mean radius of the Moon.

`GM_Moon`
:   Gravitational constant times the mass of the Moon.

`Mass_Moon`
:   Mass of the Moon.

`rho_bar_Moon`
:   Average density of the Moon.

`a_Moon`
:   Semi-major axis of the lunar orbit.

`g0_Moon`
:   Gravitational acceleration at R_Moon, not including rotation.

`R0_pot_Moon`
:   Reference radius of the gravitational-potential model.

`Omega_Moon`
:   Angular rotation rate of the Moon.

`C_MOI_Moon`
:   Polar moment of inertia, using R=1738 km.

`I_MOI_Moon`
:   Average moment of inertia, using R=1738 km.

`Gamma_MOI_Moon`
:   Gamma = (B-A)/C.

`Beta_MOI_Moon`
:   Beta = (C-A)/B.

## Mars

`R_Mars`
:   Average radius of Mars.

`GM_Mars` 
:   Gravitational constant times the mass of the Mars.

`Mass_Mars`
:   Mass of Mars.

`rho_bar_Mars`
:   Average density of Mars.

`g0_Mars`
:   Gravitational acceleration at R_Mars, not including rotation.

`R0_pot_Mars`
:   Reference radius of the gravitational-potential model.

`Omega_Mars`
:   Angular rotation rate of Mars.

`f_Mars`
:   Topographic flattening of Mars.

`a_Mars`
:   Semi-major axis radius of Mars.

`b_Mars`
:	 Semi-minor axis radius of Mars.

`W0_mars`
:	Reference potential of Mars.

## Venus

`R_Venus`
:   Average radius of Venus.

`R0_pot_Venus` 
:   Reference radius of the gravitational-potential model.

`GM_Venus`
:   Gravitational constant times the mass of the Venus.

`Mass_Venus`
:   Mass of Venus.

`Omega_Venus`
:   Angular rotation rate of Venus.

`rho_bar_Venus`
:   Average density of Venus.

`g0_Venus`
:   Gravitational acceleration at R_Venus, not including rotation.

## Earth

`GM_Earth`
:   Gravitational constant times the mass of the Earth.

`R0_pot_Earth`
:   Reference radius of the gravitational-potential model.

`Mass_Earth`
:   The mass of the Earth.

`WGS84_a`
:   The semi-major axis of the WGS84 ellipsoid.

`WGS84_b`
:   The semi-minor axis of the WGS84 ellipsoid.

`WGS84_r3`
:   The radius of a sphere of Earth's volume.

`WGS84_f`
:   The flattening of the WGS84 ellipsoid.

`WGS84_gm`
:   The adopted GM of the WGS84 model, which includes the atmosphere.

`WGS84_gma`
:   The GM of the atmosphere adopted by the WGS84 model.

`WGS84_omega`
:   The adopted angular rotation rate of the Earth of the WGS84 model.

`WGS84_U0`
:   The theoretical normal potential associated with the WGS84 model.

## Mercury

`R_Mercury`
:   Average radius of Mercury.

`GM_Mercury`
:   Gravitational constant times the mass of the Mercury.

`Mass_Mercury`
:   Mass of Mercury.

`R0_pot_Mercury`
:   Reference radius used for the spherical harmonic gravity solutions.

`Omega_Mercury_orbit`
:   Angular rotation rate of Mercury about the Sun.

`Omega_Mercury_spin`
:   Angular rotation rate of Mercury.

`rho_bar_Mercury`
:   Average density of Mercury.

`g0_Mercury`
:   Gravitational acceleration at R_Mercury, not including rotation.

# See also

shtools

# Copyright and License

Copyright 2015 Mark Wieczorek.

This is free software. You can distribute and modify it under the terms of the revised BSD license.
