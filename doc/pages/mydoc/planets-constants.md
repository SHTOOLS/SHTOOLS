---
title: PlanetsConstants
keywords: spherical harmonics software package, spherical harmonic transform, legendre functions, multitaper spectral analysis, fortran, Python, gravity, magnetic field
sidebar: mydoc_sidebar
permalink: planets-constants.html
summary:
tags: [fortran]
toc: true
---

<style>
table:nth-of-type(n) {
    display:table;
    width:100%;
}
table:nth-of-type(n) th:nth-of-type(2) {
    width:75%;
}
</style>

The `PlanetsConstants`module defines several constants that are useful when working with gravity and topography data of the terrestrial planets. All units are SI. Confer with the Python `info()` method or fortran source file for exact values and references.

## Fundamental constants

| Constant | Description |
| -------- | ----------- |
| `grav_constant` | Gravitational Constant |
| `pi_constant` | Pi |
| `mu0_constant` | Magnetic constant |

## Mercury

| Constant | Description |
| -------- | ----------- |
| `r_mercury` | Average radius of Mercury |
| `gm_mercury` | Gravitational constant times the mass of the Mercury |
| `mass_mercury` | Mass of Mercury |
| `r0_pot_mercury` | Reference radius used for the spherical harmonic gravity solutions |
| `omega_mercury_orbit` | Angular rotation rate of Mercury about the Sun |
| `omega_mercury_spin` | Angular rotation rate of Mercury |
| `rho_bar_mercury` | Average density of Mercury |
| `g0_mercury` | Gravitational acceleration at `r_mercury`, not including rotation |

## Venus

| Constant | Description |
| -------- | ----------- |
| `r_venus` | Average radius of Venus |
| `r0_pot_venus` | Reference radius of the gravitational-potential models of Venus |
| `gm_venus` | Gravitational constant times the mass of the Venus |
| `mass_venus` | Mass of Venus |
| `omega_venus` | Angular rotation rate of Venus |
| `rho_bar_venus` | Average density of Venus |
| `g0_venus` | Gravitational acceleration at `r_venus`, not including rotation |

## Earth

| Constant | Description |
| -------- | ----------- |
| `gm_earth` | Gravitational constant times the mass of Earth |
| `r0_pot_earth` | Reference radius of the terrestrial gravitational-potential models |
| `mass_earth` | The mass of Earth |
| `wgs84_a` | The semi-major axis of the WGS84 ellipsoid |
| `wgs84_b` | The semi-minor axis of the WGS84 ellipsoid |
| `wgs84_r3` | The radius of a sphere of Earth's volume |
| `wgs84_f` | The flattening of the WGS84 ellipsoid |
| `wgs84_gm` | The adopted GM of the WGS84 model, which includes the atmosphere |
| `wgs84_gma` | The GM of the atmosphere adopted by the WGS84 model |
| `wgs84_omega` | The adopted angular rotation rate of the Earth of the WGS84 model |
| `wgs84_u0` | The theoretical normal potential associated with the WGS84 model |

## The Moon

| Constant | Description |
| -------- | ----------- |
| `r_moon` | Mean radius of the Moon |
| `gm_moon` | Gravitational constant times the mass of the Moon |
| `mass_moon` | Mass of the Moon |
| `rho_bar_moon` | Average density of the Moon |
| `a_moon` | Semi-major axis of the lunar orbit |
| `g0_moon` | Mean gravitational acceleration of the Moon at the mean surface radius (`r_moon`), not including rotation |
| `r0_pot_moon` | Reference radius of the gravitational-potential models of the Moon |
| `omega_moon` | Angular rotation rate of the Moon |
| `c_moi_moon` | Polar moment of inertia of the Moon, using R=1738 km |
| `i_moi_moon` | Average moment of inertia of the Moon, using R=1738 km |
| `gamma_moi_moon` | Libration parameter of the Moon, (B-A)/C |
| `beta_moi_moon` | Libration parameter of the Moon, (C-A)/B |

## Mars

| Constant | Description |
| -------- | ----------- |
| `r_mars` |  Average radius of Mars |
| `gm_mars` | Gravitational constant times the mass of Mars |
| `mass_mars` | Mass of Mars |
| `rho_bar_mars` | Average density of Mars |
| `g0_mars` | Gravitational acceleration of Mars at `r_mars`, not including rotation |
| `r0_pot_mars` | Reference radius of the gravitational-potential models of Mars |
| `omega_mars` | Angular rotation rate of Mars |
| `f_mars` | Topographic flattening of Mars |
| `a_mars` | Semi-major axis radius of Mars |
| `b_mars` | Semi-minor axis radius of Mars |
| `w0_mars` | Reference potential of Mars |
