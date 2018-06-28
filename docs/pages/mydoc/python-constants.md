---
title: Constants
keywords: spherical harmonics software package, spherical harmonic transform, legendre functions, multitaper spectral analysis, fortran, Python, gravity, magnetic field
sidebar: mydoc_sidebar
permalink: python-constants.html
summary:
tags: [python]
toc: 
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

This subpackage defines several constants that are used in analyzing gravity,
topography, and magnetic field data of the terrestrial planets. Each object is
an astropy `Constant` that possesses the attributes name, value, unit,
uncertainty, and reference.

## Fundamental constants

| Constant | Description |
| -------- | ----------- |
| `G` | Gravitational Constant |
| `mu0` | Magnetic constant |

## Mercury

| Constant | Description |
| -------- | ----------- |
| `r_mercury` | Average radius of Mercury |
| `gm_mercury` | Gravitational constant times the mass of the Mercury |
| `mass_mercury` | Mass of Mercury |
| `omega_orbit_mercury` | Angular rotation rate of Mercury about the Sun |
| `omega_mercury` | Angular rotation rate of Mercury |
| `density_mercury` | Average density of Mercury |
| `g0_mercury` | Gravitational acceleration at `r_mercury`, not including rotation |

## Venus

| Constant | Description |
| -------- | ----------- |
| `r_venus` | Average radius of Venus |
| `gm_venus` | Gravitational constant times the mass of the Venus |
| `mass_venus` | Mass of Venus |
| `omega_venus` | Angular rotation rate of Venus |
| `density_venus` | Average density of Venus |
| `g0_venus` | Gravitational acceleration at `r_venus`, not including rotation |

## Earth

| Constant | Description |
| -------- | ----------- |
| `gm_egm2008` | Gravitational constant times the mass of Earth from the EGM2008 model |
| `mass_egm2008` | The mass of Earth from the EGM2008 model |
| `a_wgs84` | The semi-major axis of the WGS84 ellipsoid |
| `b_wgs84` | The semi-minor axis of the WGS84 ellipsoid |
| `r3_wgs84` | The radius of a sphere of Earth's volume |
| `f_wgs84` | The flattening of the WGS84 ellipsoid |
| `gm_wgs84` | The adopted GM of the WGS84 model, which includes the atmosphere |
| `mass_wgs84` | The mass of Earth from the WGS84 model, which includes the atmosphere |
| `gma_wgs84` | The GM of the atmosphere adopted by the WGS84 model |
| `omega_wgs84` | The adopted angular rotation rate of the Earth of the WGS84 model |
| `u0_wgs84` | The theoretical normal potential associated with the WGS84 model |

## The Moon

| Constant | Description |
| -------- | ----------- |
| `r_moon` | Mean radius of the Moon |
| `gm_moon` | Gravitational constant times the mass of the Moon |
| `mass_moon` | Mass of the Moon |
| `density_moon` | Average density of the Moon |
| `a_orbit_moon` | Semi-major axis of the lunar orbit |
| `g0_moon` | Mean gravitational acceleration of the Moon at the mean surface radius (`r_moon`), not including rotation |
| `omega_moon` | Angular rotation rate of the Moon |
| `i_solid_moon` | Average moment of inertial of the solid portion of the Moon |
| `gamma_moon` | Libration parameter of the Moon, (B-A)/C |
| `beta_moon` | Libration parameter of the Moon, (C-A)/B |

## Mars

| Constant | Description |
| -------- | ----------- |
| `r_mars` |  Average radius of Mars |
| `gm_mars` | Gravitational constant times the mass of Mars |
| `mass_mars` | Mass of Mars |
| `density_mars` | Average density of Mars |
| `g0_mars` | Gravitational acceleration of Mars at `r_mars`, not including rotation |
| `omega_mars` | Angular rotation rate of Mars |
| `f_mars` | Topographic flattening of Mars |
| `a_mars` | Semi-major axis radius of Mars |
| `b_mars` | Semi-minor axis radius of Mars |
| `u0_mars` | Reference potential of Mars |
