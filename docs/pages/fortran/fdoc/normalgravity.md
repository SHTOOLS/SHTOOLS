---
title: NormalGravity (Fortran)
keywords: spherical harmonics software package, spherical harmonic transform, legendre functions, multitaper spectral analysis, fortran, Python, gravity, magnetic field
sidebar: fortran_sidebar
permalink: normalgravity.html
summary:
tags: [fortran]
toc: false
editdoc: fdoc
---

Calculate the normal gravity on an ellipsoid or sphere using geocentric coordinates.

## Usage

`value` = NormalGravity (`geocentriclat`, `gm`, `omega`, `a`, `b`)

## Parameters

`value` : output, real(dp)
:   The normal gravity in SI units.

`geocentriclat`: input, real(dp)
:   Geocentric latitude in degrees.

`gm` : input, real(dp)
:   The gravitational constant multiplied by the mass of the planet.

`omega` : input, real(dp)
:   The angular rotation rate of the planet.

`a` : input, real(dp)
:   The semi-major axis of the ellipsoid on which the normal gravity is computed.

`b` : input, real(dp)
:   The semi-minor axis of the ellipsoid on which the normal gravity is computed.

## Description

`NormalGravity` will compute the magnitude of the total gravity (gravitation and centrifugal) on the surface of a rotating ellipsoid (in m/s^2).  The latitude is input in geocentric coordinates in degrees.

For a rotating ellipsoid, the surface corresponds to a constant potential, and the gravity vector is normal to the surface. The normal gravity is computed using Somigliana's formula, as taken from Physical Geodesy (Hofmann-Wellenhof and Moritz, 2nd ed., sections 2.7 and 2.8). In this routine, the geodetic latitude is computed internally from the input geocentric latitude.

For the case of a sphere (a is equal to b), the normal gravity is defined as the magnitude of the sum of the normal gravitation (GM/r^2) and the centrifugal gravity. For a rotating sphere, the surface does not correspond to a constant potential and the gravity vector is not normal to the surface.

## References

Hofmann-Wellenhof B, and H. Moritz, "Physical Geodesy," second edition, Springer, Wien, 403 pp., 2006.

## See also

[makegravgriddh](makegravgriddh.html), [makegeoidgrid](makegeoidgrid.html)
