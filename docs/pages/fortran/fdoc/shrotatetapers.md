---
title: SHRotateTapers (Fortran)
keywords: spherical harmonics software package, spherical harmonic transform, legendre functions, multitaper spectral analysis, fortran, Python, gravity, magnetic field
sidebar: fortran_sidebar
permalink: shrotatetapers.html
summary:
tags: [fortran]
toc: false
editdoc: fdoc
---

Rotate orthogonal spherical-cap Slepian functions centered at the North pole to a different location.

## Usage

call SHRotateTapers(`tapersrot`, `tapers`, `taper_order`, `lmax`, `nrot`, `x`,
                    `dj`, `exitstatus`)

## Parameters

`tapersrot` : output, real(dp), dimension ((`lmax`+1)**2, `nrot`)
:   An array containing the spherical harmonic coefficients of the rotated spherical-cap functions. Each column corresponds to a single function of which the spherical harmonic coefficients can be unpacked with `SHVectorToCilm`.

`tapers` : input, real(dp), dimension (`lmax`+1, `nrot`)
:   An array containing the eigenfunctions of the spherical-cap concentration problem obtained from `SHReturnTapers`. The functions are listed by columns, ordered from best to worst concentrated.

`taper_order` : input, integer(int32), dimension (`nrot`)
:   The angular order of the non-zero spherical harmonic coefficients in each column of `tapers`.

`lmax` : input, integer(int32)
:   The spherical harmonic bandwidth of the functions.

`nrot` : input, integer(int32)
:   The number of functions to rotate, which must be less than or equal to (`lmax`+1)**2.

`x` : input, real(dp), dimension(3)
:   The three Euler angles, alpha, beta, and gamma, in radians.

`dj` : input, real(dp), dimension (`lmax`+1, `lmax`+1, `lmax`+1)
:   The rotation matrix `dj(pi/2)`, obtained from a call to `djpi2`.

`exitstatus` : output, optional, integer(int32)
:   If present, instead of executing a STOP when an error is encountered, the variable exitstatus will be returned describing the error. 0 = No errors; 1 = Improper dimensions of input array; 2 = Improper bounds for input variable; 3 = Error allocating memory; 4 = File IO error.

## Description

`SHRotateTapers` will rotate a set of orthogonal spherical-cap Slepian functions originally centered at the North pole to a different location according to the three Euler angles in the vector `x`. The original matrix `tapers` is computed by a call to `SHReturnTapers`. Only the first `nrot` tapers are rotated, each of which is returned in a column of the output matrix `tapersrot`. The spherical harmonic coefficients are geodesy 4pi normalized, and each column of `tapersrot` can be unpacked using `SHVectorCilm`. The input rotation matrix `dj` is computed by a call to `djpi2`.

The rotation of a coordinate system or body can be viewed in two complementary ways involving three successive rotations. Both methods have the same initial and final configurations, and the angles listed in both schemes are the same.

`Scheme A:`

(I) Rotation about the z axis by alpha.
(II) Rotation about the new y axis by beta.
(III) Rotation about the new z axis by gamma.

`Scheme B:`

(I) Rotation about the z axis by gamma.
(II) Rotation about the initial y axis by beta.
(III) Rotation about the initial z axis by alpha.

The rotations can further be viewed either as a rotation of the coordinate system or the physical body. For a rotation of the coordinate system without rotation of the physical body, use 

`x(alpha, beta, gamma)`.

For a rotation of the physical body without rotation of the coordinate system, use 

`x(-gamma, -beta, -alpha)`.

To perform the inverse transform of `x(alpha, beta, gamma)`, use `x(-gamma, -beta, -alpha)`.

Note that this routine uses the "y convention", where the second rotation is with respect to the new y axis. If alpha, beta, and gamma were originally defined in terms of the "x convention", where the second rotation was with respect to the new x axis, the Euler angles according to the y convention would be `alpha_y=alpha_x-pi/2`, `beta_x=beta_y`, and `gamma_y=gamma_x+pi/2`.

This routine first converts the real coefficients to complex form using `SHrtoc`. Then the coefficients are converted to indexed form using `SHCilmToCindex`, these are sent to `SHRotateCoef`, the result if converted back to `cilm` complex form using `SHCindexToCilm`, and these are finally converted back to real form using `SHctor`.

## See also

[djpi2](djpi2.html), [shreturntapers](shreturntapers.html), [shctor](shctor.html), [shrtoc](shrtoc.html), [shcilmtocindex](shcilmtocindex.html), [shcindextocilm](shcindextocilm.html), [shvectortocilm](shvectortocilm.html)
