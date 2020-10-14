---
title: MakeMagGridDH (Fortran)
keywords: spherical harmonics software package, spherical harmonic transform, legendre functions, multitaper spectral analysis, fortran, Python, gravity, magnetic field
sidebar: fortran_sidebar
permalink: makemaggriddh.html
summary:
tags: [fortran]
toc: false
editdoc: fdoc
---

Create 2D cylindrical maps on a flattened ellipsoid of all three vector components of the magnetic field, the magnitude of the magnetic field, and the magnetic potential.

## Usage

call MakeMagGridDH (`cilm`, `lmax`, `r0`, `a`, `f`, `rad`, `theta`, `phi`, `total`, `n`, `sampling`, `lmax_calc`, `pot_grid`, `extend`, `exitstatus`)

## Parameters

`cilm` : input, real(dp), dimension (2, `lmax`+1, `lmax`+1)
:   The real Schmidt semi-normalized spherical harmonic coefficients to be expanded in the space domain. The coefficients `C1lm` and `C2lm` refer to the cosine (`Clm`) and sine (`Slm`) coefficients, respectively, with `Clm=cilm(1,l+1,m+1)` and `Slm=cilm(2,l+1,m+1)`. Alternatively, `C1lm` and `C2lm` correspond to the positive and negative order coefficients, respectively. The coefficients are assumed to have units of nT.

`lmax` : input, integer(int32)
:   The maximum spherical harmonic degree of the coefficients `cilm`. This determines the number of samples of the output grids, `n=2*lmax+2`, and the latitudinal sampling interval, `90/(lmax+1)`.

`r0` : input, real(dp)
:   The reference radius of the spherical harmonic coefficients.

`a` : input, real(dp)
:   The semi-major axis of the flattened ellipsoid on which the field is computed.

`f` : input, real(dp)
:   The flattening of the reference ellipsoid: i.e., `F=(R_equator-R_pole)/R_equator`.

`rad` : output, real(dp), dimension(nlat, nlong)
:   A 2D map of radial component of the magnetic field that conforms to the sampling theorem of Driscoll and Healy (1994). If `sampling` is 1, the grid is equally sampled and is dimensioned as (`n` by `n`), where `n` is `2lmax+2`. If sampling is 2, the grid is equally spaced and is dimensioned as (`n` by 2`n`). The first latitudinal band of the grid corresponds to 90 N, the latitudinal sampling interval is 180/`n` degrees, and the default behavior is to exclude the latitudinal band for 90 S. The first longitudinal band of the grid is 0 E, by default the longitudinal band for 360 E is not included, and the longitudinal sampling interval is 360/`n` for an equally sampled and 180/`n` for an equally spaced grid, respectively. If `extend` is 1, the longitudinal band for 360 E and the latitudinal band for 90 S will be included, which increases each of the dimensions of the grid by 1.

`theta` : output, real(dp), dimension(nlat, nlong)
:   A 2D equally sampled or equally spaced grid of the theta component of the magnetic field.

`phi` : output, real(dp), dimension(nlat, nlong)
:   A 2D equally sampled or equally spaced grid of the phi component of the magnetic field.

`total` : output, real(dp), dimension(nlat, nlong)
:   A 2D equally sampled or equally spaced grid of the total magnetic field strength.

`n` : output, integer(int32)
:   The number of samples in latitude of the output grids. This is equal to `2lmax+2`.

`sampling` : optional, input, integer(int32), default = 1
:   If 1 (default) the output grids are equally sampled (`n` by `n`). If 2, the grids are equally spaced (`n` by 2`n`).

`lmaxcalc` : optional, input, integer(int32), default = `lmax`
:   The maximum spherical harmonic degree used in evaluating the functions. This must be less than or equal to `lmax`.

`potgrid` : output, real(dp), dimension(nlat, nlong)
:   A 2D equally sampled or equaly spaced grid of the magnetic potential.

`extend` : input, optional, integer(int32), default = 0
:   If 1, compute the longitudinal band for 360 E and the latitudinal band for 90 S. This increases each of the dimensions of the grids by 1.

`exitstatus` : output, optional, integer(int32)
:   If present, instead of executing a STOP when an error is encountered, the variable exitstatus will be returned describing the error. 0 = No errors; 1 = Improper dimensions of input array; 2 = Improper bounds for input variable; 3 = Error allocating memory; 4 = File IO error.

## Description

`MakeMagGridDH` will create 2-dimensional cylindrical maps from the spherical harmonic coefficients `cilm` of all three components of the magnetic field, the total field strength, and the magnetic potential. The magnetic potential is given by

`V = R0 Sum_{l=1}^LMAX (R0/r)^{l+1} Sum_{m=-l}^l C_{lm} Y_{lm}`

and the magnetic field is

`B = - Grad V`.

The coefficients are referenced to a radius `r0`, and the function is computed on a flattened ellipsoid with semi-major axis `a` (i.e., the mean equatorial radius) and flattening `f`.

The default is to use an input grid that is equally sampled (`n` by `n`), but this can be changed to use an equally spaced grid (`n` by 2`n`) by the optional argument `sampling`. The redundant longitudinal band for 360 E and the latitudinal band for 90 S are excluded by default, but these can be computed by specifying the optional argument `extend`.

## Reference

Driscoll, J.R. and D.M. Healy, Computing Fourier transforms and convolutions on the 2-sphere, Adv. Appl. Math., 15, 202-250, 1994.
