---
title: SHMagPowerL (Fortran)
keywords: spherical harmonics software package, spherical harmonic transform, legendre functions, multitaper spectral analysis, fortran, Python, gravity, magnetic field
sidebar: mydoc_sidebar
permalink: shmagpowerl.html
summary:
tags: [fortran]
toc: false
editdoc: fdoc
---

Compute the power of the magnetic field for a single degree `l` given the Schmidt seminormalized magnetic potential spherical harmonic coefficients.

## Usage

`power` = SHMagPowerL (`c`, `a`, `r`, `l`)

## Parameters

`power` : output, real\*8
:   The power at degree `l`

`c` : input, real\*8, dimension (2, l+1, l+1)
:   The Schmidt seminormalized spherical harmonic coefficients of the magnetic potential.

`a` : input, real\*8
:   The reference radius of the magnetic potential spherical harmonic coefficients.

`r` : input, real\*8
:   The radius to evaluate the magnetic field.

`l` : input, integer
:   The spherical harmonic degree for which the power will be calculated.

## Description

`SHMagPowerL` will calculate the power of the magnetic field at radius `r` for a single degree `l` given the magnetic potential Schmidt seminormalized spherical harmonic coefficients `c` evaluated at radius `a`. This is explicitly calculated as:

`S(l) = (l+1) (a/r)**(2l+4) Sum_{m=0}^l [ c(1, l+1, m+1)**2 + c(2, l+1, m+1)**2 ].`

## See also

[shmagpowerspectrum](shmagpowerspectrum.html)
