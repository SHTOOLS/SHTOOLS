---
title: SHMagPowerSpectrum (Fortran)
keywords: spherical harmonics software package, spherical harmonic transform, legendre functions, multitaper spectral analysis, fortran, Python, gravity, magnetic field
sidebar: fortran_sidebar
permalink: shmagpowerspectrum.html
summary:
tags: [fortran]
toc: false
editdoc: fdoc
---

Compute the power spectrum of the magnetic field given the Schmidt seminormalized magnetic potential spherical harmonic coefficients.

## Usage

call SHMagPowerSpectrum (`cilm`, `a`, `r`, `lmax`, `spectrum`, `exitstatus`)

## Parameters

`cilm` : input, real(dp), dimension (2, `lmax`+1, `lmax`+1)
:   The Schmidt seminormalized spherical harmonic coefficients of the magnetic potential.

`a` : input, real(dp)
:   The reference radius of the magnetic potential spherical harmonic coefficients.

`r` : input, real(dp)
:   The radius to evaluate the magnetic field.

`lmax` : input, integer
:   The maximum spherical harmonic degree to calculate the power spectrum.

`spectrum` : output, real(dp), dimension (`lmax`+1)
:   The power spectrum of the magnetic field.

`exitstatus` : output, optional, integer
:   If present, instead of executing a STOP when an error is encountered, the variable exitstatus will be returned describing the error. 0 = No errors; 1 = Improper dimensions of input array; 2 = Improper bounds for input variable; 3 = Error allocating memory; 4 = File IO error.

## Description

`SHMagPowerSpectrum` will calculate the power spectrum of the magnetic field at radius `r` given the magnetic potential Schmidt seminormalized spherical harmonic coefficients `cilm` evaluated at radius `a`. For a given degree `l`, this is explicitly calculated as (Lowes 1966):

`S(l) = (l+1) (a/r)**(2l+4) Sum_{m=0}^l [ cilm(1, l+1, m+1)**2 + cilm(2, l+1, m+1)**2 ].`

## Reference

Lowes, F. J., Mean-square values on sphere of spherical harmonic fields, J. Geophys. Res., 71(8), 2179, 1966.

## See also

[shmagpowerl](shmagpowerl.html)
