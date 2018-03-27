---
title: SHPowerSpectrumC (Fortran)
keywords: spherical harmonics software package, spherical harmonic transform, legendre functions, multitaper spectral analysis, fortran, Python, gravity, magnetic field
sidebar: mydoc_sidebar
permalink: shpowerspectrumc.html
summary:
tags: [fortran]
toc: false
editdoc: fdoc
---

Compute the power spectrum of a complex function.

## Usage

call SHPowerSpectrumC (`cilm`, `lmax`, `pspectrum`, `exitstatus`)

## Parameters

`cilm` : input, complex\*16, dimension (2, `lmaxin`+1, `lmaxin`+1)
:   The complex function expressed in complex spherical harmonics.

`lmax` : input, integer
:   The maximum spherical harmonic degree of the power spectrum. This must be less than or equal to `lmaxin`.

`pspectrum` : output, real\*8, dimension (`lmax`+1)
:   The power spectrum of the complex function.

`exitstatus` : output, optional, integer
:   If present, instead of executing a STOP when an error is encountered, the variable exitstatus will be returned describing the error. 0 = No errors; 1 = Improper dimensions of input array; 2 = Improper bounds for input variable; 3 = Error allocating memory; 4 = File IO error.

## Description

`SHPowerSpectrumC` will calculate the power spectrum of a complex function expressed in complex 4-pi normalized spherical harmonics. For a given spherical harmonic degree `l`, this is  calculated as:

`pspectrum(l) = Sum_{i=1}^2 Sum_{m=0}^l | cilm(i, l+1, m+1) |**2`.

## See also

[shpowerlc](shpowerlc.html), [shpowerdensitylc](shpowerdensitylc.html), [shcrosspowerlc](shcrosspowerlc.html), [shcrosspowerdensitylc](shcrosspowerdensitylc.html), [shpowerspectrumdensityc](shpowerspectrumdensityc.html), [shcrosspowerspectrumc](shcrosspowerspectrumc.html), [shcrosspowerspectrumdensityc](shcrosspowerspectrumdensityc.html)
