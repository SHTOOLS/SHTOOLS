---
title: SHCrossPowerSpectrumDensityC (Fortran)
keywords: spherical harmonics software package, spherical harmonic transform, legendre functions, multitaper spectral analysis, fortran, Python, gravity, magnetic field
sidebar: mydoc_sidebar
permalink: shcrosspowerspectrumdensityc.html
summary:
tags: [fortran]
toc: false
editdoc: fdoc
---

Compute the cross-power spectral density of two complex functions.

## Usage

call SHCrossPowerSpectrumDensityC (`cilm1`, `cilm2`, `lmax`, `cspectrum`, `exitstatus`)

## Parameters

`cilm1` : input, complex\*16, dimension (2, `lmaxin1`+1, `lmaxin1`+1)
:   The complex spherical harmonics of the first complex function.

`cilm2` : input, complex\*16, dimension (2, `lmaxin2`+1, `lmaxin2`+1)
:   The complex spherical harmonics of the first complex function.

`lmax` : input, integer
:   The maximum spherical harmonic degree of the cross power spectral density. This must be less than or equal to the minimum of `lmaxin1` and `lmaxin2`.

`cspectrum` : output, complex\*16, dimension (`lmax`+1)
:   The cross-power spectral density of the two complex functions.

`exitstatus` : output, optional, integer
:   If present, instead of executing a STOP when an error is encountered, the variable exitstatus will be returned describing the error. 0 = No errors; 1 = Improper dimensions of input array; 2 = Improper bounds for input variable; 3 = Error allocating memory; 4 = File IO error.

## Description

`SHCrossPowerSpectrumDensityC` will calculate the cross-power spectral density of two complex functions expressed in complex 4-pi normalized spherical harmonics. For a given spherical harmonic degree `l`, this is calculated as:

`cspectrum(l) = Sum_{i=1}^2 Sum_{m=0}^l cilm1(i, l+1, m+1) * conjg[cilm2(i, l+1, m+1)] / (2l + 1)`.

## See also

[shpowerlc](shpowerlc.html), [shpowerdensitylc](shpowerdensitylc.html), [shcrosspowerlc](shcrosspowerlc.html), [shcrosspowerdensitylc](shcrosspowerdensitylc.html), [shpowerspectrumc](shpowerspectrumc.html), [shpowerspectrumdensityc](shpowerspectrumdensityc.html), [shcrosspowerspectrumc](shcrosspowerspectrumc.html)
