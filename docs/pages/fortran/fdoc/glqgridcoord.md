---
title: GLQGridCoord (Fortran)
keywords: spherical harmonics software package, spherical harmonic transform, legendre functions, multitaper spectral analysis, fortran, Python, gravity, magnetic field
sidebar: fortran_sidebar
permalink: glqgridcoord.html
summary:
tags: [fortran]
toc: false
editdoc: fdoc
---

Compute the latitude and longitude coordinates used in Gauss-Legendre quadrature grids.

## Usage

call GLQGridCoord (`latglq`, `longlq`, `lmax`, `nlat`, `nlong`, `extend`, `exitstatus`)

## Parameters

`latglq` : output, real(dp), dimension (`lmax`+1)
:   The latitude coordinates of a Gauss-Legendred quadrature grid in degrees.

`longlq` : output, real(dp), dimension (nlong)
:   The longitude coordinates of a Gauss-Legendre quadrature grid in degrees, dimensioned as (2\*`lmax`+1) when `extend` is 0 or (2\*`lmax`+2) when `extend` is 1.

`lmax` : input, integer(int32)
:   The maximum spherical harmonic degree that will be integrated exactly by Gauss-Legendre quadrature.

`nlat` : output, integer(int32)
:   The number of samples in latitude.

`nlong` : output, integer(int32)
:   The number of samples in longitude.

`extend` : input, optional, integer(int32), default = 0
:   If 1, include 360 E longitude.

`exitstatus` : output, optional, integer(int32)
:   If present, instead of executing a STOP when an error is encountered, the variable exitstatus will be returned describing the error. 0 = No errors; 1 = Improper dimensions of input array; 2 = Improper bounds for input variable; 3 = Error allocating memory; 4 = File IO error.

## Description

`GLQGridCoord` will compute the latitude and longitude coordinates that are used in Gauss-Legendre quadrature grids for performing spherical harmonic transforms and reconstructions. The latitudinal nodes correspond to the zeros of the Legendre polynomial of degree `lmax+1`, and the longitudinal nodes are equally spaced with an interval of `360/(2*lmax+1)` degrees.

## See also

[shglq](shglq.html), [shexpandglq](shexpandglq.html), [makegridglq](makegridglq.html), [shexpandglqc](shexpandglqc.html), [makegridglqc](makegridglqc.html), [preglq](preglq.html)
