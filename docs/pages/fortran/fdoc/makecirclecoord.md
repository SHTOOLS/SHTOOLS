---
title: MakeCircleCoord (Fortran)
keywords: spherical harmonics software package, spherical harmonic transform, legendre functions, multitaper spectral analysis, fortran, Python, gravity, magnetic field
sidebar: fortran_sidebar
permalink: makecirclecoord.html
summary:
tags: [fortran]
toc: false
editdoc: fdoc
---

Compute the coordinates of a circle placed at a given latitude and longitude.

## Usage

call MakeCircleCoord (`coord`, `lat`, `lon`, `theta0`, `cinterval`, `cnum`, `exitstatus`)

## Parameters

`coord` : output, real(dp), dimension(360/`cinterval`, 2)
:   The latitude (:,1) and longitude (:,2) coordinates of the circle in degrees. If not specified, `cintervaL` is assumed to 1.

`lat` : input, real(dp)
:   The latitude of the center of the circle in degrees.

`lon` : input, real(dp)
:   The longitude of the center of the circle in degrees.

`theta0` : input, real(dp)
:   The angular radius of the circle in degrees.

`cinterval` : optional, input, real(dp), default = 1
:   Angular spacing in degrees of the output latitude and longitude points. If not present, the default is 1.

`cnum` : optional, output, integer
:   Number of elements in the output arrays.

`exitstatus` : output, optional, integer
:   If present, instead of executing a STOP when an error is encountered, the variable exitstatus will be returned describing the error. 0 = No errors; 1 = Improper dimensions of input array; 2 = Improper bounds for input variable; 3 = Error allocating memory; 4 = File IO error.

## Description

`MakeCircleCoord` will calculate the latitude and longitude coordinates of a circle of angular radius `theta0` placed on a sphere at position (`lat`, `lon`). This is useful for plotting circles on geographic maps. The first index in the output vectors corresponds to the point directly north of the cirlce origin, and subsequent points are arranged in a clockwise manner. Input and output units are in degrees.

## See also

[makeellipsecoord](makeellipsecoord.html)
