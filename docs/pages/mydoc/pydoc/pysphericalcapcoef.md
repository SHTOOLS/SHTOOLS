---
title: SphericalCapCoef (Python)
keywords: spherical harmonics software package, spherical harmonic transform, legendre functions, multitaper spectral analysis, fortran, Python, gravity, magnetic field
sidebar: mydoc_sidebar
permalink: pysphericalcapcoef.html
summary:
tags: [python]
toc: false
editdoc: pydoc
---

Calculate the spherical harmonic coefficients of a spherical cap.

## Usage

`coef` = SphericalCapCoef (`theta`, `lmax`)

## Returns

`coef` : float, dimension(`lmax`+1)
:   The zonal spherical harmonic coefficients of a spherical cap centered over the north pole.

## Parameters

`theta` : float
:   The angular radius of the spherical cap in radians.

`lmax` : integer
:   The maximum spherical harmonic degree to calculate the spherical harmonic coefficients.

## Description

`SphericalCapCoef` will calculate the spherical harmonic coefficients of a spherical cap centered over the north pole. The zonal coefficients, returned in the array `coef`, are normalized such that the degree-0 term is 1, and are to be used with either the geodesy 4-pi normalized or orthonormalized spherical harmonics.
