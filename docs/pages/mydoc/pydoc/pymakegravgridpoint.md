---
title: MakeGravGridPoint()
keywords: spherical harmonics software package, spherical harmonic transform, legendre functions, multitaper spectral analysis, fortran, Python, gravity, magnetic field
sidebar: mydoc_sidebar
permalink: pymakegravgridpoint.html
summary:
tags: [python]
toc: false
editdoc: pydoc
---

Determine the vector components of the gravitational acceleration at a single point.

## Usage

```python
`value` = MakeGravGridPoint (`cilm`, `gm`, `r0`, `r`, `lat`, `lon`, [`lmax`, `dealloc`])
```

## Returns

`value` : float, dimension (3)
:   Vector components of the gravitational acceleration at (`r`, `lat`, `lon`) in spherical coordinates (r, theta, phi).

## Parameters

`cilm` : float, dimension (2, `lmaxin`+1, `lmaxin`+1)
:   The real 4-pi normalized spherical harmonic coefficients of the gravitational potential. The coefficients `C0lm` and `C1lm` refer to the cosine (`Clm`) and sine (`Slm`) coefficients, respectively, with `Clm=cilm[0,1,m]` and `Slm=cilm[1,l,m]`.

`gm` : float
:   The gravitational constant multiplied by the mass of the planet.

`r0`: float
:   The reference radius of the spherical harmonic coefficients.

`r`: float
:   The radius to evaluate the gravity field.

`lat` : float
:   The latitude of the point in DEGREES.

`lon` : float
:   The longitude of the point in DEGREES.

`lmax` : optional, integer, default = `lmaxin`
:   The maximum spherical harmonic degree used in evaluating the function.

`dealloc` : optional, integer, default = 0
:   0 (default) = Save variables used in the external Legendre function calls. (1) Deallocate this memory at the end of the funcion call.

## Description

`MakeGravGridPoint` will compute the three vector components of the gravitational acceleration at a single point. The input latitude and longitude are in degrees, and the output components of the gravitational acceleration are in spherical coordinates (r, theta, phi). The gravitational potential is given by

`V = GM/r Sum_{l=0}^lmax (r0/r)^l Sum_{m=-l}^l C_{lm} Y_{lm}`,

and the gravitational acceleration is

`B = Grad V`.

The coefficients are referenced to a radius `r0`, and the output accelerations are in m/s^2. To convert m/s^2 to mGals, multiply by 10^5.

