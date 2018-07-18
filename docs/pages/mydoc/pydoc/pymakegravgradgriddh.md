---
title: MakeGravGradGridDH (Python)
keywords: spherical harmonics software package, spherical harmonic transform, legendre functions, multitaper spectral analysis, fortran, Python, gravity, magnetic field
sidebar: mydoc_sidebar
permalink: pymakegravgradgriddh.html
summary:
tags: [python]
toc: false
editdoc: pydoc
---

Create 2D cylindrical maps on a flattened ellipsoid of the components of the gravity "gradient" tensor in a local north-oriented reference frame.

## Usage

`vxx`, `vyy`, `vzz`, `vxy`, `vxz`, `vyz` = MakeGravGradGridDH (`cilm`,  `gm`, `r0`, [`a`, `f`, `lmax`, `sampling`, `lmax_calc`])

## Returns

`vxx` : float, dimension (2\*`lmax`+2, `sampling`\*(2*`lmax`+2))
:   A 2D equally sampled (`n` by `n`) or equally spaced (`n` by 2`n`) grid of the `xx` component of the gravity tensor. The first latitudinal band corresponds to 90 N, the latitudinal band for 90 S is not included, and the latitudinal sampling interval is 180/`n` degrees. The first longitudinal band is 0 E, the longitudinal band for 360 E is not included, and the longitudinal sampling interval is 360/`n` for an equally sampled and 180/`n` for an equally spaced grid, respectively.

`vyy` : float, dimension (2\*`lmax`+2, `sampling`\*(2*`lmax`+2))
:   A 2D equally sampled or equally spaced grid of the `yy` component of the gravity tensor.

`vzz` : float, dimension (2\*`lmax`+2, `sampling`\*(2*`lmax`+2))
:   A 2D equally sampled or equally spaced grid of the `zz` component of the gravity tensor.

`vxy` : float, dimension (2\*`lmax`+2, `sampling`\*(2*`lmax`+2))
:   A 2D equally sampled or equally spaced grid of the `xy` component of the gravity tensor.

`vxz` : float, dimension (2\*`lmax`+2, `sampling`\*(2*`lmax`+2))
:   A 2D equally sampled or equally spaced grid of the `xz` component of the gravity tensor.

`vyz` : float, dimension (2\*`lmax`+2, `sampling`\*(2*`lmax`+2))
:   A 2D equally sampled or equally spaced grid of the YZ component of the gravity tensor.

## Parameters

`cilm` : float, dimension (2, `lmaxin`+1, `lmaxin`+1)
:   The real gravitational potential spherical harmonic coefficients. The coefficients `c1lm` and `c2lm` refer to the cosine and sine coefficients, respectively, with `c1lm=cilm[0,l,m]` and `c2lm=cilm[1,l,m]`. 

`gm` : float
:   The gravitational constant multiplied by the mass of the planet.

`r0`: float
:   The reference radius of the spherical harmonic coefficients.

`a` : float
:   The semi-major axis of the flattened ellipsoid on which the field is computed.

`f` : float
:   The flattening of the reference ellipsoid: `f=(R_equator-R_pole)/R_equator`.

`lmax` : optional, integer, default = `lmaxin`
:   The maximum spherical harmonic degree of the coefficients `cilm`. This determines the number of samples of the output grids, `n=2lmax+2`, and the latitudinal sampling interval, `90/(lmax+1)`.

`sampling` : optional, integer, default = 2
:   If 1 the output grids are equally sampled (`n` by `n`). If 2, the grids are equally spaced (`n` by 2`n`).

`lmax_calc` : optional, integer, default = `lmax`
:   The maximum spherical harmonic degree used in evaluating the functions. This must be less than or equal to `lmax`.

## Description

`MakeGravGradGridDH` will create 2-dimensional cylindrical maps from the spherical harmonic coefficients `cilm`, equally sampled (`n` by `n`) or equally spaced (`n` by 2`n`) in latitude and longitude, for six components of the gravity "gradient" tensor (all using geocentric coordinates):

`(Vxx, Vxy, Vxz)`  
`(Vyx, Vyy, Vyz)`  
`(Vzx, Vzy, Vzz)`

The reference frame is north-oriented, where `x` points north, `y` points west, and `z` points upward (all tangent or perpendicular to a sphere of radius r). The gravitational potential is defined as

`V = GM/r Sum_{l=0}^lmax (r0/r)^l Sum_{m=-l}^l C_{lm} Y_{lm}`,

where `r0` is the reference radius of the spherical harmonic coefficients `Clm`, and the gravitational acceleration is

`B = Grad V`.

The gravity tensor is symmetric, and satisfies `Vxx+Vyy+Vzz=0`, though all three diagonal elements are calculated independently in this routine.

The components of the gravity tensor are calculated according to eq. 1 in Petrovskaya and Vershkov (2006), which is based on eq. 3.28 in Reed (1973) (noting that Reed's equations are in terms of latitude and that the `y` axis points east):

`Vzz = Vrr`  
`Vxx = 1/r Vr + 1/r^2 Vtt`  
`Vyy = 1/r Vr + 1/r^2 /tan(t) Vt + 1/r^2 /sin(t)^2 Vpp`  
`Vxy = 1/r^2 /sin(t) Vtp - cos(t)/sin(t)^2 /r^2 Vp`  
`Vxz = 1/r^2 Vt - 1/r Vrt`  
`Vyz = 1/r^2 /sin(t) Vp - 1/r /sin(t) Vrp`

where `r`, `t`, `p` stand for radius, theta, and phi, respectively, and subscripts on `V` denote partial derivatives.

The output grid are in units of s^-2 and are cacluated on a flattened ellipsoid with semi-major axis `a` and flattening `f`. To obtain units of Eotvos (10^-9 s^-2), multiply the output by 10^9. The calculated values should be considered exact only when the radii on the ellipsoid are less than the maximum radius of the planet (the potential coefficients are simply downward continued in the spectral domain).

The default is to calculate grids for use in the Driscoll and Healy (1994) routines that are equally spaced (`n` by `2n`), but this can be changed to calculate equally sampled grids (`n` by `n`) by setting the optional argument `sampling` to 1. The input value of `lmax` determines the number of samples, `n=2lmax+2`, and the latitudinal sampling interval, 90/(`lmax`+1). The first latitudinal band of the grid corresponds to 90 N, the latitudinal band for 90 S is not calculated, and the latitudinal sampling interval is 180/`n` degrees. The first longitudinal band is 0 E, the longitudinal band for 360 E is not calculated, and the longitudinal sampling interval is 360/`n` for equally sampled and 180/`n` for equally spaced grids, respectively.

## References

Reed, G.B., Application of kinematical geodesy for determining
the short wave length components of the gravity field by satellite gradiometry, Ohio State University, Dept. of Geod. Sciences, Rep. No. 201, Columbus, Ohio, 1973.

Driscoll, J.R. and D.M. Healy, Computing Fourier transforms and convolutions on the 2-sphere, Adv. Appl. Math., 15, 202-250, 1994.

Petrovskaya, M.S. and A.N. Vershkov, Non-singular expressions for the gravity gradients in the local north-oriented and orbital reference frames, J. Geod., 80, 117-127, 2006.

## See also

[makegravgriddh](pymakegravgriddh.html), [makegeoidgriddh](pymakegeoidgriddh.html), [makegriddh](pymakegriddh.html)
