---
title: DownContFilterMC (Python)
keywords: spherical harmonics software package, spherical harmonic transform, legendre functions, multitaper spectral analysis, fortran, Python, gravity, magnetic field
sidebar: mydoc_sidebar
permalink: pydowncontfiltermc.html
summary:
tags: [python]
toc: false
editdoc: pydoc
---

Calculate a minimum-curvature downward continuation filter for a given spherical harmonic degree.

## Usage

`wl` = DownContFilterMC (`l`, `half`, `r`, `d`)

## Returns

`wl` : float
:   The amplitude of the downward continuation filter.

## Parameters

`l` : integer
:   The spherical harmonic degree.

`half` : integer
:   The spherical harmonic degree where the filter is equal to 0.5.

`r` : float
:   The reference radius of the gravitational field.

`d` : float
:   The radius of the surface to downward continue to.

## Description

`DownContFilterMC` will calculate a minimum-curvature downward continuation filter for a given spherical harmonic degree `l`. The input parameters include `half`, which is the degree where the filter is equal to 0.5, and `r` and `d`, which are the reference radius of the gravitational field and the radius of the surface to downward continue to, respectively. 

A simple analytic expression exists for the downward continuation filter, following the methodology of Wieczorek and Phillips (1998), only when taking the first, third, fifth, and so on, derivatives of their equation 17. For this minimum-curvature filter, which corresponds to the second derivative, the form has simply been generalized using the solutions of the odd derivatives. This may or may not turn out to be exact. In any case, the form of this filter is numerically very similar to the Cartesian minimum-curvature filter of Phipps Morgan and Blackman (1993).

## References

Phipps Morgan, J., and D. K. Blackman, Inversion of combined gravity and bathymetry data for crustal structure: A prescription for downward continuation, Earth Planet. Sci. Lett., 119, 167-179, 1993.

Wieczorek, M. A. and R. J. Phillips, Potential anomalies on a sphere: applications to the thickness of the lunar crust, J. Geophys. Res., 103, 1715-1724, 1998.

## See also

[downcontfilterma](pydowncontfilterma.html), [batohilmdh](pybatohilmdh.html) [batohilmrhohdh](pybatohilmrhohdh.html)
