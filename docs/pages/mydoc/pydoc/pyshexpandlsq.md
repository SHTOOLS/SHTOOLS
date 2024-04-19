---
title: SHExpandLSQ()
keywords: spherical harmonics software package, spherical harmonic transform, legendre functions, multitaper spectral analysis, Python, gravity, magnetic field
sidebar: mydoc_sidebar
permalink: pyshexpandlsq.html
summary:
tags: [python]
toc: false
editdoc: pydoc
---

Determine the spherical harmonic coefficients of an irregularly sampled function using a least squares inversion.

## Usage

cilm, chi2 = SHExpandLSQ (d, lat, lon, lmax, [norm,  csphase])

## Returns

cilm : float, dimension (2, lmax+1, lmax+1)
:   The real spherical harmonic coefficients of the function. The coefficients C0lm and C1lm refer to the cosine (Clm) and sine (Slm) coefficients, respectively, with Clm=cilm[0,l,m] and Slm=cilm[1,l,m].

chi2 : float
:   The residual sum of squares misfit.

## Parameters

d : float, dimension (nmax)
:   The value of the function at the coordinates (lat, lon).

lat : float, dimension (nmax)
:   The latitude in degrees corresponding to the value in d.

lon : float, dimension (nmax)
:   The longitude in degrees corresponding to the value in d.

lmax : integer
:   The maximum spherical harmonic degree of the output coefficients cilm.

norm : optional, integer, default = 1
:   1 (default) = Geodesy 4-pi normalized harmonics; 2 = Schmidt semi-normalized harmonics; 3 = unnormalized harmonics; 4 = orthonormal harmonics.

csphase : optional, integer, default = 1
:   1 (default) = do not apply the Condon-Shortley phase factor to the associated Legendre functions; -1 = append the Condon-Shortley phase factor of (-1)^m to the associated Legendre functions.

## Description

SHExpandLSQ will determine the spherical harmonic coefficients of an irregularly sampled function using a least squares inversion. When the number of data points is greater or equal to the number of spherical harmonic coefficients (i.e., nmax>=(lmax+1)**2), the solution of the overdetermined system will be determined. If there are more coefficients than data points, then the solution of the underdetermined system that minimizes the solution norm will be determined. The inversions are performed using the LAPACK routine DGELS.

The employed spherical harmonic normalization and Condon-Shortley phase convention can be set by the optional arguments norm and csphase; if not set, the default is to use geodesy 4-pi normalized harmonics that exclude the Condon-Shortley phase of (-1)^m.
