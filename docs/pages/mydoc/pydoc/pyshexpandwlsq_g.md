---
title: SHExpandWLSQ_G()
keywords: spherical harmonics software package, spherical harmonic transform, legendre functions, multitaper spectral analysis, Python, gravity, magnetic field
sidebar: mydoc_sidebar
permalink: pyshexpandwlsq_g.html
summary:
tags: [python]
toc: false
editdoc: pydoc
---

Determine the spherical harmonic coefficients of an irregularly sampled function using a weighted least squares inversion with a precomputed data kernel matrix.

## Usage

cilm, chi2 = SHExpandWLSQ_G (d, w, lat, lon, lmax, g, [norm,  csphase])

## Returns

cilm : float, dimension (2, lmax+1, lmax+1)
:   The real spherical harmonic coefficients of the function. The coefficients C0lm and C1lm refer to the cosine (Clm) and sine (Slm) coefficients, respectively, with Clm=cilm[0,l,m] and Slm=cilm[1,l,m].

chi2 : float
:   The residual weighted sum of squares misfit.

## Parameters

d : float, dimension (nmax)
:   The value of the function at the coordinates (lat, lon).

w : float, dimension (nmax)
:   The weights used in the weighted least squares inversion.

lat : float, dimension (nmax)
:   The latitude in degrees corresponding to the value in d.

lon : float, dimension (nmax)
:   The longitude in degrees corresponding to the value in d.

lmax : integer
:   The maximum spherical harmonic degree of the output coefficients cilm.

g : float, dimension(nmax, (lmax+1)**2)
:   The precomputed data kernel matrix G obtained from LSQ_G.

norm : optional, integer, default = 1
:   1 (default) = Geodesy 4-pi normalized harmonics; 2 = Schmidt semi-normalized harmonics; 3 = unnormalized harmonics; 4 = orthonormal harmonics.

csphase : optional, integer, default = 1
:   1 (default) = do not apply the Condon-Shortley phase factor to the associated Legendre functions; -1 = append the Condon-Shortley phase factor of (-1)^m to the associated Legendre functions.

## Description

SHExpandWLSQ_G will determine the spherical harmonic coefficients of an irregularly sampled function using a weighted least squares inversion with a precomputed data kernel matrix G. The matrix G can be obtained from the routine LSQ_G. The weights should be set equal to the inverse of the data variance, and it is assumed explicitly that each measurement is statistically independent (i.e., the weighting matrix is diagonal). The weighted least squares inversion must be overdetermined (i.e., nmax>(lmax+1)**2), and the inversion is performed using the LAPACK routine DGELS after scaling the data vector and inversion matrix.

The employed spherical harmonic normalization and Condon-Shortley phase convention can be set by the optional arguments norm and csphase; if not set, the default is to use geodesy 4-pi normalized harmonics that exclude the Condon-Shortley phase of (-1)^m.
