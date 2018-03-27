---
title: SHAdmitCorr (Python)
keywords: spherical harmonics software package, spherical harmonic transform, legendre functions, multitaper spectral analysis, fortran, Python, gravity, magnetic field
sidebar: mydoc_sidebar
permalink: pyshadmitcorr.html
summary:
tags: [python]
toc: false
editdoc: pydoc
---

Calculate the admittance and correlation spectra of two real functions.

## Usage

`admit`, `error`, `corr` = SHAdmitCorr (`gilm`, `tilm`, [`lmax`])

## Returns

`admit` : float, dimension (`lmax`+1)
:   The admittance function, which is equal to `Sgt/Stt`.

`error` : float, dimension (`lmax`+1)
:   The uncertainty of the admittance function, assuming that `gilm` and `tilm` are related by a linear isotropic transfer function, and that the lack of correlation is a result of uncorrelated noise.

`corr` : float, dimension (`lmax`+1)
:   The degree correlation function, which is equal to `Sgt/sqrt(Sgg Stt)`.

## Parameters

`gilm` : float, dimension (2, `lmaxg`+1, `lmaxg`+1)
:   The real spherical harmonic coefficients of the function `G`.

`tilm` : float, dimension (2, `lmaxt`+1, `lmaxt`+1)
:   The real spherical harmonic coefficients of the function `T`.

`lmax` : optional, integer, default = min(`lmaxg`, `lmaxt`)
:   The maximum spherical harmonic degree that will be calculated for the admittance and correlation spectra. This must be less than or equal to the minimum of `lmaxg` and `lmaxt`.

## Description

`SHAdmitCorr` will calculate the admittance, admittance error, and correlation spectra associated with two real functions expressed in real spherical harmonics. The admittance is defined as `Sgt/Stt`, where `Sgt` is the cross-power spectrum of two functions `G` and `T`. The degree-correlation spectrum is defined as `Sgt/sqrt(Sgg Stt)`, which can possess values between -1 and 1. The error of the admittance is calculated assuming that `G` and `T` are related by a linear isotropic transfer function:` Gilm = Ql Tilm + Nilm`, where `N` is noise that is uncorrelated with the topography. It is important to note that the relationship between two fields is often not described by such an isotropic expression.

## See also

[spectrum](spectrum.html), [cross_spectrum](cross_spectrum.html)
