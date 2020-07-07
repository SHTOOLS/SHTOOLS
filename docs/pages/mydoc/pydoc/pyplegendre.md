---
title: PLegendre()
keywords: spherical harmonics software package, spherical harmonic transform, legendre functions, multitaper spectral analysis, fortran, Python, gravity, magnetic field
sidebar: mydoc_sidebar
permalink: pyplegendre.html
summary:
tags: [python]
toc: false
editdoc: pydoc
---

Compute all the unnormalized Legendre polynomials.

## Usage

`p` = PLegendre (`lmax`, `z`)

## Returns

`p` : float, dimension (`lmax`+1)
:   An array of unnormalized Legendre polynomials up to degree `lmax`. Degree `l` corresponds to array index `l`.

## Parameters

`lmax` : integer
:   The maximum degree of the Legendre polynomials to be computed.

`z` : float
:   The argument of the Legendre polynomial.

## Description

`PLegendre` will calculate all of the unnormalized Legendre polynomials up to degree `lmax` for a given argument. These are calculated using a standard three-term recursion formula. The integral of the Legendre polynomials over the interval [-1, 1] is `2/(2l+1)`.
