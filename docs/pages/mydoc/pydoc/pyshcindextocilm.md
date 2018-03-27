---
title: SHCindexToCilm (Python)
keywords: spherical harmonics software package, spherical harmonic transform, legendre functions, multitaper spectral analysis, fortran, Python, gravity, magnetic field
sidebar: mydoc_sidebar
permalink: pyshcindextocilm.html
summary:
tags: [python]
toc: false
editdoc: pydoc
---

Convert a two-dimensional indexed array of spherical harmonic coefficients to a three-dimensional array.

## Usage

`cilm` = SHCindexToCilm (`cindex`, [`lmax`])

## Returns

`cilm` : float, dimension (2, `lmax`+1, `lmax`+1)
:   The input spherical harmonic coefficients. `cilm[0,:,:]` and `cilm[1,:,:]` correspond to either the real and imaginary components, or cosine and sine coefficients, respectively.

## Parameters

`cindex` : float, dimension (2, (`lmaxin`+1)\*(`lmaxin`+2)/2)
:   The indexed spherical harmonic coefficients.

`lmax` : optional, integer, default = `lmaxin`
:   The maximum degree of the output coefficients to convert.

## Description

`SHCindexToCilm` will convert a two-dimensional indexed array of spherical harmonic coefficients to a three-dimensional array of complex spherical harmonic coefficients.  The degree `l` and order `m` corresponds to the index `l*(l+1)/2+m`.

## See also

[shcilmtocindex](pyshcilmtocindex.html), [shcilmtovector](pyshcilmtovector.html), [shvectortocilm](pyshvectortocilm.html)
