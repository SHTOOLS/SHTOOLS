---
title: shwrite()
keywords: spherical harmonics software package, spherical harmonic transform, legendre functions, multitaper spectral analysis, fortran, Python, gravity, magnetic field
sidebar: mydoc_sidebar
permalink: shwrite.html
summary:
tags: [python]
toc: false
editdoc: pydoc
---

Write shtools-formatted spherical harmonic coefficients to a text file.

## Usage

shwrite(filename, coeffs, [errors, header, header2, lmax])

## Parameters

**filename : str**
:   File name of the shtools-formatted spherical harmonic coefficients. If filename ends with '.gz' the file will be automatically compressed with gzip.

**coeffs : ndarray, size(2, lmaxin+1, lmaxin+1)**
:   The spherical harmonic coefficients.

**errors : ndarray, size(2, lmaxin+1, lmaxin+1), optional, default = None**
:   The errors associated with the spherical harmonic coefficients.

**header : str, optional default = None**
:   A string to be written directly before the spherical harmonic coefficients.

**header2 : str, optional default = None**
:   A second string to be written directly before the spherical harmonic coefficients.

**lmax : int, optional, default = None**
:   The maximum spherical harmonic degree to write to the file.

## Notes

This function will write spherical harmonic coefficients (and optionally
the errors) to an shtools-formatted text file. If header or header2 are
specified, these strings will be written first, directly before the
spherical harmonic coefficients. Both real and complex spherical harmonic
coefficients are supported.

The spherical harmonic coefficients in the file will be formatted as

l, m, coeffs[0, l, m], coeffs[1, l, m]

where l and m are the spherical harmonic degree and order, respectively.
If the errors are included, each line will be formatted as

l, m, coeffs[0, l, m], coeffs[1, l, m], errors[0, l, m], errors[1, l, m]

For each value of increasing l, all the angular orders are listed in
inceasing order, from 0 to l.

If the filename ends with '.gz', the file will be automatically compressed
using gzip.
    