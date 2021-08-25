---
title: write_dov()
keywords: spherical harmonics software package, spherical harmonic transform, legendre functions, multitaper spectral analysis, Python, gravity, magnetic field
sidebar: mydoc_sidebar
permalink: write_dov.html
summary:
tags: [python]
toc: false
editdoc: pydoc
---

Write spherical harmonic coefficients to a text file formatted as
[degree, order, value].

## Usage

write_dov(filename, coeffs, [errors, header, header2, lmax, encoding])

## Parameters

filename : str
:   File name of the 'dov'-formatted spherical harmonic coefficients. If
    filename ends with '.gz' the file will be automatically compressed with
    gzip.

coeffs : ndarray, size(2, lmaxin+1, lmaxin+1)
:   The spherical harmonic coefficients.

errors : ndarray, size(2, lmaxin+1, lmaxin+1), optional, default = None
:   The errors associated with the spherical harmonic coefficients.

header : str, optional default = None
:   A string to be written directly before the spherical harmonic
    coefficients.

header2 : str, optional default = None
:   A second string to be written directly before the spherical harmonic
    coefficients.

lmax : int, optional, default = None
:   The maximum spherical harmonic degree to write to the file.

encoding : str, optional, default = None
:   Encoding of the output file. The default is to use the system default.

## Notes

This function will write spherical harmonic coefficients (and optionally
the errors) to a text file formatted as [degree, order, value]. If header
or header2 are specified, these strings will be written first, directly
before the spherical harmonic coefficients. Both real and complex spherical
harmonic coefficients are supported.

The spherical harmonic coefficients in the file will be formatted as pairs
of lines as

l, m, coeffs[0, l, m]
l, -m, coeffs[1, l, m]

where l and m are the spherical harmonic degree and order, respectively.
If the errors are included, each pair of lines will be formatted as

l, m, coeffs[0, l, m], errors[0, l, m]
l, -m, coeffs[1, l, m], errors[1, l, m]

For each value of increasing l, all the angular orders are listed in
inceasing order, from 0 to l.

If the filename ends with '.gz', the file will be automatically compressed
using gzip.

