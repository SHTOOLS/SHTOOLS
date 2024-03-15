---
title: write_bshc()
keywords: spherical harmonics software package, spherical harmonic transform, legendre functions, multitaper spectral analysis, Python, gravity, magnetic field
sidebar: mydoc_sidebar
permalink: write_bshc.html
summary:
tags: [python]
toc: false
editdoc: pydoc
---

Write real spherical harmonic coefficients to a binary bshc file.

## Usage

write_bshc(filename, coeffs, [lmax])

## Parameters

filename : str or pathlib.Path
:   File name of the binary 'bshc'-formatted spherical harmonic
    coefficients. If filename ends with '.gz' the file will be
    automatically compressed with gzip.

coeffs : ndarray, size(2, lmaxin+1, lmaxin+1)
:   The spherical harmonic coefficients.

lmax : int, optional, default = None
:   The maximum spherical harmonic degree to write to the file. The
    default is to write all coefficients.

## Notes

This function writes real spherical harmonic coefficients to a binary
'bshc'-formatted file as used at Curtin University. The file is composed
solely of 8-byte floats, starting with the minimum and maximum degree,
and followed by the cosine coefficients and then sine coefficients
(with all orders being listed, one degree at a time). For a 100 degree
file, the contents are

0 100
C(0,0), C(1,0), C(1,1), C(2,0), C(2,1), ... C(100,99), C(100,100)
S(0,0), S(1,0), S(1,1), S(2,0), S(2,1), ... S(100,99), S(100,100).

If the filename ends with '.gz', the file will be automatically
compressed using gzip.

