---
title: read_bshc()
keywords: spherical harmonics software package, spherical harmonic transform, legendre functions, multitaper spectral analysis, Python, gravity, magnetic field
sidebar: mydoc_sidebar
permalink: read_bshc.html
summary:
tags: [python]
toc: false
editdoc: pydoc
---

Read real spherical harmonic coefficients from a binary bshc file.

## Usage

coeffs, lmaxout = read_bshc(filename, [lmax])

## Returns

coeffs : ndarray, size(2, lmaxout+1, lmaxout+1)
:   The spherical harmonic coefficients.

lmaxout : int
:   The maximum spherical harmonic degree read from the file.

## Parameters

filename : str or pathlib.Path
:   File name or URL that contains the spherical harmonic coefficients.
    filename will be treated as a URL if it starts with 'http://',
    'https://', or 'ftp://'. If filename ends with '.gz' or '.zip', the
    file will be uncompressed before parsing.

lmax : int, optional, default = None
:   The maximum spherical harmonic degree to read from the file. The
    default is to read the entire file.

## Notes

This function reads real spherical harmonic coefficients from binary
'bshc'-formatted files as used at Curtin University. The file is composed
solely of 8-byte floats, starting with the minimum and maximum degree,
and followed by the cosine coefficients and then sine coefficients
(with all orders being listed, one degree at a time). For a 100 degree
file, the contents are

0 10800
C(0,0), C(1,0), C(1,1), C(2,0), C(2,1), ... C(100,99), C(100,100)
S(0,0), S(1,0), S(1,1), S(2,0), S(2,1), ... S(100,99), S(100,100).

If filename starts with 'http://', 'https://', or 'ftp://', the file will
be treated as a URL. In this case, the file will be downloaded in its
entirety before it is parsed.

If the filename ends with '.gz' or '.zip', the file will be automatically
uncompressed before parsing. For zip files, archives with only a single
file are supported.

