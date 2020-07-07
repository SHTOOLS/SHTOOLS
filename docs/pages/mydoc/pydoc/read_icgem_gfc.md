---
title: read_icgem_gfc()
keywords: spherical harmonics software package, spherical harmonic transform, legendre functions, multitaper spectral analysis, fortran, Python, gravity, magnetic field
sidebar: mydoc_sidebar
permalink: read_icgem_gfc.html
summary:
tags: [python]
toc: false
editdoc: pydoc
---

Read real spherical harmonic gravity coefficients from an ICGEM formatted
file.

## Returns

**cilm : ndarray, size (2, lmax + 1, lmax + 1)**
:   Array of '4pi' normalized spherical harmonic coefficients for the given epoch.

**gm : float**
:   Gravitational constant of the model, in m\*\*3/s\*\*2.

**r0 : float**
:   Reference radius of the model, in meters.

**errors : ndarray, optional, shape (2, lmax + 1, lmax + 1)**
:   Array of the spherical harmonic error coefficients for the given epoch.

## Parameters

**filename : str**
:   The filename containing the spherical harmonic ICGEM-formatted coefficients. filename will be treated as a URL if it starts with 'http://', 'https://', or 'ftp://'. If filename ends with '.gz' or '.zip' (or if the path contains '/zip/'), the file will be uncompressed before parsing.

**errors : str, optional, default = None**
:   Which errors to read. Can be 'unknown', 'calibrated', 'formal' or None.

**lmax : int, optional, default = None**
:   Maximum degree to read from the file. If lmax is None, less than 0, or greater than lmax_model, the maximum degree of the model will be used.

**epoch : str or float, optional, default = None**
:   The epoch time to calculate time-variable coefficients in YYYYMMDD.DD format. If None then the reference epoch t0 of the model will be used. If the format of the file is 'icgem2.0' then the epoch must be specified.

**encoding : str, optional**
:   Encoding of the input file. Try to use 'iso-8859-1' if the default (UTF-8) fails.
    