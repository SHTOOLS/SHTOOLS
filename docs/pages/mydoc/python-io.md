---
title: pysh.shio
keywords: spherical harmonics, file io, python, pyshtools
sidebar: mydoc_sidebar
permalink: python-io.html
summary: This module provides routines for reading and writing spherical harmonic coefficients to and from files, and for converting spherical harmonic coefficients to different normalization conventions.
toc: true
folder: mydoc
---

<style>
table:nth-of-type(n) {
    display:table;
    width:100%;
}
table:nth-of-type(n) th:nth-of-type(2) {
    width:75%;
}
</style>

## Spherical harmonic I/O

| Function name | Description |
| ------------- | ----------- |
| [shread](pyshread.html) | Read shtools-formatted spherical harmonic coefficients from a text file. |
| [shwrite](shwrite.html) | Write shtools-formatted spherical harmonic coefficients to a text file. |
| [read_dov](read_dov.html) | Read spherical harmonic coefficients from a text file formatted as [degree, order, value]. |
| [write_dov](write_dov.html) | Write spherical harmonic coefficients to a text file formatted as [degree, order, value]. |
| [read_bshc](read_bshc.html) | Read real spherical harmonic coefficients from a binary bshc-formatted file. |
| [write_bshc](write_bshc.html) | Write real spherical harmonic coefficients to a binary bshc-formatted file. |
| [read_icgem_gfc](read_icgem_gfc.html) | Read real spherical harmonic gravitational potential coefficients and associated errors from an ICGEM GFC formatted file. |
| [write_icgem_gfc](write_icgem_gfc.html) | Write real spherical harmonic gravitational potential coefficients and associated errors to an ICGEM GFC formatted file. |
| [read_igrf](read_igrf.html) | Read IGRF real spherical harmonic coefficients, and return the magnetic potential coefficients for the specified year. |
| [SHRead2](pyshread2.html) | Read spherical harmonic coefficients from a CHAMP or GRACE-like ascii-formatted file. |
| [SHRead2Error](pyshread2error.html) | Read spherical harmonic coefficients and associated errors from a CHAMP or GRACE-like ascii-formatted file. |
| [SHReadJPL](pyshreadjpl.html) | Read spherical harmonic coefficients from a JPL ascii-formatted file. |
| [SHReadJPLError](pyshreadjplerror.html) | Read spherical harmonic coefficients and associated errors from a JPL ascii-formatted file. |

## Spherical harmonic storage

| Function name | Description |
| ------------- | ----------- |
| [SHCilmToCindex](pyshcilmtocindex.html) | Convert a three-dimensional array of complex spherical harmonic coefficients to a two-dimensional indexed array. |
| [SHCindexToCilm](pyshcindextocilm.html) | Convert a two-dimensional indexed array of complex spherical harmonic coefficients to a three-dimensional array. |
| [SHCilmToVector](pyshcilmtovector.html) | Convert a 3-dimensional array of real spherical harmonic coefficients to a 1-dimensional ordered array. |
| [SHVectorToCilm](pyshvectortocilm.html) | Convert a 1-dimensional indexed vector of real spherical harmonic coefficients to a 3-dimensional array. |
| [YilmIndexVector](pyyilmindexvector.html) | Determine the index of a 1-dimensional ordered vector of spherical harmonic coefficients corresponding to *i*, *l*, and *m*.

## Spherical harmonic conversions

| Function name | Description |
| ------------- | ----------- |
| [convert](convert.html) | Convert an array of spherical harmonic coefficients to a different normalization convention. |
| [SHrtoc](pyshrtoc.html) | Convert real spherical harmonics to complex form. |
| [SHctor](pyshctor.html) | Convert complex spherical harmonics to real form. |
