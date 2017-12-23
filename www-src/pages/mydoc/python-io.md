---
title: "Spherical harmonic I/O, storage, and conversions"
keywords: spherical harmonics software package, spherical harmonic transform, legendre functions, multitaper spectral analysis, fortran, Python, gravity, magnetic field
sidebar: mydoc_sidebar
permalink: python-io.html
summary: 
toc: true
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
| [SHRead](pyshread.html) | Read spherical harmonic coefficients from an ascii-formatted file. |
| [SHReadH](pyshreadh.html) | Read spherical harmonic coefficients from an ascii-formatted file with a header line. |
| [SHReadError](pyshreaderror.html) | Read spherical harmonic coefficients and associated errors from an ascii-formatted file. |
| [SHReadErrorH](pyshreaderrorh.html) | Read spherical harmonic coefficients and associated errors from an ascii-formatted file with a header line. |
| [SHRead2](pyshread2.html) | Read spherical harmonic coefficients from a CHAMP or GRACE-like ascii-formatted file. |
| [SHRead2Error](pyshread2error.html) | Read spherical harmonic coefficients and associated errors from a CHAMP or GRACE-like ascii-formatted file. |
| [SHReadJPL](pyshreadjpl.html) | Read spherical harmonic coefficients from a JPL ascii-formatted file. |
| [SHReadJPLError](pyshreadjplerror.html) | Read spherical harmonic coefficients and associated errors from a JPL ascii-formatted file. |
| [read_icgem_gfc](read_icgem_gfc.html) | Read spherical harmonic coefficients from an ICGEM GFC ascii-formatted file. |

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
| [SHrtoc](pyshrtoc.html) | Convert real spherical harmonics to complex form. |
| [SHctor](pyshctor.html) | Convert complex spherical harmonics to real form. |
