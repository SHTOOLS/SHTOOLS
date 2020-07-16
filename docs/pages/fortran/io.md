---
title: "Spherical harmonic I/O, storage, and conversions"
keywords: spherical harmonics software package, spherical harmonic transform, legendre functions, multitaper spectral analysis, fortran, Python, gravity, magnetic field
sidebar: fortran_sidebar
permalink: io.html
summary: 
toc: true
folder: fortran
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

| Routine name | Description |
| ------------ | ----------- |
| [SHRead](shread.html) | Read spherical harmonic coefficients from an ascii-formatted file. |
| [SHRead2](shread2.html) | Read spherical harmonic coefficients from a CHAMP or GRACE-like ASCII file. |
| [SHReadJPL](shreadjpl.html) | Read spherical harmonic coefficients from a JPL ASCII file. |

## Spherical harmonic storage

| Routine name | Description |
| ------------ | ----------- |
| [SHCilmToCindex](shcilmtocindex.html) | Convert a three-dimensional array of complex spherical harmonic coefficients to a two-dimensional indexed array. |
| [SHCindexToCilm](shcindextocilm.html) | Convert a two-dimensional indexed array of complex spherical harmonic coefficients to a three-dimensional array. |
| [SHCilmToVector](shcilmtovector.html) | Convert a 3-dimensional array of real spherical harmonic coefficients to a 1-dimensional ordered array. |
| [SHVectorToCilm](shvectortocilm.html) | Convert a 1-dimensional indexed vector of real spherical harmonic coefficients to a 3-dimensional array. |
| [YilmIndexVector](yilmindexvector.html) | Determine the index of a 1D ordered vector of spherical harmonic coefficients corresponding to i, l, and m. |

## Spherical harmonic conversions

| Routine name | Description |
| ------------ | ----------- |
| [SHrtoc](shrtoc.html) | Convert real spherical harmonics to complex form. |
| [SHctor](shctor.html) | Convert complex spherical harmonics to real form. |
