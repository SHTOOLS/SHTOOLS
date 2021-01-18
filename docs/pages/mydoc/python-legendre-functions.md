---
title: "Legendre functions"
keywords: spherical harmonics software package, spherical harmonic transform, legendre functions, multitaper spectral analysis, fortran, Python, gravity, magnetic field
sidebar: mydoc_sidebar
permalink: python-legendre-functions.html
summary: 
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

## Convenience functions

| Function name | Description |
| ------------- | ----------- |
| [legendre](pylegendre.html) | Compute all the associated Legendre functions up to a maximum degree and order. |
| [legendre_lm](pylegendre_lm.html) | Compute the associated Legendre function for specific degrees l and orders m. |

## 4&pi; normalized

| Function name | Description |
| ------------- | ----------- |
| [PlmBar](pyplmbar.html) | Compute all the geodesy-normalized associated Legendre functions. |
| [PlmBar_d1](pyplmbar_d1.html) | Compute all the geodesy-normalized associated Legendre functions and first derivatives. |
| [PlBar](pyplbar.html) | Compute all the geodesy-normalized Legendre polynomials. |
| [PlBar_d1](pyplbar_d1.html) | Compute all the geodesy-normalized Legendre Polynomials and first derivatives. |

## Orthonormalized

| Function name | Description |
| ------------- | ----------- |
| [PlmON](pyplmon.html) | Compute all the orthonormalized associated Legendre functions. |
| [PlmON_d1](pyplmon_d1.html) | Compute all the orthonormalized associated Legendre functions and first derivatives. |
| [PlON](pyplon.html) | Compute all the orthonormalized Legendre polynomials. |
| [PlON_d1](pyplon_d1.html) | Compute all the orthonormalized Legendre polynomials and first derivatives. |

## Schmidt normalized

| Function name | Description |
| ------------- | ----------- |
| [PlmSchmidt](pyplmschmidt.html) | Compute all the Schmidt-normalized associated Legendre functions. |
| [PlmSchmidt_d1](pyplmschmidt_d1.html) | Compute all the Schmidt-normalized associated Legendre functions and first derivatives. |
| [PlSchmidt](pyplschmidt.html) | Compute all the Schmidt-normalized Legendre polynomials. |
| [PlSchmidt_d1](pyplschmidt_d1.html) | Compute all the Schmidt-normalized Legendre polynomials and first derivatives. |

## Unnormalized

| Function name | Description |
| ------------- | ----------- |
| [PLegendreA](pyplegendrea.html) | Compute all the unnormalized associated Legendre functions. |
| [PLegendreA_d1](pyplegendrea_d1.html) | Compute all the unnormalized associated Legendre functions and first derivatives.
| [PLegendre](pyplegendre.html) | Compute all the unnormalized Legendre polynomials. |
| [PLegendre_d1](pyplegendre_d1.html) | Compute all the unnormalized Legendre polynomials and first derivatives. |

## Utilities

| Function name | Description |
| ------------- | ----------- |
| [PlmIndex](pyplmindex.html) | Compute the index of an array of Legendre function corresponding to degree *l* and angular order *m*. |
