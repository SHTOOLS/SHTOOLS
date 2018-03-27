---
title: Wigner3j (Python)
keywords: spherical harmonics software package, spherical harmonic transform, legendre functions, multitaper spectral analysis, fortran, Python, gravity, magnetic field
sidebar: mydoc_sidebar
permalink: pywigner3j.html
summary:
tags: [python]
toc: false
editdoc: pydoc
---

Compute the Wigner-3j symbols for all allowable values of J.

## Usage

`w3j`, `jmin`, `jmax` = Wigner3j (`j2`, `j3`, `m1`, `m2`, `m3`)

## Returns

`w3j` : float, dimension (`j2`+`j3`+1)
:   An array of the Wigner-3j symbols evaluated for all allowable values of `j`. The minimum and maximum values of `j` are given by `jmin` and `jmax`.

`jmin` : integer
:   The minimum value of `j` in the array `w3j`. This corresponds to the first element of `w3j`.

`jmax` : integer
:   The maximum value of `j` in the array `w3j`. This corresponds to the last non-zero element of `w3j`.

## Parameters

`j2` : integer
:   A positive integer.

`j3` : integer
:   A positive integer.

`m1` : integer
:   An integer.

`m2` : integer
:   An integer.

`m3` : integer
:   An integer.

## Description

`Wigner3j` will calculate the Wigner 3J symbols

`/ j  j2 j3 \`  
`\ m1 m2 m3 /`

for all allowable values of `j`. The returned values in the array `w3j` are calculated only for the limits

`jmin = max(|j2-j3|, |m1|)`  and

`jmax = j2 + j3`.

To be non-zero, `m1+m2+m3` must equal 0. It is assumed that all `j`s and `m`s are integers. Returned values have a relative error less than ~1.d-8 when `j2` and `j3` are less than about 100 (see below). In practice, this routine is probably usable up to about 165.

The employed algorithm is based upon the stable non-linear recurrence relations of Luscombe and Luban (1998) for the "non classical" regions near `jmin and jmax`. The direction of the iteration starts from low values of `j` to high values, but when `abs(w3j(j+2)/w3j(j))` is less than one, the iteration will restart from high to low values. For the classical region, the standard three term recursion relationship is used (e.g., Schulten and Gordon 1975). As this three term recursion can lead to overflows, the values are rescaled by a factor "scalef" whenever the absolute value of the 3j coefficient becomes greater than unity.  More efficient algorithms probably exist for specific cases (for instance, when all Ms are zero).

The results of this routine have been verified against the same routine run in quadruple precision. For 1.e7 acceptable random values of `j2`, `j3`, `m2`, and `m3` between -200 and 200, the relative error was calculated only for those 3j coefficients that had an absolute value greater than 1.d-17 (values smaller than this are for all practical purposed zero, and can be heavily affected by machine roundoff errors or underflow). 853 combinations of parameters were found to have relative errors greater than 1.d-8. Here I list the minimum value of `max(j2,j3)` for different ranges of error, as well as the number of times this error occurred:

`max(j2,j3) = 103: 1.d-7 < error <=1.d-8 ; Number of occurrences = 483`  
`max(j2,j3) = 116: 1.d-6 < error <= 1.d-7 ; Number of occurrences = 240`  
`max(j2,j3) = 165: 1.d-5 < error <= 1.d-6 ; Number of occurrences = 93`  
`max(j2,j3) = 167: 1.d-4 < error <= 1.d-5 ; Number of occurrences = 36`

Many times, the large relative errors occur when the 3j coefficient changes sign and is very close to zero (i.e., adjacent values are about 1.e7 times greater in magnitude). Thus, if one does not need to know highly accurate values of the 3j coefficients when they are almost zero (i.e., ~1.e-10) then this routine is probably usable up to about 160.

These results have also been verified for parameter values less than 100 using a code based on the algorith of de Blanc (1987), which was originally coded by Olav van Genabeek, and modified by M. Fang. (This code was run in quadruple precision and only calculates one coefficient for each call.) Maximum relative errors between the two routines were less than 1.d-8 for a large number of values (again, only 3j coefficients greater than 1.d-17 were considered here).

## References

Luscombe, J. J., and M. Luban, Simplified recursive algorithm for Wigner 3j and 6j symbols, Phys. Rev. E, 57, 7274-7277, 1998.

Schulten, K., and R. G. Gordon, Exact recursive evaluation of 3j-coefficients
and 6j-coefficients for quantum-mechanical coupling of angular momenta, J. Math. Phys., 16, 1961-1970, 1975.
