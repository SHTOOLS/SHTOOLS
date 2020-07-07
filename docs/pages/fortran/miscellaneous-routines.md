---
title: "Miscellaneous routines"
keywords: spherical harmonics software package, spherical harmonic transform, legendre functions, multitaper spectral analysis, fortran, Python, gravity, magnetic field
sidebar: fortran_sidebar
permalink: miscellaneous-routines.html
summary: 
toc: false
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

| Routine name | Description |
| ------------ | ----------- |
| [MakeCircleCoord](makecirclecoord.html) | Compute coordinates of a circle placed at a given latitude and longitude. | 
| [MakeEllipseCoord](makeellipsecoord.html) | Compute coordinates of an ellipse placed at a given latitude and longitude. |
| [RandomN](randomn.html) | Return a pseudo uniform random deviate between 0 and 1 using the algorithm of Park and Miller with a Marsaglia shift sequence. |
| [RandomGaussian](randomgaussian.html) | Return a pseudo Gaussian deviate of zero mean and unit variance.|
| [PreGLQ](preglq.html) | Calculate the weights and nodes used in integrating a function by Gauss-Legendre quadrature. |
| [EigValVecSym](eigvalvecsym.html) | Compute the eigenvalues and eigenvectors of a real symmetric matrix. |
| [EigValVecSymTri](eigvalvecsymtri.html) | Compute the eigenvalues and eigenvectors of a real symmetric tridiagonal matrix. |
| [EigValSym](eigvalsym.html) | Compute the eigenvalues of a real symmetric matrix. |
| [Wigner3j](wigner3j.html) | Compute the Wigner-3j symbols for all allowable values of j. |
| [DHaj](dhaj.html) | Compute the latitudinal weights used in the *Driscoll and Healy* (1994) spherical harmonic transform. |
