---
title: "SHGravTensor class"
keywords: spherical harmonics software package, spherical harmonic transform, legendre functions, multitaper spectral analysis, fortran, Python, gravity, magnetic field
sidebar: mydoc_sidebar
permalink: python-shgravtensor.html
summary: 
toc: true
---

<style>
table:nth-of-type(n) {
    display:table;
    width:100%;
}
table:nth-of-type(n) th:nth-of-type(2) {
    width:70%;
}
</style>

## Initialization

| Initialization method | Description |
| --------------------- | ----------- |
| `x = SHGravCoeffs.tensor()` | Initialize using an SHGravCoeffs class instance. |

## Class attributes

| Attribute | Description |
| --------- | ----------- |
| `vxx`, `vxy`, `vxz`, `vyx`, `vyy`, `vyz`, `vzx`, `vzy`, `vzz`| The nine component of the gravity tensor. |
| `i0`, `i1`, `i2` | First, second and third invariants of the gravity tensor. |
| `i` | Derived quantity from `i1` and `i2` that is bounded between 0 and 1. |
| `gm` | Gravitational constant time the mass of the body. |
| `a` | Semimajor axis of the reference ellipsoid. |
| `f` | Flattening of the reference ellipsoid, f = (a - b) / a. |
| `lmax` | The maximum spherical harmonic degree resolvable by the grids. |
| `lmax_calc` | The maximum spherical harmonic degree of the gravitational potential used in creating the grids. |
| `nlat`, `nlon` | The number of latitude and longitude bands in the grids. |
| `sampling` | The longitudinal sampling scheme of the grids: either 1 for `nlon` = `nlat` or 2 for `nlon` = 2 * `nlat`. |

## Class methods

| Method | Description |
| ------ | ----------- |
| `plot()` | Plot all 9 components of the gravity tensor. |
| `plot_vxx()` | Plot the Vxx component of the gravity tensor. |
| `plot_vxy()` | Plot the Vxy component of the gravity tensor. |
| `plot_vxz()` | Plot the Vxz component of the gravity tensor. |
| `plot_vyx()` | Plot the Vyx component of the gravity tensor. |
| `plot_vyy()` | Plot the Vyy component of the gravity tensor. |
| `plot_vyz()` | Plot the Vyz component of the gravity tensor. |
| `plot_vzx()` | Plot the Vzx component of the gravity tensor. |
| `plot_vzy()` | Plot the Vzy component of the gravity tensor. |
| `plot_vzz()` | Plot the Vzz component of the gravity tensor. |
| `compute_invar()` | Compute the invariants of the gravity tensor. |
| `plot_i0()` | Plot the first invariant I0 of the gravity tensor. |
| `plot_i1()` | Plot the second invariant I1 of the gravity tensor. |
| `plot_i2()` | Plot the third invariant I2 of the gravity tensor. |
| `plot_i()` | Plot the derived quantity -(I2/2)\*\*2 / (I1/3)\*\*3. |
| `compute_eig()` | Compute the three eigenvalues of the gravity tensor. |
| `plot_eig()` | Plot the three eigenvalues of the gravity tensor. |
| `plot_eig1()` | Plot the first eigenvalue of the gravity tensor. |
| `plot_eig2()` | Plot the second eigenvalue of the gravity tensor. |
| `plot_eig3()` | Plot the third eigenvalue of the gravity tensor. |
| `compute_eigh()` | Compute the horizontal eigenvalues of the gravity tensor. |
| `plot_eigh()` | Plot the two horizontal eigenvalues and the combined maximum absolute eigenvalue of the gravity tensor. |
| `plot_eigh1()` | Plot the first horizontal eigenvalue of the gravity tensor.|
| `plot_eigh2()` | Plot the second horizontal eigenvalue of the gravity tensor. |
| `plot_eighh()` | Plot the combined maximum absolute eigenvalue of the gravity tensor.|
| `copy()` | Return a copy of the class instance. |
| `info()` | Print a summary of the data stored in the SHGravTensor instance. |

