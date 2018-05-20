---
title: "Gravity and Magnetics routines"
keywords: spherical harmonics software package, spherical harmonic transform, legendre functions, multitaper spectral analysis, fortran, Python, gravity, magnetic field
sidebar: mydoc_sidebar
permalink: gravity-magnetics.html
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

## Gravity routines

| Routine name | Description |
| ------------ | ----------- |
| [MakeGravGridDH](makegravgriddh.html) | Create 2D cylindrical maps on a flattened and rotating ellipsoid of all three components of the gravity field, the gravity disturbance, and the gravitational potential. |
| [MakeGravGradGridDH](makegravgradgriddh.html) | Calculate the components of the gravity "gradient" tensor on a flattened ellipsoid. |
| [MakeGeoidGrid](makegeoidgrid.html) | Create a global map of the geoid. |
| [CilmPlus](cilmplus.html) | Calculate the gravitational potential exterior to relief along a spherical interface using the finite-amplitude algorithm of *Wieczorek and Phillips* (1998). |
| [CilmMinus](cilmminus.html) | Calculate the gravitational potential interior to relief along a spherical interface using the finite-amplitude algorithm of *Wieczorek and Phillips* (1998). |
| [CilmPlusRhoH](cilmplusrhoh.html) | Calculate the gravitational potential exterior to relief along a spherical interface with laterally varying density using the finite amplitude algorithm of *Wieczorek* (2007).|
| [CilmMinusRhoH](cilmminusrhoh.html) | Calculate the gravitational potential interior to relief along a spherical interface with laterally varying density using the finite amplitude algorithm of *Wieczorek* (2007).|
| [BAtoHilm](batohilm.html) | Calculate iteratively the relief along an interface with constant density contrast that corresponds to a given Bouguer anomaly using the algorithm of *Wieczorek and Phillips* (1998). |
| [BAtoHilmRhoH](batohilmrhoh.html) | Iteratively calculate the relief along an interface with laterally varying density contrast that corresponds to a given Bouguer anomaly using the algorithm of *Wieczorek and Phillips* (1998). |
| [DownContFilterMA](downcontfilterma.html) | Compute the minimum-amplitude downward continuation filter of *Wieczorek and Phillips* (1998). |
| [DownContFilterMC](downcontfiltermc.html) | Calculate a minimum-curvature downward continuation filter for a given spherical harmonic degree. |
| [NormalGravity](normalgravity.html) | Calculate the normal gravity on a flattened ellipsoid using the formula of Somigliana. |

## Magnetics routines

| Routine name | Description |
| ------------ | ----------- |
| [MakeMagGridDH](makemaggriddh.html) | Create 2D cylindrical maps on a flattened ellipsoid of all three vector components of the magnetic field, the magnitude of the magnetic field, and the magnetic potential. |
| [SHMagPowerSpectrum](shmagpowerspectrum.html) | Compute the power spectrum of the magnetic field given the Schmidt semi-normalized magnetic potential spherical harmonic coefficients. |
| [SHMagPowerL](shmagpowerl.html) | Compute the power of the magnetic field for a single degree L given the Schmidt semi-normalized magnetic potential spherical harmonic coefficients. |

## References

* Wieczorek, M. A. and R. J. Phillips, Potential anomalies on a sphere: applications to the thickness of the lunar crust, J. Geophys. Res., 103, 1715-1724, doi:[10.1029/97JE03136](https://doi.org/10.1029/97JE03136), 1998.
* Wieczorek, M. A. Gravity and topography of the terrestrial planets, Treatise on Geophysics, 10, 165-206, doi:[10.1016/B978-044452748-6/00156-5](https://doi.org/10.1016/B978-044452748-6/00156-5), 2007.
