# CilmPlusDH

Calculate the gravitational potential exterior to relief referenced to a spherical interface using the finite-amplitude algorithm of Wieczorek and Phillips (1998).

# Usage

`cilm`, `d` = CilmPlusDH (`gridin`, `nmax`, `mass`, `rho`, [`lmax`,  `n`, `dref`, `sampling`])

# Returns

`cilm` : float, dimension (2, `lmax`+1, `lmax`+1)
:   The real spherical harmonic coefficients (geodesy normalized) of the gravitational potential corresponding to constant density relief referenced to a spherical interface of radius `d`.

`d` : float
:   The mean radius of the relief in meters.

# Parameters

`gridin` : float, dimension (`nin`, `sampling`\*`nin`)
:   The radii of the interface evaluated on a grid, determined by a call to `MakeGridDH`.

`nmax` : integer
:   The maximum order used in the Taylor-series expansion used in calculating the potential coefficients.

`mass` : float
:   The mass of the planet in kg.

`rho` : float
:   The density contrast of the relief in kg/m^3.

`lmax` : optional, integer, default = `n/2-1`
:   The maximum spherical harmonic degree of the output spherical harmonic coefficients. `lmax` must be less than or equal to `n/2-1`. 

`n` : optional, integer, default = `nin`
:   The number of samples in latitude when using Driscoll-Healy grids. For a function bandlimited to `lmax`, `n=2(lmax+1)`.

`dref` : optional, float
:   The reference radius to be used when calculating both the relief and spherical harmonic coefficients. If this is not specified, this parameter will be set equal to the mean radius of `gridin`.

`sampling` : optional, integer, default determined by dimensions of `gridin`
:   If 1 the output grids are equally sampled (`n` by `n`). If 2, the grids are equally spaced (`n` by 2`n`).

# Description

`CilmPlus` will calculate the spherical harmonic coefficients of the gravitational potential exterior to constant density relief referenced to a spherical interface. This is equation 10 of Wieczorek and Phillips (1998), where the potential is strictly valid only when the coefficients are evaluated at a radius greater than the maximum radius of the relief. The input relief `gridin` must correspond to absolute radii. The parameter `nmax` is the order of the Taylor series used in the algorithm to approximate the potential coefficients. By default, the relief and spherical harmonic coefficients will be referenced to the mean radius of `gridin`. However, if the optional parameter `dref` is specified, this will be used instead as the reference radius. 

It is important to understand that as an intermediate step, this routine calculates the spherical harmonic coefficients of the relief (referenced to the mean radius of `gridin` or `dref`) raised to the nth power, i.e., `(gridin-d)\*\*n`. As such, if the input function is bandlimited to degree `L`, the resulting function will be bandlimited to degree `L*nmax`. This subroutine implicitly assumes that the `gridin` has an effective spherical harmonic bandwidth greater or equal to this value. (The effective bandwidth is equal to `n/2-1`) If this is not the case, aliasing will occur. In practice, for accurate results, it is found that the effective bandwidth needs only to be about three times the size of `L`, though this should be verified for each application. Thus, if the input function is considered to be bandlimited to degree `L`, the function should be evaluated on a grid corresponding to a maximum degree of about `3*`L.

This routine uses geodesy 4-pi normalized spherical harmonics that exclude the Condon-Shortley phase.

# References

Wieczorek, M. A. and R. J. Phillips, Potential anomalies on a sphere: applications to the thickness of the lunar crust, J. Geophys. Res., 103, 1715-1724, 1998.

# See also

[cilmplusrhohdh](pycilmplusrhohdh.html), [cilmminusdh](pycilmminusdh.html), [cilmminusrhohdh](pycilmminusrhohdh.html), [makegriddh](pymakegriddh.html)
