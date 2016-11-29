# CilmMinus

Calculate the gravitational potential interior to relief referenced to a spherical interface using the finite-amplitude algorithm of Wieczorek and Phillips (1998).

# Usage

call CilmMinus (`cilm`, `gridin`, `lmax`, `nmax`, `mass`, `d`, `rho`, `gridtype`, `w`, `zero`, `plx`, `n`, `dref`, `exitstatus`)

# Parameters

`cilm` : output, real\*8, dimension (2, `lmax`+1, `lmax`+1)
:   The real spherical harmonic coefficients (geodesy normalized) of the gravitational potential corresponding to constant density relief referenced to a spherical interface of radius `d`.

`gridin` : input, real\*8, dimension (`lmax`+1, 2\*`lmax`+1) for gridtype 1, (`n`, `n`) for gridtype 2, (`n`, 2\*`n`) for gridtype 3
:   The radii of the interface evaluated on a grid corresponding to a function of maximum spherical harmonic degree `lmax`. This is calculated by a call to either `MakeGridGLQ` or `MakeGridDH`.

`lmax` : input, integer
:   The maximum spherical harmonic degree of the output spherical harmonic coefficients. This degree also determines the dimension of the input relief `gridin` for `gridtype` 1. For Driscoll-Healy grids, `lmax` must be less than or equal to `n/2-1`.

`nmax` : input, integer
:   The maximum order used in the Taylor-series expansion used in calculating the potential coefficients.

`mass` : input, real\*8
:   The mass of the planet in kg.

`d` : output, real\*8
:   The mean radius of the relief in meters.

`rho` : input, real\*8
:   The density contrast of the relief in kg/m^3.

`gridtype` : input, integer
:   1 = Gauss-Legendre grids, calculated using `SHGLQ` and `MakeGridGLQ>`. 2 = Equally sampled Driscoll-Healy grids, `n` by `n`, calculated using `MakeGridDH`. 3 = Equally spaced Driscoll-Healy grids, `n` by 2`n`, calculated using `MakeGridDH`.

`w` : optional, input, real\*8, dimension (`lmax`+1)
:   The weights used in the Gauss-Legendre quadrature, which are required for `gridtype` = 1. These are calculated from a call to `SHGLQ`.

`zero` : optional, input, real\*8, dimension (`lmax`+1)
:   The nodes used in the Gauss-Legendre quadrature over latitude for `gridtype` 1, calculated by a call to `SHGLQ`. One of `plx` or `zero` must be present when `gridtype=1`, but not both.

`plx` : optional, input, real\*8, dimension (`lmax`+1, (`lmax`+1)\*(`lmax`+2)/2)
:   An array of the associated Legendre functions calculated at the nodes used in the Gauss-Legendre quadrature for `gridtype` 1. These are determined from a call to `SHGLQ`. One of `plx` or `zero` must be present when `gridtype=1`, but not both.

`n` : optional, input, integer
:   The number of samples in latitude when using Driscoll-Healy grids. For a function bandlimited to `lmax`, `n=2(lmax+1)`. This is required for gridtypes 2 and 3.

`dref` : optional, input, real\*8
:   The reference radius to be used when calculating both the relief and spherical harmonic coefficients. If this is not specified, this parameter will be set equal to the mean radius of `gridin`.

`exitstatus` : output, optional, integer
:   If present, instead of executing a STOP when an error is encountered, the variable exitstatus will be returned describing the error. 0 = No errors; 1 = Improper dimensions of input array; 2 = Improper bounds for input variable; 3 = Error allocating memory; 4 = File IO error.

# Description

`CilmMinus` will calculate the spherical harmonic coefficients of the gravitational potential that correspond to constant density relief referenced to a spherical interface. This is equation 12 of Wieczorek and Phillips (1998), where the potential is strictly valid only when the coefficients are evaluated at a radius greater than the maximum radius of the relief. The relief is input as a grid, whose type is specified by `gridtype` (1 for Gauss-Legendre quadrature grids, 2 for `n` by `n` Driscoll and Healy sampled grids, and 3 for `n` by 2`n` Driscoll and Healy sampled grids). The input relief `gridin` must correspond to absolute radii. The parameter `nmax` is the order of the Taylor series used in the algorithm to approximate the potential coefficients. By default, the relief and spherical harmonic coefficients will be referenced to the mean radius of `gridin`. However, if the optional parameter `dref` is specified, this will be used instead as the reference radius.

It is important to understand that as an intermediate step, this routine calculates the spherical harmonic coefficients of the relief (referenced to the mean radius of `gridin` or `dref`) raised to the nth power, i.e., `(gridin-d)\*\*n`. As such, if the input function is bandlimited to degree `L`, the resulting function will be bandlimited to degree `L*nmax`. This subroutine implicitly assumes that the `gridin` has an effective spherical harmonic bandwidth greater or equal to this value. (The effective bandwidth is equal to `lmax` for `gridtype` 1, and is `n/2-1` for `gridtype` 2 or 3.) If this is not the case, aliasing will occur. In practice, for accurate results, it is found that the effective bandwidth needs only to be about three times the size of `L`, though this should be verified for each application. Thus, if the input function is considered to be bandlimited to degree `L`, the function should be evaluated on a grid corresponding to a maximum degree of about `3*`L. Aliasing effects can be partially mitigated by using Driscoll and Healy `n` by 2`n` grids.

If the input grid is evaluated on the Gauss-Legendre points, it is necessary to specify the optional parameters `w` and `zero`, or `w` and `plx`, which are calculated by a call to `SHGLQ`. In contast, if Driscoll-Healy grids are used, it is necessary to specify the optional parameter `n`. If memory is not an issue, the algorithm can be speeded up when using Gauss-Lengendre grids by inputing the optional array `plx` (along with `w`) of precomputed associated Legendre functions on the Gauss-Legendre nodes. Both of these variables are computed by a call to `SHGLQ`.

This routine uses geodesy 4-pi normalized spherical harmonics that exclude the Condon-Shortley phase.

# References

Wieczorek, M. A. and R. J. Phillips, Potential anomalies on a sphere: applications to the thickness of the lunar crust, J. Geophys. Res., 103, 1715-1724, 1998.

# See also

[cilmminusrhoh](cilmminusrhoh.html), [cilmplus](cilmplus.html), [cilmplusrhoh](cilmplusrhoh.html), [shexpandglq](shexpandglq.html), [makegridglq](makegridglq.html), [shglq](shglq.html), [glqgridcoord](glqgridcoord.html), [makegriddh](makegriddh.html)
