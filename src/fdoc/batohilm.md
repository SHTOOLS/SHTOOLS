# BAtoHilm

Calculate iteratively the relief along an interface of constant density contrast that corresponds to a given Bouguer anomaly using the algorithm of Wieczorek and Phillips (1998).

# Usage

call BAtoHilm (`cilm`, `ba`, `grid`, `lmax`, `nmax`, `mass`, `r0`, `rho`, `gridtype`, `w`, `plx`, `zero`, `filtertype`, `filterdeg`, `lmaxcalc`, `exitstatus`)

# Parameters

`cilm` : output, real\*8, dimension (2, `lmaxcalc`+1, `lmaxcalc`+1)
:   An estimate of the real spherical harmonic coefficients (geodesy normalized) of relief along an interface with density contrast `rho` that satisfies the Bouguer anomaly `ba`. The degree zero term corresponds to the mean radius of the relief.

`ba` : input, real\*8, dimension (2, `lmaxcalc`+1, `lmaxcalc`+1)
:   The real spherical harmonic coefficients of the Bouguer anomaly referenced to a spherical interface `r0`.

`grid` : input, real\*8, dimension (`lmax`+1, 2\*`lmax`+1) for `gridtype` 1, (2\*`lmax`+2, 2\*`lmax`+2) for `gridtype` 2, (2\*`lmax`+2, 4\*`lmax`+4) for `gridtype` 3
:   The initial estimate for the radii of the interface evaluated on a grid corresponding to a function of maximum spherical harmonic degree `lmax`. This is calculated by a call to either `MakeGridGLQ` or `MakeGridDH`. This grid must contain the degree-0 average radius of the interface.

`lmax` : input, integer
:   The spherical harmonic bandwidth of the input relief `grid`, which determines the dimensions of `grid`. If `lmaxcalc` is not set, this determines also the maximum spherical harmonic degree of the output spherical harmonic coefficients of the relief and the input spherical harmonics of the Bouguer anomaly.

`nmax` : input, integer
:   The maximum order used in the Taylor-series expansion used in calculating the potential coefficients.

`mass` : input, real\*8
:   The mass of the planet in kg.

`r0` : input, real\*8
:   The reference radius of the Bouguer anomaly `ba`.

`rho` : input, real\*8
:   The density contrast of the relief in kg/m^3.

`gridtype` : input, integer
:   1 = Gauss-Legendre grids, calculated using `SHGLQ` and `MakeGridGLQ`. 2 = Equally sampled Driscoll-Healy grids, `n` by `n`, calculated using `MakeGridDH`. 3 = Equally spaced Driscoll-Healy grids, `n` by 2`n`, calculated using `MakeGridDH`.

`w` : optional, input, real\*8, dimension (`lmax`+1)
:   The weights used in the Gauss-Legendre quadrature. These are calculated from a call to `SHGLQ`. If present, one of `plx` or `zero` must also be present.

`plx` : optional, input, real\*8, dimension (`lmax`+1, (`lmax`+1)\*(`lmax`+2)/2)
:   An array of the associated Legendre functions calculated at the nodes used in the Gauss-Legendre quadrature. These are determined from a call to `SHGLQ`.

`zero` : optional, input, real\*8, dimension (`lmax`+1)
:   The nodes used in the Gauss-Legendre quadrature over latitude, calculated by a call to `SHGLQ`.

`filtertype` : optional, input, integer, default = 0
:   Apply a filter when calculating the relief in order to minimize the destabilizing effects of downward continuation which amplify uncertainties in the Bouguer anomaly. If 0, no filtering is applied. If 1, use the minimum amplitude filter `DownContFilterMA`. If 2, use the minimum curvature filter `DownContFilterMC`.

`filterdeg` : optional, input, integer
:   The spherical harmonic degree for which the filter is 0.5.

`lmaxcalc` : optional, input, integer, default = `lmax`
:   The maximum degree that will be calculated in the spherical harmonic expansions.

`exitstatus` : output, optional, integer
:   If present, instead of executing a STOP when an error is encountered, the variable exitstatus will be returned describing the error. 0 = No errors; 1 = Improper dimensions of input array; 2 = Improper bounds for input variable; 3 = Error allocating memory; 4 = File IO error.

# Description

I`BAtoHilm` is used to solve iteratively for the relief along an interface that corresponds to a given Bouguer anomaly. This is equation 18 of Wieczorek and Phillips (1998) which implicitly takes into consideration the finite-amplitude correction. Each iteration takes as input a guess for the relief (specified by `grid`) and outputs the iteratively improved spherical harmonic coefficients of this relief. These coefficients can then be re-expanded and re-input into this routine as the next guess. For the initial guess, it is often sufficient to use the relief predicted using the first-order "mass sheet" approximation. The input relief `grid` can be of one of three type specified by `gridtype`: 1 for Gauss-Legendre grids, 2 for equally sampled Driscoll-Healy grids (`n` by `n`), and 3 for equally spaced Driscoll-Healy grid (`n` by 2`n`).

If the algorithm does not converge, one might want to try damping the initial estimate. Alternatively, iterations of the following form have proven successfulin in damping oscilations between successive iterations:

`h3 = (h2+h1)/2`
`h4 = f(h3)`

It is important to understand that as an intermediate step, this routine calculates the spherical harmonic coefficients of the relief raised to the nth power, i.e., grid\*\*n. As such, if the input function is bandlimited to degree `L`, the resulting function will thus be bandlimited to degree `L*nmax`. This subroutine implicitly assumes that `lmax` is greater than or equal to `L*nmax`. If this is not the case, aliasing will occur. In practice, for accurate results, it is found that `lmax` needs only to be about twice the size of `L`, though this should be verified for each application. Thus, if the input function is considered to be bandlimited to degree `L`, the function should be evaluated on a grid corresponding to a maximum degree of about `2L`.

If the input grid is evaluated on the Gauss-Legendre points, it is necessary to specify the optional parameters `w` and` zero`, or `w` and `plx`, which are calculated by a call to `SHGLQ`. If memory is not an issue, the algorithm can be speeded up by inputing the optional array `plx` of precomputed associated Legendre functions on the Gauss-Legendre nodes. If `plx` is not specified, then it is necessary to input the optional array `zero` that contains the latitudinal Gauss-Legendre quadrature nodes.

This routine uses geodesy 4-pi normalized spherical harmonics that exclude the Condon-Shortley phase.

# References

Wieczorek, M. A. and R. J. Phillips, Potential anomalies on a sphere: applications to the thickness of the lunar crust, J. Geophys. Res., 103, 1715-1724, 1998.

# See also

[batohilmrhoh](batohilmrhoh.html), [cilmplus](cilmplus.html), [cilmplusrhoh](cilmplusrhoh.html), [shexpandglq](shexpandglq.html), [makegridglq](makegridglq.html), [shglq](shglq.html), [shexpanddh](shexpanddh.html), [makegriddh](makegriddh.html), [glqgridcoord](glqgridcoord.html), [downcontfilterma](downcontfilterma.html), [downcontfiltermc](downcontfiltermc.html)
