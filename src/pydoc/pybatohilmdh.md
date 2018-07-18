# BAtoHilmDH

Calculate iteratively the relief along an interface of constant density contrast that corresponds to a given Bouguer anomaly using the algorithm of Wieczorek and Phillips (1998).

# Usage

`cilm` = BAtoHilmDH (`ba`, `grid`, `nmax`, `mass`, `r0`, `rho`, [`filter_type`, `filter_deg`, `lmax`, `lmax_calc`, `smapling`])

# Returns

`cilm` : float, dimension (2, `lmax_calc`+1, `lmax_calc`+1)
:   An estimate of the real spherical harmonic coefficients (geodesy normalized) of relief along an interface with density contrast `rho` that satisfies the Bouguer anomaly `ba`. The degree zero term corresponds to the mean radius of the relief.

# Parameters

`ba` : float, dimension (2, `lmax_calc`+1, `lmax_calc`+1)
:   The real spherical harmonic coefficients of the Bouguer anomaly referenced to a spherical interface `r0`.

`grid` : float, dimension (2\*`lmaxin`+2, `sampling`\*(2\*`lmaxin`+2)) 
:   The initial estimate for the radii of the interface evaluated on a grid corresponding to a function of maximum spherical harmonic degree `lmaxin`. This is calculated by a call to `MakeGridDH` and must contain the degree-0 average radius of the interface.

`nmax` : integer
:   The maximum order used in the Taylor-series expansion used in calculating the potential coefficients.

`mass` : float
:   The mass of the planet in kg.

`r0` : float
:   The reference radius of the Bouguer anomaly `ba`.

`rho` : float
:   The density contrast of the relief in kg/m^3.

`filter_type` : optional, integer, default = 0
:   Apply a filter when calculating the relief in order to minimize the destabilizing effects of downward continuation which amplify uncertainties in the Bouguer anomaly. If 0, no filtering is applied. If 1, use the minimum amplitude filter `DownContFilterMA`. If 2, use the minimum curvature filter `DownContFilterMC`. 

`filter_deg` : optional, integer, default = 0
:   The spherical harmonic degree for which the filter is 0.5.

`lmax` : optional, integer, default = `lmaxin`
:   The spherical harmonic bandwidth of the input relief `grid`, which determines the dimensions of `grid`. If `lmax_calc` is not set, this determines also the maximum spherical harmonic degree of the output spherical harmonic coefficients of the relief and the input spherical harmonics of the Bouguer anomaly.

`lmax_calc` : optional, integer, default = `lmax`
:   The maximum degree that will be calculated in the spherical harmonic expansions.

`sampling` : optional, integer, default set by dimensions of `grid`
:   If 1 the output grids are equally sampled (`n` by `n`). If 2, the grids are equally spaced (`n` by 2`n`).

# Description

`BAtoHilm` is used to solve iteratively for the relief along an interface that corresponds to a given Bouguer anomaly. This is equation 18 of Wieczorek and Phillips (1998) which implicitly takes into consideration the finite-amplitude correction. Each iteration takes as input a guess for the relief (specified by `grid`) and outputs the iteratively improved spherical harmonic coefficients of this relief. These coefficients can then be re-expanded and re-input into this routine as the next guess. For the initial guess, it is often sufficient to use the relief predicted using the first-order "mass sheet" approximation.

If the algorithm does not converge, one might want to try damping the initial estimate. Alternatively, iterations of the following form have proven successfulin in damping oscilations between successive iterations:

`h3 = (h2+h1)/2`  
`h4 = f(h3)`  

It is important to understand that as an intermediate step, this routine calculates the spherical harmonic coefficients of the relief raised to the nth power, i.e., grid\*\*n. As such, if the input function is bandlimited to degree `L`, the resulting function will thus be bandlimited to degree `L*nmax`. This subroutine implicitly assumes that `lmax` is greater than or equal to `L*nmax`. If this is not the case, aliasing will occur. In practice, for accurate results, it is found that `lmax` needs only to be about twice the size of `L`, though this should be verified for each application. Thus, if the input function is considered to be bandlimited to degree `L`, the function should be evaluated on a grid corresponding to a maximum degree of about `2L`.

This routine uses geodesy 4-pi normalized spherical harmonics that exclude the Condon-Shortley phase.

# References

Wieczorek, M. A. and R. J. Phillips, Potential anomalies on a sphere: applications to the thickness of the lunar crust, J. Geophys. Res., 103, 1715-1724, 1998.

# See also

[batohilmrhohdh](pybatohilmrhohdh.html), [cilmplusdh](pycilmplusdh.html), [cilmplusrhohdh](pycilmplusrhohdh.html), [shexpanddh](pyshexpanddh.html), [makegriddh](pymakegriddh.html), [downcontfilterma](pydowncontfilterma.html), [downcontfiltermc](pydowncontfiltermc.html)
