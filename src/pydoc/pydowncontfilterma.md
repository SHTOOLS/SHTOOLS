# DownContFilterMA

Compute the minimum-amplitude downward continuation filter of Wieczorek and Phillips (1998).

# Usage

`wl` = DownContFilterMA (`l`, `half`, `r`, `d`)

# Returns

`wl` : float
:   The amplitude of the downward continuation filter.

# Parameters

`l` : integer
:   The spherical harmonic degree.

`half` : integer
:   The spherical harmonic degree where the filter is equal to 0.5.

`r` : float
:   The reference radius of the gravitational field.

`d` : float
:   The radius of the surface to downward continue to.

# Description

`DownContFilterMA` will calculate the downward continuation filter of Wieczorek and Phillips (1998; eq. 19) for a given spherical harmonic degree `l`. The input parameters include `half`, which is the degree where the filter is equal to 0.5, and `r` and `d`, which are the reference radius of the gravitational field and the radius of the surface to downward continue to, respectively.

# References

Wieczorek, M. A. and R. J. Phillips, Potential anomalies on a sphere: applications to the thickness of the lunar crust, J. Geophys. Res., 103, 1715-1724, 1998.

# See also

[downcontfiltermc](pydowncontfiltermc.html), [batohilmdh](pybatohilmdh.html) [batohilmrhohdh](pybatohilmrhohdh.html)
