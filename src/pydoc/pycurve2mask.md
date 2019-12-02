# Curve2Mask

Given a set of latitude and longitude coordinates representing a closed curve, output a gridded binary mask.

# Usage

`mask_dh` = Curve2Mask (`n`, `profile`, `np`, [`nprofile`, `sampling`])

# Returns

`dh_mask` : integer, dimension (`n`, `n`\*`sampling`)
:   A Driscoll and Healy (1994) sampled grid describing the concentration region R. All elements on output will either be 1 (for inside the concentration region) or 0 (for outside R).

# Parameters

`n` : integer
:   The number of latitudinal samples in `dh_mask`. The effective spherical harmonic bandwidth of this grid is `L=n/2-1`.

`profile` : float, dimension (`nprofilein`, 2)
:   List of latitude [:,0] and longitude [:,1] coordinates in degrees specifying a single closed curve.

`np` : integer
:   The value of the returned mask at the North pole (90N, 0E). If the North pole is outside of the concentration region, set this to 0; if it is inside the concentration region, set this to 1.

`nprofile` : optional, integer, default = `nprofilein`
:   The number of coordinates in the curve `profile`.

`sampling` : optional, integer, default = 1
:   For 1, `dh_mask` is dimensioned as (`n`, `n`), whereass for 2, `dh_mask` is dimensssioned as (`n`, `2n`).

# Description

`Curve2Mask` will take a list of latitude and longitude coordinates that represents a single closed curve, and output a mask `mask_dh` that contains 1s and 0s where the grid nodes are inside and outside of the curve, respectively. `mask_dh` will be sampled according to the Driscoll and Healy (1994) sampling theorem with `n` samples in latitude, and either `n` (`sampling=1`) or `2n` (`sampling=2`) samples in longitude. It is necessary to specify a single point as being inside or outside of the curve, and for this the value at the North pole (90N, 0E) must be specified as either 0 or 1.

Longitudes of the curve can span the range from -360 to 720 degrees. If the longitudes of two adjacent points differ by more than 180 degrees, it will be assumed that the curve passes from 360 to 0 degrees, or from -180 to 180 degrees.

# References

Driscoll, J.R. and D.M. Healy, Computing Fourier transforms and convolutions on the 2-sphere, Adv. Appl. Math., 15, 202-250, 1994.

# See also

[shreturntapersmap](pyshreturntapersmap.html), [computedmap](pycomputedmap.html)
