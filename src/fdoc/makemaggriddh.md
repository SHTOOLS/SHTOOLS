# MakeMagGridDH

Create 2D cylindrical maps on a flattened ellipsoid of all three vector components of the magnetic field, the magnitude of the magnetic field, and the magnetic potential.

# Usage

call MakeMagGridDH (`cilm`, `lmax`, `r0`, `a`, `f`, `rad`, `theta`, `phi`, `total`, `n`, `sampling`, `lmaxcalc`, `potgrid`, `exitstatus`)

# Parameters

`cilm` : input, real\*8, dimension (2, `lmax`+1, `lmax`+1)
:   The real Schmidt semi-normalized spherical harmonic coefficients to be expanded in the space domain. The coefficients `C1lm` and `C2lm` refer to the cosine (`Clm`) and sine (`Slm`) coefficients, respectively, with `Clm=cilm(1,l+1,m+1)` and `Slm=cilm(2,l+1,m+1)`. Alternatively, `C1lm` and `C2lm` correspond to the positive and negative order coefficients, respectively.

`lmax` : input, integer
:   The maximum spherical harmonic degree of the coefficients `cilm`. This determines the number of samples of the output grids, `n=2*lmax+2`, and the latitudinal sampling interval, `90/(lmax+1)`.

`r0` : input, real\*8
:   The reference radius of the spherical harmonic coefficients.

`a` : input, real\*8 
:   The semi-major axis of the flattened ellipsoid on which the field is computed.

`f` : input, real\*8
:   The flattening of the reference ellipsoid: i.e., `F=(R_equator-R_pole)/R_equator`.

`rad` : output, real\*8, dimension(2\*`lmax`+2, `sampling`\*(2\*`lmax`+2))
:   A 2D equally sampled (`n` by `n`) or equally spaced (`n` by 2`n`) grid of the radial component of the magnetic field corresponding to the input spherical harmonic coefficients `cilm`. The first latitudinal band corresponds to 90 N, the latitudinal band for 90 S is not included, and the latitudinal sampling interval is 180/`n` degrees. The first longitudinal band is 0 E, the longitudinal band for 360 E is not included, and the longitudinal sampling interval is 360/`n` for an equally sampled and 180/`n` for an equally spaced grid, respectively.

`theta` : output, real\*8, dimension(2\*`lmax`+2, `sampling`\*(2\*`lmax`+2))
:   A 2D equally sampled or equally spaced grid of the theta component of the magnetic field.

`phi` : output, real\*8, dimension(2\*`lmax`+2, `sampling`\*(2\*`lmax`+2))
:   A 2D equally sampled or equally spaced grid of the phi component of the magnetic field. 

`total` : output, real\*8, dimension(2\*`lmax`+2, `sampling`\*(2\*`lmax`+2))
:   A 2D equally sampled or equally spaced grid of the total magnetic field strength. 

`n` : output, integer
:   The number of samples in latitude of the output grids. This is equal to `2lmax+2`.

`sampling` : optional, input, integer, default = 1
:   If 1 (default) the output grids are equally sampled (`n` by `n`). If 2, the grids are equally spaced (`n` by 2`n`).

`lmaxcalc` : optional, input, integer, default = `lmax`
:   The maximum spherical harmonic degree used in evaluating the functions. This must be less than or equal to `lmax`.

`potgrid` : output, real\*8, dimension(2\*`lmax`+2, `sampling`\*(2\*`lmax`+2))
:   A 2D equally sampled or equaly spaced grid of the magnetic potential.

`exitstatus` : output, optional, integer
:   If present, instead of executing a STOP when an error is encountered, the variable exitstatus will be returned describing the error. 0 = No errors; 1 = Improper dimensions of input array; 2 = Improper bounds for input variable; 3 = Error allocating memory; 4 = File IO error.

# Description

`MakeMagGridDH` will create 2-dimensional cylindrical maps from the spherical harmonic coefficients `cilm` of all three components of the magnetic field, the total field strength, and the magnetic potential. The magnetic potential is given by

`V = R0 Sum_{l=1}^LMAX (R0/r)^{l+1} Sum_{m=-l}^l C_{lm} Y_{lm}`

and the magnetic field is

`B = - Grad V`.

The coefficients are referenced to a radius `r0`, and the function is computed on a flattened ellipsoid with semi-major axis `a` (i.e., the mean equatorial radius) and flattening `f`.

The default is to calculate grids for use in the Driscoll and Healy routines that are equally sampled (`n` by `n`), but this can be changed to calculate equally spaced grids (`n` by 2`n`) by setting the optional argument `sampling` to 2. The input value of `lmax` determines the number of samples, `n=2lmax+2`, and the latitudinal sampling interval, `90/(lmax+1)`. The first latitudinal band of the grid corresponds to 90 N, the latitudinal band for 90 S is not calculated, and the latitudinal sampling interval is 180/`n` degrees. The first longitudinal band is 0 E, the longitudinal band for 360 E is not calculated, and the longitudinal sampling interval is 360/`n` for equally sampled and 180/`n` for equally spaced grids, respectively.
