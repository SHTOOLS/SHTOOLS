# SHRotateRealCoef

Determine the spherical harmonic coefficients of a real function rotated by three Euler angles.

# Usage

`cilmrot` = SHRotateRealCoef (`cilm`, `x`, `dj`, [`lmax`])

# Returns

`cilmrot` : float, dimension (2, `lmax`+1, `lmax`+1)
:   The spherical harmonic coefficients of the rotated function, normalized for use with the geodesy 4-pi spherical harmonics.

# Parameters

`cilm` : float, dimension (2, `lmaxin`+1, `lmaxin`+1)
:   The input real spherical harmonic coefficients. The coefficients must correspond to geodesy 4-pi normalized spherical harmonics that do not possess the Condon-Shortley phase convention.

`x` : float, dimension(3)
:   The three Euler angles, alpha, beta, and gamma, in radians.

`dj` : float, dimension (`lmaxin2`+1, `lmaxin2`+1, `lmaxin2`+1)
:   The rotation matrix `dj(pi/2)`, obtained from a call to `djpi2`.

`lmax` : optional, integer, default = `lmaxin`
:   The maximum spherical harmonic degree of the input and output coefficients.

# Description

`SHRotateRealCoef` will take the real spherical harmonic coefficients of a function, rotate it according to the three Euler anlges in `x`, and output the spherical harmonic coefficients of the rotated function. The input and output coefficients must correspond to geodesy 4-pi normalized spherical harmonics that do not possess the Condon-Shortley phase convention. The input rotation matrix `dj` is computed by a call to `djpi2`.

The rotation of a coordinate system or body can be viewed in two complementary ways involving three successive rotations. Both methods have the same initial and final configurations, and the angles listed in both schemes are the same. This routine uses the 'y convention', where the second rotation axis corresponds to the y axis.

`Scheme A:`

(I) Rotation about the z axis by alpha.
(II) Rotation about the new y axis by beta.
(III) Rotation about the new z axis by gamma.

`Scheme B:`

(I) Rotation about the z axis by gamma.
(II) Rotation about the initial y axis by beta.
(III) Rotation about the initial z axis by alpha.

The rotations can further be viewed either as a rotation of the coordinate system or the physical body. For a rotation of the coordinate system without rotation of the physical body, use

`x(alpha, beta, gamma)`.

For a rotation of the physical body without rotation of the coordinate system, use 

`x(-gamma, -beta, -alpha)`.

The inverse transform of `x(alpha, beta, gamma)` is `x(-gamma, -beta, -alpha)`.

Note that this routine uses the "y convention", where the second rotation is with respect to the new y axis. If alpha, beta, and gamma were originally defined in terms of the "x convention", where the second rotation was with respect to the new x axis, the Euler angles according to the y convention would be `alpha_y=alpha_x-pi/2`, `beta_x=beta_y`, and `gamma_y=gamma_x+pi/2`.

This routine first converts the real coefficients to complex form using `SHrtoc`. Then the coefficients are converted to indexed form using `SHCilmToCindex`, these are sent to `SHRotateCoef`, the result if converted back to `cilm` complex form using `SHCindexToCilm`, and these are finally converted back to real form using `SHctor`.

# See also

[djpi2](pydjpi2.html), [shrotatecoef](pyshrotatecoef.html), [shctor](pyshctor.html), [shrtoc](pyshrtoc.html), [shcilmtocindex](pyshcilmtocindex.html), [shcindextocilm](pyshcindextocilm.html)
