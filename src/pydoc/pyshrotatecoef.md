# SHRotateCoef

Determine the spherical harmonic coefficients of a complex function rotated by three Euler angles.

# Usage

`rcoef` = pyshtools.SHRotateCoef (`x`, `coef`, `dj`, [`lmax`])

# Returns

`rcoef` : flaot, dimension (2, (`lmax`+1)\*(`lmax`+2)/2)
:   The spherical harmonic coefficients of the rotated function in indexed form.

# Parameters

`x` : float, dimension(3)
:   The three Euler angles, alpha, beta, and gamma, in radians.

`coef` : float, dimension (2, (`lmaxin`+1)\*(`lmaxin`+2)/2)
:   The input complex spherical harmonic coefficients. This is an indexed array where the real and complex components are given by `coef[0,:]` and `coef[1,:]`, respectively. The functions `SHCilmToCindex` and `SHCindexToCilm` are used to convert to and from indexed and `cilm[:,:,:]` form. The coefficients must correspond to unit-normalized spherical harmonics that possess the Condon-Shortley phase convention.

`dj` : float, dimension (`lmaxin2`+1, `lmaxin2`+1, `lmaxin2`+1)
:   The rotation matrix `dj(pi/2)` obtained from a call to `djpi2`.

`lmax` : optional, integer, default = `lmaxin`
:   The maximum spherical harmonic degree of the input and output coefficients. 

# Description

`SHRotateCoef` will take the complex spherical harmonic coefficients of a function, rotate it according to the three Euler anlges in `x`, and output the spherical harmonic coefficients of the rotated function. The input and output coefficients are in an indexed form that can be converted to and from `cilm[:,:,:]` form by using the functions `SHCilmToCindex` and `SHCindexToCilm`. The coefficients must correspond to unit-normalized spherical harmonics that possess the Condon-Shortley phase convention. Real spherical harmonics can be converted to and from complex form using `SHrtoc` and `SHctor`. The input rotation matrix `dj` is computed by a call to `djpi2`.

The rotation of a coordinate system or body can be viewed in two complementary ways involving three successive rotations. Both methods have the same initial and final configurations, and the angles listed in both schemes are the same.

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

To perform the inverse transform of `x(alpha, beta, gamma)`, use `x(-gamma, -beta, -alpha)`.

Note that this routine uses the "y convention", where the second rotation is with respect to the new y axis. If alpha, beta, and gamma were orginally defined in terms of the "x convention", where the second rotation was with respect to the newx axis, the Euler angles according to the y convention would be `alpha_y=alpha_x-pi/2`, `beta_x=beta_y`, and `gamma_y=gamma_x+pi/2`.

# See also

[djpi2](pydjpi2.html), [shrotaterealcoef](pyshrotaterealcoef.html), [shctor](pyshctor.html), [shrtoc](pyshrtoc.html), [shcilmtocindex](pyshcilmtocindex.html), [shcindextocilm](pyshcindextocilm.html)
