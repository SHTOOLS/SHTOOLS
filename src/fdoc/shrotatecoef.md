# SHRotateCoef

Determine the spherical harmonic coefficients of a complex function rotated by three Euler angles.

# Usage

call SHRotateCoef (`x`, `coef`, `rcoef`, `dj`, `lmax`, `exitstatus`)

# Parameters

`x` : input, real\*8, dimension(3)
:   The three Euler angles, alpha, beta, and gamma, in radians.

`coef` : input, real\*8, dimension (2, (`lmax`+1)\*(`lmax`+2)/2)
:   The input complex spherical harmonic coefficients. This is an indexed array where the real and complex components are given by `coef(1,:)` and `coef(2,:)`, respectively. The functions `SHCilmToCindex` and `SHCindexToCilm` are used to convert to and from indexed and `cilm(2,:,:)` form. The coefficients must correspond to unit-normalized spherical harmonics that possess the Condon-Shortley phase convention.

`rcoef` : output, real\*8, dimension (2, (`lmax`+1)\*(`lmax`+2)/2)
:   The spherical harmonic coefficients of the rotated function in indexed form.

`dj` : input, real\*8, dimension (`lmax`+1, `lmax`+1, `lmax`+1)
:   The rotation matrix `dj(pi/2)` obtained from a call to `djpi2`.

`lmax` : input, integer
:   The maximum spherical harmonic degree of the input and output coefficients.

`exitstatus` : output, optional, integer
:   If present, instead of executing a STOP when an error is encountered, the variable exitstatus will be returned describing the error. 0 = No errors; 1 = Improper dimensions of input array; 2 = Improper bounds for input variable; 3 = Error allocating memory; 4 = File IO error.

# Description

`SHRotateCoef` will take the complex spherical harmonic coefficients of a function, rotate it according to the three Euler anlges in `x`, and output the spherical harmonic coefficients of the rotated function. The input and output coefficients are in an indexed form that can be converted to and from `cilm(2,:,:)` form by using the functions `SHCilmToCindex` and `SHCindexToCilm`. The coefficients must correspond to unit-normalized spherical harmonics that possess the Condon-Shortley phase convention. Real spherical harmonics can be converted to and from complex form using `SHrtoc` and `SHctor`. The input rotation matrix `dj` is computed by a call to `djpi2`.

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

[djpi2](djpi2.html), [shrotaterealcoef](shrotaterealcoef.html), [shctor](shctor.html), [shrtoc](shrtoc.html), [shcilmtocindex](shcilmtocindex.html), [shcindextocilm](shcindextocilm.html)
