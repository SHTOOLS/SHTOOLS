---
title: "Complex spherical harmonics"
keywords: spherical harmonics software package, spherical harmonic transform, legendre functions, multitaper spectral analysis, fortran, Python, gravity, magnetic field
sidebar: mydoc_sidebar
permalink: complex-spherical-harmonics.html
summary: SHTOOLS uses by default 4&pi;-normalized spherical harmonic functions that exclude the Condon-Shortley phase factor. Schmidt semi-normalized, orthonormalized, and unnormalized harmonics can be employed in most routines by specifying optional parameters.
toc: true
---

<style>
table:nth-of-type(n) {
    display:table;
    width:100%;
}
table:nth-of-type(n) th:nth-of-type(2) {
    width:75%;
}
</style>

## Definitions: Complex $$4\pi$$-normalized harmonics

Any real square-integrable function can be expressed as a series of spherical harmonic functions

\begin{equation}
f\left(\theta,\phi\right) = \sum_{l=0}^{\infty} \sum_{m=-l}^l f_l^m \, Y_l^m\left(\theta,\phi \right),
\label{eq:f-complex}
\end{equation}

where $$f_l^m$$ is the complex spherical harmonic coefficient, $$Y_l^m$$ is the
corresponding complex spherical harmonic function, $$\theta$$ is co-latitude, $$\phi$$ is longitude, and $$l$$ and $$m$$ are the spherical harmonic degree and order, respectively. The complex spherical harmonics are defined as
\begin{equation}
Y_l^m(\theta,\phi) = \bar{P}_l^m(\cos \theta) \, e^{im\phi},
\end{equation}

where the normalized associated Legendre functions for use with the complex $$4\pi$$-normalized spherical harmonic functions are given by

$$ \begin{eqnarray}
\bar{P}_l^m(\mu) = \sqrt{\left(2l+1\right)\frac{(l-m)!}{(l+m)!}}\, P_{lm}(\mu),
\end{eqnarray} $$

and where $$\delta_{ij}$$ is the Kronecker delta function. The unnormalized
associated Legendre functions are derived from the standard Legendre
polynomials using the relations
\begin{equation}
P_{lm}(\mu) = \left( 1-\mu^2\right)^{m/2} \frac{d^m}{d\mu^m} P_l(\mu)
\end{equation}
and

\begin{equation}
P_l(\mu) = \frac{1}{2^l l!}\frac{d^l}{d\mu^l}\left(\mu^2-1\right)^l.
\end{equation}

The normalized associated Legendre functions are orthogonal for a given value of $$m$$,

$$ \begin{equation}
\int_{-1}^{1} \bar{P}_l^m(\mu) \,\bar{P}_{l'}^m(\mu) = 2 \, \delta_{ll'}.
\end{equation} $$

The complex spherical harmonic functions possess the symmetry relationship for positive
and negative angular orders

$$ \begin{equation}
Y_l^{*m}(\theta,\phi) = (-1)^m \, Y_l^{-m}(\theta,\phi),
\end{equation} $$

where the asterisk denotes complex conjugation, and they satisfy the orthogonality relationship
\begin{equation}
\int_\Omega Y_l^m(\theta,\phi) \, Y_{l'}^{m'}(\theta,\phi) \, d\Omega = 4 \pi \, \delta_{ll'} \, \delta_{mm'},
\end{equation}
where $$d\Omega$$ is the differential surface area on the unit sphere, $$\sin
\theta \, d\theta \, d\phi$$. By multiplying the equation \eqref{eq:f-complex} by $$Y_{l'm'}$$ and integrating over all space, it is straightforward to show that the spherical harmonic coefficients of a function can be calculated by the integral
\begin{equation}
f_l^m = \frac{1}{4\pi} \int_\Omega f(\theta,\phi) \, Y_l^{*m}(\theta,\phi) \, d\Omega.
\end{equation}

Finally, it is noted that if a function defined on the sphere is entirely real, then the real and complex spherical harmonic coefficients are related by

$$ \begin{eqnarray}
f_l^m = \left \lbrace \begin{array}{ll} (f_{lm} - if_{l-m}) /
	\sqrt{2} & \mbox{if $m >0$} \\
	f_{l0} & \mbox{if $m = 0$}\\
	 (-1)^m \, f_l^{*-m} & \mbox{if $m <0$}.
	\end{array} \right.
\end{eqnarray} $$

## Power spectrum

Parseval's theorem in Cartesian geometry relates the integral of a function squared to the sum of the squares of the function's Fourier coefficients. This relation is easily extended to spherical geometry using the orthogonality properties of the spherical harmonic functions. Defining *power* to be the integral of the function squared divided by the area it spans, the total power of a function is equal to a sum over its power spectrum
\begin{equation}
\frac{1}{4\pi} \int_\Omega |f|^2(\theta,\phi) \, d\Omega
= \sum_{l=0}^{\infty} S_{ff}(l),
\end{equation}
where the power spectrum $$S$$ is related to the spherical harmonic coefficients by

$$ \begin{equation}
S_{ff}(l) = \sum\limits_{m=-l}^l |f_l^m|^2.
\end{equation} $$

Similarly, the cross power of two functions $$f$$ and $$g$$ is given by
\begin{equation}
\frac{1}{4\pi} \int_\Omega f(\theta,\phi) \, g^*(\theta,\phi)\, d\Omega
= \sum_{l=0}^{\infty} S_{fg}(l),
\end{equation}
with

$$ \begin{equation}
S_{fg}(l) = \sum\limits_{m=-l}^l f_{lm} \, g_{l}^{*m}.
\end{equation} $$

The power spectrum is unmodified by a rotation of the coordinate system. Furthermore, the numerical values of the power spectrum are independent of the normalization convention used for the spherical harmonic functions (though the mathematical formulae will be different, as given [below](#supported-normalizations)). If the functions $$f$$ and $$g$$ have a zero mean, $$S_{ff}$$ and $$S_{fg}$$ represent the contribution to the variance and covariance, respectively, as a function of degree $$l$$. It should be noted that while the power spectrum of a function is inherently real, the cross power of two functions may be a complex quantity.

$$S$$ is the total power of the function at spherical harmonic degree $$l$$, which in SHTOOLS is called the *power per degree $$l$$*. Alternatively, one can calculate the average power per coefficient at spherical harmonic degree $$l$$, which in SHTOOLS is referred to as the *power per $$lm$$*. Since there are $$(2l+1)$$ spherical harmonic coefficients at degree $$l$$, this is
\begin{equation}
\mbox{power per $lm$} = \frac{S(l)}{(2l+1)}.
\end{equation}
One can also calculate the power from all angular orders over an infinitesimal logarithmic spherical harmonic degree band $$d \log_a l$$, where $$a$$ is the logarithmic base. In SHTOOLS, this is referred to as the *power per $$d\log_a l$$*, which is given by
\begin{equation}
\mbox{power per $d\log_a l$} = S(l)\, l \, \ln a.
\end{equation}
Finally, SHTOOLS defines the *energy* of a function as the integral of its square. The energy spectrum is thus equal to the power spectrum multiplied by $$4\pi$$.

## Condon-Shortley phase factor

The above definitions of the Legendre functions and spherical harmonic functions do not include the Condon-Shortley phase factor of $$(-1)^m$$ that is often employed in the physics and seismology communities [Varshalovich et al. 1988, Dahlen and Tromp 1998]. Nevertheless, this phase can be included in most SHTOOLS routines by specifying the optional parameter

* `csphase = 0` : exclude the Condon-Shortley phase factor (default)
* `csphase = 1` : append the Condon-Shortley phase factor to the Legendre functions.

The choice of the Condon-Shortley phase factor does not affect the numerical value of the power spectrum.

## Supported normalizations

SHTOOLS supports the use of $$4\pi$$-normalized, Schmidt semi-normalized, orthonormalized, and unnormalized spherical harmonic functions. To specify which normalization should be used, it is only necessary to specify the optional parameter `norm` in the Fortran 95 routines or `normalization` in the Python routines:

* `norm = 1`, `normalization = '4pi'`: $$4\pi$$ normalized (default, unless stated otherwise)
* `norm = 2`, `normalization = 'schmidt'`: Schmidt semi-normalized
* `norm = 3`, `normalization = 'unnorm'`: Unnormalized
* `norm = 4`, `normalization = 'ortho'`: Orthonormalized.

Each of these normalizations has slightly different definitions for the normalized Legendre functions, the orthogonality conditions of the Legendre functions and spherical harmonic functions, and the power spectrum. These equations are provided below.

### $$4\pi$$ normalized

| $$ \displaystyle \bar{P}_{l}^m(\mu) =  \sqrt{\left(2l+1\right) \frac{(l-m)!}{(l+m)!}}\, P_{lm}\left(\mu\right) $$ |
| $$\displaystyle \int_{-1}^{1} \bar{P}_{l}^m(\mu) \,\bar{P}_{l'}^m(\mu) \,d\mu=  2 \, \delta_{ll'}$$ |
| $$\displaystyle \int_\Omega {Y_{l}^{m}}^*(\theta,\phi) \, Y_{l'}^{m'}(\theta,\phi) \, d\Omega = 4\pi\, \delta_{ll'}\, \delta_{mm'}$$ |
| $$\displaystyle S_{fg}\left(l\right) = \sum_{m=-l}^l f_l^m \, g_l^{*m}$$ |

### Schmidt semi-normalized

| $$\displaystyle \bar{P}_{l}^m(\mu) = \sqrt{ \frac{(l-m)!}{(l+m)!}}\, P_{lm}\left(\mu\right) $$ |
| $$\displaystyle \int_{-1}^{1} \bar{P}_{l}^m(\mu) \,\bar{P}_{l'}^m(\mu) \,d\mu=  \frac{2}{(2l+1)} \, \delta_{ll'}$$ |
| $$\displaystyle \int_\Omega Y_{l}^{*m}(\theta,\phi) \,Y_{l'}^{m'}(\theta,\phi) \, d\Omega = \frac{4\pi}{(2l+1)} \, \delta_{ll'}\, \delta_{mm'}$$ |
| $$\displaystyle S_{fg}\left(l\right) = \frac{1}{(2l+1)}\sum_{m=-l}^l f_l^m \, g_l^{*m}$$ |

### Orthonormalized

| $$\displaystyle \bar{P}_{l}^m(\mu) = \sqrt{\frac{ \left(2l+1\right)}{4 \pi} \frac{(l-m)!}{(l+m)!}}\, P_{lm}\left(\mu\right)$$ |
| $$\displaystyle \int_{-1}^{1} \bar{P}_{l}^m(\mu) \,\bar{P}_{l'}^m(\mu) \,d\mu= \frac{1}{2 \pi} \, \delta_{ll'}$$ |
| $$\displaystyle\int_\Omega Y_{l}^{*m}(\theta,\phi) \,Y_{l'}^{m'}(\theta,\phi)\, d\Omega =  \delta_{ll'}\, \delta_{mm'}$$ |
| $$\displaystyle S_{fg}\left(l\right) = \frac{1}{4\pi}\sum_{m=-l}^l f_l^m \, g_l^{*m}$$ |

### Unnormalized

| $$\displaystyle \bar{P}_{l}^m(\mu) =  P_{lm}\left(\mu\right)$$ |
| $$\displaystyle \int_{-1}^{1} \bar{P}_{l}^m(\mu) \,\bar{P}_{l'}^m(\mu)= \frac{2}{(2l+1)} \frac{(l+m)!}{(l-m)!} \, \delta_{ll'}$$ |
| $$\displaystyle \int_\Omega Y_{l}^{*m}(\theta,\phi) \,Y_{l'}^{m'}(\theta,\phi)\, d\Omega = \frac{4\pi}{(2l+1)} \frac{(l+m)!}{(l-m)!}\,  \delta_{ll'}\, \delta_{mm'}$$ |
| $$\displaystyle S_{fg}\left(l\right) = \frac{(l+m)!}{(2l+1)(l-m)!}\sum_{m=-l}^l f_l^m \, g_l^{*m}$$ |

## References

* Dahlen, F. A. and J. Tromp, "Theoretical Global Seismology," *Princeton University Press*, Princeton, New Jersey, 1025 pp., 1998.

* Varshalovich, D. A., A. N. Moskalev, and V. K. Khersonskii, "Quantum theory of angular momentum," *World Scientific*, Singapore, 1988.
