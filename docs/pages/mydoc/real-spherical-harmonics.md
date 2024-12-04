---
title: "Real spherical harmonics"
keywords: spherical harmonics software package, spherical harmonic transform, legendre functions, multitaper spectral analysis, fortran, Python, gravity, magnetic field
sidebar: mydoc_sidebar
permalink: real-spherical-harmonics.html
summary: pyshtools uses by default 4&pi;-normalized spherical harmonic functions that exclude the Condon-Shortley phase factor. Schmidt semi-normalized, orthonormalized, and unnormalized harmonics can be employed in most routines by specifying optional parameters.
toc: true
folder: mydoc
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

## Definitions: Real $$4\pi$$-normalized harmonics

Any real square-integrable function can be expressed as a series of spherical harmonic functions

\begin{equation}
f\left(\theta,\phi\right) = \sum_{l=0}^{\infty} \sum_{m=-l}^l f_{lm} \, Y_{lm}\left(\theta,\phi \right),
\label{eq:f}
\end{equation}

where $$f_{lm}$$ is the spherical harmonic coefficient, $$Y_{lm}$$ is the
corresponding spherical harmonic function, $$\theta$$ is co-latitude, $$\phi$$ is longitude, and $$l$$ and $$m$$ are the spherical harmonic degree and order, respectively. The real spherical harmonics are defined as

$$ \begin{equation}
Y_{lm}(\theta,\phi) = \left \lbrace \begin{array}{ll} \bar{P}_{lm}(\cos
    \theta) \cos m \phi & \mbox{if $m \ge 0$} \\
    \bar{P}_{l|m|}(\cos \theta) \sin |m| \phi & \mbox{if $m < 0$},
    \end{array} \right.
\end{equation} $$

where the normalized associated Legendre functions for use with the $$4\pi$$-normalized spherical harmonic functions are given by

$$ \begin{eqnarray}
\bar{P}_{lm}(\mu) = \sqrt{\left(2-\delta_{m0}\right) \left(2l+1\right)\frac{(l-m)!}{(l+m)!}}\, P_{lm}(\mu)
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
\int_{-1}^{1} \bar{P}_{lm}(\mu) \,\bar{P}_{l'm}(\mu) = 2\left(2-\delta_{0m}\right) \delta_{ll'},
\end{equation} $$

and the spherical harmonic functions are orthogonal for all degrees $$l$$ and orders $$m$$

\begin{equation}
\int_\Omega Y_{lm}(\theta,\phi) \, Y_{l'm'}(\theta,\phi) \, d\Omega = 4 \pi \, \delta_{ll'}\, \delta_{mm'},
\end{equation}
where $$d\Omega$$ is the differential surface area on the unit sphere, $$\sin
\theta \, d\theta \, d\phi$$. By multiplying equation \eqref{eq:f} by $$Y_{l'm'}$$ and integrating over all space, it is straightforward to show that the spherical harmonic coefficients of a function can be calculated by the integral
\begin{equation}
f_{lm} = \frac{1}{4\pi} \int_\Omega f(\theta,\phi) \, Y_{lm}(\theta,\phi) \, d\Omega.
\end{equation}

## Power spectrum

Parseval's theorem in Cartesian geometry relates the integral of a function squared to the sum of the squares of the function's Fourier coefficients. This relation is easily extended to spherical geometry using the orthogonality properties of the spherical harmonic functions. Defining *power* to be the integral of the function squared divided by the area it spans, the total power of a function is equal to a sum over its power spectrum
\begin{equation}
\frac{1}{4\pi} \int_\Omega f^2(\theta,\phi) \, d\Omega
= \sum_{l=0}^{\infty} S_{ff}(l),
\end{equation}
where the power spectrum $$S$$ is related to the spherical harmonic coefficients by

$$ \begin{equation}
S_{ff}(l) = \sum\limits_{m=-l}^l f^2_{lm}.
\end{equation} $$

Similarly, the cross power of two functions $$f$$ and $$g$$ is given by
\begin{equation}
\frac{1}{4\pi} \int_\Omega f(\theta,\phi) \, g(\theta,\phi)\, d\Omega
= \sum_{l=0}^{\infty} S_{fg}(l),
\end{equation}
with

$$ \begin{equation}
S_{fg}(l) = \sum\limits_{m=-l}^l f_{lm} \, g_{lm}.
\end{equation} $$

The power spectrum is unmodified by a rotation of the coordinate system. Furthermore, the numerical values of the power spectrum are independent of the normalization convention used for the spherical harmonic functions (though the mathematical formulae will be different, as given [below](#supported-normalizations)). If the functions $$f$$ and $$g$$ have a zero mean, $$S_{ff}$$ and $$S_{fg}$$ represent the contribution to the variance and covariance, respectively, as a function of degree $$l$$.

$$S$$ is the total power of the function at spherical harmonic degree $$l$$, which in pyshtools is called the *power per degree $$l$$*. Alternatively, one can calculate the average power per coefficient at spherical harmonic degree $$l$$, which in pyshtools is referred to as the *power per $$lm$$*. Since there are $$(2l+1)$$ spherical harmonic coefficients at degree $$l$$, this is
\begin{equation}
\mbox{power per $lm$} = \frac{S(l)}{(2l+1)}.
\end{equation}
One can also calculate the power from all angular orders over an infinitesimal logarithmic spherical harmonic degree band $$d \log_a l$$, where $$a$$ is the logarithmic base. In pyshtools, this is referred to as the *power per $$d\log_a l$$*, which is given by
\begin{equation}
\mbox{power per $d\log_a l$} = S(l)\, l \, \ln a.
\end{equation}
Finally, pyshtools defines the *energy* of a function as the integral of its square. The energy spectrum is thus equal to the power spectrum multiplied by $$4\pi$$.

## Condon-Shortley phase factor

The above definitions of the Legendre functions and spherical harmonic functions do not include the Condon-Shortley phase factor of $$(-1)^m$$ that is often employed in the physics and seismology communities [Varshalovich et al. 1988, Dahlen and Tromp 1998]. Nevertheless, this phase can be included in most pyshtools routines by specifying the optional parameter

* `csphase = 1` : exclude the Condon-Shortley phase factor (default)
* `csphase = -1` : append the Condon-Shortley phase factor to the Legendre functions.

The choice of the Condon-Shortley phase factor does not affect the numerical value of the power spectrum.

## Supported normalizations

pyshtools supports the use of $$4\pi$$-normalized, Schmidt semi-normalized, orthonormalized, and unnormalized spherical harmonic functions. To specify which normalization should be used, it is only necessary to specify the optional parameter `normalization` in the Python routines:

* `normalization = '4pi'`: $$4\pi$$ normalized (default, unless stated otherwise)
* `normalization = 'schmidt'`: Schmidt semi-normalized
* `normalization = 'unnorm'`: Unnormalized
* `normalization = 'ortho'`: Orthonormalized.

Each of these normalizations has slightly different definitions for the normalized Legendre functions, the orthogonality conditions of the Legendre functions and spherical harmonic functions, and the power spectrum. These equations are provided below.

### $$4\pi$$ normalized

| $$ \displaystyle \bar{P}_{lm}(\mu) = \sqrt{\left(2-\delta_{m0}\right) \left(2l+1\right)\frac{(l-m)!}{(l+m)!}}\, P_{lm}(\mu) $$ |
| $$\displaystyle \int_{-1}^{1} \bar{P}_{lm}(\mu) \,\bar{P}_{l'm}(\mu)=  2\left(2-\delta_{0m}\right) \, \delta_{ll'}$$ |
| $$\displaystyle \int_\Omega Y_{lm}(\theta,\phi)\, Y_{l'm'}(\theta,\phi) \, d\Omega = 4 \pi \, \delta_{ll'} \, \delta_{mm'}$$ |
| $$\displaystyle S_{fg}(l) = \sum\limits_{m=-l}^l f_{lm} \, g_{lm}$$ |

### Schmidt semi-normalized

| $$\displaystyle \bar{P}_{lm}(\mu) = \sqrt{\left(2-\delta_{m0}\right) \frac{(l-m)!}{(l+m)!}}\, P_{lm}(\mu) $$ |
| $$\displaystyle \int_{-1}^{1} \bar{P}_{lm}(\mu) \,\bar{P}_{l'm}(\mu)= \frac{2\left(2-\delta_{0m}\right)}{(2l+1)} \, \delta_{ll'}$$ |
| $$\displaystyle \int_\Omega {Y_{lm}}(\theta,\phi) \,Y_{l'm'}(\theta,\phi)\, d\Omega = \frac{4\pi}{(2l+1)} \, \delta_{ll'}\, \delta_{mm'}$$ |
| $$\displaystyle S_{fg}(l) = \frac{1}{(2l+1)}\sum\limits_{m=-l}^l f_{lm} \, g_{lm}$$ |

### Orthonormalized

| $$\displaystyle \bar{P}_{lm}(\mu) = \sqrt{\frac{\left(2-\delta_{0m}\right) \left(2l+1\right)}{4 \pi} \frac{(l-m)!}{(l+m)!}}\, P_{lm}\left(\mu\right)$$ |
| $$\displaystyle \int_{-1}^{1} \bar{P}_{lm}(\mu) \,\bar{P}_{l'm}(\mu)= \frac{\left(2-\delta_{0m}\right)}{2 \pi} \, \delta_{ll'}$$ |
| $$\displaystyle \int_\Omega {Y_{lm}}(\theta,\phi) \,Y_{l'm'}(\theta,\phi)\, d\Omega = \delta_{ll'}\, \delta_{mm'}$$ |
| $$\displaystyle S_{fg}(l) = \frac{1}{4\pi}\sum\limits_{m=-l}^l f_{lm} \, g_{lm}$$ |

### Unnormalized

| $$\displaystyle \bar{P}_{lm}(\mu) =  P_{lm}(\mu)$$ |
| $$\displaystyle \int_{-1}^{1} \bar{P}_{lm}(\mu) \,\bar{P}_{l'm}(\mu)= \frac{2}{(2l+1)} \frac{(l+m)!}{(l-m)!} \, \delta_{ll'}$$ |
| $$\displaystyle \int_\Omega Y_{lm}(\theta,\phi) \,Y_{l'm'}(\theta,\phi)\, d\Omega = \frac{4\pi\, (l+m)!}{(2-\delta_{0m})(2l+1) (l-m)!} \, \delta_{ll'}\, \delta_{mm'} $$ |
| $$\displaystyle S_{fg}(l) = \sum\limits_{m=-l}^l \frac{(l+m)!}{(2-\delta_{0m})(2l+1)(l-m)!} f_{lm} \, g_{lm}$$ |

## References

* Dahlen, F. A. and J. Tromp, "Theoretical Global Seismology," *Princeton University Press*, Princeton, New Jersey, 1025 pp., 1998.

* Varshalovich, D. A., A. N. Moskalev, and V. K. Khersonskii, "Quantum theory of angular momentum," *World Scientific*, Singapore, 1988.
