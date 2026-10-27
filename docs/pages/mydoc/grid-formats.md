---
title: "Grid formats"
keywords: spherical harmonics software package, spherical harmonic transform, legendre functions, multitaper spectral analysis, fortran, Python, gravity, magnetic field
sidebar: mydoc_sidebar
permalink: grid-formats.html
summary: pyshtools supports equally sampled, equally spaced, and Gauss-Legendre quadrature grids.
toc: true
folder: mydoc
---

<style>
table:nth-of-type(n) {
    display:table;
    width:100%;
}
</style>

## Supported grid formats

pyshtools uses grid formats that accommodate exact quadrature. These include regularly spaced grids that satisfy the *Driscoll and Healy* (1994) sampling theorem, endpoint-including Clenshaw-Curtis grids, and Gauss-Legendre quadrature grids [e.g., Press et al. 1992]. The rows and columns correspond to latitude and longitude nodes, respectively. The first row is closest to the north pole and the first column is at 0 degrees E.

### Gauss-Legendre Quadrature

For the case of Gauss-Legendre quadrature (`GLQ`), the quadrature is exact when the function $$f$$ is sampled in latitude at the $$(L+1)$$ zeros of the Legendre Polynomial of degree $$(L+1)$$. Since the function also needs to be sampled on $$(2L+1)$$ equally space grid nodes for the Fourier transforms in longitude, the function $$f$$ is sampled on a grid of size $$(L+1)\times(2L+1)$$. The redundant data points at 360$$^{\circ}$$ E longitude are not required by the spherical harmonic transformation routines, but can be computed by specifying the optional argument `extend`.

### Driscoll and Healy [1994]

The second type of grid is for data that are sampled on regular grids. As shown by *Driscoll and Healy* [1994], an exact quadrature exists when the function $$f$$ is sampled at $$N$$ equally spaced nodes in latitude and $$N$$ equally spaced nodes in longitude. For this sampling (`DH`), the grids make use of the longitude band at 90$$^{\circ}$$ N, but not 90$$^{\circ}$$ S, and the number of samples is $$2(L+1)$$, which is always even. Given that the sampling in latitude was imposed a priori, these grids contain almost twice as many samples in latitude as the grids used with Gauss-Legendre quadrature. It should be noted that for this quadrature, the longitude band at 90$$^{\circ}$$ N is ultimately downweighted to zero, and hence has no influence on the returned spherical harmonic coefficients.

For geographic data, it is common to work with grids that are equally spaced in degrees latitude and longitude. pyshtools provides the option of using grids of size $$N\times2N$$, and when performing the Fourier transforms for this case (`DH2`), the coefficients $$c_{lm}$$ and $$s_{lm}$$ with $$m>L$$ are discarded. The redundant data points at 360$$^{\circ}$$ E longitude and the latitudinal band at 90$$^{\circ}$$ S are not required by the spherical harmonic transformation routines, but can be computed by specifying the optional argument `extend`.

### Clenshaw-Curtis

For maximum spherical harmonic degree $$L$$, a non-extended Clenshaw-Curtis grid (`CC`) contains $$L+2$$ equally spaced latitude nodes from 90$$^{\circ}$$ N through 90$$^{\circ}$$ S and $$2L+2$$ equally spaced longitude nodes from 0 through less than 360$$^{\circ}$$ E. Both poles are included. Thus $$L=N_{lat}-2$$ is the maximum recoverable spherical harmonic degree. Setting `extend=True` appends only the redundant longitude column at 360$$^{\circ}$$ E. CC does not use the Driscoll-Healy sampling parameter.

For example, a band-limited field can be synthesized and analyzed without changing its grid identity:

```python
grid = coeffs.expand(grid='CC', backend='ducc', extend=False)
recovered = grid.expand(backend='ducc')
```

## Comparison of DH, CC, and GLQ grids

The properties of Driscoll and Healy [1994], Clenshaw-Curtis, and Gauss-Legendre quadrature grids are summarized below. Here, $$N$$ is the DH latitude count and $$N_{lat}$$ is the CC latitude count; for GLQ, $$N=L+1$$.

|                                    | DH1                | DH2                | CC                   | GLQ                        |
| ---------------------------------- | ------------------ | ------------------ | -------------------- | -------------------------- |
| Name                               | Driscoll and Healy | Driscoll and Healy | Clenshaw-Curtis      | Gauss-Legendre Quadrature  |
| Shape ($$N_{lat} \times N_{lon}$$) | $$N \times N$$     | $$N \times 2N$$    | $$(L+2)\times(2L+2)$$ | $$N \times 2N$$            |
| $$L$$                              | $$N/2-1$$          | $$N/2-1$$          | $$N_{lat}-2$$        | $$N-1$$                    |
| Latitude poles                     | North only         | North only         | North and South      | Neither                    |
| $$\Delta \theta$$                  | $$180^{\circ}/N$$  | $$180^{\circ}/N$$  | $$180^{\circ}/(L+1)$$ | Variable                   |
| $$\Delta \phi$$                    | $$360^{\circ}/N$$  | $$180^{\circ}/N$$  | $$180^{\circ}/(L+1)$$ | $$360^{\circ}/(2N-1)$$     |

The existing figure below illustrates the `DH1`, `DH2`, and `GLQ` grids only. The `DH1` and `DH2` grids have the same latitude sampling, but `DH2` has twice as many longitude nodes as `DH1`. `GLQ` is regularly sampled in longitude and irregularly sampled in latitude. The red points at the south pole and 360$$^{\circ}$$ E are optional for the illustrated DH and GLQ transforms; CC pole and extension semantics are described above.

{% include image.html file="grids.png" alt="Spherical harmonic grid formats" caption="Schematic diagram illustrating the properties of the grids used with the Gauss-Legendre quadrature and Driscoll and Healy routines. The red points are not required by the spherical harmonic transform routines, but can be computed by specifying the optional argument `extend`." %}

## References

* Driscoll, J. R. and D. M. Healy, Computing Fourier transforms and convolutions on the 2-sphere, Adv. Appl. Math., 15, 202-250, doi:[10.1006/aama.1994.1008](https://doi.org/10.1006/aama.1994.1008), 1994.

* Press, W. H., S. A. Teukolsky, W. T. Vetterling, and B. P. Flannery, "Numerical Recipes in FORTRAN: The Art of Scientific Computing," 2nd ed., Cambridge Univ. Press, Cambridge, UK, 1992.
