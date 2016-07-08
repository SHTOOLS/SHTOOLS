# PreGLQ

Calculate the weights and nodes used in integrating a function by Gauss-Legendre quadrature.

# Usage

`zero`, `w` = pyshtools.PreGLQ (`lower`, `upper`, `n`)

# Returns

`zero` : float, dimension (`n`)
:   The zeros used in the Gauss-Legendre quadrature.

`w` : float, dimension (`n`)
:   The weights used in the Gauss-Legendre quadrature.

# Parameters

`lower` : float
:   The lower bound of the integration.

`upper` : float
:   The upper bound of the integration.

`n` : integer
:   The number of integration points to use. This will integrate exactly a polynomial of degree `2n-1`.

# Description

`PreGLQ` will calculate the weights and zeros used to integrate a function using Gauss-Legendre quadrature. For `n` quadrature points, the integration will be exact if the function is a polynomial of degree `2n-1`, or less. The quadrature nodes correspond to the zeros of the Legendre polynomial of degree `n`. The number of quadrature points required to integrate a polynomial of degree `L` is `ceiling((L+1)/2)`.

To integrate a function between the bounds `lower` and `upper` it is only necessary to calculate the sum of the function evaluated at the nodes `zero` multiplied by the weights.

This is a slightly modified version of the algorithm that was published in NUMERICAL RECIPES.

# References

Press, W.H., S.A. Teukolsky, W.T. Vetterling, and B.P. Flannery, Numerical Recipes in FORTRAN: The Art of Scientific Computing, 2nd ed., Cambridge Univ. Press, Cambridge, UK, 1992.

# See also

[shglq](pyshglq.html)
