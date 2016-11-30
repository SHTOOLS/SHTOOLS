# PreGLQ

Calculate the weights and nodes used in integrating a function by Gauss-Legendre quadrature.

# Usage

call PreGLQ (`lower`, `upper`, `n`, `zero`, `w`, `exitstatus`)

# Parameters

`lower` : input, real\*8
:   The lower bound of the integration.

`upper` : input, real\*8
:   The upper bound of the integration.

`n` : input, integer
:   The number of integration points to use. This will integrate exactly a polynomial of degree `2n-1`.

`zero` : output, real\*8, dimension (`n`)
:   The zeros used in the Gauss-Legendre quadrature.

`w` : output, real\*8, dimension (`n`)
:   The weights used in the Gauss-Legendre quadrature.

`exitstatus` : output, optional, integer
:   If present, instead of executing a STOP when an error is encountered, the variable exitstatus will be returned describing the error. 0 = No errors; 1 = Improper dimensions of input array; 2 = Improper bounds for input variable; 3 = Error allocating memory; 4 = File IO error.

# Description

`PreGLQ` will calculate the weights and zeros used to integrate a function using Gauss-Legendre quadrature. For `n` quadrature points, the integration will be exact if the function is a polynomial of degree `2n-1`, or less. The quadrature nodes correspond to the zeros of the Legendre polynomial of degree `n`. The number of quadrature points required to integrate a polynomial of degree `L` is `ceiling((L+1)/2)`.

To integrate a function between the bounds `lower` and `upper` it is only necessary to calculate the sum of the function evaluated at the nodes `zero` multiplied by the weights.

This is a slightly modified version of the algorithm that was published in NUMERICAL RECIPES.

# References

Press, W.H., S.A. Teukolsky, W.T. Vetterling, and B.P. Flannery, Numerical Recipes in FORTRAN: The Art of Scientific Computing, 2nd ed., Cambridge Univ. Press, Cambridge, UK, 1992.

# See also

[shglq](shglq.html)
