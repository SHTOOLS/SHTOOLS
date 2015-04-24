# RandomGaussian

Return a pseudo-Gaussian deviate of zero mean and unit variance.

# Usage

`rg` = RandomGaussian (`seed`)

# Parameters

`rg` : output, real\*8
:   The radom Gaussian deviate.

`seed` : input/output, integer
:   Input a negative integer to (re-)initialize the random number generator. Afterwards, this argument should not be modified.

# Description

`RandomGaussian` will return a Gaussian random deviate with unit variance and zero mean. The underlying random number generator uses the algorithm of Park and Miller combined with a Marsaglia shift sequence, which is claimed to have a periodicity of about 3.1 10^18. The random number generator is intialized by calling with a negative value of `seed`, and afterwards, this variable should not be modified. To obtain a Gaussian deviate with a standard deviation of `sigma`, it is only necessary to multiply the unit variance deviate by `sigma`.

This is a slightly modified version of the algorithm that was published in NUMERICAL RECIPES as GASDEV.

# References

Press, W.H., S.A. Teukolsky, W.T. Vetterling, and B.P. Flannery, Numerical Recipes in FORTRAN: The Art of Scientific Computing, 2nd ed., Cambridge Univ. Press, Cambridge, UK, 1992. 

# See also

[randomn](randomn.html)
