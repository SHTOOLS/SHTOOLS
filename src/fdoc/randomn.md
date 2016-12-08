# RandomN

Return a pseudo uniform random deviate between 0 and 1 using the algorithm of Park and Miller with a Marsaglia shift sequence.

# Usage

`rn` = RandomN (`seed`)

# Parameters

`rn` : output, real\*8
:   The uniform random deviate.

`seed` : input/output, integer
:   Input a negative integer to (re-)initialize the random number generator. Afterwards, this argument should not be modified.

# Description

`RandomN` will return a uniform random deviate between 0 and 1 (exclusive of the endpoints) using the algorithm of Park and Miller combined with a Marsaglia shift sequence. The random number generator is intialized by calling with a negative value of `seed`, and afterwards this variable should not be modified. The period of this generator is claimed to be about 3.1 10^18 (see Press et al. 2002).

# References

Press, W.H., S.A. Teukolsky, W.T. Vetterling, and B.P. Flannery, Numerical Recipes in FORTRAN 90: Volume 2 of FORTRAN numerical recipes, 2nd ed., Cambridge Univ. Press, Cambridge, UK, 2002.

# See also

[randomgaussian](randomgaussian.html)
