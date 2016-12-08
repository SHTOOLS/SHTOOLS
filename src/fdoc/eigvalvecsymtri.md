# EigValVecSymTri

Compute the eigenvalues and eigenvectors of a real symmetric tridiagonal matrix.

# Usage

call EigValVecSymTri (`ain`, `n`, `eval`, `evec`, `ul`, `exitstatus`)

# Parameters

`ain` : input, real\*8, dimension (`n`, `n`)
:   The input real symmetric tridiagonal matrix.

`n` : input, integer
:   The rank of the matrix `ain`.

`eval` : output, real\*8, dimension (`n`)
:   The eigenvalues of `ain`, sorted from largest to smallest.

`evec` : output, real\*8, dimension (`n`, `n`)
:   The eigenvectors of `ain`, sorted from largest to smallest eigenvalues. The sign of the first element of each eigenvector is chosen to be positive.

`ul` : optional, input, character, default = `L`
:   If `U` then the upper portion of the matrix `ain` will be used. If `L` then the lower portion of the matrix `ain` will be used (default).

`exitstatus` : output, optional, integer
:   If present, instead of executing a STOP when an error is encountered, the variable exitstatus will be returned describing the error. 0 = No errors; 1 = Improper dimensions of input array; 2 = Improper bounds for input variable; 3 = Error allocating memory; 4 = File IO error.

# Description

`EigValVecSymTri` will calculate the eigenvalues and eigenvectors of a real symmetric tridiagonal matrix. By default, only the lower portion of the matrix is used, but this can be changed by the optional argument `ul`. The eigenvalues and eigenvectors are sorted from largest to smallest eigenvalues, and the sign of the first element of each eigenvector is chosen to be positive. This routine factors the matrix `ain` using the LAPACK routine `DSTEGR`.

# See also

[eigvalsym](eigvalsym.html), [eigvalvecsym](eigvalvecsym.html)
