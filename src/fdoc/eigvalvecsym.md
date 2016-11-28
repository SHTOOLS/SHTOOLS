# EigValVecSym

Compute the eigenvalues and eigenvectors of a real symmetric matrix.

# Usage

call EigValVecSym (`ain`, `n`, `eval`, `evec`, `ul`, `k`, `exitstatus`)

# Parameters

`ain` : input, real\*8, dimension (`n`, `n`)
:   The input real symmetric matrix. By default, only the upper portion of the matrix is used.

`n` : input, integer
:   The rank of the matrix `ain`.

`eval` : output, real\*8, dimension (`k`)
:   The eigenvalues of `ain`, sorted from largest to smallest.

`evec` : output, real\*8, dimension (`n`, `k`)
:   The eigenvectors of `ain`, sorted from largest to smallest eigenvalues. The sign of the first element of each eigenvector is chosen to be positive.

`ul` : input, character, optional, default = `U`
:   If `U` then the upper portion of the matrix `ain` will be used (default). If `L` then the lower portion of the matrix `ain` will be used.

`k` : input, optional, integer, default = `n`
:   The number of largest eigenvalues and corresponding eigenvectors to be output.

`exitstatus` : output, optional, integer
:   If present, instead of executing a STOP when an error is encountered, the variable exitstatus will be returned describing the error. 0 = No errors; 1 = Improper dimensions of input array; 2 = Improper bounds for input variable; 3 = Error allocating memory; 4 = File IO error.

# Description

`EigValVecSym` will calculate the eigenvalues and eigenvectors of a real symmetric matrix. By default, only the upper portion of the matrix is used, but this can be changed by the optional argument `ul`. The eigenvalues and eigenvectors are sorted from largest to smallest eigenvalues. If the optional parameter `k` is specified, then only the `k` largest eigenvalues and their corresponding eigenvectors will be output. 

The matrix `ain` is first factorized into a tridiagonal matrix using the LAPACK routine `DSYTRD`, and then the eigenvalues are calculated by calls to `DSTEGR` and `DORMTR`.

# See also

[eigvalsym](eigvalsym.html), [eigvalvecsymtri](eigvalvecsymtri.html)
