# EigValSym

Compute the eigenvalues of a real symmetric matrix.

# Usage

call EigValSym (`ain`, `n`, `eval`, `ul`)

# Parameters

`ain` : input, real(dp), dimension (`n`, `n`)
:   The input real symmetric matrix. By default, only the upper portion of the matrix is used.

`n` : input, integer(int32)
:   The rank of the matrix `ain`.

`eval` : output, real(dp), dimension (`n`)
:   The eigenvalues of `ain`, sorted from largest to smallest.

`ul` : optional, input, character, default = `U`
:   If `U` then the upper portion of the matrix `ain` will be used (default). If `L` then the lower portion of the matrix `ain` will be used.

# Description

`EigValSym` will calculate the eigenvalues of a real symmetric matrix. By default, only the upper portion of the matrix is used, but this can be changed by the optional argument `ul`. The eigenvalues are sorted from largest to smallest. The matrix `ain` is first factorized into a tridiagonal matrix using the LAPACK routine `DSYTRD`, and then the eigenvalues are calculated by a call to `DSTEGR`.

# See also

[eigvalvecsym](eigvalvecsym.html), [eigvalvecsymtri](eigvalvecsymtri.html)
