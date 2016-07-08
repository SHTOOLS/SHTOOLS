# EigValVecSymTri

Compute the eigenvalues and eigenvectors of a real symmetric tridiagonal matrix.

# Usage

`eval`, `evec` = pyshtools.EigValVecSymTri (`ain`, [`n`, `ul`])

# Returns

`eval` : float, dimension (`n`)
:   The eigenvalues of `ain`, sorted from largest to smallest.

`evec` : float, dimension (`n`, `n`)
:   The eigenvectors of `ain`, sorted from largest to smallest eigenvalues. The sign of the first element of each eigenvector is chosen to be positive.

# Parameters

`ain` : float, dimension (`nin`, `nin`)
:   The input real symmetric tridiagonal matrix. 

`n` : optional, integer, default = `nin`
:   The rank of the matrix `ain`.

`ul` : optional, character, default = `L`
:   If `U` then the upper portion of the matrix `ain` will be used. If `L` then the lower portion of the matrix `ain` will be used (default).

# Description

`EigValVecSymTri` will calculate the eigenvalues and eigenvectors of a real symmetric tridiagonal matrix. By default, only the lower portion of the matrix is used, but this can be changed by the optional argument `ul`. The eigenvalues and eigenvectors are sorted from largest to smallest eigenvalues, and the sign of the first element of each eigenvector is chosen to be positive. This routine factors the matrix `ain` using the LAPACK routine `DSTEGR`.

# See also

[eigvalsym](pyeigvalsym.html), [eigvalvecsym](pyeigvalvecsym.html)
