# EigValVecSym

Compute the eigenvalues and eigenvectors of a real symmetric matrix.

# Usage

`eval`, `evec` = pyshtools.EigValVecSym (`ain`, [`n`, `ul`, `k`])

# Returns

`eval` : float, dimension (`k`)
:   The eigenvalues of `ain`, sorted from largest to smallest.

`evec` : float, dimension (`n`, `k`)
:   The eigenvectors of `ain`, sorted from largest to smallest eigenvalues. The sign of the first element of each eigenvector is chosen to be positive.

# Parameters

`ain` : float, dimension (`nin`, `nin`)
:   The input real symmetric matrix. By default, only the upper portion of the matrix is used.
	
`n` : optional, integer, default = `nin`
:   The rank of the matrix `ain`.
	
`ul` : optional, character, default = `U`
:   If `U` then the upper portion of the matrix `ain` will be used (default). If `L` then the lower portion of the matrix `ain` will be used.

`k` : optional, optional, integer, default = `nin`
:   The number of largest eigenvalues and corresponding eigenvectors to be output.

# Description

`EigValVecSym` will calculate the eigenvalues and eigenvectors of a real symmetric matrix. By default, only the upper portion of the matrix is used, but this can be changed by the optional argument `ul`. The eigenvalues and eigenvectors are sorted from largest to smallest eigenvalues. If the optional parameter `k` is specified, then only the `k` largest eigenvalues and their corresponding eigenvectors will be output. 

The matrix `ain` is first factorized into a tridiagonal matrix using the LAPACK routine `DSYTRD`, and then the eigenvalues are calculated by calls to `DSTEGR` and `DORMTR`.

# See also

[eigvalsym](pyeigvalsym.html), [eigvalvecsymtri](pyeigvalvecsymtri.html)
