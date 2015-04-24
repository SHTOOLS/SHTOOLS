# EigValSym

Compute the eigenvalues of a real symmetric matrix.

# Usage

`eval` = pyshtools.EigValSym (`ain`, [`n`, `ul`])

# Returns

`eval` : float, dimension (`n`)
:   The eigenvalues of `ain`, sorted from largest to smallest.

# Parameters

`ain` : float, dimension (`nin`, `nin`)
:   The input real symmetric matrix. By default, only the upper portion of the matrix is used.
	
`n` : optional, integer, default = `nin`
:   The rank of the matrix `ain`.

`ul` : optional, input, character, default = `U`
:   If `U` then the upper portion of the matrix `ain` will be used (default). If `L` then the lower portion of the matrix `ain` will be used.

# Description

`EigValSym` will calculate the eigenvalues of a real symmetric matrix. By default, only the upper portion of the matrix is used, but this can be changed by the optional argument `ul`. The eigenvalues are sorted from largest to smallest. The matrix `ain` is first factorized into a tridiagonal matrix using the LAPACK routine `DSYTRD`, and then the eigenvalues are calculated by a call to `DSTEGR`.

# See also

[eigvalvecsym](pyeigvalvecsym.html), [eigvalvecsymtri](pyeigvalvecsymtri.html)
