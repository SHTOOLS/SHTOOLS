subroutine EigValSym(ain, n, eval, ul)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	This subroutine will return the eigenvalues of the symmetric square 
!	matrix Ain, ordered from greatest to least.
!
!	Calling Parameters
!		IN
!			Ain	Input symmetric matrix. By default, only the
!				upper portion is used.
!			n	Order of the matrix Ain.
!		OUT
!			eval	Vector of length n of the eigenvalues of Ain.
!		OPTIONAL
!			uplo	Use the upper 'U' or lower 'L' portion of the 
!				input symmetric matrix.
!
!	Notes:
!		
!	1.	The eigenvalues and eigenvectors are determined by reducing the matrix to
!			A = Z L Z = Q (S L S') Q' 
!		by the two operations:
!
!		(1) The real symmetric square matrix is reduced to tridiagonal form
!			A = Q T Q'
!		where Q is orthogonal, and T is symmetric tridiagonal.
!		(2) The tridiagonal matrix is reduced to 
!			T = S L S'
!
!		The eigenvalues of A correspond to L (which is a diagonal)
!		
!	Dependencies: LAPACK, BLAS
!
!	Written by Mark Wieczorek June 2004.
!
!	Copyright (c) 2005,2006 Mark A. Wieczorek
!	All rights reserved.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	implicit none
	real*8, intent(in) ::	ain(:,:)
	integer, intent(in) ::	n
	real*8, intent(out) ::	eval(:)
	character, intent(in), optional :: ul
	integer, parameter ::	nb = 80, nbl = 10
	character ::	uplo
	real*8 ::	d(n), e(n), tau(n-1), work(nb*n), vl, vu, &
			abstol, w(n)
	real*8, allocatable ::	a(:,:), z(:,:)
	integer ::	lwork, info, il, iu, m, isuppz(2*n), liwork, iwork(nbl*n), i, astat(2)
	external	dsytrd_, dstegr_
	
	
	if (size(ain(:,1)) < n .or. size(ain(1,:)) < n) then
		print*, "Error --- EigValSym"
		print*, "AIN must be dimensioned as (N, N) where N is ", n
		print*, "Input array is dimensioned as ", size(ain(:,1)), size(ain(1,:))
		stop
	elseif(size(eval) < n) then
		print*, "Error --- EigValSym"
		print*, "EVAL must be dimensioned as (N) where N is ", n
		print*, "Input array is dimensioned as ", size(eval)
		stop
	endif
		
	allocate(a(n,n), stat = astat(1))
	allocate(z(n,n), stat = astat(2))
	
	if (astat(1) /= 0 .or. astat(2) /= 0) then
		print*, "Error --- EigValSym"
		print*, "Problem allocating arrays A and Z", astat(1), astat(2)
		stop
	endif
	
	lwork = nb*n
	liwork = nbl*n
	
	eval = 0.0d0
	a(1:n,1:n) = ain(1:n,1:n)
	
	if (present(ul)) then
		uplo = ul
	else
		uplo = "U"
	endif
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!
	!	Factor A = Q T Q'
	!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	call dsytrd_(uplo, n, a, n, d, e(1:n-1), tau, work, lwork, info)

	if (info /= 0) then
		print*, "Error --- EigValSym"
		print*, "Problem tri-diagonalizing input matrix"
		stop
	else
		if ( work(1) > dble(lwork) ) then
			print*, "Warning --- EigValSym"
			print*, "Consider changing value of nb to ", work(1)/n, " and recompile the SHTOOLS archive."
		endif
	endif
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!
	!	Factor T = S L S'
	!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	abstol = 0.0d0
	
	call dstegr_('n','a', n, d, e, vl, vu, il, iu, abstol, m,  w, &
			z, n, isuppz, work, lwork, iwork, liwork, info)

       	if (info /= 0) then
		print*, "Error --- EigValSym"
		print*, "Problem determining eigenvalues and eigenvectors of tridiagonal matrix."
		if (info==1) print*, "Internal error  in  DLARRE"
		if (info==2) print*, "Internal error in DLARRV"
		stop
	else
		if (work(1) > dble(lwork) ) then
			print*, "Warning --- EigValSym"
			print*, "Consider changing value of nb to ", work(1)/n, " and recompile the SHTOOLS archive."
		endif
		
		if (iwork(1) > liwork ) then
			print*, "Warning --- Eigsym"
			print*, "Consider changing value of nb to ", iwork(1)/n, " and recompile the SHTOOLS archive."
		endif
	endif

	! Reorder eigenvalues from greatest to least.
	
	do i=1, n
		eval(i) = w(n+1-i)
	enddo
	
	deallocate(a)
	deallocate(z)

end subroutine EigValSym

