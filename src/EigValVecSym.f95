subroutine EigValVecSym(ain, n, eig, evec, ul, K)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	This subroutine will return the eigenvalues and eigenvectors
!	of the symmetric square matrix Ain. The output eigenvectors
!	are ordered from greatest to least, and the norm of the eigenvectors
!	is unity. If the optional parameter K is specified, only the K largest
!	eigenvalues and corresponding vectors will be output.
!
!	Calling Parameters
!		IN
!			Ain	Input symmetric matrix. By default, only the
!				upper portion is used.
!			n	Order of the matrix Ain.
!		OUT
!			eig	Vector of length n of the eigenvalues of Ain.
!			evec	Matrix of dimension n of the 
!				eigenvectors of Ain.
!		OPTIONAL
!			ul	Use the upper 'U' or lower 'L' portion of the 
!				input symmetric matrix.
!			K	The K largest eigenvalues and corresponding eigenvectors
!				to calculate and output.
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
!		The eigenvalues of A correspond to the L (which is a diagonal), and the
!		eigenvectors correspond to Z = Q S.
!	2. 	The sign convention for the eigenvalues might want to be changed.
!		For instance, IMSL chooses the sign such that the value with 
!		max(abs(evec)) is positive.
!		
!	Dependencies: LAPACK, BLAS
!
!	Written by Mark Wieczorek June 2004.
!	Modified August 15:	Option to only calcuate K largest eigenvectors and corresponding eigenvalues. Also,
!				with LAPACK3.1, the DSTEGR routine no longer uses ABSTOL, so this parameter does
!				not need to be specified.
!
!	Copyright (c) 2005-2009, Mark A. Wieczorek
!	All rights reserved.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	implicit none
	real*8, intent(in) ::		ain(:,:)
	integer, intent(in) ::		n
	real*8, intent(out) ::		eig(:), evec(:,:)
	character, intent(in), optional ::	ul
	integer, intent(in), optional ::	K
	integer, parameter ::			nb = 80, nbl = 10
	character ::			uplo
	real*8 ::			d(n), e(n), tau(n-1), work(nb*n), vl, vu, &
					abstol, w(n)
	real*8, allocatable :: 		a(:,:), z(:,:)
	integer ::			lwork, info, il, iu, m, isuppz(2*n), liwork, iwork(nbl*n), &
					i, astat(2)
	external	 		dsytrd, dstegr, dormtr
	
	
	if (size(ain(:,1)) < n .or. size(ain(1,:)) < n) then
		print*, "Error --- EigValVecSym"
		print*, "AIN must be dimensioned as (N, N) where N is ", n
		print*, "Input array is dimensioned as ", size(ain(:,1)), size(ain(1,:))
		stop
	endif
	
	if (present(K)) then
		if (K > n .or. K < 1) then
			print*, "Error --- EigValVecSym"
			print*, "The number of eigenvalues to output must be between 0 and N."
			print*, "N = ", n
			print*, "K = ", k
			stop
		endif
		if(size(eig) < K) then
			print*, "Error --- EigValVecSym"
			print*, "EIG must be dimensioned as (K) where K is ", K
			print*, "Input array is dimensioned as ", size(eig)
			stop
		elseif(size(evec(:,1)) < n .or. size(evec(1,:)) < K) then
			print*, "Error --- EigValVecSym"
			print*, "EVEC must be dimensioned as (N, K)."
			print*, "N = ", n
			print*, "K = ", k
			print*, "Input array is dimensioned as ", size(evec(:,1)), size(evec(1,:))
			stop
		endif
	else
		if(size(eig) < n) then
			print*, "Error --- EigValVecSym"
			print*, "EIG must be dimensioned as (N) where N is ", n
			print*, "Input array is dimensioned as ", size(eig)
			stop
		elseif(size(evec(:,1)) < n .or. size(evec(1,:)) < n) then
			print*, "Error --- EigValVecSym"
			print*, "EVEC must be dimensioned as (N, N) where N is ", n
			print*, "Input array is dimensioned as ", size(evec(:,1)), size(evec(1,:))
			stop
		endif
	endif

	allocate(a(n,n), stat = astat(1))
	allocate(z(n,n), stat = astat(2))
	
	if (astat(1) /= 0 .or. astat(2) /= 0) then
		print*, "Error --- EigValVecSym2"
		print*, "Problem allocating arrays A and Z", astat(1), astat(2)
		stop
	endif

	
	lwork = nb*n
	liwork = nbl*n
	
	eig = 0.0d0
	evec = 0.0d0
	a(1:n,1:n) = ain(1:n,1:n)
	
	if (present(ul)) then
		uplo = ul
	else
		uplo = "U"
	endif
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!
	!	Factor A to Q T Q' where T is a tridiagonal matrix.
	!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	call dsytrd(uplo, n, a, n, d, e(1:n-1), tau, work, lwork, info)

	if (info /= 0) then
		print*, "Error --- EigValVecSym"
		print*, "Problem tri-diagonalizing input matrix"
		print*, "DSYTRD info = ", info
		stop
	else
		if ( work(1) > dble(lwork) ) then
			print*, "Warning --- EigValVecSym"
			print*, "Consider changing value of nb to ", work(1)/n, " and recompile."
		endif
	endif
		
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!
	!	Factor T to S L S' where L is a diagonal matrix.
	!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	abstol = 0.0d0
	
	if (present(K)) then
		call dstegr('v','I', n, d, e, vl, vu, n-K+1, n, abstol, m,  w, &
			z, n, isuppz, work, lwork, iwork, liwork, info)
	else
		call dstegr('v','I', n, d, e, vl, vu, il, iu, abstol, m,  w, &
			z, n, isuppz, work, lwork, iwork, liwork, info)
	endif
	
       	if (info /= 0) then
		print*, "Error --- EigValVecSym"
		print*, "Problem determining eigenvalues and eigenvectors of tridiagonal matrix."
		if (info==1) print*, "Internal error in DLARRE"
		if (info==2) print*, "Internal error in DLARRV"
		print*, "DSTEGR info = ", info
		stop
	else
		if (work(1) > dble(lwork) ) then
			print*, "Warning --- EigValVecSym"
			print*, "Consider changing value of nb to ", work(1)/n, " and recompile SHTOOLS archive."
		endif
		
		if (iwork(1) > liwork ) then
			print*, "Warning --- EigValVecSym"
			print*, "Consider changing value of nb to ", iwork(1)/n, " and recompile SHTOOLS archive."
		endif
	endif

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!
	!	Determine eigenvalues Z = Q S (note that Q is stored in a
	!	bizarre manner, see LAPACK notes), and reorder eigenvalues and
	!	eigenvectors from greatest to least.
	!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	call dormtr('L', uplo, 'N', n, n, a, n, tau, z, n, work, lwork, info)

	if (info /= 0) then
		print*, "Error --- EigValVecSym"
		print*, "Problem multiplying matrices."
		print*, "DORMTR info = ", info
		stop
	else
		if ( work(1) > dble(lwork) ) then
			print*, "Warning --- EigValVecSym"
			print*, "Consider changing value of nb to ", work(1)/n, " and recompile."
		endif
	endif

	if (present(k)) then
		do i=n-K+1, n
			eig(i-n+k) = w(n+1-i)
			evec(1:n,i-n+k) = z(1:n,n+1-i)
		enddo
	else
		do i=1, n
			eig(i) = w(n+1-i)
			evec(1:n,i) = z(1:n,n+1-i)
		enddo
	endif
	
	deallocate(a)
	deallocate(z)
	
end subroutine EigValVecSym

