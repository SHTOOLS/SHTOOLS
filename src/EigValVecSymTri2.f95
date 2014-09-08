subroutine EigValVecSymTri(ain, n, eig, evec, ul)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	This subroutine will return the eigenvalues and eigenvectors
!	of the symmetric square tridiagonal matrix Ain. The output eigenvectors
!	are ordered from greatest to least, and the norm of the eigenvectors
!	is unity. The sign convention for the eigenvectors is that the 
!	first value of each eigenvector is positive.
!
!	Calling Parameters
!		IN
!			Ain	Input symmetric tridiagonal matrix.
!			n	Order of the matrix Ain.
!		OUT
!			eig	Vector of length n of the eigenvalues of Ain.
!			evec	Matrix of dimension n of the 
!				eigenvectors of Ain.
!		IN (OPTIONAL)
!			ul	Use the upper 'U' or lower 'L' portion of the 
!				input symmetric matrix. By default, the lower portion 
!				of the matrix will be used.
!
!	Notes:
!		
!	1.	The tridiagonal matrix is reduced to  A = S L S'.
!	2.	If accuracy is a problem, consider changing the value of ABSTOL.
!		I have arbitrarily set this to zero, which forces the routines to 
!		pick a default value (which may or may not be good).
!	3. 	The sign convention for the eigenvalues might want to be changed.
!		For instance, IMSL chooses the sign such that the value with 
!		max(abs(evec)) is positive.
!		
!	Dependencies: LAPACK, BLAS
!
!	Written by Mark Wieczorek June 2004.
!
!	Copyright (c) 2005, Mark A. Wieczorek
!	All rights reserved.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	implicit none
	real*8, intent(in) ::	ain(:,:)
	integer, intent(in) ::	n
	real*8, intent(out) ::	eig(:), evec(:,:)
	character, intent(in), optional :: ul
	integer, parameter ::	nb = 80, nbl = 10
	real*8 ::	d(n), e(n), work(nb*n), vl, vu, abstol, w(n)
	real*8, allocatable ::	z(:, :)
	integer ::	lwork, info, il, iu, m, isuppz(2*n), liwork, iwork(nbl*n), i, astat
	external 	dstegr_
	
	
	if (size(ain(:,1)) < n .or. size(ain(1,:)) < n) then
		print*, "Error --- EigValVecSymTri"
		print*, "AIN must be dimensioned as (N, N) where N is ", n
		print*, "Input array is dimensioned as ", size(ain(:,1)), size(ain(1,:))
		stop
	elseif(size(eig) < n) then
		print*, "Error --- EigValVecSymTri"
		print*, "EIG must be dimensioned as (N) where N is ", n
		print*, "Input array is dimensioned as ", size(eig)
		stop
	elseif(size(evec(:,1)) < n .or. size(evec(1,:)) < n) then
		print*, "Error --- EigValVecSymTri"
		print*, "EVEC must be dimensioned as (N, N) where N is ", n
		print*, "Input array is dimensioned as ", size(evec(:,1)), size(evec(1,:))
		stop

	endif
	
	allocate(z(n, n), stat = astat)
	if (astat /= 0) then
		print*, "Error --- EigValVecSymTri"
		print*, "Problem allocating arrays Z", astat
		stop
	endif
	
	
	lwork = nb*n
	liwork = nbl*n
	
	eig = 0.0d0
	evec = 0.0d0
	
	d(1) = ain(1,1)

	if (present(ul)) then
		if (ul == "U" .or. ul == "u") then
			do i=2, n, 1
				d(i) = ain(i,i)
				e(i-1) = ain(i-1, i)
			enddo
		elseif(ul =="L" .or. ul == "l") then
			do i=2, n, 1
				d(i) = ain(i,i)
				e(i-1) = ain(i, i-1)
			enddo
		else
			print*, "Error --- EigValVecSym"
			print*, "UL must be either U, u, L, or l"
			print*, "Input value is ", ul
			stop
		endif
	else
	
		do i=2, n, 1
			d(i) = ain(i,i)
			e(i-1) = ain(i, i-1)
		enddo
	endif
	
	e(n) = 0.0d0
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!
	!	Factor tridiagonal matric A = S L S', and re-order
	!	eignevalues and vectors from greatest to least.
	!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	abstol = 0.0d0
	
	call dstegr_('v','a', n, d, e, vl, vu, il, iu, abstol, m,  w, &
			z, n, isuppz, work, lwork, iwork, liwork, info)
	
       	if (info /= 0) then
		print*, "Error --- EigValVecSymTri"
		print*, "Problem determining eigenvalues and eigenvectors of tridiagonal matrix."
		if (info==1) print*, "Internal error in DLARRE"
		if (info==2) print*, "Internal error in DLARRV"
		stop
	else
		if (work(1) > dble(lwork) ) then
			print*, "Warning --- EigValVecSymTri"
			print*, "Consider changing value of nb to ", work(1)/n, " and recompile SHTOOLS archive."
		endif
		
		if (iwork(1) > liwork ) then
			print*, "Warning --- EigValVecSymTri"
			print*, "Consider changing value of nb to ", iwork(1)/n, " and recompile SHTOOLS archive."
		endif
	endif

	do i=1, n
		eig(i) = w(n+1-i)
		evec(1:n,i) = z(1:n,n+1-i)
		if (evec(1,i) < 0.0d0) evec(1:n,i) = -evec(1:n,i)
	enddo
	
	deallocate(z)
	
end subroutine EigValVecSymTri

