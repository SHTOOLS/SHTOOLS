subroutine SHReturnTapersM(theta0, lmax, m, tapers, eigenvalues, shannon)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	This subroutine will return all the eigenvalues and 
!	eigenfunctions for the space concentration problem of
!	a spherical cap of angular radius theta0. The returned
!	eigenfunctions correspond to "geodesy" normalized spherical 
!	harmonic coefficients, and the eigenfunctions are further
!	normalized such that they have unit power (i.e., the integral 
!	of the function squared over the sphere divided by 4 pi is 1, and the sum of 
!	the squares of their coefficients is 1).
!
!	Note that the eigenfunctions are calculated using the
!	kernel of Grunbaum et al. 1982, and that the eigenvalues are
!	then calculated using the definition of the the space concentration
!	problem with its corresponding space concentration kerel. This 
!	is done because the eigenfunctions of the former are unreliable 
!	when the there are several eigenvalues
!	identical (to machine precision) to either 1 or zero.
!
!
!	Calling Parameters
!		IN
!			theta0		Angular radius of spherical cap
!					in RADIANS.
!			lmax		Maximum spherical harmonic degree
!					for the concentration problem.
!			m		Angular order of the concentration
!					problem (m=0 corresponds to isotropic
!					case).
!		OUT
!			tapers		An (lmax+1) by (lmax+1) array containing
!					all the eigenfunctions of the space-
!					concentration kernel. Eigenfunctions
!					are listed by columns in decreasing order
!					corresponding to value of their eigenvalue.
!			eigenvalues	A vector of length lmax+1 containing the
!					eigenvalued corresponding to the individual
!					eigenfunctions.
!		OPTIONAL	
!			shannon		Shannon number as calculated from the trace of the 
!					kernel.
!
!	Dependencies: LAPACK, BLAS, ComputeDG82, EigValVecSymTri, PreGLQ, PlmBar
!
!	Written by Mark Wieczorek 2003
!
!	Copyright (c) 2005, Mark A. Wieczorek
!	All rights reserved.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	use SHTOOLS, only: ComputeDG82, EigValVecSymTri, PreGLQ, PlmBar, PlmIndex

	implicit none
		
	real*8, intent(in) ::	theta0
	integer, intent(in) ::	lmax, m
	real*8, intent(out) ::	tapers(:,:), eigenvalues(:)
	real*8, intent(out), optional ::	shannon
	integer ::	l, n, n_int, j, i, astat(3)
	real*8 :: 	eval(lmax+1), pi, upper, lower, zero(lmax+1), w(lmax+1), h
	real*8, allocatable ::	evec(:,:), dllmtri(:,:), p(:)
	
	
	if (size(tapers(:,1)) < (lmax+1) .or. size(tapers(1,:)) < (lmax+1) ) then
		print*, "Error --- SHReturnTapersM"
		print*, "TAPERS must be dimensioned as ( LMAX+1, LMAX+1) where LMAX is ", lmax
		print*, "Input array is dimensioned as ", size(tapers(:,1)), size(tapers(1,:))
		stop
	elseif(size(eigenvalues) < (lmax+1) ) then
		print*, "Error --- SHReturnTapersM"
		print*, "EIGENVALUES must be dimensioned as (LMAX+1) where LMAX is ", lmax
		print*, "Input array is dimensioned as ", size(eigenvalues)
		stop
	elseif (m>lmax) then
		print*, "Error --- SHReturnTapersM"
		print*, "M must be less than or equal to LMAX."
		print*, "M = ", m
		print*, "LMAX = ", lmax
		stop
	endif
	
	allocate(evec(lmax+1, lmax+1), stat = astat(1))
	allocate(dllmtri(lmax+1, lmax+1), stat = astat(2))
	allocate(p((lmax+1)*(lmax+2)/2), stat = astat(3))
	if (astat(1) /= 0 .or. astat(2) /= 0 .or. astat(3) /= 0) then
		print*, "Error --- SHReturnTapersM"
		print*, "Problem allocating arrays EVEC, DLLMTRI, and P", &
			astat(1), astat(2), astat(3)
		stop
	endif

	
	pi = acos(-1.0d0)
	
	tapers = 0.0d0
	eigenvalues = 0.0d0
	eval = 0.0d0
	evec = 0.0d0
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!
	!	Calculate space concentration Kernel, and the
	!	corresponding eigenfunctions of the Grunbaum et al. kernel.
	!	Calculate eigenvalues using concentration criter.
	!	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	n = lmax + 1 - abs(m)
	
	call ComputeDG82(dllmtri(1:n,1:n), lmax, m, theta0)
	
	call EigValVecSymTri(dllmtri(1:n,1:n), n, eval(1:n), evec(1+abs(m):lmax+1,1:n))

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!
	!	Calculate true eigenvalues
	!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		
	upper = 1.0d0
	lower = cos(theta0)
	n_int = lmax+1
	
	call PreGLQ(lower, upper, n_int, zero, w)
	
	do i=1, n_int
		
		call PlmBar(p, lmax, zero(i))
		
		do j=1, lmax+1
		
			h = 0.0d0
			
			do l=abs(m), lmax
				h = h + p(PlmIndex(l, abs(m))) * evec(l+1, j)
			enddo
			
			eigenvalues(j) = eigenvalues(j) + w(i) * h**2
		enddo
	enddo
	
	if (m==0) then
		eigenvalues(1:lmax+1) = eigenvalues(1:lmax+1) / 2.0d0
	else
		eigenvalues(1:lmax+1) = eigenvalues(1:lmax+1) / 4.0d0
	endif
	
	if (present(shannon)) then
		shannon = sum(eigenvalues(1:lmax+1))
	endif
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!
	! 	Normalize eigenvectors. By default, the eigenvectors have 
	!	an L2 norm of 1, which corresponds to a total unit power 
	!	(function^2 integreated over all space / 4pi = 1) when
	!	using the geodesy spherical harmonic normalization convenction.
	!	Modify sign convention of the eigenvectors such that the eigenfunction
	!	has a positive value at the north pole (form m=0 only). 
	!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	if (m == 0) then
	
		do l=0, lmax
			p(l+1) = sqrt(dble(2*l+1))	! values of normalized Legendre Polynomials at north pole
		enddo
	
		do l=1, lmax+1
			!evec(:,l) =  evec(:,l) / sqrt(4.0d0*pi) 	! The integral of the taper squared is 1
			if (dot_product(evec(1:lmax+1,l), p(1:lmax+1)) < 0.0d0) evec(:,l) = -evec(:,l) 
		enddo
		
	endif
		
	tapers(1:lmax+1,1:lmax+1) = evec(1:lmax+1,1:lmax+1)
	
	! deallocate memory
	call PlmBar(p, -1, zero(1))
	deallocate(evec)
	deallocate(dllmtri)
	deallocate(p)
	
end subroutine SHReturnTapersM

