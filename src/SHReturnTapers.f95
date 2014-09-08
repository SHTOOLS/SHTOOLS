subroutine SHReturnTapers(theta0, lmax, tapers, eigenvalues, taper_order)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	This subroutine will return all the eigenvalues and 
!	eigenfunctions for the space concentration problem of
!	a spherical cap of angular radius theta0. The returned
!	eigenfunctions correspond to "geodesy" normalized spherical 
!	harmonic coefficients, and the eigenfunctions are further
!	normalized such that they have unit power (i.e., the integral 
!	of the function squared over the sphere divided by 4 pi is 1, 
!	and the sum of the squares of their coefficients is 1). 
!	The eigenfunctions are calculated using the
!	kernel of Grunbaum et al. 1982, and the eigenvalues are 
!	calculated by numerical integration. 
!
!
!	Calling Parameters
!		IN
!			theta0		Angular radius of spherical cap
!					in RADIANS.
!			lmax		Maximum spherical harmonic degree
!					for the concentration problem.
!		OUT
!			tapers		An (lmax+1) by (lmax+1) array containing
!					all the eigenfunctions of the space-
!					concentration kernel. Eigenfunctions
!					are listed by columns in decreasing order
!					corresponding to value of their eigenvalue.
!			eigenvalues	A vector of length lmax+1 containing the
!					eigenvalued corresponding to the individual
!					eigenfunctions.
!			taper_order	A vector of dimension X denoting which order m
!					corresponds to the column of tapers and 
!					eigenvalues.
!
!	Dependencies: LAPACK, BLAS, ComputeDG82, EigValVecSymTri, PreGLQ, PlmBar, PlmIndex
!
!	Written by Mark Wieczorek 2006
!
!	Copyright (c) 2005, Mark A. Wieczorek
!	All rights reserved.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	use SHTOOLS, only: ComputeDG82, EigValVecSymTri, PreGLQ, PlmBar, PlmIndex

	implicit none
		
	real*8, intent(in) ::	theta0
	integer, intent(in) ::	lmax
	real*8, intent(out) ::	tapers(:,:), eigenvalues(:)
	integer, intent(out) ::	taper_order(:)
	integer ::		l, m, nt, nt2, n, i, j, jj(1), astat(8), n_int
	real*8, allocatable :: 	eval(:), evec(:, :), tapers_unordered(:,:), &
				dllmtri(:,:), eval_unordered(:), p(:)
	real*8 :: 		pi, upper, lower, w(lmax+1), zero(lmax+1), h
	integer, allocatable ::	m_unordered(:)
	
	
	if (size(tapers(:,1)) < (lmax+1) .or. size(tapers(1,:)) < (lmax+1)**2) then
		print*, "Error --- SHReturnTapers"
		print*, "TAPERS must be dimensioned as ( LMAX+1, (LMAX+1)**2 ) where LMAX is ", lmax
		print*, "Input array is dimensioned as ", size(tapers(:,1)), size(tapers(1,:))
		stop
	elseif(size(eigenvalues) < (lmax+1)**2) then
		print*, "Error --- SHReturnTapers"
		print*, "EIGENVALUES must be dimensioned as (LMAX+1)**2 where LMAX is ", lmax
		print*, "Input array is dimensioned as ", size(eigenvalues)
		stop
	elseif(size(taper_order) < (lmax+1)**2) then
		print*, "Error --- SHReturnTapers"
		print*, "TAPER_ORDER must be dimensioned as (LMAX+1)**2 where LMAX is ", lmax
		print*, "Input array is dimensioned as ", size(taper_order)
		stop
	endif

	
	allocate(eval(lmax+1), stat = astat(1))
	allocate(evec(lmax+1, lmax+1), stat = astat(2))
	allocate(tapers_unordered(lmax+1, (lmax+1)**2), stat = astat(3))
	allocate(dllmtri(lmax+1, lmax+1), stat = astat(4))
	allocate(eval_unordered((lmax+1)**2), stat = astat(5))
	allocate(m_unordered((lmax+1)**2), stat = astat(6))
	allocate(p((lmax+1)*(lmax+2)/2), stat = astat(7))
	
	if (astat(1) /= 0 .or. astat(2) /= 0 .or. astat(3) /= 0 .or. astat(4) /=0 &
		.or. astat(5) /= 0 .or. astat(6) /=0 .or. astat(7) /=0 ) then
		print*, "ERROR ---- SHReturnTapers"
		print*, "Problem allocating memory for temporary arrays", astat
		stop
	endif
	
	pi = acos(-1.0d0)
	
	tapers = 0.0d0
	tapers_unordered = 0.0d0
	eigenvalues = 0.0d0
	eval_unordered = 0.0d0
	m_unordered = 0
	taper_order = 0
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!
	!	Calculate eigenfunctions of the Grunbaum et al. kernel.
	!	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	nt = 0
	
	do m=0, lmax
		
		n = lmax + 1 - m
		
		call ComputeDG82(dllmtri(1:n,1:n), lmax, m, theta0)
		call EigValVecSymTri(dllmtri(1:n, 1:n), n, eval(1:n), evec(1:n,1:n))
		
		tapers_unordered(m+1:lmax+1, nt+1:nt+n) = evec(1:n,1:n)
		m_unordered(nt+1:nt+n) = m
		
		nt = nt + n
	enddo
		
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!
	!	Calculate true eigenvalues
	!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	upper = 1.0d0
	lower = cos(theta0)
	n_int = lmax+1
	
	call PreGLQ(lower, upper, n_int, zero, w)
	
	do i=1, n_int
	
		call PlmBar(p, lmax, zero(i))
		
		do j=1, nt
		
			h = 0.0d0
			
			do l=m_unordered(j), lmax
				h = h + p(PlmIndex(l, m_unordered(j))) * tapers_unordered(l+1, j)
			enddo
			
			eval_unordered(j) = eval_unordered(j) + w(i) * h**2
		enddo
	enddo
	
	do j=1, nt
		if (m_unordered(j) == 0) then
			eval_unordered(j) = eval_unordered(j) / 2.0d0
		else
			eval_unordered(j) = eval_unordered(j) / 4.0d0
		endif
	enddo
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!
	!	Reorder tapers and eigenvalues
	!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	nt2 = 0
	
	do i=1, nt
	
		jj = maxloc(eval_unordered(1:nt))
		j = jj(1)
		
		if (m_unordered(j) == 0) then
			nt2 = nt2 + 1
			taper_order(nt2) = m_unordered(j)
			eigenvalues(nt2) = eval_unordered(j)
			tapers(1:lmax+1,nt2) = tapers_unordered(1:lmax+1,j)
		else
			nt2 = nt2 + 1
			taper_order(nt2) = -m_unordered(j)
			eigenvalues(nt2) = eval_unordered(j)
			tapers(1:lmax+1,nt2) = tapers_unordered(1:lmax+1,j)
			nt2 = nt2 + 1
			taper_order(nt2) = m_unordered(j)
			eigenvalues(nt2) = eval_unordered(j)
			tapers(1:lmax+1,nt2) = tapers_unordered(1:lmax+1,j)
		endif
		
		eval_unordered(j) = -1.0d25
	
	enddo
	
	call PlmBar(p, -1, zero(1))
	deallocate(eval) ;  deallocate(evec) ; deallocate(tapers_unordered)
	deallocate(dllmtri) ; deallocate(eval_unordered)
	deallocate(m_unordered) ; deallocate(p)
					
end subroutine SHReturnTapers
