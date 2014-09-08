subroutine SHMTDebias (mtdebias, mtspectra, lmax, tapers, lwin, K, nl, lmid, n, taper_wt)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	This routine will invert for the global power spectrum given a 
!	multitaper spectrum estimate, its associated uncertiainty, the 
!	coefficients of the windows used in the mutitaper analysis, and
!	the wieghts used with each window. This routine will only work using
!	tapers obtained from the spherical cap concentration problem. It is assumed
!	that the global power spectrum is constant in bins of nl. It is furthermore
!	assumed that the global spectrum is constant beyond lmax. The inverse problem
!	is solved by SVD as described in Numerical Recipes (pp. 670-672)
!
!	Written by Mark Wieczorek, September 13, 2006
!
!	Copyright (c) 2006, Mark A. Wieczorek
!	All rights reserved.
!	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	use SHTOOLS, only: wigner3j
	
	implicit none
	real*8, intent(out) ::	mtdebias(:,:), lmid(:)
	real*8, intent(in) :: 	mtspectra(:,:), tapers(:,:)
	real*8, intent(in), optional :: taper_wt(:)
	integer, intent(in) :: 	lmax, K, lwin, nl
	integer, intent(out) :: n
	real*8 ::	w3j(lwin+2*lmax+1), sum1, y(lmax+1), ss(lmax+1)
	integer :: 	i, j, l, wmin, wmax, nstart, nstop, info, lwork, m, astat(5), iwork(8*(lmax+lwin+1))
	real*8, allocatable :: work(:), Mmt(:, :), a(:, :), vt(:, :), uu(:, :)


	n = ceiling( dble(lmax+1)/dble(nl) )
	m = lmax+1
			
	if (size(mtspectra(:,1)) /= 2 .or. size(mtspectra(1,:)) < lmax+1) then
		print*, "Error --- SHMTDebias"
		print*, "MTSPECTRA must be dimensioned as (2, LMAX+1) where LMAX is ", lmax
		print*, "Input array is dimensioned as ", size(mtspectra(:,1)), size(mtspectra(1,:))
		stop
	elseif (size(mtdebias(:,1)) /= 2 .or. size(mtdebias(1,:)) < lmax+1) then
		print*, "Error --- SHMTDebias"
		print*, "MTDEBIAS must be dimensioned as (2, LMAX+1) where LMAX is ", lmax
		print*, "Input array is dimensioned as ", size(mtdebias(:,1)), size(mtdebias(1,:))
		stop
	elseif (size(tapers(:,1)) < lwin+1 .or. size(tapers(1,:)) < K) then
		print*, "Error --- SHMTDebias"
		print*, "TAPERS must be dimensioned as (LWIN+1, K) where LWIN and K are ", lwin, k
		print*, "Input array is dimensioned as ", size(tapers(:,1)), size(tapers(1,:))
		stop
	elseif (size(lmid) < n) then
		print*, "Error --- SHMTDebias"
		print*, "LMID must be dimensioned as N=ceiling(dble(lmax+1)/dble(nl), where N is ", n
		print*, "Input array is dimensioned as ", size(lmid)
		stop
	endif
	
	if (present(taper_wt)) then
		if (size(taper_wt) < K) then
			print*, "Error ---  SHMTDebias"
			print*, "TAPER_WT must be dimensioned as (K) where K is ", k 
			print*, "Input array is dimensioned as ", size(taper_wt)
			stop
		endif
	endif
	
	lwork = 3*min(m,n)*min(m,n) + max(max(m,n),4*min(m,n)*min(m,n)+4*min(m,n))
	
	allocate(Mmt(lmax+1, lmax+lwin+1), stat = astat(1))
	allocate(a(lmax+1, (lmax+2)/nl), stat = astat(2))
	allocate(vt(lmax+lwin+1, lmax+lwin+1), stat = astat(3))
	allocate(uu(lmax+1, lmax+1), stat = astat(4))
	allocate(work(lwork), stat=astat(5))
	
	if (astat(1) /= 0 .or. astat(2) /= 0 .or. astat(3) /= 0 .or. astat(4) /= 0 .or. astat(5) /= 0) then
		print*, "Error --- SHMTDebias"
		print*, "Problem allocating arrays MMT, A, VT, UU and WORK", astat(1), astat(2), astat(3), &
			astat(4), astat(5)
		stop
	endif

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!
	!	Compute coupling matrix, M^(mt)
	!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		
	if (present(taper_wt)) then
	
		do i=0, lmax
			do j=0, lmax+lwin
				call Wigner3j(w3j, wmin, wmax, i, j, 0, 0, 0)
				sum1 = 0.0d0
				do l=wmin, min(wmax,lwin), 2
					sum1 = sum1 + dot_product(taper_wt(1:K), tapers(l+1,1:K)**2) * w3j(l-wmin+1)**2
				enddo
				Mmt(i+1,j+1) = sum1 * dble(2*i+1)
			enddo
		enddo	
		
	else

		do i=0, lmax
			do j=0, lmax+lwin
				call Wigner3j(w3j, wmin, wmax, i, j, 0, 0, 0)
				sum1 = 0.0d0
				do l=wmin, min(wmax,lwin), 2
					sum1 = sum1 + sum(tapers(l+1,1:K)**2) * w3j(l-wmin+1)**2
				enddo
				Mmt(i+1,j+1) = sum1 * dble(2*i+1) / dble(K)
			enddo
		enddo	
	
	endif
	
	! Divide linear equations by their uncertainty.
	
	do i=1, m
		Mmt(i,:) = Mmt(i,:) / mtspectra(2,i)
		y(i) = mtspectra(1,i) / mtspectra(2,i)
	enddo
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!
	!	Compute matrix A by assuming that the global power spectrum is constant
	!	in intervals of deltal.
	!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	a = 0.0d0
	
	do j=1, n
	
		nstart =  1 + (j-1)*nl
		
		if (j==n) then
			nstop = lmax+lwin+1
		else
			nstop =  nstart + nl - 1
		endif
		
		lmid(j) = dble(nstart+nstop)/2.0d0 - 1
		
		do i=1, m
			a(i,j) = sum( Mmt(i,nstart:nstop) )
		enddo
		
	enddo

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!
	!	Do least squares inversion by SVD
	!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		
	call dgesdd_('a', m, n, a, m, ss, uu, m, vt, lmax+lwin+1, work, lwork, iwork, info)
	
	if (info /=0) then
		print*, "Error --- SHMTDebias"
		print*, "Problem with DGESDD, INFO = ", info
		print*, "if INFO = -i, the i-th argument had an illegal value."
		print*, "if INFO > 0:  DBDSDC did not converge, updating process failed."
		stop
	endif
		
	mtdebias = 0.0d0
	
	do i=1, min(m,n)
		mtdebias(1,1:n) = mtdebias(1,1:n) + dot_product(uu(1:m,i), y(1:m)) * vt(i, 1:n) / ss(i)
	enddo
	
	do j=1, n
		do i=1, min(m,n)
			mtdebias(2,j) = mtdebias(2,j) + (vt(i,j)/ss(i))**2
		enddo
	enddo
	
	mtdebias(2,1:n) = sqrt(mtdebias(2,1:n))
	 
	deallocate(Mmt)
	deallocate(a)
	deallocate(vt)
	deallocate(uu)
	deallocate(work)
		
end subroutine SHMTDebias