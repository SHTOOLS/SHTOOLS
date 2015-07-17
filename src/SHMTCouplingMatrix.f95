subroutine SHMTCouplingMatrix(Mmt, lmax, tapers, lwin, K, taper_wt)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       This routine returns the multitaper coupling matrix, weighted by
!       the standard deviation
!
!	Written by Mark A. Wieczorek and Matthias Meschede
!
!	Copyright (c) 2015, Mark A. Wieczorek
!	All rights reserved.
!	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	use SHTOOLS, only: wigner3j
	
	implicit none
	real*8, intent(out) :: 	Mmt(:,:)
	real*8, intent(in) :: 	tapers(:,:)
	real*8, intent(in), optional :: taper_wt(:)
	integer, intent(in) :: 	lmax, K, lwin
	real*8 ::	w3j(lwin+2*lmax+1), sum1, ss(lmax+1)
	integer :: 	i, j, l, wmin, wmax, nstart, nstop, info, lwork, m, astat(5), iwork(8*(lmax+lwin+1))


	m = lmax+1
			
	if (size(tapers(:,1)) < lwin+1 .or. size(tapers(1,:)) < K) then
		print*, "Error --- SHMTDebias"
		print*, "TAPERS must be dimensioned as (LWIN+1, K) where LWIN and K are ", lwin, k
		print*, "Input array is dimensioned as ", size(tapers(:,1)), size(tapers(1,:))
		stop
	endif
	
	if (present(taper_wt)) then
		if (size(taper_wt) < K) then
			print*, "Error ---  Coupling Matrix"
			print*, "TAPER_WT must be dimensioned as (K) where K is ", k 
			print*, "Input array is dimensioned as ", size(taper_wt)
			stop
		endif
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
	
end subroutine SHMTCouplingMatrix
