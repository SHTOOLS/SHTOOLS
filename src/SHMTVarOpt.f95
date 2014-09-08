Subroutine SHMTVarOpt(l, tapers, taper_order, lwin, kmax, Sff, var_opt, var_unit, weight_opt, unweighted_covar, nocross)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	Given the first Kmax tapers of a matrix TAPERS, and an input global power spectrum
!	Sff, this subroutine will compute the variance of the multitaper spectral estimate
!	assuming that the weights are equal to 1/kmax, and by using the optimal weights that
!	minimize the variance. This routine only works using the tapers of the spherical cap
!	concentration problem.
!
!	Calling Parameters
!		IN
!			l		Spherical harmonic degree to compute variances.
!			tapers		An array of (lwin+1, kmax) tapers (arranged in columns).
!			taper_order	An array of dimension kmax containing the REAL angular order of the tapers.
!			lwin		Maximum spherical harmonic degree of the bandlimited tapers.
!			kmax 		Maximum number of tapers to be used in making the spectral estimate.
!			Sff		Known global power spectrum of the unwindowed field.
!		OUT
!			var_opt		Array of dimension kmax containing the variance computed using the
!					optimal weights when using the first 1 to kmax tapers.
!			var_unit	Array of dimension kmax containg the variance computed using
!					equal weights when using the first 1 to kmax tapers.
!		OPTIONAL(out)
!			weight_opt	Matrix of dimension (kmax, kmax) containing the numerical values of the
!					optimal weights to be applied to each spectral estimate (arranged by columns).
!			unweighted_covar	Unweighted covariance matrix, Fij
!		OPTIONAL(in)
!			nocross		If present and equal to 1, then the off-diagonal terms of the covariance
!					matrix will be assumed to be zero.
!
!	Dependencies: SHSjkPG, SSYSV (LAPACK)
!
!	Copyright (c) 2006, Mark A. Wieczorek
!	All rights reserved.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	use SHTOOLS, only: SHSjkPG

	implicit none

	real*8, intent(in) ::	tapers(:,:), Sff(:)
	real*8, intent(out) ::	var_opt(:), var_unit(:)
	integer, intent(in) ::	l, lwin, kmax, taper_order(:)
	real*8, intent(out), optional ::	weight_opt(:,:), unweighted_covar(:,:)
	integer, intent(in), optional ::	nocross
	integer, parameter ::	nb = 64
	real*8 ::		Fij(kmax, kmax), ww(kmax), bb(kmax+1), MM(kmax+1, kmax+1), work((kmax+1)*nb)
	integer ::	i,j, m, mp, k, ipiv(kmax+1), info, lwork
	external :: 	DSYSV
	complex*16 ::	temp1
	
	if (size(Sff) < l+lwin + 1) then
		print*, "Error --- SHMTVarOpt"
		print*, "Sff must be dimensioned (L+LWIN+1) where L and LWIN are ", l, lwin
		print*, "Input array is dimensioned ", size(Sff)
		stop
	elseif (size(tapers(:,1)) < lwin+1 .or. size(tapers(1,:)) < kmax) then
		print*, "Error --- SHMTVarOpt"
		print*, "TAPERS must be dimensioned as (LWIN+1, KMAX) where LWIN and KMAX are ", lwin, kmax
		print*, "Input array is dimensioned ", size(tapers(:,1)), size(tapers(1,:))
		stop
	elseif ( size(var_opt) < kmax) then
		print*, "Error --- SHMTVarOpt"
		print*, "VAR_OPT must be dimensioned (KMAX) where KMAX is ",  kmax
		print*, "Input array is dimensioned ",  size(var_opt)
		stop
	elseif ( size(var_unit) < kmax) then
		print*, "Error --- SHMTVarOpt"
		print*, "VAR_UNIT must be dimensioned (KMAX) where KMAX is ",  kmax
		print*, "Input array is dimensioned ",  size(var_unit)
		stop
	elseif (size(taper_order) < kmax) then
		print*, "Error --- SHMTVarOpt"
		print*, "TAPER_ORDER must be dimensioned as (KMAX) where KMAX is ", kmax
		print*, "Input array is dimensioned ", size(taper_order)
		stop
	endif

	if (present(weight_opt)) then
		if (size(weight_opt(:,1)) < kmax .or. size(weight_opt(1,:)) < kmax) then
			print*, "Error --- SHMTVarOpt"
			print*, "WEIGHT_OPT must be dimensioned (KMAX, KMAX) where KMAX is ", kmax
			stop
		endif
	endif
	
	if (present(unweighted_covar)) then
		if (size(unweighted_covar(:,1)) < kmax .or. size(unweighted_covar(1,:)) < kmax) then
			print*, "Error --- SHMTVarOpt"
			print*, "UNWEIGHTED_COVAR must be dimensioned (KMAX, KMAX) where KMAX is ", kmax
			print*, "Input array is dimensioned ", size(unweighted_covar(:,1)), size(unweighted_covar(1,:))
			stop
		endif
	endif
	
	if (present(nocross)) then
		if (nocross .ne. 1 .and. nocross .ne. 0) then
			print*, "Error --- SHMTVarOpt0"
			print*, "NOCROSS must be either 0 (use all elements of covariance matrix) or"
			print*, "1 (set off-diagonal elements of covariance matrix to zero)."
			print*, "Input value is ", nocross
			stop
		endif
	endif

	
	Fij = 0.0d0
	
	lwork = (kmax+1)*nb
	
	var_opt = 0.0d0
	var_unit = 0.0d0
	if (present(weight_opt)) then
		weight_opt = 0.0d0
	endif
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!
	!	Calculate matrix Fij, bb
	!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	if (present(nocross)) then
		
		if (nocross==1) then

			do i=1, kmax, 1
				!print*, i
				do m=-l, l
					do mp = -l, l
						temp1 = SHSjkPG(Sff, l, m, mp, tapers(1:lwin+1,i), tapers(1:lwin+1,i),&
								taper_order(i), taper_order(i), lwin, 1)
						Fij(i,i) = Fij(i,i) + 2.0d0 * temp1*dconjg(temp1)
					enddo
				enddo
			enddo
		
		else
		
			do i=1, kmax, 1
				!print*, i
				do j=i, kmax, 1
					do m=-l, l
						do mp = -l, l
							temp1 = SHSjkPG(Sff, l, m, mp, tapers(1:lwin+1,i), tapers(1:lwin+1,j),&
									taper_order(i), taper_order(j), lwin, 1)
							Fij(i,j) = Fij(i,j) + 2.0d0 * temp1*dconjg(temp1)
						enddo
					enddo
					if (i/=j) Fij(j,i) = Fij(i,j)
				enddo
			enddo
		
		endif
		
	else
	
		do i=1, kmax, 1
			!print*, i
			do j=i, kmax, 1
				do m=-l, l
					do mp = -l, l
						temp1 = SHSjkPG(Sff, l, m, mp, tapers(1:lwin+1,i), tapers(1:lwin+1,j),&
								taper_order(i), taper_order(j), lwin, 1)
						Fij(i,j) = Fij(i,j) + 2.0d0 * temp1*dconjg(temp1)
					enddo
				enddo
				if (i/=j) Fij(j,i) = Fij(i,j)
			enddo
		enddo
		
	endif
	
	if (present(unweighted_covar)) then
		unweighted_covar = 0.0d0
		unweighted_covar(1:kmax, 1:kmax) = Fij(1:kmax,1:kmax)
	endif
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!
	!	Calculate variance using equal weights
	!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	do k=1, kmax
		var_unit(k) = sum(Fij(1:k,1:k)) / dble(k)**2
	enddo

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!
	!	Calculate optimal weights and variance using optimal weights
	!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	do k=1, kmax
		MM(1:k,1:k) = 2.0d0 * Fij(1:k,1:k)
		MM(k+1,1:k) = 1.0d0
		MM(1:k,k+1) = 1.0d0
		MM(k+1,k+1) = 0.0d0
		bb(1:k) = 0.0d0
		bb(k+1) = 1.0d0
		
		call DSYSV("u", k+1, 1, MM, kmax+1, ipiv,  bb,  kmax+1,  work,  lwork, info)
		
		if (info /=0) then
			print*, "Error --- SHMTVarOpt"
			print*, "Problem with call to DSYSV." 
			print*, "DSYSV info = ", info
			stop
		endif
		
		if (work(1)/(k+1) > nb) then
			print*, "Warning --- SHMTVarOpt0"
			print*, "The optimal size of nb is ", work(1)/(k+1)
			print*, "Present size is ", nb
			print*, "Please consider changing this parameter and recompiling SHTOOLS"
		endif
		
		do i=1, k
			ww(i) = bb(i)
		enddo
		
		if (present(weight_opt)) then
			weight_opt(1:k,k) = ww(1:k)
		endif
		
		do i=1, k
			do j=1, k
				var_opt(k) = var_opt(k) + ww(i)*ww(j)*Fij(i,j)
			enddo
		enddo
	
	enddo

end subroutine SHMTVarOpt
