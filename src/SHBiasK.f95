subroutine SHBiasK(tapers, lwin, numk, incspectra, ldata, outcspectra, taper_wt, save_cg)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	This subroutine will compute the expected multitaper windowed (cross-)power spectra
!	using the first K tpaers of the spherical cap concentration problem for a given input 
!	(cross-)power spectrum. Note that this routine makes the assumption that 
!	the known spectra can be described as a random variable. If the input
!	spectra is a cross-power spectra, then it is assumed that the two
!	fields are linearly related. The default is to apply an equal weight to each expected
!	spectrum, but this can be changed by specifying the optional arguement taper_wt.
!	The sum of taper_wt will be normalized to unity.
!
!	Calling Parameteters
!		IN
!			tapers		The coefficients of the tapers (lwin+1, >=numk) corresponding to the
!					spherical cap concentration problem, where each column corresponds
!					to the coefficients for a given value of m. Note the the exact value
!					of m is not important for this routine.
!			lwin		Maximum spherical harmonic degree of the window coefficients.
!			incspectra	Knonw input (cross-)power spectrum as a function of degree.
!			ldata		Maximum degree of incspectra. Beyong this degree
!					incspectra is assumed to be zero.
!		IN (Optional)
!			taper_wt	Vector of length numk corresponding to the weights applied
!					to each spectal estimate. The sum of taper_wt will be normalized
!					to unity.
!			save_cg		If 1, the Clebsch-Gordon coefficients will be calculated
!					and saved, and then used in all subsequent calls (if lwin
!					and ldata are not changed).
!					If -1, the allocated memory for these terms will be deallocated.
!		OUT	
!			outcspectra	Biassed estimate of the windowed 
!					power spectra. Maximum degree calculated is equal
!					to lwin + ldata, or the dimension of outspectra.
!
!	Dependencies: Wigner3j
!
!	Written by Mark Wieczorek, June 2006.
!
!	Modified May 5 2009 in oder to allow the Clebsch-Gordon coefficients to be precomputed and saved.
!	November 5 (2009). Fixed a major bug that would give wrong answers when SAVE_CG=1.
!
!	Copyright (c) 2006-2009, Mark A. Wieczorek
!	All rights reserved.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!			
	use SHTOOLS, only: Wigner3j
	implicit none
	
	real*8, intent(in) ::	tapers(:,:), incspectra(:)
	real*8, intent(out) ::	outcspectra(:)
	integer, intent(in) ::	lwin, ldata, numk
	real*8, intent(in), optional :: taper_wt(:)
	integer, intent(in), optional :: save_cg
	integer ::	l, i, j, lmax, imin, imax, k, astat
	real*8 :: 	wig(2*lwin+ldata+1)
	real*8, allocatable, save ::	cg2(:,:,:)
	
	lmax = ldata + lwin
	outcspectra = 0.0d0
	
	if (size(tapers(:,1)) < lwin+1 .or. size(tapers(1,:)) < numk) then
		print*, "Error --- SHBiasK"
		print*, "TAPERS must be dimensioned as (LWIN+1, NUMK) where LWIN and NUMK are ", LWIN, NUMK
		print*, "Input array is dimensioned ", size(tapers(:,1)), size(tapers(1,:))
		stop
	elseif (size(incspectra) < ldata +1) then
		print*, "Error --- SHBiasK"
		print*, "INCSPECTRA must be dimensioned as (LDATA+1) where LDATA is ", ldata
		print*, "Input array is dimensioned ", size(incspectra)
		stop
	endif
	
	if (present(taper_wt)) then
		if (size(taper_wt) < numk) then	
			print*, "Error --- SHBiasK"
			print*, "TAPER_WT must be dimensioned as (NUMK) where NUMK is ", numk
			print*, "Input array is dimensioned as ", size(taper_wt)
			stop
		elseif (sum(taper_wt(1:numk)) /= 1.0d0) then
			print*, "Error --- SHBiasK"
			print*, "TAPER_WT must sum to unity."
			print*, "Input array sums to ", sum(taper_wt(1:numk))
			stop
		endif
	endif
	
	if (present(save_cg)) then
		if (save_cg /= 1 .and. save_cg /= -1 .and. save_cg /= 0) then
			print*, "Error --- SHBiasK"
			print*, "SAVE_CG must be 1 (to save the Clebsch-Gordon coefficients), -1 (to deallocate the memory)" 
			print*, "or 0 (to do nothing)."
			print*, "Input value is ", save_cg
			stop
		endif
	endif
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!
	!	Calculate the biassed power spectrum
	!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	if (present(save_cg)) then
		
		if (save_cg == -1) then
		
			if (allocated(cg2)) deallocate(cg2)
			return
			
		elseif (save_cg == 0) then
		
			if (present(taper_wt)) then
				do l=0, min(lmax, size(outcspectra)-1)
					do j=0, lwin
						call Wigner3j(wig, imin, imax, j, l, 0, 0, 0)
						do i=imin, min(imax,ldata), 2		
							do k=1, numk	
								outcspectra(l+1) = outcspectra(l+1) + taper_wt(k) * tapers(j+1,k)**2 * &
									incspectra(i+1) * (2.0d0*l+1.0d0) * wig(i-imin+1)**2	
							enddo	
						enddo
					enddo
				enddo

			else
		
				do l=0, min(lmax, size(outcspectra)-1)
					do j=0, lwin
						call Wigner3j(wig, imin, imax, j, l, 0, 0, 0)
						do i=imin, min(imax,ldata), 2		
							do k=1, numk					
								outcspectra(l+1) = outcspectra(l+1) + tapers(j+1,k)**2 * &
									incspectra(i+1) * (2.0d0*l+1.0d0) * wig(i-imin+1)**2
							enddo	
						enddo
					enddo
				enddo

				outcspectra = outcspectra / dble(numk)
		
			endif
			
			return
			
		endif
		
		if (allocated(cg2) .and. (size(cg2(:,1,1)) /= lmax+1 .or. size(cg2(1,:,1)) /= lwin+1 .or. &
			size(cg2(1,1,:)) /= lmax+lwin+1) ) deallocate(cg2)	
		
		if (.not.allocated(cg2)) then
			allocate(cg2(lmax+1, lwin+1, lmax+lwin+1), stat = astat)
			if (astat /= 0) then
				print*, "Error --- SHBiasK"
				print*, "Problem allocating internal array CG2"
				stop
			endif
			
			cg2 = 0.0d0
			
			do l=0, lmax
				do j=0, lwin
					call Wigner3j(wig, imin, imax, j, l, 0, 0, 0)
					cg2(l+1,j+1,1:imax-imin+1) = (2.0d0*l+1.0d0) * wig(1:imax-imin+1)**2
				enddo
			enddo

			
		endif
		
		if (present(taper_wt)) then
			do l=0, min(lmax, size(outcspectra)-1)
				do j=0, lwin
					imin = abs(j-l)
					imax = j+l
					do i=imin, min(imax,ldata), 2		
						do k=1, numk	
							outcspectra(l+1) = outcspectra(l+1) + taper_wt(k) * tapers(j+1,k)**2 * &
								incspectra(i+1) * cg2(l+1,j+1,i-imin+1)	
						enddo	
					enddo
				enddo
			enddo

		else
	
			do l=0, min(lmax, size(outcspectra)-1)
				do j=0, lwin
					imin = abs(j-l)
					imax = j+l
					do i=imin, min(imax,ldata), 2		
						do k=1, numk					
							outcspectra(l+1) = outcspectra(l+1) + tapers(j+1,k)**2 * &
								incspectra(i+1) * cg2(l+1,j+1,i-imin+1)
						enddo	
					enddo
				enddo
			enddo
		
			outcspectra = outcspectra / dble(numk)
		
		endif
			
	else
	
		if (present(taper_wt)) then
			do l=0, min(lmax, size(outcspectra)-1)
				do j=0, lwin
					call Wigner3j(wig, imin, imax, j, l, 0, 0, 0)
					do i=imin, min(imax,ldata), 2		
						do k=1, numk	
							outcspectra(l+1) = outcspectra(l+1) + taper_wt(k) * tapers(j+1,k)**2 * &
								incspectra(i+1) * (2.0d0*l+1.0d0) * wig(i-imin+1)**2	
						enddo	
					enddo
				enddo
			enddo

		else
	
			do l=0, min(lmax, size(outcspectra)-1)
				do j=0, lwin
					call Wigner3j(wig, imin, imax, j, l, 0, 0, 0)
					do i=imin, min(imax,ldata), 2		
						do k=1, numk					
							outcspectra(l+1) = outcspectra(l+1) + tapers(j+1,k)**2 * &
								incspectra(i+1) * (2.0d0*l+1.0d0) * wig(i-imin+1)**2
						enddo	
					enddo
				enddo
			enddo

			outcspectra = outcspectra / dble(numk)
		
		endif
	
	endif

	
end subroutine SHBiasK


