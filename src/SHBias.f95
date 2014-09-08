subroutine SHBias(Shh, lwin, incspectra, ldata, outcspectra, save_cg)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	This subroutine will compute the expected windowed (cross-)power spectra
!	given the power spectrum of an arbitrary window, and the "known" input (cross-)power
!	spectra. Note that this routine makes the assumption that 
!	the known spectra can be described as a random variable. 
!
!	Calling Parameteters
!		IN
!			Shh		Window power spectrum.
!			lwin		Maximum spherical harmonic degree of the window.
!			incspectra	Knonw input (cross-)power spectrum as a function of degree.
!			ldata		Maximum degree of incspectra. Beyong this degree
!					incspectra is assumed to be zero.
!		OUT	
!			outcspectra	Biassed estimate of the windowed 
!					power spectra. Maximum degree calculated is equal
!					to lwin + ldata, or the dimension of outcspectra.
!		IN, OPTIONAL
!			save_cg		If 1, the Clebsch-Gordon coefficients will be calculated
!					and saved, and then used in all subsequent calls (if lwin
!					and ldata are not changed).
!					If -1, the allocated memory for these terms will be deallocated.
!
!	Dependencies: Wigner3j
!
!	Written by Mark Wieczorek, August 2004.
!
!	Modified May 28, 2006 so that one inputs the window's power spectrum instead
!	of the window's coefficients. Also, the routine now assumes that incspectra is zero
!	beyond ldata.
!
!	Modified May 5 2009 in oder to allow the Clebsch-Gordon coefficients to be precomputed and saved.
!	
!
!	Copyright (c) 2005-2009, Mark A. Wieczorek
!	All rights reserved.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!			
	use SHTOOLS, only: Wigner3j
	implicit none
	
	real*8, intent(in) ::	Shh(:), incspectra(:)
	real*8, intent(out) ::	outcspectra(:)
	integer, intent(in) ::	lwin, ldata
	integer, intent(in), optional :: save_cg
	integer ::	l, i, j, lmax, imin, imax, astat
	real*8 ::	wig(2*lwin+ldata+1)
	real*8, allocatable, save ::	cg2(:,:,:)
	
	lmax = ldata + lwin
	outcspectra = 0.0d0
	
	if (size(Shh) < lwin+1) then
		print*, "Error --- SHBias"
		print*, "SHH must be dimensioned as (LWIN+1) where LWIN is ", lwin
		print*, "Input array is dimensioned as ", size(Shh)
		stop
	elseif (size(incspectra) < ldata+1) then
		print*, "Error --- SHBias"
		print*, "INCSPECTRA must be dimensioned as (LDATA+1) where LDATA is ", ldata
		print*, "Input array is dimensioned as ", size(incspectra)
		stop
	endif
	
	if (present(save_cg)) then
		if (save_cg /= 1 .and. save_cg /= -1 .and. save_cg /= 0) then
			print*, "Error --- SHBias"
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
		
			do l=0, min(lmax, size(outcspectra)-1)
				do j=0, lwin
					call Wigner3j(wig, imin, imax, j, l, 0, 0, 0)
					do i=imin, min(imax, ldata), 2								
						outcspectra(l+1) = outcspectra(l+1) + Shh(j+1) * &
							incspectra(i+1) * (2.0d0*l+1.0d0) * wig(i-imin+1)**2
					enddo
				enddo
			enddo
			
			return
			
		endif
		
		if (allocated(cg2) .and. (size(cg2(:,1,1)) /= lmax+1 .or. size(cg2(1,:,1)) /= lwin+1 .or. &
			size(cg2(1,1,:)) /= lmax+lwin+1) ) deallocate(cg2)	
		
		if (.not.allocated(cg2)) then
			allocate(cg2(lmax+1, lwin+1, lmax+lwin+1), stat = astat)
			if (astat /= 0) then
				print*, "Error --- SHBias"
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
		
		do l=0, min(lmax, size(outcspectra)-1)
			do j=0, lwin	
				imin = abs(j-l)
				imax = j+l
				do i=imin, min(imax, ldata), 2								
					outcspectra(l+1) = outcspectra(l+1) + Shh(j+1) * &
						incspectra(i+1) * cg2(l+1,j+1,i-imin+1)
				enddo
			enddo
		enddo

			
	else
	
		do l=0, min(lmax, size(outcspectra)-1)
			do j=0, lwin
				call Wigner3j(wig, imin, imax, j, l, 0, 0, 0)
				do i=imin, min(imax, ldata), 2								
					outcspectra(l+1) = outcspectra(l+1) + Shh(j+1) * &
						incspectra(i+1) * (2.0d0*l+1.0d0) * wig(i-imin+1)**2
				enddo
			enddo
		enddo
	
	endif
	
end subroutine SHBias


