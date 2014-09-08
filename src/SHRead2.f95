subroutine SHRead2(filename, cilm, lmax, gm, r0_pot, error, dot, doystart, doyend, epoch)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	This program will read the file format of spherical harmonic
!	coefficients used by the CHAMP/GRACE group into standard arrays for
!	of Cilm and the error. Not all data for all fields are returned.
!
!	Calling Parameters
!		IN
!			filename: 	The name of the file.
!		OUT
!			lmax:		Maximum spherical harmonic degree.
!			cilm:		An array of the spherical harmonic coefficients.
!			gm:		GM
!			R0_pot		Reference radius of gravity field
!		OPTIONAL
!			error:		An array containing the error of the cilm coefficients.
!			dot:		An array containing the assumed derivatives in time.
!					of the low degree spherical harmonic coefficients.
!			doy1:		Beginning date of the gravity solution.
!			doy2:		Ending date of the gravity solution.
!			epoch:		Epochal date for the time derivatives.
!
!	Written by Mark Wieczorek (September 2005). Tweaked August (2009) to read "dot" terms, among others
!
!	Copyright (c) 2005-2009, Mark A. Wieczorek
!	All rights reserved.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	implicit none
	
	character(*), intent(in) ::	filename
	integer, intent(out) ::	lmax
	real*8, intent(out) ::	cilm(:,:,:), gm, r0_pot
	real*8, intent(out  ), optional :: error(:,:,:), dot(:,:,:), doystart, doyend, epoch
	integer ::	temp_max, mmax, maxl_cilm, maxl_error, maxl_dot
	real*8 ::	error_scale, doy1, doy2, epoch_temp
	character*114 ::	comment
		
	integer:: l, m, stat
	character :: c*6
	real*8 ::	clm, slm, clmsig, slmsig
	
	cilm=0.0d0
	
	doy1 = 0.0d0
	doy2 = 0.0d0
	epoch_temp = 0.0d0
	
	if (size(cilm(:,1,1)) < 2) then
		print*, "Error --- SHRead2"
		print*, "CILM must be dimensioned (2, *, *)."
		print*, "Input array is dimensioned ", size(cilm(:,1,1)), size(cilm(1,:,1)), size(cilm(1,1,:))
		stop
	endif
	
	maxl_cilm = min(size(cilm(1,:,1)) - 1, size(cilm(1,1,:)) - 1)
	
	if (present(error)) then
		if (size(error(:,1,1)) < 2) then
			print*, "Error --- SHRead2"
			print*, "ERROR must be dimensioned (2, *, *)."
			print*, "Input array is dimensioned ", size(error(:,1,1)), size(error(1,:,1)), size(error(1,1,:))
			stop
		endif
		
		error=0.0d0
		maxl_error = min(size(error(1,:,1)) - 1, size(error(1,1,:)) - 1)
	endif
	
	if (present(dot)) then
		if (size(dot(:,1,1)) < 2) then
			print*, "Error --- SHRead2"
			print*, "DOT must be dimensioned (2, *, *)."
			print*, "Input array is dimensioned ", size(dot(:,1,1)), size(dot(1,:,1)), size(dot(1,1,:))
			stop
		endif
		
		dot=0.0d0
		maxl_dot = min(size(dot(1,:,1)) - 1, size(dot(1,1,:)) - 1)
	endif

	
	open(13,file=filename)
	
	lmax = 0
	temp_max = 0
	gm = 0.0d0
	r0_pot = 0.0d0
	error_scale = 1.0d0
	
	do 
		read(13,"(a6)", advance="no", iostat=stat) c
		if(stat /= 0) exit
		
		if (c=="EARTH " .or. c(1:3) =="GGM") then
		
			read(13,*) gm, r0_pot
			
		elseif(c=="SHM   ") then 
		
			read(13,*) lmax, mmax, error_scale
			if (mmax /= lmax) then
				print*, "Warning --- SHREAD2"
				print*, "Maximum order is not equal to the maximum degree."
				print*, "LMAX = ", lmax
				print*, "MMax = ", mmax
			endif
			if (present(error) .and. error_scale /= 1.0d0) then
				print*, "Warning --- SHREAD2"
				print*, "Standard deviations will be scale by ", error_scale
			endif
			if (size(cilm(1,:,1)) < lmax+1 .or. size(cilm(1,1,:)) < lmax+1) then
				print*, "Error ---SHRead2"
				print*, "CILM must be dimensioned (2, LMAX+1, LMAX+1) where LMAX is ", lmax
				print*, "Input array is dimensioned ", size(cilm(1,:,1)),  size(cilm(1,1,:)) 
				stop
			elseif (present(error)) then
				if (size(error(1,:,1)) < lmax+1 .or. size(error(1,1,:)) < lmax+1) then
					print*, "Error ---SHRead2"
					print*, "ERROR must be dimensioned (2, LMAX+1, LMAX+1) where LMAX is ", lmax
					print*, "Input array is dimensioned ", size(error(1,:,1)),  size(error(1,1,:)) 
					stop
				endif
			endif
			
		elseif(c=="GRCOF2") then
		
			read(13,*) l, m, clm, slm, clmsig, slmsig, doy1, doy2
			if (l > maxl_cilm) then
				print*, "Error ---SHRead2"
				print*, "CILM must be dimensioned (2, LMAX+1, LMAX+1) where LMAX is ", l
				print*, "Input array is dimensioned ", size(cilm(1,:,1)),  size(cilm(1,1,:)) 
				stop
			endif
			
			cilm(1,l+1,m+1) = clm
			cilm(2,l+1,m+1) = slm
			
			if (present(error)) then
				if(l > maxl_error) then
					print*, "Error ---SHRead2"
					print*, "ERROR must be dimensioned (2, LMAX+1, LMAX+1) where LMAX is ", l
					print*, "Input array is dimensioned ", size(error(1,:,1)),  size(error(1,1,:)) 
					stop
				endif
				
				if (error_scale /= 1.0d0) then
					error(1,l+1,m+1) = clmsig * error_scale
					error(2,l+1,m+1) = slmsig * error_scale
				else
					error(1,l+1,m+1) = clmsig
					error(2,l+1,m+1) = slmsig
				endif
			endif
			
			if (l > temp_max) temp_max = l
			
		elseif(c=="GRDOTA") then
			read(13,*) l, m, clm, slm, clmsig, slmsig, epoch_temp
			if (present(dot)) then
				if(l > maxl_dot) then
					print*, "Error ---SHRead2"
					print*, "DOT must be dimensioned (2, LMAX+1, LMAX+1) where LMAX is ", l
					print*, "Input array is dimensioned ", size(dot(1,:,1)),  size(dot(1,1,:)) 
					stop
				endif
			
			dot(1,l+1, m+1) = clm
			dot(2,l+1, m+1) = slm
			
			endif
			
			if (l > temp_max) temp_max = l
			
		elseif(c=="CALSDV") then
		
			read(13,*) l, m, clm, slm
			if (l > maxl_cilm) then
				print*, "Error ---SHRead2"
				print*, "CILM must be dimensioned (2, LMAX+1, LMAX+1) where LMAX is ", l
				print*, "Input array is dimensioned ", size(cilm(1,:,1)),  size(cilm(1,1,:)) 
				stop
			endif
				cilm(1,l+1,m+1) = clm
				cilm(2,l+1,m+1) = slm
			
			if (l > temp_max) temp_max = l
			
		elseif(c=="gfc") then
		
			read(13,*) l, m, clm, slm, clmsig, slmsig
			
			if (l > maxl_cilm) then
				print*, "Error ---SHRead2"
				print*, "CILM must be dimensioned (2, LMAX+1, LMAX+1) where LMAX is ", l
				print*, "Input array is dimensioned ", size(cilm(1,:,1)),  size(cilm(1,1,:)) 
				stop
			endif
			
			cilm(1,l+1,m+1) = clm
			cilm(2,l+1,m+1) = slm
			
			if (present(error)) then
				if(l > maxl_error) then
					print*, "Error ---SHRead2"
					print*, "ERROR must be dimensioned (2, LMAX+1, LMAX+1) where LMAX is ", l
					print*, "Input array is dimensioned ", size(error(1,:,1)),  size(error(1,1,:)) 
					stop
				endif
				
				if (error_scale /= 1.0d0) then
					error(1,l+1,m+1) = clmsig * error_scale
					error(2,l+1,m+1) = slmsig * error_scale
				else
					error(1,l+1,m+1) = clmsig
					error(2,l+1,m+1) = slmsig
				endif
			endif
			
			if (l > temp_max) temp_max = l
			
		elseif(c=="gfct") then
		
			read(13,*) l, m, clm, slm, clmsig, slmsig, epoch_temp
			
			if (l > maxl_cilm) then
				print*, "Error ---SHRead2"
				print*, "CILM must be dimensioned (2, LMAX+1, LMAX+1) where LMAX is ", l
				print*, "Input array is dimensioned ", size(cilm(1,:,1)),  size(cilm(1,1,:)) 
				stop
			endif
			
			cilm(1,l+1,m+1) = clm
			cilm(2,l+1,m+1) = slm
			
			if (present(error)) then
				if(l > maxl_error) then
					print*, "Error ---SHRead2"
					print*, "ERROR must be dimensioned (2, LMAX+1, LMAX+1) where LMAX is ", l
					print*, "Input array is dimensioned ", size(error(1,:,1)),  size(error(1,1,:)) 
					stop
				endif
				
				if (error_scale /= 1.0d0) then
					error(1,l+1,m+1) = clmsig * error_scale
					error(2,l+1,m+1) = slmsig * error_scale
				else
					error(1,l+1,m+1) = clmsig
					error(2,l+1,m+1) = slmsig
				endif
			endif
			
			if (l > temp_max) temp_max = l
		
		elseif(c=="dot") then
			read(13,*) l, m, clm, slm, clmsig, slmsig
			if (present(dot)) then
				if(l > maxl_dot) then
					print*, "Error ---SHRead2"
					print*, "DOT must be dimensioned (2, LMAX+1, LMAX+1) where LMAX is ", l
					print*, "Input array is dimensioned ", size(dot(1,:,1)),  size(dot(1,1,:)) 
					stop
				endif
			
			dot(1,l+1, m+1) = clm
			dot(2,l+1, m+1) = slm
			
			endif
			
			if (l > temp_max) temp_max = l

		elseif(c=="CMMNT") then
		
			read(13,"(A114)") comment
			print*, "SHRead2 --- Ignoring comment : ", comment
			
		else
			read(13,*)
			print*, "SHRead2 --- Ignoring line ", c
		endif
	enddo
	
	if (lmax == 0) lmax = temp_max
	
	if (present(doystart)) then
		doystart = doy1
	endif
	if (present(doyend)) then
		doyend = doy2
	endif
	if (present(epoch)) then
		epoch = epoch_temp
	endif

	close(13)
	
end subroutine SHRead2
