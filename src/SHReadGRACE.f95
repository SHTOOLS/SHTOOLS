subroutine SHReadGRACE(filename, cilm, lmax, gm, r0_pot, error)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	This program will read the file format of spherical harmonic
!	coefficients used by the GRACE group into standard arrays for
!	of Cilm and the error.
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
!			error:	An array containing the error coefficients.
!
!	Written by Mark Wieczorek (September 2005)
!
!	Copyright (c) 2005, Mark A. Wieczorek
!	All rights reserved.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	implicit none
	
	character(*), intent(in) ::	filename
	integer, intent(out) ::	lmax
	real*8, intent(out) ::	cilm(:,:,:), gm, r0_pot
	real*8, intent(out), optional :: error(:,:,:)
		
	integer:: l, m, stat
	character :: c*20, junk*20
	real*8 ::	clm, slm, clmsig, slmsig
	character*40 ::	formatstring
	
	cilm=0.0d0
	
	if (size(cilm(:,1,1)) < 2) then
		print*, "Error --- SHReadGRACE"
		print*, "CILM must be dimensioned (2, *, *)."
		print*, "Input array is dimensioned ", size(cilm(:,1,1)), size(cilm(1,:,1)), size(cilm(1,1,:))
		stop
	endif
	
	if (present(error)) then
		if (size(error(:,1,1)) < 2) then
			print*, "Error --- SHReadGRACE"
			print*, "ERROR must be dimensioned (2, *, *)."
			print*, "Input array is dimensioned ", size(error(:,1,1)), size(error(1,:,1)), size(error(1,1,:))
			stop
		endif
		
		error=0.0d0
	endif

	
	open(13,file=filename)
	
	read(13,"(a40)", iostat=stat) formatstring
	if(stat /=0) then
		print*, "Error --- SHReadGRACE"
		print*, "Problem reading first FORMATSTRING from file ", filename
		stop
	endif
	
	read(13,formatstring, iostat=stat) c, junk, gm, r0_pot
	if(stat /=0) then
		print*, "Error --- SHReadGRACE"
		print*, "Problem reading second line of file ", filename
		stop
	endif
	
	read(13,"(a40)", iostat=stat) formatstring
	if(stat /=0) then
		print*, "Error --- SHReadGRACE"
		print*, "Problem reading second FORMATSTRING from file ", filename
		stop
	endif
	
	lmax = -1
	
	do 
		read(13,formatstring, iostat=stat) c, l, m, clm, slm, clmsig, slmsig
		if(stat /= 0) exit
		
		if (l > lmax) lmax = l

		if (size(cilm(1,:,1)) < lmax+1 .or. size(cilm(1,1,:)) < lmax+1) then
			print*, "Error ---SHReadGRACE"
			print*, "CILM must be dimensioned (2, LMAX+1, LMAX+1) where LMAX is ", lmax
			print*, "Input array is dimensioned ", size(cilm(1,:,1)),  size(cilm(1,1,:)) 
			stop
		elseif (present(error)) then
			if (size(error(1,:,1)) < lmax+1 .or. size(error(1,1,:)) < lmax+1) then
				print*, "Error ---SHReadGRACE"
				print*, "ERROR must be dimensioned (2, LMAX+1, LMAX+1) where LMAX is ", lmax
				print*, "Input array is dimensioned ", size(error(1,:,1)),  size(error(1,1,:)) 
				stop
			endif
		endif

		if (c/="RECOEF") cycle
		
		if (present(error)) then
			cilm(1,l+1,m+1) = clm
			cilm(2,l+1,m+1) = slm
			error(1,l+1,m+1) = clmsig
			error(2,l+1,m+1) = slmsig
		else
			cilm(1,l+1,m+1) = clm
			cilm(2,l+1,m+1) = slm
		endif

	enddo

	close(13)
	
end subroutine SHReadGRACE
