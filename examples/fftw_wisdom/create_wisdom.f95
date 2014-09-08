program create_wisdom
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	This program will demonstrate how to create an fftw file for a specific sized 
!	transform using the unix command fftw-wisdom. For a given spherical harmonic 
!	degree l, there will be 2l+1 data points used in the longitudinal FFTs employed 
!	in the Gauss-Legendre integrations, and 2l+2 for the Driscoll and Healy intergrations.
!	When using an equally SPACED DH grid there are 4l+4 points in longitude.
!
!	In order to compute a wisdom file, use the following unix command
!		
!		fftw-wisdom -v -x -o wisdom_file [SIZES]
!
!	The -x command performs an exhaustive search.
!	
!	SHTOOLS will try to read the wisdom file in the default location
!		/etc/fftw/wisdom
!	Otherwise, one can specify its location using the optional argument wisdom_file 
!	with the routine PreCompute_FFTW.
!	
!	Notes:
!		1.	Only real data FFTs are used in SHTOOLS
!		2.	The FFTs are "out-of-place" in that the length of the 
!			input and output arrays are different.
!		3.	Both forward and backward FFTs are used.
!		4.	I have found that my code runs at exactly the same speed 
!			when using the above created wisdom file on a Mac G4 laptop. Results may
!			be different on differnt machine architectures.
!
!	Written by Mark Wieczorek June 2004
!
!	Copyright (c) 2005-2006, Mark A. Wieczorek
!	All rights reserved.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	implicit none
	integer ::	lmax, n, i
	character*9 ::	str(6)
	
	print*, "Spherical harmonic degree to compute transforms > "
	read(*,*) lmax

	if (lmax >= 2800) then
		print*, "Warning --- SHTOOLS is not accurate beyond ~2800."
	endif
	
	do i=1, 3
	
		if (i==1) then
			n = 2*lmax +1
		elseif (i==2) then
			n = 2*(lmax+1)
		elseif (i==3) then
			n = 4*(lmax+1)
		endif
	
		if (n<10) then
			write(str(i) ,fmt="('rof', i1)") n
			write(str(3+i), fmt="('rob', i1)") n
		elseif (n<100) then
			write(str(i),fmt="('rof', i2)") n
			write(str(3+i),fmt="('rob', i2)") n
		elseif(n<1000 .and. n>=100) then
			write(str(i),fmt="('rof', i3)") n
			write(str(3+i),fmt="('rob', i3)") n
		elseif(n<10000 .and. n>=1000) then
			write(str(i),fmt="('rof', i4)") n
			write(str(3+i),fmt="('rob', i4)") n
		elseif(n<100000 .and. n>=10000) then
			write(str(i),fmt="('rof', i5)") n
			write(str(3+i),fmt="('rob', i5)") n
		endif
		
	enddo	

	print*, "Create wisdom file by using the following command:"
	write(*,*)
	print*,  "fftw-wisdom -v -x -o wisdom ", trim(str(1)), " ", trim(str(2)), " ", trim(str(3)), " ", trim(str(4)), &
		" ", trim(str(5)), " ", trim(str(6))
	write(*,*)
	print*, "If a system wisdom file does not already exist, or if you do not want "
	print*, "to read the system wisdom file, add the option -n."
	print*, "To read a specific wisdom file, add the option -w wisdom_file_name."
	print*, "The default location for the system wisdom file is"
	print*, "/etc/fftw/wisdom (copying to this directory will require root privileges)."

end program create_wisdom 

