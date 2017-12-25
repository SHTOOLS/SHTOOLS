program TestSHRotate
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	This program will read in a spherical harmonic file, then rotate 
!	the body by the three euler angles.
!
!	Written by Mark Wieczorek (August 2003)
!
!	Copyright (c) 2005, Mark A. Wieczorek
!	All rights reserved.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	use SHTOOLS

	implicit none
	real*8 ::		x(3), pi, interval, alpha, beta, gamma
	real*8, allocatable ::	cilm1(:,:,:), cilm2(:,:,:), dj(:,:,:), grid(:,:)
	character*120 ::	outsh, outgrid, infile
	integer ::		lmax, nlat, nlong, i, j, l, m, bodycoord, astat(4)
	
	pi = acos(-1.0d0)
	
	lmax = 719
	allocate(cilm1(2, lmax+1, lmax+1), stat = astat(1))
	if (astat(1) /= 0) then
		print*, "Problem allocating arrays"
		stop
	endif
	
	infile = "../../ExampleDataFiles/MarsTopo719.shape"
	print*, "Input file = ", infile
	call SHRead(infile, cilm1, lmax)
	print*, "Lmax of input file = ", lmax
	
	print*, "Maximum degree to be used with rotations > "
	read(*,*) lmax
	
	print*, "Input Euler Angles"
	print*, "Alpha (deg) = "
	read(*,*) alpha
	print*, "Beta (deg) = "
	read(*,*) beta
	print*, "Gamma (deg) = "
	read(*,*) gamma
	
	print*, "Do these correspond to "
	print*, "(1) Rotation of the coordinate system without rotation of the physical body, or"
	print*, "(2) Rotation of the physical body without rotation of the coordinate system."
	read(*,*) bodycoord
	
	if (bodycoord ==1) then
		x(1) = alpha ; x(2) = beta ; x(3) = gamma
	elseif (bodycoord ==2) then
		x(1) = -gamma ; x(2) = -beta ; x(3) = -alpha
	else
		print*, "Must enter 1 or 2"
		stop
	endif
	
	x = x * pi/180.0d0

	print*, "Output file name of rotated spherical harmonics > "
	read(*,*) outsh
	
	print*, "Grid spacing for output grid (deg) > "
	read(*,*) interval
	print*, "File name of gridded rotated field > "
	read(*,*) outgrid
	
	allocate(cilm2(2, lmax+1, lmax+1), stat = astat(2))
	allocate(dj(lmax+1, lmax+1, lmax+1), stat = astat(3))
	nlat = 180.0d0/interval + 1
	nlong = 360.0d0/interval + 1
	allocate(grid(nlat,nlong), stat = astat(4))
	
	if (sum(astat(2:4)) /= 0) then
		print*, "Problem allocating arrays"
		stop
	endif
		
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!
	!	Rotate coefficients to a specific coordinate
	!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	print*, "Computing rotation matrix"
	call djpi2(dj, lmax)	! Create rotation matrix used in the rotation routine.
	
	print*, "Rotating spherical harmonics"
	call SHRotateRealCoef(cilm2, cilm1, lmax, x, dj)	
	
	print*, "Making 2D grid"
	call MakeGrid2d(grid, cilm2, lmax, interval, nlat, nlong)		
				
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!
	!	Output data to disk
	!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	open(12,file=outsh)
	do l=0, lmax
		do m=0, l
			write(12,*) l, m, cilm2(1,l+1,m+1), cilm2(2,l+1,m+1)
		enddo
	enddo
	close(12)
	
	open(13, file=outgrid)
	write(13,*) nlat, nlong
	
	do i=1, nlat
		do j=1, nlong
			write(13,*) grid(i,j)
		enddo
	enddo		
	close(13)
	
	deallocate(cilm1)
	deallocate(cilm2)
	deallocate(dj)
	deallocate(grid)
		
end program TestSHRotate
