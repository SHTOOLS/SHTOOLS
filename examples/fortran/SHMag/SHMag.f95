program SHMag
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	This program will read in a set of magnetic spherical harmonic coefficients
!	and make 2D equally spaced grid of the three vector components of the magnetic
!	field, the total magnetic field strength, and the magnetic potential. The field 
!	is calculated on a speroid with mean radius r and flattening f.
!
!	The included spherical harmonic file FSU_mars90.sh is the martian magnetic
!	field model of Cain et al., 2003.
!
!	Dependencies:	SHRead, MakeMagGridDH
!
!	Written by Mark Wieczorek, February 2004
!	
!	October 2012:	Modified to use the new routine MakeMagGridDH, instead of the old routine MakeMagGrid2D.
!
!	Copyright (c) 2012, Mark A. Wieczorek
!	All rights reserved.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	use SHTOOLS

	implicit none
	
	character*200 ::	infile, radf, thetaf, phif, totalf, potf	
	real*8 ::		header(4), interval, r0, a, f, mpr, z, timein, timeout, temp
	real*8, allocatable ::	glm(:,:,:), rad(:,:), phi(:,:), theta(:,:), total(:,:), pot(:,:)		
	integer ::		lmax, lmaxp, n, n_out, nlong, nlat, i, j, astat(6), sampling
	
	infile = "../../ExampleDataFiles/FSU_mars90.sh"
	f = 1.0d0/169.864881d0	! Mars flattening = (R_eq - R_p)/R_eq
	mpr = 3389.508d3	! Mean radius of mars
	z = 145.d3		! mean altitude to calculate field at.
	
	lmax = 359	! Maximum spherical harmonic degree that is resolved in the output grids.
	N = 2*lmax+2	! Number of latitudinal samples in the output grids.
	interval = 180.0d0 / dble(n)	! Sampling interval in latitude.
	
	sampling = 2
	nlat = n
	nlong = 2*n
	
	print*, "LMAX of output grids = ", lmax
	print*, "N = ", n
	print*, "Latitudinal sampling interval = ", interval
	
	allocate(glm(2,lmax+1,lmax+1), stat = astat(1))
	allocate(rad(nlat,nlong), stat = astat(2))
	allocate(theta(nlat,nlong), stat = astat(3))
	allocate(phi(nlat,nlong), stat = astat(4))
	allocate(total(nlat,nlong), stat = astat(5))
	allocate(pot(nlat,nlong), stat = astat(6))
	
	if (sum(astat(1:6)) /= 0) then
		print*, "Problem allocating GLM, RAD, THETA, PHI, TOTAL, and POT", astat(1), astat(2), astat(3), astat(4), astat(5), astat(6)
		stop
	endif
	
	call SHRead(infile, glm, lmaxp, header=header(1:4), skip=1)
	r0 = header(1)*1.d3
	print*, "Reference radius of spherical harmonic coefficients R0 (km) = ", r0/1.d3
	print*, "Lmax of spherical harmonic model = ", lmaxp
	
	a = mpr + z
	print*, "Field calculated on a flattened ellipsoid with semi-major axis A (km) = ", a/1.d3
	print*, "Flattening of the ellipsoid = ", f
		
	radf = "radial_145f.dat"
	thetaf = "theta_145f.dat"
	phif = "phi_145f.dat"
	totalf = "total_145f.dat"
	potf = "pot_145f.dat"
		
	call cpu_time(timein)
	
	call MakeMagGridDH(glm, lmax, r0, a, f, rad, theta, phi, total, n_out, sampling = sampling, lmax_calc = lmaxp, pot_grid = pot)
	
	call cpu_time(timeout)
	
	print*, "Elapsed time (sec) = ", timeout-timein
	
	if (n /= n_out) then 
		print*, "Problem. N is not equal to N_OUT."
		print*, "N_OUT = ", n_out
		stop
	endif
	
	print*, "Maximum and minimum field intensity (nT) = ", maxval(total(1:nlat,1:nlong)), minval(total(1:nlat,1:nlong))
	print*, "Maximum and minimum Br (nT) = ", maxval(rad(1:nlat,1:nlong)), minval(rad(1:nlat,1:nlong))
	print*, "Maximum and minimum Btheta (nT) = ", maxval(theta(1:nlat,1:nlong)), minval(theta(1:nlat,1:nlong))
	print*, "Maximum and minimum Bphi (nT) = ", maxval(phi(1:nlat,1:nlong)), minval(phi(1:nlat,1:nlong))
	print*, "Maximum and minimum potential (nT m) = ", maxval(pot(1:nlat,1:nlong)), minval(pot(1:nlat,1:nlong))

	
	open(12,file=radf)
	open(13,file=phif)
	open(14,file=thetaf)
	open(15,file=totalf)
	open(16,file=potf)
	
	write(12,*) nlat+1, nlong+1
	write(13,*) nlat+1, nlong+1
	write(14,*) nlat+1, nlong+1
	write(15,*) nlat+1, nlong+1
	write(16,*) nlat+1, nlong+1
	
	do i=1, nlat
		do j=1, nlong
			write(12,*) rad(i,j)
			write(13,*) phi(i,j)
			write(14,*) theta(i,j)
			write(15,*) total(i,j)
			write(16,*) pot(i,j)
		enddo
		! write out value at 360 degrees longitude, which is redundant
		write(12,*) rad(i,1)
		write(13,*) phi(i,1)
		write(14,*) theta(i,1)
		write(15,*) total(i,1)
		write(16,*) pot(i,1)
	enddo
	
	! write out values at 90 S, which are not calculated in the Driscoll and Healy routine, and which are just the average of the values 
	! of the last row.
	temp = sum(rad(nlat,1:nlong)) / dble(nlong)
	do j=1, nlong+1
		write(12,*) temp
	enddo
	temp = sum(phi(nlat,1:nlong)) / dble(nlong)
	do j=1, nlong+1
		write(13,*) temp
	enddo
	temp = sum(theta(nlat,1:nlong)) / dble(nlong)
	do j=1, nlong+1
		write(14,*) temp
	enddo
	temp = sum(total(nlat,1:nlong)) / dble(nlong)
	do j=1, nlong+1
		write(15,*) temp
	enddo
	temp = sum(pot(nlat,1:nlong)) / dble(nlong)
	do j=1, nlong+1
		write(16,*) temp
	enddo
	
	close(12)
	close(13)
	close(14)
	close(15)
	close(16)
	
	deallocate(rad)
	deallocate(theta)
	deallocate(phi)
	deallocate(total)
	deallocate(pot)
	deallocate(glm)
	
end program SHMag
