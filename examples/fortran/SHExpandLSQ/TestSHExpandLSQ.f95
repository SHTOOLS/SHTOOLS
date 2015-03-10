program TestSHExpandLSQ
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	This program will (1) create a set of unevenly spaced
!	data points from a spherical harmonic file, and then
!	(2) expand these into spherical harmonics using a least-squares 
!	inversion if there are more data points than spherical harmonic
!	coefficients (overdetermined case), or using a minimum norm
!	solution if there are more coefficients than data points
!	(underdetermined case).
!
!	Written by Mark Wieczorek (June 2004)
!
!	Copyright (c) 2005, Mark A. Wieczorek
!	All rights reserved.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	use SHTOOLS
	
	implicit none
	integer, parameter ::	dmax = 10000, degmax = 120
	real*8 ::		d(dmax), lat(dmax), lon(dmax), cilm(2,degmax+1, degmax+1), &
				x, y, z, pi, maxerror, cilm1(2,degmax+1, degmax+1), &
				dd(dmax), misfit
	integer ::		nmax, lmax, l, i, seed, lmaxfile
	character*80 ::		infile 

	
	d = 0.0d0
	dd = 0.0d0
	lon = 0.0d0
	lat = 0.0d0
	pi = acos(-1.0d0)
	
	infile = "../../ExampleDataFiles/MarsTopo719.shape"
	
	call SHRead(infile, cilm, lmaxfile)
	print*, "Maximum degree of spherical harmonic file = ", lmaxfile
	
	print*, "Number of random data points to use > "
	read(*,*) nmax
	lmax = floor(sqrt(dble(nmax)) - 1)
	print*, "Maximum spherical harmonic degree for overdetermined least-squares inversion = ", lmax
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!
	!	Create synthetic data from the known spherical harmonic
	!	coeffients. Data points are located at random points
	!	on the sphere.
	!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	seed = -13453
	
	do i = 1, nmax
		x = 2.0d0*RandomN(seed)-1.0d0
		y = 2.0d0*RandomN(seed)-1.0d0
		z = 2.0d0*RandomN(seed)-1.0d0
		lat(i) = atan2(z, sqrt(x**2 + y**2)) * 180.0d0/pi
		lon(i) = atan2(y, x) * 180.0d0/pi
		d(i) = MakeGridPoint(cilm, lmaxfile, lat(i), lon(i), norm=1)
	enddo
	
	print*, "maximum (km) = ", maxval(d(1:nmax))/1.d3
	print*, "minimum (km) = ", minval(d(1:nmax))/1.d3
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!
	!	Do least squares inversion for increasing values of l.
	!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	do l=1, lmax, 1
	
		print*, "l = ", l
	
		call SHExpandLSQ(cilm1, d(1:nmax), lat(1:nmax), lon(1:nmax), nmax, l, norm=1, chi2=misfit)
	
		do i=1, nmax
			dd(i) = MakeGridPoint(cilm1, l, lat(i), lon(i), norm=1)
		enddo
		
		maxerror = maxval(abs(d(1:nmax) - dd(1:nmax)))
	
		print*, "Maximum error between input and output data points (km) = ", maxerror/1.d3
		print*, "Sum of squares residuals (km^2) = ", misfit/1.d6
		
	enddo
	
end program TestSHExpandLSQ
	
	
	