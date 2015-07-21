program TestCilmPlus
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	This program will 
!
!	(1)	Read in a spherical harmonic file of maximum degree Lmax, expand it in the 
!		spatial domain, and then re-expand it into spherical harmonics using
!		Guass-Legendre quadrature. The relativive errors in the space and spectral 
!		domaines will be determined. 
!	(2) 	Using these routines, the potential coefficients will be computed according 
!		to the algorithm of Wieczorek and Phillips (1998). The maximum order of
!		the exspansion will be determined as well.
!	(3)	The difference between using nmax*lmax coefficients when calculating the 
!		gravity coefficients vs. n*lmax will be determined.
!	(4)	A 2d gridded data file of the input coefficients with equal latitude and 
!		longitude spacings) will be output. 
!
!	Notes: The only two parameters in this program that might want
!	to be changed is MAX and MAX2DGRID. The first sets the maximum 
!	degree that can be computed in the spherical harmonic expansions,
!	while the second sets the maximum size of the output 2d grid.
!
!	Comments on compiling and modifying this code:
!	
!	!. This example code uses a lot of memory, and runs slowly in order to demonstrate
!	the accuracy of the routines. On OS X, it was necessary for me to compile it with
!	the -s option (static memory). Also, if memory is not an issue, then the arrays plx and plxn
!	should be precomputed.
!	2. As is seen in the output, for Mars (and also the Moon) nmax should be about 5
!	in order to achieve an accuracy comparable to the spherical harmonic model.
!	3. In order to exactly calculate the gravity field using the Wieczorek and Phillips (1998)
!	method, if the original field is defined up to lmax, the final output field should
!	be nmax*lmax. This is achieved in practice by zero padding the input field from 
!	lmax+1 to nmax*lmax+1. However, as is seen in the output, in practice, it is only 
!	necessary to increase the size of the input field by about 50%. At the same time, since the topography
!	spherical harmonic coefficients are NOT really zero beyond lmax, zero padding between 
!	(lmax+1, and nmax*lmax+1) still only yields an approximation to the gravity field. 
!	This would only be strictly correct if the topography was bandlimited to degree lmax.
!
!	Written by Mark Wieczorek (October, 2003)
!
!	Copyright (c) 2005, Mark A. Wieczorek
!	All rights reserved.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	use SHTOOLS
	use PlanetsConstants
	
	implicit none

	integer, parameter ::	degmax = 650, max2dgrid = 181
	real*8 :: 		zero(degmax+1), w(degmax+1), c1ilm(2, degmax+1, degmax+1), &
				zeron(degmax+1), wn(degmax+1), &
				cilm(2, degmax+1, degmax+1), c2ilm(2,degmax+1, degmax+1), grid1glq(degmax+1, 2*degmax+1), &
				gridglq(degmax+1, 2*degmax+1), grid2glq(degmax+1, 2*degmax+1), &
				grid3glq(degmax+1, 2*degmax+1), pi, maxerror, err1, err2, &
				d, rho, mass, grav, gm, grid2d(max2dgrid, 2*max2dgrid), interval, &
				timein, timeout
	real*8, allocatable ::		plxn(:,:),  plx(:,:)
	integer :: 		lmax, i, j, l, m, nlat, nlong, nmax, n, astat, &
				lmaxfile, lmaxn, nlatn, nlongn, precomputeplx
	character :: 		infile*80, outfile*80
	
	call cpu_time(timein)
	
	pi = acos(-1.0d0)
	
	rho = 2900.0d0			! assumed crustal density
	grav = Grav_Constant
	gm = GM_Mars
	mass = Mass_Mars	
	
	precomputeplx = 0			! set to 1 precompute arrays PLXN and PLX
	
	
	!print*, "Spherical Harmonic file > "
	!read(*,*) infile
	infile = "../../ExampleDataFiles/MarsTopo719.shape"
	
	!print*, "Maximum degree to read > "
	!read(*,*) lmax
	lmax = 90
	
	if (precomputeplx == 1) then
		allocate (plx(lmax+1, (lmax+1)*(lmax+2)/2), stat = astat)
		if (astat /=0) then
			print*, "Problem allocating plx"
			stop
		endif
	endif
	
	if (lmax > degmax) then
		print*, "lmax must be less than or equal to , ", degmax
		print*, "Change max and recompile"
		stop
	endif
	
	nlat = NGLQSH(lmax)	! These are the sizes of the grids used in  (lmax+1)
	nlong = 2*lmax + 1	! the gauss-legendre quadtrature routines
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!		
	!	Read Spherical harmonic coefficients
	!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	call SHRead(infile, cilm, lmaxfile)

	print*, "Input file name = ", infile
	print*, "Maximum degree of input file = ", lmaxfile
	print*, "Maximum degree used in the calculations = ", lmax
	if (lmaxfile > lmax) then
		cilm(:,lmax+2:lmaxfile+1,lmax+2:lmaxfile+1) = 0.0d0	! zero coefficients greater than lmax
	endif
		
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!
	! 	Call a routine which will compute a number of arrays that will be used later on 
	!	(these include the guass points, weights, and plx).
	!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	print*, "Precomputing data structures..."
	
	if (precomputeplx == 1) then
		call SHGLQ(lmax, zero, w, plx = plx)
	else	
		call SHGLQ(lmax, zero, w)
	endif

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!
	!	Compute a grid from the spherical harmonic coefficients.
	!	Note that this grid is equally spaced in longitude,
	!	but unevenly space in latitude according to the gauss points.
	!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		
	print*, "Making Grid..."
	
	if (precomputeplx == 1) then
		call MakeGridGLQ(gridglq, cilm, lmax, plx = plx)
	else
		call MakeGridGLQ(gridglq, cilm, lmax, zero = zero)
	endif
		
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!
	!	This routine expands the grid into spherical harmonics
	!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	print*, "Expanding into spherical harmonics"
	
	if (precomputeplx == 1) then
		call SHExpandGLQ(c2ilm, lmax, gridglq, w, plx = plx)
	else
		call SHExpandGLQ(c2ilm, lmax, gridglq, w, zero = zero)
	endif
		
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!
	!	Next, compute relative errors between the original and
	!	computed spherical harmonic coefficients
	!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	maxerror = 0.0d0
	do l=0, lmax
		do m=0, l
			if (m==0) then
				err1 = abs(cilm(1,l+1,m+1)-c2ilm(1,l+1,m+1))/abs(cilm(1,l+1,m+1))
				if (err1 >= maxerror) maxerror = err1
			else
				err1 = abs(cilm(1,l+1,m+1)-c2ilm(1,l+1,m+1))/abs(cilm(1,l+1,m+1))
				err2 = abs(cilm(2,l+1,m+1)-c2ilm(2,l+1,m+1))/abs(cilm(2,l+1,m+1))
				if (err1 >= maxerror) maxerror = err1
				if (err2 >= maxerror) maxerror = err2
			endif
		enddo
	enddo
	
	print*, "Maximum relative error between spherical harmonic coefficients = ", maxerror
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!
	!	Next, compute the maximum error in the space domain 
	!	at the gauss points.
	!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	if (precomputeplx == 1) then
		call MakeGridGLQ(grid2glq, cilm(:,1:lmax+1,1:lmax+1)-c2ilm(:,1:lmax+1,1:lmax+1), lmax, plx = plx)
	else
		call MakeGridGLQ(grid2glq, cilm(:,1:lmax+1,1:lmax+1)-c2ilm(:,1:lmax+1,1:lmax+1), lmax, zero = zero)
	endif
	
	print*, "Maximum error (meters) at the gauss points = ", maxval(abs(grid2glq(1:nlat,1:nlong)))
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!
	!	Compute the potential coefficients due to the topography
	!	using nmax, and then determine the maximum difference between 
	!	this gravity field and that calculated at a smaller n
	!	in the space domain. Note, even though the gravity field
	!	is calculated up to nmax*lmax, the differences that are
	!	output only correspond to the first lmax coefficients.
	!
	!	Note: the minimum gravity anomaly of Lemoine et al. (2001) is
	!	about 6 mgals. I have also determined that the the relative
	!	error for the Moho coefficients is smaller than that of the
	!	topography (The opposite is the case for the Moon). Thus the 
	!	value of n used for the topography is sufficient for the moho as well. 
	!
	!	The reason that the maximum degree of the spherical harmonic
	!	fields needs to be increased beyond lmax is that raising a field to the nth
	!	power results in a field with coefficients up to n*lmax.
	!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	nmax = 7
	lmaxn = lmax*nmax
	nlatn = NGLQSH(lmaxn)	! These are the sizes of the grids used in  (lmax+1)
	nlongn = 2*lmaxn + 1	! the gauss-legendre quadtrature routines
	
	print*, "Maximum spherical harmonic degree corresponding to n = ", nmax, " = ", lmaxn
	if (lmaxn > degmax ) then
		print*, "lmaxn is greater than max. Change max and recompile."
		stop
	endif
	
	if (precomputeplx == 1) then
		! This must be called again since the value of lmax is different than above
		allocate (plxn(lmaxn+1, (lmaxn+1)*(lmaxn+2)/2), stat = astat)
		if (astat /=0) then
			print*, "Problem allocating plxn"
			stop
		endif
	endif
	
	
	print*, "Computing gravity coefficients..."
	if (precomputeplx == 1) then
		call SHGLQ(lmaxn, zeron, wn, plx = plxn)
		call MakeGridGLQ(gridglq, cilm, lmaxn, plx = plxn)
		call CilmPlus(c1ilm, gridglq, lmaxn, nmax, mass, d, rho, gridtype = 1, w = wn, plx = plxn)
	else
		call SHGLQ(lmaxn, zeron, wn)
		call MakeGridGLQ(gridglq, cilm, lmaxn, zero = zeron)
		call CilmPlus(c1ilm, gridglq, lmaxn, nmax, mass, d, rho, gridtype = 1, w = wn, zero = zeron)
	endif
	
	do l=0, lmaxn
		c1ilm(:,l+1,1:l+1) = c1ilm(:,l+1,1:l+1)*(l+1.0d0)*gm*(1.0d5)/d**2	! convert to gravity coefficients in mgals
	enddo
	
	do n=1, nmax-1
	
		if (precomputeplx == 1) then
			call CilmPlus(c2ilm, gridglq, lmaxn, n, mass, d, rho, gridtype = 1, w = wn, plx = plxn)
		else
			call CilmPlus(c2ilm, gridglq, lmaxn, n, mass, d, rho, gridtype = 1, w = wn, zero = zeron)
		endif
		
		do l=0, lmaxn
			c2ilm(:,l+1,1:l+1) = c2ilm(:,l+1,1:l+1)*(l+1.0d0)*gm*(1.0d5)/d**2	
		enddo
		
		! Compute difference only up to lmax, not lmaxn
		
		if (precomputeplx == 1) then
			call MakeGridGLQ(grid3glq, c1ilm(:,1:lmax+1,1:lmax+1) - c2ilm(:,1:lmax+1,1:lmax+1), lmax, plx = plx)
		else
			call MakeGridGLQ(grid3glq, c1ilm(:,1:lmax+1,1:lmax+1) - c2ilm(:,1:lmax+1,1:lmax+1), lmax, zero = zero)
		endif
		
		print*, "Maximum difference (mgals) between nmax = ", nmax, " and ", n, " = ", &
			maxval(abs(grid3glq(1:nlat,1:nlong)))

	enddo
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!
	!	Next calculate the difference in the gravity field that arrises from 
	!	doing the GLQ integrations with a field zero-padded up to lmaxn as opposed
	!	to nmax*lmax.
	!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	do n=1, 10
		!lmaxn = lmax*n
		lmaxn = lmax + (n-1)*10
		nlatn = NGLQSH(lmaxn)	
		nlongn = 2*lmaxn + 1	
		
		if (precomputeplx == 1) then
			call SHGLQ(lmaxn, zeron, wn, plx = plxn)
			call MakeGridGLQ(gridglq, cilm, lmaxn, plx = plxn)
			call CilmPlus(c2ilm, gridglq, lmaxn, nmax, mass, d, rho, gridtype = 1, w = wn, plx = plxn)
		else
			call SHGLQ(lmaxn, zeron, wn)
			call MakeGridGLQ(gridglq, cilm, lmaxn, zero = zeron)
			call CilmPlus(c2ilm, gridglq, lmaxn, nmax, mass, d, rho, gridtype = 1, w = wn, zero = zeron)
		endif
		
		do l=0, lmaxn
			c2ilm(:,l+1,1:l+1) = c2ilm(:,l+1,1:l+1)*(l+1.0d0)*gm*(1.0d5)/d**2	! convert to gravity coefficients in mgals
		enddo
	
		! Compute difference only up to lmax, not lmaxn
		
		if (precomputeplx == 1) then
			call MakeGridGLQ(grid1glq, c1ilm(:,1:lmax+1,1:lmax+1)-c2ilm(:,1:lmax+1,1:lmax+1), lmax, plx = plx)
		else
			call MakeGridGLQ(grid1glq, c1ilm(:,1:lmax+1,1:lmax+1)-c2ilm(:,1:lmax+1,1:lmax+1), lmax, zero = zero)
		endif
	
		print*, "Maximum difference (mgals) between using Lmax of expansion equal to ", lmax*nmax, " and LMAX = ", lmaxn, &
			 maxval(abs(grid1glq(1:nlat,1:nlong)))
	
	enddo
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!
	!	Calculate a 2d grid with equal longitude and latitude spacings. 
	!	Note that if you want a grid finer than 1 degree, that you must
	!	increase the parameter maxgrid2d.
	!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	print*, "Making 2d grid with equal latitude and longitude spacings"
	interval = 2
	cilm(1,1,1) = 0.0d0		! set degree-zero term to zero
	cilm(1,3,1) = 0.0d0		! set j2 term to zero	

	outfile = "grid2d.dat"
	open(12,file = outfile)

	call MakeGrid2d(grid2d, cilm, lmax, interval, nlat, nlong)

	write(12,*) nlat, nlong 	! header line giving the dimensions of the following
					! raster image
	do i=1, nlat
		do j=1, nlong
			write(12,*) grid2d(i,j)
		enddo
	enddo			
	
	close(12)
	
	if (precomputeplx == 1) then
		deallocate(plx)
		deallocate(plxn)
	endif

	call cpu_time(timeout)
	print*, "time (sec) = ", timeout-timein

	
end program TestCilmPlus

