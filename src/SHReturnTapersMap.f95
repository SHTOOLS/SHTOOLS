subroutine SHReturnTapersMap(tapers, eigenvalues, dh_mask, n_dh, sampling, lmax, Ntapers)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	This subroutine will calculate the eigenvalues and eigenfunctions of the generalized
!	concentration problem where the region of interest R is given by the grid DH_MASK. This
!	matrix is sampled according to the Driscoll and Healy sampling theorem, and possesses a value 
!	of 1 inside of R, and 0 elsewhere. Returned tapers and eigenvalues are ordered from the largest
!	to smallest eigenvalue, and the spectral coefficients are packed into a 1D column vector
!	according to the scheme described in YilmIndex.
!
!	The elements Dij are calculated approximately using spherical harmonic transforms.
!	The effective bandwidth of the grid DH_MASK should in general be larger than LMAX by about a factor 
!	of 4 or so for accurate results. In addition, as a result of the approximate manner in which
!	the space-concentration kernel is calculated, it is preferable to use SAMPLING=1.
!	As one test of the accuracy, the area of the input grid DH_MASK, is compared to the calculated 
!	area given in the elements D(1,1).
!
!	Calling Parameters
!		IN
!			dh_mask		Integer grid sampled according to the Driscoll and Healy sampling
!					theorem. A value of 1 indicates the the grid node is in the concentration
!					domain, and a value of 0 indicates that it is outside. Dimensioned as
!					(n_dh, n_dh) for SAMPLING = 1 or (n_dh, 2*n_dh) for SAMPLING = 2.
!			n_dh		The number of latitude samples in the Driscoll and Healy sampled grid.
!			SAMPLING	1 corresponds to equal sampling (n_dh, n_dh), whereas 2 corresponds
!					to equal spaced grids (n_dh, 2*n_dh).
!			lmax		Maximum spherical harmonic degree of the outpt spherical harmonic coefficients.
!		IN, OPTIONAL
!			ntapers		Number of tapers and eigenvalues to output
!		OUT
!			Tapers		Column vectors contain the spherical harmonic coefficients, packed according
!					to the scheme described in YilmIndex. The dimension of this array is 
!					(lmax+1)**2 by (lmax+1)**2, or (lmax+1)**2 by NTapers if NTapers is present.
!			Eigenvalues	A 1-dimensional vector containing the eigenvalues corresponding to the columns
!					of Tapers, dimensioned as (lmax+1)**2.
!
!	Written by Mark Wieczorek (August 2009)
!
!	Copyright (c) 2009, Mark A. Wieczorek
!	All rights reserved.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	use SHTOOLS, only : EigValVecSym, ComputeDMap
	
	implicit none
		
	real*8, intent(out) ::	tapers(:,:), eigenvalues(:)
	integer, intent(in) ::	dh_mask(:,:), n_dh, sampling, lmax
	integer, intent(in), optional ::	ntapers
	
	real*8, allocatable ::	dij(:,:)
	integer ::		nlat, nlong, lmax_dh, astat, i, j
	real*8 ::		area, areaf, pi, colat, long_int, lat_int, da
	
	pi = acos(-1.0d0)
	
	if (present(ntapers)) then
		if (ntapers > (lmax+1)**2) then
			print*, "Error --- SHRetrunTapersMap"
			print*, "The number of output tapers must be less than or equal to (LMAX+1)**2."
			print*, "LMAX = ", lmax
			print*, "NTAPERS = ", ntapers
			stop
		elseif (size(tapers(:,1)) < (lmax+1)**2 .or. size(tapers(1,:)) < ntapers) then
				print*, "Error --- SHReturnTapersMap"
				print*, "TAPERS must be dimensioned as ((LMAX+1)**2, NTAPERS)."
				print*, "Dimension of TAPERS = ", size(tapers(:,1)), size(tapers(1,:))
				print*, "LMAX = ", lmax
				print*, "NTAPERS = ", ntapers
				stop
		elseif (size(eigenvalues(:)) < ntapers) then
			print*, "Error --- SHReturnTapersMap"
			print*, "EIGENVALUES must be dimensioned as NTAPERS."
			print*, "Dimension of EIGENVALUES = ", size(eigenvalues)
			print*, "NTAPERS = ", ntapers
			stop
		endif
	else
		if (size(tapers(:,1)) < (lmax+1)**2 .or. size(tapers(1,:)) < (lmax+1)**2) then
			print*, "Error --- SHReturnTapersMap"
			print*, "TAPERS must be dimensioned as ((LMAX+1)**2, (LMAX+1)**2)."
			print*, "Dimension of TAPERS = ", size(tapers(:,1)), size(tapers(1,:))
			print*, "LMAX = ", lmax
			stop
		elseif (size(eigenvalues(:)) < (lmax+1)**2) then
			print*, "Error --- SHReturnTapersMap"
			print*, "EIGENVALUES must be dimensioned as (LMAX+1)**2"
			print*, "Dimension of EIGENVALUES = ", size(eigenvalues)
			print*, "LMAX = ", lmax
			stop
		endif
	endif
	
	if (mod(n_dh,2) /= 0) then
		print*, "Error --- SHReturnTapersMap"
		print*, "Number of samples in latitude must be even for the Driscoll and Healy sampling theorem."
		print*, "N_DH = ", n_dh
		stop
	endif
	
	nlat = n_dh
	lat_int = pi/dble(nlat)
	if (sampling == 1) then
		nlong = nlat
		long_int = 2.0d0 * lat_int
	elseif (sampling == 2) then
		nlong = 2*nlat
		long_int = lat_int
	else 
		print*, "Error --- SHReturnTapersMap"
		print*, "SAMPLING must be either 1 (equally sampled) or 2 (equally spaced)."
		print*, "SAMPLING = ", sampling
		stop
	endif
	
	if (size(dh_mask(:,1)) < nlat .or. size(dh_mask(1,:)) < nlong) then
		print*, "Error --- SHReturnTapersMap"
		print*, "DH_MASK must be dimensioned as ", nlat, nlong
		print*, "Dimensions of DH_MASK are ", size(dh_mask(:,1)), size(dh_mask(1,:))
		stop
	endif
	
	lmax_dh = n_dh/2 - 1

	if (lmax_dh < lmax) then
		print*, "Error --- SHReturnTapersMap"
		print*, "The effective bandwith of the input grid DH_MASK must be greater or equal than LMAX."
		print*, "Experience suggests that this should be about 4 times LMAX."
		print*, "LMAX = ", lmax
		print*, "Effective bandwidth of DH_MASK = (N/2 -1) = ", lmax_dh
		stop
	endif
	
	allocate(dij((lmax+1)**2, (lmax+1)**2), stat = astat)
	if (astat /= 0) then
		print*, "Error --- SHReturnTapersMap"
		print*, "Problem allocating DIJ((LMAX+1)**2,(LMAX+1)**2)"
		print*, "LMAX = ", lmax
		stop
	endif
	
	call ComputeDMap(Dij, dh_mask, n_dh, sampling, lmax)
	
	if (present(ntapers)) then
		call EigValVecSym(Dij, (lmax+1)**2, eigenvalues, tapers, k = ntapers)
	else
		call EigValVecSym(Dij, (lmax+1)**2, eigenvalues, tapers)
	endif
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!
	! 	Check and see how good the solutions are.
	!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	area = 0.0d0
	areaf = 0.0d0
	
	do i=2, nlat
		colat = dble(i-1)*lat_int
		do j=1, nlong
			da = -( cos(colat + lat_int/2.0d0) - cos(colat - lat_int/2.0d0) ) * long_int
			area = area + da
			if (dh_mask(i,j) == 1) areaf = areaf + da
		enddo
	enddo

	da = 2.0d0 * pi * (1.0d0 - cos(lat_int/2.0d0) )
	area = area + 2.0d0 * da
	if (dh_mask(1,1) == 1) areaf = areaf +  da
	if (dh_mask(nlat,nlong) == 1) areaf = areaf + da
	
	print*, "SHReturnTapersMap --- (Area of entire grid) / 4pi = ", area / 4.0d0 / pi
	print*, "SHReturnTapersMap --- (Area of DH_MASK) / 4pi = ", areaf / 4.0d0 / pi
	print*, "SHReturnTapersMap --- (Area of DH_MASK) / 4pi / D00 = ", areaf / 4.0d0 / pi / Dij(1,1)

	deallocate(dij)
		
end subroutine SHReturnTapersMap
