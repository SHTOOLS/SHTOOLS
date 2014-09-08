subroutine Curve2Mask(dhgrid, n, sampling, profile, nprofile, NP)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	Given a list of coordinates of a SINGLE CLOSED CURVE, this routine 
!	will fill the interior and exterior with 0s and 1s. The value at the 
!	north pole (either 0 or 1) is specified by the input parameter NP.
!
!	Calling Parameters
!		IN
!			profile		Latitude (:,1) and longitude (:,2) coordinates of a 
!					single close curve having dimension (nprofile, 2).
!			nprofile 	Number of coordinates in the vector profile.
!			np		Value of the output function at the North Pole, which can
!					be either 0 or 1.
!			n		Number of latitude bands in the output Driscoll and Healy sampled
!					grid.
!			sampling	1 sets the number of longitude bands equal to 1, whereas 2
!					sets the number to twice that of n.
!		OUT
!			dhgrid		A Driscoll and Healy sampled grid specifiying whether the point
!					is in the curve (1), or outside of it (0).
!
!
!	Written by Mark Wieczorek (August 2009)
!
!	Copyright (c) 2009, Mark A. Wieczorek
!	All rights reserved.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	implicit none
	integer, intent(out) ::		dhgrid(:,:)
	real*8, intent(in) ::		profile(:,:)
	integer, intent(in) ::		n, sampling, nprofile, np
	integer, parameter ::		maxcross = 1000
	integer ::			i, j, k, k_loc, nlat, nlong, numcross, next, ind1, ind2
	real*8 ::			lat_int, long_int, lon, cross(maxcross), cross_sort(maxcross)
	
	nlat = n
	lat_int = 180.0d0 / dble(nlat)
	dhgrid = 0
	
	if (sampling == 1) then 
		nlong = nlat
		long_int = 2.0d0 * lat_int
	elseif(sampling == 2) then
		nlong = 2 * nlat
		long_int = lat_int
	else
		print*, "Error --- Curve2Mask"
		print*, "SAMPLING of DHGRID must be 1 (equally sampled) or 2 (equally spaced)."
		print*, "SAMPLING = ", sampling
		stop
	endif
	
	if (NP /=1 .and. NP /= 0) then
		print*, "Error --- Curve2Mask"
		print*, "NP must be 0 if the North pole is outside of curve,"
		print*, "or 1 if the North pole is inside of the curve."
		print*, "NP = ", np
		stop
	endif
	
	if (size(dhgrid(1,:)) < nlong .or. size(dhgrid(:,1)) < nlat ) then
		print*, "Error --- Curve2Mask"
		print*, "DHGRID must be dimensioned as (NLAT, NLONG)."
		print*, "NLAT = ", nlat
		print*, "NLONG = ", nlong
		print*, "Size of GRID = ", size(dhgrid(:,1)), size(dhgrid(1,:))
	endif

	if (size(profile(:,1)) < nprofile .or. size(profile(1,:)) < 2) then
		print*, "Error --- Curve2Mask"
		print*, "PROFILE must be dimensioned as (NPROFILE, 2)."
		print*, "Dimension of NPROFILE = ", size(profile(:,1)), size(profile(1,:))
		stop
	endif
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!
	!	Start at 90N and 0E. Determine where the curve crosses this longitude band,
	!	sort the values, and then set the pixels between zero crossings to either 0
	!	or 1.
	!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	do j=1, nlong
		lon = dble(j-1)*long_int
	
		numcross = 0
		do i=1, nprofile-1
			if (profile(i,2) <= lon .and. profile(i+1,2) > lon) then
				numcross = numcross + 1
				if (numcross > maxcross) then
					print*, "Error --- Curve2Mask"
					print*, "Internal variable MAXCROSS needs to be increased."
					print*, "MAXCROSS = ", maxcross
					stop
				endif
				cross(numcross) = profile(i,1) + (profile(i+1,1)-profile(i,1)) / &
						(profile(i+1,2)-profile(i,2)) * (lon - profile(i,2))
			elseif (profile(i,2) >= lon .and. profile(i+1,2) < lon) then
				numcross = numcross + 1
				if (numcross > maxcross) then
					print*, "Error --- Curve2Mask"
					print*, "Internal variable MAXCROSS needs to be increased."
					print*, "MAXCROSS = ", maxcross
					stop
				endif
				cross(numcross) = profile(i+1,1) + (profile(i,1)-profile(i+1,1)) / &
						(profile(i,2)-profile(i+1,2)) * (lon - profile(i+1,2))
			endif
		enddo
		
		! do first and last points
		
		if (profile(nprofile,2) <= lon .and. profile(1,2) > lon) then
			numcross = numcross + 1
			if (numcross > maxcross) then
				print*, "Error --- Curve2Mask"
				print*, "Internal variable MAXCROSS needs to be increased."
				print*, "MAXCROSS = ", maxcross
				stop
			endif
			cross(numcross) = profile(nprofile,1) + (profile(1,1)-profile(nprofile,1)) / &
					(profile(1,2)-profile(nprofile,2)) * (lon - profile(nprofile,2))
		elseif (profile(nprofile,2) >= lon .and. profile(1,2) < lon) then
			numcross = numcross + 1
			if (numcross > maxcross) then
				print*, "Error --- Curve2Mask"
				print*, "Internal variable MAXCROSS needs to be increased."
				print*, "MAXCROSS = ", maxcross
				stop
			endif
			cross(numcross) = profile(1,1) + (profile(nprofile,1)-profile(1,1)) / &
					(profile(nprofile,2)-profile(1,2)) * (lon - profile(1,2))
		endif


		if (numcross > 0) then	! sort crossings by decreasing latitude
			do k=1, numcross
				k_loc = maxloc(cross(1:numcross), 1)
				cross_sort(k) = cross(k_loc)
				cross(k_loc) = -999.0d0
			enddo
		endif
			
			
		if (numcross == 0) then
			dhgrid(1:nlat, j) = np
			
		elseif (numcross == 1) then
			ind1 = int( (90.0d0 - cross_sort(1)) / lat_int) + 1
			dhgrid(1:ind1, j) = np
			if (ind1 == nlat) then
				cycle
			elseif (np == 0) then
				dhgrid(ind1+1:nlat, j) = 1
			else
				dhgrid(ind1+1:nlat, j) = 0
			endif
			
		else
			ind1 = 1
			next = np
			do k=1, numcross
				ind2 = int( (90.0d0 - cross_sort(k)) / lat_int) + 1
				if (ind2 >= ind1) dhgrid(ind1:ind2, j) = next
				if (next == 0) then
					next = 1
				else 
					next = 0
				endif
				ind1 = ind2 + 1
			enddo
			
			if (ind1 <= nlat) dhgrid(ind1:nlat, j) = next
			
		endif
		
	enddo

end subroutine Curve2Mask
