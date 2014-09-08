subroutine ComputeDMap(Dij, dh_mask, n_dh, sampling, lmax)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	This subroutine will compute the matrix D(lm,l'm') for the spatiospectral concentration
!	problem using the input grid DH_MASK to define the concentration region. The
!	input grid must be sampled according to the Driscoll and Healy sampling theorem
!	and must possess values of 1 inside the concentration region, and 0 elsewhere. The grid
!	possess N_DH samples in latitude, and either N_DH samples in longitude for SAMPLING=1
!	or 2*N_DH samples in longitude for SAMPLING=2. Elements of the matrix D
!	(which are ordered accoring to YilmIndex) are computed up to a maximum degree LMAX.
!
!	It should be noted that the elements of D are computed approximately. These are explicitly
!	calculated using
!
!		Dlm,l'm' = 1/(4 pi) Int_Omega F Yl'm' dOmega
!
!	where
!
!		F = Ylm * DH_MASK. 
!
!	This integral would be computed exactly if F had a bandwith of
!	L = N_DH/2 - 1. However, since DH_MASK has an infinite bandwidth, this will not be true. However,
!	if L is larger than LMAX, then we might expect that the spherical harmonic expansions will
!	be good approximations. This needs to be verified by chosing different values of N_DH and
!	comparing the results. A simple comparison is to see how good the D(1,1) element approximates
!	the area of the concentration region.
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
!			lmax		Maximum spherical harmonic degree of elements in the matrix D.
!		OUT
!			Dij		Elements of the matrix D, which is symmetric, and whose elements
!					are packed into a 1D array according to YilmINDEX. Dimesions of 
!					Dij are ((lmax+1)**2, (lmax+1)**2).
!
!	Written by Mark Wieczorek (August 2009)
!
!	Copyright (c) 2009, Mark A. Wieczorek
!	All rights reserved.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	use SHTOOLS, only : SHExpandDH, PlmIndex, PlmBar, SHCilmToVector, YilmIndex
	
	implicit none
		
	real*8, intent(out) ::	Dij(:,:)
	integer, intent(in) ::	dh_mask(:,:), n_dh, sampling, lmax
	integer ::		nlat, nlong, lmax_dh, astat, i, j, k, l, m, max_mask, min_mask
	real*8, allocatable ::	f(:,:), plm(:), clm(:,:,:), vec(:)
	real*8 ::		colat, lat_int, temp(2*n_dh), lon, lon_int

	
	if (size(Dij(:,1)) < (lmax+1)**2 .or. size(Dij(1,:)) < (lmax+1)**2) then
		print*, "Error --- ComputeDMap"
		print*, "Dij must be dimesioned as ((LMAX+1)**2, (LMAX+1)**2)."
		print*, "Dij is dimensioned as ", size(Dij(:,1)), size(Dij(1,:))
		stop
	endif
	
	if (mod(n_dh,2) /= 0) then
		print*, "Error --- ComputeDMap"
		print*, "Number of samples in latitude must be even for the Driscoll and Healy sampling theorem."
		print*, "N_DH = ", n_dh
		stop
	endif
	
	nlat = n_dh
	lat_int = acos(-1.0d0) / dble(nlat)
	if (sampling == 1) then
		nlong = nlat
		lon_int = 2.0d0*lat_int
	elseif (sampling == 2) then
		nlong = 2*nlat
		lon_int = lat_int
	else 
		print*, "Error --- ComputeDMap"
		print*, "SAMPLING must be either 1 (equally sampled) or 2 (equally spaced)."
		print*, "SAMPLING = ", sampling
		stop
	endif
	
	if (size(dh_mask(:,1)) < nlat .or. size(dh_mask(1,:)) < nlong) then
		print*, "Error --- ComputeDMap"
		print*, "DH_MASK must be dimensioned as ", nlat, nlong
		print*, "Dimensions of DH_MASK are ", size(dh_mask(:,1)), size(dh_mask(1,:))
		stop
	endif
	
	lmax_dh = n_dh/2 - 1
	
	if (lmax_dh < lmax) then
		print*, "Error --- ComputeDMap"
		print*, "The effective bandwith of the input grid DH_MASK must be greater or equal than LMAX."
		print*, "LMAX = ", lmax
		print*, "Effective bandwidth of DH_MASK = ", lmax_dh
		stop
	endif
	
	allocate(f(nlat,nlong), stat = astat)
	if (astat /= 0) then
		print*, "Error --- ComputeDMap"
		print*, "Problem allocating memory for grid F(NLAT, NLONG)."
		print*, "NLAT = ", nlat
		print*, "NLONG = ", nlong
		stop
	endif
	
	allocate(plm((lmax+1)*(lmax+2)/2), stat = astat)
	if (astat /= 0) then
		print*, "Error --- ComputeDMap"
		print*, "Problem allocating memory for grid PLM((LMAX+1)*(LMAX+2)/2)."
		print*, "LMAX = ", lmax
		stop
	endif
	
	allocate(clm(2, lmax+1, lmax+1), stat = astat)
	if (astat /= 0) then
		print*, "Error --- ComputeDMap"
		print*, "Problem allocating memory for CLM(2, LMAX+1, LMAX+1)."
		stop
	endif
	
	allocate(vec((lmax+1)**2), stat = astat)
	if (astat /= 0) then
		print*, "Error --- ComputeDMap"
		print*, "Problem allocating memory for VEC((LMAX+1)**2)."
		print*, "LMAX = ", lmax
		stop
	endif
	
	min_mask = minval(dh_mask(1:nlat, 1:nlong))
	max_mask = maxval(dh_mask(1:nlat, 1:nlong))
	if (min_mask < 0 .or. max_mask > 1) then
		print*, "Error --- ComputeDMap"
		print*, "DH_MASK must consist of 0s exterior to the concentration region"
		print*, "and 1s inside of the concentration region."
		print*, "Minimum value of DH_MASK = ", min_mask
		print*, "Maximum value of DH_MASK = ", max_mask
		stop
	endif

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!
	!	Compute elements of Dij. Start at high degrees, and work towards
	!	low degrees. Dij is symmetric.
	!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	dij = 0.0d0
	
	do l=lmax, 0, -1
		print*, "ComputeDMap --- Working on degree", l
		do m=0, l	
			do j=1, 2
			
				if (m==0 .and. j==2) cycle
				
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				!
				!	Make Map of Yilm * mask
				!
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				
				do k=1, nlong
					lon = dble(k-1)*lon_int
					if (j==1) then
						temp(k) = cos(dble(m)*lon)
					else
						temp(k) = sin(dble(m)*lon)
					endif
				enddo
				
				do k=1, nlat
					colat = dble(k-1)*lat_int
					call plmbar(plm, l, cos(colat))
					f(k, 1:nlong) = plm(plmindex(l,m)) * temp(1:nlong) &
							* dble(dh_mask(k, 1:nlong))
				enddo
				
				i = YilmIndex(j, l, m)
				
				call SHExpandDH(f, n_dh, clm, lmax_dh, sampling = sampling, lmax_calc = l)
				
				call SHCilmToVector(clm, vec, l)
				
				dij(i, 1:(l+1)**2) = vec(1:(l+1)**2)
				
				if (l /= 0) dij(1:(l+1)**2 - 1, i) = vec(1:(l+1)**2 - 1) 
				
			enddo
		enddo
	enddo
	
	
	call plmbar(plm, -1, cos(colat))

	deallocate(f)
	deallocate(plm)
	deallocate(clm)
	deallocate(vec)
		
end subroutine ComputeDMap
