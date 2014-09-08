subroutine SHMultiTaperCSE(mtse, sd, sh1, lmax1, sh2, lmax2, tapers, taper_order, lmaxt, K, alpha, &
	lat, lon, taper_wt, norm, csphase)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	This subroutine will calculate the multitaper spectrum estimate utilizing
!	the first K tapers. The standard error is calculated using an unbiased
!	estimate of the sample variance.
!
!	Calling Parameters
!
!		IN
!			sh1		Input spherical harmonic file.
!			lmax1		Maximum degree of sh1.
!			sh2		Input spherical harmonic file.
!			lmax2		Maximum degree of sh2.
!			tapers		The eigenvector matrix returned from a program
!					such as SHReturnAllTapers, where each column corresponds
!					to the coefficients for a window for a single non-zero
!					angular order.
!			taper_order	Angular order of the tapers. 
!			lmaxt		Maximum degree of the eigentapers.
!			K		Number of tapers to use in the multitaper spectral estimation.
!		OUT
!			mtse		Multitaper spectrum estimate, valid up to and including a 
!					maximum degree lmax-lmaxt.
!			sd		Standard error of the multitaper spectrum estimate.
!		OPTIONAL (IN)		
!			alpha		Euler angles used to rotate the localizing windows.
!			lat 		Latitude to perform localized analysis (degrees).
!			lon		Longitude to perform localized analysis (degrees).
!			taper_wt	Weight to be applied to each direct spectral estimate. This
!					should sum to unity.
!			csphase:	1: Do not include the phase factor of (-1)^m
!					-1: Apply the phase factor of (-1)^m.
!			norm:		Normalization to be used when calculating Legendre functions
!						(1) "geodesy" (default)
!						(2) Schmidt
!						(3) unnormalized
!						(4) orthonormalized
!
!		If the optional parameter alpha (or lat and lon) is input, then the spherical harmonic coefficients of
!		the localizing windows will be rotated accordingly. To rotate a window originally centered
!		at the north pole to (lat, long) given in degrees, use
!
!			alpha(1) = 0.0
!			alpha(2) = -(90.0d0 - lat)*pi/180.0d0
!			alpha(3) =  -long*pi/180.0d0
!
!		See documentation in file ShRotateRealCoef for further information on spherical
!		harmonic rotations.
!
!	Dependencies: SHMultiply, SHCrossPowerSpectrum, SHRotateRealCoef, djpi2, CSPHASE_DEFAULT
!
!	Written by Mark Wieczorek (April 2005)
!		April 2009: Modified routine so that the rotation matrices are only calculated once, instead of
!		for all K tapers.
!
!	Copyright (c) 2009, Mark A. Wieczorek
!	All rights reserved.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	use SHTOOLS, only: SHMultiply, SHCrossPowerSpectrum, SHRotateRealCoef, djpi2, CSPHASE_DEFAULT
	
	implicit none
	
	real*8, intent(out) ::	mtse(:), sd(:)
	real*8, intent(in) ::	sh1(:,:,:), sh2(:,:,:), tapers(:,:)
	integer, intent(in) ::	lmax1, lmax2, lmaxt, K, taper_order(:)
	real*8, intent(in), optional :: alpha(:), lat, lon, taper_wt(:)
	integer, intent(in), optional ::	csphase, norm
	integer ::	i, l, lmax, phase, mnorm, astat(5)
	real*8 ::	se(lmax1-lmaxt+1,K), x(3), pi, factor
	real*8, allocatable ::	shwin(:,:,:), shloc1(:,:,:),  shloc2(:,:,:), dj(:,:,:), shwinrot(:,:,:)
	pi = acos(-1.0d0)
	
	lmax = min(lmax1, lmax2)
	
	if (size(sh1(:,1,1)) < 2 .or. size(sh1(1,:,1)) < lmax+1 .or. size(sh1(1,1,:)) < lmax+1) then
		print*, "ERROR - SHMultiTaperCSE"
		print*, "SH1 must be dimensioned (2,LMAX+1, LMAX+1) where LMAX is ", lmax
		print*, "Input array is dimensioned ", size(sh1(:,1,1)), size(sh1(1,:,1)), size(sh1(1,1,:))
		stop
	elseif (size(sh2(:,1,1)) < 2 .or. size(sh2(1,:,1)) < lmax+1 .or. size(sh2(1,1,:)) < lmax+1) then
		print*, "ERROR - SHMultiTaperCSE"
		print*, "SH2 must be dimensioned (2,LMAX+1, LMAX+1) where lmax is ", lmax
		print*, "Input array is dimensioned ", size(sh2(:,1,1)), size(sh2(1,:,1)), size(sh2(1,1,:))
		stop
	elseif(size(tapers(:,1)) < lmaxt+1 .or. size(tapers(1,:)) < K) then
		print*, "ERROR - SHMultiTaperCSE"
		print*, "TAPERS must be dimensioned (LMAXT+1, K) where LMAXT and K are, ",lmaxt, K
		print*, "Input array is dimensioned ", size(tapers(:,1)), size(tapers(1,:))
		stop
	elseif(size(taper_order) < K) then
		print*, "ERROR --- SHMutltiTaperCSE"
		print*, "TAPER_ORDER must be dimensioned as (K) where K is ", K
		print*, "Input dimension of array is ", size(taper_order)
		stop
	elseif(size(mtse) < lmax-lmaxt+1) then
		print*, "Error --- SHMultiTaperCSE"
		print*, "MTSE must be dimensioned as (LMAX-LMAXT+1) where LMAX and LMAXT are ", lmax, lmaxt
		print*, "Input dimension of array is ", size(mtse)
		stop
	elseif(size(sd) < lmax-lmaxt+1) then
		print*, "Error --- SHMultiTaperCSE"
		print*, "SD must be dimensioned as (LMAX-LMAXT1) where LMAX and LMAXT are ", lmax, lmaxt
		print*, "Input dimension of array is ", size(sd)
		stop
	elseif (lmax < lmaxt) then
		print*, "Error --- SHMultiTaperCSE"
		print*, "LMAX must be larger than LMAXT."
		print*, "Input valuse of LMAX and LMAXT are ", lmax, lmaxt
		stop
	endif
	
	if(present(taper_wt)) then
		if (size(taper_wt) < K) then
			print*, "ERROR --- SHMultiTaperCSE"
			print*, "TAPER_WT must be dimensioned as (K) where K is ", K
			print*, "Input dimension of array is ", size(taper_wt)
			stop
		endif
	endif
     	
	if (present(norm)) then
		if (norm > 4 .or. norm < 1) then
			print*, "Error - SHMultiTaperCSE"
			print*, "Parameter NORM must be 1 (geodesy), 2 (Schmidt), 3 (unnormalized), or 4 (orthonormalized)."
			print*, "Input value is ", norm
			stop
		endif
		mnorm = norm
	else
		mnorm = 1
	endif

	if (present(csphase)) then
     		if (csphase /= -1 .and. csphase /= 1) then
     			print*, "SHMultiTaperCSE --- Error"
     			print*, "CSPHASE must be 1 (exclude) or -1 (include)."
     			print*, "Input value is ", csphase
     			stop
     		else
     			phase = csphase
     		endif
     	else
     		phase = CSPHASE_DEFAULT
     	endif

	
	if (present(lat) .and. present(lon)) then
		x(1) = 0.0d0
		x(2) = -(90.0d0 - lat)*pi/180.0d0
		x(3) =  -lon*pi/180.0d0
		
	elseif (present(alpha)) then
		if (size(alpha) < 3) then
			print*, "Error --- SHMultiTaperCSE"
			print*, "ALPHA must be dimensioned as (3)."
			print*, "Input array is dimensioned as ", size(alpha)
			stop
		endif

		x = alpha
		
	endif
	
	allocate(shwin(2,lmaxt+1,lmaxt+1), stat = astat(1))
	allocate(shloc1(2, lmax1+lmaxt+1, lmax1+lmaxt+1), stat= astat(2))
	allocate(shloc2(2, lmax2+lmaxt+1, lmax2+lmaxt+1), stat = astat(3))
	allocate(dj(lmaxt+1,lmaxt+1,lmaxt+1), stat = astat(4))
	allocate(shwinrot(2,lmaxt+1,lmaxt+1), stat = astat(5))
	if (sum(astat(1:5)) /= 0) then
		print*, "Error --- SHMultiTaperCSE"
		print*, "Problem allocating arrays SHWIN, SHLOC1, SHLOC2, DJ and SHWINROT", &
			astat(1), astat(2), astat(3), astat(4), astat(5)
		stop
	endif

	mtse = 0.0d0
	sd = 0.0d0

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!
	!	Calculate localized power spectra
	!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	if (present(alpha) .or. (present(lat) .and. present(lon)) ) then
		call djpi2(dj, lmaxt)	! Create rotation matrix used in the rotation routine.
	endif

	do i=1, K
		shwin = 0.0d0
		if (taper_order(i) < 0) then
			shwin(2,1:lmaxt+1,abs(taper_order(i))+1) = tapers(1:lmaxt+1,i)
		else
			shwin(1,1:lmaxt+1,taper_order(i)+1) = tapers(1:lmaxt+1,i)
		endif
		
		if (present(alpha) .or. (present(lat) .and. present(lon)) ) then
			call SHRotateRealCoef(shwinrot, shwin, lmaxt, x, dj)
			shwin = shwinrot
		endif
		
		call SHMultiply(shloc1, sh1, lmax1, shwin, lmaxt, csphase = phase, norm = mnorm, precomp = 0)
		call SHMultiply(shloc2, sh2, lmax2, shwin, lmaxt, csphase = phase, norm = mnorm, precomp = 0)
		call SHCrossPowerSpectrum(shloc1, shloc2, lmax-lmaxt, se(:,i))
		
	enddo
	
	
	if (present(taper_wt)) then
	
		factor = sum(taper_wt(1:K))**2 - sum(taper_wt(1:K)**2)
		factor = factor * sum(taper_wt(1:K))
		factor= sum(taper_wt(1:K)**2) / factor
	
		do l=0, lmax-lmaxt, 1
			mtse(l+1) = dot_product(se(l+1,1:K), taper_wt(1:K)) / sum(taper_wt(1:K))
		
			if (K > 1) then 
				sd(l+1) = dot_product( (se(l+1,1:K) - mtse(l+1) )**2, taper_wt(1:K) ) * factor
			endif
		
		enddo
	
	else
	
		do l=0, lmax-lmaxt, 1
			mtse(l+1) = sum(se(l+1,1:K))/dble(K)
		
			if (K > 1) then
				sd(l+1) = sum( ( se(l+1,1:K) - mtse(l+1) )**2 ) / dble(K-1) /dble(K) ! standard error !
			endif
		
		enddo
		
	endif
	
	if (K > 1) sd(1:lmax-lmaxt+1) = sqrt(sd(1:lmax-lmaxt+1) )
	
	deallocate(shwin)
	deallocate(shloc1)
	deallocate(shloc2)
	deallocate(dj)
	deallocate(shwinrot)
	
end subroutine SHMultiTaperCSE

