subroutine SHLocalizedAdmitCorr(tapers, taper_order, lwin, lat, lon, g, t, lmax, admit, corr, K, &
	admit_error, corr_error, taper_wt, mtdef, k1linsig)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	Given two spherical harmonic fields (G and T), this routine will calculate
!	the localized admittance and correlation using the first space-concentrated
!	window of Wieczorek and Simons (2005). All functions must be 4-pi normalized, 
!	and exclude the Condon-Shortley phase factor. Two manners of calculating the 
!	localized admittance and correlation are possible according to the optional 
!	parameter MTDEF. In one case, the multitaper cross-power spectra are calculated,
!	and from these, the admittance and correlation. In the second, the admittance and 
!	correlation are calculated for each taper, and these are then averaged.
!
!	Calling Parameters
!		IN
!			tapers		A matrix of tapers obtained from SHReturnTapers.
!			taper_order	A vector continaing the angular order of each column of taper.
!			lwin		Spectral bandwidth of the localizing window.
!			lat, lon	Latitude and longitude that the window will be
!					rotated to, in DEGREES.
!			G, T		Input spherical harmonic fields.
!			K		Number of tapers to use in Multitaper spectral estimations.
!			lmax		Maximum spherical harmonic degree of the intput fields.
!		OUT
!			admit		Admittance between the localized G and T assuming 
!					that G = Z T.
!			
!			corr		Correlation of the two fields.
!		OPTIONAL (OUT)
!			admit_error	Error of the admittance (only when K>1)
!			corr_error	Error of the admittance (only when K>1)
!		OPTIONAL (IN)
!			mtdef		1 (default): Calculate multitaper cross-spectral estimates, and use these
!					to calculate a single admittance and correlations.
!					2: Calculate the admittance and correlation using each individual taper, and
!					then average these to get the admittance and correlation.
!			taper_wt	Weights to be applied to the spectral estimates. This can only be used
!					when MTDEF is 1.
!			k1linsig:	If present and equal to 1, the uncertainty in the addmittance
!					will be calculated by assuming the gravity and topography coefficients
!					are linearly correlated, and that any lack of correlation is the 
!					result of uncorrelated noise. This should only be used when one expects
!					the gravity and topography to be linearly correlated and when only a 
!					single taper is being used. This should not be used with a Forsyth 
!					type model that predicts a less than 1correlation coefficient. 
!					This is the square root of eq. 33 of Simons et al. 1997.
!
!
!	Notes:
!		1. The units of the output admittance will correspond to the units of the
!		spherical harmonic coefficients. If gravity/topography admittances are
!		desired, then either the gravity coefficients should be multiplied by
!		G M (l+1) * (r0/r)**(l+2) / r0**2 before calling this routine, or
!		the admittances should be multiplied by this factor afterwards.
!
!		2. The correlation is defined as Sgt / sqrt(Sgg Stt), which varies between -1 and 1.
!		To obtain the "coherence" (which is also sometimes referred to as the "coherence squared"),
!		just square this number.
!
!	Dependencies: 	SHMultitaperSE, SHMultiTaperCSE, CSPHASE_DEFAULT,  djpi2, SHRotateRealCoef,
!			SHMultiply, SHCrossPowerSpectrum, SHPowerSpectrum
!
!	Written by Mark Wieczorek July 2004
!		Completely rewritten June 1 2006 to take into account non-zonal windows.
!		April 2009: added the optional parameter MTDEF for defining two manners of calculating
!			mutlitapered admittances and correlations.
!		July 25, 2012. Replaced the calls to SHMultiTaperSE and SHMultiTaperCSE with their consituent
!			calls to reduce redundant overhead. Save array DJ for use in multiple calls to the routine.
!
!	Copyright (c) 2009, Mark A. Wieczorek
!	All rights reserved.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	use SHTOOLS, only:  djpi2, SHRotateRealCoef, SHMultiply, SHCrossPowerSpectrum, SHPowerSpectrum
			
	implicit none

	real*8, intent(in) ::	tapers(:,:), lat, lon, g(:,:,:), t(:,:,:)
	integer, intent(in) ::	lwin, lmax, K, taper_order(:)
	real*8, intent(out) ::	admit(:), corr(:)
	real*8, intent(out), optional ::	admit_error(:), corr_error(:)
	integer, intent(in), optional ::	mtdef, k1linsig
	real*8, intent(in), optional ::	taper_wt(:)
	integer ::		lmaxwin, l, def, astat(5), phase, norm, i
	integer, save :: 	first = 1, lwin_last = 0
	real*8 ::		pi, g_power(2,lwin+lmax+1), t_power(2,lwin+lmax+1), gt_power(2,lwin+lmax+1), x(3), &
				sgt(lmax-lwin+1, K), sgg(lmax-lwin+1, K), stt(lmax-lwin+1, K), admit_k(lmax-lwin+1, K), &
				corr_k(lmax-lwin+1, K), factor
	real*8, allocatable ::	shwin(:,:,:), shwinrot(:,:,:), shloc_g(:,:,:), shloc_t(:,:,:)	
	real*8, allocatable, save ::	dj(:,:,:)

	phase = 1
	norm = 1

	pi = acos(-1.0d0)
	lmaxwin = lmax + lwin

	if (present(k1linsig) .and. K /= 1) then
		if (k1linsig == 1) then
			print*, "SHlocalizedAdmitCorr --- Error"
			print*, "If K1LINSIG is present and equal to 1, K must be equal to 1."
			print*, "Input value of K1LINSIG is ", k1linsig
			stop
		endif
	endif

	if (size(admit) < lmax-lwin+1) then
		print*, "Error --- SHLocalizedAdmitCorr"
		print*, "ADMIT must be dimensioned as (LMAX-LWIN+1) where LMAX and LWIN are ", lmax, lwin
		print*, "Input array is dimensioned ", size(admit)
		stop
	elseif (size(corr) < lmax-lwin+1) then
		print*, "Error --- SHLocalizedAdmitCorr"
		print*, "CORR must be dimensioned as (LMAX-LWIN+1) where LMAX and LWIN are ", lmax, lwin
		print*, "Input array is dimensioned ", size(corr)
		stop
	elseif (size(tapers(:,1)) < lwin+1 .or. size(tapers(1,:)) < K) then
		print*, "Error --- SHLocalizedAdmitCorr"
		print*, "TAPERS must be dimensioned as (LWIN+1, K) where lwin and K are ", lwin, K
		print*, "Iinput array is dimensioned as ", size(tapers(:,1)), size(tapers(1,:))
		stop
	elseif (size(taper_order) < K ) then
		print*, "Error --- SHLocalizedAdmitCorr"
		print*, "TAPER_ORDER must be dimensioned as (K) where K is ", K
		print*, "Input array is dimensioned ", size(taper_order)
		stop
	elseif (size(g(:,1,1)) < 2 .or. size(g(1,:,1)) < lmax+1 .or. size(g(1,1,:)) < lmax+1) then
		print*, "Error --- SHLocalizedAdmitCorr"
		print*, "G must be dimensioned as (2, LMAX+1, LMAX+1) where lmax is ", lmax
		print*, "Input array is dimensioned ", size(g(:,1,1)), size(g(1,:,1)), size(g(1,1,:))
		stop
	elseif (size(t(:,1,1)) < 2 .or. size(t(1,:,1)) < lmax+1 .or. size(t(1,1,:)) < lmax+1) then
		print*, "Error --- SHLocalizedAdmitCorr"
		print*, "T must be dimensioned as (2, LMAX+1, LMAX+1) where lmax is ", lmax
		print*, "Input array is dimensioned ", size(t(:,1,1)), size(t(1,:,1)), size(t(1,1,:))
		stop
	endif
	
	if (present(admit_error)) then
		if (size(admit_error) < lmax - lwin+1) then
			print*, "Error ---SHLocalizedAdmitCorr"
			print*, "ADMIT_ERROR must be dimensioned as (LMAX-LWIN+1) where lmax and lmaxt are ", lmax, lwin
			print*, "Input array is dimensioned ", size(admit_error)
			stop	
		endif
	endif
	
	if (present(corr_error)) then
		if (size(corr_error) < lmax - lwin + 1) then
			print*, "Error ---SHLocalizedAdmitCorr"
			print*, "CORR_ERROR  must be dimensioned as (LMAX-LWIN+1) where lmax and lmaxt are ", lmax, lwin
			print*, "Input array is dimensioedn ", size(corr_error)
			stop	
		endif	
	endif
	
	if(present(taper_wt)) then
		if (size(taper_wt) < K) then
			print*, "Error --- SHLocalizedAdmitCorr"
			print*, "TAPER_WT must be dimensioned as (K) where K is ", K
			print*, "Input array has dimension ", size(taper_wt)
			stop
		endif
	endif
	
	if (present(mtdef)) then
		if (mtdef == 2 .and. present(taper_wt)) then
			print*, "Error --- SHLocalizedAdmitCorr"
			print*, "TAPER_WT can only be used when MTDEF is 1."
			stop
		endif
	endif
	
	if (present(mtdef)) then
		if (mtdef /= 1 .and. mtdef /= 2) then
			print*, "SHLocalizedAdmitCorr --- Error"
     			print*, "MTDEF must be 1 or 2."
     			print*, "Input value is ", mtdef
     			stop
     		else
     			def = mtdef
     		endif
     	else
     		def = 1
     	endif
     	
     	
   	admit = 0.0d0
   	corr = 0.0d0
   	
     	if (present(admit_error)) then
		admit_error = 0.0d0
	endif
	
	if (present(corr_error)) then
		corr_error = 0.0d0
	endif

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!
	!	Determine multitaper cross spectra estimates of G and T, and then calculate 
	!	the admittance and correlation. Errors for the latter are calculated by adding
	!	the error sources in quadrature. Taper weights can be specified in order
	!	to minimize the variance of the multitaper cross-spectral estimates.
	!
	!	Note that only the admittances with degrees greater than lwin and less than
	!	Lmax - lwin should be interpretted.
	!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	x(1) = 0.0d0
	x(2) = -(90.0d0 - lat)*pi/180.0d0
	x(3) =  -lon*pi/180.0d0
		
	if (first == 1) then
		lwin_last = lwin
		first = 0
		allocate(dj(lwin+1,lwin+1,lwin+1), stat = astat(1))
		if (astat(1) /= 0) then
			print*, "Error --- SHLocalizedAdmitCorr"
			print*, "Problem allocating array DJ", astat(1)
			stop
		endif
		dj = 0.0d0
		call djpi2(dj, lwin)
	endif
	
	if (lwin > lwin_last) then
		lwin_last = lwin
		deallocate(dj)
		allocate(dj(lwin+1,lwin+1,lwin+1), stat = astat(1))
		if (astat(1) /= 0) then
			print*, "Error --- SHLocalizedAdmitCorr"
			print*, "Problem allocating array DJ", astat(1)
			stop
		endif
		dj = 0.0d0
		call djpi2(dj, lwin)
	endif
	
	allocate(shwin(2,lwin+1,lwin+1), stat = astat(1))
	allocate(shwinrot(2,lwin+1,lwin+1), stat = astat(2))
	allocate(shloc_g(2, lmaxwin+1, lmaxwin+1), stat= astat(3))
	allocate(shloc_t(2, lmaxwin+1, lmaxwin+1), stat= astat(4))
		
	if (sum(astat(1:4)) /= 0) then
		print*, "Error --- SHLocalizedAdmitCorr"
		print*, "Problem allocating arrays SHWIN, SHWINROT, SHLOC_G, SHLOC_T", &
			astat(1), astat(2), astat(3), astat(4)
		stop
	endif
		
	if (def == 1) then
	
		do i=1, K
		
			shwin = 0.0d0
			if (taper_order(i) < 0) then
				shwin(2,1:lwin+1,abs(taper_order(i))+1) = tapers(1:lwin+1,i)
			else
				shwin(1,1:lwin+1,taper_order(i)+1) = tapers(1:lwin+1,i)
			endif
			
			call SHRotateRealCoef(shwinrot, shwin, lwin, x, dj)
		
			call SHMultiply(shloc_g, g, lmax, shwinrot, lwin, csphase = phase, norm = norm)
			call SHMultiply(shloc_t, t, lmax, shwinrot, lwin, csphase = phase, norm = norm)
			
			call SHCrossPowerSpectrum(shloc_g, shloc_t, lmax-lwin, sgt(:,i))
			call SHPowerSpectrum(shloc_g, lmax-lwin, sgg(:,i))
			call SHPowerSpectrum(shloc_t, lmax-lwin, stt(:,i))
			
			if (present(taper_wt)) then

				factor = sum(taper_wt(1:K))**2 - sum(taper_wt(1:K)**2)
				factor = factor * sum(taper_wt(1:K))
				factor= sum(taper_wt(1:K)**2) / factor
	
				do l=0, lmax-lwin, 1
					g_power(1,l+1) = dot_product(sgg(l+1,1:K), taper_wt(1:K)) / sum(taper_wt(1:K))
					t_power(1,l+1) = dot_product(stt(l+1,1:K), taper_wt(1:K)) / sum(taper_wt(1:K))
					gt_power(1,l+1) = dot_product(sgt(l+1,1:K), taper_wt(1:K)) / sum(taper_wt(1:K))
		
					if (K > 1) then 
						g_power(2,l+1) = dot_product( (sgg(l+1,1:K) - g_power(1,l+1) )**2, taper_wt(1:K) ) * factor
						t_power(2,l+1) = dot_product( (stt(l+1,1:K) - t_power(1,l+1) )**2, taper_wt(1:K) ) * factor
						gt_power(2,l+1) = dot_product( (sgt(l+1,1:K) - gt_power(1,l+1) )**2, taper_wt(1:K) ) * factor
					endif
		
				enddo
		
			else
			
				do l=0, lmax-lwin, 1
					g_power(1,l+1) = sum(sgg(l+1,1:K))/dble(K)
					t_power(1,l+1) = sum(stt(l+1,1:K))/dble(K)
					gt_power(1,l+1) = sum(sgt(l+1,1:K))/dble(K)
	
					if (K > 1) then 
						g_power(2,l+1) = sum( ( sgg(l+1,1:K) - g_power(1,l+1) )**2 ) / dble(K-1) / dble(K) ! standard error!
						t_power(2,l+1) = sum( ( stt(l+1,1:K) - t_power(1,l+1) )**2 ) / dble(K-1) / dble(K) ! standard error!
						gt_power(2,l+1) = sum( ( sgt(l+1,1:K) - gt_power(1,l+1) )**2 ) / dble(K-1) / dble(K) ! standard error!
					endif	
				
				enddo
	
			endif
			
			if (K > 1) then
				g_power(2,1:lmax-lwin+1) = sqrt(g_power(2,1:lmax-lwin+1) )
				t_power(2,1:lmax-lwin+1) = sqrt(t_power(2,1:lmax-lwin+1) )
				gt_power(2,1:lmax-lwin+1) = sqrt(gt_power(2,1:lmax-lwin+1) )
			endif
	
		enddo
		
		admit(1:lmax-lwin+1) = gt_power(1,1:lmax-lwin+1) / t_power(1,1:lmax-lwin+1)
		corr(1:lmax-lwin+1) =  gt_power(1,1:lmax-lwin+1) / sqrt(t_power(1,1:lmax-lwin+1)*g_power(1,1:lmax-lwin+1))

		if (K > 1 .and. present(admit_error)) then
			admit_error(1:lmax-lwin+1) = ( gt_power(2,1:lmax-lwin+1) / t_power(1,1:lmax-lwin+1) )**2 + &
				( gt_power(1,1:lmax-lwin+1) / t_power(1,1:lmax-lwin+1)**2 * t_power(2,1:lmax-lwin+1) )**2
			admit_error(1:lmax-lwin+1) = sqrt(admit_error(1:lmax-lwin+1))
		endif
	
		if (K > 1 .and. present(corr_error)) then
			corr_error(1:lmax-lwin+1) = gt_power(2,1:lmax-lwin+1)**2 / t_power(1,1:lmax-lwin+1) / g_power(1,1:lmax-lwin+1) + &
				( gt_power(1,1:lwin-lmax+1) * t_power(2,1:lwin-lmax+1) / sqrt(g_power(1,1:lmax-lwin+1)) &
				/ 2.0d0 / t_power(1,1:lmax-lwin+1)**(3.0d0/2.0d0) )**2 + &
				( gt_power(1,1:lwin-lmax+1) * g_power(2,1:lwin-lmax+1) / sqrt(t_power(1,1:lmax-lwin+1)) &
				/ 2.0d0 / g_power(1,1:lmax-lwin+1)**(3.0d0/2.0d0) )**2
			corr_error(1:lmax-lwin+1) = sqrt(corr_error(1:lmax-lwin+1))
		endif
	
		if (K==1 .and. present(k1linsig) .and. present(admit_error)) then
			admit_error = 0.0d0
			if (k1linsig == 1) then
				do l=1, lmax-lwin
					admit_error(l+1) = g_power(1,l+1)*(1.0d0 - corr(l+1)**2) / ( t_power(1,l+1) * dble(2*l) )
				enddo
				admit_error(1:lmax-lwin+1) = sqrt(admit_error(1:lmax-lwin+1))
			endif
		endif
		
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!
	!	Calculate the admittance and correlation for each individual taper. Then average
	!	these in order to get the multitaper estimates and uncertainties.
	!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	else
				
		do i=1, K
			shwin = 0.0d0
			if (taper_order(i) < 0) then
				shwin(2,1:lwin+1,abs(taper_order(i))+1) = tapers(1:lwin+1,i)
			else
				shwin(1,1:lwin+1,taper_order(i)+1) = tapers(1:lwin+1,i)
			endif
		
			call SHRotateRealCoef(shwinrot, shwin, lwin, x, dj)
		
			call SHMultiply(shloc_g, g, lmax, shwinrot, lwin, csphase = phase, norm=norm)
			call SHMultiply(shloc_t, t, lmax, shwinrot, lwin, csphase = phase, norm=norm)
			
			call SHCrossPowerSpectrum(shloc_g, shloc_t, lmax-lwin, sgt(:,i))
			call SHPowerSpectrum(shloc_g, lmax-lwin, sgg(:,i))
			call SHPowerSpectrum(shloc_t, lmax-lwin, stt(:,i))
			
			admit_k(1:lmax-lwin+1, i) =  sgt(1:lmax-lwin+1, i) / stt(1:lmax-lwin+1, i)
			corr_k(1:lmax-lwin+1, i) =  sgt(1:lmax-lwin+1, i) / sqrt(stt(1:lmax-lwin+1, i)) / sqrt(sgg(1:lmax-lwin+1, i))
			
		enddo
		
		do l=0, lmax-lwin, 1
			admit(l+1) = sum(admit_k(l+1,1:K)) / dble(K)
			corr(l+1) = sum(corr_k(l+1,1:K)) / dble(K)
		enddo
		
		if (present(admit_error) .or. present(corr_error)) then
		
			do l=0, lmax-lwin, 1
				
				if (present(admit_error)) then
				
					if (K > 1) then 
					
						admit_error(l+1) = sum( ( admit_k(l+1,1:K) - admit(l+1) )**2 ) / dble(K-1) / dble(K) ! standard error!
						admit_error(l+1) = sqrt(admit_error(l+1))

					elseif (K == 1 .and. present(k1linsig)) then
					
						if (k1linsig == 1) then
						
							if (l==0) then
								admit_error(1) = 0.0d0
							else
								admit_error(l+1) = sgg(l+1,1)*(1.0d0 - corr(l+1)**2) / ( stt(l+1,1) * dble(2*l) )
								admit_error(l+1) = sqrt(admit_error(l+1))
							endif
							
						endif
						
					endif
					
				endif
				
				if (present(corr_error)) then
				
					if (K > 1) then 
					
						corr_error(l+1) = sum( ( corr_k(l+1,1:K) - corr(l+1) )**2 ) / dble(K-1) / dble(K) ! standard error!
						corr_error(l+1) = sqrt(corr_error(l+1))
					endif
					
				endif
				
			enddo
			
		endif
		
	endif	
	
	deallocate(shwin)
	deallocate(shwinrot)
	deallocate(shloc_g)
	deallocate(shloc_t)
	
end subroutine SHLocalizedAdmitCorr

