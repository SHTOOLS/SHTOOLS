subroutine SHBiasAdmitCorr(Sgt, Sgg, Stt, lmax, tapers, lwin, K, admit, corr, mtdef, taper_wt)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	This subroutine will take as input the three global (unbiassed) cross-power spectra
!	of two functions G and T (Sgt, Sgg, and Stt), bias these using the first K localization
!	windows, and then calculated the expected biassed admittance (Sgt/Stt) and correlation
!	(Sgt/sqrt(Sgg Stt)) spectra. Two manners of calculating the localized admittance and correlation 
!	are possible according to the optional parameter MTDEF. In one case, the multitaper cross-power 
!	spectra are calculated, and from these, the admittance and correlation. In the second, the 
!	admittance and correlation are calculated for each taper, and these are then averaged.
!
!	The admittance and correlation functions are only calculated up to a maximum degree of 
!	lmax - lwin. This implicitly assumes that the cross-power spectra beyond lmax are simply
!	unknown, and not zero.
!	
!	All cross-power spectra must be constructed using 4-pi normalized or orthonormalized
!	spherical harmonic functions.
!
!	Calling Parameters
!		INPUT
!			Sgt		Global cross-power spectrum of G and T.
!			Sgg		Global power spectrum of G.
!			Stt		Global power spectrum of T.
!			lmax		Maximum spherical harmonic degree of the intput cross-power spectra.
!			tapers		A matrix of tapers obtained from SHReturnTapers.
!			lwin		Spectral bandwidth of the localizing window.
!			K		Number of tapers to use in the multitaper spectral estimations.
!		OUTPUT
!			admit		Biassed admittance between the localized G and T: Sgt/Stt. 
!			corr		Biassed correlation of the two functions: Sgt / sqrt(Sgg Stt)
!		OPTIONAL (IN)
!			mtdef		1 (default): Calculate multitaper cross-spectral estimates, and use these
!					to calculate a single admittance and correlations.
!					2: Calculate the admittance and correlation using each individual taper, and
!					then average these to get the admittance and correlation.
!			taper_wt	Vector of length numk corresponding to the weights applied
!					to each spectal estimate. The sum of taper_wt will be normalized
!					to unity.
!
!	Dependencies: SHBias, SHBiasK
!
!	Written by Mark Wieczorek April 2009
!
!	Copyright (c) 2009, Mark A. Wieczorek
!	All rights reserved.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	use SHTOOLS, only: SHBias, SHBiasK
	
	implicit none
	
	real*8, intent(in) ::	sgt(:), sgg(:), stt(:), tapers(:,:)
	integer, intent(in) ::	lmax, lwin, K
	real*8, intent(out) ::	admit(:), corr(:)
	integer, intent(in), optional ::	mtdef
	real*8, intent(in), optional ::	taper_wt(:)
	real*8 ::	sgt_bias(lmax-lwin+1), sgg_bias(lmax-lwin+1), stt_bias(lmax-lwin+1), shh(lwin+1)
	integer ::	lmax_calc, def, i
	
	
	if (size(sgt) < lmax+1) then
		print*, "Error --- SHBiasAdmitCorr"
		print*, "SGT must be dimensioned as (LMAX+1), or more, where LMAX is ", lmax
		print*, "Input array is dimensioned ", size(sgt)
		stop
	elseif (size(sgg) < lmax+1) then
		print*, "Error --- SHBiasAdmitCorr"
		print*, "SGG must be dimensioned as (LMAX+1), or more, where LMAX is ", lmax
		print*, "Input array is dimensioned ", size(sgg)
		stop
	elseif (size(stt) < lmax+1) then
		print*, "Error --- SHBiasAdmitCorr"
		print*, "STT must be dimensioned as (LMAX+1), or more, where LMAX is ", lmax
		print*, "Input array is dimensioned ", size(stt)
		stop
	elseif (size(admit) < lmax-lwin+1) then
		print*, "Error --- SHBiasAdmitCorr"
		print*, "ADMIT must be dimensioned as (LMAX-LWIN+1) where LMAX and LWIN are ", lmax, lwin
		print*, "Input array is dimensioned ", size(admit)
		stop
	elseif (size(corr) < lmax-lwin+1) then
		print*, "Error --- SHBiasAdmitCorr"
		print*, "CORR must be dimensioned as (LMAX-LWIN+1) where LMAX and LWIN are ", lmax, lwin
		print*, "Input array is dimensioned ", size(corr)
		stop
	elseif (size(tapers(:,1)) < lwin+1 .or. size(tapers(1,:)) < K) then
		print*, "Error --- SHBiasAdmitCorr"
		print*, "TAPERS must be dimensioned as (LWIN+1, K) where lwin and K are ", lwin, K
		print*, "Iinput array is dimensioned as ", size(tapers(:,1)), size(tapers(1,:))
		stop
	endif
	
	if (present(taper_wt)) then
		if (size(taper_wt) < K) then	
			print*, "Error --- SHBiasAdmitCorr"
			print*, "TAPER_WT must be dimensioned as K where K is ", K
			print*, "Input array is dimensioned as ", size(taper_wt)
			stop
		elseif (sum(taper_wt(1:K)) /= 1.0d0) then
			print*, "Error --- SHBiasAdmitCorr"
			print*, "TAPER_WT must sum to unity."
			print*, "Input array sums to ", sum(taper_wt(1:K))
			stop
		endif
	endif
	
	if (present(mtdef)) then
		if (mtdef /= 1 .and. mtdef /= 2) then
			print*, "SHBiasAdmitCorr --- Error"
     			print*, "MTDEF must be 1 or 2."
     			print*, "Input value is ", mtdef
     			stop
     		else
     			def = mtdef
     		endif
     	else
     		def = 1
     	endif
     	
     	lmax_calc = lmax - lwin
     	
   	admit = 0.0d0
   	corr = 0.0d0
   	
   	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   	!
   	!	Calculate baissed cross-power spectra and biassed admittances and correlations.
   	!
   	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   	
   	if (def == 1) then
   	
   		if (present(taper_wt)) then
  	
   			call SHBiasK (tapers, lwin, K, sgt, lmax, sgt_bias(1:lmax_calc+1), taper_wt = taper_wt, save_cg = 1)
   			call SHBiasK (tapers, lwin, K, sgg, lmax, sgg_bias(1:lmax_calc+1), taper_wt = taper_wt, save_cg = 1)
   			call SHBiasK (tapers, lwin, K, stt, lmax, stt_bias(1:lmax_calc+1), taper_wt = taper_wt, save_cg = 1)

   		else
   	
   			call SHBiasK (tapers, lwin, K, sgt, lmax, sgt_bias(1:lmax_calc+1), save_cg = 1)
   			call SHBiasK (tapers, lwin, K, sgg, lmax, sgg_bias(1:lmax_calc+1), save_cg = 1)
   			call SHBiasK (tapers, lwin, K, stt, lmax, stt_bias(1:lmax_calc+1), save_cg = 1)
   			
   		endif
   		
   		admit(1:lmax_calc+1) = sgt_bias(1:lmax_calc+1) / stt_bias(1:lmax_calc+1)
   		corr(1:lmax_calc+1) = sgt_bias(1:lmax_calc+1) / sqrt(stt_bias(1:lmax_calc+1)) / sqrt(sgg_bias(1:lmax_calc+1))
   	
   	elseif (def == 2) then
   	
   		do i=1, K
   		
   			shh(1:lwin+1) = tapers(1:lwin+1, i)**2
   			
   			call SHBias (shh, lwin, sgt, lmax, sgt_bias(1:lmax_calc+1), save_cg = 1)
   			call SHBias (shh, lwin, sgg, lmax, sgg_bias(1:lmax_calc+1), save_cg = 1)
   			call SHBias (shh, lwin, stt, lmax, stt_bias(1:lmax_calc+1), save_cg = 1)
   			
   			admit(1:lmax_calc+1) = admit(1:lmax_calc+1) + sgt_bias(1:lmax_calc+1) / stt_bias(1:lmax_calc+1)
   			corr(1:lmax_calc+1) = corr(1:lmax_calc+1) + sgt_bias(1:lmax_calc+1) / &
   				sqrt(stt_bias(1:lmax_calc+1)) / sqrt(sgg_bias(1:lmax_calc+1))
   				
   		enddo
   		
   		admit(1:lmax_calc+1) = admit(1:lmax_calc+1) / dble(K)
   		corr(1:lmax_calc+1) = corr(1:lmax_calc+1) / dble(K)
   	
   	endif

end subroutine SHBiasAdmitCorr