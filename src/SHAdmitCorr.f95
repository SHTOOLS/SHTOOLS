subroutine SHAdmitCorr(G, T, lmax, admit, corr, admit_error)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	This subroutine will compute the admittance and correlation between the two input
!	spherical harmonic coefficients. The optional "error" on the admittance estimates 
!	assumes that the two fields are perfectly correlated, which may not necessarily true
!	for certain models (such as Forsyth's loading model). 
!	
!	Assuming that the two input fields are related by 
!
!	G = Z T + N
!
!	where N is noise, the admittance is calcalated as 
!
!		Z = < G T > / < T T >,
!
!	the correlation is 
!
!		correlation = < G T > / sqrt(< G G >) / sqrt(< T T >).
!
!	Assuming that the two fields are perfectly correlated, and that the lack of coherence is
!	a result of noise, the uncertainty on the admittance can be calculaed from
!
!		sigma^2 = ( < G G > / < T T > ) ( 1 - gamma^2) / (2 l).
!
!	< > signifies an average of over all m. Since only ratios 
!	are being used, < > is implemented only as sums over all m (i.e., by the cross 
!	spectra of the fields). Note that the correlation can possess values between -1 and 1, 
!	in contrast to the "coherence" (or coherence-squared; the terminology is ambiguous) 
!	which only possesses values between 0 and 1.
!
!	Note that in order for the magnitude of the admittance to be correct, the 
!	input coefficients must be multiplied by all constants. For instance, in the 
!	case of gravity-topography admittances, the potential coefficients must first be 
!	multiplied by GM (l+1) (R/r)**(l+2) where R is the refernce radius of the coefficients,
!	and r is the radius at which the field is evaluated.
!
!	Calling Parameters
!		IN
!			G		First input spherical harmonic coefficients.
!			T		Second input spherical harmonic coefficients.
!			lmax		Maximum spherical harmonic degree of the input fields.
!		OUT
!			admit		G/T admittance spectra (length lmax+1).
!			corr		G/T correlation spectra (length lmax+1).
!		OPTIONAL
!			admit_error	G/T addittance error spectra (length lmax+1)
!					assuming that the two fields perfectly correlated in 
!					the abscence of noise.
!
!	Dependencies: SHPowerSpectrum, SHCrossPowerSpectrum
!
!	Written by Mark Wieczorek (June 2004)
!
!	Copyright (c) 2005-2006, Mark A. Wieczorek
!	All rights reserved.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	use SHTOOLS, only: SHPowerSpectrum, SHCrossPowerSpectrum

	implicit none
	
	real*8, intent(in) ::	G(:,:,:), T(:,:,:)
	integer, intent(in) ::	lmax
	real*8, intent(out) ::	admit(:), corr(:)
	real*8, intent(out), optional ::	admit_error(:)
	real*8 ::	gt(lmax+1), gg(lmax+1), tt(lmax+1)
	integer ::	l, l1
	
	if (size(G(:,1,1)) < 2 .or. size(G(1,:,1)) < lmax+1 .or. size(G(1,1,:)) < lmax+1) then
		print*, "Error --- SHAdmitCorr"
		print*, "G must be dimensioned as (2, LMAX+1, LMAX+1) where LMAX is ", lmax
		print*, "Input dimension is ", size(G(:,1,1)), size(G(1,:,1)), size(G(1,1,:))
		stop
	elseif (size(T(:,1,1)) < 2 .or. size(T(1,:,1)) < lmax+1 .or. size(T(1,1,:)) < lmax+1) then
		print*, "Error --- SHAdmitCorr"
		print*, "T must be dimensioned as (2, LMAX+1, LMAX+1) where LMAX is ", lmax
		print*, "Input dimension is ", size(T(:,1,1)), size(T(1,:,1)), size(T(1,1,:))
		stop
	elseif (size(admit) < lmax+1) then
		print*, "Error --- SHAdmitCorr"
		print*, "ADMIT must be dimensioned as (LMAX+1) where LMAX is ", lmax
		print*, "Input dimension is ", size(admit)
		stop
	elseif (size(corr) < lmax+1) then
		print*, "Error --- SHAdmitCorr"
		print*, "CORR must be dimensioned as (LMAX+1) where LMAX is ", lmax
		print*, "Input dimension is ", size(corr)
		stop
	endif
	
	if (present(admit_error)) then
		if (size(admit_error) < lmax+1) then
			print*, "Error --- SHAdmitCorr"
			print*, "ADMIT_ERROR must be dimensioned as (LMAX+1) where LMAX is ", lmax
			print*, "Input dimension is ", size(admit_error)
			stop
		endif
	endif
	
	admit = 0.0d0
	corr = 0.0d0

	call SHCrossPowerSpectrum(G, T, lmax, gt) ! < G T > 

	call SHPowerSpectrum(G, lmax, gg) ! < G G >
	
	call SHPowerSpectrum(T, lmax, tt) ! < T T >
	
	admit(1:lmax+1) = gt(1:lmax+1) / tt(1:lmax+1)
	
	corr(1:lmax+1) = gt(1:lmax+1) / sqrt( tt(1:lmax+1) * gg(1:lmax+1) )
	
	if (present(admit_error)) then
		admit_error = 0.0d0
		do l = 1, lmax
			l1 = l+1
			admit_error(l1) = gg(l1)*(1.0d0 - corr(l1)**2) / ( tt(l1) * dble(2*l) )
		enddo
		admit_error(1:lmax+1) = sqrt(admit_error(1:lmax+1))
	endif
	
end subroutine SHAdmitCorr

