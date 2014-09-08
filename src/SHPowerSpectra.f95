real*8 function SHPowerL(c, l)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	This function will compute the dimensionless power at
!	degree l corresponding to the spherical harmonic coefficients c(i, l, m)
!
!	Power = Sum_{m=0}^l ( C1lm**2 + C2lm**2 )
!
!	Calling Parameters
!		c	Spherical harmonic coefficients, dimensioned as 
!			(2, lmax+1, lmax+1).
!		l	Spherical harmonic degree to compute power.
!
!	Written by Mark Wieczorek 2004.
!
!	Copyright (c) 2005, Mark A. Wieczorek
!	All rights reserved.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	implicit none
	
	real*8, intent(in) :: c(:,:,:)
	integer, intent(in) :: l
	integer i, m, l1, m1
	
	l1 = l+1
	
	if (size(c(:,1,1)) < 2 .or. size(c(1,:,1)) < l1 .or. size(c(1,1,:)) < l1) then
		print*, "SHPowerL --- Error"
		print*, "C must be dimensioned as (2, L+1, L+1) where L is ", l
		print*, "Input array is dimensioned ", size(c(:,1,1)),  size(c(1,:,1)), size(c(1,1,:)) 
		stop
	endif

	SHPowerL =  c(1, l1, 1)**2	! m=0 term
	
	do m = 1, l, 1
		m1 = m+1
		do i=1, 2
			SHPowerL = SHPowerL + c(i, l1, m1)**2
		enddo
	enddo

end function SHPowerL


real*8 function SHPowerDensityL(c, l)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	This function will compute the dimensionless power per coefficient
!	(density) at degree l of the spherical harmonic coefficients c(i, l, m)
!
!	PowerSpectralDensity = Sum_{m=0}^l ( C1lm**2 + C2lm**2 ) / (2l+1)
!
!	Calling Parameters
!		c	Spherical harmonic coefficients, dimensioned as 
!			(2, lmax+1, lmax+1).
!		l	Spherical harmonic degree to compute power.
!
!	Written by Mark Wieczorek 2004.
!
!	Copyright (c) 2005, Mark A. Wieczorek
!	All rights reserved.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	implicit none
	real*8, intent(in) :: c(:,:,:)
	integer, intent(in) :: l
	integer i, m, l1, m1	
	
	l1 = l+1
	
	if (size(c(:,1,1)) < 2 .or. size(c(1,:,1)) < l1 .or. size(c(1,1,:)) < l1) then
		print*, "SHPowerDensityL --- Error"
		print*, "C must be dimensioned as (2, L+1, L+1) where L is ", l
		print*, "Input array is dimensioned ", size(c(:,1,1)),  size(c(1,:,1)), size(c(1,1,:)) 
		stop
	endif

	SHPowerDensityL = c(1, l1, 1)**2	! m=0 term

	do m = 1, l, 1
		m1 = m+1
		do i=1, 2
			SHPowerDensityL = SHPowerDensityL + c(i, l1, m1)**2
		enddo
	enddo
	
	SHPowerDensityL = SHPowerDensityL/dble(2*l+1)

end function SHPowerDensityL


real*8 function SHCrossPowerL(c1, c2, l)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	This function will compute the dimensionless cross power at
!	degree l for the spherical harmonic coefficients c1(i, l, m)
!	and c2(i,l,m).
!
!	CrossPower =  Sum_{m=0}^l ( A1lm*B1lm + A2lm*B2lm )
!
!	Calling Parameters
!		c1	Spherical harmonic coefficients, dimensioned as 
!			(2, lmax+1, lmax+1).
!		c2	Spherical harmonic coefficients, dimensioned as 
!			(2, lmax+1, lmax+1).
!		l	Spherical harmonic degree to compute power.
!
!	Written by Mark Wieczorek 2004.
!
!	Copyright (c) 2005, Mark A. Wieczorek
!	All rights reserved.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	implicit none
	real*8, intent(in) :: c1(:,:,:), c2(:,:,:)
	integer, intent(in) :: l
	integer i, m, l1, m1	
	
	l1 = l+1
	
	if (size(c1(:,1,1)) < 2 .or. size(c1(1,:,1)) < l1 .or. size(c1(1,1,:)) < l1) then
		print*, "SHCrossPowerL --- Error"
		print*, "C1 must be dimensioned as (2, L+1, L+1) where L is ", l
		print*, "Input array is dimensioned ", size(c1(:,1,1)),  size(c1(1,:,1)), size(c1(1,1,:)) 
		stop
	elseif (size(c2(:,1,1)) < 2 .or. size(c2(1,:,1)) < l1 .or. size(c2(1,1,:)) < l1) then
		print*, "SHCrossPowerL --- Error"
		print*, "C2 must be dimensioned as (2, L+1, L+1) where L is ", l
		print*, "Input array is dimensioned ", size(c2(:,1,1)),  size(c2(1,:,1)), size(c2(1,1,:)) 
		stop
	endif

	SHCrossPowerL = c1(1, l1, 1)*c2(1,l1,1)

	do m = 1, l, 1
		m1 = m+1
		do i=1, 2
			SHCrossPowerL = SHCrossPowerL + c1(i, l1, m1)*c2(i,l1,m1)
		enddo
	enddo

end function SHCrossPowerL


real*8 function SHCrossPowerDensityL(c1, c2, l)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	This function will compute the dimensionless cross power
!	(density) at degree l of the spherical harmonic coefficients c1(i, l, m)
!	and c2(i,l,m).
!
!	CrossPower =  Sum_{m=0}^l ( A1lm*B1lm + A2lm*B2lm ) / (2l+1)
!
!	Calling Parameters
!		c1	Spherical harmonic coefficients, dimensioned as 
!			(2, lmax+1, lmax+1).
!		c2	Spherical harmonic coefficients, dimensioned as 
!			(2, lmax+1, lmax+1).
!		l	Spherical harmonic degree to compute power.
!
!	Written by Mark Wieczorek 2004.
!
!	Copyright (c) 2005, Mark A. Wieczorek
!	All rights reserved.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	implicit none
	real*8, intent(in) :: c1(:,:,:), c2(:,:,:)
	integer, intent(in) :: l
	integer i, m, l1, m1	
	
	l1 = l+1
	
	if (size(c1(:,1,1)) < 2 .or. size(c1(1,:,1)) < l1 .or. size(c1(1,1,:)) < l1) then
		print*, "SHCrossPowerDensityL --- Error"
		print*, "C1 must be dimensioned as (2, L+1, L+1) where L is ", l
		print*, "Input array is dimensioned ", size(c1(:,1,1)),  size(c1(1,:,1)), size(c1(1,1,:)) 
		stop
	elseif (size(c2(:,1,1)) < 2 .or. size(c2(1,:,1)) < l1 .or. size(c2(1,1,:)) < l1) then
		print*, "SHCrossPowerDensityL --- Error"
		print*, "C2 must be dimensioned as (2, L+1, L+1) where L is ", l
		print*, "Input array is dimensioned ", size(c2(:,1,1)),  size(c2(1,:,1)), size(c2(1,1,:)) 
		stop
	endif

	SHCrossPowerDensityL = c1(1, l1, 1)*c2(1,l1,1)
	
	do m = 1, l, 1
		m1 = m+1
		do i=1, 2
			SHCrossPowerDensityL = SHCrossPowerDensityL + c1(i, l1, m1)*c2(i,l1,m1)
		enddo
	enddo
	
	SHCrossPowerDensityL = SHCrossPowerDensityL/dble(2*l+1)

end function SHCrossPowerDensityL


subroutine SHPowerSpectrum(c, lmax, spectra)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	This function will compute the dimensionless power spectrum
!	of the spherical harmonic coefficients c(i, l, m).
!
!	Power(l) = Sum_{m=0}^l ( C1lm**2 + C2lm**2 )
!
!
!	Calling Parameters
!		IN
!			c	Spherical harmonic coefficients, dimensioned as 
!				(2, lmax+1, lmax+1).
!			lmax	Spherical harmonic degree to compute power.
!		OUT
!			spectra	Array of length (lmax+1) containing the power
!				spectra of c.
!
!	Written by Mark Wieczorek 2004.
!
!	Copyright (c) 2005, Mark A. Wieczorek
!	All rights reserved.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	implicit none
	real*8, intent(in) :: c(:,:,:)
	integer, intent(in) :: lmax
	real*8, intent(out) ::	spectra(:)
	integer i, m, l1, m1, l
	
	if (size(c(:,1,1)) < 2 .or. size(c(1,:,1)) < lmax+1 .or. size(c(1,1,:)) < lmax+1) then
		print*, "SHPowerSpectrum --- Error"
		print*, "C must be dimensioned as (2, LMAX+1, LMAX+1) where LMAX is ", lmax
		print*, "Input array is dimensioned ", size(c(:,1,1)),  size(c(1,:,1)), size(c(1,1,:)) 
		stop
	elseif (size(spectra) < lmax+1) then
		print*, "SHPowerSpectrum --- Error"
		print*, "SPECTRA must be dimensioned as (LMAX+1) where LMAX is ", lmax
		print*, "Input vector has dimension ", size(spectra)
		stop
	endif
	
	spectra = 0.0d0
	
	do l=0, lmax
		l1 = l+1
		
		spectra(l1) = c(1, l1, 1)**2
		
		do m = 1, l, 1
			m1 = m+1
			
			do i=1, 2
				spectra(l1) = spectra(l1) + c(i, l1, m1)**2
			enddo
		enddo
	enddo

end subroutine SHPowerSpectrum


subroutine SHPowerSpectrumDensity(c, lmax, spectra)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	This function will compute the dimensionless power spectrum
!	density of the spherical harmonic coefficients c(i, l, m)
!
!	PowerSpectralDensity = Sum_{m=0}^l ( C1lm**2 + C2lm**2 ) / (2l+1)
!
!	Calling Parameters
!		IN
!			c	Spherical harmonic coefficients, dimensioned as 
!				(2, lmax+1, lmax+1).
!			lmax	Spherical harmonic degree to compute power.
!		OUT
!			spectra	Array of length (lmax+1) containing the power
!				spectra density of c.
!
!	Written by Mark Wieczorek 2004.
!
!	Copyright (c) 2005, Mark A. Wieczorek
!	All rights reserved.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	implicit none
	real*8, intent(in) :: c(:,:,:)
	integer, intent(in) :: lmax
	real*8, intent(out) ::	spectra(:)
	integer i, m, l1, m1, l
	
	if (size(c(:,1,1)) < 2 .or. size(c(1,:,1)) < lmax+1 .or. size(c(1,1,:)) < lmax+1) then
		print*, "SHPowerSpectrumDensity --- Error"
		print*, "C must be dimensioned as (2, LMAX+1, LMAX+1) where LMAX is ", lmax
		print*, "Input array is dimensioned ", size(c(:,1,1)),  size(c(1,:,1)), size(c(1,1,:)) 
		stop
	elseif (size(spectra) < lmax+1) then
		print*, "SHPowerSpectrumDensity --- Error"
		print*, "SPECTRA must be dimensioned as (LMAX+1) where LMAX is ", lmax
		print*, "Input vector has dimension ", size(spectra)
		stop
	endif

	spectra = 0.0d0
	
	do l=0, lmax
		l1 = l+1
		
		spectra(l1) = c(1, l1, 1)**2
		
		do m = 1, l, 1
			m1 = m+1
			
			do i=1, 2
				spectra(l1) = spectra(l1) + c(i, l1, m1)**2
			enddo
		enddo
		
		spectra(l1) = spectra(l1)/dble(2*l+1)
	enddo

end subroutine SHPowerSpectrumDensity


subroutine SHCrossPowerSpectrum(c1, c2, lmax, cspectra)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	This function will compute the dimensionless cross power spectrum
!	of the spherical harmonic coefficients c1(i, l, m) and c2(1,l,m).
!
!	CrossPower(l) =  Sum_{m=0}^l ( A1lm*B1lm + A2lm*B2lm )
!
!	Calling Parameters
!		IN
!			c1		Spherical harmonic coefficients, dimensioned as 
!					(2, lmax+1, lmax+1).
!			c2		Spherical harmonic coefficients, dimensioned as 
!					(2, lmax+1, lmax+1).
!			lmax		Spherical harmonic degree to compute power.
!		OUT
!			cspectra	Array of length (lmax+1) containing the cross power
!					spectra of c.
!
!	Written by Mark Wieczorek 2004.
!
!	Copyright (c) 2005, Mark A. Wieczorek
!	All rights reserved.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	implicit none
	real*8, intent(in) :: c1(:,:,:), c2(:,:,:)
	integer, intent(in) :: lmax
	real*8, intent(out) ::	cspectra(:)
	integer i, m, l1, m1, l
	
	if (size(c1(:,1,1)) < 2 .or. size(c1(1,:,1)) < lmax+1 .or. size(c1(1,1,:)) < lmax+1) then
		print*, "SHCrossPowerSpectrum --- Error"
		print*, "C1 must be dimensioned as (2, LMAX+1, LMAX+1) where lmax is", lmax
		print*, "Input array is dimensioned ", size(c1(:,1,1)),  size(c1(1,:,1)), size(c1(1,1,:)) 
		stop
	elseif (size(c2(:,1,1)) < 2 .or. size(c2(1,:,1)) < lmax+1 .or. size(c2(1,1,:)) < lmax+1) then
		print*, "SHCrossPowerSpectrum --- Error"
		print*, "C2 must be dimensioned as (2, LMAX+1, LMAX+1)"
		print*, "Input array is dimensioned ", size(c2(:,1,1)),  size(c2(1,:,1)), size(c2(1,1,:)) 
		stop
	elseif (size(cspectra) < lmax+1) then
		print*, "SHCrossPowerSpectrum --- Error"
		print*, "CSPECTRA must be dimensioned as (LMAX+1) where lmax is ", lmax
		print*, "Input vector has dimension ", size(cspectra)
		stop
	endif

	cspectra = 0.0d0
	
	do l=0, lmax
		l1 = l+1
		
		 cspectra(l1) = c1(1, l1, 1)*c2(1, l1, 1)

		do m = 1, l, 1
			m1 = m+1
			
			do i=1, 2
				cspectra(l1) = cspectra(l1) + c1(i, l1, m1)*c2(i, l1, m1)
			enddo
		enddo
	enddo

end subroutine SHCrossPowerSpectrum


subroutine SHCrossPowerSpectrumDensity(c1, c2, lmax, cspectra)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	This function will compute the dimensionless cross power spectrum
!	density of the spherical harmonic coefficients c1(i, l, m) and c2(i,l,m).
!
!	CrossPower =  Sum_{m=0}^l ( A1lm*B1lm + A2lm*B2lm ) / (2l+1)
!
!	Calling Parameters
!		IN
!			c1		Spherical harmonic coefficients, dimensioned as 
!					(2, lmax+1, lmax+1).
!			c2		Spherical harmonic coefficients, dimensioned as 
!					(2, lmax+1, lmax+1).
!			lmax		Spherical harmonic degree to compute power.
!		OUT
!			cspectra	Array of length (lmax+1) containing the cross power
!					spectral density of c.
!
!	Written by Mark Wieczorek 2004.
!
!	Copyright (c) 2005, Mark A. Wieczorek
!	All rights reserved.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	implicit none
	real*8, intent(in) :: c1(:,:,:), c2(:,:,:)
	integer, intent(in) :: lmax
	real*8, intent(out) ::	cspectra(:)
	integer i, m, l1, m1, l
	
	if (size(c1(:,1,1)) < 2 .or. size(c1(1,:,1)) < lmax+1 .or. size(c1(1,1,:)) < lmax+1) then
		print*, "SHCrossPowerSpectrumDensity --- Error"
		print*, "C1 must be dimensioned as (2, LMAX+1, LMAX+1) where lmax is", lmax
		print*, "Input array is dimensioned ", size(c1(:,1,1)),  size(c1(1,:,1)), size(c1(1,1,:)) 
		stop
	elseif (size(c2(:,1,1)) < 2 .or. size(c2(1,:,1)) < lmax+1 .or. size(c2(1,1,:)) < lmax+1) then
		print*, "SHCrossPowerSpectrumDensity --- Error"
		print*, "C2 must be dimensioned as (2, LMAX+1, LMAX+1)"
		print*, "Input array is dimensioned ", size(c2(:,1,1)),  size(c2(1,:,1)), size(c2(1,1,:)) 
		stop
	elseif (size(cspectra) < lmax+1) then
		print*, "SHCrossPowerSpectrumDensity --- Error"
		print*, "CSPECTRA must be dimensioned as (LMAX+1) where lmax is ", lmax
		print*, "Input vector has dimension ", size(cspectra)
		stop
	endif

	cspectra = 0.0d0
	
	do l=0, lmax
		l1 = l+1
		
		cspectra(l1) = c1(1, l1, 1)*c2(1, l1, 1)

		do m = 1, l, 1
			m1 = m+1
			
			do i=1, 2
				cspectra(l1) = cspectra(l1) + c1(i, l1, m1)*c2(i, l1, m1)
			enddo
			
		enddo
		
		cspectra(l1) = cspectra(l1)/dble(2*l+1)
		
	enddo

end subroutine SHCrossPowerSpectrumDensity

