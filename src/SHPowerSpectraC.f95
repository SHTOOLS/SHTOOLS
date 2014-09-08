real*8 function SHPowerLC(c, l)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	This function will compute the dimensionless power at
!	degree l corresponding to the complex spherical harmonic coefficients c(i, l, m)
!
!	Power(l) = Sum_{m=0}^l | C1lm | **2 + | C2lm | **2 )
!
!	Calling Parameters
!		c	Complex spherical harmonic coefficients, dimensioned as 
!			(2, lmax+1, lmax+1).
!		l	Spherical harmonic degree to compute power.
!
!	Written by Mark Wieczorek (May 2008).
!
!	Copyright (c) 2008, Mark A. Wieczorek
!	All rights reserved.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	implicit none
	
	complex*16, intent(in) :: c(:,:,:)
	integer, intent(in) :: l
	integer i, m, l1, m1
	
	l1 = l+1
	
	if (size(c(:,1,1)) < 2 .or. size(c(1,:,1)) < l1 .or. size(c(1,1,:)) < l1) then
		print*, "SHPowerLC --- Error"
		print*, "C must be dimensioned as (2, L+1, L+1) where L is ", l
		print*, "Input array is dimensioned ", size(c(:,1,1)),  size(c(1,:,1)), size(c(1,1,:)) 
		stop
	endif

	SHPowerLC = dble(c(1, l1, 1))**2 + aimag(c(1, l1, 1))**2	! m=0 term

	do m = 1, l, 1
		m1 = m+1
		do i=1, 2
			SHPowerLC = SHPowerLC + dble(c(i, l1, m1))**2 + aimag(c(i, l1, m1))**2
		enddo
	enddo

end function SHPowerLC


real*8 function SHPowerDensityLC(c, l)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	This function will compute the dimensionless power per coefficient
!	(density) at degree l of the complex spherical harmonic coefficients c(i, l, m)
!
!	PowerSpectralDensity(l) = Sum_{m=0}^l ( | C1lm |**2 + | C2lm |**2 ) / (2l+1)
!
!	Calling Parameters
!		c	Complex spherical harmonic coefficients, dimensioned as 
!			(2, lmax+1, lmax+1).
!		l	Spherical harmonic degree to compute power.
!
!	Written by Mark Wieczorek (May 2008).
!
!	Copyright (c) 2008, Mark A. Wieczorek
!	All rights reserved.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	implicit none
	complex*16, intent(in) :: c(:,:,:)
	integer, intent(in) :: l
	integer i, m, l1, m1	
	
	l1 = l+1
	
	if (size(c(:,1,1)) < 2 .or. size(c(1,:,1)) < l1 .or. size(c(1,1,:)) < l1) then
		print*, "SHPowerDensityLC --- Error"
		print*, "C must be dimensioned as (2, L+1, L+1) where L is ", l
		print*, "Input array is dimensioned ", size(c(:,1,1)),  size(c(1,:,1)), size(c(1,1,:)) 
		stop
	endif

	SHPowerDensityLC = dble(c(1, l1, 1))**2 + aimag(c(1, l1, 1))**2

	do m = 1, l, 1
		m1 = m+1
		do i=1, 2
			SHPowerDensityLC = SHPowerDensityLC + dble(c(i, l1, m1))**2 + aimag(c(i, l1, m1))**2
		enddo
	enddo
	
	SHPowerDensityLC = SHPowerDensityLC/dble(2*l+1)

end function SHPowerDensityLC


complex*16 function SHCrossPowerLC(c1, c2, l)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	This function will compute the dimensionless cross power at
!	degree l for the complex spherical harmonic coefficients c1(i, l, m)
!	and c2(i,l,m).
!
!	CrossPower =  Sum_{m=0}^l ( C11lm*C21lm^* + C12lm*C22lm^* )
!
!	Calling Parameters
!		c1	Complex spherical harmonic coefficients, dimensioned as 
!			(2, lmax+1, lmax+1).
!		c2	Complex spherical harmonic coefficients, dimensioned as 
!			(2, lmax+1, lmax+1).
!		l	Spherical harmonic degree to compute power.
!
!	Written by Mark Wieczorek (May 2008).
!
!	Copyright (c) 2008, Mark A. Wieczorek
!	All rights reserved.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	implicit none
	complex*16, intent(in) :: c1(:,:,:), c2(:,:,:)
	integer, intent(in) :: l
	integer i, m, l1, m1	
	
	l1 = l+1
	
	if (size(c1(:,1,1)) < 2 .or. size(c1(1,:,1)) < l1 .or. size(c1(1,1,:)) < l1) then
		print*, "SHCrossPowerLC --- Error"
		print*, "C1 must be dimensioned as (2, L+1, L+1) where L is ", l
		print*, "Input array is dimensioned ", size(c1(:,1,1)),  size(c1(1,:,1)), size(c1(1,1,:)) 
		stop
	elseif (size(c2(:,1,1)) < 2 .or. size(c2(1,:,1)) < l1 .or. size(c2(1,1,:)) < l1) then
		print*, "SHCrossPowerLC --- Error"
		print*, "C2 must be dimensioned as (2, L+1, L+1) where L is ", l
		print*, "Input array is dimensioned ", size(c2(:,1,1)),  size(c2(1,:,1)), size(c2(1,1,:)) 
		stop
	endif

	SHCrossPowerLC = c1(1, l1, 1)*dconjg(c2(1,l1,1))

	do m = 1, l, 1
		m1 = m+1
		do i=1, 2
			SHCrossPowerLC = SHCrossPowerLC + c1(i, l1, m1)*dconjg(c2(i,l1,m1))
		enddo
	enddo

end function SHCrossPowerLC


complex*16 function SHCrossPowerDensityLC(c1, c2, l)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	This function will compute the dimensionless cross power
!	(density) at degree l of the complext spherical harmonic coefficients c1(i, l, m)
!	and c2(i,l,m).
!
!	CrossPower(l) =  Sum_{m=0}^l ( C11lm*C21lm^* + C12lm*C22lm^* ) / (2l+1)
!
!	Calling Parameters
!		c1	Complex spherical harmonic coefficients, dimensioned as 
!			(2, lmax+1, lmax+1).
!		c2	Complex spherical harmonic coefficients, dimensioned as 
!			(2, lmax+1, lmax+1).
!		l	Spherical harmonic degree to compute power.
!
!	Written by Mark Wieczorek (May 2008).
!
!	Copyright (c) 2008, Mark A. Wieczorek
!	All rights reserved.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	implicit none
	complex*16, intent(in) :: c1(:,:,:), c2(:,:,:)
	integer, intent(in) :: l
	integer i, m, l1, m1	
	
	l1 = l+1
	
	if (size(c1(:,1,1)) < 2 .or. size(c1(1,:,1)) < l1 .or. size(c1(1,1,:)) < l1) then
		print*, "SHCrossPowerDensityLC --- Error"
		print*, "C1 must be dimensioned as (2, L+1, L+1) where L is ", l
		print*, "Input array is dimensioned ", size(c1(:,1,1)),  size(c1(1,:,1)), size(c1(1,1,:)) 
		stop
	elseif (size(c2(:,1,1)) < 2 .or. size(c2(1,:,1)) < l1 .or. size(c2(1,1,:)) < l1) then
		print*, "SHCrossPowerDensityLC --- Error"
		print*, "C2 must be dimensioned as (2, L+1, L+1) where L is ", l
		print*, "Input array is dimensioned ", size(c2(:,1,1)),  size(c2(1,:,1)), size(c2(1,1,:)) 
		stop
	endif

	SHCrossPowerDensityLC =  c1(1, l1, 1)*dconjg(c2(1,l1,1))
	
	do m = 1, l, 1
		m1 = m+1
		do i=1, 2
			SHCrossPowerDensityLC = SHCrossPowerDensityLC + c1(i, l1, m1)*dconjg(c2(i,l1,m1))
		enddo
	enddo
	
	SHCrossPowerDensityLC = SHCrossPowerDensityLC/dble(2*l+1)

end function SHCrossPowerDensityLC


subroutine SHPowerSpectrumC(c, lmax, spectra)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	This function will compute the dimensionless power spectrum
!	of the complex spherical harmonic coefficients c(i, l, m).
!
!	Power(l) = Sum_{m=0}^l ( | C1lm |**2 + | C2lm |**2 )
!
!
!	Calling Parameters
!		IN
!			c	Complex pherical harmonic coefficients, dimensioned as 
!				(2, lmax+1, lmax+1).
!			lmax	Spherical harmonic degree to compute power.
!		OUT
!			spectra	Array of length (lmax+1) containing the power
!				spectra of c.
!
!	Written by Mark Wieczorek (May 2008).
!
!	Copyright (c) 2008, Mark A. Wieczorek
!	All rights reserved.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	implicit none
	complex*16, intent(in) :: c(:,:,:)
	integer, intent(in) :: lmax
	real*8, intent(out) ::	spectra(:)
	integer i, m, l1, m1, l
	
	if (size(c(:,1,1)) < 2 .or. size(c(1,:,1)) < lmax+1 .or. size(c(1,1,:)) < lmax+1) then
		print*, "SHPowerSpectrumC --- Error"
		print*, "C must be dimensioned as (2, LMAX+1, LMAX+1) where LMAX is ", lmax
		print*, "Input array is dimensioned ", size(c(:,1,1)),  size(c(1,:,1)), size(c(1,1,:)) 
		stop
	elseif (size(spectra) < lmax+1) then
		print*, "SHPowerSpectrumC --- Error"
		print*, "SPECTRA must be dimensioned as (LMAX+1) where LMAX is ", lmax
		print*, "Input vector has dimension ", size(spectra)
		stop
	endif
	
	spectra = 0.0d0
	
	do l=0, lmax
		l1 = l+1
		
		spectra(l1) = dble(c(1, l1, 1))**2 + aimag(c(1, l1, 1))**2
		
		do m = 1, l, 1
			m1 = m+1
			
			do i=1, 2
				spectra(l1) = spectra(l1) + dble(c(i, l1, m1))**2 + aimag(c(i, l1, m1))**2
			enddo
		enddo
	enddo

end subroutine SHPowerSpectrumC


subroutine SHPowerSpectrumDensityC(c, lmax, spectra)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	This function will compute the dimensionless power spectrum
!	density of the complex spherical harmonic coefficients c(i, l, m)
!
!	PowerSpectralDensityC(l) = Sum_{m=0}^l ( C1lm**2 + C2lm**2 ) / (2l+1)
!
!	Calling Parameters
!		IN
!			c	Complex spherical harmonic coefficients, dimensioned as 
!				(2, lmax+1, lmax+1).
!			lmax	Spherical harmonic degree to compute power.
!		OUT
!			spectra	Array of length (lmax+1) containing the power
!				spectra density of c.
!
!	Written by Mark Wieczorek (May 2008).
!
!	Copyright (c) 2008, Mark A. Wieczorek
!	All rights reserved.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	implicit none
	complex*16, intent(in) :: c(:,:,:)
	integer, intent(in) :: lmax
	real*8, intent(out) ::	spectra(:)
	integer i, m, l1, m1, l
	
	if (size(c(:,1,1)) < 2 .or. size(c(1,:,1)) < lmax+1 .or. size(c(1,1,:)) < lmax+1) then
		print*, "SHPowerSpectrumDensityC --- Error"
		print*, "C must be dimensioned as (2, LMAX+1, LMAX+1) where LMAX is ", lmax
		print*, "Input array is dimensioned ", size(c(:,1,1)),  size(c(1,:,1)), size(c(1,1,:)) 
		stop
	elseif (size(spectra) < lmax+1) then
		print*, "SHPowerSpectrumDensityC --- Error"
		print*, "SPECTRA must be dimensioned as (LMAX+1) where LMAX is ", lmax
		print*, "Input vector has dimension ", size(spectra)
		stop
	endif

	spectra = 0.0d0
	
	do l=0, lmax
		l1 = l+1
		
		spectra(l1) = dble(c(1, l1, 1))**2 + aimag(c(1, l1, 1))**2
		
		do m = 1, l, 1
			m1 = m+1
			
			do i=1, 2
				spectra(l1) = spectra(l1) + dble(c(i, l1, m1))**2 + aimag(c(i, l1, m1))**2
			enddo
		enddo
		
		spectra(l1) = spectra(l1)/dble(2*l+1)
	enddo

end subroutine SHPowerSpectrumDensityC


subroutine SHCrossPowerSpectrumC(c1, c2, lmax, cspectra)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	This function will compute the dimensionless cross power spectrum
!	of the complex spherical harmonic coefficients c1(i, l, m) and c2(1,l,m).
!
!	CrossPower(l) =  Sum_{m=0}^l ( C11lm*C21lm^* + C12lm*C22lm^* )
!
!	Calling Parameters
!		IN
!			c1		Complex spherical harmonic coefficients, dimensioned as 
!					(2, lmax+1, lmax+1).
!			c2		Complex spherical harmonic coefficients, dimensioned as 
!					(2, lmax+1, lmax+1).
!			lmax		Spherical harmonic degree to compute power.
!		OUT
!			cspectra	Array of length (lmax+1) containing the complex cross power
!					spectra of c.
!
!	Written by Mark Wieczorek (May 2008).
!
!	Copyright (c) 2008, Mark A. Wieczorek
!	All rights reserved.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	implicit none
	complex*16, intent(in) :: c1(:,:,:), c2(:,:,:)
	integer, intent(in) :: lmax
	complex*16, intent(out) ::	cspectra(:)
	integer i, m, l1, m1, l
	
	if (size(c1(:,1,1)) < 2 .or. size(c1(1,:,1)) < lmax+1 .or. size(c1(1,1,:)) < lmax+1) then
		print*, "SHCrossPowerSpectrumC --- Error"
		print*, "C1 must be dimensioned as (2, LMAX+1, LMAX+1) where lmax is", lmax
		print*, "Input array is dimensioned ", size(c1(:,1,1)),  size(c1(1,:,1)), size(c1(1,1,:)) 
		stop
	elseif (size(c2(:,1,1)) < 2 .or. size(c2(1,:,1)) < lmax+1 .or. size(c2(1,1,:)) < lmax+1) then
		print*, "SHCrossPowerSpectrumC --- Error"
		print*, "C2 must be dimensioned as (2, LMAX+1, LMAX+1)"
		print*, "Input array is dimensioned ", size(c2(:,1,1)),  size(c2(1,:,1)), size(c2(1,1,:)) 
		stop
	elseif (size(cspectra) < lmax+1) then
		print*, "SHCrossPowerSpectrumC --- Error"
		print*, "CSPECTRA must be dimensioned as (LMAX+1) where lmax is ", lmax
		print*, "Input vector has dimension ", size(cspectra)
		stop
	endif

	cspectra = 0.0d0
	
	do l=0, lmax
		l1 = l+1
		
		cspectra(l1) = c1(1, l1, 1)*dconjg(c2(1, l1, 1))

		do m = 1, l, 1
			m1 = m+1
			
			do i=1, 2
				cspectra(l1) = cspectra(l1) + c1(i, l1, m1)*dconjg(c2(i, l1, m1))
			enddo
		enddo
	enddo

end subroutine SHCrossPowerSpectrumC


subroutine SHCrossPowerSpectrumDensityC(c1, c2, lmax, cspectra)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	This function will compute the dimensionless cross power spectrum
!	density of the complex spherical harmonic coefficients c1(i, l, m) and c2(i,l,m).
!
!	CrossPower(l) =  Sum_{m=0}^l ( C11lm*C21lm^* + C12lm*C22lm^* ) / (2l+1)
!
!	Calling Parameters
!		IN
!			c1		Complex spherical harmonic coefficients, dimensioned as 
!					(2, lmax+1, lmax+1).
!			c2		Complex spherical harmonic coefficients, dimensioned as 
!					(2, lmax+1, lmax+1).
!			lmax		Spherical harmonic degree to compute power.
!		OUT
!			cspectra	Array of length (lmax+1) containing the complex cross power
!					spectral density of c.
!
!	Written by Mark Wieczorek (May 2008).
!
!	Copyright (c) 2008, Mark A. Wieczorek
!	All rights reserved.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	implicit none
	complex*16, intent(in) :: c1(:,:,:), c2(:,:,:)
	integer, intent(in) :: lmax
	complex*16, intent(out) ::	cspectra(:)
	integer i, m, l1, m1, l
	
	if (size(c1(:,1,1)) < 2 .or. size(c1(1,:,1)) < lmax+1 .or. size(c1(1,1,:)) < lmax+1) then
		print*, "SHCrossPowerSpectrumDensityC --- Error"
		print*, "C1 must be dimensioned as (2, LMAX+1, LMAX+1) where lmax is", lmax
		print*, "Input array is dimensioned ", size(c1(:,1,1)),  size(c1(1,:,1)), size(c1(1,1,:)) 
		stop
	elseif (size(c2(:,1,1)) < 2 .or. size(c2(1,:,1)) < lmax+1 .or. size(c2(1,1,:)) < lmax+1) then
		print*, "SHCrossPowerSpectrumDensityC --- Error"
		print*, "C2 must be dimensioned as (2, LMAX+1, LMAX+1)"
		print*, "Input array is dimensioned ", size(c2(:,1,1)),  size(c2(1,:,1)), size(c2(1,1,:)) 
		stop
	elseif (size(cspectra) < lmax+1) then
		print*, "SHCrossPowerSpectrumDensityC --- Error"
		print*, "CSPECTRA must be dimensioned as (LMAX+1) where lmax is ", lmax
		print*, "Input vector has dimension ", size(cspectra)
		stop
	endif

	cspectra = 0.0d0
	
	do l=0, lmax
		l1 = l+1
		
		cspectra(l1) = c1(1, l1, 1)*dconjg(c2(1, l1, 1))

		do m = 1, l, 1
			m1 = m+1
			
			do i=1, 2
				cspectra(l1) = cspectra(l1) + c1(i, l1, m1)*dconjg(c2(i, l1, m1))
			enddo
			
		enddo
		
		cspectra(l1) = cspectra(l1)/dble(2*l+1)
		
	enddo

end subroutine SHCrossPowerSpectrumDensityC

