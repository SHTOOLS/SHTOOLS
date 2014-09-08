subroutine SphericalCapCoef(coef, theta, lmax)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	This routine will return the coefficients for a spherical
!	cap (located at the north pole) with angular radius theta,
!	using the geodesy 4-pi normalization.
!
!	Calling Parameters:
!		IN
!			theta	Angular radius of spherical cap in RADIANS.
!		OUT
!			coef	m=0 coefficients normalized according to the 
!				geodesy convention, further normalized such
!				that the degree-0 term is 1. (i.e., function 
!				has an average value of one over the entire
!				sphere).
!		OPTIONAL
!			lmax	Maximum spherical harmonic degree.
!
!	Dependencies: PlBar
!
!	Written by Mark Wieczorek
!
!	Copyright (c) 2005, Mark A. Wieczorek
!	All rights reserved.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	use SHTOOLS, only: PlBar

	implicit none
	real*8, intent(out) ::	coef(:)
	real*8, intent(in) ::	theta
	integer, intent(in), optional ::	lmax
	real*8 ::	x, top, bot, pi
	real*8, allocatable ::	pl(:)
	integer ::	l, lmax2, astat
	
	
	pi = acos(-1.0d0)
	coef = 0.0d0
	x = cos(theta)
	
	if (present(lmax) ) then
		lmax2 = lmax
	else
		lmax2 = size(coef) - 1
	endif
	
	allocate(pl(lmax2+3), stat = astat)
	
	if(astat /=0) then
		print*, "Error --- SphericalCapCoef"
		print*, "Unable to allocate array pl", astat
		stop
	endif
	
	call PlBar(pl, lmax2+2, x)
	
	coef(1) = 1.0d0
	
	bot = pl(1) - pl(2) / sqrt(3.0d0)
	
	do l=1, lmax2
		top =  pl(l) / sqrt(dble(2*l-1)) -pl(l+2) / sqrt(dble(2*l+3))
		coef(l+1) = top / (bot * sqrt(dble(2*l+1)) )
	enddo
	
	deallocate (pl)
		
end subroutine SphericalCapCoef

