subroutine ComputeD0(D0, lmax, theta0)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	This routine will compute the kernel D0 for the isotropic space
!	concentration problem of a spherical cap. All terms are computed
!	using exact expressions. The eigenfunctions of this kernel are
!	spherical harmonic coefficients normalized according to the "geodesy"
!	convention. The diagonal elements of D0 approach unity as theta
!	approaches 180 degrees.
!
!	Calling Parameters
!		IN
!			lmax:		Maximum spherical harmonic degree.
!			theta0:		Angular radius of spherical cap IN RADIANS.
!		OUT
!			D0:		Symmetric kernel of size lmax+1 by lmax+1.
!
!	Dependencies:	PLegendreL_d1
!
!	Written by Mark Wieczorek (June 2004)
!
!	Copyright (c) 2005, Mark A. Wieczorek
!	All rights reserved.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	use SHTOOLS, only: PLegendre_d1

	implicit none

	real*8, intent(out) ::	D0(:,:)
	real*8, intent(in) ::	theta0
	integer, intent(in) :: 	lmax
	real*8 ::		x, wll, p(2*lmax+1), dp(2*lmax+1)
	integer ::		l, lp, j

	if (size(D0(1,:)) < lmax+1 .or. size(D0(:,1)) < lmax + 1) then
		print*, "Error --- ComputeD0"
		print*, "D0 must be dimensioned as (LMAX+1,LMAX+1) where LMAX is ", lmax
		print*, "Input array is dimensioned as ", size(D0(1,:)), size(D0(:,1))
		stop
	endif

	D0 = 0.0d0

	x = cos(theta0)
	
	call PLegendre_d1(p, dp, 2*lmax, x)
					
	! Compute diagonal terms using an exact sum involving Wigner 3j symbols
	
	do l=0, lmax
	
		wll = (-1.0d0)**l / sqrt(2.0d0*l+1.0d0)		! Wig3j0(l,l,0)
		D0(l+1,l+1) = wll**2 * (1.0d0 - x)		! Fist term inside sum for l_prime = 0
		
		do j = 2, 2*l, 2
			wll = wll * dble(1-j) * sqrt( dble(2*l+j)*dble(2*l-j+2) )  &	! Wig3j0(l,l,j)
				 / dble(j) / sqrt( dble(2*l+j+1)*dble(2*l-j+1) ) 
			D0(l+1,l+1) = D0(l+1,l+1) +  wll**2 * dble(2*j+1)  &
				* (p(j) - x * p(j+1) ) / dble(j+1)
		enddo
		
		D0(l+1,l+1) = D0(l+1,l+1)*(2.0d0*l+1.0d0)
			
	enddo
	
	! compute off diagonal terms using exact expressions
	
	do l=0, lmax
	
		do lp = l+1, lmax, 1
		
			D0(l+1, lp+1) = sqrt( dble(2*l+1)*dble(2*lp+1) ) * (1.0d0-x**2) * &
				(p(lp+1) * dp(l+1) - p(l+1) * dp(lp+1) ) / &
				(dble(l**2+l) - dble(lp**2+lp))	
			D0(lp+1,l+1) = D0(l+1,lp+1)
			
		enddo
	enddo
	
	D0 = D0 / 2.0d0
	
end subroutine ComputeD0

