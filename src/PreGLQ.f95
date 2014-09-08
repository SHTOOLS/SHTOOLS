subroutine PreGLQ(x1, x2, n, zero, w)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	This routine will find the zeros and weights that are
!	used in Gauss-Legendre quadrature routines. (Based on routines
!	in Numerical Recipes).
!
!	Calling Parameters:
!		IN
!			x1: 	Lower bound of integration.
!			x2:	Upper bound of integration.
!			n:	Number of points used in the quadrature. n points
!				will integrate a polynomial of degree 2n-1 exactly.
!		OUT
!			zero:	Array of n Gauss points, which correspond to the zeros
!				of P(n,0).
!			w:	Array of n weights used in the quadrature.
!
!
!	Note 
!		1.	If EPS is less than what is defined, then the do 
!			loop for finding the roots might never terminate for some
!			values of lmax. If the algorithm doesn't converge, consider
!			increasing itermax, or decreasing eps.
!
!	Dependencies:	None
!
!	Written by Mark Wieczorek 2003
!
!	Copyright (c) 2005-2006, Mark A. Wieczorek
!	All rights reserved.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	implicit none
	real*8, intent(in) :: 	x1, x2
	real*8, intent(out) ::	zero(:), w(:)
	integer, intent(in) ::	n
	integer ::		i, j, m, iter
	integer, parameter ::	itermax = 1000
	real*8, parameter ::	eps=1.0d-15
	real*8 ::		p1, p2, p3, pp, z, z1, xm, xu
	
	
	if (size(zero) < n) then
		print*, "Error --- PreGLQ"
		print*, "ZERO must be dimensioned as (N) where N is ", n
		print*, "Input array is dimensioned ", size(zero)
		stop
	elseif (size(w) < n) then
		print*, "Error --- PreGLQ"
		print*, "W must be dimensioned as (N) where N is ", n
		print*, "Input array is dimensioned ", size(w)
		stop
	endif
	
	
	zero = 0.0d0
	w = 0.0d0
	
	m=(n+1)/2 		! The roots are symmetric in the interval, so we
				! only have to find half of them. 
	xm = (x2 + x1) / 2.0d0	! midpoint of integration
	xu = (x2 - x1) / 2.0d0	! Scaling factor between interval of integration, and that of 
				! the -1 to 1 interval for the Gauss-Legendre interval
	
	! Compute roots and weights
				
	do i=1,m 
		iter = 0
		z=cos(3.141592654d0*(i-.25d0)/(n+.5d0))	! Approximation for the ith root
							
		! Find the true value using newtons method
		
		do
		
			iter = iter +1
			
			p1=1.0d0
			p2=0.0d0
			do j=1, n 	! determine the legendre polynomial evaluated at z (p1) using 
					! recurrence relationships
				p3=p2
				p2=p1
				p1=(dble(2*j-1)*z*p2-dble(j-1)*p3)/dble(j)
			enddo
		
			pp=dble(n)*(z*p1-p2)/(z*z-1.0d0)	! This is the derivative of the legendre polynomial
								! using recurrence relationships
		
			z1=z					! This is newtons method here
			z=z1-p1/pp 
		
			if (abs(z-z1) <= eps) exit
			
			if (iter >itermax) then
				print*, "ERROR --- PreGLQ"
				print*, "Root Finding of PreGLQ not converging."
				print*, "m , n = ", m, n
				stop
			endif
			
		enddo
		
		zero(i) = xm + xu*z
		zero(n+1-i) = xm - xu*z
		w(i)= 2.0d0 * xu / ((1.0d0-z*z)*pp*pp)
		w(n+1-i)=w(i) 
		
	enddo
	
end subroutine PreGLQ


integer function NGLQ(degree)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	For a polynomial of order degree, this simple function
!	will determine how many gauss-legendre quadrature points
!	are needed in order to integrate the function exactly.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	implicit none
	integer, intent(in) ::	degree
	
	nglq = ceiling((degree+1.0d0)/2.0d0) 	
	
end function NGLQ


integer function NGLQSH(degree)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	This function returns the number of gauss-legendre points that
!	are needed to exactly integrate a spherical harmonic field of 
!	Lmax = degree
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	implicit none
	integer, intent(in) ::	degree
	
	nglqsh = degree + 1.0d0

end function NGLQSH


integer function NGLQSHN(degree, n)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	This function returns the number of gauss-legendre points that
!	are needed to exactly integrate a spherical harmonic field of 
!	Lmax = degree raised to the nth power. Here, the maximum degree
!	of the integrand is n*lmax + lmax, or (n+1)*lmax
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	implicit none
	integer, intent(in) ::	degree, n
	
	nglqshn = ceiling( ((n+1.0d0)*degree + 1.0d0)/2.0d0)

end function NGLQSHN

