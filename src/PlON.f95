subroutine PlON(p,lmax,z)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	This function evalutates all of the "ortho-normalized legendre 
!	polynomials up to degree lmax. 
!
!	Calling Parameters:
!		Out
!			p:	A vector of all associated Legendgre polynomials evaluated at 
!				z up to lmax. The lenght must by greater or equal to (lmax+1).
!		IN
!			lmax:	Maximum degree to compute.
!			z:	cos(colatitude) or sin(latitude).
!
!	Notes:
!	
!	1.	The employed normalization is the "orthonormalization convention." 
!	2.	The integral of PlBar**2 over all space on the sphere is 1. 
!		The integral of PlBar**2 over (-1,1) is 1/2pi.
!
!	Dependencies:	None
!
!	Written by Mark Wieczorek June 2004
!
!	Copyright (c) 2005, Mark A. Wieczorek
!	All rights reserved.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
	implicit none
	integer, intent(in) ::	lmax
	real*8, intent(out) ::	p(:)
       	real*8, intent(in) ::	z
       	real*8 ::	pm2, pm1, pl, pi
      	integer ::	l


	if (size(p) < lmax+1) then
		print*, "Error --- PlBar"
     		print*, "P must be dimensioned as (LMAX+1) where LMAX is ", lmax 
     		print*, "Input array is dimensioned ", size(p)
     		stop
     	elseif (lmax < 0) then 
     		print*, "Error --- PlBar"
     		print*, "LMAX must be greater than or equal to 0."
     		print*, "Input value is ", lmax
     		stop
     	elseif(abs(z) > 1.0d0) then
     		print*, "Error --- PlBar"
     		print*, "ABS(Z) must be less than or equal to 1."
     		print*, "Input value is ", z
     		stop
     	endif
      	
      	pi = acos(-1.0d0)
      	
   	pm2  = 1.d0/sqrt(4*pi)
      	p(1) = pm2
      	
      	pm1  = sqrt(3.0d0) * z / sqrt(4*pi)
      	p(2) = pm1
      	
      	do l = 2, lmax
         	pl = ( sqrt(dble(2*l-1)) * z * pm1 - &
         		(l-1) * pm2 / sqrt(dble(2*l-3)) ) * &
         		sqrt(dble(2*l+1))/dble(l)
         	p(l+1) = pl
         	pm2  = pm1
         	pm1  = pl
      	enddo

end subroutine PlON

