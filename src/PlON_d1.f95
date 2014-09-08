subroutine PlON_d1(p, dp, lmax, z)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	This function evalutates all of the "ortho-normalized legendre 
!	polynomials, as well as their first derivatives,  up to degree lmax. 
!
!	Calling Parameters:
!		Out
!			p:	A vector of all normalized Legendgre polynomials evaluated at 
!				z up to lmax with dimension (lmax+1).
!			dp:	A vector of all first derivatives of the normalized Legendgre polynomials evaluated at 
!				z up to lmax with dimension (lmax+1).
!		IN
!			lmax:	Maximum degree to compute.
!			z:	[-1, 1], cos(colatitude), or sin(latitude).
!
!	Notes:
!	
!	1.	The employed normalization is the "orthonomalization convention."
!	2.	The integral of PlBar**2 over all space on the sphere is 1. 
!	3.	The integral of PlBar**2 over (-1,1) is 1/2pi.	
!	4. 	The derivative is evaluated with respecte to z, and NOT cos(colatitude) or sin(latitude).
!	5.	Derivatives are calculated according to the geodesy-normalized relationships
!			P'_0(z) = 0.0, P'_1(z) = 1.0, and
!			P'_l(z) = l * (P_{l-1}(z) - z * P_l(z) ) / (1.0d0 - z**2)
!			At z = 1, Pl(1) = 1, and P'l(1) = l (l+1) / 2 	(Boyd 2001)
!			At z = -1 Pl(-1) = (-1)**l, and P'l(-1) = (-1)**(l-1) l (l+1) / 2
!
!	Dependencies:	None
!
!	Written by Mark Wieczorek June 2004
!
!	Copyright (c) 2005-2006, Mark A. Wieczorek
!	All rights reserved.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
	implicit none
	integer, intent(in) ::	lmax
	real*8, intent(out) ::	p(:), dp(:)
       	real*8, intent(in) ::	z
       	real*8 ::	pm2, pm1, pl, sinsq, pi
      	integer ::	l

	pi = acos(-1.0d0)

	if (size(p) < lmax+1) then
		print*, "Error --- PlBar_d1"
     		print*, "P must be dimensioned as (LMAX+1) where LMAX is ", lmax 
     		print*, "Input array is dimensioned ", size(p)
     		stop
     	elseif (size(dp) < lmax+1) then
		print*, "Error --- PlBar_d1"
     		print*, "DP must be dimensioned as (LMAX+1) where LMAX is ", lmax 
     		print*, "Input array is dimensioned ", size(dp)
     		stop
     	elseif (lmax < 0) then 
     		print*, "Error --- PlBar_d1"
     		print*, "LMAX must be greater than or equal to 0."
     		print*, "Input value is ", lmax
     		stop
     	elseif(abs(z) > 1.0d0) then
     		print*, "Error --- PlBar_d1"
     		print*, "ABS(Z) must be less than or equal to 1."
     		print*, "Input value is ", z
     		stop
     	endif
     	      	
      	if (z == 1.0d0) then
      	
      		do l=0, lmax
      			p(1:lmax+1) = sqrt( dble(2*l+1))
      			dp(l+1) = sqrt( dble(2*l+1)) * dble(l)*dble(l+1)/2.0d0
      		enddo
      		
      	elseif (z == -1.0d0) then
      		
      		do l=0, lmax
      			p(l+1) = sqrt( dble(2*l+1) ) * dble((-1)**l)
      			dp(l+1) =  sqrt( dble(2*l+1)) *  dble(l)*dble(l+1) * dble((-1)**(l-1)) / 2.0d0
      		enddo
      	
      	else

		sinsq = (1.0d0 - z**2)
      	
   		pm2  = 1.d0/sqrt(4*pi)
      		p(1) = pm2
      		dp(1) = 0.0d0
      	
      		pm1  = sqrt(3.0d0) * z / sqrt(4*pi)
      		p(2) = pm1
      		dp(2) = sqrt(3.0d0) / sqrt(4*pi)
      	      	
      		do l = 2, lmax
         		pl = ( sqrt(dble(2*l-1))  * z * pm1 - &
         			(l-1) * pm2 / sqrt(dble(2*l-3)) ) * &
         			sqrt(dble(2*l+1))/dble(l)
         		p(l+1) = pl
         		dp(l+1) = l *( sqrt( dble(2*l+1)/dble(2*l-1) ) * &
      				p(l) - z * pl ) / sinsq
         		pm2  = pm1
         		pm1  = pl
      		enddo     
      		
      	endif 	
      
end subroutine PlON_d1

