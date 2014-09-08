subroutine PlSchmidt_d1(p, dp, lmax, z)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	This function evalutates all of the Schmidt normalized legendre 
!	polynomials and their first derivatives up to degree lmax. 
!
!	Calling Parameters:
!		Out
!			p:	A vector of all Schmidt normalized Legendgre polynomials evaluated at 
!				z up to lmax. The lenght must by greater or equal to (lmax+1).
!			dp:	A vector of all associated Legendgre polynomials evaluated at 
!				z up to lmax. The lenght must by greater or equal to (lmax+1).
!		IN
!			lmax:	Maximum degree to compute.
!			z:	[-1, 1], cos(colatitude) or sin(latitude).
!
!	Notes:
!	
!	1.	The integral of plm**2 over (-1,1) is 2 * / (2l+1).
!	2.	The integral of Plm**2 over all space is 4 pi / (2l+1).
!	3.	Derivatives are calculated according to the unnormalized relationships:
!			P'_0(z) = 0.0, P'_1(z) = 1.0, and
!			P'_l(z) = l * (P'_{l-1}(z) - z * P_l(z) ) / (1.0d0 - z**2)
!			At z = 1, Pl(1) = 1, and P'l(1) = l (l+1) / 2 	(Boyd 2001)
!			At z = -1 Pl(-1) = (-1)**l, and P'l(-1) = (-1)**(l-1) l (l+1) / 2
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
	real*8, intent(out) ::	p(:), dp(:)
       	real*8, intent(in) ::	z
       	real*8 ::	pm2, pm1, pl, sinsq
      	integer ::	l

	if (size(p) < lmax+1) then
		print*, "Error --- PlSchmidt_d1"
     		print*, "P must be dimensioned as (LMAX+1) where LMAX is ", lmax 
     		print*, "Input array is dimensioned ", size(p)
     		stop
     	elseif (size(dp) < lmax+1) then
		print*, "Error --- PlSchmidt_d1"
     		print*, "DP must be dimensioned as (LMAX+1) where LMAX is ", lmax 
     		print*, "Input array is dimensioned ", size(dp)
     		stop
     	elseif (lmax < 0) then 
     		print*, "Error --- PlSchmidt_d1"
     		print*, "LMAX must be greater than or equal to 0."
     		print*, "Input value is ", lmax
     		stop
     	elseif(abs(z) > 1.0d0) then
     		print*, "Error --- PlSchmidt_d1"
     		print*, "ABS(Z) must be less than or equal to 1."
     		print*, "Input value is ", z
     		stop
     	endif
     	
	if (z == 1.0d0) then
      	
      		p(1:lmax+1) = 1.0d0
      		do l=0, lmax
      			dp(l+1) = dble(l)*dble(l+1)/2.0d0
      		enddo
      		
      	elseif (z == -1.0d0) then
      		
      		do l=0, lmax
      			p(l+1) = dble((-1)**l)
      			dp(l+1) = dble(l)*dble(l+1) * dble((-1)**(l-1)) / 2.0d0
      		enddo
      	
      	else	

		sinsq = (1.0d0-z)*(1.0d0+z)

   		pm2  = 1.0d0
      		p(1) = 1.0d0
      		dp(1) = 0.0d0
      	
      		pm1  = z
      		p(2) = pm1
      		dp(2) = 1.0d0
      	
      		do l = 2, lmax
         		pl = ( (2*l-1) * z * pm1 - (l-1) * pm2 ) / dble(l)
         		p(l+1) = pl
         		dp(l+1) =  dble(l) * (pm1 - z * pl) / sinsq
         		pm2  = pm1
         		pm1  = pl
      		enddo
      		
      	endif

end subroutine PlSchmidt_d1

