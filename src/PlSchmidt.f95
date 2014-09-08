subroutine PlSchmidt(p,lmax,z)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	This function evalutates all of the Schmidt normalized legendre 
!	polynomials up to degree lmax. 
!
!	Calling Parameters:
!		Out
!			p:	A vector of all Schmidt normalized Legendgre polynomials evaluated at 
!				z up to lmax. The lenght must by greater or equal to (lmax+1).
!		IN
!			lmax:	Maximum degree to compute.
!			z:	[-1, 1], cos(colatitude) or sin(latitude).
!
!	Notes:
!	
!	1.	The integral of plm**2 over (-1,1) is 2 * / (2l+1).
!	2.	The integral of Plm**2 over all space is 4 pi / (2l+1).
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
       	real*8 ::	pm2, pm1, pl
      	integer ::	l


	if (size(p) < lmax+1) then
		print*, "Error --- PlSchmidt"
     		print*, "P must be dimensioned as (LMAX+1) where LMAX is ", lmax 
     		print*, "Input array is dimensioned ", size(p)
     		stop
     	elseif (lmax < 0) then 
     		print*, "Error --- PlSchmidt"
     		print*, "LMAX must be greater than or equal to 0."
     		print*, "Input value is ", lmax
     		stop
     	elseif(abs(z) > 1.0d0) then
     		print*, "Error --- PlSchmidt"
     		print*, "ABS(Z) must be less than or equal to 1."
     		print*, "Input value is ", z
     		stop
     	endif
     	
   	pm2  = 1.d0
      	p(1) = 1.d0
      	
      	pm1  = z
      	p(2) = pm1
      	
      	do l = 2, lmax
         	pl = ( dble(2*l-1)  * z * pm1 - dble(l-1) * pm2 )  / dble(l)
         	p(l+1) = pl
         	pm2  = pm1
         	pm1  = pl
      	enddo

end subroutine PlSchmidt

