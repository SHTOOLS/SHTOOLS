subroutine DHaj(n, aj)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	This subroutine will compute the weights a_j that are used to
!	expand an equally sampled grid into spherical harmonics by using
!	the sampling theorem presented in Driscoll and Healy (1994).
!
!	Note that the number of samples, n, must be even! Also, a_j(1) = 0
!
!	Calling parameters:
!		IN
!			n	Number of samples in longitude and latitude.
!		OUT
!			aj	Vector of length n containing the weights.
!
!	Dependencies:	None
!
!	Written by Mark Wieczorek (May, 2004)
!
!	Copyright (c) 2005, 2006, Mark A. Wieczorek
!	All rights reserved.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	implicit none
	integer, intent(in) ::	n
	real*8, intent(out) ::	aj(:)
	integer ::		j, l
	real*8 ::		sum1, pi
	
	pi = acos(-1.0d0)
	
	aj = 0.0d0
	
	if (mod(n,2) /= 0) then
		print*, "Error --- DH_aj"
		print*, "The number of samples in the equi-dimensional grid must be even for use with SHExpandDH"
		print*, "Input value of N is ", n
		stop
	elseif (size(aj) < n) then
		print*, "Error --- DH_aj"
		print*, "The size of AJ must be greater than or equal to N where N is ", n
		print*, "Input array is dimensioned as ", size(aj)
		stop
	endif
	
	do j=0, n-1
		sum1 = 0.0d0
		
		do l=0, n/2 -1
			sum1 = sum1 + sin( dble(2*l+1) * pi * dble(j) / dble(n) ) / dble(2*l+1)
		enddo
	
		aj(j+1) = sum1  * sin(pi * dble(j)/dble(n)) * sqrt(8.0d0) / dble(n)
	
	enddo
	
end subroutine DHaj

