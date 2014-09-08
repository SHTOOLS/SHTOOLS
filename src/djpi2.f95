subroutine djpi2(dj, lmax)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	This subroutine computes 
!	  j
!	d     (pi/2)
!	  m N
!
!	for all posible values of m (>=0) and N for 0<l<lmax.
!	The output array corresponds to dj(N,m,j). (I'm not positive
!	about the N and m ordering).
!
!	Note that this algorithm is easy to modify to work on a
!	single l (just need a loop to compute f1) since there are
!	no l recursions involved.
!
!	Calling Parameters
!		IN
!			lmax	Maximum spherical harmonic degree to be computed.
!		OUT
!			dj	Rotation matrix with dimensions (lmax+1, lmax+1, lmax+1).
!
!	Dependencies:	None
!
!	History
!
!		1. Based on routine from Guy Masters (July 16, 1993)
!		2. Modified by Mark Simons (July 19, 1993)
!		3. Turned into readable f95 code by Mark Wieczorek (August, 2003).
!
!	Copyright (c) 2005, Mark A. Wieczorek
!	All rights reserved.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	implicit none
	integer, intent(in) ::	lmax
	real*8, intent(out) ::	dj(:,:,:)
	integer ::		i, l, n, lp1, j, m, isn, np
	real*8 ::		f((lmax+1)*8), f1, f2, g1, g2, en2, fl2p1
	
	if (size(dj(:,1,1)) < lmax+1 .or. size(dj(1,:,1)) < lmax + 1 .or. size(dj(1,1,:)) < lmax +1) then
		print*, "Error --- djpi2"
		print*, "DJ must be dimensioned (LMAX+1, LMAX+1, LMAX+1) where LMAX is ", lmax
		print*, "Input array is dimensioned as ", size(dj(:,1,1)), size(dj(1,:,1)), size(dj(1,1,:))
		stop
	endif
	
	dj = 0.0d0
	
	dj(1,1,1) = 1.d0
	dj(1,1,2) = 0.d0
	dj(1,2,2) = -1.0d0/sqrt(2.0d0)
	dj(2,1,2) = -dj(1,2,2)
	dj(2,2,2) = 0.5d0
	f1        = 0.5d0

	do l = 2, lmax

		lp1   = l+1
		fl2p1 = l+lp1
		
		do i = 1, l
			f(i) = dsqrt(i*(fl2p1-i))
        	enddo
        	
        	f1 = f1*(l+l-1.0d0)/(l+l)

		! Do N = 0 terms
		
		dj(lp1,1,lp1) = -dsqrt(f1)
		dj(l,1,lp1)   = 0.0d0
		
		do i = 2, l
			j = lp1-i
			dj(j,1,lp1) = -f(i-1)*dj(j+2,1,lp1)/f(i)
		enddo

		! Do positive N terms (bottom triangle)
		
		f2 = f1
		g1 = l
		g2 = lp1
		
		do n = 1, l
			np = n+1
			en2 = n+n
			g1 = g1 + 1.0d0
			g2 = g2 - 1.0d0
			f2 = f2*g2/g1
			dj(lp1,np,lp1) = -dsqrt(f2)
			dj(l,np,lp1) = dj(lp1,np,lp1)*en2/f(1)

			do i = 2, l-n
				j = lp1-i
				dj(j,np,lp1) = ( en2*dj(j+1,np,lp1) &
					- f(i-1)*dj(j+2,np,lp1) ) / f(i)
			enddo
			
		enddo

		! Fill upper triangle and fix signs

		do j = 1, l
			do m = j,l
				dj(j,m+1,lp1) = dj(m+1,j,lp1)
			enddo
		enddo

		isn = 1 + mod(l,2)
		
		do np = 1, lp1
			do i = isn, lp1, 2
				dj(i,np,lp1) = -dj(i,np,lp1)
			enddo
		enddo
		
	enddo
	
end subroutine djpi2

