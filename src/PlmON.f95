subroutine PlmON(p, lmax, z, csphase, cnorm)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	This function evalutates all of the normalized associated Legendre
!	functions up to degree lmax. The functions are initially scaled by 
!	10^280 sin^m in order to minimize the effects of underflow at large m 
!	near the poles (see Holmes and Featherstone 2002, J. Geodesy, 76, 279-299). 
!	On a Mac OSX system with a maximum allowable double precision value of 
!	2.225073858507203E-308 the scaled portion of the algorithm will not overflow 
!	for degrees less than or equal to 2800.
!
!	For each value of m, the rescaling factor is computed as rescalem=rescalem*sin(theta), 
!	with the intial value of rescalem being equal to 1/scalef (which is here equal 
!	to 10^280). This will gradually reduce this huge number to a tiny number, and will 
!	ultimately underflow. In order to prevent this underflow, when rescalem becomes less than
!	10^(-280), the subsequent rescaling factors of sin(theta) will be directly applied to Plm, and then this
!	number will be multipled by the old value of rescalem.
!
!	Temporary variables in saved in an allocated array. In order to explicitly deallocate this
!	memory, call this routine with a spherical harmonic degree of -1.
!
!	Calling Parameters:
!		OUT
!			p:		A vector of all associated Legendgre polynomials evaluated at 
!					z up to lmax. The length must by greater or equal to (lmax+1)*(lmax+2)/2.
!		OPTIONAL (IN)
!			csphase:	1: Do not include the phase factor of (-1)^m
!					-1: Apply the phase factor of (-1)^m.
!			cnorm:		0: Use real normalization.
!					1: Use complex normalization.
!		IN
!			lmax:		Maximum spherical harmonic degree to compute.
!			z:		cos(colatitude) or sin(latitude).
!
!	Notes:
!	
!	1.	The employed normalization is the "orthonormalized convention." The integral of
!		(plm*cos(m theta))**2 or (plm*sin (m theta))**2 over all space is 1.
!	2.	The integral of plm**2 over (-1,1) is (2 - delta(0,m))/2pi. If CNORM=1, then
!		this is equal to 1/2pi.
!	3.	The index of the array p corresponds to l*(l+1)/2 + m + 1. As such
!		the array p should be dimensioned as (lmax+1)*(lmax+2)/2 in the 
!		calling routine.
!	4. 	The default is to exclude the Condon-Shortley phase of (-1)^m.
!
!
!	Dependencies:	CSPHASE_DEFAULT
!
!	Written by Mark Wieczorek September 25, 2005.
!
!	April 19, 2008: Added CNORM optional parameter compute complex normalized functions.
!	August 14, 2012: Modified to save auxilliary variables whenever lmax <= lmax_old (instead of
!			reinitializing whenever lmax /= lmax_old). Converted PHASE from REAL*8 to 
!			INTEGER*1.
!
!	Copyright (c) 2008, Mark A. Wieczorek
!	All rights reserved.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	

	use SHTOOLS, only: CSPHASE_DEFAULT

	implicit none
	integer, intent(in) ::	lmax
	real*8, intent(out) ::	p(:)
       	real*8, intent(in) ::	z
       	integer, intent(in), optional :: csphase, cnorm
       	real*8 ::	pm2, pm1, pmm, plm, rescalem, pi, u, scalef
      	real*8, save, allocatable ::	f1(:), f2(:), sqr(:)
      	integer ::	k, kstart, m, l, astat(3)
      	integer, save ::	lmax_old = 0
      	integer*1 ::	phase

	if (lmax == -1) then
		if (allocated(sqr)) deallocate(sqr)
		if (allocated(f1)) deallocate(f1)
		if (allocated(f2)) deallocate(f2)
		lmax_old = 0
		return
	endif

	pi = acos(-1.0d0)
	
	if (size(p) < (lmax+1)*(lmax+2)/2) then 
		print*, "Error --- PlmBar"
     		print*, "P must be dimensioned as (LMAX+1)*(LMAX+2)/2 where LMAX is ", lmax
     		print*, "Input array is dimensioned ", size(p)
     		stop
     	elseif (lmax < 0) then 
     		print*, "Error --- PlmBar"
     		print*, "LMAX must be greater than or equal to 0."
     		print*, "Input value is ", lmax
     		stop
     	elseif(abs(z) > 1.0d0) then
     		print*, "Error --- PlmBar"
     		print*, "ABS(Z) must be less than or equal to 1."
     		print*, "Input value is ", z
     		stop
     	endif     	
     	
     	if (present(csphase)) then
     		if (csphase == -1) then
     			phase = -1
     		elseif (csphase == 1) then
     			phase = 1
     		else
     			print*, "PlmBar --- Error"
     			print*, "CSPHASE must be 1 (exclude) or -1 (include)."
     			print*, "Input value is ", csphase
     			stop
     		endif
     	else
     		phase = CSPHASE_DEFAULT
     	endif
     		
	scalef = 1.0d-280
	
	
	if (lmax > lmax_old) then
		
		if (allocated(sqr)) deallocate(sqr)
		if (allocated(f1)) deallocate(f1)
		if (allocated(f2)) deallocate(f2)
		
		allocate(sqr(2*lmax+1), stat=astat(1))
		allocate(f1((lmax+1)*(lmax+2)/2), stat=astat(2))
		allocate(f2((lmax+1)*(lmax+2)/2), stat=astat(3))
		
		if (astat(1) /= 0 .or. astat(2) /= 0 .or. astat(3) /= 0) then
			print*, "PlmON --- Error"
			print*, "Problem allocating arrays SQR, F1 and F2", astat(1), astat(2), astat(3)
			stop
		endif
	
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!
		!	Precompute square roots of integers that are used several times.
		!
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
		do l=1, 2*lmax+1
			sqr(l) = sqrt(dble(l))
		enddo

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!
		! 	Precompute multiplicative factors used in recursion relationships
		! 		Plmbar(l,m) = x*f1(l,m)*Plmbar(l-1,m) - Plmbar(l-2,m)*f2(l,m)
		!		k = l*(l+1)/2 + m + 1
		!	Note that prefactors are not used for the case when m=l and m=l-1,
		!	as a different recursion is used for these two values.
		!
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
		k = 3
	
		do l=2, lmax, 1
			k = k + 1
			f1(k) = sqr(2*l-1) * sqr(2*l+1) / dble(l)
			f2(k) = dble(l-1) * sqr(2*l+1) / sqr(2*l-3) / dble(l)
			do m=1, l-2
				k = k+1
				f1(k) = sqr(2*l+1) * sqr(2*l-1) / sqr(l+m) / sqr(l-m)
                		f2(k) = sqr(2*l+1) * sqr(l-m-1) * sqr(l+m-1) &
                  			 / sqr(2*l-3) / sqr(l+m) / sqr(l-m) 
			enddo
			k = k + 2
		enddo
		
		lmax_old = lmax
		
	endif
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!	
	!	Calculate P(l,0). These are not scaled.
	!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	u = sqrt((1.0d0-z)*(1.0d0+z)) ! sin(theta)

      	pm2  = 1.0d0/sqrt(4.0d0*pi)
      	p(1) = pm2
      	
      	if (lmax == 0) return
      	
      	pm1  = sqr(3)*z/sqrt(4*pi)
      	p(2) = pm1
      		
	k = 2

      	do l = 2, lmax, 1
         	k = k+l
         	plm = f1(k)*z*pm1-f2(k)*pm2
         	p(k) = plm
         	pm2  = pm1
         	pm1  = plm
      	enddo

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!
	!	Calculate P(m,m), P(m+1,m), and P(l,m)
	!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
	
 	if (present(cnorm)) then
		if (cnorm == 1) then
			pmm  = scalef/sqrt(4*pi)
		else
			pmm  = sqr(2)*scalef/sqrt(4*pi)
		endif
	else
      		pmm  = sqr(2)*scalef/sqrt(4*pi)
      	endif

      	rescalem = 1.0d0/scalef
      	kstart = 1

      	do m = 1, lmax - 1, 1
      		
      		rescalem = rescalem*u

		! Calculate P(m,m)
        	kstart = kstart+m+1
         
         	pmm = phase * pmm * sqr(2*m+1) / sqr(2*m)
        	p(kstart) = pmm*rescalem
        	pm2 = pmm

		! Calculate P(m+1,m)
		k = kstart+m+1
	   	pm1 = z * sqr(2*m+3) * pmm
	    	p(k) = pm1*rescalem

		! Calculate P(l,m)
               	do l = m+2, lmax, 1
               		k = k+l
                  	plm  = z*f1(k)*pm1-f2(k)*pm2
                  	p(k) = plm*rescalem
                  	pm2  = pm1
                  	pm1  = plm
               	enddo
              
      	enddo
      	
      	! Calculate P(lmax,lmax)
      	
      	rescalem = rescalem*u
      	
        kstart = kstart+m+1
        pmm = phase * pmm * sqr(2*lmax+1) / sqr(2*lmax)
        p(kstart) = pmm*rescalem
      		
end subroutine PlmON

