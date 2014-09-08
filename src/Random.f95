real*8 function RandomN(idum)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	Minimal random number generator of Park and Miller combined with a Marsaglia 
!	shift sequence. Returns a uniform random deviate between 0.0 and 1.0 (exclusive 
!	of the endpoint values). This fully portable scalar generator has the 
!	"traditional" (not Fortran 90) calling sequence with a random deviate as the 
!	returned function value. Call with idum a negative integer to initialize; 
!	thereafter, do not alter idum except to reinitialize. The period of this
!	generator is about 3.1 10^18.
!
!	Taken from numerical recipes, f90
!
!	Copyright (c) 2005, Mark A. Wieczorek
!	All rights reserved.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	implicit none
	integer, parameter ::	K4B=selected_int_kind(9)
	integer(K4B), intent(inout) ::	idum
	integer(K4B), parameter ::	IA=16807,IM=2147483647,IQ=127773,IR=2836, &
					one = 1, param1 = 888889999, param2 = 777755555
	real*8, save ::	am
	integer(K4B), save ::	ix=-1, iy=-1, k

	if (idum <= 0 .or.iy < 0) then 	! Initialize.
		am=nearest(1.0,-1.0)/IM
		iy=ior(ieor(param1, abs(idum)), one)
		ix=ieor(param2, abs(idum))
		idum=abs(idum)+1 	! Set idum positive.
	endif
	
	ix=ieor(ix,ishft(ix,13))	! Marsaglia shift sequence with period 2^32 -1
	ix=ieor(ix,ishft(ix,-17))
	ix=ieor(ix,ishft(ix,5))
	k=iy/IQ 			! Park-Miller sequence by Schrage's method, period 2^31 -2
	iy=IA * (iy - k * IQ)-IR * k
	if (iy < 0) iy = iy + IM
	RandomN=am * ior(iand(IM, ieor(ix,iy)), one) 	! Combine the two generators with masking to
						! ensure nonzero value
						
end function RandomN


real*8 function RandomGaussian(idum)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	Returns a normally distributed deviate with 
!	zero mean and unit variance, using RandomN(idum)
!	as the source of uniform deviates. To convert to a 
!	normal distribution with standard deviation sigma,
!	multiply this function by sigma.
!
!	Taken from numerical recipes f77.
!
!	Copyright (c) 2005, Mark A. Wieczorek
!	All rights reserved.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	implicit none
	integer, parameter ::	K4B=selected_int_kind(9)
	integer(K4B), intent(inout) :: 	idum
	real*8 ::	rsq, v1, v2, RandomN
	real*8, save :: gset
	external	RandomN
	logical, save :: gaus_stored = .false.
	
	if (idum < 0) gaus_stored = .false. 	! Reinitialize.
	
	if (gaus_stored) then
		RandomGaussian = gset
		gaus_stored = .false.
	else					! We don't have an extra deviate handy, so
		do				! pick two uniform numbers in the square 
			v1 = RandomN(idum)		! extending from -1 to +1 in each direction,
			v2 = RandomN(idum)		! see if they are in the unit circle,
			v1 = 2.0d0*v1 -1.0d0	! and if they are not, try again.
			v2 = 2.0d0*v2 -1.0d0
			rsq = v1**2 + v2**2
			if (rsq > 0.0d0 .and. rsq < 1.0d0) exit
		enddo
		
		rsq = sqrt(-2.0d0*log(rsq)/rsq)		! Now make the Box-Muller transformation to get
							! two normal deviates. Return one and save
							! the other for next time.
		RandomGaussian = v1*rsq
		gset = v2*rsq
		gaus_stored = .true.
	endif
	
end function RandomGaussian

