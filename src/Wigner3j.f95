subroutine Wigner3j(w3j, jmin, jmax, j2, j3, m1, m2, m3)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	This subroutine will calculate the Wigner 3j symbols
!
!		j  j2 j3
!		m1 m2 m3
!
!	for all allowable values of j. The returned values in the array j are 
!	calculated only for the limits
!
!		jmin = max(|j2-j3|, |m1|)
!		jmax = j2 + j3
!
!	To be non-zero, m1 + m2 + m3 = 0. In addition, it is assumed that all j and m are 
!	integers. Returned values have a relative error less than ~1.d-8 when j2 and j3 
!	are less than 103 (see below). In practice, this routine is probably usable up to 165.
!
!	This routine is based upon the stable non-linear recurence relations of Luscombe and 
!	Luban (1998) for the "non classical" regions near jmin and jmax. For the classical 
!	region, the standard three term recursion relationship is used (Schulten and Gordon 1975). 
!	Note that this three term recursion can be unstable and can also lead to overflows. Thus 
!	the values are rescaled by a factor "scalef" whenever the absolute value of the 3j coefficient 
!	becomes greater than unity. Also, the direction of the iteration starts from low values of j
!	to high values, but when abs(w3j(j+2)/w3j(j)) is less than one, the iteration will restart 
!	from high to low values. More efficient algorithms might be found for specific cases 
!	(for instance, when all m's are zero).
!
!	Verification: 
!
!	The results have been verified against this routine run in quadruple precision.
!	For 1.e7 acceptable random values of j2, j3, m2, and m3 between -200 and 200, the relative error
!	was calculated only for those 3j coefficients that had an absolute value greater than 
!	1.d-17 (values smaller than this are for all practical purposed zero, and can be heavily 
!	affected by machine roundoff errors or underflow). 853 combinations of parameters were found
!	to have relative errors greater than 1.d-8. Here I list the minimum value of max(j2,j3) for
!	different ranges of error, as well as the number of times this occured
!	
!	1.d-7 < error  <=1.d-8 = 103	# = 483
!	1.d-6 < error <= 1.d-7 =  116	# = 240
!	1.d-5 < error <= 1.d-6 =  165	# = 93
!	1.d-4 < error <= 1.d-5 = 167	# = 36
!
!	Many times (maybe always), the large relative errors occur when the 3j coefficient 
!	changes sign and is close to zero. (I.e., adjacent values are about 10.e7 times greater 
!	in magnitude.) Thus, if one does not need to know highly accurate values of the 3j coefficients
!	when they are almost zero (i.e., ~1.d-10) then this routine is probably usable up to about 160.
!
!	These results have also been verified for parameter values less than 100 using a code
!	based on the algorith of de Blanc (1987), which was originally coded by Olav van Genabeek, 
!	and modified by M. Fang (note that this code was run in quadruple precision, and
!	only calculates one coefficient for each call. I also have no idea if this code
!	was verified.) Maximum relative errors in this case were less than 1.d-8 for a large number
!	of values (again, only 3j coefficients greater than 1.d-17 were considered here).
!	
!	The biggest improvement that could be made in this routine is to determine when one should
!	stop iterating in the forward direction, and start iterating from high to low values. 
!
!	Calling parameters
!		IN	
!			j2, j3, m1, m2, m3 	Integer values.
!		OUT	
!			w3j			Array of length jmax - jmin + 1.
!			jmin, jmax		Minimum and maximum values
!						out output array.
!	Dependencies: None
!	
!	Written by Mark Wieczorek August (2004)
!
!	August 2009: Based on the suggestions of Roelof Rietbroek, the calculation of RS has been slightly
!	modified so that division by zero will not cause a run time crash (this behavior depends on how the 
!	compiler treats IEEE floating point exceptions). These values were never used in the original code 
!	when this did occur.
!
!	Copyright (c) 2005-2009, Mark A. Wieczorek
!	All rights reserved.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	implicit none
	integer, intent(in) ::	j2, j3, m1, m2, m3
	integer, intent(out) ::	jmin, jmax
	real*8, intent(out) ::	w3j(:)
	real*8 ::		wnmid, wpmid, scalef, denom, rs(j2+j3+1), &
				wl(j2+j3+1), wu(j2+j3+1), xjmin, yjmin, yjmax, zjmax, xj, zj
	integer :: 		j, jnum, jp, jn, k, flag1, flag2, jmid
	
	
	if (size(w3j) < j2+j3+1) then
		print*, "Error --- Wigner3j"
		print*, "W3J must be dimensioned (J2+J3+1) where J2 and J3 are ", j2, j3
		print*, "Input array is dimensioned ", size(w3j)
		stop
	endif
	
	w3j = 0.0d0
	
	flag1 = 0
	flag2 = 0
	
	scalef = 1.0d3
	
	jmin = max(abs(j2-j3), abs(m1))
	jmax = j2 + j3
	jnum = jmax - jmin + 1
	
	if (abs(m2) > j2 .or. abs(m3) > j3) then
		return
	elseif (m1 + m2 + m3 /= 0) then
		return
	elseif (jmax < jmin) then
		return
	endif
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!
	!	Only one term is present
	!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	if (jnum == 1) then
		w3j(1) = 1.0d0 / sqrt(2.0d0*jmin+1.0d0)
		if ( (w3j(1) < 0.0d0 .and. (-1)**(j2-j3+m2+m3) > 0) .or. &
			(w3j(1) > 0.0d0 .and. (-1)**(j2-j3+m2+m3) < 0) ) &
			w3j(1) = -w3j(1)
		return	
	endif
		
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!
	! 	Calculate lower non-classical values for [jmin, jn]. If the second term
	!	can not be calculated because the recursion relationsips give rise to a
	!	1/0, then set flag1 to 1.  If all m's are zero, then this is not a problem 
	!	as all odd terms must be zero.
	!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	rs = 0.0d0
	wl = 0.0d0
	
	xjmin = x(jmin)
	yjmin = y(jmin)
	
	if (m1 == 0 .and. m2 == 0 .and. m3 == 0) then		! All m's are zero
	
		wl(jindex(jmin)) = 1.0d0
		wl(jindex(jmin+1)) = 0.0d0
		jn = jmin+1
		
	elseif (yjmin == 0.0d0) then				! The second terms is either zero
	
		if (xjmin == 0.0d0) then			! or undefined
			flag1 = 1
			jn = jmin
		else
			wl(jindex(jmin)) = 1.0d0
			wl(jindex(jmin+1)) = 0.0d0
			jn = jmin+1
		endif
		
	elseif ( xjmin * yjmin >= 0.0d0) then			! The second term is outside of the 
								! non-classical region 
		wl(jindex(jmin)) = 1.0d0
		wl(jindex(jmin+1)) = -yjmin / xjmin
		jn = jmin+1
		
	else							! Calculate terms in the non-classical region
	
		rs(jindex(jmin)) = -xjmin / yjmin
		
		jn = jmax
		do j=jmin + 1, jmax-1, 1
			denom =  y(j) + z(j)*rs(jindex(j-1))
			xj = x(j)
			if (abs(xj) > abs(denom) .or. xj * denom >= 0.0d0 .or. denom == 0.0d0) then
				jn = j-1
				exit
			else
				rs(jindex(j)) = -xj / denom
			endif
				
		enddo
		
		wl(jindex(jn)) = 1.0d0
		
		do k=1, jn - jmin, 1
			wl(jindex(jn-k)) = wl(jindex(jn-k+1)) * rs(jindex(jn-k))
		enddo
		
		if (jn == jmin) then					! Calculate at least two terms so that
			wl(jindex(jmin+1)) = -yjmin / xjmin		! these can be used in three term
			jn = jmin+1					! recursion
			
		endif

	endif
	
	if (jn == jmax) then					! All terms are calculated
	
		w3j(1:jnum) = wl(1:jnum)
		call normw3j
		call fixsign
		
		return

	endif

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!
	! 	Calculate upper non-classical values for [jp, jmax].
	!	If the second last term can not be calculated because the
	!	recursion relations give a 1/0, then set flag2 to 1. 
	!	(Note, I don't think that this ever happens).
	!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	wu = 0.0d0
	
	yjmax = y(jmax)
	zjmax = z(jmax)
	
	if (m1 == 0 .and. m2 == 0 .and. m3 == 0) then
	
		wu(jindex(jmax)) = 1.0d0
		wu(jindex(jmax-1)) = 0.0d0
		jp = jmax-1
		
	elseif (yjmax == 0.0d0) then
	
		if (zjmax == 0.0d0) then
			flag2 = 1
			jp = jmax
		else
			wu(jindex(jmax)) = 1.0d0
			wu(jindex(jmax-1)) = - yjmax / zjmax
			jp = jmax-1
		endif
		
	elseif (yjmax * zjmax >= 0.0d0) then
	
		wu(jindex(jmax)) = 1.0d0
		wu(jindex(jmax-1)) = - yjmax / zjmax
		jp = jmax-1

	else
		rs(jindex(jmax)) = -zjmax / yjmax

		jp = jmin
		do j=jmax-1, jn, -1
			denom = y(j) + x(j)*rs(jindex(j+1))
			zj = z(j)
			if (abs(zj) > abs(denom) .or. zj * denom >= 0.0d0 .or. denom == 0.0d0) then
				jp = j+1
				exit
			else
				rs(jindex(j)) = -zj / denom
			endif
		enddo	
		
		wu(jindex(jp)) = 1.0d0
		
		do k=1, jmax - jp, 1
			wu(jindex(jp+k)) = wu(jindex(jp+k-1))*rs(jindex(jp+k))
		enddo
		
		if (jp == jmax) then
			wu(jindex(jmax-1)) = - yjmax / zjmax
			jp = jmax-1
		endif
		
	endif
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!
	! 	Calculate classical terms for [jn+1, jp-1] using standard three
	! 	term rercusion relationship. Start from both jn and jp and stop at the
	! 	midpoint. If flag1 is set, then perform the recursion solely from high to
	! 	low values. If flag2 is set, then perform the recursion solely from low to high.
	!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	if (flag1 == 0) then
	
		jmid = (jn + jp)/2
		
		do j=jn, jmid - 1, 1			
			wl(jindex(j+1)) = - (z(j)*wl(jindex(j-1)) +y(j)*wl(jindex(j))) / x(j)
			
			if (abs(wl(jindex(j+1))) > 1.0d0) then				! watch out for overflows.
				wl(jindex(jmin):jindex(j+1)) = wl(jindex(jmin):jindex(j+1)) / scalef
			endif
			
			if (abs(wl(jindex(j+1)) / wl(jindex(j-1))) < 1.0d0 .and. &	! if values are decreasing
				wl(jindex(j+1)) /= 0.0d0) then				! then stop upward iteration
				jmid = j+1						! and start with the downward
				exit							! iteration.
			endif
		enddo
		
		wnmid = wl(jindex(jmid))
		
		if (abs(wnmid/wl(jindex(jmid-1))) < 1.d-6 .and. &
			wl(jindex(jmid-1)) /= 0.0d0) then				! Make sure that the stopping
			wnmid = wl(jindex(jmid-1))					! midpoint value is not a zero,
			jmid = jmid - 1							! or close to it!
		endif
		
		
		do j=jp, jmid+1, -1
			wu(jindex(j-1)) = - (x(j)*wu(jindex(j+1)) + y(j)*wu(jindex(j)) ) / z(j)
			if (abs(wu(jindex(j-1))) > 1.0d0) then
				wu(jindex(j-1):jindex(jmax)) = wu(jindex(j-1):jindex(jmax)) / scalef
			endif
	
		enddo
		
		wpmid = wu(jindex(jmid))
		
		! rescale two sequences to common midpoint
		
		if (jmid == jmax) then
			w3j(1:jnum) = wl(1:jnum)
		elseif (jmid == jmin) then
			w3j(1:jnum) = wu(1:jnum)
		else
			w3j(1:jindex(jmid)) = wl(1:jindex(jmid)) * wpmid / wnmid 
			w3j(jindex(jmid+1):jindex(jmax)) = wu(jindex(jmid+1):jindex(jmax))
		endif
		
	elseif (flag1 == 1 .and. flag2 == 0) then	! iterature in downward direction only
		
		do j=jp, jmin+1, -1
			wu(jindex(j-1)) = - (x(j)*wu(jindex(j+1)) + y(j)*wu(jindex(j)) ) / z(j)
			if (abs(wu(jindex(j-1))) > 1) then
				wu(jindex(j-1):jindex(jmax)) = wu(jindex(j-1):jindex(jmax)) / scalef
			endif
		enddo
		
		w3j(1:jnum) = wu(1:jnum)
		
	elseif (flag2 == 1 .and. flag1 == 0) then	! iterature in upward direction only
		
		do j=jn, jp-1, 1
			wl(jindex(j+1)) = - (z(j)*wl(jindex(j-1)) +y(j)*wl(jindex(j))) / x(j)
			if (abs(wl(jindex(j+1))) > 1) then
				wl(jindex(jmin):jindex(j+1)) = wl(jindex(jmin):jindex(j+1))/ scalef
			endif
		enddo
		
		w3j(1:jnum) = wl(1:jnum)
		
	elseif (flag1 == 1 .and. flag2 == 1) then

		print*, "Fatal Error --- Wigner3j"
		print*, "Can not calculate function for input values, both flag1 and flag 2 are set."
		stop
	endif

	
	call normw3j
	call fixsign
		
	
	contains
	
		integer function jindex(j)
			integer :: j
			jindex = j-jmin+1
		end function jindex
	
		real*8 function a(j)
			integer :: j
			a = (dble(j)**2 - dble(j2-j3)**2) * (dble(j2+j3+1)**2 - dble(j)**2) * (dble(j)**2-dble(m1)**2)
			a = sqrt(a)
		end function a
		
		real*8 function y(j)
			integer :: j
			y = -dble(2*j+1) * &
				( dble(m1) * (dble(j2)*dble(j2+1) - dble(j3)*dble(j3+1) ) - dble(m3-m2)*dble(j)*dble(j+1) )
		end function y

		real*8 function x(j)	
			integer :: j
			x = dble(j) * a(j+1)
		end function x
		
		real*8 function z(j)
			integer :: j
			z = dble(j+1)*a(j)
		end function z
		
		subroutine normw3j
			real*8:: norm
			integer j
			
			norm = 0.0d0
			do j = jmin, jmax
				norm = norm + dble(2*j+1) * w3j(jindex(j))**2
			enddo
			
			w3j(1:jnum) = w3j(1:jnum) / sqrt(norm)
			
		end subroutine normw3j
		
		subroutine fixsign
		
			if ( (w3j(jindex(jmax)) < 0.0d0 .and. (-1)**(j2-j3+m2+m3) > 0) .or. &
				(w3j(jindex(jmax)) > 0.0d0 .and. (-1)**(j2-j3+m2+m3) < 0) ) then
				w3j(1:jnum) = -w3j(1:jnum)
			endif
			
		end subroutine fixsign

		
end subroutine Wigner3j


