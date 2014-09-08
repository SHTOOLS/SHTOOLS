real*8 function Wl(l, half, r, d)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	This function will compute the minimum amplitude downward continuation 
!	filter of Wieczorek and Phillips 1998 for degree l, where the filter is assumed
!	to be equal to 0.5 at degree half.
!
!	Calling Parameters:
!		l	spherical harmonic degree
!		half 	spherical harmonic degree where the filter is 0.5
!		r 	reference radius for surface gravity field
!		d	mean radius of downward continuation
!
!	Dependencies:	None
!
!	Written by Mark Wieczorek (May 2004)
!
!	Copyright (c) 2005-2006, Mark A. Wieczorek
!	All rights reserved.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	implicit none
	integer, intent(in) ::	l, half
	real*8, intent(in) ::	r, d
	real*8 ::		const
	
	if (half == 0) then
	
		wl=1.0d0
		
	else
	
		const = ( dble(2*half+1) * (r/d)**half )**2
		const = 1.0d0/const
	
		wl = 1.0d0 + const * ( dble(2*l+1)*(r/d)**l )**2
		wl = 1.0d0/wl
		
	endif
	
end function Wl
	
	
real*8 function WlCurv(l, half, r, d)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	This function will compute a minimum curvature downward continuation 
!	filter for degree l, where the filter is assumed to be equal to 0.5 at degree half.
!
!	Calling Parameters:
!		l	spherical harmonic degree
!		half 	spherical harmonic degree where the filter is 0.5
!		r 	reference radius for surface gravity field
!		d	mean radius of downward continuation
!
!	Dependencies:	None
!
!	Note: This filter is analogous to (and numerically very similar to) 
!		the minimum curvature filter in Phipps Morgan and Blackman (1993)
!
!	Written by Mark Wieczorek (May 2004)
!
!	Copyright (c) 2005-2006, Mark A. Wieczorek
!	All rights reserved.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	implicit none
	integer, intent(in) ::	l, half
	real*8, intent(in) ::	r, d
	real*8 ::	const
	
	if (half == 0) then
		WlCurv=1.0d0
	else
		const = dble(half*half+half) * ( dble(2*half+1) * (r/d)**half )**2
		const = 1.0d0/const
	
		WlCurv = 1.0d0 + const * dble(l*l+l) * ( dble(2*l+1)*(r/d)**l )**2
		WlCurv = 1.0d0/WlCurv
	endif
	
end function WlCurv

	