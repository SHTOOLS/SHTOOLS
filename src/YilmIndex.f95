integer function YilmIndex(i, l, m)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!	This function will give the index in a 1-dimensional array of the spherical 
!	harmonic coefficient corresponding to the element Cilm. The elements of the 
!	1D vector array are packed according to (where positive m corresponds 
!	to i=1 (the cosine coefficients), and negative m corresponds to i=1 
!	(the sine coefficients))
!
!	0,0
!	1, 0; 1, 1; 1, -1
!	2,0 ; 2, 1; 2, 2; 2, -1; 2, -2
!
!	This mapping is given by the function:
!
!		YilmIndex = 1 + l**2 + (i-1)*l + m
!	
!
!	Written by Mark Wieczorek (August 2009)
!
!	Copyright (c) 2009, Mark A. Wieczorek
!	All rights reserved.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	implicit none
	integer, intent(in) ::	i, l, m
	
	if (i /= 1 .and. i /= 2) then
		print*, "Error --- YilmIndex"
		print*, "I must be 1 (for cosine terms) or 2 (for sine terms)."
		print*, "I = ", i
		stop
	endif
	
	if (l < 0) then
		print*, "Error --- YilmIndex"
		print*, "L must be positive."
		print*, "L = ", l
		stop
	endif
	
	if (m < 0 .or. m > l) then
		print*, "Error --- YilmIndex"
		print*, "M must be positive and less than L."
		print*, "M = ", m
		print*, "L = ", l
		stop
	endif
	
	if (m == 0 .and. i == 2) then
		print*, "Error --- YilmIndex"
		print*, "When M = 0, I must be 1."
		print*, "I = ", i
		print*, "M = ", m
		stop
	endif
	
	yilmindex = 1 + l**2 + (i-1)*l + m	
	
end function YilmIndex
