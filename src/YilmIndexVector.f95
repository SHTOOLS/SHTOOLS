integer function YilmIndexVector(i, l, m)
!-------------------------------------------------------------------------------
!
!   This function will give the index in a 1-dimensional array of the spherical 
!   harmonic coefficient corresponding to the element Cilm. The elements of the 
!   1D vector array are packed according to (where positive m corresponds 
!   to i=1 (the cosine coefficients), and negative m corresponds to i=1 
!   (the sine coefficients))
!
!   0,0
!   1, 0; 1, 1; 1, -1
!   2,0 ; 2, 1; 2, 2; 2, -1; 2, -2
!
!   This mapping is given by the function:
!
!       YilmIndex = l**2 + (i-1)*l + m + 1
!
!   Copyright (c) 2016, SHTOOLS
!   All rights reserved.
!
!-------------------------------------------------------------------------------
    implicit none
    
    integer, intent(in) :: i, l, m
    
    if (i /= 1 .and. i /= 2) then
        print*, "Error --- YilmIndexVector"
        print*, "I must be 1 (for cosine terms) or 2 (for sine terms)."
        print*, "I = ", i
        stop
        
    end if
    
    if (l < 0) then
        print*, "Error --- YilmIndexVector"
        print*, "L must be positive."
        print*, "L = ", l
        stop
        
    end if
    
    if (m < 0 .or. m > l) then
        print*, "Error --- YilmIndexVector"
        print*, "M must be positive and less than L."
        print*, "M = ", m
        print*, "L = ", l
        stop
        
    end if
    
    if (m == 0 .and. i == 2) then
        print*, "Error --- YilmIndexVector"
        print*, "When M = 0, I must be 1."
        print*, "I = ", i
        print*, "M = ", m
        stop
        
    end if
    
    yilmindexvector = l**2 + (i-1)*l + m + 1
    
end function YilmIndexVector
