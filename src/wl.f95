real*8 function Wl(l, half, r, d)
!-------------------------------------------------------------------------------
!
!   This function will compute the minimum amplitude downward continuation 
!   filter of Wieczorek and Phillips 1998 for degree l, where the filter is 
!   assumed to be equal to 0.5 at degree half.
!
!   Calling Parameters
!
!       l       Spherical harmonic degree
!       half    Spherical harmonic degree where the filter is 0.5
!       r       Reference radius for surface gravity field
!       d       Mean radius of downward continuation
!
!   Dependencies:   None
!
!   Copyright (c) 2015, Mark A. Wieczorek
!   All rights reserved.
!
!-------------------------------------------------------------------------------
    implicit none
    
    integer, intent(in) :: l, half
    real*8, intent(in) ::  r, d
    real*8 :: const
    
    if (l < 0) then
        print*, "Error --- Wl"
        print*, "L must be greater or equal to zero."
        print*, "Input value = ", l
        stop
    end if
    
    if (half == 0) then
        wl=1.0d0
        
    else
        const = ( dble(2*half+1) * (r/d)**half )**2
        const = 1.0d0/const
    
        wl = 1.0d0 + const * ( dble(2*l+1)*(r/d)**l )**2
        wl = 1.0d0 / wl
        
    end if
    
end function Wl
    
    
real*8 function WlCurv(l, half, r, d)
!-------------------------------------------------------------------------------
!
!   This function will compute a minimum curvature downward continuation 
!   filter for degree l, where the filter is assumed to be equal to 0.5 at 
!   degree half.
!
!   Calling Parameters
!
!       l       Spherical harmonic degree
!       half    Spherical harmonic degree where the filter is 0.5
!       r       Reference radius for surface gravity field
!       d       Mean radius of downward continuation
!
!   Dependencies:   None
!
!   Note: This filter is analogous to (and numerically very similar to) 
!       the minimum curvature filter in Phipps Morgan and Blackman (1993)
!
!   Copyright (c) 2015, Mark A. Wieczorek
!   All rights reserved.
!
!-------------------------------------------------------------------------------
    implicit none
    
    integer, intent(in) ::  l, half
    real*8, intent(in) ::   r, d
    real*8 ::   const
    
    if (l < 0) then
        print*, "Error --- WlCurv"
        print*, "L must be greater or equal to zero."
        print*, "Input value = ", l
        stop
    end if

    if (half == 0) then
        WlCurv=1.0d0
        
    else
        const = dble(half*half+half) * ( dble(2*half+1) * (r/d)**half )**2
        const = 1.0d0/const
    
        WlCurv = 1.0d0 + const * dble(l*l+l) * ( dble(2*l+1)*(r/d)**l )**2
        WlCurv = 1.0d0 / WlCurv
        
    endif
    
end function WlCurv
