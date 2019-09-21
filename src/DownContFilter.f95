function DownContFilterMA(l, half, r, d)
!------------------------------------------------------------------------------
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
!   Copyright (c) 2005-2019, SHTOOLS
!   All rights reserved.
!
!------------------------------------------------------------------------------
    use ftypes

    implicit none

    real(dp) :: DownContFilterMA
    integer, intent(in) :: l, half
    real(dp), intent(in) :: r, d
    real(dp) :: const

    if (l < 0) then
        print*, "Error --- DownContFilterMA"
        print*, "L must be greater or equal to zero."
        print*, "Input value = ", l
        stop
    end if

    if (half == 0) then
        DownContFilterMA = 1.0_dp

    else
        const = (dble(2*half+1) * (r / d)**half )**2
        const = 1.0_dp / const

        DownContFilterMA = 1.0_dp + const * ( dble(2*l+1) * (r / d)**l )**2
        DownContFilterMA = 1.0_dp / DownContFilterMA

    end if

end function DownContFilterMA


function DownContFilterMC(l, half, r, d)
!------------------------------------------------------------------------------
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
!   Note: This filter is analogous to (and numerically very similar to)
!       the minimum curvature filter in Phipps Morgan and Blackman (1993)
!
!   Copyright (c) 2005-2019, SHTOOLS
!   All rights reserved.
!
!------------------------------------------------------------------------------
    use ftypes

    implicit none

    real(dp) :: DownContFilterMC
    integer, intent(in) :: l, half
    real(dp), intent(in) :: r, d
    real(dp) :: const

    if (l < 0) then
        print*, "Error --- DownContFilterMC"
        print*, "L must be greater or equal to zero."
        print*, "Input value = ", l
        stop
    end if

    if (half == 0) then
        DownContFilterMC = 1.0_dp

    else
        const = dble(half * half+half) * ( dble(2*half+1) * (r / d)**half)**2
        const = 1.0_dp / const

        DownContFilterMC = 1.0_dp + const * dble(l*l+l) * (dble(2*l+1) * (r/d)**l)**2
        DownContFilterMC = 1.0_dp / DownContFilterMC

    end if

end function DownContFilterMC
