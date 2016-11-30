subroutine SphericalCapCoef(coef, theta, lmax, exitstatus)
!------------------------------------------------------------------------------
!
!   This routine will return the coefficients for a spherical
!   cap (located at the north pole) with angular radius theta,
!   using the geodesy 4-pi normalization.
!
!   Calling Parameters
!
!       IN
!           theta       Angular radius of spherical cap in RADIANS.
!
!       OUT
!           coef        m=0 coefficients normalized according to the
!                       geodesy convention, further normalized such
!                       that the degree-0 term is 1. (i.e., function
!                       has an average value of one over the entire
!                       sphere).
!
!       OPTIONAL
!           lmax        Maximum spherical harmonic degree.
!
!       OPTIONAL (OUT)
!           exitstatus  If present, instead of executing a STOP when an error
!                       is encountered, the variable exitstatus will be
!                       returned describing the error.
!                       0 = No errors;
!                       1 = Improper dimensions of input array;
!                       2 = Improper bounds for input variable;
!                       3 = Error allocating memory;
!                       4 = File IO error.
!
!   Dependencies: PlBar
!
!   Copyright (c) 2016, SHTOOLS
!   All rights reserved.
!
!-------------------------------------------------------------------------------
    use SHTOOLS, only: PlBar

    implicit none

    real*8, intent(out) :: coef(:)
    real*8, intent(in) :: theta
    integer, intent(in), optional :: lmax
    integer, intent(out), optional :: exitstatus
    real*8 :: x, top, bot, pi
    real*8, allocatable :: pl(:)
    integer :: l, lmax2, astat

    if (present(exitstatus)) exitstatus = 0

    pi = acos(-1.0d0)
    coef = 0.0d0
    x = cos(theta)
    
    if (present(lmax) ) then
        lmax2 = lmax

    else
        lmax2 = size(coef) - 1

    end if

    allocate (pl(lmax2+3), stat = astat)

    if(astat /= 0) then
        print*, "Error --- SphericalCapCoef"
        print*, "Unable to allocate array pl", astat
        if (present(exitstatus)) then
            exitstatus = 3
            return
        else
            stop
        end if

    end if

    if (present(exitstatus)) then
        call PlBar(pl, lmax2+2, x, exitstatus = exitstatus)
        if (exitstatus /= 0) return
    else
        call PlBar(pl, lmax2+2, x)
    end if

    coef(1) = 1.0d0

    bot = pl(1) - pl(2) / sqrt(3.0d0)

    do l = 1, lmax2
        top = pl(l) / sqrt(dble(2*l-1)) -pl(l+2) / sqrt(dble(2*l+3))
        coef(l+1) = top / (bot * sqrt(dble(2*l+1)))

    end do

    deallocate (pl)

end subroutine SphericalCapCoef
