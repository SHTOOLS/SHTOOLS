subroutine GLQGridCoord(latglq, longlq, lmax, nlat, nlong, extend, exitstatus)
!------------------------------------------------------------------------------
!
!   Given a maximum spherical harmonic degree lmax, this routine
!   will determine the latitude and longitude coordinates associated with
!   grids that are used in the Gauss-Legendre quadratue spherical harmonic
!   expansion routines. The coordinates are output in DEGREES.
!
!   Calling Parameters
!       IN
!           lmax        Maximum spherical harmonic degree of the expansion.
!
!       OUT
!           latglq      Array of latitude points used in Gauss-Legendre
!                       grids, in degrees.
!           longlq      Array of longitude points used in Gauss-Legendre
!                       grids, in degrees.
!           nlat        Number of latitude points.
!           nlong       Number of longitude points.
!
!       OPTIONAL (IN)
!           extend      If 1, return a grid that contains an additional column
!                       for 360 E longitude.
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
!   Copyright (c) 2005-2019, SHTOOLS
!   All rights reserved.
!
!------------------------------------------------------------------------------
    use SHTOOLS, only: PreGLQ
    use ftypes

    implicit none

    integer, intent(in) :: lmax
    integer, intent(out) :: nlat, nlong
    real(dp), intent(out) :: latglq(:), longlq(:)
    integer, intent(in), optional :: extend
    integer, intent(out), optional :: exitstatus
    real(dp) :: pi, upper, lower, zero(lmax+1), w(lmax+1)
    integer :: i, nlong_out

    if (present(exitstatus)) exitstatus = 0

    nlat = lmax + 1
    nlong = 2 * lmax + 1

    if (present(extend)) then
        if (extend == 0) then
            nlong_out = nlong
        else if (extend == 1) then
            nlong_out = nlong + 1
        else
            print*, "Error --- GLQGridCoord"
            print*, "Optional parameter EXTEND must be 0 or 1."
            print*, "Input value is ", extend
            if (present(exitstatus)) then
                exitstatus = 2
                return
            else
                stop
            end if
        end if
    else
        nlong_out = nlong
    end if

    if (size(latglq) < nlat) then
        print*, "Error --- GLQGridCoord"
        print*, "LATGLQ must be dimensioned as (LMAX+1) where LMAX is ", lmax
        print*, "Input array is dimensioned as ", size(latglq)
        if (present(exitstatus)) then
            exitstatus = 1
            return
        else
            stop
        end if

    else if (size(longlq) < nlong_out) then
        print*, "Error --- GLQGridCoord"
        print*, "LONGLQ must be dimensioned as ", nlong_out
        print*, "Input array is dimensioned as ", size(longlq)
        if (present(exitstatus)) then
            exitstatus = 1
            return
        else
            stop
        end if

    end if

    pi = acos(-1.0_dp)

    upper = 1.0_dp
    lower = -1.0_dp

    if (present(exitstatus)) then
        call PreGLQ(lower, upper, nlat, zero, w, exitstatus = exitstatus)
        if (exitstatus /= 0) return
    else
        call PreGLQ(lower, upper, nlat, zero, w)
    end if

    do i = 1, nlong_out
        longlq(i) = 360.0_dp * dble(i-1) / dble(nlong)
    end do

    do i = 1, nlat
        latglq(i) = asin(zero(i)) * 180.0_dp / pi
    end do

end subroutine GLQGridCoord
