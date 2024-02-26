subroutine MakeCircleCoord(coord, lat, lon, theta0, cinterval, cnum, &
                           exitstatus)
!------------------------------------------------------------------------------
!
!   This subroutine will return the latitude and longitude
!   coordinates of a circle of radius theta from a given
!   point. The first index in the output vectors corresponds
!   to the point directly north of the central circle coordinates,
!   and subsequent points arranged in a clockwise manner.
!
!   Calling Parameters
!
!       IN
!           lat         Latitude in degrees.
!           lon         Longitude in degrees.
!           theta0      Angular radius of the circle in degrees.
!
!       OUT
!           coord       360/interval (latitude, longitude) coordinates.
!
!       OPTIONAL, IN
!           cinterval   Angular spacing of latitude and longitude points
!                       in degrees (default=1).
!
!       OPTIONAL, OUT
!           cnum        Number of points in the output vectors.
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
    use ftypes

    implicit none

    real(dp), intent(in) :: lat, lon, theta0
    real(dp), intent(out) :: coord(:,:)
    real(dp), intent(in), optional :: cinterval
    integer(int32), intent(out), optional :: cnum, exitstatus
    real(dp) :: pi, interval, xold, yold, zold, x, y, z, x1, phi
    integer(int32) :: k, num

    if (present(exitstatus)) exitstatus = 0

    if (present(cinterval)) then
        interval = cinterval
    else
        interval = 1.0_dp
    end if

    num = int(360.0_dp / interval)

    if (present(cnum)) then
        cnum = num
    end if

    if (size(coord(:,1)) < num .or. size(coord(1,:)) < 2) then
        print*, "Error --- MakeCircleCoord"
        print*, "COORD must be dimensioned as (NUM, 2) where NUM is ", num
        print*, "Input array is dimensioned as ", size(coord(:,1)), &
                size(coord(1,:))
        if (present(exitstatus)) then
            exitstatus = 1
            return
        else
            stop
        end if
    end if

    pi = acos(-1.0_dp)

    !--------------------------------------------------------------------------
    !
    !   Treat the case where theta0 = 0 separately.
    !
    !--------------------------------------------------------------------------

    if (theta0 == 0.0_dp) then
        coord(1:cnum,1) = lat
        coord(1:cnum,2) = lon

        return

    end if

    !--------------------------------------------------------------------------
    !
    !   Calculate grid points. First create a circle, then rotate these
    !   points.
    !
    !--------------------------------------------------------------------------

    zold = cos(theta0 * pi / 180.0_dp)

    do k = 1, num
        phi = pi - dble(k-1) * (2.0_dp * pi / dble(num))
        xold = sin(theta0 * pi / 180.0_dp) * cos(phi)
        yold = sin(theta0 * pi / 180.0_dp) * sin(phi)

        ! rotate coordinate system 90-lat degrees about y axis

        x1 = xold * cos(pi / 2.0_dp - lat * pi / 180.0_dp) + &
             zold * sin(pi / 2.0_dp - lat * pi / 180.0_dp)
        z = -xold * sin(pi / 2.0_dp - lat * pi / 180.0_dp) + &
             zold * cos(pi / 2.0_dp - lat * pi / 180.0_dp)

        ! rotate coordinate system lon degrees about z axis

        x = x1 * cos(lon * pi / 180.0_dp) - yold * sin(lon * pi / 180.0_dp)
        y = x1 * sin(lon * pi / 180.0_dp) + yold * cos(lon * pi / 180.0_dp)

        coord(k,1) = (pi / 2.0_dp - acos(z / sqrt(x**2 + y**2 + z**2))) * &
                     180.0_dp / pi
        coord(k,2) = atan2(y, x) * 180.0_dp / pi

    end do

end subroutine MakeCircleCoord
