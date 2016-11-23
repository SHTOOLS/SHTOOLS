subroutine MakeEllipseCoord(coord, lat, lon, dec, A_theta, B_theta, &
                            cinterval, cnum, exitstatus)
!------------------------------------------------------------------------------
!
!   This subroutine will return the latitude and longitude
!   coordinates of an ellipse with semi-major and semi-minor axes
!   A_theta and B_theta (in degrees) at postition LAT and LON (in degrees)
!   and with a clockwise rotation of the semi-major axis DEC with respect
!   to local north (in degrees).
!
!   Calling Parameters
!
!       IN
!           lat         Latitude and longitude in degrees.
!           lon         Latitude and longitude in degrees.
!           A_theta     Angular radius of the semi-major axis in degrees.
!           B_theta     Angular radius of the semi-minor axis in degrees.
!           dec         Clockwise rotation of the semi-major axis with
!                       respect to local north in degrees.
!
!       OUT
!           coord       360/interval (latitude, longitude) coordinates.
!
!       OPTIONAL, IN
!           cinterval   Angular spacing of latitude and longitude points
!                       in degrees (default=1).
!           exitstatus  If present, instead of executing a STOP when an error
!                       is encountered, the variable exitstatus will be
!                       returned describing the error.
!                       0 = No errors;
!                       1 = Improper dimensions of input array;
!                       2 = Improper bounds for input variable;
!                       3 = Error allocating memory;
!                       4 = File IO error.
!
!   Dependencies: None
!
!   Copyright (c) 2016, SHTOOLS
!   All rights reserved.
!
!-------------------------------------------------------------------------------
    implicit none

    real*8, intent(in) :: lat, lon, A_theta, B_theta, dec
    real*8, intent(out) :: coord(:,:)
    real*8, intent(in), optional :: cinterval
    integer, intent(out), optional :: cnum, exitstatus
    real*8 :: pi, interval, xold, yold, zold, x, y, z, x1, phi, r
    integer ::  k, num

    if (present(exitstatus)) exitstatus = 0

    if (present(cinterval)) then
        interval = cinterval
    else
        interval = 1.0d0
    end if

    num = int(360.0d0 / interval)

    if (present(cnum)) then
        cnum = num
    end if

    if (size(coord(:,1)) < num .or. size(coord(1,:)) < 2) then
        print*, "Error --- MakeEllipseCoord"
        print*, "COORD must be dimensioned as (NUM, 2) where NUM is ", NUM
        print*, "Input array is dimensioned as ", size(coord(:,1)), &
                size(coord(1,:))
        if (present(exitstatus)) then
            exitstatus = 1
            return
        else
            stop
        end if
    end if

    pi = acos(-1.0d0)

    !--------------------------------------------------------------------------
    !
    !   Calculate grid points. First create a cirlce, then rotate these
    !   points.
    !
    !--------------------------------------------------------------------------

    do k = 1, num

        phi = pi - dble(k-1) * (2.0d0 * pi / dble(num))
        r = a_theta * b_theta / sqrt( (b_theta*cos(phi))**2 &
            + (a_theta * sin(phi))**2)
        xold = sin(r * pi / 180.0d0) * cos(phi - dec * pi / 180.0d0)
        yold = sin(r * pi / 180.0d0) * sin(phi - dec * pi / 180.0d0)
        zold = cos(r * pi / 180.0d0)

        ! rotate coordinate system 90-lat degrees about y axis

        x1 = xold * cos(pi / 2.0 - lat * pi / 180.0d0) + &
             zold * sin(pi / 2.0 - lat * pi / 180.0d0)
        z = -xold * sin(pi / 2.0 - lat * pi / 180.0d0) + &
            zold * cos(pi / 2.0 - lat * pi / 180.0d0)

        ! rotate coordinate system lon degrees about z axis

        x = x1 * cos(lon * pi / 180.0d0) - yold * sin(lon * pi / 180.0d0)
        y = x1 * sin(lon * pi / 180.0d0) + yold * cos(lon * pi/180.0d0)

        coord(k, 1) = (pi / 2.0d0 - acos(z / sqrt(x**2 + y**2 + z**2))) * &
                      180.0d0 / pi
        coord(k, 2) = atan2(y, x) * 180.0d0 / pi

    enddo

end subroutine MakeEllipseCoord
