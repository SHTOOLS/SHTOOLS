subroutine Curve2Mask(dhgrid, n, sampling, profile, nprofile, NP, extend, &
                      exitstatus)
!------------------------------------------------------------------------------
!
!   Given a list of coordinates of a single closed curve, this routine
!   will output a Driscoll and Healy sampled grid where the interior and
!   exterior of the curve are filled with zeros and ones. The value at the
!   north pole (either 0 or 1) is specified by the input parameter NP.
!
!   Longitudes can span the range from -360 to 720 degrees. If the longitudes
!   of two adjacent points differ by more than 180 degrees, it will be assumed
!   that the curve passes from 360 to 0 degrees, or from -180 to 180 degrees.
!
!   Calling Parameters
!
!       IN
!           profile     Latitude (:,1) and longitude (:,2) coordinates of a
!                       single close curve having dimension (nprofile, 2).
!           nprofile    Number of coordinates in the vector profile.
!           np          Value of the output function at the North Pole, which
!                       can be either 0 or 1.
!           n           Number of latitude bands in the output Driscoll and
!                       Healy sampled grid.
!           sampling    1 sets the number of longitude bands equal to 1,
!                       whereas 2 sets the number to twice that of n.
!
!       OUT
!           dhgrid      A Driscoll and Healy sampled grid specifiying whether
!                       the point is in the curve (1), or outside of it (0).
!
!       OPTIONAL (IN)
!           extend      If 1, return a grid that contains an additional column
!                       and row corresponding to 360 E longitude and 90 S
!                       latitude, respectively.
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
    use ftypes

    implicit none

    integer, intent(out) :: dhgrid(:,:)
    real(dp), intent(in) :: profile(:,:)
    integer, intent(in) :: n, sampling, nprofile, np
    integer, intent(in), optional :: extend
    integer, intent(out), optional :: exitstatus
    integer, parameter :: maxcross = 2000
    integer :: i, j, k, k_loc, nlat, nlong, numcross, next, ind1, ind2, &
               nlat_out, nlong_out, extend_grid
    real(dp) :: lat_int, long_int, lon, cross(maxcross), &
                cross_sort(maxcross), lat1, lat2, lon1, lon2

    if (present(exitstatus)) exitstatus = 0

    nlat = n
    lat_int = 180.0_dp / dble(nlat)
    dhgrid = 0

    if (mod(n,2) /= 0) then
        print*, "Error --- Curve2Mask"
        print*, "N must be even"
        print*, "N = ", n
        if (present(exitstatus)) then
            exitstatus = 2
            return
        else
            stop
        end if
    end if

    if (sampling == 1) then
        nlong = nlat
        long_int = 2.0_dp * lat_int
    else if (sampling == 2) then
        nlong = 2 * nlat
        long_int = lat_int
    else
        print*, "Error --- Curve2Mask"
        print*, "SAMPLING of DHGRID must be 1 (equally sampled) or 2 (equally spaced)."
        print*, "SAMPLING = ", sampling
        if (present(exitstatus)) then
            exitstatus = 2
            return
        else
            stop
        end if
    end if

    if (present(extend)) then
        if (extend == 0) then
            extend_grid = 0
            nlat_out = nlat
            nlong_out = nlong
        else if (extend == 1) then
            extend_grid = 1
            nlat_out = nlat + 1
            nlong_out = nlong + 1
        else
            print*, "Error --- Curve2Mask"
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
        extend_grid = 0
        nlat_out = nlat
        nlong_out = nlong
    end if

    if (NP /= 1 .and. NP /= 0) then
        print*, "Error --- Curve2Mask"
        print*, "NP must be 0 if the North pole is outside of curve,"
        print*, "or 1 if the North pole is inside of the curve."
        print*, "NP = ", np
        if (present(exitstatus)) then
            exitstatus = 2
            return
        else
            stop
        end if

    end if

    if (size(dhgrid(1,:)) < nlong_out .or. size(dhgrid(:,1)) < nlat_out) then
        print*, "Error --- Curve2Mask"
        print*, "DHGRID must be dimensioned as (NLAT, NLONG)."
        print*, "NLAT = ", nlat_out
        print*, "NLONG = ", nlong_out
        print*, "Size of GRID = ", size(dhgrid(:,1)), size(dhgrid(1,:))
        if (present(exitstatus)) then
            exitstatus = 1
            return
        else
            stop
        end if
    end if

    if (size(profile(:,1)) < nprofile .or. size(profile(1,:)) < 2) then
        print*, "Error --- Curve2Mask"
        print*, "PROFILE must be dimensioned as (NPROFILE, 2)."
        print*, "Dimension of NPROFILE = ", size(profile(:,1)), &
                size(profile(1,:))
        if (present(exitstatus)) then
            exitstatus = 1
            return
        else
            stop
        end if

    end if

    !--------------------------------------------------------------------------
    !
    !   Start at 90N and 0E. Determine where the curve crosses in this
    !   longitude band, sort the values, and then set the pixels between zero
    !   crossings to either 0 or 1.
    !
    !--------------------------------------------------------------------------

    do j = 1, nlong

        lon = dble(j-1) * long_int
        numcross = 0

        do i = 1, nprofile

            lat1 = profile(i, 1)
            lon1 = profile(i, 2)
            if (i /= nprofile) then
                lat2 = profile(i+1, 1)
                lon2 = profile(i+1, 2)
            else
                lat2 = profile(1, 1)
                lon2 = profile(1, 2)
            end if

            ! The output grid will be for longitudes between 0 and 360 degrees.
            ! Check if both longitudes are less than 0 or greater than 360,
            ! and then convert to this range.
            if (lon1 < 0._dp .and. lon2 < 0._dp) then
                lon1 = lon1 + 360._dp
                lon2 = lon2 + 360._dp
            end if

            if (lon1 > 360._dp .and. lon2 > 360._dp) then
                lon1 = lon1 - 360._dp
                lon2 = lon2 - 360._dp
            end if

            ! Depending on how the longitude coordinates are defined, adjacent
            ! points may differ by almost 360 degrees. Check for such jumps by
            ! looking for differences between adjacent points of more than
            ! 180 degrees.
            if (abs(lon1 - lon2) > 180._dp) then

                if (lon1 < 0._dp .or. lon2 < 0._dp) then
                    ! There is a jump from -180 to 180
                    if (lon1 < 0._dp) lon1 = lon1 + 360._dp
                    if (lon2 < 0._dp) lon2 = lon2 + 360._dp
                elseif (lon1 > 0._dp .and. lon2 > 0._dp) then
                    ! There is a jump from 360 to 0. Set the largest value
                    ! to a negative longitude.
                    if (lon1 > lon2) then
                        lon1 = lon1 - 360._dp
                    else
                        lon2 = lon2 - 360._dp
                    end if
                end if

            end if

            ! Find the latitude crossing by interpolating between the two
            ! points.
            if (lon1 <= lon .and. lon2 > lon) then
                numcross = numcross + 1

                if (numcross > maxcross) then
                    print*, "Error --- Curve2Mask"
                    print*, "Internal variable MAXCROSS needs to be increased."
                    print*, "MAXCROSS = ", maxcross
                    if (present(exitstatus)) then
                        exitstatus = 5
                        return
                    else
                        stop
                    end if

                end if

                cross(numcross) = lat1 + (lat2 - lat1) / (lon2 - lon1) &
                                  * (lon - lon1)

            else if (lon1 > lon .and. lon2 <= lon) then
                numcross = numcross + 1

                if (numcross > maxcross) then
                    print*, "Error --- Curve2Mask"
                    print*, "Internal variable MAXCROSS needs to be increased."
                    print*, "MAXCROSS = ", maxcross
                    if (present(exitstatus)) then
                        exitstatus = 5
                        return
                    else
                        stop
                    end if

                end if

                cross(numcross) = lat2 + (lat1 - lat2) / (lon1 - lon2) &
                                  * (lon - lon2)

            end if

        end do

        if (numcross == 0) then
            dhgrid(1:nlat_out, j) = np
            cycle

        else  ! sort crossings by decreasing latitude
            do k = 1, numcross
                k_loc = maxloc(cross(1:numcross), 1)
                cross_sort(k) = cross(k_loc)
                cross(k_loc) = -999.0_dp
            end do

        end if

        ind1 = 1
        next = np
        do k = 1, numcross
            ind2 = nint( (90.0_dp - cross_sort(k)) / lat_int) + 1
            dhgrid(ind1:ind2, j) = next

            if (next == 0) then
                next = 1
            else
                next = 0
            end if
            ind1 = ind2 + 1

        end do

        if (ind1 <= nlat_out) dhgrid(ind1:nlat_out, j) = next

    end do

    if (extend_grid == 1) then
        dhgrid(1:nlat_out, nlong_out) = dhgrid(1:nlat_out, 1)
    end if

end subroutine Curve2Mask
