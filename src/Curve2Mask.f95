subroutine Curve2Mask(dhgrid, n, sampling, profile, nprofile, NP, &
                      centralmeridian, exitstatus)
!------------------------------------------------------------------------------
!
!   Given a list of coordinates of a SINGLE CLOSED CURVE, this routine
!   will fill the interior and exterior with 0s and 1s. The value at the
!   north pole (either 0 or 1) is specified by the input parameter NP.
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
!           centralmeridian     If 1, the curve is assumed to pass through the
!                               central meridian: passing from < 360 degrees
!                               to > 0 degrees. The curve makes a complete
!                               circle about the planet in longitude. default
!                               is zero.
!
!       OUT
!           dhgrid      A Driscoll and Healy sampled grid specifiying whether
!                       the point is in the curve (1), or outside of it (0).
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
!   Copyright (c) 2016, SHTOOLS
!   All rights reserved.
!
!------------------------------------------------------------------------------
    implicit none

    integer, intent(out) :: dhgrid(:,:)
    real*8, intent(in) ::   profile(:,:)
    integer, intent(in) ::  n, sampling, nprofile, np
    integer, intent(in), optional :: centralmeridian
    integer, intent(out), optional :: exitstatus
    integer, parameter ::   maxcross = 2000
    integer ::  i, j, k, k_loc, nlat, nlong, numcross, next, ind1, ind2, cm
    real*8 :: lat_int, long_int, lon, cross(maxcross), cross_sort(maxcross)

    if (present(exitstatus)) exitstatus = 0

    nlat = n
    lat_int = 180.0d0 / dble(nlat)
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
        long_int = 2.0d0 * lat_int
        
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
    
    if (NP /=1 .and. NP /= 0) then
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

    if (size(dhgrid(1,:)) < nlong .or. size(dhgrid(:,1)) < nlat ) then
        print*, "Error --- Curve2Mask"
        print*, "DHGRID must be dimensioned as (NLAT, NLONG)."
        print*, "NLAT = ", nlat
        print*, "NLONG = ", nlong
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
    
    if (present(centralmeridian)) then
        if (centralmeridian /= 0 .and. centralmeridian /= 1) then
             print*, "Error --- Curve2Mask"
             print*, "CENTRALMERIDIAN must be either 0 or 1."
             print*, "Input value is ", centralmeridian
            if (present(exitstatus)) then
                exitstatus = 2
                return
            else
                stop
            end if

        endif
        cm = centralmeridian

    else
        cm = 0

    endif

    !--------------------------------------------------------------------------
    !
    !   Start at 90N and 0E. Determine where the curve crosses this longitude
    !   band, sort the values, and then set the pixels between zero crossings
    !   to either 0 or 1.
    !
    !--------------------------------------------------------------------------

    do j = 1, nlong

        lon = dble(j-1) * long_int
        numcross = 0

        do i = 1, nprofile - 1

            if (profile(i,2) <= lon .and. profile(i+1,2) > lon) then
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

                cross(numcross) = profile(i,1) + (profile(i+1,1)-profile(i,1)) &
                                  / (profile(i+1,2)-profile(i,2)) &
                                  * (lon - profile(i,2))

            else if (profile(i,2) > lon .and. profile(i+1,2) <= lon) then
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

                cross(numcross) = profile(i+1,1) + &
                                  (profile(i,1)-profile(i+1,1)) / &
                                  (profile(i,2)-profile(i+1,2)) &
                                  * (lon - profile(i+1,2))

            end if

        end do

        ! do first and last points
        if (cm == 0) then

            if (profile(nprofile,2) <= lon .and. profile(1,2) > lon) then
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

                cross(numcross) = profile(nprofile,1) + &
                                  (profile(1,1)-profile(nprofile,1)) / &
                                  (profile(1,2)-profile(nprofile,2)) * &
                                  (lon - profile(nprofile,2))

            else if (profile(nprofile,2) > lon .and. profile(1,2) <= lon) then
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

                cross(numcross) = profile(1,1) + &
                                  (profile(nprofile,1)-profile(1,1)) / &
                                  (profile(nprofile,2)-profile(1,2)) &
                                  * (lon - profile(1,2))

            end if

        end if

        if (numcross > 0) then  ! sort crossings by decreasing latitude
            do k = 1, numcross
                k_loc = maxloc(cross(1:numcross), 1)
                cross_sort(k) = cross(k_loc)
                cross(k_loc) = -999.0d0
            end do
        end if

        if (numcross == 0) then
            dhgrid(1:nlat, j) = np

        else if (numcross == 1) then
            ind1 = int( (90.0d0 - cross_sort(1)) / lat_int) + 1
            dhgrid(1:ind1, j) = np

            if (ind1 == nlat) then
                cycle
            else if (np == 0) then
                dhgrid(ind1+1:nlat, j) = 1
            else
                dhgrid(ind1+1:nlat, j) = 0
            end if

        else
            ind1 = 1
            next = np
            do k = 1, numcross
                ind2 = int( (90.0d0 - cross_sort(k)) / lat_int) + 1
                if (ind2 >= ind1) dhgrid(ind1:ind2, j) = next

                if (next == 0) then
                    next = 1
                else 
                    next = 0
                end if
                ind1 = ind2 + 1

            end do

            if (ind1 <= nlat) dhgrid(ind1:nlat, j) = next

        end if

    end do

end subroutine Curve2Mask
