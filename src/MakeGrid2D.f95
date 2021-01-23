subroutine MakeGrid2D(grid, cilm, lmax, interval, nlat, nlong, norm, csphase, &
                      f, a, north, south, east, west, dealloc, exitstatus)
!------------------------------------------------------------------------------
!
!   Given the Spherical Harmonic coefficients cilm, this subroutine
!   will compute a 2D grid with equal latitude and longitude spacings.
!   Note that since this is NOT done using FFTs, this routine is therefore
!   relatively SLOW! You might want to instead use MakeGridDH.
!
!   The value at grid(1,1) correspons to 90 degrees latitude
!   and 0 degrees longitude, and the longitude spacing goes from 0 to 360,
!   with both points being calculated. If the optional
!   parameters NORTH, SOUTH, EAST and WEST are specified, the upper-left and
!   lower right coordinates of the output grid are (NORTH, WEST) and
!   (SOUTH, EAST), respectively.
!
!   Calling Parameters
!
!       IN
!           cilm        Spherical harmonic coefficients.
!           lmax        Maximum degree of expansions to be performed.
!           interval    Spacing of output grid in DEGREES.
!
!       OUT
!           grid        Gridded expansion of spherical harmonic coefficients.
!           nlat        Number of latitude points for the grid.
!           nlong       Number of longitude points for the grid.
!
!       OPTIONAL (IN)
!           norm        Spherical harmonic normalization:
!                           (1) "geodesy" (default)
!                           (2) Schmidt
!                           (3) unnormalized
!                           (4) orthonormalized
!           csphase     1: Do not include the phase factor of (-1)^m (default).
!                       -1: Apply the phase factor of (-1)^m.
!           f           Flattening of the reference ellipsoid (a-c)/a. If
!                       included, this ellipsoid will be subtracted from
!                       the data.
!           a           Semimajor axis of the reference ellipsoid. If included,
!                       an ellipsoid with these parameters will be subtracted
!                       from the data.
!           north       Maximum latitude to compute, in degrees.
!           south       Minimum latitude to compute, in degrees.
!           east        Maximum longitude to compute, in degrees.
!           west        Minimum latitude to compute, in degrees.
!           dealloc     0 (default): Do not deallocate saved variables in the
!                       Legendre function routines. (1) Dellocate this memory
!                       at the end of the routine.
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
!   Notes:
!       1.  If lmax is greater than the the maximum spherical harmonic
!           degree of the input file, then this file will be ZERO PADDED!
!           (i.e., those degrees after lmax are assumed to be zero).
!
!   Copyright (c) 2005-2019, SHTOOLS
!   All rights reserved.
!
!------------------------------------------------------------------------------
    use SHTOOLS, only:  PlmBar, PlBar, PlmSchmidt, PlSchmidt, PLegendreA, &
                        PLegendre, PlmON, PlON
    use ftypes

    implicit none

    real(dp), intent(in) :: cilm(:,:,:), interval
    real(dp), intent(out) :: grid(:,:)
    integer(int32), intent(in) :: lmax
    integer(int32), intent(out) :: nlat, nlong
    integer(int32), intent(in), optional :: norm, csphase, dealloc
    real(dp), intent(in), optional :: f, a, north, south, east, west
    integer(int32), intent(out), optional :: exitstatus
    integer(int32) :: l, m, j, k, index, l1, m1, lmax_comp, phase, lnorm, &
                      temp, astat(4)
    real(dp) :: pi, latmax, latmin, longmin, longmax, lat, longitude, &
                x, intervalrad, r_ex
    real(dp), allocatable :: pl(:), cosm(:, :), sinm(:, :), cilm2(:,:, :)

    if (present(exitstatus)) exitstatus = 0

    temp = 0

    if (present(north)) temp = temp + 1
    if (present(south)) temp = temp + 1
    if (present(east)) temp = temp + 1
    if (present(west)) temp = temp + 1

    if (temp /= 0 .and. temp /= 4) then
        print*, "Error --- MakeGrid2d"
        print*, "The optional parameters NORTH, SOUTH, EAST, and WEST " // &
                "must all be specified", present(north), present(south), &
                present(east), present(west)
        if (present(exitstatus)) then
            exitstatus = 5
            return
        else
            stop
        end if

    end if

    if (temp == 4) then
        latmax = north
        latmin = south
        longmin = west
        longmax = east

        if (latmax < latmin) then
            print*, "Error --- MakeGrid2d"
            print*, "NORTH must be larger than SOUTH."
            print*, "NORTH = ", latmax
            print*, "SOUTH = ", latmin
            if (present(exitstatus)) then
                exitstatus = 2
                return
            else
                stop
            end if
        end if
        
        if (latmax > 90.0_dp .or. latmin < -90.0_dp) then
            print*, "Error --- MakeGrid2d"
            print*, "NORTH and SOUTH must lie between 90 and -90."
            print*, "NORTH = ", latmax
            print*, "SOUTH = ", latmin
            if (present(exitstatus)) then
                exitstatus = 2
                return
            else
                stop
            end if

        end if

        if (longmin > longmax) longmax = longmax + 360.0_dp

    else
        latmax = 90.0_dp
        latmin = -90.0_dp
        longmin = 0.0_dp
        longmax = 360.d0

    end if

    nlat = int((latmax - latmin) / interval + 1)
    nlong = int((longmax - longmin) / interval + 1)

    if (present(norm)) then
        if (norm > 4 .or. norm < 1) then
            print*, "Error --- MakeGrid2d"
            print*, "Parameter NORM must be 1 (geodesy), " // &
                    "2 (Schmidt), 3 (unnormalized), or 4 (orthonormalized)."
            print*, "Input value is ", norm
            if (present(exitstatus)) then
                exitstatus = 2
                return
            else
                stop
            end if
        end if

        lnorm = norm

    else 
        lnorm = 1

    end if

    if (present(csphase)) then
            if (csphase /= -1 .and. csphase /= 1) then
                print*, "Error --- MakeGrid2D"
                print*, "CSPHASE must be 1 (exclude) or -1 (include)"
                print*, "Input valuse is ", csphase
                if (present(exitstatus)) then
                    exitstatus = 2
                    return
                else
                    stop
                end if

            else
                phase = csphase

            end if
        else
            phase = 1

        end if

    if (size(grid(:,1)) < nlat .or. size(grid(1,:)) < nlong ) then
        print*, "Error --- MakeGrid2D"
        print*, "GRID must be dimensioned ( (LATMAX-LATMIN)/INTERVAL+1, " // &
                "(LONGMAX-LONGMIN)/INTERVAL+1 ) where"
        print*, "INTERVAL = ", interval
        print*, "LATMAX = ", latmax
        print*, "LATMIN = ", latmin
        print*, "LONGMIN = ", longmin
        print*, "LONGMAX = ", longmax
        print*, "Input array is dimensioned ", size(grid(:,1)), size(grid(1,:))
        if (present(exitstatus)) then
            exitstatus = 1
            return
        else
            stop
        end if

    else if (size(cilm(:,1,1)) < 2 .or. size(cilm(1,:,1)) < lmax+1 .or. &
                size(cilm(1,1,:)) < lmax+1) then
        print*, "Error --- MakeGrid2D"
        print*, "CILM must be dimensioned as (2, LMAX+1, LMAX+1) " // &
                "where LMAX is ", lmax
        print*, "Input dimension is ", size(cilm(:,1,1)), size(cilm(1,:,1)), &
                size(cilm(1,1,:))
        if (present(exitstatus)) then
            exitstatus = 1
            return
        else
            stop
        end if

    end if

    if ((present(f) .and. .not. present(a)) .or. (present(a) .and. &
            .not. present(f)) ) then
        print*, "Error --- MakeGrid2D"
        print*, "Both F and A must be present"
        print*, "F ", present(f)
        print*, "A ", present(a)
        if (present(exitstatus)) then
            exitstatus = 5
            return
        else
            stop
        end if

    end if

    allocate (pl(((lmax+1) * (lmax+2)) / 2), stat = astat(1))
    allocate (cosm(lmax+1, int(360.0_dp/interval + 1)), stat = astat(2))
    allocate (sinm(lmax+1, int(360.0_dp/interval + 1)), stat = astat(3))
    allocate (cilm2(2,lmax+1,lmax+1), stat = astat(4))

    if (astat(1) /= 0 .or. astat(2) /= 0 .or. astat(3) /=0 &
            .or. astat(4) /= 0) then
        print*, "Error --- MakeGrid2D"
        print*, "Problem allocating arrays PL, COSM, SINM, and CILM2", &
            astat(1), astat(2), astat(3), astat(4)
        if (present(exitstatus)) then
            exitstatus = 3
            return
        else
            stop
        end if

    end if

    pi = acos(-1.0_dp)
    grid = 0.0_dp

    lmax_comp = min(lmax, size(cilm(1,1,:))-1)

    intervalrad = interval*pi/180.0_dp

    cilm2(1:2,1:lmax+1,1:lmax+1) = cilm(1:2,1:lmax+1,1:lmax+1)

    !--------------------------------------------------------------------------
    !
    !   Precomputing sines and cosines leads to an increase in speed by a
    !   factor of almost 4 with no optimization, and by a factor of about 15
    !   with normal optimizations.
    !
    !--------------------------------------------------------------------------

    do k = 1, nlong
        longitude = longmin * pi / 180.0_dp + dble(k-1) * intervalrad
        sinm(1,k) = 0.0_dp
        cosm(1,k) = 1.0_dp

        if (lmax > 0) then
            sinm(2,k) = sin(longitude)
            cosm(2,k) = cos(longitude)
        end if

        do m = 2, lmax, 1
            sinm(m+1,k) = 2 * sinm(m,k)*cosm(2,k) - sinm(m-1,k)
            cosm(m+1,k) = 2 * cosm(m,k)*cosm(2,k) - cosm(m-1,k)
        end do

    end do

    do j = 1, nlat
        lat = latmax - (j-1)*interval
        x = sin(lat*pi/180.0_dp)

        if (lat == 90.0_dp .or. lat == -90.0_dp) then

            if (present(exitstatus)) then
                select case (lnorm)
                    case (1); call PlBar(pl, lmax_comp, x, &
                                         exitstatus = exitstatus)
                    case (2); call PlSchmidt(pl, lmax_comp, x, &
                                             exitstatus = exitstatus)
                    case (3); call PLegendre(pl, lmax_comp, x, &
                                             exitstatus = exitstatus)
                    case (4); call PlON(pl, lmax_comp, x, &
                                        exitstatus = exitstatus)
                end select
                if (exitstatus /= 0) return

            else
                select case (lnorm)
                    case (1); call PlBar(pl, lmax_comp, x)
                    case (2); call PlSchmidt(pl, lmax_comp, x)
                    case (3); call PLegendre(pl, lmax_comp, x)
                    case (4); call PlON(pl, lmax_comp, x)
                end select
            end if

            do l = lmax_comp, 0, -1
                l1 = l + 1
                grid(j,1) = grid(j,1) + cilm2(1,l1,1) * pl(l1)
            end do

            grid(j, 2:nlong) = grid(j,1)

            if (present(f)) grid(j,1:nlong) = grid(j,1:nlong) - a * (1.0_dp - f)

        else

            if (present(exitstatus)) then
                select case (lnorm)
                    case (1); call PlmBar(pl, lmax_comp, x, csphase = phase, &
                                          exitstatus = exitstatus)
                    case (2); call PlmSchmidt(pl, lmax_comp, x, &
                                              csphase = phase, &
                                              exitstatus = exitstatus)
                    case (3); call PLegendreA(pl, lmax_comp, x, &
                                              csphase = phase, &
                                              exitstatus = exitstatus)
                    case (4); call PlmON(pl, lmax_comp, x, csphase = phase, &
                                         exitstatus = exitstatus)
                end select
                if (exitstatus /= 0) return

            else
                select case (lnorm)
                    case (1); call PlmBar(pl, lmax_comp, x, csphase = phase)
                    case (2); call PlmSchmidt(pl, lmax_comp, x, csphase = phase)
                    case (3); call PLegendreA(pl, lmax_comp, x, csphase = phase)
                    case (4); call PlmON(pl, lmax_comp, x, csphase = phase)
                end select
            end if

            do k = 1, nlong

                ! do m = 0 term first
                m1 = 1
                do l = 0, lmax_comp, 1
                    l1 = l + 1
                    index = ((l+1) * l) / 2 + 1
                    grid(j,k) = grid(j,k) + cilm2(1,l1,1) * cosm(1,k) * &
                                pl(index)
                end do

                do m = 1, lmax_comp, 1
                    m1 = m + 1
                    do l = m, lmax_comp, 1
                        l1 = l + 1
                        index = ((l+1) * l) / 2 + m + 1
                        grid(j,k) = grid(j,k) + (cilm2(1,l1,m1) * cosm(m1,k) +&
                                                 cilm2(2,l1,m1) * sinm(m1,k)) &
                                                 * pl(index)
                    end do
                end do

            end do

            if (present(f)) then
                lat = lat * pi / 180.0_dp
                r_ex = cos(lat)**2 + sin(lat)**2 / (1.0_dp - f)**2
                r_ex = a * sqrt(1.0_dp / r_ex)

                grid(j,1:nlong) = grid(j,1:nlong) - r_ex
            end if

        end if

    end do

    ! deallocate memory
    if (present(dealloc)) then
        if (dealloc == 1) then
            select case (lnorm)
                case (1); call PlmBar(pl, -1, x, csphase = phase)
                case (2); call PlmSchmidt(pl, -1, x, csphase = phase)
                case (4); call PlmON(pl, -1, x, csphase = phase)
            end select
        end if
    end if

    deallocate (pl)
    deallocate (cosm)
    deallocate (sinm)
    deallocate (cilm2)

end subroutine MakeGrid2D
