subroutine MakeGeoidGrid(geoid, cilm, lmax, r0pot, GM, PotRef, omega, r, &
                         gridtype, order, nlat, nlong, interval, lmax_calc,&
                         a, f, extend, exitstatus)
!------------------------------------------------------------------------------
!
!   This subroutine will calculate the height of an equipotential surface with
!   respect to a spherical reference radius R for a body rotating with angular
!   rotation rate OMEGA. This method uses a first-, second-, or third-order
!   Taylor series expansion of the potential, evaluated at radius R. The output
!   gridype is specificied by GRIDTYPE, which can be either Cartesian 2D,
!   Guass-Legendre quadrature, or Driscoll and Healy. If the optional
!   parameters A and F are specified, the geoid heights will be with respect
!   to a flattened ellipsoid with semimajor axis A and flattening F.
!
!   Note that this routine is only strictly valid when the geoid is above the
!   surface. To calculated the height of the geoid when it is below the
!   surface, one would need to know the density structure of the planet.
!   Furthermore, the calculation of the potential (and its derivatives) is
!   only strictly valid when R lies above the maximum radius of the planet.
!
!   Latitude is geocentric spherical latitude.
!
!   Calling Parameters
!
!       IN
!           cilm        Spherical harmonic potential coefficients.
!           lmax        Maximum spherical harmonic degree. For gridtype=1, this
!                       is simply the maximum degree used in evaluating the
!                       spherical harmonic coefficients. For the other grids,
!                       this sets the spacing of the output grid.
!           r0pot       Reference radius of the potential coefficients.
!           GM          GM associated with the potential coefficients.
!           PotRef      Reference potential of the geoid.
!           omega       Angular rotation rate used when calculating the geoid
!                       on in the rotating coordinate system.
!           r           Reference sperical radius to calculate the geoid
!                       heights.
!           gridtype    1 = Gauss-Legendre quadrature grid corresponding
!                           to LMAX.
!                       2 = N by N Driscoll and Healy grid corresponding
!                           to LMAX.
!                       3 = N by 2N Driscoll and Healy grid corresponding
!                           to LMAX.
!                       4 = 2D Cartesian using MakeGrid2D.
!           order       Order of the Taylor expansion. Either 1, 2, or 3.
!
!       OUT
!           geoid       Gridded values of the geoid, in meters, referenced to
!                       the spherical radius R.
!           nlat        Number of latitude points for the grid.
!           nlong       Number of longitude points for the grid.
!
!       OPTIONAL, IN
!           interval    Grid spacing of the output grid in DEGREES. Used only
!                       when GRIDTYPE = 1.
!           lmax_calc   For GRIDTYPE 2, 3, and 4, this specifies the maximum
!                       spherical harmonic degree to evaluate the function to.
!           a           Semimajor axis of the reference ellipsoid.
!           f           Flattening of the reference ellipsoid (a-c)/a.
!           extend      If 1, return a grid that contains an additional column
!                       corresponding to 360 E longitude (for GRIDTYPE 1, 2,
!                       and 3) and an additional row corresponding 90 S
!                       latitude (for GRIDTYPE 1 and 2).
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
    use SHTOOLS, only: MakeGrid2D, MakeGridGLQ, MakeGridDH, SHGLQ
    use ftypes

    implicit none

    real(dp), intent(out) :: geoid(:,:)
    real(dp), intent(in) :: cilm(:,:,:), r0pot, GM, r, PotRef, omega
    integer(int32), intent(in) :: lmax, order, gridtype
    integer(int32), intent(in), optional :: lmax_calc, extend
    integer(int32), intent(out) :: nlat, nlong
    real(dp), intent(in), optional :: interval, a, f
    integer(int32), intent(out), optional :: exitstatus
    real(dp) :: pi, r_ex, lat
    integer(int32) :: l, nlat1, nlong1, lmax_comp, astat, n, i, astat1, &
                      astat2, extend_grid
    real(dp), allocatable :: grida(:,:), gridb(:,:), gridc(:,:), gridd(:,:), &
                             zero(:), w(:), qq(:,:), pp(:,:), uu(:,:), &
                             cilm1(:,:,:), cilm2(:,:,:)

    if (present(exitstatus)) exitstatus = 0

    if ( (present(f) .and. .not. present(a)) .or. (present(a) .and. &
            .not. present(f)) ) then
        print*, "Error --- MakeGeoidGrid"
        print*, "Both F and A must be specified."
        print*, "A, ", present(a)
        print*, "F, ", present(f)
        if (present(exitstatus)) then
            exitstatus = 5
            return
        else
            stop
        end if

    end if

    if (present(lmax_calc)) then
        if (lmax_calc > lmax) then
            print*, "Error --- MakeGeoidGrid"
            print*, "LMAX_CALC must be less than or equal to LMAX."
            print*, "LMAX_CALC = ", lmax_calc
            print*, "LMAX = ", lmax
            if (present(exitstatus)) then
                exitstatus = 2
                return
            else
                stop
            end if

        end if
    end if

    if (gridtype < 1 .or. gridtype > 4) then
        print*, "Error --- MakeGeoidGrid"
        print*, "GRIDTYPE must be either (1) Cartesian 2D, " // &
                "(2) Gauss-Legendre quadrature,"
        print*, "(3) N by N Driscoll and Healy, or " // &
                "(4) N by 2N Driscoll and Healy."
        print*, "Input value = ", gridtype
        if (present(exitstatus)) then
            exitstatus = 2
            return
        else
            stop
        end if

    end if

    if (gridtype == 4 .and. .not. present(interval)) then
        print*, "Error --- MakeGeoidGrid"
        print*, "If GRIDTYPE = 4 (2D Cartesian), the optional " // &
                "parameter INTERVAL must be specified."
        if (present(exitstatus)) then
            exitstatus = 5
            return
        else
            stop
        end if

    end if

    if (order > 3 .or. order < 1) then
        print*, "Error --- MakeGeoidGrid"
        print*, "ORDER must be 1, 2, or 3."
        if (present(exitstatus)) then
            exitstatus = 2
            return
        else
            stop
        end if

    end if

    if (gridtype == 1) then
        nlong = 2 * lmax + 1
        nlat = lmax + 1

    else if (gridtype == 2) then
        nlat = 2 * lmax + 2
        nlong = nlat

    else if (gridtype == 3) then
        nlat = 2 * lmax + 2
        nlong = 2 * nlat

    else if (gridtype == 4) then
        nlat = int(180.0_dp / interval + 1)
        nlong = int(360.0_dp / interval + 1)

    end if

    if (present(extend)) then
        extend_grid = extend
        if (extend == 0) then
            continue
        else if (extend == 1) then
            if (gridtype == 1) then
                nlong = nlong + 1
            else if (gridtype == 2 .or. gridtype == 3) then
                nlong = nlong + 1
                nlat = nlat + 1
            else
                print*, "Error --- MakeGeoidGrid"
                print*, "EXTEND can only be set to 1 when GRIDTYPE " // &
                        "is 2 or 3."
                if (present(exitstatus)) then
                    exitstatus = 2
                    return
                else
                    stop
                end if
            end if

        else
            print*, "Error --- MakeGeoidGrid"
            print*, "EXTEND must be 0 or 1. Input value is ", extend
            if (present(exitstatus)) then
                exitstatus = 2
                return
            else
                stop
            end if

        end if
    end if

    if (size(geoid(:,1)) < nlat .or. size(geoid(1,:)) < nlong) then
        print*, "Error --- MakeGeoidGrid"
        print*, "GEOID must be dimensioned as (180/INTERVAL+1, " // &
                "360/INTERVAL+1) where INTERVAL is ", interval
        print*, "Input array has dimension ", size(geoid(:,1)), size(geoid(1,:))
        if (present(exitstatus)) then
            exitstatus = 1
            return
        else
            stop
        end if

    else if (size(cilm(:,1,1)) < 2) then
        print*, "Error --- MakeGeoidGrid"
        print*, "CILM must be dimensioned as (2, *, *)."
        print*, "Input dimension is ", size(cilm(:,1,1)), size(cilm(1,:,1)), &
                size(cilm(1,1,:))
        if (present(exitstatus)) then
            exitstatus = 1
            return
        else
            stop
        end if

    end if

    if (present(lmax_calc)) then
        lmax_comp = min(lmax, size(cilm(1,1,:))-1, size(cilm(1,:,1))-1, &
                        lmax_calc)

    else
        lmax_comp = min(lmax, size(cilm(1,1,:))-1, size(cilm(1,:,1))-1)

    end if

    if (gridtype == 1) then
        allocate (zero(lmax+1), stat = astat)

        if (astat /= 0) then
            print*, "Error --- MakeGeoidGrid"
            print*, "Problem allocating array ZERO, ", astat
            if (present(exitstatus)) then
                exitstatus = 3
                return
            else
                stop
            end if

        end if

        allocate (w(lmax+1), stat = astat)

        if (astat /= 0) then
            print*, "Error --- MakeGeoidGrid"
            print*, "Problem allocating array W, ", astat
            if (present(exitstatus)) then
                exitstatus = 3
                return
            else
                stop
            end if

        end if

        if (present(exitstatus)) then
            call SHGLQ(lmax, zero, w, norm=1, csphase=1, &
                       exitstatus=exitstatus)
            if (exitstatus /= 0) return
        else
            call SHGLQ(lmax, zero, w, norm=1, csphase=1)
        end if

    end if

    allocate (cilm1(2,lmax+1,lmax+1), stat = astat1)
    allocate (cilm2(2,lmax+1,lmax+1), stat = astat2)
    if (astat1 /= 0 .or. astat2 /= 0) then
        print*, "Error --- MakeGeoidGrid"
        print*, "Problem allocating arrays CILM1 and CILM2", astat1, astat2
        if (present(exitstatus)) then
             exitstatus = 3
             return
        else
             stop
        end if

    end if

    pi = acos(-1.0_dp)

    cilm1(1:2,1:lmax_comp+1,1:lmax_comp+1) = cilm(1:2,1:lmax_comp+1, &
                                                  1:lmax_comp+1)
    cilm1(1,1,1) = 1.0_dp    ! Make sure that the degree-0 term is included

    !--------------------------------------------------------------------------
    !
    !   Solve the equation a + b x + c x**2 + d x**3 = 0
    !   for each grid point.
    !
    !--------------------------------------------------------------------------
    ! Create grid A.
    do l = 0, lmax_comp
        cilm2(1:2,l+1,1:l+1) = (GM / r) * cilm1(1:2,l+1,1:l+1) * (r0pot / r)**l
    end do

    ! add rotation and subtract reference field
    cilm2(1,1,1) = cilm2(1,1,1) + (omega * r)**2 / 3.0_dp - PotRef
    cilm2(1,3,1) = cilm2(1,3,1) - (omega * r)**2 / (3.0_dp * sqrt(5.0_dp))

    allocate (grida(nlat, nlong), stat = astat)
    if (astat /= 0) then
        print*, "Error --- MakeGeoidGrid"
        print*, "Problem allocating array GRIDA", astat
        if (present(exitstatus)) then
            exitstatus = 3
            return
        else
            stop
         end if

    end if

    select case (gridtype)
        case (1)
            if (present(exitstatus)) then
                call MakeGridGLQ(grida, cilm2, lmax, zero=zero, norm=1, &
                                 csphase=1, lmax_calc=lmax_comp, &
                                 extend=extend_grid, exitstatus=exitstatus)
                if (exitstatus /= 0) return
            else
                call MakeGridGLQ(grida, cilm2, lmax, zero=zero, norm=1, &
                                 csphase=1, lmax_calc=lmax_comp, &
                                 extend=extend_grid)
            end if

        case (2)
            if (present(exitstatus)) then
                call MakeGridDH(grida, n, cilm2, lmax, norm=1, csphase=1, &
                                sampling=1, lmax_calc=lmax_comp, &
                                extend=extend_grid, exitstatus=exitstatus)
                if (exitstatus /= 0) return
            else
                call MakeGridDH(grida, n, cilm2, lmax, norm=1, csphase=1, &
                                sampling=1, lmax_calc=lmax_comp, &
                                extend=extend_grid)
            end if

        case (3)
            if (present(exitstatus)) then
                call MakeGridDH(grida, n, cilm2, lmax, norm=1, csphase=1, &
                                sampling=2, lmax_calc=lmax_comp, &
                                extend=extend_grid, exitstatus=exitstatus)
                if (exitstatus /= 0) return
            else
                call MakeGridDH(grida, n, cilm2, lmax, norm=1, csphase=1, &
                                sampling=2, lmax_calc=lmax_comp, &
                                extend=extend_grid)
            end if

        case (4)
            if (present(exitstatus)) then
                call MakeGrid2D(grida, cilm2, lmax_comp, interval, nlat1, &
                                nlong1, norm=1, csphase=1, &
                                exitstatus=exitstatus)
                if (exitstatus /= 0) return
            else
                call MakeGrid2D(grida, cilm2, lmax_comp, interval, nlat1, &
                                nlong1, norm=1, csphase=1)
            end if

    end select

    ! Create Grid B
    do l = 0, lmax_comp
        cilm2(1:2,l+1,1:l+1) = -(GM / r**2) * dble(l+1) &
                               * cilm1(1:2,l+1,1:l+1) * (r0pot / r)**l
    end do

    ! add rotational terms
    cilm2(1,1,1) = cilm2(1,1,1) + 2.0_dp * r * omega**2 / 3.0_dp
    cilm2(1,3,1) = cilm2(1,3,1) - 2.0_dp * r * omega**2 / &
                   (3.0_dp * sqrt(5.0_dp))

    allocate (gridb(nlat, nlong), stat = astat)
    if (astat /= 0) then
        print*, "Error --- MakeGeoidGrid"
        print*, "Problem allocating array GRIDB", astat
        if (present(exitstatus)) then
            exitstatus = 3
            return
        else
            stop
        end if

    end if

    select case (gridtype)
        case (1)
            if (present(exitstatus)) then
                call MakeGridGLQ(gridb, cilm2, lmax, zero=zero, norm=1, &
                                 csphase=1, lmax_calc=lmax_comp, &
                                 extend=extend_grid, exitstatus=exitstatus)
                if (exitstatus /= 0) return
            else
                call MakeGridGLQ(gridb, cilm2, lmax, zero=zero, norm=1, &
                                 csphase=1, lmax_calc=lmax_comp, &
                                 extend=extend_grid)
            end if

        case (2)
            if (present(exitstatus)) then
                call MakeGridDH(gridb, n, cilm2, lmax, norm=1, csphase=1, &
                                sampling=1, lmax_calc=lmax_comp, &
                                extend=extend_grid, exitstatus=exitstatus)
                if (exitstatus /= 0) return
            else
                call MakeGridDH(gridb, n, cilm2, lmax, norm=1, csphase=1, &
                                sampling=1, lmax_calc=lmax_comp, &
                                extend=extend_grid)
            end if

        case (3)
            if (present(exitstatus)) then
                call MakeGridDH(gridb, n, cilm2, lmax, norm=1, csphase=1, &
                                sampling=2, lmax_calc=lmax_comp, &
                                extend=extend_grid, exitstatus=exitstatus)
                if (exitstatus /= 0) return
            else
                call MakeGridDH(gridb, n, cilm2, lmax, norm=1, csphase=1, &
                                sampling=2, lmax_calc=lmax_comp, &
                                extend=extend_grid)
            end if

        case (4)
            if (present(exitstatus)) then
                call MakeGrid2D(gridb, cilm2, lmax_comp, interval, nlat1, &
                                nlong1, norm=1, csphase=1, &
                                exitstatus=exitstatus)
                if (exitstatus /= 0) return
            else
                call MakeGrid2D(gridb, cilm2, lmax_comp, interval, nlat1, &
                                nlong1, norm=1, csphase=1)
            end if

    end select

    ! Create Grid C
    if (order == 2 .or. order == 3) then

        do l=0, lmax_comp
            cilm2(1:2,l+1,1:l+1) = GM / (2.0_dp * r**3) * dble(l+1) &
                                   * dble(l+2) * cilm1(1:2,l+1,1:l+1) &
                                   * (r0pot / r)**l
        end do

        ! add rotational terms
        cilm2(1,1,1) = cilm2(1,1,1) + omega**2 / 3.0_dp
        cilm2(1,3,1) = cilm2(1,3,1) - omega**2 / (3.0_dp * sqrt(5.0_dp))

        allocate (gridc(nlat, nlong), stat = astat)
        if (astat /= 0) then
            print*, "Error --- MakeGeoidGrid"
            print*, "Problem allocating array GRIDC", astat
            if (present(exitstatus)) then
                exitstatus = 3
                return
            else
                stop
            end if

        end if

    select case (gridtype)
        case (1)
            if (present(exitstatus)) then
                call MakeGridGLQ(gridc, cilm2, lmax, zero=zero, norm=1, &
                                 csphase=1, lmax_calc=lmax_comp, &
                                 extend=extend_grid, exitstatus=exitstatus)
                if (exitstatus /= 0) return
            else
                call MakeGridGLQ(gridc, cilm2, lmax, zero=zero, norm=1, &
                                 csphase=1, lmax_calc=lmax_comp, &
                                 extend=extend_grid)
            end if

        case (2)
            if (present(exitstatus)) then
                call MakeGridDH(gridc, n, cilm2, lmax, norm=1, csphase=1, &
                                sampling=1, lmax_calc=lmax_comp, &
                                extend=extend_grid, exitstatus=exitstatus)
                if (exitstatus /= 0) return
            else
                call MakeGridDH(gridc, n, cilm2, lmax, norm=1, csphase=1, &
                                sampling=1, lmax_calc=lmax_comp, &
                                extend=extend_grid)
            end if

        case (3)
            if (present(exitstatus)) then
                call MakeGridDH(gridc, n, cilm2, lmax, norm=1, csphase=1, &
                                sampling=2, lmax_calc=lmax_comp, &
                                extend=extend_grid, exitstatus=exitstatus)
                if (exitstatus /= 0) return
            else
                call MakeGridDH(gridc, n, cilm2, lmax, norm=1, csphase=1, &
                                sampling=2, lmax_calc=lmax_comp, &
                                extend=extend_grid)
            end if

        case (4)
            if (present(exitstatus)) then
                call MakeGrid2D(gridc, cilm2, lmax_comp, interval, nlat1, &
                                nlong1, norm=1, csphase=1, &
                                exitstatus=exitstatus)
                if (exitstatus /= 0) return
            else
                call MakeGrid2D(gridc, cilm2, lmax_comp, interval, nlat1, &
                                nlong1, norm=1, csphase=1)
            end if

    end select

    end if

    ! Create Grid D
    if (order == 3) then

        do l=0, lmax_comp
            cilm2(1:2,l+1,1:l+1) = -GM / (6.0_dp * r**4) * dble(l+1) &
                                   * dble(l+2) * dble(l+3) * &
                                    cilm1(1:2,l+1,1:l+1) * (r0pot / r)**l
        end do

        allocate (gridd(nlat, nlong), stat = astat)
        if (astat /= 0) then
            print*, "Error --- MakeGeoidGrid"
            print*, "Problem allocating array GRIDD", astat
            if (present(exitstatus)) then
                exitstatus = 3
                return
            else
                stop
            end if

        end if

    select case (gridtype)
        case (1)
            if (present(exitstatus)) then
                call MakeGridGLQ(gridd, cilm2, lmax, zero=zero, norm=1, &
                                 csphase=1, lmax_calc=lmax_comp, &
                                 extend=extend_grid, exitstatus=exitstatus)
                if (exitstatus /= 0) return
            else
                call MakeGridGLQ(gridd, cilm2, lmax, zero=zero, norm=1, &
                                 csphase=1, lmax_calc=lmax_comp, &
                                 extend=extend_grid)
            end if

        case (2)
            if (present(exitstatus)) then
                call MakeGridDH(gridd, n, cilm2, lmax, norm=1, csphase=1, &
                                sampling=1, lmax_calc = lmax_comp, &
                                extend=extend_grid, exitstatus=exitstatus)
                if (exitstatus /= 0) return
            else
                call MakeGridDH(gridd, n, cilm2, lmax, norm=1, csphase=1, &
                                sampling=1, lmax_calc=lmax_comp, &
                                extend=extend_grid)
            end if

        case (3)
            if (present(exitstatus)) then
                call MakeGridDH(gridd, n, cilm2, lmax, norm=1, csphase=1, &
                                sampling=2, lmax_calc=lmax_comp, &
                                extend=extend_grid, exitstatus=exitstatus)
                if (exitstatus /= 0) return
            else
                call MakeGridDH(gridd, n, cilm2, lmax, norm=1, csphase=1, &
                                sampling=2, lmax_calc=lmax_comp, &
                                extend=extend_grid)
            end if

        case (4)
            if (present(exitstatus)) then
                call MakeGrid2D(gridd, cilm2, lmax_comp, interval, nlat1, &
                                nlong1, norm=1, csphase=1, &
                                exitstatus=exitstatus)
                if (exitstatus /= 0) return
            else
                call MakeGrid2D(gridd, cilm2, lmax_comp, interval, nlat1, &
                                nlong1, norm=1, csphase=1)
            end if

    end select

    end if

    ! solve equation
    if (order == 1) then
        geoid(1:nlat,1:nlong) = -grida(1:nlat,1:nlong) / gridb(1:nlat,1:nlong)

    else if (order == 2) then
        geoid(1:nlat, 1:nlong) = (-gridb(1:nlat, 1:nlong) - &
                                  sqrt(gridb(1:nlat, 1:nlong)**2 &
                                  - 4.0_dp * gridc(1:nlat,1:nlong) * &
                                  grida(1:nlat, 1:nlong)) ) &
                                  / (2.0_dp * gridc(1:nlat,1:nlong))

    else if (order == 3) then
        allocate (qq(nlat, nlong), stat = astat)
        if (astat /= 0) then
            print*, "Error --- MakeGeoidGrid"
            print*, "Problem allocating array QQ", astat
            if (present(exitstatus)) then
                exitstatus = 3
                return
            else
                stop
            end if

        end if

        allocate (pp(nlat, nlong), stat = astat)
        if (astat /= 0) then
            print*, "Error --- MakeGeoidGrid"
            print*, "Problem allocating array PP", astat
            if (present(exitstatus)) then
                exitstatus = 3
                return
            else
                stop
            end if

        end if

        allocate (uu(nlat, nlong), stat = astat)
        if (astat /= 0) then
            print*, "Error --- MakeGeoidGrid"
            print*, "Problem allocating array UU", astat
            if (present(exitstatus)) then
                exitstatus = 3
                return
            else
                stop
            end if

        end if

        pp(1:nlat,1:nlong) = gridb(1:nlat,1:nlong) / gridd(1:nlat,1:nlong) &
                             - ((gridc(1:nlat,1:nlong) &
                             / gridd(1:nlat,1:nlong))**2) / 3.0_dp
        qq(1:nlat,1:nlong) = grida(1:nlat,1:nlong) / gridd(1:nlat,1:nlong) &
             + 2.0_dp * ((gridc(1:nlat,1:nlong) / gridd(1:nlat,1:nlong))**3) &
             / 27.0_dp - 9.0_dp * gridc(1:nlat,1:nlong) * &
             gridb(1:nlat,1:nlong) / (gridd(1:nlat,1:nlong)**2) / 27.0_dp
        uu(1:nlat,1:nlong) = (qq(1:nlat,1:nlong) / 2.0_dp &
                             + sqrt((qq(1:nlat,1:nlong)**2) / 4.0_dp &
                             + (pp(1:nlat,1:nlong)**3) / 27.0_dp) &
                               )**(1.0_dp / 3.0_dp)
        geoid(1:nlat,1:nlong) = pp(1:nlat,1:nlong) / 3.0_dp &
                                / uu(1:nlat,1:nlong) - uu(1:nlat,1:nlong) &
                                - gridc(1:nlat,1:nlong) &
                                / gridd(1:nlat,1:nlong) / 3.0_dp

        deallocate (qq)
        deallocate (pp)
        deallocate (uu)

    end if

    !--------------------------------------------------------------------------
    !
    !   Reference geoid to a flattened ellipsoid
    !
    !--------------------------------------------------------------------------
    if (present(a) .and. present(f)) then

        do i = 1, nlat
            if (gridtype == 4) then
                lat = 90.0_dp - dble(i-1) * interval

            else if (gridtype == 1) then
                lat = asin(zero(i)) * 180.0_dp / pi

            else if (gridtype == 2 .or. gridtype == 3) then
                lat = 90.0_dp - 180.0_dp * dble(i-1) / dble(2 * lmax + 2)

            end if

            lat = lat * pi / 180.0_dp
            r_ex = cos(lat)**2 + sin(lat)**2 / (1.0_dp - f)**2
            r_ex = a * sqrt(1.0_dp / r_ex)

            geoid(i,1:nlong) = geoid(i,1:nlong) + r - r_ex

        end do

    end if

    !--------------------------------------------------------------------------
    !
    !   Clean up
    !
    !--------------------------------------------------------------------------
    if (gridtype == 2) then
        deallocate(zero)
        deallocate(w)
    end if

    deallocate (grida)
    deallocate (gridb)
    if (order == 2) deallocate (gridc)

    if (order == 3) then
        deallocate (gridc)
        deallocate (gridd)
    end if

    deallocate (cilm1)
    deallocate (cilm2)

end subroutine MakeGeoidGrid
