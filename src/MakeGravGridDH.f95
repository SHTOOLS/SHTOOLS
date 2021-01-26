subroutine MakeGravGridDH(cilm, lmax, gm, r0, a, f, rad_grid, theta_grid, &
                        phi_grid, total_grid, n, sampling, lmax_calc, omega, &
                        normal_gravity, pot_grid, extend, exitstatus)
!------------------------------------------------------------------------------
!
!   Given the spherical harmonic coefficients CILM of the gravitational
!   potential, this subroutine will compute 2D Driscol and Healy sampled grids
!   of the three vector components of the gravity vector (gravitational force +
!   centrifugal force), and the magnitude of the gravity vector in SI units.
!   The output grids are calculated on a flattened ellipsoid with semi-major
!   axis A and flattening F. In order to calculate the entire gravity vector,
!   it is necessary that the degree-0 term be set equal to 1. Radial gravity is
!   positive when directed UPWARDS.
!
!   If the optional parameter OMEGA is specified, the gravity vector will be
!   calculated in the reference frame of a rotating body and will include the
!   contribution of the centrifugal force. If the parameter NORMAL_GRAVITY is
!   set to 1, the normal gravity predicted for a flattened ellipsoid with A, F,
!   and OMEGA will be removed from the magnitude of the total gravity, yielding
!   the gravity disturbance.
!
!   The gravitational potential is defined as
!
!       V = GM/r Sum_{l=0}^LMAX (R0/r)^l Sum_{m=-l}^l C_{lm} Y_{lm},
!
!   and the gravitational acceleration is
!
!       B = Grad V.
!
!   The output grids contain N samples in latitude and longitude by default,
!   but if the optional parameter SAMPLING is set to 2, the grids will contain
!   N samples in latitude and 2N samples in longitude. When SAMPLING = 1, the
!   output grids contain N samples in latitude from 90 to -90 + interval and N
!   samples in longitude from 0 to 360-2*interval, where N=2*(LMAX+1) and
!   interval=180/N. When SAMPLING = 2, the grids are equally spaced in degrees
!   latitude and longitude with dimension (N x 2N). If the optional parameter
!   EXTEND is set to 1, the output grids will contain an extra column
!   corresponding to 360 E and an extra row corresponding to 90 S, which
!   increases each of the dimensions of the grid by one.
!
!   This routine assumes that the spherical harmonic coefficients are
!   4-pi normalized and exclude the Condon-Shortley phase factor.
!
!   Calling Parameters
!
!       IN
!           cilm        Gravitational spherical harmonic coefficients.
!           lmax        The maximum spherical harmonic degree of the function,
!                       used to determine the number of samples N.
!           GM          Product of the gravitatonal constant and the planet's
!                       mass.
!           r0          Reference radius of the potential coefficients.
!           a           The semimajor axis of the flattened ellipsoid.
!           f           Flattening of the planet.
!
!       IN, OPTIONAL
!           sampling    (1) The output grids are N by N (default).
!                       (2) The output grids are N by 2N.
!           lmax_calc   The maximum spherical harmonic degree to evaluate
!                       the coefficients up to.
!           omega       Angular rotation rate of the planet.
!           normal_gravity
!                       If 1, the magnitude of the normal gravity on the
!                       ellipsoid will be removed from the magnitude of the
!                       total gravity vector. This is the "gravity
!                       disturbance."
!           extend      If 1, return a grid that contains an additional column
!                       and row corresponding to 360 E longitude and 90 S
!                       latitude, respectively.
!
!       OUT
!           rad_grid    Gridded expansion of the radial component of the
!                       gravity field.
!           theta_grid  Gridded expansaion of the theta component of the
!                       gravity field.
!           phi_grid    Gridded expansaion of the phi component of the
!                       gravity field.
!           total_grid  Gridded expansaion of the the magnitude of the
!                       gravity field.
!           N           Number of samples in latitude. Number of samples in
!                       longitude is N when sampling is 1 (default), and is 2N
!                       when sampling is 2.
!
!       OUT, OPTIONAL
!           pot_grid    Gravity potential on the ellipsoid in SI units.
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
!           degree of the input coefficients, then the coefficients will be
!           zero padded.
!       2.  Latitude is geocentric latitude.
!
!   Copyright (c) 2005-2019, SHTOOLS
!   All rights reserved.
!
!------------------------------------------------------------------------------
    use FFTW3
    use SHTOOLS, only: normalgravity
    use ftypes
    use, intrinsic :: iso_c_binding

    implicit none

    real(dp), intent(in) :: cilm(:,:,:), gm, r0, a, f
    real(dp), intent(out) :: rad_grid(:,:), theta_grid(:,:), phi_grid(:,:), &
                             total_grid(:,:)
    real(dp), intent(in), optional :: omega
    real(dp), intent(out), optional :: pot_grid(:,:)
    integer(int32), intent(in) :: lmax
    integer(int32), intent(out) :: n
    integer(int32), intent(in), optional :: sampling, lmax_calc, &
                                            normal_gravity, extend
    integer(int32), intent(out), optional :: exitstatus
    integer(int32) :: l, m, i, l1, m1, lmax_comp, i_eq, i_s, astat(4), nlong, &
                      nlat_out, nlong_out, extend_grid
    real(dp) :: grid(4*lmax+4), pi, theta, scalef, rescalem, u, p, dpl, pmm, &
                pm1, pm2, z, tempr, r_ex, lat, prefactor(lmax), coefr0, &
                coefu0, coefrs0, coeft0, coefts0, coefp0, coefps0, coefus0, b
    complex(dp) :: coef(2*lmax+3), coefr(2*lmax+3), coefrs(2*lmax+3), &
                   coeft(2*lmax+3), coefts(2*lmax+3), coefp(2*lmax+3), &
                   coefps(2*lmax+3), coefu(2*lmax+3), coefus(2*lmax+3), tempc
    type(C_PTR) :: plan
    real(dp), save, allocatable :: ff1(:,:), ff2(:,:), sqr(:)
    integer(int8), save, allocatable :: fsymsign(:,:)
    integer(int32), save :: lmax_old = 0
    logical :: calcu

!$OMP   threadprivate(ff1, ff2, sqr, fsymsign, lmax_old)

    if (present(exitstatus)) exitstatus = 0

    n = 2 * lmax + 2

    if (present(sampling)) then
        if (sampling == 1) then
            nlong = n
        else if (sampling == 2) then
            nlong = 2 * n
        else
            print*, "Error --- MakeGravGridDH"
            print*, "Optional parameter SAMPLING must be 1 (N by N) " // &
                    "or 2 (N by 2N)."
            print*, "Input value is ", sampling
            if (present(exitstatus)) then
                exitstatus = 2
                return
            else
                stop
            end if
        end if
    else
        nlong = n
    end if

    if (present(extend)) then
        if (extend == 0) then
            extend_grid = 0
            nlat_out = n
            nlong_out = nlong
        else if (extend == 1) then
            extend_grid = 1
            nlat_out = n + 1
            nlong_out = nlong + 1
        else
            print*, "Error --- MakeGravGridDH"
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
        nlat_out = n
        nlong_out = nlong
    end if

    if (size(cilm(:,1,1)) < 2) then
        print*, "Error --- MakeGravGridDH"
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

    if (size(rad_grid(:,1)) < nlat_out .or. size(rad_grid(1,:)) < nlong_out &
        .or. size(theta_grid(:,1)) < nlat_out .or. size(theta_grid(1,:)) < &
        nlong_out .or. size(phi_grid(:,1)) < nlat_out .or. size(phi_grid(1,:)) &
        < nlong_out .or. size(total_grid(:,1)) < nlat_out .or. &
        size(total_grid(1,:)) < nlong_out) then
        print*, "Error --- MakeGravGridDH"
        print*, "RAD_GRID, THETA_GRID, PHI_GRID, and TOTAL_GRID " // &
                "must be dimensioned as: ", nlat_out, nlong_out
        print*, "Input dimensions are ", size(rad_grid(:,1)), &
                size(rad_grid(1,:)), size(theta_grid(:,1)), &
                size(theta_grid(1,:)), size(phi_grid(:,1)),  &
                size(phi_grid(1,:)), size(total_grid(:,1)), &
                size(total_grid(1,:))
        if (present(exitstatus)) then
            exitstatus = 1
            return
        else
            stop
        end if

    end if

    if (present(normal_gravity)) then
        if (normal_gravity /= 0 .and. normal_gravity /= 1) then
            print*, "Error --- MakeGravGridDH"
            print*, "NORMAL_GRAVITY must be either 1 (remove normal gravity)"
            print*, "or 0 (do not remove normal gravity)."
            print*, "Input value of NORMAL_GRAVITY is ", normal_gravity
            if (present(exitstatus)) then
                exitstatus = 2
                return
            else
                stop
            end if

        else if (.not. present(omega) .and. normal_gravity == 1) then
            print*, "Error --- MakeGravGridDH"
            print*, "OMEGA must be specified when removing the normal gravity."
            if (present(exitstatus)) then
                exitstatus = 5
                return
            else
                stop
            end if

        end if

    end if

    if (cilm(1,1,1) == 0.0_dp .and. f /= 0.0_dp) then
        print*, "Warning --- MakeGravGridDH"
        print*, "The degree-0 term of the spherical harmonic " // &
                "coefficients is equal to zero."
        print*, "The variation in gravity resulting from variations in " // &
                "radius of the flattened ellipsoid will not be " // &
                "taken into account."
        print*, "C00 = ", cilm(1,1,1)
        print*, "F = ", f

    end if

    if (present(pot_grid)) then
        calcu = .true.

        if (size(pot_grid(:,1)) < nlat_out .or. size(pot_grid(1,:)) &
            < nlong_out) then
            print*, "Error --- MakeGravGridDH"
            print*, "POT_GRID must be dimensioned as: ", nlat_out, nlong_out
            print*, "Input dimensions are ", size(pot_grid(:,1)), &
                    size(pot_grid(1,:))
            if (present(exitstatus)) then
                exitstatus = 1
                return
            else
                stop
            end if

        end if

    else
        calcu = .false.

    end if

    pi = acos(-1.0_dp)

    scalef = 1.0e-280_dp

    if (present(lmax_calc)) then
        if (lmax_calc > lmax) then
            print*, "Error --- MakeGravGridDH"
            print*, "LMAX_CALC must be less than or equal to LMAX."
            print*, "LMAX = ", lmax
            print*, "LMAX_CALC = ", lmax_calc
            if (present(exitstatus)) then
                exitstatus = 2
                return
            else
                stop
            end if

        else
            lmax_comp = min(lmax, size(cilm(1,1,:))-1, size(cilm(1,:,1))-1, &
                            lmax_calc)

        end if

    else
        lmax_comp = min(lmax, size(cilm(1,1,:))-1, size(cilm(1,:,1))-1)

    end if

    !--------------------------------------------------------------------------
    !
    !   Calculate recursion constants used in computing Legendre functions.
    !
    !--------------------------------------------------------------------------
    if (lmax_comp /= lmax_old) then

        if (allocated (sqr)) deallocate (sqr)
        if (allocated (ff1)) deallocate (ff1)
        if (allocated (ff2)) deallocate (ff2)
        if (allocated (fsymsign)) deallocate (fsymsign)

        allocate (sqr(2*lmax_comp+1), stat=astat(1))
        allocate (ff1(lmax_comp+1,lmax_comp+1), stat=astat(2))
        allocate (ff2(lmax_comp+1,lmax_comp+1), stat=astat(3))
        allocate (fsymsign(lmax_comp+1,lmax_comp+1), stat=astat(4))

        if (sum(astat(1:4)) /= 0) then
            print*, "Error --- MakeGravGridDH"
            print*, "Problem allocating arrays SQR, FF1, FF2, or FSYMSIGN", &
                    astat(1), astat(2), astat(3), astat(4)
            if (present(exitstatus)) then
                exitstatus = 3
                return
            else
                stop
            end if

        end if

        !----------------------------------------------------------------------
        !
        !   Calculate signs used for symmetry of Legendre functions about
        !   equator. For the first derivative in theta, these signs are
        !   reversed.
        !
        !----------------------------------------------------------------------
        do l = 0, lmax_comp, 1
            do m = 0, l, 1
                if (mod(l-m,2) == 0) then
                    fsymsign(l+1,m+1) = 1

                else
                    fsymsign(l+1,m+1) = -1

                end if

            end do

        end do

        !----------------------------------------------------------------------
        !
        !   Precompute square roots of integers that are used several times.
        !
        !----------------------------------------------------------------------
        do l = 1, 2*lmax_comp+1
            sqr(l) = sqrt(dble(l))
        end do

        !----------------------------------------------------------------------
        !
        !   Precompute multiplicative factors used in recursion relationships
        !       P(l,m) = x*f1(l,m)*P(l-1,m) - P(l-2,m)*f2(l,m)
        !       k = l*(l+1)/2 + m + 1
        !   Note that prefactors are not used for the case when m=l as a
        !   different recursion is used. Furthermore, for m=l-1, Plmbar(l-2,m)
        !   is assumed to be zero.
        !
        !----------------------------------------------------------------------
        if (lmax_comp /= 0) then
            ff1(2,1) = sqr(3)
            ff2(2,1) = 0.0_dp
        end if

        do l = 2, lmax_comp, 1
            ff1(l+1,1) = sqr(2*l-1) * sqr(2*l+1) / dble(l)
            ff2(l+1,1) = dble(l-1) * sqr(2*l+1) / sqr(2*l-3) / dble(l)

            do m = 1, l-2, 1
                ff1(l+1,m+1) = sqr(2*l+1) * sqr(2*l-1) / sqr(l+m) / sqr(l-m)
                ff2(l+1,m+1) = sqr(2*l+1) * sqr(l-m-1) * sqr(l+m-1) &
                               / sqr(2*l-3) / sqr(l+m) / sqr(l-m)
            end do

            ff1(l+1,l) = sqr(2*l+1) * sqr(2*l-1) / sqr(l+m) / sqr(l-m)
            ff2(l+1,l) = 0.0_dp

        end do

        lmax_old = lmax_comp

    end if

    !--------------------------------------------------------------------------
    !
    !   Create generic plan for grid.
    !
    !--------------------------------------------------------------------------
    plan = fftw_plan_dft_c2r_1d(nlong, coef(1:nlong/2+1), grid(1:nlong), &
                                FFTW_MEASURE)

    !--------------------------------------------------------------------------
    !
    !   Determine Clms one l at a time by intergrating over latitude.
    !
    !--------------------------------------------------------------------------
    i_eq = n / 2 + 1  ! Index correspondong to zero latitude

    ! First do equator
    r_ex = a
    theta = pi / 2.0_dp
    z = 0.0_dp
    u = 1.0_dp

    coefr(1:lmax+2) = cmplx(0.0_dp, 0.0_dp, dp)
    coefr0 = 0.0_dp

    if (calcu) then
        coefu(1:lmax+2) = cmplx(0.0_dp, 0.0_dp, dp)
        coefu0 = 0.0_dp
    end if

    coeft(1:lmax+2) = cmplx(0.0_dp, 0.0_dp, dp)
    coeft0 = 0.0_dp

    coefp(1:lmax+2) = cmplx(0.0_dp, 0.0_dp, dp)
    coefp0 = 0.0_dp

    pm2 = 1.0_dp

    coefr0 = coefr0 - cilm(1,1,1) * pm2

    if (calcu) coefu0 = coefu0 + cilm(1,1,1) * pm2

    ! derivative of l=0 term is 0, so no need to calculate this

    if (lmax_comp /= 0) then    ! l = 1
        prefactor(1) = r0 / r_ex

        do l = 2, lmax_comp,1
            prefactor(l) = prefactor(l-1) * r0 / r_ex
        end do

        pm1 = 0.0_dp

        if (calcu) coefu0 = coefu0 + cilm(1,2,1) * pm1 * prefactor(1)

        dpl = ff1(2,1)
        coeft0 = coeft0 + cilm(1,2,1) * dpl * prefactor(1)

    end if

    do l = 2, lmax_comp, 1
        l1 = l + 1
        p = - ff2(l1,1) * pm2
        coefr0 = coefr0 + cilm(1,l1,1) * p * (-l1) * prefactor(l)

        if (calcu) coefu0 = coefu0 + cilm(1,l1,1) * p * prefactor(l)

        dpl = l * (sqr(2*l+1) / sqr(2*l-1) * pm1)
        coeft0 = coeft0 + cilm(1,l1,1) * dpl * prefactor(l)

        pm2 = pm1
        pm1 = p

    end do

    pmm = sqr(2) * scalef

    rescalem = 1.0_dp / scalef

    do m = 1, lmax_comp-1, 1
        m1 = m + 1

        pmm = pmm * sqr(2*m+1) / sqr(2*m)
        pm2 = pmm

        coefr(m1) = coefr(m1) + cmplx(cilm(1,m1,m1), - cilm(2,m1,m1), dp) &
                                * pm2 * (-m-1) * prefactor(m)

        if (calcu) coefu(m1) = coefu(m1) + cmplx(cilm(1,m1,m1), &
                                - cilm(2,m1,m1), dp) * pm2 * prefactor(m)

        coefp(m1) = coefp(m1) + cmplx(cilm(2,m1,m1), cilm(1,m1,m1), dp) &
                                * pm2 * prefactor(m) * m

        pm1 = 0.0_dp

        dpl = (sqr(2*m+3) * pmm)
        coeft(m1) = coeft(m1) + cmplx(cilm(1,m1+1,m1), &
                                - cilm(2,m1+1,m1), dp) * dpl * prefactor(m+1)

        do l = m+2, lmax_comp, 1
            l1 = l + 1
            p = - ff2(l1,m1) * pm2
            coefr(m1) = coefr(m1) + cmplx(cilm(1,l1,m1), - cilm(2,l1,m1), dp) &
                                    * p * (-l1) * prefactor(l)

            if (calcu) coefu(m1) = coefu(m1) + cmplx(cilm(1,l1,m1), &
                                    - cilm(2,l1,m1), dp) * p * prefactor(l)

            coefp(m1) = coefp(m1) + cmplx(cilm(2,l1,m1), cilm(1,l1,m1), dp) &
                                    * p * prefactor(l) * m

            dpl = ( sqr(2*l+1) * sqr(l-m) * sqr(l+m) / sqr(2*l-1) * pm1)
            coeft(m1) = coeft(m1) + cmplx(cilm(1,l1,m1), - cilm(2,l1,m1), dp) &
                                    * dpl * prefactor(l)

            pm2 = pm1
            pm1 = p

        end do

        coefr(m1) = coefr(m1) * rescalem
        if (calcu) coefu(m1) = coefu(m1) * rescalem
        coefp(m1) = coefp(m1) * rescalem
        coeft(m1) = coeft(m1) * rescalem

    end do

    if (lmax_comp /= 0) then

        pmm = pmm * sqr(2*lmax_comp+1) / sqr(2*lmax_comp) * rescalem
        coefr(lmax_comp+1) = coefr(lmax_comp+1) + &
                            cmplx(cilm(1,lmax_comp+1,lmax_comp+1), &
                            - cilm(2,lmax_comp+1,lmax_comp+1), dp) &
                            * pmm * (-lmax_comp-1) * prefactor(lmax_comp)

        if (calcu) coefu(lmax_comp+1) = coefu(lmax_comp+1) + &
                            cmplx(cilm(1,lmax_comp+1,lmax_comp+1), &
                            - cilm(2,lmax_comp+1,lmax_comp+1), dp) &
                            * pmm * prefactor(lmax_comp)

        coefp(lmax_comp+1) = coefp(lmax_comp+1) &
                            + cmplx(cilm(2,lmax_comp+1,lmax_comp+1), &
                            cilm(1,lmax_comp+1,lmax_comp+1), dp) &
                            * pmm * prefactor(lmax_comp) * lmax_comp

        dpl = -lmax_comp * z * pmm / u**2
        coeft(lmax_comp+1) = coeft(lmax_comp+1) + &
                            cmplx(cilm(1,lmax_comp+1,lmax_comp+1), &
                            - cilm(2,lmax_comp+1,lmax_comp+1), dp) &
                            * dpl * prefactor(lmax_comp)

    end if

    coef(1) = cmplx(coefr0, 0.0_dp, dp)
    coef(2:lmax+1) = coefr(2:lmax+1) / 2.0_dp

    if (present(sampling)) then
        if (sampling == 2) then
            coef(lmax+2:2*lmax+3) = cmplx(0.0_dp, 0.0_dp, dp)
        end if
    end if

    call fftw_execute_dft_c2r(plan, coef, grid)
    rad_grid(i_eq,1:nlong) = grid(1:nlong) * gm / r_ex**2

    if (calcu) then
        coef(1) = cmplx(coefu0, 0.0_dp, dp)
        coef(2:lmax+1) = coefu(2:lmax+1) / 2.0_dp

        if (present(sampling)) then
            if (sampling == 2) then
                coef(lmax+2:2*lmax+3) = cmplx(0.0_dp, 0.0_dp, dp)
            end if
        end if

        call fftw_execute_dft_c2r(plan, coef, grid)
        pot_grid(i_eq,1:nlong) = grid(1:nlong) * gm / r_ex

    end if

    coef(1) = cmplx(coeft0, 0.0_dp, dp)
    coef(2:lmax+1) = coeft(2:lmax+1) / 2.0_dp

    if (present(sampling)) then
        if (sampling == 2) then
            coef(lmax+2:2*lmax+3) = cmplx(0.0_dp, 0.0_dp, dp)
        end if
    end if

    call fftw_execute_dft_c2r(plan, coef, grid)
    theta_grid(i_eq,1:nlong) = -sin(theta) * grid(1:nlong) * gm / r_ex**2

    coef(1) = cmplx(coefp0, 0.0_dp, dp)
    coef(2:lmax+1) = coefp(2:lmax+1) / 2.0_dp

    if (present(sampling)) then
        if (sampling == 2) then
            coef(lmax+2:2*lmax+3) = cmplx(0.0_dp, 0.0_dp, dp)
        end if
    end if

    call fftw_execute_dft_c2r(plan, coef, grid)
    phi_grid(i_eq,1:nlong) = grid(1:nlong) * (gm / r_ex**2) / sin(theta)

    do i=1, i_eq-1, 1

        i_s = 2 * i_eq - i

        theta = pi * dble(i-1) / dble(n)
        z = cos(theta)
        u = sqrt( (1.0_dp-z) * (1.0_dp+z) )

        lat = pi / 2.0_dp - theta

        if (i == 1) then      ! Reference ellipsoid radius
            r_ex = a * (1.0_dp - f)

        else
            r_ex = cos(lat)**2 + sin(lat)**2 / (1.0_dp - f)**2
            r_ex = a * sqrt(1.0_dp / r_ex)

        end if

        coefr(1:lmax+2) = cmplx(0.0_dp, 0.0_dp, dp)
        coefr0 = 0.0_dp
        coefrs(1:lmax+2) = cmplx(0.0_dp, 0.0_dp, dp)
        coefrs0 = 0.0_dp

        coeft(1:lmax+2) = cmplx(0.0_dp, 0.0_dp, dp)
        coeft0 = 0.0_dp
        coefts(1:lmax+2) = cmplx(0.0_dp, 0.0_dp, dp)
        coefts0 = 0.0_dp

        coefp(1:lmax+2) = cmplx(0.0_dp, 0.0_dp, dp)
        coefp0 = 0.0_dp
        coefps(1:lmax+2) = cmplx(0.0_dp, 0.0_dp, dp)
        coefps0 = 0.0_dp

        if (calcu) then
            coefu(1:lmax+2) = cmplx(0.0_dp, 0.0_dp, dp)
            coefu0 = 0.0_dp
            coefus(1:lmax+2) = cmplx(0.0_dp, 0.0_dp, dp)
            coefus0 = 0.0_dp
        end if

        pm2 = 1.0_dp

        tempr = -cilm(1,1,1) * pm2 ! l = 0 
        coefr0 = coefr0 + tempr
        coefrs0 = coefrs0 + tempr   ! fsymsign is always 1 for l=m=0

        if (calcu) then
            tempr = cilm(1,1,1) * pm2  ! l = 0 
            coefu0 = coefu0 + tempr
            coefus0 = coefus0 + tempr   ! fsymsign is always 1 for l=m=0
        end if

        ! derivative of l=0 term is 0, so no need to calculate this

        if (lmax_comp /= 0) then    ! l = 1
            prefactor(1) = r0 / r_ex

            do l = 2, lmax_comp, 1
                prefactor(l) = prefactor(l-1) * r0 / r_ex
            end do

            pm1 = ff1(2,1) * z 
            ! -2 = (l+1) prefactor
            tempr = cilm(1,2,1) * pm1 * (-2) * prefactor(1)
            coefr0 = coefr0 + tempr
            coefrs0 = coefrs0 - tempr   ! fsymsign = -1

            if (calcu) then
                tempr = cilm(1,2,1) * pm1 * prefactor(1)
                coefu0 = coefu0 + tempr
                coefus0 = coefus0 - tempr   ! fsymsign = -1
            end if

            dpl = ff1(2,1)
            tempr = cilm(1,2,1) * dpl * prefactor(1)
            coeft0 = coeft0 + tempr
            coefts0 = coefts0 + tempr   ! reverse fsymsign

        end if

        do l = 2, lmax_comp, 1
            l1 = l + 1
            p = ff1(l1,1) * z * pm1 - ff2(l1,1) * pm2
            tempr = cilm(1,l1,1) * p * (-l1) * prefactor(l)
            coefr0 = coefr0 + tempr
            coefrs0 = coefrs0 + tempr * fsymsign(l1,1)

            if (calcu) then
                tempr = cilm(1,l1,1) * p * prefactor(l)
                coefu0 = coefu0 + tempr
                coefus0 = coefus0 + tempr * fsymsign(l1,1)
            end if

            dpl = l * ( sqr(2*l+1) / sqr(2*l-1) * pm1 - z * p ) / u**2
            tempr = cilm(1,l1,1) * dpl * prefactor(l)
            coeft0 = coeft0 + tempr
            coefts0 = coefts0 - tempr * fsymsign(l1,1)  ! reverse fsymsign

            pm2 = pm1
            pm1 = p

        end do

        pmm = sqr(2) * scalef

        rescalem = 1.0_dp / scalef

        do m = 1, lmax_comp-1, 1

            m1 = m + 1
            rescalem = rescalem * u

            pmm = pmm * sqr(2*m+1) / sqr(2*m)
            pm2 = pmm

            tempc = cmplx(cilm(1,m1,m1), - cilm(2,m1,m1), dp) * pm2 * (-m-1) &
                    * prefactor(m)    ! (m,m)
            coefr(m1) = coefr(m1) + tempc
            coefrs(m1) = coefrs(m1) + tempc
            ! fsymsign = 1

            if (calcu) then
                tempc = cmplx(cilm(1,m1,m1), - cilm(2,m1,m1), dp) * pm2 &
                        * prefactor(m) ! (m,m)
                coefu(m1) = coefu(m1) + tempc
                coefus(m1) = coefus(m1) + tempc
                ! fsymsign = 1
            end if

            tempc = cmplx(cilm(2,m1,m1), cilm(1,m1,m1), dp) * pm2 &
                    * prefactor(m) * m ! (m,m)
            coefp(m1) = coefp(m1) + tempc
            coefps(m1) = coefps(m1) + tempc
            ! fsymsign = 1

            dpl = -m * z * pm2 / u**2
            tempc = cmplx(cilm(1,m1,m1), - cilm(2,m1,m1), dp) * dpl &
                    * prefactor(m)  ! (m,m)
            coeft(m1) = coeft(m1) + tempc
            coefts(m1) = coefts(m1) - tempc ! reverse fsymsign

            pm1 = z * ff1(m1+1,m1) * pm2
            tempc = cmplx(cilm(1,m1+1,m1), - cilm(2,m1+1,m1), dp) * pm1 &
                    * (-m-2) * prefactor(m+1)  ! (m+1,m)
            coefr(m1) = coefr(m1) + tempc
            coefrs(m1) = coefrs(m1) - tempc
            ! fsymsign = -1

            if (calcu) then
                tempc = cmplx(cilm(1,m1+1,m1), - cilm(2,m1+1,m1), dp) * pm1 &
                        * prefactor(m+1)   ! (m+1,m)
                coefu(m1) = coefu(m1) + tempc
                coefus(m1) = coefus(m1) - tempc
            end if

            tempc = cmplx(cilm(2,m1+1,m1), cilm(1,m1+1,m1), dp) * pm1 &
                    * prefactor(m+1) * m ! (m+1,m)
            coefp(m1) = coefp(m1) + tempc
            coefps(m1) = coefps(m1) - tempc
            ! fsymsign = -1

            dpl = ( sqr(2*m+3) * pmm - z * (m+1) * pm1) / u**2
            tempc = cmplx(cilm(1,m1+1,m1), - cilm(2,m1+1,m1), dp) * dpl &
                    * prefactor(m+1)    ! (m+1,m)
            coeft(m1) = coeft(m1) + tempc
            coefts(m1) = coefts(m1) + tempc ! reverse fsymsign

            do l=m+2, lmax_comp, 1
                l1 = l + 1
                p = z * ff1(l1,m1) * pm1 - ff2(l1,m1) * pm2
                tempc = cmplx(cilm(1,l1,m1), - cilm(2,l1,m1), dp) * p &
                        * (-l1) * prefactor(l)
                coefr(m1) = coefr(m1) + tempc
                coefrs(m1) = coefrs(m1) + tempc * fsymsign(l1,m1)

                if (calcu) then
                    tempc = cmplx(cilm(1,l1,m1), - cilm(2,l1,m1), dp) * p &
                            * prefactor(l)
                    coefu(m1) = coefu(m1) + tempc
                    coefus(m1) = coefus(m1) + tempc * fsymsign(l1,m1)
                end if

                tempc = cmplx(cilm(2,l1,m1), cilm(1,l1,m1), dp) * p &
                        * prefactor(l) * m
                coefp(m1) = coefp(m1) + tempc
                coefps(m1) = coefps(m1) + tempc * fsymsign(l1,m1)

                dpl = ( sqr(2*l+1) * sqr(l-m) * sqr(l+m) / sqr(2*l-1) * pm1 &
                       - z * l * p) / u**2
                tempc = cmplx(cilm(1,l1,m1), - cilm(2,l1,m1), dp) * dpl &
                        * prefactor(l)
                coeft(m1) = coeft(m1) + tempc
                ! reverse fsymsign
                coefts(m1) = coefts(m1) - tempc * fsymsign(l1,m1)

                pm2 = pm1
                pm1 = p

            end do

            coefr(m1) = coefr(m1) * rescalem
            coefrs(m1) = coefrs(m1) * rescalem

            if (calcu) then
                coefu(m1) = coefu(m1) * rescalem
                coefus(m1) = coefus(m1) * rescalem
            end if

            coeft(m1) = coeft(m1) * rescalem
            coefts(m1) = coefts(m1) * rescalem

            coefp(m1) = coefp(m1) * rescalem
            coefps(m1) = coefps(m1) * rescalem

        end do

        if (lmax_comp /= 0) then

            rescalem = rescalem * u

            pmm = pmm * sqr(2*lmax_comp+1) / sqr(2*lmax_comp) * rescalem
            tempc = cmplx(cilm(1,lmax_comp+1,lmax_comp+1), &
                            - cilm(2,lmax_comp+1,lmax_comp+1), dp) * pmm &
                            * (-lmax_comp-1) * prefactor(lmax_comp)
            coefr(lmax_comp+1) = coefr(lmax_comp+1) + tempc
            coefrs(lmax_comp+1) = coefrs(lmax_comp+1) + tempc
            ! fsymsign = 1

            if (calcu) then
                tempc = cmplx(cilm(1,lmax_comp+1,lmax_comp+1), &
                                - cilm(2,lmax_comp+1,lmax_comp+1), dp) * pmm &
                                * prefactor(lmax_comp)
                coefu(lmax_comp+1) = coefu(lmax_comp+1) + tempc
                coefus(lmax_comp+1) = coefus(lmax_comp+1) + tempc
            end if

            tempc = cmplx(cilm(2,lmax_comp+1,lmax_comp+1), &
                            cilm(1,lmax_comp+1,lmax_comp+1), dp) * pmm &
                            * prefactor(lmax_comp) * lmax_comp
            coefp(lmax_comp+1) = coefp(lmax_comp+1) + tempc
            coefps(lmax_comp+1) = coefps(lmax_comp+1) + tempc
            ! fsymsign = 1

            dpl = -lmax_comp * z * pmm / u**2
            tempc = cmplx(cilm(1,lmax_comp+1,lmax_comp+1), &
                            - cilm(2,lmax_comp+1,lmax_comp+1), dp) * dpl &
                            * prefactor(lmax_comp)
            coeft(lmax_comp+1) = coeft(lmax_comp+1) + tempc
            ! reverse fsymsign
            coefts(lmax_comp+1) = coefts(lmax_comp+1) - tempc

        end if

        coef(1) = cmplx(coefr0, 0.0_dp, dp)
        coef(2:lmax+1) = coefr(2:lmax+1) / 2.0_dp

        if (present(sampling)) then
            if (sampling == 2) then
                coef(lmax+2:2*lmax+3) = cmplx(0.0_dp, 0.0_dp, dp)
            end if
        end if

        call fftw_execute_dft_c2r(plan, coef, grid)
        rad_grid(i,1:nlong) = grid(1:nlong) * gm / r_ex**2

        if (calcu) then
            coef(1) = cmplx(coefu0, 0.0_dp, dp)
            coef(2:lmax+1) = coefu(2:lmax+1) / 2.0_dp

            if (present(sampling)) then
                if (sampling == 2) then
                    coef(lmax+2:2*lmax+3) = cmplx(0.0_dp, 0.0_dp, dp)
                end if
            end if

            call fftw_execute_dft_c2r(plan, coef, grid)
            pot_grid(i,1:nlong) = grid(1:nlong) * gm / r_ex

        end if

        if (i == 1) then
            ! These two derivatives are undefined at the pole
            theta_grid(1,1:nlong) = 0.0_dp
            phi_grid(1,1:nlong) = 0.0_dp

        else
            coef(1) = cmplx(coeft0, 0.0_dp, dp)
            coef(2:lmax+1) = coeft(2:lmax+1) / 2.0_dp

            if (present(sampling)) then
                if (sampling == 2) then
                    coef(lmax+2:2*lmax+3) = cmplx(0.0_dp, 0.0_dp, dp)
                end if
            end if

            call fftw_execute_dft_c2r(plan, coef, grid)
            theta_grid(i,1:nlong) = -sin(theta) * grid(1:nlong) * gm / r_ex**2

            coef(1) = cmplx(coefp0, 0.0_dp, dp)
            coef(2:lmax+1) = coefp(2:lmax+1) / 2.0_dp

            if (present(sampling)) then
                if (sampling == 2) then
                    coef(lmax+2:2*lmax+3) = cmplx(0.0_dp, 0.0_dp, dp)
                end if
            end if

            call fftw_execute_dft_c2r(plan, coef, grid)
            phi_grid(i,1:nlong) = grid(1:nlong) * (gm/r_ex**2) / sin(theta)

        end if

        ! don't compute value for south pole when extend = 0.
        if (.not. (i == 1 .and. extend_grid == 0) ) then
            coef(1) = cmplx(coefrs0, 0.0_dp, dp)
            coef(2:lmax+1) = coefrs(2:lmax+1) / 2.0_dp

            if (present(sampling)) then
                if (sampling == 2) then
                    coef(lmax+2:2*lmax+3) = cmplx(0.0_dp, 0.0_dp, dp)
                end if
            end if

            call fftw_execute_dft_c2r(plan, coef, grid)
            rad_grid(i_s,1:nlong) = grid(1:nlong) * gm / r_ex**2

            if (calcu) then
                coef(1) = cmplx(coefus0, 0.0_dp, dp)
                coef(2:lmax+1) = coefus(2:lmax+1) / 2.0_dp

                if (present(sampling)) then
                    if (sampling == 2) then
                        coef(lmax+2:2*lmax+3) = cmplx(0.0_dp, 0.0_dp, dp)
                    end if
                end if

                call fftw_execute_dft_c2r(plan, coef, grid)
                pot_grid(i_s,1:nlong) = grid(1:nlong) * gm / r_ex
            end if

            if (i == 1) then
                ! These two derivatives are undefined at the pole
                theta_grid(i_s,1:nlong) = 0.0_dp
                phi_grid(i_s,1:nlong) = 0.0_dp

            else
                coef(1) = cmplx(coefts0, 0.0_dp, dp)
                coef(2:lmax+1) = coefts(2:lmax+1) / 2.0_dp

                if (present(sampling)) then
                    if (sampling == 2) then
                        coef(lmax+2:2*lmax+3) = cmplx(0.0_dp, 0.0_dp, dp)
                    end if
                end if

                call fftw_execute_dft_c2r(plan, coef, grid)
                theta_grid(i_s,1:nlong) = -sin(theta) * grid(1:nlong) &
                                          * gm / r_ex**2

                coef(1) = cmplx(coefps0, 0.0_dp, dp)
                coef(2:lmax+1) = coefps(2:lmax+1) / 2.0_dp
                if (present(sampling)) then
                    if (sampling == 2) then
                        coef(lmax+2:2*lmax+3) = cmplx(0.0_dp, 0.0_dp, dp)
                    end if
                end if

                call fftw_execute_dft_c2r(plan, coef, grid)
                phi_grid(i_s,1:nlong) = grid(1:nlong) * (gm/r_ex**2) / &
                                        sin(theta)

            end if

        end if

    end do

    call fftw_destroy_plan(plan)

    !--------------------------------------------------------------------------
    !
    !   Add rotational effects
    !
    !--------------------------------------------------------------------------

    if (present(omega)) then

        do i = 1, nlat_out

            theta = pi * dble(i-1) / dble(n)
            lat = (pi / 2.0_dp - theta)

            r_ex = (1.0_dp + tan(lat)**2) / &
                   (1.0_dp  + tan(lat)**2 / (1.0_dp - f)**2)
            r_ex = a * sqrt(r_ex)

            rad_grid(i,1:nlong) = rad_grid(i,1:nlong)  &
                                  + r_ex * ( sin(theta) * omega )**2

            theta_grid(i,1:nlong) = theta_grid(i,1:nlong) &
                            + sin(theta) * cos(theta) * r_ex * omega**2

            if (calcu) then
                pot_grid(i,1:nlong) = pot_grid(i,1:nlong) &
                                      + 0.50_dp &
                                      * ( r_ex * sin(theta) * omega )**2
            end if

        end do

    end if

    total_grid(1:nlat_out, 1:nlong) = sqrt(rad_grid(1:nlat_out,1:nlong)**2 &
                                           + phi_grid(1:nlat_out,1:nlong)**2 &
                                           + theta_grid(1:nlat_out,1:nlong)**2)

    ! remove normal gravity from total gravitational acceleration
    if (present(normal_gravity)) then
        if (normal_gravity == 1) then
            b = a * (1.0_dp - f)

            do i = 1, nlat_out
                theta = pi * dble(i-1) / dble(n)
                lat = (pi / 2.0_dp - theta) * 180.0_dp / pi
                total_grid(i,1:nlong) = total_grid(i,1:nlong) - &
                                        NormalGravity(lat, GM, omega, a, b)
            end do

        end if

    end if

    if (extend_grid == 1) then
        rad_grid(1:nlat_out, nlong_out) = rad_grid(1:nlat_out, 1)
        theta_grid(1:nlat_out, nlong_out) = theta_grid(1:nlat_out, 1)
        phi_grid(1:nlat_out, nlong_out) = phi_grid(1:nlat_out, 1)
        total_grid(1:nlat_out, nlong_out) = total_grid(1:nlat_out, 1)

        if (calcu) then
            pot_grid(1:nlat_out, nlong_out) = pot_grid(1:nlat_out, 1)
        end if

    end if

end subroutine MakeGravGridDH
