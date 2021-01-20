subroutine MakeGradientDH(cilm, lmax, theta_grid, phi_grid, n, sampling, &
                          lmax_calc, extend, exitstatus)
!------------------------------------------------------------------------------
!
!   Given the spherical harmonic coefficients CILM of a scalar function defined
!   on the sphere, this subroutine will compute the horizontal gradient of the
!   function, and return 2D Driscol and Healy sampled grids for the theta and
!   phi vector components. The gradient is given by the formula:
!
!       Grad F = 1/r dF/theta theta-hat + 1/(r sin theta) dF/dphi phi-hat
!
!   where theta is colatitude and phi is longitude. The radius r is taken from
!   the degree zero coefficient of the function CILM(1,1,1).
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
!           cilm        Spherical harmonic coefficients of the scalar function.
!           lmax        The maximum spherical harmonic degree of the function,
!                       used to determine the number of samples N.
!
!       IN, OPTIONAL
!           sampling    (1) The output grids are N by N (default).
!                       (2) The output grids are N by 2N.
!           lmax_calc   The maximum spherical harmonic degree to evaluate
!                       the coefficients up to.
!           extend      If 1, return a grid that contains an additional column
!                       and row corresponding to 360 E longitude and 90 S
!                       latitude, respectively.
!
!       OUT
!           theta_grid  Gridded expansaion of the theta component of the
!                       gradient.
!           phi_grid    Gridded expansaion of the phi component of the
!                       gradient.
!           N           Number of samples in latitude. Number of samples in
!                       longitude is N when sampling is 1 (default), and is 2N
!                       when sampling is 2.
!
!       OUT, OPTIONAL
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
!   Copyright (c) 2021, SHTOOLS
!   All rights reserved.
!
!------------------------------------------------------------------------------
    use FFTW3
    use ftypes
    use, intrinsic :: iso_c_binding

    implicit none

    real(dp), intent(in) :: cilm(:,:,:)
    real(dp), intent(out) :: theta_grid(:,:), phi_grid(:,:)
    integer(int32), intent(in) :: lmax
    integer(int32), intent(out) :: n
    integer(int32), intent(in), optional :: sampling, lmax_calc, extend
    integer(int32), intent(out), optional :: exitstatus
    integer(int32) :: l, m, i, l1, m1, lmax_comp, i_eq, i_s, astat(4), nlong, &
                      nlat_out, nlong_out, extend_grid
    real(dp) :: grid(4*lmax+4), pi, theta, scalef, rescalem, u, p, dpl, pmm, &
                pm1, pm2, z, tempr, lat, coeft0, coefts0, coefp0, &
                coefps0, r
    complex(dp) :: coef(2*lmax+3), coeft(2*lmax+3), coefts(2*lmax+3), &
                   coefp(2*lmax+3), coefps(2*lmax+3), tempc
    type(C_PTR) :: plan
    real(dp), save, allocatable :: ff1(:,:), ff2(:,:), sqr(:)
    integer(int8), save, allocatable :: fsymsign(:,:)
    integer(int32), save :: lmax_old = 0

!$OMP   threadprivate(ff1, ff2, sqr, fsymsign, lmax_old)

    if (present(exitstatus)) exitstatus = 0

    n = 2 * lmax + 2

    if (present(sampling)) then
        if (sampling == 1) then
            nlong = n
        else if (sampling == 2) then
            nlong = 2 * n
        else
            print*, "Error --- MakeGradientDH"
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
            print*, "Error --- MakeGradientDH"
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
        print*, "Error --- MakeGradientDH"
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

    if (size(theta_grid(:,1)) < nlat_out .or. size(theta_grid(1,:)) < &
        nlong_out .or. size(phi_grid(:,1)) < nlat_out .or. &
        size(phi_grid(1,:)) < nlong_out) then
        print*, "Error --- MakeGradientDH"
        print*, "THETA_GRID and PHI_GRID " // &
                "must be dimensioned as: ", nlat_out, nlong_out
        print*, "Input dimensions are ", &
                size(theta_grid(:,1)), size(theta_grid(1,:)), &
                size(phi_grid(:,1)), size(phi_grid(1,:))
        if (present(exitstatus)) then
            exitstatus = 1
            return
        else
            stop
        end if

    end if

    pi = acos(-1.0_dp)

    scalef = 1.0e-280_dp

    if (present(lmax_calc)) then
        if (lmax_calc > lmax) then
            print*, "Error --- MakeGradientDH"
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

    r = cilm(1,1,1)  ! set the radius equal to the degree 0 term of the
                     ! function.
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
            print*, "Error --- MakeGradientDH"
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
    theta = pi / 2.0_dp
    z = 0.0_dp
    u = 1.0_dp

    coeft(1:lmax+2) = cmplx(0.0_dp, 0.0_dp, dp)
    coeft0 = 0.0_dp

    coefp(1:lmax+2) = cmplx(0.0_dp, 0.0_dp, dp)
    coefp0 = 0.0_dp

    pm2 = 1.0_dp

    ! derivative of l=0 term is 0, so no need to calculate this

    if (lmax_comp /= 0) then    ! l = 1
        pm1 = 0.0_dp
        dpl = ff1(2,1)
        coeft0 = coeft0 + cilm(1,2,1) * dpl

    end if

    do l = 2, lmax_comp, 1
        l1 = l + 1
        p = - ff2(l1,1) * pm2

        dpl = l * (sqr(2*l+1) / sqr(2*l-1) * pm1)
        coeft0 = coeft0 + cilm(1,l1,1) * dpl

        pm2 = pm1
        pm1 = p

    end do

    pmm = sqr(2) * scalef

    rescalem = 1.0_dp / scalef

    do m = 1, lmax_comp-1, 1
        m1 = m + 1

        pmm = pmm * sqr(2*m+1) / sqr(2*m)
        pm2 = pmm

        coefp(m1) = coefp(m1) + cmplx(cilm(2,m1,m1), cilm(1,m1,m1), dp) &
                                * pm2 * m

        pm1 = 0.0_dp

        dpl = (sqr(2*m+3) * pmm)
        coeft(m1) = coeft(m1) + cmplx(cilm(1,m1+1,m1), &
                                - cilm(2,m1+1,m1), dp) * dpl

        do l = m+2, lmax_comp, 1
            l1 = l + 1
            p = - ff2(l1,m1) * pm2

            coefp(m1) = coefp(m1) + cmplx(cilm(2,l1,m1), cilm(1,l1,m1), dp) &
                                    * p * m

            dpl = ( sqr(2*l+1) * sqr(l-m) * sqr(l+m) / sqr(2*l-1) * pm1)
            coeft(m1) = coeft(m1) + cmplx(cilm(1,l1,m1), - cilm(2,l1,m1), dp) &
                                    * dpl

            pm2 = pm1
            pm1 = p

        end do

        coefp(m1) = coefp(m1) * rescalem
        coeft(m1) = coeft(m1) * rescalem

    end do

    if (lmax_comp /= 0) then

        pmm = pmm * sqr(2*lmax_comp+1) / sqr(2*lmax_comp) * rescalem

        coefp(lmax_comp+1) = coefp(lmax_comp+1) &
                            + cmplx(cilm(2,lmax_comp+1,lmax_comp+1), &
                            cilm(1,lmax_comp+1,lmax_comp+1), dp) &
                            * pmm * lmax_comp

        dpl = -lmax_comp * z * pmm / u**2
        coeft(lmax_comp+1) = coeft(lmax_comp+1) + &
                            cmplx(cilm(1,lmax_comp+1,lmax_comp+1), &
                            - cilm(2,lmax_comp+1,lmax_comp+1), dp) &
                            * dpl

    end if

    coef(1) = cmplx(coeft0, 0.0_dp, dp)
    coef(2:lmax+1) = coeft(2:lmax+1) / 2.0_dp

    if (present(sampling)) then
        if (sampling == 2) then
            coef(lmax+2:2*lmax+3) = cmplx(0.0_dp, 0.0_dp, dp)
        end if
    end if

    call fftw_execute_dft_c2r(plan, coef, grid)
    theta_grid(i_eq,1:nlong) = -sin(theta) * grid(1:nlong) / r

    coef(1) = cmplx(coefp0, 0.0_dp, dp)
    coef(2:lmax+1) = coefp(2:lmax+1) / 2.0_dp

    if (present(sampling)) then
        if (sampling == 2) then
            coef(lmax+2:2*lmax+3) = cmplx(0.0_dp, 0.0_dp, dp)
        end if
    end if

    call fftw_execute_dft_c2r(plan, coef, grid)
    phi_grid(i_eq,1:nlong) = grid(1:nlong) / r / sin(theta)

    do i=1, i_eq-1, 1

        i_s = 2 * i_eq - i

        theta = pi * dble(i-1) / dble(n)
        z = cos(theta)
        u = sqrt( (1.0_dp-z) * (1.0_dp+z) )

        lat = pi / 2.0_dp - theta

        coeft(1:lmax+2) = cmplx(0.0_dp, 0.0_dp, dp)
        coeft0 = 0.0_dp
        coefts(1:lmax+2) = cmplx(0.0_dp, 0.0_dp, dp)
        coefts0 = 0.0_dp

        coefp(1:lmax+2) = cmplx(0.0_dp, 0.0_dp, dp)
        coefp0 = 0.0_dp
        coefps(1:lmax+2) = cmplx(0.0_dp, 0.0_dp, dp)
        coefps0 = 0.0_dp

        pm2 = 1.0_dp

        ! derivative of l=0 term is 0, so no need to calculate this

        if (lmax_comp /= 0) then    ! l = 1
            pm1 = ff1(2,1) * z

            dpl = ff1(2,1)
            tempr = cilm(1,2,1) * dpl
            coeft0 = coeft0 + tempr
            coefts0 = coefts0 + tempr   ! reverse fsymsign

        end if

        do l = 2, lmax_comp, 1
            l1 = l + 1
            p = ff1(l1,1) * z * pm1 - ff2(l1,1) * pm2
            dpl = l * ( sqr(2*l+1) / sqr(2*l-1) * pm1 - z * p ) / u**2
            tempr = cilm(1,l1,1) * dpl
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

            ! (m,m)
            tempc = cmplx(cilm(2,m1,m1), cilm(1,m1,m1), dp) * pm2 * m
            coefp(m1) = coefp(m1) + tempc
            coefps(m1) = coefps(m1) + tempc
            ! fsymsign = 1

            dpl = -m * z * pm2 / u**2
            tempc = cmplx(cilm(1,m1,m1), - cilm(2,m1,m1), dp) * dpl  ! (m,m)
            coeft(m1) = coeft(m1) + tempc
            coefts(m1) = coefts(m1) - tempc ! reverse fsymsign

            pm1 = z * ff1(m1+1,m1) * pm2

            tempc = cmplx(cilm(2,m1+1,m1), cilm(1,m1+1,m1), dp) * pm1 &
                    * m ! (m+1,m)
            coefp(m1) = coefp(m1) + tempc
            coefps(m1) = coefps(m1) - tempc
            ! fsymsign = -1

            ! (m+1,m)
            dpl = ( sqr(2*m+3) * pmm - z * (m+1) * pm1) / u**2
            tempc = cmplx(cilm(1,m1+1,m1), - cilm(2,m1+1,m1), dp) * dpl
            coeft(m1) = coeft(m1) + tempc
            coefts(m1) = coefts(m1) + tempc ! reverse fsymsign

            do l=m+2, lmax_comp, 1
                l1 = l + 1
                p = z * ff1(l1,m1) * pm1 - ff2(l1,m1) * pm2

                tempc = cmplx(cilm(2,l1,m1), cilm(1,l1,m1), dp) * p * m
                coefp(m1) = coefp(m1) + tempc
                coefps(m1) = coefps(m1) + tempc * fsymsign(l1,m1)

                dpl = ( sqr(2*l+1) * sqr(l-m) * sqr(l+m) / sqr(2*l-1) * pm1 &
                       - z * l * p) / u**2
                tempc = cmplx(cilm(1,l1,m1), - cilm(2,l1,m1), dp) * dpl
                coeft(m1) = coeft(m1) + tempc
                ! reverse fsymsign
                coefts(m1) = coefts(m1) - tempc * fsymsign(l1,m1)

                pm2 = pm1
                pm1 = p

            end do

            coeft(m1) = coeft(m1) * rescalem
            coefts(m1) = coefts(m1) * rescalem

            coefp(m1) = coefp(m1) * rescalem
            coefps(m1) = coefps(m1) * rescalem

        end do

        if (lmax_comp /= 0) then

            rescalem = rescalem * u

            pmm = pmm * sqr(2*lmax_comp+1) / sqr(2*lmax_comp) * rescalem

            tempc = cmplx(cilm(2,lmax_comp+1,lmax_comp+1), &
                            cilm(1,lmax_comp+1,lmax_comp+1), dp) * pmm &
                            * lmax_comp
            coefp(lmax_comp+1) = coefp(lmax_comp+1) + tempc
            coefps(lmax_comp+1) = coefps(lmax_comp+1) + tempc
            ! fsymsign = 1

            dpl = -lmax_comp * z * pmm / u**2
            tempc = cmplx(cilm(1,lmax_comp+1,lmax_comp+1), &
                            - cilm(2,lmax_comp+1,lmax_comp+1), dp) * dpl
            coeft(lmax_comp+1) = coeft(lmax_comp+1) + tempc
            ! reverse fsymsign
            coefts(lmax_comp+1) = coefts(lmax_comp+1) - tempc

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
            theta_grid(i,1:nlong) = -sin(theta) * grid(1:nlong) / r

            coef(1) = cmplx(coefp0, 0.0_dp, dp)
            coef(2:lmax+1) = coefp(2:lmax+1) / 2.0_dp

            if (present(sampling)) then
                if (sampling == 2) then
                    coef(lmax+2:2*lmax+3) = cmplx(0.0_dp, 0.0_dp, dp)
                end if
            end if

            call fftw_execute_dft_c2r(plan, coef, grid)
            phi_grid(i,1:nlong) = grid(1:nlong) / r / sin(theta)

        end if

        ! don't compute value for south pole when extend = 0.
        if (.not. (i == 1 .and. extend_grid == 0) ) then

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
                theta_grid(i_s,1:nlong) = -sin(theta) * grid(1:nlong) / r

                coef(1) = cmplx(coefps0, 0.0_dp, dp)
                coef(2:lmax+1) = coefps(2:lmax+1) / 2.0_dp
                if (present(sampling)) then
                    if (sampling == 2) then
                        coef(lmax+2:2*lmax+3) = cmplx(0.0_dp, 0.0_dp, dp)
                    end if
                end if

                call fftw_execute_dft_c2r(plan, coef, grid)
                phi_grid(i_s,1:nlong) = grid(1:nlong) / r / sin(theta)

            end if

        end if

    end do

    call fftw_destroy_plan(plan)

    if (extend_grid == 1) then
        theta_grid(1:nlat_out, nlong_out) = theta_grid(1:nlat_out, 1)
        phi_grid(1:nlat_out, nlong_out) = phi_grid(1:nlat_out, 1)

    end if

end subroutine MakeGradientDH
