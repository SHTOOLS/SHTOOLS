subroutine MakeMagGradGridDH(cilm, lmax, r0, a, f, vxx, vyy, vzz, vxy, &
                             vxz, vyz, n, sampling, lmax_calc, extend, &
                             exitstatus)
!------------------------------------------------------------------------------
!
!   Given the magnetic potential spherical harmonic coefficients CILM, this
!   subroutine will compute 2D Driscol and Healy sampled grids of the six
!   components of the magnetic field tensor in a local north-oriented
!   reference frame:
!
!       (Vxx,   Vxy,    Vxz)
!       (Vyx,   Vyy,    Vyz)
!       (Vzx,   Vzy,    Vzz)
!
!   where X points NORTH, Y points WEST, and Z points UPWARD. The magnetic
!   potential is defined as
!
!       V = R0 Sum_{l=0}^LMAX (R0/r)^(l+1) Sum_{m=-l}^l C_{lm} Y_{lm},
!
!   where the Gauss coefficients are in units of nT, and the spherical
!   harmonic functions are Schmidt semi-normalized.
!
!   Laplace's equation implies that Vxx + Vyy + Vzz = 0, and the tensor
!   is symmetric. The components are calculated according to eq. 1 in
!   Petrovskaya and Vershkov (2006, J. Geod, 80, 117-127), which is based on
!   eq. 3.28 in Reed (1973, Ohio State Univ., Dept. Geod. Sci., Rep. 201,
!   Columbus, OH). Note that Reed's equations are in terms of latitude, and
!   that the Y axis points East:
!
!       Vzz = Vrr
!       Vxx = 1/r Vr + 1/r^2 Vtt
!       Vyy = 1/r Vr + 1/r^2 /tan(t) Vt + 1/r^2 /sin(t)^2 Vpp
!       Vxy = 1/r^2 /sin(t) Vtp - cos(t)/sin(t)^2 /r^2 Vp
!       Vxz = 1/r^2 Vt - 1/r Vrt
!       Vyz = 1/r^2 /sin(t) Vp - 1/r /sin(t) Vrp
!
!   where r, t, p stand for radius, theta, and phi, and subscripts on V denote
!   partial derivatives.
!
!   The output grids are in units of nT / m and are cacluated on a flattened
!   ellipsoid with semi-major axis A and flattening F.
!
!   When SAMPLING = 1, the output grids contain N samples in latitude from
!   90 to -90 + interval and N samples in longitude from 0 to 360-2*interval,
!   where N=2*(LMAX+1) and interval=180/N. When SAMPLING = 2, the grids are
!   equally spaced in degrees latitude and longitude with dimension (N x 2N).
!   If the optional parameter EXTEND is set to 1, the output grids will contain
!   an extra column corresponding to 360 E and an extra row corresponding to
!   90 S, which increases each of the dimensions of the grid by one.
!
!   Calling Parameters
!
!       IN
!           cilm        Schmidth seminormalized magnetic potential spherical
!                       harmonic coefficients.
!           lmax        The maximum spherical harmonic degree of the function,
!                       used to determine the number of samples N.
!           r0          Reference radius of potential coefficients.
!           a           The semimajor axis of the flattened ellipsoid.
!           f           Flattening of the planet. (a-c)/a.
!
!       IN, OPTIONAL
!           sampling    (1) The output grids are N by N longitudes (default).
!                       (2) The output grids are N by 2N.
!           lmax_calc   The maximum spherical harmonic degree to evaluate
!                       the coefficients up to.
!           extend      If 1, return a grid that contains an additional column
!                       and row corresponding to 360 E longitude and 90 S
!                       latitude, respectively.
!
!       OUT
!           Vxx         x-x component of the magnetic field tensor.
!           Vyy         y-y component of the magnetic field tensor.
!           Vzz         z-z component of the magnetic field tensor.
!           Vxy         x-y component of the magnetic field tensor.
!           Vxz         x-z component of the magnetic field tensor.
!           Vyz         y-z component of the magnetic field tensor.
!           N           Number of samples in latitude. Number of samples in
!                       longitude is N when sampling is 1 (default), and is
!                       2N when sampling is 2.
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
!           degree of the input coefficients, then the coefficients will be
!           zero padded.
!       2.  Latitude is geocentric latitude.
!
!   Copyright (c) 2005-2019, SHTOOLS
!   All rights reserved.
!
!------------------------------------------------------------------------------
    use FFTW3
    use ftypes
    use, intrinsic :: iso_c_binding

    implicit none

    real(dp), intent(in) :: cilm(:,:,:), r0, a, f
    real(dp), intent(out) :: vxx(:,:), vyy(:,:), vzz(:,:), vxy(:,:), &
                             vxz(:,:), vyz(:,:)
    integer(int32), intent(in) :: lmax
    integer(int32), intent(out) :: n
    integer(int32), intent(in), optional :: sampling, lmax_calc, extend
    integer(int32), intent(out), optional :: exitstatus
    integer(int32) :: l, m, i, l1, m1, lmax_comp, i_eq, i_s, astat(4), nlong, &
                      nlat_out, nlong_out, extend_grid
    real(dp) :: grid(4*lmax+4), pi, theta, scalef, rescalem, u, p, dpl, dpl2, &
                dpl2s, pmm, sint, pm1, pm2, z, tempr, r_ex, lat, &
                prefactor(lmax), coefr0, coefrs0, coeft0, coefts0, coefp0, &
                coefps0, coefrr0, coefrrs0, coefrt0, coefrts0, coefrp0, &
                coefrps0, coeftp0, coeftps0, coefpp0, coefpps0, coeftt0, &
                coeftts0
    complex(dp) :: coef(2*lmax+3), coefr(2*lmax+3), coefrs(2*lmax+3), &
                coeft(2*lmax+3), coefts(2*lmax+3), coefp(2*lmax+3), &
                coefps(2*lmax+3), tempc, coefrr(2*lmax+3), coefrrs(2*lmax+3), &
                coefrt(2*lmax+3), coefrts(2*lmax+3), coefrp(2*lmax+3), &
                coefrps(2*lmax+3), coeftp(2*lmax+3), coeftps(2*lmax+3), &
                coefpp(2*lmax+3), coefpps(2*lmax+3), coeftt(2*lmax+3), &
                coeftts(2*lmax+3)
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
            print*, "Error --- MakeMagGradGridDH"
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
            print*, "Error --- MakeMagGradGridDH"
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
        print*, "Error --- MakeMagGradGridDH"
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

    if (size(vxx(:,1)) < nlat_out .or. size(vxx(1,:)) < nlong_out .or. &
            size(vyy(:,1)) < nlat_out .or. size(vyy(1,:)) < nlong_out .or. &
            size(vzz(:,1)) < nlat_out .or. size(vzz(1,:)) < nlong_out .or. &
            size(vxy(:,1)) < nlat_out .or. size(vxy(1,:)) < nlong_out .or. &
            size(vxz(:,1)) < nlat_out .or. size(vxz(1,:)) < nlong_out .or. &
            size(vyz(:,1)) < nlat_out .or. size(vyz(1,:)) < nlong_out) then
        print*, "Error --- MakeMagGradGridDH"
        print*, "VXX, VYY, VZZ, VXY, VXZ, and VYZ must be dimensioned " // &
                "as: ", nlat_out, nlong_out
        print*, "Input dimensions are ", size(vxx(:,1)), size(vxx(1,:)), &
            size(vyy(:,1)), size(vyy(1,:)), size(vzz(:,1)), size(vzz(1,:)), &
            size(vxy(:,1)), size(vxy(1,:)), size(vxz(:,1)), size(vxz(1,:)), &
            size(vyz(:,1)), size(vyz(1,:))
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
            print*, "Error --- MakeMagGradGridDH"
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

        allocate (sqr(2 * lmax_comp + 1), stat=astat(1))
        allocate (ff1(lmax_comp+1,lmax_comp+1), stat=astat(2))
        allocate (ff2(lmax_comp+1,lmax_comp+1), stat=astat(3))
        allocate (fsymsign(lmax_comp+1,lmax_comp+1), stat=astat(4))

        if (sum(astat(1:4)) /= 0) then
            print*, "Error --- MakeMagGradGridDH"
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
        do l=1, 2*lmax_comp+1
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
            ff1(2,1) = 1.0_dp
            ff2(2,1) = 0.0_dp

        end if

        do l = 2, lmax_comp, 1
            ff1(l+1,1) = dble(2*l-1) / dble(l)
            ff2(l+1,1) = dble(l-1) / dble(l)

            do m = 1, l-2, 1
                ff1(l+1,m+1) = dble(2*l-1) / sqr(l+m) / sqr(l-m)
                ff2(l+1,m+1) = sqr(l-m-1) * sqr(l+m-1) / sqr(l+m) / sqr(l-m)
            end do

            ff1(l+1,l)= dble(2*l-1) / sqr(l+m) / sqr(l-m)
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
    lat = 0.0_dp

    coefr(1:lmax+2) = cmplx(0.0_dp, 0.0_dp, dp)
    coefr0 = 0.0_dp

    coefrr(1:lmax+2) = cmplx(0.0_dp, 0.0_dp, dp)
    coefrr0 = 0.0_dp

    coeft(1:lmax+2) = cmplx(0.0_dp, 0.0_dp, dp)
    coeft0 = 0.0_dp

    coeftp(1:lmax+2) = cmplx(0.0_dp, 0.0_dp, dp)
    coeftp0 = 0.0_dp

    coefp(1:lmax+2) = cmplx(0.0_dp, 0.0_dp, dp)
    coefp0 = 0.0_dp

    coefrt(1:lmax+2) = cmplx(0.0_dp, 0.0_dp, dp)
    coefrt0 = 0.0_dp

    coefrp(1:lmax+2) = cmplx(0.0_dp, 0.0_dp, dp)
    coefrp0 = 0.0_dp

    coefpp(1:lmax+2) = cmplx(0.0_dp, 0.0_dp, dp)
    coefpp0 = 0.0_dp

    coeftt(1:lmax+2) = cmplx(0.0_dp, 0.0_dp, dp)
    coeftt0 = 0.0_dp

    pm2 = 1.0_dp

    tempr = -cilm(1,1,1) * pm2 ! l = 0
    coefr0 = coefr0 + tempr

    tempr = 2 * cilm(1,1,1) * pm2  ! l = 0
    coefrr0 = coefrr0 + tempr

    ! derivative in theta and phi of l=0 term is 0. No need to calculate

    if (lmax_comp /= 0) then    ! l = 1
        prefactor(1) = r0 / r_ex

        do l = 2, lmax_comp, 1
            prefactor(l) = prefactor(l-1) * r0 / r_ex
        end do

        pm1 = 0.0_dp

        dpl = ff1(2,1)
        tempr = cilm(1,2,1) * dpl * prefactor(1)
        coeft0 = coeft0 + tempr

        tempr = cilm(1,2,1) * dpl * (-2) * prefactor(1)  ! -2 = (l+1) prefactor
        coefrt0 = coefrt0 + tempr

    end if

    do l = 2, lmax_comp, 1
        l1 = l + 1
        p = - ff2(l1,1) * pm2
        tempr = cilm(1,l1,1) * p * (-l1) * prefactor(l)
        coefr0 = coefr0 + tempr

        tempr = cilm(1,l1,1) * p * (l1) * (l1+1) * prefactor(l)
        coefrr0 = coefrr0 + tempr

        dpl = l * (pm1)
        tempr = cilm(1,l1,1) * dpl * prefactor(l)
        coeft0 = coeft0 + tempr

        tempr = cilm(1,l1,1) * dpl * (-l1) * prefactor(l)
        coefrt0 = coefrt0 + tempr

        dpl2 = -l*l1 * p
        tempr = cilm(1,l1,1) * dpl2 * prefactor(l)
        coeftt0 = coeftt0 + tempr

        pm2 = pm1
        pm1 = p

    end do

    pmm = sqr(2) * scalef

    rescalem = 1.0_dp / scalef

    do m = 1, lmax_comp-1, 1

        m1 = m + 1

        pmm = pmm * sqr(2*m+1) / sqr(2*m)
        pm2 = pmm / sqr(2*m+1)

        tempc = cmplx(cilm(1,m1,m1), - cilm(2,m1,m1), dp) * pm2 * (-m-1) &
                * prefactor(m)    ! (m,m)
        coefr(m1) = coefr(m1) + tempc

        tempc = cmplx(cilm(1,m1,m1), - cilm(2,m1,m1), dp) * pm2 * (m+1) &
                * (m+2) * prefactor(m) ! (m,m)
        coefrr(m1) = coefrr(m1) + tempc

        tempc = cmplx(cilm(2,m1,m1), cilm(1,m1,m1), dp) * pm2 &
                * prefactor(m) * m ! (m,m)
        coefp(m1) = coefp(m1) + tempc

        tempc = - cmplx(cilm(1,m1,m1), - cilm(2,m1,m1), dp) * pm2 &
                * prefactor(m) * m**2 ! (m,m)
        coefpp(m1) = coefpp(m1) + tempc

        tempc = cmplx(cilm(2,m1,m1), cilm(1,m1,m1), dp) * pm2 *(-m-1) &
                * prefactor(m) * m ! (m,m)
        coefrp(m1) = coefrp(m1) + tempc

        dpl2 = -(m*m1 -(m**2)) * pm2
        tempc = cmplx(cilm(1,m1,m1), - cilm(2,m1,m1), dp) * dpl2 &
                * prefactor(m) ! (m,m)
        coeftt(m1) = coeftt(m1) + tempc

        pm1 = 0.0_dp

        dpl = pm2 * sqr(2*m+1)
        tempc = cmplx(cilm(1,m1+1,m1), - cilm(2,m1+1,m1), dp) * dpl &
                * prefactor(m+1)    ! (m+1,m)
        coeft(m1) = coeft(m1) + tempc

        tempc = cmplx(cilm(1,m1+1,m1), - cilm(2,m1+1,m1), dp) * dpl * (-m-2) &
                * prefactor(m+1)   ! (m+1,m)
        coefrt(m1) = coefrt(m1) + tempc

        tempc = cmplx(cilm(2,m1+1,m1), cilm(1,m1+1,m1), dp) * dpl &
                * prefactor(m+1) * m  ! (m+1,m)
        coeftp(m1) = coeftp(m1) + tempc

        do l = m + 2, lmax_comp, 1
            l1 = l + 1
            p = - ff2(l1,m1) * pm2
            tempc = cmplx(cilm(1,l1,m1), - cilm(2,l1,m1), dp) * p * (-l1) &
                    * prefactor(l)
            coefr(m1) = coefr(m1) + tempc

            tempc = cmplx(cilm(1,l1,m1), - cilm(2,l1,m1), dp) * p * (l1) &
                    * (l1+1) * prefactor(l)
            coefrr(m1) = coefrr(m1) + tempc

            tempc = cmplx(cilm(2,l1,m1), cilm(1,l1,m1), dp) * p &
                    * prefactor(l) * m
            coefp(m1) = coefp(m1) + tempc

            tempc = - cmplx(cilm(1,l1,m1), - cilm(2,l1,m1), dp) * p &
                    * prefactor(l) * m**2
            coefpp(m1) = coefpp(m1) + tempc

            tempc = cmplx(cilm(2,l1,m1), cilm(1,l1,m1), dp) * p * (-l1) &
                    * prefactor(l) * m
            coefrp(m1) = coefrp(m1) + tempc

            dpl = sqr(l+m) * sqr(l-m) * pm1
            tempc = cmplx(cilm(1,l1,m1), - cilm(2,l1,m1), dp) * dpl &
                    * prefactor(l)
            coeft(m1) = coeft(m1) + tempc

            tempc = cmplx(cilm(1,l1,m1), - cilm(2,l1,m1), dp) * dpl * (-l1) &
                    * prefactor(l)
            coefrt(m1) = coefrt(m1) + tempc

            tempc = cmplx(cilm(2,l1,m1), cilm(1,l1,m1), dp) * dpl * &
                    prefactor(l) * m
            coeftp(m1) = coeftp(m1) + tempc

            dpl2 = - (l * l1 -(m**2) / u**2) * p
            tempc = cmplx(cilm(1,l1,m1), - cilm(2,l1,m1), dp) * dpl2 &
                    * prefactor(l)
            coeftt(m1) = coeftt(m1) + tempc

            pm2 = pm1
            pm1 = p

        end do

        coefr(m1) = coefr(m1) * rescalem
        coefrr(m1) = coefrr(m1) * rescalem
        coeft(m1) = coeft(m1) * rescalem
        coeftt(m1) = coeftt(m1) * rescalem
        coefrt(m1) = coefrt(m1) * rescalem
        coefp(m1) = coefp(m1) * rescalem
        coefpp(m1) = coefpp(m1) * rescalem
        coeftp(m1) = coeftp(m1) * rescalem
        coefrp(m1) = coefrp(m1) * rescalem

    end do

    if (lmax_comp /= 0) then

        pmm = pmm / sqr(2*lmax_comp) * rescalem
        tempc = cmplx(cilm(1,lmax_comp+1,lmax_comp+1), &
                        - cilm(2,lmax_comp+1,lmax_comp+1), dp) * pmm &
                        * (-lmax_comp-1) * prefactor(lmax_comp)
        coefr(lmax_comp+1) = coefr(lmax_comp+1) + tempc

        tempc = cmplx(cilm(1,lmax_comp+1,lmax_comp+1), &
                        - cilm(2,lmax_comp+1,lmax_comp+1), dp) * pmm &
                        * (lmax_comp+1) * (lmax_comp+2) * prefactor(lmax_comp)
        coefrr(lmax_comp+1) = coefrr(lmax_comp+1) + tempc

        tempc = cmplx(cilm(2,lmax_comp+1,lmax_comp+1), &
                        cilm(1,lmax_comp+1,lmax_comp+1), dp) * pmm &
                        * prefactor(lmax_comp) * lmax_comp
        coefp(lmax_comp+1) = coefp(lmax_comp+1) + tempc

        tempc = - cmplx(cilm(1,lmax_comp+1,lmax_comp+1), &
                        - cilm(2,lmax_comp+1,lmax_comp+1), dp) * pmm &
                        * prefactor(lmax_comp) * lmax_comp**2
        coefpp(lmax_comp+1) = coefpp(lmax_comp+1) + tempc

        tempc = cmplx(cilm(2,lmax_comp+1,lmax_comp+1), &
                        cilm(1,lmax_comp+1,lmax_comp+1), dp) * pmm &
                        * (-lmax_comp-1) * prefactor(lmax_comp) * lmax_comp
        coefrp(lmax_comp+1) = coefrp(lmax_comp+1) + tempc

        dpl2 = -(lmax_comp*(lmax_comp+1)-(lmax_comp**2)) * pmm
        tempc = cmplx(cilm(1,lmax_comp+1,lmax_comp+1), &
                        - cilm(2,lmax_comp+1,lmax_comp+1), dp) * dpl2 &
                        * prefactor(lmax_comp)
        coeftt(lmax_comp+1) = coeftt(lmax_comp+1) + tempc

    end if

    coefr0 = coefr0 * r0 / r_ex
    coefr(2:lmax+1) = coefr(2:lmax+1) * r0 / r_ex

    coefrt0 = - coefrt0 * r0 / r_ex
    coefrt(2:lmax+1) = - coefrt(2:lmax+1) * r0 / r_ex

    coefrr0 = coefrr0 * r0 / r_ex**2
    coefrr(2:lmax+1) = coefrr(2:lmax+1) * r0 / r_ex**2

    coeft0 = - coeft0 * r0
    coeft(2:lmax+1) = - coeft(2:lmax+1) * r0

    coeftt0 = coeftt0 * r0
    coeftt(2:lmax+1) = coeftt(2:lmax+1) * r0

    coeftp0 = - coeftp0 * r0
    coeftp(2:lmax+1) = - coeftp(2:lmax+1) * r0

    coefp0 = coefp0 * r0
    coefp(2:lmax+1) = coefp(2:lmax+1) * r0

    coefpp0 = coefpp0 * r0
    coefpp(2:lmax+1) = coefpp(2:lmax+1) * r0

    coefrp0 = coefrp0 * r0 / r_ex
    coefrp(2:lmax+1) = coefrp(2:lmax+1) * r0 / r_ex

    ! Vzz = Vrr
    coef(1) = cmplx(coefrr0, 0.0_dp, dp)
    coef(2:lmax+1) = coefrr(2:lmax+1) / 2.0_dp

    if (present(sampling)) then
        if (sampling == 2) then
            coef(lmax+2:2*lmax+3) = cmplx(0.0_dp, 0.0_dp, dp)
        end if
    end if

    call fftw_execute_dft_c2r(plan, coef, grid)
    vzz(i_eq,1:nlong) = grid(1:nlong)

    ! Vxx = 1/r Vr + 1/r^2 Vtt
    coef(1) = cmplx(coefr0/r_ex + coeftt0/r_ex**2, 0.0_dp, dp)
    coef(2:lmax+1) = (coefr(2:lmax+1)/r_ex + coeftt(2:lmax+1)/r_ex**2 ) / 2.0_dp

    if (present(sampling)) then
        if (sampling == 2) then
            coef(lmax+2:2*lmax+3) = cmplx(0.0_dp, 0.0_dp, dp)
        end if
    end if

    call fftw_execute_dft_c2r(plan, coef, grid)
    vxx(i_eq,1:nlong) = grid(1:nlong)

    ! Vyy = 1/r Vr + 1/r^2 /tan(t) Vt + 1/r^2 /sin(t)^2 Vpp
    coef(1) = cmplx(coefr0/r_ex + coefpp0/(r_ex**2), 0.0_dp, dp)
    coef(2:lmax+1) = (coefr(2:lmax+1)/r_ex &
                     + coefpp(2:lmax+1)/(r_ex**2) ) / 2.0_dp

    if (present(sampling)) then
        if (sampling == 2) then
            coef(lmax+2:2*lmax+3) = cmplx(0.0_dp, 0.0_dp, dp)
        end if
    end if

    call fftw_execute_dft_c2r(plan, coef, grid)
    vyy(i_eq,1:nlong) = grid(1:nlong)

    ! Vxy = 1/r^2 /sin(t) Vtp - cos(t)/sin(t)^2 /r^2 Vp
    coef(1) = cmplx(coeftp0/r_ex**2, 0.0_dp, dp)
    coef(2:lmax+1) = (coeftp(2:lmax+1)/r_ex**2 ) / 2.0_dp

    if (present(sampling)) then
        if (sampling == 2) then
            coef(lmax+2:2*lmax+3) = cmplx(0.0_dp, 0.0_dp, dp) 
        end if
    end if

    call fftw_execute_dft_c2r(plan, coef, grid)
    vxy(i_eq,1:nlong) = grid(1:nlong)

    ! Vxz = 1/r^2 Vt - 1/r Vrt
    coef(1) = cmplx(coeft0/r_ex**2 - coefrt0/r_ex, 0.0_dp, dp)
    coef(2:lmax+1) = (coeft(2:lmax+1) / r_ex**2 &
                    - coefrt(2:lmax+1)/r_ex ) / 2.0_dp

    if (present(sampling)) then
        if (sampling == 2) then
            coef(lmax+2:2*lmax+3) = cmplx(0.0_dp, 0.0_dp, dp)
        end if
    end if

    call fftw_execute_dft_c2r(plan, coef, grid)
    vxz(i_eq,1:nlong) = grid(1:nlong)

    ! Vyz = 1/r^2 /sin(t) Vp - 1/r /sin(t) Vrp
    coef(1) = cmplx(coefp0/r_ex**2 - coefrp0/r_ex, 0.0_dp, dp)
    coef(2:lmax+1) = (coefp(2:lmax+1) / r_ex**2 &
                    - coefrp(2:lmax+1)/r_ex ) / 2.0_dp

    if (present(sampling)) then
        if (sampling == 2) then
            coef(lmax+2:2*lmax+3) = cmplx(0.0_dp, 0.0_dp, dp)
        end if
    end if

    call fftw_execute_dft_c2r(plan, coef, grid)
    vyz(i_eq,1:nlong) = grid(1:nlong)

    do i = 1, i_eq - 1, 1

        i_s = 2 * i_eq - i

        theta = pi * dble(i-1) / dble(n)
        z = cos(theta)
        u = sqrt( (1.0_dp-z) * (1.0_dp+z) )
        sint = sin(theta)
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

        coefrr(1:lmax+2) = cmplx(0.0_dp, 0.0_dp, dp)
        coefrr0 = 0.0_dp
        coefrrs(1:lmax+2) = cmplx(0.0_dp, 0.0_dp, dp)
        coefrrs0 = 0.0_dp

        coeft(1:lmax+2) = cmplx(0.0_dp, 0.0_dp, dp)
        coeft0 = 0.0_dp
        coefts(1:lmax+2) = cmplx(0.0_dp, 0.0_dp, dp)
        coefts0 = 0.0_dp

        coeftp(1:lmax+2) = cmplx(0.0_dp, 0.0_dp, dp)
        coeftp0 = 0.0_dp
        coeftps(1:lmax+2) = cmplx(0.0_dp, 0.0_dp, dp)
        coeftps0 = 0.0_dp

        coefp(1:lmax+2) = cmplx(0.0_dp, 0.0_dp, dp)
        coefp0 = 0.0_dp
        coefps(1:lmax+2) = cmplx(0.0_dp, 0.0_dp, dp)
        coefps0 = 0.0_dp

        coefrt(1:lmax+2) = cmplx(0.0_dp, 0.0_dp, dp)
        coefrt0 = 0.0_dp
        coefrts(1:lmax+2) = cmplx(0.0_dp, 0.0_dp, dp)
        coefrts0 = 0.0_dp

        coefrp(1:lmax+2) = cmplx(0.0_dp, 0.0_dp, dp)
        coefrp0 = 0.0_dp
        coefrps(1:lmax+2) = cmplx(0.0_dp, 0.0_dp, dp)
        coefrps0 = 0.0_dp

        coefpp(1:lmax+2) = cmplx(0.0_dp, 0.0_dp, dp)
        coefpp0 = 0.0_dp
        coefpps(1:lmax+2) = cmplx(0.0_dp, 0.0_dp, dp)
        coefpps0 = 0.0_dp

        coeftt(1:lmax+2) = cmplx(0.0_dp, 0.0_dp, dp)
        coeftt0 = 0.0_dp
        coeftts(1:lmax+2) = cmplx(0.0_dp, 0.0_dp, dp)
        coeftts0 = 0.0_dp

        pm2 = 1.0_dp

        ! l = 0 terms are zero

        ! derivative in theta and phi of l=0 term is 0, so no need to
        ! calculate this

        if (lmax_comp /= 0) then    ! l = 1
            prefactor(1) = (r0 / r_ex)**2

            do l = 2, lmax_comp, 1
                prefactor(l) = prefactor(l-1) * r0 / r_ex
            end do

            pm1 = ff1(2,1) * z

            ! -2 = (l+1) prefactor
            tempr = cilm(1,2,1) * pm1 * (-2) * prefactor(1)
            coefr0 = coefr0 + tempr
            coefrs0 = coefrs0 - tempr   ! fsymsign = -1

            ! 6 = (l+1)*(l+2) prefactor
            tempr = cilm(1,2,1) * pm1 * (6) * prefactor(1)
            coefrr0 = coefrr0 + tempr
            coefrrs0 = coefrrs0 - tempr     ! fsymsign = -1

            ! dpl is the first derivative with respect to Z
            dpl = ff1(2,1)
            tempr = cilm(1,2,1) * dpl * prefactor(1)
            coeft0 = coeft0 + tempr
            coefts0 = coefts0 + tempr   ! reverse fsymsign

            ! -2 = (l+1) prefactor
            tempr = cilm(1,2,1) * dpl * (-2) * prefactor(1)
            coefrt0 = coefrt0 + tempr
            coefrts0 = coefrts0 + tempr     ! reverse fsymsign

            ! dpl2 is the second derivative with respect to THETA. Must
            ! multiply dpl by -sin(theta) in the recurrence relation
            ! Note that dpl2 is symmetric about the equator according
            ! the fsymsign.
            dpl2 = -2 * pm1 + z * dpl
            dpl2s = - dpl2
            tempr = cilm(1,2,1) * prefactor(1)
            coeftt0 = coeftt0 + tempr * dpl2
            coeftts0 = coeftts0 + tempr * dpl2s

        end if

        do l = 2, lmax_comp, 1
            l1 = l + 1
            p = ff1(l1,1) * z * pm1 - ff2(l1,1) * pm2
            tempr = cilm(1,l1,1) * p * (-l1) * prefactor(l)
            coefr0 = coefr0 + tempr
            coefrs0 = coefrs0 + tempr * fsymsign(l1,1)

            tempr = cilm(1,l1,1) * p * (l1) * (l1+1) * prefactor(l)
            coefrr0 = coefrr0 + tempr
            coefrrs0 = coefrrs0 + tempr * fsymsign(l1,1)

            dpl = l * (pm1 - z * p) / u**2
            tempr = cilm(1,l1,1) * dpl * prefactor(l)
            coeft0 = coeft0 + tempr
            coefts0 = coefts0 - tempr * fsymsign(l1,1)  ! reverse fsymsign

            tempr = cilm(1,l1,1) * dpl * (-l1) * prefactor(l)
            coefrt0 = coefrt0 + tempr
            coefrts0 = coefrts0 - tempr * fsymsign(l1,1)

            dpl2 = -l * l1 * p + z * dpl
            dpl2s = dpl2 * fsymsign(l1,1)
            tempr = cilm(1,l1,1) * prefactor(l)
            coeftt0 = coeftt0 + tempr * dpl2
            coeftts0 = coeftts0 + tempr * dpl2s

            pm2 = pm1
            pm1 = p

        end do

        pmm = sqr(2) * scalef

        rescalem = 1.0_dp / scalef

        do m = 1, lmax_comp-1, 1

            m1 = m + 1
            rescalem = rescalem * u

            pmm = pmm * sqr(2*m+1) / sqr(2*m)
            pm2 = pmm / sqr(2*m+1)

            tempc = cmplx(cilm(1,m1,m1), - cilm(2,m1,m1), dp) * pm2 * (-m-1) &
                    * prefactor(m)    ! (m,m)
            coefr(m1) = coefr(m1) + tempc
            coefrs(m1) = coefrs(m1) + tempc ! fsymsign = 1

            tempc = cmplx(cilm(1,m1,m1), - cilm(2,m1,m1), dp) * pm2 * (m+1) &
                    * (m+2) * prefactor(m) ! (m,m)
            coefrr(m1) = coefrr(m1) + tempc
            coefrrs(m1) = coefrrs(m1) + tempc ! fsymsign = 1

            tempc = cmplx(cilm(2,m1,m1), cilm(1,m1,m1), dp) * pm2 &
                    * prefactor(m) * m ! (m,m)
            coefp(m1) = coefp(m1) + tempc
            coefps(m1) = coefps(m1) + tempc ! fsymsign = 1

            tempc = - cmplx(cilm(1,m1,m1), - cilm(2,m1,m1), dp) * pm2 &
                    * prefactor(m) * m**2 ! (m,m)
            coefpp(m1) = coefpp(m1) + tempc
            coefpps(m1) = coefpps(m1) + tempc ! fsymsign = 1

            tempc = cmplx(cilm(2,m1,m1), cilm(1,m1,m1), dp) * pm2 *(-m-1) &
                    * prefactor(m) * m ! (m,m)
            coefrp(m1) = coefrp(m1) + tempc
            coefrps(m1) = coefrps(m1) + tempc ! fsymsign = 1

            dpl = -m * z * pm2 / u**2
            tempc = cmplx(cilm(1,m1,m1), - cilm(2,m1,m1), dp) * dpl &
                    * prefactor(m)  ! (m,m)
            coeft(m1) = coeft(m1) + tempc
            coefts(m1) = coefts(m1) - tempc ! reverse fsymsign

            tempc = cmplx(cilm(1,m1,m1), - cilm(2,m1,m1), dp) * dpl * (-m-1) &
                    * prefactor(m) ! (m,m)
            coefrt(m1) = coefrt(m1) + tempc
            coefrts(m1) = coefrts(m1) - tempc ! reverse fsymsign

            tempc = cmplx(cilm(2,m1,m1), cilm(1,m1,m1), dp) * dpl &
                    * prefactor(m) * m ! (m,m)
            coeftp(m1) = coeftp(m1) + tempc
            coeftps(m1) = coeftps(m1) - tempc ! reverse fsymsign

            dpl2 = -(m*m1 - (m**2)/u**2) * pm2 + z * dpl
            dpl2s = dpl2
            tempc = cmplx(cilm(1,m1,m1), - cilm(2,m1,m1), dp) &
                    * prefactor(m)   ! (m,m)
            coeftt(m1) = coeftt(m1) + tempc * dpl2
            coeftts(m1) = coeftts(m1) + tempc * dpl2s

            pm1 = z * ff1(m1+1,m1) * pm2
            tempc = cmplx(cilm(1,m1+1,m1), - cilm(2,m1+1,m1), dp) * pm1 &
                    * (-m-2) * prefactor(m+1)  ! (m+1,m)
            coefr(m1) = coefr(m1) + tempc 
            coefrs(m1) = coefrs(m1) - tempc ! fsymsign = -1

            tempc = cmplx(cilm(1,m1+1,m1), - cilm(2,m1+1,m1), dp) * pm1 &
                    * (m+2) * (m+3) * prefactor(m+1)   ! (m+1,m)
            coefrr(m1) = coefrr(m1) + tempc
            coefrrs(m1) = coefrrs(m1) - tempc ! fsymsign = -1

            tempc = cmplx(cilm(2,m1+1,m1), cilm(1,m1+1,m1), dp) * pm1 &
                    * prefactor(m+1) * m ! (m+1,m)
            coefp(m1) = coefp(m1) + tempc
            coefps(m1) = coefps(m1) - tempc ! fsymsign = -1

            tempc = - cmplx(cilm(1,m1+1,m1), - cilm(2,m1+1,m1), dp) * pm1 &
                    * prefactor(m+1) * m**2   ! (m+1,m)
            coefpp(m1) = coefpp(m1) + tempc
            coefpps(m1) = coefpps(m1) - tempc ! fsymsign = -1

            tempc = cmplx(cilm(2,m1+1,m1), cilm(1,m1+1,m1), dp) * pm1 &
                    * (-m-2) * prefactor(m+1) * m    ! (m+1,m)
            coefrp(m1) = coefrp(m1) + tempc
            coefrps(m1) = coefrps(m1) - tempc ! fsymsign = -1

            dpl = (pm2 * sqr(2*m+1) - z * (m+1) * pm1) / u**2
            tempc = cmplx(cilm(1,m1+1,m1), - cilm(2,m1+1,m1), dp) * dpl &
                    * prefactor(m+1)    ! (m+1,m)
            coeft(m1) = coeft(m1) + tempc
            coefts(m1) = coefts(m1) + tempc ! reverse fsymsign

            tempc = cmplx(cilm(1,m1+1,m1), - cilm(2,m1+1,m1), dp) * dpl &
                    * (-m-2) * prefactor(m+1)   ! (m+1,m)
            coefrt(m1) = coefrt(m1) + tempc
            coefrts(m1) = coefrts(m1) + tempc   ! reverse fsymsign

            tempc = cmplx(cilm(2,m1+1,m1), cilm(1,m1+1,m1), dp) * dpl &
                    * prefactor(m+1) * m  ! (m+1,m)
            coeftp(m1) = coeftp(m1) + tempc
            coeftps(m1) = coeftps(m1) + tempc ! reverse fsymsign

            dpl2 = -(m1*(m1+1) - (m**2)/u**2) * pm1 + z * dpl
            dpl2s = - dpl2
            tempc = cmplx(cilm(1,m1+1,m1), - cilm(2,m1+1,m1), dp) &
                    * prefactor(m+1) ! (m+1,m)
            coeftt(m1) = coeftt(m1) + tempc  * dpl2
            coeftts(m1) = coeftts(m1) + tempc * dpl2s

            do l=m+2, lmax_comp, 1
                l1 = l + 1
                p = z * ff1(l1,m1) * pm1 - ff2(l1,m1) * pm2
                tempc = cmplx(cilm(1,l1,m1), - cilm(2,l1,m1), dp) * p * (-l1) &
                        * prefactor(l)
                coefr(m1) = coefr(m1) + tempc
                coefrs(m1) = coefrs(m1) + tempc * fsymsign(l1,m1)

                tempc = cmplx(cilm(1,l1,m1), - cilm(2,l1,m1), dp) * p * (l1) &
                        * (l1+1) * prefactor(l)
                coefrr(m1) = coefrr(m1) + tempc
                coefrrs(m1) = coefrrs(m1) + tempc * fsymsign(l1,m1)

                tempc = cmplx(cilm(2,l1,m1), cilm(1,l1,m1), dp) * p &
                        * prefactor(l) * m
                coefp(m1) = coefp(m1) + tempc
                coefps(m1) = coefps(m1) + tempc * fsymsign(l1,m1)

                tempc = - cmplx(cilm(1,l1,m1), - cilm(2,l1,m1), dp) * p &
                        * prefactor(l) * m**2
                coefpp(m1) = coefpp(m1) + tempc
                coefpps(m1) = coefpps(m1) + tempc * fsymsign(l1,m1)

                tempc = cmplx(cilm(2,l1,m1), cilm(1,l1,m1), dp) * p * (-l1) &
                        * prefactor(l) * m
                coefrp(m1) = coefrp(m1) + tempc
                coefrps(m1) = coefrps(m1) + tempc * fsymsign(l1,m1)

                dpl = ( sqr(l+m) * sqr(l-m) * pm1 - l * z * p ) / u**2
                tempc = cmplx(cilm(1,l1,m1), - cilm(2,l1,m1), dp) * dpl &
                        * prefactor(l)
                coeft(m1) = coeft(m1) + tempc
                ! reverse fsymsign
                coefts(m1) = coefts(m1) - tempc * fsymsign(l1,m1)

                tempc = cmplx(cilm(1,l1,m1), - cilm(2,l1,m1), dp) * dpl &
                        * (-l1) * prefactor(l)
                coefrt(m1) = coefrt(m1) + tempc
                ! reverse fsymsign
                coefrts(m1) = coefrts(m1) - tempc * fsymsign(l1,m1)

                tempc = cmplx(cilm(2,l1,m1), cilm(1,l1,m1), dp) * dpl &
                        * prefactor(l) * m
                coeftp(m1) = coeftp(m1) + tempc
                ! reverse fsymsign
                coeftps(m1) = coeftps(m1) - tempc * fsymsign(l1,m1)

                dpl2 = -(l * l1 -(m**2)/u**2) * p + z * dpl
                dpl2s = dpl2 * fsymsign(l1,m1) 
                tempc = cmplx(cilm(1,l1,m1), - cilm(2,l1,m1), dp) &
                        * prefactor(l)
                coeftt(m1) = coeftt(m1) + tempc * dpl2
                coeftts(m1) = coeftts(m1) + tempc * dpl2s

                pm2 = pm1
                pm1 = p

            end do

            coefr(m1) = coefr(m1) * rescalem
            coefrs(m1) = coefrs(m1) * rescalem

            coefrr(m1) = coefrr(m1) * rescalem
            coefrrs(m1) = coefrrs(m1) * rescalem

            coeft(m1) = coeft(m1) * rescalem
            coefts(m1) = coefts(m1) * rescalem

            coeftt(m1) = coeftt(m1) * rescalem
            coeftts(m1) = coeftts(m1) * rescalem

            coefrt(m1) = coefrt(m1) * rescalem
            coefrts(m1) = coefrts(m1) * rescalem

            coefp(m1) = coefp(m1) * rescalem
            coefps(m1) = coefps(m1) * rescalem

            coefpp(m1) = coefpp(m1) * rescalem
            coefpps(m1) = coefpps(m1) * rescalem

            coeftp(m1) = coeftp(m1) * rescalem
            coeftps(m1) = coeftps(m1) * rescalem

            coefrp(m1) = coefrp(m1) * rescalem
            coefrps(m1) = coefrps(m1) * rescalem

        end do

        if (lmax_comp /= 0) then

            rescalem = rescalem * u

            pmm = pmm / sqr(2*lmax_comp) * rescalem
            tempc = cmplx(cilm(1,lmax_comp+1,lmax_comp+1), &
                            - cilm(2,lmax_comp+1,lmax_comp+1), dp) * pmm &
                            * (-lmax_comp-1) * prefactor(lmax_comp)
            coefr(lmax_comp+1) = coefr(lmax_comp+1) + tempc
            coefrs(lmax_comp+1) = coefrs(lmax_comp+1) + tempc ! fsymsign = 1

            tempc = cmplx(cilm(1,lmax_comp+1,lmax_comp+1), &
                            - cilm(2,lmax_comp+1,lmax_comp+1), dp) * pmm &
                            * (lmax_comp+1) * (lmax_comp+2) &
                            * prefactor(lmax_comp)
            coefrr(lmax_comp+1) = coefrr(lmax_comp+1) + tempc
            coefrrs(lmax_comp+1) = coefrrs(lmax_comp+1) + tempc ! fsymsign = 1

            tempc = cmplx(cilm(2,lmax_comp+1,lmax_comp+1), &
                            cilm(1,lmax_comp+1,lmax_comp+1), dp) * pmm &
                            * prefactor(lmax_comp) * lmax_comp
            coefp(lmax_comp+1) = coefp(lmax_comp+1) + tempc
            coefps(lmax_comp+1) = coefps(lmax_comp+1) + tempc ! fsymsign = 1

            tempc = - cmplx(cilm(1,lmax_comp+1,lmax_comp+1), &
                            - cilm(2,lmax_comp+1,lmax_comp+1), dp) * pmm &
                            * prefactor(lmax_comp) * lmax_comp**2
            coefpp(lmax_comp+1) = coefpp(lmax_comp+1) + tempc
            coefpps(lmax_comp+1) = coefpps(lmax_comp+1) + tempc ! fsymsign = 1

            tempc = cmplx(cilm(2,lmax_comp+1,lmax_comp+1), &
                            cilm(1,lmax_comp+1,lmax_comp+1), dp) * pmm &
                            * (-lmax_comp-1) * prefactor(lmax_comp) * lmax_comp
            coefrp(lmax_comp+1) = coefrp(lmax_comp+1) + tempc
            coefrps(lmax_comp+1) = coefrps(lmax_comp+1) + tempc ! fsymsign = 1

            dpl = -lmax_comp * z * pmm / u**2
            tempc = cmplx(cilm(1,lmax_comp+1,lmax_comp+1), &
                            - cilm(2,lmax_comp+1,lmax_comp+1), dp) * dpl &
                            * prefactor(lmax_comp)
            coeft(lmax_comp+1) = coeft(lmax_comp+1) + tempc
            ! reverse fsymsign
            coefts(lmax_comp+1) = coefts(lmax_comp+1) - tempc

            tempc = cmplx(cilm(1,lmax_comp+1,lmax_comp+1), &
                            - cilm(2,lmax_comp+1,lmax_comp+1), dp) * dpl &
                            * (-lmax_comp-1) * prefactor(lmax_comp)
            coefrt(lmax_comp+1) = coefrt(lmax_comp+1) + tempc
            ! reverse fsymsign
            coefrts(lmax_comp+1) = coefrts(lmax_comp+1) - tempc

            tempc = cmplx(cilm(2,lmax_comp+1,lmax_comp+1), &
                            cilm(1,lmax_comp+1,lmax_comp+1), dp) * dpl &
                            * prefactor(lmax_comp) * lmax_comp
            coeftp(lmax_comp+1) = coeftp(lmax_comp+1) + tempc
            coeftps(lmax_comp+1) = coeftps(lmax_comp+1) - tempc

            dpl2 = -(lmax_comp*(lmax_comp+1)-(lmax_comp**2)/u**2) * pmm + z &
                   * dpl
            dpl2s = dpl2
            tempc = cmplx(cilm(1,lmax_comp+1,lmax_comp+1), &
                            - cilm(2,lmax_comp+1,lmax_comp+1), dp) &
                            * prefactor(lmax_comp)
            coeftt(lmax_comp+1) = coeftt(lmax_comp+1) + tempc * dpl2
            coeftts(lmax_comp+1) = coeftts(lmax_comp+1) + tempc * dpl2s

        end if

        ! Note that the first angular derivatives are with repsect to z,
        ! but that the second is with respect to theta.

        coefr0 = coefr0 * r0 / r_ex
        coefr(2:lmax+1) = coefr(2:lmax+1) * r0 / r_ex

        coefrs0 = coefrs0 * r0 / r_ex
        coefrs(2:lmax+1) = coefrs(2:lmax+1) * r0 / r_ex

        coefrt0 = -sint * coefrt0 * r0 / r_ex
        coefrt(2:lmax+1) = -sint * coefrt(2:lmax+1) * r0 / r_ex

        coefrts0 = -sint * coefrts0 * r0 / r_ex
        coefrts(2:lmax+1) = -sint * coefrts(2:lmax+1) * r0 / r_ex

        coefrr0 = coefrr0 * r0 / r_ex**2
        coefrr(2:lmax+1) = coefrr(2:lmax+1) * r0 / r_ex**2

        coefrrs0 = coefrrs0 * r0 / r_ex**2
        coefrrs(2:lmax+1) = coefrrs(2:lmax+1) * r0 / r_ex**2

        coeft0 = -sint * coeft0 * r0
        coeft(2:lmax+1) = -sint * coeft(2:lmax+1) * r0

        coefts0 = -sint * coefts0 * r0
        coefts(2:lmax+1) = -sint * coefts(2:lmax+1) * r0

        coeftt0 = coeftt0 * r0
        coeftt(2:lmax+1) = coeftt(2:lmax+1) * r0

        coeftts0 = coeftts0 * r0
        coeftts(2:lmax+1) = coeftts(2:lmax+1) * r0

        coeftp0 = -sint * coeftp0 * r0
        coeftp(2:lmax+1) = -sint * coeftp(2:lmax+1) * r0

        coeftps0 = -sint * coeftps0 * r0
        coeftps(2:lmax+1) = -sint * coeftps(2:lmax+1) * r0

        coefp0 = coefp0 * r0
        coefp(2:lmax+1) = coefp(2:lmax+1) * r0

        coefps0 = coefps0 * r0
        coefps(2:lmax+1) = coefps(2:lmax+1) * r0

        coefpp0 = coefpp0 * r0
        coefpp(2:lmax+1) = coefpp(2:lmax+1) * r0

        coefpps0 = coefpps0 * r0
        coefpps(2:lmax+1) = coefpps(2:lmax+1) * r0

        coefrp0 = coefrp0 * r0 / r_ex
        coefrp(2:lmax+1) = coefrp(2:lmax+1) * r0 / r_ex

        coefrps0 = coefrps0 * r0 / r_ex
        coefrps(2:lmax+1) = coefrps(2:lmax+1) * r0 / r_ex

        ! Vzz = Vrr
        coef(1) = cmplx(coefrr0, 0.0_dp, dp)
        coef(2:lmax+1) = coefrr(2:lmax+1) / 2.0_dp

        if (present(sampling)) then
            if (sampling == 2) then
                coef(lmax+2:2*lmax+3) = cmplx(0.0_dp, 0.0_dp, dp)
            end if

        end if

        call fftw_execute_dft_c2r(plan, coef, grid)
        vzz(i,1:nlong) = grid(1:nlong)

        if (i == 1) then
            ! These derivatives are undefined at the pole
            vxx(1,1:nlong) = 0.0_dp
            vyy(1,1:nlong) = 0.0_dp
            vxy(1,1:nlong) = 0.0_dp
            vxz(1,1:nlong) = 0.0_dp
            vyz(1,1:nlong) = 0.0_dp

        else
            ! Vxx = 1/r Vr + 1/r^2 Vtt
            coef(1) = cmplx(coefr0/r_ex + coeftt0/r_ex**2, 0.0_dp, dp)
            coef(2:lmax+1) = (coefr(2:lmax+1)/r_ex &
                             + coeftt(2:lmax+1)/r_ex**2 ) / 2.0_dp

            if (present(sampling)) then
                if (sampling == 2) then
                    coef(lmax+2:2*lmax+3) = cmplx(0.0_dp, 0.0_dp, dp)
                end if
            end if

            call fftw_execute_dft_c2r(plan, coef, grid)
            vxx(i,1:nlong) = grid(1:nlong)

            ! Vyy = 1/r Vr + 1/r^2 /tan(t) Vt + 1/r^2 /sin(t)^2 Vpp
            coef(1) = cmplx(coefr0/r_ex + coeft0/(r_ex**2)/tan(theta) &
                            + coefpp0/(r_ex**2)/u**2, 0.0_dp, dp)
            coef(2:lmax+1) = (coefr(2:lmax+1)/r_ex &
                            + coeft(2:lmax+1)/(r_ex**2)/tan(theta) + &
                            coefpp(2:lmax+1)/(r_ex**2)/u**2 ) / 2.0_dp

            if (present(sampling)) then
                if (sampling == 2) then
                    coef(lmax+2:2*lmax+3) = cmplx(0.0_dp, 0.0_dp, dp)
                end if
            end if

            call fftw_execute_dft_c2r(plan, coef, grid)
            vyy(i,1:nlong) = grid(1:nlong)

            ! Vxy = 1/r^2 /sin(t) Vtp - cos(t)/sin(t)^2 /r^2 Vp
            coef(1) = cmplx(coeftp0/sint/r_ex**2 &
                            - coefp0/(r_ex**2)*z/u**2, 0.0_dp, dp)
            coef(2:lmax+1) = (coeftp(2:lmax+1)/sint/r_ex**2 &
                            - coefp(2:lmax+1)/(r_ex**2)*z/u**2 ) / 2.0_dp

            if (present(sampling)) then
                if (sampling == 2) then
                    coef(lmax+2:2*lmax+3) = cmplx(0.0_dp, 0.0_dp, dp)
                end if
            end if

            call fftw_execute_dft_c2r(plan, coef, grid)
            vxy(i,1:nlong) = grid(1:nlong)

            ! Vxz = 1/r^2 Vt - 1/r Vrt
            coef(1) = cmplx(coeft0/r_ex**2 - coefrt0/r_ex, 0.0_dp, dp)
            coef(2:lmax+1) = (coeft(2:lmax+1)/r_ex**2 &
                            - coefrt(2:lmax+1)/r_ex ) / 2.0_dp

            if (present(sampling)) then
                if (sampling == 2) then
                    coef(lmax+2:2*lmax+3) = cmplx(0.0_dp, 0.0_dp, dp)
                end if
            end if

            call fftw_execute_dft_c2r(plan, coef, grid)
            vxz(i,1:nlong) = grid(1:nlong)

            ! Vyz = 1/r^2 /sin(t) Vp - 1/r /sin(t) Vrp
            coef(1) = cmplx(coefp0/sint/r_ex**2 - coefrp0/sint/r_ex, 0.0_dp, dp)
            coef(2:lmax+1) = (coefp(2:lmax+1)/sint/r_ex**2 &
                            - coefrp(2:lmax+1)/sint/r_ex ) / 2.0_dp

            if (present(sampling)) then
                if (sampling == 2) then
                    coef(lmax+2:2*lmax+3) = cmplx(0.0_dp, 0.0_dp, dp)
                end if
            end if

            call fftw_execute_dft_c2r(plan, coef, grid)
            vyz(i,1:nlong) = grid(1:nlong)

        end if

        ! don't compute value for south pole when extend = 0.
        if (.not. (i == 1 .and. extend_grid == 0) ) then
            ! Vzz = Vrr
            coef(1) = cmplx(coefrrs0, 0.0_dp, dp)
            coef(2:lmax+1) = coefrrs(2:lmax+1) / 2.0_dp

            if (present(sampling)) then
                if (sampling == 2) then
                    coef(lmax+2:2*lmax+3) = cmplx(0.0_dp, 0.0_dp, dp)
                end if
            end if

            call fftw_execute_dft_c2r(plan, coef, grid)
            vzz(i_s,1:nlong) = grid(1:nlong)

            if (i == 1) then
                ! These derivatives are undefined at the pole
                vxx(i_s,1:nlong) = 0.0_dp
                vyy(i_s,1:nlong) = 0.0_dp
                vxy(i_s,1:nlong) = 0.0_dp
                vxz(i_s,1:nlong) = 0.0_dp
                vyz(i_s,1:nlong) = 0.0_dp

            else
                ! Vxx = 1/r Vr + 1/r^2 Vtt
                coef(1) = cmplx(coefrs0/r_ex + coeftts0/r_ex**2, 0.0_dp, dp)
                coef(2:lmax+1) = (coefrs(2:lmax+1)/r_ex &
                                  + coeftts(2:lmax+1)/r_ex**2 ) / 2.0_dp

                if (present(sampling)) then
                    if (sampling == 2) then
                        coef(lmax+2:2*lmax+3) = cmplx(0.0_dp, 0.0_dp, dp)
                    end if
                end if

                call fftw_execute_dft_c2r(plan, coef, grid)
                vxx(i_s,1:nlong) = grid(1:nlong)

                ! Vyy = 1/r Vr + 1/r^2 /tan(t) Vt + 1/r^2 /sin(t)^2 Vpp
                ! Note that tan(t) changes sign in the southern hemisphere
                coef(1) = cmplx(coefrs0/r_ex - coefts0/(r_ex**2)/tan(theta) &
                                + coefpps0/(r_ex**2)/u**2, 0.0_dp, dp)
                coef(2:lmax+1) = (coefrs(2:lmax+1)/r_ex &
                                  - coefts(2:lmax+1)/(r_ex**2)/tan(theta) + &
                                  coefpps(2:lmax+1)/(r_ex**2)/u**2 ) / 2.0_dp

                if (present(sampling)) then
                    if (sampling == 2) then
                        coef(lmax+2:2*lmax+3) = cmplx(0.0_dp, 0.0_dp, dp)
                    end if
                end if

                call fftw_execute_dft_c2r(plan, coef, grid)
                vyy(i_s,1:nlong) = grid(1:nlong)

                ! Vxy = 1/r^2 /sin(t) Vtp - cos(t)/sin(t)^2 /r^2 Vp
                ! Note that cos(t) changes sign in the southern hemisphere
                coef(1) = cmplx(coeftps0/sint/r_ex**2 &
                                + coefps0/(r_ex**2)*z/u**2, 0.0_dp, dp)
                coef(2:lmax+1) = (coeftps(2:lmax+1)/sint/r_ex**2 &
                                  + coefps(2:lmax+1)/(r_ex**2)*z/u**2 ) &
                                  / 2.0_dp

                if (present(sampling)) then
                    if (sampling == 2) then
                        coef(lmax+2:2*lmax+3) = cmplx(0.0_dp, 0.0_dp, dp)
                    end if
                end if

                call fftw_execute_dft_c2r(plan, coef, grid)
                vxy(i_s,1:nlong) = grid(1:nlong)

                ! Vxz = 1/r^2 Vt - 1/r Vrt
                coef(1) = cmplx(coefts0/r_ex**2 - coefrts0/r_ex, 0.0_dp, dp)
                coef(2:lmax+1) = (coefts(2:lmax+1)/r_ex**2 &
                                  - coefrts(2:lmax+1)/r_ex ) / 2.0_dp

                if (present(sampling)) then
                    if (sampling == 2) then
                        coef(lmax+2:2*lmax+3) = cmplx(0.0_dp, 0.0_dp, dp)
                    end if
                end if

                call fftw_execute_dft_c2r(plan, coef, grid)
                vxz(i_s,1:nlong) = grid(1:nlong)

                ! Vyz = 1/r^2 /sin(t) Vp - 1/r /sin(t) Vrp
                coef(1) = cmplx(coefps0/sint/r_ex**2 - coefrps0/sint/r_ex, &
                                0.0_dp, dp)
                coef(2:lmax+1) = (coefps(2:lmax+1)/sint/r_ex**2 &
                                  - coefrps(2:lmax+1)/sint/r_ex ) / 2.0_dp

                if (present(sampling)) then
                    if (sampling == 2) then
                        coef(lmax+2:2*lmax+3) = cmplx(0.0_dp, 0.0_dp, dp)
                    end if
                end if

                call fftw_execute_dft_c2r(plan, coef, grid)
                vyz(i_s,1:nlong) = grid(1:nlong)

            end if

        end if

    end do

    if (extend_grid == 1) then
        vxx(1:nlat_out, nlong_out) = vxx(1:nlat_out, 1)
        vyy(1:nlat_out, nlong_out) = vyy(1:nlat_out, 1)
        vzz(1:nlat_out, nlong_out) = vzz(1:nlat_out, 1)
        vxy(1:nlat_out, nlong_out) = vxy(1:nlat_out, 1)
        vxz(1:nlat_out, nlong_out) = vxz(1:nlat_out, 1)
        vyz(1:nlat_out, nlong_out) = vyz(1:nlat_out, 1)
    end if

    call fftw_destroy_plan(plan)

end subroutine MakeMagGradGridDH
