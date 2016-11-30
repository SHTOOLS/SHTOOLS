subroutine MakeGridGLQC(gridglq, cilm, lmax, plx, zero, norm, csphase, &
                        lmax_calc, exitstatus)
!------------------------------------------------------------------------------
!
!   Given the complex spherical harmonic coefficients CILM, this subroutine
!   will evalate these coefficients on a grid with equal spacing in
!   longitude and with latitude points appropriate for Gauss-Legendre
!   quadrature integrations. Note that this is done using FFTs for each degree
!   of each latitude band. The grid spacing is determined by the spherical
!   harmonic bandwidth LMAX. Nevertheless, the coefficients can be evaluated
!   up to smaller spherical harmonic degree by specifying the optional parameter
!   LMAX_CALC.
!
!   The optional array PLX contains precomputed associated legendre functions
!   evaluated on the Gauss-Legendre quadrature nodes (obtained from SHGLQ with
!   CNORM=1) and should not be precomputed when memory is an issue (i.e.,
!   lmax>360). If PLX is not present, the Legendre functions are computed on
!   the fly using the scaling methodolgy presented in Holmes and Featherston
!   (2002). When NORM=1,2 or 4, these are accurate to degree 2800. When NORM=3,
!   the routine is only stable to about degree 15!
!
!   Calling Parameters
!
!       IN
!           cilm        Input spherical harmonic coefficients with
!                       dimensions (2, LMAX+1, LMAX+1).
!           lmax        Maximum spherical harmonic degree used in the expansion.
!                       This value determines the grid spacing of the output
!                       function.
!
!       OUT
!           gridglq     Gridded data of the spherical harmonic
!                       coefficients CILM with dimensions (LMAX+1 , 2*LMAX+1).
!                       The first index (latitude) corresponds to the
!                       Gauss points, and the second index corresponds to
!                       360*(k-1)/nlong = 360*(k-1)/(2*LMAX +1).
!
!       OPTIONAL (IN)
!           plx         Input array of Associated Legendre Polnomials computed
!                       at the Gauss points (determined from a call to
!                       SHGLQ). If this is not included, then the optional
!                       array ZERO MUST be inlcuded. PLX must be computed using
!                       CNORM=1.
!           zero        Array of dimension lmax+1 that contains the latitudinal
!                       gridpoints used in the Gauss-Legendre quadrature
!                       integration
!                       scheme. Only needed if PLX is not included.
!           norm        Normalization to be used when calculating Legendre
!                       functions
!                           (1) "geodesy" (default)
!                           (2) Schmidt
!                           (3) unnormalized
!                           (4) orthonormalized
!           csphase 1   Do not include the phase factor of (-1)^m
!                       -1: Apply the phase factor of (-1)^m.
!           lmax_calc   The maximum spherical harmonic degree to evaluate the
!                       coefficients up to.
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
!       2.  Latitudes are geocentric latitude.
!
!   Dependencies:   FFTW3, CSPHASE_DEFAULT
!
!   Copyright (c) 2016, SHTOOLS
!   All rights reserved.
!
!------------------------------------------------------------------------------
    use FFTW3
    use SHTOOLS, only: CSPHASE_DEFAULT
#ifdef FFTW3_UNDERSCORE
#define dfftw_plan_dft_1d dfftw_plan_dft_1d_
#define dfftw_execute dfftw_execute_
#define dfftw_destroy_plan dfftw_destroy_plan_
#endif

    implicit none

    complex*16, intent(in) :: cilm(:,:,:)
    real*8, intent(in), optional :: plx(:,:), zero(:)
    complex*16, intent(out) :: gridglq(:,:)
    integer, intent(in) :: lmax
    integer, intent(in), optional :: norm, csphase, lmax_calc
    integer, intent(out), optional :: exitstatus
    integer :: l, m, i, nlat, nlong, l1, m1, lmax_comp, i_s, astat(4), lnorm, k
    real*8 :: pi, scalef, rescalem, u, p, pmm, pm1, pm2, z
    complex*16 :: coef(2*lmax+1), coefs(2*lmax+1), grid(2*lmax+1)
    integer*8 :: plan
    real*8, save, allocatable :: ff1(:,:), ff2(:,:), sqr(:)
    integer*1, save, allocatable :: fsymsign(:,:)
    integer, save :: lmax_old = 0, norm_old = 0
    integer :: phase
    external :: dfftw_plan_dft_1d, dfftw_execute, dfftw_destroy_plan

!$OMP   threadprivate(ff1, ff2, sqr, fsymsign, lmax_old, norm_old)

    if (present(exitstatus)) exitstatus = 0

    if (size(cilm(:,1,1)) < 2) then
        print*, "Error --- MakeGridGLQC"
        print*, "CILM must be dimensioned as (2, *, *)."
        print*, "Input dimension is ", size(cilm(:,1,1)), size(cilm(1,:,1)), &
            size(cilm(1,1,:))
        if (present(exitstatus)) then
            exitstatus = 1
            return
        else
            stop
        end if

    else if (size(gridglq(1,:)) < 2*lmax+1 .or. &
                size(gridglq(:,1)) < lmax+1 ) then
        print*, "Error --- MakeGridGLQC"
        print*, "GRIDGLQ must be dimensioned as (LMAX+1, 2*LMAX+1) " // &
                "where LMAX is ", lmax
        print*, "Input array is dimensioned ", size(gridglq(:,1)), &
                size(gridglq(1,:))
        if (present(exitstatus)) then
            exitstatus = 1
            return
        else
            stop
        end if

    end if

    if (present(plx)) then
        if (size(plx(:,1)) < lmax+1 .or. &
            size(plx(1,:)) < ((lmax+1)*(lmax+2))/2) then
            print*, "Error --- MakeGridGLQC"
            print*, "PLX must be dimensioned as (LMAX+1, " // &
                    "(LMAX+1)*(LMAX+2)/2) where LMAX is ", lmax
            print*, "Input array is dimensioned as ", size(plx(:,1)), &
                    size(plx(1,:))
            if (present(exitstatus)) then
                exitstatus = 1
                return
            else
                stop
            end if
        end if

    else if (present(zero)) then
        if (size(zero) < lmax+1) then
            print*, "Error --- MakeGridGLQC"
            print*, "ZERO must be dimensioned as (LMAX+1) where LMAX is ", lmax
            print*, "Input array is dimensioned ", size(zero)
            if (present(exitstatus)) then
                exitstatus = 1
                return
            else
                stop
            end if
        end if

    else
        print*, "Error --- MakeGridGLQC"
        print*, "Either PLX or ZERO must be specified."
        if (present(exitstatus)) then
            exitstatus = 5
            return
        else
            stop
        end if

    end if

    if (present(norm)) then
        if (norm > 4 .or. norm < 1) then
            print*, "Error - MakeGridGLQC"
            print*, "Parameter NORM must be 1 (geodesy), 2 (Schmidt), " // &
                    "3 (unnormalized), or 4 (orthonormalized)."
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
            print*, "Error --- MakeGridGLQC"
            print*, "CSPHASE must be 1 (exclude) or -1 (include)."
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
        phase = CSPHASE_DEFAULT

    end if

    pi = acos(-1.0d0)

    nlong = 2 * lmax + 1

    nlat = lmax+1

    scalef = 1.0d-280

    if (present(lmax_calc)) then
        if (lmax_calc > lmax) then
            print*, "Error --- MakeGridGLQC"
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
    !   Calculate recursion constants used in computing the Legendre polynomials
    !
    !--------------------------------------------------------------------------
    if ((lmax_comp /= lmax_old .or. lnorm /= norm_old) .and. &
            .not. present(plx)) then

        if (allocated (sqr)) deallocate (sqr)
        if (allocated (ff1)) deallocate (ff1)
        if (allocated (ff2)) deallocate (ff2)
        if (allocated (fsymsign)) deallocate (fsymsign)

        allocate (sqr(2*lmax_comp+1), stat=astat(1))
        allocate (ff1(lmax_comp+1,lmax_comp+1), stat=astat(2))
        allocate (ff2(lmax_comp+1,lmax_comp+1), stat=astat(3))
        allocate (fsymsign(lmax_comp+1,lmax_comp+1), stat=astat(4))

        if (sum(astat(1:4)) /= 0) then
            print*, "MakeGridGLQ --- Error"
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
        !   equator.
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
        do l = 1, 2 * lmax_comp+1
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
        select case (lnorm)
            case (1,4)
                if (lmax_comp /= 0) then
                    ff1(2,1) = sqr(3)
                    ff2(2,1) = 0.0d0
                end if

                do l = 2, lmax_comp, 1
                    ff1(l+1,1) = sqr(2*l-1) * sqr(2*l+1) / dble(l)
                    ff2(l+1,1) = dble(l-1) * sqr(2*l+1) / sqr(2*l-3) / dble(l)

                    do m = 1, l-2, 1
                        ff1(l+1,m+1) = sqr(2*l+1) * sqr(2*l-1) / sqr(l+m) &
                                       / sqr(l-m)
                        ff2(l+1,m+1) = sqr(2*l+1) * sqr(l-m-1) * sqr(l+m-1) &
                                       / sqr(2*l-3) / sqr(l+m) / sqr(l-m) 
                    end do

                    ff1(l+1,l) = sqr(2*l+1) * sqr(2*l-1) / sqr(l+m) / sqr(l-m)
                    ff2(l+1,l) = 0.0d0

                end do

            case (2)
                if (lmax_comp /= 0) then
                    ff1(2,1) = 1.0d0
                    ff2(2,1) = 0.0d0
                end if

                do l = 2, lmax_comp, 1
                    ff1(l+1,1) = dble(2*l-1) / dble(l)
                    ff2(l+1,1) = dble(l-1) / dble(l)

                    do m = 1, l - 2, 1
                        ff1(l+1,m+1) = dble(2*l-1) / sqr(l+m) / sqr(l-m)
                        ff2(l+1,m+1) = sqr(l-m-1) * sqr(l+m-1) / sqr(l+m) &
                                       / sqr(l-m)
                    enddo

                    ff1(l+1,l)= dble(2*l-1) / sqr(l+m) / sqr(l-m)
                    ff2(l+1,l) = 0.0d0

                end do

            case (3)
                do l = 1, lmax_comp, 1
                    ff1(l+1,1) = dble(2*l-1) / dble(l)
                    ff2(l+1,1) = dble(l-1) / dble(l)

                    do m = 1, l-1, 1
                        ff1(l+1,m+1) = dble(2*l-1) / dble(l-m)
                        ff2(l+1,m+1) = dble(l+m-1) / dble(l-m)
                    end do

                end do

        end select

        lmax_old = lmax_comp
        norm_old = lnorm

    end if

    !--------------------------------------------------------------------------
    !
    !   Do special case of lmax_comp = 0
    !
    !--------------------------------------------------------------------------
    if (lmax_comp == 0) then
        select case (lnorm)
            case (1,2,3); pm2 = 1.0d0
            case (4); pm2 = 1.0d0 / sqrt(4 * pi)
        end select

        gridglq(1:nlat, 1:nlong) = cilm(1,1,1) * pm2

        return

    end if

    !--------------------------------------------------------------------------
    !
    !   Determine Cilms, one l at a time I by integrating over all
    !   latitudes using Gauss-Legendre Quadrature. When PLX is not
    !   present, the Legendre functions are computed on the fly
    !   during the summations over l and m. These are scaled using
    !   the methodology of Holmesand Featherstone (2002), with the
    !   exception of the m=0 terms that do not need to be scaled
    !
    !--------------------------------------------------------------------------
    call dfftw_plan_dft_1d(plan, nlong, coef(1:nlong), grid(1:nlong), &
                           FFTW_BACKWARD, FFTW_MEASURE)  ! create generic plan

    if (present(plx)) then
        do i = 1, nlat
            coef = dcmplx(0.0d0,0.0d0)
            ! This summation order is intended to add the smallest terms first

            do l = lmax_comp, 0, -1
                l1 = l + 1
                k =  (l1 * l) / 2 + 1   ! m=0
                coef(1) = coef(1) + cilm(1,l1,1) * plx(i,k)

                do m = 1, l, 1
                    m1 = m + 1
                    k = (l1 * l) / 2 + m1
                    coef(m1) = coef(m1) + cilm(1,l1,m1) * plx(i,k)
                    coef(nlong-(m-1)) = coef(nlong-(m-1)) + &
                                        cilm(2,l1,m1) * plx(i,k) * &
                                        ((-1)**mod(m,2))
                end do

            end do

            call dfftw_execute(plan)    ! take fourier transform
            gridglq(i,1:nlong) = grid(1:nlong)

        end do

    else
        do i = 1, (nlat+1) / 2
            coef = dcmplx(0.0d0,0.0d0)

            if ( i == (nlat+1)/2 .and. mod(nlat,2) /= 0) then
                ! This latitude is the equator; z=0, u=1
                u = 1.0d0

                select case (lnorm)
                    case (1,2,3); pm2 = 1.0d0
                    case (4); pm2 = 1.0d0 / sqrt(4 * pi)
                end select

                coef(1) = coef(1) + cilm(1,1,1) * pm2

                do l = 2, lmax_comp, 2
                    l1 = l + 1
                    p = - ff2(l1,1) * pm2
                    pm2 = p
                    coef(1) = coef(1) + cilm(1,l1,1) * p
                end do

                select case (lnorm)
                    case (1,2,3); pmm = scalef
                    case (4); pmm = scalef / sqrt(4 * pi)
                end select

                rescalem = 1.0d0 / scalef

                do m = 1, lmax_comp-1, 1
                    m1 = m + 1

                    select case (lnorm)
                        case (1,4)
                            pmm = phase * pmm * sqr(2*m+1) / sqr(2*m)
                            pm2 = pmm
                        case (2)
                            pmm = phase * pmm * sqr(2*m+1) / sqr(2*m)
                            pm2 = pmm / sqr(2*m+1)
                        case (3)
                            pmm = phase * pmm * (2*m-1)
                            pm2 = pmm
                    end select

                    coef(m1) = coef(m1) + cilm(1,m1,m1) * pm2
                    coef(nlong-(m-1)) = coef(nlong-(m-1)) + cilm(2,m1,m1) * pm2

                    do l = m + 2, lmax_comp, 2
                        l1 = l + 1
                        p = - ff2(l1,m1) * pm2
                        coef(m1) = coef(m1) + cilm(1,l1,m1) * p
                        coef(nlong-(m-1)) = coef(nlong-(m-1)) &
                                            + cilm(2,l1,m1) * p
                        pm2 = p
                    end do

                    coef(m1) = coef(m1) * rescalem
                    coef(nlong-(m-1)) = coef(nlong-(m-1)) * rescalem &
                                        * ((-1)**mod(m,2))

                end do

                select case(lnorm)
                    case(1,4)
                        pmm = phase * pmm * sqr(2*lmax_comp+1) &
                              / sqr(2*lmax_comp) * rescalem
                    case(2)
                        pmm = phase * pmm / sqr(2*lmax_comp) * rescalem
                    case(3)
                        pmm = phase * pmm * (2*lmax_comp-1) * rescalem
                end select

                coef(lmax_comp+1) = coef(lmax_comp+1) &
                                    + cilm(1,lmax_comp+1,lmax_comp+1) * pmm
                coef(nlong-(lmax_comp-1)) = coef(nlong-(lmax_comp-1)) + &
                                            cilm(2,lmax_comp+1,lmax_comp+1) * &
                                            pmm * ((-1)**mod(lmax_comp,2))

                call dfftw_execute(plan)    ! take fourier transform 
                gridglq(i,1:nlong) = grid(1:nlong)

            else
                z = zero(i)
                u = sqrt( (1.0d0-z) * (1.0d0+z) )

                i_s = nlat + 1 - i

                coefs = dcmplx(0.0d0,0.0d0)

                select case (lnorm)
                    case (1,2,3);    pm2 = 1.0d0
                    case (4);    pm2 = 1.0d0 / sqrt(4 * pi)
                end select

                coef(1) = coef(1) + cilm(1,1,1) * pm2
                coefs(1) = coefs(1) + cilm(1,1,1) * pm2
                ! fsymsign is always 1 for l=m=0

                pm1 =  ff1(2,1) * z * pm2
                coef(1) = coef(1) + cilm(1,2,1) * pm1
                coefs(1) = coefs(1) - cilm(1,2,1) * pm1 ! fsymsign = -1

                do l = 2, lmax_comp, 1
                    l1 = l + 1
                    p = ff1(l1,1) * z * pm1 - ff2(l1,1) * pm2
                    coef(1) = coef(1) + cilm(1,l1,1) * p
                    coefs(1) = coefs(1) + cilm(1,l1,1) * p * fsymsign(l1,1)
                    pm2 = pm1
                    pm1 = p
                end do

                select case(lnorm)
                    case(1,2,3); pmm = scalef
                    case(4); pmm = scalef / sqrt(4 * pi)
                end select

                rescalem = 1.0d0 / scalef

                do m = 1, lmax_comp-1, 1
                    m1 = m + 1
                    rescalem = rescalem * u

                    select case (lnorm)
                        case (1,4)
                            pmm = phase * pmm * sqr(2*m+1) / sqr(2*m)
                            pm2 = pmm
                        case (2)
                            pmm = phase * pmm * sqr(2*m+1) / sqr(2*m)
                            pm2 = pmm / sqr(2*m+1)
                        case (3)
                            pmm = phase * pmm * (2*m-1)
                            pm2 = pmm
                    end select

                    coef(m1) = coef(m1) + cilm(1,m1,m1) * pm2
                    coef(nlong-(m-1)) = coef(nlong-(m-1)) &
                                        + cilm(2,m1,m1) * pm2
                    coefs(m1) = coefs(m1) + cilm(1,m1,m1) * pm2
                    coefs(nlong-(m-1)) = coefs(nlong-(m-1)) &
                                         + cilm(2,m1,m1) * pm2
                    ! fsymsign = 1

                    pm1 = z * ff1(m1+1,m1) * pm2

                    coef(m1) = coef(m1) + cilm(1,m1+1,m1) * pm1
                    coef(nlong-(m-1)) = coef(nlong-(m-1)) &
                                        + cilm(2,m1+1,m1) * pm1
                    coefs(m1) = coefs(m1) - cilm(1,m1+1,m1) * pm1
                    coefs(nlong-(m-1)) = coefs(nlong-(m-1)) &
                                        - cilm(2,m1+1,m1) * pm1
                    ! fsymsign = -1

                    do l = m + 2, lmax_comp, 1
                        l1 = l + 1
                        p = z * ff1(l1,m1) * pm1 - ff2(l1,m1) * pm2
                        pm2 = pm1
                        pm1 = p
                        coef(m1) = coef(m1) + cilm(1,l1,m1) * p
                        coef(nlong-(m-1)) = coef(nlong-(m-1)) &
                                            + cilm(2,l1,m1) * p
                        coefs(m1) = coefs(m1) + cilm(1,l1,m1) &
                                    * p * fsymsign(l1,m1)
                        coefs(nlong-(m-1)) = coefs(nlong-(m-1)) &
                                             + cilm(2,l1,m1) * p * &
                                             fsymsign(l1,m1)
                    end do

                    coef(m1) = coef(m1) * rescalem
                    coef(nlong-(m-1)) = coef(nlong-(m-1)) * rescalem &
                                        * ((-1)**mod(m,2))
                    coefs(m1) = coefs(m1) * rescalem
                    coefs(nlong-(m-1)) = coefs(nlong-(m-1)) * rescalem &
                                         * ((-1)**mod(m,2))

                end do

                rescalem = rescalem * u

                select case(lnorm)
                    case(1,4)
                        pmm = phase * pmm * sqr(2*lmax_comp+1) &
                              / sqr(2*lmax_comp) * rescalem
                    case(2)
                        pmm = phase * pmm / sqr(2*lmax_comp) * rescalem
                    case(3)
                        pmm = phase * pmm * (2*lmax_comp-1) * rescalem
                end select

                coef(lmax_comp+1) = coef(lmax_comp+1) &
                                    + cilm(1,lmax_comp+1,lmax_comp+1) * pmm
                coef(nlong-(lmax_comp-1)) = coef(nlong-(lmax_comp-1)) &
                                            + cilm(2,lmax_comp+1,lmax_comp+1) &
                                            * pmm * ((-1)**mod(lmax_comp,2))
                coefs(lmax_comp+1) = coefs(lmax_comp+1) &
                                     + cilm(1,lmax_comp+1,lmax_comp+1) * pmm
                coefs(nlong-(lmax_comp-1)) = coefs(nlong-(lmax_comp-1)) &
                                             + cilm(2,lmax_comp+1,lmax_comp+1) &
                                             * pmm * ((-1)**mod(lmax_comp,2))
                ! fsymsign = 1

                call dfftw_execute(plan)    ! take fourier transform
                gridglq(i,1:nlong) = grid(1:nlong)

                coef(1:nlong) = coefs(1:nlong)

                call dfftw_execute(plan)    ! take fourier transform
                gridglq(i_s,1:nlong) = grid(1:nlong)

            end if

        end do

    end if

    call dfftw_destroy_plan(plan)

end subroutine MakeGridGLQC
