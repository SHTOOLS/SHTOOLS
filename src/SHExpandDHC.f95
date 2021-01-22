subroutine SHExpandDHC(grid, n, cilm, lmax, norm, sampling, csphase, &
                       lmax_calc, exitstatus)
!------------------------------------------------------------------------------
!
!   This routine will expand a grid containing N samples in both longitude
!   and latitude (or N x 2N, see below) into spherical harmonics. This routine
!   makes use of the sampling theorem presented in Driscoll and Healy (1994)
!   and employs FFTs when calculating the exponential terms. The number of
!   samples, N, must be EVEN for this routine to work, and the spherical
!   harmonic expansion is exact if the function is bandlimited to degree N/2-1.
!   Legendre functions are computed on the fly using the scaling methodology
!   presented in Holmes and Featherston (2002). When NORM is 1, 2 or 4, these
!   are accurate to about degree 2800. When NORM is 3, the routine is only
!   stable to about degree 15! If the optional parameter LMAX_CALC is
!   specified, the spherical harmonic coefficients will only be calculated up
!   to this degree.
!
!   If SAMPLING is 1 (default), the input grid contains N samples in latitude
!   from 90 to -90+interval, and N samples in longitude from 0 to
!   360-2*interval, where interval is the latitudinal sampling interval 180/N.
!   Note that the datum at 90 degees North latitude is ultimately downweighted
!   to zero, so this point does not contribute to the spherical harmonic
!   coefficients. If SAMPLING is 2, the input grid must contain N samples in
!   latitude and 2N samples in longitude. In this case, the sampling intervals
!   in latitude and longitude are 180/N and 360/N respectively. When performing
!   the FFTs in longitude, the frequencies greater than N/2-1 are simply
!   discarded to prevent aliasing.
!
!   The complex spherical harmonics are output in the array cilm. Cilm(1,,)
!   contains the positive m term, wheras cilm(2,,) contains the negative m
!   term. The negative order Legendre functions are calculated making use of
!   the identity Y_{lm}^* = (-1)^m Y_{l,-m}.
!
!   Calling Parameters
!
!       IN
!           grid        Equally sampled grid in latitude and longitude of
!                       dimension (1:N, 1:N) or and equally spaced grid of
!                       dimension (1:N,2N).
!           N           Number of samples in latitude and longitude (for
!                       SAMPLING=1), or the number of samples in latitude (for
!                       SAMPLING=2).
!
!       OUT
!           cilm        Array of spherical harmonic coefficients with dimension
!                       (2, LMAX+1, LMAX+1), or, if LMAX_CALC is present
!                       (2, LMAX_CALC+1, LMAX_CALC+1).
!           lmax        Spherical harmonic bandwidth of the grid. This
!                       corresponds to the maximum spherical harmonic degree of
!                       the expansion if the optional parameter LMAX_CALC is
!                       not specified.
!
!       OPTIONAL (IN)
!           norm        Normalization to be used when calculating Legendre
!                       functions
!                           (1) "geodesy" (default)
!                           (2) Schmidt
!                           (3) unnormalized
!                           (4) orthonormalized
!           sampling    (1) Grid is N latitudes by N longitudes (default).
!                       (2) Grid is N by 2N. The higher frequencies resulting
!                       from this oversampling are discarded, and hence not
!                       aliased into lower frequencies.
!           csphase     1: Do not include the Condon-Shortley phase factor of
!                       (-1)^m (default). -1: Apply the Condon-Shortley phase
!                       factor of (-1)^m.
!           lmax_calc   The maximum spherical harmonic degree calculated in the
!                       spherical harmonic expansion.
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
!       1.  This routine does not use the fast legendre transforms that
!           are presented in Driscoll and Heally (1994).
!       2.  Use of a N by 2N grid is implemented because many geographic grids
!           are sampled this way. When taking the Fourier transforms in
!           longitude, all of the higher frequencies are ultimately discarded.
!           If, instead, every other column of the grid were discarded to form
!           a NxN grid, higher frequencies could be aliased into lower
!           frequencies.
!
!   Copyright (c) 2005-2019, SHTOOLS
!   All rights reserved.
!
!------------------------------------------------------------------------------
    use FFTW3
    use SHTOOLS, only: DHaj
    use ftypes
    use, intrinsic :: iso_c_binding

    implicit none

    complex(dp), intent(in) :: grid(:,:)
    complex(dp), intent(out) :: cilm(:,:,:)
    integer(int32), intent(in) :: n
    integer(int32), intent(out) :: lmax
    integer(int32), intent(in), optional :: norm, sampling, csphase, lmax_calc
    integer(int32), intent(out), optional :: exitstatus
    complex(dp) :: cc(2*n), ccs(2*n), gridl(2*n), gridls(2*n), fcoef1(2*n), &
                   fcoef2(2*n), ffc1(-1:1), ffc2(-1:1)
    integer(int32) :: l, m, i, l1, m1, i_eq, i_s, lnorm, astat(5), lmax_comp, &
                      nlong
    type(C_PTR) :: plan, plans
    real(dp) :: pi, aj(n), theta, prod, scalef, rescalem, u, p, pmm, pm1, &
                pm2, z
    real(dp), save, allocatable :: ff1(:,:), ff2(:,:), sqr(:)
    integer(int8), save, allocatable :: fsymsign(:,:)
    integer(int32), save :: lmax_old = 0, norm_old = 0
    integer(int32) :: phase

!$OMP   threadprivate(sqr, ff1, ff2, fsymsign, lmax_old, norm_old)

    if (present(exitstatus)) exitstatus = 0

    lmax = n / 2 - 1

    if (present(lmax_calc)) then
        if (lmax_calc > lmax) then
            print*, "Error --- SHExpandDHC"
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
            lmax_comp = min(lmax, lmax_calc)

        end if

    else
        lmax_comp = lmax

    end if

    if (present(sampling)) then
        if (sampling /= 1 .and. sampling /= 2) then
            print*, "Error --- SHExpandDHC"
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

    end if

    if (mod(n,2) /= 0) then
        print*, "Error --- SHExpandDHC"
        print*, "The number of samples in latitude and longitude, " // &
                "n, must be even."
        print*, "Input value is ", n
        if (present(exitstatus)) then
            exitstatus = 2
            return
        else
            stop
        end if

    else if (size(cilm(:,1,1)) < 2 .or. size(cilm(1,:,1)) < lmax_comp+1 .or. &
             size(cilm(1,1,:)) < lmax_comp+1) then
        print*, "Error --- SHExpandDHC"
        print*, "CILM must be dimensioned as (2, LMAX_COMP+1, LMAX_COMP+1)"
        print*, "where LMAX_COMP = MIN(N/2, LMAX_CALC+1)"
        print*, "N = ", n
        if (present(lmax_calc)) print*, "LMAX_CALC = ", lmax_calc
        print*, "Input dimension is ", size(cilm(:,1,1)), size(cilm(1,:,1)), &
                size(cilm(1,1,:))
        if (present(exitstatus)) then
            exitstatus = 1
            return
        else
            stop
        end if

    end if

    if (present(sampling)) then
        if (sampling == 1) then
            if (size(grid(:,1)) < n .or. size(grid(1,:)) < n) then
                print*, "Error --- SHExpandDHC"
                print*, "GRIDDH must be dimensioned as (N, N) where N is ", n
                print*, "Input dimension is ", size(grid(:,1)), size(grid(1,:))
                if (present(exitstatus)) then
                    exitstatus = 1
                    return
                else
                    stop
                end if

            end if

        else if (sampling == 2) then
            if (size(grid(:,1)) < n .or. size(grid(1,:)) < 2*n) then
                print*, "Error --- SHExpandDHC"
                print*, "GRIDDH must be dimensioned as (N, 2*N) where N is ", n
                print*, "Input dimension is ", size(grid(:,1)), size(grid(1,:))
                if (present(exitstatus)) then
                    exitstatus = 1
                    return
                else
                    stop
                end if

            end if

        end if

    else
        if (size(grid(:,1)) < n .or. size(grid(1,:)) < n) then
            print*, "Error --- SHExpandDHC"
            print*, "GRIDDH must be dimensioned as (N, N) where N is ", n
            print*, "Input dimension is ", size(grid(:,1)), size(grid(1,:))
            if (present(exitstatus)) then
                exitstatus = 1
                return
            else
                stop
            end if

        end if

    end if

    if (present(csphase)) then
        if (csphase /= -1 .and. csphase /= 1) then
            print*, "SHExpandDHC --- Error"
            print*, "CSPHASE must be 1 (exclude) or -1 (include)."
            print*, "Input value is ", csphase
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

    if (present(norm)) then
        if (norm > 4 .or. norm < 1) then
            print*, "Error --- SHExpandDHC"
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

    pi = acos(-1.0_dp)

    cilm = cmplx(0.0_dp, 0.0_dp, dp)

    scalef = 1.0e-280_dp

    if (present(exitstatus)) then
        call DHaj(n, aj, exitstatus=exitstatus)
        if (exitstatus /= 0) return
    else
        call DHaj(n, aj)
    end if

    aj(1:n) = aj(1:n) * sqrt(4.0_dp * pi)
    ! Driscoll and Heally use unity normalized spherical harmonics

    if (present(sampling)) then
        if (sampling == 1) then
            nlong = n

        else
            nlong = 2*n

        end if

    else
        nlong = n

    end if

    !--------------------------------------------------------------------------
    !
    !   Calculate recursion constants used in computing the Legendre functions.
    !
    !--------------------------------------------------------------------------
    if (lmax_comp /= lmax_old .or. lnorm /= norm_old) then
        if (allocated (sqr)) deallocate (sqr)
        if (allocated (ff1)) deallocate (ff1)
        if (allocated (ff2)) deallocate (ff2)
        if (allocated (fsymsign)) deallocate (fsymsign)

        allocate (sqr(2*lmax_comp+1), stat=astat(1))
        allocate (ff1(lmax_comp+1,lmax_comp+1), stat=astat(2))
        allocate (ff2(lmax_comp+1,lmax_comp+1), stat=astat(3))
        allocate (fsymsign(lmax_comp+1,lmax_comp+1), stat=astat(4))

        if (sum(astat(1:4)) /= 0) then
            print*, "Error --- SHExpandDHC"
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
                    ff2(2,1) = 0.0_dp

                end if

                do l = 2, lmax_comp, 1
                    ff1(l+1,1) = sqr(2*l-1) * sqr(2*l+1) / dble(l)
                    ff2(l+1,1) = dble(l-1) * sqr(2*l+1) / sqr(2*l-3) / dble(l)

                    do m = 1, l - 2, 1
                        ff1(l+1,m+1) = sqr(2*l+1) * sqr(2*l-1) / sqr(l+m) &
                                       / sqr(l-m)
                        ff2(l+1,m+1) = sqr(2*l+1) * sqr(l-m-1) * sqr(l+m-1) &
                                       / sqr(2*l-3) / sqr(l+m) / sqr(l-m)
                    end do

                    ff1(l+1,l) = sqr(2*l+1) * sqr(2*l-1) / sqr(l+m) / sqr(l-m)
                    ff2(l+1,l) = 0.0_dp

                end do

            case(2)
                if (lmax_comp /= 0) then
                    ff1(2,1) = 1.0_dp
                    ff2(2,1) = 0.0_dp

                end if

                do l = 2, lmax_comp, 1
                    ff1(l+1,1) = dble(2*l-1) / dble(l)
                    ff2(l+1,1) = dble(l-1) / dble(l)

                    do m = 1, l - 2, 1
                        ff1(l+1,m+1) = dble(2*l-1) / sqr(l+m) / sqr(l-m)
                        ff2(l+1,m+1) = sqr(l-m-1) * sqr(l+m-1) / sqr(l+m) &
                                       / sqr(l-m)

                    end do

                    ff1(l+1,l) = dble(2*l-1) / sqr(l+m) / sqr(l-m)
                    ff2(l+1,l) = 0.0_dp

                end do

            case(3)
                do l = 1, lmax_comp, 1
                    ff1(l+1,1) = dble(2*l-1) / dble(l)
                    ff2(l+1,1) = dble(l-1) / dble(l)

                    do m = 1, l - 1, 1
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
    !   Create generic plan for gridl and gridls.
    !
    !--------------------------------------------------------------------------
    plan = fftw_plan_dft_1d(nlong, gridl(1:nlong), cc(1:nlong), &
                            FFTW_FORWARD, FFTW_MEASURE)
    plans = fftw_plan_dft_1d(nlong, gridls(1:nlong), ccs(1:nlong), &
                             FFTW_FORWARD, FFTW_MEASURE)

    !--------------------------------------------------------------------------
    !
    !   Integrate over all latitudes. Take into account symmetry of the
    !   Plms about the equator.
    !
    !--------------------------------------------------------------------------
    i_eq = n / 2 + 1  ! Index correspondong to the equator

    ! First do equator
    i = i_eq

    z = 0.0_dp
    u = 1.0_dp

    gridl(1:nlong) = grid(i,1:nlong)
    call fftw_execute_dft(plan, gridl, cc)    ! take fourier transform
    fcoef1(1:nlong) = cc(1:nlong) * sqrt(2*pi) * aj(i) / dble(nlong)

    select case (lnorm)
        case (1,2,3); pm2 = 1.0_dp
        case (4);     pm2 = 1.0_dp / sqrt(4 * pi)
    end select

    cilm(1,1,1) = cilm(1,1,1) + pm2 * fcoef1(1)

    if (lmax_comp /= 0) then
        do l = 2, lmax_comp, 2
            l1 = l + 1
            p = - ff2(l1,1) * pm2
            pm2 = p
            cilm(1,l1,1) = cilm(1,l1,1) + p * fcoef1(1)

        end do

        select case (lnorm)
            case (1,2,3); pmm = scalef
            case (4);     pmm = scalef / sqrt(4 * pi)
        end select

        rescalem = 1.0_dp / scalef

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

            fcoef1(m1) = fcoef1(m1) * rescalem
            fcoef1(nlong-(m-1)) = fcoef1(nlong-(m-1)) * rescalem &
                                  * ((-1)**mod(m,2))

            cilm(1,m1,m1) = cilm(1,m1,m1) + pm2 * fcoef1(m1)
            cilm(2,m1,m1) = cilm(2,m1,m1) + pm2 * fcoef1(nlong-(m-1))

            do l = m + 2, lmax_comp, 2
                l1 = l + 1
                p = - ff2(l1,m1) * pm2
                pm2 = p
                cilm(1,l1,m1) = cilm(1,l1,m1) + p * fcoef1(m1)
                cilm(2,l1,m1) = cilm(2,l1,m1) + p * fcoef1(nlong-(m-1))

            end do

        end do

        select case (lnorm)
            case (1,4)
                pmm = phase * pmm * sqr(2*lmax_comp+1) &
                      / sqr(2*lmax_comp) * rescalem

            case(2)
                pmm = phase * pmm / sqr(2*lmax_comp) * rescalem

            case(3)
                pmm = phase * pmm * (2*lmax_comp-1) * rescalem

        end select

        cilm(1,lmax_comp+1,lmax_comp+1) = cilm(1,lmax_comp+1,lmax_comp+1) &
                                          + pmm * fcoef1(lmax_comp+1)
        cilm(2,lmax_comp+1,lmax_comp+1) = cilm(2,lmax_comp+1,lmax_comp+1) &
                                          + ((-1)**mod(lmax_comp,2)) * pmm &
                                          * fcoef1(nlong-(lmax_comp-1))

    end if

    do i = 2, i_eq - 1, 1
        theta = (i-1) * pi /dble(n)
        z = cos(theta)
        u = sqrt( (1.0_dp-z) * (1.0_dp+z) )

        gridl(1:nlong) = grid(i,1:nlong)
        call fftw_execute_dft(plan, gridl, cc)    ! take fourier transform
        fcoef1(1:nlong) = cc(1:nlong) * sqrt(2*pi) * aj(i) / dble(nlong)
        ! positive frequencies up to n/2,
        ! negative frequencies beyond in oposite order

        i_s = 2 * i_eq - i

        gridls(1:nlong) = grid(i_s,1:nlong)
        call fftw_execute_dft(plans, gridls, ccs)    ! take fourier transform
        fcoef2(1:nlong) = ccs(1:nlong) * sqrt(2*pi) * aj(i_s) / dble(nlong)

        select case (lnorm)
            case (1,2,3); pm2 = 1.0_dp
            case (4); pm2 = 1.0_dp / sqrt(4 * pi)
        end select

        cilm(1,1,1) = cilm(1,1,1) + pm2 * (fcoef1(1) + fcoef2(1))
        ! fsymsign = 1

        if (lmax_comp == 0) cycle

        pm1 = ff1(2,1) * z * pm2
        cilm(1,2,1) = cilm(1,2,1) + pm1 * ( fcoef1(1) - fcoef2(1) )
        ! fsymsign = -1

        ffc1(-1) = fcoef1(1) - fcoef2(1)
        ffc1(1) = fcoef1(1) + fcoef2(1)

        do l = 2, lmax_comp, 1
            l1 = l + 1
            p = ff1(l1,1) * z * pm1 - ff2(l1,1) * pm2
            pm2 = pm1
            pm1 = p
            cilm(1,l1,1) = cilm(1,l1,1) + p * ffc1(fsymsign(l1,1))

        end do

        select case (lnorm)
            case(1,2,3); pmm = scalef
            case(4);     pmm = scalef / sqrt(4*pi)
        end select

        rescalem = 1.0_dp / scalef

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

            fcoef1(m1) = fcoef1(m1) * rescalem
            fcoef1(nlong-(m-1)) = fcoef1(nlong-(m-1)) * rescalem &
                                  * ((-1)**mod(m,2))
                                  ! multiply by (-1)^m for P_{l,-m}
            fcoef2(m1) = fcoef2(m1) * rescalem
            fcoef2(nlong-(m-1)) = fcoef2(nlong-(m-1)) * rescalem &
                                  * ((-1)**mod(m,2))

            cilm(1,m1,m1) = cilm(1,m1,m1) + pm2 * ( fcoef1(m1) + fcoef2(m1) )
            cilm(2,m1,m1) = cilm(2,m1,m1) + pm2 * ( fcoef1(nlong-(m-1)) &
                            + fcoef2(nlong-(m-1)) )
                            ! fsymsign = 1

            pm1 = z * ff1(m1+1,m1) * pm2

            cilm(1,m1+1,m1) = cilm(1,m1+1,m1) + pm1 * (fcoef1(m1) - fcoef2(m1))
            cilm(2,m1+1,m1) = cilm(2,m1+1,m1) +  pm1 &
                              * ( fcoef1(nlong-(m-1)) - fcoef2(nlong-(m-1)) )
                              ! fsymsign = -1

            ffc1(-1) = fcoef1(m1) - fcoef2(m1)
            ffc1(1) = fcoef1(m1) + fcoef2(m1)
            ffc2(-1) = fcoef1(nlong-(m-1)) - fcoef2(nlong-(m-1))
            ffc2(1) = fcoef1(nlong-(m-1)) + fcoef2(nlong-(m-1))

            do l = m + 2, lmax_comp, 1
                l1 = l + 1
                p = z * ff1(l1,m1) * pm1-ff2(l1,m1) * pm2
                pm2 = pm1
                pm1 = p
                cilm(1,l1,m1) = cilm(1,l1,m1) + p * ffc1(fsymsign(l1,m1))
                cilm(2,l1,m1) = cilm(2,l1,m1) + p * ffc2(fsymsign(l1,m1))

            end do

        end do

        rescalem = rescalem * u

        select case(lnorm)
            case(1,4)
                pmm = phase * pmm * sqr(2*lmax_comp+1) &
                      / sqr(2*lmax_comp) * rescalem
            case(2);    pmm = phase * pmm / sqr(2*lmax_comp) * rescalem
            case(3);    pmm = phase * pmm * (2*lmax_comp-1) * rescalem
    
        end select

        cilm(1,lmax_comp+1,lmax_comp+1) = cilm(1,lmax_comp+1,lmax_comp+1) &
                                          + pmm * (fcoef1(lmax_comp+1) + &
                                          fcoef2(lmax_comp+1))
        cilm(2,lmax_comp+1,lmax_comp+1) = cilm(2,lmax_comp+1,lmax_comp+1) &
                                          + pmm * (fcoef1(nlong-(lmax_comp-1))&
                                          + fcoef2(nlong-(lmax_comp-1)))
                                          ! fsymsign = 1

    end do

    call fftw_destroy_plan(plan)
    call fftw_destroy_plan(plans)

    !--------------------------------------------------------------------------
    !
    !   Divide by integral of Ylm*Ylm 
    !
    !--------------------------------------------------------------------------
    select case (lnorm)
        case(1)
            do l = 0, lmax_comp, 1
                cilm(1:2,l+1, 1:l+1) = cilm(1:2,l+1, 1:l+1) / (4*pi)

            end do

        case(2)
            do l = 0, lmax_comp, 1
                cilm(1:2,l+1, 1:l+1) = cilm(1:2,l+1, 1:l+1) * (2*l+1) / (4*pi)

            end do

        case (3)
            do l = 0, lmax_comp, 1
                prod = 4 * pi / dble(2*l+1)
                cilm(1,l+1,1) = cilm(1,l+1,1) / prod

                do m = 1, l-1, 1
                    prod = prod * (l+m) * (l-m+1)
                    cilm(1:2,l+1,m+1) = cilm(1:2,l+1,m+1) / prod

                end do

                !do m=l case
                if (l /= 0) cilm(1:2,l+1,l+1) = cilm(1:2,l+1, l+1) / (prod*2*l)

            end do

    end select

end subroutine SHExpandDHC
