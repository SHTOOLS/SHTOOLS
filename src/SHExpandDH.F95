subroutine SHExpandDH(grid, n, cilm, lmax, norm, sampling, csphase, &
                      lmax_calc, exitstatus)
!------------------------------------------------------------------------------
!
!   This routine will expand a grid containing n samples in both longitude
!   and latitude (or n x 2n, see below) into spherical harmonics. This routine
!   makes use of the sampling theorem presented in driscoll and healy (1994)
!   and employs ffts when calculating the sin and cos terms. The number of
!   samples, n, must be even for this routine to work, and the spherical
!   harmonic expansion is exact if the function is bandlimited to degree n/2-1.
!   Legendre functions are computed on the fly using the scaling methodology
!   presented in holmes and featherston (2002). When norm is 1,2 or 4, these are
!   accurate to about degree 2800. When norm is 3, the routine is only stable
!   to about degree 15. If the optional parameter lmax_calc is specified, the
!   spherical harmonic coefficients will only be calculated up to this degree.
!
!   If sampling is 1 (default), the input grid contains n samples in latitude
!   from 90 to -90+interval, and n samples in longitude from 0 to
!   360-2*interval, where interval is the latitudinal sampling interval 180/n.
!   Note that the datum at 90 degees north latitude is ultimately downweighted
!   to zero, so this point does not contribute to the spherical harmonic
!   coefficients. If sampling is 2, the input grid must contain n samples in
!   latitude and 2n samples in longitude. In this case, the sampling intervals
!   in latitude and longitude are 180/n and 360/n respectively. when performing
!   the ffts in longitude, the frequencies greater than n/2-1 are simply
!   discarded to prevent aliasing.
!
!   Calling Parameters
!
!       IN
!           grid        Equally sampled grid in latitude and longitude of
!                       dimension (1:n, 1:n) or and equally spaced grid of
!                       dimension (1:n,2n).
!           n           Number of samples in latitude and longitude (for
!                       sampling=1), or the number of samples in latitude (for
!                       sampling=2).
!
!       OUT
!           cilm        Array of spherical harmonic coefficients with dimension
!                       (2, lmax+1, lmax+1), or, if lmax_calc is present
!                       (2, lmax_calc+1, lmax_calc+1).
!           lmax        Spherical harmonic bandwidth of the grid. This
!                       corresponds to the maximum spherical harmonic degree of
!                       the expansion if the optional parameter lmax_calc is not
!                       specified.
!
!       OPTIONAL (IN)
!           norm        Normalization to be used when calculating legendre
!                       functions
!                           (1) "geodesy" (default)
!                           (2) schmidt
!                           (3) unnormalized
!                           (4) orthonormalized
!           sampling    (1) Grid is n latitudes by n longitudes (default).
!                       (2) Grid is n by 2n. the higher frequencies resulting
!                       from this oversampling are discarded, and hence not
!                       aliased into lower frequencies.
!           csphase     1: Do not include the condon-shortley phase factor of
!                       (-1)^m. -1: Apply the condon-shortley phase factor of
!                       (-1)^m.
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
!           are presented in driscoll and heally (1994).
!       2.  Use of a n by 2n grid is implemented because many geographic grids
!           are sampled this way. when taking the fourier transforms in
!           longitude, all of the higher frequencies are ultimately discarded.
!           If, instead, every other column of the grid were discarded to form
!           a nxn grid, higher frequencies could be aliased into lower
!           frequencies.
!
!   Dependencies:       dhaj, fftw3, csphase_default
!
!   Copyright (c) 2016, SHTOOLS
!   All rights reserved.
!
!------------------------------------------------------------------------------
    use FFTW3
    use SHTOOLS, only: dhaj, csphase_default
#ifdef FFTW3_UNDERSCORE
#define dfftw_plan_dft_r2c_1d dfftw_plan_dft_r2c_1d_
#define dfftw_execute dfftw_execute_
#define dfftw_destroy_plan dfftw_destroy_plan_
#endif

    implicit none

    real*8, intent(in) :: grid(:,:)
    real*8, intent(out) :: cilm(:,:,:)
    integer, intent(in) :: n
    integer, intent(out) :: lmax
    integer, intent(in), optional :: norm, sampling, csphase, lmax_calc
    integer, intent(out), optional :: exitstatus
    complex*16 :: cc(n+1)
    integer :: l, m, i, l1, m1, i_eq, i_s, lnorm, astat(4), lmax_comp, nlong
    integer*8 :: plan
    real*8 :: pi, gridl(2*n), aj(n), fcoef1(2, n/2+1), fcoef2(2, n/2+1), &
              theta, prod, scalef, rescalem, u, p, pmm, pm1, pm2, z, &
              ffc(1:2,-1:1)
    real*8, save, allocatable :: sqr(:), ff1(:,:), ff2(:,:)
    integer*1, save, allocatable :: fsymsign(:,:)
    integer, save :: lmax_old = 0, norm_old = 0
    integer :: phase
    external :: dfftw_plan_dft_r2c_1d, dfftw_execute, dfftw_destroy_plan

!$OMP   threadprivate(sqr, ff1, ff2, fsymsign, lmax_old, norm_old)

    if (present(exitstatus)) exitstatus = 0

    lmax = n / 2 - 1

    if (present(lmax_calc)) then
        if (lmax_calc > lmax) then
            print*, "Error --- SHExpandDH"
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
            print*, "Error --- SHExpandDH"
            print*, "Optional parameter sampling must be " // &
                    "1 (N by N) or 2 (N by 2N)."
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
        print*, "Error --- SHExpandDH"
        print*, "The number of samples in latitude and longitude, " // &
                "N, must be even."
        print*, "Input value is ", n
        if (present(exitstatus)) then
            exitstatus = 2
            return
        else
            stop
        end if

    else if (size(cilm(:,1,1)) < 2 .or. size(cilm(1,:,1)) < lmax_comp+1 .or. &
             size(cilm(1,1,:)) < lmax_comp+1) then
        print*, "Error --- SHExpandDH"
        print*, "CILM must be dimensioned as (2, LMAX_COMP+1, " // &
                "LMAX_COMP+1) where LMAX_COMP = MIN(N/2, LMAX_CALC+1)"
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
                print*, "Error --- SHExpandDH"
                print*, "GRIDDH must be dimensioned as (N, N) where N is ", n
                print*, "Input dimension is ", size(grid(:,1)), size(grid(1,:))
                if (present(exitstatus)) then
                    exitstatus = 1
                    return
                else
                    stop
                end if

            end if

        elseif (sampling == 2) then
            if (size(grid(:,1)) < n .or. size(grid(1,:)) < 2*n) then
                print*, "Error --- SHExpandDH"
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
            print*, "Error --- SHExpandDH"
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
            print*, "Error --- SHExpandDH"
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
        phase = csphase_default

    end if

    if (present(norm)) then
        if (norm > 4 .or. norm < 1) then
            print*, "Error --- SHExpandDH"
            print*, "Parameter norm must be 1 (geodesy), 2 (schmidt), " // &
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

    pi = acos(-1.0d0)

    cilm = 0.0d0

    scalef = 1.0d-280

    if (present(exitstatus)) then
        call DHaj(n, aj, exitstatus=exitstatus)
        if (exitstatus /= 0) return
    else
        call DHaj(n, aj)
    endif

    aj(1:n) = aj(1:n) * sqrt(4.0d0*pi)
    ! Driscoll and Heally use unity normalized spherical harmonics

    if (present(sampling)) then
        if (sampling == 1) then
            nlong = n

        else
            nlong = 2 * n

        endif

    else
        nlong = n

    end if

    !--------------------------------------------------------------------------
    !
    !   Calculate recursion constants used in computing the legendre polynomials
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
            print*, "Error --- SHExpandDH"
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
        !   Calculate signs used for symmetry of legendre functions about
        !   equator
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
        do l = 1, 2 * lmax_comp + 1
            sqr(l) = sqrt(dble(l))
        end do

        !----------------------------------------------------------------------
        !
        !   Precompute multiplicative factors used in recursion relationships
        !       p(l,m) = x*f1(l,m)*p(l-1,m) - p(l-2,m)*f2(l,m)
        !       k = l*(l+1)/2 + m + 1
        !   Note that prefactors are not used for the case when m=l as a
        !   different recursion is used. furthermore, for m=l-1, plmbar(l-2,m)
        !   is assumed to be zero.
        !
        !----------------------------------------------------------------------
        select case (lnorm)
            case(1,4)
                if (lmax_comp /= 0) then
                    ff1(2,1) = sqr(3)
                    ff2(2,1) = 0.0d0

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
                    ff2(l+1,l) = 0.0d0

                end do

            case(2)
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

                    end do

                    ff1(l+1,l)= dble(2*l-1) / sqr(l+m) / sqr(l-m)
                    ff2(l+1,l) = 0.0d0

                end do

            case(3)
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
    !   Create generic plan for grid
    !
    !--------------------------------------------------------------------------

    call dfftw_plan_dft_r2c_1d(plan, nlong, gridl(1:nlong), cc, fftw_measure)

    !--------------------------------------------------------------------------
    !
    !   Integrate over all latitudes. Take into account symmetry of the
    !   plms about the equator.
    !
    !--------------------------------------------------------------------------
    i_eq = n / 2 + 1  ! index correspondong to the equator

    do i = 2, i_eq - 1, 1
        theta = (i-1) * pi / dble(n)
        z = cos(theta)
        u = sqrt( (1.0d0-z) * (1.0d0+z) )

        gridl(1:nlong) = grid(i,1:nlong)
        call dfftw_execute(plan)    ! take fourier transform
        fcoef1(1,1:n/2) = sqrt(2*pi) * aj(i) * dble(cc(1:n/2)) / dble(nlong)
        fcoef1(2,1:n/2) = -sqrt(2*pi) * aj(i) * dimag(cc(1:n/2)) / dble(nlong)

        i_s = 2 * i_eq - i

        gridl(1:nlong) = grid(i_s,1:nlong)
        call dfftw_execute(plan)    ! take fourier transform
        fcoef2(1,1:n/2) = sqrt(2*pi) * aj(i_s) * dble(cc(1:n/2)) / dble(nlong)
        fcoef2(2,1:n/2) = -sqrt(2*pi) * aj(i_s) * dimag(cc(1:n/2)) / &
                          dble(nlong)

        select case (lnorm)
            case (1,2,3);    pm2 = 1.0d0
            case (4);        pm2 = 1.0d0 / sqrt(4*pi)

        end select

        cilm(1,1,1) = cilm(1,1,1) + pm2 * (fcoef1(1,1) + fcoef2(1,1))
        ! fsymsign = 1

        if (lmax_comp == 0) cycle

        pm1 = ff1(2,1) * z * pm2
        cilm(1,2,1) = cilm(1,2,1) + pm1 * (fcoef1(1,1) - fcoef2(1,1))
        ! fsymsign = -1

        ffc(1,-1) = fcoef1(1,1) - fcoef2(1,1)
        ffc(1, 1) = fcoef1(1,1) + fcoef2(1,1)

        do l = 2, lmax_comp, 1
            l1 = l + 1
            p = ff1(l1,1) * z * pm1 - ff2(l1,1) * pm2
            pm2 = pm1
            pm1 = p
            cilm(1,l1,1) = cilm(1,l1,1) + p * ffc(1,fsymsign(l1,1))

        end do

        select case (lnorm)
            case (1,2);  pmm = sqr(2) * scalef
            case (3);    pmm = scalef
            case (4);    pmm = sqr(2) * scalef / sqrt(4*pi)

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

            fcoef1(1:2,m1) = fcoef1(1:2,m1) * rescalem
            fcoef2(1:2,m1) = fcoef2(1:2,m1) * rescalem

            cilm(1:2,m1,m1) = cilm(1:2,m1,m1) + pm2 * &
                              (fcoef1(1:2,m1) + fcoef2(1:2,m1))
            ! fsymsign = 1

            pm1 = z * ff1(m1+1,m1) * pm2

            cilm(1:2,m1+1,m1) = cilm(1:2,m1+1,m1) + pm1 * &
                                (fcoef1(1:2,m1) - fcoef2(1:2,m1))
            ! fsymsign = -1

            ffc(1:2,-1) = fcoef1(1:2,m1) - fcoef2(1:2,m1)
            ffc(1:2, 1) = fcoef1(1:2,m1) + fcoef2(1:2,m1)

            do l = m + 2, lmax_comp, 1
                l1 = l + 1
                p = z * ff1(l1,m1) * pm1-ff2(l1,m1) * pm2
                pm2 = pm1
                pm1 = p
                cilm(1:2,l1,m1) = cilm(1:2,l1,m1) + p * &
                                  ffc(1:2,fsymsign(l1,m1))

            end do

        end do

        rescalem = rescalem * u

        select case (lnorm)
            case (1,4)
                pmm = phase * pmm * sqr(2*lmax_comp+1) &
                      / sqr(2*lmax_comp) * rescalem
            case (2)
                pmm = phase * pmm / sqr(2*lmax_comp) * rescalem
            case (3)
                pmm = phase * pmm * (2*lmax_comp-1) * rescalem

        end select

        cilm(1:2,lmax_comp+1,lmax_comp+1) = cilm(1:2,lmax_comp+1,lmax_comp+1) &
                                            + pmm * ( fcoef1(1:2,lmax_comp+1) &
                                            + fcoef2(1:2,lmax_comp+1) )
                        ! fsymsign = 1
    end do

    ! finally, do equator
    i = i_eq

    z = 0.0d0
    u = 1.0d0

    gridl(1:nlong) = grid(i,1:nlong)
    call dfftw_execute(plan)    ! take fourier transform
    fcoef1(1,1:n/2) = sqrt(2*pi) * aj(i) * dble(cc(1:n/2)) / dble(nlong)
    fcoef1(2,1:n/2) = -sqrt(2*pi) * aj(i) * dimag(cc(1:n/2)) / dble(nlong)

    select case (lnorm)
        case (1,2,3); pm2 = 1.0d0
        case (4);     pm2 = 1.0d0 / sqrt(4*pi)
    end select

    cilm(1,1,1) = cilm(1,1,1) + pm2 * fcoef1(1,1)

    if (lmax_comp /= 0) then
        do l = 2, lmax_comp, 2
            l1 = l + 1
            p = - ff2(l1,1) * pm2
            pm2 = p
            cilm(1,l1,1) = cilm(1,l1,1) + p * fcoef1(1,1)

        end do

        select case (lnorm)
            case (1,2);  pmm = sqr(2) * scalef
            case (3);    pmm = scalef
            case (4);    pmm = sqr(2) * scalef / sqrt(4 * pi)
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

            fcoef1(1:2,m1) = fcoef1(1:2,m1) * rescalem

            cilm(1:2,m1,m1) = cilm(1:2,m1,m1) + pm2 * fcoef1(1:2,m1)

            do l = m + 2, lmax_comp, 2
                l1 = l + 1
                p = - ff2(l1,m1) * pm2
                pm2 = p
                cilm(1:2,l1,m1) = cilm(1:2,l1,m1) + p * fcoef1(1:2,m1)

            end do

        end do

            select case (lnorm)
                    case(1,4)
                        pmm = phase * pmm * sqr(2*lmax_comp+1) &
                              / sqr(2*lmax_comp) * rescalem
                    case(2)
                        pmm = phase * pmm / sqr(2*lmax_comp) * rescalem
                    case(3)
                        pmm = phase * pmm * (2*lmax_comp-1) * rescalem
            end select

            cilm(1:2,lmax_comp+1,lmax_comp+1) = &
                            cilm(1:2,lmax_comp+1,lmax_comp+1) &
                            + pmm * fcoef1(1:2,lmax_comp+1)

        end if

    call dfftw_destroy_plan(plan)

    !--------------------------------------------------------------------------
    !
    !   Divide by integral of Ylm*Ylm 
    !
    !--------------------------------------------------------------------------
    select case(lnorm)
        case(1)
            do l = 0, lmax_comp, 1
                cilm(1:2,l+1, 1:l+1) = cilm(1:2,l+1, 1:l+1) / (4*pi)
            end do

        case (2)
            do l = 0, lmax_comp, 1
                cilm(1:2,l+1, 1:l+1) = cilm(1:2,l+1, 1:l+1) * (2*l+1) / (4*pi)
            end do

        case(3)
            do l = 0, lmax_comp, 1
                prod = 4 * pi / dble(2*l+1)
                cilm(1,l+1,1) = cilm(1,l+1,1) / prod
                prod = prod / 2.0d0

                do m = 1, l-1, 1
                    prod = prod * (l+m) * (l-m+1)
                    cilm(1:2,l+1,m+1) = cilm(1:2,l+1,m+1) / prod

                enddo

                !do m=l case
                if (l /= 0) cilm(1:2,l+1,l+1) = cilm(1:2,l+1, l+1) / (prod*2*l)

            end do

    end select

end subroutine SHExpandDH
