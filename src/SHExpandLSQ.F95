subroutine SHExpandLSQ(cilm, d, lat, lon, nmax, lmax, norm, chi2, &
                       csphase, exitstatus)
!------------------------------------------------------------------------------
!
!   This subroutine will expand a set of discrete data points into
!   spherical harmonics using a least squares inversion. When there are
!   more data points than spherical harmonic coefficients (nmax > (lmax+1)**2)
!   the solution of the overdetermined system will be determined by least
!   squares. If there are more coefficients than data points, then the solution
!   of the underdetermined system will be determined by minimizing the solution
!   norm. (See LAPACK DGELS documentation).
!
!   The default normalization convention for the output spherical harmonics
!   (and the calculation of the matrix G) is the "geodesy" normalization, though
!   this can be modified by supplying the optional argument norm.
!
!   Note that this routine takes lots of memory (~8*nmax*(lmax+1)**2 bytes) and
!   is very slow for large lmax
!
!   Calling Parameters
!
!       IN
!           d       Vector of length nmax of the raw data points.
!           lat     Vector of length nmax of the corresponding latitude points
!                   (in degrees).
!           lon     Vector of length nmax of the corresponding longitude points
!                   (in degrees).
!           nmax    Number of data points.
!           lmax    Maximum degree of spherical harmonic expansion.
!
!       OUT 
!           cilm    Spherical harmonic coefficients.
!
!       OPTIONAL (IN)
!           norm    Spherical harmonic normalizaton for output coefficients and
!                   calculation of matrix G:
!                       1. PlmBar (geodesy)
!                       2. PlmSchmidt
!                       3. PLegendreA (unnormalized)
!                       4. PlmBar/sqrt(4 pi) (orthonormalized)
!           csphase:    1: Do not include the phase factor of (-1)^m
!                       -1: Apply the phase factor of (-1)^m.
!
!       OPTIONAL (OUT)
!           chi2        This is the residual sum of squares misfit for an
!                       overdetermined inversion.
!           exitstatus  If present, instead of executing a STOP when an error
!                       is encountered, the variable exitstatus will be
!                       returned describing the error.
!                       0 = No errors;
!                       1 = Improper dimensions of input array;
!                       2 = Improper bounds for input variable;
!                       3 = Error allocating memory;
!                       4 = File IO error.
!
!   Dependencies:   LAPACK, PlmBar, PLegendreA, PlmSchmidt, PlmON,
!                   CSPHASE_DEFAULT
!
!   Copyright (c) 2016, SHTOOLS
!   All rights reserved.
!
!------------------------------------------------------------------------------
    use SHTOOLS, only: PlmBar, PLegendreA, PlmSchmidt, PlmON, CSPHASE_DEFAULT

    implicit none

    real*8, intent(in) :: d(:), lat(:), lon(:)
    real*8, intent(out) :: cilm(:,:,:)
    integer, intent(in) :: nmax, lmax
    integer, intent(in), optional :: norm, csphase
    real*8, intent(out), optional :: chi2
    integer, intent(out), optional :: exitstatus
    integer, parameter :: opt = 80
    integer :: ncoef, i, l, m, ind1, ind2, info, lwork, opt1, phase, &
               astat(4), lnorm
    real*8 :: pi, lonr
    real*8, allocatable :: mm(:), gg(:, :), p(:), work(:)
#ifdef LAPACK_UNDERSCORE
#define dgels dgels_
#endif
    external :: dgels

    if (present(exitstatus)) exitstatus = 0

    if (size(cilm(:,1,1)) < 2 .or. size(cilm(1,:,1)) < lmax+1 &
            .or. size(cilm(1,1,:)) < lmax+1) then
        print*, "Error --- SHExpandLSQ"
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

    else if (size(d) < nmax) then
        print*, "Error --- SHExpandLSQ"
        print*, "D must be dimensioned as (NMAX) where NMAX is ", nmax
        print*, "Input array is dimensioned ", size(d)
        if (present(exitstatus)) then
            exitstatus = 1
            return
        else
            stop
        end if

    else if (size(lat) < nmax) then
        print*, "Error --- SHExpandLSQ"
        print*, "LAT must be dimensioned as (NMAX) where NMAX is ", nmax
        print*, "Input array is dimensioned ", size(lat)
        if (present(exitstatus)) then
            exitstatus = 1
            return
        else
            stop
        end if

    else if (size(lon) < nmax) then
        print*, "Error --- SHExpandLSQ"
        print*, "LON must be dimensioned as (NMAX) where NMAX is ", nmax
        print*, "Input array is dimensioned ", size(lon)
        if (present(exitstatus)) then
            exitstatus = 1
            return
        else
            stop
        end if

    end if

    if (present(norm)) then
        if (norm > 4 .or. norm < 1) then
            print*, "Error - SHExpandLSQ"
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
            print*, "Error ---- SHExpandLSQ"
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
            phase = CSPHASE_DEFAULT

    end if

    allocate (mm(max((lmax+1)**2, nmax)), stat = astat(1))
    allocate (gg(nmax, (lmax+1)**2), stat = astat(2))
    allocate (p((lmax+1)*(lmax+2)/2), stat = astat(3))
    allocate (work(min((lmax+1)**2, nmax)*(1+opt)), stat = astat(4))

    if (astat(1) /= 0 .or. astat(2) /= 0 .or. astat(3) /= 0 &
            .or. astat(4) /= 0) then
        print*, "Error --- SHExpandLSQ"
        print*, "Problem allocating arrays MM, GG, P, and WORK", astat(1), &
                astat(2), astat(3), astat(4)
        if (present(exitstatus)) then
            exitstatus = 3
            return
        else
            stop
        end if

    end if

    lwork = min((lmax+1)**2, nmax)*(1+opt)
    pi = acos(-1.0d0)
    mm = 0.0d0
    gg = 0.0d0

    ncoef = (lmax+1)**2

    !--------------------------------------------------------------------------
    !
    !   Calculate matrix G (nmax by ncoef)
    !
    !--------------------------------------------------------------------------
    do i = 1, nmax
        if (present(exitstatus)) then
            select case(lnorm)
                case(1)
                    call PlmBar(p, lmax, sin(lat(i) * pi / 180.0d0), &
                                csphase = phase, exitstatus = exitstatus)
                case(2)
                    call PlmSchmidt(p, lmax, sin(lat(i) * pi / 180.0d0), &
                                    csphase = phase, exitstatus = exitstatus)
                case(3)
                    call PLegendreA(p, lmax, sin(lat(i) * pi / 180.0d0), &
                                    csphase = phase, exitstatus = exitstatus)
                case(4)
                    call PlmON(p, lmax, sin(lat(i) * pi / 180.0d0), &
                               csphase = phase, exitstatus = exitstatus)
            end select
            if (exitstatus /= 0) return

        else
            select case(lnorm)
                case(1)
                    call PlmBar(p, lmax, sin(lat(i) * pi / 180.0d0), &
                                csphase = phase)
                case(2)
                    call PlmSchmidt(p, lmax, sin(lat(i) * pi / 180.0d0), &
                                    csphase = phase)
                case(3)
                    call PLegendreA(p, lmax, sin(lat(i) * pi / 180.0d0), &
                                    csphase = phase)
                case(4)
                    call PlmON(p, lmax, sin(lat(i) * pi / 180.0d0), &
                               csphase = phase)
            end select

        end if

        lonr = lon(i) * pi / 180.0d0
        ind1 = 0

        do l = 0, lmax
            ! do cos terms

            do m = 0, l
                ind1 = ind1 + 1
                ind2 = l * (l + 1) / 2 + m + 1
                gg(i,ind1) = p(ind2) * cos(m*lonr)

            end do

            ! do sin terms
            do m = 1, l, 1
                ind1 = ind1 + 1
                ind2 = l*(l+1)/2 + m + 1
                gg(i,ind1) = p(ind2) * sin(m*lonr)

            end do

        end do

    end do

    mm(1:nmax) = d(1:nmax)

    !--------------------------------------------------------------------------
    !
    !   Do least squares inversion, i.e.,
    !   m = [G' G]^-1 G' d
    !   using LAPACK routine DGELS.
    !
    !--------------------------------------------------------------------------
    call dgels('N', nmax, ncoef, 1, gg, nmax, mm, max((lmax+1)**2, nmax), &
                work, lwork, info)

    if (info /= 0) then
        print*, "Error --- SHExpandLSQ"
        print*, "DGELS: Problem performing least squares inversion."
        print*, "DGELS INFO = ", info
        if (present(exitstatus)) then
            exitstatus = 5
            return
        else
            stop
        end if

    end if

    if (work(1) >  dble(lwork)) then
        opt1 = int(work(1) / min((lmax+1)**2, nmax) - 1)
        print*, "Warning --- SHExpandLSQ"
        print*, "Consider changing parameter value of OPT to ", opt1, &
                " and recompiling the SHTOOLS archive."
    end if

    !--------------------------------------------------------------------------
    !
    !   Convert mm into cilm
    !
    !--------------------------------------------------------------------------
    ind1 = 0

    do l = 0, lmax
        ! do cos terms
        do m = 0, l
            ind1 = ind1 + 1
            cilm(1,l+1, m+1) = mm(ind1)

        end do

        ! do sin terms
        do m = 1, l, 1
            ind1 = ind1 + 1
            cilm(2, l+1, m+1) = mm(ind1)

        end do

    end do

    !--------------------------------------------------------------------------
    !
    !   Compute residual sum of sqaures misfit for the overdetermined case.
    !
    !--------------------------------------------------------------------------
    if (present(chi2) .and. nmax > ncoef) then
        chi2 = 0.0d0

        do i = ncoef + 1, nmax
            chi2 = chi2 + mm(i)**2

        end do

    end if

    ! deallocate memory
    select case(lnorm)
        case(1)
            call PlmBar(p, -1, sin(lat(1) * pi / 180.0d0), csphase = phase)
        case(2)
            call PlmSchmidt(p, -1, sin(lat(1) * pi / 180.0d0), csphase = phase)
        case(4)
            call PlmON(p, -1, sin(lat(1) * pi / 180.0d0), csphase = phase)

    end select

    deallocate (mm)
    deallocate (gg)
    deallocate (p)
    deallocate (work)

end subroutine SHExpandLSQ
