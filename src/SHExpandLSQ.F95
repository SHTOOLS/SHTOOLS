subroutine SHExpandLSQ(cilm, d, lat, lon, nmax, lmax, norm, chi2, &
                       csphase, weights, g, exitstatus)
!------------------------------------------------------------------------------
!
!   This subroutine will determine the spherical harmonic coefficients of a
!   function using a least squares inversion:
!
!       m = [G' W G]^-1 G' W d
!
!   When there are more data points than spherical harmonic coefficients
!   (nmax > (lmax+1)**2), the solution of the overdetermined system will be
!   determined by least squares. If there are more coefficients than data
!   points, then the solution of the underdetermined system will be determined
!   by minimizing the solution norm. The inversion is performed using the
!   LAPACK routine DGELS.
!
!   A weighted least squares inversion will be performed if the optional vector
!   weights is specified. In this case, the problem must be overdetermined, and
!   it is assumed that each measurement is statistically independent. The
!   inversion is performed using the LAPACK routine DGELS be scaling the data
!   vector and data kernel matrix G.
!
!   The default normalization convention for the output spherical harmonics
!   (and the calculation of the matrix G) is the "geodesy" 4pi normalization,
!   though this can be modified by supplying the optional argument norm.
!
!   Note that this routine requires a lot of memory (~8*nmax*(lmax+1)**2 bytes)
!   and is very slow for large lmax. If the subroutine is to be called several
!   times using the same lat and lon coordinates, the data kernel matrix G can
!   be precomputed using the routine LSQ_G.
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
!           cilm    The spherical harmonic coefficients.
!
!       OPTIONAL (IN)
!           norm    Spherical harmonic normalizaton for output coefficients and
!                   calculation of matrix G:
!                       1. PlmBar (geodesy)
!                       2. PlmSchmidt
!                       3. PLegendreA (unnormalized)
!                       4. PlmBar/sqrt(4 pi) (orthonormalized)
!           csphase     1: Do not include the phase factor of (-1)^m (default).
!                       -1: Apply the phase factor of (-1)^m.
!           weights  The weights to be used when performing the weighted least
!                    squares inversion.
!           G        The data kernel matrix used in the least squares inversion.
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
!   Copyright (c) 2005-2019, SHTOOLS
!   All rights reserved.
!
!------------------------------------------------------------------------------
    use SHTOOLS, only: PlmBar, PLegendreA, PlmSchmidt, PlmON
    use ftypes

    implicit none

    real(dp), intent(in) :: d(:), lat(:), lon(:)
    real(dp), intent(out) :: cilm(:,:,:)
    integer(int32), intent(in) :: nmax, lmax
    integer(int32), intent(in), optional :: norm, csphase
    real(dp), intent(out), optional :: chi2
    real(dp), intent(in), optional :: weights(:), g(:,:)
    integer(int32), intent(out), optional :: exitstatus
    integer(int32), parameter :: opt = 80
    integer(int32) :: ncoef, i, l, m, ind1, ind2, info, lwork, opt1, phase, &
                      astat(7), lnorm
    real(dp) :: pi, lonr
    real(dp), allocatable :: mm(:), gg(:, :), p(:), work(:), y(:)
#ifdef LAPACK_UNDERSCORE
#define dgels dgels_
#endif
    external :: dgels

    if (present(exitstatus)) exitstatus = 0

    ncoef = (lmax+1)**2

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
            phase = 1

    end if

    if (present(weights)) then
        if (size(weights) < nmax) then
            print*, "Error --- SHExpandLSQ"
            print*, "WEIGHTS must be dimensioned as (NMAX) where NMAX is ", &
                nmax
            print*, "Input array is dimensioned ", size(weights)
            if (present(exitstatus)) then
                exitstatus = 1
                return
            else
                stop
            end if

        else if (nmax < ncoef) then
            print*, "Error --- SHExpandLSQ"
            print*, "When performing a weighted least-squares inversion, "
            print*, "NMAX must be greater or equal to (LMAX+1)**2, where"
            print*, "NMAX = ", nmax
            print*, "LMAX = ", lmax
            print*, "Input array is dimensioned ", size(weights)
            if (present(exitstatus)) then
                exitstatus = 1
                return
            else
                stop
            end if

        end if

    end if

    if (present(g)) then
        if (size(g(:,1)) < nmax .or. size(g(1,:)) < ncoef) then
            print*, "Error --- SHExpandLSQ"
            print*, "G must be dimensioned as (NMAX, (LMAX+1)**2)."
            print*, "NMAX = ", nmax
            print*, "LMAX = ", lmax
            print*, "Input array is dimensioned ", size(g(:,1)), size(g(1,:))
            if (present(exitstatus)) then
                exitstatus = 1
                return
            else
                stop
            end if
        end if
    end if

    lwork = min(ncoef, nmax)*(1+opt)

    astat = 0
    allocate (mm(max(ncoef, nmax)), stat = astat(1))
    allocate (gg(nmax, ncoef), stat = astat(2))
    allocate (p((lmax+1)*(lmax+2)/2), stat = astat(3))
    allocate (work(lwork), stat = astat(4))

    if (astat(1) /= 0 .or. astat(2) /= 0 .or. astat(3) /= 0 &
            .or. astat(4) /= 0) then
        print*, "Error --- SHExpandLSQ"
        print*, "Problem allocating arrays MM, GG, P, WORK: ", astat(1), &
            astat(2), astat(3), astat(4)
        if (present(exitstatus)) then
            exitstatus = 3
            return
        else
            stop
        end if

    end if

    pi = acos(-1.0_dp)
    mm = 0.0_dp

    !--------------------------------------------------------------------------
    !
    !   Compute the matrix G (nmax by ncoef)
    !
    !--------------------------------------------------------------------------
    if (present(g)) then
      gg = g

    else
        do i = 1, nmax
            if (present(exitstatus)) then
                select case(lnorm)
                    case(1)
                        call PlmBar(p, lmax, sin(lat(i) * pi / 180.0_dp), &
                                    csphase = phase, exitstatus = exitstatus)
                    case(2)
                        call PlmSchmidt(p, lmax, sin(lat(i) * pi / 180.0_dp), &
                                        csphase = phase, &
                                        exitstatus = exitstatus)
                    case(3)
                        call PLegendreA(p, lmax, sin(lat(i) * pi / 180.0_dp), &
                                        csphase = phase, &
                                        exitstatus = exitstatus)
                    case(4)
                        call PlmON(p, lmax, sin(lat(i) * pi / 180.0_dp), &
                                   csphase = phase, exitstatus = exitstatus)
                end select
                if (exitstatus /= 0) return

            else
                select case(lnorm)
                    case(1)
                        call PlmBar(p, lmax, sin(lat(i) * pi / 180.0_dp), &
                                    csphase = phase)
                    case(2)
                        call PlmSchmidt(p, lmax, sin(lat(i) * pi / 180.0_dp), &
                                        csphase = phase)
                    case(3)
                        call PLegendreA(p, lmax, sin(lat(i) * pi / 180.0_dp), &
                                        csphase = phase)
                    case(4)
                        call PlmON(p, lmax, sin(lat(i) * pi / 180.0_dp), &
                                   csphase = phase)
                end select

            end if

            lonr = lon(i) * pi / 180.0_dp
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

    end if
    !--------------------------------------------------------------------------
    !
    !   Do least squares inversion with no weights
    !
    !       m = [G' G]^-1 G' d
    !
    !   or
    !
    !       m = [G' W G]^-1 G' W d
    !
    !   with weights using the LAPACK routine DGELS.
    !
    !--------------------------------------------------------------------------
    if (present(weights)) then
        do i = 1, nmax
            mm(i) = d(i) * sqrt(weights(i))
            gg(i, :) = gg(i, :) * sqrt(weights(i))
        enddo

    else
        mm(1:nmax) = d(1:nmax)

    end if

    call dgels('N', nmax, ncoef, 1, gg, nmax, mm, max((lmax+1)**2, nmax), &
               work, lwork, info)

    if (info /= 0) then
        print*, "Error --- SHExpandLSQ"
        print*, "Problem performing the least squares inversion."
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
    if (present(chi2) .and. nmax >= ncoef) then
        chi2 = 0.0_dp

        do i = ncoef + 1, nmax, 1
            chi2 = chi2 + mm(i)**2
        end do

    end if

    ! deallocate memory
    select case(lnorm)
        case(1)
            call PlmBar(p, -1, 0.0_dp)
        case(2)
            call PlmSchmidt(p, -1, 0.0_dp)
        case(4)
            call PlmON(p, -1, 0.0_dp)

    end select

    deallocate (mm)
    deallocate (gg)
    deallocate (p)
    deallocate (work)

end subroutine SHExpandLSQ
