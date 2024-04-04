subroutine LSQ_G(G, lat, lon, nmax, lmax, norm, csphase, exitstatus)
!------------------------------------------------------------------------------
!
!   This subroutine will compute the matrix G that is used when computing the
!   spherical harmonic coefficients of an irregularly sampled function by
!   least squares inversion, as used with SHExpandLSQ. The matrix G has
!   dimensions (nmax, (lmax+1)**2) where nmax is the number of data points and
!   lmax is the maximum spherical harmonic degree of the expansion. Each
!   element in a given row corresponds to the values of the spherical harmonic
!   functions for a given latitude and longitude. The elements in each row are
!   ordered by increasing degree, with all cosine terms for a given degree
!   followed by all sin terms for the same degree (with increasing order).
!
!   The default normalization convention for the matrix G is the "geodesy" 4pi
!   normalization, though this can be modified by supplying the optional
!   argument norm.
!
!   Calling Parameters
!
!       IN
!           lat     Vector of length nmax of the corresponding latitude points
!                   (in degrees).
!           lon     Vector of length nmax of the corresponding longitude points
!                   (in degrees).
!           nmax    Number of data points.
!           lmax    Maximum degree of spherical harmonic expansion.
!
!       OUT
!           G       The matrix G.
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
!   Copyright (c) 2024, SHTOOLS
!   All rights reserved.
!
!------------------------------------------------------------------------------
    use SHTOOLS, only: PlmBar, PLegendreA, PlmSchmidt, PlmON
    use ftypes

    implicit none

    real(dp), intent(in) :: lat(:), lon(:)
    real(dp), intent(out) :: G(:, :)
    integer(int32), intent(in) :: nmax, lmax
    integer(int32), intent(in), optional :: norm, csphase
    integer(int32), intent(out), optional :: exitstatus
    integer(int32) :: ncoef, i, l, m, ind1, ind2, phase, astat, lnorm
    real(dp) :: pi, lonr
    real(dp), allocatable :: p(:)

    if (present(exitstatus)) exitstatus = 0

    ncoef = (lmax+1)**2

    if (size(G(1, :)) < ncoef .or. size(G(:,1)) < nmax) then
        print*, "Error --- LSQ_G"
        print*, "G must be dimensioned as (NMAX, (LMAX+1)**2) with "
        print*, "NMAX = ", nmax
        print*, "LMAX = ", lmax
        print*, "Input dimension is ", size(G(:,1)), size(G(1,:))
        if (present(exitstatus)) then
            exitstatus = 1
            return
        else
            stop
        end if

    else if (size(lat) < nmax) then
        print*, "Error --- LSQ_G"
        print*, "LAT must be dimensioned as (NMAX) where NMAX is ", nmax
        print*, "Input array is dimensioned ", size(lat)
        if (present(exitstatus)) then
            exitstatus = 1
            return
        else
            stop
        end if

    else if (size(lon) < nmax) then
        print*, "Error --- LSQ_G"
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
            print*, "Error --- LSQ_G"
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
            print*, "Error ---- LSQ_G"
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

    allocate (p((lmax+1)*(lmax+2)/2), stat = astat)

    if (astat /= 0) then
        print*, "Error --- LSQ_G"
        print*, "Problem allocating array P: ", astat
        if (present(exitstatus)) then
            exitstatus = 3
            return
        else
            stop
        end if

    end if

    pi = acos(-1.0_dp)
    G = 0.0_dp

    !--------------------------------------------------------------------------
    !
    !   Calculate matrix G (nmax by ncoef)
    !
    !--------------------------------------------------------------------------
    do i = 1, nmax
        if (present(exitstatus)) then
            select case(lnorm)
                case(1)
                    call PlmBar(p, lmax, sin(lat(i) * pi / 180.0_dp), &
                                csphase = phase, exitstatus = exitstatus)
                case(2)
                    call PlmSchmidt(p, lmax, sin(lat(i) * pi / 180.0_dp), &
                                    csphase = phase, exitstatus = exitstatus)
                case(3)
                    call PLegendreA(p, lmax, sin(lat(i) * pi / 180.0_dp), &
                                    csphase = phase, exitstatus = exitstatus)
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
                G(i,ind1) = p(ind2) * cos(m*lonr)

            end do

            ! do sin terms
            do m = 1, l, 1
                ind1 = ind1 + 1
                ind2 = l*(l+1)/2 + m + 1
                G(i,ind1) = p(ind2) * sin(m*lonr)

            end do

        end do

    end do

    ! deallocate memory
    select case(lnorm)
        case(1)
            call PlmBar(p, -1, 0.0_dp)
        case(2)
            call PlmSchmidt(p, -1, 0.0_dp)
        case(4)
            call PlmON(p, -1, 0.0_dp)

    end select

end subroutine LSQ_G
