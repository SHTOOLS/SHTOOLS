real*8 function MakeGridPoint(cilm, lmax, lat, longitude, norm, &
                              csphase, dealloc)
!------------------------------------------------------------------------------
!
!   This function will determine the value at a given latitude and
!   longitude corresponding to the given set of spherical harmonics.
!   Latitude and Longitude are assumed to be in DEGREES!
!
!   Calling Parameters
!
!       IN
!           cilm        Spherical harmonic coefficients, with dimensions
!                       (2, lmax+1, lmax+1).
!           lmax        Maximum degree used in the expansion.
!           lat         latitude (degrees).
!           long        longitude (degrees).
!
!       OPTIONAL (IN)
!           norm        Spherical harmonic normalization:
!                           (1) "geodesy" (default)
!                           (2) Schmidt
!                           (3) unnormalized
!                           (4) orthonormalized
!           csphase     1: Do not include the phase factor of (-1)^m
!                       -1: Apply the phase factor of (-1)^m.
!           dealloc     If (1) Deallocate saved memory in Legendre function
!                       routines. Default (0) is not to deallocate memory.
!
!   Dependencies:       PlmBar, PlBar, PlmSchmidt, PlmON, CSPHASE_DEFAULT
!
!   Notes:
!       1.  If lmax is greater than the the maximum spherical harmonic
!           degree of the input file, then this file will be ZERO PADDED!
!           (i.e., those degrees after lmax are assumed to be zero).
!
!   Copyright (c) 2016, SHTOOLS
!   All rights reserved.
!
!------------------------------------------------------------------------------
    use SHTOOLS, only: PlmBar, PLegendreA, PlmSchmidt, PlmON, CSPHASE_DEFAULT

    implicit none

    real*8, intent(in):: cilm(:,:,:), lat, longitude
    integer, intent(in) :: lmax
    integer, intent(in), optional :: norm, csphase, dealloc
    real*8 :: pi, x, expand, lon
    integer :: index, l, m, l1, m1, lmax_comp, phase, astat(3)
    real*8, allocatable :: pl(:), cosm(:), sinm(:)

    if (size(cilm(:,1,1)) < 2 .or. size(cilm(1,:,1)) < lmax+1 .or. &
            size(cilm(1,1,:)) < lmax+1) then
        print*, "Error --- MakeGridPoint"
        print*, "CILM must be dimensioned as (2, LMAX+1, LMAX+1) " // &
                "where LMAX is ", lmax
        print*, "Input dimension is ", size(cilm(:,1,1)), size(cilm(1,:,1)), &
                size(cilm(1,1,:))
        stop
    end if

    if (present(norm)) then
        if (norm >4 .or. norm < 1) then
            print*, "Error - MakeGridPoint"
            print*, "Parameter NORM must be 1, 2, 3, or 4"
            stop
        end if
    end if

    if (present(csphase)) then
        if (csphase == -1) then
             phase = -1.0d0
        else if (csphase == 1) then
                phase = 1.0d0
        else
            print*, "Error --- MakeGridPoint"
            print*, "CSPHASE must be 1 (exclude) or -1 (include)."
            print*, "Input value is ", csphase
            stop
        end if

    else
        phase = dble(CSPHASE_DEFAULT)

    end if

    allocate (pl(((lmax+1) * (lmax+2)) / 2), stat = astat(1))
    allocate (cosm(lmax+1), stat = astat(2))
    allocate (sinm(lmax+1), stat = astat(3))

    if (sum(astat(1:3)) /= 0) then
        print*, "Error --- MakeGridPoint"
        print*, "Cannot allocate memory for arrays PL, MCOS and MSIN", &
                astat(1), astat(2), astat(3)
        stop
    end if

    pi = acos(-1.0d0)
    x = sin(lat * pi / 180.0d0)
    lon = longitude * pi / 180.0d0

    lmax_comp = min(lmax, size(cilm(1,1,:)) - 1)

    if (present(norm)) then
        if (norm == 1) call PlmBar(pl, lmax_comp, x, csphase = phase)
        if (norm == 2) call PlmSchmidt(pl, lmax_comp, x, csphase = phase)
        if (norm == 3) call PLegendreA(pl, lmax_comp, x, csphase = phase)
        if (norm == 4) call PlmON(pl, lmax_comp, x, csphase = phase)

    else
        call PlmBar(pl, lmax_comp, x, csphase = phase)

    end if

    expand = 0.0d0

    ! Precompute sines and cosines. Use multiple angle identity to minimize
    ! number of calls to SIN and COS.
    sinm(1) = 0.0d0
    cosm(1) = 1.0d0

    if (lmax_comp > 0) then
        sinm(2) = sin(lon)
        cosm(2) = cos(lon)
    end if

    do m = 2, lmax_comp, 1
        m1 = m + 1
        sinm(m1) = 2 * sinm(m) * cosm(2) - sinm(m-1)
        cosm(m1) = 2 * cosm(m) * cosm(2) - cosm(m-1)
    end do

    do l = lmax_comp, 0, -1
        l1 = l + 1
        index = (l+1) * l / 2 + 1
        expand = expand + cilm(1,l1,1) * pl(index)

        do m = 1, l, 1
            m1 = m + 1
            index = index + 1
            expand = expand + (cilm(1,l1,m1) * cosm(m1) + &
                               cilm(2,l1,m1) * sinm(m1)) * pl(index)
        end do

    end do

    MakeGridPoint = expand

    ! deallocate memory
    if (present(dealloc)) then
        if (dealloc == 1) then
            if (present(norm)) then
                if (norm == 1) call PlmBar(pl, -1, x, csphase = phase)
                if (norm == 2) call PlmSchmidt(pl, -1, x, csphase = phase)
                if (norm == 4) call PlmON(pl, -1, x, csphase = phase)
            else
                call PlmBar(pl, -1, x, csphase = phase)
            end if
        end if
    endif

    deallocate (pl)
    deallocate (cosm)
    deallocate (sinm)

end function MakeGridPoint
