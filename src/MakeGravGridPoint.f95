function MakeGravGridPoint(cilm, lmax, gm, r0, r, lat, lon, omega, dealloc)
!------------------------------------------------------------------------------
!
!   This function will determine the 3 components (r-hat, theta-hat, phi-hat)
!   of the gravity vector (gravitational force + centrifugal force) at a given
!   latitude, longitude and radius. Latitude and longitude must be input in
!   degrees. The input coefficients must be 4-pi normalized, excluding the
!   Condon-Shortley phase factor.
!
!   Calling Parameters
!
!       IN
!           cilm        Spherical harmonic coefficients, with dimensions
!                       (2, lmax+1, lmax+1).
!           lmax        Maximum degree used in the expansion.
!           gm          Product of the gravitatonal constant and the planet's
!                       mass.
!           r0          Reference radius of the potential coefficients.
!           r           Radius where the gravity vector is computed (meters).
!           lat         Latitude where the gravity vector is computed, in
!                       degrees.
!           lon         Longitude where the gravity vector is computed, in
!                       degrees.
!
!       OPTIONAL (IN)
!           omega       The angular rotation rate of the planet. If present,
!                       compute the gravity vector in a body-fixed rotating
!                       frame.
!           dealloc     If (1) deallocate saved memory in the Legendre function
!                       routines. Default (0) is not to deallocate memory.
!
!   Copyright (c) 2005-2020, SHTOOLS
!   All rights reserved.
!
!------------------------------------------------------------------------------
    use SHTOOLS, only: PlmBar_d1
    use ftypes

    implicit none

    real(dp), dimension(3) :: MakeGravGridPoint
    real(dp), intent(in) :: cilm(:,:,:), gm, r0, r, lat, lon
    real(dp), intent(in), optional :: omega
    integer(int32), intent(in) :: lmax
    integer(int32), intent(in), optional :: dealloc
    real(dp) :: pi, x, expand(3), lon_rad, prefactor(lmax)
    integer(int32) :: index, l, m, l1, m1, lmax_comp, astat(4)
    real(dp), allocatable :: pl(:), dpl(:), cosm(:), sinm(:)

    if (size(cilm(:,1,1)) < 2 .or. size(cilm(1,:,1)) < lmax+1 .or. &
            size(cilm(1,1,:)) < lmax+1) then
        print*, "Error --- MakeGravGridPoint"
        print*, "CILM must be dimensioned as (2, LMAX+1, LMAX+1) " // &
                "where LMAX is ", lmax
        print*, "Input dimension is ", size(cilm(:,1,1)), size(cilm(1,:,1)), &
                size(cilm(1,1,:))
        stop
    end if

    lmax_comp = min(lmax, size(cilm(1,1,:)) - 1)

    allocate (pl(((lmax_comp+1) * (lmax_comp+2)) / 2), stat = astat(1))
    allocate (dpl(((lmax_comp+1) * (lmax_comp+2)) / 2), stat = astat(2))
    allocate (cosm(lmax_comp+1), stat = astat(3))
    allocate (sinm(lmax_comp+1), stat = astat(4))

    if (sum(astat(1:4)) /= 0) then
        print*, "Error --- MakeGravGridPoint"
        print*, "Cannot allocate memory for arrays PL, DPL, COSM and SINM", &
                astat(1), astat(2), astat(3), astat(4)
        stop
    end if

    pi = acos(-1.0_dp)
    x = sin(lat * pi / 180.0_dp)
    lon_rad = lon * pi / 180.0_dp

    call PlmBar_d1(pl, dpl, lmax_comp, x, csphase = 1)
    dpl = -dpl * cos(lat * pi / 180.0_dp)

    ! Precompute sines and cosines. Use multiple angle identity to minimize
    ! number of calls to SIN and COS.
    sinm(1) = 0.0_dp
    cosm(1) = 1.0_dp

    if (lmax_comp > 0) then
        sinm(2) = sin(lon_rad)
        cosm(2) = cos(lon_rad)
    end if

    do m = 2, lmax_comp, 1
        m1 = m + 1
        sinm(m1) = 2 * sinm(m) * cosm(2) - sinm(m-1)
        cosm(m1) = 2 * cosm(m) * cosm(2) - cosm(m-1)
    end do

    prefactor(1) = r0 / r  ! l=1
    do l = 2, lmax_comp, 1
        prefactor(l) = prefactor(l-1) * (r0 / r)
    end do

    expand(1) = -cilm(1,1,1)
    expand(2:3) = 0.0_dp

    do l = 1, lmax_comp, 1
        l1 = l + 1
        index = (l+1) * l / 2 + 1
        expand(1) = expand(1) - prefactor(l) * l1 * cilm(1,l1,1) * pl(index)
        expand(2) = expand(2) + prefactor(l) * cilm(1,l1,1) * dpl(index)

        do m = 1, l, 1
            m1 = m + 1
            index = index + 1
            expand(1) = expand(1) - prefactor(l) * l1 * pl(index) * &
                        (cilm(1,l1,m1) * cosm(m1) + cilm(2,l1,m1) * sinm(m1))
            expand(2) = expand(2) + prefactor(l) * dpl(index) * &
                        (cilm(1,l1,m1) * cosm(m1) + cilm(2,l1,m1) * sinm(m1))
            expand(3) = expand(3) + prefactor(l) * pl(index) * &
                        (-m * cilm(1,l1,m1) * sinm(m1) &
                         + m * cilm(2,l1,m1) * cosm(m1))
        end do

    end do

    expand(1:3) = expand(1:3) * gm / r**2
    if (abs(lat) /= 90.0_dp) then
        expand(3) = expand(3) / cos(lat * pi / 180.0_dp)
    else
        expand(3) = 0.0_dp
    end if

    ! Add rotational effects
    if (present(omega)) then
        expand(1) = expand(1) + r * ( cos(lat * pi / 180.0_dp) * omega )**2
        expand(2) = expand(2) + r * cos(lat * pi / 180.0_dp) * &
                                sin(lat * pi / 180.0_dp) * omega**2

    end if

    MakeGravGridPoint(1:3) = expand(1:3)

    ! deallocate memory
    if (present(dealloc)) then
        if (dealloc == 1) then
            call PlmBar_d1(pl, dpl, -1, x)
        end if
    end if

    deallocate (pl)
    deallocate (dpl)
    deallocate (cosm)
    deallocate (sinm)

end function MakeGravGridPoint
