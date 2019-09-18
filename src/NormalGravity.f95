function NormalGravity(geocentric_lat, GM, omega, a, b)
!------------------------------------------------------------------------------
!
!   Compute the total predicted gravity normal to the ellipsoid at a given
!   GEOCENTRIC latitude (input in DEGREES), using Somigliana's formula.
!   This is taken from Physical Geodesy (Hofmann-Wellenhof and Moritz,
!   sec. ed.), sections 2.7 and 2.8.
!
!   Copyright (c) 2016, SHTOOLS
!   All rights reserved.
!
!------------------------------------------------------------------------------
    use ftypes

    implicit none

    real(dp) :: NormalGravity
    real(dp), intent(in) :: geocentric_lat, gm, omega, a, b
    real(dp) :: geodetic_lat, pi, ga, gb, m, ep, bigE, q0, q0p

    if (a < b) then
        print*, "Warning --- NormalGravity"
        print*, "The semimajor axis A should be greater than the " // &
                "semiminor axis B."
    end if

    if (a == b .and. omega == 0.0_dp) then
        NormalGravity = GM / a**2

    else if (a == b .and. omega /= 0.0_dp) then
        print*, "Warning --- NormalGravity"
        print*, "A can not be equal to B when OMEGA is non zero."
        print*, "Setting OMEGA equal to zero."
        NormalGravity = GM / a**2

    else
        pi = acos(-1.0_dp)

        m = omega**2 * a**2 * b / GM

        bigE = sqrt(a**2 - b**2) ! linear eccentricity

        ep = bigE / b ! second eccentricity

        q0 = 0.50_dp * ((1.0_dp + 3.0_dp * (b / bigE)**2) * atan(bigE / b) &
                        - 3.0_dp * b / bigE)

        q0p = 3.0_dp * (1.0_dp + (b / bigE)**2) &
              * (1.0_dp - b / bigE * atan(bigE / b) ) - 1.0_dp

        ! gravity on equator
        ga = GM / (a*b) * (1.0_dp - m - m * ep * q0p / 6.0_dp / q0)

        gb = GM / a**2 * (1.0_dp + m * ep * q0p / 3.0_dp / q0) ! gravity at poles

        geodetic_lat = atan((a / b)**2 * tan(geocentric_lat * pi / 180.0_dp))

        NormalGravity = a * ga * cos(geodetic_lat)**2 + b * gb &
                        * sin(geodetic_lat)**2

        NormalGravity = NormalGravity / sqrt(a**2 * cos(geodetic_lat)**2 &
                        + b**2 * sin(geodetic_lat)**2 )
    end if

end function NormalGravity
