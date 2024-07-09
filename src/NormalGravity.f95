function NormalGravity(geocentric_lat, GM, omega, a, b)
!------------------------------------------------------------------------------
!
!   Compute the normal gravity at a given GEOCENTRIC latitude (input in
!   degrees). The normal gravity is the norm of the total gravity (gravitation
!   and centrifugal) on the surface of a rotating ellipsoid.
!
!   For a rotating ellipsoid, the surface corresponds to a constant potential,
!   and the gravity vector is normal to the surface. The normal gravity is
!   computed using Somigliana's formula, as taken from Physical Geodesy
!   (Hofmann-Wellenhof and Moritz, 2nd ed., sections 2.7 and 2.8). In this
!   routine, the geodetic latitude is computed from the input geocentric
!   latitude.
!
!   For the case of a sphere (a is equal to b), the normal gravity is defined
!   as the magnitude of the sum of the normal gravitation (GM/r**2) and
!   centrifugal gravity. For a rotating sphere, the surface does not correspond
!   to a constant potential and the gravity vector is not normal to the
!   surface.
!
!   Copyright (c) 2005-2024, SHTOOLS
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

    pi = acos(-1.0_dp)

    if (a == b .and. omega == 0.0_dp) then
        NormalGravity = GM / a**2

    else if (a == b .and. omega /= 0.0_dp) then
        NormalGravity = GM**2 / a**4 &
          + a * omega**2 * cos(geocentric_lat * pi / 180.0_dp)**2 &
          * (a * omega**2 - 2.0_dp * GM / a**2)
        NormalGravity = sqrt(NormalGravity)

    else
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
