program TestSHExpandLSQ
!------------------------------------------------------------------------------
!
!    This program will (1) create a set of unevenly spaced data points from a
!    spherical harmonic file, and then (2) expand these into spherical
!    harmonics using a least-squares inversion if there are more data points
!    than spherical harmonic coefficients (overdetermined case), or using a
!    minimum norm solution if there are more coefficients than data points
!    (underdetermined case).
!
!    Copyright (c) 2005, SHTOOLS
!    All rights reserved.
!
!------------------------------------------------------------------------------
    use SHTOOLS
    use ftypes

    implicit none

    integer, parameter :: dmax = 10000, degmax = 120
    real(dp) :: d(dmax), lat(dmax), lon(dmax), cilm(2,degmax+1, degmax+1), &
                x, y, z, pi, maxerror, cilm1(2,degmax+1, degmax+1), &
                dd(dmax), misfit
    integer(int32) :: nmax, lmax, l, i, seed, lmaxfile
    character(80) :: infile

    d = 0.0_dp
    dd = 0.0_dp
    lon = 0.0_dp
    lat = 0.0_dp
    pi = acos(-1.0_dp)

    infile = "../../ExampleDataFiles/MarsTopo719.shape"

    call SHRead(infile, cilm, lmaxfile)
    print*, "Maximum degree of spherical harmonic file = ", lmaxfile

    print*, "Number of random data points to use > "
    read(*,*) nmax
    lmax = floor(sqrt(dble(nmax)) - 1)
    print*, "Maximum spherical harmonic degree for overdetermined least-squares inversion = ", lmax

    !--------------------------------------------------------------------------
    !
    !    Create synthetic data from the known spherical harmonic coefficients.
    !    Data points are located at random points on the sphere.
    !
    !--------------------------------------------------------------------------

    seed = -13453

    do i = 1, nmax
        x = 2.0_dp * RandomN(seed) - 1.0_dp
        y = 2.0_dp * RandomN(seed) - 1.0_dp
        z = 2.0_dp * RandomN(seed) - 1.0_dp
        lat(i) = atan2(z, sqrt(x**2 + y**2)) * 180.0_dp  /pi
        lon(i) = atan2(y, x) * 180.0_dp / pi
        d(i) = MakeGridPoint(cilm, lmaxfile, lat(i), lon(i), norm=1)
    end do

    print*, "maximum (km) = ", maxval(d(1:nmax)) / 1.e3_dp
    print*, "minimum (km) = ", minval(d(1:nmax)) / 1.e3_dp

    !--------------------------------------------------------------------------
    !
    !    Do least squares inversion for increasing values of l.
    !
    !--------------------------------------------------------------------------

    do l = 1, lmax, 1

        print*, "l = ", l

        call SHExpandLSQ(cilm1, d(1:nmax), lat(1:nmax), lon(1:nmax), nmax, &
                         l, norm=1, chi2=misfit)

        do i = 1, nmax
            dd(i) = MakeGridPoint(cilm1, l, lat(i), lon(i), norm=1)
        end do

        maxerror = maxval(abs(d(1:nmax) - dd(1:nmax)))

        print*, "Maximum error between input and output data points (km) = ", maxerror / 1.e3_dp
        print*, "Sum of squares residuals (km^2) = ", misfit / 1.e6_dp

    end do

end program TestSHExpandLSQ