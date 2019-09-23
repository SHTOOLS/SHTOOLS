module PlanetsConstants
!-------------------------------------------------------------------------------
!
!   Constants related to the planets that are useful for gravity and topography
!   analyses.
!
!   Copyright (c) 2005-2019, SHTOOLS
!   All rights reserved.
!
!-------------------------------------------------------------------------------
    use ftypes

    implicit none

    real(dp), parameter :: &
        Grav_constant = 6.67408e-11_dp, & ! CODATA 2014, Mohr et al. (2016)
        pi_constant = 3.14159265358979_dp
    real(dp), parameter :: &
        mu0_constant = pi_constant * 4.0e-7_dp ! magnetic constant

    ! The Moon
    real(dp), parameter :: &
        R_Moon = 1737151.0_dp, & ! Wieczorek (2015)
        GM_Moon = 4902.80007e9_dp, & ! DE430, Williams et al. (2014)
        Mass_Moon = GM_Moon / Grav_constant, & ! Wieczorek et al. (2005)
        a_orbit_Moon = 384399.014e3_dp , & ! DE430, Williams et al. (2013)
        Omega_Moon = 2.0_dp * pi_constant / 27.3215820_dp / 24.0_dp / 60.0_dp / 60.0_dp, & ! Yoder (1995)
        I_Solid_Moon = 0.3931120_dp, & ! Williams et al. (2016)
        gamma_Moon = 227.7317e-6_dp, & ! (B-A)/C, Williams et al. (2014)
        beta_Moon = 631.0213e-6_dp ! (C-A)/B, Williams et al. (2014)
    real(dp), parameter :: &
        density_Moon = Mass_Moon * 3.0_dp / 4.0_dp / pi_constant / R_Moon**3, &
        g0_Moon = GM_Moon / R_Moon**2 ! gravitational acceleration at R_Moon, not including rotation.

    ! Mars
    real(dp), parameter :: &
        R_Mars = 3389500.12207057_dp, & ! MarsTopo2600, Wieczorek (2015)
        GM_Mars = 0.4282837581575610e+14_dp, & ! jgmro_120d, Konopliv et al. (2016)
        Mass_Mars = GM_Mars / Grav_constant, & ! GM / Grav_constant
        Omega_Mars = 350.891985307_dp * pi_constant / 180.0_dp / 24.0_dp / 60.0_dp / 60.0_dp, & ! Konopliv et al. (2016)
        f_Mars = 0.0052276178_dp, & ! Ardalan et al. (2010)
        a_Mars = 3395.42800e3_dp, & ! Ardalan et al. (2010)
        b_Mars = 3377.6780e3_dp, & ! Ardalan et al. (2010)
        U0_mars = 12654875.0_dp ! Ardalan et al. (2010)
    real(dp), parameter :: &
        g0_Mars = GM_Mars / R_Mars**2, & ! gravitational acceleration at R_Mars, not including rotation.&
        density_Mars = Mass_Mars * 3.0_dp / 4.0_dp / pi_constant / R_Mars**3

    ! Venus
    real(dp), parameter :: &
        R_Venus = 6051.877e3_dp, & ! VenusTopo719, Wieczorek (2015)
        GM_Venus = 324858.592079e9_dp, & ! MGNP180U, Konopliv et al. (1999)
        Mass_Venus = GM_Venus / Grav_constant, & ! GM / Grav_constant
        Omega_Venus = -2.0_dp * pi_constant / 243.020_dp / 24.0_dp / 60.0_dp / 60.0_dp ! Konopliv et al. (1999)
    real(dp), parameter :: &
        density_Venus = Mass_Venus * 3.0_dp / 4.0_dp / pi_constant / R_Venus**3, &
        g0_Venus = GM_Venus / R_Venus**2

    ! Earth
    real(dp), parameter :: &
        GM_egm2008 = 3986004.415e8_dp, & ! EGM2008, Pavlis et al. (2012)
        Mass_egm2008 = GM_egm2008 / Grav_constant, &
        a_WGS84 = 6378137.0_dp, & ! WGS84 ellipsoid semi-major axis
        f_WGS84 = 1.0_dp / 298.2572235630_dp, & ! WGS84 ellipsoid flattening
        gm_WGS84 = 3986004.418e8_dp, & ! WGS84 gm, including the atmosphere
        Mass_WGS84 = GM_WGS84 / Grav_constant, &
        omega_WGS84 = 7292115.0e-11_dp, & ! WGS84 angular velocity
        gma_WGS84 = 3.50e8_dp, & ! WGS84 GM of the atmosphere
        b_WGS84 = 6356752.31420_dp, & ! WGS84 semi-minor axis
        U0_WGS84 = 62636851.7146_dp, & ! WGS84 Theoretical normal potential
        r3_WGS84 = 6371000.7900_dp ! WGS84 Radius of sphere of equal volume

    ! Mercury
    real(dp), parameter :: &
        GM_Mercury = 2.2031815411154894e+13_dp, & ! ggmes_100v07, Mazarico et al. (2014)
        Mass_Mercury = GM_Mercury / Grav_constant, &
        R_Mercury = 2439.40197456433e3_dp, & ! gtmes_150v05, Smith et al. (2012)
        Omega_orbit_Mercury = 2.0_dp * pi_constant / 87.969216879_dp / 24.0_dp / 60.0_dp / 60.0_dp, & ! ggmes_100v07, Mazarico et al. (2014)
        Omega_Mercury = 6.1385108_dp * 2.0_dp * pi_constant / 360.0_dp / (24.0_dp * 60.0_dp * 60.0_dp) ! ggmes_100v07, Mazarico et al. (2014)
    real(dp), parameter :: &
        density_Mercury = Mass_Mercury * 3.0_dp / 4.0_dp / pi_constant / R_Mercury**3, &
        g0_Mercury = GM_Mercury / R_Mercury**2

end module PlanetsConstants
