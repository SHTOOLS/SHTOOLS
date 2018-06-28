module PlanetsConstants
!-------------------------------------------------------------------------------
!
!   Constants related to the planets that are useful for gravity and topography
!   analyses.
!
!   Copyright (c) 2016, SHTOOLS
!   All rights reserved.
!
!-------------------------------------------------------------------------------
    implicit none

    real*8, parameter :: &
        Grav_constant = 6.67384d-11, & ! CODATA, Mohm and Taylor (2012)
        pi_constant =  3.14159265358979
    real*8, parameter :: &
        mu0_constant = pi_constant * 4.0d-7 ! magnetic constant

    ! The Moon
    real*8, parameter :: &
        R_Moon = 1737151.0d0, & ! Wieczorek (2015)
        GM_Moon = 4902.80007d9, & ! DE430, Williams et al. (2014)
        Mass_Moon = GM_Moon / Grav_constant, & ! Wieczorek et al. (2005)
        a_orbit_Moon = 384399.014d3 , & ! DE430, Williams et al. (2013)
        Omega_Moon = 2.0d0 * pi_constant / 27.3215820d0 / 24.0d0 / 60.0d0 / 60.0d0, & ! Yoder (1995)
        I_Solid_Moon = 0.3931120d0, & ! Williams et al. (2016)
        gamma_Moon = 227.7317d-6, & ! (B-A)/C, Williams et al. (2014)
        beta_Moon = 631.0213d-6 ! (C-A)/B, Williams et al. (2014)
    real*8, parameter ::  &
        density_Moon = Mass_Moon * 3.0d0 / 4.0d0 / pi_constant / R_Moon**3, &
        g0_Moon = GM_Moon / R_Moon**2 ! gravitational acceleration at R_Moon, not including rotation.

    ! Mars
    real*8, parameter :: &
        R_Mars = 3389500.12207057, & ! MarsTopo2600, Wieczorek (2015)
        GM_Mars = 0.4282837581575610d+14, & ! jgmro_120d, Konopliv et al. (2016)
        Mass_Mars = GM_Mars / Grav_constant, & ! GM / Grav_constant
        Omega_Mars = 350.891985307d0 * pi_constant / 180.0d0 / 24.0d0 / 60.0d0 / 60.0d0, & ! Konopliv et al. (2016)
        f_Mars = 0.0052276178, & ! Ardalan et al. (2010)
        a_Mars = 3395.42800d3, & ! Ardalan et al. (2010)
        b_Mars = 3377.6780d3, & ! Ardalan et al. (2010)
        U0_mars = 12654875.0d0 ! Ardalan et al. (2010)
    real*8, parameter :: &
        g0_Mars = GM_Mars / R_Mars**2, & ! gravitational acceleration at R_Mars, not including rotation.&
        density_Mars = Mass_Mars * 3.0d0 / 4.0d0 / pi_constant / R_Mars**3

    ! Venus
    real*8, parameter :: &
        R_Venus = 6051.877d3, & ! VenusTopo719, Wieczorek (2015)
        GM_Venus = 324858.592079d9, & ! MGNP180U, Konopliv et al. (1999)
        Mass_Venus = GM_Venus / Grav_constant, & ! GM / Grav_constant
        Omega_Venus = -2.0d0 * pi_constant / 243.020d0 / 24.0d0 / 60.0d0 / 60.0d0 ! Konopliv et al. (1999)
    real*8, parameter :: &
        density_Venus = Mass_Venus * 3.0d0 / 4.0d0 / pi_constant / R_Venus**3, &
        g0_Venus = GM_Venus / R_Venus**2

    ! Earth
    real*8, parameter :: &
        GM_egm2008 = 3986004.415d8, & ! EGM2008, Pavlis et al. (2012)
        Mass_egm2008 = GM_egm2008 / Grav_constant, &
        a_WGS84 = 6378137.0d0, & ! WGS84 ellipsoid semi-major axis
        f_WGS84 = 1.0d0 / 298.2572235630d0, & ! WGS84 ellipsoid flattening
        gm_WGS84 = 3986004.418d8, & ! WGS84 gm, including the atmosphere
        Mass_WGS84 = GM_WGS84 / Grav_constant, &
        omega_WGS84 = 7292115.0d-11, & ! WGS84 angular velocity
        gma_WGS84 = 3.50d8, & ! WGS84 GM of the atmosphere
        b_WGS84 = 6356752.31420d0, & ! WGS84 semi-minor axis
        U0_WGS84 = 62636851.7146d0, & ! WGS84 Theoretical normal potential
        r3_WGS84 = 6371000.7900d0 ! WGS84 Radius of sphere of equal volume

    ! Mercury
    real*8, parameter :: &
        GM_Mercury = 2.2031815411154894d+13, & ! ggmes_100v07, Mazarico et al. (2014)
        Mass_Mercury = GM_Mercury / Grav_constant, &
        R_Mercury = 2439.40197456433d3, & ! gtmes_150v05, Smith et al. (2012)
        Omega_orbit_Mercury = 2.0d0 * pi_constant / 87.969216879d0 / 24.0d0 / 60.0d0 / 60.0d0, & ! ggmes_100v07, Mazarico et al. (2014)
        Omega_Mercury = 6.1385108d0 * 2.0d0 * pi_constant / 360.0d0 / (24.0d0 * 60.0d0 * 60.0d0) ! ggmes_100v07, Mazarico et al. (2014)
    real*8, parameter :: &
        density_Mercury = Mass_Mercury * 3.0d0 / 4.0d0 / pi_constant / R_Mercury**3, &
        g0_Mercury = GM_Mercury / R_Mercury**2

end module PlanetsConstants
