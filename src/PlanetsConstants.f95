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
        R_Moon = 1737151.198265080d0, & ! degree-0 term of the LOLA model LOLA2600p.sh
        GM_Moon = 4902.801076d9, & ! Konopliv et al. (2001), LP150Q 
        Mass_Moon = GM_Moon / Grav_constant, & ! Wieczorek et al. (2005)
        a_Moon = 384399.0d3, & ! Semi-major axis (meters), Williams et al. (2001)
        R0_pot_Moon = 1738.d3, & ! Reference radius of LP150Q gravity model.
        Omega_Moon = 2.0d0*pi_constant/27.3215820d0/24.0d0/60.0d0/60.0d0, & ! Lunar rotation rate, Yoder 1995
        c_moi_Moon = 0.39320d0, & ! normalized polar moment of inertia (R=1738), Konopliv et al. (1998)
        i_moi_Moon = 0.39310d0, & ! average moment of inertia (R=1738), Konopliv et al. (1998)
        gamma_moi_Moon = 227.7317d-6, & ! (B-A)/C from Williams et al. 2014
        beta_moi_Moon = 631.0213d-6 ! (C-A)/B) from Williams et al. 2014
    real*8, parameter ::  &
        rho_bar_Moon =  Mass_Moon * 3.0d0 / 4.0d0 / pi_constant / R_Moon**3, &
        g0_Moon = GM_Moon / R_Moon**2 ! gravitational acceleration at R_Moon, not including rotation.
    
    ! Mars
    real*8, parameter :: &
        R_Mars = 3389500.12207057, & ! MarsTopo2600
        GM_Mars = 42828.3756410d9, & ! JGMRO_110c, Konopliv et al. 2011
        Mass_Mars = GM_Mars / Grav_constant, & ! GM/Grav_constant   
        R0_pot_Mars = 3396.d3, & ! reference radius of JGMRO_110c gravity model.
        Omega_Mars = 350.8919830d0*pi_constant/180.0d0/24.0d0/60.0d0/60.0d0, & ! Rotation rate (radians/sec), from Yuan et al. 2001
        f_Mars = 0.0052276178, & ! Ardalan et al. 2010
        a_Mars = 3395.42800d3, & ! mean equatorial radius; Ardalan et al. 2010
        b_Mars = 3377.6780d3, & ! Ardalan et al. 2010
        W0_mars = 12654875.0d0 ! reference potential, Ardalan et al. 2010
    real*8, parameter :: &
        g0_Mars = GM_Mars/R_Mars**2, & ! gravitational acceleration at R_Mars, not including rotation.&
        rho_bar_Mars = Mass_Mars * 3.0d0 / 4.0d0 / pi_constant / R_Mars**3

    ! Venus 
    real*8, parameter :: &
        R_Venus = 6051.877d3, & ! VenusTopo719
        R0_pot_Venus = 6051.0d3, & ! Reference radius of MGNP180U gravity model
        GM_Venus = 324858.592079d9, & ! Konopliv et al. (1999), MGNP180U
        Mass_Venus = GM_Venus / Grav_constant, & ! GM/Grav_constant
        Omega_Venus = -2.0d0*pi_constant/243.020d0/24.0d0/60.0d0/60.0d0 ! Rotation rate, Konopliv et al. 1999
    real*8, parameter :: &
        rho_bar_Venus = Mass_Venus * 3.0d0 / 4.0d0 / pi_constant / R_Venus**3, &
        g0_Venus = GM_Venus/R_Venus**2
                
    ! Earth
    real*8, parameter :: &
        GM_Earth = 3986004.415d8, & ! EGM2008        
        R0_pot_Earth = 6378136.30d0, & ! EGM2008
        Mass_Earth = GM_Earth / Grav_constant, &
        WGS84_a = 6378137.0d0, & ! WGS84 ellipsoid semi-major axis
        WGS84_f = 1.0d0 /298.2572235630d0, & ! WGS84 ellipsoid flattening
        WGS84_gm = 3986004.418d8, & ! WGS84 gm, includes the atmosphere
        WGS84_omega = 7292115.0d-11, & ! WGS84 angular velocity, radians/sec
        WGS84_gma = 3.50d8, & ! WGS84 gm of the atmosphere
        WGS84_b = 6356752.31420d0, & ! WGS84 semi-minor axis
        WGS84_U0 = 62636851.7146d0, & ! WGS84 Theoretical normal potential
        WGS84_r3 = 6371000.7900d0 ! WGS84 Radius of sphere of equal volume
                
    ! Mercury
    real*8, parameter ::  &
        GM_Mercury = 22031.839224134801e+09, & ! model ggmes_20v04_sha available on the PDS
        Mass_Mercury = GM_Mercury / Grav_constant, &
        R0_pot_Mercury = 2440.0000000000000e+03, & ! for model ggmes_20v04_sha
        R_Mercury = 2.4393717480908799e+06, & ! degree 0 term of the model gtmes_24v01_sha available on the PDS
        Omega_Mercury_orbit = 2.0d0 *pi_constant/87.9690d0/24.0d0/60.0d0/60.0d0, & ! Rambaux and Bois (2004)
        Omega_Mercury_spin = 3.0d0 / 2.0d0 * Omega_Mercury_orbit 
    real*8, parameter :: &
        rho_bar_Mercury = Mass_Mercury * 3.0d0 / 4.0d0 / pi_constant / R_Mercury**3, &
        g0_Mercury = GM_Mercury/R_Mercury**2

end module PlanetsConstants
