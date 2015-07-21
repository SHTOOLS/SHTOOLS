subroutine MakeGravGradGridDH(cilm, lmax, gm, r0, a, f, vxx, vyy, vzz, vxy, &
                                vxz, vyz, n, sampling, lmax_calc)
!-------------------------------------------------------------------------------
!
!   Given the gravitational spherical harmonic coefficients CILM, this 
!   subroutine will compute 2D Driscol and Healy sampled grids of the six 
!   components of the gravity "gradient" tensor in a local north-oriented 
!   reference frame:
!
!       (Vxx,   Vxy,    Vxz)
!       (Vyx,   Vyy,    Vyz)
!       (Vzx,   Vzy,    Vzz)
!   
!   where X points NORTH, Y points WEST, and Z points UPWARD. The gravitational 
!   potential is defined as
!
!       V = GM/r Sum_{l=0}^LMAX (R0/r)^l Sum_{m=-l}^l C_{lm} Y_{lm}, 
!
!   Laplaces equation implies that Vxx + Vyy + Vzz = 0, and the gravity tensor 
!   is symmetric. The components are calculated according to eq. 1 in 
!   Petrovskaya and Vershkov (2006, J. Geod, 80, 117-127), which is based on the 
!   eq. 3.28 in Reed (1973, Ohio State Univ., Dept. Geod. Sci., Rep. 201, 
!   Columbus, OH). Note that Reed's equations are in terms of latitude, and 
!   that the Y axis points East.
!
!       Vzz = Vrr
!       Vxx = 1/r Vr + 1/r^2 Vtt
!       Vyy = 1/r Vr + 1/r^2 /tan(t) Vt + 1/r^2 /sin(t)^2 Vpp
!       Vxy = 1/r^2 /sin(t) Vtp - cos(t)/sin(t)^2 /r^2 Vp
!       Vxz = 1/r^2 Vt - 1/r Vrt
!       Vyz = 1/r^2 /sin(t) Vp - 1/r /sin(t) Vrp
!
!   where r, t, p stand for radius, theta, and phi, and subscripts on V denote 
!   partial derivatives.
!
!   The output grid are in units of 1/s and are cacluated on a flattened 
!   ellipsoid with semi-major axis A and flattening F. To obtain units of 
!   Eotvos (10^-9 s^-1), multiply by 10^9. The output grids contain N samples 
!   in latitude and longitude by default, but if the optional parameter SAMPLING
!   is set to 2, the grids will contain N samples in latitude and 2N samples 
!   in longitude. The first latitudinal band of the grid corresponds to 90 N, 
!   the latitudinal band for 90 S is not calculated, and the latitudinal 
!   sampling interval is 180/N degrees. The first longitudinal band is 0 E, 
!   the longitudinal band for 360 E is not calculated, and the longitudinal 
!   sampling interval is 360/N for equally sampled and 180/N for equally spaced
!   grids, respectively.
!
!   This routine assumes that the spherical harmonic coefficients of geodesy 
!   4-pi normalized.
!
!   Calling Parameters
! 
!       IN
!           cilm        Gravitational spherical harmonic coefficients.
!           lmax        The maximum spherical harmonic degree of the function, 
!                       used to determine the number of samples N.
!           GM          Product of the gravitatonal constant and the planet's 
!                       mass.
!           r0          Reference radius of potential coefficients.
!           a           The semimajor axis of the flattened ellipsoid.
!           f           Flattening of the planet.

!       IN, OPTIONAL
!           sampling    (1) Grid is N latitudes by N longitudes (default).
!                       (2) Grid is N by 2N. The higher frequencies resulting
!                       from this oversampling in longitude are discarded, and 
!                       hence not aliased into lower frequencies.
!           lmax_calc   The maximum spherical harmonic degree to evaluate
!                       the coefficients up to.
!
!       OUT
!           Vxx         x-x component of the gravity gradient tensor.
!           Vyy         y-y component of the gravity gradient tensor.
!           Vzz         z-z component of the gravity gradient tensor.
!           Vxy         x-y component of the gravity gradient tensor.
!           Vxz         x-z component of the gravity gradient tensor.
!           Vyz         y-z component of the gravity gradient tensor.
!           N           Number of samples in latitude. Number of samples in 
!                       longitude is N when sampling is 1 (default), and is 
!                       2N when sampling is 2.
!
!   Notes:
!       1.  If lmax is greater than the the maximum spherical harmonic
!           degree of the input file, Cilm will be ZERO PADDED!
!           (i.e., those degrees after lmax are assumed to be zero).
!       2.  Latitude is geocentric latitude.
!
!   Dependencies:   FFTW3
!
!   Copyright (c) 2015, Mark A. Wieczorek
!   All rights reserved.
!
!-------------------------------------------------------------------------------
    use FFTW3
    
    implicit none
    
    real*8, intent(in) :: cilm(:,:,:), gm, r0, a, f
    real*8, intent(out) :: vxx(:,:), vyy(:,:), vzz(:,:), vxy(:,:), vxz(:,:), &
                            vyz(:,:)
    integer, intent(in) :: lmax
    integer, intent(out) :: n
    integer, intent(in), optional :: sampling, lmax_calc
    integer :: l, m, i, l1, m1, lmax_comp, i_eq, i_s, astat(4), nlong
    real*8 ::  grid(4*lmax+4), pi, theta, scalef, rescalem, u, p, dp, dp2, &
                dp2s, pmm, sint, pm1, pm2, z, tempr, r_ex, lat, &
                prefactor(lmax), coefr0, coefrs0, coeft0, coefts0, coefp0, &
                coefps0, coefrr0, coefrrs0, coefrt0, coefrts0, coefrp0, &
                coefrps0, coeftp0, coeftps0, coefpp0, coefpps0, coeftt0, &
                coeftts0
    complex*16 :: coef(2*lmax+3), coefr(2*lmax+3), coefrs(2*lmax+3), &
                coeft(2*lmax+3), coefts(2*lmax+3), coefp(2*lmax+3), &
                coefps(2*lmax+3), tempc, coefrr(2*lmax+3), coefrrs(2*lmax+3), &
                coefrt(2*lmax+3), coefrts(2*lmax+3), coefrp(2*lmax+3), &
                coefrps(2*lmax+3), coeftp(2*lmax+3), coeftps(2*lmax+3), &
                coefpp(2*lmax+3), coefpps(2*lmax+3), coeftt(2*lmax+3), &
                coeftts(2*lmax+3)
    integer*8 :: plan
    real*8, save, allocatable :: ff1(:,:), ff2(:,:), sqr(:)
    integer*1, save, allocatable :: fsymsign(:,:)
    integer, save :: lmax_old = 0
    
    n = 2 * lmax+2
    
    if (present(sampling)) then
        if (sampling /= 1 .and. sampling /=2) then
            print*, "Error --- MakeGravGradGridDH"
            print*, "Optional parameter SAMPLING must be 1 (N by N) " // &
                    "or 2 (N by 2N)."
            print*, "Input value is ", sampling
            stop
        end if
    end if
    
    if (size(cilm(:,1,1)) < 2) then
        print*, "Error --- MakeGravGradGridDH"
        print*, "CILM must be dimensioned as (2, *, *)."
        print*, "Input dimension is ", size(cilm(:,1,1)), size(cilm(1,:,1)), &
                size(cilm(1,1,:))
        stop
    end if 
    
    if (present(sampling)) then
    
        if (sampling == 1) then
            nlong = n
            
        else
            nlong = 2 * n
            
        end if
    else
        nlong = n
        
    end if
    
    if (size(vxx(:,1)) < n .or. size(vxx(1,:)) < nlong .or. size(vyy(:,1)) < n &
            .or. size(vyy(1,:)) < nlong .or. size(vzz(:,1)) < n .or. &
            size(vzz(1,:)) < nlong .or. size(vxy(:,1)) < n .or. &
            size(vxy(1,:)) < nlong .or. size(vxz(:,1)) < n .or. &
            size(vxz(1,:)) < nlong  .or. size(vyz(:,1)) < n .or. &
            size(vyz(1,:)) < nlong) then
        print*, "Error --- MakeGravGradGridDH"
        
        if (present(sampling)) then
            if (sampling == 1) then
                print*, "VXX, VYY, VZZ, VXY, VXZ, and VYZ must be " // &
                        "dimensioned as (N, N) where N is ", n
                        
            else if (sampling == 2) then
                print*, "VXX, VYY, VZZ, VXY, VXZ, and VYZ must be " // &
                        "dimensioned as (N, 2N) where N is ", n
                        
            end if
            
        else
            print*, "VXX, VYY, VZZ, VXY, VXZ, and VYZ must be dimensioned " // &
                    "as (N, N) where N is ", n
                    
        end if
        
        print*, "Input dimensions are ", size(vxx(:,1)), size(vxx(1,:)), &
            size(vyy(:,1)), size(vyy(1,:)), size(vzz(:,1)), size(vzz(1,:)), &
            size(vxy(:,1)), size(vxy(1,:)), size(vxz(:,1)), size(vxz(1,:)), &
            size(vyz(:,1)), size(vyz(1,:))
        stop
        
    end if
    
    if (cilm(1,1,1) /= 1.0d0 .and. f /= 0.0d0) then
        print*, "Warning --- MakeGravGradGridDH"
        print*, "The degree-0 term of the spherical harmonic coefficients is not equal to 1."
        print*, "The variation in gravity resulting from variations in radius of the flattened"
        print*, "ellipsoid will not be taken into account."
        print*, "C00 = ", cilm(1,1,1)
        print*, "F = ", f
    end if
            
    pi = acos(-1.0d0)
    
    scalef = 1.0d-280
    
    if (present(lmax_calc)) then
        if (lmax_calc > lmax) then
            print*, "Error --- MakeGravGradGridDH"
            print*, "LMAX_CALC must be less than or equal to LMAX."
            print*, "LMAX = ", lmax
            print*, "LMAX_CALC = ", lmax_calc
            stop
            
        else
            lmax_comp = min(lmax, size(cilm(1,1,:))-1, size(cilm(1,:,1))-1, &
                            lmax_calc)
                            
        end if
        
    else
        lmax_comp = min(lmax, size(cilm(1,1,:))-1, size(cilm(1,:,1))-1)
        
    end if

    !---------------------------------------------------------------------------
    !
    !   Calculate recursion constants used in computing the Legendre polynomials
    !
    !---------------------------------------------------------------------------
    if (lmax_comp /= lmax_old) then
        
        if (allocated (sqr)) deallocate (sqr)
        if (allocated (ff1)) deallocate (ff1)
        if (allocated (ff2)) deallocate (ff2)
        if (allocated (fsymsign)) deallocate (fsymsign)
        
        allocate (sqr(2 * lmax_comp + 1), stat=astat(1))
        allocate (ff1(lmax_comp+1,lmax_comp+1), stat=astat(2))
        allocate (ff2(lmax_comp+1,lmax_comp+1), stat=astat(3))
        allocate (fsymsign(lmax_comp+1,lmax_comp+1), stat=astat(4))
        
        if (sum(astat(1:4)) /= 0) then
            print*, "Error --- MakeGravGradGridDH"
            print*, "Problem allocating arrays SQR, FF1, FF2, or FSYMSIGN", &
                    astat(1), astat(2), astat(3), astat(4)
            stop
        end if
        
        !-----------------------------------------------------------------------
        !
        !   Calculate signs used for symmetry of Legendre functions about 
        !   equator. For the first derivative in theta, these signs are 
        !   reversed.
        !
        !-----------------------------------------------------------------------
        do l = 0, lmax_comp, 1
            do m = 0, l, 1
                if (mod(l-m,2) == 0) then
                    fsymsign(l+1,m+1) = 1
                    
                else
                    fsymsign(l+1,m+1) = -1
                    
                end if
                
            end do
            
        end do
            
        !-----------------------------------------------------------------------
        !
        !   Precompute square roots of integers that are used several times.
        !
        !-----------------------------------------------------------------------
    
        do l=1, 2 * lmax_comp + 1
            sqr(l) = sqrt(dble(l))
        end do

        !-----------------------------------------------------------------------
        !
        !   Precompute multiplicative factors used in recursion relationships
        !       P(l,m) = x*f1(l,m)*P(l-1,m) - P(l-2,m)*f2(l,m)
        !       k = l*(l+1)/2 + m + 1
        !   Note that prefactors are not used for the case when m=l as a 
        !   different recursion is used. Furthermore, for m=l-1, Plmbar(l-2,m) 
        !   is assumed to be zero.
        !
        !-----------------------------------------------------------------------
    
        if (lmax_comp /= 0) then
            ff1(2,1) = sqr(3)
            ff2(2,1) = 0.0d0
            
        end if
                
        do l=2, lmax_comp, 1
            ff1(l+1,1) = sqr(2*l-1) * sqr(2*l+1) / dble(l)
            ff2(l+1,1) = dble(l-1) * sqr(2*l+1) / sqr(2*l-3) / dble(l)
            
            do m = 1, l-2, 1
                ff1(l+1,m+1) = sqr(2*l+1) * sqr(2*l-1) / sqr(l+m) / sqr(l-m)
                ff2(l+1,m+1) = sqr(2*l+1) * sqr(l-m-1) * sqr(l+m-1) &
                            / sqr(2*l-3) / sqr(l+m) / sqr(l-m) 
            end do
            
            ff1(l+1,l) = sqr(2*l+1) * sqr(2*l-1) / sqr(l+m) / sqr(l-m)
            ff2(l+1,l) = 0.0d0
            
        end do
    
        lmax_old = lmax_comp
    
    end if
    
    !---------------------------------------------------------------------------
    !
    !   Determine Clms one l at a time by intergrating over latitude.
    !   
    !---------------------------------------------------------------------------
    
    call dfftw_plan_dft_c2r_1d(plan, nlong, coef(1:nlong/2+1), grid(1:nlong), &
                                FFTW_MEASURE)

    i_eq = n/2 + 1  ! Index correspondong to zero latitude
    
    do i = 1, i_eq - 1, 1
    
        i_s = 2 * i_eq -i
    
        theta = pi * dble(i-1)/dble(n)
        z = cos(theta)
        u = sqrt( (1.0d0-z) * (1.0d0+z) )
        sint = sin(theta)
        lat = pi / 2.0d0 - theta
        
        if (i == 1) then      ! Reference ellipsoid radius
            r_ex = a * (1.0d0 - f)  
            
        else
            r_ex = (1.0d0 + tan(lat)**2) / &
                    (1.0d0  + tan(lat)**2 / (1.0d0 - f)**2 )
            r_ex = a * sqrt(r_ex)
            
        end if

        coefr(1:lmax+1) = dcmplx(0.0d0,0.0d0)
        coefr0 = 0.0d0
        coefrs(1:lmax+1) = dcmplx(0.0d0,0.0d0)
        coefrs0 = 0.0d0
        
        coefrr(1:lmax+1) = dcmplx(0.0d0,0.0d0)
        coefrr0 = 0.0d0
        coefrrs(1:lmax+1) = dcmplx(0.0d0,0.0d0)
        coefrrs0 = 0.0d0
        
        coeft(1:lmax+1) = dcmplx(0.0d0,0.0d0)
        coeft0 = 0.0d0
        coefts(1:lmax+1) = dcmplx(0.0d0,0.0d0)
        coefts0 = 0.0d0
        
        coeftp(1:lmax+1) = dcmplx(0.0d0,0.0d0)
        coeftp0 = 0.0d0
        coeftps(1:lmax+1) = dcmplx(0.0d0,0.0d0)
        coeftps0 = 0.0d0
        
        coefp(1:lmax+1) = dcmplx(0.0d0,0.0d0)
        coefp0 = 0.0d0
        coefps(1:lmax+1) = dcmplx(0.0d0,0.0d0)
        coefps0 = 0.0d0
        
        coefrt(1:lmax+1) = dcmplx(0.0d0,0.0d0)
        coefrt0 = 0.0d0
        coefrts(1:lmax+1) = dcmplx(0.0d0,0.0d0)
        coefrts0 = 0.0d0
        
        coefrp(1:lmax+1) = dcmplx(0.0d0,0.0d0)
        coefrp0 = 0.0d0
        coefrps(1:lmax+1) = dcmplx(0.0d0,0.0d0)
        coefrps0 = 0.0d0
        
        coefpp(1:lmax+1) = dcmplx(0.0d0,0.0d0)
        coefpp0 = 0.0d0
        coefpps(1:lmax+1) = dcmplx(0.0d0,0.0d0)
        coefpps0 = 0.0d0
        
        coeftt(1:lmax+1) = dcmplx(0.0d0,0.0d0)
        coeftt0 = 0.0d0
        coeftts(1:lmax+1) = dcmplx(0.0d0,0.0d0)
        coeftts0 = 0.0d0
                
        pm2 = 1.0d0

        tempr =  -cilm(1,1,1) * pm2 ! l = 0 
        coefr0 = coefr0 + tempr
        coefrs0 = coefrs0 + tempr   ! fsymsign is always 1 for l=m=0
        
        tempr =  2 * cilm(1,1,1) * pm2  ! l = 0 
        coefrr0 = coefrr0 + tempr
        coefrrs0 = coefrrs0 + tempr
        
        ! derivative in theta and phi of l=0 term is 0, so no need to 
        ! calculate this
                
        if (lmax_comp /= 0) then    ! l = 1
            prefactor(1) = r0/r_ex
            
            do l = 2, lmax_comp, 1
                prefactor(l) = prefactor(l-1) * r0/r_ex
            end do
            
            pm1 =  ff1(2,1) * z 
            
            ! -2 = (l+1) prefactor
            tempr = cilm(1,2,1) * pm1 * (-2) * prefactor(1) 
            coefr0 = coefr0 + tempr
            coefrs0 = coefrs0 - tempr   ! fsymsign = -1
            
            ! 6 = (l+1)*(l+2) prefactor
            tempr = cilm(1,2,1) * pm1 * (6) * prefactor(1)  
            coefrr0 = coefrr0 + tempr
            coefrrs0 = coefrrs0 - tempr     ! fsymsign = -1
            
            ! dp is the first derivative with respect to Z
            dp = ff1(2,1)
            tempr = cilm(1,2,1) * dp * prefactor(1)
            coeft0 = coeft0 + tempr
            coefts0 = coefts0 + tempr   ! reverse fsymsign
            
            ! -2 = (l+1) prefactor
            tempr = cilm(1,2,1) * dp * (-2) * prefactor(1)  
            coefrt0 = coefrt0 + tempr
            coefrts0 = coefrts0 + tempr     ! reverse fsymsign
            
            ! dp2 is the second derivative with respect to THETA. Must 
            ! multiply dp by -sin(theta) in the recurrence relation
            ! note that the first term is symmetric about the equator accoring 
            ! the fsymsign, whereas the second term differs
            ! by -1.
            dp2 = -2 * pm1 + z * dp
            dp2s = 2 * pm1 + z * dp
            tempr = cilm(1,2,1) * prefactor(1)
            coeftt0 = coeftt0 + tempr * dp2
            coeftts0 = coeftts0 + tempr * dp2s
            
        end if
                
        do l = 2, lmax_comp, 1
            l1 = l+1
            p = ff1(l1,1) * z * pm1 - ff2(l1,1) * pm2
            tempr = cilm(1,l1,1) * p * (-l1) * prefactor(l)
            coefr0 = coefr0 + tempr
            coefrs0 = coefrs0 + tempr * fsymsign(l1,1)
            
            tempr = cilm(1,l1,1) * p * (l1) * (l1+1) * prefactor(l)
            coefrr0 = coefrr0 + tempr
            coefrrs0 = coefrrs0 + tempr * fsymsign(l1,1)
            
            dp = l * ( sqr(2*l+1) / sqr(2*l-1) * pm1 - z * p ) / u**2
            tempr = cilm(1,l1,1) * dp * prefactor(l)
            coeft0 = coeft0 + tempr
            coefts0 = coefts0 - tempr * fsymsign(l1,1)  ! reverse fsymsign
            
            tempr = cilm(1,l1,1) * dp * (-l1) * prefactor(l)
            coefrt0 = coefrt0 + tempr
            coefrts0 = coefrts0 - tempr * fsymsign(l1,1)
            
            dp2 = -l*l1 * p + z * dp
            dp2s = -l*l1 * p * fsymsign(l1,1) - z * dp * fsymsign(l1,1)
            tempr = cilm(1,l1,1) * prefactor(l)
            coeftt0 = coeftt0 + tempr * dp2
            coeftts0 = coeftts0 + tempr * dp2s
            
            pm2 = pm1
            pm1 = p
            
        end do
                
        pmm = sqr(2) * scalef
                
        rescalem = 1.0d0/scalef
            
        do m = 1, lmax_comp-1, 1
            
            m1 = m+1
            rescalem = rescalem * u
                    
            pmm = pmm * sqr(2*m+1) / sqr(2*m)
            pm2 = pmm
            
            tempc = dcmplx(cilm(1,m1,m1), - cilm(2,m1,m1)) * pm2 * (-m-1) &
                    * prefactor(m)    ! (m,m)
            coefr(m1) = coefr(m1) + tempc
            coefrs(m1) = coefrs(m1) + tempc ! fsymsign = 1
            
            tempc = dcmplx(cilm(1,m1,m1), - cilm(2,m1,m1)) * pm2 * (m+1) &
                            * (m+2) * prefactor(m) ! (m,m)
            coefrr(m1) = coefrr(m1) + tempc
            coefrrs(m1) = coefrrs(m1) + tempc ! fsymsign = 1
            
            tempc = dcmplx(cilm(2,m1,m1), cilm(1,m1,m1)) * pm2 &
                            * prefactor(m) * m ! (m,m)
            coefp(m1) = coefp(m1) + tempc
            coefps(m1) = coefps(m1) + tempc ! fsymsign = 1
            
            tempc = -dcmplx(cilm(1,m1,m1), - cilm(2,m1,m1)) * pm2 &
                            * prefactor(m) * m**2 ! (m,m)
            coefpp(m1) = coefpp(m1) + tempc
            coefpps(m1) = coefpps(m1) + tempc ! fsymsign = 1
            
            tempc = dcmplx(cilm(2,m1,m1), cilm(1,m1,m1)) * pm2 *(-m-1) &
                            * prefactor(m) * m ! (m,m)
            coefrp(m1) = coefrp(m1) + tempc
            coefrps(m1) = coefrps(m1) + tempc ! fsymsign = 1
            
            dp = -m * z * pm2 / u**2
            tempc = dcmplx(cilm(1,m1,m1), - cilm(2,m1,m1)) * dp &
                            * prefactor(m)  ! (m,m)
            coeft(m1) = coeft(m1) + tempc
            coefts(m1) = coefts(m1) - tempc ! reverse fsymsign
            
            tempc = dcmplx(cilm(1,m1,m1), - cilm(2,m1,m1)) * dp * (-m-1) &
                            * prefactor(m) ! (m,m)
            coefrt(m1) = coefrt(m1) + tempc
            coefrts(m1) = coefrts(m1) - tempc ! reverse fsymsign
            
            tempc = dcmplx(cilm(2,m1,m1), cilm(1,m1,m1)) * dp &
                            * prefactor(m) * m ! (m,m)
            coeftp(m1) = coeftp(m1) + tempc
            coeftps(m1) = coeftps(m1) - tempc ! reverse fsymsign
            
            dp2 = -(m*m1 - (m**2)/u**2) * pm2 + z * dp  
            dp2s = -(m*m1 - (m**2)/u**2) * pm2 - z * dp 
            tempc = dcmplx(cilm(1,m1,m1), - cilm(2,m1,m1)) &
                            * prefactor(m)   ! (m,m)
            coeftt(m1) = coeftt(m1) + tempc * dp2
            coeftts(m1) = coeftts(m1) + tempc * dp2s
                                        
            pm1 = z * ff1(m1+1,m1) * pm2        
            tempc = dcmplx(cilm(1,m1+1,m1), - cilm(2,m1+1,m1)) * pm1 * (-m-2) &
                            * prefactor(m+1)  ! (m+1,m)
            coefr(m1) = coefr(m1) + tempc   
            coefrs(m1) = coefrs(m1) - tempc ! fsymsign = -1
            
            tempc = dcmplx(cilm(1,m1+1,m1), - cilm(2,m1+1,m1)) * pm1 * (m+2) &
                            * (m+3) * prefactor(m+1)   ! (m+1,m)
            coefrr(m1) = coefrr(m1) + tempc 
            coefrrs(m1) = coefrrs(m1) - tempc ! fsymsign = -1
            
            tempc = dcmplx(cilm(2,m1+1,m1), cilm(1,m1+1,m1)) * pm1 &
                            * prefactor(m+1) * m ! (m+1,m)
            coefp(m1) = coefp(m1) + tempc   
            coefps(m1) = coefps(m1) - tempc ! fsymsign = -1
            
            tempc = -dcmplx(cilm(1,m1+1,m1), - cilm(2,m1+1,m1)) * pm1 &
                            * prefactor(m+1) * m**2   ! (m+1,m)
            coefpp(m1) = coefpp(m1) + tempc 
            coefpps(m1) = coefpps(m1) - tempc ! fsymsign = -1
            
            tempc = dcmplx(cilm(2,m1+1,m1), cilm(1,m1+1,m1)) * pm1 * (-m-2) &
                            * prefactor(m+1) * m    ! (m+1,m)
            coefrp(m1) = coefrp(m1) + tempc 
            coefrps(m1) = coefrps(m1) - tempc ! fsymsign = -1
            
            dp =  ( sqr(2*m+3) * pmm - z * (m+1) * pm1) / u**2
            tempc = dcmplx(cilm(1,m1+1,m1), - cilm(2,m1+1,m1)) * dp &
                            * prefactor(m+1)    ! (m+1,m)
            coeft(m1) = coeft(m1) + tempc   
            coefts(m1) = coefts(m1) + tempc ! reverse fsymsign
            
            tempc = dcmplx(cilm(1,m1+1,m1), - cilm(2,m1+1,m1)) * dp &
                            * (-m-2) * prefactor(m+1)   ! (m+1,m)
            coefrt(m1) = coefrt(m1) + tempc 
            coefrts(m1) = coefrts(m1) + tempc   ! reverse fsymsign
            
            tempc = dcmplx(cilm(2,m1+1,m1), cilm(1,m1+1,m1)) * dp &
                            * prefactor(m+1) * m  ! (m+1,m)
            coeftp(m1) = coeftp(m1) + tempc 
            coeftps(m1) = coeftps(m1) + tempc ! reverse fsymsign
            
            dp2 = -(m1*(m1+1) - (m**2)/u**2) * pm1 + z * dp
            dp2s = (m1*(m1+1) - (m**2)/u**2) * pm1 + z * dp
            tempc = dcmplx(cilm(1,m1+1,m1), - cilm(2,m1+1,m1)) &
                            * prefactor(m+1) ! (m+1,m)
            coeftt(m1) = coeftt(m1) + tempc  * dp2
            coeftts(m1) = coeftts(m1) + tempc * dp2s
                    
            do l = m + 2, lmax_comp, 1
                l1 = l+1
                p = z * ff1(l1,m1) * pm1 - ff2(l1,m1) * pm2
                tempc = dcmplx(cilm(1,l1,m1), - cilm(2,l1,m1)) * p * (-l1) &
                                * prefactor(l)
                coefr(m1) = coefr(m1) + tempc
                coefrs(m1) = coefrs(m1) + tempc * fsymsign(l1,m1)
                
                tempc = dcmplx(cilm(1,l1,m1), - cilm(2,l1,m1)) * p * (l1) &
                                * (l1+1) * prefactor(l)
                coefrr(m1) = coefrr(m1) + tempc
                coefrrs(m1) = coefrrs(m1) + tempc * fsymsign(l1,m1)
                
                tempc = dcmplx(cilm(2,l1,m1), cilm(1,l1,m1)) * p &
                                * prefactor(l) * m
                coefp(m1) = coefp(m1) + tempc
                coefps(m1) = coefps(m1) + tempc * fsymsign(l1,m1)
                
                tempc = -dcmplx(cilm(1,l1,m1), - cilm(2,l1,m1)) * p &
                                * prefactor(l) * m**2
                coefpp(m1) = coefpp(m1) + tempc
                coefpps(m1) = coefpps(m1) + tempc * fsymsign(l1,m1)
                
                tempc = dcmplx(cilm(2,l1,m1), cilm(1,l1,m1)) * p * (-l1) &
                                * prefactor(l) * m
                coefrp(m1) = coefrp(m1) + tempc
                coefrps(m1) = coefrps(m1) + tempc * fsymsign(l1,m1)
                
                dp = ( sqr(2*l+1) * sqr(l-m) * sqr(l+m) / sqr(2*l-1) &
                                * pm1 - z * l * p) / u**2
                tempc = dcmplx(cilm(1,l1,m1), - cilm(2,l1,m1)) * dp &
                                * prefactor(l)
                coeft(m1) = coeft(m1) + tempc
                ! reverse fsymsign
                coefts(m1) = coefts(m1) - tempc * fsymsign(l1,m1) 
                
                tempc = dcmplx(cilm(1,l1,m1), - cilm(2,l1,m1)) * dp * (-l1) &
                                * prefactor(l)
                coefrt(m1) = coefrt(m1) + tempc
                ! reverse fsymsign
                coefrts(m1) = coefrts(m1) - tempc * fsymsign(l1,m1) 
                
                tempc = dcmplx(cilm(2,l1,m1), cilm(1,l1,m1)) * dp &
                                * prefactor(l) * m
                coeftp(m1) = coeftp(m1) + tempc
                ! reverse fsymsign
                coeftps(m1) = coeftps(m1) - tempc * fsymsign(l1,m1) 
                
                dp2 = -(l*l1 -(m**2)/u**2) * p + z * dp
                dp2s = -(l*l1 -(m**2)/u**2) * p * fsymsign(l1,m1) &
                                - z * dp * fsymsign(l1,m1) 
                tempc = dcmplx(cilm(1,l1,m1), - cilm(2,l1,m1)) * prefactor(l)
                coeftt(m1) = coeftt(m1) + tempc * dp2
                coeftts(m1) = coeftts(m1) +  tempc * dp2s

                pm2 = pm1
                pm1 = p
                    
            end do
                    
            coefr(m1) = coefr(m1) * rescalem
            coefrs(m1) = coefrs(m1) * rescalem
                        
            coefrr(m1) = coefrr(m1) * rescalem
            coefrrs(m1) = coefrrs(m1) * rescalem
            
            coeft(m1) = coeft(m1) * rescalem
            coefts(m1) = coefts(m1) * rescalem
            
            coeftt(m1) = coeftt(m1) * rescalem
            coeftts(m1) = coeftts(m1) * rescalem
            
            coefrt(m1) = coefrt(m1) * rescalem
            coefrts(m1) = coefrts(m1) * rescalem
            
            coefp(m1) = coefp(m1) * rescalem
            coefps(m1) = coefps(m1) * rescalem
            
            coefpp(m1) = coefpp(m1) * rescalem
            coefpps(m1) = coefpps(m1) * rescalem
            
            coeftp(m1) = coeftp(m1) * rescalem
            coeftps(m1) = coeftps(m1) * rescalem
            
            coefrp(m1) = coefrp(m1) * rescalem
            coefrps(m1) = coefrps(m1) * rescalem
                    
        end do           
                                
        if (lmax_comp /= 0) then

            rescalem = rescalem * u
                
            pmm = pmm * sqr(2*lmax_comp+1) / sqr(2*lmax_comp) * rescalem                    
            tempc = dcmplx(cilm(1,lmax_comp+1,lmax_comp+1), &
                            - cilm(2,lmax_comp+1,lmax_comp+1)) * pmm &
                            * (-lmax_comp-1) * prefactor(lmax_comp)
            coefr(lmax_comp+1) = coefr(lmax_comp+1) + tempc
            coefrs(lmax_comp+1) = coefrs(lmax_comp+1) + tempc ! fsymsign = 1
            
            tempc = dcmplx(cilm(1,lmax_comp+1,lmax_comp+1), &
                            - cilm(2,lmax_comp+1,lmax_comp+1)) * pmm &
                            * (lmax_comp+1) * (lmax_comp+2) &
                            * prefactor(lmax_comp)
            coefrr(lmax_comp+1) = coefrr(lmax_comp+1) + tempc
            coefrrs(lmax_comp+1) = coefrrs(lmax_comp+1) + tempc ! fsymsign = 1
        
            tempc = dcmplx(cilm(2,lmax_comp+1,lmax_comp+1), &
                            cilm(1,lmax_comp+1,lmax_comp+1)) * pmm &
                            * prefactor(lmax_comp) * lmax_comp
            coefp(lmax_comp+1) = coefp(lmax_comp+1) + tempc
            coefps(lmax_comp+1) = coefps(lmax_comp+1) + tempc   ! fsymsign = 1
            
            tempc = -dcmplx(cilm(1,lmax_comp+1,lmax_comp+1), &
                            - cilm(2,lmax_comp+1,lmax_comp+1)) * pmm &
                            * prefactor(lmax_comp) * lmax_comp**2
            coefpp(lmax_comp+1) = coefpp(lmax_comp+1) + tempc
            coefpps(lmax_comp+1) = coefpps(lmax_comp+1) + tempc ! fsymsign = 1
            
            tempc = dcmplx(cilm(2,lmax_comp+1,lmax_comp+1), &
                            cilm(1,lmax_comp+1,lmax_comp+1)) * pmm &
                            * (-lmax_comp-1) * prefactor(lmax_comp) * lmax_comp
            coefrp(lmax_comp+1) = coefrp(lmax_comp+1) + tempc
            coefrps(lmax_comp+1) = coefrps(lmax_comp+1) + tempc ! fsymsign = 1
        
            dp = -lmax_comp * z * pmm / u**2
            tempc = dcmplx(cilm(1,lmax_comp+1,lmax_comp+1), &
                            - cilm(2,lmax_comp+1,lmax_comp+1)) * dp &
                            * prefactor(lmax_comp)
            coeft(lmax_comp+1) = coeft(lmax_comp+1) + tempc
            ! reverse fsymsign
            coefts(lmax_comp+1) = coefts(lmax_comp+1) - tempc   
            
            tempc = dcmplx(cilm(1,lmax_comp+1,lmax_comp+1), &
                            - cilm(2,lmax_comp+1,lmax_comp+1)) * dp &
                            * (-lmax_comp-1) * prefactor(lmax_comp)
            coefrt(lmax_comp+1) = coefrt(lmax_comp+1) + tempc
            ! reverse fsymsign
            coefrts(lmax_comp+1) = coefrts(lmax_comp+1) - tempc 
            
            tempc = dcmplx(cilm(2,lmax_comp+1,lmax_comp+1), &
                            cilm(1,lmax_comp+1,lmax_comp+1)) * dp &
                            * prefactor(lmax_comp) * lmax_comp
            coeftp(lmax_comp+1) = coeftp(lmax_comp+1) + tempc
            coeftps(lmax_comp+1) = coeftps(lmax_comp+1) - tempc
            
            dp2 = -(lmax_comp*(lmax_comp+1)-(lmax_comp**2)/u**2) * pmm + z * dp
            dp2s = -(lmax_comp*(lmax_comp+1)-(lmax_comp**2)/u**2) * pmm - z * dp
            tempc = dcmplx(cilm(1,lmax_comp+1,lmax_comp+1), &
                            - cilm(2,lmax_comp+1,lmax_comp+1)) &
                            * prefactor(lmax_comp)
            coeftt(lmax_comp+1) = coeftt(lmax_comp+1) + tempc * dp2
            coeftts(lmax_comp+1) = coeftts(lmax_comp+1) + tempc * dp2s
        
        end if
        
        ! Note that the first angular derivatives are with repsect to z, 
        ! but that the second is with respect to theta.
        
        coefr0 = coefr0 * gm/r_ex**2
        coefr(2:lmax+1) = coefr(2:lmax+1) * gm/r_ex**2
        
        coefrs0 = coefrs0 * gm/r_ex**2
        coefrs(2:lmax+1) = coefrs(2:lmax+1) * gm/r_ex**2
        
        coefrt0 =  -sint * coefrt0 * gm/r_ex**2 
        coefrt(2:lmax+1) =  -sint * coefrt(2:lmax+1) * gm/r_ex**2
        
        coefrts0 =  -sint * coefrts0 * gm/r_ex**2
        coefrts(2:lmax+1) =  -sint * coefrts(2:lmax+1) * gm/r_ex**2
        
        coefrr0 = coefrr0 * gm/r_ex**3
        coefrr(2:lmax+1) = coefrr(2:lmax+1) * gm/r_ex**3
        
        coefrrs0 = coefrrs0 * gm/r_ex**3
        coefrrs(2:lmax+1) = coefrrs(2:lmax+1) * gm/r_ex**3
        
        coeft0 = -sint * coeft0 * gm/r_ex 
        coeft(2:lmax+1) = -sint * coeft(2:lmax+1) * gm/r_ex 
        
        coefts0 = -sint * coefts0 * gm/r_ex 
        coefts(2:lmax+1) = -sint * coefts(2:lmax+1) * gm/r_ex
        
        coeftt0 = coeftt0 * gm/r_ex 
        coeftt(2:lmax+1) = coeftt(2:lmax+1) * gm/r_ex 
        
        coeftts0 = coeftts0 * gm/r_ex 
        coeftts(2:lmax+1) = coeftts(2:lmax+1) * gm/r_ex
        
        coeftp0 = -sint * coeftp0 * gm/r_ex 
        coeftp(2:lmax+1) = -sint * coeftp(2:lmax+1) * gm/r_ex 
        
        coeftps0 = -sint * coeftps0 * gm/r_ex 
        coeftps(2:lmax+1) = -sint * coeftps(2:lmax+1) * gm/r_ex
        
        coefp0 = coefp0 * gm/r_ex
        coefp(2:lmax+1) = coefp(2:lmax+1) * gm/r_ex
        
        coefps0 = coefps0 * gm/r_ex
        coefps(2:lmax+1) = coefps(2:lmax+1) * gm/r_ex
        
        coefpp0 = coefpp0 * gm/r_ex
        coefpp(2:lmax+1) = coefpp(2:lmax+1) * gm/r_ex
        
        coefpps0 = coefpps0 * gm/r_ex
        coefpps(2:lmax+1) = coefpps(2:lmax+1) * gm/r_ex
        
        coefrp0 = coefrp0 * gm/r_ex**2
        coefrp(2:lmax+1) = coefrp(2:lmax+1) * gm/r_ex**2
        
        coefrps0 = coefrps0 * gm/r_ex**2
        coefrps(2:lmax+1) = coefrps(2:lmax+1) * gm/r_ex**2
        
        ! Vzz = Vrr
        coef(1) = dcmplx(coefrr0,0.0d0)
        coef(2:lmax+1) = coefrr(2:lmax+1) / 2.0d0
        
        if (present(sampling)) then
            if (sampling == 2) then
                coef(lmax+2:2*lmax+3) = dcmplx(0.0d0,0.0d0) 
            end if
            
        end if       
        
        call dfftw_execute(plan)    ! take fourier transform
        vzz(i,1:nlong) = grid(1:nlong)
                                
        if (i == 1) then
            vyy(1,1:nlong) = 0.0d0  ! These derivatives are 
            vzz(1,1:nlong) = 0.0d0  ! undefined at the pole
            vxy(1,1:nlong) = 0.0d0  
            vxz(1,1:nlong) = 0.0d0  
            vyz(1,1:nlong) = 0.0d0
                    
        else
            ! Vxx = 1/r Vr + 1/r^2 Vtt
            coef(1) = dcmplx(coefr0/r_ex + coeftt0/r_ex**2,0.0d0)
            coef(2:lmax+1) = (coefr(2:lmax+1)/r_ex &
                            + coeftt(2:lmax+1)/r_ex**2 ) / 2.0d0
                            
            if (present(sampling)) then
                if (sampling == 2) then
                    coef(lmax+2:2*lmax+3) = dcmplx(0.0d0,0.0d0) 
                end if
            end if 
                  
            call dfftw_execute(plan)    ! take fourier transform
            vxx(i,1:nlong) = grid(1:nlong)
                    
            ! Vyy = 1/r Vr + 1/r^2 /tan(t) Vt + 1/r^2 /sin(t)^2 Vpp
            coef(1) = dcmplx(coefr0/r_ex + coeft0/(r_ex**2)/tan(theta) &
                            + coefpp0/(r_ex**2)/u**2 ,0.0d0)
            coef(2:lmax+1) = (coefr(2:lmax+1)/r_ex &
                            + coeft(2:lmax+1)/(r_ex**2)/tan(theta) + &
                            coefpp(2:lmax+1)/(r_ex**2)/u**2 ) / 2.0d0
                            
            if (present(sampling)) then
                if (sampling == 2) then
                    coef(lmax+2:2*lmax+3) = dcmplx(0.0d0,0.0d0) 
                end if
            end if       
            
            call dfftw_execute(plan)    ! take fourier transform
            vyy(i,1:nlong) = grid(1:nlong)
                    
            ! Vxy = 1/r^2 /sin(t) Vtp - cos(t)/sin(t)^2 /r^2 Vp
            coef(1) = dcmplx(coeftp0/sint/r_ex**2 &
                            - coefp0/(r_ex**2)*z/u**2, 0.0d0)
            coef(2:lmax+1) = (coeftp(2:lmax+1)/sint/r_ex**2 &
                            - coefp(2:lmax+1)/(r_ex**2)*z/u**2 ) / 2.0d0
                            
            if (present(sampling)) then
                if (sampling == 2) then
                    coef(lmax+2:2*lmax+3) = dcmplx(0.0d0,0.0d0) 
                end if
            end if
                   
            call dfftw_execute(plan)    ! take fourier transform
            vxy(i,1:nlong) = grid(1:nlong)
                    
            ! Vxz = 1/r^2 Vt - 1/r Vrt
            coef(1) = dcmplx(coeft0/r_ex**2 - coefrt0/r_ex, 0.0d0)
            coef(2:lmax+1) = (coeft(2:lmax+1)/r_ex**2 &
                            - coefrt(2:lmax+1)/r_ex ) / 2.0d0
                            
            if (present(sampling)) then
                if (sampling == 2) then
                    coef(lmax+2:2*lmax+3) = dcmplx(0.0d0,0.0d0) 
                end if
            end if
                  
            call dfftw_execute(plan)    ! take fourier transform
            vxz(i,1:nlong) = grid(1:nlong)
                    
            ! Vyz = 1/r^2 /sin(t) Vp - 1/r /sin(t) Vrp
            coef(1) = dcmplx(coefp0/sint/r_ex**2 - coefrp0/sint/r_ex, 0.0d0)
            coef(2:lmax+1) = (coefp(2:lmax+1)/sint/r_ex**2 &
                            - coefrp(2:lmax+1)/sint/r_ex ) / 2.0d0
                            
            if (present(sampling)) then
                if (sampling == 2) then
                    coef(lmax+2:2*lmax+3) = dcmplx(0.0d0,0.0d0) 
                end if
            end if       
            
            call dfftw_execute(plan)    ! take fourier transform
            vyz(i,1:nlong) = grid(1:nlong) 
        
        end if
                
        if (i /= 1) then    ! don't compute value for south pole.
            ! Vzz = Vrr
            coef(1) = dcmplx(coefrrs0,0.0d0)
            coef(2:lmax+1) = coefrrs(2:lmax+1) / 2.0d0
            
            if (present(sampling)) then
                if (sampling == 2) then
                    coef(lmax+2:2*lmax+3) = dcmplx(0.0d0,0.0d0) 
                end if
            end if       
            
            call dfftw_execute(plan)    ! take fourier transform
            vzz(i_s,1:nlong) = grid(1:nlong) 

            ! Vxx = 1/r Vr + 1/r^2 Vtt
            coef(1) = dcmplx(coefrs0/r_ex + coeftts0/r_ex**2,0.0d0)
            coef(2:lmax+1) = (coefrs(2:lmax+1)/r_ex &
                                + coeftts(2:lmax+1)/r_ex**2 ) / 2.0d0
                                
            if (present(sampling)) then
                if (sampling == 2) then
                    coef(lmax+2:2*lmax+3) = dcmplx(0.0d0,0.0d0) 
                end if
            end if       
            
            call dfftw_execute(plan)    ! take fourier transform
            vxx(i_s,1:nlong) = grid(1:nlong)
                    
            ! Vyy = 1/r Vr + 1/r^2 /tan(t) Vt + 1/r^2 /sin(t)^2 Vpp
            coef(1) = dcmplx(coefrs0/r_ex + coefts0/(r_ex**2)/tan(theta) &
                            + coefpps0/(r_ex**2)/u**2 ,0.0d0)
            coef(2:lmax+1) = (coefrs(2:lmax+1)/r_ex &
                            + coefts(2:lmax+1)/(r_ex**2)/tan(theta) + &
                            coefpps(2:lmax+1)/(r_ex**2)/u**2 ) / 2.0d0
                            
            if (present(sampling)) then
                if (sampling == 2) then
                    coef(lmax+2:2*lmax+3) = dcmplx(0.0d0,0.0d0) 
                end if
            end if       
            
            call dfftw_execute(plan)    ! take fourier transform
            vyy(i_s,1:nlong) = grid(1:nlong)
                    
            ! Vxy = 1/r^2 /sin(t) Vtp - cos(t)/sin(t)^2 /r^2 Vp
            coef(1) = dcmplx(coeftps0/sint/r_ex**2 &
                            - coefps0/(r_ex**2)*z/u**2, 0.0d0)
            coef(2:lmax+1) = (coeftps(2:lmax+1)/sint/r_ex**2 &
                            - coefps(2:lmax+1)/(r_ex**2)*z/u**2 ) / 2.0d0
                            
            if (present(sampling)) then
                if (sampling == 2) then
                    coef(lmax+2:2*lmax+3) = dcmplx(0.0d0,0.0d0) 
                end if
            end if       
            
            call dfftw_execute(plan)    ! take fourier transform
            vxy(i_s,1:nlong) = grid(1:nlong)
                    
            ! Vxz = 1/r^2 Vt - 1/r Vrt
            coef(1) = dcmplx(coefts0/r_ex**2 - coefrts0/r_ex, 0.0d0)
            coef(2:lmax+1) = (coefts(2:lmax+1)/r_ex**2 &
                            - coefrts(2:lmax+1)/r_ex ) / 2.0d0
                            
            if (present(sampling)) then
                if (sampling == 2) then
                    coef(lmax+2:2*lmax+3) = dcmplx(0.0d0,0.0d0) 
                end if
            end if       
            
            call dfftw_execute(plan)    ! take fourier transform
            vxz(i_s,1:nlong) = grid(1:nlong)
                    
            ! Vyz = 1/r^2 /sin(t) Vp - 1/r /sin(t) Vrp
            coef(1) = dcmplx(coefps0/sint/r_ex**2 - coefrps0/sint/r_ex, 0.0d0)
            coef(2:lmax+1) = (coefps(2:lmax+1)/sint/r_ex**2 &
                            - coefrps(2:lmax+1)/sint/r_ex ) / 2.0d0
                            
            if (present(sampling)) then
                if (sampling == 2) then
                    coef(lmax+2:2*lmax+3) = dcmplx(0.0d0,0.0d0) 
                end if
            end if       
            
            call dfftw_execute(plan)    ! take fourier transform
            vyz(i_s,1:nlong) = grid(1:nlong) 
        
        end if
            
    end do
    
    ! Finally, do equator
    
    r_ex = a
    theta = pi/2.0d0
    z = 0.0d0
    u = 1.0d0
    lat = 0.0d0

    coefr(1:lmax+1) = dcmplx(0.0d0,0.0d0)
    coefr0 = 0.0d0
        
    coefrr(1:lmax+1) = dcmplx(0.0d0,0.0d0)
    coefrr0 = 0.0d0
        
    coeft(1:lmax+1) = dcmplx(0.0d0,0.0d0)
    coeft0 = 0.0d0
        
    coeftp(1:lmax+1) = dcmplx(0.0d0,0.0d0)
    coeftp0 = 0.0d0
        
    coefp(1:lmax+1) = dcmplx(0.0d0,0.0d0)
    coefp0 = 0.0d0
        
    coefrt(1:lmax+1) = dcmplx(0.0d0,0.0d0)
    coefrt0 = 0.0d0
        
    coefrp(1:lmax+1) = dcmplx(0.0d0,0.0d0)
    coefrp0 = 0.0d0

    coefpp(1:lmax+1) = dcmplx(0.0d0,0.0d0)
    coefpp0 = 0.0d0
        
    coeftt(1:lmax+1) = dcmplx(0.0d0,0.0d0)
    coeftt0 = 0.0d0
                
    pm2 = 1.0d0

    tempr =  -cilm(1,1,1) * pm2 ! l = 0 
    coefr0 = coefr0 + tempr
        
    tempr =  2 * cilm(1,1,1) * pm2  ! l = 0 
    coefrr0 = coefrr0 + tempr
        
    ! derivative in theta and phi of l=0 term is 0, so no need to calculate this
                
    if (lmax_comp /= 0) then    ! l = 1
        prefactor(1) = r0 / r_ex
        
        do l = 2, lmax_comp, 1
            prefactor(l) = prefactor(l-1) * r0 / r_ex
        end do

        pm1 = 0.0d0
            
        dp = ff1(2,1)
        tempr = cilm(1,2,1) * dp * prefactor(1)
        coeft0 = coeft0 + tempr
            
        tempr = cilm(1,2,1) * dp * (-2) * prefactor(1)  ! -2 = (l+1) prefactor
        coefrt0 = coefrt0 + tempr
            
    end if
                
    do l = 2, lmax_comp, 1
        l1 = l+1
        p =  - ff2(l1,1) * pm2
        tempr = cilm(1,l1,1) * p * (-l1) * prefactor(l)
        coefr0 = coefr0 + tempr
            
        tempr = cilm(1,l1,1) * p * (l1) * (l1+1) * prefactor(l)
        coefrr0 = coefrr0 + tempr
            
        dp = l * ( sqr(2*l+1) / sqr(2*l-1) * pm1 )
        tempr = cilm(1,l1,1) * dp * prefactor(l)
        coeft0 = coeft0 + tempr
            
        tempr = cilm(1,l1,1) * dp * (-l1) * prefactor(l)
        coefrt0 = coefrt0 + tempr
            
        dp2 = -l*l1 * p 
        tempr = cilm(1,l1,1) * dp2 * prefactor(l)
        coeftt0 = coeftt0 + tempr
            
        pm2 = pm1
        pm1 = p
        
    end do
                
    pmm = sqr(2) * scalef
                
    rescalem = 1.0d0/scalef
            
    do m = 1, lmax_comp-1, 1
            
        m1 = m+1
                    
        pmm = pmm * sqr(2*m+1) / sqr(2*m)
        pm2 = pmm
            
        tempc = dcmplx(cilm(1,m1,m1), - cilm(2,m1,m1)) * pm2 * (-m-1) &
                        * prefactor(m)    ! (m,m)
        coefr(m1) = coefr(m1) + tempc
            
        tempc = dcmplx(cilm(1,m1,m1), - cilm(2,m1,m1)) * pm2 * (m+1) * (m+2) &
                        * prefactor(m) ! (m,m)
        coefrr(m1) = coefrr(m1) + tempc
            
        tempc = dcmplx(cilm(2,m1,m1), cilm(1,m1,m1)) * pm2 &
                        * prefactor(m) * m ! (m,m)
        coefp(m1) = coefp(m1) + tempc
            
        tempc = -dcmplx(cilm(1,m1,m1), - cilm(2,m1,m1)) * pm2 &
                        * prefactor(m) * m**2 ! (m,m)
        coefpp(m1) = coefpp(m1) + tempc
            
        tempc = dcmplx(cilm(2,m1,m1), cilm(1,m1,m1)) * pm2 *(-m-1) &
                        * prefactor(m) * m ! (m,m)
        coefrp(m1) = coefrp(m1) + tempc

        dp2 = -(m*m1 -(m**2)) * pm2 
        tempc = dcmplx(cilm(1,m1,m1), - cilm(2,m1,m1)) * dp2 &
                        * prefactor(m) ! (m,m)
        coeftt(m1) = coeftt(m1) + tempc 
        
        pm1 = 0.0d0
        
        dp =  ( sqr(2*m+3) * pmm ) 
        tempc = dcmplx(cilm(1,m1+1,m1), - cilm(2,m1+1,m1)) * dp &
                        * prefactor(m+1)    ! (m+1,m)
        coeft(m1) = coeft(m1) + tempc   
            
        tempc = dcmplx(cilm(1,m1+1,m1), - cilm(2,m1+1,m1)) * dp * (-m-2) &
                        * prefactor(m+1)   ! (m+1,m)
        coefrt(m1) = coefrt(m1) + tempc 
            
        tempc = dcmplx(cilm(2,m1+1,m1), cilm(1,m1+1,m1)) * dp &
                        * prefactor(m+1) * m  ! (m+1,m)
        coeftp(m1) = coeftp(m1) + tempc 
                    
        do l = m + 2, lmax_comp, 1
            l1 = l+1
            p =  - ff2(l1,m1) * pm2
            tempc = dcmplx(cilm(1,l1,m1), - cilm(2,l1,m1)) * p * (-l1) &
                            * prefactor(l)
            coefr(m1) = coefr(m1) + tempc
                
            tempc = dcmplx(cilm(1,l1,m1), - cilm(2,l1,m1)) * p * (l1) &
                            * (l1+1) * prefactor(l)
            coefrr(m1) = coefrr(m1) + tempc
                
            tempc = dcmplx(cilm(2,l1,m1), cilm(1,l1,m1)) * p &
                            * prefactor(l) * m
            coefp(m1) = coefp(m1) + tempc
                
            tempc = -dcmplx(cilm(1,l1,m1), - cilm(2,l1,m1)) * p &
                            * prefactor(l) * m**2
            coefpp(m1) = coefpp(m1) + tempc
                
            tempc = dcmplx(cilm(2,l1,m1), cilm(1,l1,m1)) * p * (-l1) &
                            * prefactor(l) * m
            coefrp(m1) = coefrp(m1) + tempc
                
            dp = ( sqr(2*l+1) * sqr(l-m) * sqr(l+m) / sqr(2*l-1) * pm1) 
            tempc = dcmplx(cilm(1,l1,m1), - cilm(2,l1,m1)) * dp * prefactor(l)
            coeft(m1) = coeft(m1) + tempc
                
            tempc = dcmplx(cilm(1,l1,m1), - cilm(2,l1,m1)) * dp * (-l1) &
                            * prefactor(l)
            coefrt(m1) = coefrt(m1) + tempc
                
            tempc = dcmplx(cilm(2,l1,m1), cilm(1,l1,m1)) * dp * prefactor(l) * m
            coeftp(m1) = coeftp(m1) + tempc
                
            dp2 = -(l*l1 -(m**2)/u**2) * p 
            tempc = dcmplx(cilm(1,l1,m1), - cilm(2,l1,m1)) * dp2 * prefactor(l)
            coeftt(m1) = coeftt(m1) + tempc

            pm2 = pm1
            pm1 = p
                    
        end do
                    
        coefr(m1) = coefr(m1) * rescalem                  
        coefrr(m1) = coefrr(m1) * rescalem        
        coeft(m1) = coeft(m1) * rescalem   
        coeftt(m1) = coeftt(m1) * rescalem
        coefrt(m1) = coefrt(m1) * rescalem
        coefp(m1) = coefp(m1) * rescalem
        coefpp(m1) = coefpp(m1) * rescalem
        coeftp(m1) = coeftp(m1) * rescalem
        coefrp(m1) = coefrp(m1) * rescalem
                    
    end do           
                                
    if (lmax_comp /= 0) then

        pmm = pmm * sqr(2*lmax_comp+1) / sqr(2*lmax_comp) * rescalem                    
        tempc = dcmplx(cilm(1,lmax_comp+1,lmax_comp+1), &
                        - cilm(2,lmax_comp+1,lmax_comp+1)) * pmm &
                        * (-lmax_comp-1) * prefactor(lmax_comp)
        coefr(lmax_comp+1) = coefr(lmax_comp+1) + tempc
            
        tempc = dcmplx(cilm(1,lmax_comp+1,lmax_comp+1), &
                        - cilm(2,lmax_comp+1,lmax_comp+1)) * pmm &
                        * (lmax_comp+1) * (lmax_comp+2) * prefactor(lmax_comp)
        coefrr(lmax_comp+1) = coefrr(lmax_comp+1) + tempc
        
        tempc = dcmplx(cilm(2,lmax_comp+1,lmax_comp+1), &
                        cilm(1,lmax_comp+1,lmax_comp+1)) * pmm &
                        * prefactor(lmax_comp) * lmax_comp
        coefp(lmax_comp+1) = coefp(lmax_comp+1) + tempc
            
        tempc = -dcmplx(cilm(1,lmax_comp+1,lmax_comp+1), &
                        - cilm(2,lmax_comp+1,lmax_comp+1)) * pmm &
                        * prefactor(lmax_comp) * lmax_comp**2
        coefpp(lmax_comp+1) = coefpp(lmax_comp+1) + tempc
            
        tempc = dcmplx(cilm(2,lmax_comp+1,lmax_comp+1), &
                        cilm(1,lmax_comp+1,lmax_comp+1)) * pmm &
                        * (-lmax_comp-1) * prefactor(lmax_comp) * lmax_comp
        coefrp(lmax_comp+1) = coefrp(lmax_comp+1) + tempc
            
        dp2 = -(lmax_comp*(lmax_comp+1)-(lmax_comp**2)) * pmm 
        tempc = dcmplx(cilm(1,lmax_comp+1,lmax_comp+1), &
                        - cilm(2,lmax_comp+1,lmax_comp+1)) * dp2 &
                        * prefactor(lmax_comp)
        coeftt(lmax_comp+1) = coeftt(lmax_comp+1) + tempc
        
    end if
        
    coefr0 = coefr0 * gm / r_ex**2
    coefr(2:lmax+1) = coefr(2:lmax+1) * gm / r_ex**2
            
    coefrt0 = - coefrt0 * gm / r_ex**2
    coefrt(2:lmax+1) = - coefrt(2:lmax+1) * gm / r_ex**2
        
    coefrr0 = coefrr0 * gm / r_ex**3
    coefrr(2:lmax+1) = coefrr(2:lmax+1) * gm / r_ex**3
        
    coeft0 = - coeft0 * gm / r_ex 
    coeft(2:lmax+1) = - coeft(2:lmax+1) * gm / r_ex 
        
    coeftt0 = coeftt0 * gm / r_ex 
    coeftt(2:lmax+1) = coeftt(2:lmax+1) * gm / r_ex 
        
    coeftp0 = - coeftp0 * gm / r_ex 
    coeftp(2:lmax+1) = - coeftp(2:lmax+1) * gm / r_ex 
        
    coefp0 = coefp0 * gm / r_ex
    coefp(2:lmax+1) = coefp(2:lmax+1) * gm / r_ex
        
    coefpp0 = coefpp0 * gm / r_ex
    coefpp(2:lmax+1) = coefpp(2:lmax+1) * gm / r_ex
        
    coefrp0 = coefrp0 * gm / r_ex**2
    coefrp(2:lmax+1) = coefrp(2:lmax+1) * gm / r_ex**2
        
    ! Vzz = Vrr
    coef(1) = dcmplx(coefrr0,0.0d0)
    coef(2:lmax+1) = coefrr(2:lmax+1) / 2.0d0
    
    if (present(sampling)) then
        if (sampling == 2) then
            coef(lmax+2:2*lmax+3) = dcmplx(0.0d0,0.0d0) 
        end if
    end if 
          
    call dfftw_execute(plan)    ! take fourier transform
    vzz(i_eq,1:nlong) = grid(1:nlong)
                                
    ! Vxx = 1/r Vr + 1/r^2 Vtt
    coef(1) = dcmplx(coefr0/r_ex + coeftt0/r_ex**2,0.0d0)
    coef(2:lmax+1) = (coefr(2:lmax+1)/r_ex + coeftt(2:lmax+1)/r_ex**2 ) / 2.0d0
    
    if (present(sampling)) then
        if (sampling == 2) then
            coef(lmax+2:2*lmax+3) = dcmplx(0.0d0,0.0d0) 
        end if
    end if
           
    call dfftw_execute(plan)    ! take fourier transform
    vxx(i_eq,1:nlong) = grid(1:nlong)
                    
    ! Vyy = 1/r Vr + 1/r^2 /tan(t) Vt + 1/r^2 /sin(t)^2 Vpp
    coef(1) = dcmplx(coefr0/r_ex  + coefpp0/(r_ex**2) ,0.0d0)
    coef(2:lmax+1) = (coefr(2:lmax+1)/r_ex &
                    + coefpp(2:lmax+1)/(r_ex**2) ) / 2.0d0
                    
    if (present(sampling)) then
        if (sampling == 2) then
            coef(lmax+2:2*lmax+3) = dcmplx(0.0d0,0.0d0) 
        end if
    end if 
          
    call dfftw_execute(plan)    ! take fourier transform
    vyy(i_eq,1:nlong) = grid(1:nlong)
                    
    ! Vxy = 1/r^2 /sin(t) Vtp - cos(t)/sin(t)^2 /r^2 Vp
    coef(1) = dcmplx(coeftp0/r_ex**2 , 0.0d0)
    coef(2:lmax+1) = (coeftp(2:lmax+1)/r_ex**2 ) / 2.0d0
    
    if (present(sampling)) then
        if (sampling == 2) then
            coef(lmax+2:2*lmax+3) = dcmplx(0.0d0,0.0d0) 
        end if
    end if       
    
    call dfftw_execute(plan)    ! take fourier transform
    vxy(i_eq,1:nlong) = grid(1:nlong)
                    
    ! Vxz = 1/r^2 Vt - 1/r Vrt
    coef(1) = dcmplx(coeft0/r_ex**2 - coefrt0/r_ex, 0.0d0)
    coef(2:lmax+1) = (coeft(2:lmax+1) / r_ex**2 &
                    - coefrt(2:lmax+1)/r_ex ) / 2.0d0
                    
    if (present(sampling)) then
        if (sampling == 2) then
            coef(lmax+2:2*lmax+3) = dcmplx(0.0d0,0.0d0) 
        endif
    end if
          
    call dfftw_execute(plan)    ! take fourier transform
    vxz(i_eq,1:nlong) = grid(1:nlong)
                    
    ! Vyz = 1/r^2 /sin(t) Vp - 1/r /sin(t) Vrp
    coef(1) = dcmplx(coefp0/r_ex**2 - coefrp0/r_ex, 0.0d0)
    coef(2:lmax+1) = (coefp(2:lmax+1) / r_ex**2 &
                    - coefrp(2:lmax+1)/r_ex ) / 2.0d0
                    
    if (present(sampling)) then
        if (sampling == 2) then
            coef(lmax+2:2*lmax+3) = dcmplx(0.0d0,0.0d0) 
        end if
    end if
        
    call dfftw_execute(plan)    ! take fourier transform
    vyz(i_eq,1:nlong) = grid(1:nlong)
                    
    call dfftw_destroy_plan(plan)
                
end subroutine MakeGravGradGridDH
