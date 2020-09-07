subroutine SHLocalizedAdmitCorr(tapers, taper_order, lwin, lat, lon, gilm, &
                                tilm, lmax, admit, corr, K, admit_error, &
                                corr_error, taper_wt, mtdef, k1linsig, &
                                exitstatus)
!------------------------------------------------------------------------------
!
!   Given two spherical harmonic fields (gilm and tilm), this routine will
!   calculate the localized admittance and correlation using the first
!   space-concentrated window of Wieczorek and Simons (2005). All functions
!   must be 4-pi normalized, and exclude the Condon-Shortley phase factor. Two
!   manners of calculating the localized admittance and correlation are
!   possible according to the optional parameter MTDEF. In one case, the
!   multitaper cross-power spectra are calculated, and from these, the
!   admittance and correlation. In the second, the admittance and correlation
!   are calculated for each taper, and these are then averaged.
!
!   Calling Parameters
!
!       IN
!
!           tapers          A matrix of tapers obtained from SHReturnTapers.
!           taper_order     A vector continaing the angular order of each column
!                           of taper.
!           lwin            Spectral bandwidth of the localizing window.
!           lat, lon        Latitude and longitude that the window will be
!                           rotated to, in DEGREES.
!           gilm, tilm      Input spherical harmonic fields.
!           K               Number of tapers to use in Multitaper spectral
!                           estimations.
!           lmax            Maximum spherical harmonic degree of the intput
!                           fields.
!
!       OUT
!           admit           Admittance between the localized gilm and tilm
!                           assuming that gilm = Z tilm.
!           corr            Correlation of the two fields.
!
!       OPTIONAL (OUT)
!           admit_error     Error of the admittance (only when K>1)
!           corr_error      Error of the admittance (only when K>1)
!
!       OPTIONAL (IN)
!           mtdef           1 (default): Calculate multitaper cross-spectral
!                           estimates, and use these to calculate a single
!                           admittance and correlations.
!                           2: Calculate the admittance and correlation using
!                           each individual taper, and then average these to
!                           get the admittance and correlation.
!           taper_wt        Weights to be applied to the spectral estimates.
!                           This can only be used when MTDEF is 1.
!           k1linsig:       If present and equal to 1, the uncertainty in the
!                           admittance will be calculated by assuming the
!                           gravity and topography coefficients are linearly
!                           correlated, and that any lack of correlation is the
!                           result of uncorrelated noise. This should only be
!                           used when one expects the gravity and topography
!                           to be linearly correlated and when only a single
!                           taper is being used. This should not be used with a
!                           Forsyth type model that predicts a less than 1
!                           correlation coefficient. This is the square root
!                           of eq. 33 of Simons et al. 1997.
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
!   Notes:
!       1. The units of the output admittance will correspond to the units of
!       the spherical harmonic coefficients. If gravity/topography admittances
!       are desired, then either the gravity coefficients should be multiplied
!       by G M (l+1) * (r0/r)**(l+2) / r0**2 before calling this routine, or
!       the admittances should be multiplied by this factor afterwards.
!
!       2. The correlation is defined as Sgt / sqrt(Sgg Stt), which varies
!       between -1 and 1. To obtain the "coherence" (which is also sometimes
!       referred to as the "coherence squared"), just square this number.
!
!   Copyright (c) 2005-2019, SHTOOLS
!   All rights reserved.
!
!------------------------------------------------------------------------------
    use SHTOOLS, only:  djpi2, SHRotateRealCoef, &
                        SHCrossPowerSpectrum, SHPowerSpectrum, MakeGridGLQ, &
                        SHGLQ, SHExpandGLQ
    use ftypes

    implicit none

    real(dp), intent(in) :: tapers(:,:), lat, lon, gilm(:,:,:), tilm(:,:,:)
    integer, intent(in) :: lwin, lmax, K, taper_order(:)
    real(dp), intent(out) :: admit(:), corr(:)
    real(dp), intent(out), optional :: admit_error(:), corr_error(:)
    integer, intent(in), optional :: mtdef, k1linsig
    real(dp), intent(in), optional :: taper_wt(:)
    integer, intent(out), optional :: exitstatus
    integer :: lmaxwin, l, def, astat(9), phase, norm, i, nlat, nlong
    integer, save :: first = 1, lmaxwin_last = -1, lwin_last = -1
    real(dp) :: pi, g_power(2,lwin+lmax+1), t_power(2,lwin+lmax+1), &
                gt_power(2,lwin+lmax+1), x(3), sgt(lmax-lwin+1, K), &
                sgg(lmax-lwin+1, K), stt(lmax-lwin+1, K), &
                admit_k(lmax-lwin+1, K), corr_k(lmax-lwin+1, K), factor
    real(dp), allocatable :: shwin(:,:,:), shwinrot(:,:,:), shloc_g(:,:,:), &
                             shloc_t(:,:,:), gridtglq(:,:), gridgglq(:,:), &
                             gridwinglq(:,:), temp(:,:)
    real(dp), allocatable, save :: dj(:,:,:), zero(:), w(:)

!$OMP   threadprivate(first, lmaxwin_last, lwin_last, dj)

    if (present(exitstatus)) exitstatus = 0

    phase = 1
    norm = 1

    pi = acos(-1.0_dp)
    lmaxwin = lmax + lwin

    if (present(k1linsig) .and. K /= 1) then
        if (k1linsig == 1) then
            print*, "Error --- SHlocalizedAdmitCorr"
            print*, "If K1LINSIG is present and equal to 1, K must be " // &
                    "equal to 1."
            print*, "Input value of K1LINSIG is ", k1linsig
            if (present(exitstatus)) then
                exitstatus = 2
                return
            else
                stop
            end if

        end if

    end if

    if (size(admit) < lmax-lwin+1) then
        print*, "Error --- SHLocalizedAdmitCorr"
        print*, "ADMIT must be dimensioned as (LMAX-LWIN+1) where " // &
                "LMAX and LWIN are ", lmax, lwin
        print*, "Input array is dimensioned ", size(admit)
        if (present(exitstatus)) then
            exitstatus = 1
            return
        else
            stop
        end if

    else if (size(corr) < lmax-lwin+1) then
        print*, "Error --- SHLocalizedAdmitCorr"
        print*, "CORR must be dimensioned as (LMAX-LWIN+1) where " // &
                "LMAX and LWIN are ", lmax, lwin
        print*, "Input array is dimensioned ", size(corr)
        if (present(exitstatus)) then
            exitstatus = 1
            return
        else
            stop
        end if

    else if (size(tapers(:,1)) < lwin+1 .or. size(tapers(1,:)) < K) then
        print*, "Error --- SHLocalizedAdmitCorr"
        print*, "TAPERS must be dimensioned as (LWIN+1, K) where " // &
                "LWIN and K are ", lwin, K
        print*, "Iinput array is dimensioned as ", size(tapers(:,1)), &
                size(tapers(1,:))
        if (present(exitstatus)) then
            exitstatus = 1
            return
        else
            stop
        end if

    else if (size(taper_order) < K ) then
        print*, "Error --- SHLocalizedAdmitCorr"
        print*, "TAPER_ORDER must be dimensioned as (K) where K is ", K
        print*, "Input array is dimensioned ", size(taper_order)
        if (present(exitstatus)) then
            exitstatus = 1
            return
        else
            stop
        end if

    else if (size(gilm(:,1,1)) < 2 .or. size(gilm(1,:,1)) < lmax+1 .or. &
             size(gilm(1,1,:)) < lmax+1) then
        print*, "Error --- SHLocalizedAdmitCorr"
        print*, "GILM must be dimensioned as (2, LMAX+1, LMAX+1) " // &
                "where LMAX is ", lmax
        print*, "Input array is dimensioned ", size(gilm(:,1,1)), &
                size(gilm(1,:,1)), size(gilm(1,1,:))
        if (present(exitstatus)) then
            exitstatus = 1
            return
        else
            stop
        end if

    else if (size(tilm(:,1,1)) < 2 .or. size(tilm(1,:,1)) < lmax+1 .or. &
             size(tilm(1,1,:)) < lmax+1) then
        print*, "Error --- SHLocalizedAdmitCorr"
        print*, "TILM must be dimensioned as (2, LMAX+1, LMAX+1) " // & 
                "where LMAX is ", lmax
        print*, "Input array is dimensioned ", size(tilm(:,1,1)), &
                size(tilm(1,:,1)), size(tilm(1,1,:))
        if (present(exitstatus)) then
            exitstatus = 1
            return
        else
            stop
        end if

    end if

    if (present(admit_error)) then
        if (size(admit_error) < lmax - lwin+1) then
            print*, "Error ---SHLocalizedAdmitCorr"
            print*, "ADMIT_ERROR must be dimensioned as (LMAX-LWIN+1) " // &
                    "where LMAX and LMAXT are ", lmax, lwin
            print*, "Input array is dimensioned ", size(admit_error)
            if (present(exitstatus)) then
                exitstatus = 1
                return
            else
                stop
            end if

        end if

    end if

    if (present(corr_error)) then
        if (size(corr_error) < lmax - lwin + 1) then
            print*, "Error ---SHLocalizedAdmitCorr"
            print*, "CORR_ERROR  must be dimensioned as (LMAX-LWIN+1) " // &
                    "where LMAX and LMAXT are ", lmax, lwin
            print*, "Input array is dimensioned ", size(corr_error)
            if (present(exitstatus)) then
                exitstatus = 1
                return
            else
                stop
            end if

        end if

    end if

    if(present(taper_wt)) then
        if (size(taper_wt) < K) then
            print*, "Error --- SHLocalizedAdmitCorr"
            print*, "TAPER_WT must be dimensioned as (K) where K is ", K
            print*, "Input array has dimension ", size(taper_wt)
            if (present(exitstatus)) then
                exitstatus = 1
                return
            else
                stop
            end if

        end if

    end if

    if (present(mtdef)) then
        if (mtdef == 2 .and. present(taper_wt)) then
            print*, "Error --- SHLocalizedAdmitCorr"
            print*, "TAPER_WT can only be used when MTDEF is 1."
            if (present(exitstatus)) then
                exitstatus = 2
                return
            else
                stop
            end if

        end if

    end if

    if (present(mtdef)) then
        if (mtdef /= 1 .and. mtdef /= 2) then
            print*, "SHLocalizedAdmitCorr --- Error"
            print*, "MTDEF must be 1 or 2."
            print*, "Input value is ", mtdef
            if (present(exitstatus)) then
                exitstatus = 2
                return
            else
                stop
            end if

        else
            def = mtdef

        end if

    else
        def = 1

    end if


    admit = 0.0_dp
    corr = 0.0_dp

    if (present(admit_error)) then
        admit_error = 0.0_dp
    end if

    if (present(corr_error)) then
        corr_error = 0.0_dp
    end if

    !--------------------------------------------------------------------------
    !
    !   Determine multitaper cross spectra estimates of gilm and tilm, and then
    !   calculate the admittance and correlation. Errors for the latter are
    !   calculated by adding the error sources in quadrature. Taper weights can
    !   be specified in order to minimize the variance of the multitaper
    !   cross-spectral estimates.
    !
    !   Note that only the admittances with degrees greater than lwin and less
    !   than Lmax - lwin should be interpretted.
    !
    !--------------------------------------------------------------------------
    x(1) = 0.0_dp
    x(2) = -(90.0_dp - lat) * pi / 180.0_dp
    x(3) = -lon * pi / 180.0_dp
        
    if (first == 1) then
        lwin_last = lwin
        lmaxwin_last = lmaxwin
        first = 0
 
        allocate (zero(lmaxwin+1), stat = astat(1))
        allocate (w(lmaxwin+1), stat = astat(2))
        allocate (dj(lwin+1,lwin+1,lwin+1), stat = astat(3))
        
        if (sum(astat(1:3)) /= 0) then
            print*, "Error --- SHLocalizedAdmitCorr"
            print*, "Problem allocating arrays ZERO, W and DJ", &
                    astat(1), astat(2), astat(3)
            if (present(exitstatus)) then
                exitstatus = 3
                return
            else
                stop
            end if

        end if

        if (present(exitstatus)) then
            call SHGLQ(lmaxwin, zero, w, csphase = phase, norm = 1, &
                       exitstatus = exitstatus)
            if (exitstatus /= 0) return
            dj = 0.0_dp
            call djpi2(dj, lwin, exitstatus = exitstatus)
            if (exitstatus /= 0) return
        else
            call SHGLQ(lmaxwin, zero, w, csphase = phase, norm = 1)
            dj = 0.0_dp
            call djpi2(dj, lwin)
        end if

    end if

    if (lwin > lwin_last) then
        lwin_last = lwin
        deallocate (dj)

        allocate (dj(lwin+1,lwin+1,lwin+1), stat = astat(1))

        if (astat(1) /= 0) then
            print*, "Error --- SHLocalizedAdmitCorr"
            print*, "Problem allocating array DJ", astat(1)
            if (present(exitstatus)) then
                exitstatus = 3
                return
            else
                stop
            end if

        end if

        dj = 0.0_dp

        if (present(exitstatus)) then
            call djpi2(dj, lwin, exitstatus = exitstatus)
            if (exitstatus /= 0) return
        else
            call djpi2(dj, lwin)
        end if

    end if

    if (lmaxwin /= lmaxwin_last) then
        lmaxwin_last = lmaxwin

        deallocate (zero)
        deallocate (w)
        allocate (zero(lmaxwin+1), stat = astat(1))
        allocate (w(lmaxwin+1), stat = astat(2))
        
        if (sum(astat(1:2)) /= 0) then
            print*, "Error --- SHLocalizedAdmitCorr"
            print*, "Problem allocating arrays ZERO and W", astat(1), astat(2)
            if (present(exitstatus)) then
                exitstatus = 3
                return
            else
                stop
            end if

        end if

        if (present(exitstatus)) then
            call SHGLQ(lmaxwin, zero, w, csphase = phase, norm = 1, &
                       exitstatus = exitstatus)
            if (exitstatus /= 0) return
        else
            call SHGLQ(lmaxwin, zero, w, csphase = phase, norm = 1)
        end if

    end if

    nlat = lmax + lwin + 1
    nlong = 2 * (lmax + lwin) + 1

    allocate (shwin(2,lwin+1,lwin+1), stat = astat(1))
    allocate (shwinrot(2,lwin+1,lwin+1), stat = astat(2))
    allocate (shloc_g(2, lmaxwin+1, lmaxwin+1), stat= astat(3))
    allocate (shloc_t(2, lmaxwin+1, lmaxwin+1), stat= astat(4))
    allocate (gridtglq(nlat,nlong), stat = astat(5))
    allocate (gridgglq(nlat,nlong), stat = astat(6))
    allocate (gridwinglq(nlat,nlong), stat = astat(7))
    allocate (temp(nlat,nlong), stat = astat(8))

    if (sum(astat(1:8)) /= 0) then
        print*, "Error --- SHLocalizedAdmitCorr"
        print*, "Problem allocating arrays SHWIN, SHWINROT, SHLOC_G, " // &
            "SHLOC_T, GRIDTGLQ, GRIDGGLQ, GRIDWINGLQ, and TEMP", &
            astat(1), astat(2), astat(3), astat(4), astat(5), astat(6), &
            astat(7), astat(8)
            if (present(exitstatus)) then
                exitstatus = 3
                return
            else
                stop
            end if

    end if

    if (present(exitstatus)) then
        call MakeGridGLQ(gridtglq, tilm(1:2,1:lmax+1, 1:lmax+1), &
                         lmaxwin, zero = zero, csphase = phase, norm = 1, &
                         exitstatus = exitstatus)
        if (exitstatus /= 0) return
        call MakeGridGLQ(gridgglq, gilm(1:2,1:lmax+1, 1:lmax+1), &
                         lmaxwin, zero = zero, csphase = phase, norm = 1, &
                         exitstatus = exitstatus)
        if (exitstatus /= 0) return
    else
        call MakeGridGLQ(gridtglq, tilm(1:2,1:lmax+1, 1:lmax+1), &
                         lmaxwin, zero = zero, csphase = phase, norm = 1)
        call MakeGridGLQ(gridgglq, gilm(1:2,1:lmax+1, 1:lmax+1), &
                         lmaxwin, zero = zero, csphase = phase, norm = 1)
    end if

    if (def == 1) then
        do i = 1, K
            shwin = 0.0_dp

            if (taper_order(i) < 0) then
                shwin(2,1:lwin+1,abs(taper_order(i))+1) = tapers(1:lwin+1,i)

            else
                shwin(1,1:lwin+1,taper_order(i)+1) = tapers(1:lwin+1,i)

            end if

            if (present(exitstatus)) then
                call SHRotateRealCoef(shwinrot, shwin, lwin, x, dj, &
                                      exitstatus = exitstatus)
                if (exitstatus /= 0) return
                call MakeGridGLQ(gridwinglq, shwinrot(1:2,1:lwin+1, 1:lwin+1),&
                                 lmaxwin, zero = zero, csphase = phase, &
                                 norm = 1, exitstatus = exitstatus)
                if (exitstatus /= 0) return
                temp(1:nlat,1:nlong) = gridtglq(1:nlat,1:nlong) &
                                       * gridwinglq(1:nlat,1:nlong)
                call SHExpandGLQ(shloc_t, lmaxwin, temp, w, zero = zero, &
                                 csphase = phase, norm = 1, &
                                 exitstatus = exitstatus)
                if (exitstatus /= 0) return
                temp(1:nlat,1:nlong) = gridgglq(1:nlat,1:nlong) &
                                       * gridwinglq(1:nlat,1:nlong)
                call SHExpandGLQ(shloc_g, lmaxwin, temp, w, zero = zero, &
                                csphase = phase, norm = 1, &
                                exitstatus = exitstatus)
                if (exitstatus /= 0) return
                call SHCrossPowerSpectrum(shloc_g, shloc_t, lmax-lwin, &
                                          sgt(:,i), exitstatus = exitstatus)
                if (exitstatus /= 0) return
                call SHPowerSpectrum(shloc_g, lmax-lwin, sgg(:,i), &
                                     exitstatus = exitstatus)
                if (exitstatus /= 0) return
                call SHPowerSpectrum(shloc_t, lmax-lwin, stt(:,i), &
                                     exitstatus = exitstatus)
                if (exitstatus /= 0) return
            else
                call SHRotateRealCoef(shwinrot, shwin, lwin, x, dj)
                call MakeGridGLQ(gridwinglq, shwinrot(1:2,1:lwin+1, 1:lwin+1),&
                                 lmaxwin, zero = zero, csphase = phase, &
                                 norm = 1)
                temp(1:nlat,1:nlong) = gridtglq(1:nlat,1:nlong) &
                                       * gridwinglq(1:nlat,1:nlong)
                call SHExpandGLQ(shloc_t, lmaxwin, temp, w, zero = zero, &
                                 csphase = phase, norm = 1)
                temp(1:nlat,1:nlong) = gridgglq(1:nlat,1:nlong) &
                                       * gridwinglq(1:nlat,1:nlong)
                call SHExpandGLQ(shloc_g, lmaxwin, temp, w, zero = zero, &
                                csphase = phase, norm = 1)
                call SHCrossPowerSpectrum(shloc_g, shloc_t, lmax-lwin, &
                                          sgt(:,i))
                call SHPowerSpectrum(shloc_g, lmax-lwin, sgg(:,i))
                call SHPowerSpectrum(shloc_t, lmax-lwin, stt(:,i))
        end if

            if (present(taper_wt)) then
                factor = sum(taper_wt(1:K))**2 - sum(taper_wt(1:K)**2)
                factor = factor * sum(taper_wt(1:K))
                factor = sum(taper_wt(1:K)**2) / factor

                do  l= 0, lmax-lwin, 1
                    g_power(1,l+1) = dot_product(sgg(l+1,1:K), taper_wt(1:K)) &
                                     / sum(taper_wt(1:K))
                    t_power(1,l+1) = dot_product(stt(l+1,1:K), taper_wt(1:K)) &
                                     / sum(taper_wt(1:K))
                    gt_power(1,l+1) = dot_product(sgt(l+1,1:K), &
                                      taper_wt(1:K)) / sum(taper_wt(1:K))

                    if (K > 1) then
                        g_power(2,l+1) = dot_product( (sgg(l+1,1:K) &
                                         - g_power(1,l+1) )**2, &
                                         taper_wt(1:K) ) * factor
                        t_power(2,l+1) = dot_product( (stt(l+1,1:K) &
                                         - t_power(1,l+1) )**2, &
                                         taper_wt(1:K) ) * factor
                        gt_power(2,l+1) = dot_product( (sgt(l+1,1:K) &
                                          - gt_power(1,l+1) )**2, &
                                          taper_wt(1:K)) * factor

                    end if

                end do

            else
                do l = 0, lmax-lwin, 1
                    g_power(1,l+1) = sum(sgg(l+1,1:K)) / dble(K)
                    t_power(1,l+1) = sum(stt(l+1,1:K)) / dble(K)
                    gt_power(1,l+1) = sum(sgt(l+1,1:K)) / dble(K)

                    if (K > 1) then
                        g_power(2,l+1) = sum( ( sgg(l+1,1:K) &
                                         - g_power(1,l+1) )**2 ) / dble(K-1) &
                                         / dble(K) ! standard error!
                        t_power(2,l+1) = sum( ( stt(l+1,1:K) &
                                         - t_power(1,l+1) )**2 ) / dble(K-1) &
                                         / dble(K) ! standard error!
                        gt_power(2,l+1) = sum( ( sgt(l+1,1:K) &
                                          - gt_power(1,l+1) )**2 ) / dble(K-1) &
                                          / dble(K) ! standard error!

                    end if

                end do

            end if

            if (K > 1) then
                g_power(2,1:lmax-lwin+1) = sqrt(g_power(2,1:lmax-lwin+1))
                t_power(2,1:lmax-lwin+1) = sqrt(t_power(2,1:lmax-lwin+1))
                gt_power(2,1:lmax-lwin+1) = sqrt(gt_power(2,1:lmax-lwin+1))

            end if

        end do

        admit(1:lmax-lwin+1) = gt_power(1,1:lmax-lwin+1) &
                               / t_power(1,1:lmax-lwin+1)
        corr(1:lmax-lwin+1) = gt_power(1,1:lmax-lwin+1) &
                              / sqrt(t_power(1,1:lmax-lwin+1) &
                              * g_power(1,1:lmax-lwin+1))

        if (K > 1 .and. present(admit_error)) then
            admit_error(1:lmax-lwin+1) = ( gt_power(2,1:lmax-lwin+1) &
                                           / t_power(1,1:lmax-lwin+1) )**2 + &
                                           ( gt_power(1,1:lmax-lwin+1) &
                                           / t_power(1,1:lmax-lwin+1)**2 &
                                           * t_power(2,1:lmax-lwin+1) )**2
            admit_error(1:lmax-lwin+1) = sqrt(admit_error(1:lmax-lwin+1))

        end if

        if (K > 1 .and. present(corr_error)) then
            corr_error(1:lmax-lwin+1) = gt_power(2,1:lmax-lwin+1)**2 &
                            / t_power(1,1:lmax-lwin+1) &
                            / g_power(1,1:lmax-lwin+1) + &
                            ( gt_power(1,1:lwin-lmax+1) &
                            * t_power(2,1:lwin-lmax+1) / &
                            sqrt(g_power(1,1:lmax-lwin+1)) / 2.0_dp &
                            / t_power(1,1:lmax-lwin+1)**(3.0_dp/2.0_dp))**2 + &
                            ( gt_power(1,1:lwin-lmax+1) &
                            * g_power(2,1:lwin-lmax+1) &
                            / sqrt(t_power(1,1:lmax-lwin+1)) / 2.0_dp &
                            / g_power(1,1:lmax-lwin+1)**(3.0_dp/2.0_dp))**2
            corr_error(1:lmax-lwin+1) = sqrt(corr_error(1:lmax-lwin+1))

        end if

        if (K == 1 .and. present(k1linsig) .and. present(admit_error)) then
            admit_error = 0.0_dp

            if (k1linsig == 1) then
                do l = 1, lmax-lwin
                    admit_error(l+1) = g_power(1,l+1)*(1.0_dp - corr(l+1)**2) &
                                       / ( t_power(1,l+1) * dble(2*l) )
                end do

                admit_error(1:lmax-lwin+1) = sqrt(admit_error(1:lmax-lwin+1))

            end if

        end if

    !--------------------------------------------------------------------------
    !
    !   Calculate the admittance and correlation for each individual taper.
    !   Then average these in order to get the multitaper estimates and
    !   uncertainties.
    !
    !--------------------------------------------------------------------------
    else
        do i = 1, K
            shwin = 0.0_dp

            if (taper_order(i) < 0) then
                shwin(2,1:lwin+1,abs(taper_order(i))+1) = tapers(1:lwin+1,i)

            else
                shwin(1,1:lwin+1,taper_order(i)+1) = tapers(1:lwin+1,i)

            end if

            if (present(exitstatus)) then
                call SHRotateRealCoef(shwinrot, shwin, lwin, x, dj, &
                                      exitstatus = exitstatus)
                if (exitstatus /= 0) return
                call MakeGridGLQ(gridwinglq, shwinrot(1:2,1:lwin+1, 1:lwin+1),&
                                 lmaxwin, zero = zero, csphase = phase, &
                                 norm = 1, exitstatus = exitstatus)
                if (exitstatus /= 0) return
                temp(1:nlat,1:nlong) = gridtglq(1:nlat,1:nlong) &
                                       * gridwinglq(1:nlat,1:nlong)
                call SHExpandGLQ(shloc_t, lmaxwin, temp, w, zero = zero, &
                                csphase = phase, norm = 1, &
                                exitstatus = exitstatus)
                if (exitstatus /= 0) return
                temp(1:nlat,1:nlong) = gridgglq(1:nlat,1:nlong) &
                                       * gridwinglq(1:nlat,1:nlong)
                call SHExpandGLQ(shloc_g, lmaxwin, temp, w, zero = zero, &
                                 csphase = phase, norm = 1, &
                                 exitstatus = exitstatus)
                if (exitstatus /= 0) return
                call SHCrossPowerSpectrum(shloc_g, shloc_t, lmax-lwin, &
                                          sgt(:,i), exitstatus = exitstatus)
                if (exitstatus /= 0) return
                call SHPowerSpectrum(shloc_g, lmax-lwin, sgg(:,i), &
                                     exitstatus = exitstatus)
                if (exitstatus /= 0) return
                call SHPowerSpectrum(shloc_t, lmax-lwin, stt(:,i), &
                                     exitstatus = exitstatus)
                if (exitstatus /= 0) return
            else
                call SHRotateRealCoef(shwinrot, shwin, lwin, x, dj)
                call MakeGridGLQ(gridwinglq, shwinrot(1:2,1:lwin+1, 1:lwin+1),&
                                 lmaxwin, zero = zero, csphase = phase, &
                                 norm = 1)
                temp(1:nlat,1:nlong) = gridtglq(1:nlat,1:nlong) &
                                       * gridwinglq(1:nlat,1:nlong)
                call SHExpandGLQ(shloc_t, lmaxwin, temp, w, zero = zero, &
                                csphase = phase, norm = 1)
                temp(1:nlat,1:nlong) = gridgglq(1:nlat,1:nlong) &
                                       * gridwinglq(1:nlat,1:nlong)
                call SHExpandGLQ(shloc_g, lmaxwin, temp, w, zero = zero, &
                                csphase = phase, norm = 1)
                call SHCrossPowerSpectrum(shloc_g, shloc_t, lmax-lwin, &
                                          sgt(:,i))
                call SHPowerSpectrum(shloc_g, lmax-lwin, sgg(:,i))
                call SHPowerSpectrum(shloc_t, lmax-lwin, stt(:,i))
            end if

            admit_k(1:lmax-lwin+1, i) = sgt(1:lmax-lwin+1, i) &
                                         / stt(1:lmax-lwin+1, i)
            corr_k(1:lmax-lwin+1, i) = sgt(1:lmax-lwin+1, i) &
                                        / sqrt(stt(1:lmax-lwin+1, i)) &
                                        / sqrt(sgg(1:lmax-lwin+1, i))

        end do

        do l = 0, lmax-lwin, 1
            admit(l+1) = sum(admit_k(l+1,1:K)) / dble(K)
            corr(l+1) = sum(corr_k(l+1,1:K)) / dble(K)

        end do

        if (present(admit_error) .or. present(corr_error)) then
            do l = 0, lmax-lwin, 1
                if (present(admit_error)) then
                    if (K > 1) then 
                        admit_error(l+1) = sum( ( admit_k(l+1,1:K) &
                                           - admit(l+1) )**2 ) / dble(K-1) &
                                           / dble(K) ! standard error!
                        admit_error(l+1) = sqrt(admit_error(l+1))

                    else if (K == 1 .and. present(k1linsig)) then
                        if (k1linsig == 1) then
                            if (l == 0) then
                                admit_error(1) = 0.0_dp

                            else
                                admit_error(l+1) = sgg(l+1,1)*(1.0_dp &
                                                   - corr(l+1)**2) &
                                                   / (stt(l+1,1) * dble(2*l))
                                admit_error(l+1) = sqrt(admit_error(l+1))

                            end if

                        end if

                    end if

                end if

                if (present(corr_error)) then
                    if (K > 1) then
                        ! standard error!
                        corr_error(l+1) = sum((corr_k(l+1,1:K) &
                                               - corr(l+1))**2 ) &
                                          / dble(K-1) / dble(K)
                        corr_error(l+1) = sqrt(corr_error(l+1))

                    end if

                end if

            end do

        end if

    end if

    deallocate (shwin)
    deallocate (shwinrot)
    deallocate (shloc_g)
    deallocate (shloc_t)
    deallocate (gridtglq)
    deallocate (gridgglq)
    deallocate (gridwinglq)
    deallocate (temp)

end subroutine SHLocalizedAdmitCorr
