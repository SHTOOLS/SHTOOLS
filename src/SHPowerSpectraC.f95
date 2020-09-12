function SHPowerLC(cilm, l)
!------------------------------------------------------------------------------
!
!   This function will compute the dimensionless power at
!   degree l corresponding to the complex 4-pi normalized spherical harmonic
!   coefficients cilm(i, l, m).
!
!   Power(l) = Sum_{m=0}^l | c1lm | **2 + | c2lm | **2 )
!
!   Calling Parameters
!
!       cilm    Complex spherical harmonic coefficients, dimensioned as
!               (2, lmax+1, lmax+1).
!       l       Spherical harmonic degree to compute power.
!
!   Copyright (c) 2016-2019, SHTOOLS
!   All rights reserved.
!
!------------------------------------------------------------------------------
    use ftypes

    implicit none

    real(dp) :: SHPowerLC
    complex(dp), intent(in) :: cilm(:,:,:)
    integer, intent(in) :: l
    integer :: i, m, l1, m1

    l1 = l + 1

    if (size(cilm(:,1,1)) < 2 .or. size(cilm(1,:,1)) < l1 &
            .or. size(cilm(1,1,:)) < l1) then
        print*, "Error --- SHPowerLC"
        print*, "CILM must be dimensioned as (2, L+1, L+1) where L is ", l
        print*, "Input array is dimensioned ", size(cilm(:,1,1)), &
                size(cilm(1,:,1)), size(cilm(1,1,:))
        stop

    end if

    SHPowerLC = dble(cilm(1, l1, 1))**2 + aimag(cilm(1, l1, 1))**2

    do m = 1, l, 1
        m1 = m + 1

        do i = 1, 2
            SHPowerLC = SHPowerLC + dble(cilm(i, l1, m1))**2 &
                                  + aimag(cilm(i, l1, m1))**2

        end do

    end do

end function SHPowerLC


function SHPowerDensityLC(cilm, l)
!------------------------------------------------------------------------------
!
!   This function will compute the dimensionless power per coefficient
!   (density) at degree l of the complex 4-pi normalized spherical harmonic
!   coefficients cilm(i, l, m).
!
!   PowerSpectralDensity(l) = Sum_{m=0}^l ( | c1lm |**2 + | c2lm |**2 ) /(2l+1)
!
!   Calling Parameters
!
!       cilm    Complex spherical harmonic coefficients, dimensioned as
!               (2, lmax+1, lmax+1).
!       l       Spherical harmonic degree to compute power.
!
!   Copyright (c) 2016-2019, SHTOOLS
!   All rights reserved.
!
!------------------------------------------------------------------------------
    use ftypes

    implicit none

    real(dp) :: SHPowerDensityLC
    complex(dp), intent(in) :: cilm(:,:,:)
    integer, intent(in) :: l
    integer :: i, m, l1, m1

    l1 = l + 1

    if (size(cilm(:,1,1)) < 2 .or. size(cilm(1,:,1)) < l1 &
            .or. size(cilm(1,1,:)) < l1) then
        print*, "Error --- SHPowerDensityLC"
        print*, "CILM must be dimensioned as (2, L+1, L+1) where L is ", l
        print*, "Input array is dimensioned ", size(cilm(:,1,1)), &
                size(cilm(1,:,1)), size(cilm(1,1,:))
        stop

    end if

    SHPowerDensityLC = dble(cilm(1, l1, 1))**2 + aimag(cilm(1, l1, 1))**2

    do m = 1, l, 1
        m1 = m + 1

        do i = 1, 2
            SHPowerDensityLC = SHPowerDensityLC + dble(cilm(i, l1, m1))**2 &
                                                + aimag(cilm(i, l1, m1))**2

        end do
    end do

    SHPowerDensityLC = SHPowerDensityLC / dble(2*l+1)

end function SHPowerDensityLC


function SHCrossPowerLC(cilm1, cilm2, l)
!------------------------------------------------------------------------------
!
!   This function will compute the dimensionless cross power at
!   degree l for the complex 4-pi normalized spherical harmonic
!   coefficients cilm1(i, l, m) and cilm2(i,l,m).
!
!   CrossPower = Sum_{m=0}^l ( c11lm * c21lm^* + c12lm * c22lm^* )
!
!   Calling Parameters
!
!       cilm1   Complex spherical harmonic coefficients, dimensioned as
!               (2, lmax+1, lmax+1).
!       cilm2   Complex spherical harmonic coefficients, dimensioned as
!               (2, lmax+1, lmax+1).
!       l       Spherical harmonic degree to compute power.
!
!   Copyright (c) 2016-2019, SHTOOLS
!   All rights reserved.
!
!------------------------------------------------------------------------------
    use ftypes

    implicit none

    complex(dp) :: SHCrossPowerLC
    complex(dp), intent(in) :: cilm1(:,:,:), cilm2(:,:,:)
    integer, intent(in) :: l
    integer :: i, m, l1, m1

    l1 = l + 1

    if (size(cilm1(:,1,1)) < 2 .or. size(cilm1(1,:,1)) < l1 &
            .or. size(cilm1(1,1,:)) < l1) then
        print*, "Error --- SHCrossPowerLC"
        print*, "CILM1 must be dimensioned as (2, L+1, L+1) where L is ", l
        print*, "Input array is dimensioned ", size(cilm1(:,1,1)), &
                size(cilm1(1,:,1)), size(cilm1(1,1,:))
        stop

    else if (size(cilm2(:,1,1)) < 2 .or. size(cilm2(1,:,1)) < l1 .or. &
             size(cilm2(1,1,:)) < l1) then
        print*, "SHCrossPowerLC --- Error"
        print*, "CILM2 must be dimensioned as (2, L+1, L+1) where L is ", l
        print*, "Input array is dimensioned ", size(cilm2(:,1,1)), &
                size(cilm2(1,:,1)), size(cilm2(1,1,:))
        stop

    end if

    SHCrossPowerLC = cilm1(1, l1, 1) * conjg(cilm2(1,l1,1))

    do m = 1, l, 1
        m1 = m + 1

        do i = 1, 2
            SHCrossPowerLC = SHCrossPowerLC + cilm1(i,l1,m1) &
                                              * conjg(cilm2(i,l1,m1))

        end do

    end do

end function SHCrossPowerLC


function SHCrossPowerDensityLC(cilm1, cilm2, l)
!------------------------------------------------------------------------------
!
!   This function will compute the dimensionless cross power
!   (density) at degree l of the complex 4-pi normalized spherical harmonic
!   coefficients cilm1(i, l, m) and cilm2(i,l,m).
!
!   CrossPower(l) = Sum_{m=0}^l ( c11lm * C21lm^* + c12lm * C22lm^* ) / (2l+1)
!
!   Calling Parameters
!
!       cilm1   Complex spherical harmonic coefficients, dimensioned as
!               (2, lmax+1, lmax+1).
!       cilm2   Complex spherical harmonic coefficients, dimensioned as
!               (2, lmax+1, lmax+1).
!       l       Spherical harmonic degree to compute power.
!
!   Copyright (c) 2016-2019, SHTOOLS
!   All rights reserved.
!
!------------------------------------------------------------------------------
    use ftypes

    implicit none

    complex(dp) :: SHCrossPowerDensityLC
    complex(dp), intent(in) :: cilm1(:,:,:), cilm2(:,:,:)
    integer, intent(in) :: l
    integer :: i, m, l1, m1

    l1 = l + 1

    if (size(cilm1(:,1,1)) < 2 .or. size(cilm1(1,:,1)) < l1 &
            .or. size(cilm1(1,1,:)) < l1) then
        print*, "Error --- SHCrossPowerDensityLC"
        print*, "CILM1 must be dimensioned as (2, L+1, L+1) where L is ", l
        print*, "Input array is dimensioned ", size(cilm1(:,1,1)), &
                size(cilm1(1,:,1)), size(cilm1(1,1,:))
        stop

    else if (size(cilm2(:,1,1)) < 2 .or. size(cilm2(1,:,1)) < l1 &
            .or. size(cilm2(1,1,:)) < l1) then
        print*, "Error --- SHCrossPowerDensityLC"
        print*, "CILM2 must be dimensioned as (2, L+1, L+1) where L is ", l
        print*, "Input array is dimensioned ", size(cilm2(:,1,1)), &
                size(cilm2(1,:,1)), size(cilm2(1,1,:))
        stop

    end if

    SHCrossPowerDensityLC = cilm1(1, l1, 1) * conjg(cilm2(1,l1,1))

    do m = 1, l, 1
        m1 = m + 1

        do i = 1, 2
            SHCrossPowerDensityLC = SHCrossPowerDensityLC &
                                    + cilm1(i, l1, m1) * conjg(cilm2(i,l1,m1))

        end do

    end do

    SHCrossPowerDensityLC = SHCrossPowerDensityLC / dble(2*l+1)

end function SHCrossPowerDensityLC


subroutine SHPowerSpectrumC(cilm, lmax, spectra, exitstatus)
!------------------------------------------------------------------------------
!
!   This function will compute the dimensionless power spectrum
!   of the complex 4-pi normalized spherical harmonic coefficients
!   cilm(i, l, m).
!
!   Power(l) = Sum_{m=0}^l ( | c1lm |**2 + | c2lm |**2 )
!
!   Calling Parameters
!
!       IN
!           cilm        Complex pherical harmonic coefficients, dimensioned as
!                       (2, lmax+1, lmax+1).
!           lmax        Spherical harmonic degree to compute power.
!
!       OUT
!           spectra     Array of length (lmax+1) containing the power
!                       spectra of cilm.
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
!   Copyright (c) 2016-2019, SHTOOLS
!   All rights reserved.
!
!------------------------------------------------------------------------------
    use ftypes

    implicit none

    complex(dp), intent(in) :: cilm(:,:,:)
    integer, intent(in) :: lmax
    real(dp), intent(out) :: spectra(:)
    integer, intent(out), optional :: exitstatus
    integer :: i, m, l1, m1, l

    if (present(exitstatus)) exitstatus = 0

    if (size(cilm(:,1,1)) < 2 .or. size(cilm(1,:,1)) < lmax+1 &
            .or. size(cilm(1,1,:)) < lmax+1) then
        print*, "Error --- SHPowerSpectrumC"
        print*, "CILM must be dimensioned as (2, LMAX+1, LMAX+1) " // &
                "where LMAX is ", lmax
        print*, "Input array is dimensioned ", size(cilm(:,1,1)), &
                size(cilm(1,:,1)), size(cilm(1,1,:))
        if (present(exitstatus)) then
            exitstatus = 1
            return
        else
            stop
        end if

    else if (size(spectra) < lmax+1) then
        print*, "Error --- SHPowerSpectrumC"
        print*, "SPECTRA must be dimensioned as (LMAX+1) where LMAX is ", lmax
        print*, "Input vector has dimension ", size(spectra)
        if (present(exitstatus)) then
            exitstatus = 1
            return
        else
            stop
        end if

    end if

    spectra = 0.0_dp

    do l = 0, lmax
        l1 = l + 1

        spectra(l1) = dble(cilm(1, l1, 1))**2 + aimag(cilm(1, l1, 1))**2

        do m = 1, l, 1
            m1 = m + 1

            do i = 1, 2
                spectra(l1) = spectra(l1) + dble(cilm(i, l1, m1))**2 &
                                          + aimag(cilm(i, l1, m1))**2

            end do
        end do
    end do

end subroutine SHPowerSpectrumC


subroutine SHPowerSpectrumDensityC(cilm, lmax, spectra, exitstatus)
!------------------------------------------------------------------------------
!
!   This function will compute the dimensionless power spectrum
!   density of the complex 4-pi normalized spherical harmonic coefficients
!   cilm(i, l, m).
!
!   PowerSpectralDensityC(l) = Sum_{m=0}^l ( c1lm**2 + c2lm**2 ) / (2l+1)
!
!   Calling Parameters
!
!       IN
!           cilm        Complex spherical harmonic coefficients, dimensioned as
!                       (2, lmax+1, lmax+1).
!           lmax        Spherical harmonic degree to compute power.
!
!       OUT
!           spectra     Array of length (lmax+1) containing the power
!                       spectra density of cilm.
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
!   Copyright (c) 2016-2019, SHTOOLS
!   All rights reserved.
!
!-------------------------------------------------------------------------------
    use ftypes

    implicit none

    complex(dp), intent(in) :: cilm(:,:,:)
    integer, intent(in) :: lmax
    real(dp), intent(out) :: spectra(:)
    integer, intent(out), optional :: exitstatus
    integer :: i, m, l1, m1, l

    if (present(exitstatus)) exitstatus = 0

    if (size(cilm(:,1,1)) < 2 .or. size(cilm(1,:,1)) < lmax+1 &
            .or. size(cilm(1,1,:)) < lmax+1) then
        print*, "Error --- SHPowerSpectrumDensityC"
        print*, "CILM must be dimensioned as (2, LMAX+1, LMAX+1) " // &
                "where LMAX is ", lmax
        print*, "Input array is dimensioned ", size(cilm(:,1,1)), &
                size(cilm(1,:,1)), size(cilm(1,1,:))
        if (present(exitstatus)) then
            exitstatus = 1
            return
        else
            stop
        end if

    else if (size(spectra) < lmax+1) then
        print*, "Error --- SHPowerSpectrumDensityC"
        print*, "SPECTRA must be dimensioned as (LMAX+1) where LMAX is ", lmax
        print*, "Input vector has dimension ", size(spectra)
        if (present(exitstatus)) then
            exitstatus = 1
            return
        else
            stop
        end if

    end if

    spectra = 0.0_dp

    do l = 0, lmax
        l1 = l + 1

        spectra(l1) = dble(cilm(1, l1, 1))**2 + aimag(cilm(1, l1, 1))**2

        do m = 1, l, 1
            m1 = m + 1

            do i = 1, 2
                spectra(l1) = spectra(l1) + dble(cilm(i, l1, m1))**2 &
                                          + aimag(cilm(i, l1, m1))**2

            end do

        end do

        spectra(l1) = spectra(l1) / dble(2*l+1)

    end do

end subroutine SHPowerSpectrumDensityC


subroutine SHCrossPowerSpectrumC(cilm1, cilm2, lmax, cspectra, exitstatus)
!------------------------------------------------------------------------------
!
!   This function will compute the dimensionless cross power spectrum
!   of the complex 4-pi normalized spherical harmonic coefficients
!   cilm1(i, l, m) and cilm2(1,l,m).
!
!   CrossPower(l) = Sum_{m=0}^l ( c11lm * c21lm^* + c12lm * c22lm^* )
!
!   Calling Parameters
!
!       IN
!           cilm1       Complex spherical harmonic coefficients, dimensioned as
!                       (2, lmax+1, lmax+1).
!           cilm2       Complex spherical harmonic coefficients, dimensioned as
!                       (2, lmax+1, lmax+1).
!           lmax        Spherical harmonic degree to compute power.
!       OUT
!           cspectra    Array of length (lmax+1) containing the complex cross
!                       power spectra of cilm1 and cilm2.
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
!   Copyright (c) 2016-2019, SHTOOLS
!   All rights reserved.
!
!------------------------------------------------------------------------------
    use ftypes

    implicit none

    complex(dp), intent(in) :: cilm1(:,:,:), cilm2(:,:,:)
    integer, intent(in) :: lmax
    complex(dp), intent(out) :: cspectra(:)
    integer, intent(out), optional :: exitstatus
    integer :: i, m, l1, m1, l

    if (present(exitstatus)) exitstatus = 0

    if (size(cilm1(:,1,1)) < 2 .or. size(cilm1(1,:,1)) < lmax+1 &
            .or. size(cilm1(1,1,:)) < lmax+1) then
        print*, "Error --- SHCrossPowerSpectrumC"
        print*, "CILM1 must be dimensioned as (2, LMAX+1, LMAX+1) " // &
                "where lmax is", lmax
        print*, "Input array is dimensioned ", size(cilm1(:,1,1)), &
                size(cilm1(1,:,1)), size(cilm1(1,1,:))
        if (present(exitstatus)) then
            exitstatus = 1
            return
        else
            stop
        end if

    else if (size(cilm2(:,1,1)) < 2 .or. size(cilm2(1,:,1)) < lmax+1 &
            .or. size(cilm2(1,1,:)) < lmax+1) then
        print*, "Error --- SHCrossPowerSpectrumC"
        print*, "CILM2 must be dimensioned as (2, LMAX+1, LMAX+1)"
        print*, "Input array is dimensioned ", size(cilm2(:,1,1)), &
                size(cilm2(1,:,1)), size(cilm2(1,1,:))
        if (present(exitstatus)) then
            exitstatus = 1
            return
        else
            stop
        end if

    else if (size(cspectra) < lmax+1) then
        print*, "Error --- SHCrossPowerSpectrumC"
        print*, "CSPECTRA must be dimensioned as (LMAX+1) where lmax is ", lmax
        print*, "Input vector has dimension ", size(cspectra)
        if (present(exitstatus)) then
            exitstatus = 1
            return
        else
            stop
        end if

    end if

    cspectra = 0.0_dp

    do l = 0, lmax
        l1 = l + 1

        cspectra(l1) = cilm1(1, l1, 1) * conjg(cilm2(1, l1, 1))

        do m = 1, l, 1
            m1 = m + 1

            do i = 1, 2
                cspectra(l1) = cspectra(l1) &
                               + cilm1(i, l1, m1) * conjg(cilm2(i, l1, m1))

            end do

        end do

    end do

end subroutine SHCrossPowerSpectrumC


subroutine SHCrossPowerSpectrumDensityC(cilm1, cilm2, lmax, cspectra, &
                                        exitstatus)
!------------------------------------------------------------------------------
!
!   This function will compute the dimensionless cross power spectrum
!   density of the complex 4-pi normalized spherical harmonic coefficients
!   cilm1(i, l, m) and cilm2(i,l,m).
!
!   CrossPower(l) = Sum_{m=0}^l ( c11lm * c21lm^* + c12lm * c22lm^* ) / (2l+1)
!
!   Calling Parameters
!
!       IN
!           cilm1       Complex spherical harmonic coefficients, dimensioned as
!                       (2, lmax+1, lmax+1).
!           cilm2       Complex spherical harmonic coefficients, dimensioned as
!                       (2, lmax+1, lmax+1).
!           lmax        Spherical harmonic degree to compute power.
!
!       OUT
!           cspectra    Array of length (lmax+1) containing the complex cross
!                       power spectral density of cilm1 and cilm2.
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
!   Copyright (c) 2005-2019, SHTOOLS
!   All rights reserved.
!
!------------------------------------------------------------------------------
    use ftypes

    implicit none

    complex(dp), intent(in) :: cilm1(:,:,:), cilm2(:,:,:)
    integer, intent(in) :: lmax
    complex(dp), intent(out) :: cspectra(:)
    integer, intent(out), optional :: exitstatus
    integer :: i, m, l1, m1, l

    if (present(exitstatus)) exitstatus = 0

    if (size(cilm1(:,1,1)) < 2 .or. size(cilm1(1,:,1)) < lmax+1 &
            .or. size(cilm1(1,1,:)) < lmax+1) then
        print*, "Error --- SHCrossPowerSpectrumDensityC"
        print*, "CILM1 must be dimensioned as (2, LMAX+1, LMAX+1) " // &
                "where lmax is", lmax
        print*, "Input array is dimensioned ", size(cilm1(:,1,1)), &
                size(cilm1(1,:,1)), size(cilm1(1,1,:))
        if (present(exitstatus)) then
            exitstatus = 1
            return
        else
            stop
        end if

    else if (size(cilm2(:,1,1)) < 2 .or. size(cilm2(1,:,1)) < lmax+1 &
            .or. size(cilm2(1,1,:)) < lmax+1) then
        print*, "Error --- SHCrossPowerSpectrumDensityC"
        print*, "CILM2 must be dimensioned as (2, LMAX+1, LMAX+1)"
        print*, "Input array is dimensioned ", size(cilm2(:,1,1)), &
                size(cilm2(1,:,1)), size(cilm2(1,1,:))
        if (present(exitstatus)) then
            exitstatus = 1
            return
        else
            stop
        end if

    else if (size(cspectra) < lmax+1) then
        print*, "Error --- SHCrossPowerSpectrumDensityC"
        print*, "CSPECTRA must be dimensioned as (LMAX+1) where lmax is ", lmax
        print*, "Input vector has dimension ", size(cspectra)
        if (present(exitstatus)) then
            exitstatus = 1
            return
        else
            stop
        end if

    end if

    cspectra = 0.0_dp

    do l = 0, lmax
        l1 = l + 1

        cspectra(l1) = cilm1(1, l1, 1) * conjg(cilm2(1, l1, 1))

        do m = 1, l, 1
            m1 = m + 1

            do i = 1, 2
                cspectra(l1) = cspectra(l1) &
                               + cilm1(i, l1, m1) * conjg(cilm2(i, l1, m1))

            end do

        end do

        cspectra(l1) = cspectra(l1) / dble(2*l+1)

    end do

end subroutine SHCrossPowerSpectrumDensityC
