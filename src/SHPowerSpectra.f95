function SHPowerL(cilm, l)
!------------------------------------------------------------------------------
!
!   This function will compute the dimensionless power at
!   degree l corresponding to the 4pi spherical harmonic coefficients \
!   cilm(i, l, m)
!
!   Power = Sum_{m=0}^l ( c1lm**2 + c2lm**2 )
!
!   Calling Parameters
!
!       cilm    Spherical harmonic coefficients, dimensioned as
!               (2, lmax+1, lmax+1).
!       l       Spherical harmonic degree to compute power.
!
!   Copyright (c) 2005-2019, SHTOOLS
!   All rights reserved.
!
!------------------------------------------------------------------------------
    use ftypes

    implicit none

    real(dp) :: SHPowerL
    real(dp), intent(in) :: cilm(:,:,:)
    integer, intent(in) :: l
    integer i, m, l1, m1

    l1 = l + 1

    if (size(cilm(:,1,1)) < 2 .or. size(cilm(1,:,1)) < l1 .or. &
            size(cilm(1,1,:)) < l1) then
        print*, "Error --- SHPowerL"
        print*, "CILM must be dimensioned as (2, L+1, L+1) where L is ", l
        print*, "Input array is dimensioned ", size(cilm(:,1,1)), &
                size(cilm(1,:,1)), size(cilm(1,1,:))
        stop

    end if

    SHPowerL = cilm(1, l1, 1)**2  ! m=0 term

    do m = 1, l, 1
        m1 = m + 1

        do i = 1, 2
            SHPowerL = SHPowerL + cilm(i, l1, m1)**2

        end do

    end do

end function SHPowerL


function SHPowerDensityL(cilm, l)
!------------------------------------------------------------------------------
!
!   This function will compute the dimensionless power per coefficient
!   (density) at degree l of the 4pi spherical harmonic coefficients
!   cilm(i, l, m)
!
!   PowerSpectralDensity = Sum_{m=0}^l ( c1lm**2 + c2lm**2 ) / (2l+1)
!
!   Calling Parameters
!
!       cilm    Spherical harmonic coefficients, dimensioned as
!               (2, lmax+1, lmax+1).
!       l       Spherical harmonic degree to compute power.
!
!   Copyright (c) 2016, SHTOOLS
!   All rights reserved.
!
!------------------------------------------------------------------------------
    use ftypes

    implicit none

    real(dp) :: SHPowerDensityL
    real(dp), intent(in) :: cilm(:,:,:)
    integer, intent(in) :: l
    integer :: i, m, l1, m1

    l1 = l + 1

    if (size(cilm(:,1,1)) < 2 .or. size(cilm(1,:,1)) < l1 &
            .or. size(cilm(1,1,:)) < l1) then
        print*, "Error --- SHPowerDensityL"
        print*, "CILM must be dimensioned as (2, L+1, L+1) where L is ", l
        print*, "Input array is dimensioned ", size(cilm(:,1,1)), &
                size(cilm(1,:,1)), size(cilm(1,1,:))
        stop

    end if

    SHPowerDensityL = cilm(1, l1, 1)**2
    ! m=0 term

    do m = 1, l, 1
        m1 = m + 1

        do i = 1, 2
            SHPowerDensityL = SHPowerDensityL + cilm(i, l1, m1)**2

        end do

    end do

    SHPowerDensityL = SHPowerDensityL / dble(2*l+1)

end function SHPowerDensityL


function SHCrossPowerL(cilm1, cilm2, l)
!------------------------------------------------------------------------------
!
!   This function will compute the dimensionless cross power at
!   degree l for the 4pi spherical harmonic coefficients cilm1(i, l, m)
!   and cilm2(i,l,m).
!
!   CrossPower = Sum_{m=0}^l ( c1lm1*c1lm2 + c2lm1*c2lm2 )
!
!   Calling Parameters
!
!       cilm1   Spherical harmonic coefficients, dimensioned as
!               (2, lmax+1, lmax+1).
!       cilm2   Spherical harmonic coefficients, dimensioned as
!               (2, lmax+1, lmax+1).
!       l       Spherical harmonic degree to compute power.
!
!   Copyright (c) 2016, SHTOOLS
!   All rights reserved.
!
!------------------------------------------------------------------------------
    use ftypes

    implicit none

    real(dp) :: SHCrossPowerL
    real(dp), intent(in) :: cilm1(:,:,:), cilm2(:,:,:)
    integer, intent(in) :: l
    integer :: i, m, l1, m1

    l1 = l + 1

    if (size(cilm1(:,1,1)) < 2 .or. size(cilm1(1,:,1)) < l1 &
            .or. size(cilm1(1,1,:)) < l1) then
        print*, "Error --- SHCrossPowerL"
        print*, "CILM1 must be dimensioned as (2, L+1, L+1) where L is ", l
        print*, "Input array is dimensioned ", size(cilm1(:,1,1)), &
                size(cilm1(1,:,1)), size(cilm1(1,1,:))
        stop

    else if (size(cilm2(:,1,1)) < 2 .or. size(cilm2(1,:,1)) < l1 &
            .or. size(cilm2(1,1,:)) < l1) then
        print*, "Error --- SHCrossPowerL"
        print*, "CILM2 must be dimensioned as (2, L+1, L+1) where L is ", l
        print*, "Input array is dimensioned ", size(cilm2(:,1,1)), &
                size(cilm2(1,:,1)), size(cilm2(1,1,:))
        stop

    end if

    SHCrossPowerL = cilm1(1, l1, 1) * cilm2(1,l1,1)

    do m = 1, l, 1
        m1 = m + 1

        do i = 1, 2
            SHCrossPowerL = SHCrossPowerL + cilm1(i,l1,m1) * cilm2(i,l1,m1)

        end do

    end do

end function SHCrossPowerL


function SHCrossPowerDensityL(cilm1, cilm2, l)
!------------------------------------------------------------------------------
!
!   This function will compute the dimensionless cross power
!   (density) at degree l of the 4pi spherical harmonic coefficients 
!   cilm1(i, l, m) and cilm2(i,l,m).
!
!   CrossPower = Sum_{m=0}^l ( c1lm1*c1lm2 + c2lm1*c2lm2 ) / (2l+1)
!
!   Calling Parameters
!
!       cilm1   Spherical harmonic coefficients, dimensioned as
!               (2, lmax+1, lmax+1).
!       cilm2   Spherical harmonic coefficients, dimensioned as
!               (2, lmax+1, lmax+1).
!       l       Spherical harmonic degree to compute power.
!
!   Copyright (c) 2016, SHTOOLS
!   All rights reserved.
!
!------------------------------------------------------------------------------
    use ftypes

    implicit none

    real(dp) :: SHCrossPowerDensityL
    real(dp), intent(in) :: cilm1(:,:,:), cilm2(:,:,:)
    integer, intent(in) :: l
    integer :: i, m, l1, m1

    l1 = l + 1

    if (size(cilm1(:,1,1)) < 2 .or. size(cilm1(1,:,1)) < l1 &
            .or. size(cilm1(1,1,:)) < l1) then
        print*, "Error --- SHCrossPowerDensityL"
        print*, "CILM1 must be dimensioned as (2, L+1, L+1) where L is ", l
        print*, "Input array is dimensioned ", size(cilm1(:,1,1)), &
                size(cilm1(1,:,1)), size(cilm1(1,1,:))
        stop

    else if (size(cilm2(:,1,1)) < 2 .or. size(cilm2(1,:,1)) < l1 &
            .or. size(cilm2(1,1,:)) < l1) then
        print*, "Error --- SHCrossPowerDensityL"
        print*, "CILM2 must be dimensioned as (2, L+1, L+1) where L is ", l
        print*, "Input array is dimensioned ", size(cilm2(:,1,1)), &
                size(cilm2(1,:,1)), size(cilm2(1,1,:))
        stop

    end if

    SHCrossPowerDensityL = cilm1(1, l1, 1) * cilm2(1,l1,1)

    do m = 1, l, 1
        m1 = m + 1
        do i = 1, 2
            SHCrossPowerDensityL = SHCrossPowerDensityL &
                                   + cilm1(i,l1,m1) * cilm2(i,l1,m1)

        end do

    end do

    SHCrossPowerDensityL = SHCrossPowerDensityL / dble(2*l+1)

end function SHCrossPowerDensityL


subroutine SHPowerSpectrum(cilm, lmax, spectra, exitstatus)
!------------------------------------------------------------------------------
!
!   This function will compute the dimensionless power spectrum
!   of the 4pi spherical harmonic coefficients cilm(i, l, m).
!
!   Power(l) = Sum_{m=0}^l ( c1lm**2 + c2lm**2 )
!
!   Calling Parameters
!
!       IN
!           cilm        Spherical harmonic coefficients, dimensioned as
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
!   Copyright (c) 2016, SHTOOLS
!   All rights reserved.
!
!------------------------------------------------------------------------------
    use ftypes

    implicit none

    real(dp), intent(in) :: cilm(:,:,:)
    integer, intent(in) :: lmax
    real(dp), intent(out) :: spectra(:)
    integer, intent(out), optional :: exitstatus
    integer :: i, m, l1, m1, l

    if (present(exitstatus)) exitstatus = 0

    if (size(cilm(:,1,1)) < 2 .or. size(cilm(1,:,1)) < lmax+1 &
            .or. size(cilm(1,1,:)) < lmax+1) then
        print*, "Error --- SHPowerSpectrum"
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
        print*, "Error --- SHPowerSpectrum"
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

        spectra(l1) = cilm(1, l1, 1)**2

        do m = 1, l, 1
            m1 = m + 1

            do i = 1, 2
                spectra(l1) = spectra(l1) + cilm(i, l1, m1)**2

            end do

        end do

    end do

end subroutine SHPowerSpectrum


subroutine SHPowerSpectrumDensity(cilm, lmax, spectra, exitstatus)
!------------------------------------------------------------------------------
!
!   This function will compute the dimensionless power spectrum
!   density of the 4pi spherical harmonic coefficients cilm(i, l, m).
!
!   PowerSpectralDensity = Sum_{m=0}^l ( c1lm**2 + c2lm**2 ) / (2l+1)
!
!   Calling Parameters
!
!       IN
!           cilm        Spherical harmonic coefficients, dimensioned as
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
!   Copyright (c) 2016, SHTOOLS
!   All rights reserved.
!
!------------------------------------------------------------------------------
    use ftypes

    implicit none

    real(dp), intent(in) :: cilm(:,:,:)
    integer, intent(in) :: lmax
    real(dp), intent(out) :: spectra(:)
    integer, intent(out), optional :: exitstatus
    integer :: i, m, l1, m1, l

    if (present(exitstatus)) exitstatus = 0

    if (size(cilm(:,1,1)) < 2 .or. size(cilm(1,:,1)) < lmax+1 &
            .or. size(cilm(1,1,:)) < lmax+1) then
        print*, "Error --- SHPowerSpectrumDensity"
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
        print*, "Error --- SHPowerSpectrumDensity"
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

        spectra(l1) = cilm(1, l1, 1)**2

        do m = 1, l, 1
            m1 = m + 1

            do i = 1, 2
                spectra(l1) = spectra(l1) + cilm(i, l1, m1)**2

            end do

        end do

        spectra(l1) = spectra(l1) / dble(2*l+1)

    end do

end subroutine SHPowerSpectrumDensity


subroutine SHCrossPowerSpectrum(cilm1, cilm2, lmax, cspectra, exitstatus)
!-------------------------------------------------------------------------------
!
!   This function will compute the dimensionless cross power spectrum
!   of the 4pi spherical harmonic coefficients c1(i, l, m) and c2(1,l,m).
!
!   CrossPower(l) = Sum_{m=0}^l ( c1lm1*c1lm2 + c2lm1*c2lm2 )
!
!   Calling Parameters
!
!       IN
!           cilm1       Spherical harmonic coefficients, dimensioned as
!                       (2, lmax+1, lmax+1).
!           cilm2       Spherical harmonic coefficients, dimensioned as
!                       (2, lmax+1, lmax+1).
!           lmax        Spherical harmonic degree to compute power.
!       OUT
!           cspectra    Array of length (lmax+1) containing the cross power
!                       spectra of cilm1 and cilm2.
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
!   Copyright (c) 2016, SHTOOLS
!   All rights reserved.
!
!-------------------------------------------------------------------------------
    use ftypes

    implicit none

    real(dp), intent(in) :: cilm1(:,:,:), cilm2(:,:,:)
    integer, intent(in) :: lmax
    real(dp), intent(out) :: cspectra(:)
    integer, intent(out), optional :: exitstatus
    integer :: i, m, l1, m1, l

    if (present(exitstatus)) exitstatus = 0

    if (size(cilm1(:,1,1)) < 2 .or. size(cilm1(1,:,1)) < lmax+1 &
            .or. size(cilm1(1,1,:)) < lmax+1) then
        print*, "Error --- SHCrossPowerSpectrum"
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
        print*, "Error --- SHCrossPowerSpectrum"
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
        print*, "Error --- SHCrossPowerSpectrum"
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

        cspectra(l1) = cilm1(1, l1, 1) * cilm2(1, l1, 1)

        do m = 1, l, 1
            m1 = m + 1

            do i = 1, 2
                cspectra(l1) = cspectra(l1) + cilm1(i,l1,m1) * cilm2(i, l1, m1)

            end do

        end do

    end do

end subroutine SHCrossPowerSpectrum


subroutine SHCrossPowerSpectrumDensity(cilm1, cilm2, lmax, cspectra, &
                                       exitstatus)
!------------------------------------------------------------------------------
!
!   This function will compute the dimensionless cross power spectrum
!   density of the 4pi spherical harmonic coefficients cilm1(i, l, m) and
!   cilm2(i,l,m).
!
!   CrossPower = Sum_{m=0}^l ( c1lm1*c1lm2 + c2lm1*c2lm2 ) / (2l+1)
!
!   Calling Parameters
!
!       IN
!           cilm1       Spherical harmonic coefficients, dimensioned as
!                       (2, lmax+1, lmax+1).
!           cilm2       Spherical harmonic coefficients, dimensioned as
!                       (2, lmax+1, lmax+1).
!           lmax        Spherical harmonic degree to compute power.
!
!       OUT
!           cspectra    Array of length (lmax+1) containing the cross power
!                       spectral density of cilm1 and cilm2.
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
!   Copyright (c) 2016, SHTOOLS
!   All rights reserved.
!
!------------------------------------------------------------------------------
    use ftypes

    implicit none

    real(dp), intent(in) :: cilm1(:,:,:), cilm2(:,:,:)
    integer, intent(in) :: lmax
    real(dp), intent(out) :: cspectra(:)
    integer, intent(out), optional :: exitstatus
    integer :: i, m, l1, m1, l

    if (size(cilm1(:,1,1)) < 2 .or. size(cilm1(1,:,1)) < lmax+1 &
            .or. size(cilm1(1,1,:)) < lmax+1) then
        print*, "Error --- SHCrossPowerSpectrumDensity"
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
        print*, "Error --- SHCrossPowerSpectrumDensity"
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
        print*, "Error --- SHCrossPowerSpectrumDensity"
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
        
        cspectra(l1) = cilm1(1, l1, 1) * cilm2(1, l1, 1)

        do m = 1, l, 1
            m1 = m + 1

            do i = 1, 2
                cspectra(l1) = cspectra(l1) + cilm1(i,l1,m1) * cilm2(i, l1, m1)

            end do

        end do

        cspectra(l1) = cspectra(l1) / dble(2*l+1)

    end do

end subroutine SHCrossPowerSpectrumDensity
