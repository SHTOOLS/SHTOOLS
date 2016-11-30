real*8 function SHPowerLC(c, l)
!------------------------------------------------------------------------------
!
!   This function will compute the dimensionless power at
!   degree l corresponding to the complex 4-pi normalized spherical harmonic
!   coefficients c(i, l, m).
!
!   Power(l) = Sum_{m=0}^l | C1lm | **2 + | C2lm | **2 )
!
!   Calling Parameters
!
!       c       Complex spherical harmonic coefficients, dimensioned as
!               (2, lmax+1, lmax+1).
!       l       Spherical harmonic degree to compute power.
!
!   Copyright (c) 2016, SHTOOLS
!   All rights reserved.
!
!------------------------------------------------------------------------------
    implicit none

    complex*16, intent(in) :: c(:,:,:)
    integer, intent(in) :: l
    integer i, m, l1, m1

    l1 = l + 1

    if (size(c(:,1,1)) < 2 .or. size(c(1,:,1)) < l1 &
            .or. size(c(1,1,:)) < l1) then
        print*, "Error --- SHPowerLC"
        print*, "C must be dimensioned as (2, L+1, L+1) where L is ", l
        print*, "Input array is dimensioned ", size(c(:,1,1)), &
                size(c(1,:,1)), size(c(1,1,:)) 
        stop

    end if

    SHPowerLC = dble(c(1, l1, 1))**2 + aimag(c(1, l1, 1))**2    ! m=0 term

    do m = 1, l, 1
        m1 = m + 1

        do i = 1, 2
            SHPowerLC = SHPowerLC + dble(c(i, l1, m1))**2 &
                                  + aimag(c(i, l1, m1))**2

        end do

    end do

end function SHPowerLC


real*8 function SHPowerDensityLC(c, l)
!------------------------------------------------------------------------------
!
!   This function will compute the dimensionless power per coefficient
!   (density) at degree l of the complex 4-pi normalized spherical harmonic
!   coefficients c(i, l, m).
!
!   PowerSpectralDensity(l) = Sum_{m=0}^l ( | C1lm |**2 + | C2lm |**2 ) / (2l+1)
!
!   Calling Parameters
!
!       c       Complex spherical harmonic coefficients, dimensioned as
!               (2, lmax+1, lmax+1).
!       l       Spherical harmonic degree to compute power.
!
!   Copyright (c) 2016, SHTOOLS
!   All rights reserved.
!
!------------------------------------------------------------------------------
    implicit none

    complex*16, intent(in) :: c(:,:,:)
    integer, intent(in) :: l
    integer i, m, l1, m1

    l1 = l + 1

    if (size(c(:,1,1)) < 2 .or. size(c(1,:,1)) < l1 &
            .or. size(c(1,1,:)) < l1) then
        print*, "Error --- SHPowerDensityLC"
        print*, "C must be dimensioned as (2, L+1, L+1) where L is ", l
        print*, "Input array is dimensioned ", size(c(:,1,1)), &
                size(c(1,:,1)), size(c(1,1,:)) 
        stop

    end if

    SHPowerDensityLC = dble(c(1, l1, 1))**2 + aimag(c(1, l1, 1))**2

    do m = 1, l, 1
        m1 = m + 1

        do i = 1, 2
            SHPowerDensityLC = SHPowerDensityLC + dble(c(i, l1, m1))**2 &
                                                + aimag(c(i, l1, m1))**2

        end do
    end do

    SHPowerDensityLC = SHPowerDensityLC / dble(2*l+1)

end function SHPowerDensityLC


complex*16 function SHCrossPowerLC(c1, c2, l)
!------------------------------------------------------------------------------
!
!   This function will compute the dimensionless cross power at
!   degree l for the complex 4-pi normalized spherical harmonic
!   coefficients c1(i, l, m) and c2(i,l,m).
!
!   CrossPower =  Sum_{m=0}^l ( C11lm*C21lm^* + C12lm*C22lm^* )
!
!   Calling Parameters
!
!       c1      Complex spherical harmonic coefficients, dimensioned as
!               (2, lmax+1, lmax+1).
!       c2      Complex spherical harmonic coefficients, dimensioned as
!               (2, lmax+1, lmax+1).
!       l       Spherical harmonic degree to compute power.
!
!   Copyright (c) 2016, SHTOOLS
!   All rights reserved.
!
!------------------------------------------------------------------------------
    implicit none

    complex*16, intent(in) :: c1(:,:,:), c2(:,:,:)
    integer, intent(in) :: l
    integer i, m, l1, m1

    l1 = l + 1

    if (size(c1(:,1,1)) < 2 .or. size(c1(1,:,1)) < l1 &
            .or. size(c1(1,1,:)) < l1) then
        print*, "Error --- SHCrossPowerLC"
        print*, "C1 must be dimensioned as (2, L+1, L+1) where L is ", l
        print*, "Input array is dimensioned ", size(c1(:,1,1)), &
                size(c1(1,:,1)), size(c1(1,1,:))
        stop

    else if (size(c2(:,1,1)) < 2 .or. size(c2(1,:,1)) < l1 .or. &
             size(c2(1,1,:)) < l1) then
        print*, "SHCrossPowerLC --- Error"
        print*, "C2 must be dimensioned as (2, L+1, L+1) where L is ", l
        print*, "Input array is dimensioned ", size(c2(:,1,1)), &
                size(c2(1,:,1)), size(c2(1,1,:))
        stop

    end if

    SHCrossPowerLC = c1(1, l1, 1) * dconjg(c2(1,l1,1))

    do m = 1, l, 1
        m1 = m + 1

        do i = 1, 2
            SHCrossPowerLC = SHCrossPowerLC + c1(i,l1,m1) * dconjg(c2(i,l1,m1))

        end do

    end do

end function SHCrossPowerLC


complex*16 function SHCrossPowerDensityLC(c1, c2, l)
!------------------------------------------------------------------------------
!
!   This function will compute the dimensionless cross power
!   (density) at degree l of the complex 4-pi normalized spherical harmonic
!   coefficients c1(i, l, m) and c2(i,l,m).
!
!   CrossPower(l) =  Sum_{m=0}^l ( C11lm*C21lm^* + C12lm*C22lm^* ) / (2l+1)
!
!   Calling Parameters
!
!       c1      Complex spherical harmonic coefficients, dimensioned as
!               (2, lmax+1, lmax+1).
!       c2      Complex spherical harmonic coefficients, dimensioned as
!               (2, lmax+1, lmax+1).
!       l       Spherical harmonic degree to compute power.
!
!   Copyright (c) 2016, SHTOOLS
!   All rights reserved.
!
!------------------------------------------------------------------------------
    implicit none

    complex*16, intent(in) :: c1(:,:,:), c2(:,:,:)
    integer, intent(in) :: l
    integer i, m, l1, m1

    l1 = l + 1

    if (size(c1(:,1,1)) < 2 .or. size(c1(1,:,1)) < l1 &
            .or. size(c1(1,1,:)) < l1) then
        print*, "Error --- SHCrossPowerDensityLC"
        print*, "C1 must be dimensioned as (2, L+1, L+1) where L is ", l
        print*, "Input array is dimensioned ", size(c1(:,1,1)), &
                size(c1(1,:,1)), size(c1(1,1,:))
        stop

    else if (size(c2(:,1,1)) < 2 .or. size(c2(1,:,1)) < l1 &
            .or. size(c2(1,1,:)) < l1) then
        print*, "Error --- SHCrossPowerDensityLC"
        print*, "C2 must be dimensioned as (2, L+1, L+1) where L is ", l
        print*, "Input array is dimensioned ", size(c2(:,1,1)), &
                size(c2(1,:,1)), size(c2(1,1,:))
        stop

    end if

    SHCrossPowerDensityLC =  c1(1, l1, 1) * dconjg(c2(1,l1,1))

    do m = 1, l, 1
        m1 = m + 1

        do i = 1, 2
            SHCrossPowerDensityLC = SHCrossPowerDensityLC &
                                    + c1(i, l1, m1)*dconjg(c2(i,l1,m1))

        end do

    end do

    SHCrossPowerDensityLC = SHCrossPowerDensityLC / dble(2*l+1)

end function SHCrossPowerDensityLC


subroutine SHPowerSpectrumC(c, lmax, spectra, exitstatus)
!------------------------------------------------------------------------------
!
!   This function will compute the dimensionless power spectrum
!   of the complex 4-pi normalized spherical harmonic coefficients c(i, l, m).
!
!   Power(l) = Sum_{m=0}^l ( | C1lm |**2 + | C2lm |**2 )
!
!   Calling Parameters
!
!       IN
!           c           Complex pherical harmonic coefficients, dimensioned as
!                       (2, lmax+1, lmax+1).
!           lmax        Spherical harmonic degree to compute power.
!
!       OUT
!           spectra     Array of length (lmax+1) containing the power
!                       spectra of c.
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
    implicit none

    complex*16, intent(in) :: c(:,:,:)
    integer, intent(in) :: lmax
    real*8, intent(out) ::  spectra(:)
    integer, intent(out), optional :: exitstatus
    integer i, m, l1, m1, l

    if (present(exitstatus)) exitstatus = 0

    if (size(c(:,1,1)) < 2 .or. size(c(1,:,1)) < lmax+1 &
            .or. size(c(1,1,:)) < lmax+1) then
        print*, "Error --- SHPowerSpectrumC"
        print*, "C must be dimensioned as (2, LMAX+1, LMAX+1) " // &
                "where LMAX is ", lmax
        print*, "Input array is dimensioned ", size(c(:,1,1)), &
                size(c(1,:,1)), size(c(1,1,:))
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

    spectra = 0.0d0

    do l = 0, lmax
        l1 = l + 1

        spectra(l1) = dble(c(1, l1, 1))**2 + aimag(c(1, l1, 1))**2

        do m = 1, l, 1
            m1 = m + 1

            do i = 1, 2
                spectra(l1) = spectra(l1) + dble(c(i, l1, m1))**2 &
                                          + aimag(c(i, l1, m1))**2

            end do
        end do
    end do

end subroutine SHPowerSpectrumC


subroutine SHPowerSpectrumDensityC(c, lmax, spectra, exitstatus)
!------------------------------------------------------------------------------
!
!   This function will compute the dimensionless power spectrum
!   density of the complex 4-pi normalized spherical harmonic coefficients
!   c(i, l, m).
!
!   PowerSpectralDensityC(l) = Sum_{m=0}^l ( C1lm**2 + C2lm**2 ) / (2l+1)
!
!   Calling Parameters
!
!       IN
!           c           Complex spherical harmonic coefficients, dimensioned as
!                       (2, lmax+1, lmax+1).
!           lmax        Spherical harmonic degree to compute power.
!
!       OUT
!           spectra     Array of length (lmax+1) containing the power
!                       spectra density of c.
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
    implicit none
    
    complex*16, intent(in) :: c(:,:,:)
    integer, intent(in) :: lmax
    real*8, intent(out) ::  spectra(:)
    integer, intent(out), optional :: exitstatus
    integer i, m, l1, m1, l

    if (present(exitstatus)) exitstatus = 0

    if (size(c(:,1,1)) < 2 .or. size(c(1,:,1)) < lmax+1 &
            .or. size(c(1,1,:)) < lmax+1) then
        print*, "Error --- SHPowerSpectrumDensityC"
        print*, "C must be dimensioned as (2, LMAX+1, LMAX+1) " // &
                "where LMAX is ", lmax
        print*, "Input array is dimensioned ", size(c(:,1,1)), &
                size(c(1,:,1)), size(c(1,1,:))
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

    spectra = 0.0d0

    do l = 0, lmax
        l1 = l + 1

        spectra(l1) = dble(c(1, l1, 1))**2 + aimag(c(1, l1, 1))**2

        do m = 1, l, 1
            m1 = m + 1

            do i = 1, 2
                spectra(l1) = spectra(l1) + dble(c(i, l1, m1))**2 &
                                          + aimag(c(i, l1, m1))**2

            end do

        end do

        spectra(l1) = spectra(l1) / dble(2*l+1)

    end do

end subroutine SHPowerSpectrumDensityC


subroutine SHCrossPowerSpectrumC(c1, c2, lmax, cspectra, exitstatus)
!------------------------------------------------------------------------------
!
!   This function will compute the dimensionless cross power spectrum
!   of the complex 4-pi normalized spherical harmonic coefficients c1(i, l, m)
!   and c2(1,l,m).
!
!   CrossPower(l) =  Sum_{m=0}^l ( C11lm*C21lm^* + C12lm*C22lm^* )
!
!   Calling Parameters
!
!       IN
!           c1      C   omplex spherical harmonic coefficients, dimensioned as
!                       (2, lmax+1, lmax+1).
!           c2          Complex spherical harmonic coefficients, dimensioned as
!                       (2, lmax+1, lmax+1).
!           lmax        Spherical harmonic degree to compute power.
!       OUT
!           cspectra    Array of length (lmax+1) containing the complex cross
!                       power spectra of c.
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
    implicit none
    
    complex*16, intent(in) :: c1(:,:,:), c2(:,:,:)
    integer, intent(in) :: lmax
    complex*16, intent(out) ::  cspectra(:)
    integer, intent(out), optional :: exitstatus
    integer i, m, l1, m1, l

    if (present(exitstatus)) exitstatus = 0

    if (size(c1(:,1,1)) < 2 .or. size(c1(1,:,1)) < lmax+1 &
            .or. size(c1(1,1,:)) < lmax+1) then
        print*, "Error --- SHCrossPowerSpectrumC"
        print*, "C1 must be dimensioned as (2, LMAX+1, LMAX+1) " // &
                "where lmax is", lmax
        print*, "Input array is dimensioned ", size(c1(:,1,1)), &
                size(c1(1,:,1)), size(c1(1,1,:))
        if (present(exitstatus)) then
            exitstatus = 1
            return
        else
            stop
        end if

    else if (size(c2(:,1,1)) < 2 .or. size(c2(1,:,1)) < lmax+1 &
            .or. size(c2(1,1,:)) < lmax+1) then
        print*, "Error --- SHCrossPowerSpectrumC"
        print*, "C2 must be dimensioned as (2, LMAX+1, LMAX+1)"
        print*, "Input array is dimensioned ", size(c2(:,1,1)), &
                size(c2(1,:,1)), size(c2(1,1,:))
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

    cspectra = 0.0d0

    do l = 0, lmax
        l1 = l + 1

        cspectra(l1) = c1(1, l1, 1) * dconjg(c2(1, l1, 1))

        do m = 1, l, 1
            m1 = m + 1

            do i = 1, 2
                cspectra(l1) = cspectra(l1) &
                               + c1(i, l1, m1) * dconjg(c2(i, l1, m1))

            end do

        end do

    end do

end subroutine SHCrossPowerSpectrumC


subroutine SHCrossPowerSpectrumDensityC(c1, c2, lmax, cspectra, exitstatus)
!------------------------------------------------------------------------------
!
!   This function will compute the dimensionless cross power spectrum
!   density of the complex 4-pi normalized spherical harmonic coefficients
!   c1(i, l, m) and c2(i,l,m).
!
!   CrossPower(l) =  Sum_{m=0}^l ( C11lm*C21lm^* + C12lm*C22lm^* ) / (2l+1)
!
!   Calling Parameters
!
!       IN
!           c1          Complex spherical harmonic coefficients, dimensioned as
!                       (2, lmax+1, lmax+1).
!           c2          Complex spherical harmonic coefficients, dimensioned as
!                       (2, lmax+1, lmax+1).
!           lmax        Spherical harmonic degree to compute power.
!
!       OUT
!           cspectra    Array of length (lmax+1) containing the complex cross
!                       power spectral density of c.
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
    implicit none

    complex*16, intent(in) :: c1(:,:,:), c2(:,:,:)
    integer, intent(in) :: lmax
    complex*16, intent(out) ::  cspectra(:)
    integer, intent(out), optional :: exitstatus
    integer i, m, l1, m1, l

    if (present(exitstatus)) exitstatus = 0

    if (size(c1(:,1,1)) < 2 .or. size(c1(1,:,1)) < lmax+1 &
            .or. size(c1(1,1,:)) < lmax+1) then
        print*, "Error --- SHCrossPowerSpectrumDensityC"
        print*, "C1 must be dimensioned as (2, LMAX+1, LMAX+1) " // &
                "where lmax is", lmax
        print*, "Input array is dimensioned ", size(c1(:,1,1)), &
                size(c1(1,:,1)), size(c1(1,1,:))
        if (present(exitstatus)) then
            exitstatus = 1
            return
        else
            stop
        end if

    else if (size(c2(:,1,1)) < 2 .or. size(c2(1,:,1)) < lmax+1 &
            .or. size(c2(1,1,:)) < lmax+1) then
        print*, "Error --- SHCrossPowerSpectrumDensityC"
        print*, "C2 must be dimensioned as (2, LMAX+1, LMAX+1)"
        print*, "Input array is dimensioned ", size(c2(:,1,1)), &
                size(c2(1,:,1)), size(c2(1,1,:))
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

    cspectra = 0.0d0

    do l = 0, lmax
        l1 = l + 1

        cspectra(l1) = c1(1, l1, 1) * dconjg(c2(1, l1, 1))

        do m = 1, l, 1
            m1 = m + 1

            do i = 1, 2
                cspectra(l1) = cspectra(l1) &
                               + c1(i, l1, m1) * dconjg(c2(i, l1, m1))

            end do

        end do

        cspectra(l1) = cspectra(l1) / dble(2*l+1)

    end do

end subroutine SHCrossPowerSpectrumDensityC
