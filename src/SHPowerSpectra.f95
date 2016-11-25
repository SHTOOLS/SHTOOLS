real*8 function SHPowerL(c, l)
!------------------------------------------------------------------------------
!
!   This function will compute the dimensionless power at
!   degree l corresponding to the 4pi spherical harmonic coefficients \
!   c(i, l, m)
!
!   Power = Sum_{m=0}^l ( C1lm**2 + C2lm**2 )
!
!   Calling Parameters
!
!       c       Spherical harmonic coefficients, dimensioned as
!               (2, lmax+1, lmax+1).
!       l       Spherical harmonic degree to compute power.
!
!
!   Copyright (c) 2016, SHTOOLS
!   All rights reserved.
!
!------------------------------------------------------------------------------
    implicit none

    real*8, intent(in) :: c(:,:,:)
    integer, intent(in) :: l
    integer i, m, l1, m1

    l1 = l + 1

    if (size(c(:,1,1)) < 2 .or. size(c(1,:,1)) < l1 .or. &
            size(c(1,1,:)) < l1) then
        print*, "Error --- SHPowerL"
        print*, "C must be dimensioned as (2, L+1, L+1) where L is ", l
        print*, "Input array is dimensioned ", size(c(:,1,1)), &
                size(c(1,:,1)), size(c(1,1,:)) 
        stop

    end if

    SHPowerL = c(1, l1, 1)**2  ! m=0 term

    do m = 1, l, 1
        m1 = m + 1

        do i = 1, 2
            SHPowerL = SHPowerL + c(i, l1, m1)**2

        end do

    end do

end function SHPowerL


real*8 function SHPowerDensityL(c, l)
!------------------------------------------------------------------------------
!
!   This function will compute the dimensionless power per coefficient
!   (density) at degree l of the 4pi spherical harmonic coefficients c(i, l, m)
!
!   PowerSpectralDensity = Sum_{m=0}^l ( C1lm**2 + C2lm**2 ) / (2l+1)
!
!   Calling Parameters
!
!       c       Spherical harmonic coefficients, dimensioned as
!               (2, lmax+1, lmax+1).
!       l       Spherical harmonic degree to compute power.
!
!   Copyright (c) 2016, SHTOOLS
!   All rights reserved.
!
!------------------------------------------------------------------------------
    implicit none

    real*8, intent(in) :: c(:,:,:)
    integer, intent(in) :: l
    integer i, m, l1, m1

    l1 = l + 1

    if (size(c(:,1,1)) < 2 .or. size(c(1,:,1)) < l1 &
            .or. size(c(1,1,:)) < l1) then
        print*, "Error --- SHPowerDensityL"
        print*, "C must be dimensioned as (2, L+1, L+1) where L is ", l
        print*, "Input array is dimensioned ", size(c(:,1,1)), &
                size(c(1,:,1)), size(c(1,1,:)) 
        stop

    end if

    SHPowerDensityL = c(1, l1, 1)**2
    ! m=0 term

    do m = 1, l, 1
        m1 = m + 1

        do i = 1, 2
            SHPowerDensityL = SHPowerDensityL + c(i, l1, m1)**2

        end do

    end do

    SHPowerDensityL = SHPowerDensityL/dble(2*l+1)

end function SHPowerDensityL


real*8 function SHCrossPowerL(c1, c2, l)
!------------------------------------------------------------------------------
!
!   This function will compute the dimensionless cross power at
!   degree l for the 4pi spherical harmonic coefficients c1(i, l, m)
!   and c2(i,l,m).
!
!   CrossPower =  Sum_{m=0}^l ( A1lm*B1lm + A2lm*B2lm )
!
!   Calling Parameters
!
!       c1      Spherical harmonic coefficients, dimensioned as
!               (2, lmax+1, lmax+1).
!       c2      Spherical harmonic coefficients, dimensioned as
!               (2, lmax+1, lmax+1).
!       l       Spherical harmonic degree to compute power.
!
!   Copyright (c) 2016, SHTOOLS
!   All rights reserved.
!
!------------------------------------------------------------------------------
    implicit none

    real*8, intent(in) :: c1(:,:,:), c2(:,:,:)
    integer, intent(in) :: l
    integer i, m, l1, m1

    l1 = l + 1

    if (size(c1(:,1,1)) < 2 .or. size(c1(1,:,1)) < l1 &
            .or. size(c1(1,1,:)) < l1) then
        print*, "Error --- SHCrossPowerL"
        print*, "C1 must be dimensioned as (2, L+1, L+1) where L is ", l
        print*, "Input array is dimensioned ", size(c1(:,1,1)), &
                size(c1(1,:,1)), size(c1(1,1,:)) 
        stop

    else if (size(c2(:,1,1)) < 2 .or. size(c2(1,:,1)) < l1 &
            .or. size(c2(1,1,:)) < l1) then
        print*, "Error --- SHCrossPowerL"
        print*, "C2 must be dimensioned as (2, L+1, L+1) where L is ", l
        print*, "Input array is dimensioned ", size(c2(:,1,1)), &
                size(c2(1,:,1)), size(c2(1,1,:)) 
        stop

    end if

    SHCrossPowerL = c1(1, l1, 1) * c2(1,l1,1)

    do m = 1, l, 1
        m1 = m + 1

        do i = 1, 2
            SHCrossPowerL = SHCrossPowerL + c1(i,l1,m1) * c2(i,l1,m1)

        end do

    end do

end function SHCrossPowerL


real*8 function SHCrossPowerDensityL(c1, c2, l)
!------------------------------------------------------------------------------
!
!   This function will compute the dimensionless cross power
!   (density) at degree l of the 4pi spherical harmonic coefficients 
!   c1(i, l, m) and c2(i,l,m).
!
!   CrossPower =  Sum_{m=0}^l ( A1lm*B1lm + A2lm*B2lm ) / (2l+1)
!
!   Calling Parameters
!
!       c1      Spherical harmonic coefficients, dimensioned as
!               (2, lmax+1, lmax+1).
!       c2      Spherical harmonic coefficients, dimensioned as
!               (2, lmax+1, lmax+1).
!       l       Spherical harmonic degree to compute power.
!
!   Copyright (c) 2016, SHTOOLS
!   All rights reserved.
!
!------------------------------------------------------------------------------
    implicit none

    real*8, intent(in) :: c1(:,:,:), c2(:,:,:)
    integer, intent(in) :: l
    integer i, m, l1, m1

    l1 = l + 1

    if (size(c1(:,1,1)) < 2 .or. size(c1(1,:,1)) < l1 &
            .or. size(c1(1,1,:)) < l1) then
        print*, "Error --- SHCrossPowerDensityL"
        print*, "C1 must be dimensioned as (2, L+1, L+1) where L is ", l
        print*, "Input array is dimensioned ", size(c1(:,1,1)), &
                size(c1(1,:,1)), size(c1(1,1,:)) 
        stop

    else if (size(c2(:,1,1)) < 2 .or. size(c2(1,:,1)) < l1 &
            .or. size(c2(1,1,:)) < l1) then
        print*, "Error --- SHCrossPowerDensityL"
        print*, "C2 must be dimensioned as (2, L+1, L+1) where L is ", l
        print*, "Input array is dimensioned ", size(c2(:,1,1)), &
                size(c2(1,:,1)), size(c2(1,1,:)) 
        stop

    end if

    SHCrossPowerDensityL = c1(1, l1, 1) * c2(1,l1,1)

    do m = 1, l, 1
        m1 = m + 1
        do i = 1, 2
            SHCrossPowerDensityL = SHCrossPowerDensityL &
                                   + c1(i,l1,m1)*c2(i,l1,m1)

        end do

    end do

    SHCrossPowerDensityL = SHCrossPowerDensityL / dble(2*l+1)

end function SHCrossPowerDensityL


subroutine SHPowerSpectrum(c, lmax, spectra, exitstatus)
!------------------------------------------------------------------------------
!
!   This function will compute the dimensionless power spectrum
!   of the 4pi spherical harmonic coefficients c(i, l, m).
!
!   Power(l) = Sum_{m=0}^l ( C1lm**2 + C2lm**2 )
!
!   Calling Parameters
!
!       IN
!           c           Spherical harmonic coefficients, dimensioned as
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

    real*8, intent(in) :: c(:,:,:)
    integer, intent(in) :: lmax
    real*8, intent(out) ::  spectra(:)
    integer, intent(out), optional :: exitstatus
    integer i, m, l1, m1, l

    if (present(exitstatus)) exitstatus = 0

    if (size(c(:,1,1)) < 2 .or. size(c(1,:,1)) < lmax+1 &
            .or. size(c(1,1,:)) < lmax+1) then
        print*, "Error --- SHPowerSpectrum"
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

    spectra = 0.0d0

    do l = 0, lmax
        l1 = l + 1

        spectra(l1) = c(1, l1, 1)**2

        do m = 1, l, 1
            m1 = m + 1

            do i = 1, 2
                spectra(l1) = spectra(l1) + c(i, l1, m1)**2

            end do

        end do

    end do

end subroutine SHPowerSpectrum


subroutine SHPowerSpectrumDensity(c, lmax, spectra, exitstatus)
!------------------------------------------------------------------------------
!
!   This function will compute the dimensionless power spectrum
!   density of the 4pi spherical harmonic coefficients c(i, l, m).
!
!   PowerSpectralDensity = Sum_{m=0}^l ( C1lm**2 + C2lm**2 ) / (2l+1)
!
!   Calling Parameters
!
!       IN
!           c           Spherical harmonic coefficients, dimensioned as
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
!------------------------------------------------------------------------------
    implicit none

    real*8, intent(in) :: c(:,:,:)
    integer, intent(in) :: lmax
    real*8, intent(out) ::  spectra(:)
    integer, intent(out), optional :: exitstatus
    integer i, m, l1, m1, l

    if (present(exitstatus)) exitstatus = 0

    if (size(c(:,1,1)) < 2 .or. size(c(1,:,1)) < lmax+1 &
            .or. size(c(1,1,:)) < lmax+1) then
        print*, "Error --- SHPowerSpectrumDensity"
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

    spectra = 0.0d0

    do l = 0, lmax
        l1 = l + 1

        spectra(l1) = c(1, l1, 1)**2

        do m = 1, l, 1
            m1 = m + 1

            do i = 1, 2
                spectra(l1) = spectra(l1) + c(i, l1, m1)**2

            end do

        end do

        spectra(l1) = spectra(l1) / dble(2*l+1)

    enddo

end subroutine SHPowerSpectrumDensity


subroutine SHCrossPowerSpectrum(c1, c2, lmax, cspectra, exitstatus)
!-------------------------------------------------------------------------------
!
!   This function will compute the dimensionless cross power spectrum
!   of the 4pi spherical harmonic coefficients c1(i, l, m) and c2(1,l,m).
!
!   CrossPower(l) =  Sum_{m=0}^l ( A1lm*B1lm + A2lm*B2lm )
!
!   Calling Parameters
!
!       IN
!           c1          Spherical harmonic coefficients, dimensioned as
!                       (2, lmax+1, lmax+1).
!           c2          Spherical harmonic coefficients, dimensioned as
!                       (2, lmax+1, lmax+1).
!           lmax        Spherical harmonic degree to compute power.
!       OUT
!           cspectra    Array of length (lmax+1) containing the cross power
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
!-------------------------------------------------------------------------------
    implicit none

    real*8, intent(in) :: c1(:,:,:), c2(:,:,:)
    integer, intent(in) :: lmax
    real*8, intent(out) ::  cspectra(:)
    integer, intent(out), optional :: exitstatus
    integer i, m, l1, m1, l

    if (present(exitstatus)) exitstatus = 0

    if (size(c1(:,1,1)) < 2 .or. size(c1(1,:,1)) < lmax+1 &
            .or. size(c1(1,1,:)) < lmax+1) then
        print*, "Error --- SHCrossPowerSpectrum"
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
        print*, "Error --- SHCrossPowerSpectrum"
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

    cspectra = 0.0d0

    do l = 0, lmax
        l1 = l + 1

        cspectra(l1) = c1(1, l1, 1) * c2(1, l1, 1)

        do m = 1, l, 1
            m1 = m + 1

            do i = 1, 2
                cspectra(l1) = cspectra(l1) + c1(i,l1,m1) * c2(i, l1, m1)

            end do

        end do

    end do

end subroutine SHCrossPowerSpectrum


subroutine SHCrossPowerSpectrumDensity(c1, c2, lmax, cspectra, exitstatus)
!------------------------------------------------------------------------------
!
!   This function will compute the dimensionless cross power spectrum
!   density of the 4pi spherical harmonic coefficients c1(i, l, m) and
!   c2(i,l,m).
!
!   CrossPower =  Sum_{m=0}^l ( A1lm*B1lm + A2lm*B2lm ) / (2l+1)
!
!   Calling Parameters
!
!       IN
!           c1          Spherical harmonic coefficients, dimensioned as
!                       (2, lmax+1, lmax+1).
!           c2          Spherical harmonic coefficients, dimensioned as
!                       (2, lmax+1, lmax+1).
!           lmax        Spherical harmonic degree to compute power.
!
!       OUT
!           cspectra    Array of length (lmax+1) containing the cross power
!                       spectral density of c.
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

    real*8, intent(in) :: c1(:,:,:), c2(:,:,:)
    integer, intent(in) :: lmax
    real*8, intent(out) ::  cspectra(:)
    integer, intent(out), optional :: exitstatus
    integer i, m, l1, m1, l

    if (size(c1(:,1,1)) < 2 .or. size(c1(1,:,1)) < lmax+1 &
            .or. size(c1(1,1,:)) < lmax+1) then
        print*, "Error --- SHCrossPowerSpectrumDensity"
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
        print*, "Error --- SHCrossPowerSpectrumDensity"
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

    cspectra = 0.0d0
    
    do l = 0, lmax
        l1 = l + 1
        
        cspectra(l1) = c1(1, l1, 1) * c2(1, l1, 1)

        do m = 1, l, 1
            m1 = m + 1

            do i = 1, 2
                cspectra(l1) = cspectra(l1) + c1(i,l1,m1) * c2(i, l1, m1)

            end do

        end do

        cspectra(l1) = cspectra(l1) / dble(2*l+1)

    end do

end subroutine SHCrossPowerSpectrumDensity
