function SHMagPowerL(cilm, a, r, l)
!------------------------------------------------------------------------------
!
!   This function will compute the dimensionless power for
!   degree l of the magnetic field evaluated at radius r,
!   given the Schmidt seminormalized spherical harmonic
!   coefficients cilm(i, l, m) of the magnetic potential evaluated
!   at radius a.
!
!   power(l) = (a/r)^(2l+4) (l+1) sum_{m=0}^l ( c1lm**2 + c2lm**2 )
!
!   Calling Parameters
!
!       cilm   Schmidt-seminormalized spherical harmonic coefficients of the
!              magnetic potential evaluated at radius a, dimensioned as
!              (2, lmax+1, lmax+1).
!       a      Reference radius of the magnetic potential spherical harmonic
!              coefficients.
!       r      Radius to evaluate the magnetic field.
!       l      Spherical harmonic degree to compute power.
!
!   Copyright (c) 2016, SHTOOLS
!   All rights reserved.
!
!------------------------------------------------------------------------------
    use ftypes

    implicit none

    real(dp) :: SHMagPowerL
    real(dp), intent(in) :: cilm(:,:,:)
    real(dp), intent(in) :: a, r
    integer, intent(in) :: l
    integer i, m, l1, m1

    l1 = l + 1

    if (size(cilm(:,1,1)) < 2 .or. size(cilm(1,:,1)) < l1 .or. &
        size(cilm(1,1,:)) < l1) then
        print*, "Error --- SHMagPowerL"
        print*, "CILM must be dimensioned as (2, L+1, L+1) where L is ", l
        print*, "Input array is dimensioned ", size(cilm(:,1,1)), &
                                               size(cilm(1,:,1)), &
                                               size(cilm(1,1,:))
        stop

    end if

    SHMagPowerL = 0.0_dp

    do m = 0, l
        m1 = m + 1
        do i = 1, 2
            SHMagPowerL = SHMagPowerL + cilm(i, l1, m1)**2
        end do
    end do

    SHMagPowerL = SHMagPowerL * dble(l+1) * (a / r)**(2*l+4)

end function SHMagPowerL


subroutine SHMagPowerSpectrum(cilm, a, r, lmax, spectra, exitstatus)
!------------------------------------------------------------------------------
!
!   This function will compute the dimensionless power spectrum
!   of the magnetic field evaluated at radius r,
!   given the Schmidt seminormalized spherical harmonic
!   coefficients cilm(i, l, m) of the magnetic potential evaluated
!   at radius a.
!
!   power(l) = (a/r)^(2l+4) (l+1) sum_{m=0}^l ( c1lm**2 + c2lm**2 )
!
!   Calling Parameters
!
!       IN
!           cilm    Schmidt seminormalized spherical harmonic coefficients of
!                   the magnetic potential evaluated at radius a, dimensioned
!                   as (2, lmax+1, lmax+1).
!           a       Reference radius of the magnetic potential spherical
!                   harmonic coefficients.
!           r       Radius to evaluate the magnetic field.
!           lmax    Maximum spherical harmonic degree to compute power of.
!
!       OUT
!           spectra Array of length (lmax+1) containing the power spectra of c.
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

    real(dp), intent(in) :: cilm(:,:,:)
    real(dp), intent(in) :: a, r
    integer, intent(in) :: lmax
    real(dp), intent(out) ::spectra(:)
    integer, intent(out), optional :: exitstatus
    integer i, m, l1, m1, l

    if (present(exitstatus)) exitstatus = 0

    if (size(cilm(:,1,1)) < 2 .or. size(cilm(1,:,1)) < lmax+1 &
        .or. size(cilm(1,1,:)) < lmax+1) then
        print*, "Error --- SHMagPowerSpectrum"
        print*, "CILM must be dimensioned as (2, LMAX+1, LMAX+1) " // &
                "where LMAX is ", lmax
        print*, "Input array is dimensioned ", size(cilm(:,1,1)), &
                                               size(cilm(1,:,1)), &
                                               size(cilm(1,1,:))
        if (present(exitstatus)) then
            exitstatus = 1
            return
        else
            stop
        end if

    else if (size(spectra) < lmax+1) then
        print*, "Error --- SHMagPowerSpectrum"
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

        do m = 0, l
            m1 = m + 1

            do i = 1, 2
                spectra(l1) = spectra(l1) + cilm(i, l1, m1)**2
            end do
        end do

        spectra(l1) = spectra(l1) * dble(l+1) * (a / r)**(2*l+4)

    end do

end subroutine SHMagPowerSpectrum
