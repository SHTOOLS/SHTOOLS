subroutine SHAdmitCorr(gilm, tilm, lmax, admit, corr, admit_error, exitstatus)
!------------------------------------------------------------------------------
!
!   This subroutine will compute the admittance and correlation between the two
!   input spherical harmonic coefficients. The optional "error" on the
!   admittance estimates assumes that the two fields are perfectly correlated,
!   which may not necessarily true for certain models (such as Forsyth's
!   loading model).
!
!   Assuming that the two input fields are related by
!
!   gilm = Z tilm + N
!
!   where N is noise, the admittance is calcalated as
!
!       Z = < gilm tilm > / < tilm tilm >,
!
!   the correlation is
!
!       correlation = < gilm tilm > / sqrt(< gilm gilm > < tilm tilm >).
!
!   Assuming that the two fields are perfectly correlated, and that the lack of
!   coherence is a result of noise, the uncertainty on the admittance can be
!   calculaed from
!
!       sigma^2 = ( < gilm gilm > / < tilm tilm > ) ( 1 - gamma^2) / (2 l).
!
!   < > signifies an average of over all m. Since only ratios are being used,
!   < > is implemented only as sums over all m (i.e., by the cross spectra of
!   the fields). Note that the correlation can possess values between -1 and 1,
!   in contrast to the "coherence" (or coherence-squared; the terminology is
!   ambiguous) which only possesses values between 0 and 1.
!
!   Note that in order for the magnitude of the admittance to be correct, the
!   input coefficients must be multiplied by all constants. For instance, in the
!   case of gravity-topography admittances, the potential coefficients must
!   first be multiplied by GM (l+1) (R/r)**(l+2) where R is the refernce radius
!   of the coefficients, and r is the radius at which the field is evaluated.
!
!   Calling Parameters
!
!       IN
!           gilm        First input spherical harmonic coefficients.
!           tilm        Second input spherical harmonic coefficients.
!           lmax        Maximum spherical harmonic degree of the input fields.
!
!       OUT
!           admit       Admittance spectra (length lmax+1).
!           corr        Correlation spectra (length lmax+1).
!
!       OPTIONAL
!           admit_error Addittance error spectra (length lmax+1)
!                       assuming that the two fields perfectly correlated in
!                       the abscence of noise.
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
    use SHTOOLS, only: SHPowerSpectrum, SHCrossPowerSpectrum
    use ftypes

    implicit none

    real(dp), intent(in) :: gilm(:,:,:), tilm(:,:,:)
    integer, intent(in) :: lmax
    real(dp), intent(out) :: admit(:), corr(:)
    real(dp), intent(out), optional :: admit_error(:)
    integer, intent(out), optional :: exitstatus
    real(dp) :: gt(lmax+1), gg(lmax+1), tt(lmax+1)
    integer :: l, l1

    if (present(exitstatus)) exitstatus = 0

    if (size(gilm(:,1,1)) < 2 .or. size(gilm(1,:,1)) < lmax+1 .or. &
            size(gilm(1,1,:)) < lmax+1) then
        print*, "Error --- SHAdmitCorr"
        print*, "GILM must be dimensioned as (2, LMAX+1, LMAX+1) where " // &
                "LMAX is ", lmax
        print*, "Input dimension is ", size(gilm(:,1,1)), size(gilm(1,:,1)), &
                size(gilm(1,1,:))
        if (present(exitstatus)) then
            exitstatus = 1
            return
        else
            stop
        end if

    else if (size(tilm(:,1,1)) < 2 .or. size(tilm(1,:,1)) < lmax+1 .or. &
            size(tilm(1,1,:)) < lmax+1) then
        print*, "Error --- SHAdmitCorr"
        print*, "TILM must be dimensioned as (2, LMAX+1, LMAX+1) " // &
                "where LMAX is ", lmax
        print*, "Input dimension is ", size(tilm(:,1,1)), size(tilm(1,:,1)), &
                size(tilm(1,1,:))
        if (present(exitstatus)) then
            exitstatus = 1
            return
        else
            stop
        end if

    else if (size(admit) < lmax+1) then
        print*, "Error --- SHAdmitCorr"
        print*, "ADMIT must be dimensioned as (LMAX+1) where LMAX is ", lmax
        print*, "Input dimension is ", size(admit)
        if (present(exitstatus)) then
            exitstatus = 1
            return
        else
            stop
        end if

    else if (size(corr) < lmax+1) then
        print*, "Error --- SHAdmitCorr"
        print*, "CORR must be dimensioned as (LMAX+1) where LMAX is ", lmax
        print*, "Input dimension is ", size(corr)
        if (present(exitstatus)) then
            exitstatus = 1
            return
        else
            stop
        end if

    end if
    
    if (present(admit_error)) then
        if (size(admit_error) < lmax+1) then
            print*, "Error --- SHAdmitCorr"
            print*, "ADMIT_ERROR must be dimensioned as (LMAX+1) " // &
                    "where LMAX is ", lmax
            print*, "Input dimension is ", size(admit_error)
            if (present(exitstatus)) then
                exitstatus = 1
                return
            else
                stop
            end if

        end if

    end if

    admit = 0.0_dp
    corr = 0.0_dp

    if (present(exitstatus)) then
        call SHCrossPowerSpectrum(gilm, tilm, lmax, gt, exitstatus=exitstatus)
        if (exitstatus /= 0) return
        call SHPowerSpectrum(gilm, lmax, gg, exitstatus=exitstatus)
        if (exitstatus /= 0) return
        call SHPowerSpectrum(tilm, lmax, tt, exitstatus=exitstatus)
        if (exitstatus /= 0) return
    else
        call SHCrossPowerSpectrum(gilm, tilm, lmax, gt) ! < gilm tilm >
        call SHPowerSpectrum(gilm, lmax, gg) ! < gilm gilm >
        call SHPowerSpectrum(tilm, lmax, tt) ! < tilm tilm >
    end if

    admit(1:lmax+1) = gt(1:lmax+1) / tt(1:lmax+1)
    corr(1:lmax+1) = gt(1:lmax+1) / sqrt( tt(1:lmax+1) * gg(1:lmax+1) )

    if (present(admit_error)) then
        admit_error = 0.0_dp
        do l = 1, lmax, 1
            l1 = l + 1
            admit_error(l1) = gg(l1)*(1.0_dp - corr(l1)**2) &
                              / ( tt(l1) * dble(2*l) )
        end do

        admit_error(1:lmax+1) = sqrt(admit_error(1:lmax+1))

    end if

end subroutine SHAdmitCorr
