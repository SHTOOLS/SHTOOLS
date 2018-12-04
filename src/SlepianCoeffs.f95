subroutine SlepianCoeffs(salpha, galpha, flm, lmax, nalpha, exitstatus)
!------------------------------------------------------------------------------
!
!   This subroutine will compute the Slepian coefficients of an input function
!   FLM given the Slepian functions GALPHA. The Slepian functions are
!   determined by a call to SHReturnTapers and then SHRotateTapers, or
!   SHReturnTapersMap. Each row of GALPHA contains the (LMAX+1)**2 spherical
!   harmonic coefficients of the Slepian function ordered according to the
!   subroutine SHCilmToVector. The Slepian functions must be normalized to
!   have unit power (that is the sum of the coefficients squared is 1), and the
!   Slepian coefficients are calculated as
!
!       s(alpha) = sum_{lm}^{lmax} f_lm g_lm(alpha)
!
!   Calling Parameters
!
!       IN
!           galpha  Matrix of dimension ( (LMAX+1)**2, NALPHA) containing the
!                   spherical harmonic coefficients of the Slepian functions.
!                   Each column corresponds to a Slepian function ordered from
!                   best to worst concentrated.
!           flm     Input spherical harmonic coefficients with dimension
!                      (2, LMAX+1, LMAX+1).
!           lmax    Maximum spherical harmonic degree of the Slepian functions.
!           nalpha  Maximum number of Slepian coefficients to return.
!
!       OUT
!           salpha  1D vector of dimension NALPHA containing the Slepian
!                   coefficients of the function FLM.
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
!   Copyright (c) 2018, SHTOOLS
!   All rights reserved.
!
!------------------------------------------------------------------------------
    use SHTOOLS, only: SHCilmToVector

    implicit none
    real*8, intent(out) :: salpha(:)
    real*8, intent(in) :: galpha(:,:), flm(:,:,:)
    integer, intent(in) :: lmax, nalpha
    integer, intent(out), optional :: exitstatus
    real*8, allocatable :: f(:)
    integer :: i, astat

    if (present(exitstatus)) exitstatus = 0

    if (size(salpha) < nalpha) then
        print*, "Error --- SlepianCoeffs"
        print*, "SALPHA must be dimensioned as (NALPHA)."
        print*, "NALPHA = ", NALPHA
        print*, "Dimension of SALPHA = ", size(salpha)
        if (present(exitstatus)) then
            exitstatus = 1
            return
        else
            stop
        end if

    else if (size(galpha(:,1)) < (lmax+1)**2 .or. &
            size(galpha(1,:)) < nalpha) then
        print*, "Error --- SlepianCoeffs"
        print*, "GALPHA must be dimensioned as ( (LMAX+1)**2, " // &
                "NALPHA ), where LMAX = ", lmax, "and NALPHA = ", nalpha
        print*, "Input array is dimensioned as ", size(galpha(:,1)), &
                size(galpha(1,:))
        if (present(exitstatus)) then
            exitstatus = 1
            return
        else
            stop
        end if

    else if (size(flm(:,1,1)) < 2 .or. size(flm(1,:,1)) < lmax+1 .or. &
            size(flm(1,1,:)) < lmax+1) then
        print*, "Error --- SlepianCoeffs"
        print*, "FLM must be dimensioned as (2, LMAX+1, LMAX + 1)."
        print*, "LMAX = ", lmax
        print*, "Dimension of FLM = ", size(flm(:,1,1)), size(flm(1,:,1)), &
            size(flm(1,1,:))
        if (present(exitstatus)) then
            exitstatus = 1
            return
        else
            stop
        end if

    end if

    allocate(f((lmax+1)**2), stat = astat)
    if (astat /= 0) then
        print*, "Error --- SlepianCoeffs"
        print*, "Problem allocating array f", astat
        if (present(exitstatus)) then
            exitstatus = 3
            return
        else
            stop
        end if
    end if

    if (present(exitstatus)) then
        call SHCilmToVector(flm, f, lmax, exitstatus=exitstatus)
        if (exitstatus /= 0) return
    else
        call SHCilmToVector(flm, f, lmax)
    end if

    salpha = 0.0d0

    do i=1, nalpha, 1
        salpha(i) = dot_product(f(1:(lmax+1)**2), galpha(1:(lmax+1)**2, i))
    end do

    deallocate(f)

end subroutine SlepianCoeffs