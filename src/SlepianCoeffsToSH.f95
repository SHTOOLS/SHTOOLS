subroutine SlepianCoeffsToSH(flm, falpha, galpha, lmax, nmax, exitstatus)
!------------------------------------------------------------------------------
!
!   This subroutine will compute the spherical harmonic coefficients of a
!   function FLM given the input Slepian functions GALPHA and associated
!   coefficients FALPHA. The Slepian functions are determined by a call to
!   SHReturnTapers and then SHRotateTapers, or SHReturnTapersMap. Each row of
!   GALPHA contains the (LMAX+1)**2 spherical harmonic coefficients of the
!   Slepian function ordered according to the subroutine SHCilmToVector.
!   The Slepian functions must be normalized to have unit power (that is the
!   sum of the coefficients squared is 1), and the spherical harmonic
!   coefficients are calculated as
!
!       f_lm = sum_{alpha}^{nmax} f_alpha g(alpha)_lm
!
!   Calling Parameters
!
!       IN
!           falpha  1D vector of dimension nmax containing the Slepian
!                   coefficients of the function FLM.
!           galpha  Matrix of dimension ( (LMAX+1)**2, nmax) containing the
!                   spherical harmonic coefficients of the Slepian functions.
!                   Each column corresponds to a Slepian function ordered from
!                   best to worst concentrated.
!           lmax    Maximum spherical harmonic degree of the Slepian functions.
!           nmax    Maximum number of Slepian coefficients used to construct
!                   the function.
!
!       OUT
!           flm     Spherical harmonic coefficients of the function FLM with
!                      dimension (2, LMAX+1, LMAX+1).
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
    use SHTOOLS, only: SHVectorToCilm

    implicit none
    real*8, intent(out) :: flm(:,:,:)
    real*8, intent(in) :: falpha(:), galpha(:,:)
    integer, intent(in) :: lmax, nmax
    integer, intent(out), optional :: exitstatus
    real*8, allocatable :: f(:)
    integer :: i, astat

    if (present(exitstatus)) exitstatus = 0

    if (size(flm(:,1,1)) < 2 .or. size(flm(1,:,1)) < lmax+1 .or. &
            size(flm(1,1,:)) < lmax+1) then
        print*, "Error --- SlepianCoeffsToSH"
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

    else if (size(falpha) < nmax) then
        print*, "Error --- SlepianCoeffsToSH"
        print*, "FALPHA must be dimensioned as (NMAX)."
        print*, "NMAX = ", nmax
        print*, "Dimension of FALPHA = ", size(falpha)
        if (present(exitstatus)) then
            exitstatus = 1
            return
        else
            stop
        end if

    else if (size(galpha(:,1)) < (lmax+1)**2 .or. &
            size(galpha(1,:)) < nmax) then
        print*, "Error --- SlepianCoeffsToSH"
        print*, "GALPHA must be dimensioned as ( (LMAX+1)**2, " // &
                "nmax ), where LMAX = ", lmax, "and NMAX = ", nmax
        print*, "Input array is dimensioned as ", size(galpha(:,1)), &
                size(galpha(1,:))
        if (present(exitstatus)) then
            exitstatus = 1
            return
        else
            stop
        end if

    end if

    allocate(f((lmax+1)**2), stat = astat)
    if (astat /= 0) then
        print*, "Error --- SlepianCoeffsToSH"
        print*, "Problem allocating array f", astat
        if (present(exitstatus)) then
            exitstatus = 3
            return
        else
            stop
        end if
    endif

    f = 0.0d0

    do i=1, nmax, 1
        f(1:(lmax+1)**2) = f(1:(lmax+1)**2) + &
                           falpha(i) * galpha(1:(lmax+1)**2, i)
    end do

    if (present(exitstatus)) then
        call SHVectorToCilm(f, flm, lmax, exitstatus=exitstatus)
        if (exitstatus /= 0) return
    else
        call SHVectorToCilm(f, flm, lmax)
    endif

    deallocate(f)

end subroutine SlepianCoeffsToSH