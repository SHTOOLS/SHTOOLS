subroutine SlepianCoeffs(falpha, galpha, film, lmax, nmax, exitstatus)
!------------------------------------------------------------------------------
!
!   This subroutine will compute the Slepian coefficients of an input function
!   FIiILM given the Slepian functions GALPHA. The Slepian functions are
!   determined by a call to SHReturnTapers and then SHRotateTapers, or
!   SHReturnTapersMap. Each row of GALPHA contains the (LMAX+1)**2 spherical
!   harmonic coefficients of the Slepian function ordered according to the
!   subroutine SHCilmToVector. The Slepian functions must be normalized to
!   have unit power (that is the sum of the coefficients squared is 1), and the
!   Slepian coefficients are calculated as
!
!       f_alpha = sum_{lm}^{lmax} f_lm g(alpha)_lm
!
!   Calling Parameters
!
!       IN
!           galpha  Matrix of dimension ( (LMAX+1)**2, nmax) containing the
!                   spherical harmonic coefficients of the Slepian functions.
!                   Each column corresponds to a Slepian function ordered from
!                   best to worst concentrated.
!           film     Input spherical harmonic coefficients with dimension
!                      (2, LMAX+1, LMAX+1).
!           lmax    Maximum spherical harmonic degree of the Slepian functions.
!           nmax    Maximum number of Slepian coefficients to return.
!
!       OUT
!           falpha  1D vector of dimension nmax containing the Slepian
!                   coefficients of the function FILM.
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
    use SHTOOLS, only: SHCilmToVector
    use ftypes

    implicit none

    real(dp), intent(out) :: falpha(:)
    real(dp), intent(in) :: galpha(:,:), film(:,:,:)
    integer, intent(in) :: lmax, nmax
    integer, intent(out), optional :: exitstatus
    real(dp), allocatable :: f(:)
    integer :: i, astat

    if (present(exitstatus)) exitstatus = 0

    if (size(falpha) < nmax) then
        print*, "Error --- SlepianCoeffs"
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
        print*, "Error --- SlepianCoeffs"
        print*, "GALPHA must be dimensioned as ( (LMAX+1)**2, " // &
                "NMAX ), where LMAX = ", lmax, "and NMAX = ", nmax
        print*, "Input array is dimensioned as ", size(galpha(:,1)), &
                size(galpha(1,:))
        if (present(exitstatus)) then
            exitstatus = 1
            return
        else
            stop
        end if

    else if (size(film(:,1,1)) < 2 .or. size(film(1,:,1)) < lmax+1 .or. &
            size(film(1,1,:)) < lmax+1) then
        print*, "Error --- SlepianCoeffs"
        print*, "FILM must be dimensioned as (2, LMAX+1, LMAX + 1)."
        print*, "LMAX = ", lmax
        print*, "Dimension of FILM = ", size(film(:,1,1)), size(film(1,:,1)), &
            size(film(1,1,:))
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
        call SHCilmToVector(film, f, lmax, exitstatus=exitstatus)
        if (exitstatus /= 0) return
    else
        call SHCilmToVector(film, f, lmax)
    end if

    falpha = 0.0_dp

    do i=1, nmax, 1
        falpha(i) = dot_product(f(1:(lmax+1)**2), galpha(1:(lmax+1)**2, i))
    end do

    deallocate(f)

end subroutine SlepianCoeffs