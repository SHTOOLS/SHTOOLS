subroutine SHSlepianVar(l, galpha, galpha_order, lmax, kmax, Sff, variance, &
                        exitstatus)
!------------------------------------------------------------------------------
!
!   Given the first KMAX Slepian functions in the matrix GALPHA, and an input
!   global power spectrum Sff, this subroutine will compute the theoretical
!   variance of the power of a function expanded in the first KMAX Slepian
!   functions at spherical harmonic degree l. This routine works using only the
!   Slepian functions of the spherical-cap concentration problem obtained from
!   SHReturnTapers.
!
!   Calling Parameters
!
!       IN
!           l               The single spherical harmonic degree to compute.
!           galpha          An array of (lmax+1, kmax) Slepian functions
!                           (arranged in columns) obtained from
!                           SHReturnTapers.
!           galpha_order    An array of dimension kmax containing the real
!                           angular order of the functions in galpha, as
!                           obtained from SHReturnTapers.
!           lmax            Maximum spherical harmonic degree of the
!                           bandlimited Slepian functions.
!           kmax            Maximum number of Slepian functions to be used in
!                           making the spectral estimate.
!           Sff             Known global power spectrum of the function.
!
!       OUT
!           variance        The variance of the function expanded in Slepian
!                           functions at spherical harmonic degree l.
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

    real(dp), intent(in) :: galpha(:,:), Sff(:)
    real(dp), intent(out) :: variance
    integer, intent(in) :: l, lmax, kmax, galpha_order(:)
    integer, intent(out), optional :: exitstatus
    integer :: m, lp, alpha, beta
    real(dp) :: fmm

    if (present(exitstatus)) exitstatus = 0

    variance = 0.0_dp

    if (size(Sff) < lmax + 1) then
        print*, "Error --- SHSlepianVar"
        print*, "Sff must be dimensioned (LMAX+1) where LMAX is ", lmax
        print*, "Input array is dimensioned ", size(Sff)
        if (present(exitstatus)) then
            exitstatus = 1
            return
        else
            stop
        end if

    else if (size(galpha(:,1)) < lmax+1 .or. size(galpha(1,:)) < kmax) then
        print*, "Error ---  SHSlepianVar"
        print*, "GALPHA must be dimensioned as (LMAX+1, KMAX) where " // &
                "LMAX and KMAX are ", lmax, kmax
        print*, "Input array is dimensioned ", size(galpha(:,1)), &
                size(galpha(1,:))
        if (present(exitstatus)) then
            exitstatus = 1
            return
        else
            stop
        end if

    else if (size(galpha_order) < kmax) then
        print*, "Error ---  SHSlepianVar"
        print*, "GALPHA_ORDER must be dimensioned as (KMAX) where KMAX is ", &
                kmax
        print*, "Input array is dimensioned ", size(galpha_order)
        if (present(exitstatus)) then
            exitstatus = 1
            return
        else
            stop
        end if

    end if

    !--------------------------------------------------------------------------
    !
    !   Calculate the variance,
    !
    !     var S_{\tilde{f}}(l) =
    !         \sum_{m=-l}^l 2 < \tilde{f}_{lm} \tilde{f}_{lm} >^2
    !
    !     < \tilde{f}_{lm} \tilde{f}_{lm} > =
    !         \sum_{l'=0}^L \frac{S_{f}(l')}{(2l'+1)} \sum_{\alpha=1}^N
    !         \sum_{\beta=1}^N g^{(\alpha)}_{l'm} g^{(\alpha)}_{lm}
    !         g^{(\beta)}_{l'm} g^{(\beta)}_{lm}.
    !
    !--------------------------------------------------------------------------
    do m = -l, l, 1

        fmm = 0.0_dp

        do lp = 0, lmax

            do alpha = 1, kmax
                if (galpha_order(alpha) /= m) cycle

                do beta = 1, kmax
                    if (galpha_order(beta) /= m) cycle

                    fmm = fmm + Sff(lp+1) / dble(2*lp+1) * &
                        galpha(lp+1, alpha) * galpha(l+1, alpha) * &
                        galpha(lp+1, beta) * galpha(l+1, beta)

                end do

            end do

        end do

        variance = variance + 2.0_dp * fmm**2

    end do

end subroutine SHSlepianVar
