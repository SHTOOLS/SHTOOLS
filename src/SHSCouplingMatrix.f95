subroutine SHSCouplingMatrix(kij, galpha, lmax, nmax, exitstatus)
!------------------------------------------------------------------------------
!
!   SHSCouplingMatrix returns the spherical harmonic coupling matrix that
!   relates the expectation of the function expressed in a subset of the best-
!   localized Slepian functions to the expectation of the global power spectrum
!   (assumed to be stationary). The Slepian functions are determined by a call
!   to either (1) `SHReturnTapers` and then `SHRotateTapers`, or
!   (2) `SHReturnTapersMap`. Each row of galpha contains the (lmax+1)**2
!   spherical harmonic coefficients of a Slepian function that can be unpacked
!   using `SHVectorToCilm`. The Slepian functions must be normalized to have
!   unit power (that is the sum of the coefficients squared is 1).
!
!   The relationship between the global and localized power spectra is given by
!
!       < S_{\tilde{f}}(l) > = \sum_{l'=0}^lmax K_{ll'} S_{f}(l')
!
!   where S_{\tilde{f}} is the expectation of the power spectrum at degree l of
!   the function expressed in Slepian functions, S_{f}(l') is the expectation
!   of the global power spectrum, and < ... > is the expectation operator.
!   The coupling matrix is given explicitly by
!
!   K_{ij} = \frac{1}{2l'+1} Sum_{m=-l}^l Sum_{m'=-l'}^l'
!            ( Sum_{alpha=1}^nmax g_{l'm'}(alpha) g_{lm}(alpha) )**2
!
!   Calling Parameters
!
!       IN
!           galpha  Matrix of dimension ( (LMAX+1)**2, nmax) containing the
!                   spherical harmonic coefficients of the Slepian functions.
!                   Each column corresponds to a Slepian function ordered from
!                   best to worst concentrated.
!           lmax    Maximum spherical harmonic degree of the Slepian functions.
!           nmax    Maximum number of Slepian coefficients used to construct
!                   the function.
!
!       OUT
!           kij     Spherical harmonic coupling matrix of dimension
!                   ((LMAX+1)**2, nmax).
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

    real(dp), intent(out) :: kij(:,:)
    real(dp), intent(in) :: galpha(:,:)
    integer(int32), intent(in) :: lmax, nmax
    integer(int32), intent(out), optional :: exitstatus
    integer(int32) :: l, lp, m, mp, alpha, ind, indp
    real(dp) :: temp

    if (present(exitstatus)) exitstatus = 0

    if (size(kij(:,1)) < (lmax+1) .or. size(kij(1,:)) < (lmax+1)) then
        print*, "Error --- SHSCouplingMatrix"
        print*, "KIJ must be dimensioned as (LMAX+1, LMAX+1)."
        print*, "LMAX = ", lmax
        print*, "Dimension of KIJ = ", size(kij(:,1)), size(kij(1,:))
        if (present(exitstatus)) then
            exitstatus = 1
            return
        else
            stop
        end if

    else if (size(galpha(:,1)) < (lmax+1)**2 .or. &
            size(galpha(1,:)) < nmax) then
        print*, "Error --- SHSCouplingMatrix"
        print*, "GALPHA must be dimensioned as ( (LMAX+1)**2, NMAX )." 
        print*, "LMAX = ", lmax
        print*, "NMAX = ", nmax
        print*, "Input array is dimensioned as ", size(galpha(:,1)), &
                size(galpha(1,:))
        if (present(exitstatus)) then
            exitstatus = 1
            return
        else
            stop
        end if

    end if

    kij = 0.0_dp

    ! Calculate k(l, lp) * (2lp +1 ), which is symmetric
    do l = 0, lmax
        do lp = l, lmax, 1
            do m = -l, l, 1

                if (m >= 0) then
                    ind = l**2+m+1
                else
                    ind = l**2+l+abs(m)+1
                end if

                do mp = -lp, lp, 1

                    if (mp >= 0) then
                        indp = lp**2+mp+1
                    else
                        indp = lp**2+lp+abs(mp)+1
                    end if

                    temp = 0.0_dp
                    do alpha = 1, nmax, 1

                       temp = temp + galpha(indp, alpha) * galpha(ind, alpha)

                    end do

                    kij(l+1, lp+1) = kij(l+1, lp+1) + temp**2

                end do

            end do

        end do

    end do

    ! Fill lower half of matrix
    do l = 1, lmax, 1

        do lp = 0, l-1, 1

            kij(l+1, lp+1) = kij(lp+1, l+1)

        end do

    end do

    ! Multiply by 1/(2lp+1)
    do lp = 0, lmax

        kij(1:lmax+1, lp+1) = kij(1:lmax+1, lp+1) / dble(2*lp+1)

    end do

end subroutine SHSCouplingMatrix