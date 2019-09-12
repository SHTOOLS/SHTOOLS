subroutine SHMTVar(l, tapers, taper_order, lwin, kmax, Sff, variance, &
                   taper_wt, unweighted_covar, nocross, exitstatus)
!------------------------------------------------------------------------------
!
!   Given the first Kmax tapers of a matrix TAPERS, and an input global power
!   spectrum Sff, this subroutine will compute the theoretical variance of a
!   multitaper spectral estimate at a given degree l, and for a given set of
!   optional input taper weights. This routine works only using the tapers of
!   the spherical-cap concentration problem.
!
!   Calling Parameters
!
!       IN
!           l               The single spherical harmonic degree to compute.
!           tapers          An array of (lwin+1, kmax) tapers (arranged in
!                           columns).
!           taper_order     An array of dimension kmax containing the REAL
!                           angular order of the tapers.
!           lwin            Maximum spherical harmonic degree of the
!                           bandlimited tapers.
!           kmax            Maximum number of tapers to be used in making the
!                           spectral estimate.
!           Sff             Known global power spectrum of the unwindowed
!                           field.
!
!       OUT
!           variance        The variance at spherical harmonic degree l.
!
!       OPTIONAL (IN)
!           taper_wt        Vector of dimension (kmax) containing the numerical
!                           values of the weights to be applied to each
!                           spectral estimate.
!           nocross         If present and equal to 1, then the off-diagonal
!                           terms of the covariance matrix will be assumed to
!                           be zero.
!
!       OPTIONAL (OUT)
!           unweighted_covar   Unweighted covariance matrix, Fij
!           exitstatus  If present, instead of executing a STOP when an error
!                       is encountered, the variable exitstatus will be
!                       returned describing the error.
!                       0 = No errors;
!                       1 = Improper dimensions of input array;
!                       2 = Improper bounds for input variable;
!                       3 = Error allocating memory;
!                       4 = File IO error.
!
!   Copyright (c) 2019, SHTOOLS
!   All rights reserved.
!
!------------------------------------------------------------------------------
    use SHTOOLS, only: SHSjkPG
    use ftypes

    implicit none

    real(dp), intent(in) :: tapers(:,:), Sff(:)
    real(dp), intent(out) :: variance
    integer, intent(in) :: l, lwin, kmax, taper_order(:)
    real(dp), intent(in), optional :: taper_wt(:)
    real(dp), intent(out), optional :: unweighted_covar(:, :)
    integer, intent(in), optional :: nocross
    integer, intent(out), optional :: exitstatus
    real(dp) :: Fij(kmax, kmax)
    integer :: i, j, m, mp
    complex(dp) :: temp1

    if (present(exitstatus)) exitstatus = 0

    if (size(Sff) < l+lwin + 1) then
        print*, "Error --- SHMTVar"
        print*, "Sff must be dimensioned (L+LWIN+1) where L and LWIN " // &
                "are ", l, lwin
        print*, "Input array is dimensioned ", size(Sff)
        if (present(exitstatus)) then
            exitstatus = 1
            return
        else
            stop
        end if

    else if (size(tapers(:,1)) < lwin+1 .or. size(tapers(1,:)) < kmax) then
        print*, "Error --- SHMTVar"
        print*, "TAPERS must be dimensioned as (LWIN+1, KMAX) where " // &
                "LWIN and KMAX are ", lwin, kmax
        print*, "Input array is dimensioned ", size(tapers(:,1)), &
                size(tapers(1,:))
        if (present(exitstatus)) then
            exitstatus = 1
            return
        else
            stop
        end if

    else if (size(taper_order) < kmax) then
        print*, "Error --- SHMTVar"
        print*, "TAPER_ORDER must be dimensioned as (KMAX) where KMAX is ", &
                kmax
        print*, "Input array is dimensioned ", size(taper_order)
        if (present(exitstatus)) then
            exitstatus = 1
            return
        else
            stop
        end if

    end if

    if (present(taper_wt)) then
        if (size(taper_wt(:)) < kmax) then
            print*, "Error --- SHMTVar"
            print*, "TAPER_WT must be dimensioned (KMAX) " // &
                    "where KMAX is ", kmax
            if (present(exitstatus)) then
                exitstatus = 1
                return
            else
                stop
            end if

        end if

    end if

    if (present(unweighted_covar)) then
        if (size(unweighted_covar(:,1)) < kmax .or. &
                size(unweighted_covar(1,:)) < kmax) then
            print*, "Error --- SHMTVar"
            print*, "UNWEIGHTED_COVAR must be dimensioned (KMAX, KMAX) " // &
                    "where KMAX is ", kmax
            print*, "Input array is dimensioned ", &
                    size(unweighted_covar(:,1)), size(unweighted_covar(1,:))
            if (present(exitstatus)) then
                exitstatus = 1
                return
            else
                stop
            end if

        end if

    end if

    if (present(nocross)) then
        if (nocross /= 1 .and. nocross /= 0) then
            print*, "Error --- SHMTVar"
            print*, "NOCROSS must be either 0 (use all elements of " // &
                    "covariance matrix) or"
            print*, "1 (set off-diagonal elements of " // &
                    "covariance matrix to zero)."
            print*, "Input value is ", nocross
            if (present(exitstatus)) then
                exitstatus = 2
                return
            else
                stop
            end if

        end if

    end if

    Fij = 0.0_dp
    variance = 0.0_dp

    !--------------------------------------------------------------------------
    !
    !   Calculate matrix Fij
    !
    !--------------------------------------------------------------------------
    if (present(nocross)) then
        if (nocross == 1) then
            do i = 1, kmax, 1
                do m = -l, l
                    do mp = -l, l
                        temp1 = SHSjkPG(Sff, l, m, mp, tapers(1:lwin+1,i), &
                                        tapers(1:lwin+1,i), taper_order(i), &
                                        taper_order(i), lwin, 1)
                        Fij(i,i) = Fij(i,i) + 2.0_dp * dble(temp1*dconjg(temp1))

                    end do

                end do

            end do

        else
            do i = 1, kmax, 1
                do j = i, kmax, 1
                    do m = -l, l
                        do mp = -l, l
                            temp1 = SHSjkPG(Sff, l, m, mp, &
                                            tapers(1:lwin+1,i), &
                                            tapers(1:lwin+1,j), &
                                            taper_order(i), taper_order(j), &
                                            lwin, 1)
                            Fij(i,j) = Fij(i,j) + &
                                       2.0_dp * dble(temp1*dconjg(temp1))

                        end do

                    end do

                    if (i /= j) Fij(j,i) = Fij(i,j)

                end do

            end do

        endif

    else
        do i = 1, kmax, 1
            do j = i, kmax, 1
                do m = -l, l
                    do mp = -l, l
                        temp1 = SHSjkPG(Sff, l, m, mp, tapers(1:lwin+1,i), &
                                        tapers(1:lwin+1,j), taper_order(i), &
                                        taper_order(j), lwin, 1)
                        Fij(i,j) = Fij(i,j) + 2.0_dp * dble(temp1*dconjg(temp1))

                    end do

                end do

                if (i /= j) Fij(j,i) = Fij(i,j)

            end do

        end do

    endif

    if (present(unweighted_covar)) then
        unweighted_covar = 0.0_dp
        unweighted_covar(1:kmax, 1:kmax) = Fij(1:kmax,1:kmax)

    end if

    !----------------------------------------------------------------------
    !
    !   Calculate variances either using input weights, or equal weights.
    !
    !----------------------------------------------------------------------
    if (present(taper_wt)) then
        do i = 1, kmax
            do j = 1, kmax
                variance = variance + Fij(i, j) * taper_wt(i) * taper_wt(j)

            end do
        end do

    else
        variance = sum(Fij(1:kmax,1:kmax)) / dble(kmax)**2

    end if


end subroutine SHMTVar
