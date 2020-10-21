subroutine SHMTCouplingMatrix(Mmt, lmax, tapers_power, lwin, K, taper_wt, &
                              exitstatus)
!------------------------------------------------------------------------------
!
!   This routine returns the multitaper coupling matrix, which relates the
!   global input spectrum to the expectation of the localized multitaper
!   spectrum.
!
!   < S_{Phi Phi}^(mt) > = M^(mt) S_{ff}
!
!   This is given by eqs 4.5 and 4.6 of Wieczorek and Simons (2007).
!   The input tapers_power is a matrix containing the power spectra (in
!   columns) of all the localization windows. When using spherical cap
!   localization windows, this is simply
!
!   tapers_power[:,:] = tapers[:,:]**2
!
!   Note that this routine returns the "full" coupling matrix of dimension
!   (lmax+lwin+1, lmax+1). When multiplied by a global input power spectrum
!   with bandwidth lmax, it returns the output power spectrum with a bandwidth
!   of lmax+lwin. In doing so, it is implicitly assumed that input power
!   spectrum is exactly zero for all degrees greater than lmax. If this is not
!   the case, the ouput power spectrum should be considered valid only
!   for the degrees up to and including lmax-lwin.
!
!   Copyright (c) 2005-2019, SHTOOLS
!   All rights reserved.
!
!------------------------------------------------------------------------------
    use SHTOOLS, only: wigner3j
    use ftypes

    implicit none
    real(dp), intent(out) :: Mmt(:,:)
    real(dp), intent(in) :: tapers_power(:,:)
    real(dp), intent(in), optional :: taper_wt(:)
    integer(int32), intent(in) :: lmax, K, lwin
    integer(int32), intent(out), optional :: exitstatus
    real(dp) :: w3j(lwin+2*lmax+1), sum1
    integer(int32) :: i, j, l, wmin, wmax

    if (present(exitstatus)) exitstatus = 0

    if  (size(Mmt(:,1)) < lmax+lwin+1 .or. size(Mmt(1,:)) < lmax+1) then
        print*, "Error --- SHMTCouplingMatrix"
        print*, "MMT must be dimensioned as (LMAX+LWIN+1, LMAX+1) where "// &
                "LMAX and LWIN are ", lmax, lwin
        print*, "Input array is dimensioned as ", size(Mmt(:,1)), size(Mmt(1,:))
        if (present(exitstatus)) then
            exitstatus = 1
            return
        else
            stop
        end if

    end if

    if (size(tapers_power(:,1)) < lwin+1 .or. size(tapers_power(1,:)) < K) then
        print*, "Error --- SHMTCouplingMatrix"
        print*, "TAPERS_POWER must be dimensioned as (LWIN+1, K) where LWIN "// &
                "and K are ", lwin, k
        print*, "Input array is dimensioned as ", size(tapers_power(:,1)), &
            size(tapers_power(1,:))
        if (present(exitstatus)) then
            exitstatus = 1
            return
        else
            stop
        end if

    end if

    if (present(taper_wt)) then
        if (size(taper_wt) < K) then
            print*, "Error --- SHMTCouplingMatrix"
            print*, "TAPER_WT must be dimensioned as (K) where K is ", k
            print*, "Input array is dimensioned as ", size(taper_wt)
            if (present(exitstatus)) then
                exitstatus = 1
                return
            else
                stop
            end if
        end if

    end if

    !--------------------------------------------------------------------------
    !
    !   Compute coupling matrix, M^(mt)
    !
    !--------------------------------------------------------------------------

    if (present(taper_wt)) then

        do i=0, lmax+lwin
            do j=0, lmax
                if (present(exitstatus)) then
                    call Wigner3j(w3j, wmin, wmax, i, j, 0, 0, 0, &
                                  exitstatus = exitstatus)
                    if (exitstatus /= 0) return
                else
                    call Wigner3j(w3j, wmin, wmax, i, j, 0, 0, 0)
                end if
                sum1 = 0.0_dp

                do l = wmin, min(wmax,lwin), 2
                    sum1 = sum1 + dot_product(taper_wt(1:K), &
                           tapers_power(l+1,1:K)) * w3j(l-wmin+1)**2
                end do

                Mmt(i+1,j+1) = sum1 * dble(2*i+1)

            end do
        end do

    else

        do i = 0, lmax+lwin
            do j = 0, lmax
                if (present(exitstatus)) then
                    call Wigner3j(w3j, wmin, wmax, i, j, 0, 0, 0, &
                                  exitstatus = exitstatus)
                    if (exitstatus /= 0) return
                else
                    call Wigner3j(w3j, wmin, wmax, i, j, 0, 0, 0)
                end if
                sum1 = 0.0_dp

                do l = wmin, min(wmax,lwin), 2
                    sum1 = sum1 + sum(tapers_power(l+1,1:K)) * w3j(l-wmin+1)**2
                end do

                Mmt(i+1,j+1) = sum1 * dble(2*i+1) / dble(K)

            end do
        end do

    end if

end subroutine SHMTCouplingMatrix
