subroutine PLegendreA_d1(p, dp, lmax, z, csphase, exitstatus)
!------------------------------------------------------------------------------
!
!   This function evalutates all of the unnormalized associated legendre
!   polynomials and their first derivatives up to degree lmax.
!
!   Calling Parameters
!
!       IN
!           lmax        Maximum spherical harmonic degree to compute.
!           z           [-1, 1], cos(colatitude) or sin(latitude).
!
!       OPTIONAL (IN)
!           csphase     1: Do not include the phase factor of (-1)^m
!                       -1: Apply the phase factor of (-1)^m.
!
!       OUT
!           p           A vector of all associated Legendgre polynomials
!                       evaluated at z up to lmax. The length must by greater
!                       or equal to (lmax+1)*(lmax+2)/2.
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
!   Notes:
!
!   1.  The integral of plm**2 over (-1,1) is 2 * (l+m)! / (2l+1) / (l-m)!.
!   2.  The index of the array p corresponds to l*(l+1)/2 + m + 1.
!   3.  The index of the array p corresponds to l*(l+1)/2 + m + 1.
!   4.  The derivative is evaluated with respecte to z, and NOT cos(colatitude)
!       or sin(latitude).
!   5.  Derivatives are calculated using the unnormalized identities
!           P'l,m = ( (l+m) Pl-1,m - l z Plm ) / (1-z**2)   (for l>m), and
!           P'll = - l z Pll / (1-z**2) (for l=m).
!   6.  The derivative is not defined at z=+-1 for all m>0, and is therefore not
!       calculated here.
!   7.  The default is to exlude the Condon-Shortley phase of (-1)^m.
!
!   Dependencies:   CSPHASE_DEFAULT
!
!   Copyright (c) 2016, SHTOOLS
!   All rights reserved.
!
!------------------------------------------------------------------------------
    use SHTOOLS, only: CSPHASE_DEFAULT

    implicit none

    integer, intent(in) :: lmax
    real*8, intent(out) :: p(:), dp(:)
    real*8, intent(in) :: z
    integer, intent(in), optional :: csphase
    integer, intent(out), optional :: exitstatus
    real*8 :: pm2, pm1, pmm, sinsq, sinsqr, fact, plm
    integer :: k, kstart, m, l, sdim
    integer*1 :: phase

    if (present(exitstatus)) exitstatus = 0

    sdim = (lmax+1)*(lmax+2)/2

    if (size(p) < sdim) then
        print*, "Error --- PLegendreA_d1"
        print*, "P must be dimensioned as (LMAX+1)*(LMAX+2)/2 " // &
                "where LMAX is ", lmax
        print*, "Input array is dimensioned ", size(p)
        if (present(exitstatus)) then
            exitstatus = 1
            return
        else
            stop
        endif

    else if (size(dp) < sdim) then
        print*, "Error --- PLegendreA_d1"
        print*, "DP must be dimensioned as (LMAX+1)*(LMAX+2)/2 " // &
                "where LMAX is ", lmax
        print*, "Input array is dimensioned ", size(dp)
        if (present(exitstatus)) then
            exitstatus = 1
            return
        else
            stop
        endif

    else if (lmax < 0) then
        print*, "Error --- PLegendreA_d1"
        print*, "LMAX must be greater than or equal to 0."
        print*, "Input value is ", lmax
        if (present(exitstatus)) then
            exitstatus = 2
            return
        else
            stop
        endif

    elseif(abs(z) > 1.0d0) then
        print*, "Error --- PLegendreA_d1"
        print*, "ABS(Z) must be less than or equal to 1."
        print*, "Input value is ", z
        if (present(exitstatus)) then
            exitstatus = 2
            return
        else
            stop
        endif

    else if (abs(z) == 1.0d0) then
        print*, "Error --- PLegendreA_d1"
        print*, "Derivative can not be calculated at Z = 1 or -1."
        print*, "Input value is ", z
        if (present(exitstatus)) then
            exitstatus = 2
            return
        else
            stop
        endif

    end if

    if (present(csphase)) then
        if (csphase == -1) then
            phase = -1

        else if (csphase == 1) then
            phase = 1

        else
            print*, "Error --- PLegendreA_d1"
            print*, "CSPHASE must be 1 (exclude) or -1 (include)."
            print*, "Input value is ", csphase
            if (present(exitstatus)) then
                exitstatus = 2
                return
            else
                stop
            endif

        end if

    else
        phase = CSPHASE_DEFAULT

    end if

    !--------------------------------------------------------------------------
    !   
    !   Calculate P(l,0)
    !
    !--------------------------------------------------------------------------

    sinsq = (1.0d0-z) * (1.0d0+z)
    sinsqr = sqrt(sinsq)

    pm2 = 1.d0
    p(1) = 1.d0
    dp(1) = 0.0d0

    if (lmax == 0) return

    pm1 = z
    p(2) = pm1
    dp(2) = 1.0d0

    k = 2

    do l = 2, lmax, 1
        k = k + l
        plm = ( z * (2*l-1) * pm1 - (l-1) * pm2 ) / dble(l)
        p(k) = plm
        dp(k) = l * (pm1 -z*plm) / sinsq
        pm2 = pm1
        pm1 = plm

    end do

    !--------------------------------------------------------------------------
    !
    !   Calculate P(m,m), P(m+1,m), and P(l,m)
    !
    !--------------------------------------------------------------------------
    pmm = 1.0d0
    fact = -1.0d0
    kstart = 1

    do m = 1, lmax - 1, 1

        ! Calculate P(m,m)
        kstart = kstart + m + 1
        fact = fact + 2.0d0
        pmm = phase * pmm * sinsqr * fact
        p(kstart) = pmm
        dp(kstart) = -m * z * pmm / sinsq
        pm2 = pmm

        ! Calculate P(m+1,m)
        k = kstart + m + 1
        pm1 = z * pmm * (2 * m + 1)
        p(k) = pm1
        dp(k) = ( (2*m+1) * pmm - (m+1) * z * pm1) / sinsq

        ! Calculate P(l,m)
        do l = m + 2, lmax, 1
            k = k + l
            plm  = ( z * (2*l-1) * pm1 - (l+m-1) * pm2 ) / dble(l-m)
            p(k) = plm
            dp(k) = ( (l+m) * pm1 - l * z * plm) / sinsq
            pm2 = pm1
            pm1 = plm
        end do

    end do

    ! P(lmax, lmax)
    kstart = kstart + m + 1
    fact = fact + 2.0d0
    pmm = phase * pmm * sinsqr * fact
    p(kstart) = pmm
    dp(kstart) = -lmax*z*pmm / sinsq

end subroutine PLegendreA_d1
