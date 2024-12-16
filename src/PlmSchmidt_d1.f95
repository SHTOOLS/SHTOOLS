subroutine PlmSchmidt_d1(p, dp1, lmax, z, csphase, cnorm, exitstatus)
!------------------------------------------------------------------------------
!
!   This function evalutates all of the normalized associated Legendre
!   functions up to degree lmax. The functions are initially scaled by
!   10^280 sin^m in order to minimize the effects of underflow at large m
!   near the poles (see Holmes and Featherstone 2002, J. Geodesy, 76, 279-299).
!   On a macOS system with a maximum allowable double precision value of
!   2.225073858507203E-308 the scaled portion of the algorithm will not overflow
!   for degrees less than or equal to 2800.
!
!   For each value of m, the rescaling factor is computed as
!   rescalem=rescalem*sin(theta), with the intial value of rescalem being equal
!   to 1/scalef (which is here equal to 10^280). This will gradually reduce this
!   huge number to a tiny number, and will ultimately underflow. In order to
!   prevent this underflow, when rescalem becomes less than 10^(-280), the
!   subsequent rescaling factors of sin(theta) will be directly applied to Plm,
!   and then this number will be multipled by the old value of rescalem.
!
!   Temporary variables are saved in an allocated array. In order to explicitly
!   deallocate this memory, call this routine with a spherical harmonic degree
!   of -1.
!
!   Calling Parameters
!
!       IN
!           lmax        Maximum spherical harmonic degree to compute.
!           z           cos(colatitude) or sin(latitude).
!
!       OPTIONAL (IN)
!           csphase     1: Do not include the phase factor of (-1)^m (default).
!                       -1: Apply the phase factor of (-1)^m.
!           cnorm       0: Use real normalization.
!                       1: Use complex normalization.
!
!       OUT
!           p           A vector of all associated Legendgre polynomials
!                       evaluated at z up to lmax. The lenght must by greater
!                       or equal to (lmax+1)*(lmax+2)/2.
!           dp1         A vector of all first derivatives of the normalized
!                       Legendgre polynomials evaluated at z up to lmax with
!                       dimension (lmax+1).
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
!   1.  The employed normalization is the "Schmidt semi-normalized convention."
!       The integral of (PlmSchmidt*cos(m theta))**2 or
!       (PlmSchmidt*sin (m theta))**2 over all space is 4 pi/(2l+1).
!   2.  The integral of PlmSchmidt**2 over (-1,1) is 2 * (2 - delta(0,m)) /
!       (2l+1). If CNORM = 1, then this is equal to 2/(2l+1).
!   3.  The index of the array p corresponds to l*(l+1)/2 + m + 1. As such
!       the array p should be dimensioned as (lmax+1)*(lmax+2)/2 in the
!       calling routine.
!   4.  Derivatives are calculated using the unnormalized identities
!           P'l,m = ( (l+m) Pl-1,m - l z Plm ) / (1-z**2)   (for l>m), and
!           P'll = - l z Pll / (1-z**2) (for l=m).
!   5.  The derivative is not defined at z=+-1 for all m>0, and is therefore not
!       calculated here.
!   6.  The default is to exlude the Condon-Shortley phase of (-1)^m.
!
!   Copyright (c) 2005-2019, SHTOOLS
!   All rights reserved.
!
!------------------------------------------------------------------------------
    use ftypes

    implicit none

    integer(int32), intent(in) :: lmax
    real(dp), intent(out) :: p(:), dp1(:)
    real(dp), intent(in) :: z
    integer(int32), intent(in), optional :: csphase, cnorm
    integer(int32), intent(out), optional :: exitstatus
    real(dp) :: pm2, pm1, pmm, plm, rescalem, u, scalef
    real(dp), save, allocatable :: f1(:), f2(:), sqr(:)
    integer(int32) :: k, kstart, m, l, sdim, astat(3)
    integer(int32), save :: lmax_old = 0
    integer(int32) :: phase

!$OMP    threadprivate(f1, f2, sqr, lmax_old)

    if (present(exitstatus)) exitstatus = 0

    if (lmax == -1) then
        if (allocated (sqr)) deallocate (sqr)
        if (allocated (f1)) deallocate (f1)
        if (allocated (f2)) deallocate (f2)
        lmax_old = 0
        return
    end if

    sdim = (lmax+1)*(lmax+2)/2

    if (size(p) < sdim) then
        print*, "Error --- PlmSchmidt_d1"
        print*, "P must be dimensioned as (LMAX+1)*(LMAX+2)/2 " // &
                "where LMAX is ", lmax
        print*, "Input array is dimensioned ", size(p)
        if (present(exitstatus)) then
            exitstatus = 1
            return
        else
            stop
        end if

    else if (size(dp1) < sdim) then
        print*, "Error --- PlmSchmidt_d1"
        print*, "DP1 must be dimensioned as (LMAX+1)*(LMAX+2)/2 " // &
                "where LMAX is ", lmax
        print*, "Input array is dimensioned ", size(dp1)
        if (present(exitstatus)) then
            exitstatus = 1
            return
        else
            stop
        end if

    else if (lmax < 0) then
        print*, "Error --- PlmSchmidt_d1"
        print*, "LMAX must be greater than or equal to 0."
        print*, "Input value is ", lmax
        if (present(exitstatus)) then
            exitstatus = 2
            return
        else
            stop
        end if

    else if (abs(z) > 1.0_dp) then
        print*, "Error --- PlmSchmidt_d1"
        print*, "ABS(Z) must be less than or equal to 1."
        print*, "Input value is ", z
        if (present(exitstatus)) then
            exitstatus = 2
            return
        else
            stop
        end if

    else if (abs(z) == 1.0_dp) then
        print*, "Error --- PlmSchmidt_d1"
        print*, "Derivative can not be calculated at Z = 1 or -1."
        print*, "Input value is ", z
        if (present(exitstatus)) then
            exitstatus = 2
            return
        else
            stop
        end if

    end if

    if (present(csphase)) then
        if (csphase == -1) then
            phase = -1

        else if (csphase == 1) then
                phase = 1

        else
            print*, "Error --- PlmSchmidt_d1"
            print*, "CSPHASE must be 1 (exclude) or -1 (include)."
            print*, "Input value is ", csphase
            if (present(exitstatus)) then
                exitstatus = 2
                return
            else
                stop
            end if

        end if

    else
        phase = 1

    end if

    scalef = 1.0e-280_dp

    if (lmax > lmax_old) then
        if (allocated (sqr)) deallocate (sqr)
        if (allocated (f1)) deallocate (f1)
        if (allocated (f2)) deallocate (f2)

        allocate (sqr(2*lmax+1), stat=astat(1))
        allocate (f1((lmax+1)*(lmax+2)/2), stat=astat(2))
        allocate (f2((lmax+1)*(lmax+2)/2), stat=astat(3))

        if (astat(1) /= 0 .or. astat(2) /= 0 .or. astat(3) /= 0) then
            print*, "Error --- PlmSchmidt_d1"
            print*, "Problem allocating arrays SQR, F1 and F2", astat(1), &
                    astat(2), astat(3)
            if (present(exitstatus)) then
                exitstatus = 3
                return
            else
                stop
            end if
        end if

        !----------------------------------------------------------------------
        !
        !   Precompute square roots of integers that are used several times.
        !
        !----------------------------------------------------------------------
        do l = 1, 2 * lmax+1
            sqr(l) = sqrt(dble(l))
        end do

        !----------------------------------------------------------------------
        !
        !   Precompute multiplicative factors used in recursion relationships
        !       PlmSchmidt(l,m) = x*f1(l,m)*PlmSchmidt(l-1,m) -
        !       PlmSchmidt(l-2,m)*f2(l,m)
        !       k = l*(l+1)/2 + m + 1
        !   Note that prefactors are not used for the case when m=l and m=l-1,
        !   as a different recursion is used for these two values.
        !
        !----------------------------------------------------------------------
        k = 3

        do l = 2, lmax, 1
            k = k + 1
            f1(k) = dble(2*l-1) / dble(l)
            f2(k) = dble(l-1) / dble(l)

            do m = 1, l-2
                k = k + 1
                f1(k) = dble(2*l-1) / sqr(l+m) / sqr(l-m)
                f2(k) = sqr(l-m-1) * sqr(l+m-1) / sqr(l+m) / sqr(l-m)
            end do

            k = k + 2

        end do

        lmax_old = lmax

    end if

    !--------------------------------------------------------------------------
    !
    !   Calculate P(l,0). These are not scaled.
    !
    !--------------------------------------------------------------------------
    u = sqrt((1.0_dp - z) * (1.0_dp + z)) ! sin(theta)

    pm2 = 1.0_dp
    p(1) = 1.0_dp
    dp1(1) = 0.0_dp

    if (lmax == 0) return

    pm1 = z
    p(2) = pm1
    dp1(2) = 1.0_dp

    k = 2

    do l = 2, lmax, 1
        k = k + l
        plm = f1(k) * z * pm1 - f2(k) * pm2
        p(k) = plm
        dp1(k) = l * (pm1 - z * plm) / u**2
        pm2 = pm1
        pm1 = plm
    end do

    !--------------------------------------------------------------------------
    !
    !   Calculate P(m,m), P(m+1,m), and P(l,m)
    !
    !--------------------------------------------------------------------------
    if (present(cnorm)) then
        if (cnorm == 1) then
            pmm = scalef

        else
            pmm = sqr(2) * scalef

        end if

    else
        pmm = sqr(2) * scalef

    end if

    rescalem = 1.0_dp / scalef
    kstart = 1

    do m = 1, lmax - 1, 1

        rescalem = rescalem * u

        ! Calculate P(m,m)
        kstart = kstart + m + 1

        pmm = phase * pmm * sqr(2*m+1) / sqr(2*m)
        p(kstart) = pmm * rescalem / sqr(2*m+1)
        dp1(kstart) = -m * z * p(kstart) / u**2
        pm2 = pmm / sqr(2*m+1)

        ! Calculate P(m+1,m)
        k = kstart + m + 1
        pm1 = z * sqr(2*m+1) * pm2
        p(k) = pm1 * rescalem
        dp1(k) = (p(k-m-1) * sqr(2*m+1) - z * (m+1) * p(k)) / u**2

        ! Calculate P(l,m)
        do l = m + 2, lmax, 1
            k = k + l
            plm = z * f1(k) * pm1 - f2(k) * pm2
            p(k) = plm * rescalem
            dp1(k) = ( sqr(l+m) * sqr(l-m) * p(k-l) - l * z * p(k) ) / u**2
            pm2 = pm1
            pm1 = plm
        end do

    end do

    ! Calculate P(lmax,lmax)
    rescalem = rescalem * u

    kstart = kstart + m + 1
    pmm = phase * pmm  / sqr(2*lmax)
    p(kstart) = pmm*rescalem
    dp1(kstart) = -lmax * z * p(kstart) / u**2

end subroutine PlmSchmidt_d1
