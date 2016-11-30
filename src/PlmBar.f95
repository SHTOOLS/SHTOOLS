subroutine PlmBar(p, lmax, z, csphase, cnorm, exitstatus)
!------------------------------------------------------------------------------
!
!   This function evalutates all of the normalized associated Legendre
!   functions up to degree lmax. The functions are initially scaled by
!   10^280 sin^m in order to minimize the effects of underflow at large m
!   near the poles (see Holmes and Featherstone 2002, J. Geodesy, 76, 279-299).
!   On a Mac OSX system with a maximum allowable double precision value of
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
!   Temporary variables in saved in an allocated array. In order to explicitly
!   deallocate this memory, call this routine with a spherical harmonic degree
!   of -1.
!
!   Calling Parameters
!
!       OUT
!           p           A vector of all associated Legendgre polynomials.
!                       evaluated at z up to lmax. The length must by greater
!                       or equal to (lmax+1)*(lmax+2)/2.
!
!       OPTIONAL (IN)
!           csphase     1: Do not include the phase factor of (-1)^m,
!                       -1: Apply the phase factor of (-1)^m.
!           cnorm       0: Use real normalization.
!                       1: Use complex normalization.
!
!       IN
!           lmax        Maximum spherical harmonic degree to compute.
!           z           cos(colatitude) or sin(latitude).
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
!   1.  The employed normalization is the "geophysical convention." The integral
!       of (plm*cos(m theta))**2 or (plm*sin (m theta))**2 over all space is
!       4 pi.
!   2.  The integral of plm**2 over (-1,1) is 2 * (2 - delta(0,m)). If CNORM=1,
!       then this is equal to 2.
!   3.  The index of the array p corresponds to l*(l+1)/2 + m + 1. As such
!       the array p should be dimensioned as (lmax+1)*(lmax+2)/2 in the
!       calling routine.
!   4.  The default is to exclude the Condon-Shortley phase of (-1)^m.
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
    real*8, intent(out) :: p(:)
    real*8, intent(in) :: z
    integer, intent(in), optional :: csphase, cnorm
    integer, intent(out), optional :: exitstatus
    real*8 :: pmm, rescalem, u, scalef
    real*8, save, allocatable :: f1(:), f2(:), sqr(:)
    integer :: k, kstart, m, l, astat(3)
    integer, save :: lmax_old = 0
    integer*1 :: phase

!$OMP    threadprivate(f1, f2, sqr, lmax_old)

    if (present(exitstatus)) exitstatus = 0

    if (lmax == -1) then
        if (allocated (sqr)) deallocate (sqr)
        if (allocated (f1)) deallocate (f1)
        if (allocated (f2)) deallocate (f2)
        lmax_old = 0
        return
    end if

    if (size(p) < (lmax+1)*(lmax+2)/2) then
        print*, "Error --- PlmBar"
        print*, "P must be dimensioned as (LMAX+1)*(LMAX+2)/2 " // &
                "where LMAX is ", lmax
        print*, "Input array is dimensioned ", size(p)
        if (present(exitstatus)) then
            exitstatus = 1
            return
        else
            stop
        end if

    else if (lmax < 0) then
        print*, "Error --- PlmBar"
        print*, "LMAX must be greater than or equal to 0."
        print*, "Input value is ", lmax
        if (present(exitstatus)) then
            exitstatus = 2
            return
        else
            stop
        end if

    else if(abs(z) > 1.0d0) then
        print*, "Error --- PlmBar"
        print*, "ABS(Z) must be less than or equal to 1."
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
        elseif (csphase == 1) then
                phase = 1
        else
            print*, "Error --- PlmBar"
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
            phase = CSPHASE_DEFAULT

    end if

    scalef = 1.0d-280

    if (lmax > lmax_old) then
        if (allocated (sqr)) deallocate (sqr)
        if (allocated (f1)) deallocate (f1)
        if (allocated (f2)) deallocate (f2)

        allocate (sqr(2*lmax+1), stat=astat(1))
        allocate (f1((lmax+1)*(lmax+2)/2), stat=astat(2))
        allocate (f2((lmax+1)*(lmax+2)/2), stat=astat(3))

        if (astat(1) /= 0 .or. astat(2) /= 0 .or. astat(3) /= 0) then
            print*, "Error --- PlmBar"
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
        !       Plmbar(l,m) = x*f1(l,m)*Plmbar(l-1,m) - Plmbar(l-2,m)*f2(l,m)
        !       k = l*(l+1)/2 + m + 1
        !   Note that prefactors are not used for the case when m=l and m=l-1,
        !   as a different recursion is used for these two values.
        !
        !----------------------------------------------------------------------
        k = 3

        do l = 2, lmax, 1
            k = k + 1
            f1(k) = sqr(2*l-1) * sqr(2*l+1) / dble(l)
            f2(k) = dble(l-1) * sqr(2*l+1) / sqr(2*l-3) / dble(l)

            do m = 1, l-2
                k = k + 1
                f1(k) = sqr(2*l+1) * sqr(2*l-1) / sqr(l+m) / sqr(l-m)
                f2(k) = sqr(2*l+1) * sqr(l-m-1) * sqr(l+m-1) &
                             / sqr(2*l-3) / sqr(l+m) / sqr(l-m)
            end do

            k = k + 2

        enddo

        lmax_old = lmax

    end if

    !--------------------------------------------------------------------------
    !
    !   Calculate P(l,0). These are not scaled.
    !
    !--------------------------------------------------------------------------
    u = sqrt((1.0d0 - z) * (1.0d0 + z)) ! sin(theta)

    p(1) = 1.0d0

    if (lmax == 0) return

    p(2) = sqr(3) * z

    k = 2

    do l = 2, lmax, 1
        k = k + l
        p(k) = f1(k) * z * p(k-l) - f2(k) * p(k-2*l+1)
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

    rescalem = 1.0d0/scalef
    kstart = 1

    do m = 1, lmax - 1, 1
        rescalem = rescalem * u

        ! Calculate P(m,m)
        kstart = kstart + m + 1

        pmm = phase * pmm * sqr(2*m+1) / sqr(2*m)
        p(kstart) = pmm

        ! Calculate P(m+1,m)
        k = kstart + m + 1
        p(k) = z * sqr(2*m+3) * pmm

        ! Calculate P(l,m)
        do l = m + 2, lmax, 1
            k = k + l
                p(k) = z * f1(k) * p(k-l) - f2(k) * p(k-2*l+1)
                p(k-2*l+1) = p(k-2*l+1) * rescalem
        end do

        p(k) = p(k) * rescalem
        p(k-lmax) = p(k-lmax) * rescalem

        end do

        ! Calculate P(lmax,lmax)
        rescalem = rescalem * u

        kstart = kstart + m + 1
        p(kstart) = phase * pmm * sqr(2*lmax+1) / sqr(2*lmax) * rescalem

end subroutine PlmBar
