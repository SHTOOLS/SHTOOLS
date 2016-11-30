subroutine PLegendre(p, lmax, z, exitstatus)
!------------------------------------------------------------------------------
!
!   This subroutine evalutates all of the unnormalized Legendre polynomials
!   up to degree lmax.
!
!   Calling Parameters
!
!       Out
!           p       A vector of all unnormalized Legendgre polynomials
!                   evaluated at z up to lmax. The lenght must by greater or
!                   equal to (lmax+1).
!
!       IN
!           lmax    Maximum degree to compute.
!           z       Value within [-1, 1], cos(colatitude) or sin(latitude).
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
!   1.  The integral of Pl**2 over (-1,1) is 2/(2l+1).
!   2.  Values are calculated accoring to the following recursion scheme:
!           P_0(z) = 1.0, P_1(z) = z, and 
!           P_l(z) = (2l-1) * z * P_{l-1}(z) / l - (l-1) * P_{l-2}(z) / l
!
!   Dependencies:   None
!
!   Copyright (c) 2016, SHTOOLS
!   All rights reserved.
!
!------------------------------------------------------------------------------
    implicit none

    integer, intent(in) :: lmax
    real*8, intent(out) :: p(:)
    real*8, intent(in) :: z
    integer, intent(out), optional :: exitstatus
    real*8 :: pm2, pm1, pl
    integer :: l

    if (present(exitstatus)) exitstatus = 0

    if (size(p) < lmax+1) then
        print*, "Error --- PlegendreL"
        print*, "P must be dimensioned as (LMAX+1) where LMAX is ", lmax
        print*, "Input array is dimensioned ", size(p)
        if (present(exitstatus)) then
            exitstatus = 1
            return
        else
            stop
        endif

    else if (lmax < 0) then 
        print*, "Error --- PlegendreL"
        print*, "LMAX must be greater than or equal to 0."
        print*, "Input value is ", lmax
        if (present(exitstatus)) then
            exitstatus = 2
            return
        else
            stop
        endif

    else if (abs(z) > 1.0d0) then
        print*, "Error --- PlegendreL"
        print*, "ABS(Z) must be less than or equal to 1."
        print*, "Input value is ", z
        if (present(exitstatus)) then
            exitstatus = 2
            return
        else
            stop
        endif

    end if

    pm2 = 1.d0
    p(1) = 1.d0

    pm1 = z
    p(2) = pm1

    do l = 2, lmax, 1
        pl = ((2*l-1) * z * pm1 - (l-1) * pm2) / dble(l)
        p(l+1) = pl
        pm2 = pm1
        pm1 = pl

    end do

end subroutine PLegendre
