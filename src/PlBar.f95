subroutine PlBar(p, lmax, z, exitstatus)
!------------------------------------------------------------------------------
!
!   This function evalutates all of the "geophysical normalized legendre
!   polynomials up to degree lmax.
!
!   Calling Parameters
!
!       Out
!           p       A vector of all associated Legendgre polynomials evaluated
!                   at z up to lmax. The lenght must by greater or equal to
!                   (lmax+1).
!
!       IN
!           lmax    Maximum degree to compute.
!           z       cos(colatitude) or sin(latitude).
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
!   1.  The employed normalization is the "geophysical convention."
!   2.  The integral of PlBar**2 over all space on the sphere is 4 pi.
!       The integral of PlBar**2 over (-1,1) is 2.
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
        print*, "Error --- PlBar"
        print*, "P must be dimensioned as (LMAX+1) where LMAX is ", lmax
        print*, "Input array is dimensioned ", size(p)
        if (present(exitstatus)) then
            exitstatus = 1
            return
        else
            stop
        endif

    else if (lmax < 0) then 
        print*, "Error --- PlBar"
        print*, "LMAX must be greater than or equal to 0."
        print*, "Input value is ", lmax
        if (present(exitstatus)) then
            exitstatus = 2
            return
        else
            stop
        endif

    else if (abs(z) > 1.0d0) then
        print*, "Error --- PlBar"
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

    pm1 = sqrt(3.0d0) * z
    p(2) = pm1

    do l = 2, lmax, 1
        pl = ( sqrt(dble(2*l-1)) * z * pm1 - &
                (l-1) * pm2 / sqrt(dble(2*l-3)) ) * &
                sqrt(dble(2*l+1)) / dble(l)
        p(l+1) = pl
        pm2 = pm1
        pm1 = pl

    end do

end subroutine PlBar
