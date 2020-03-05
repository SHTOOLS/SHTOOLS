subroutine PlSchmidt(p, lmax, z, exitstatus)
!------------------------------------------------------------------------------
!
!   This function evalutates all of the Schmidt normalized legendre
!   polynomials up to degree lmax.
!
!   Calling Parameters
!
!       Out
!           p       A vector of all Schmidt normalized Legendgre polynomials
!                   evaluated at z up to lmax. The lenght must by greater or
!                   equal to (lmax+1).
!
!       IN
!           lmax    Maximum degree to compute.
!           z       [-1, 1], cos(colatitude) or sin(latitude).
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
!   1.  The integral of Pl**2 over (-1,1) is 2 * / (2l+1).
!   2.  The integral of Pl**2 over all space is 4 pi / (2l+1).
!
!   Copyright (c) 2005-2019, SHTOOLS
!   All rights reserved.
!
!------------------------------------------------------------------------------
    use ftypes

    implicit none

    integer, intent(in) :: lmax
    real(dp), intent(out) :: p(:)
    real(dp), intent(in) :: z
    integer, intent(out), optional :: exitstatus
    real(dp) :: pm2, pm1, pl
    integer :: l

    if (present(exitstatus)) exitstatus = 0

    if (size(p) < lmax+1) then
        print*, "Error --- PlSchmidt"
        print*, "P must be dimensioned as (LMAX+1) where LMAX is ", lmax
        print*, "Input array is dimensioned ", size(p)
        if (present(exitstatus)) then
            exitstatus = 1
            return
        else
            stop
        end if

    else if (lmax < 0) then
        print*, "Error --- PlSchmidt"
        print*, "LMAX must be greater than or equal to 0."
        print*, "Input value is ", lmax
        if (present(exitstatus)) then
            exitstatus = 2
            return
        else
            stop
        end if

    else if(abs(z) > 1.0_dp) then
        print*, "Error --- PlSchmidt"
        print*, "ABS(Z) must be less than or equal to 1."
        print*, "Input value is ", z
        if (present(exitstatus)) then
            exitstatus = 2
            return
        else
            stop
        end if

    end if

    pm2 = 1.0_dp
    p(1) = 1.0_dp

    pm1 = z
    p(2) = pm1

    do l = 2, lmax, 1
        pl = ( dble(2*l-1)  * z * pm1 - dble(l-1) * pm2 ) / dble(l)
        p(l+1) = pl
        pm2 = pm1
        pm1 = pl

    end do

end subroutine PlSchmidt
