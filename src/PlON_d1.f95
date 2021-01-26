subroutine PlON_d1(p, dp1, lmax, z, exitstatus)
!------------------------------------------------------------------------------
!
!   This function evalutates all of the "ortho-normalized legendre
!   polynomials, as well as their first derivatives,  up to degree lmax.
!
!   Calling Parameters
!
!       Out
!           p       A vector of all normalized Legendgre polynomials evaluated
!                   at z up to lmax with dimension (lmax+1).
!           dp1     A vector of all first derivatives of the normalized
!                   Legendgre polynomials evaluated at z up to lmax with
!                   dimension (lmax+1).
!
!       IN
!           lmax    Maximum degree to compute.
!           z       [-1, 1], cos(colatitude), or sin(latitude).
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
!   1.  The employed normalization is the "orthonomalization convention."
!   2.  The integral of PlON**2 over all space on the sphere is 1.
!   3.  The integral of PlON**2 over (-1,1) is 1/2pi.
!   4.  The derivative is evaluated with respect to z, and NOT cos(colatitude)
!       or sin(latitude).
!   5.  Derivatives are calculated according to the geodesy-normalized
!       relationships
!           P'_0(z) = 0.0, P'_1(z) = 1.0, and
!           P'_l(z) = l * (P_{l-1}(z) - z * P_l(z) ) / (1.0_dp - z**2)
!           At z = 1, Pl(1) = 1, and P'l(1) = l (l+1) / 2   (Boyd 2001)
!           At z = -1 Pl(-1) = (-1)**l, and P'l(-1) = (-1)**(l-1) l (l+1) / 2
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
    integer(int32), intent(out), optional :: exitstatus
    real(dp) :: pm2, pm1, pl, sinsq, pi
    integer(int32) :: l

    if (present(exitstatus)) exitstatus = 0

    pi = acos(-1.0_dp)

    if (size(p) < lmax+1) then
        print*, "Error --- PlON_d1"
        print*, "P must be dimensioned as (LMAX+1) where LMAX is ", lmax
        print*, "Input array is dimensioned ", size(p)
        if (present(exitstatus)) then
            exitstatus = 1
            return
        else
            stop
        end if

    else if (size(dp1) < lmax+1) then
        print*, "Error --- PlON_d1"
        print*, "DP1 must be dimensioned as (LMAX+1) where LMAX is ", lmax
        print*, "Input array is dimensioned ", size(dp1)
        if (present(exitstatus)) then
            exitstatus = 1
            return
        else
            stop
        end if

     else if (lmax < 0) then
        print*, "Error --- PlON_d1"
        print*, "LMAX must be greater than or equal to 0."
        print*, "Input value is ", lmax
        if (present(exitstatus)) then
            exitstatus = 2
            return
        else
            stop
        end if

    else if (abs(z) > 1.0_dp) then
        print*, "Error --- PlON_d1"
        print*, "ABS(Z) must be less than or equal to 1."
        print*, "Input value is ", z
        if (present(exitstatus)) then
            exitstatus = 2
            return
        else
            stop
        end if

    end if

    if (z == 1.0_dp) then
        do l = 0, lmax
            p(l+1) = sqrt( dble(2*l+1)) / sqrt(4*pi)
            dp1(l+1) = sqrt( dble(2*l+1) ) * dble(l) * dble(l+1) / 2.0_dp &
                        / sqrt(4*pi)
        end do

    else if (z == -1.0_dp) then
        do l = 0, lmax
            p(l+1) = sqrt( dble(2*l+1)) * dble((-1)**l) / sqrt(4*pi)
            dp1(l+1) = sqrt( dble(2*l+1)) * dble(l) * dble(l+1) &
                        * dble((-1)**(l-1)) / 2.0_dp / sqrt(4*pi)
        end do

    else
        sinsq = (1.0_dp - z**2)

        pm2 = 1.0_dp / sqrt(4*pi)
        p(1) = pm2
        dp1(1) = 0.0_dp

        pm1 = sqrt(3.0_dp) * z / sqrt(4*pi)
        p(2) = pm1
        dp1(2) = sqrt(3.0_dp) / sqrt(4*pi)

        do l = 2, lmax, 1
            pl = ( sqrt(dble(2*l-1)) * z * pm1 - &
                    (l-1) * pm2 / sqrt(dble(2*l-3)) ) * &
                    sqrt(dble(2*l+1)) / dble(l)
            p(l+1) = pl
            dp1(l+1) = l * ( sqrt( dble(2*l+1) / dble(2*l-1) ) * &
                    p(l) - z * pl ) / sinsq
            pm2 = pm1
            pm1 = pl
        end do

    end if

end subroutine PlON_d1
