subroutine DHaj(n, aj, extend, exitstatus)
!------------------------------------------------------------------------------
!
!   This subroutine will compute the weights a_j that are used to
!   expand an equally sampled grid into spherical harmonics by using
!   the sampling theorem presented in Driscoll and Healy (1994).
!
!   Note that the number of samples, n, must be even! Also, a_j(1) = 0
!
!   Calling parameters
!
!       IN
!           n       Number of samples in longitude and latitude.
!
!       IN, optional
!           extend  If 1, include the latitudinal band for 90 S, which
!                   increases the dimension of aj by 1.
!
!       OUT, optional
!           aj      Vector of length n or n+1 containing the weights.
!
!   Copyright (c) 2005-2019, SHTOOLS
!   All rights reserved.
!
!------------------------------------------------------------------------------
    use ftypes

    implicit none

    integer(int32), intent(in) :: n
    real(dp), intent(out) :: aj(:)
    integer(int32), intent(in), optional :: extend
    integer(int32), intent(out), optional :: exitstatus
    integer(int32) :: j, l, n_out, extend_grid
    real(dp) :: sum1, pi

    if (present(exitstatus)) exitstatus = 0

    pi = acos(-1.0_dp)

    if (mod(n,2) /= 0) then
        print*, "Error --- DH_aj"
        print*, "The number of samples in the equi-dimensional grid must " // &
                "be even for use with SHExpandDH"
        print*, "Input value of N is ", n
        if (present(exitstatus)) then
            exitstatus = 2
            return
        else
            stop
        end if
    end if

    if (present(extend)) then
        if (extend == 0) then
            extend_grid = 0
            n_out = n
        else if (extend == 1) then
            extend_grid = 1
            n_out = n + 1
        else
            print*, "Error --- DHaj"
            print*, "Optional parameter EXTEND must be 0 or 1."
            print*, "Input value is ", extend
            if (present(exitstatus)) then
                exitstatus = 2
                return
            else
                stop
            end if
        end if
    else
        extend_grid = 0
        n_out = n
    end if

    if (size(aj) < n_out) then
        print*, "Error --- DH_aj"
        print*, "The size of AJ must be greater than or equal " // &
                "to ", n_out
        print*, "Input array is dimensioned as ", size(aj)
        if (present(exitstatus)) then
            exitstatus = 1
            return
        else
            stop
        end if

    end if

    do j = 0, n-1
        sum1 = 0.0_dp

        do l = 0, n/2 -1
            sum1 = sum1 + sin( dble(2*l+1) * pi * dble(j) / dble(n) ) &
                          / dble(2*l+1)
        end do

        aj(j+1) = sum1 * sin(pi * dble(j) / dble(n)) * sqrt(8.0_dp) / dble(n)

    end do

    if (extend_grid == 1) then
        aj(n_out) = 0.0_dp
    end if

end subroutine DHaj
