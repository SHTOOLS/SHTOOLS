subroutine DHaj(n, aj, exitstatus)
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
!           n   Number of samples in longitude and latitude.
!
!       OUT
!           aj  Vector of length n containing the weights.
!
!   Copyright (c) 2016, SHTOOLS
!   All rights reserved.
!
!------------------------------------------------------------------------------
    implicit none

    integer, intent(in) ::  n
    real*8, intent(out) ::  aj(:)
    integer, intent(out), optional :: exitstatus
    integer :: j, l
    real*8 ::  sum1, pi

    if (present(exitstatus)) exitstatus = 0

    pi = acos(-1.0d0)

    aj = 0.0d0

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

    else if (size(aj) < n) then
        print*, "Error --- DH_aj"
        print*, "The size of AJ must be greater than or equal " // &
                "to N where N is ", n
        print*, "Input array is dimensioned as ", size(aj)
        if (present(exitstatus)) then
            exitstatus = 1
            return
        else
            stop
        end if

    end if

    do j = 0, n-1
        sum1 = 0.0d0

        do l = 0, n/2 -1
            sum1 = sum1 + sin( dble(2*l+1) * pi * dble(j) / dble(n) ) &
                          / dble(2*l+1)
        end do

        aj(j+1) = sum1 * sin(pi * dble(j) / dble(n)) * sqrt(8.0d0) / dble(n)

    end do

end subroutine DHaj
