subroutine ComputeDG82(DG82, lmax, m, theta0, exitstatus)
!-------------------------------------------------------------------------------
!
!   This routine will compute the kernel of Grunbaum et al. 1982
!   that commutes with the space concentration kernel. This
!   kernel is tridiagonal and has simple expressions for
!   its elements. Note that this kernel is multiplied by -1 in
!   comparison to Grunbaum et al. in order that the eigenvectors
!   will be listed in the same order as the eigenvectors of the
!   equivalent space concentration matrix Dllm. While the eigenfunctions
!   of this kernel correspond to the eigenvalues of Dllm, the eigenvalues
!   do NOT!
!
!   Calling Parameters
!       IN
!           lmax        Maximum spherical harmonic degree.
!           theta0      Angular radius of spherical cap IN RADIANS.
!           m           Angular order to concentration problem.
!
!       OUT
!           D0G82       Symmetric tridiagonal kernel with a maximum size of
!                       lmax+1 by lmax+1.
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
!   Dependencies:   none
!
!   Copyright (c) 2016, SHTOOLS
!   All rights reserved.
!
!------------------------------------------------------------------------------
    implicit none

    real*8, intent(out) :: DG82(:,:)
    real*8, intent(in) :: theta0
    integer, intent(in) :: lmax, m
    integer, intent(out), optional :: exitstatus
    real*8 :: x
    integer :: i, n

    if (present(exitstatus)) exitstatus = 0

    n = lmax + 1 - abs(m)

    if (n < 1) then
        print*, "Error --- ComputeDG82"
        print*, "abs(M) must be less than or equal to LMAX."
        print*, "Input values of l and m are ", lmax, m
        if (present(exitstatus)) then
            exitstatus = 2
            return
        else
            stop
        end if

    else if (size(DG82(1,:)) < n .or. size(DG82(:,1)) < n) then
        print*, "Error --- ComputeDG82"
        print*, "DG82 must be dimensioned as (LMAX-abs(M)+1," // &
                "LMAX-abs(M)+1) where LMAX and M are ", lmax, m
        print*, "Input array is dimensioned as ", size(DG82(1,:)), &
                size(DG82(:,1))
        if (present(exitstatus)) then
            exitstatus = 1
            return
        else
            stop
        end if

    end if

    DG82 = 0.0d0

    x = cos(theta0)

    DG82(1,1) = x * dble(1+m) * dble(m)

    do i = 2, n, 1
        DG82(i,i) = x * dble(i+m) * dble(i+m-1)

        DG82(i,i-1) = - sqrt(dble(i-1+m)**2 - dble(m)**2) * &
                      ( dble(i-1+m)**2 - dble(lmax+1)**2 ) / &
                      sqrt(4.0d0*dble(i-1+m)**2 - 1.0d0)

        DG82(i-1,i) = DG82(i,i-1)

    end do

end subroutine ComputeDG82
