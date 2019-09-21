subroutine ComputeDm(dllm, lmax, m, theta0, degrees, exitstatus)
!------------------------------------------------------------------------------
!
!   This routine will compute the kernel D for the space-concentration
!   problem of a spherical cap for a given spherical harmonic order m. All
!   terms are computed exactly using Gauss-Legendre quadrature. The
!   eigenfunctions of this kernel are 4pi-normalized spherical harmonic
!   coefficients that exclude the Condon-Shortley phase. The diagonal elements
!   of Dllm approaches unity when theta0 appoaches 180 degrees. If the optional
!   vector DEGREES is specified, then the matrix will be computed only for
!   degrees l where DEGREES(l+1) is not zero.
!
!   Calling Parameters
!
!       IN
!           lmax        Maximum spherical harmonic degree.
!           theta0      Angular radius of spherical cap IN RADIANS.
!           m           Angular order used in the concentration problem.
!
!       OUT
!           dllm        Symmetric kernel of size lmax+1 by lmax+1.
!
!       OPTIONAL (IN)
!           degrees     Specify those degrees of the coupling matrix to
!                       compute. If degrees(l+1) is zero, degree l will not
!                       be computed.
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
!   Copyright (c) 2005-2019, SHTOOLS
!   All rights reserved.
!
!------------------------------------------------------------------------------ 
    use SHTOOLS, only: PreGLQ, PlmBar, NGLQ, PlmIndex
    use ftypes

    implicit none

    real(dp), intent(out) :: dllm(:,:)
    real(dp), intent(in) :: theta0
    integer, intent(in) :: lmax, m
    integer, intent(in), optional :: degrees(:)
    integer, intent(out), optional :: exitstatus
    real(dp) :: upper, zero(2*lmax+1), w(2*lmax+1), x
    real(dp), allocatable :: plm(:)
    integer :: l, lp, i, n, astat

    if (present(exitstatus)) exitstatus = 0

    if (abs(m) > lmax) then
        print*, "Error --- ComputeDm"
        print*, "M must be less than or equal to LMAX."
        print*, "Input values of LMAX and M are", lmax, m
        if (present(exitstatus)) then
            exitstatus = 2
            return
        else
            stop
        end if

    else if (size(dllm(1,:)) < lmax+1 .or. size(dllm(:,1)) < lmax+1) then
        print*, "Error --- ComputeDm"
        print*, "DLLM must be dimensioned as (LMAX+1, LMAX+1) where " // &
                "LMAX is ", lmax
        print*, "Input array is dimensioned as ", size(dllm(1,:)), & 
                size(dllm(:,1))
        if (present(exitstatus)) then
            exitstatus = 1
            return
        else
            stop
        end if
    end if

    if (present(degrees)) then
        if (size(degrees) < lmax+1) then
            print*, "Error --- ComputeDm"
            print*, "DEGREES must have dimension LMAX+1, where LMAX is ", lmax
            print*, "Input array is dimensioned as ", size(degrees)
            if (present(exitstatus)) then
                exitstatus = 1
                return
            else
                stop
            end if
        end if
    end if

    allocate (plm( (lmax+1)*(lmax+2)/2 ), stat = astat)

    if (astat /= 0) then
        print*, "Error --- ComputeDM"
        print*, "Problem allocating array PLM", astat
        if (present(exitstatus)) then
            exitstatus = 3
            return
        else
            stop
        end if
    end if

    dllm = 0.0_dp

    n = NGLQ(2*lmax)

    upper = 1.0_dp
    x = cos(theta0)

    if (present(exitstatus)) then
        call PreGLQ(x, upper, n, zero, w, exitstatus = exitstatus)
        if (exitstatus /= 0) return
    else
        call PreGLQ(x, upper, n, zero, w)
    end if

    do i = 1, n
        if (present(exitstatus)) then
            call PlmBar(plm, lmax, zero(i), exitstatus = exitstatus)
            if (exitstatus /= 0) return
        else
            call PlmBar(plm, lmax, zero(i))
        end if

        do l = abs(m), lmax
            if (present(degrees)) then
                if (degrees(l+1) == 0) cycle
            end if

            do lp = l, lmax
                if (present(degrees)) then
                    if (degrees(lp+1) == 0) cycle
                end if

                dllm(l+1, lp+1) = dllm(l+1,lp+1) + w(i) * &
                                  plm(PlmIndex(l,abs(m))) * &
                                  plm(PlmIndex(lp,abs(m)))
            end do

        end do

    end do

    do l = abs(m), lmax
        do lp = l + 1, lmax, 1
            dllm(lp+1,l+1) = dllm(l+1,lp+1)
        end do
    end do

    if (m == 0) then
        dllm = dllm / 2.0_dp
    else
        dllm = dllm / 4.0_dp
    end if

    call PlmBar(plm, -1, zero(1))   ! deallocate memory

    deallocate (plm)

end subroutine ComputeDm
