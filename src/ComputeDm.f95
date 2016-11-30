subroutine ComputeDm(dllm, lmax, m, theta0, exitstatus)
!------------------------------------------------------------------------------
!
!   This routine will compute the kernel D for the space-concentration
!   problem of a spherical cap for a given spherical harmonic order.
!   (When m /= 0, the eigenfunctions will not be isotropic). All terms
!   are computed exactly using Gauss-Legendre quadrature. The eigenfunctions
!   of this kernel are spherical harmonic coefficients normalized according
!   to the "geodesy" convention. The diagonal elements of Dllm approaches
!   unity when theta0 appoaches 180 degrees.
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
!   Copyright (c) 2016, SHTOOLS
!   All rights reserved.
!
!------------------------------------------------------------------------------ 
    use SHTOOLS, only: PreGLQ, PlmBar, NGLQ, PlmIndex

    implicit none

    real*8, intent(out) ::  dllm(:,:)
    real*8, intent(in) ::   theta0
    integer, intent(in) ::  lmax, m
    integer, intent(out), optional :: exitstatus
    real*8 ::   upper, zero(2*lmax+1), w(2*lmax+1), x
    real*8, allocatable ::  plm(:)
    integer ::  l, lp, i, n, astat

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

    dllm = 0.0d0

    n = NGLQ(2*lmax)

    upper = 1.0d0
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
            do lp = l, lmax
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
        dllm = dllm / 2.0d0
    else
        dllm = dllm / 4.0d0
    end if

    call PlmBar(plm, -1, zero(1))   ! deallocate memory

    deallocate (plm)

end subroutine ComputeDm
