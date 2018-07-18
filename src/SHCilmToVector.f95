subroutine SHCilmToVector(cilm, vector, lmax, exitstatus)
!------------------------------------------------------------------------------
!
!   This subroutine will convert an array of spherical harmonic
!   coefficients Cilm to an ordered 1D vector.
!
!   Calling Parameters
!
!       IN
!           cilm    Input spherical harmonic coefficients with dimension
!                      cilm(2, lmax+1, lmax+1).
!           lmax    Maximum spherical harmonic degree of the input coefficients:
!
!       OUT
!           vector  1D vector of ordered spherical harmonic coefficients with
!                   dimension (lmax+1)**2. The ordering is described in 
!                   YilmIndexVector.
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
    implicit none

    real*8, intent(in) :: cilm(:,:,:)
    real*8, intent(out) :: vector(:)
    integer, intent(in) :: lmax
    integer, intent(out), optional :: exitstatus
    integer :: i, l, m

    if (present(exitstatus)) exitstatus = 0

    if (size(cilm(:,1,1)) < 2 .or. size(cilm(1,:,1)) < lmax+1 .or. &
            size(cilm(1,1,:)) < lmax+1) then
        print*, "Error --- SHCilmToVector"
        print*, "CILM must be dimensioned as (2, LMAX+1, LMAX + 1)."
        print*, "LMAX = ", lmax
        print*, "Dimension of CILM = ", size(cilm(:,1,1)),  size(cilm(1,:,1)), &
            size(cilm(1,1,:))
        if (present(exitstatus)) then
            exitstatus = 1
            return
        else
            stop
        end if

    end if

    if (size(vector(:)) < (lmax+1)**2) then
        print*, "Error --- SHCilmToVector"
        print*, "VECTOR must be have dimension (LMAX+1)**2."
        print*, "LMAX = ", lmax
        print*, "Dimension of VECTOR = ", size(vector(:))
        if (present(exitstatus)) then
            exitstatus = 1
            return
        else
            stop
        end if

    end if

    if (lmax < 0) then
        print*, "Error --- SHCilmToVector"
        print*, "LMAX must be positive."
        print*, "LMAX = ", lmax
        if (present(exitstatus)) then
            exitstatus = 2
            return
        else
            stop
        end if

    end if

    i = 0

    do l = 0, lmax
        do m = 0, l
            i = i + 1
            vector(i) = cilm(1, l+1, m+1)

        end do

        do m = 1, l, 1
            i = i + 1
            vector(i) = cilm(2, l+1, m+1)
        end do

    end do

end subroutine SHCilmToVector


subroutine SHVectorToCilm(vector, cilm, lmax, exitstatus)
!------------------------------------------------------------------------------
!
!   This subroutine will convert a 1D ordered vector of spherical harmonic
!   coefficients to a 3D array Cilm(i, l, m).
!
!   Calling Parameters
!
!       IN
!           vector  1D vector of ordered spherical harmonic coefficients with
!                   dimension (lmax+1)**2. The ordering is described in 
!                   YilmIndexVector.
!           lmax    Maximum spherical harmonic degree of the input coefficients.
!
!       OUT
!           cilm    Output spherical harmonic coefficients with dimension
!                   cilm(2, lmax+1, lmax+1).
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
    implicit none

    real*8, intent(out) :: cilm(:,:,:)
    real*8, intent(in) :: vector(:)
    integer, intent(in) :: lmax
    integer, intent(out), optional :: exitstatus
    integer :: k, i, l, m

    if (present(exitstatus)) exitstatus = 0

    if (size(cilm(:,1,1)) < 2 .or. size(cilm(1,:,1)) < lmax+1 &
            .or. size(cilm(1,1,:)) < lmax+1) then
        print*, "Error --- SHVectorToCilm"
        print*, "CILM must be dimensioned as (2, LMAX+1, LMAX + 1)."
        print*, "LMAX = ", lmax
        print*, "Dimension of CILM = ", size(cilm(:,1,1)),  size(cilm(1,:,1)), &
            size(cilm(1,1,:))
        if (present(exitstatus)) then
            exitstatus = 1
            return
        else
            stop
        end if

    end if

    if (size(vector(:)) < (lmax+1)**2) then
        print*, "Error --- SHVectorToCilm"
        print*, "VECTOR must be have dimension (LMAX+1)**2."
        print*, "LMAX = ", lmax
        print*, "Dimension of VECTOR = ", size(vector(:))
        if (present(exitstatus)) then
            exitstatus = 1
            return
        else
            stop
        end if

    end if

    if (lmax < 0) then
        print*, "Error --- SHVectorToCilm"
        print*, "LMAX must be positive."
        print*, "LMAX = ", lmax
        if (present(exitstatus)) then
            exitstatus = 2
            return
        else
            stop
        end if

    end if

    l = 0
    m = 0
    i = 1

    cilm(1,1,1) = vector(1)

    i = 2

    do k=2, (lmax+1)**2, 1
        m = m + 1
    
        if (m > l .and. i == 1) then
            i = 2
            m = 1

        else if (m > l .and. i == 2) then
            l = l + 1
            m = 0
            i = 1

        end if

        cilm(i,l+1,m+1) = vector(k)

    enddo

end subroutine SHVectorToCilm
