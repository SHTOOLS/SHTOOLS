subroutine SHRotateRealCoef(cilmrot, cilm, lmax, x, dj, exitstatus)
!------------------------------------------------------------------------------
!
!   This subroutine will rotate a set of real spherical harmonic coefficients
!   corresponding to the angles listed in the input array x.
!
!   The rotation of a coordinate system or body can be viewed in two
!   complementary ways involving three successive rotations. Both methods have
!   the same initial and final configurations, and the angles listed in both
!   schemes are the same.
!
!   Scheme A:   (I) Rotation about the z axis by alpha.
!               (II) Rotation about the new y axis by beta.
!               (III) Rotation about the new z axis by gamma.
!
!   Scheme B:   (I) Rotation about the z axis by gamma.
!               (II) Rotation about the initial y axis by beta.
!               (III) Rotation about the initial z axis by alpha.
!
!   The rotations can further be viewed either as either a rotation of the
!   coordinate system or the physical body.
!
!   1. Rotation of the coordinate system without rotation of the physical body,
!       use x(alpha, beta, gamma).
!
!   2. Rotation of the physical body without rotation of the coordinate system,
!       use x(-gamma, -beta, -alpha).
!
!   To perform the inverse trasform of x(alpha, beta, gamma), use
!   x(-gamma, -beta, -alpha).
!
!   This routine uses the "y-convention" were rotations are about the y-axis
!   instead of the x-axis.
!
!   Calling Parameters
!
!       IN
!           cilm        Real "geodesy" normalized spherical harmonic
!                       coefficients with dimension (2, lmax+1, lmax+1).
!           x           Array or rotation angles in radians.
!           lmax        Maximum spherical harmonic degree.
!           dj          Rotation matrix with dimension (lmax+1, lmax+1, lmax+1).
!
!       OUT
!           cilmrot     Rotated real "geodesy" normalized spherical harmonic
!                       coefficients with dimension (2, lmax+1, lmax+1).
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
!
!   Dependencies:   SHrtoc, SHctor, SHcilmtocindex, SHcindextocilm,
!                   SHRotateCoef, CSPHASE_DEFAULT
!
!   Note: Before using this routine, please verify that the input Euler
!   angles and signs give the expected results. Some people define the angle
!   beta as a rotation with respect to the x axis.
!
!   Copyright (c) 2016, SHTOOLS
!   All rights reserved.
!
!-------------------------------------------------------------------------------
    use SHTOOLS, only: SHrtoc, SHctor, SHcilmtocindex, SHcindextocilm, &
                       SHRotateCoef, CSPHASE_DEFAULT

    implicit none

    real*8, intent(in) :: cilm(:,:,:), x(:), dj(:,:,:)
    real*8, intent(out) :: cilmrot(:,:,:)
    integer, intent(in) :: lmax
    integer, intent(out), optional :: exitstatus
    integer :: astat(3)
    real*8, allocatable :: ccilm(:,:,:), cof(:,:), rcof(:,:)

    if (present(exitstatus)) exitstatus = 0

    if (size(cilm(:,1,1)) < 2 .or. size(cilm(1,:,1)) < lmax+1 &
            .or. size(cilm(1,1,:)) < lmax+1) then
        print*, "Error --- SHRotateRealCoef"
        print*, "CILM must be dimensioned as (2, LMAX+1, LMAX+1) " // &
                "where LMAX is", lmax
        print*, "Input array is dimensioned ", size(cilm(:,1,1)), &
                size(cilm(1,:,1)), size(cilm(1,1,:)) 
        if (present(exitstatus)) then
            exitstatus = 1
            return
        else
            stop
        end if

    elseif (size(cilmrot(:,1,1)) < 2 .or. size(cilmrot(1,:,1)) < lmax+1 &
            .or. size(cilmrot(1,1,:)) < lmax+1) then
        print*, "Error --- SHRotateRealCoef"
        print*, "CILMROT must be dimensioned as (2, LMAX+1, LMAX+1) " // &
                "where LMAX is", lmax
        print*, "Input array is dimensioned ", size(cilmrot(:,1,1)), &
                size(cilmrot(1,:,1)), size(cilmrot(1,1,:)) 
        if (present(exitstatus)) then
            exitstatus = 1
            return
        else
            stop
        end if

    else if (size(dj(:,1,1)) < lmax+1 .or. size(dj(1,:,1)) < lmax+1 &
            .or. size(dj(1,1,:)) < lmax+1) then
        print*, "Error --- SHRotateRealCoef"
        print*, "DJ must be dimensioned as (LMAX+1, LMAX+1, LMAX+1) " // &
                "where LMAX is", lmax
        print*, "Input array is dimensioned ", size(dj(:,1,1)), &
                size(dj(1,:,1)), size(dj(1,1,:)) 
        if (present(exitstatus)) then
            exitstatus = 1
            return
        else
            stop
        end if

    else if (size(x) < 3) then
        print*, "Error --- SHRotateRealCoef"
        print*, "X must be dimensioned as (3)"
        print*, "Input array is dimensioned ", size(x)
        if (present(exitstatus)) then
            exitstatus = 1
            return
        else
            stop
        end if

    end if

    allocate (ccilm(2,lmax+1,lmax+1), stat = astat(1))
    allocate (cof(2,((lmax+1)*(lmax+2))/2), stat = astat(2))
    allocate (rcof(2,((lmax+1)*(lmax+2))/2), stat = astat(3))

    if (sum(astat(1:3)) /= 0) then
        print*, "Error --- SHRotateRealCoef"
        print*, "Problem allocating arrays CCILM, COF, and RCOF", &
            astat(1), astat(2), astat(3)
        if (present(exitstatus)) then
            exitstatus = 3
            return
        else
            stop
        end if

    end if

    if (CSPHASE_DEFAULT == 1) then
        ! Convert geodesy coefficients to Varshalovich et al. complex form
        if (present(exitstatus)) then
            call SHrtoc(cilm, ccilm, degmax=lmax, convention=2, switchcs=1, &
                        exitstatus=exitstatus)
            if (exitstatus /= 0) return
        else
            call SHrtoc(cilm, ccilm, degmax=lmax, convention=2, switchcs=1)
        endif

    else
        if (present(exitstatus)) then
            call SHrtoc(cilm, ccilm, degmax=lmax, convention=2, switchcs=0, &
                        exitstatus=exitstatus)
            if (exitstatus /= 0) return
        else
            call SHrtoc(cilm, ccilm, degmax=lmax, convention=2, switchcs=0)
        endif

    end if

    if (present(exitstatus)) then
        ! Re-order complex coefficients to form a vector
        call SHcilmtocindex(ccilm, cof, lmax, exitstatus)
        if (exitstatus /= 0) return

        ! Rotate complex re-ordered coefficients
        call SHRotateCoef(x, cof, rcof, dj, lmax, exitstatus)
        if (exitstatus /= 0) return

        ! Convert ordered coefficients back to an array
        call SHcindextocilm(rcof, ccilm, lmax, exitstatus)
        if (exitstatus /= 0) return

    else
        ! Re-order complex coefficients to form a vector
        call SHcilmtocindex(ccilm, cof, lmax)

        ! Rotate complex re-ordered coefficients
        call SHRotateCoef(x, cof, rcof, dj, lmax)

        ! Convert ordered coefficients back to an array
        call SHcindextocilm(rcof, ccilm, lmax)

    endif

    if (CSPHASE_DEFAULT == 1) then
        ! Convert Varshalovich et al complex coefficients back to geodesy form
        if (present(exitstatus)) then
            call SHctor(ccilm, cilmrot, degmax=lmax, convention=2, switchcs=1, &
                        exitstatus=exitstatus)
            if (exitstatus /= 0) return
        else
            call SHctor(ccilm, cilmrot, degmax=lmax, convention=2, switchcs=1)

        endif

    else
        if (present(exitstatus)) then
            call SHctor(ccilm, cilmrot, degmax=lmax, convention=2, switchcs=0, &
                        exitstatus = exitstatus)
            if (exitstatus /= 0) return
        else
            call SHctor(ccilm, cilmrot, degmax=lmax, convention=2, switchcs=0)
        endif

    end if

    deallocate (ccilm)
    deallocate (cof)
    deallocate (rcof)

end subroutine SHRotateRealCoef
