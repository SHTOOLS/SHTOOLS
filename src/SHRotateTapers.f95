subroutine SHRotateTapers(tapersrot, tapers, taper_order, lmax, nrot, x, dj, &
                          exitstatus)
!------------------------------------------------------------------------------
!
!   This subroutine will rotate a set of real spherical cap tapers (originally
!   centered at the North pole) corresponding to the angles listed in the
!   input array X. Only the first NROT tapers are rotated, each of which is
!   returned in a column of the output matrix TAPERSROT with spherical harmonic
!   coefficients ordered according to SHCilmToVector.
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
!           tapers          An (lmax+1) by (lmax+1)**2 array containing
!                           all the eigenfunctions of the space-
!                           concentration kernel. Eigenfunctions
!                           are listed by columns in decreasing order
!                           corresponding to value of their eigenvalue.
!           taper_order     A vector of dimension (lmax+1)**2 denoting which
!                           order m corresponds to the column of tapers and
!                           eigenvalues.
!           lmax            Maximum spherical harmonic degree of the tapers.
!           nrot            Number of tapers to rotate.
!           x               Array or rotation angles in radians.
!           dj              Rotation matrix with dimension (lmax+1, lmax+1,
!                                                           lmax+1).
!
!       OUT
!           tapersrot       Rotated real "geodesy" normalized spherical
!                           harmonic coefficients with dimension (lmax+1)**2
!                           by nrot.
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
!                   SHRotateCoef, SHCilmToVector, CSPHASE_DEFAULT
!
!   Note: Before using this routine, please verify that the input Euler
!   angles and signs give the expected results. Some people define the angle
!   beta as a rotation with respect to the x axis.
!
!   Copyright (c) 2018, SHTOOLS
!   All rights reserved.
!
!-------------------------------------------------------------------------------
    use SHTOOLS, only: SHrtoc, SHctor, SHcilmtocindex, SHcindextocilm, &
                       SHRotateCoef, SHCilmToVector, CSPHASE_DEFAULT

    implicit none

    real*8, intent(in) :: tapers(:,:), x(:), dj(:,:,:)
    real*8, intent(out) :: tapersrot(:,:)
    integer, intent(in) :: taper_order(:), lmax, nrot
    integer, intent(out), optional :: exitstatus
    integer :: astat(5), i
    real*8, allocatable :: ccilm(:,:,:), cilm(:,:,:), cof(:,:), rcof(:,:), &
                           vec(:)

    if (present(exitstatus)) exitstatus = 0

    if (size(tapers(:,1)) < (lmax+1) .or. size(tapers(1,:)) < nrot) then
        print*, "Error --- SHRotateTapers"
        print*, "TAPERS must be dimensioned as ( LMAX+1, NROT ) " // &
                "where LMAX = ", lmax, "and NROT = ", nrot
        print*, "Input array is dimensioned as ", size(tapers(:,1)), &
                size(tapers(1,:))
        if (present(exitstatus)) then
            exitstatus = 1
            return
        else
            stop
        end if

    else if (size(taper_order) < nrot) then
        print*, "Error --- SHRotateTapers"
        print*, "TAPER_ORDER must be dimensioned as NROT " // &
                "where NROT = ", nrot
        print*, "Input array is dimensioned as ", size(taper_order)
        if (present(exitstatus)) then
            exitstatus = 1
            return
        else
            stop
        end if

    else if (size(tapersrot(:,1)) < (lmax+1)**2 .or. &
        size(tapersrot(1,:)) < nrot) then
        print*, "Error --- SHRotateTapers"
        print*, "TAPERSROT must be dimensioned as ( (LMAX+1)**2, " // &
                "NROT ), where LMAX = ", lmax, "and NROT = ", nrot
        print*, "Input array is dimensioned as ", size(tapersrot(:,1)), &
                size(tapersrot(1,:))
        if (present(exitstatus)) then
            exitstatus = 1
            return
        else
            stop
        end if

    else if (size(dj(:,1,1)) < lmax+1 .or. size(dj(1,:,1)) < lmax+1 &
            .or. size(dj(1,1,:)) < lmax+1) then
        print*, "Error --- SHRotateTapers"
        print*, "DJ must be dimensioned as (LMAX+1, LMAX+1, LMAX+1) " // &
                "where LMAX = ", lmax
        print*, "Input array is dimensioned ", size(dj(:,1,1)), &
                size(dj(1,:,1)), size(dj(1,1,:)) 
        if (present(exitstatus)) then
            exitstatus = 1
            return
        else
            stop
        end if

    else if (size(x) < 3) then
        print*, "Error --- SHRotateTapers"
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
    allocate (cilm(2,lmax+1,lmax+1), stat = astat(2))
    allocate (cof(2,((lmax+1)*(lmax+2))/2), stat = astat(3))
    allocate (rcof(2,((lmax+1)*(lmax+2))/2), stat = astat(4))
    allocate (vec((lmax+1)**2), stat = astat(5))

    if (sum(astat(1:5)) /= 0) then
        print*, "Error --- SHRotateTapers"
        print*, "Problem allocating arrays CCILM, CILM, COF, RCOF and VEC", &
            astat(1), astat(2), astat(3), astat(4), astat(5)
        if (present(exitstatus)) then
            exitstatus = 3
            return
        else
            stop
        end if

    end if

    do i=1, nrot, 1

        cilm = 0.0d0
        if (taper_order(i) >= 0) then
            cilm(1, 1:lmax+1, taper_order(i)+1) = tapers(1:lmax+1, i)
        else
            cilm(2, 1:lmax+1, abs(taper_order(i))+1) = tapers(1:lmax+1, i)
        endif

        if (CSPHASE_DEFAULT == 1) then
            ! Convert geodesy coefficients to Varshalovich et al. complex form
            if (present(exitstatus)) then
                call SHrtoc(cilm, ccilm, degmax=lmax, convention=2, &
                            switchcs=1, exitstatus=exitstatus)
                if (exitstatus /= 0) return
            else
                call SHrtoc(cilm, ccilm, degmax=lmax, convention=2, switchcs=1)
            endif

        else
            if (present(exitstatus)) then
                call SHrtoc(cilm, ccilm, degmax=lmax, convention=2, &
                            switchcs=0, exitstatus=exitstatus)
                if (exitstatus /= 0) return
            else
                call SHrtoc(cilm, ccilm, degmax=lmax, convention=2, switchcs=0)
            endif

        end if

        if (present(exitstatus)) then
            ! Re-order complex coefficients to form a 2D vector
            call SHcilmtocindex(ccilm, cof, lmax, exitstatus=exitstatus)
            if (exitstatus /= 0) return

            ! Rotate complex re-ordered coefficients
            call SHRotateCoef(x, cof, rcof, dj, lmax, exitstatus=exitstatus)
            if (exitstatus /= 0) return

            ! Convert ordered complex coefficients back to a 3D array
            call SHcindextocilm(rcof, ccilm, lmax, exitstatus=exitstatus)
            if (exitstatus /= 0) return

        else
            ! Re-order complex coefficients to form a 2D vector
            call SHcilmtocindex(ccilm, cof, lmax)

            ! Rotate complex re-ordered coefficients
            call SHRotateCoef(x, cof, rcof, dj, lmax)

            ! Convert ordered complex coefficients back to a 3D array
            call SHcindextocilm(rcof, ccilm, lmax)

        endif

        if (CSPHASE_DEFAULT == 1) then
            ! Convert Varshalovich et al complex coefficients back to 4pi
            ! geodesy form
            if (present(exitstatus)) then
                call SHctor(ccilm, cilm, degmax=lmax, convention=2, &
                            switchcs=1, exitstatus=exitstatus)
                if (exitstatus /= 0) return
            else
                call SHctor(ccilm, cilm, degmax=lmax, convention=2, &
                            switchcs=1)

            endif

        else
            if (present(exitstatus)) then
                call SHctor(ccilm, cilm, degmax=lmax, convention=2, &
                            switchcs=0, exitstatus=exitstatus)
                if (exitstatus /= 0) return
            else
                call SHctor(ccilm, cilm, degmax=lmax, convention=2, &
                            switchcs=0)
            endif

        end if

        if (present(exitstatus)) then
            ! Re-order real coefficients to form a 1-D vector
            call SHCilmToVector(cilm, vec, lmax, exitstatus=exitstatus)
            if (exitstatus /= 0) return

        else
            ! Re-order real coefficients to form a 1-D vector
            call SHCilmToVector(cilm, vec, lmax)

        endif

        tapersrot(1:(lmax+1)**2, i) = vec(1:(lmax+1)**2)

    enddo

    deallocate (ccilm)
    deallocate (cof)
    deallocate (rcof)
    deallocate (cilm)
    deallocate (vec)

end subroutine SHRotateTapers
