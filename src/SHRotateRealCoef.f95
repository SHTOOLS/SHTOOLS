subroutine SHRotateRealCoef(cilmrot, cilm, lmax, x, dj)
!-------------------------------------------------------------------------------
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
!   Dependencies:   SHrtoc, SHctor, SHcilmtocindex, SHcindextocilm, 
!                   SHRotateCoef, CSPHASE_DEFAULT
!
!   Note: Before using this routine, I would verify that the input euler 
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
    integer :: astat(3)
    real*8, allocatable :: ccilm(:,:,:), cof(:,:), rcof(:,:)
    
    if (size(cilm(:,1,1)) < 2 .or. size(cilm(1,:,1)) < lmax+1 &
            .or. size(cilm(1,1,:)) < lmax+1) then
        print*, "Error --- SHRotateRealCoef"
        print*, "CILM must be dimensioned as (2, LMAX+1, LMAX+1) " // &
                "where LMAX is", lmax
        print*, "Input array is dimensioned ", size(cilm(:,1,1)), &
                size(cilm(1,:,1)), size(cilm(1,1,:)) 
        stop
        
    elseif (size(cilmrot(:,1,1)) < 2 .or. size(cilmrot(1,:,1)) < lmax+1 &
            .or. size(cilmrot(1,1,:)) < lmax+1) then
        print*, "Error --- SHRotateRealCoef"
        print*, "CILMROT must be dimensioned as (2, LMAX+1, LMAX+1) " // &
                "where LMAX is", lmax
        print*, "Input array is dimensioned ", size(cilmrot(:,1,1)), &
                size(cilmrot(1,:,1)), size(cilmrot(1,1,:)) 
        stop
        
    else if (size(dj(:,1,1)) < lmax+1 .or. size(dj(1,:,1)) < lmax+1 &
            .or. size(dj(1,1,:)) < lmax+1) then
        print*, "Error --- SHRotateRealCoef"
        print*, "DJ must be dimensioned as (LMAX+1, LMAX+1, LMAX+1) " // &
                "where LMAX is", lmax
        print*, "Input array is dimensioned ", size(dj(:,1,1)), &
                size(dj(1,:,1)), size(dj(1,1,:)) 
        stop
        
    else if (size(x) < 3) then
        print*, "Error --- SHRotateRealCoef"
        print*, "X must be dimensioned as (3)"
        print*, "Input array is dimensioned ", size(x)
        stop
        
    end if
    
    allocate (ccilm(2,lmax+1,lmax+1), stat = astat(1))
    allocate (cof(2,((lmax+1)*(lmax+2))/2), stat = astat(2))
    allocate (rcof(2,((lmax+1)*(lmax+2))/2), stat = astat(3))
    
    if (sum(astat(1:3)) /= 0) then
        print*, "Error --- SHRotateRealCoef"
        print*, "Problem allocating arrays CCILM, COF, and RCOF", &
            astat(1), astat(2), astat(3)
        stop
        
    end if
    
    if (CSPHASE_DEFAULT == 1) then
        call SHrtoc(cilm, ccilm, degmax=lmax, convention=2, switchcs=1)
        ! Convert geodesy coefficients to Varshalovich et al. complex form
        
    else
        call SHrtoc(cilm, ccilm, degmax=lmax, convention=2, switchcs=0) 
        
    end if
    
    call SHcilmtocindex(ccilm, cof, lmax)
    ! Re-order complex coefficients to form a vector
                            
    call SHRotateCoef(x, cof, rcof, dj, lmax)
    ! Rotate complex re-ordered coefficients
            
    call SHcindextocilm(rcof, ccilm, lmax)
    ! Convert ordered coefficients back to an array
    
    if (CSPHASE_DEFAULT == 1) then
        call SHctor(ccilm, cilmrot, degmax=lmax, convention=2, switchcs=1)
        ! Convert Varshalovich et al complex coefficients back to geodesy form
        
    else
        call SHctor(ccilm, cilmrot, degmax=lmax, convention=2, switchcs=0)
        
    end if
    
    deallocate (ccilm)
    deallocate (cof)
    deallocate (rcof)
    
end subroutine SHRotateRealCoef
