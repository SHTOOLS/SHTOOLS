subroutine SHRotateCoef(x, cof, rcof, dj, lmax)
!-------------------------------------------------------------------------------
!
!   This subroutine will rotate a set of spherical harmonic coefficients (see 
!   documentation in SHRotateRealCoef for a description of how rotation angles 
!   are defined.) All angles are measured in radians, and spherical harmonic
!   coefficients are fully normalized as in Edmonds, 1960. Note that euler 
!   angles are different in Goldstein and Edmonds.
!
!   This routine uses the "y-convention" where rotations are about the y axis
!   instead of the x axis.
!
!   djpi2 must be called before using this routine.
!
!   Calling Parameters
!
!       IN
!           x       Array of dimension 3 containing the three Euler angles.
!           dj      Roation matrix with dimension (lmax+1, lmax+1, lmax+1).
!           cof     Indexed spherical harmonic coefficients with dimensions
!                   (2, (lmax+1)*(lmax+2)/2).
!           lmax    Maximum spherical harmonic degree.
!
!       OUT
!           rcof    Indexed rotated spherical harmonic coefficients with 
!                   dimensions (2, (lmax+1)*(lmax+2)/2).
!
!   History:
!       1.  Based on routine from Guy Masters (July16, 1993)
!       2.  Modified by Mark Simons (July 25, 1993)
!       3.  Turned into readable f95 code by Mark Wieczorek (August, 2003)
!
!   Copyright (c) 2015, Mark A. Wieczorek
!   All rights reserved.
!
!-------------------------------------------------------------------------------
    implicit none
    
    real*8, intent(in) :: cof(:,:), dj(:,:,:), x(3)
    real*8, intent(out) :: rcof(:,:)
    integer, intent(in) :: lmax
    real*8 :: sum(2), temp(2,lmax+1), temp2(2,lmax+1), cgam(lmax+1), &
                sgam(lmax+1), calf(lmax+1), salf(lmax+1), cbet(lmax+1), &
                sbet(lmax+1), pi2, alpha, beta, gamma
    integer ::  ind, lp1, l, mp1, jp1, isgn, ii, indx
    
    if (size(cof(:,1)) < 2 .or. size(cof(1,:)) < ((lmax+1)*(lmax+2))/2) then
        print*, "Error --- SHRotateCoef"
        print*, "COEF must be dimensioned (2, (LMAX+1)*(LMAX+2)/2) " // &
                "where LMAX is ", lmax
        print*, "Input array is dimensioned ", size(cof(:,1)),  size(cof(1,:))
        stop
        
    else if (size(rcof(:,1)) < 2 .or. size(rcof(1,:)) &
            < ((lmax+1)*(lmax+2))/2) then
        print*, "Error --- SHRotateCoef"
        print*, "RCOEF must be dimensioned (2, (LMAX+1)*(LMAX+2)/2) " // &
                "where LMAX is ", lmax
        print*, "Input array is dimensioned ", size(rcof(:,1)),  size(rcof(1,:))
        stop
        
    else if (size(dj(:,1,1)) < lmax+1 .or. size(dj(1,:,1)) < lmax+1 &
            .or. size(dj(1,1,:)) < lmax+1 ) then
        print*, "Error --- SHRotateCoef"
        print*, "DJ must be dimensioned (LMAX+1, LMAX+1, LMAX+1) " // &
                "where LMAX is ", lmax
        print*, "Input array is dimensioned ", size(dj(:,1,1)), &
                size(dj(1,:,1)), size(dj(1,1,:))
        stop
        
    end if
    
    rcof = 0.0d0

    pi2 = 1.570796326794895d0
    
    alpha = x(1)
    beta  = x(2)
    gamma = x(3)
    
    alpha = alpha-pi2
    gamma = gamma+pi2
    beta = -beta
    
    ind = 0
    
    ! Loop over degrees
    do lp1 = 1, lmax+1 
        l = lp1-1
        cbet(lp1) = cos(l*beta)
        sbet(lp1) = sin(l*beta)
        cgam(lp1) = cos(l*gamma)
        sgam(lp1) = sin(l*gamma)
        calf(lp1) = cos(l*alpha)
        salf(lp1) = sin(l*alpha)

        ! Alpha rotation
        do mp1 = 1, lp1
            indx = ind+mp1
            temp(1,mp1) = cof(1,indx) * calf(mp1) - cof(2,indx) * salf(mp1)
            temp(2,mp1) = cof(2,indx) * calf(mp1) + cof(1,indx) * salf(mp1)
            
        end do
        
        ! B rotation and beta rotation
        do jp1 = 1, lp1
            sum(1) = dj(jp1,1,lp1)*temp(1,1)
            sum(2) = 0.d0
            isgn = 1-2*mod((lp1-jp1),2)
            
            do mp1 = 2, lp1
                isgn = -isgn
                ii = (3-isgn)/2
                sum(ii) = sum(ii)+2.d0*dj(jp1,mp1,lp1)*temp(ii,mp1)
            end do
            
            temp2(1,jp1) = sum(1)*cbet(jp1)-sum(2)*sbet(jp1)
            temp2(2,jp1) = sum(2)*cbet(jp1)+sum(1)*sbet(jp1)
        
        end do

        ! Inverse B rotation and gamma rotation
        do jp1 = 1, lp1
            sum(1) = dj(1,jp1,lp1)*temp2(1,1)
            sum(2) = 0.d0
            isgn = 1-2*mod((lp1-jp1),2)
            
            do mp1 = 2,lp1
                isgn = -isgn
                ii = (3-isgn)/2
                sum(ii) = sum(ii)+2.d0*dj(mp1,jp1,lp1)*temp2(ii,mp1)
                
            end do
            
            indx = ind+jp1
            rcof(1,indx) = sum(1)*cgam(jp1)-sum(2)*sgam(jp1)
            rcof(2,indx) = sum(2)*cgam(jp1)+sum(1)*sgam(jp1)
            
        end do
        
        ind = ind+lp1
        
    end do
    
end subroutine SHRotateCoef
