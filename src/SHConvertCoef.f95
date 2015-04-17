subroutine SHrtoc(rcilm, ccilm, degmax, convention, switchcs)
!-------------------------------------------------------------------------------
!
!   This subroutine will convert the real "geodesy 4-pi" spherical harmonic 
!   coefficients into complex form with either a 4 pi or unit normalization. 
!   The complex coefficients are only calculated for positive m's (the negative 
!   m's are given by 
!
!       c_l,m* = (-1)^m c_l,-m
!
!   If degmax is not specified, then the maximum degree of the conversion is 
!   taken from the size of the input arrays.
!
!   Calling Parameters
!
!       IN
!           rcilm       Real "geodesy" spherical harmonic coefficients with 
!                       dimensions (2, lmax+1, lmax+1).
!  
!       OUT
!           ccilm       Complex unity-normalized spherical harmonic 
!                       coefficients, dimensioned as (2, lmax+1, lmax+1). The 
!                       first index corresponds to the real and complex 
!                       coefficients, respectively.
!
!       OPTIONAL
!           degmax      Maximum degree of conversion to be performed.
!           convention  1=output 4-pi normalized coefficients
!                       2=output  Varshalovich et al. normalized coefficients
!           switchcs    If 1, Change between different Condon-Shortley phase 
!                       conventions. If 0, use consistent phase convention.
!
!   Copyright (c) 2015, Mark A. Wieczorek
!   All rights reserved.
!
!-------------------------------------------------------------------------------
    implicit none
    
    real*8, intent(in) :: rcilm(:,:,:)
    real*8, intent(out) :: ccilm(:,:,:) 
    integer, intent(in), optional :: degmax, convention, switchcs
    integer :: lmax, l, m, convention_flag, switchcs_flag
    real*8 :: pi
    
    switchcs_flag = 0
    
    if (present(switchcs)) then
        if (switchcs /= 1 .and. switchcs /= 0) then
            print*, "Error --- SHrtoc"
            print*, "switchcs must be equal to either 0 (keep same " // &
                    "convention) of 1 (change Condon-Shortley phase)"
            print*, "Input value is ", switchcs
            stop
            
        end if
        
        switchcs_flag = switchcs
        
    end if
    
    convention_flag = 1
    
    if (present(convention) ) then
        if (convention /= 1 .and. convention /= 2) then
            print*, "Error --- SHrtoc"
            print*, "CONVENTION must be 1 or 2."
            print*, "Input valuse is ", convention
            stop
            
        end if
        
        convention_flag = convention
        
    end if

    if (present(degmax)) then
        lmax = degmax
        
        if (size(rcilm(:,1,1)) < 2 .or. size(rcilm(1,:,1)) < lmax +1 &
                .or. size(rcilm(1,1,:)) < lmax +1) then
            print*, "Error --- SHrtoc"
            print*, "RCILM must be dimensioned as (2, DEGMAX+1, DEGMAX+1) " // &
                    "where DEGMAX is ", degmax
            print*, "Input array is dimensioned as ", size(rcilm(:,1,1)), &
                    size(rcilm(1,:,1)),  size(rcilm(1,1,:))
            stop
            
        else if (size(ccilm(:,1,1)) < 2 .or. size(ccilm(1,:,1)) < lmax +1 &
                .or. size(ccilm(1,1,:)) < lmax +1) then
            print*, "Error --- SHrtoc"
            print*, "CCILM must be dimensioned as (2, DEGMAX+1, DEGMAX+1) " // &
                    "where DEGMAX is ", degmax
            print*, "Input array is dimensioned as ", size(ccilm(:,1,1)), &
                    size(ccilm(1,:,1)),  size(ccilm(1,1,:))
            stop

        end if
        
    else
        if (size(rcilm(:,1,1)) < 2) then
            print*, "Error --- SHrtoc"
            print*, "RCILM must be dimensioned as (2,*,*)."
            print*, "Input array is dimensioned as ",  size(rcilm(:,1,1)), &
                    size(rcilm(1,:,1)), size(rcilm(1,1,:))
            stop
            
        else if (size(ccilm(:,1,1)) < 2) then
            print*, "Error --- SHrtoc"
            print*, "CCILM must be dimensioned as (2,*,*)."
            print*, "Input array is dimensioned as ",  size(ccilm(:,1,1)), &
                    size(ccilm(1,:,1)), size(ccilm(1,1,:))
            stop
            
        end if
        
        lmax = min(size(rcilm(1,1,:)) -1, size(ccilm(1,1,:)) -1, &
                size(rcilm(1,:,1)) -1, size(ccilm(1,:,1)) -1)
                    
    end if

    pi = acos(-1.0d0)
    ccilm = 0.0d0
        
    do l = 0, lmax, 1
        if (convention_flag == 2) then
            ccilm(1,l+1, 1) = sqrt(4.0d0*pi)*rcilm(1,l+1,1)
            ccilm(2,l+1, 1) = 0.0d0
            
            do m = 1, l, 1
                if (switchcs_flag == 1) then
                    ccilm(1, l+1, m+1) = sqrt(2.0d0*pi) * rcilm(1,l+1,m+1) &
                                    * (-1.0d0)**m
                    ccilm(2, l+1, m+1) = -sqrt(2.0d0*pi) * rcilm(2,l+1,m+1) &
                                    * (-1.0d0)**m
                                    
                else
                    ccilm(1, l+1, m+1) = sqrt(2.0d0*pi) * rcilm(1,l+1,m+1) 
                    ccilm(2, l+1, m+1) = -sqrt(2.0d0*pi) * rcilm(2,l+1,m+1)
                    
                end if  
                 
            enddo
            
        else
            ccilm(1,l+1, 1) = rcilm(1,l+1,1)
            ccilm(2,l+1, 1) = 0.0d0
            
            do m = 1, l, 1
                if (switchcs_flag == 1) then
                    ccilm(1, l+1, m+1) = rcilm(1,l+1,m+1) / sqrt(2.0d0) &
                                    * (-1.0d0)**m
                    ccilm(2, l+1, m+1) = -rcilm(2,l+1,m+1) / sqrt(2.0d0) &
                                    * (-1.0d0)**m
                                    
                else
                    ccilm(1, l+1, m+1) = rcilm(1,l+1,m+1) / sqrt(2.0d0)
                    ccilm(2, l+1, m+1) = -rcilm(2,l+1,m+1) / sqrt(2.0d0)
                    
                end if
                
            end do
            
        end if
                
    enddo
    
end subroutine SHrtoc


subroutine SHctor(ccilm, rcilm, degmax, convention, switchcs)
!-------------------------------------------------------------------------------
!
!   This subroutine will convert either "geodesy 4-pi" or "Varshalovich et al." 
!   complex spherical harmonic coefficients into real "geodesy 4-pi" spherical 
!   harmonic coefficients.
!
!   If degmax is not specified, then the maximum degree of the 
!   conversion is taken from the size of the input arrays.
!
!   Calling Parameters
!
!       IN
!           ccilm       Complex unity-normalized spherical harmonic 
!                       coefficients, dimensioned as (2, lmax+1, lmax+1). The 
!                       first index corresponds to the real and complex 
!                       coefficients, respectively.
!       OUT
!           rcilm       Real "geodesy" spherical harmonic coefficients with 
!                       dimensions (2, lmax+1, lmax+1).
!
!       OPTIONAL
!           degmax      Maximum degree of conversion to be performed.
!           convention  1=input coefficients are 4-pi normalized
!                       2=input coefficients are Varshalovich et al. normalized.
!           switchcs    If 1, Change between different Condon-Shortley phase 
!                       conventions. If 0, use consistent phase convention.
!
!   Copyright (c) 2015, Mark A. Wieczorek
!   All rights reserved.
!
!-------------------------------------------------------------------------------
    implicit none
    
    real*8, intent(in) :: ccilm(:,:,:)
    real*8, intent(out) :: rcilm(:,:,:)
    integer, intent(in), optional :: degmax, convention, switchcs
    integer :: lmax, l, m, convention_flag, switchcs_flag
    real*8 :: pi
    
    switchcs_flag = 0
    
    if (present(switchcs)) then
        if (switchcs /= 1 .and. switchcs /=0) then
            print*, "Error --- SHrtoc"
            print*, "switchcs must be equal to either 0 (keep same " // &
                    "convention) of 1 (change Condon-Shortley phase)"
            print*, "Input value is ", switchcs
            stop
            
        end if
        
        switchcs_flag = switchcs
        
    end if
    
    convention_flag = 1
    
    if (present(convention) ) then
        if (convention /=1 .and. convention /=2) then
            print*, "Error --- SHrtoc"
            print*, "CONVENTION must be 1 or 2."
            print*, "Input valuse is ", convention
            stop
            
        end if
        
        convention_flag = convention
        
    end if

    if (present(degmax)) then
        lmax = degmax
        if (size(rcilm(:,1,1)) < 2 .or. size(rcilm(1,:,1)) < lmax +1 .or. &
                size(rcilm(1,1,:)) < lmax +1) then
            print*, "Error --- SHrtoc"
            print*, "RCILM must be dimensioned as (2, DEGMAX+1, DEGMAX+1) " // &
                    "where DEGMAX is ", degmax
            print*, "Input array is dimensioned as ", size(rcilm(:,1,1)), &
                    size(rcilm(1,:,1)),  size(rcilm(1,1,:))
            stop
            
        else if (size(ccilm(:,1,1)) < 2 .or. size(ccilm(1,:,1)) < lmax +1 &
                .or. size(ccilm(1,1,:)) < lmax +1) then
            print*, "Error --- SHrtoc"
            print*, "CCILM must be dimensioned as (2, DEGMAX+1, DEGMAX+1) " // &
                    "where DEGMAX is ", degmax
            print*, "Input array is dimensioned as ", size(ccilm(:,1,1)), &
                    size(ccilm(1,:,1)),  size(ccilm(1,1,:))
            stop

        end if
        
    else
        if (size(rcilm(:,1,1)) < 2) then
            print*, "Error --- SHrtoc"
            print*, "RCILM must be dimensioned as (2,*,*)."
            print*, "Input array is dimensioned as ",  size(rcilm(:,1,1)), &
                    size(rcilm(1,:,1)), size(rcilm(1,1,:))
            stop
            
        else if (size(ccilm(:,1,1)) < 2) then
            print*, "Error --- SHrtoc"
            print*, "CCILM must be dimensioned as (2,*,*)."
            print*, "Input array is dimensioned as ",  size(ccilm(:,1,1)), &
                    size(ccilm(1,:,1)), size(ccilm(1,1,:))
            stop
            
        end if
        
        lmax = min(size(rcilm(1,1,:)) -1, size(ccilm(1,1,:)) -1, &
                size(rcilm(1,:,1)) -1, size(ccilm(1,:,1)) -1)
                    
    end if
    
    pi = acos(-1.0d0)
    rcilm = 0.0d0
        
    do l = 0, lmax, 1
        if (convention == 2) then
            rcilm(1,l+1,1) = ccilm(1,l+1,1) / sqrt(4.0d0*pi)
            rcilm(2, l+1, 1) = 0.0d0
            
            do m = 1, l, 1
                if (switchcs == 1) then
                    rcilm(1,l+1, m+1) =  ccilm(1, l+1, m+1) / sqrt(2.0d0*pi) &
                                    * (-1.0d0)**m
                    rcilm(2,l+1, m+1) =  - ccilm(2,l+1,m+1) / sqrt(2.0d0*pi) &
                                    * (-1.0d0)**m
                                    
                else
                    rcilm(1,l+1, m+1) =  ccilm(1, l+1, m+1) / sqrt(2.0d0*pi)
                    rcilm(2,l+1, m+1) =  - ccilm(2,l+1,m+1) / sqrt(2.0d0*pi)
                    
                end if
                
            end do
        
        else
            rcilm(1,l+1,1) = ccilm(1,l+1,1)
            rcilm(2, l+1, 1) = 0.0d0
            
            do m = 1, l, 1
                if(switchcs == 1) then
                    rcilm(1,l+1, m+1) = sqrt(2.0d0) * ccilm(1, l+1, m+1) &
                                    * (-1.0d0)**m
                    rcilm(2,l+1, m+1) = -sqrt(2.0d0) * ccilm(2, l+1, m+1) &
                                    * (-1.0d0)**m
                else
                    rcilm(1,l+1, m+1) = sqrt(2.0d0) * ccilm(1, l+1, m+1)
                    rcilm(2,l+1, m+1) = -sqrt(2.0d0) * ccilm(2, l+1, m+1)
                    
                end if
                
            end do
        
        end if
        
    end do
    
end subroutine SHctor


subroutine SHCilmToCindex(cilm, cindex, degmax)
!-------------------------------------------------------------------------------
!
!   This routine will convert a 3D matrix of spherical harmonics indexed as 
!   (i, l+1, m+1) into a 2D matrix that is indexed as (i, index) where 
!   index = l(l+1)/2+m+1.
!
!   Calling Parameters
!
!       IN
!           cilm    Array of spherical harmonic coefficients with dimensions 
!                   (2, lmax+1, lmax+1).
!       OUT
!           cindex  Array of indexed spherical harmonic coefficnets with 
!                   dimensions (2, (lmax+1)*(lmax+2)/2).
!
!   Copyright (c) 2015, Mark A. Wieczorek
!   All rights reserved.
!
!-------------------------------------------------------------------------------
    implicit none
    
    real*8, intent(in) :: cilm(:,:,:)
    real*8, intent(out) :: cindex(:,:)
    integer, intent(in), optional :: degmax
    integer :: lmax, l, m, index
    
    if (present(degmax)) then
        lmax = degmax
        
        if (size(cilm(:,1,1)) < 2 .or. size(cilm(1,:,1)) < lmax + 1 .or. &
                size(cilm(1,1,:)) < lmax+1 ) then
            print*, "Error --- SHcilmtocindex"
            print*, "CILM must be dimensioned as (2, DEGMAX+1, DEGMAX+1) " // &
                    "where DEGMAX is ", degmax
            print*, "Input array is dimensioned ", size(cilm(:,1,1)), &
                    size(cilm(1,:,1)), size(cilm(1,1,:))
            stop
            
        else if (size(cindex(:,1)) < 2 .or. size(cindex(1,:)) &
                < ((lmax+1)*(lmax+2))/2) then
            print*, "Error --- SHcilmtocindex"
            print*, "CINDEX must be dimensioned as " // &
                    "(2, (DEGMAX+1)*(DEGMAX+2)/2) where DEGMAX is ", degmax
            print*, "Input array is dimensioned ", size(cindex(:,1)), &
                    size(cindex(1,:))
            stop
            
        end if
        
    else
        lmax = min(size(cilm(1,1,:)) - 1, size(cilm(1,:,1)) - 1)
        
        if (size(cilm(:,1,1)) < 2) then
            print*, "Error --- SHcilmtocindex"
            print*, "CILM must be dimensioned as (2, *, *)."
            print*, "Input array is dimensioned ", size(cilm(:,1,1)), &
                size(cilm(1,:,1)), size(cilm(1,1,:))
            stop
            
        else if (size(cindex(:,1)) < 2 .or. size(cindex(1,:)) &
                < ((lmax+1)*(lmax+2))/2) then
            print*, "Error --- SHcilmtocindex"
            print*, "CINDEX must be dimensioned as " // &
                    "(2, (LMAX+1)*(LMAX+2)/2) where LMAX is ", lmax
            stop
            
        end if
        
    end if
    
    cindex = 0.0d0
        
    do l = 0, lmax
        do m  =0, l
            index = (l*(l+1))/2+m+1
            cindex(1, index) = cilm(1,l+1,m+1)
            cindex(2, index) = cilm(2,l+1,m+1)
        end do
    end do

end subroutine SHCilmToCindex
        
        
subroutine SHCindexToCilm(cindex, cilm, degmax)
!-------------------------------------------------------------------------------
!
!   This subroutine will convert a 2D matrix of spherical harmonics indexed as 
!   (i, index) into a 2D matrix that is index as (i, l+1, m+1) where 
!   index = l(l+1)/2+m+1.
!
!   Calling Parameters
!
!       IN
!           cindex  Array of indexed spherical harmonic coefficnets with 
!                   dimensions (2, (lmax+1)*(lmax+2)/2).
!       OUT
!           cilm    Array of spherical harmonic coefficients with dimensions 
!                   (2, lmax+1, lmax+1).
!
!   Copyright (c) 2015, Mark A. Wieczorek
!   All rights reserved.
!
!-------------------------------------------------------------------------------
    implicit none
    
    real*8, intent(out) :: cilm(:,:,:)
    real*8, intent(in) :: cindex(:,:)
    integer, intent(in), optional :: degmax
    integer :: lmax, l, m, index, n
    
    n = size(cindex(1,:))

    if (present(degmax)) then
        lmax = degmax
        
        if (lmax > nint((-3.0d0 + sqrt(1.0d0 + 8.0d0*n) ) / 2.0d0) ) then
            print*, "Error - SHcindextocilm"
            print*, "The output spherical harmonic degree DEGMAX is larger " // &
                    "than the input coefficients."
            print*, "Input value of DEGMAX ", degmax
            print*, "Maximum spherical harmonic degree of CINDEX ", &
                    nint((-3.0d0 + sqrt(1.0d0 + 8.0d0*n) )/2.0d0)
            stop
            
        elseif (size(cilm(:,1,1)) < 2 .or. size(cilm(1,:,1)) < lmax + 1 &
                .or. size(cilm(1,1,:)) < lmax+1 ) then
            print*, "Error --- SHcindextocilm"
            print*, "CILM must be dimensioned as (2, DEGMAX+1, DEGMAX+1) " // &
                    "where DEGMAX is ", degmax
            print*, "Input array is dimensioned ", size(cilm(:,1,1)), &
                    size(cilm(1,:,1)), size(cilm(1,1,:))
            stop
            
        end if
        
    else    
        lmax = nint((-3.0d0 + sqrt(1.0d0 + 8.0d0*n) )/2.0d0)
        
        if (size(cilm(:,1,1)) < 2 .or. size(cilm(1,:,1)) < lmax + 1 .or. &
                size(cilm(1,1,:)) < lmax+1 ) then
            print*, "Error --- SHcindextocilm"
            print*, "CILM must be dimensioned as (2, DEGMAX+1, DEGMAX+1) " // &
                    "where DEGMAX is ", degmax
            print*, "Input array is dimensioned ", size(cilm(:,1,1)), &
                    size(cilm(1,:,1)), size(cilm(1,1,:))
            stop
            
        end if
        
    end if
    
    cilm = 0.0d0
        
    do l = 0, lmax
        do m = 0, l
            index = (l*(l+1))/2+m+1
            cilm(1,l+1,m+1) = cindex(1, index)
            cilm(2,l+1,m+1) = cindex(2, index)  
        enddo
    enddo   

end subroutine SHCindexToCilm
