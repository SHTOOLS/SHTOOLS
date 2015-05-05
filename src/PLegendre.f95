subroutine PLegendre(p, lmax, z)
!-------------------------------------------------------------------------------
!
!   This subroutine evalutates all of the unnormalized Legendre polynomials 
!   up to degree lmax. 
!
!   Calling Parameters
!
!       Out
!           p       A vector of all unnormalized Legendgre polynomials 
!                   evaluated at z up to lmax. The lenght must by greater or 
!                   equal to (lmax+1).
!
!       IN
!           lmax    Maximum degree to compute.
!           z       Value within [-1, 1], cos(colatitude) or sin(latitude).
!
!   Notes:
!   
!   1.  The integral of Pl**2 over (-1,1) is 2/(2l+1).
!   2.  Values are calculated accoring to the following recursion scheme:
!           P_0(z) = 1.0, P_1(z) = z, and 
!           P_l(z) = (2l-1) * z * P_{l-1}(z) / l - (l-1) * P_{l-2}(z) / l
!
!   Dependencies:   None
!
!   Copyright (c) 2015, Mark A. Wieczorek
!   All rights reserved.
!
!-------------------------------------------------------------------------------
    implicit none
    
    integer, intent(in) :: lmax
    real*8, intent(out) :: p(:)
    real*8, intent(in) :: z
    real*8 :: pm2, pm1, pl
    integer :: l

    if (size(p) < lmax+1) then
        print*, "Error --- PlegendreL"
        print*, "P must be dimensioned as (LMAX+1) where LMAX is ", lmax 
        print*, "Input array is dimensioned ", size(p)
        stop
        
    else if (lmax < 0) then 
        print*, "Error --- PlegendreL"
        print*, "LMAX must be greater than or equal to 0."
        print*, "Input value is ", lmax
        stop
        
    else if (abs(z) > 1.0d0) then
        print*, "Error --- PlegendreL"
        print*, "ABS(Z) must be less than or equal to 1."
        print*, "Input value is ", z
        stop
        
    end if
        
    pm2 = 1.d0
    p(1) = 1.d0
        
    pm1  = z
    p(2) = pm1
        
    do l = 2, lmax, 1
        pl = ((2*l-1) * z * pm1 - (l-1) * pm2) / dble(l)
        p(l+1) = pl
        pm2  = pm1
        pm1  = pl
        
    end do

end subroutine PLegendre


subroutine PLegendre_d1(p, dp, lmax, z)
!-------------------------------------------------------------------------------
!
!   This function evalutates all of the unnormalized Legendre polynomials 
!   and their first derivatives up to degree lmax. 
!
!   Calling Parameters
!
!       Out
!           p       A vector of all associated Legendgre polynomials evaluated 
!                   at z up to lmax. The lenght must by greater or equal to 
!                   (lmax+1)
!           dp      A vector of all first derivatives of the associated 
!                   Legendgre polynomials evaluated at z up to lmax. The length 
!                   must by greater or equal to (lmax+1).
!
!       IN
!           lmax    Maximum degree to compute.
!           z       Value within [-1, 1], cos(colatitude) or sin(latitude)
!
!   Notes:
!
!   1.  The integral of Pl**2 over (-1,1) is 2/(2l+1)
!   2.  Values are calculated accoring to the following recursion scheme:
!           P_0(z) = 1.0, P_1(z) = z, and 
!           P_l(z) = (2l-1) * z * P_{l-1}(z) / l - (l-1) * P_{l-2}(z) / l
!   3.  Derivatives are calculated according to
!           P'_0(z) = 0.0, P'_1(z) = 1.0, and
!           P'_l(z) = l * (P_{l-1}(z) - z * P_l(z) ) / (1.0d0 - z**2)
!           At z = 1, Pl(1) = 1, and P'l(1) = l (l+1) / 2   (Boyd 2001)
!           At z = -1 Pl(-1) = (-1)**l, and P'l(-1) = (-1)**(l-1) l (l+1) / 2
!
!   Dependencies:   None
!
!   Copyright (c) 2015, Mark A. Wieczorek
!   All rights reserved.
!
!-------------------------------------------------------------------------------
    implicit none
    
    integer, intent(in) :: lmax
    real*8, intent(out) :: p(:), dp(:)
    real*8, intent(in) :: z
    real*8 ::   pm2, pm1, pl, sinsq
    integer :: l

    if (size(p) < lmax+1) then
        print*, "Error --- PLegendre_d1"
        print*, "P must be dimensioned as (LMAX+1) where LMAX is ", lmax 
        print*, "Input array is dimensioned ", size(p)
        stop
        
    else if (size(dp) < lmax+1) then
        print*, "Error --- PLegendre_d1"
        print*, "DP must be dimensioned as (LMAX+1) where LMAX is ", lmax 
        print*, "Input array is dimensioned ", size(dp)
        stop
        
    else if (lmax < 0) then 
        print*, "Error --- PLegendre_d1"
        print*, "LMAX must be greater than or equal to 0."
        print*, "Input value is ", lmax
        stop
        
    else if (abs(z) > 1.0d0) then
        print*, "Error --- PLegendre_d1"
        print*, "ABS(Z) must be less than or equal to 1."
        print*, "Input value is ", z
        stop
        
    end if
                
    if (z == 1.0d0) then
        p(1:lmax+1) = 1.0d0
        
        do l = 0, lmax
            dp(l+1) = dble(l)*dble(l+1)/2.0d0
        end do
            
    else if (z == -1.0d0) then
        do l =0, lmax
            p(l+1) = dble((-1)**l)
            dp(l+1) = dble(l)*dble(l+1) * dble((-1)**(l-1)) / 2.0d0
        end do
        
    else
        sinsq = (1.0d0 - z**2)
                
        pm2  = 1.d0
        p(1) = 1.d0
        dp(1) = 0.0d0
        
        pm1  = z
        p(2) = pm1
        dp(2) = 1.0d0
        
        do l = 2, lmax, 1
            pl = ((2*l-1) * z * pm1 - (l-1) * pm2) / dble(l)
            p(l+1) = pl
            dp(l+1) = l * (pm1 - z * pl) / sinsq
            pm2  = pm1
            pm1  = pl
        end do
            
    end if

end subroutine PLegendre_d1
