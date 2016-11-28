integer function SHFindLWin(theta0, m, alpha, taper_number)
!------------------------------------------------------------------------------
!
!   This function will output the spherical harmonic bandwidth for
!   the space concentrated taper that yields a concetration factor
!   equal to or greater than the input value of alpha. By default, the
!   concentration factor of the first taper is computed, but this can be
!   modified by inputing the optional arguement taper_number.
!
!   Calling Parameters
!
!       IN
!           theta0          Angular radius of the concentration spherical
!                           cap in RADIANS.
!           alpha           Concentration factor of the window.
!           m               Angular order of the space concetrated tapers.
!
!       IN (OPTIONAL)
!           taper_number    Taper number used to calculate concentration
!                           factors. 
!
!   Dependencies: LAPACK, BLAS, ComputeDm, EigValSym
!
!   Copyright (c) 2016, SHTOOLS
!   All rights reserved.
!
!------------------------------------------------------------------------------
    use SHTOOLS, only: ComputeDm, EigvalSym

    implicit none

    real*8, intent(in) :: theta0, alpha
    integer, intent(in) ::  m
    integer, intent(in), optional :: taper_number
    real*8, allocatable ::  dllm(:,:), eval(:)
    real*8 ::   pi, alpha1
    integer ::  l, astat(2), tn

    if (alpha < 0.0d0 .or. alpha > 1.0d0) then
        print*, "Error --- SHFindLWin"
        print*, "The concentration factor alpha must be between 0 and 1."
        print*, "Input value is ", alpha
        stop
    end if

    pi = acos(-1.0d0)

    if (present(taper_number)) then
        if (taper_number < 1) then
            print*, "Error --- SHFindLWin"
            print*, "TAPER_NUMBER must be greater than 0."
            print*, "Input value is ", taper_number
            stop
        end if

        tn = taper_number

    else
        tn = 1

    end if

    l = tn

    do
        l = l + 1

        allocate (dllm(l+1,l+1), stat = astat(1)) ; dllm = 0.0d0
        allocate (eval(l+1), stat = astat(2)) ; eval = 0.0d0

        if (astat(1) /=0 .or. astat(2) /=0) then 
            print*, "Error --- SHFindLWin"
            print*, "Probelm allocating arrays."
            stop
        end if

        call ComputeDm(dllm, l, abs(m), theta0)
        call EigValSym(dllm, l+1, eval(1:l+1))

        alpha1 = eval(tn)

        deallocate (dllm)
        deallocate (eval)

        if (alpha1 >= alpha) then
            SHFindLWin = l
            return
        end if
    end do

end function SHFindLWin
