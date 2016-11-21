subroutine SHGLQ(lmax, zero, w, plx, norm, csphase, cnorm, exitstatus)
!------------------------------------------------------------------------------
!
!   This routine will precompute the zeros and weights that are
!   necessary for the Gauss-Legendre quadrature sherical harmonic
!   expansion routines. This routine can also calculate a matrix of
!   legendre functions on the grid nodes to speed up computations
!   in the case where more than one expansion is performed.
!   Note, if memory is an issue because of large spherical harmonic
!   degrees, DO NOT PRECOMPUTE THE ARRAY PLX!
!
!   Calling Parameters
!
!       IN
!           lmax        Maximum spherical harmonic degree to
!                       be used in all following expansions.
!
!       OUT
!           zero        An array of length (lmax + 1) which contains the Gauss
!                       points (i.e., zeros of P(l+1,m=0)).
!           w           An array of length lmax + 1 which contains the weights
!                       used in the Gauss-Legendre quadrature routines.
!
!       OPTIONAL
!           plx         Array of normalized associated Legendre Polnomials
!                       computed at the Gauss points. The size of this array
!                       must be (1:lmax+1, 1:(lmax+1)*(lmax+2)/2).
!           norm        Normalization to be used when calculating Legendre
!                       functions
!                           (1) "geodesy" (default)
!                           (2) Schmidt
!                           (3) unnormalized
!                           (4) orthonormalized
!           csphase     1: Do not include the phase factor of (-1)^m
!                       -1: Apply the phase factor of (-1)^m.
!           cnorm:      1: compute complex normalized Legendre functions.
!                       0 (default): compute real normalized Legendre functions.
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
!   Dependencies:   PreGLQ, PlmBar, PlmSchmidt, PLegendreA, PlmON,
!                   CSPHASE_DEFAULT, PlmIndex
!
!   Copyright (c) 2016, SHTOOLS
!   All rights reserved.
!
!------------------------------------------------------------------------------
    use SHTOOLS, only: PreGLQ, PlmBar, PlmSchmidt, PLegendreA, PlmON, &
                       CSPHASE_DEFAULT, PlmIndex

    implicit none

    integer, intent(in) :: lmax
    real*8, intent(out) :: zero(:), w(:)
    real*8, intent(out), optional :: plx(:,:)
    integer, intent(in), optional :: norm, csphase, cnorm
    integer, intent(out), optional :: exitstatus
    integer :: n, i, astat, phase, l, m, i_s, cnormin, lnorm
    real*8 :: upper, lower, pi
    real*8, allocatable ::  pl(:)

    if (present(exitstatus)) exitstatus = 0

    if (size(zero) < lmax+1) then
        print*, "Error --- SHGLQ"
        print*, "ZERO must be dimensioned as (LMAX+1) where LMAX is ", lmax
        print*, "Input array is dimensioned ", size(zero)
        if (present(exitstatus)) then
            exitstatus = 1
            return
        else
            stop
        end if

    else if (size(w) < lmax+1) then
        print*, "Error --- SHGLQ"
        print*, "W must be dimensioned as (LMAX+1) where LMAX is ", lmax
        print*, "Input array is dimensioned ", size(w)
        if (present(exitstatus)) then
            exitstatus = 1
            return
        else
            stop
        end if

    end if

    if (present(plx)) then
        if (size(plx(:,1)) < lmax +1 .or. size(plx(1,:)) &
                < (lmax+1)*(lmax+2)/2) then
            print*, "Error --- SHGLQ"
            print*, "PLX must be dimensioned as (LMAX+1, " // &
                    "(LMAX+1)*(LMAX+2)/2) where LMAX is ", lmax
            print*, "Input array is dimensioned as ", size(plx(:,1)), &
                    size(plx(1,:))
            if (present(exitstatus)) then
                exitstatus = 1
                return
            else
                stop
            end if

        end if

    end if
    
    if (present(norm)) then
        if (norm > 4 .or. norm < 1) then
            print*, "Error --- SHGLQ"
            print*, "Parameter NORM must be 1 (geodesy), 2 (Schmidt), " // &
                    "3 (unnormalized), or 4 (orthonormalized)."
            print*, "Input value is ", norm
            if (present(exitstatus)) then
                exitstatus = 2
                return
            else
                stop
            end if

        end if

        lnorm = norm

    else
        lnorm = 1
        
    end if

    if (present(csphase)) then
        if (csphase /= -1 .and. csphase /= 1) then
            print*, "Error --- SHGLQ"
            print*, "CSPHASE must be 1 (exclude) or -1 (include)."
            print*, "Input value is ", csphase
            if (present(exitstatus)) then
                exitstatus = 2
                return
            else
                stop
            end if

        else
            phase = csphase

        end if

    else
        phase = CSPHASE_DEFAULT

    end if

    allocate (pl(((lmax+2)*(lmax+1))/2), stat = astat)

    if (astat /= 0) then
        print*, "Error --- SHGLQ"
        print*, "Problem allocating array PL", astat
        if (present(exitstatus)) then
            exitstatus = 3
            return
        else
            stop
        end if
    end if

    if (present(cnorm)) then
        cnormin = cnorm

    else
        cnormin = 0

    end if

    pi = acos(-1.0d0)

    upper = 1.0d0
    lower = -1.0d0

    n = lmax + 1
    ! lmax is the maxium degree of the equation we
    ! are going to integrate.

    ! Determine Gauss Points and Weights.
    if (present(exitstatus)) then
        call PreGLQ(lower, upper, n, zero, w, exitstatus = exitstatus)
        if (exitstatus /= 0) return
    else
        call PreGLQ(lower, upper, n, zero, w)
    endif
    
    !--------------------------------------------------------------------------
    !
    !   Next, fill the array plx with all of the associated legendre polynomials
    !   evaluated at the Gauss points. The first index of this array corresponds
    !   to the Gauss point, whereas the (lmax+2)*(lmax+1)/2 elements in the
    !   second column correspond to all of the assiated legendre polynomials,
    !   ordered according to index = (l+1)*l/2 + m + 1.
    !
    !--------------------------------------------------------------------------
    if (present(plx)) then
        do i = 1, (n + 1) / 2
            if (present(exitstatus)) then
                select case(lnorm)
                    case(1)
                        call PlmBar(pl, lmax, zero(i), csphase = phase, &
                                    cnorm=cnormin, exitstatus = exitstatus)
                    case(2)
                        call PlmSchmidt(pl, lmax, zero(i), csphase = phase, &
                                        cnorm=cnormin, exitstatus = exitstatus)
                    case(3)
                        call PLegendreA(pl, lmax, zero(i), csphase = phase, &
                                        exitstatus = exitstatus)
                    case(4)
                        call PlmON(pl, lmax, zero(i), csphase = phase, &
                                   cnorm=cnormin, exitstatus = exitstatus)
                    if (exitstatus /= 0) return
                    
                end select

            else
                select case(lnorm)
                    case(1)
                        call PlmBar(pl, lmax, zero(i), csphase = phase, &
                                    cnorm=cnormin)
                    case(2)
                        call PlmSchmidt(pl, lmax, zero(i), csphase = phase, &
                                        cnorm=cnormin)
                    case(3)
                        call PLegendreA(pl, lmax, zero(i), csphase = phase)
                    case(4)
                        call PlmON(pl, lmax, zero(i), csphase = phase, &
                                   cnorm=cnormin)

                end select

            end if

            if (i == (n + 1) / 2 .and. mod(n,2) /= 0) then
                plx(i, 1:((lmax+1)*(lmax+2))/2) = pl(1:((lmax+1)*(lmax+2))/2)

            else
                i_s = n + 1 - i
                plx(i, 1:((lmax+1)*(lmax+2))/2) = pl(1:((lmax+1)*(lmax+2))/2)

                do l = 0, lmax
                    do m = 0, l
                        plx(i_s,PlmIndex(l,m)) = plx(i,PlmIndex(l,m)) &
                                                 * (-1.0d0)**(l-m)

                    end do

                end do

            end if

        end do

        select case(lnorm)
            case(1)
                call PlmBar(pl, -1, zero(1), csphase = phase, cnorm=cnormin)
            case(2)
                call PlmSchmidt(pl, -1, zero(1), csphase = phase, cnorm=cnormin)
            case(4)
                call PlmON(pl, -1, zero(1), csphase = phase, cnorm=cnormin)

        end select

    end if

    deallocate (pl)

end subroutine SHGLQ
