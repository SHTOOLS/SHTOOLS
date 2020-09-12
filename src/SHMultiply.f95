subroutine SHMultiply(cilmout, cilm1, lmax1, cilm2, lmax2, precomp, norm, &
                      csphase, exitstatus)
!------------------------------------------------------------------------------
!
!   This subroutine will multiply two spherical harmonic fields which are
!   expressed up to maximum spherical harmonic degrees lmax1 and lmax2. The
!   output spherical harmonic coefficients will have a maximum spherical
!   harmonic degree equal to lmax1 + lmax2.
!
!   Calling Parameters:
!
!       IN
!           cilm1       Spherical harmonic field with maximum spherical
!                       harmonic degree lmax1.
!           cilm2       Spherical harmonic field with maximun spherical
!                       harmonic degree lmax2.
!           lmax1       Maximum spherical harmonic degree of cilm1.
!           lmax2       Maximum spherical harmonic degree of cilm2.
!
!       OUT 
!           cilmout     Spherical harmonic expansion of spatial
!                       multiplication of cilm1 and cilm2, with a maximum
!                       spherical harmonic degree of lmax1 + lmax2.
!
!       OPTIONAL (IN)
!           precomp     If 1, the array plx will be precomputed when calling
!                       the subroutine SHGLQ. If 0 (default), then this array
!                       will not be precomputed.
!           csphase     1: Do not include the phase factor of (-1)^m
!                       -1: Apply the phase factor of (-1)^m.
!           norm:       Normalization to be used when calculating Legendre
!                       functions
!                           (1) "geodesy" (default)
!                           (2) Schmidt
!                           (3) unnormalized
!                           (4) orthonormalized
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
!   Copyright (c) 2005-2019, SHTOOLS
!   All rights reserved.
!
!------------------------------------------------------------------------------
    use SHTOOLS, only: SHGLQ, MakeGridGLQ, SHExpandGLQ, CSPHASE_DEFAULT
    use ftypes

    implicit none

    real(dp), intent(out) :: cilmout(:,:,:)
    real(dp), intent(in) :: cilm1(:,:,:), cilm2(:,:,:)
    integer, intent(in) :: lmax1, lmax2
    integer, intent(in), optional :: precomp, norm, csphase
    integer, intent(out), optional :: exitstatus
    integer :: lmaxout, phase, mnorm, astat(2), nlat, nlong
    real(dp), allocatable, save :: zero(:), w(:)
    integer, save :: first = 1, lmaxout_last = -1
    real(dp), allocatable :: grid1glq(:,:), grid2glq(:,:), plx(:,:)

!$OMP   threadprivate(zero, w, first, lmaxout_last)

    if (present(exitstatus)) exitstatus = 0

    if (size(cilm1(:,1,1)) < 2 .or. size(cilm1(1,:,1)) < lmax1+1 .or. &
            size(cilm1(1,1,:)) < lmax1+1) then
        print*, "Error --- SHMultiply"
        print*, "CILM1 must be dimensioned as (2, LMAX1+1, LMAX1+1) " // &
                "where LMAX1 is", lmax1
        print*, "Input array is dimensioned ", size(cilm1(:,1,1)), &
                size(cilm1(1,:,1)), size(cilm1(1,1,:))
        if (present(exitstatus)) then
            exitstatus = 1
            return
        else
            stop
        end if

    else if (size(cilm2(:,1,1)) < 2 .or. size(cilm2(1,:,1)) < lmax2+1 &
            .or. size(cilm2(1,1,:)) < lmax2+1) then
        print*, "Error --- SHMultiply"
        print*, "CILM2 must be dimensioned as (2,LMAX2+1, LMAX2+1) " // &
                "where LMAX2 is", lmax2
        print*, "Input array is dimensioned ", size(cilm2(:,1,1)), &
                size(cilm2(1,:,1)), size(cilm2(1,1,:))
        if (present(exitstatus)) then
            exitstatus = 1
            return
        else
            stop
        end if

    else if (size(cilmout(:,1,1)) < 2 &
                .or. size(cilmout(1,:,1)) < lmax1+lmax2+1 &
                .or. size(cilmout(1,1,:)) < lmax1+lmax2+1) then
        print*, "Error --- SHMultiply"
        print*, "CILMOUT must be dimensioned as (2, LMAX1+LMAX2+1, " // &
                "LMAX1+LMAX2+1) where LMAX1 and LMAX2 are", lmax1, lmax2
        print*, "Input array is dimensioned ", size(cilmout(:,1,1)), &
                size(cilmout(1,:,1)), size(cilmout(1,1,:)) 
        if (present(exitstatus)) then
            exitstatus = 1
            return
        else
            stop
        end if

    end if

    if (present(csphase)) then
        if (csphase /= -1 .and. csphase /= 1) then
            print*, "Error --- SHMultiply"
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

    if (present(precomp)) then
        if (precomp /= 1 .and. precomp /= 0) then
            print*, "Error --- SHMultiply"
            print*, "PRECOMP must be either 0 (do not precompute PLX) " // &
                    "or 1 (precompute PLX)."
            print*, "Input value is ", precomp
            if (present(exitstatus)) then
                exitstatus = 2
                return
            else
                stop
            end if

        end if

    end if

    if (present(norm)) then
        if (norm > 4 .or. norm < 1) then
            print*, "Error - SHMultiply"
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

        mnorm = norm

    else
        mnorm = 1

    end if

    lmaxout = lmax1 + lmax2
    nlat = lmax1 + lmax2 + 1
    nlong = 2 * (lmax1 + lmax2) + 1

    if (first == 1) then
        first = 0
        lmaxout_last = lmaxout

        allocate (zero(lmaxout+1), stat = astat(1))
        allocate (w(lmaxout+1), stat = astat(2))

        if (sum(astat(1:2)) /= 0) then
            print*, "Error --- SHMultiply"
            print*, "Problem allocating arrays ZERO and W", astat(1), astat(2)
            if (present(exitstatus)) then
                exitstatus = 3
                return
            else
                stop
            end if

        end if

        if (present(exitstatus)) then
            call SHGLQ(lmaxout, zero, w, csphase = phase, norm = mnorm, &
                       exitstatus = exitstatus)
            if (exitstatus /= 0) return
        else
            call SHGLQ(lmaxout, zero, w, csphase = phase, norm = mnorm)
        end if

    end if

    if (lmaxout /= lmaxout_last) then
        lmaxout_last = lmaxout

        deallocate (zero)
        deallocate (w)
        allocate (zero(lmaxout+1), stat = astat(1))
        allocate (w(lmaxout+1), stat = astat(2))

        if (sum(astat(1:2)) /= 0) then
            print*, "Error --- SHMultiply"
            print*, "Problem allocating arrays ZERO and W", astat(1), astat(2)
            if (present(exitstatus)) then
                exitstatus = 3
                return
            else
                stop
            end if

        end if

        if (present(exitstatus)) then
            call SHGLQ(lmaxout, zero, w, csphase = phase, norm = mnorm, &
                       exitstatus = exitstatus)
            if (exitstatus /= 0) return
        else
            call SHGLQ(lmaxout, zero, w, csphase = phase, norm = mnorm)
        end if

    end if

    allocate (grid1glq(nlat, nlong), stat = astat(1))
    allocate (grid2glq(nlat, nlong), stat = astat(2))

    if (sum(astat(1:2)) /= 0) then
        print*, "Error --- SHMultiply"
        print*, "Problem allocating arrays GRID1GLQ and GRID2GLQ", &
                astat(1), astat(2)
        if (present(exitstatus)) then
            exitstatus = 3
            return
        else
            stop
        end if

    end if

    if (present(precomp)) then
        if (precomp == 0) then
            if (present(exitstatus)) then
                call MakeGridGLQ(grid1glq, cilm1(1:2,1:lmax1+1, 1:lmax1+1), &
                                 lmaxout, zero = zero, csphase = phase, &
                                 norm = mnorm, exitstatus = exitstatus)
                if (exitstatus /= 0) return
                call MakeGridGLQ(grid2glq, cilm2(1:2,1:lmax2+1, 1:lmax2+1), &
                                 lmaxout, zero = zero, csphase = phase, &
                                 norm = mnorm, exitstatus = exitstatus)
                if (exitstatus /= 0) return
                grid1glq(1:nlat,1:nlong) = grid1glq(1:nlat,1:nlong) &
                                           * grid2glq(1:nlat,1:nlong)
                call SHExpandGLQ(cilmout, lmaxout, grid1glq, w, zero = zero, &
                                 csphase = phase, norm = mnorm, &
                                 exitstatus = exitstatus)
                if (exitstatus /= 0) return

            else
                call MakeGridGLQ(grid1glq, cilm1(1:2,1:lmax1+1, 1:lmax1+1), &
                                 lmaxout, zero = zero, csphase = phase, &
                                 norm = mnorm)
                call MakeGridGLQ(grid2glq, cilm2(1:2,1:lmax2+1, 1:lmax2+1), &
                                 lmaxout, zero = zero, csphase = phase, &
                                 norm = mnorm)
                grid1glq(1:nlat,1:nlong) = grid1glq(1:nlat,1:nlong) &
                                           * grid2glq(1:nlat,1:nlong)
                call SHExpandGLQ(cilmout, lmaxout, grid1glq, w, zero = zero, &
                                 csphase = phase, norm = mnorm)
            end if

        else
            allocate (plx(lmax1+lmax2+1, (lmax1+lmax2+1)*(lmax1+lmax2+2)/2 ), &
                      stat=astat(1))
            if (astat(1) /= 0) then
                print*, "Error --- SHMultiply"
                print*, "Problem allocating array PLX", astat(1)
                if (present(exitstatus)) then
                    exitstatus = 3
                    return
                else
                    stop
                end if

            end if

            if (present(exitstatus)) then
                call SHGLQ(lmaxout, zero, w, plx = plx, csphase = phase, &
                           norm = mnorm, exitstatus = exitstatus)
                if (exitstatus /= 0) return
                call MakeGridGLQ(grid1glq, cilm1(1:2,1:lmax1+1, 1:lmax1+1), &
                                 lmaxout, plx = plx, csphase = phase, &
                                 norm = mnorm, exitstatus = exitstatus)
                if (exitstatus /= 0) return
                call MakeGridGLQ(grid2glq, cilm2(1:2,1:lmax2+1, 1:lmax2+1), &
                                 lmaxout, plx = plx, csphase = phase, &
                                 norm = mnorm, exitstatus = exitstatus)
                if (exitstatus /= 0) return
                grid1glq(1:nlat,1:nlong) = grid1glq(1:nlat,1:nlong) &
                                           * grid2glq(1:nlat,1:nlong)
                call SHExpandGLQ(cilmout, lmaxout, grid1glq, w, plx = plx, &
                                 csphase = phase, norm = mnorm, &
                                 exitstatus = exitstatus)
                if (exitstatus /= 0) return

            else
                call SHGLQ(lmaxout, zero, w, plx = plx, csphase = phase, &
                           norm = mnorm)
                call MakeGridGLQ(grid1glq, cilm1(1:2,1:lmax1+1, 1:lmax1+1), &
                                 lmaxout, plx = plx, csphase = phase, &
                                 norm = mnorm)
                call MakeGridGLQ(grid2glq, cilm2(1:2,1:lmax2+1, 1:lmax2+1), &
                                 lmaxout, plx = plx, csphase = phase, &
                                 norm = mnorm)
                grid1glq(1:nlat,1:nlong) = grid1glq(1:nlat,1:nlong) &
                                           * grid2glq(1:nlat,1:nlong)
                call SHExpandGLQ(cilmout, lmaxout, grid1glq, w, plx = plx, &
                                 csphase = phase, norm = mnorm)
            end if

            deallocate (plx)

        end if

    else
        if (present(exitstatus)) then
            call MakeGridGLQ(grid1glq, cilm1(1:2,1:lmax1+1, 1:lmax1+1), &
                             lmaxout, zero = zero, csphase = phase, &
                             norm = mnorm, exitstatus = exitstatus)
            if (exitstatus /= 0) return
            call MakeGridGLQ(grid2glq, cilm2(1:2,1:lmax2+1, 1:lmax2+1), &
                             lmaxout, zero = zero, csphase = phase, &
                             norm = mnorm, exitstatus = exitstatus)
            if (exitstatus /= 0) return
            grid1glq(1:nlat,1:nlong) = grid1glq(1:nlat,1:nlong) &
                                       * grid2glq(1:nlat,1:nlong)
            call SHExpandGLQ(cilmout, lmaxout, grid1glq, w, zero = zero, &
                             csphase = phase, norm = mnorm, &
                             exitstatus = exitstatus)
            if (exitstatus /= 0) return

        else
            call MakeGridGLQ(grid1glq, cilm1(1:2,1:lmax1+1, 1:lmax1+1), &
                             lmaxout, zero = zero, csphase = phase, &
                             norm = mnorm)
            call MakeGridGLQ(grid2glq, cilm2(1:2,1:lmax2+1, 1:lmax2+1), &
                             lmaxout, zero = zero, csphase = phase, &
                             norm = mnorm)
            grid1glq(1:nlat,1:nlong) = grid1glq(1:nlat,1:nlong) &
                                       * grid2glq(1:nlat,1:nlong)
            call SHExpandGLQ(cilmout, lmaxout, grid1glq, w, zero = zero, &
                             csphase = phase, norm = mnorm)
        end if

    end if

    deallocate(grid1glq)
    deallocate(grid2glq)

end subroutine SHMultiply
