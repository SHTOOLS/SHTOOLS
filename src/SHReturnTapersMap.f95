subroutine SHReturnTapersMap(tapers, eigenvalues, dh_mask, n_dh, lmax, &
                             sampling, ntapers, degrees, exitstatus)
!------------------------------------------------------------------------------
!
!   This subroutine will calculate the eigenvalues and eigenfunctions of the
!   generalized concentration problem where the region of interest R is given
!   by the grid DH_MASK. This matrix is sampled according to the Driscoll and
!   Healy sampling theorem, and possesses a value of 1 inside of R, and 0
!   elsewhere. Returned tapers and eigenvalues are ordered from the largest to
!   smallest eigenvalue, and the spectral coefficients are packed into a 1D
!   column vector according to the scheme described in YilmIndexVector.
!
!   The elements Dij are calculated approximately using spherical harmonic
!   transforms. The effective bandwidth of the grid DH_MASK should in general
!   be larger than LMAX by about a factor of 3 or so for accurate results. In
!   addition, as a result of the approximate manner in which the space-
!   concentration kernel is calculated, it is preferable to use SAMPLING=2.
!
!   Calling Parameters
!
!       IN
!           dh_mask     Integer grid sampled according to the Driscoll and
!                       Healy sampling theorem. A value of 1 indicates that the
!                       grid node is in the concentration domain, and a value
!                       of 0 indicates that it is outside. Dimensioned as
!                       (n_dh, n_dh) for SAMPLING = 1 or (n_dh, 2*n_dh) for
!                       SAMPLING = 2.
!           n_dh        The number of latitude samples in the Driscoll and
!                       Healy sampled grid.
!           lmax        Maximum spherical harmonic degree of the outpt
!                       spherical harmonic coefficients.
!           sampling    1 corresponds to equal sampling (n_dh, n_dh), whereas
!                       2 corresponds to equal spaced grids (n_dh, 2*n_dh).
!
!       IN, OPTIONAL
!           ntapers     Number of tapers and eigenvalues to output.
!           degrees     Specify those degrees of the coupling matrix to
!                       compute. If degrees(l+1) is zero, degree l will not
!                       be used.
!
!       OUT
!           tapers      Column vectors contain the spherical harmonic
!                       coefficients, packed according to the scheme described
!                       in YilmIndexVector. The dimension of this array is
!                       (lmax+1)**2 by (lmax+1)**2, or (lmax+1)**2 by ntapers
!                       if ntapers is present.
!           eigenvalues A 1-dimensional vector containing the eigenvalues
!                       corresponding to the columns of tapers, dimensioned as
!                       (lmax+1)**2.
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
!   Copyright (c) 2019, SHTOOLS
!   All rights reserved.
!
!------------------------------------------------------------------------------
    use SHTOOLS, only : EigValVecSym, ComputeDMap

    implicit none

    real*8, intent(out) :: tapers(:,:), eigenvalues(:)
    integer, intent(in) :: dh_mask(:,:), n_dh, lmax, sampling
    integer, intent(in), optional :: ntapers, degrees(:)
    integer, intent(out), optional :: exitstatus
    real*8, allocatable :: dij(:,:), dijex(:, :), evec(:, :)
    integer :: nlat, nlong, lmax_dh, astat(2), i, j, l, m, exclude, n, &
               ind((lmax+1)**2), numk

    if (present(exitstatus)) exitstatus = 0

    if (present(ntapers)) then
        if (ntapers > (lmax+1)**2) then
            print*, "Error --- SHRetrunTapersMap"
            print*, "The number of output tapers must be less than or " // &
                    "equal to (LMAX+1)**2."
            print*, "LMAX = ", lmax
            print*, "NTAPERS = ", ntapers
            if (present(exitstatus)) then
                exitstatus = 2
                return
            else
                stop
            end if

        else if (size(tapers(:,1)) < (lmax+1)**2 .or. &
            size(tapers(1,:)) < ntapers) then
                print*, "Error --- SHReturnTapersMap"
                print*, "TAPERS must be dimensioned as ((LMAX+1)**2, NTAPERS)."
                print*, "Dimension of TAPERS = ", size(tapers(:,1)), &
                        size(tapers(1,:))
                print*, "LMAX = ", lmax
                print*, "NTAPERS = ", ntapers
                if (present(exitstatus)) then
                    exitstatus = 1
                    return
                else
                    stop
                end if

        else if (size(eigenvalues(:)) < ntapers) then
            print*, "Error --- SHReturnTapersMap"
            print*, "EIGENVALUES must be dimensioned as NTAPERS."
            print*, "Dimension of EIGENVALUES = ", size(eigenvalues)
            print*, "NTAPERS = ", ntapers
            if (present(exitstatus)) then
                exitstatus = 1
                return
            else
                stop
            end if

        end if

    else
        if (size(tapers(:,1)) < (lmax+1)**2 .or. &
                size(tapers(1,:)) < (lmax+1)**2) then
            print*, "Error --- SHReturnTapersMap"
            print*, "TAPERS must be dimensioned as ((LMAX+1)**2, (LMAX+1)**2)."
            print*, "Dimension of TAPERS = ", size(tapers(:,1)), &
                    size(tapers(1,:))
            print*, "LMAX = ", lmax
            if (present(exitstatus)) then
                exitstatus = 1
                return
            else
                stop
            end if

        else if (size(eigenvalues(:)) < (lmax+1)**2) then
            print*, "Error --- SHReturnTapersMap"
            print*, "EIGENVALUES must be dimensioned as (LMAX+1)**2"
            print*, "Dimension of EIGENVALUES = ", size(eigenvalues)
            print*, "LMAX = ", lmax
            if (present(exitstatus)) then
                exitstatus = 1
                return
            else
                stop
            end if

        end if

    endif

    if (mod(n_dh,2) /= 0) then
        print*, "Error --- SHReturnTapersMap"
        print*, "Number of samples in latitude must be even for the " // &
                "Driscoll and Healy sampling theorem."
        print*, "N_DH = ", n_dh
        if (present(exitstatus)) then
             exitstatus = 2
             return
        else
             stop
        end if

    end if

    if (present(degrees)) then
        if (size(degrees) < lmax+1) then
            print*, "Error --- SHReturnTapersMap"
            print*, "DEGREES must have dimension LMAX+1, where LMAX is ", lmax
            print*, "Input array is dimensioned as ", size(degrees)
            if (present(exitstatus)) then
                exitstatus = 1
                return
            else
                stop
            end if
        end if
    end if

    nlat = n_dh

    if (sampling == 1) then
        nlong = nlat

    else if (sampling == 2) then
        nlong = 2 * nlat

    else
        print*, "Error --- SHReturnTapersMap"
        print*, "SAMPLING must be either 1 (equally sampled) " // &
                "or 2 (equally spaced)."
        print*, "SAMPLING = ", sampling
        if (present(exitstatus)) then
            exitstatus = 2
            return
        else
            stop
        end if

    end if

    if (size(dh_mask(:,1)) < nlat .or. size(dh_mask(1,:)) < nlong) then
        print*, "Error --- SHReturnTapersMap"
        print*, "DH_MASK must be dimensioned as ", nlat, nlong
        print*, "Dimensions of DH_MASK are ", size(dh_mask(:,1)), &
                size(dh_mask(1,:))
        if (present(exitstatus)) then
            exitstatus = 1
            return
        else
            stop
        end if

    end if

    lmax_dh = n_dh/2 - 1

    if (lmax_dh < lmax) then
        print*, "Error --- SHReturnTapersMap"
        print*, "The effective bandwith of the input grid DH_MASK must " // &
                "be greater or equal than LMAX."
        print*, "In practice, this should be about 3 * LMAX."
        print*, "LMAX = ", lmax
        print*, "Effective bandwidth of DH_MASK = (N/2 -1) = ", lmax_dh
        if (present(exitstatus)) then
            exitstatus = 2
            return
        else
            stop
        end if

    end if

    eigenvalues = 0.0d0
    tapers = 0.0d0

    allocate (dij((lmax+1)**2, (lmax+1)**2), stat = astat(1))
    allocate (evec((lmax+1)**2, (lmax+1)**2), stat = astat(2))

    if (astat(1) /= 0 .or. astat(2) /= 0) then
        print*, "Error --- SHReturnTapersMap"
        print*, "Problem allocating arrays DIJ and EVEC."
        if (present(exitstatus)) then
            exitstatus = 3
            return
        else
            stop
        end if

    end if

    if (present(exitstatus)) then
        call ComputeDMap(Dij, dh_mask, n_dh, lmax, sampling = sampling, &
                         exitstatus = exitstatus)
        if (exitstatus /= 0) return

    else
        call ComputeDMap(Dij, dh_mask, n_dh, lmax, sampling = sampling)

    end if

    exclude = 0
    n = (lmax+1)**2

    if (present(degrees)) then
        if (sum(degrees(1:lmax+1)) /= lmax + 1 .and. &
                sum(degrees(1:lmax+1)) /= 0) then
            exclude = 1
            ind = 0
            i = 0

            do l = 0, lmax
                if (degrees(l+1) /= 0) then
                    do m = 0, l
                        i = i + 1
                        ind(i) = l**2 + m + 1
                    end do

                    do m = 1, l, 1
                        i = i + 1
                        ind(i) = l**2 + l + m + 1
                    end do

                end if
            end do

            n = i

            allocate (dijex(n, n), stat = astat(1))
            if (astat(1) /= 0) then
                print*, "Error --- SHReturnTapersMap"
                print*, "Problem allocating array DIJEX.", astat(1)
                if (present(exitstatus)) then
                    exitstatus = 3
                    return
                else
                    stop
                end if
            end if

            dijex = 0.0d0
            do i = 1, n
                do j = 1, n
                    dijex(i, j) = dij(ind(i), ind(j))
                end do
            end do

            dij(1:n, 1:n) = dijex(1:n, 1:n)

        end if

    end if

    if (present(exitstatus)) then
        if (present(ntapers)) then
            call EigValVecSym(Dij(1:n, 1:n), n, eigenvalues, evec, &
                              k = min(ntapers, n), exitstatus = exitstatus)
            if (exitstatus /= 0) return
        else
            call EigValVecSym(Dij(1:n, 1:n), n, eigenvalues, evec, &
                              exitstatus = exitstatus)
            if (exitstatus /= 0) return
        endif

    else
        if (present(ntapers)) then
            call EigValVecSym(Dij(1:n, 1:n), n, eigenvalues, evec, &
                              k = min(ntapers, n))
        else
            call EigValVecSym(Dij(1:n, 1:n), n, eigenvalues, evec)
        endif

    end if

    if (present(ntapers)) then
        numk = min(ntapers, n)
    else
        numk = n
    end if

    if (exclude == 0) then
        tapers(1:n, 1:numk) = evec(1:n, 1:numk)

    else
        do i = 1, n
            tapers(ind(i), 1:numk) = evec(i, 1:numk)

        end do

    end if

    deallocate (dij)
    deallocate (evec)
    if (exclude == 1) deallocate (dijex)

end subroutine SHReturnTapersMap
