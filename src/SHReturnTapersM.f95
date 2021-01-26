subroutine SHReturnTapersM(theta0, lmax, m, tapers, eigenvalues, shannon, &
                           degrees, ntapers, exitstatus)
!------------------------------------------------------------------------------
!
!   This subroutine will return all the eigenvalues and eigenfunctions for the
!   space-concentration problem of a spherical cap of angular radius theta0.
!   The returned eigenfunctions correspond to "geodesy" normalized spherical
!   harmonic coefficients, and the eigenfunctions are further normalized such
!   that they have unit power (i.e., the integral of the function squared over
!   the sphere divided by 4 pi is 1, and the sum of the squares of their
!   coefficients is 1). If the optional vector DEGREES is specified, then the
!   eigenfunctions will be computed using only those degrees where DEGREES(l+1)
!   is not zero.
!
!   When possible, the eigenfunctions are calculated using the kernel of
!   Grunbaum et al. 1982 and the eigenvalues are then calculated by integration
!   using the definition of the space-concentration problem. Use of the
!   Grunbaum et al. kernel is prefered over the space-concentration kernel as
!   the eigenfunctions of the later are unreliable when there are several
!   eigenvalues identical (within machine precision) to either 1 or zero. If,
!   the optional parameter DEGREES is specified, and at least one element is
!   zero for degrees greater or equal to abs(m), then the eigenfunctions and
!   eigenvalues will instead be computed directly using the space-concentration
!   kernel.
!
!   Calling Parameters
!
!       IN
!           theta0          Angular radius of spherical cap in RADIANS.
!           lmax            Maximum spherical harmonic degree
!                           for the concentration problem.
!           m               Angular order of the concentration
!                           problem (m=0 corresponds to isotropic case).
!
!       OUT
!           tapers          An (lmax+1) by (lmax+1-abs(m)) array containing
!                           all the eigenfunctions of the space-
!                           concentration kernel. Eigenfunctions
!                           are listed by columns in decreasing order
!                           corresponding to value of their eigenvalue.
!                           Only the first ntapers columns are non-zero.
!           eigenvalues     A vector of length lmax+1-abs(m) containing the
!                           eigenvalued corresponding to the individual
!                           eigenfunctions.
!
!       OPTIONAL (IN)
!           degrees     Specify those degrees of the coupling matrix to
!                       compute. If degrees(l+1) is zero, degree l will not
!                       be computed.
!
!       OPTIONAL (OUT)
!           shannon     Shannon number as calculated from the trace of the
!                       kernel.
!           ntapers     The number of non-zero tapers.
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
    use SHTOOLS, only: ComputeDG82, ComputeDm, EigValVecSymTri, EigValVecSym, &
                       PreGLQ, PlmBar, PlmIndex
    use ftypes

    implicit none

    real(dp), intent(in) :: theta0
    integer(int32), intent(in) :: lmax, m
    real(dp), intent(out) :: tapers(:,:), eigenvalues(:)
    real(dp), intent(out), optional :: shannon
    integer(int32), intent(in), optional :: degrees(:)
    integer(int32), intent(out), optional :: ntapers
    integer(int32), intent(out), optional :: exitstatus
    integer(int32) :: l, n, n_int, j, i, astat(3), ind(lmax+1), use_dg82
    real(dp) :: eval(lmax+1), pi, upper, lower, zero(lmax+1), w(lmax+1), h
    real(dp), allocatable :: evec(:,:), dllmtri(:,:), p(:), dllm(:,:)

    use_dg82 = 1

    if (present(exitstatus)) exitstatus = 0

    if (size(tapers(:,1)) < (lmax+1) .or. &
            size(tapers(1,:)) < (lmax+1-abs(m)) ) then
        print*, "Error --- SHReturnTapersM"
        print*, "TAPERS must be dimensioned as (LMAX+1, LMAX+1-ABS(M)) " // &
                "where LMAX is ", lmax
        print*, "Input array is dimensioned as ", size(tapers(:,1)), &
                size(tapers(1,:))
        if (present(exitstatus)) then
            exitstatus = 1
            return
        else
            stop
        end if

    else if(size(eigenvalues) < (lmax+1-abs(m)) ) then
        print*, "Error --- SHReturnTapersM"
        print*, "EIGENVALUES must be dimensioned as (LMAX+1-ABS(m)) " // &
                "where LMAX is ", lmax
        print*, "Input array is dimensioned as ", size(eigenvalues)
        if (present(exitstatus)) then
            exitstatus = 1
            return
        else
            stop
        end if

    else if (m > lmax) then
        print*, "Error --- SHReturnTapersM"
        print*, "M must be less than or equal to LMAX."
        print*, "M = ", m
        print*, "LMAX = ", lmax
        if (present(exitstatus)) then
            exitstatus = 2
            return
        else
            stop
        end if

    end if

    if (present(degrees)) then
        if (size(degrees) < lmax+1) then
            print*, "Error --- SHReturnTapersM"
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

    if (present(degrees)) then
        if (sum(degrees(1+abs(m):lmax+1)) /= lmax + 1 - abs(m)) use_dg82 = 0
    end if

    allocate (evec(lmax+1, lmax+1-abs(m)), stat = astat(1))
    allocate (dllmtri(lmax+1, lmax+1), stat = astat(2))
    allocate (p((lmax+1)*(lmax+2)/2), stat = astat(3))

    if (astat(1) /= 0 .or. astat(2) /= 0 .or. astat(3) /= 0) then
        print*, "Error --- SHReturnTapersM"
        print*, "Problem allocating arrays EVEC, DLLMTRI, and P", &
            astat(1), astat(2), astat(3)
        if (present(exitstatus)) then
            exitstatus = 3
            return
        else
            stop
        end if

    end if

    pi = acos(-1.0_dp)

    tapers = 0.0_dp
    eigenvalues = 0.0_dp
    eval = 0.0_dp
    evec = 0.0_dp

    if (use_dg82 == 1) then
        !----------------------------------------------------------------------
        !
        !   Calculate space-concentration kernel and the corresponding
        !   eigenfunctions of the Grunbaum et al. kernel.
        !
        !----------------------------------------------------------------------
        n = lmax + 1 - abs(m)

        if (present(exitstatus)) then
            call ComputeDG82(dllmtri(1:n,1:n), lmax, m, theta0, &
                             exitstatus = exitstatus)
        if (exitstatus /= 0) return
            call EigValVecSymTri(dllmtri(1:n,1:n), n, eval(1:n), &
                                 tapers(1+abs(m):lmax+1,1:n), &
                                 exitstatus = exitstatus)
            if (exitstatus /= 0) return
        else
            call ComputeDG82(dllmtri(1:n,1:n), lmax, m, theta0)
            call EigValVecSymTri(dllmtri(1:n,1:n), n, eval(1:n), &
                                 tapers(1+abs(m):lmax+1,1:n))
        end if

        !----------------------------------------------------------------------
        !
        !   Calculate true eigenvalues
        !
        !----------------------------------------------------------------------
        upper = 1.0_dp
        lower = cos(theta0)
        n_int = lmax + 1

        if (present(exitstatus)) then
            call PreGLQ(lower, upper, n_int, zero, w, exitstatus = exitstatus)
            if (exitstatus /= 0) return
        else
            call PreGLQ(lower, upper, n_int, zero, w)
        end if

        do i=1, n_int
            if (present(exitstatus)) then
                call PlmBar(p, lmax, zero(i), exitstatus = exitstatus)
                if (exitstatus /= 0) return
            else
                call PlmBar(p, lmax, zero(i))
            end if

            do j = 1, n
                h = 0.0_dp

                do l = abs(m), lmax
                    h = h + p(PlmIndex(l, abs(m))) * tapers(l+1, j)

                end do

                eigenvalues(j) = eigenvalues(j) + w(i) * h**2

            end do
        end do

        if (m == 0) then
            eigenvalues(1:n) = eigenvalues(1:n) / 2.0_dp

        else
            eigenvalues(1:n) = eigenvalues(1:n) / 4.0_dp

        end if

    else
        !----------------------------------------------------------------------
        !
        !   Calculate space-concentration kernel and the corresponding
        !   eigenfunctions.
        !
        !----------------------------------------------------------------------
        ind = 0
        i = 0

        do l=abs(m), lmax, 1
            if (degrees(l+1) /= 0) then
                i = i + 1
                ind(i) = l + 1
            end if
        end do

        n = i

        if (n /= 0) then

            allocate (dllm(n, n), stat = astat(1))
            if (astat(1) /= 0) then
                print*, "Error --- SHReturnTapersM"
                print*, "Problem allocating array DLLM ", astat(1)
                if (present(exitstatus)) then
                    exitstatus = 3
                    return
                else
                    stop
                end if
            end if

            dllm = 0.0_dp

            if (present(exitstatus)) then
                call ComputeDm(dllmtri, lmax, m, theta0, degrees=degrees, &
                               exitstatus=exitstatus)
                if (exitstatus /= 0) return
            else
                call ComputeDm(dllmtri, lmax, m, theta0, degrees=degrees)
            end if

            do i=1, n
                do j=1, n
                    dllm(i, j) = dllmtri(ind(i), ind(j))
                end do
            end do

            if (present(exitstatus)) then
                call EigValVecSym(dllm(1:n, 1:n), n, &
                                  eval(1:n), evec(1:n,1:n), &
                                exitstatus = exitstatus)
                if (exitstatus /= 0) return
            else
                call EigValVecSym(dllm, n, eval(1:n),&
                                  evec(1:n,1:n))
            end if

            do i=1, n
                tapers(ind(i), 1:n) = evec(i, 1:n)
            end do

            eigenvalues(1:n) = eval(1:n)

        end if

    end if

    if (present(shannon)) then
        shannon = sum(eigenvalues(1:lmax+1-abs(m)))
    end if

    if (present(ntapers)) then
        ntapers = n
    end if

    ! deallocate memory
    call PlmBar (p, -1, zero(1))
    deallocate (evec)
    deallocate (dllmtri)
    deallocate (p)
    if (use_dg82 == 0 .and. n /= 0) deallocate (dllm)

end subroutine SHReturnTapersM
