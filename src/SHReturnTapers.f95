subroutine SHReturnTapers(theta0, lmax, tapers, eigenvalues, taper_order, &
                          degrees, exitstatus)
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
!
!       OUT
!           tapers          An (lmax+1) by (lmax+1)**2 array containing
!                           all the eigenfunctions of the space-
!                           concentration kernel. Eigenfunctions
!                           are listed by columns in decreasing order
!                           corresponding to value of their eigenvalue.
!           eigenvalues     A vector of length lmax+1 containing the
!                           eigenvalued corresponding to the individual
!                           eigenfunctions.
!           taper_order     A vector of dimension (lmax+1)**2 denoting which
!                           order m corresponds to the column of tapers and
!                           eigenvalues.
!
!       OPTIONAL (IN)
!           degrees     Specify those degrees of the coupling matrix to
!                       compute. If degrees(l+1) is zero, degree l will not
!                       be computed.
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
    use SHTOOLS, only: SHReturnTapersM
    use ftypes

    implicit none

    real(dp), intent(in) :: theta0
    integer(int32), intent(in) :: lmax
    real(dp), intent(out) :: tapers(:,:), eigenvalues(:)
    integer(int32), intent(out) :: taper_order(:)
    integer(int32), intent(in), optional :: degrees(:)
    integer(int32), intent(out), optional :: exitstatus
    integer(int32) :: m, nt, nt2, n, i, j, jj(1), astat(8)
    real(dp), allocatable :: eval(:), evec(:, :), tapers_unordered(:,:), &
                             eval_unordered(:)
    integer(int32), allocatable :: m_unordered(:)

    if (present(exitstatus)) exitstatus = 0

    if (size(tapers(:,1)) < (lmax+1) .or. size(tapers(1,:)) < (lmax+1)**2) then
        print*, "Error --- SHReturnTapers"
        print*, "TAPERS must be dimensioned as ( LMAX+1, (LMAX+1)**2 ) " // &
                "where LMAX is ", lmax
        print*, "Input array is dimensioned as ", size(tapers(:,1)), &
                size(tapers(1,:))
        if (present(exitstatus)) then
            exitstatus = 1
            return
        else
            stop
        end if

    else if (size(eigenvalues) < (lmax+1)**2) then
        print*, "Error --- SHReturnTapers"
        print*, "EIGENVALUES must be dimensioned as (LMAX+1)**2 " // &
                "where LMAX is ", lmax
        print*, "Input array is dimensioned as ", size(eigenvalues)
        if (present(exitstatus)) then
            exitstatus = 1
            return
        else
            stop
        end if

    else if (size(taper_order) < (lmax+1)**2) then
        print*, "Error --- SHReturnTapers"
        print*, "TAPER_ORDER must be dimensioned as (LMAX+1)**2 " // &
                "where LMAX is ", lmax
        print*, "Input array is dimensioned as ", size(taper_order)
        if (present(exitstatus)) then
            exitstatus = 1
            return
        else
            stop
        end if

    end if

    if (present(degrees)) then
        if (size(degrees) < lmax+1) then
            print*, "Error --- SHReturnTapers"
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

    allocate (eval(lmax+1), stat = astat(1))
    allocate (evec(lmax+1, lmax+1), stat = astat(2))
    allocate (tapers_unordered(lmax+1, (lmax+1)**2), stat = astat(3))
    allocate (eval_unordered((lmax+1)**2), stat = astat(4))
    allocate (m_unordered((lmax+1)**2), stat = astat(5))

    if (astat(1) /= 0 .or. astat(2) /= 0 .or. astat(3) /= 0 .or. &
            astat(4) /= 0 .or. astat(5) /= 0) then
        print*, "Error ---- SHReturnTapers"
        print*, "Problem allocating memory for temporary arrays: ", astat(1:5)
        if (present(exitstatus)) then
            exitstatus = 3
            return
        else
            stop
        end if

    end if

    tapers = 0.0_dp
    tapers_unordered = 0.0_dp
    eigenvalues = 0.0_dp
    eval_unordered = 0.0_dp
    m_unordered = 0
    taper_order = 0

    !----------------------------------------------------------------------
    !
    !   Calculate eigenfunctions and eigenvalues using SHReturnTapersM.
    !   When possible, this routine will use the kernel of Grunbaum et al.
    !   (1982).
    !
    !----------------------------------------------------------------------
    nt = 0

    do m = 0, lmax

        if (present(exitstatus)) then
            if (present(degrees)) then
                call SHReturnTapersM(theta0, lmax, m, evec, eval, ntapers=n, &
                                     degrees=degrees, exitstatus=exitstatus)
            else
                call SHReturnTapersM(theta0, lmax, m, evec, eval, ntapers=n, &
                                     exitstatus=exitstatus)
            end if

            if (exitstatus /= 0) return

        else
            if (present(degrees)) then
                call SHReturnTapersM(theta0, lmax, m, evec, eval, ntapers=n, &
                                     degrees=degrees)
            else
                call SHReturnTapersM(theta0, lmax, m, evec, eval, ntapers=n)
            end if

        end if

        if (n /= 0) then
            tapers_unordered(1:lmax+1, nt+1:nt+n) = evec(1:lmax+1,1:n)
            m_unordered(nt+1:nt+n) = m
            eval_unordered(nt+1:nt+n) = eval(1:n)

            nt = nt + n

        end if

    end do

    !--------------------------------------------------------------------------
    !
    !   Reorder tapers and eigenvalues, and add negative angular orders tapers.
    !
    !--------------------------------------------------------------------------
    nt2 = 0

    do i = 1, nt
        jj = maxloc(eval_unordered(1:nt))
        j = jj(1)

        if (m_unordered(j) == 0) then
            nt2 = nt2 + 1
            taper_order(nt2) = m_unordered(j)
            eigenvalues(nt2) = eval_unordered(j)
            tapers(1:lmax+1,nt2) = tapers_unordered(1:lmax+1,j)

        else
            nt2 = nt2 + 1
            taper_order(nt2) = -m_unordered(j)
            eigenvalues(nt2) = eval_unordered(j)
            tapers(1:lmax+1,nt2) = tapers_unordered(1:lmax+1,j)
            nt2 = nt2 + 1
            taper_order(nt2) = m_unordered(j)
            eigenvalues(nt2) = eval_unordered(j)
            tapers(1:lmax+1,nt2) = tapers_unordered(1:lmax+1,j)

        end if

        eval_unordered(j) = -1.0d25

    end do

    deallocate (eval)
    deallocate (evec)
    deallocate (tapers_unordered)
    deallocate (eval_unordered)
    deallocate (m_unordered)

end subroutine SHReturnTapers
