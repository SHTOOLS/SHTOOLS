subroutine SHReadJPL(filename, cilm, lmax, error, gm, formatstring, exitstatus)
!------------------------------------------------------------------------------
!
!   This program will read the file format of spherical harmonic
!   coefficients sometimes used by the JPL group into standard arrays for
!   of Cilm and the error. In order to for this to work, you need to know
!   a priori the maximum spherical harmonic degree of the file.
!
!   The input file contains
!       1. comment lines starting with "#"
!       2. GM (optional)
!       3. A list of J_l, which is -C(1,l+1,1)
!       4. A list of the cosine and sine terms
!       5. Starting from 2, the same thing over for the errors.
!
!   Calling Parameters
!
!       IN
!           filename        The name of the file.
!           lmax            Maximum spherical harmonic degree.
!
!       OUT
!           cilm            An array of the spherical harmonic coefficients.
!
!       OPTIONAL
!           error           An array containing the error coefficients.
!           formatstring    This is a string containing an I/O specification
!                           for the numbers of the spherical harmonic 
!                           coefficients.
!           exitstatus  If present, instead of executing a STOP when an error
!                       is encountered, the variable exitstatus will be
!                       returned describing the error.
!                       0 = No errors;
!                       1 = Improper dimensions of input array;
!                       2 = Improper bounds for input variable;
!                       3 = Error allocating memory;
!                       4 = File IO error.
!
!   Copyright (c) 2016, SHTOOLS
!   All rights reserved.
!
!------------------------------------------------------------------------------
    implicit none

    character(*), intent(in) :: filename
    integer, intent(in) :: lmax
    real*8, intent(out) :: cilm(:,:,:)
    real*8, intent(out), optional :: error(:,:,:), gm(2)
    integer, intent(out), optional :: exitstatus
    character, intent(in), optional :: formatstring*6
    real*8 :: gm1, gm2
    logical ::  gmpresent
    integer l, m, stat, i, ll1, mm1, ll2, mm2, skip
    character :: c*4, s*4, j*4, dumb*14, dumb2*14, js*4, cs*4, ss*4, inum*2

    if (present(exitstatus)) exitstatus = 0

    if (size(cilm(:,1,1)) < 2 .or. size(cilm(1,:,1)) < lmax+1 .or. &
            size(cilm(1,1,:)) < lmax+1) then
        print*, "Error --- SHReadJPL"
        print*, "CILM must be dimensioned (2, LMAX+1, LMAX+1) " // &
                "where LMAX is ", lmax
        print*, "Input array is dimensioned ", size(cilm(:,1,1)), &
                size(cilm(1,:,1)), size(cilm(1,1,:))
        if (present(exitstatus)) then
            exitstatus = 1
            return
        else
            stop
        end if

    end if 

    if (present(error)) then
        if (size(error(:,1,1)) < 2 .or. size(error(1,:,1)) < lmax+1 .or. &
                size(error(1,1,:)) < lmax+1) then
            print*, "Error --- SHReadJPL"
            print*, "error must be dimensioned (2, LMAX+1, LMAX+1) " // &
                    "where LMAX is ", lmax
            print*, "Input array is dimensioned ", size(error(:,1,1)), &
                    size(error(1,:,1)), size(error(1,1,:))
            if (present(exitstatus)) then
                exitstatus = 1
                return
            else
                stop
            end if

        end if

    end if

    gmpresent = .false.

    c = "OBAC"
    cs = "SIGC"
    s = "OBAS"
    ss = "SIGS"
    j = "OBAJ"
    js = "SIGJ"

    if (lmax <= 99) then
         inum = "i2"

    else if (lmax >= 100 .and. lmax <=999) then
        inum = "i3"

    else
        inum = "i4"

    end if

    open(13,file=filename)

    ! Skip past header lines
    skip = 0

    do
        read(13,"(a)", iostat=stat) dumb
        if (stat /= 0) then
            print*, "Error --- SHReadJPL"
            print*, "Problem reading header lines of ", filename
            if (present(exitstatus)) then
                exitstatus = 4
                return
            else
                stop
            end if

        end if

        if (dumb(1:1) == "#") then
            skip = skip + 1

        else
            exit

        end if

    end do

    ! Try to read GM, if present
    rewind(13)

    do i = 1, skip, 1
        read(13, *, iostat=stat)

    end do

    read(13, "(a12, f14.6)", iostat=stat) dumb, gm1

    if (stat /= 0) then
        print*, "Error --- SHReadJPL"
        print*, "Problem reading GM(1)"
        if (present(exitstatus)) then
            exitstatus = 4
            return
        else
            stop
        end if

    end if

    if (dumb(2:3) == "GM") then
        if ( present(gm) ) then
            gm(1) = gm1
            gmpresent = .true.

        end if

    else
        rewind(13)

        do i = 1, skip, 1
            read(13,*)
        end do

    end if

    ! Start reading J_l coefficients
    do i = 1, lmax
        if (present(formatstring)) then
            read(13,"(1x, a5, "//inum//", 5x,"//formatstring//")", &
                 iostat=stat) dumb, l, cilm(1,i+1,1)

        else
            read(13,"(1x, a5, "//inum//", 5x, e19.12)", iostat=stat) dumb, &
                l, cilm(1,i+1,1)

        end if

        if (stat /= 0) then
            print*, "Error --- SHReadJPL"
            print*, "Problem reading file", i
            if (present(exitstatus)) then
                exitstatus = 4
                return
            else
                stop
            end if

        else if (dumb(1:4) /= j) then
            print*, "Error --- SHReadJPL"
            print*, "Problem with line specifiers in file ", filename
            print*, "Read sring = ", dumb(1:4)
            print*, "Expected string = ", j
            if (present(exitstatus)) then
                exitstatus = 4
                return
            else
                stop
            end if

        else if (i /= l) then
            print*, "Error --- SHReadJPL"
            print*, "Problem with indices in file ", filename
            print*, "Read indice = ", l
            print*, "Expected indice = ", i
            if (present(exitstatus)) then
                exitstatus = 4
                return
            else
                stop
            end if

        end if

        cilm(1,i+1,1) = -cilm(1,i+1,1)

    end do

    ! Read the Clm and Slm coefficients
    do l = 1, lmax
        do m = 1, l
            if (present(formatstring)) then
                read(13, "(1x, a5, "//inum//", 1x, "//inum//", 5x," &
                    //formatstring//", 5x, a5, "//inum//", 1x, " &
                    //inum//", 5x,"//formatstring//")", &
                    iostat=stat) dumb, ll1, mm1, cilm(1,l+1,m+1), &
                                 dumb2, ll2, mm2, cilm(2,l+1,m+1)

            else
                read(13,"(1x, a5, "//inum//", 1x, " &
                     //inum//", 5x, e19.12, 5x, a5, "//inum//", 1x, " &
                     //inum//", 5x, e19.12)", iostat=stat) dumb, ll1, mm1, &
                                                           cilm(1,l+1,m+1), &
                                                           dumb2, ll2, mm2,&
                                                           cilm(2,l+1,m+1)

            end if

            if (stat /= 0) then
                print*, "Error --- SHReadJPL"
                print*, "Problem reading file", l, m
                if (present(exitstatus)) then
                    exitstatus = 4
                    return
                else
                    stop
                end if

            else if (dumb(1:4) /= c .or. dumb2(1:4) /= s) then
                print*, "Error --- SHReadJPL"
                print*, "Problem with line specifiers in file ", filename
                print*, "Read srings = ", dumb(1:4), dumb2(1:4)
                print*, "Expected strings = ", c, s
                if (present(exitstatus)) then
                    exitstatus = 4
                    return
                else
                    stop
                end if

            else if (ll1 /= l .or. ll2 /= l .or. mm1 /= m .or. mm2 /= m ) then
                print*, "Error --- SHReadJPL"
                print*, "Problem with indices in file ", filename
                print*, "Read indices (l1, m1), (l2, m2) = ", ll1, mm1, ll2, mm2
                print*, "Expected indices (l, m) = ", l, m
                if (present(exitstatus)) then
                    exitstatus = 4
                    return
                else
                    stop
                end if

            end if

        end do

    end do

    ! Next read uncertainties
    if (present(error)) then
        if (gmpresent) then
            read(13, "(a12, f17.6)", iostat=stat) dumb, gm2
            gm(2) = gm2

            if (stat /= 0) then
                print*, "Error --- SHReadJPL"
                print*, "Problem reading GM(2)"
                if (present(exitstatus)) then
                    exitstatus = 4
                    return
                else
                    stop
                end if

            end if

        end if

        do i = 1, lmax
            if (present(formatstring)) then
                read(13,"(1x, a5, "//inum//", 5x,"//formatstring//")", &
                     iostat=stat) dumb, l, error(1,i+1,1)

            else
                read(13,"(1x, a5, "//inum//", 5x, e19.12)", iostat=stat) &
                    dumb, l, error(1,i+1,1)

            end if

            if (stat /= 0) then
                print*, "Error --- SHReadJPL"
                print*, "Problem reading file during J2 error coefficients", i
                if (present(exitstatus)) then
                    exitstatus = 4
                    return
                else
                    stop
                end if

            else if (dumb(1:4) /= js) then
                print*, "Error --- SHReadJPL"
                print*, "Problem with line specifiers in file ", filename
                print*, "Read sring = ", dumb(1:4)
                print*, "Expected string = ", js
                if (present(exitstatus)) then
                    exitstatus = 4
                    return
                else
                    stop
                end if

            else if (i /= l) then
                print*, "Error --- SHReadJPL"
                print*, "Problem with indices in file ", filename
                print*, "Read indice = ", l
                print*, "Expected indice = ", i
                if (present(exitstatus)) then
                    exitstatus = 4
                    return
                else
                    stop
                end if

            end if

        enddo

        do l = 1, lmax
            do m = 1, l
                if (present(formatstring)) then
                    read(13, "(1x, a5, "//inum//", 1x, "//inum//", 5x," &
                         //formatstring//", 5x, a5, "//inum//", 1x, " &
                         //inum//", 5x,"//formatstring//")", iostat=stat) &
                             dumb, ll1, mm1, error(1,l+1,m+1), dumb2, ll2, &
                             mm2, error(2,l+1,m+1)

                else
                    read(13,"(1x, a5, "//inum//", 1x, " &
                         //inum//", 5x, e19.12, 5x, a5, "//inum//", 1x, " &
                         //inum//", 5x, e19.12)", iostat=stat) dumb, ll1, mm1,&
                             error(1,l+1,m+1), dumb2, ll2, mm2, error(2,l+1,m+1)

                end if

                if (stat /= 0) then
                    print*, "Error --- SHReadJPL"
                    print*, "Problem reading file during error coefficients", l, m
                    if (present(exitstatus)) then
                        exitstatus = 4
                        return
                    else
                        stop
                    end if

                else if (dumb(1:4) /= cs .or. dumb2(1:4) /= ss) then
                    print*, "Error --- SHReadJPL"
                    print*, "Problem with line specifiers in file ", filename
                    print*, "Read srings = ", dumb(1:4), dumb2(1:4)
                    print*, "Expected strings = ", cs, ss
                    if (present(exitstatus)) then
                        exitstatus = 4
                        return
                    else
                        stop
                    end if

                else if (ll1 /= l .or. ll2 /= l .or. mm1 /= m .or. mm2 /= m ) then
                    print*, "Error --- SHReadJPL"
                    print*, "Problem with indices in file ", filename
                    print*, "Read indices (l1, m1), (l2, m2) = ", ll1, mm1, ll2, mm2
                    print*, "Expected indices (l, m) = ", l, m
                    if (present(exitstatus)) then
                        exitstatus = 4
                        return
                    else
                        stop
                    end if

                end if

            end do

        end do

    end if

    close(13)

end subroutine SHReadJPL
