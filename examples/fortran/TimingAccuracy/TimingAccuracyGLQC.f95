Program TimingAccuracyGLQC
!------------------------------------------------------------------------------
!
!    This program will test the accuracy of the spherical harmonic GLQ
!    tranformation routines by expanding a field in the space domain,
!    transforming this to spherical harmonics, and comparing the relative error
!    of the coefficients.
!
!    Copyright (c) 2005, SHTOOLS
!    All rights reserved.
!
!------------------------------------------------------------------------------
    use SHTOOLS
    use ftypes

    implicit none

    integer(int32), parameter :: maxdeg = 2800
    character(200) :: outfile1, outfile2, outfile3, outfile4, outfile, infile
    complex(dp), allocatable :: cilm(:,:,:), cilm2(:,:,:), gridglq(:,:)
    real(dp), allocatable :: zero(:), w(:)
    real(dp) :: maxerror, err1, err2, beta, rms, timein(3), timeout(3)
    integer(int32) :: lmax, l, m, seed, astat(5)

    allocate(cilm(2,maxdeg+1,maxdeg+1), stat=astat(1))
    allocate(cilm2(2,maxdeg+1,maxdeg+1), stat=astat(2))
    allocate(zero(maxdeg+1), stat=astat(3))
    allocate(gridglq(maxdeg+1,2*maxdeg+1), stat=astat(4))
    allocate(w(maxdeg+1), stat=astat(5))

    if (sum(astat(1:5)) /= 0) then
        print*, "Error --- TimingAccuracyGLQ"
        print*, "Problem allocating arrays CILM, CILM2, ZERO, GRIDGLQ, W"
        stop
    end if

    ! A data input file may be passed as the first argument. Otherwise prompt for required settings.
    if (command_argument_count() > 0) then
        call get_command_argument(1, infile)
        open(unit=20, file=infile, action="read")
        read(20,*) beta, outfile
        close(20)
    else
        print*, "Value of beta for power law Sff = l^(-beta) > "
        read(*,*) beta

        print*, "output file name > "
        read(*,*) outfile
    end if

    outfile1 = trim(outfile) // ".timef"
    outfile2 = trim(outfile) // ".timei"
    outfile3 = trim(outfile) // ".maxerror"
    outfile4 = trim(outfile) // ".rmserror"

    seed = -1053253

    cilm = cmplx(0.0_dp, 0.0_dp, dp)

    do l = 1, maxdeg

        do m = 0, l
            if (m == 0) then
                cilm(1,l+1,m+1) = cmplx(RandomGaussian(seed), &
                                        RandomGaussian(seed), dp)
            else
                cilm(1,l+1,m+1) = cmplx(RandomGaussian(seed), &
                                        RandomGaussian(seed), dp)
                cilm(2,l+1,m+1) = cmplx(RandomGaussian(seed), &
                                        RandomGaussian(seed), dp)
            end if
        end do
        cilm(1:2, l+1, 1:l+1) = cilm(1:2, l+1, 1:l+1) * sqrt(dble(l)**beta) &
                                / sqrt(2.0_dp * l + 1)
    end do

    print*, "Lmax, Maximum rel. error of Cilm, RMS relative error, Precompute time (sec), Time inverse (sec), Time forward (sec)"

    lmax = 1

    do

        lmax = lmax * 2
        if (lmax > maxdeg) lmax = maxdeg

        if (lmax == 2) then
            open(12,file=outfile1, status="replace")
            open(13,file=outfile2, status="replace")
            open(14,file=outfile3, status="replace")
            open(15,file=outfile4, status="replace")
        else
            open(12,file=outfile1, position="append")
            open(13,file=outfile2, position="append")
            open(14,file=outfile3, position="append")
            open(15,file=outfile4, position="append")
        end if

        call cpu_time(timein(1))
        call SHGLQ(lmax, zero(1:lmax+1), w(1:lmax+1))
        call cpu_time(timeout(1))

        call cpu_time(timein(2))
        call MakeGridGLQC(gridglq(1:lmax+1, 1:2*lmax+1), &
                          cilm(1:2,1:lmax+1, 1:lmax+1), lmax, &
                          zero=zero(1:lmax+1), norm=1)
        call cpu_time(timeout(2))

        call cpu_time(timein(3))
        call SHExpandGLQC(cilm2(1:2, 1:lmax+1, 1:lmax+1), lmax, &
                          gridglq(1:lmax+1, 1:2*lmax+1), w, &
                          zero=zero(1:lmax+1), norm=1)
        call cpu_time(timeout(3))

        maxerror = 0.0_dp
        rms = 0.0_dp

        do l = 1, lmax

            do m = 0, l
                if (m == 0) then
                    err1 = abs( (cilm(1,l+1,m+1) - cilm2(1,l+1,m+1)) ) / abs( cilm(1,l+1,m+1) )
                    if (err1 >= maxerror) maxerror = err1
                    rms = rms + err1**2
                else
                    err1 = abs( (cilm(1,l+1,m+1) - cilm2(1,l+1,m+1)) ) / abs( cilm(1,l+1,m+1) )
                    err2 = abs( (cilm(2,l+1,m+1) - cilm2(2,l+1,m+1)) ) / abs( cilm(2,l+1,m+1) )
                    if (err1 >= maxerror) maxerror = err1
                    if (err2 >= maxerror) maxerror = err2
                    rms = rms + err1**2 + err2**2
                end if
            end do
        end do
        rms = sqrt(rms / dble(l+1)**2)

        ! elasped time in seconds!
        print*, lmax, maxerror, rms, timeout(1)-timein(1), timeout(2)-timein(2), timeout(3)-timein(3)
        write(12,*) lmax, timeout(2)-timein(2)
        write(13,*) lmax, timeout(3)-timein(3)
        write(14,*) lmax, maxerror
        write(15,*) lmax, rms

        if (maxerror > 100.0_dp) then
            print*, "TESTS FAILED"
            print*, "Degree = ", lmax
            print*, "Maximum relative error = ", maxerror
            close(12)
            close(13)
            close(14)
            close(15)
            stop
        end if

        close(12)
        close(13)
        close(14)
        close(15)

        if (lmax == maxdeg) exit

    enddo

end program TimingAccuracyGLQC
