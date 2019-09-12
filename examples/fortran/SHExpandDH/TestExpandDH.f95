Program TestExpandDH
!------------------------------------------------------------------------------
!
!    This program will expand a grid containing n samples in both longitude and
!    latitude into spherical harmonics using the sampling theorems given in
!    Driscoll and Heally (1994).
!
!    Copyright (c) 2005, SHTOOLS
!    All rights reserved.
!
!------------------------------------------------------------------------------
    use SHTOOLS
    use ftypes

    implicit none
    integer, parameter :: maxdeg = 120, maxgrid = 800
    character(80) :: infile
    real(dp) :: cilm(2, maxdeg+1, maxdeg+1), grid(maxgrid, maxgrid), &
                griddh(maxgrid, maxgrid), cilm2(2, maxdeg+1, maxdeg+1), &
                interval, error, maxerror
    integer :: lmax, n, nlat, nlong, lmax2, l, m, i

    infile = "../../ExampleDataFiles/MarsTopo719.shape"

    call SHRead(infile, cilm, lmax)
    print*, "Lmax of data file = ", lmax

    !--------------------------------------------------------------------------
    !
    !     Create an equally sampled grid
    !
    !--------------------------------------------------------------------------

    lmax = 89
    print*, "For this test of SHExpandDH, the initial grid will be created with lmax = ", lmax

    interval = 1.0_dp

    print*, "Making grid using MakeGrid2D..."

    call MakeGrid2D(grid, cilm, lmax, interval, nlat, nlong)
    print*, "nlat, nlong = ", nlat, nlong

    print*, "Minimum and maximum values of the input grid = ", minval(grid(1:nlat, 1:nlong)), maxval(grid(1:nlat,1:nlong))

    n = nlat -1

    print*, "Number of samples in latitude used with SHExpandDH = ", n

    griddh(1:n,1:2*n) = grid(1:n, 1:2*n)

    !--------------------------------------------------------------------------
    !
    !    Expand equally spaced grid using Driscoll and Heally routine, and
    !    compute relative error between input and output spherical harmonic
    !    coefficients.
    !
    !--------------------------------------------------------------------------

    call SHExpandDH(griddh, n, cilm2, lmax2, sampling = 2)
    print*, "Lmax of Driscoll Heally expansion = ", lmax2

    maxerror = 0.0_dp

    do l = 0, lmax
        error = (cilm2(1,l+1,1)-cilm(1,l+1,1)) / cilm(1,l+1,1)
        if (error > maxerror) maxerror = error
        do m = 1, l, 1
            error = (cilm2(1,l+1,m+1)-cilm(1,l+1,m+1)) / cilm(1,l+1,m+1)
            if (error > maxerror) maxerror = error
            error = (cilm2(2,l+1,m+1)-cilm(2,l+1,m+1)) / cilm(2,l+1,m+1)
            if (error > maxerror) maxerror = error
        end do
    end do

    print*, "Maximum relative error between initial and output spherical harmonic coefficients = ", maxerror

    !--------------------------------------------------------------------------
    !
    !    Re-expand coefficients in the space domain, and compare
    !    with input grid.
    !
    !--------------------------------------------------------------------------

    griddh = 0.0_dp
    call MakeGridDH(griddh, n, cilm2, lmax2)
    print*, "Number of samples in latitude and longitude = ", n

    maxerror = 0.0_dp

    do i = 1, n
        error = maxval( abs( griddh(1:n,i) - grid(1:n, 2*i-1) ) )
        if (error > maxerror) maxerror = error
    end do

    print*, "Maximum error (meters) in the space domain = ", maxerror

end program TestExpandDH
