!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!      This is a simple Python wrapper for the shtools library written by Mark
!      Wieczorek (http://shtools.ipgp.fr/). It uses the f2py compiler to create
!      the shtools.so shared object that can be opened from python. I tried to
!      make the wrapper guess the array shapes and other optional
!      parameters automatically from the input arrays. This wrapper is far from complete!
!      After installation and setup of "f2py" (usually comes with
!      numpy/scipy)(http://www.scipy.org/F2py) compile with: f2py
!      -I<SHTOOLS_INCLUDE_DIR> -L<SHTOOLS_LIB_DIR> -lSHTOOLS2.8 \ -lfftw3 -lm
!      --f90flags="-m64 -fPIC" --f77flags="-m64 -fPIC" -c -m shtools wrapper.f90
!
!      in python type: >>> import shtools. Examine content with shtools.<TAB>
!      and use ipythons ? operator with shtools.<subroutine>? . You can also use
!      pythons dir(shtools) as well as print shtools.pymakegriddh.__doc__ to get
!      info on input and output variables. Also check out SHTOOLS documentary
!      for more information. We use geodesy + no Condon-Shortley Phase. To
!      change it, you have to change the subroutine below and recompile.
!
!      Authors: Paul Anton Letnes (2011),  Matthias Meschede (2013)
!
!      Functions and subroutines that do not need wrapping:
!          Random.f95: Replace by 'import random' or 'import numpy.random' in python
!          FFTW3.f95:  Seems to be 'internal use only', no need to expose it
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module params
    implicit none

    !configuration variables:
    integer :: norm = 1
    integer :: csphase = 1
end module

!========  LEGENDRE FUNCTIONS ========
    !== Geodesy 4pi normalized ==
    subroutine pyPlmBar(p, lmax, z)
        use shtools, only: PlmBar
        use params
        implicit none
        integer, intent(in) :: lmax
        real(8), intent(out) :: p((lmax + 1)*(lmax + 2) / 2)
        real(8), intent(in) :: z
        !f2py depend(lmax) p
    
        call PlmBar(p, lmax, z, csphase, norm)
    end subroutine
    !missing: PlmBar_d1
    subroutine pyPlBar(p, lmax, z)
        use shtools, only: PlBar
        use params
        implicit none
        integer, intent(in) :: lmax
        real(8), intent(out) :: p(lmax + 1)
        real(8), intent(in) :: z
        !f2py depend(lmax) p
    
        call PlBar(p, lmax, z)
    end subroutine
    
    subroutine pyPlBar_d1(p, dp, lmax, z)
        use shtools, only: PlBar_d1
        use params
        implicit none
        integer, intent(in) :: lmax
        real(8), intent(out), dimension(lmax + 1) :: p, dp
        real(8), intent(in) :: z
        !f2py depend(lmax) p, dp
    
        call PlBar_d1(p, dp, lmax, z)
    end subroutine
    
    !== Orthonormalized ==
    subroutine pyPlmON(p, lmax, z)
        use shtools, only: PlmON
        use params
        implicit none
        integer, intent(in) :: lmax
        real(8), intent(out) :: p((lmax + 1)*(lmax + 2) / 2)
        real(8), intent(in) :: z
        !f2py depend(lmax) p
    
        call PlmON(p, lmax, z, csphase, norm)
    end subroutine
    !PlmON_d1
    subroutine pyPlON(p, lmax, z)
        use shtools, only: PlON
        use params
        implicit none
        integer, intent(in) :: lmax
        real(8), intent(out) :: p(lmax + 1)
        real(8), intent(in) :: z
        !f2py depend(lmax) p
        
        call PlON(p, lmax, z)
    end subroutine
    
    subroutine pyPlON_d1(p, dp, lmax, z)
        use shtools, only: PlON_d1
        use params
        implicit none
        integer, intent(in) :: lmax
        real(8), intent(out), dimension(lmax + 1) :: p, dp
        real(8), intent(in) :: z
        !f2py depend(lmax) p, dp
    
        call PlON_d1(p, dp, lmax, z)
    end subroutine
    
    !== Schmidt normalized ==
    subroutine pyPlmSchmidt(p, lmax, z)
        use shtools, only: PlmSchmidt
        use params
        implicit none
        integer, intent(in) :: lmax
        real(8), intent(out) :: p((lmax + 1)*(lmax + 2) / 2)
        real(8), intent(in) :: z
        !f2py depend(lmax) p
    
        call PlmSchmidt(p, lmax, z, csphase, norm)
    end subroutine
    !missing: PlmSchmidt_d1
    subroutine pyPlSchmidt(p, lmax, z)
        use shtools, only: PlSchmidt
        use params
        implicit none
        integer, intent(in) :: lmax
        real(8), intent(out) :: p(lmax + 1)
        real(8), intent(in) :: z
        !f2py depend(lmax) p
        
        call PlSchmidt(p, lmax, z)
    end subroutine
    
    subroutine pyPlSchmidt_d1(p, dp, lmax, z)
        use shtools, only: PlSchmidt_d1
        use params
        implicit none
        integer, intent(in) :: lmax
        real(8), intent(out), dimension(lmax + 1) :: p, dp
        real(8), intent(in) :: z
        !f2py depend(lmax) p, dp
    
        call PlSchmidt_d1(p, dp, lmax, z)
    end subroutine
    
    !== Unnormalized ==
    !missing: PLegendreA
    subroutine pyPLegendreA_d1(p, dp, lmax, z)
        use shtools, only: PLegendreA_d1
        use params
        implicit none
        integer, intent(in) :: lmax
        real(8), intent(out), dimension((lmax + 1)*(lmax + 2) / 2) :: p, dp
        real(8), intent(in) :: z
        !f2py depend(lmax) p
        !f2py depend(lmax) dp
    
        call PLegendreA_d1(p, dp, lmax, z, csphase)
    end subroutine
    
    subroutine pyPLegendre(p, lmax, z)
        use shtools, only: PLegendre
        use params
        implicit none
        integer, intent(in) :: lmax
        real(8), intent(out) :: p(lmax + 1)
        real(8), intent(in) :: z
        !f2py depend(lmax) p
    
        call PLegendre(p, lmax, z)
    end subroutine
    !missing: Plegendre_d1
    
    !Other
    function pyPlmIndex(l, m)
        use shtools, only: PlmIndex
        use params
        implicit none
        integer :: pyPlmIndex
        integer, intent(in) :: l, m
    
        ! TODO check fortran -> python index issue.
        ! Fortran: 1-based. Python: 0-based
        pyPlmIndex = PlmIndex(l, m) - 1
    end function

!======== SPHERICAL HARMONICS TRANSFORMATIONS ========
    subroutine pySHExpandDH(grid,nlat,nlon,coeffs,lmax_coeffs,sampling)
        use shtools, only: SHExpandDH
        use params
        implicit none
        double precision grid(nlat,nlon)
        double precision coeffs(2,lmax_coeffs+1,lmax_coeffs+1)
        integer nlat,lmax_coeffs,lmax_grid, nlon, sampling
        !f2py intent(out) coeffs
        !f2py intent(in) grid
        !f2py integer intent(hide), optional, depend(grid) :: nlat = shape(grid,0)
        !f2py integer intent(hide), optional, depend(grid) :: nlon = shape(grid,1)
        !f2py integer intent(in), optional, depend(nlat):: lmax_coeffs = nlat/2-1
        !f2py integer optional :: sampling = 2

        call SHExpandDH(grid, nlat, coeffs, lmax_grid, norm, sampling, csphase, lmax_coeffs)
    end subroutine

    subroutine pyMakeGridDH(grid, nlat, coeffs, lmax, sampling, lmax_calc)
        use shtools, only: MakeGridDH
        use params
        implicit none
        double precision :: interval
        integer :: n, nlat, lmax, lmax_calc, sampling
        double precision grid(nlat,sampling*nlat)
        double precision coeffs(2,lmax+1,lmax+1)
        !f2py intent(out) grid
        !f2py intent(in) coeffs
        !f2py integer intent(in), depend(coeffs), optional :: lmax = shape(coeffs,1)-1
        !f2py integer intent(in), depend(lmax),   optional :: nlat = 2*(lmax +1)
        !f2py integer optional :: lmax_calc=lmax
        !f2py integer optional :: sampling = 2
    
        call MakeGridDH(grid, n, coeffs, lmax, norm, sampling, csphase, lmax_calc)
    end subroutine
