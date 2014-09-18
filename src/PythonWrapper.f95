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

module pyshtools
    !can contain global module variables 
    !(e.g. csphase and cnorm could be moved here)
    implicit none
contains

!========  LEGENDRE FUNCTIONS ========

    !== Geodesy 4pi normalized ==
    subroutine PlmBar(p, lmax, z, csphase, cnorm)
        use shtools, only: PlmBar_f => PlmBar
        implicit none
        integer, intent(in) :: csphase, cnorm
        integer, intent(in) :: lmax
        real(8), intent(out) :: p((lmax + 1)*(lmax + 2) / 2)
        real(8), intent(in) :: z
        !f2py depend(lmax) p
        !f2py integer, optional :: csphase=1
        !f2py integer, optional :: cnorm=0
    
        call PlmBar_f(p, lmax, z, csphase, cnorm)
    end subroutine PlmBar
    !PlmBar_d1
    subroutine PlBar(p, lmax, z)
        use shtools, only: PlBar_f => PlBar
        implicit none
        integer, intent(in) :: lmax
        real(8), intent(out) :: p(lmax + 1)
        real(8), intent(in) :: z
        !f2py depend(lmax) p
    
        call PlBar_f(p, lmax, z)
    end subroutine PlBar
    
    subroutine PlBar_d1(p, dp, lmax, z)
        use shtools, only: PlBar_d1_f => PlBar_d1
        implicit none
        integer, intent(in) :: lmax
        real(8), intent(out), dimension(lmax + 1) :: p, dp
        real(8), intent(in) :: z
        !f2py depend(lmax) p, dp
    
        call PlBar_d1_f(p, dp, lmax, z)
    end subroutine PlBar_d1
    
    !== Orthonormalized ==
    subroutine PlmON(p, lmax, z, csphase, cnorm)
        use shtools, only: PlmON_f => PlmON
        implicit none
        integer, intent(in) :: lmax
        real(8), intent(out) :: p((lmax + 1)*(lmax + 2) / 2)
        real(8), intent(in) :: z
        integer, intent(in) :: csphase, cnorm
        !f2py depend(lmax) p
        !f2py integer, optional :: csphase=1
        !f2py integer, optional :: cnorm=0
    
        call PlmON_f(p, lmax, z, csphase, cnorm)
    end subroutine PlmOn
    !PlmON_d1
    subroutine PlON(p, lmax, z)
        use shtools, only: PlON_f => PlON
        implicit none
        integer, intent(in) :: lmax
        real(8), intent(out) :: p(lmax + 1)
        real(8), intent(in) :: z
        !f2py depend(lmax) p
        
        call PlON_f(p, lmax, z)
    end subroutine PlON
    
    subroutine PlON_d1(p, dp, lmax, z)
        use shtools, only: PlON_d1_f => PlON_d1
        implicit none
        integer, intent(in) :: lmax
        real(8), intent(out), dimension(lmax + 1) :: p, dp
        real(8), intent(in) :: z
        !f2py depend(lmax) p, dp
    
        call PlON_d1_f(p, dp, lmax, z)
    end subroutine PlON_d1
    
    !== Schmidt normalized ==
    subroutine PlmSchmidt(p, lmax, z, csphase, cnorm)
        use shtools, only: PlmSchmidt_f => PlmSchmidt
        implicit none
        integer, intent(in) :: csphase, cnorm
        integer, intent(in) :: lmax
        real(8), intent(out) :: p((lmax + 1)*(lmax + 2) / 2)
        real(8), intent(in) :: z
        !f2py depend(lmax) p
        !f2py integer, optional :: csphase=1
        !f2py integer, optional :: cnorm=0
    
        call PlmSchmidt_f(p, lmax, z, csphase, cnorm)
    end subroutine PlmSchmidt
    !PlmSchmidt_d1
    subroutine PlSchmidt(p, lmax, z)
        use shtools, only: PlSchmidt_f => PlSchmidt
        implicit none
        integer, intent(in) :: lmax
        real(8), intent(out) :: p(lmax + 1)
        real(8), intent(in) :: z
        !f2py depend(lmax) p
        
        call PlSchmidt_f(p, lmax, z)
    end subroutine PlSchmidt
    
    subroutine PlSchmidt_d1(p, dp, lmax, z)
        use shtools, only: PlSchmidt_d1_f => PlSchmidt_d1
        implicit none
        integer, intent(in) :: lmax
        real(8), intent(out), dimension(lmax + 1) :: p, dp
        real(8), intent(in) :: z
        !f2py depend(lmax) p, dp
    
        call PlSchmidt_d1_f(p, dp, lmax, z)
    end subroutine PlSchmidt_d1
    
    !== Unnormalized ==
    !PLegendreA
    subroutine PLegendreA_d1(p, dp, lmax, z, csphase)
        use shtools, only: PLegendreA_d1_f => PLegendreA_d1
        implicit none
        integer, intent(in) :: lmax
        real(8), intent(out), dimension((lmax + 1)*(lmax + 2) / 2) :: p, dp
        !real(8), intent(out), dimension((lmax + 1)*(lmax + 2) / 2) :: dp
        real(8), intent(in) :: z
        integer, intent(in) :: csphase
        !f2py depend(lmax) p
        !f2py depend(lmax) dp
        !f2py integer, optional :: csphase=1
    
        call PLegendreA_d1_f(p, dp, lmax, z, csphase)
    end subroutine PLegendreA_d1
    
    subroutine PLegendre(p, lmax, z)
        use shtools, only: PLegendre_f => PLegendre
        implicit none
        integer, intent(in) :: lmax
        real(8), intent(out) :: p(lmax + 1)
        real(8), intent(in) :: z
        !f2py depend(lmax) p
    
        call PLegendre_f(p, lmax, z)
    end subroutine PLegendre
    !Plegendre_d1
    
    !Other
    function PlmIndex(l, m)
        use shtools, only: PlmIndex_f => PlmIndex
        implicit none
        integer :: PlmIndex
        integer, intent(in) :: l, m
    
        ! TODO check fortran -> python index issue.
        ! Fortran: 1-based. Python: 0-based
        PlmIndex = PlmIndex_f(l, m) - 1
    end function PlmIndex

!======== SPHERICAL HARMONICS TRANSFORMATIONS ========



end module
