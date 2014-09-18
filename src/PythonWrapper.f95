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

!==== GLOBAL CONFIGURATION VARIABLES ====
module params
    implicit none
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

!missing: PlmON_d1

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
!== Equally Sampled Grids ==
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

!missing: SHExpandDHC
!missing: MakeGridDHC
!missing: DHaj

!== Gauss-Legendre quadrature grids ==

!missing: PreCompute
!missing: SHExpandGLQ
!missing: MakeGridGLQ
!missing: SHExpandGLQC
!missing: MakeGridGLQC
!missing: GLQGridCoord
!missing: PreGLQ
!missing: NGLQ
!missing: NGLQSH
!missing: NGLQSHN

!== Other ==

!missing: SHExpandLSQ
!missing: MakeGrid2D
subroutine pyMakeGridPoint(coeffs,lmax_calc,lmax_coeffs,lat,lon,return_value)
    use shtools, only: MakeGridPoint
    implicit none
    integer lmax_calc,lmax_coeffs
    double precision coeffs(2,lmax_coeffs+1,lmax_coeffs+1)
    double precision lat, lon
    double precision return_value
    !f2py intent(in) coeffs
    !f2py intent(in) lat
    !f2py intent(in) lon
    !f2py intent(out) return_value
    !f2py integer intent(in),depend(coeffs),optional :: lmax_coeffs = shape(coeffs,1)-1
    !f2py integer intent(in),depend(lmax_coeffs),optional :: lmax_calc = lmax_coeffs
    return_value = MakeGridPoint(coeffs, lmax_calc, lat, lon)
end subroutine

subroutine pySHMultiply(shout, sh1, lmax1, sh2, lmax2, precomp)
    use shtools, only: SHMultiply
    use params
    implicit none
    double precision shout(2,lmax1+lmax2+1, lmax1+lmax2+1)
    double precision sh1(2,lmax1+1,lmax1+1)
    double precision sh2(2,lmax2+1,lmax2+1)
    integer lmax1,lmax2,precomp
    !f2py intent(out) shout
    !f2py intent(in) sh1, sh2
    !f2py integer intent(in), depend(sh1), optional :: lmax1 = shape(sh1,1)-1
    !f2py integer intent(in), depend(sh2), optional :: lmax2 = shape(sh2,1)-1
    !f2py integer, optional :: precomp = 0
    call SHMultiply(shout, sh1, lmax1, sh2, lmax2, precomp, norm, csphase)
end subroutine

!======== SPHERICAL HARMONIC I/O, STORAGE, AND CONVERSIONS ========

!== Spherical harmonic I/O ==
!missing: SHRead
!missing: SHRead2
!missing: SHReadJPL

!== Spherical harmonic storage ==
!missing: SHCilmToCindex
!missing: SHCindexToCilm
!missing: SHCilmToVector
subroutine pySHVectorToCilm(vector, cilm, lmax)
    use shtools, only: SHVectorToCilm
    implicit none
    double precision vector((lmax+1)*(lmax+1))
    double precision cilm(2,lmax+1,lmax+1)
    integer lmax
    !f2py intent(in)  lmax,vector
    !f2py intent(out) cilm
    call SHVectorToCilm(vector, cilm, lmax)
end subroutine
!missing: YilmIndex

!== Spherical harmonic conversions ==
!missing: SHrtoc
!missing: SHctor

!======== GLOBAL SPECTRAL ANALYSIS ========

!== Real spectral analysis ==
!missing: SHPowerL
!missing: SHPowerDensityL
!missing: SHCrossPowerL
!missing: SHCrossPowerDensityL
subroutine pySHPowerSpectrum(coeffs,lmax,spectrum)
    use shtools, only: SHPowerSpectrum
    implicit none
    double precision coeffs(2,lmax+1,lmax+1)
    double precision spectrum(lmax+1)
    integer lmax
    !f2py intent(out) spectrum
    !f2py intent(in)  coeffs
    !f2py integer intent(hide), depend(coeffs) :: lmax = shape(coeffs,1) - 1

    call SHPowerSpectrum(coeffs, lmax, spectrum)
end subroutine
!missing: SHPowerSpectrumDensity
subroutine pySHCrossPowerSpectrum(coeffs1, coeffs2, lmax, cspectrum)
    use shtools, only: SHCrossPowerSpectrum
    implicit none
    double precision coeffs1(2,lmax+1,lmax+1)
    double precision coeffs2(2,lmax+1,lmax+1)
    double precision cspectrum(lmax+1)
    integer lmax
    !f2py intent(in) coeffs1
    !f2py intent(in) coeffs2
    !f2py integer intent(hide),depend(coeffs1) :: lmax = shape(coeffs1,1) - 1
    !f2py intent(out) cspectrum

    call SHCrossPowerSpectrum(coeffs1,coeffs2,lmax,cspectrum)
end subroutine

subroutine pySHCrossPowerSpectrumDensity(coeffs1, coeffs2, lmax, cspectrum)
    use shtools, only: SHCrossPowerSpectrum
    implicit none
    double precision coeffs1(2,lmax+1,lmax+1)
    double precision coeffs2(2,lmax+1,lmax+1)
    double precision cspectrum(lmax+1)
    integer lmax
    !f2py intent(in) coeffs1
    !f2py intent(in) coeffs2
    !f2py integer intent(hide),depend(coeffs1) :: lmax = shape(coeffs1,1) - 1
    !f2py intent(out) cspectrum

    call SHCrossPowerSpectrumDensity(coeffs1,coeffs2,lmax,cspectrum)
end subroutine

subroutine pySHAdmitCorr(coeffs1,coeffs2,lmax,admit,corr, admit_error)
    use shtools, only: SHAdmitCorr
    implicit none
    double precision coeffs1(2,lmax+1,lmax+1)
    double precision coeffs2(2,lmax+1,lmax+1)
    double precision admit(lmax+1)
    double precision corr(lmax+1)
    double precision admit_error(lmax+1)
    integer lmax
    !f2py intent(in) coeffs1
    !f2py intent(in) coeffs2
    !f2py integer intent(hide),depend(coeffs1) :: lmax = shape(coeffs1,1) - 1
    !f2py intent(out) admit
    !f2py intent(out) corr
    !f2py intent(out) admit_error

    call SHAdmitCorr(coeffs1, coeffs2, lmax, admit, corr, admit_error=admit_error)
end subroutine

subroutine pySHConfidence(l, corr, return_value)
    use shtools, only: SHConfidence
    implicit none
    integer l
    double precision corr, return_value
    !f2py intent(in) l
    !f2py intent(in) corr
    !f2py intent(out) return_value

    return_value = SHConfidence(l,corr)
end subroutine

!== Complex spectral analysis ==
!missing: SHPowerLC
!missing: SHPowerDensityLC
!missing: SHCrossPowerLC
!missing: SHCrossPowerDensityLC
!missing: SHPowerSpectrumC
!missing: SHPowerSpectrumDensityC
!missing: SHCrossPowerSpectrumC
!missing: SHCrossPowerSpectrumDenistyC

!======== LOCALIZED SPECTRAL ANALYSIS ========
!== Multitaper spectral estimation (spherical cap domain) ==
subroutine pySHMultiTaperSE(mtse, mtse_dev, coeffs, lmax_coeffs, tapers, taper_order, &
                                lmaxt, ntapers, k, angles, clat, clon, taper_wt)
    use shtools, only: SHMultiTaperSE
    implicit none
    integer lmax_coeffs, lmaxt, k, ntapers
    double precision mtse(lmax_coeffs - lmaxt + 1)
    double precision mtse_dev(lmax_coeffs - lmaxt + 1)
    double precision coeffs(2,lmax_coeffs+1,lmax_coeffs+1)
    double precision tapers(lmaxt+1, ntapers), taper_wt(k)
    integer taper_order( ntapers )
    double precision clat, clon
    double precision angles(3)
    !f2py intent(out) mtse,mtse_dev
    !f2py integer intent(in),depend(tapers),optional :: lmaxt = shape(tapers,0)-1
    !f2py integer intent(in),depend(tapers),optional :: ntapers = shape(tapers,1)
    !f2py intent(in) coeffs
    !f2py integer intent(in),depend(coeffs),optional :: lmax_coeffs = shape(coeffs,1)-1
    !f2py intent(in) k
    !f2py intent(in) tapers
    !f2py intent(in) taper_order
    !f2py intent(in), optional :: clat=0
    !f2py intent(in), optional :: clon=0
    !f2py intent(in), optional :: angles = {0} 
    !f2py intent(in), optional :: taper_wt = {-1}
    if (taper_wt(1) > 0) then
        call SHMultiTaperSE(mtse, mtse_dev, coeffs, lmax_coeffs, tapers, taper_order, lmaxt, k, &
                                             angles, clat, clon, taper_wt=taper_wt)
    else
        call SHMultiTaperSE(mtse, mtse_dev, coeffs, lmax_coeffs, tapers, taper_order, lmaxt, k, &
                                                                  angles, clat, clon)
    endif
end subroutine

!missing: SHMultiTaperCSE
subroutine pySHLocalizedAdmitCorr(tapers, taper_order, lwin, clat, clon,&
    coeffs1, coeffs2, lmax,admit,corr, k, admit_error, corr_error)
    use shtools, only: SHLocalizedAdmitCorr
    implicit none
    integer lwin, lmax, k, i, j 
    integer taper_order(k)
    double precision tapers(lwin+1, k)
    double precision coeffs1(2,lmax+1,lmax+1)
    double precision coeffs2(2,lmax+1,lmax+1)
    double precision clat,clon
    double precision, dimension(lmax-lwin+1) :: admit,corr,admit_error,corr_error
    !f2py intent(in) tapers, taper_order
    !f2py intent(in) coeffs1, coeffs2
    !f2py intent(out) admit, corr, admit_error, corr_error
    !f2py integer intent(in),depend(tapers),optional :: k = shape(tapers,1)
    !f2py integer intent(in),depend(tapers),optional :: lwin = shape(tapers,0)-1
    !f2py integer intent(in),depend(coeffs1),optional :: lmax = shape(coeffs1,1)-1
    call SHLocalizedAdmitCorr(tapers, taper_order, lwin, clat, clon, coeffs1,&
                             coeffs2, lmax, admit, corr, k, admit_error, corr_error)
end subroutine

subroutine pySHReturnTapers(theta0, lmax, tapers, eigenvalues, taper_order)
    use shtools, only: SHReturnTapers
    implicit none
    integer lmax
    double precision theta0
    double precision tapers(lmax+1, (lmax+1)*(lmax+1))
    double precision eigenvalues((lmax+1)*(lmax+1))
    integer taper_order( (lmax+1)*(lmax+1) )
    !f2py intent(in) theta0, lmax
    !f2py intent(out) tapers, eigenvalues, taper_order
    call SHReturnTapers(theta0, lmax, tapers, eigenvalues, taper_order)
end subroutine

subroutine pySHReturnTapersM(theta0, lmax, m, tapers, eigenvalues)
    use shtools, only: SHReturnTapersM
    implicit none
    integer m, lmax
    double precision theta0
    double precision tapers(lmax+1, (lmax+1))
    double precision eigenvalues((lmax+1))
    !f2py intent(in) theta0, lmax, m
    !f2py intent(out) tapers, eigenvalues
    call SHReturnTapersM(theta0, lmax, m, tapers, eigenvalues)
end subroutine

subroutine pyComputeD0(d0, lmax, theta0)
    use shtools, only: ComputeD0
    implicit none
    integer lmax
    double precision theta0
    double precision d0(lmax+1, lmax+1)
    !f2py intent(in) lmax,theta0
    !f2py intent(out) d0
    call ComputeD0(d0, lmax, theta0)
end subroutine

subroutine pyComputeDm(dllm, lmax, m, theta0)
    use shtools, only: ComputeDm
    implicit none
    integer lmax,m
    double precision theta0
    double precision dllm(lmax+1, lmax+1)
    !f2py intent(in) lmax,theta0,m
    !f2py intent(out) dllm
    call ComputeDm(dllm, lmax, m, theta0)
end subroutine
!missing: ComputeDG82
!missing: SHFindLWin
subroutine pySHBiasK(tapers, lmaxt, k, inspec, lmax_in, outspec)
    use shtools, only: SHBiasK
    implicit none
    double precision tapers(lmaxt+1, k)
    double precision inspec(lmax_in+1)
    double precision outspec(lmax_in+1)
    integer lmaxt, k, lmax_in
    !f2py intent(in) tapers, inspec
    !f2py integer intent(in), depend(tapers), optional :: lmaxt   = shape(tapers,0)-1
    !f2py integer intent(in), depend(tapers), optional :: k   = shape(tapers,1)
    !f2py integer intent(in), depend(inspec), optional :: lmax_in = shape(inspec,0)-1
    !f2py intent(out) outspec
    call SHBiasK(tapers, lmaxt, k, inspec, lmax_in, outspec)
end subroutine

!missing: SHBiasAdmitCorr

subroutine pySHMTDebias(mtdebias, mtspectra, lmax, nl, tapers, lwin, k, lmid,n )
    use shtools, only: SHMTDebias
    implicit none
    double precision mtdebias(2,lmax+1 )
    double precision mtspectra(2, lmax+1)
    double precision tapers(lwin+1, k)
    double precision lmid(n)
    integer lmax, nl, lwin, k, n
    !!integer n should be: np.ceiling(lmax+1/nl)
    !f2py intent(in) mtspectra
    !f2py intent(in) nl
    !f2py intent(in) tapers
    !f2py intent(in) n
    !f2py integer intent(in), depend(mtspectra), optional :: lmax = shape(mtspectra,1)-1
    !f2py integer intent(in), depend(tapers), optional :: lwin = shape(tapers,0)-1
    !f2py integer intent(in), depend(tapers), optional :: k    = shape(tapers,1)
    !f2py intent(out) mtdebias
    !f2py intent(out) lmid
    call SHMTDebias(mtdebias, mtspectra, lmax, tapers, lwin, k, nl, lmid, n)
end subroutine

!missing: SHMTVarOpt0
!missing: SHMTVarOpt
!missing: SHSjkPG0
!missing: SHSjkPG

!== Localization windows (arbitrary domain) ==
subroutine pySHReturnTapersMap(tapers, eigenvalues, dh_mask, n, lmaxt, ntapers)
    use shtools, only: SHReturnTapersMap
    implicit none
    double precision tapers((lmaxt+1)*(lmaxt+1), ntapers)
    double precision eigenvalues(ntapers)
    integer lmaxt,n,ntapers
    integer dh_mask(N,2*N)
    !f2py intent(out) tapers, eigenvalues
    !f2py intent(in) dh_mask, n, lmaxt, ntapers
    !f2py integer intent(in), depend(dh_mask), optional :: n = shape(dh_mask,0)
    call SHReturnTapersMap(tapers, eigenvalues, dh_mask, n, 2, lmaxt, ntapers)
end subroutine
!missing: SHComputeDMap
!missing: SHCurve2Mask

!== Localization Bias (General) ==
!missing: SHBias

!== Other ==
!missing: SphericalCapCoef

!======== SPHERICAL HARMONIC ROTATIONS ========
subroutine pydjpi2(DJ, lmax)
    use shtools, only: djpi2
    implicit none
    integer lmax
    double precision DJ(lmax+1,lmax+1,lmax+1)
    !f2py intent(in) lmax
    !f2py intent(out) DJ
    call djpi2(DJ,lmax)
end subroutine

!missing: SHRotateCoef

subroutine pySHRotateRealCoef(coeffs_rot, coeffs, lmax, angles, DJ)
    use shtools, only: SHRotateRealCoef
    implicit none
    double precision coeffs_rot(2,lmax+1,lmax+1)
    double precision coeffs(2,lmax+1,lmax+1)
    double precision angles(3)
    double precision DJ(lmax+1,lmax+1,lmax+1)
    integer lmax
    !f2py intent(in) coeffs
    !f2py intent(in) angles
    !f2py intent(in) DJ
    !f2py integer intent(hide),depend(coeffs) :: lmax = shape(coeffs,1) - 1
    !f2py intent(out) coeffs_rot

    call SHRotateRealCoef(coeffs_rot, coeffs, lmax, angles, DJ)
end subroutine

subroutine pySHRotateRealCoef_direct(coeffs_rot, coeffs, lmax, angles)
    !----convenience function----
    use shtools, only: SHRotateRealCoef, djpi2
    implicit none
    double precision coeffs_rot(2,lmax+1,lmax+1)
    double precision coeffs(2,lmax+1,lmax+1)
    double precision angles(3)
    double precision DJ(lmax+1,lmax+1,lmax+1)
    integer lmax
    !f2py intent(in) coeffs
    !f2py intent(in) angles
    !f2py intent(in) DJ
    !f2py integer intent(hide),depend(coeffs) :: lmax = shape(coeffs,1) - 1
    !f2py intent(out) coeffs_rot

    call djpi2(DJ,lmax)
    call SHRotateRealCoef(coeffs_rot, coeffs, lmax, angles, DJ)
end subroutine

!======== GRAVITY AND MAGNETICS ========
!== Gravity routines ==
!== Magnetics routines ==

!======== OTHER ROUTINES ========
