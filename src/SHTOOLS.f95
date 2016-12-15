module SHTOOLS
!------------------------------------------------------------------------------
!
!   This module contains an interface block defining all the routines
!   used in the archive SHTOOLS. These are necessary in order to use
!   implicitly shaped arrays with most subroutines.
!
!   Copyright (c) 2016, SHTOOLS
!   All rights reserved.
!
!------------------------------------------------------------------------------
    integer, parameter ::   CSPHASE_DEFAULT = 1
                            ! The default for all routines is to EXCLUDE
                            ! the Condon-Shortley phase of (-1)^m
                            ! in front of the Legendre functions.

    interface

        subroutine PlmBar(p, lmax, z, csphase, cnorm, exitstatus)
            integer, intent(in) ::  lmax
            real*8, intent(out) ::  p(:)
            real*8, intent(in) ::   z
            integer, intent(in), optional :: csphase, cnorm
            integer, intent(out), optional :: exitstatus
        end subroutine PlmBar

        subroutine PlmBar_d1(p, dp, lmax, z, csphase, cnorm, exitstatus)
            integer, intent(in) ::  lmax
            real*8, intent(out) ::  p(:), dp(:)
            real*8, intent(in) ::   z
            integer, intent(in), optional :: csphase, cnorm
            integer, intent(out), optional :: exitstatus
        end subroutine PlmBar_d1

        subroutine PlBar(p, lmax, z, exitstatus)
            integer, intent(in) ::  lmax
            real*8, intent(out) ::  p(:)
            real*8, intent(in) ::   z
            integer, intent(out), optional :: exitstatus
        end subroutine PlBar

        subroutine PlBar_d1(p, dp, lmax, z, exitstatus)
            integer, intent(in) ::  lmax
            real*8, intent(out) ::  p(:), dp(:)
            real*8, intent(in) ::   z
            integer, intent(out), optional :: exitstatus
        end subroutine PlBar_d1

        subroutine PlmON(p, lmax, z, csphase, cnorm, exitstatus)
            integer, intent(in) ::  lmax
            real*8, intent(out) ::  p(:)
            real*8, intent(in) ::   z
            integer, intent(in), optional :: csphase, cnorm
            integer, intent(out), optional :: exitstatus
        end subroutine PlmON

        subroutine PlmON_d1(p, dp, lmax, z, csphase, cnorm, exitstatus)
            integer, intent(in) ::  lmax
            real*8, intent(out) ::  p(:), dp(:)
            real*8, intent(in) ::   z
            integer, intent(in), optional :: csphase, cnorm
            integer, intent(out), optional :: exitstatus
        end subroutine PlmON_d1

        subroutine PlON(p, lmax, z, exitstatus)
            integer, intent(in) ::  lmax
            real*8, intent(out) ::  p(:)
            real*8, intent(in) ::   z
            integer, intent(out), optional :: exitstatus
        end subroutine PlON

        subroutine PlON_d1(p, dp, lmax, z, exitstatus)
            integer, intent(in) ::  lmax
            real*8, intent(out) ::  p(:), dp(:)
            real*8, intent(in) ::   z
            integer, intent(out), optional :: exitstatus
            end subroutine PlON_d1

        subroutine PlmSchmidt(p,lmax,z, csphase, cnorm, exitstatus)
            integer, intent(in) ::  lmax
            real*8, intent(out) ::  p(:)
            real*8, intent(in) ::   z
            integer, intent(in), optional :: csphase, cnorm
            integer, intent(out), optional :: exitstatus
        end subroutine PlmSchmidt

        subroutine PlmSchmidt_d1(p, dp, lmax, z, csphase, cnorm, exitstatus)
            integer, intent(in) ::  lmax
            real*8, intent(out) ::  p(:), dp(:)
            real*8, intent(in) ::   z
            integer, intent(in), optional :: csphase, cnorm
            integer, intent(out), optional :: exitstatus
        end subroutine PlmSchmidt_d1

        subroutine PlSchmidt(p,lmax,z, exitstatus)
            integer, intent(in) ::  lmax
            real*8, intent(out) ::  p(:)
            real*8, intent(in) ::   z
            integer, intent(out), optional :: exitstatus
        end subroutine PlSchmidt

        subroutine PlSchmidt_d1(p, dp, lmax, z, exitstatus)
            integer, intent(in) ::  lmax
            real*8, intent(out) ::  p(:), dp(:)
            real*8, intent(in) ::   z
            integer, intent(out), optional :: exitstatus
        end subroutine PlSchmidt_d1

        subroutine PLegendreA(p,lmax,z, csphase, exitstatus)
            integer, intent(in) ::  lmax
            real*8, intent(out) ::  p(:)
            real*8, intent(in) ::   z
            integer, intent(in), optional :: csphase
            integer, intent(out), optional :: exitstatus
        end subroutine PLegendreA

        subroutine PLegendreA_d1(p, dp, lmax, z, csphase, exitstatus)
            integer, intent(in) ::  lmax
            real*8, intent(out) ::  p(:), dp(:)
            real*8, intent(in) ::   z
            integer, intent(in), optional :: csphase
            integer, intent(out), optional :: exitstatus
        end subroutine PLegendreA_d1

        subroutine PLegendre(p,lmax,z, exitstatus)
            integer, intent(in) ::  lmax
            real*8, intent(out) ::  p(:)
            real*8, intent(in) ::   z
            integer, intent(out), optional :: exitstatus
        end subroutine PLegendre

        subroutine PLegendre_d1(p, dp, lmax, z, exitstatus)
            integer, intent(in) ::  lmax
            real*8, intent(out) ::  p(:), dp(:)
            real*8, intent(in) ::   z
            integer, intent(out), optional :: exitstatus
        end subroutine PLegendre_d1

        integer function PlmIndex(l,m)
            integer, intent(in) :: l, m
        end function PlmIndex

        subroutine SHExpandDH(grid, n, cilm, lmax, norm, sampling, &
                              csphase, lmax_calc, exitstatus)
            real*8, intent(in) ::   grid(:,:)
            real*8, intent(out) ::  cilm(:,:,:)
            integer, intent(in) ::  n
            integer, intent(out) :: lmax
            integer, intent(in), optional :: norm, sampling, csphase, lmax_calc
            integer, intent(out), optional :: exitstatus
        end subroutine SHExpandDH

        subroutine MakeGridDH(griddh, n, cilm, lmax, norm, sampling, &
                              csphase, lmax_calc, exitstatus)
            real*8, intent(in) ::   cilm(:,:,:)
            real*8, intent(out) ::  griddh(:,:)
            integer, intent(in) ::  lmax
            integer, intent(out) :: n
            integer, intent(in), optional :: norm, sampling, csphase, lmax_calc
            integer, intent(out), optional :: exitstatus
        end subroutine MakeGridDH

        subroutine SHExpandDHC(grid, n, cilm, lmax, norm, sampling, &
                               csphase, lmax_calc, exitstatus)
            complex*16, intent(in) ::   grid(:,:)
            complex*16, intent(out) ::  cilm(:,:,:)
            integer, intent(in) ::  n
            integer, intent(out) :: lmax
            integer, intent(in), optional :: norm, sampling, csphase, lmax_calc
            integer, intent(out), optional :: exitstatus
        end subroutine SHExpandDHC

        subroutine MakeGridDHC(griddh, n, cilm, lmax, norm, sampling, &
                               csphase, lmax_calc, exitstatus)
            complex*16, intent(in) ::   cilm(:,:,:)
            complex*16, intent(out) ::  griddh(:,:)
            integer, intent(in) ::  lmax
            integer, intent(out) :: n
            integer, intent(in), optional :: norm, sampling, csphase, lmax_calc
            integer, intent(out), optional :: exitstatus
        end subroutine MakeGridDHC

        subroutine SHGLQ(lmax, zero, w, plx, norm, csphase, cnorm, exitstatus)
            integer, intent(in) ::  lmax
            real*8, intent(out) ::  zero(:), w(:)
            real*8, intent(out), optional :: plx(:,:)
            integer, intent(in), optional :: norm, csphase, cnorm
            integer, intent(out), optional :: exitstatus
        end subroutine SHGLQ

        subroutine SHExpandGLQ(cilm, lmax, gridglq, w, plx, zero, norm, &
                               csphase, lmax_calc, exitstatus)
            real*8, intent(in) ::   w(:), gridglq(:,:)
            real*8, intent(in), optional :: plx(:,:), zero(:)
            real*8, intent(out) ::  cilm(:,:,:)
            integer, intent(in) ::  lmax
            integer, intent(in), optional :: norm, csphase, lmax_calc
            integer, intent(out), optional :: exitstatus
        end subroutine SHExpandGLQ

        subroutine MakeGridGLQ(gridglq, cilm, lmax, plx, zero, norm, &
                               csphase, lmax_calc, exitstatus)
            real*8, intent(in) ::   cilm(:,:,:)
            real*8, intent(in), optional :: plx(:,:), zero(:)
            real*8, intent(out) ::  gridglq(:,:)
            integer, intent(in) ::  lmax
            integer, intent(in), optional :: norm, csphase, lmax_calc
            integer, intent(out), optional :: exitstatus
        end subroutine MakeGridGLQ

        subroutine SHExpandGLQC(cilm, lmax, gridglq, w, plx, zero, norm, &
                                csphase, lmax_calc, exitstatus)
            real*8, intent(in) ::   w(:)
            complex*16, intent(in) ::   gridglq(:,:)
            real*8, intent(in), optional :: plx(:,:), zero(:)
            complex*16, intent(out) ::  cilm(:,:,:)
            integer, intent(in) ::  lmax
            integer, intent(in), optional :: norm, csphase, lmax_calc
            integer, intent(out), optional :: exitstatus
        end subroutine SHExpandGLQC

        subroutine MakeGridGLQC(gridglq, cilm, lmax, plx, zero, norm, &
                                csphase, lmax_calc, exitstatus)
            complex*16, intent(in) ::   cilm(:,:,:)
            real*8, intent(in), optional :: plx(:,:), zero(:)
            complex*16, intent(out) ::  gridglq(:,:)
            integer, intent(in) ::  lmax
            integer, intent(in), optional :: norm, csphase, lmax_calc
            integer, intent(out), optional :: exitstatus
        end subroutine MakeGridGLQC

        subroutine GLQGridCoord(latglq, longlq, lmax, nlat, nlong, exitstatus)
            integer, intent(in) ::  lmax
            integer, intent(out) :: nlat, nlong
            real*8, intent(out) ::  latglq(:), longlq(:)
            integer, intent(out), optional :: exitstatus
        end subroutine GLQGridCoord

        subroutine SHExpandLSQ(cilm, d, lat, lon, nmax, lmax, norm, &
                               chi2, csphase, exitstatus)
            real*8, intent(in) ::   d(:), lat(:), lon(:)
            real*8, intent(out) ::  cilm(:,:,:)
            integer, intent(in) ::  nmax, lmax
            integer, intent(in), optional ::  norm, csphase
            real*8, intent(out), optional ::  chi2
            integer, intent(out), optional :: exitstatus
        end subroutine SHExpandLSQ

        subroutine MakeGrid2d(grid, cilm, lmax, interval, nlat, nlong, &
                              norm, csphase, f, a, north, south, east, west, &
                              dealloc, exitstatus)
            real*8, intent(in) ::   cilm(:,:,:), interval
            real*8, intent(out) ::  grid(:,:)
            integer, intent(in) ::  lmax
            integer, intent(out) :: nlat, nlong
            integer, intent(in), optional ::    norm, csphase, dealloc
            real*8, intent(in), optional ::     f, a, north, south, east, west
            integer, intent(out), optional :: exitstatus
        end subroutine MakeGrid2D

        real*8 function MakeGridPoint(cilm, lmax, lat, lon, norm, &
                                      csphase, dealloc)
            real*8, intent(in) ::   cilm(:,:,:), lat, lon
            integer, intent(in) ::  lmax
            integer, intent(in), optional :: norm, csphase, dealloc
        end function MakeGridPoint

        complex*16 function MakeGridPointC(cilm, lmax, lat, lon, norm, &
                                           csphase, dealloc)
            complex*16, intent(in) :: cilm(:,:,:)
            real*8, intent(in) :: lat, lon
            integer, intent(in) :: lmax
            integer, intent(in), optional :: norm, csphase, dealloc
        end function MakeGridPointC

        subroutine SHMultiply(shout, sh1, lmax1, sh2, lmax2, precomp, &
                              norm, csphase, exitstatus)
            real*8, intent(out) ::  shout(:,:,:)
            real*8, intent(in) ::   sh1(:,:,:), sh2(:,:,:)
            integer, intent(in) ::  lmax1, lmax2
            integer, intent(in), optional ::  precomp, norm, csphase
            integer, intent(out), optional :: exitstatus
        end subroutine SHMultiply

        subroutine SHRead(filename, cilm, lmax, skip, header, error, &
                          exitstatus)
            character(*), intent(in) :: filename
            integer, intent(out) :: lmax
            real*8, intent(out) ::  cilm(:,:,:)
            real*8, intent(out), optional :: header(:), error(:,:,:)
            integer, intent(in), optional :: skip
            integer, intent(out), optional :: exitstatus
        end subroutine SHRead

        subroutine SHRead2(filename, cilm, lmax, gm, r0_pot, error, dot, &
                           doystart, doyend, epoch, exitstatus)
            character(*), intent(in) :: filename
            integer, intent(out) :: lmax
            real*8, intent(out) ::  cilm(:,:,:), gm, r0_pot
            real*8, intent(out), optional :: error(:,:,:), dot(:,:,:), &
                                             doystart, doyend, epoch
            integer, intent(out), optional :: exitstatus
        end subroutine SHRead2

        subroutine SHReadJPL(filename, cilm, lmax, error, gm, formatstring, &
                             exitstatus)
            character(*), intent(in) :: filename
            integer, intent(in) ::  lmax
            real*8, intent(out) ::  cilm(:,:,:)
            real*8, intent(out), optional :: error(:,:,:), gm(2)
            character, intent(in), optional :: formatstring*6
            integer, intent(out), optional :: exitstatus
        end subroutine SHReadJPL

        subroutine SHCilmToVector(cilm, vector, lmax, exitstatus)
            real*8, intent(in) ::   cilm(:,:,:)
            real*8, intent(out) ::  vector(:)
            integer, intent(in) ::  lmax
            integer, intent(out), optional :: exitstatus
        end subroutine SHCilmToVector

        subroutine SHVectorToCilm(vector, cilm, lmax, exitstatus)
            real*8, intent(out) ::  cilm(:,:,:)
            real*8, intent(in) ::   vector(:)
            integer, intent(in) ::  lmax
            integer, intent(out), optional :: exitstatus
        end subroutine SHVectorToCilm

        subroutine SHCilmToCindex(cilm, cindex, degmax, exitstatus)
            real*8, intent(in) ::   cilm(:,:,:)
            real*8, intent(out) ::  cindex(:,:)
            integer, intent(in), optional :: degmax
            integer, intent(out), optional :: exitstatus
        end subroutine SHCilmToCindex

        subroutine SHCindexToCilm(cindex, cilm, degmax, exitstatus)
            real*8, intent(out) ::  cilm(:,:,:)
            real*8, intent(in) ::   cindex(:,:)
            integer, intent(in), optional :: degmax
            integer, intent(out), optional :: exitstatus
        end subroutine SHCindexToCilm

        subroutine SHrtoc(rcilm, ccilm, degmax, convention, switchcs, &
                          exitstatus)
            real*8, intent(in) ::   rcilm(:,:,:)
            real*8, intent(out) ::  ccilm(:,:,:) 
            integer, intent(in), optional :: degmax, convention, switchcs
            integer, intent(out), optional :: exitstatus
        end subroutine SHrtoc

        subroutine SHctor(ccilm, rcilm, degmax, convention, switchcs, &
                          exitstatus)
            real*8, intent(in) ::   ccilm(:,:,:)
            real*8, intent(out) ::  rcilm(:,:,:)
            integer, intent(in), optional :: degmax, convention, switchcs
            integer, intent(out), optional :: exitstatus
        end subroutine SHctor

        subroutine djpi2(dj, lmax, exitstatus)
            integer, intent(in) ::  lmax
            real*8, intent(out) ::  dj(:,:,:)
            integer, intent(out), optional :: exitstatus
        end subroutine djpi2

        subroutine SHRotateCoef(x, cof, rcof, dj, lmax, exitstatus)
            real*8, intent(in) ::   cof(:,:), dj(:,:,:), x(3)
            real*8, intent(out) ::  rcof(:,:)
            integer, intent(in) ::  lmax
            integer, intent(out), optional :: exitstatus
        end subroutine SHRotateCoef

        subroutine SHRotateRealCoef(cilmrot, cilm, lmax, x, dj, exitstatus)
            real*8, intent(in) ::   cilm(:,:,:), x(:), dj(:,:,:)
            real*8, intent(out) ::  cilmrot(:,:,:)
            integer, intent(in) ::  lmax
            integer, intent(out), optional :: exitstatus
        end subroutine SHRotateRealCoef

        real*8 function SHPowerL(c, l)
            real*8, intent(in) :: c(:,:,:)
            integer, intent(in) :: l
        end function SHPowerL

        real*8 function SHPowerDensityL(c, l)
            real*8, intent(in) ::   c(:,:,:)
            integer, intent(in) ::  l
        end function SHPowerDensityL

        real*8 function SHCrossPowerL(c1, c2, l)
            real*8, intent(in) ::   c1(:,:,:), c2(:,:,:)
            integer, intent(in) ::  l
        end function SHCrossPowerL

        real*8 function SHCrossPowerDensityL(c1, c2, l)
            real*8, intent(in) ::   c1(:,:,:), c2(:,:,:)
            integer, intent(in) ::  l
        end function SHCrossPowerDensityL

        subroutine SHPowerSpectrum(c, lmax, spectra, exitstatus)
            real*8, intent(in) ::   c(:,:,:)
            integer, intent(in) ::  lmax
            real*8, intent(out) ::  spectra(:)
            integer, intent(out), optional :: exitstatus
        end subroutine SHPowerSpectrum

        subroutine SHPowerSpectrumDensity(c, lmax, spectra, exitstatus)
            real*8, intent(in) ::   c(:,:,:)
            integer, intent(in) ::  lmax
            real*8, intent(out) ::  spectra(:)
            integer, intent(out), optional :: exitstatus
        end subroutine SHPowerSpectrumDensity

        subroutine SHCrossPowerSpectrum(c1, c2, lmax, cspectra, exitstatus)
            real*8, intent(in) ::   c1(:,:,:), c2(:,:,:)
            integer, intent(in) ::  lmax
            real*8, intent(out) ::  cspectra(:)
            integer, intent(out), optional :: exitstatus
        end subroutine SHCrossPowerSpectrum

        subroutine SHCrossPowerSpectrumDensity(c1, c2, lmax, cspectra, &
                                               exitstatus)
            real*8, intent(in) ::   c1(:,:,:), c2(:,:,:)
            integer, intent(in) ::  lmax
            real*8, intent(out) ::  cspectra(:)
            integer, intent(out), optional :: exitstatus
        end subroutine SHCrossPowerSpectrumDensity

        subroutine SHAdmitCorr(G, T, lmax, admit, corr, admit_error, exitstatus)
            real*8, intent(in) ::   G(:,:,:), T(:,:,:)
            integer, intent(in) ::  lmax
            real*8, intent(out) ::  admit(:), corr(:)
            real*8, intent(out), optional ::  admit_error(:)
            integer, intent(out), optional :: exitstatus
        end subroutine SHAdmitCorr

        real*8 function SHConfidence(l_conf, r)
            real*8, intent(in) :: r
            integer, intent(in) :: l_conf
        end function SHConfidence

        real*8 function SHPowerLC(c, l)
            complex*16, intent(in) :: c(:,:,:)
            integer, intent(in) :: l
        end function SHPowerLC

        real*8 function SHPowerDensityLC(c, l)
            complex*16, intent(in) ::   c(:,:,:)
            integer, intent(in) ::  l
        end function SHPowerDensityLC

        complex*16 function SHCrossPowerLC(c1, c2, l)
            complex*16, intent(in) :: c1(:,:,:), c2(:,:,:)
            integer, intent(in) ::  l
        end function SHCrossPowerLC

        Complex*16 function SHCrossPowerDensityLC(c1, c2, l)
            complex*16, intent(in) :: c1(:,:,:), c2(:,:,:)
            integer, intent(in) ::  l
        end function SHCrossPowerDensityLC

        subroutine SHPowerSpectrumC(c, lmax, spectra, exitstatus)
            complex*16, intent(in) ::   c(:,:,:)
            integer, intent(in) ::  lmax
            real*8, intent(out) ::  spectra(:)
            integer, intent(out), optional :: exitstatus
        end subroutine SHPowerSpectrumC

        subroutine SHPowerSpectrumDensityC(c, lmax, spectra, exitstatus)
            complex*16, intent(in) ::   c(:,:,:)
            integer, intent(in) ::  lmax
            real*8, intent(out) ::  spectra(:)
            integer, intent(out), optional :: exitstatus
        end subroutine SHPowerSpectrumDensityC

        subroutine SHCrossPowerSpectrumC(c1, c2, lmax, cspectra, exitstatus)
            complex*16, intent(in) ::   c1(:,:,:), c2(:,:,:)
            integer, intent(in) ::  lmax
            complex*16, intent(out) ::  cspectra(:)
            integer, intent(out), optional :: exitstatus
        end subroutine SHCrossPowerSpectrumC

        subroutine SHCrossPowerSpectrumDensityC(c1, c2, lmax, cspectra,&
                                                exitstatus)
            complex*16, intent(in) ::   c1(:,:,:), c2(:,:,:)
            integer, intent(in) ::  lmax
            complex*16, intent(out) ::  cspectra(:)
            integer, intent(out), optional :: exitstatus
        end subroutine SHCrossPowerSpectrumDensityC

        subroutine SHMultiTaperSE(mtse, sd, sh, lmax, tapers, taper_order, &
                                  lmaxt, k, alpha, lat, lon, taper_wt, norm, &
                                  csphase, exitstatus)
            real*8, intent(out) ::  mtse(:), sd(:)
            real*8, intent(in) ::   sh(:,:,:), tapers(:,:)
            integer, intent(in) ::  lmax, lmaxt, k, taper_order(:)
            real*8, intent(in), optional :: alpha(:), lat, lon, taper_wt(:)
            integer, intent(in), optional :: csphase, norm
            integer, intent(out), optional :: exitstatus
        end subroutine SHMultiTaperSE

        subroutine SHMultiTaperCSE(mtse, sd, sh1, lmax1, sh2, lmax2, tapers, &
                                   taper_order, lmaxt, k, alpha, lat, lon, &
                                   taper_wt, norm, csphase, exitstatus)
            real*8, intent(out) ::  mtse(:), sd(:)
            real*8, intent(in) ::   sh1(:,:,:), sh2(:,:,:), tapers(:,:)
            integer, intent(in) ::  lmax1, lmax2, lmaxt, k, taper_order(:)
            real*8, intent(in), optional :: alpha(:), lat, lon, taper_wt(:)
            integer, intent(in), optional :: csphase, norm
            integer, intent(out), optional :: exitstatus
        end subroutine SHMultiTaperCSE

        subroutine SHLocalizedAdmitCorr(tapers, taper_order, lwin, lat, lon, &
                                        g, t, lmax, admit, corr, k, &
                                        admit_error, corr_error, taper_wt, &
                                        mtdef, k1linsig, exitstatus)
            real*8, intent(in) ::   tapers(:,:), lat, lon, g(:,:,:), t(:,:,:)
            integer, intent(in) ::  lwin, lmax, k, taper_order(:)
            real*8, intent(out) ::  admit(:), corr(:)
            real*8, intent(out), optional ::    admit_error(:), corr_error(:)
            integer, intent(in), optional ::    mtdef, k1linsig
            real*8, intent(in), optional :: taper_wt(:)
            integer, intent(out), optional :: exitstatus
        end subroutine SHLocalizedAdmitCorr

        subroutine SHReturnTapers(theta0, lmax, tapers, eigenvalues, &
                                  taper_order, exitstatus)
            real*8, intent(in) ::   theta0
            integer, intent(in) ::  lmax
            real*8, intent(out) ::  tapers(:,:), eigenvalues(:)
            integer, intent(out) :: taper_order(:)
            integer, intent(out), optional :: exitstatus
        end subroutine SHReturnTapers

        subroutine SHReturnTapersM(theta0, lmax, m, tapers, &
                                   eigenvalues, shannon, exitstatus)
            real*8, intent(in) ::   theta0
            integer, intent(in) ::  lmax, m
            real*8, intent(out) ::  tapers(:,:), eigenvalues(:)
            real*8, intent(out), optional :: shannon
            integer, intent(out), optional :: exitstatus
        end subroutine SHReturnTapersM

        subroutine ComputeDm(dllm, lmax, m, theta0, exitstatus)
            real*8, intent(out) ::  dllm(:,:)
            real*8, intent(in) ::   theta0
            integer, intent(in) ::  lmax, m
            integer, intent(out), optional :: exitstatus
        end subroutine ComputeDm

        subroutine ComputeDG82(dG82, lmax, m, theta0, exitstatus)
            real*8, intent(out) ::  dG82(:,:)
            real*8, intent(in) ::   theta0
            integer, intent(in) ::  lmax, m
            integer, intent(out), optional :: exitstatus
        end subroutine ComputeDG82

        integer function SHFindLWin(theta0, m, alpha, taper_number)
            real*8, intent(in) ::   theta0, alpha
            integer, intent(in) ::  m
            integer, intent(in), optional :: taper_number
        end function SHFindLWin

        subroutine SHBiasK(tapers, lwin, k, incspectra, ldata, &
                           outcspectra, taper_wt, save_cg, exitstatus)
            real*8, intent(in) ::   tapers(:,:), incspectra(:)
            real*8, intent(out) ::  outcspectra(:)
            integer, intent(in) ::  lwin, ldata, k
            real*8, intent(in), optional :: taper_wt(:)
            integer, intent(in), optional :: save_cg
            integer, intent(out), optional :: exitstatus
        end subroutine SHBiasK

        subroutine SHMTCouplingMatrix(Mmt, lmax, tapers_power, lwin, k, &
                                      taper_wt, exitstatus)
            real*8, intent(out) :: Mmt(:,:)
            real*8, intent(in)  :: tapers_power(:,:)
            real*8, intent(in), optional :: taper_wt(:)
            integer, intent(in) :: lmax, k, lwin
            integer, intent(out), optional :: exitstatus
        end subroutine

        subroutine SHBiasAdmitCorr(sgt, sgg, stt, lmax, tapers, lwin, k, &
                                   admit, corr, mtdef, taper_wt, exitstatus)
            real*8, intent(in) ::   sgt(:), sgg(:), stt(:), tapers(:,:)
            integer, intent(in) ::  lmax, lwin, k
            real*8, intent(out) ::  admit(:), corr(:)
            integer, intent(in), optional ::    mtdef
            real*8, intent(in), optional :: taper_wt(:)
            integer, intent(out), optional :: exitstatus
        end subroutine SHBiasAdmitCorr

        subroutine SHMTDebias (mtdebias, mtspectra, lmax, tapers, lwin, k, nl, &
                               lmid, n, taper_wt, exitstatus)
            real*8, intent(out) ::  mtdebias(:,:), lmid(:)
            real*8, intent(in) :: mtspectra(:,:), tapers(:,:)
            real*8, intent(in), optional :: taper_wt(:)
            integer, intent(in) :: lmax, k, lwin, nl
            integer, intent(out) :: n
            integer, intent(out), optional :: exitstatus
        end subroutine SHMTDebias

        subroutine SHMTVarOpt(l, tapers, taper_order, lwin, kmax, Sff, &
                              var_opt, var_unit, weight_opt, &
                              unweighted_covar, nocross, exitstatus)
            real*8, intent(in) ::   tapers(:,:), Sff(:)
            real*8, intent(out) ::  var_opt(:), var_unit(:)
            integer, intent(in) ::  l, lwin, kmax, taper_order(:)
            real*8, intent(out), optional :: weight_opt(:,:), unweighted_covar(:,:)
            integer, intent(in), optional :: nocross
            integer, intent(out), optional :: exitstatus
        end subroutine SHMTVarOpt

        complex*16 function SHSjkPG(incspectra, l, m, mprime, hj_real, &
                                    hk_real, mj, mk, lwin, hkcc)
            real*8, intent(in) ::   incspectra(:), hj_real(:), hk_real(:)
            integer, intent(in) ::  lwin, l, m, mprime, mj, mk, hkcc
        end function SHSjkPG

        subroutine SHReturnTapersMap(tapers, eigenvalues, dh_mask, n_dh, &
                                        lmax, sampling, ntapers, exitstatus)
            real*8, intent(out) ::  tapers(:,:), eigenvalues(:)
            integer, intent(in) ::  dh_mask(:,:), n_dh, lmax
            integer, intent(in), optional ::    sampling, ntapers
            integer, intent(out), optional :: exitstatus
        end subroutine SHReturnTapersMap

        subroutine SHBiasKMask(tapers, lwin, k, incspectra, ldata, &
                           outcspectra, taper_wt, save_cg, exitstatus)
            real*8, intent(in) ::   tapers(:,:), incspectra(:)
            real*8, intent(out) ::  outcspectra(:)
            integer, intent(in) ::  lwin, ldata, k
            real*8, intent(in), optional :: taper_wt(:)
            integer, intent(in), optional :: save_cg
            integer, intent(out), optional :: exitstatus
        end subroutine SHBiasKMask

        subroutine SHMultiTaperMaskSE(mtse, sd, sh, lmax, tapers, &
                                      lmaxt, k, taper_wt, norm, csphase, &
                                      exitstatus)
            real*8, intent(out) ::  mtse(:), sd(:)
            real*8, intent(in) ::   sh(:,:,:), tapers(:,:)
            integer, intent(in) ::  lmax, lmaxt, k
            real*8, intent(in), optional :: taper_wt(:)
            integer, intent(in), optional :: csphase, norm
            integer, intent(out), optional :: exitstatus
        end subroutine SHMultiTaperMaskSE

        subroutine SHMultiTaperMaskCSE(mtse, sd, sh1, lmax1, sh2, lmax2, &
                                       tapers, lmaxt, k, taper_wt, norm, &
                                       csphase, exitstatus)
            real*8, intent(out) ::  mtse(:), sd(:)
            real*8, intent(in) ::   sh1(:,:,:), sh2(:,:,:), tapers(:,:)
            integer, intent(in) ::  lmax1, lmax2, lmaxt, k
            real*8, intent(in), optional :: taper_wt(:)
            integer, intent(in), optional ::    csphase, norm
            integer, intent(out), optional :: exitstatus
        end subroutine SHMultiTaperMaskCSE

        subroutine ComputeDMap(Dij, dh_mask, n_dh, lmax, sampling, exitstatus)
            real*8, intent(out) ::  Dij(:,:)
            integer, intent(in) ::  dh_mask(:,:), n_dh, lmax
            integer, intent(in), optional :: sampling
            integer, intent(out), optional :: exitstatus
        end subroutine ComputeDMap

        subroutine Curve2Mask(dhgrid, n, sampling, profile, nprofile, NP, &
                              centralmeridian, exitstatus)
            integer, intent(out) :: dhgrid(:,:)
            real*8, intent(in) ::   profile(:,:)
            integer, intent(in) ::  n, sampling, nprofile, np
            integer, intent(in), optional :: centralmeridian
            integer, intent(out), optional :: exitstatus
        end subroutine Curve2Mask

        subroutine SHBias(Shh, lwin, incspectra, ldata, outcspectra, save_cg, &
                          exitstatus)
            real*8, intent(in) ::   Shh(:), incspectra(:)
            real*8, intent(out) ::  outcspectra(:)
            integer, intent(in) ::  lwin, ldata
            integer, intent(in), optional :: save_cg
            integer, intent(out), optional :: exitstatus
        end subroutine SHBias

        subroutine SphericalCapCoef(coef, theta, lmax, exitstatus)
            real*8, intent(out) ::  coef(:)
            real*8, intent(in) ::   theta
            integer, intent(in), optional ::  lmax
            integer, intent(out), optional :: exitstatus
        end subroutine SphericalCapCoef

        subroutine MakeGravGridDH(cilm, lmax, gm, r0, a, f, rad, theta, phi, &
                                  total, n, sampling, lmax_calc, omega, &
                                  normal_gravity, pot, exitstatus)
            real*8, intent(in) ::   cilm(:,:,:), gm, r0, a, f
            real*8, intent(out) ::  rad(:,:), theta(:,:), phi(:,:), total(:,:)
            real*8, intent(in), optional :: omega
            real*8, intent(out), optional :: pot(:,:)
            integer, intent(in) ::  lmax
            integer, intent(out) :: n
            integer, intent(in), optional :: sampling, lmax_calc, normal_gravity
            integer, intent(out), optional :: exitstatus
        end subroutine MakeGravGridDH

        subroutine MakeGravGradGridDH(cilm, lmax, gm, r0, a, f, vxx, vyy, &
                                      vzz, vxy, vxz, vyz, n, sampling, &
                                      lmax_calc, exitstatus)
            real*8, intent(in) ::   cilm(:,:,:), gm, r0, a, f
            real*8, intent(out) ::  vxx(:,:), vyy(:,:), vzz(:,:), vxy(:,:), &
                                    vxz(:,:), vyz(:,:)
            integer, intent(in) ::  lmax
            integer, intent(out) :: n
            integer, intent(in), optional :: sampling, lmax_calc
            integer, intent(out), optional :: exitstatus
        end subroutine MakeGravGradGridDH

        subroutine MakeGeoidGrid(geoid, cilm, lmax, r0pot, GM, PotRef, omega, &
                                 r, gridtype, order, nlat, nlong, interval, &
                                 lmax_calc, a, f, exitstatus)
            real*8, intent(out) ::  geoid(:,:)
            real*8, intent(in) ::   cilm(:,:,:), r0pot, GM, r, PotRef, omega
            integer, intent(in) :: lmax, order, gridtype
            integer, intent(in), optional :: lmax_calc
            integer, intent(out) :: nlat, nlong
            real*8, intent(in), optional :: interval, a, f
            integer, intent(out), optional :: exitstatus
        end subroutine MakeGeoidGrid

        subroutine CilmPlus(cilm, gridin, lmax, nmax, mass, d, rho, gridtype, &
                            w, zero, plx, n, dref, exitstatus)
            real*8, intent(in) ::   gridin(:,:), mass, rho
            real*8, intent(in), optional :: w(:), zero(:), plx(:,:), dref
            real*8, intent(out) ::  cilm(:,:,:), d
            integer, intent(in) ::  lmax, nmax, gridtype
            integer, intent(in), optional :: n
            integer, intent(out), optional :: exitstatus
        end subroutine CilmPlus

        subroutine CilmMinus(cilm, gridin, lmax, nmax, mass, d, rho, &
                             gridtype, w, zero, plx, n, dref, exitstatus)
            real*8, intent(in) ::   gridin(:,:), mass, rho
            real*8, intent(in), optional :: w(:), zero(:), plx(:,:), dref
            real*8, intent(out) ::  cilm(:,:,:), d
            integer, intent(in) ::  lmax, nmax, gridtype
            integer, intent(in), optional :: n
            integer, intent(out), optional :: exitstatus
        end subroutine CilmMinus

        subroutine CilmPlusRhoH(cilm, gridin, lmax, nmax, mass, d, rho, &
                                gridtype, w, zero, plx, n, dref, exitstatus)
            real*8, intent(in) ::   gridin(:,:), mass, rho(:,:)
            real*8, intent(in), optional :: w(:), zero(:), plx(:,:), dref
            real*8, intent(out) ::  cilm(:,:,:), d
            integer, intent(in) ::  lmax, nmax, gridtype
            integer, intent(in), optional :: n
            integer, intent(out), optional :: exitstatus
        end subroutine CilmPlusRhoH

        subroutine CilmMinusRhoH(cilm, gridin, lmax, nmax, mass, d, rho, &
                                 gridtype, w, zero, plx, n, dref, exitstatus)
            real*8, intent(in) ::   gridin(:,:), mass, rho(:,:)
            real*8, intent(in), optional :: w(:), zero(:), plx(:,:), dref
            real*8, intent(out) ::  cilm(:,:,:), d
            integer, intent(in) ::  lmax, nmax, gridtype
            integer, intent(in), optional :: n
            integer, intent(out), optional :: exitstatus
        end subroutine CilmMinusRhoH

        subroutine BAtoHilm(cilm, ba, gridglq, lmax, nmax, mass, r0, rho, &
                            gridtype, w, plx, zero, filter_type, filter_deg, &
                            lmax_calc, exitstatus)
            real*8, intent(out) ::  cilm(:,:,:)
            real*8, intent(in) ::   ba(:,:,:), gridglq(:,:), mass, r0, rho
            real*8, intent(in), optional :: plx(:,:), zero(:), w(:)
            integer, intent(in) ::  lmax, nmax, gridtype
            integer, intent(in), optional :: filter_type, filter_deg, lmax_calc
            integer, intent(out), optional :: exitstatus
        end subroutine BAtoHilm

        subroutine BAtoHilmRhoH(cilm, ba, gridglq, lmax, nmax, mass, r0, &
                                rho, gridtype, w, plx, zero, filter_type, &
                                filter_deg, lmax_calc, exitstatus)
            real*8, intent(out) ::  cilm(:,:,:)
            real*8, intent(in) ::   ba(:,:,:), gridglq(:,:), mass, r0, rho(:,:)
            real*8, intent(in), optional :: plx(:,:), zero(:), w(:)
            integer, intent(in) ::  lmax, nmax, gridtype
            integer, intent(in), optional :: filter_type, filter_deg, lmax_calc
            integer, intent(out), optional :: exitstatus
        end subroutine BAtoHilmRhoH

        real*8 function DownContFilterMA(l, half, r, d)
            integer, intent(in) ::  l, half
            real*8, intent(in) ::   r, d
        end function DownContFilterMA

        real*8 function DownContFilterMC(l, half, r, d)
            integer, intent(in) ::  l, half
            real*8, intent(in) ::  r, d
        end function DownContFilterMC

        real*8 function NormalGravity(geocentric_lat, gm, omega, a, b)
            real*8, intent(in) ::   geocentric_lat, gm, omega, a, b
        end function NormalGravity

        subroutine MakeMagGridDH(cilm, lmax, r0, a, f, rad_grid, theta_grid, &
                                 phi_grid, total_grid, n, sampling, lmax_calc, &
                                 pot_grid, exitstatus)
            real*8, intent(in) ::   cilm(:,:,:), r0, a, f
            real*8, intent(out) ::  rad_grid(:,:), theta_grid(:,:), &
                                    phi_grid(:,:), total_grid(:,:)
            real*8, intent(out), optional :: pot_grid(:,:)
            integer, intent(in) ::  lmax
            integer, intent(out) :: n
            integer, intent(in), optional :: sampling, lmax_calc
            integer, intent(out), optional :: exitstatus
        end subroutine MakeMagGridDH

        subroutine SHMagPowerSpectrum(c, a, r, lmax, spectra, exitstatus)
            real*8, intent(in) :: c(:,:,:)
            real*8, intent(in) :: a, r
            integer, intent(in) :: lmax
            real*8, intent(out) ::  spectra(:)
            integer, intent(out), optional :: exitstatus
        end subroutine SHMagPowerSpectrum

        real*8 function SHMagPowerL(c, a, r, l)
            real*8, intent(in) :: c(:,:,:)
            real*8, intent(in) :: a, r
            integer, intent(in) :: l
        end function SHMagPowerL

        subroutine MakeCircleCoord(coord, lat, lon, theta0, cinterval, cnum, &
                                   exitstatus)
            real*8, intent(in) :: lat, lon, theta0
            real*8, intent(out) :: coord(:,:)
            real*8, intent(in), optional :: cinterval
            integer, intent(out), optional :: cnum, exitstatus
        end subroutine MakeCircleCoord

        subroutine MakeEllipseCoord(coord, lat, lon, dec, A_theta, B_theta, &
                                    cinterval, cnum, exitstatus)
            real*8, intent(in) :: lat, lon, A_theta, B_theta, dec
            real*8, intent(out) :: coord(:,:)
            real*8, intent(in), optional :: cinterval
            integer, intent(out), optional :: cnum, exitstatus
        end subroutine MakeEllipseCoord

        subroutine Wigner3j(w3j, jmin, jmax, j2, j3, m1, m2, m3, exitstatus)
            integer, intent(in) :: j2, j3, m1, m2, m3
            integer, intent(out) :: jmin, jmax
            real*8, intent(out) :: w3j(:)
            integer, intent(out), optional :: exitstatus
        end subroutine Wigner3j

        real*8 function RandomN(idum)
            integer, parameter :: K4B=selected_int_kind(9)
            integer(K4B), intent(inout) ::  idum
        end function RandomN

        real*8 function RandomGaussian(idum)
            integer, parameter ::   K4B=selected_int_kind(9)
            integer(K4B), intent(inout) ::  idum
        end function RandomGaussian

        subroutine PreGLQ(x1, x2, n, zero, w, exitstatus)
            real*8, intent(in) ::   x1, x2
            real*8, intent(out) ::  zero(:), w(:)
            integer, intent(in) ::  n
            integer, intent(out), optional :: exitstatus
        end subroutine PreGLQ

        integer function NGLQ(degree)
            integer, intent(in) ::  degree
        end function NGLQ

        integer function NGLQSH(degree)
            integer, intent(in) ::  degree
        end function NGLQSH

        integer function NGLQSHN(degree, n)
            integer, intent(in) ::  degree, n
        end function NGLQSHN

        subroutine DHaj(n, aj, exitstatus)
            integer, intent(in) ::  n
            real*8, intent(out) ::  aj(:)
            integer, intent(out), optional :: exitstatus
        end subroutine DHaj

        integer function YilmIndexVector(i, l, m)
            integer, intent(in) ::  i, l, m
        end function YilmIndexVector

        subroutine EigValVecSym(ain, n, eig, evec, ul, k, exitstatus)
            real*8, intent(in) ::   ain(:,:)
            integer, intent(in) ::  n
            real*8, intent(out) ::  eig(:), evec(:,:)
            character, intent(in), optional ::  ul
            integer, intent(in), optional ::    k
            integer, intent(out), optional :: exitstatus
        end subroutine EigValVecSym

        subroutine EigValVecSymTri(ain, n, eig, evec, ul, exitstatus)
            real*8, intent(in) ::   ain(:,:)
            integer, intent(in) ::  n
            real*8, intent(out) ::  eig(:), evec(:,:)
            character, intent(in), optional :: ul
            integer, intent(out), optional :: exitstatus
        end subroutine EigValVecSymTri

        subroutine EigValSym(ain, n, eval, ul)
            real*8, intent(in) ::   ain(:,:)
            integer, intent(in) ::  n
            real*8, intent(out) ::  eval(:)
            character, intent(in), optional :: ul
        end subroutine EigValSym

    end interface

end module SHTOOLS
