module SHTOOLS
!------------------------------------------------------------------------------
!
!   This module contains an interface block defining all the routines
!   used in the archive SHTOOLS. These are necessary in order to use
!   implicitly shaped arrays with most subroutines.
!
!------------------------------------------------------------------------------
    implicit none

    integer, parameter :: CSPHASE_DEFAULT = 1
                          ! The default for all routines is to EXCLUDE
                          ! the Condon-Shortley phase of (-1)^m
                          ! in front of the Legendre functions.

    interface

        subroutine PlmBar(p, lmax, z, csphase, cnorm, exitstatus)
            integer, parameter :: dp = selected_real_kind(p=15)
            integer, intent(in) :: lmax
            real(dp), intent(out) :: p(:)
            real(dp), intent(in) :: z
            integer, intent(in), optional :: csphase, cnorm
            integer, intent(out), optional :: exitstatus
        end subroutine PlmBar

        subroutine PlmBar_d1(p, dp1, lmax, z, csphase, cnorm, exitstatus)
            integer, parameter :: dp = selected_real_kind(p=15)
            integer, intent(in) :: lmax
            real(dp), intent(out) :: p(:), dp1(:)
            real(dp), intent(in) :: z
            integer, intent(in), optional :: csphase, cnorm
            integer, intent(out), optional :: exitstatus
        end subroutine PlmBar_d1

        subroutine PlBar(p, lmax, z, exitstatus)
            integer, parameter :: dp = selected_real_kind(p=15)
            integer, intent(in) :: lmax
            real(dp), intent(out) :: p(:)
            real(dp), intent(in) :: z
            integer, intent(out), optional :: exitstatus
        end subroutine PlBar

        subroutine PlBar_d1(p, dp1, lmax, z, exitstatus)
            integer, parameter :: dp = selected_real_kind(p=15)
            integer, intent(in) :: lmax
            real(dp), intent(out) :: p(:), dp1(:)
            real(dp), intent(in) :: z
            integer, intent(out), optional :: exitstatus
        end subroutine PlBar_d1

        subroutine PlmON(p, lmax, z, csphase, cnorm, exitstatus)
            integer, parameter :: dp = selected_real_kind(p=15)
            integer, intent(in) :: lmax
            real(dp), intent(out) :: p(:)
            real(dp), intent(in) :: z
            integer, intent(in), optional :: csphase, cnorm
            integer, intent(out), optional :: exitstatus
        end subroutine PlmON

        subroutine PlmON_d1(p, dp1, lmax, z, csphase, cnorm, exitstatus)
            integer, parameter :: dp = selected_real_kind(p=15)
            integer, intent(in) :: lmax
            real(dp), intent(out) :: p(:), dp1(:)
            real(dp), intent(in) :: z
            integer, intent(in), optional :: csphase, cnorm
            integer, intent(out), optional :: exitstatus
        end subroutine PlmON_d1

        subroutine PlON(p, lmax, z, exitstatus)
            integer, parameter :: dp = selected_real_kind(p=15)
            integer, intent(in) :: lmax
            real(dp), intent(out) :: p(:)
            real(dp), intent(in) :: z
            integer, intent(out), optional :: exitstatus
        end subroutine PlON

        subroutine PlON_d1(p, dp1, lmax, z, exitstatus)
            integer, parameter :: dp = selected_real_kind(p=15)
            integer, intent(in) :: lmax
            real(dp), intent(out) :: p(:), dp1(:)
            real(dp), intent(in) :: z
            integer, intent(out), optional :: exitstatus
            end subroutine PlON_d1

        subroutine PlmSchmidt(p, lmax, z, csphase, cnorm, exitstatus)
            integer, parameter :: dp = selected_real_kind(p=15)
            integer, intent(in) :: lmax
            real(dp), intent(out) :: p(:)
            real(dp), intent(in) :: z
            integer, intent(in), optional :: csphase, cnorm
            integer, intent(out), optional :: exitstatus
        end subroutine PlmSchmidt

        subroutine PlmSchmidt_d1(p, dp1, lmax, z, csphase, cnorm, exitstatus)
            integer, parameter :: dp = selected_real_kind(p=15)
            integer, intent(in) :: lmax
            real(dp), intent(out) :: p(:), dp1(:)
            real(dp), intent(in) :: z
            integer, intent(in), optional :: csphase, cnorm
            integer, intent(out), optional :: exitstatus
        end subroutine PlmSchmidt_d1

        subroutine PlSchmidt(p, lmax, z, exitstatus)
            integer, parameter :: dp = selected_real_kind(p=15)
            integer, intent(in) :: lmax
            real(dp), intent(out) :: p(:)
            real(dp), intent(in) :: z
            integer, intent(out), optional :: exitstatus
        end subroutine PlSchmidt

        subroutine PlSchmidt_d1(p, dp1, lmax, z, exitstatus)
            integer, parameter :: dp = selected_real_kind(p=15)
            integer, intent(in) :: lmax
            real(dp), intent(out) :: p(:), dp1(:)
            real(dp), intent(in) :: z
            integer, intent(out), optional :: exitstatus
        end subroutine PlSchmidt_d1

        subroutine PLegendreA(p,lmax,z, csphase, exitstatus)
            integer, parameter :: dp = selected_real_kind(p=15)
            integer, intent(in) :: lmax
            real(dp), intent(out) :: p(:)
            real(dp), intent(in) :: z
            integer, intent(in), optional :: csphase
            integer, intent(out), optional :: exitstatus
        end subroutine PLegendreA

        subroutine PLegendreA_d1(p, dp1, lmax, z, csphase, exitstatus)
            integer, parameter :: dp = selected_real_kind(p=15)
            integer, intent(in) :: lmax
            real(dp), intent(out) :: p(:), dp1(:)
            real(dp), intent(in) :: z
            integer, intent(in), optional :: csphase
            integer, intent(out), optional :: exitstatus
        end subroutine PLegendreA_d1

        subroutine PLegendre(p,lmax,z, exitstatus)
            integer, parameter :: dp = selected_real_kind(p=15)
            integer, intent(in) :: lmax
            real(dp), intent(out) :: p(:)
            real(dp), intent(in) :: z
            integer, intent(out), optional :: exitstatus
        end subroutine PLegendre

        subroutine PLegendre_d1(p, dp1, lmax, z, exitstatus)
            integer, parameter :: dp = selected_real_kind(p=15)
            integer, intent(in) :: lmax
            real(dp), intent(out) :: p(:), dp1(:)
            real(dp), intent(in) :: z
            integer, intent(out), optional :: exitstatus
        end subroutine PLegendre_d1

        integer function PlmIndex(l,m)
            integer, intent(in) :: l, m
        end function PlmIndex

        subroutine SHExpandDH(grid, n, cilm, lmax, norm, sampling, &
                              csphase, lmax_calc, exitstatus)
            integer, parameter :: dp = selected_real_kind(p=15)
            real(dp), intent(in) :: grid(:,:)
            real(dp), intent(out) :: cilm(:,:,:)
            integer, intent(in) :: n
            integer, intent(out) :: lmax
            integer, intent(in), optional :: norm, sampling, csphase, lmax_calc
            integer, intent(out), optional :: exitstatus
        end subroutine SHExpandDH

        subroutine MakeGridDH(griddh, n, cilm, lmax, norm, sampling, &
                              csphase, lmax_calc, extend, exitstatus)
            integer, parameter :: dp = selected_real_kind(p=15)
            real(dp), intent(in) :: cilm(:,:,:)
            real(dp), intent(out) :: griddh(:,:)
            integer, intent(in) :: lmax
            integer, intent(out) :: n
            integer, intent(in), optional :: norm, sampling, csphase, &
                                             lmax_calc, extend
            integer, intent(out), optional :: exitstatus
        end subroutine MakeGridDH

        subroutine SHExpandDHC(grid, n, cilm, lmax, norm, sampling, &
                               csphase, lmax_calc, exitstatus)
            integer, parameter :: dp = selected_real_kind(p=15)
            complex(dp), intent(in) :: grid(:,:)
            complex(dp), intent(out) :: cilm(:,:,:)
            integer, intent(in) :: n
            integer, intent(out) :: lmax
            integer, intent(in), optional :: norm, sampling, csphase, lmax_calc
            integer, intent(out), optional :: exitstatus
        end subroutine SHExpandDHC

        subroutine MakeGridDHC(griddh, n, cilm, lmax, norm, sampling, &
                               csphase, lmax_calc, extend, exitstatus)
            integer, parameter :: dp = selected_real_kind(p=15)
            complex(dp), intent(in) :: cilm(:,:,:)
            complex(dp), intent(out) :: griddh(:,:)
            integer, intent(in) :: lmax
            integer, intent(out) :: n
            integer, intent(in), optional :: norm, sampling, csphase, &
                                             lmax_calc, extend
            integer, intent(out), optional :: exitstatus
        end subroutine MakeGridDHC

        subroutine SHGLQ(lmax, zero, w, plx, norm, csphase, cnorm, exitstatus)
            integer, parameter :: dp = selected_real_kind(p=15)
            integer, intent(in) :: lmax
            real(dp), intent(out) :: zero(:), w(:)
            real(dp), intent(out), optional :: plx(:,:)
            integer, intent(in), optional :: norm, csphase, cnorm
            integer, intent(out), optional :: exitstatus
        end subroutine SHGLQ

        subroutine SHExpandGLQ(cilm, lmax, gridglq, w, plx, zero, norm, &
                               csphase, lmax_calc, exitstatus)
            integer, parameter :: dp = selected_real_kind(p=15)
            real(dp), intent(in) :: w(:), gridglq(:,:)
            real(dp), intent(in), optional :: plx(:,:), zero(:)
            real(dp), intent(out) :: cilm(:,:,:)
            integer, intent(in) :: lmax
            integer, intent(in), optional :: norm, csphase, lmax_calc
            integer, intent(out), optional :: exitstatus
        end subroutine SHExpandGLQ

        subroutine MakeGridGLQ(gridglq, cilm, lmax, plx, zero, norm, &
                               csphase, lmax_calc, extend, exitstatus)
            integer, parameter :: dp = selected_real_kind(p=15)
            real(dp), intent(in) :: cilm(:,:,:)
            real(dp), intent(in), optional :: plx(:,:), zero(:)
            real(dp), intent(out) :: gridglq(:,:)
            integer, intent(in) :: lmax
            integer, intent(in), optional :: norm, csphase, lmax_calc, extend
            integer, intent(out), optional :: exitstatus
        end subroutine MakeGridGLQ

        subroutine SHExpandGLQC(cilm, lmax, gridglq, w, plx, zero, norm, &
                                csphase, lmax_calc, exitstatus)
            integer, parameter :: dp = selected_real_kind(p=15)
            real(dp), intent(in) :: w(:)
            complex(dp), intent(in) :: gridglq(:,:)
            real(dp), intent(in), optional :: plx(:,:), zero(:)
            complex(dp), intent(out) :: cilm(:,:,:)
            integer, intent(in) :: lmax
            integer, intent(in), optional :: norm, csphase, lmax_calc
            integer, intent(out), optional :: exitstatus
        end subroutine SHExpandGLQC

        subroutine MakeGridGLQC(gridglq, cilm, lmax, plx, zero, norm, &
                                csphase, lmax_calc, extend, exitstatus)
            integer, parameter :: dp = selected_real_kind(p=15)
            complex(dp), intent(in) :: cilm(:,:,:)
            real(dp), intent(in), optional :: plx(:,:), zero(:)
            complex(dp), intent(out) :: gridglq(:,:)
            integer, intent(in) :: lmax
            integer, intent(in), optional :: norm, csphase, lmax_calc, extend
            integer, intent(out), optional :: exitstatus
        end subroutine MakeGridGLQC

        subroutine GLQGridCoord(latglq, longlq, lmax, nlat, nlong, extend, &
                                exitstatus)
            integer, parameter :: dp = selected_real_kind(p=15)
            integer, intent(in) :: lmax
            integer, intent(out) :: nlat, nlong
            real(dp), intent(out) :: latglq(:), longlq(:)
            integer, intent(in), optional :: extend
            integer, intent(out), optional :: exitstatus
        end subroutine GLQGridCoord

        subroutine SHExpandLSQ(cilm, d, lat, lon, nmax, lmax, norm, &
                               chi2, csphase, weights, exitstatus)
            integer, parameter :: dp = selected_real_kind(p=15)
            real(dp), intent(in) :: d(:), lat(:), lon(:)
            real(dp), intent(out) :: cilm(:,:,:)
            integer, intent(in) :: nmax, lmax
            integer, intent(in), optional :: norm, csphase
            real(dp), intent(out), optional :: chi2
            real(dp), intent(in), optional :: weights(:)
            integer, intent(out), optional :: exitstatus
        end subroutine SHExpandLSQ

        subroutine MakeGrid2d(grid, cilm, lmax, interval, nlat, nlong, &
                              norm, csphase, f, a, north, south, east, west, &
                              dealloc, exitstatus)
            integer, parameter :: dp = selected_real_kind(p=15)
            real(dp), intent(in) :: cilm(:,:,:), interval
            real(dp), intent(out) :: grid(:,:)
            integer, intent(in) :: lmax
            integer, intent(out) :: nlat, nlong
            integer, intent(in), optional :: norm, csphase, dealloc
            real(dp), intent(in), optional :: f, a, north, south, east, west
            integer, intent(out), optional :: exitstatus
        end subroutine MakeGrid2D

        function MakeGridPoint(cilm, lmax, lat, lon, norm, csphase, dealloc)
            integer, parameter :: dp = selected_real_kind(p=15)
            real(dp) :: MakeGridPoint
            real(dp), intent(in) :: cilm(:,:,:), lat, lon
            integer, intent(in) :: lmax
            integer, intent(in), optional :: norm, csphase, dealloc
        end function MakeGridPoint

        function MakeGridPointC(cilm, lmax, lat, lon, norm, csphase, dealloc)
            integer, parameter :: dp = selected_real_kind(p=15)
            complex(dp) :: MakeGridPointC
            complex(dp), intent(in) :: cilm(:,:,:)
            real(dp), intent(in) :: lat, lon
            integer, intent(in) :: lmax
            integer, intent(in), optional :: norm, csphase, dealloc
        end function MakeGridPointC

        subroutine SHMultiply(cilmout, cilm1, lmax1, cilm2, lmax2, precomp, &
                              norm, csphase, exitstatus)
            integer, parameter :: dp = selected_real_kind(p=15)
            real(dp), intent(out) :: cilmout(:,:,:)
            real(dp), intent(in) :: cilm1(:,:,:), cilm2(:,:,:)
            integer, intent(in) :: lmax1, lmax2
            integer, intent(in), optional :: precomp, norm, csphase
            integer, intent(out), optional :: exitstatus
        end subroutine SHMultiply

        subroutine SHRead(filename, cilm, lmax, skip, header, error, &
                          exitstatus)
            integer, parameter :: dp = selected_real_kind(p=15)
            character(*), intent(in) :: filename
            integer, intent(out) :: lmax
            real(dp), intent(out) :: cilm(:,:,:)
            real(dp), intent(out), optional :: header(:), error(:,:,:)
            integer, intent(in), optional :: skip
            integer, intent(out), optional :: exitstatus
        end subroutine SHRead

        subroutine SHRead2(filename, cilm, lmax, gm, r0_pot, error, dot, &
                           doystart, doyend, epoch, exitstatus)
            integer, parameter :: dp = selected_real_kind(p=15)
            character(*), intent(in) :: filename
            integer, intent(out) :: lmax
            real(dp), intent(out) :: cilm(:,:,:), gm, r0_pot
            real(dp), intent(out), optional :: error(:,:,:), dot(:,:,:), &
                                               doystart, doyend, epoch
            integer, intent(out), optional :: exitstatus
        end subroutine SHRead2

        subroutine SHReadJPL(filename, cilm, lmax, error, gm, formatstring, &
                             exitstatus)
            integer, parameter :: dp = selected_real_kind(p=15)
            character(*), intent(in) :: filename
            integer, intent(in) :: lmax
            real(dp), intent(out) :: cilm(:,:,:)
            real(dp), intent(out), optional :: error(:,:,:), gm(2)
            character(6), intent(in), optional :: formatstring
            integer, intent(out), optional :: exitstatus
        end subroutine SHReadJPL

        subroutine SHCilmToVector(cilm, vector, lmax, exitstatus)
            integer, parameter :: dp = selected_real_kind(p=15)
            real(dp), intent(in) :: cilm(:,:,:)
            real(dp), intent(out) :: vector(:)
            integer, intent(in) :: lmax
            integer, intent(out), optional :: exitstatus
        end subroutine SHCilmToVector

        subroutine SHVectorToCilm(vector, cilm, lmax, exitstatus)
            integer, parameter :: dp = selected_real_kind(p=15)
            real(dp), intent(out) :: cilm(:,:,:)
            real(dp), intent(in) :: vector(:)
            integer, intent(in) :: lmax
            integer, intent(out), optional :: exitstatus
        end subroutine SHVectorToCilm

        subroutine SHCilmToCindex(cilm, cindex, degmax, exitstatus)
            integer, parameter :: dp = selected_real_kind(p=15)
            real(dp), intent(in) :: cilm(:,:,:)
            real(dp), intent(out) :: cindex(:,:)
            integer, intent(in), optional :: degmax
            integer, intent(out), optional :: exitstatus
        end subroutine SHCilmToCindex

        subroutine SHCindexToCilm(cindex, cilm, degmax, exitstatus)
            integer, parameter :: dp = selected_real_kind(p=15)
            real(dp), intent(out) :: cilm(:,:,:)
            real(dp), intent(in) :: cindex(:,:)
            integer, intent(in), optional :: degmax
            integer, intent(out), optional :: exitstatus
        end subroutine SHCindexToCilm

        subroutine SHrtoc(rcilm, ccilm, degmax, convention, switchcs, &
                          exitstatus)
            integer, parameter :: dp = selected_real_kind(p=15)
            real(dp), intent(in) :: rcilm(:,:,:)
            real(dp), intent(out) :: ccilm(:,:,:)
            integer, intent(in), optional :: degmax, convention, switchcs
            integer, intent(out), optional :: exitstatus
        end subroutine SHrtoc

        subroutine SHctor(ccilm, rcilm, degmax, convention, switchcs, &
                          exitstatus)
            integer, parameter :: dp = selected_real_kind(p=15)
            real(dp), intent(in) :: ccilm(:,:,:)
            real(dp), intent(out) :: rcilm(:,:,:)
            integer, intent(in), optional :: degmax, convention, switchcs
            integer, intent(out), optional :: exitstatus
        end subroutine SHctor

        subroutine djpi2(dj, lmax, exitstatus)
            integer, parameter :: dp = selected_real_kind(p=15)
            integer, intent(in) :: lmax
            real(dp), intent(out) :: dj(:,:,:)
            integer, intent(out), optional :: exitstatus
        end subroutine djpi2

        subroutine SHRotateCoef(x, cof, rcof, dj, lmax, exitstatus)
            integer, parameter :: dp = selected_real_kind(p=15)
            real(dp), intent(in) :: cof(:,:), dj(:,:,:), x(3)
            real(dp), intent(out) :: rcof(:,:)
            integer, intent(in) :: lmax
            integer, intent(out), optional :: exitstatus
        end subroutine SHRotateCoef

        subroutine SHRotateRealCoef(cilmrot, cilm, lmax, x, dj, exitstatus)
            integer, parameter :: dp = selected_real_kind(p=15)
            real(dp), intent(in) :: cilm(:,:,:), x(:), dj(:,:,:)
            real(dp), intent(out) :: cilmrot(:,:,:)
            integer, intent(in) :: lmax
            integer, intent(out), optional :: exitstatus
        end subroutine SHRotateRealCoef

        function SHPowerL(cilm, l)
            integer, parameter :: dp = selected_real_kind(p=15)
            real(dp) :: SHPowerL
            real(dp), intent(in) :: cilm(:,:,:)
            integer, intent(in) :: l
        end function SHPowerL

        function SHPowerDensityL(cilm, l)
            integer, parameter :: dp = selected_real_kind(p=15)
            real(dp) :: SHPowerDensityL
            real(dp), intent(in) :: cilm(:,:,:)
            integer, intent(in) :: l
        end function SHPowerDensityL

        function SHCrossPowerL(cilm1, cilm2, l)
            integer, parameter :: dp = selected_real_kind(p=15)
            real(dp) :: SHCrossPowerL
            real(dp), intent(in) :: cilm1(:,:,:), cilm2(:,:,:)
            integer, intent(in) :: l
        end function SHCrossPowerL

        function SHCrossPowerDensityL(cilm1, cilm2, l)
            integer, parameter :: dp = selected_real_kind(p=15)
            real(dp) :: SHCrossPowerDensityL
            real(dp), intent(in) :: cilm1(:,:,:), cilm2(:,:,:)
            integer, intent(in) :: l
        end function SHCrossPowerDensityL

        subroutine SHPowerSpectrum(cilm, lmax, spectra, exitstatus)
            integer, parameter :: dp = selected_real_kind(p=15)
            real(dp), intent(in) :: cilm(:,:,:)
            integer, intent(in) :: lmax
            real(dp), intent(out) :: spectra(:)
            integer, intent(out), optional :: exitstatus
        end subroutine SHPowerSpectrum

        subroutine SHPowerSpectrumDensity(cilm, lmax, spectra, exitstatus)
            integer, parameter :: dp = selected_real_kind(p=15)
            real(dp), intent(in) :: cilm(:,:,:)
            integer, intent(in) :: lmax
            real(dp), intent(out) :: spectra(:)
            integer, intent(out), optional :: exitstatus
        end subroutine SHPowerSpectrumDensity

        subroutine SHCrossPowerSpectrum(cilm1, cilm2, lmax, cspectra, &
                                        exitstatus)
            integer, parameter :: dp = selected_real_kind(p=15)
            real(dp), intent(in) :: cilm1(:,:,:), cilm2(:,:,:)
            integer, intent(in) :: lmax
            real(dp), intent(out) :: cspectra(:)
            integer, intent(out), optional :: exitstatus
        end subroutine SHCrossPowerSpectrum

        subroutine SHCrossPowerSpectrumDensity(cilm1, cilm2, lmax, cspectra, &
                                               exitstatus)
            integer, parameter :: dp = selected_real_kind(p=15)
            real(dp), intent(in) :: cilm1(:,:,:), cilm2(:,:,:)
            integer, intent(in) :: lmax
            real(dp), intent(out) :: cspectra(:)
            integer, intent(out), optional :: exitstatus
        end subroutine SHCrossPowerSpectrumDensity

        subroutine SHAdmitCorr(gilm, tilm, lmax, admit, corr, admit_error, &
                               exitstatus)
            integer, parameter :: dp = selected_real_kind(p=15)
            real(dp), intent(in) :: gilm(:,:,:), tilm(:,:,:)
            integer, intent(in) :: lmax
            real(dp), intent(out) :: admit(:), corr(:)
            real(dp), intent(out), optional :: admit_error(:)
            integer, intent(out), optional :: exitstatus
        end subroutine SHAdmitCorr

        function SHConfidence(l_conf, r)
            integer, parameter :: dp = selected_real_kind(p=15)
            real(dp) :: SHConfidence
            real(dp), intent(in) :: r
            integer, intent(in) :: l_conf
        end function SHConfidence

        function SHPowerLC(cilm, l)
            integer, parameter :: dp = selected_real_kind(p=15)
            real(dp) :: SHPowerLC
            complex(dp), intent(in) :: cilm(:,:,:)
            integer, intent(in) :: l
        end function SHPowerLC

        function SHPowerDensityLC(cilm, l)
            integer, parameter :: dp = selected_real_kind(p=15)
            real(dp) :: SHPowerDensityLC
            complex(dp), intent(in) :: cilm(:,:,:)
            integer, intent(in) :: l
        end function SHPowerDensityLC

        function SHCrossPowerLC(cilm1, cilm2, l)
            integer, parameter :: dp = selected_real_kind(p=15)
            complex(dp) :: SHCrossPowerLC
            complex(dp), intent(in) :: cilm1(:,:,:), cilm2(:,:,:)
            integer, intent(in) :: l
        end function SHCrossPowerLC

        function SHCrossPowerDensityLC(cilm1, cilm2, l)
            integer, parameter :: dp = selected_real_kind(p=15)
            complex(dp) :: SHCrossPowerDensityLC
            complex(dp), intent(in) :: cilm1(:,:,:), cilm2(:,:,:)
            integer, intent(in) :: l
        end function SHCrossPowerDensityLC

        subroutine SHPowerSpectrumC(cilm, lmax, spectra, exitstatus)
            integer, parameter :: dp = selected_real_kind(p=15)
            complex(dp), intent(in) :: cilm(:,:,:)
            integer, intent(in) :: lmax
            real(dp), intent(out) :: spectra(:)
            integer, intent(out), optional :: exitstatus
        end subroutine SHPowerSpectrumC

        subroutine SHPowerSpectrumDensityC(cilm, lmax, spectra, exitstatus)
            integer, parameter :: dp = selected_real_kind(p=15)
            complex(dp), intent(in) :: cilm(:,:,:)
            integer, intent(in) :: lmax
            real(dp), intent(out) :: spectra(:)
            integer, intent(out), optional :: exitstatus
        end subroutine SHPowerSpectrumDensityC

        subroutine SHCrossPowerSpectrumC(cilm1, cilm2, lmax, cspectra, &
                                         exitstatus)
            integer, parameter :: dp = selected_real_kind(p=15)
            complex(dp), intent(in) :: cilm1(:,:,:), cilm2(:,:,:)
            integer, intent(in) :: lmax
            complex(dp), intent(out) :: cspectra(:)
            integer, intent(out), optional :: exitstatus
        end subroutine SHCrossPowerSpectrumC

        subroutine SHCrossPowerSpectrumDensityC(cilm1, cilm2, lmax, cspectra, &
                                                exitstatus)
            integer, parameter :: dp = selected_real_kind(p=15)
            complex(dp), intent(in) :: cilm1(:,:,:), cilm2(:,:,:)
            integer, intent(in) :: lmax
            complex(dp), intent(out) :: cspectra(:)
            integer, intent(out), optional :: exitstatus
        end subroutine SHCrossPowerSpectrumDensityC

        subroutine SHMultiTaperSE(mtse, sd, cilm, lmax, tapers, taper_order, &
                                  lmaxt, k, alpha, lat, lon, taper_wt, norm, &
                                  csphase, exitstatus)
            integer, parameter :: dp = selected_real_kind(p=15)
            real(dp), intent(out) :: mtse(:), sd(:)
            real(dp), intent(in) :: cilm(:,:,:), tapers(:,:)
            integer, intent(in) :: lmax, lmaxt, k, taper_order(:)
            real(dp), intent(in), optional :: alpha(:), lat, lon, taper_wt(:)
            integer, intent(in), optional :: csphase, norm
            integer, intent(out), optional :: exitstatus
        end subroutine SHMultiTaperSE

        subroutine SHMultiTaperCSE(mtse, sd, cilm1, lmax1, cilm2, lmax2, &
                                   tapers, taper_order, lmaxt, k, alpha, lat, &
                                   lon, taper_wt, norm, csphase, exitstatus)
            integer, parameter :: dp = selected_real_kind(p=15)
            real(dp), intent(out) :: mtse(:), sd(:)
            real(dp), intent(in) :: cilm1(:,:,:), cilm2(:,:,:), tapers(:,:)
            integer, intent(in) :: lmax1, lmax2, lmaxt, k, taper_order(:)
            real(dp), intent(in), optional :: alpha(:), lat, lon, taper_wt(:)
            integer, intent(in), optional :: csphase, norm
            integer, intent(out), optional :: exitstatus
        end subroutine SHMultiTaperCSE

        subroutine SHLocalizedAdmitCorr(tapers, taper_order, lwin, lat, lon, &
                                        gilm, tilm, lmax, admit, corr, k, &
                                        admit_error, corr_error, taper_wt, &
                                        mtdef, k1linsig, exitstatus)
            integer, parameter :: dp = selected_real_kind(p=15)
            real(dp), intent(in) :: tapers(:,:), lat, lon, gilm(:,:,:), &
                                    tilm(:,:,:)
            integer, intent(in) :: lwin, lmax, k, taper_order(:)
            real(dp), intent(out) :: admit(:), corr(:)
            real(dp), intent(out), optional :: admit_error(:), corr_error(:)
            integer, intent(in), optional :: mtdef, k1linsig
            real(dp), intent(in), optional :: taper_wt(:)
            integer, intent(out), optional :: exitstatus
        end subroutine SHLocalizedAdmitCorr

        subroutine SHReturnTapers(theta0, lmax, tapers, eigenvalues, &
                                  taper_order, degrees, exitstatus)
            integer, parameter :: dp = selected_real_kind(p=15)
            real(dp), intent(in) :: theta0
            integer, intent(in) :: lmax
            real(dp), intent(out) :: tapers(:,:), eigenvalues(:)
            integer, intent(out) :: taper_order(:)
            integer, intent(in), optional :: degrees(:)
            integer, intent(out), optional :: exitstatus
        end subroutine SHReturnTapers

        subroutine SHReturnTapersM(theta0, lmax, m, tapers, eigenvalues, &
                                   shannon, degrees, ntapers, exitstatus)
            integer, parameter :: dp = selected_real_kind(p=15)
            real(dp), intent(in) :: theta0
            integer, intent(in) :: lmax, m
            real(dp), intent(out) :: tapers(:,:), eigenvalues(:)
            real(dp), intent(out), optional :: shannon
            integer, intent(in), optional :: degrees(:)
            integer, intent(out), optional :: ntapers
            integer, intent(out), optional :: exitstatus
        end subroutine SHReturnTapersM

        subroutine ComputeDm(dllm, lmax, m, theta0, degrees, exitstatus)
            integer, parameter :: dp = selected_real_kind(p=15)
            real(dp), intent(out) :: dllm(:,:)
            real(dp), intent(in) :: theta0
            integer, intent(in) :: lmax, m
            integer, intent(in), optional :: degrees(:)
            integer, intent(out), optional :: exitstatus
        end subroutine ComputeDm

        subroutine ComputeDG82(dG82, lmax, m, theta0, exitstatus)
            integer, parameter :: dp = selected_real_kind(p=15)
            real(dp), intent(out) :: dG82(:,:)
            real(dp), intent(in) :: theta0
            integer, intent(in) :: lmax, m
            integer, intent(out), optional :: exitstatus
        end subroutine ComputeDG82

        function SHFindLWin(theta0, m, alpha, taper_number)
            integer, parameter :: dp = selected_real_kind(p=15)
            integer :: SHFindLWin
            real(dp), intent(in) :: theta0, alpha
            integer, intent(in) :: m
            integer, intent(in), optional :: taper_number
        end function SHFindLWin

        subroutine SHBiasK(tapers, lwin, k, incspectra, ldata, &
                           outcspectra, taper_wt, save_cg, exitstatus)
            integer, parameter :: dp = selected_real_kind(p=15)
            real(dp), intent(in) :: tapers(:,:), incspectra(:)
            real(dp), intent(out) :: outcspectra(:)
            integer, intent(in) :: lwin, ldata, k
            real(dp), intent(in), optional :: taper_wt(:)
            integer, intent(in), optional :: save_cg
            integer, intent(out), optional :: exitstatus
        end subroutine SHBiasK

        subroutine SHMTCouplingMatrix(Mmt, lmax, tapers_power, lwin, k, &
                                      taper_wt, exitstatus)
            integer, parameter :: dp = selected_real_kind(p=15)
            real(dp), intent(out) :: Mmt(:,:)
            real(dp), intent(in) :: tapers_power(:,:)
            real(dp), intent(in), optional :: taper_wt(:)
            integer, intent(in) :: lmax, k, lwin
            integer, intent(out), optional :: exitstatus
        end subroutine

        subroutine SHBiasAdmitCorr(sgt, sgg, stt, lmax, tapers, lwin, k, &
                                   admit, corr, mtdef, taper_wt, exitstatus)
            integer, parameter :: dp = selected_real_kind(p=15)
            real(dp), intent(in) :: sgt(:), sgg(:), stt(:), tapers(:,:)
            integer, intent(in) :: lmax, lwin, k
            real(dp), intent(out) :: admit(:), corr(:)
            integer, intent(in), optional :: mtdef
            real(dp), intent(in), optional :: taper_wt(:)
            integer, intent(out), optional :: exitstatus
        end subroutine SHBiasAdmitCorr

        subroutine SHMTDebias (mtdebias, mtspectra, lmax, tapers, lwin, k, &
                               nl, lmid, n, taper_wt, exitstatus)
            integer, parameter :: dp = selected_real_kind(p=15)
            real(dp), intent(out) :: mtdebias(:,:), lmid(:)
            real(dp), intent(in) :: mtspectra(:,:), tapers(:,:)
            real(dp), intent(in), optional :: taper_wt(:)
            integer, intent(in) :: lmax, k, lwin, nl
            integer, intent(out) :: n
            integer, intent(out), optional :: exitstatus
        end subroutine SHMTDebias

        subroutine SHMTVarOpt(l, tapers, taper_order, lwin, kmax, Sff, &
                              var_opt, var_unit, weight_opt, &
                              unweighted_covar, nocross, exitstatus)
            integer, parameter :: dp = selected_real_kind(p=15)
            real(dp), intent(in) :: tapers(:,:), Sff(:)
            real(dp), intent(out) :: var_opt(:), var_unit(:)
            integer, intent(in) :: l, lwin, kmax, taper_order(:)
            real(dp), intent(out), optional :: weight_opt(:,:), &
                                               unweighted_covar(:,:)
            integer, intent(in), optional :: nocross
            integer, intent(out), optional :: exitstatus
        end subroutine SHMTVarOpt

        subroutine SHMTVar(l, tapers, taper_order, lwin, kmax, Sff, &
                              variance, taper_wt, unweighted_covar, nocross, &
                              exitstatus)
            integer, parameter :: dp = selected_real_kind(p=15)
            real(dp), intent(in) :: tapers(:,:), Sff(:)
            real(dp), intent(out) :: variance
            integer, intent(in) :: l, lwin, kmax, taper_order(:)
            real(dp), intent(in), optional :: taper_wt(:)
            real(dp), intent(out), optional :: unweighted_covar(:,:)
            integer, intent(in), optional :: nocross
            integer, intent(out), optional :: exitstatus
        end subroutine SHMTVar

        function SHSjkPG(incspectra, l, m, mprime, hj_real, hk_real, mj, mk, &
                         lwin, hkcc)
            integer, parameter :: dp = selected_real_kind(p=15)
            complex(dp) :: SHSjkPG
            real(dp), intent(in) :: incspectra(:), hj_real(:), hk_real(:)
            integer, intent(in) :: lwin, l, m, mprime, mj, mk, hkcc
        end function SHSjkPG

        subroutine SHReturnTapersMap(tapers, eigenvalues, dh_mask, n_dh, &
                                     lmax, sampling, ntapers, degrees, &
                                     exitstatus)
            integer, parameter :: dp = selected_real_kind(p=15)
            real(dp), intent(out) :: tapers(:,:), eigenvalues(:)
            integer, intent(in) :: dh_mask(:,:), n_dh, lmax, sampling
            integer, intent(in), optional :: ntapers, degrees(:)
            integer, intent(out), optional :: exitstatus
        end subroutine SHReturnTapersMap

        subroutine SHBiasKMask(tapers, lwin, k, incspectra, ldata, &
                           outcspectra, taper_wt, save_cg, exitstatus)
            integer, parameter :: dp = selected_real_kind(p=15)
            real(dp), intent(in) :: tapers(:,:), incspectra(:)
            real(dp), intent(out) :: outcspectra(:)
            integer, intent(in) :: lwin, ldata, k
            real(dp), intent(in), optional :: taper_wt(:)
            integer, intent(in), optional :: save_cg
            integer, intent(out), optional :: exitstatus
        end subroutine SHBiasKMask

        subroutine SHMultiTaperMaskSE(mtse, sd, cilm, lmax, tapers, &
                                      lmaxt, k, taper_wt, norm, csphase, &
                                      exitstatus)
            integer, parameter :: dp = selected_real_kind(p=15)
            real(dp), intent(out) :: mtse(:), sd(:)
            real(dp), intent(in) :: cilm(:,:,:), tapers(:,:)
            integer, intent(in) :: lmax, lmaxt, k
            real(dp), intent(in), optional :: taper_wt(:)
            integer, intent(in), optional :: csphase, norm
            integer, intent(out), optional :: exitstatus
        end subroutine SHMultiTaperMaskSE

        subroutine SHMultiTaperMaskCSE(mtse, sd, cilm1, lmax1, cilm2, lmax2, &
                                       tapers, lmaxt, k, taper_wt, norm, &
                                       csphase, exitstatus)
            integer, parameter :: dp = selected_real_kind(p=15)
            real(dp), intent(out) :: mtse(:), sd(:)
            real(dp), intent(in) :: cilm1(:,:,:), cilm2(:,:,:), tapers(:,:)
            integer, intent(in) :: lmax1, lmax2, lmaxt, k
            real(dp), intent(in), optional :: taper_wt(:)
            integer, intent(in), optional :: csphase, norm
            integer, intent(out), optional :: exitstatus
        end subroutine SHMultiTaperMaskCSE

        subroutine ComputeDMap(Dij, dh_mask, n_dh, lmax, sampling, degrees, &
                               exitstatus)
            integer, parameter :: dp = selected_real_kind(p=15)
            real(dp), intent(out) :: Dij(:,:)
            integer, intent(in) :: dh_mask(:,:), n_dh, lmax
            integer, intent(in), optional :: sampling, degrees(:)
            integer, intent(out), optional :: exitstatus
        end subroutine ComputeDMap

        subroutine Curve2Mask(dhgrid, n, sampling, profile, nprofile, NP, &
                              extend, exitstatus)
            integer, parameter :: dp = selected_real_kind(p=15)
            integer, intent(out) :: dhgrid(:,:)
            real(dp), intent(in) :: profile(:,:)
            integer, intent(in) :: n, sampling, nprofile, np
            integer, intent(in), optional :: extend
            integer, intent(out), optional :: exitstatus
        end subroutine Curve2Mask

        subroutine SHBias(Shh, lwin, incspectra, ldata, outcspectra, save_cg, &
                          exitstatus)
            integer, parameter :: dp = selected_real_kind(p=15)
            real(dp), intent(in) :: Shh(:), incspectra(:)
            real(dp), intent(out) :: outcspectra(:)
            integer, intent(in) :: lwin, ldata
            integer, intent(in), optional :: save_cg
            integer, intent(out), optional :: exitstatus
        end subroutine SHBias

        subroutine SphericalCapCoef(coef, theta, lmax, exitstatus)
            integer, parameter :: dp = selected_real_kind(p=15)
            real(dp), intent(out) :: coef(:)
            real(dp), intent(in) :: theta
            integer, intent(in), optional :: lmax
            integer, intent(out), optional :: exitstatus
        end subroutine SphericalCapCoef

        subroutine MakeGravGridDH(cilm, lmax, gm, r0, a, f, rad, theta, phi, &
                                  total, n, sampling, lmax_calc, omega, &
                                  normal_gravity, pot, extend, exitstatus)
            integer, parameter :: dp = selected_real_kind(p=15)
            real(dp), intent(in) :: cilm(:,:,:), gm, r0, a, f
            real(dp), intent(out) :: rad(:,:), theta(:,:), phi(:,:), total(:,:)
            real(dp), intent(in), optional :: omega
            real(dp), intent(out), optional :: pot(:,:)
            integer, intent(in) :: lmax
            integer, intent(out) :: n
            integer, intent(in), optional :: sampling, lmax_calc, &
                                             normal_gravity, extend
            integer, intent(out), optional :: exitstatus
        end subroutine MakeGravGridDH

        subroutine MakeGravGradGridDH(cilm, lmax, gm, r0, a, f, vxx, vyy, &
                                      vzz, vxy, vxz, vyz, n, sampling, &
                                      lmax_calc, extend, exitstatus)
            integer, parameter :: dp = selected_real_kind(p=15)
            real(dp), intent(in) :: cilm(:,:,:), gm, r0, a, f
            real(dp), intent(out) :: vxx(:,:), vyy(:,:), vzz(:,:), vxy(:,:), &
                                     vxz(:,:), vyz(:,:)
            integer, intent(in) :: lmax
            integer, intent(out) :: n
            integer, intent(in), optional :: sampling, lmax_calc, extend
            integer, intent(out), optional :: exitstatus
        end subroutine MakeGravGradGridDH

        subroutine MakeMagGradGridDH(cilm, lmax, r0, a, f, vxx, vyy, &
                                     vzz, vxy, vxz, vyz, n, sampling, &
                                     lmax_calc, extend, exitstatus)
            integer, parameter :: dp = selected_real_kind(p=15)
            real(dp), intent(in) :: cilm(:,:,:), r0, a, f
            real(dp), intent(out) :: vxx(:,:), vyy(:,:), vzz(:,:), vxy(:,:), &
                                     vxz(:,:), vyz(:,:)
            integer, intent(in) :: lmax
            integer, intent(out) :: n
            integer, intent(in), optional :: sampling, lmax_calc, extend
            integer, intent(out), optional :: exitstatus
        end subroutine MakeMagGradGridDH

        subroutine MakeGeoidGrid(geoid, cilm, lmax, r0pot, GM, PotRef, omega, &
                                 r, gridtype, order, nlat, nlong, interval, &
                                 lmax_calc, a, f, extend, exitstatus)
            integer, parameter :: dp = selected_real_kind(p=15)
            real(dp), intent(out) :: geoid(:,:)
            real(dp), intent(in) :: cilm(:,:,:), r0pot, GM, r, PotRef, omega
            integer, intent(in) :: lmax, order, gridtype
            integer, intent(in), optional :: lmax_calc, extend
            integer, intent(out) :: nlat, nlong
            real(dp), intent(in), optional :: interval, a, f
            integer, intent(out), optional :: exitstatus
        end subroutine MakeGeoidGrid

        subroutine CilmPlus(cilm, gridin, lmax, nmax, mass, d, rho, gridtype, &
                            w, zero, plx, n, dref, exitstatus)
            integer, parameter :: dp = selected_real_kind(p=15)
            real(dp), intent(in) :: gridin(:,:), mass, rho
            real(dp), intent(in), optional :: w(:), zero(:), plx(:,:), dref
            real(dp), intent(out) :: cilm(:,:,:), d
            integer, intent(in) :: lmax, nmax, gridtype
            integer, intent(in), optional :: n
            integer, intent(out), optional :: exitstatus
        end subroutine CilmPlus

        subroutine CilmMinus(cilm, gridin, lmax, nmax, mass, d, rho, &
                             gridtype, w, zero, plx, n, dref, exitstatus)
            integer, parameter :: dp = selected_real_kind(p=15)
            real(dp), intent(in) :: gridin(:,:), mass, rho
            real(dp), intent(in), optional :: w(:), zero(:), plx(:,:), dref
            real(dp), intent(out) :: cilm(:,:,:), d
            integer, intent(in) :: lmax, nmax, gridtype
            integer, intent(in), optional :: n
            integer, intent(out), optional :: exitstatus
        end subroutine CilmMinus

        subroutine CilmPlusRhoH(cilm, gridin, lmax, nmax, mass, d, rho, &
                                gridtype, w, zero, plx, n, dref, exitstatus)
            integer, parameter :: dp = selected_real_kind(p=15)
            real(dp), intent(in) :: gridin(:,:), mass, rho(:,:)
            real(dp), intent(in), optional :: w(:), zero(:), plx(:,:), dref
            real(dp), intent(out) :: cilm(:,:,:), d
            integer, intent(in) :: lmax, nmax, gridtype
            integer, intent(in), optional :: n
            integer, intent(out), optional :: exitstatus
        end subroutine CilmPlusRhoH

        subroutine CilmMinusRhoH(cilm, gridin, lmax, nmax, mass, d, rho, &
                                 gridtype, w, zero, plx, n, dref, exitstatus)
            integer, parameter :: dp = selected_real_kind(p=15)
            real(dp), intent(in) :: gridin(:,:), mass, rho(:,:)
            real(dp), intent(in), optional :: w(:), zero(:), plx(:,:), dref
            real(dp), intent(out) :: cilm(:,:,:), d
            integer, intent(in) :: lmax, nmax, gridtype
            integer, intent(in), optional :: n
            integer, intent(out), optional :: exitstatus
        end subroutine CilmMinusRhoH

        subroutine BAtoHilm(cilm, ba, gridglq, lmax, nmax, mass, r0, rho, &
                            gridtype, w, plx, zero, filter_type, filter_deg, &
                            lmax_calc, exitstatus)
            integer, parameter :: dp = selected_real_kind(p=15)
            real(dp), intent(out) :: cilm(:,:,:)
            real(dp), intent(in) :: ba(:,:,:), gridglq(:,:), mass, r0, rho
            real(dp), intent(in), optional :: plx(:,:), zero(:), w(:)
            integer, intent(in) :: lmax, nmax, gridtype
            integer, intent(in), optional :: filter_type, filter_deg, lmax_calc
            integer, intent(out), optional :: exitstatus
        end subroutine BAtoHilm

        subroutine BAtoHilmRhoH(cilm, ba, gridglq, lmax, nmax, mass, r0, &
                                rho, gridtype, w, plx, zero, filter_type, &
                                filter_deg, lmax_calc, exitstatus)
            integer, parameter :: dp = selected_real_kind(p=15)
            real(dp), intent(out) :: cilm(:,:,:)
            real(dp), intent(in) :: ba(:,:,:), gridglq(:,:), mass, r0, rho(:,:)
            real(dp), intent(in), optional :: plx(:,:), zero(:), w(:)
            integer, intent(in) :: lmax, nmax, gridtype
            integer, intent(in), optional :: filter_type, filter_deg, lmax_calc
            integer, intent(out), optional :: exitstatus
        end subroutine BAtoHilmRhoH

        function DownContFilterMA(l, half, r, d)
            integer, parameter :: dp = selected_real_kind(p=15)
            real(dp) :: DownContFilterMA
            integer, intent(in) :: l, half
            real(dp), intent(in) :: r, d
        end function DownContFilterMA

        function DownContFilterMC(l, half, r, d)
            integer, parameter :: dp = selected_real_kind(p=15)
            real(dp) :: DownContFilterMC
            integer, intent(in) :: l, half
            real(dp), intent(in) :: r, d
        end function DownContFilterMC

        function NormalGravity(geocentric_lat, gm, omega, a, b)
            integer, parameter :: dp = selected_real_kind(p=15)
            real(dp) :: NormalGravity
            real(dp), intent(in) :: geocentric_lat, gm, omega, a, b
        end function NormalGravity

        subroutine MakeMagGridDH(cilm, lmax, r0, a, f, rad_grid, theta_grid, &
                                 phi_grid, total_grid, n, sampling, &
                                 lmax_calc, pot_grid, extend, exitstatus)
            integer, parameter :: dp = selected_real_kind(p=15)
            real(dp), intent(in) :: cilm(:,:,:), r0, a, f
            real(dp), intent(out) :: rad_grid(:,:), theta_grid(:,:), &
                                     phi_grid(:,:), total_grid(:,:)
            real(dp), intent(out), optional :: pot_grid(:,:)
            integer, intent(in) :: lmax
            integer, intent(out) :: n
            integer, intent(in), optional :: sampling, lmax_calc, extend
            integer, intent(out), optional :: exitstatus
        end subroutine MakeMagGridDH

        subroutine SHMagPowerSpectrum(cilm, a, r, lmax, spectra, exitstatus)
            integer, parameter :: dp = selected_real_kind(p=15)
            real(dp), intent(in) :: cilm(:,:,:)
            real(dp), intent(in) :: a, r
            integer, intent(in) :: lmax
            real(dp), intent(out) :: spectra(:)
            integer, intent(out), optional :: exitstatus
        end subroutine SHMagPowerSpectrum

        function SHMagPowerL(cilm, a, r, l)
            integer, parameter :: dp = selected_real_kind(p=15)
            real(dp) :: SHMagPowerL
            real(dp), intent(in) :: cilm(:,:,:)
            real(dp), intent(in) :: a, r
            integer, intent(in) :: l
        end function SHMagPowerL

        subroutine MakeCircleCoord(coord, lat, lon, theta0, cinterval, cnum, &
                                   exitstatus)
            integer, parameter :: dp = selected_real_kind(p=15)
            real(dp), intent(in) :: lat, lon, theta0
            real(dp), intent(out) :: coord(:,:)
            real(dp), intent(in), optional :: cinterval
            integer, intent(out), optional :: cnum, exitstatus
        end subroutine MakeCircleCoord

        subroutine MakeEllipseCoord(coord, lat, lon, dec, A_theta, B_theta, &
                                    cinterval, cnum, exitstatus)
            integer, parameter :: dp = selected_real_kind(p=15)
            real(dp), intent(in) :: lat, lon, A_theta, B_theta, dec
            real(dp), intent(out) :: coord(:,:)
            real(dp), intent(in), optional :: cinterval
            integer, intent(out), optional :: cnum, exitstatus
        end subroutine MakeEllipseCoord

        subroutine Wigner3j(w3j, jmin, jmax, j2, j3, m1, m2, m3, exitstatus)
            integer, parameter :: dp = selected_real_kind(p=15)
            integer, intent(in) :: j2, j3, m1, m2, m3
            integer, intent(out) :: jmin, jmax
            real(dp), intent(out) :: w3j(:)
            integer, intent(out), optional :: exitstatus
        end subroutine Wigner3j

        function RandomN(idum)
            integer, parameter :: dp = selected_real_kind(p=15)
            integer, parameter :: int4 = selected_int_kind(9)
            real(dp) :: RandomN
            integer(int4), intent(inout) :: idum
        end function RandomN

        function RandomGaussian(idum)
            integer, parameter :: dp = selected_real_kind(p=15)
            integer, parameter :: int4 = selected_int_kind(9)
            real(dp) :: RandomGaussian
            integer(int4), intent(inout) :: idum
        end function RandomGaussian

        subroutine PreGLQ(x1, x2, n, zero, w, exitstatus)
            integer, parameter :: dp = selected_real_kind(p=15)
            real(dp), intent(in) :: x1, x2
            real(dp), intent(out) :: zero(:), w(:)
            integer, intent(in) :: n
            integer, intent(out), optional :: exitstatus
        end subroutine PreGLQ

        integer function NGLQ(degree)
            integer, intent(in) :: degree
        end function NGLQ

        integer function NGLQSH(degree)
            integer, intent(in) :: degree
        end function NGLQSH

        integer function NGLQSHN(degree, n)
            integer, intent(in) :: degree, n
        end function NGLQSHN

        subroutine DHaj(n, aj, extend, exitstatus)
            integer, parameter :: dp = selected_real_kind(p=15)
            integer, intent(in) :: n
            real(dp), intent(out) :: aj(:)
            integer, intent(in), optional :: extend
            integer, intent(out), optional :: exitstatus
        end subroutine DHaj

        integer function YilmIndexVector(i, l, m)
            integer, intent(in) :: i, l, m
        end function YilmIndexVector

        subroutine EigValVecSym(ain, n, eig, evec, ul, K, exitstatus)
            integer, parameter :: dp = selected_real_kind(p=15)
            real(dp), intent(in) :: ain(:,:)
            integer, intent(in) :: n
            real(dp), intent(out) :: eig(:), evec(:,:)
            character, intent(in), optional :: ul
            integer, intent(in), optional :: K
            integer, intent(out), optional :: exitstatus
        end subroutine EigValVecSym

        subroutine EigValVecSymTri(ain, n, eig, evec, ul, exitstatus)
            integer, parameter :: dp = selected_real_kind(p=15)
            real(dp), intent(in) :: ain(:,:)
            integer, intent(in) :: n
            real(dp), intent(out) :: eig(:), evec(:,:)
            character, intent(in), optional :: ul
            integer, intent(out), optional :: exitstatus
        end subroutine EigValVecSymTri

        subroutine EigValSym(ain, n, eval, ul)
            integer, parameter :: dp = selected_real_kind(p=15)
            real(dp), intent(in) :: ain(:,:)
            integer, intent(in) :: n
            real(dp), intent(out) :: eval(:)
            character, intent(in), optional :: ul
        end subroutine EigValSym

        subroutine SHRotateTapers(tapersrot, tapers, taper_order, lmax, nrot, &
                                  x, dj, exitstatus)
            integer, parameter :: dp = selected_real_kind(p=15)
            real(dp), intent(in) :: tapers(:,:), x(:), dj(:,:,:)
            real(dp), intent(out) :: tapersrot(:,:)
            integer, intent(in) :: taper_order(:), lmax, nrot
            integer, intent(out), optional :: exitstatus
        end subroutine SHRotateTapers

        subroutine SlepianCoeffs(falpha, galpha, film, lmax, nmax, &
                                 exitstatus)
            integer, parameter :: dp = selected_real_kind(p=15)
            real(dp), intent(out) :: falpha(:)
            real(dp), intent(in) :: galpha(:,:), film(:,:,:)
            integer, intent(in) :: lmax, nmax
            integer, intent(out), optional :: exitstatus
        end subroutine SlepianCoeffs

        subroutine SlepianCoeffsToSH(film, falpha, galpha, lmax, nmax, &
                                     exitstatus)
            integer, parameter :: dp = selected_real_kind(p=15)
            real(dp), intent(out) :: film(:,:,:)
            real(dp), intent(in) :: falpha(:), galpha(:,:)
            integer, intent(in) :: lmax, nmax
            integer, intent(out), optional :: exitstatus
        end subroutine SlepianCoeffsToSH

        subroutine SHSCouplingMatrix(kij, galpha, lmax, nmax, exitstatus)
            integer, parameter :: dp = selected_real_kind(p=15)
            real(dp), intent(out) :: kij(:,:)
            real(dp), intent(in) :: galpha(:,:)
            integer, intent(in) :: lmax, nmax
            integer, intent(out), optional :: exitstatus
        end subroutine SHSCouplingMatrix

        subroutine SHSlepianVar(l, galpha, galpha_order, lmax, kmax, Sff, &
                                variance, exitstatus)
            integer, parameter :: dp = selected_real_kind(p=15)
            real(dp), intent(in) :: galpha(:,:), Sff(:)
            real(dp), intent(out) :: variance
            integer, intent(in) :: l, lmax, kmax, galpha_order(:)
            integer, intent(out), optional :: exitstatus
        end subroutine SHSlepianVar

        subroutine SHSCouplingMatrixCap(kij, galpha, galpha_order, lmax, &
                                        nmax, exitstatus)
            integer, parameter :: dp = selected_real_kind(p=15)
            real(dp), intent(out) :: kij(:,:)
            real(dp), intent(in) :: galpha(:,:)
            integer, intent(in) :: galpha_order(:), lmax, nmax
            integer, intent(out), optional :: exitstatus
        end subroutine SHSCouplingMatrixCap

    end interface

end module SHTOOLS
