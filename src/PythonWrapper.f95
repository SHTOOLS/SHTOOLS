    subroutine pyPlmBar(exitstatus,p,lmax,z,csphase,cnorm,p_d0)
        use shtools, only: PlmBar
        use ftypes
        implicit none
        integer(int32),intent(in) :: p_d0
        integer(int32),intent(out) :: exitstatus
        real(dp),dimension(p_d0),intent(out) :: p
        integer(int32),intent(in) :: lmax
        real(dp),intent(in) :: z
        integer(int32),intent(in) :: csphase
        integer(int32),intent(in) :: cnorm
        call PlmBar(p,lmax,z,csphase=csphase,cnorm=cnorm,exitstatus=exitstatus)
    end subroutine pyPlmBar

    subroutine pyPlmBar_d1(exitstatus,p,dp1,lmax,z,csphase,cnorm,p_d0,dp1_d0)
        use shtools, only: PlmBar_d1
        use ftypes
        implicit none
        integer(int32),intent(in) :: p_d0
        integer(int32),intent(in) :: dp1_d0
        integer(int32),intent(out) :: exitstatus
        real(dp),dimension(p_d0),intent(out) :: p
        real(dp),dimension(dp1_d0),intent(out) :: dp1
        integer(int32),intent(in) :: lmax
        real(dp),intent(in) :: z
        integer(int32),intent(in) :: csphase
        integer(int32),intent(in) :: cnorm
        call PlmBar_d1(p,dp1,lmax,z,csphase=csphase,cnorm=cnorm,&
                       exitstatus=exitstatus)
    end subroutine pyPlmBar_d1

    subroutine pyPlBar(exitstatus,p,lmax,z,p_d0)
        use shtools, only: PlBar
        use ftypes
        implicit none
        integer(int32),intent(in) :: p_d0
        integer(int32),intent(out) :: exitstatus
        real(dp),dimension(p_d0),intent(out) :: p
        integer(int32),intent(in) :: lmax
        real(dp),intent(in) :: z
        call PlBar(p,lmax,z,exitstatus=exitstatus)
    end subroutine pyPlBar

    subroutine pyPlBar_d1(exitstatus,p,dp1,lmax,z,p_d0,dp1_d0)
        use shtools, only: PlBar_d1
        use ftypes
        implicit none
        integer(int32),intent(in) :: p_d0
        integer(int32),intent(in) :: dp1_d0
        integer(int32),intent(out) :: exitstatus
        real(dp),dimension(p_d0),intent(out) :: p
        real(dp),dimension(dp1_d0),intent(out) :: dp1
        integer(int32),intent(in) :: lmax
        real(dp),intent(in) :: z
        call PlBar_d1(p,dp1,lmax,z,exitstatus=exitstatus)
    end subroutine pyPlBar_d1

    subroutine pyPlmON(exitstatus,p,lmax,z,csphase,cnorm,p_d0)
        use shtools, only: PlmON
        use ftypes
        implicit none
        integer(int32),intent(in) :: p_d0
        integer(int32),intent(out) :: exitstatus
        real(dp),dimension(p_d0),intent(out) :: p
        integer(int32),intent(in) :: lmax
        real(dp),intent(in) :: z
        integer(int32),intent(in) :: csphase
        integer(int32),intent(in) :: cnorm
        call PlmON(p,lmax,z,csphase=csphase,cnorm=cnorm,exitstatus=exitstatus)
    end subroutine pyPlmON

    subroutine pyPlmON_d1(exitstatus,p,dp1,lmax,z,csphase,cnorm,p_d0,dp1_d0)
        use shtools, only: PlmON_d1
        use ftypes
        implicit none
        integer(int32),intent(in) :: p_d0
        integer(int32),intent(in) :: dp1_d0
        integer(int32),intent(out) :: exitstatus
        real(dp),dimension(p_d0),intent(out) :: p
        real(dp),dimension(dp1_d0),intent(out) :: dp1
        integer(int32),intent(in) :: lmax
        real(dp),intent(in) :: z
        integer(int32),intent(in) :: csphase
        integer(int32),intent(in) :: cnorm
        call PlmON_d1(p,dp1,lmax,z,csphase=csphase,cnorm=cnorm,&
                      exitstatus=exitstatus)
    end subroutine pyPlmON_d1

    subroutine pyPlON(exitstatus,p,lmax,z,p_d0)
        use shtools, only: PlON
        use ftypes
        implicit none
        integer(int32),intent(in) :: p_d0
        integer(int32),intent(out) :: exitstatus
        real(dp),dimension(p_d0),intent(out) :: p
        integer(int32),intent(in) :: lmax
        real(dp),intent(in) :: z
        call PlON(p,lmax,z,exitstatus=exitstatus)
    end subroutine pyPlON

    subroutine pyPlON_d1(exitstatus,p,dp1,lmax,z,p_d0,dp1_d0)
        use shtools, only: PlON_d1
        use ftypes
        implicit none
        integer(int32),intent(in) :: p_d0
        integer(int32),intent(in) :: dp1_d0
        integer(int32),intent(out) :: exitstatus
        real(dp),dimension(p_d0),intent(out) :: p
        real(dp),dimension(dp1_d0),intent(out) :: dp1
        integer(int32),intent(in) :: lmax
        real(dp),intent(in) :: z
        call PlON_d1(p,dp1,lmax,z,exitstatus=exitstatus)
    end subroutine pyPlON_d1

    subroutine pyPlmSchmidt(exitstatus,p,lmax,z,csphase,cnorm,p_d0)
        use shtools, only: PlmSchmidt
        use ftypes
        implicit none
        integer(int32),intent(in) :: p_d0
        integer(int32),intent(out) :: exitstatus
        real(dp),dimension(p_d0),intent(out) :: p
        integer(int32),intent(in) :: lmax
        real(dp),intent(in) :: z
        integer(int32),intent(in) :: csphase
        integer(int32),intent(in) :: cnorm
        call PlmSchmidt(p,lmax,z,csphase=csphase,cnorm=cnorm,&
                        exitstatus=exitstatus)
    end subroutine pyPlmSchmidt

    subroutine pyPlmSchmidt_d1(exitstatus,p,dp1,lmax,z,csphase,cnorm,p_d0,dp1_d0)
        use shtools, only: PlmSchmidt_d1
        use ftypes
        implicit none
        integer(int32),intent(in) :: p_d0
        integer(int32),intent(in) :: dp1_d0
        integer(int32),intent(out) :: exitstatus
        real(dp),dimension(p_d0),intent(out) :: p
        real(dp),dimension(dp1_d0),intent(out) :: dp1
        integer(int32),intent(in) :: lmax
        real(dp),intent(in) :: z
        integer(int32),intent(in) :: csphase
        integer(int32),intent(in) :: cnorm
        call PlmSchmidt_d1(p,dp1,lmax,z,csphase=csphase,cnorm=cnorm,&
                           exitstatus=exitstatus)
    end subroutine pyPlmSchmidt_d1

    subroutine pyPlSchmidt(exitstatus,p,lmax,z,p_d0)
        use shtools, only: PlSchmidt
        use ftypes
        implicit none
        integer(int32),intent(in) :: p_d0
        integer(int32),intent(out) :: exitstatus
        real(dp),dimension(p_d0),intent(out) :: p
        integer(int32),intent(in) :: lmax
        real(dp),intent(in) :: z
        call PlSchmidt(p,lmax,z,exitstatus=exitstatus)
    end subroutine pyPlSchmidt

    subroutine pyPlSchmidt_d1(exitstatus,p,dp1,lmax,z,p_d0,dp1_d0)
        use shtools, only: PlSchmidt_d1
        use ftypes
        implicit none
        integer(int32),intent(in) :: p_d0
        integer(int32),intent(in) :: dp1_d0
        integer(int32),intent(out) :: exitstatus
        real(dp),dimension(p_d0),intent(out) :: p
        real(dp),dimension(dp1_d0),intent(out) :: dp1
        integer(int32),intent(in) :: lmax
        real(dp),intent(in) :: z
        call PlSchmidt_d1(p,dp1,lmax,z,exitstatus=exitstatus)
    end subroutine pyPlSchmidt_d1

    subroutine pyPLegendreA(exitstatus,p,lmax,z,csphase,p_d0)
        use shtools, only: PLegendreA
        use ftypes
        implicit none
        integer(int32),intent(in) :: p_d0
        integer(int32),intent(out) :: exitstatus
        real(dp),dimension(p_d0),intent(out) :: p
        integer(int32),intent(in) :: lmax
        real(dp),intent(in) :: z
        integer(int32),intent(in) :: csphase
        call PLegendreA(p,lmax,z,csphase=csphase,exitstatus=exitstatus)
    end subroutine pyPLegendreA

    subroutine pyPLegendreA_d1(exitstatus,p,dp1,lmax,z,csphase,p_d0,dp1_d0)
        use shtools, only: PLegendreA_d1
        use ftypes
        implicit none
        integer(int32),intent(in) :: p_d0
        integer(int32),intent(in) :: dp1_d0
        integer(int32),intent(out) :: exitstatus
        real(dp),dimension(p_d0),intent(out) :: p
        real(dp),dimension(dp1_d0),intent(out) :: dp1
        integer(int32),intent(in) :: lmax
        real(dp),intent(in) :: z
        integer(int32),intent(in) :: csphase
        call PLegendreA_d1(p,dp1,lmax,z,csphase=csphase,exitstatus=exitstatus)
    end subroutine pyPLegendreA_d1

    subroutine pyPLegendre(exitstatus,p,lmax,z,p_d0)
        use shtools, only: PLegendre
        use ftypes
        implicit none
        integer(int32),intent(in) :: p_d0
        integer(int32),intent(out) :: exitstatus
        real(dp),dimension(p_d0),intent(out) :: p
        integer(int32),intent(in) :: lmax
        real(dp),intent(in) :: z
        call PLegendre(p,lmax,z,exitstatus=exitstatus)
    end subroutine pyPLegendre

    subroutine pyPLegendre_d1(exitstatus,p,dp1,lmax,z,p_d0,dp1_d0)
        use shtools, only: PLegendre_d1
        use ftypes
        implicit none
        integer(int32),intent(in) :: p_d0
        integer(int32),intent(in) :: dp1_d0
        integer(int32),intent(out) :: exitstatus
        real(dp),dimension(p_d0),intent(out) :: p
        real(dp),dimension(dp1_d0),intent(out) :: dp1
        integer(int32),intent(in) :: lmax
        real(dp),intent(in) :: z
        call PLegendre_d1(p,dp1,lmax,z,exitstatus=exitstatus)
    end subroutine pyPLegendre_d1

    subroutine pySHExpandDH(exitstatus,grid,n,cilm,lmax,norm,sampling,csphase,&
                            lmax_calc,cilm_d0,cilm_d1,cilm_d2,grid_d0,grid_d1)
        use shtools, only: SHExpandDH
        use ftypes
        implicit none
        integer(int32),intent(in) :: cilm_d0
        integer(int32),intent(in) :: cilm_d1
        integer(int32),intent(in) :: cilm_d2
        integer(int32),intent(in) :: grid_d0
        integer(int32),intent(in) :: grid_d1
        integer(int32),intent(out) :: exitstatus
        real(dp),dimension(grid_d0,grid_d1),intent(in) :: grid
        integer(int32),intent(in) :: n
        real(dp),dimension(cilm_d0,cilm_d1,cilm_d2),intent(out) :: cilm
        integer(int32),intent(out) :: lmax
        integer(int32),intent(in) :: norm
        integer(int32),intent(in) :: sampling
        integer(int32),intent(in) :: csphase
        integer(int32),intent(in) :: lmax_calc
        call SHExpandDH(grid,n,cilm,lmax,norm=norm,sampling=sampling,&
                        csphase=csphase,lmax_calc=lmax_calc,&
                        exitstatus=exitstatus)
    end subroutine pySHExpandDH

    subroutine pyMakeGridDH(exitstatus,griddh,n,cilm,lmax,norm,sampling,&
                            csphase,lmax_calc,extend,cilm_d0,cilm_d1,cilm_d2,&
                            griddh_d0,griddh_d1)
        use shtools, only: MakeGridDH
        use ftypes
        implicit none
        integer(int32),intent(in) :: cilm_d0
        integer(int32),intent(in) :: cilm_d1
        integer(int32),intent(in) :: cilm_d2
        integer(int32),intent(in) :: griddh_d0
        integer(int32),intent(in) :: griddh_d1
        integer(int32),intent(out) :: exitstatus
        real(dp),dimension(griddh_d0,griddh_d1),intent(out) :: griddh
        integer(int32),intent(out) :: n
        real(dp),dimension(cilm_d0,cilm_d1,cilm_d2),intent(in) :: cilm
        integer(int32),intent(in) :: lmax
        integer(int32),intent(in) :: norm
        integer(int32),intent(in) :: sampling
        integer(int32),intent(in) :: csphase
        integer(int32),intent(in) :: lmax_calc
        integer(int32),intent(in) :: extend
        call MakeGridDH(griddh,n,cilm,lmax,norm=norm,sampling=sampling,&
                        csphase=csphase,lmax_calc=lmax_calc,extend=extend,&
                        exitstatus=exitstatus)
    end subroutine pyMakeGridDH

    subroutine pySHExpandDHC(exitstatus,grid,n,cilm,lmax,norm,sampling,&
                             csphase,lmax_calc,cilm_d0,cilm_d1,cilm_d2,&
                             grid_d0,grid_d1)
        use shtools, only: SHExpandDHC
        use ftypes
        implicit none
        integer(int32),intent(in) :: cilm_d0
        integer(int32),intent(in) :: cilm_d1
        integer(int32),intent(in) :: cilm_d2
        integer(int32),intent(in) :: grid_d0
        integer(int32),intent(in) :: grid_d1
        integer(int32),intent(out) :: exitstatus
        complex(dp),dimension(grid_d0,grid_d1),intent(in) :: grid
        integer(int32),intent(in) :: n
        complex(dp),dimension(cilm_d0,cilm_d1,cilm_d2),intent(out) :: cilm
        integer(int32),intent(out) :: lmax
        integer(int32),intent(in) :: norm
        integer(int32),intent(in) :: sampling
        integer(int32),intent(in) :: csphase
        integer(int32),intent(in) :: lmax_calc
        call SHExpandDHC(grid,n,cilm,lmax,norm=norm,sampling=sampling,&
                         csphase=csphase,lmax_calc=lmax_calc,&
                         exitstatus=exitstatus)
    end subroutine pySHExpandDHC

    subroutine pyMakeGridDHC(exitstatus,griddh,n,cilm,lmax,norm,sampling,&
                             csphase,lmax_calc,extend,cilm_d0,cilm_d1,cilm_d2,&
                             griddh_d0,griddh_d1)
        use shtools, only: MakeGridDHC
        use ftypes
        implicit none
        integer(int32),intent(in) :: cilm_d0
        integer(int32),intent(in) :: cilm_d1
        integer(int32),intent(in) :: cilm_d2
        integer(int32),intent(in) :: griddh_d0
        integer(int32),intent(in) :: griddh_d1
        integer(int32),intent(out) :: exitstatus
        complex(dp),dimension(griddh_d0,griddh_d1),intent(out) :: griddh
        integer(int32),intent(out) :: n
        complex(dp),dimension(cilm_d0,cilm_d1,cilm_d2),intent(in) :: cilm
        integer(int32),intent(in) :: lmax
        integer(int32),intent(in) :: norm
        integer(int32),intent(in) :: sampling
        integer(int32),intent(in) :: csphase
        integer(int32),intent(in) :: lmax_calc
        integer(int32),intent(in) :: extend
        call MakeGridDHC(griddh,n,cilm,lmax,norm=norm,sampling=sampling,&
                         csphase=csphase,lmax_calc=lmax_calc,extend=extend,&
                         exitstatus=exitstatus)
    end subroutine pyMakeGridDHC

    subroutine pyshglq(exitstatus,lmax,zero,w,zero_d0,w_d0)
        use shtools, only: SHGLQ
        use ftypes
        implicit none
        integer(int32),intent(in) :: zero_d0
        integer(int32),intent(in) :: w_d0
        integer(int32),intent(out) :: exitstatus
        integer(int32),intent(in) :: lmax
        real(dp),dimension(zero_d0),intent(out) :: zero
        real(dp),dimension(w_d0),intent(out) :: w
        call SHGLQ(lmax,zero,w,exitstatus=exitstatus)
    end subroutine pySHGLQ

    subroutine pySHExpandGLQ(exitstatus,cilm,lmax,gridglq,w,zero,norm,csphase,&
                             lmax_calc,cilm_d0,cilm_d1,cilm_d2,gridglq_d0,&
                             gridglq_d1,zero_d0,w_d0)
        use shtools, only: SHExpandGLQ
        use ftypes
        implicit none
        integer(int32),intent(in) :: cilm_d0
        integer(int32),intent(in) :: cilm_d1
        integer(int32),intent(in) :: cilm_d2
        integer(int32),intent(in) :: gridglq_d0
        integer(int32),intent(in) :: gridglq_d1
        integer(int32),intent(in) :: zero_d0
        integer(int32),intent(in) :: w_d0
        integer(int32),intent(out) :: exitstatus
        real(dp),dimension(cilm_d0,cilm_d1,cilm_d2),intent(out) :: cilm
        integer(int32),intent(in) :: lmax
        real(dp),dimension(gridglq_d0,gridglq_d1),intent(in) :: gridglq
        real(dp),dimension(w_d0),intent(in) :: w
        real(dp),dimension(zero_d0),intent(in) :: zero
        integer(int32),intent(in) :: norm
        integer(int32),intent(in) :: csphase
        integer(int32),intent(in) :: lmax_calc
        call SHExpandGLQ(cilm,lmax,gridglq,w,zero=zero,norm=norm,&
                         csphase=csphase,lmax_calc=lmax_calc,&
                         exitstatus=exitstatus)
    end subroutine pySHExpandGLQ

    subroutine pyMakeGridGLQ(exitstatus,gridglq,cilm,lmax,zero,norm,csphase,&
                             lmax_calc,extend,gridglq_d0,gridglq_d1,cilm_d0,&
                             cilm_d1,cilm_d2,zero_d0)
        use shtools, only: MakeGridGLQ
        use ftypes
        implicit none
        integer(int32),intent(in) :: gridglq_d0
        integer(int32),intent(in) :: gridglq_d1
        integer(int32),intent(in) :: cilm_d0
        integer(int32),intent(in) :: cilm_d1
        integer(int32),intent(in) :: cilm_d2
        integer(int32),intent(in) :: zero_d0
        integer(int32),intent(out) :: exitstatus
        real(dp),dimension(gridglq_d0,gridglq_d1),intent(out) :: gridglq
        real(dp),dimension(cilm_d0,cilm_d1,cilm_d2),intent(in) :: cilm
        integer(int32),intent(in) :: lmax
        real(dp),dimension(zero_d0),intent(in) :: zero
        integer(int32),intent(in) :: norm
        integer(int32),intent(in) :: csphase
        integer(int32),intent(in) :: lmax_calc
        integer(int32),intent(in) :: extend
        call MakeGridGLQ(gridglq,cilm,lmax,zero=zero,norm=norm,&
                         csphase=csphase,lmax_calc=lmax_calc,extend=extend,&
                         exitstatus=exitstatus)
    end subroutine pyMakeGridGLQ

    subroutine pySHExpandGLQC(exitstatus,cilm,lmax,gridglq,w,zero,norm,&
                              csphase,lmax_calc,cilm_d0,cilm_d1,cilm_d2,&
                              gridglq_d0,gridglq_d1,zero_d0,w_d0)
        use shtools, only: SHExpandGLQC
        use ftypes
        implicit none
        integer(int32),intent(in) :: cilm_d0
        integer(int32),intent(in) :: cilm_d1
        integer(int32),intent(in) :: cilm_d2
        integer(int32),intent(in) :: gridglq_d0
        integer(int32),intent(in) :: gridglq_d1
        integer(int32),intent(in) :: zero_d0
        integer(int32),intent(in) :: w_d0
        integer(int32),intent(out) :: exitstatus
        complex(dp),dimension(cilm_d0,cilm_d1,cilm_d2),intent(out) :: cilm
        integer(int32),intent(in) :: lmax
        complex(dp),dimension(gridglq_d0,gridglq_d1),intent(in) :: gridglq
        real(dp),dimension(w_d0),intent(in) :: w
        real(dp),dimension(zero_d0),intent(in) :: zero
        integer(int32),intent(in) :: norm
        integer(int32),intent(in) :: csphase
        integer(int32),intent(in) :: lmax_calc
        call SHExpandGLQC(cilm,lmax,gridglq,w,zero=zero,norm=norm,&
                          csphase=csphase,lmax_calc=lmax_calc,&
                          exitstatus=exitstatus)
    end subroutine pySHExpandGLQC

    subroutine pyMakeGridGLQC(exitstatus,gridglq,cilm,lmax,zero,norm,csphase,&
                              lmax_calc,extend,gridglq_d0,gridglq_d1,cilm_d0,&
                              cilm_d1,cilm_d2,zero_d0)
        use shtools, only: MakeGridGLQC
        use ftypes
        implicit none
        integer(int32),intent(in) :: gridglq_d0
        integer(int32),intent(in) :: gridglq_d1
        integer(int32),intent(in) :: cilm_d0
        integer(int32),intent(in) :: cilm_d1
        integer(int32),intent(in) :: cilm_d2
        integer(int32),intent(in) :: zero_d0
        integer(int32),intent(out) :: exitstatus
        complex(dp),dimension(gridglq_d0,gridglq_d1),intent(out) :: gridglq
        complex(dp),dimension(cilm_d0,cilm_d1,cilm_d2),intent(in) :: cilm
        integer(int32),intent(in) :: lmax
        real(dp),dimension(zero_d0),intent(in) :: zero
        integer(int32),intent(in) :: norm
        integer(int32),intent(in) :: csphase
        integer(int32),intent(in) :: lmax_calc
        integer(int32),intent(in) :: extend
        call MakeGridGLQC(gridglq,cilm,lmax,zero=zero,norm=norm,&
                          csphase=csphase,lmax_calc=lmax_calc,extend=extend,&
                          exitstatus=exitstatus)
    end subroutine pyMakeGridGLQC

    subroutine pyGLQGridCoord(exitstatus,latglq,longlq,lmax,nlat,nlong,extend,&
                              latglq_d0,longlq_d0)
        use shtools, only: GLQGridCoord
        use ftypes
        implicit none
        integer(int32),intent(in) :: latglq_d0
        integer(int32),intent(in) :: longlq_d0
        integer(int32),intent(out) :: exitstatus
        real(dp),dimension(latglq_d0),intent(out) :: latglq
        real(dp),dimension(longlq_d0),intent(out) :: longlq
        integer(int32),intent(in) :: lmax
        integer(int32),intent(out) :: nlat
        integer(int32),intent(out) :: nlong
        integer(int32),intent(in) :: extend
        call GLQGridCoord(latglq,longlq,lmax,nlat,nlong,extend=extend,&
                          exitstatus=exitstatus)
    end subroutine pyGLQGridCoord

    subroutine pySHExpandLSQ(exitstatus,cilm,d,lat,lon,nmax,lmax,norm,chi2,&
                             csphase,cilm_d0,cilm_d1,cilm_d2,d_d0,lat_d0,&
                             lon_d0)
        use shtools, only: SHExpandLSQ
        use ftypes
        implicit none
        integer(int32),intent(in) :: cilm_d0
        integer(int32),intent(in) :: cilm_d1
        integer(int32),intent(in) :: cilm_d2
        integer(int32),intent(in) :: d_d0
        integer(int32),intent(in) :: lat_d0
        integer(int32),intent(in) :: lon_d0
        integer(int32),intent(out) :: exitstatus
        real(dp),dimension(cilm_d0,cilm_d1,cilm_d2),intent(out) :: cilm
        real(dp),dimension(d_d0),intent(in) :: d
        real(dp),dimension(lat_d0),intent(in) :: lat
        real(dp),dimension(lon_d0),intent(in) :: lon
        integer(int32),intent(in) :: nmax
        integer(int32),intent(in) :: lmax
        integer(int32),intent(in) :: norm
        real(dp),intent(out) :: chi2
        integer(int32),intent(in) :: csphase
        call SHExpandLSQ(cilm,d,lat,lon,nmax,lmax,norm=norm,chi2=chi2,&
                         csphase=csphase,exitstatus=exitstatus)
    end subroutine pySHExpandLSQ

    subroutine pySHExpandWLSQ(exitstatus,cilm,d,w,lat,lon,nmax,lmax,norm,chi2,&
                             csphase,cilm_d0,cilm_d1,cilm_d2,d_d0,w_d0,lat_d0,&
                             lon_d0)
        use shtools, only: SHExpandLSQ
        use ftypes
        implicit none
        integer(int32),intent(in) :: cilm_d0
        integer(int32),intent(in) :: cilm_d1
        integer(int32),intent(in) :: cilm_d2
        integer(int32),intent(in) :: d_d0
        integer(int32),intent(in) :: w_d0
        integer(int32),intent(in) :: lat_d0
        integer(int32),intent(in) :: lon_d0
        integer(int32),intent(out) :: exitstatus
        real(dp),dimension(cilm_d0,cilm_d1,cilm_d2),intent(out) :: cilm
        real(dp),dimension(d_d0),intent(in) :: d
        real(dp),dimension(w_d0),intent(in) :: w
        real(dp),dimension(lat_d0),intent(in) :: lat
        real(dp),dimension(lon_d0),intent(in) :: lon
        integer(int32),intent(in) :: nmax
        integer(int32),intent(in) :: lmax
        integer(int32),intent(in) :: norm
        real(dp),intent(out) :: chi2
        integer(int32),intent(in) :: csphase
        call SHExpandLSQ(cilm,d,lat,lon,nmax,lmax,norm=norm,chi2=chi2,&
                         csphase=csphase,weights=w,exitstatus=exitstatus)
    end subroutine pySHExpandWLSQ

    subroutine pyMakeGrid2d(exitstatus,grid,cilm,lmax,interval,nlat,nlong,&
                            norm,csphase,f,a,north,south,east,west,dealloc,&
                            cilm_d0,cilm_d1,cilm_d2,grid_d0,grid_d1)
        use shtools, only: MakeGrid2d
        use ftypes
        implicit none
        integer(int32),intent(in) :: cilm_d0
        integer(int32),intent(in) :: cilm_d1
        integer(int32),intent(in) :: cilm_d2
        integer(int32),intent(in) :: grid_d0
        integer(int32),intent(in) :: grid_d1
        integer(int32),intent(out) :: exitstatus
        real(dp),dimension(grid_d0,grid_d1),intent(out) :: grid
        real(dp),dimension(cilm_d0,cilm_d1,cilm_d2),intent(in) :: cilm
        integer(int32),intent(in) :: lmax
        real(dp),intent(in) :: interval
        integer(int32),intent(out) :: nlat
        integer(int32),intent(out) :: nlong
        integer(int32),intent(in) :: norm
        integer(int32),intent(in) :: csphase
        real(dp),intent(in) :: f
        real(dp),intent(in) :: a
        real(dp),intent(in) :: north
        real(dp),intent(in) :: south
        real(dp),intent(in) :: east
        real(dp),intent(in) :: west
        integer(int32),intent(in) :: dealloc
        if (f<0.0_dp .and. a<0.0_dp) then
            call MakeGrid2d(grid,cilm,lmax,interval,nlat,nlong,norm=norm,&
                            csphase=csphase,north=north,south=south,east=east,&
                            west=west,dealloc=dealloc,exitstatus=exitstatus)
        else
            call MakeGrid2d(grid,cilm,lmax,interval,nlat,nlong,norm=norm,&
                            csphase=csphase,f=f,a=a,north=north,&
                            south=south,east=east,west=west,dealloc=dealloc,&
                            exitstatus=exitstatus)
        end if
    end subroutine pyMakeGrid2d

    function pyMakeGridPoint(cilm,lmax,lat,lon,norm,csphase,dealloc,&
                             cilm_d0,cilm_d1,cilm_d2)
        use shtools, only: MakeGridPoint
        use ftypes
        implicit none
        integer(int32),intent(in) :: cilm_d0
        integer(int32),intent(in) :: cilm_d1
        integer(int32),intent(in) :: cilm_d2
        real(dp),dimension(cilm_d0,cilm_d1,cilm_d2),intent(in) :: cilm
        integer(int32),intent(in) :: lmax
        real(dp),intent(in) :: lat
        real(dp),intent(in) :: lon
        integer(int32),intent(in) :: norm
        integer(int32),intent(in) :: csphase
        integer(int32),intent(in) :: dealloc
        real(dp) :: pyMakeGridPoint
        pyMakeGridPoint=MakeGridPoint(cilm,lmax,lat,lon,norm=norm,&
                                      csphase=csphase,dealloc=dealloc)
    end function pyMakeGridPoint

    function pyMakeGridPointC(cilm,lmax,lat,lon,norm,csphase,dealloc,&
                              cilm_d0,cilm_d1,cilm_d2)
        use shtools, only: MakeGridPointC
        use ftypes
        implicit none
        integer(int32),intent(in) :: cilm_d0
        integer(int32),intent(in) :: cilm_d1
        integer(int32),intent(in) :: cilm_d2
        complex(dp),dimension(cilm_d0,cilm_d1,cilm_d2),intent(in) :: cilm
        integer(int32),intent(in) :: lmax
        real(dp),intent(in) :: lat
        real(dp),intent(in) :: lon
        integer(int32),intent(in) :: norm
        integer(int32),intent(in) :: csphase
        integer(int32),intent(in) :: dealloc
        complex(dp) :: pyMakeGridPointC
        pyMakeGridPointC=MakeGridPointC(cilm,lmax,lat,lon,norm=norm,&
                                        csphase=csphase,dealloc=dealloc)
    end function pyMakeGridPointC

    subroutine pySHMultiply(exitstatus,shout,sh1,lmax1,sh2,lmax2,precomp,norm,&
                            csphase,sh1_d0,sh1_d1,sh1_d2,sh2_d0,sh2_d1,sh2_d2,&
                            shout_d0,shout_d1,shout_d2)
        use shtools, only: SHMultiply
        use ftypes
        implicit none
        integer(int32),intent(in) :: sh1_d0
        integer(int32),intent(in) :: sh1_d1
        integer(int32),intent(in) :: sh1_d2
        integer(int32),intent(in) :: sh2_d0
        integer(int32),intent(in) :: sh2_d1
        integer(int32),intent(in) :: sh2_d2
        integer(int32),intent(in) :: shout_d0
        integer(int32),intent(in) :: shout_d1
        integer(int32),intent(in) :: shout_d2
        integer(int32),intent(out) :: exitstatus
        real(dp),dimension(shout_d0,shout_d1,shout_d2),intent(out) :: shout
        real(dp),dimension(sh1_d0,sh1_d1,sh1_d2),intent(in) :: sh1
        integer(int32),intent(in) :: lmax1
        real(dp),dimension(sh2_d0,sh2_d1,sh2_d2),intent(in) :: sh2
        integer(int32),intent(in) :: lmax2
        integer(int32),intent(in) :: precomp
        integer(int32),intent(in) :: norm
        integer(int32),intent(in) :: csphase
        call SHMultiply(shout,sh1,lmax1,sh2,lmax2,precomp=precomp,&
                        norm=norm,csphase=csphase,exitstatus=exitstatus)
    end subroutine pySHMultiply

    subroutine pySHRead2(exitstatus,filename,cilm,lmax,lmax_in,gm,r0_pot,dot,&
                         doystart,doyend,epoch,cilm_d0,cilm_d1,cilm_d2,dot_d0,&
                         dot_d1,dot_d2)
        use shtools, only: SHRead2
        use ftypes
        implicit none
        integer(int32),intent(in) :: cilm_d0
        integer(int32),intent(in) :: cilm_d1
        integer(int32),intent(in) :: cilm_d2
        integer(int32),intent(in) :: dot_d0
        integer(int32),intent(in) :: dot_d1
        integer(int32),intent(in) :: dot_d2
        integer(int32),intent(out) :: exitstatus
        character(*),intent(in) :: filename
        real(dp),dimension(cilm_d0,cilm_d1,cilm_d2),intent(out) :: cilm
        integer(int32),intent(out) :: lmax
        integer(int32),intent(in) :: lmax_in
        real(dp),intent(out) :: gm
        real(dp),intent(out) :: r0_pot
        real(dp),dimension(dot_d0,dot_d1,dot_d2),intent(out) :: dot
        real(dp),intent(out) :: doystart
        real(dp),intent(out) :: doyend
        real(dp),intent(out) :: epoch
        call SHRead2(filename,cilm,lmax,gm,r0_pot,dot=dot,doystart=doystart,&
                     doyend=doyend,epoch=epoch,exitstatus=exitstatus)
    end subroutine pySHRead2

    subroutine pySHRead2Error(exitstatus,filename,cilm,error,lmax,lmax_in,gm,&
                              r0_pot,dot,doystart,doyend,epoch,cilm_d0,&
                              cilm_d1,cilm_d2,error_d0,error_d1,error_d2,&
                              dot_d0,dot_d1,dot_d2)
        use shtools, only: SHRead2
        use ftypes
        implicit none
        integer(int32),intent(in) :: cilm_d0
        integer(int32),intent(in) :: cilm_d1
        integer(int32),intent(in) :: cilm_d2
        integer(int32),intent(in) :: error_d0
        integer(int32),intent(in) :: error_d1
        integer(int32),intent(in) :: error_d2
        integer(int32),intent(in) :: dot_d0
        integer(int32),intent(in) :: dot_d1
        integer(int32),intent(in) :: dot_d2
        integer(int32),intent(out) :: exitstatus
        character(*),intent(in) :: filename
        real(dp),dimension(cilm_d0,cilm_d1,cilm_d2),intent(out) :: cilm
        real(dp),dimension(error_d0,error_d1,error_d2),intent(out) :: error
        integer(int32),intent(out) :: lmax
        integer(int32),intent(in) :: lmax_in
        real(dp),intent(out) :: gm
        real(dp),intent(out) :: r0_pot
        real(dp),dimension(dot_d0,dot_d1,dot_d2),intent(out) :: dot
        real(dp),intent(out) :: doystart
        real(dp),intent(out) :: doyend
        real(dp),intent(out) :: epoch
        call SHRead2(filename,cilm,lmax,gm,r0_pot,error=error,dot=dot,&
                     doystart=doystart,doyend=doyend,epoch=epoch,&
                     exitstatus=exitstatus)
    end subroutine pySHRead2Error

    subroutine pySHReadJPL(exitstatus,filename,cilm,lmax,lmax_in,gm,&
                           formatstring,cilm_d0,cilm_d1,cilm_d2)
        use shtools, only: SHReadJPL
        use ftypes
        implicit none
        integer(int32),intent(in) :: cilm_d0
        integer(int32),intent(in) :: cilm_d1
        integer(int32),intent(in) :: cilm_d2
        integer(int32),intent(out) :: exitstatus
        character(*),intent(in) :: filename
        real(dp),dimension(cilm_d0,cilm_d1,cilm_d2),intent(out) :: cilm
        integer(int32),intent(out) :: lmax
        integer(int32),intent(in) :: lmax_in
        real(dp),dimension(2),intent(out) :: gm
        character(6),intent(in) :: formatstring
        call SHReadJPL(filename,cilm,lmax,gm=gm,formatstring=formatstring,&
                       exitstatus=exitstatus)
    end subroutine pySHReadJPL

    subroutine pySHReadJPLError(exitstatus,filename,cilm,error,lmax,lmax_in,&
                                gm,formatstring,cilm_d0,cilm_d1,cilm_d2,&
                                error_d0,error_d1,error_d2)
        use shtools, only: SHReadJPL
        use ftypes
        implicit none
        integer(int32),intent(in) :: cilm_d0
        integer(int32),intent(in) :: cilm_d1
        integer(int32),intent(in) :: cilm_d2
        integer(int32),intent(in) :: error_d0
        integer(int32),intent(in) :: error_d1
        integer(int32),intent(in) :: error_d2
        integer(int32),intent(out) :: exitstatus
        character(*),intent(in) :: filename
        real(dp),dimension(cilm_d0,cilm_d1,cilm_d2),intent(out) :: cilm
        real(dp),dimension(error_d0,error_d1,error_d2),intent(out) :: error
        integer(int32),intent(in) :: lmax
        integer(int32),intent(in) :: lmax_in
        real(dp),dimension(2),intent(out) :: gm
        character(6),intent(in) :: formatstring
        call SHReadJPL(filename,cilm,lmax,error=error,gm=gm,&
                       formatstring=formatstring,exitstatus=exitstatus)
    end subroutine pySHReadJPLError

    subroutine pySHCilmToVector(exitstatus,cilm,vector,lmax,vector_d0,cilm_d0,&
                                cilm_d1,cilm_d2)
        use shtools, only: SHCilmToVector
        use ftypes
        implicit none
        integer(int32),intent(in) :: vector_d0
        integer(int32),intent(in) :: cilm_d0
        integer(int32),intent(in) :: cilm_d1
        integer(int32),intent(in) :: cilm_d2
        integer(int32),intent(out) :: exitstatus
        real(dp),dimension(cilm_d0,cilm_d1,cilm_d2),intent(in) :: cilm
        real(dp),dimension(vector_d0),intent(out) :: vector
        integer(int32),intent(in) :: lmax
        call SHCilmToVector(cilm,vector,lmax,exitstatus=exitstatus)
    end subroutine pySHCilmToVector

    subroutine pySHVectorToCilm(exitstatus,vector,cilm,lmax,vector_d0,cilm_d0,&
                                cilm_d1,cilm_d2)
        use shtools, only: SHVectorToCilm
        use ftypes
        implicit none
        integer(int32),intent(in) :: vector_d0
        integer(int32),intent(in) :: cilm_d0
        integer(int32),intent(in) :: cilm_d1
        integer(int32),intent(in) :: cilm_d2
        integer(int32),intent(out) :: exitstatus
        real(dp),dimension(vector_d0),intent(in) :: vector
        real(dp),dimension(cilm_d0,cilm_d1,cilm_d2),intent(out) :: cilm
        integer(int32),intent(in) :: lmax
        call SHVectorToCilm(vector,cilm,lmax,exitstatus=exitstatus)
    end subroutine pySHVectorToCilm

    subroutine pySHCilmToCindex(exitstatus,cilm,cindex,degmax,cindex_d0,&
                                cindex_d1,cilm_d0,cilm_d1,cilm_d2)
        use shtools, only: SHCilmToCindex
        use ftypes
        implicit none
        integer(int32),intent(in) :: cindex_d0
        integer(int32),intent(in) :: cindex_d1
        integer(int32),intent(in) :: cilm_d0
        integer(int32),intent(in) :: cilm_d1
        integer(int32),intent(in) :: cilm_d2
        integer(int32),intent(out) :: exitstatus
        real(dp),dimension(cilm_d0,cilm_d1,cilm_d2),intent(in) :: cilm
        real(dp),dimension(cindex_d0,cindex_d1),intent(out) :: cindex
        integer(int32),intent(in) :: degmax
        call SHCilmToCindex(cilm,cindex,degmax=degmax,exitstatus=exitstatus)
    end subroutine pySHCilmToCindex

    subroutine pySHCindexToCilm(exitstatus,cindex,cilm,degmax,cindex_d0,&
                                cindex_d1,cilm_d0,cilm_d1,cilm_d2)
        use shtools, only: SHCindexToCilm
        use ftypes
        implicit none
        integer(int32),intent(in) :: cindex_d0
        integer(int32),intent(in) :: cindex_d1
        integer(int32),intent(in) :: cilm_d0
        integer(int32),intent(in) :: cilm_d1
        integer(int32),intent(in) :: cilm_d2
        integer(int32),intent(out) :: exitstatus
        real(dp),dimension(cindex_d0,cindex_d1),intent(in) :: cindex
        real(dp),dimension(cilm_d0,cilm_d1,cilm_d2),intent(out) :: cilm
        integer(int32),intent(in) :: degmax
        call SHCindexToCilm(cindex,cilm,degmax=degmax,exitstatus=exitstatus)
    end subroutine pySHCindexToCilm

    subroutine pySHrtoc(exitstatus,rcilm,ccilm,degmax,convention,switchcs,&
                        rcilm_d0,rcilm_d1,rcilm_d2,ccilm_d0,ccilm_d1,ccilm_d2)
        use shtools, only: SHrtoc
        use ftypes
        implicit none
        integer(int32),intent(in) :: rcilm_d0
        integer(int32),intent(in) :: rcilm_d1
        integer(int32),intent(in) :: rcilm_d2
        integer(int32),intent(in) :: ccilm_d0
        integer(int32),intent(in) :: ccilm_d1
        integer(int32),intent(in) :: ccilm_d2
        integer(int32),intent(out) :: exitstatus
        real(dp),dimension(rcilm_d0,rcilm_d1,rcilm_d2),intent(in) :: rcilm
        real(dp),dimension(ccilm_d0,ccilm_d1,ccilm_d2),intent(out) :: ccilm
        integer(int32),intent(in) :: degmax
        integer(int32),intent(in) :: convention
        integer(int32),intent(in) :: switchcs
        call SHrtoc(rcilm,ccilm,degmax=degmax,convention=convention,&
                    switchcs=switchcs,exitstatus=exitstatus)
    end subroutine pySHrtoc

    subroutine pySHctor(exitstatus,ccilm,rcilm,degmax,convention,switchcs,&
                        rcilm_d0,rcilm_d1,rcilm_d2,ccilm_d0,ccilm_d1,ccilm_d2)
        use shtools, only: SHctor
        use ftypes
        implicit none
        integer(int32),intent(in) :: rcilm_d0
        integer(int32),intent(in) :: rcilm_d1
        integer(int32),intent(in) :: rcilm_d2
        integer(int32),intent(in) :: ccilm_d0
        integer(int32),intent(in) :: ccilm_d1
        integer(int32),intent(in) :: ccilm_d2
        integer(int32),intent(out) :: exitstatus
        real(dp),dimension(ccilm_d0,ccilm_d1,ccilm_d2),intent(in) :: ccilm
        real(dp),dimension(rcilm_d0,rcilm_d1,rcilm_d2),intent(out) :: rcilm
        integer(int32),intent(in) :: degmax
        integer(int32),intent(in) :: convention
        integer(int32),intent(in) :: switchcs
        call SHctor(ccilm,rcilm,degmax=degmax,convention=convention,&
                    switchcs=switchcs,exitstatus=exitstatus)
    end subroutine pySHctor

    subroutine pydjpi2(exitstatus,dj,lmax,dj_d0,dj_d1,dj_d2)
        use shtools, only: djpi2
        use ftypes
        implicit none
        integer(int32),intent(in) :: dj_d0
        integer(int32),intent(in) :: dj_d1
        integer(int32),intent(in) :: dj_d2
        integer(int32),intent(out) :: exitstatus
        real(dp),dimension(dj_d0,dj_d1,dj_d2),intent(out) :: dj
        integer(int32),intent(in) :: lmax
        call djpi2(dj,lmax,exitstatus=exitstatus)
    end subroutine pydjpi2

    subroutine pySHRotateCoef(exitstatus,x,cof,rcof,dj,lmax,rcof_d0,rcof_d1,&
                              dj_d0,dj_d1,dj_d2,cof_d0,cof_d1)
        use shtools, only: SHRotateCoef
        use ftypes
        implicit none
        integer(int32),intent(in) :: rcof_d0
        integer(int32),intent(in) :: rcof_d1
        integer(int32),intent(in) :: dj_d0
        integer(int32),intent(in) :: dj_d1
        integer(int32),intent(in) :: dj_d2
        integer(int32),intent(in) :: cof_d0
        integer(int32),intent(in) :: cof_d1
        integer(int32),intent(out) :: exitstatus
        real(dp),dimension(3),intent(in) :: x
        real(dp),dimension(cof_d0,cof_d1),intent(in) :: cof
        real(dp),dimension(rcof_d0,rcof_d1),intent(out) :: rcof
        real(dp),dimension(dj_d0,dj_d1,dj_d2),intent(in) :: dj
        integer(int32),intent(in) :: lmax
        call SHRotateCoef(x,cof,rcof,dj,lmax,exitstatus=exitstatus)
    end subroutine pySHRotateCoef

    subroutine pySHRotateRealCoef(exitstatus,cilmrot,cilm,lmax,x,dj,x_d0,&
                                  dj_d0,dj_d1,dj_d2,cilm_d0,cilm_d1,cilm_d2,&
                                  cilmrot_d0,cilmrot_d1,cilmrot_d2)
        use shtools, only: SHRotateRealCoef
        use ftypes
        implicit none
        integer(int32),intent(in) :: x_d0
        integer(int32),intent(in) :: dj_d0
        integer(int32),intent(in) :: dj_d1
        integer(int32),intent(in) :: dj_d2
        integer(int32),intent(in) :: cilm_d0
        integer(int32),intent(in) :: cilm_d1
        integer(int32),intent(in) :: cilm_d2
        integer(int32),intent(in) :: cilmrot_d0
        integer(int32),intent(in) :: cilmrot_d1
        integer(int32),intent(in) :: cilmrot_d2
        integer(int32),intent(out) :: exitstatus
        real(dp),dimension(cilmrot_d0,cilmrot_d1,cilmrot_d2),intent(out) ::&
                                                                        cilmrot
        real(dp),dimension(cilm_d0,cilm_d1,cilm_d2),intent(in) :: cilm
        integer(int32),intent(in) :: lmax
        real(dp),dimension(x_d0),intent(in) :: x
        real(dp),dimension(dj_d0,dj_d1,dj_d2),intent(in) :: dj
        call SHRotateRealCoef(cilmrot,cilm,lmax,x,dj,exitstatus=exitstatus)
    end subroutine pySHRotateRealCoef

    subroutine pySHAdmitCorr(exitstatus,G,T,lmax,admit,admit_error,corr,G_d0,&
                             G_d1,G_d2,admit_d0,admit_error_d0,T_d0,T_d1,T_d2,&
                             corr_d0)
        use shtools, only: SHAdmitCorr
        use ftypes
        implicit none
        integer(int32),intent(in) :: G_d0
        integer(int32),intent(in) :: G_d1
        integer(int32),intent(in) :: G_d2
        integer(int32),intent(in) :: admit_d0
        integer(int32),intent(in) :: admit_error_d0
        integer(int32),intent(in) :: T_d0
        integer(int32),intent(in) :: T_d1
        integer(int32),intent(in) :: T_d2
        integer(int32),intent(in) :: corr_d0
        integer(int32),intent(out) :: exitstatus
        real(dp),dimension(G_d0,G_d1,G_d2),intent(in) :: G
        real(dp),dimension(T_d0,T_d1,T_d2),intent(in) :: T
        integer(int32),intent(in) :: lmax
        real(dp),dimension(admit_d0),intent(out) :: admit
        real(dp),dimension(admit_error_d0),intent(out) :: admit_error
        real(dp),dimension(corr_d0),intent(out) :: corr
        call SHAdmitCorr(G,T,lmax,admit,corr,admit_error=admit_error,&
                         exitstatus=exitstatus)
    end subroutine pySHAdmitCorr

    function pySHConfidence(l_conf,r)
        use shtools, only: SHConfidence
        use ftypes
        implicit none
        integer(int32),intent(in) :: l_conf
        real(dp),intent(in) :: r
        real(dp) :: pySHConfidence
        pySHConfidence=SHConfidence(l_conf,r)
    end function pySHConfidence

    subroutine pySHMultiTaperSE(exitstatus,mtse,sd,sh,lmax,tapers,taper_order,&
                                lmaxt,k,lat,lon,taper_wt,norm,csphase,&
                                taper_order_d0,taper_wt_d0,sh_d0,sh_d1,sh_d2,&
                                tapers_d0,tapers_d1,mtse_d0,sd_d0)
        use shtools, only: SHMultiTaperSE
        use ftypes
        implicit none
        integer(int32),intent(in) :: taper_order_d0
        integer(int32),intent(in) :: taper_wt_d0
        integer(int32),intent(in) :: sh_d0
        integer(int32),intent(in) :: sh_d1
        integer(int32),intent(in) :: sh_d2
        integer(int32),intent(in) :: tapers_d0
        integer(int32),intent(in) :: tapers_d1
        integer(int32),intent(in) :: mtse_d0
        integer(int32),intent(in) :: sd_d0
        integer(int32),intent(out) :: exitstatus
        real(dp),dimension(mtse_d0),intent(out) :: mtse
        real(dp),dimension(sd_d0),intent(out) :: sd
        real(dp),dimension(sh_d0,sh_d1,sh_d2),intent(in) :: sh
        integer(int32),intent(in) :: lmax
        real(dp),dimension(tapers_d0,tapers_d1),intent(in) :: tapers
        integer(int32),dimension(taper_order_d0),intent(in) :: taper_order
        integer(int32),intent(in) :: lmaxt
        integer(int32),intent(in) :: k
        real(dp),intent(in) :: lat
        real(dp),intent(in) :: lon
        real(dp),dimension(taper_wt_d0),intent(in) :: taper_wt
        integer(int32),intent(in) :: norm
        integer(int32),intent(in) :: csphase
        if (taper_wt(1) < 0.0_dp) then
            call SHMultiTaperSE(mtse,sd,sh,lmax,tapers,taper_order,&
                                lmaxt,k,lat=lat,lon=lon,norm=norm,&
                                csphase=csphase,exitstatus=exitstatus)
        else
            call SHMultiTaperSE(mtse,sd,sh,lmax,tapers,taper_order,lmaxt,&
                                k,lat=lat,lon=lon,taper_wt=taper_wt,norm=norm,&
                                csphase=csphase,exitstatus=exitstatus)
        endif
    end subroutine pySHMultiTaperSE

    subroutine pySHMultiTaperCSE(exitstatus,mtse,sd,sh1,lmax1,sh2,lmax2,&
                                 tapers,taper_order,lmaxt,k,lat,lon,taper_wt,&
                                 norm,csphase,sh1_d0,sh1_d1,sh1_d2,sh2_d0,&
                                 sh2_d1,sh2_d2,taper_order_d0,taper_wt_d0,&
                                 tapers_d0,tapers_d1,sd_d0,mtse_d0)
        use shtools, only: SHMultiTaperCSE
        use ftypes
        implicit none
        integer(int32),intent(in) :: sh1_d0
        integer(int32),intent(in) :: sh1_d1
        integer(int32),intent(in) :: sh1_d2
        integer(int32),intent(in) :: sh2_d0
        integer(int32),intent(in) :: sh2_d1
        integer(int32),intent(in) :: sh2_d2
        integer(int32),intent(in) :: taper_order_d0
        integer(int32),intent(in) :: taper_wt_d0
        integer(int32),intent(in) :: tapers_d0
        integer(int32),intent(in) :: tapers_d1
        integer(int32),intent(in) :: sd_d0
        integer(int32),intent(in) :: mtse_d0
        integer(int32),intent(out) :: exitstatus
        real(dp),dimension(mtse_d0),intent(out) :: mtse
        real(dp),dimension(sd_d0),intent(out) :: sd
        real(dp),dimension(sh1_d0,sh1_d1,sh1_d2),intent(in) :: sh1
        integer(int32),intent(in) :: lmax1
        real(dp),dimension(sh2_d0,sh2_d1,sh2_d2),intent(in) :: sh2
        integer(int32),intent(in) :: lmax2
        real(dp),dimension(tapers_d0,tapers_d1),intent(in) :: tapers
        integer(int32),dimension(taper_order_d0),intent(in) :: taper_order
        integer(int32),intent(in) :: lmaxt
        integer(int32),intent(in) :: k
        real(dp),intent(in) :: lat
        real(dp),intent(in) :: lon
        real(dp),dimension(taper_wt_d0),intent(in) :: taper_wt
        integer(int32),intent(in) :: norm
        integer(int32),intent(in) :: csphase
        if(taper_wt(1) < 0.0_dp) then
            call SHMultiTaperCSE(mtse,sd,sh1,lmax1,sh2,lmax2,tapers,&
                                 taper_order,lmaxt,k,lat=lat,lon=lon,&
                                 norm=norm,csphase=csphase,&
                                 exitstatus=exitstatus)
        else
            call SHMultiTaperCSE(mtse,sd,sh1,lmax1,sh2,lmax2,tapers,&
                                 taper_order,lmaxt,k,lat=lat,lon=lon,&
                                 taper_wt=taper_wt,norm=norm,csphase=csphase,&
                                 exitstatus=exitstatus)
        end if
    end subroutine pySHMultiTaperCSE

    subroutine pySHLocalizedAdmitCorr(exitstatus,g,t,tapers,taper_order,k,lat,&
                                      lon,lwin,lmax,admit,corr,admit_error,&
                                      corr_error,taper_wt,mtdef,k1linsig,&
                                      taper_order_d0,g_d0,g_d1,g_d2,&
                                      taper_wt_d0,corr_error_d0,admit_d0,&
                                      admit_error_d0,corr_d0,tapers_d0,&
                                      tapers_d1,t_d0,t_d1,t_d2)
        use shtools, only: SHLocalizedAdmitCorr
        use ftypes
        implicit none
        integer(int32),intent(in) :: taper_order_d0
        integer(int32),intent(in) :: g_d0
        integer(int32),intent(in) :: g_d1
        integer(int32),intent(in) :: g_d2
        integer(int32),intent(in) :: taper_wt_d0
        integer(int32),intent(in) :: corr_error_d0
        integer(int32),intent(in) :: admit_d0
        integer(int32),intent(in) :: admit_error_d0
        integer(int32),intent(in) :: corr_d0
        integer(int32),intent(in) :: tapers_d0
        integer(int32),intent(in) :: tapers_d1
        integer(int32),intent(in) :: t_d0
        integer(int32),intent(in) :: t_d1
        integer(int32),intent(in) :: t_d2
        integer(int32),intent(out) :: exitstatus
        real(dp),dimension(g_d0,g_d1,g_d2),intent(in) :: g
        real(dp),dimension(t_d0,t_d1,t_d2),intent(in) :: t
        real(dp),dimension(tapers_d0,tapers_d1),intent(in) :: tapers
        integer(int32),dimension(taper_order_d0),intent(in) :: taper_order
        real(dp),intent(in) :: lat
        real(dp),intent(in) :: lon
        integer(int32),intent(in) :: lwin
        integer(int32),intent(in) :: lmax
        real(dp),dimension(admit_d0),intent(out) :: admit
        real(dp),dimension(corr_d0),intent(out) :: corr
        integer(int32),intent(in) :: k
        real(dp),dimension(admit_error_d0),intent(out) :: admit_error
        real(dp),dimension(corr_error_d0),intent(out) :: corr_error
        real(dp),dimension(taper_wt_d0),intent(in) :: taper_wt
        integer(int32),intent(in) :: mtdef
        integer(int32),intent(in) :: k1linsig
        if(taper_wt(1) < 0.0_dp) then
            if (k1linsig<0) then
                call SHLocalizedAdmitCorr(tapers,taper_order,lwin,lat,lon,g,t,&
                                          lmax,admit,corr,k,&
                                          admit_error=admit_error,&
                                          corr_error=corr_error,&
                                          mtdef=mtdef,exitstatus=exitstatus)
            else
                call SHLocalizedAdmitCorr(tapers,taper_order,lwin,lat,lon,g,t,&
                                          lmax,admit,corr,k,&
                                          admit_error=admit_error,&
                                          corr_error=corr_error,&
                                          mtdef=mtdef,k1linsig=k1linsig,&
                                          exitstatus=exitstatus)
            end if
        else
            if (k1linsig<0) then
                call SHLocalizedAdmitCorr(tapers,taper_order,lwin,lat,lon,g,t,&
                                          lmax,admit,corr,k,&
                                          admit_error=admit_error,&
                                          corr_error=corr_error,&
                                          taper_wt=taper_wt,mtdef=mtdef,&
                                          exitstatus=exitstatus)
            else
                call SHLocalizedAdmitCorr(tapers,taper_order,lwin,lat,lon,g,t,&
                                          lmax,admit,corr,k,&
                                          admit_error=admit_error,&
                                          corr_error=corr_error,&
                                          taper_wt=taper_wt,mtdef=mtdef,&
                                          k1linsig=k1linsig,&
                                          exitstatus=exitstatus)
            end if
        end if
    end subroutine pySHLocalizedAdmitCorr

    subroutine pySHReturnTapers(exitstatus,theta0,lmax,tapers,eigenvalues,&
                                taper_order,degrees,eigenvalues_d0,tapers_d0,&
                                tapers_d1,taper_order_d0,degrees_d0)
        use shtools, only: SHReturnTapers
        use ftypes
        implicit none
        integer(int32),intent(in) :: eigenvalues_d0
        integer(int32),intent(in) :: tapers_d0
        integer(int32),intent(in) :: tapers_d1
        integer(int32),intent(in) :: taper_order_d0
        integer(int32),intent(in) :: degrees_d0
        integer(int32),intent(out) :: exitstatus
        real(dp),intent(in) :: theta0
        integer(int32),intent(in) :: lmax
        real(dp),dimension(tapers_d0,tapers_d1),intent(out) :: tapers
        real(dp),dimension(eigenvalues_d0),intent(out) :: eigenvalues
        integer(int32),dimension(taper_order_d0),intent(out) :: taper_order
        integer(int32),dimension(degrees_d0),intent(in) :: degrees
        call SHReturnTapers(theta0,lmax,tapers,eigenvalues,taper_order,&
                            degrees=degrees,exitstatus=exitstatus)
    end subroutine pySHReturnTapers

    subroutine pySHReturnTapersM(exitstatus,theta0,lmax,m,tapers,eigenvalues,&
                                 degrees,tapers_d0,tapers_d1,eigenvalues_d0,&
                                 degrees_d0)
        use shtools, only: SHReturnTapersM
        use ftypes
        implicit none
        integer(int32),intent(in) :: tapers_d0
        integer(int32),intent(in) :: tapers_d1
        integer(int32),intent(in) :: eigenvalues_d0
        integer(int32),intent(in) :: degrees_d0
        integer(int32),intent(out) :: exitstatus
        real(dp),intent(in) :: theta0
        integer(int32),intent(in) :: lmax
        integer(int32),intent(in) :: m
        real(dp),dimension(tapers_d0,tapers_d1),intent(out) :: tapers
        real(dp),dimension(eigenvalues_d0),intent(out) :: eigenvalues
        integer(int32),dimension(degrees_d0),intent(in) :: degrees
        call SHReturnTapersM(theta0,lmax,m,tapers,eigenvalues,&
                             degrees=degrees,exitstatus=exitstatus)
    end subroutine pySHReturnTapersM

    subroutine pyComputeDm(exitstatus,dllm,lmax,m,theta0,degrees,dllm_d0,dllm_d1,degrees_d0)
        use shtools, only: ComputeDm
        use ftypes
        implicit none
        integer(int32),intent(in) :: dllm_d0
        integer(int32),intent(in) :: dllm_d1
        integer(int32),intent(in) :: degrees_d0
        integer(int32),intent(out) :: exitstatus
        real(dp),dimension(dllm_d0,dllm_d1),intent(out) :: dllm
        integer(int32),intent(in) :: lmax
        integer(int32),intent(in) :: m
        real(dp),intent(in) :: theta0
        integer(int32),dimension(degrees_d0),intent(in) :: degrees
        call ComputeDm(dllm,lmax,m,theta0,degrees=degrees,exitstatus=exitstatus)
    end subroutine pyComputeDm

    subroutine pyComputeDG82(exitstatus,dG82,lmax,m,theta0,dG82_d0,dG82_d1)
        use shtools, only: ComputeDG82
        use ftypes
        implicit none
        integer(int32),intent(in) :: dG82_d0
        integer(int32),intent(in) :: dG82_d1
        integer(int32),intent(out) :: exitstatus
        real(dp),dimension(dG82_d0,dG82_d1),intent(out) :: dG82
        integer(int32),intent(in) :: lmax
        integer(int32),intent(in) :: m
        real(dp),intent(in) :: theta0
        call ComputeDG82(dG82,lmax,m,theta0,exitstatus=exitstatus)
    end subroutine pyComputeDG82

    function pySHFindLWin(theta0,m,alpha,taper_number)
        use shtools, only: SHFindLWin
        use ftypes
        implicit none
        real(dp),intent(in) :: theta0
        integer(int32),intent(in) :: m
        real(dp),intent(in) :: alpha
        integer(int32),intent(in) :: taper_number
        integer(int32) :: pySHFindLWin
        pySHFindLWin=SHFindLWin(theta0,m,alpha,taper_number=taper_number)
    end function pySHFindLWin

    subroutine pySHBiasK(exitstatus,tapers,lwin,k,incspectra,ldata,&
                         outcspectra,taper_wt,save_cg,taper_wt_d0,tapers_d0,&
                         tapers_d1,incspectra_d0,outcspectra_d0)
        use shtools, only: SHBiasK
        use ftypes
        implicit none
        integer(int32),intent(in) :: taper_wt_d0
        integer(int32),intent(in) :: tapers_d0
        integer(int32),intent(in) :: tapers_d1
        integer(int32),intent(in) :: incspectra_d0
        integer(int32),intent(in) :: outcspectra_d0
        integer(int32),intent(out) :: exitstatus
        real(dp),dimension(tapers_d0,tapers_d1),intent(in) :: tapers
        integer(int32),intent(in) :: lwin
        integer(int32),intent(in) :: k
        real(dp),dimension(incspectra_d0),intent(in) :: incspectra
        integer(int32),intent(in) :: ldata
        real(dp),dimension(outcspectra_d0),intent(out) :: outcspectra
        real(dp),dimension(taper_wt_d0),intent(in) :: taper_wt
        integer(int32),intent(in) :: save_cg
        if (taper_wt(1) < 0.0_dp) then
            call SHBiasK(tapers,lwin,k,incspectra,ldata,outcspectra,&
                         save_cg=save_cg,exitstatus=exitstatus)
        else
            call SHBiasK(tapers,lwin,k,incspectra,ldata,outcspectra,&
                         taper_wt=taper_wt,save_cg=save_cg,&
                         exitstatus=exitstatus)
        end if
    end subroutine pySHBiasK

    subroutine pySHMTCouplingMatrix(exitstatus,Mmt,lmax,tapers_power,lwin,&
                                    k,taper_wt,Mmt_d0,Mmt_d1,&
                                    tapers_power_d0,tapers_power_d1,&
                                    taper_wt_d0)
        use shtools, only: SHMTCouplingMatrix
        use ftypes
        implicit none
        integer(int32),intent(in) :: Mmt_d0
        integer(int32),intent(in) :: Mmt_d1
        integer(int32),intent(in) :: tapers_power_d0
        integer(int32),intent(in) :: tapers_power_d1
        integer(int32),intent(in) :: taper_wt_d0
        integer(int32),intent(out) :: exitstatus
        integer(int32),intent(in) :: lmax,k,lwin
        real(dp),intent(out) :: Mmt(Mmt_d0,Mmt_d1)
        real(dp),intent(in) :: tapers_power(tapers_power_d0,tapers_power_d1)
        real(dp),dimension(taper_wt_d0),intent(in) :: taper_wt
        if (taper_wt(1) < 0.0_dp) then
            call SHMTCouplingMatrix(Mmt,lmax,tapers_power,lwin,k,&
                                    exitstatus=exitstatus)
        else
            call SHMTCouplingMatrix(Mmt,lmax,tapers_power,lwin,k,&
                                    taper_wt=taper_wt,&
                                    exitstatus=exitstatus)
        end if
    end subroutine

    subroutine pySHBiasAdmitCorr(exitstatus,sgt,sgg,stt,lmax,tapers,lwin,k,&
                                 admit,corr,mtdef,taper_wt,taper_wt_d0,sgt_d0,&
                                 stt_d0,admit_d0,tapers_d0,tapers_d1,corr_d0,&
                                 sgg_d0)
        use shtools, only: SHBiasAdmitCorr
        use ftypes
        implicit none
        integer(int32),intent(in) :: taper_wt_d0
        integer(int32),intent(in) :: sgt_d0
        integer(int32),intent(in) :: stt_d0
        integer(int32),intent(in) :: admit_d0
        integer(int32),intent(in) :: tapers_d0
        integer(int32),intent(in) :: tapers_d1
        integer(int32),intent(in) :: corr_d0
        integer(int32),intent(in) :: sgg_d0
        integer(int32),intent(out) :: exitstatus
        real(dp),dimension(sgt_d0),intent(in) :: sgt
        real(dp),dimension(sgg_d0),intent(in) :: sgg
        real(dp),dimension(stt_d0),intent(in) :: stt
        integer(int32),intent(in) :: lmax
        real(dp),dimension(tapers_d0,tapers_d1),intent(in) :: tapers
        integer(int32),intent(in) :: lwin
        integer(int32),intent(in) :: k
        real(dp),dimension(admit_d0),intent(out) :: admit
        real(dp),dimension(corr_d0),intent(out) :: corr
        integer(int32),intent(in) :: mtdef
        real(dp),dimension(taper_wt_d0),intent(in) :: taper_wt
        if (taper_wt(1) < 0.0_dp) then
            call SHBiasAdmitCorr(sgt,sgg,stt,lmax,tapers,lwin,k,admit,corr,&
                                 mtdef=mtdef,exitstatus=exitstatus)
        else
            call SHBiasAdmitCorr(sgt,sgg,stt,lmax,tapers,lwin,k,admit,corr,&
                                 mtdef=mtdef,taper_wt=taper_wt,&
                                 exitstatus=exitstatus)
        end if
    end subroutine pySHBiasAdmitCorr

    subroutine pySHMTDebias(exitstatus,mtdebias,mtspectra,lmax,tapers,lwin,k,&
                            nl,lmid,n,taper_wt,mtdebias_d0,mtdebias_d1,&
                            taper_wt_d0,mtspectra_d0,mtspectra_d1,tapers_d0,&
                            tapers_d1,lmid_d0)
        use shtools, only: SHMTDebias
        use ftypes
        implicit none
        integer(int32),intent(in) :: mtdebias_d0
        integer(int32),intent(in) :: mtdebias_d1
        integer(int32),intent(in) :: taper_wt_d0
        integer(int32),intent(in) :: mtspectra_d0
        integer(int32),intent(in) :: mtspectra_d1
        integer(int32),intent(in) :: tapers_d0
        integer(int32),intent(in) :: tapers_d1
        integer(int32),intent(in) :: lmid_d0
        integer(int32),intent(out) :: exitstatus
        real(dp),dimension(mtdebias_d0,mtdebias_d1),intent(out) :: mtdebias
        real(dp),dimension(mtspectra_d0,mtspectra_d1),intent(in) :: mtspectra
        integer(int32),intent(in) :: lmax
        real(dp),dimension(tapers_d0,tapers_d1),intent(in) :: tapers
        integer(int32),intent(in) :: lwin
        integer(int32),intent(in) :: k
        integer(int32),intent(in) :: nl
        real(dp),dimension(lmid_d0),intent(out) :: lmid
        integer(int32),intent(out) :: n
        real(dp),dimension(taper_wt_d0),intent(in) :: taper_wt
        if (taper_wt(1) < 0.0_dp) then
            call SHMTDebias(mtdebias,mtspectra,lmax,tapers,lwin,k,nl,lmid,n,&
                            exitstatus=exitstatus)
        else
            call SHMTDebias(mtdebias,mtspectra,lmax,tapers,lwin,k,nl,lmid,n,&
                            taper_wt=taper_wt,exitstatus=exitstatus)
        end if
    end subroutine pySHMTDebias

    subroutine pySHMTVarOpt(exitstatus,l,tapers,taper_order,lwin,kmax,Sff,&
                            var_opt,var_unit,weight_opt,nocross,&
                            taper_order_d0,weight_opt_d0,weight_opt_d1,&
                            var_unit_d0,var_opt_d0,Sff_d0,tapers_d0,tapers_d1)
        use shtools, only: SHMTVarOpt
        use ftypes
        implicit none
        integer(int32),intent(in) :: taper_order_d0
        integer(int32),intent(in) :: weight_opt_d0
        integer(int32),intent(in) :: weight_opt_d1
        integer(int32),intent(in) :: var_unit_d0
        integer(int32),intent(in) :: var_opt_d0
        integer(int32),intent(in) :: Sff_d0
        integer(int32),intent(in) :: tapers_d0
        integer(int32),intent(in) :: tapers_d1
        integer(int32),intent(out) :: exitstatus
        integer(int32),intent(in) :: l
        real(dp),dimension(tapers_d0,tapers_d1),intent(in) :: tapers
        integer(int32),dimension(taper_order_d0),intent(in) :: taper_order
        integer(int32),intent(in) :: lwin
        integer(int32),intent(in) :: kmax
        real(dp),dimension(Sff_d0),intent(in) :: Sff
        real(dp),dimension(var_opt_d0),intent(out) :: var_opt
        real(dp),dimension(var_unit_d0),intent(out) :: var_unit
        real(dp),dimension(weight_opt_d0,weight_opt_d1),intent(out) ::&
                                                                     weight_opt
        integer(int32),intent(in) :: nocross
        call SHMTVarOpt(l,tapers,taper_order,lwin,kmax,Sff,var_opt,var_unit,&
                        weight_opt=weight_opt,nocross=nocross,&
                        exitstatus=exitstatus)
    end subroutine pySHMTVarOpt

    subroutine pySHMTVar(exitstatus,l,tapers,taper_order,Sff,kmax,lwin,&
                            variance,taper_wt,nocross,taper_order_d0,&
                            taper_wt_d0,Sff_d0,tapers_d0,tapers_d1)
        use shtools, only: SHMTVar
        use ftypes
        implicit none
        integer(int32),intent(in) :: taper_order_d0
        integer(int32),intent(in) :: taper_wt_d0
        integer(int32),intent(in) :: Sff_d0
        integer(int32),intent(in) :: tapers_d0
        integer(int32),intent(in) :: tapers_d1
        integer(int32),intent(out) :: exitstatus
        integer(int32),intent(in) :: l
        real(dp),dimension(tapers_d0,tapers_d1),intent(in) :: tapers
        integer(int32),dimension(taper_order_d0),intent(in) :: taper_order
        real(dp),dimension(Sff_d0),intent(in) :: Sff
        integer(int32),intent(in) :: kmax
        integer(int32),intent(in) :: lwin
        real(dp),intent(out) :: variance
        real(dp),dimension(taper_wt_d0),intent(in) :: taper_wt
        integer(int32),intent(in) :: nocross
        if (taper_wt(1) < 0.0_dp) then
            call SHMTVar(l,tapers,taper_order,lwin,kmax,Sff,variance,&
                         nocross=nocross,exitstatus=exitstatus)
        else
            call SHMTVar(l,tapers,taper_order,lwin,kmax,Sff,variance,&
                         taper_wt=taper_wt,nocross=nocross,&
                         exitstatus=exitstatus)
        end if
    end subroutine pySHMTVar

    function pySHSjkPG(incspectra,l,m,mprime,hj_real,hk_real,mj,mk,lwin,hkcc,&
                       hk_real_d0,incspectra_d0,hj_real_d0)
        use shtools, only: SHSjkPG
        use ftypes
        implicit none
        integer(int32),intent(in) :: hk_real_d0
        integer(int32),intent(in) :: incspectra_d0
        integer(int32),intent(in) :: hj_real_d0
        real(dp),dimension(incspectra_d0),intent(in) :: incspectra
        integer(int32),intent(in) :: l
        integer(int32),intent(in) :: m
        integer(int32),intent(in) :: mprime
        real(dp),dimension(hj_real_d0),intent(in) :: hj_real
        real(dp),dimension(hk_real_d0),intent(in) :: hk_real
        integer(int32),intent(in) :: mj
        integer(int32),intent(in) :: mk
        integer(int32),intent(in) :: lwin
        integer(int32),intent(in) :: hkcc
        complex(dp) :: pySHSjkPG
        pySHSjkPG=SHSjkPG(incspectra,l,m,mprime,hj_real,hk_real,mj,mk,lwin,&
                          hkcc)
    end function pySHSjkPG

    subroutine pySHReturnTapersMap(exitstatus,tapers,eigenvalues,dh_mask,n_dh,&
                                   lmax,sampling,ntapers,degrees,dh_mask_d0,&
                                   dh_mask_d1,tapers_d0,tapers_d1,&
                                   eigenvalues_d0,degrees_d0)
        use shtools, only: SHReturnTapersMap
        use ftypes
        implicit none
        integer(int32),intent(in) :: dh_mask_d0
        integer(int32),intent(in) :: dh_mask_d1
        integer(int32),intent(in) :: tapers_d0
        integer(int32),intent(in) :: tapers_d1
        integer(int32),intent(in) :: eigenvalues_d0
        integer(int32),intent(in) :: degrees_d0
        integer(int32),intent(out) :: exitstatus
        real(dp),dimension(tapers_d0,tapers_d1),intent(out) :: tapers
        real(dp),dimension(eigenvalues_d0),intent(out) :: eigenvalues
        integer(int32),dimension(dh_mask_d0,dh_mask_d1),intent(in) :: dh_mask
        integer(int32),intent(in) :: n_dh
        integer(int32),intent(in) :: sampling
        integer(int32),intent(in) :: lmax
        integer(int32),intent(in) :: ntapers
        integer(int32),dimension(degrees_d0),intent(in) :: degrees
        call SHReturnTapersMap(tapers,eigenvalues,dh_mask,n_dh,lmax,&
                               sampling=sampling,ntapers=ntapers,&
                               degrees=degrees,exitstatus=exitstatus)
    end subroutine pySHReturnTapersMap

    subroutine pySHBiasKMask(exitstatus,tapers,lwin,k,incspectra,ldata,&
                             outcspectra,taper_wt,save_cg,taper_wt_d0,&
                             tapers_d0,tapers_d1,incspectra_d0,outcspectra_d0)
        use shtools, only: SHBiasKMask
        use ftypes
        implicit none
        integer(int32),intent(in) :: taper_wt_d0
        integer(int32),intent(in) :: tapers_d0
        integer(int32),intent(in) :: tapers_d1
        integer(int32),intent(in) :: incspectra_d0
        integer(int32),intent(in) :: outcspectra_d0
        integer(int32),intent(out) :: exitstatus
        real(dp),dimension(tapers_d0,tapers_d1),intent(in) :: tapers
        integer(int32),intent(in) :: lwin
        integer(int32),intent(in) :: k
        real(dp),dimension(incspectra_d0),intent(in) :: incspectra
        integer(int32),intent(in) :: ldata
        real(dp),dimension(outcspectra_d0),intent(out) :: outcspectra
        real(dp),dimension(taper_wt_d0),intent(in) :: taper_wt
        integer(int32),intent(in) :: save_cg
        if (taper_wt(1) < 0.0_dp) then
            call SHBiasKMask(tapers,lwin,k,incspectra,ldata,outcspectra,&
                             save_cg=save_cg,exitstatus=exitstatus)
        else
            call SHBiasKMask(tapers,lwin,k,incspectra,ldata,outcspectra,&
                             taper_wt=taper_wt,save_cg=save_cg,&
                             exitstatus=exitstatus)
        end if
    end subroutine pySHBiasKMask

    subroutine pySHMultiTaperMaskSE(exitstatus,mtse,sd,sh,lmax,tapers,lmaxt,k,&
                                    taper_wt,norm,csphase,&
                                    taper_wt_d0,sh_d0,sh_d1,sh_d2,&
                                    tapers_d0,tapers_d1,mtse_d0,sd_d0)
        use shtools, only: SHMultiTaperMaskSE
        use ftypes
        implicit none
        integer(int32),intent(in) :: taper_wt_d0
        integer(int32),intent(in) :: sh_d0
        integer(int32),intent(in) :: sh_d1
        integer(int32),intent(in) :: sh_d2
        integer(int32),intent(in) :: tapers_d0
        integer(int32),intent(in) :: tapers_d1
        integer(int32),intent(in) :: mtse_d0
        integer(int32),intent(in) :: sd_d0
        integer(int32),intent(out) :: exitstatus
        real(dp),dimension(mtse_d0),intent(out) :: mtse
        real(dp),dimension(sd_d0),intent(out) :: sd
        real(dp),dimension(sh_d0,sh_d1,sh_d2),intent(in) :: sh
        integer(int32),intent(in) :: lmax
        real(dp),dimension(tapers_d0,tapers_d1),intent(in) :: tapers
        integer(int32),intent(in) :: lmaxt
        integer(int32),intent(in) :: k
        real(dp),dimension(taper_wt_d0),intent(in) :: taper_wt
        integer(int32),intent(in) :: norm
        integer(int32),intent(in) :: csphase
        if (taper_wt(1) < 0.0_dp) then
            call SHMultiTaperMaskSE(mtse,sd,sh,lmax,tapers,lmaxt,k,norm=norm,&
                                    csphase=csphase,exitstatus=exitstatus)
        else
            call SHMultiTaperMaskSE(mtse,sd,sh,lmax,tapers,lmaxt,k,&
                                    taper_wt=taper_wt,norm=norm,&
                                    csphase=csphase,exitstatus=exitstatus)
        end if
    end subroutine pySHMultiTaperMaskSE

    subroutine pySHMultiTaperMaskCSE(exitstatus,mtse,sd,sh1,lmax1,sh2,lmax2,&
                                     tapers,lmaxt,k,taper_wt,norm,csphase,&
                                     sh1_d0,sh1_d1,sh1_d2,sh2_d0,sh2_d1,&
                                     sh2_d2,taper_wt_d0,tapers_d0,tapers_d1,&
                                     sd_d0,mtse_d0)
        use shtools, only: SHMultiTaperMaskCSE
        use ftypes
        implicit none
        integer(int32),intent(in) :: sh1_d0
        integer(int32),intent(in) :: sh1_d1
        integer(int32),intent(in) :: sh1_d2
        integer(int32),intent(in) :: sh2_d0
        integer(int32),intent(in) :: sh2_d1
        integer(int32),intent(in) :: sh2_d2
        integer(int32),intent(in) :: taper_wt_d0
        integer(int32),intent(in) :: tapers_d0
        integer(int32),intent(in) :: tapers_d1
        integer(int32),intent(in) :: sd_d0
        integer(int32),intent(in) :: mtse_d0
        integer(int32),intent(out) :: exitstatus
        real(dp),dimension(mtse_d0),intent(out) :: mtse
        real(dp),dimension(sd_d0),intent(out) :: sd
        real(dp),dimension(sh1_d0,sh1_d1,sh1_d2),intent(in) :: sh1
        integer(int32),intent(in) :: lmax1
        real(dp),dimension(sh2_d0,sh2_d1,sh2_d2),intent(in) :: sh2
        integer(int32),intent(in) :: lmax2
        real(dp),dimension(tapers_d0,tapers_d1),intent(in) :: tapers
        integer(int32),intent(in) :: lmaxt
        integer(int32),intent(in) :: k
        real(dp),dimension(taper_wt_d0),intent(in) :: taper_wt
        integer(int32),intent(in) :: norm
        integer(int32),intent(in) :: csphase
        if(taper_wt(1) < 0.0_dp) then
            call SHMultiTaperMaskCSE(mtse,sd,sh1,lmax1,sh2,lmax2,tapers,&
                                     lmaxt,k,norm=norm,csphase=csphase,&
                                     exitstatus=exitstatus)
        else
            call SHMultiTaperMaskCSE(mtse,sd,sh1,lmax1,sh2,lmax2,tapers,&
                                     lmaxt,k,taper_wt=taper_wt,norm=norm,&
                                     csphase=csphase,exitstatus=exitstatus)
        end if
    end subroutine pySHMultiTaperMaskCSE

    subroutine pyComputeDMap(exitstatus,Dij,dh_mask,n_dh,lmax,sampling,&
                             degrees,dh_mask_d0,dh_mask_d1,Dij_d0,Dij_d1,&
                             degrees_d0)
        use shtools, only: ComputeDMap
        use ftypes
        implicit none
        integer(int32),intent(in) :: dh_mask_d0
        integer(int32),intent(in) :: dh_mask_d1
        integer(int32),intent(in) :: Dij_d0
        integer(int32),intent(in) :: Dij_d1
        integer(int32),intent(in) :: degrees_d0
        integer(int32),intent(out) :: exitstatus
        real(dp),dimension(Dij_d0,Dij_d1),intent(out) :: Dij
        integer(int32),dimension(dh_mask_d0,dh_mask_d1),intent(in) :: dh_mask
        integer(int32),intent(in) :: n_dh
        integer(int32),intent(in) :: lmax
        integer(int32),intent(in) :: sampling
        integer(int32),intent(in),dimension(degrees_d0) :: degrees
        call ComputeDMap(Dij,dh_mask,n_dh,lmax,sampling=sampling,&
                         degrees=degrees,exitstatus=exitstatus)
    end subroutine pyComputeDMap

    subroutine pyCurve2Mask(exitstatus,dhgrid,n,sampling,profile,nprofile,NP,&
                            extend,profile_d0,profile_d1,dhgrid_d0,dhgrid_d1)
        use shtools, only: Curve2Mask
        use ftypes
        implicit none
        integer(int32),intent(in) :: profile_d0
        integer(int32),intent(in) :: profile_d1
        integer(int32),intent(in) :: dhgrid_d0
        integer(int32),intent(in) :: dhgrid_d1
        integer(int32),intent(out) :: exitstatus
        integer(int32),dimension(dhgrid_d0,dhgrid_d1),intent(out) :: dhgrid
        integer(int32),intent(in) :: n
        integer(int32),intent(in) :: sampling
        real(dp),dimension(profile_d0,profile_d1),intent(in) :: profile
        integer(int32),intent(in) :: nprofile
        integer(int32),intent(in) :: NP
        integer(int32),intent(in) :: extend
        call Curve2Mask(dhgrid,n,sampling,profile,nprofile,NP,extend=extend,&
                        exitstatus=exitstatus)
    end subroutine pyCurve2Mask

    subroutine pySHBias(exitstatus,Shh,lwin,incspectra,ldata,outcspectra,&
                        save_cg,Shh_d0,incspectra_d0,outcspectra_d0)
        use shtools, only: SHBias
        use ftypes
        implicit none
        integer(int32),intent(in) :: Shh_d0
        integer(int32),intent(in) :: incspectra_d0
        integer(int32),intent(in) :: outcspectra_d0
        integer(int32),intent(out) :: exitstatus
        real(dp),dimension(Shh_d0),intent(in) :: Shh
        integer(int32),intent(in) :: lwin
        real(dp),dimension(incspectra_d0),intent(in) :: incspectra
        integer(int32),intent(in) :: ldata
        real(dp),dimension(outcspectra_d0),intent(out) :: outcspectra
        integer(int32),intent(in) :: save_cg
        call SHBias(Shh,lwin,incspectra,ldata,outcspectra,save_cg=save_cg,&
                    exitstatus=exitstatus)
    end subroutine pySHBias

    subroutine pySphericalCapCoef(exitstatus,coef,theta,lmax,coef_d0)
        use shtools, only: SphericalCapCoef
        use ftypes
        implicit none
        integer(int32),intent(in) :: coef_d0
        integer(int32),intent(out) :: exitstatus
        real(dp),dimension(coef_d0),intent(out) :: coef
        real(dp),intent(in) :: theta
        integer(int32),intent(in) :: lmax
        call SphericalCapCoef(coef,theta,lmax=lmax,exitstatus=exitstatus)
    end subroutine pySphericalCapCoef

    subroutine pyMakeGravGridDH(exitstatus,cilm,lmax,gm,r0,a,f,rad,theta,phi,&
                                total,pot,n,sampling,lmax_calc,omega,&
                                normal_gravity,extend,phi_d0,phi_d1,total_d0,&
                                total_d1,rad_d0,rad_d1,cilm_d0,cilm_d1,&
                                cilm_d2,theta_d0,theta_d1,pot_d0,pot_d1)
        use shtools, only: MakeGravGridDH
        use ftypes
        implicit none
        integer(int32),intent(in) :: phi_d0
        integer(int32),intent(in) :: phi_d1
        integer(int32),intent(in) :: total_d0
        integer(int32),intent(in) :: total_d1
        integer(int32),intent(in) :: rad_d0
        integer(int32),intent(in) :: rad_d1
        integer(int32),intent(in) :: cilm_d0
        integer(int32),intent(in) :: cilm_d1
        integer(int32),intent(in) :: cilm_d2
        integer(int32),intent(in) :: theta_d0
        integer(int32),intent(in) :: theta_d1
        integer(int32),intent(in) :: pot_d0
        integer(int32),intent(in) :: pot_d1
        integer(int32),intent(out) :: exitstatus
        real(dp),dimension(cilm_d0,cilm_d1,cilm_d2),intent(in) :: cilm
        integer(int32),intent(in) :: lmax
        real(dp),intent(in) :: gm
        real(dp),intent(in) :: r0
        real(dp),intent(in) :: a
        real(dp),intent(in) :: f
        real(dp),dimension(rad_d0,rad_d1),intent(out) :: rad
        real(dp),dimension(theta_d0,theta_d1),intent(out) :: theta
        real(dp),dimension(phi_d0,phi_d1),intent(out) :: phi
        real(dp),dimension(total_d0,total_d1),intent(out) :: total
        real(dp),dimension(pot_d0,pot_d1),intent(out) :: pot
        integer(int32),intent(out) :: n
        integer(int32),intent(in) :: sampling
        integer(int32),intent(in) :: lmax_calc
        real(dp),intent(in) :: omega
        integer(int32),intent(in) :: normal_gravity
        integer(int32),intent(in) :: extend
        call MakeGravGridDH(cilm,lmax,gm,r0,a,f,rad,theta,phi,total,n,&
                            sampling=sampling,lmax_calc=lmax_calc,omega=omega,&
                            normal_gravity=normal_gravity,pot=pot,&
                            extend=extend,exitstatus=exitstatus)
    end subroutine pyMakeGravGridDH

    subroutine pyMakeGravGridPoint(vec,cilm,lmax,gm,r0,r,lat,lon,omega,&
                                   dealloc,cilm_d0,cilm_d1,cilm_d2)
        use shtools, only: MakeGravGridPoint
        use ftypes
        implicit none
        integer(int32),intent(in) :: cilm_d0
        integer(int32),intent(in) :: cilm_d1
        integer(int32),intent(in) :: cilm_d2
        real(dp),dimension(cilm_d0,cilm_d1,cilm_d2),intent(in) :: cilm
        integer(int32),intent(in) :: lmax
        real(dp),intent(in) :: gm
        real(dp),intent(in) :: r0
        real(dp),intent(in) :: r
        real(dp),intent(in) :: lat
        real(dp),intent(in) :: lon
        real(dp),intent(in) :: omega
        integer(int32),intent(in) :: dealloc
        real(dp),dimension(3),intent(out) :: vec
        vec=MakeGravGridPoint(cilm,lmax,gm,r0,r,lat,lon,omega=omega,&
                              dealloc=dealloc)
    end subroutine pyMakeGravGridPoint

    subroutine pyMakeGravGradGridDH(exitstatus,cilm,lmax,gm,r0,a,f,vxx,vyy,&
                                    vzz,vxy,vxz,vyz,n,sampling,lmax_calc,&
                                    extend,vyz_d0,vyz_d1,vyy_d0,vyy_d1,&
                                    cilm_d0,cilm_d1,cilm_d2,vzz_d0,vzz_d1,&
                                    vxy_d0,vxy_d1,vxx_d0,vxx_d1,vxz_d0,vxz_d1)
        use shtools, only: MakeGravGradGridDH
        use ftypes
        implicit none
        integer(int32),intent(in) :: vyz_d0
        integer(int32),intent(in) :: vyz_d1
        integer(int32),intent(in) :: vyy_d0
        integer(int32),intent(in) :: vyy_d1
        integer(int32),intent(in) :: cilm_d0
        integer(int32),intent(in) :: cilm_d1
        integer(int32),intent(in) :: cilm_d2
        integer(int32),intent(in) :: vzz_d0
        integer(int32),intent(in) :: vzz_d1
        integer(int32),intent(in) :: vxy_d0
        integer(int32),intent(in) :: vxy_d1
        integer(int32),intent(in) :: vxx_d0
        integer(int32),intent(in) :: vxx_d1
        integer(int32),intent(in) :: vxz_d0
        integer(int32),intent(in) :: vxz_d1
        integer(int32),intent(out) :: exitstatus
        real(dp),dimension(cilm_d0,cilm_d1,cilm_d2),intent(in) :: cilm
        integer(int32),intent(in) :: lmax
        real(dp),intent(in) :: gm
        real(dp),intent(in) :: r0
        real(dp),intent(in) :: a
        real(dp),intent(in) :: f
        real(dp),dimension(vxx_d0,vxx_d1),intent(out) :: vxx
        real(dp),dimension(vyy_d0,vyy_d1),intent(out) :: vyy
        real(dp),dimension(vzz_d0,vzz_d1),intent(out) :: vzz
        real(dp),dimension(vxy_d0,vxy_d1),intent(out) :: vxy
        real(dp),dimension(vxz_d0,vxz_d1),intent(out) :: vxz
        real(dp),dimension(vyz_d0,vyz_d1),intent(out) :: vyz
        integer(int32),intent(out) :: n
        integer(int32),intent(in) :: sampling
        integer(int32),intent(in) :: lmax_calc
        integer(int32),intent(in) :: extend
        call MakeGravGradGridDH(cilm,lmax,gm,r0,a,f,vxx,vyy,vzz,vxy,vxz,vyz,n,&
                                sampling=sampling,lmax_calc=lmax_calc,&
                                extend=extend,exitstatus=exitstatus)
    end subroutine pyMakeGravGradGridDH

    subroutine pyMakeMagGradGridDH(exitstatus,cilm,lmax,r0,a,f,vxx,vyy,vzz,&
                                    vxy,vxz,vyz,n,sampling,lmax_calc,extend,&
                                    vyz_d0,vyz_d1,vyy_d0,vyy_d1,cilm_d0,&
                                    cilm_d1,cilm_d2,vzz_d0,vzz_d1,vxy_d0,&
                                    vxy_d1,vxx_d0,vxx_d1,vxz_d0,vxz_d1)
        use shtools, only: MakeMagGradGridDH
        use ftypes
        implicit none
        integer(int32),intent(in) :: vyz_d0
        integer(int32),intent(in) :: vyz_d1
        integer(int32),intent(in) :: vyy_d0
        integer(int32),intent(in) :: vyy_d1
        integer(int32),intent(in) :: cilm_d0
        integer(int32),intent(in) :: cilm_d1
        integer(int32),intent(in) :: cilm_d2
        integer(int32),intent(in) :: vzz_d0
        integer(int32),intent(in) :: vzz_d1
        integer(int32),intent(in) :: vxy_d0
        integer(int32),intent(in) :: vxy_d1
        integer(int32),intent(in) :: vxx_d0
        integer(int32),intent(in) :: vxx_d1
        integer(int32),intent(in) :: vxz_d0
        integer(int32),intent(in) :: vxz_d1
        integer(int32),intent(out) :: exitstatus
        real(dp),dimension(cilm_d0,cilm_d1,cilm_d2),intent(in) :: cilm
        integer(int32),intent(in) :: lmax
        real(dp),intent(in) :: r0
        real(dp),intent(in) :: a
        real(dp),intent(in) :: f
        real(dp),dimension(vxx_d0,vxx_d1),intent(out) :: vxx
        real(dp),dimension(vyy_d0,vyy_d1),intent(out) :: vyy
        real(dp),dimension(vzz_d0,vzz_d1),intent(out) :: vzz
        real(dp),dimension(vxy_d0,vxy_d1),intent(out) :: vxy
        real(dp),dimension(vxz_d0,vxz_d1),intent(out) :: vxz
        real(dp),dimension(vyz_d0,vyz_d1),intent(out) :: vyz
        integer(int32),intent(out) :: n
        integer(int32),intent(in) :: sampling
        integer(int32),intent(in) :: lmax_calc
        integer(int32),intent(in) :: extend
        call MakeMagGradGridDH(cilm,lmax,r0,a,f,vxx,vyy,vzz,vxy,vxz,vyz,n,&
                                sampling=sampling,lmax_calc=lmax_calc,&
                                extend=extend,exitstatus=exitstatus)
    end subroutine pyMakeMagGradGridDH

    subroutine pyMakeMagGridPoint(vec,cilm,lmax,a,r,lat,lon,dealloc,cilm_d0,&
                                  cilm_d1,cilm_d2)
        use shtools, only: MakeMagGridPoint
        use ftypes
        implicit none
        integer(int32),intent(in) :: cilm_d0
        integer(int32),intent(in) :: cilm_d1
        integer(int32),intent(in) :: cilm_d2
        real(dp),dimension(cilm_d0,cilm_d1,cilm_d2),intent(in) :: cilm
        integer(int32),intent(in) :: lmax
        real(dp),intent(in) :: a
        real(dp),intent(in) :: r
        real(dp),intent(in) :: lat
        real(dp),intent(in) :: lon
        integer(int32),intent(in) :: dealloc
        real(dp),dimension(3),intent(out) :: vec
        vec=MakeMagGridPoint(cilm,lmax,a,r,lat,lon,dealloc=dealloc)
    end subroutine pyMakeMagGridPoint

    subroutine pyMakeGeoidGridDH(exitstatus,geoid,cilm,lmax,r0pot,GM,PotRef,&
                                 omega,r,sampling,order,nlat,nlong,lmax_calc,&
                                 a,f,extend,cilm_d0,cilm_d1,cilm_d2,geoid_d0,&
                                 geoid_d1)
        use shtools, only: MakeGeoidGrid
        use ftypes
        implicit none
        integer(int32),intent(in) :: cilm_d0
        integer(int32),intent(in) :: cilm_d1
        integer(int32),intent(in) :: cilm_d2
        integer(int32),intent(in) :: geoid_d0
        integer(int32),intent(in) :: geoid_d1
        integer(int32),intent(out) :: exitstatus
        real(dp),dimension(geoid_d0,geoid_d1),intent(out) :: geoid
        real(dp),dimension(cilm_d0,cilm_d1,cilm_d2),intent(in) :: cilm
        integer(int32),intent(in) :: lmax
        real(dp),intent(in) :: r0pot
        real(dp),intent(in) :: GM
        real(dp),intent(in) :: PotRef
        real(dp),intent(in) :: omega
        real(dp),intent(in) :: r
        integer(int32),intent(in) :: sampling
        integer(int32),intent(in) :: order
        integer(int32),intent(out) :: nlat
        integer(int32),intent(out) :: nlong
        integer(int32),intent(in) :: lmax_calc
        real(dp),intent(in) :: a
        real(dp),intent(in) :: f
        integer(int32),intent(in) :: extend
        if (sampling == 1) then
            call MakeGeoidGrid(geoid,cilm,lmax,r0pot,GM,PotRef,omega,r,2,&
                               order,nlat,nlong,lmax_calc=lmax_calc,a=a,f=f,&
                               extend=extend,exitstatus=exitstatus)
        else if (sampling == 2) then
            call MakeGeoidGrid(geoid,cilm,lmax,r0pot,GM,PotRef,omega,r,3,&
                               order,nlat,nlong,lmax_calc=lmax_calc,a=a,f=f,&
                               extend=extend,exitstatus=exitstatus)
        end if
    end subroutine pyMakeGeoidGridDH

    subroutine pyCilmPlusDH(exitstatus,cilm,gridin,lmax,nmax,mass,d,rho,&
                            sampling,n,gridin_d0,gridin_d1,cilm_d0,cilm_d1,&
                            cilm_d2)
        use shtools, only: CilmPlus
        use ftypes
        implicit none
        integer(int32),intent(in) :: gridin_d0
        integer(int32),intent(in) :: gridin_d1
        integer(int32),intent(in) :: cilm_d0
        integer(int32),intent(in) :: cilm_d1
        integer(int32),intent(in) :: cilm_d2
        integer(int32),intent(out) :: exitstatus
        real(dp),dimension(cilm_d0,cilm_d1,cilm_d2),intent(out) :: cilm
        real(dp),dimension(gridin_d0,gridin_d1),intent(in) :: gridin
        integer(int32),intent(in) :: lmax
        integer(int32),intent(in) :: nmax
        real(dp),intent(in) :: mass
        real(dp),intent(out) :: d
        real(dp),intent(in) :: rho
        integer(int32),intent(in) :: sampling
        integer(int32),intent(in) :: n
        if (sampling == 1) then
            call CilmPlus(cilm,gridin,lmax,nmax,mass,d,rho,2,n=n,&
                          exitstatus=exitstatus)
        else 
            call CilmPlus(cilm,gridin,lmax,nmax,mass,d,rho,3,n=n,&
                          exitstatus=exitstatus)
        end if
    end subroutine pyCilmPlusDH

    subroutine pyCilmMinusDH(exitstatus,cilm,gridin,lmax,nmax,mass,d,rho,&
                             sampling,n,gridin_d0,gridin_d1,cilm_d0,cilm_d1,&
                             cilm_d2)
        use shtools, only: CilmMinus
        use ftypes
        implicit none
        integer(int32),intent(in) :: gridin_d0
        integer(int32),intent(in) :: gridin_d1
        integer(int32),intent(in) :: cilm_d0
        integer(int32),intent(in) :: cilm_d1
        integer(int32),intent(in) :: cilm_d2
        integer(int32),intent(out) :: exitstatus
        real(dp),dimension(cilm_d0,cilm_d1,cilm_d2),intent(out) :: cilm
        real(dp),dimension(gridin_d0,gridin_d1),intent(in) :: gridin
        integer(int32),intent(in) :: lmax
        integer(int32),intent(in) :: nmax
        real(dp),intent(in) :: mass
        real(dp),intent(out) :: d
        real(dp),intent(in) :: rho
        integer(int32),intent(in) :: sampling
        integer(int32),intent(in) :: n
        if (sampling == 1) then
            call CilmMinus(cilm,gridin,lmax,nmax,mass,d,rho,2,n=n,&
                           exitstatus=exitstatus)
        else 
            call CilmMinus(cilm,gridin,lmax,nmax,mass,d,rho,3,n=n,&
                           exitstatus=exitstatus)
        end if
    end subroutine pyCilmMinusDH

    subroutine pyCilmPlusRhoHDH(exitstatus,cilm,gridin,lmax,nmax,mass,d,rho,&
                                sampling,n,gridin_d0,gridin_d1,cilm_d0,&
                                cilm_d1,cilm_d2,rho_d0,rho_d1)
        use shtools, only: CilmPlusRhoH
        use ftypes
        implicit none
        integer(int32),intent(in) :: gridin_d0
        integer(int32),intent(in) :: gridin_d1
        integer(int32),intent(in) :: cilm_d0
        integer(int32),intent(in) :: cilm_d1
        integer(int32),intent(in) :: cilm_d2
        integer(int32),intent(in) :: rho_d0
        integer(int32),intent(in) :: rho_d1
        integer(int32),intent(out) :: exitstatus
        real(dp),dimension(cilm_d0,cilm_d1,cilm_d2),intent(out) :: cilm
        real(dp),dimension(gridin_d0,gridin_d1),intent(in) :: gridin
        integer(int32),intent(in) :: lmax
        integer(int32),intent(in) :: nmax
        real(dp),intent(in) :: mass
        real(dp),intent(out) :: d
        real(dp),dimension(rho_d0,rho_d1),intent(in) :: rho
        integer(int32),intent(in) :: sampling
        integer(int32),intent(in) :: n
        if (sampling == 1) then
            call CilmPlusRhoH(cilm,gridin,lmax,nmax,mass,d,rho,2,n=n,&
                              exitstatus=exitstatus)
        else
            call CilmPlusRhoH(cilm,gridin,lmax,nmax,mass,d,rho,3,n=n,&
                              exitstatus=exitstatus)
        end if
    end subroutine pyCilmPlusRhoHDH

    subroutine pyCilmMinusRhoHDH(exitstatus,cilm,gridin,lmax,nmax,mass,d,rho,&
                                 sampling,n,gridin_d0,gridin_d1,cilm_d0,&
                                 cilm_d1,cilm_d2,rho_d0,rho_d1)
        use shtools, only: CilmMinusRhoH
        use ftypes
        implicit none
        integer(int32),intent(in) :: gridin_d0
        integer(int32),intent(in) :: gridin_d1
        integer(int32),intent(in) :: cilm_d0
        integer(int32),intent(in) :: cilm_d1
        integer(int32),intent(in) :: cilm_d2
        integer(int32),intent(in) :: rho_d0
        integer(int32),intent(in) :: rho_d1
        integer(int32),intent(out) :: exitstatus
        real(dp),dimension(cilm_d0,cilm_d1,cilm_d2),intent(out) :: cilm
        real(dp),dimension(gridin_d0,gridin_d1),intent(in) :: gridin
        integer(int32),intent(in) :: lmax
        integer(int32),intent(in) :: nmax
        real(dp),intent(in) :: mass
        real(dp),intent(out) :: d
        real(dp),dimension(rho_d0,rho_d1),intent(in) :: rho
        integer(int32),intent(in) :: sampling
        integer(int32),intent(in) :: n
        if (sampling == 1) then
            call CilmMinusRhoH(cilm,gridin,lmax,nmax,mass,d,rho,2,n=n,&
                               exitstatus=exitstatus)
        else
            call CilmMinusRhoH(cilm,gridin,lmax,nmax,mass,d,rho,3,n=n,&
                               exitstatus=exitstatus)
        end if
    end subroutine pyCilmMinusRhoHDH

    subroutine pyBAtoHilmDH(exitstatus,cilm,ba,griddh,lmax,nmax,mass,r0,rho,&
                            sampling,filter_type,filter_deg,lmax_calc,ba_d0,&
                            ba_d1,ba_d2,griddh_d0,griddh_d1,cilm_d0,cilm_d1,&
                            cilm_d2)
        use shtools, only: BAtoHilm
        use ftypes
        implicit none
        integer(int32),intent(in) :: ba_d0
        integer(int32),intent(in) :: ba_d1
        integer(int32),intent(in) :: ba_d2
        integer(int32),intent(in) :: griddh_d0
        integer(int32),intent(in) :: griddh_d1
        integer(int32),intent(in) :: cilm_d0
        integer(int32),intent(in) :: cilm_d1
        integer(int32),intent(in) :: cilm_d2
        integer(int32),intent(out) :: exitstatus
        real(dp),dimension(cilm_d0,cilm_d1,cilm_d2),intent(out) :: cilm
        real(dp),dimension(ba_d0,ba_d1,ba_d2),intent(in) :: ba
        real(dp),dimension(griddh_d0,griddh_d1),intent(in) :: griddh
        integer(int32),intent(in) :: lmax
        integer(int32),intent(in) :: nmax
        real(dp),intent(in) :: mass
        real(dp),intent(in) :: r0
        real(dp),intent(in) :: rho
        integer(int32),intent(in) :: sampling
        integer(int32),intent(in) :: filter_type
        integer(int32),intent(in) :: filter_deg
        integer(int32),intent(in) :: lmax_calc
        if (sampling == 1) then
            call BAtoHilm(cilm,ba,griddh,lmax,nmax,mass,r0,rho,2,&
                          filter_type=filter_type,filter_deg=filter_deg,&
                          lmax_calc=lmax_calc,exitstatus=exitstatus)
        else
            call BAtoHilm(cilm,ba,griddh,lmax,nmax,mass,r0,rho,3,&
                          filter_type=filter_type,filter_deg=filter_deg,&
                          lmax_calc=lmax_calc,exitstatus=exitstatus)
        end if
    end subroutine pyBAtoHilmDH
    
    subroutine pyBAtoHilmRhoHDH(exitstatus,cilm,ba,griddh,rho,lmax,nmax,mass,&
                                r0,sampling,filter_type,filter_deg,lmax_calc,&
                                ba_d0,ba_d1,ba_d2,griddh_d0,griddh_d1,cilm_d0,&
                                cilm_d1,cilm_d2,rho_d0,rho_d1)
        use shtools, only: BAtoHilmRhoH
        use ftypes
        implicit none
        integer(int32),intent(in) :: ba_d0
        integer(int32),intent(in) :: ba_d1
        integer(int32),intent(in) :: ba_d2
        integer(int32),intent(in) :: griddh_d0
        integer(int32),intent(in) :: griddh_d1
        integer(int32),intent(in) :: cilm_d0
        integer(int32),intent(in) :: cilm_d1
        integer(int32),intent(in) :: cilm_d2
        integer(int32),intent(in) :: rho_d0
        integer(int32),intent(in) :: rho_d1
        integer(int32),intent(out) :: exitstatus
        real(dp),dimension(cilm_d0,cilm_d1,cilm_d2),intent(out) :: cilm
        real(dp),dimension(ba_d0,ba_d1,ba_d2),intent(in) :: ba
        real(dp),dimension(griddh_d0,griddh_d1),intent(in) :: griddh
        integer(int32),intent(in) :: lmax
        integer(int32),intent(in) :: nmax
        real(dp),intent(in) :: mass
        real(dp),intent(in) :: r0
        real(dp),dimension(rho_d0,rho_d1),intent(in) :: rho
        integer(int32),intent(in) :: sampling
        integer(int32),intent(in) :: filter_type
        integer(int32),intent(in) :: filter_deg
        integer(int32),intent(in) :: lmax_calc
        if (sampling == 1) then
            call BAtoHilmRhoH(cilm,ba,griddh,lmax,nmax,mass,r0,rho,2,&
                              filter_type=filter_type,filter_deg=filter_deg,&
                              lmax_calc=lmax_calc,exitstatus=exitstatus)
        else
            call BAtoHilmRhoH(cilm,ba,griddh,lmax,nmax,mass,r0,rho,3,&
                              filter_type=filter_type,filter_deg=filter_deg,&
                              lmax_calc=lmax_calc,exitstatus=exitstatus)
        end if
    end subroutine pyBAtoHilmRhoHDH

    function pyDownContFilterMA(l,half,r,d)
        use shtools, only: DownContFilterMA
        use ftypes
        implicit none
        integer(int32),intent(in) :: l
        integer(int32),intent(in) :: half
        real(dp),intent(in) :: r
        real(dp),intent(in) :: d
        real(dp) :: pyDownContFilterMA
        pyDownContFilterMA=DownContFilterMA(l,half,r,d)
    end function pyDownContFilterMA

    function pyDownContFilterMC(l,half,r,d)
        use shtools, only: DownContFilterMC
        use ftypes
        implicit none
        integer(int32),intent(in) :: l
        integer(int32),intent(in) :: half
        real(dp),intent(in) :: r
        real(dp),intent(in) :: d
        real(dp) :: pyDownContFilterMC
        pyDownContFilterMC=DownContFilterMC(l,half,r,d)
    end function pyDownContFilterMC

    function pyNormalGravity(geocentric_lat,gm,omega,a,b)
        use shtools, only: NormalGravity
        use ftypes
        implicit none
        real(dp),intent(in) :: geocentric_lat
        real(dp),intent(in) :: gm
        real(dp),intent(in) :: omega
        real(dp),intent(in) :: a
        real(dp),intent(in) :: b
        real(dp) :: pyNormalGravity
        pyNormalGravity=NormalGravity(geocentric_lat,gm,omega,a,b)
    end function pyNormalGravity

    subroutine pyMakeMagGridDH(exitstatus,cilm,lmax,r0,a,f,rad_grid,&
                               theta_grid,phi_grid,total_grid,pot_grid,n,&
                               sampling,lmax_calc,extend,total_grid_d0,&
                               total_grid_d1,cilm_d0,cilm_d1,cilm_d2,&
                               rad_grid_d0,rad_grid_d1,theta_grid_d0,&
                               theta_grid_d1,phi_grid_d0,phi_grid_d1,&
                               pot_grid_d0,pot_grid_d1)
        use shtools, only: MakeMagGridDH
        use ftypes
        implicit none
        integer(int32),intent(in) :: total_grid_d0
        integer(int32),intent(in) :: total_grid_d1
        integer(int32),intent(in) :: cilm_d0
        integer(int32),intent(in) :: cilm_d1
        integer(int32),intent(in) :: cilm_d2
        integer(int32),intent(in) :: rad_grid_d0
        integer(int32),intent(in) :: rad_grid_d1
        integer(int32),intent(in) :: theta_grid_d0
        integer(int32),intent(in) :: theta_grid_d1
        integer(int32),intent(in) :: phi_grid_d0
        integer(int32),intent(in) :: phi_grid_d1
        integer(int32),intent(in) :: pot_grid_d0
        integer(int32),intent(in) :: pot_grid_d1
        integer(int32),intent(out) :: exitstatus
        real(dp),dimension(cilm_d0,cilm_d1,cilm_d2),intent(in) :: cilm
        integer(int32),intent(in) :: lmax
        real(dp),intent(in) :: r0
        real(dp),intent(in) :: a
        real(dp),intent(in) :: f
        real(dp),dimension(rad_grid_d0,rad_grid_d1),intent(out) :: rad_grid
        real(dp),dimension(theta_grid_d0,theta_grid_d1),intent(out) :: &
                                                                     theta_grid
        real(dp),dimension(phi_grid_d0,phi_grid_d1),intent(out) :: phi_grid
        real(dp),dimension(total_grid_d0,total_grid_d1),intent(out) :: &
                                                                     total_grid
        real(dp),dimension(pot_grid_d0,pot_grid_d1),intent(out) :: pot_grid
        integer(int32),intent(out) :: n
        integer(int32),intent(in) :: sampling
        integer(int32),intent(in) :: lmax_calc
        integer(int32),intent(in) :: extend
        call MakeMagGridDH(cilm,lmax,r0,a,f,rad_grid,theta_grid,phi_grid,&
                           total_grid,n,sampling=sampling,lmax_calc=lmax_calc,&
                           pot_grid=pot_grid,extend=extend,&
                           exitstatus=exitstatus)
    end subroutine pyMakeMagGridDH

    subroutine pyMakeCircleCoord(exitstatus,coord,lat,lon,theta0,cinterval,&
                                 cnum,coord_d0,coord_d1)
        use shtools, only: MakeCircleCoord
        use ftypes
        implicit none
        integer(int32),intent(in) :: coord_d0
        integer(int32),intent(in) :: coord_d1
        integer(int32),intent(out) :: exitstatus
        real(dp),dimension(coord_d0,coord_d1),intent(out) :: coord
        real(dp),intent(in) :: lat
        real(dp),intent(in) :: lon
        real(dp),intent(in) :: theta0
        real(dp),intent(in) :: cinterval
        integer(int32),intent(out) :: cnum
        call MakeCircleCoord(coord,lat,lon,theta0,cinterval=cinterval,&
                             cnum=cnum,exitstatus=exitstatus)
    end subroutine pyMakeCircleCoord

    subroutine pyMakeEllipseCoord(exitstatus,coord,lat,lon,dec,A_theta,&
                                  B_theta,cinterval,cnum,coord_d0,coord_d1)
        use shtools, only: MakeEllipseCoord
        use ftypes
        implicit none
        integer(int32),intent(in) :: coord_d0
        integer(int32),intent(in) :: coord_d1
        integer(int32),intent(out) :: exitstatus
        real(dp),dimension(coord_d0,coord_d1),intent(out) :: coord
        real(dp),intent(in) :: lat
        real(dp),intent(in) :: lon
        real(dp),intent(in) :: dec
        real(dp),intent(in) :: A_theta
        real(dp),intent(in) :: B_theta
        real(dp),intent(in) :: cinterval
        integer(int32),intent(out) :: cnum
        call MakeEllipseCoord(coord,lat,lon,dec,A_theta,B_theta,&
                                cinterval=cinterval,cnum=cnum,&
                                exitstatus=exitstatus)
    end subroutine pyMakeEllipseCoord

    subroutine pyWigner3j(exitstatus,w3j,jmin,jmax,j2,j3,m1,m2,m3,w3j_d0)
        use shtools, only: Wigner3j
        use ftypes
        implicit none
        integer(int32),intent(in) :: w3j_d0
        integer(int32),intent(out) :: exitstatus
        real(dp),dimension(w3j_d0),intent(out) :: w3j
        integer(int32),intent(out) :: jmin
        integer(int32),intent(out) :: jmax
        integer(int32),intent(in) :: j2
        integer(int32),intent(in) :: j3
        integer(int32),intent(in) :: m1
        integer(int32),intent(in) :: m2
        integer(int32),intent(in) :: m3
        call Wigner3j(w3j,jmin,jmax,j2,j3,m1,m2,m3,exitstatus=exitstatus)
    end subroutine pyWigner3j

    subroutine pyDHaj(exitstatus,n,aj,extend,aj_d0)
        use shtools, only: DHaj
        use ftypes
        implicit none
        integer(int32),intent(in) :: aj_d0
        integer(int32),intent(out) :: exitstatus
        real(dp),dimension(aj_d0),intent(out) :: aj
        integer(int32),intent(in) :: n
        integer(int32),intent(in) :: extend
        call DHaj(n,aj,extend=extend,exitstatus=exitstatus)
    end subroutine pyDHaj

    subroutine pySHRotateTapers(exitstatus,tapersrot,tapers,taper_order,&
                                lmax,nrot,x,dj,tapersrot_d0,tapersrot_d1,&
                                tapers_d0,tapers_d1,taper_order_d0,&
                                dj_d0,dj_d1,dj_d2)
        use shtools, only: SHRotateTapers
        use ftypes
        implicit none
        integer(int32),intent(in) :: tapersrot_d0
        integer(int32),intent(in) :: tapersrot_d1
        integer(int32),intent(in) :: tapers_d0
        integer(int32),intent(in) :: tapers_d1
        integer(int32),intent(in) :: taper_order_d0
        integer(int32),intent(in) :: dj_d0
        integer(int32),intent(in) :: dj_d1
        integer(int32),intent(in) :: dj_d2
        integer(int32),intent(out) :: exitstatus
        real(dp),dimension(tapersrot_d0,tapersrot_d1),intent(out) :: tapersrot
        real(dp),dimension(tapers_d0,tapers_d1),intent(in) :: tapers
        integer(int32),dimension(taper_order_d0),intent(in) :: taper_order
        integer(int32),intent(in) :: lmax
        integer(int32),intent(in) :: nrot
        real(dp),intent(in) :: x(3)
        real(dp),dimension(dj_d0,dj_d1,dj_d2),intent(in) ::dj
        call SHRotateTapers(tapersrot,tapers,taper_order,lmax,nrot,&
                            x,dj,exitstatus=exitstatus)
    end subroutine pySHRotateTapers

    subroutine pySlepianCoeffs(exitstatus,falpha,galpha,flm,lmax,nmax,&
                               falpha_d0,galpha_d0,galpha_d1,flm_d0,flm_d1,&
                               flm_d2)
        use shtools, only: SlepianCoeffs
        use ftypes
        implicit none
        integer(int32),intent(in) :: falpha_d0
        integer(int32),intent(in) :: galpha_d0
        integer(int32),intent(in) :: galpha_d1
        integer(int32),intent(in) :: flm_d0
        integer(int32),intent(in) :: flm_d1
        integer(int32),intent(in) :: flm_d2
        integer(int32),intent(out) :: exitstatus
        real(dp),dimension(falpha_d0),intent(out) :: falpha
        real(dp),dimension(galpha_d0,galpha_d1),intent(in) :: galpha
        real(dp),dimension(flm_d0,flm_d1,flm_d2),intent(in) :: flm
        integer(int32),intent(in) :: lmax
        integer(int32),intent(in) :: nmax
        call SlepianCoeffs(falpha,galpha,flm,lmax,nmax,&
                           exitstatus=exitstatus)
    end subroutine pySlepianCoeffs

    subroutine pySlepianCoeffsToSH(exitstatus,flm,falpha,galpha,lmax,nmax,&
                                   flm_d0,flm_d1,flm_d2,falpha_d0,galpha_d0,&
                                   galpha_d1)
        use shtools, only: SlepianCoeffsToSH
        use ftypes
        implicit none
        integer(int32),intent(in) :: flm_d0
        integer(int32),intent(in) :: flm_d1
        integer(int32),intent(in) :: flm_d2
        integer(int32),intent(in) :: falpha_d0
        integer(int32),intent(in) :: galpha_d0
        integer(int32),intent(in) :: galpha_d1
        integer(int32),intent(out) :: exitstatus
        real(dp),dimension(flm_d0,flm_d1,flm_d2),intent(out) :: flm
        real(dp),dimension(galpha_d0,galpha_d1),intent(in) :: galpha
        real(dp),dimension(falpha_d0),intent(in) :: falpha
        integer(int32),intent(in) :: lmax
        integer(int32),intent(in) :: nmax
        call SlepianCoeffsToSH(flm,falpha,galpha,lmax,nmax,&
                               exitstatus=exitstatus)
    end subroutine pySlepianCoeffsToSH

    subroutine pySHSCouplingMatrix(exitstatus,kij,galpha,lmax,nmax,kij_d0,&
                                   kij_d1,galpha_d0,galpha_d1)
        use shtools, only: SHSCouplingMatrix
        use ftypes
        implicit none
        integer(int32),intent(in) :: kij_d0
        integer(int32),intent(in) :: kij_d1
        integer(int32),intent(in) :: galpha_d0
        integer(int32),intent(in) :: galpha_d1
        integer(int32),intent(out) :: exitstatus
        real(dp),dimension(kij_d0,kij_d1),intent(out) :: kij
        real(dp),dimension(galpha_d0,galpha_d1),intent(in) :: galpha
        integer(int32),intent(in) :: lmax
        integer(int32),intent(in) :: nmax
        call SHSCouplingMatrix(kij,galpha,lmax,nmax,exitstatus=exitstatus)
    end subroutine pySHSCouplingMatrix

    subroutine pySHSlepianVar(exitstatus,l,galpha,galpha_order,Sff,kmax,lmax,&
                              variance,galpha_order_d0,Sff_d0,galpha_d0,&
                              galpha_d1)
        use shtools, only: SHSlepianVar
        use ftypes
        implicit none
        integer(int32),intent(in) :: galpha_order_d0
        integer(int32),intent(in) :: Sff_d0
        integer(int32),intent(in) :: galpha_d0
        integer(int32),intent(in) :: galpha_d1
        integer(int32),intent(out) :: exitstatus
        integer(int32),intent(in) :: l
        real(dp),dimension(galpha_d0,galpha_d1),intent(in) :: galpha
        integer(int32),dimension(galpha_order_d0),intent(in) :: galpha_order
        real(dp),dimension(Sff_d0),intent(in) :: Sff
        integer(int32),intent(in) :: kmax
        integer(int32),intent(in) :: lmax
        real(dp),intent(out) :: variance
        call SHSlepianVar(l,galpha,galpha_order,lmax,kmax,Sff,variance,&
                          exitstatus=exitstatus)
    end subroutine pySHSlepianVar

    subroutine pySHSCouplingMatrixCap(exitstatus,kij,galpha,galpha_order,lmax,&
                                      nmax,kij_d0,kij_d1,galpha_d0,galpha_d1,&
                                      galpha_order_d0)
        use shtools, only: SHSCouplingMatrixCap
        use ftypes
        implicit none
        integer(int32),intent(in) :: kij_d0
        integer(int32),intent(in) :: kij_d1
        integer(int32),intent(in) :: galpha_d0
        integer(int32),intent(in) :: galpha_d1
        integer(int32),intent(in) :: galpha_order_d0
        integer(int32),intent(out) :: exitstatus
        real(dp),dimension(kij_d0,kij_d1),intent(out) :: kij
        real(dp),dimension(galpha_d0,galpha_d1),intent(in) :: galpha
        integer(int32),dimension(galpha_order_d0),intent(in) :: galpha_order
        integer(int32),intent(in) :: lmax
        integer(int32),intent(in) :: nmax
        call SHSCouplingMatrixCap(kij,galpha,galpha_order,lmax,nmax,&
                                  exitstatus=exitstatus)
    end subroutine pySHSCouplingMatrixCap

    subroutine pyMakeGradientDH(exitstatus,cilm,lmax,theta,phi,n,sampling,&
                                lmax_calc,extend,radius,phi_d0,phi_d1,cilm_d0,&
                                cilm_d1,cilm_d2,theta_d0,theta_d1)
        use shtools, only: MakeGradientDH
        use ftypes
        implicit none
        integer(int32),intent(in) :: phi_d0
        integer(int32),intent(in) :: phi_d1
        integer(int32),intent(in) :: cilm_d0
        integer(int32),intent(in) :: cilm_d1
        integer(int32),intent(in) :: cilm_d2
        integer(int32),intent(in) :: theta_d0
        integer(int32),intent(in) :: theta_d1
        integer(int32),intent(out) :: exitstatus
        real(dp),dimension(cilm_d0,cilm_d1,cilm_d2),intent(in) :: cilm
        integer(int32),intent(in) :: lmax
        real(dp),dimension(theta_d0,theta_d1),intent(out) :: theta
        real(dp),dimension(phi_d0,phi_d1),intent(out) :: phi
        integer(int32),intent(out) :: n
        integer(int32),intent(in) :: sampling
        integer(int32),intent(in) :: lmax_calc
        integer(int32),intent(in) :: extend
        real(dp),intent(in) :: radius
        call MakeGradientDH(cilm,lmax,theta,phi,n,sampling=sampling,&
                            lmax_calc=lmax_calc,extend=extend,radius=radius,&
                            exitstatus=exitstatus)
    end subroutine pyMakeGradientDH
