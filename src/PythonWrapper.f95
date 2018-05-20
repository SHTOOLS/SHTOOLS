    subroutine pyPlmBar(exitstatus,p,lmax,z,csphase,cnorm,p_d0)
        use shtools, only: PlmBar
        implicit none
        integer, intent(out) :: exitstatus
        real*8, dimension(p_d0),intent(out) :: p
        integer, intent(in) :: lmax
        real*8, intent(in) :: z
        integer, optional,intent(in) :: csphase
        integer, optional,intent(in) :: cnorm
        integer, intent(in) :: p_d0
        call PlmBar(p,lmax,z,csphase=csphase,cnorm=cnorm,exitstatus=exitstatus)
    end subroutine pyPlmBar

    subroutine pyPlmBar_d1(exitstatus,p,dp,lmax,z,csphase,cnorm,p_d0,dp_d0)
        use shtools, only: PlmBar_d1
        implicit none
        integer, intent(out) :: exitstatus
        real*8, dimension(p_d0),intent(out) :: p
        real*8, dimension(dp_d0),intent(out) :: dp
        integer, intent(in) :: lmax
        real*8, intent(in) :: z
        integer, optional,intent(in) :: csphase
        integer, optional,intent(in) :: cnorm
        integer, intent(in) :: p_d0
        integer, intent(in) :: dp_d0
        call PlmBar_d1(p,dp,lmax,z,csphase=csphase,cnorm=cnorm,exitstatus=exitstatus)
    end subroutine pyPlmBar_d1

    subroutine pyPlBar(exitstatus,p,lmax,z,p_d0)
        use shtools, only: PlBar
        implicit none
        integer, intent(out) :: exitstatus
        real*8, dimension(p_d0),intent(out) :: p
        integer, intent(in) :: lmax
        real*8, intent(in) :: z
        integer, intent(in) :: p_d0
        call PlBar(p,lmax,z,exitstatus=exitstatus)
    end subroutine pyPlBar

    subroutine pyPlBar_d1(exitstatus,p,dp,lmax,z,p_d0,dp_d0)
        use shtools, only: PlBar_d1
        implicit none
        integer, intent(out) :: exitstatus
        real*8, dimension(p_d0),intent(out) :: p
        real*8, dimension(dp_d0),intent(out) :: dp
        integer, intent(in) :: lmax
        real*8, intent(in) :: z
        integer, intent(in) :: p_d0
        integer, intent(in) :: dp_d0
        call PlBar_d1(p,dp,lmax,z,exitstatus=exitstatus)
    end subroutine pyPlBar_d1

    subroutine pyPlmON(exitstatus,p,lmax,z,csphase,cnorm,p_d0)
        use shtools, only: PlmON
        implicit none
        integer, intent(out) :: exitstatus
        real*8, dimension(p_d0),intent(out) :: p
        integer, intent(in) :: lmax
        real*8, intent(in) :: z
        integer, optional,intent(in) :: csphase
        integer, optional,intent(in) :: cnorm
        integer, intent(in) :: p_d0
        call PlmON(p,lmax,z,csphase=csphase,cnorm=cnorm,exitstatus=exitstatus)
    end subroutine pyPlmON

    subroutine pyPlmON_d1(exitstatus,p,dp,lmax,z,csphase,cnorm,p_d0,dp_d0)
        use shtools, only: PlmON_d1
        implicit none
        integer, intent(out) :: exitstatus
        real*8, dimension(p_d0),intent(out) :: p
        real*8, dimension(dp_d0),intent(out) :: dp
        integer, intent(in) :: lmax
        real*8, intent(in) :: z
        integer, optional,intent(in) :: csphase
        integer, optional,intent(in) :: cnorm
        integer, intent(in) :: p_d0
        integer, intent(in) :: dp_d0
        call PlmON_d1(p,dp,lmax,z,csphase=csphase,cnorm=cnorm,exitstatus=exitstatus)
    end subroutine pyPlmON_d1

    subroutine pyPlON(exitstatus,p,lmax,z,p_d0)
        use shtools, only: PlON
        implicit none
        integer, intent(out) :: exitstatus
        real*8, dimension(p_d0),intent(out) :: p
        integer, intent(in) :: lmax
        real*8, intent(in) :: z
        integer, intent(in) :: p_d0
        call PlON(p,lmax,z,exitstatus=exitstatus)
    end subroutine pyPlON

    subroutine pyPlON_d1(exitstatus,p,dp,lmax,z,p_d0,dp_d0)
        use shtools, only: PlON_d1
        implicit none
        integer, intent(out) :: exitstatus
        real*8, dimension(p_d0),intent(out) :: p
        real*8, dimension(dp_d0),intent(out) :: dp
        integer, intent(in) :: lmax
        real*8, intent(in) :: z
        integer, intent(in) :: p_d0
        integer, intent(in) :: dp_d0
        call PlON_d1(p,dp,lmax,z,exitstatus=exitstatus)
    end subroutine pyPlON_d1

    subroutine pyPlmSchmidt(exitstatus,p,lmax,z,csphase,cnorm,p_d0)
        use shtools, only: PlmSchmidt
        implicit none
        integer, intent(out) :: exitstatus
        real*8, dimension(p_d0),intent(out) :: p
        integer, intent(in) :: lmax
        real*8, intent(in) :: z
        integer, optional,intent(in) :: csphase
        integer, optional,intent(in) :: cnorm
        integer, intent(in) :: p_d0
        call PlmSchmidt(p,lmax,z,csphase=csphase,cnorm=cnorm,exitstatus=exitstatus)
    end subroutine pyPlmSchmidt

    subroutine pyPlmSchmidt_d1(exitstatus,p,dp,lmax,z,csphase,cnorm,p_d0,dp_d0)
        use shtools, only: PlmSchmidt_d1
        implicit none
        integer, intent(out) :: exitstatus
        real*8, dimension(p_d0),intent(out) :: p
        real*8, dimension(dp_d0),intent(out) :: dp
        integer, intent(in) :: lmax
        real*8, intent(in) :: z
        integer, optional,intent(in) :: csphase
        integer, optional,intent(in) :: cnorm
        integer, intent(in) :: p_d0
        integer, intent(in) :: dp_d0
        call PlmSchmidt_d1(p,dp,lmax,z,csphase=csphase,cnorm=cnorm,exitstatus=exitstatus)
    end subroutine pyPlmSchmidt_d1

    subroutine pyPlSchmidt(exitstatus,p,lmax,z,p_d0)
        use shtools, only: PlSchmidt
        implicit none
        integer, intent(out) :: exitstatus
        real*8, dimension(p_d0),intent(out) :: p
        integer, intent(in) :: lmax
        real*8, intent(in) :: z
        integer, intent(in) :: p_d0
        call PlSchmidt(p,lmax,z,exitstatus=exitstatus)
    end subroutine pyPlSchmidt

    subroutine pyPlSchmidt_d1(exitstatus,p,dp,lmax,z,p_d0,dp_d0)
        use shtools, only: PlSchmidt_d1
        implicit none
        integer, intent(out) :: exitstatus
        real*8, dimension(p_d0),intent(out) :: p
        real*8, dimension(dp_d0),intent(out) :: dp
        integer, intent(in) :: lmax
        real*8, intent(in) :: z
        integer, intent(in) :: p_d0
        integer, intent(in) :: dp_d0
        call PlSchmidt_d1(p,dp,lmax,z,exitstatus=exitstatus)
    end subroutine pyPlSchmidt_d1

    subroutine pyPLegendreA(exitstatus,p,lmax,z,csphase,p_d0)
        use shtools, only: PLegendreA
        implicit none
        integer, intent(out) :: exitstatus
        real*8, dimension(p_d0),intent(out) :: p
        integer, intent(in) :: lmax
        real*8, intent(in) :: z
        integer, optional,intent(in) :: csphase
        integer, intent(in) :: p_d0
        call PLegendreA(p,lmax,z,csphase=csphase,exitstatus=exitstatus)
    end subroutine pyPLegendreA

    subroutine pyPLegendreA_d1(exitstatus,p,dp,lmax,z,csphase,p_d0,dp_d0)
        use shtools, only: PLegendreA_d1
        implicit none
        integer, intent(out) :: exitstatus
        real*8, dimension(p_d0),intent(out) :: p
        real*8, dimension(dp_d0),intent(out) :: dp
        integer, intent(in) :: lmax
        real*8, intent(in) :: z
        integer, optional,intent(in) :: csphase
        integer, intent(in) :: p_d0
        integer, intent(in) :: dp_d0
        call PLegendreA_d1(p,dp,lmax,z,csphase=csphase,exitstatus=exitstatus)
    end subroutine pyPLegendreA_d1

    subroutine pyPLegendre(exitstatus,p,lmax,z,p_d0)
        use shtools, only: PLegendre
        implicit none
        integer, intent(out) :: exitstatus
        real*8, dimension(p_d0),intent(out) :: p
        integer, intent(in) :: lmax
        real*8, intent(in) :: z
        integer, intent(in) :: p_d0
        call PLegendre(p,lmax,z,exitstatus=exitstatus)
    end subroutine pyPLegendre

    subroutine pyPLegendre_d1(exitstatus,p,dp,lmax,z,p_d0,dp_d0)
        use shtools, only: PLegendre_d1
        implicit none
        integer, intent(out) :: exitstatus
        real*8, dimension(p_d0),intent(out) :: p
        real*8, dimension(dp_d0),intent(out) :: dp
        integer, intent(in) :: lmax
        real*8, intent(in) :: z
        integer, intent(in) :: p_d0
        integer, intent(in) :: dp_d0
        call PLegendre_d1(p,dp,lmax,z,exitstatus=exitstatus)
    end subroutine pyPLegendre_d1

    subroutine pySHExpandDH(exitstatus,grid,n,cilm,lmax,norm,sampling,csphase,&
                            lmax_calc, cilm_d0,cilm_d1,cilm_d2,grid_d0,grid_d1)
        use shtools, only: SHExpandDH
        implicit none
        integer, intent(out) :: exitstatus
        real*8, dimension(grid_d0,grid_d1),intent(in) :: grid
        integer, intent(in) :: n
        real*8, dimension(cilm_d0,cilm_d1,cilm_d2),intent(out) :: cilm
        integer, intent(out) :: lmax
        integer, optional,intent(in) :: norm
        integer, optional,intent(in) :: sampling
        integer, optional,intent(in) :: csphase
        integer, optional,intent(in) :: lmax_calc
        integer, intent(in) :: cilm_d0
        integer, intent(in) :: cilm_d1
        integer, intent(in) :: cilm_d2
        integer, intent(in) :: grid_d0
        integer, intent(in) :: grid_d1
        call SHExpandDH(grid,n,cilm,lmax,norm=norm,sampling=sampling, &
                        csphase=csphase,lmax_calc=lmax_calc,&
                        exitstatus=exitstatus)
    end subroutine pySHExpandDH

    subroutine pyMakeGridDH(exitstatus,griddh,n,cilm,lmax,norm,sampling,&
                            csphase,lmax_calc,cilm_d0,cilm_d1,cilm_d2,&
                            griddh_d0,griddh_d1)
        use shtools, only: MakeGridDH
        implicit none
        integer, intent(out) :: exitstatus
        real*8, dimension(griddh_d0,griddh_d1),intent(out) :: griddh
        integer, intent(out) :: n
        real*8, dimension(cilm_d0,cilm_d1,cilm_d2),intent(in) :: cilm
        integer, intent(in) :: lmax
        integer, optional,intent(in) :: norm
        integer, optional,intent(in) :: sampling
        integer, optional,intent(in) :: csphase
        integer, optional,intent(in) :: lmax_calc
        integer, intent(in) :: cilm_d0
        integer, intent(in) :: cilm_d1
        integer, intent(in) :: cilm_d2
        integer, intent(in) :: griddh_d0
        integer, intent(in) :: griddh_d1
        call MakeGridDH(griddh,n,cilm,lmax,norm=norm,sampling=sampling, &
                        csphase=csphase,lmax_calc=lmax_calc,&
                        exitstatus=exitstatus)
    end subroutine pyMakeGridDH

    subroutine pySHExpandDHC(exitstatus,grid,n,cilm,lmax,norm,sampling,csphase,&
                             lmax_calc,cilm_d0,cilm_d1,cilm_d2,grid_d0,grid_d1)
        use shtools, only: SHExpandDHC
        implicit none
        integer, intent(out) :: exitstatus
        complex*16, dimension(grid_d0,grid_d1),intent(in) :: grid
        integer, intent(in) :: n
        complex*16, dimension(cilm_d0,cilm_d1,cilm_d2),intent(out) :: cilm
        integer, intent(out) :: lmax
        integer, optional,intent(in) :: norm
        integer, optional,intent(in) :: sampling
        integer, optional,intent(in) :: csphase
        integer, optional,intent(in) :: lmax_calc
        integer, intent(in) :: cilm_d0
        integer, intent(in) :: cilm_d1
        integer, intent(in) :: cilm_d2
        integer, intent(in) :: grid_d0
        integer, intent(in) :: grid_d1
        call SHExpandDHC(grid,n,cilm,lmax,norm=norm,sampling=sampling, &
                         csphase=csphase,lmax_calc=lmax_calc,&
                         exitstatus=exitstatus)
    end subroutine pySHExpandDHC

    subroutine pyMakeGridDHC(exitstatus,griddh,n,cilm,lmax,norm,sampling,&
                             csphase,lmax_calc,cilm_d0,cilm_d1,cilm_d2,&
                             griddh_d0,griddh_d1)
        use shtools, only: MakeGridDHC
        implicit none
        integer, intent(out) :: exitstatus
        complex*16, dimension(griddh_d0,griddh_d1),intent(out) :: griddh
        integer, intent(out) :: n
        complex*16, dimension(cilm_d0,cilm_d1,cilm_d2),intent(in) :: cilm
        integer, intent(in) :: lmax
        integer, optional,intent(in) :: norm
        integer, optional,intent(in) :: sampling
        integer, optional,intent(in) :: csphase
        integer, optional,intent(in) :: lmax_calc
        integer, intent(in) :: cilm_d0
        integer, intent(in) :: cilm_d1
        integer, intent(in) :: cilm_d2
        integer, intent(in) :: griddh_d0
        integer, intent(in) :: griddh_d1
        call MakeGridDHC(griddh,n,cilm,lmax,norm=norm,sampling=sampling, &
                         csphase=csphase,lmax_calc=lmax_calc,&
                         exitstatus=exitstatus)
    end subroutine pyMakeGridDHC

    subroutine pyshglq(exitstatus,lmax,zero,w,zero_d0,w_d0)
        use shtools, only: SHGLQ
        implicit none
        integer, intent(out) :: exitstatus
        integer, intent(in) :: lmax
        real*8, dimension(zero_d0),intent(out) :: zero
        real*8, dimension(w_d0),intent(out) :: w
        integer, intent(in) :: zero_d0
        integer, intent(in) :: w_d0
        call SHGLQ(lmax,zero,w,exitstatus=exitstatus)
    end subroutine pySHGLQ

    subroutine pySHExpandGLQ(exitstatus,cilm,lmax,gridglq,w,zero,norm,csphase,&
                             lmax_calc,cilm_d0,cilm_d1,cilm_d2,gridglq_d0,&
                             gridglq_d1,zero_d0,w_d0)
        use shtools, only: SHExpandGLQ
        implicit none
        integer, intent(out) :: exitstatus
        real*8, dimension(cilm_d0,cilm_d1,cilm_d2),intent(out) :: cilm
        integer, intent(in) :: lmax
        real*8, dimension(gridglq_d0,gridglq_d1),intent(in) :: gridglq
        real*8, dimension(w_d0),intent(in) :: w
        real*8, optional,dimension(zero_d0),intent(in) :: zero
        integer, optional,intent(in) :: norm
        integer, optional,intent(in) :: csphase
        integer, optional,intent(in) :: lmax_calc
        integer, intent(in) :: cilm_d0
        integer, intent(in) :: cilm_d1
        integer, intent(in) :: cilm_d2
        integer, intent(in) :: gridglq_d0
        integer, intent(in) :: gridglq_d1
        integer, intent(in) :: zero_d0
        integer, intent(in) :: w_d0
        call SHExpandGLQ(cilm,lmax,gridglq,w,zero=zero,norm=norm, &
                         csphase=csphase,lmax_calc=lmax_calc,&
                         exitstatus=exitstatus)
    end subroutine pySHExpandGLQ

    subroutine pyMakeGridGLQ(exitstatus,gridglq,cilm,lmax,zero,norm,csphase,&
                             lmax_calc,gridglq_d0,gridglq_d1,cilm_d0,cilm_d1,&
                             cilm_d2,zero_d0)
        use shtools, only: MakeGridGLQ
        implicit none
        integer, intent(out) :: exitstatus
        real*8, dimension(gridglq_d0,gridglq_d1),intent(out) :: gridglq
        real*8, dimension(cilm_d0,cilm_d1,cilm_d2),intent(in) :: cilm
        integer, intent(in) :: lmax
        real*8, optional,dimension(zero_d0),intent(in) :: zero
        integer, optional,intent(in) :: norm
        integer, optional,intent(in) :: csphase
        integer, optional,intent(in) :: lmax_calc
        integer, intent(in) :: gridglq_d0
        integer, intent(in) :: gridglq_d1
        integer, intent(in) :: cilm_d0
        integer, intent(in) :: cilm_d1
        integer, intent(in) :: cilm_d2
        integer, intent(in) :: zero_d0
        call MakeGridGLQ(gridglq,cilm,lmax,zero=zero,norm=norm, &
                         csphase=csphase,lmax_calc=lmax_calc,&
                         exitstatus=exitstatus)
    end subroutine pyMakeGridGLQ

    subroutine pySHExpandGLQC(exitstatus,cilm,lmax,gridglq,w,zero,norm,csphase,&
                              lmax_calc,cilm_d0,cilm_d1,cilm_d2,gridglq_d0,&
                              gridglq_d1,zero_d0,w_d0)
        use shtools, only: SHExpandGLQC
        implicit none
        integer, intent(out) :: exitstatus
        complex*16, dimension(cilm_d0,cilm_d1,cilm_d2),intent(out) :: cilm
        integer, intent(in) :: lmax
        complex*16, dimension(gridglq_d0,gridglq_d1),intent(in) :: gridglq
        real*8, dimension(w_d0),intent(in) :: w
        real*8, optional,dimension(zero_d0),intent(in) :: zero
        integer, optional,intent(in) :: norm
        integer, optional,intent(in) :: csphase
        integer, optional,intent(in) :: lmax_calc
        integer, intent(in) :: cilm_d0
        integer, intent(in) :: cilm_d1
        integer, intent(in) :: cilm_d2
        integer, intent(in) :: gridglq_d0
        integer, intent(in) :: gridglq_d1
        integer, intent(in) :: zero_d0
        integer, intent(in) :: w_d0
        call SHExpandGLQC(cilm,lmax,gridglq,w,zero=zero,norm=norm, &
                          csphase=csphase,lmax_calc=lmax_calc,&
                          exitstatus=exitstatus)
    end subroutine pySHExpandGLQC

    subroutine pyMakeGridGLQC(exitstatus,gridglq,cilm,lmax,zero,norm,csphase,&
                              lmax_calc,gridglq_d0,gridglq_d1,cilm_d0,cilm_d1,&
                              cilm_d2,zero_d0)
        use shtools, only: MakeGridGLQC
        implicit none
        integer, intent(out) :: exitstatus
        complex*16, dimension(gridglq_d0,gridglq_d1),intent(out) :: gridglq
        complex*16, dimension(cilm_d0,cilm_d1,cilm_d2),intent(in) :: cilm
        integer, intent(in) :: lmax
        real*8, optional,dimension(zero_d0),intent(in) :: zero
        integer, optional,intent(in) :: norm
        integer, optional,intent(in) :: csphase
        integer, optional,intent(in) :: lmax_calc
        integer, intent(in) :: gridglq_d0
        integer, intent(in) :: gridglq_d1
        integer, intent(in) :: cilm_d0
        integer, intent(in) :: cilm_d1
        integer, intent(in) :: cilm_d2
        integer, intent(in) :: zero_d0
        call MakeGridGLQC(gridglq,cilm,lmax,zero=zero,norm=norm, &
                          csphase=csphase,lmax_calc=lmax_calc,&
                          exitstatus=exitstatus)
    end subroutine pyMakeGridGLQC

    subroutine pyGLQGridCoord(exitstatus,latglq,longlq,lmax,nlat,nlong,&
                              latglq_d0,longlq_d0)
        use shtools, only: GLQGridCoord
        implicit none
        integer, intent(out) :: exitstatus
        real*8, dimension(latglq_d0),intent(out) :: latglq
        real*8, dimension(longlq_d0),intent(out) :: longlq
        integer, intent(in) :: lmax
        integer, intent(out) :: nlat
        integer, intent(out) :: nlong
        integer, intent(in) :: latglq_d0
        integer, intent(in) :: longlq_d0
        call GLQGridCoord(latglq,longlq,lmax,nlat,nlong,exitstatus=exitstatus)
    end subroutine pyGLQGridCoord

    subroutine pySHExpandLSQ(exitstatus,cilm,d,lat,lon,nmax,lmax,norm,chi2,&
                             csphase,d_d0,lon_d0,cilm_d0,cilm_d1,cilm_d2,lat_d0)
        use shtools, only: SHExpandLSQ
        implicit none
        integer, intent(out) :: exitstatus
        real*8, dimension(cilm_d0,cilm_d1,cilm_d2),intent(out) :: cilm
        real*8, dimension(d_d0),intent(in) :: d
        real*8, dimension(lat_d0),intent(in) :: lat
        real*8, dimension(lon_d0),intent(in) :: lon
        integer, intent(in) :: nmax
        integer, intent(in) :: lmax
        integer, optional,intent(in) :: norm
        real*8, optional,intent(out) :: chi2
        integer, optional,intent(in) :: csphase
        integer, intent(in) :: d_d0
        integer, intent(in) :: lon_d0
        integer, intent(in) :: cilm_d0
        integer, intent(in) :: cilm_d1
        integer, intent(in) :: cilm_d2
        integer, intent(in) :: lat_d0
        call SHExpandLSQ(cilm,d,lat,lon,nmax,lmax,norm=norm,chi2=chi2, &
                         csphase=csphase,exitstatus=exitstatus)
    end subroutine pySHExpandLSQ

    subroutine pyMakeGrid2d(exitstatus,grid,cilm,lmax,interval,nlat,nlong,norm,&
                            csphase,f,a,north,south,east,west,dealloc,cilm_d0,&
                            cilm_d1,cilm_d2,grid_d0,grid_d1)
        use shtools, only: MakeGrid2d
        implicit none
        integer, intent(out) :: exitstatus
        real*8, dimension(grid_d0,grid_d1),intent(out) :: grid
        real*8, dimension(cilm_d0,cilm_d1,cilm_d2),intent(in) :: cilm
        integer, intent(in) :: lmax
        real*8, intent(in) :: interval
        integer, intent(out) :: nlat
        integer, intent(out) :: nlong
        integer, optional,intent(in) :: norm
        integer, optional,intent(in) :: csphase
        real*8, optional,intent(in) :: f
        real*8, optional,intent(in) :: a
        real*8, optional,intent(in) :: north
        real*8, optional,intent(in) :: south
        real*8, optional,intent(in) :: east
        real*8, optional,intent(in) :: west
        integer, optional,intent(in) :: dealloc
        integer, intent(in) :: cilm_d0
        integer, intent(in) :: cilm_d1
        integer, intent(in) :: cilm_d2
        integer, intent(in) :: grid_d0
        integer, intent(in) :: grid_d1
        if (f<0.0d0 .and. a<0.0d0) then
            call MakeGrid2d(grid,cilm,lmax,interval,nlat,nlong,norm=norm, &
                            csphase=csphase,north=north,south=south,east=east,&
                            west=west,dealloc=dealloc,exitstatus=exitstatus)
        else
            call MakeGrid2d(grid,cilm,lmax,interval,nlat,nlong,norm=norm, &
                            csphase=csphase,f=f,a=a,north=north, &
                            south=south,east=east,west=west,dealloc=dealloc,&
                            exitstatus=exitstatus)
        endif
    end subroutine pyMakeGrid2d

    function pyMakeGridPoint(cilm,lmax,lat,lon,norm,csphase,dealloc, &
                             cilm_d0,cilm_d1,cilm_d2)
        use shtools, only: MakeGridPoint
        implicit none
        real*8, dimension(cilm_d0,cilm_d1,cilm_d2),intent(in) :: cilm
        integer, intent(in) :: lmax
        real*8, intent(in) :: lat
        real*8, intent(in) :: lon
        integer, optional,intent(in) :: norm
        integer, optional,intent(in) :: csphase
        integer, optional,intent(in) :: dealloc
        integer, intent(in) :: cilm_d0
        integer, intent(in) :: cilm_d1
        integer, intent(in) :: cilm_d2
        real*8 :: pyMakeGridPoint
        pyMakeGridPoint=MakeGridPoint(cilm,lmax,lat,lon,norm=norm, &
                                      csphase=csphase,dealloc=dealloc)
    end function pyMakeGridPoint

    function pyMakeGridPointC(cilm,lmax,lat,lon,norm,csphase,dealloc, &
                              cilm_d0,cilm_d1,cilm_d2)
        use shtools, only: MakeGridPointC
        implicit none
        complex*16, dimension(cilm_d0,cilm_d1,cilm_d2),intent(in) :: cilm
        integer, intent(in) :: lmax
        real*8, intent(in) :: lat
        real*8, intent(in) :: lon
        integer, optional,intent(in) :: norm
        integer, optional,intent(in) :: csphase
        integer, optional,intent(in) :: dealloc
        integer, intent(in) :: cilm_d0
        integer, intent(in) :: cilm_d1
        integer, intent(in) :: cilm_d2
        complex*16 :: pyMakeGridPointC
        pyMakeGridPointC=MakeGridPointC(cilm,lmax,lat,lon,norm=norm, &
                                        csphase=csphase,dealloc=dealloc)
    end function pyMakeGridPointC

    subroutine pySHMultiply(exitstatus,shout,sh1,lmax1,sh2,lmax2,precomp,norm,&
                            csphase,sh1_d0,sh1_d1,sh1_d2,sh2_d0,sh2_d1,sh2_d2,&
                            shout_d0,shout_d1,shout_d2)
        use shtools, only: SHMultiply
        implicit none
        integer, intent(out) :: exitstatus
        real*8, dimension(shout_d0,shout_d1,shout_d2),intent(out) :: shout
        real*8, dimension(sh1_d0,sh1_d1,sh1_d2),intent(in) :: sh1
        integer, intent(in) :: lmax1
        real*8, dimension(sh2_d0,sh2_d1,sh2_d2),intent(in) :: sh2
        integer, intent(in) :: lmax2
        integer, optional,intent(in) :: precomp
        integer, optional,intent(in) :: norm
        integer, optional,intent(in) :: csphase
        integer, intent(in) :: sh1_d0
        integer, intent(in) :: sh1_d1
        integer, intent(in) :: sh1_d2
        integer, intent(in) :: sh2_d0
        integer, intent(in) :: sh2_d1
        integer, intent(in) :: sh2_d2
        integer, intent(in) :: shout_d0
        integer, intent(in) :: shout_d1
        integer, intent(in) :: shout_d2
        call SHMultiply(shout,sh1,lmax1,sh2,lmax2,precomp=precomp, &
                        norm=norm,csphase=csphase,exitstatus=exitstatus)
    end subroutine pySHMultiply

    subroutine pySHRead2(exitstatus,filename,cilm,lmax,lmax_in,gm,r0_pot,dot,&
                         doystart,doyend,epoch,cilm_d0,cilm_d1,cilm_d2,dot_d0,&
                         dot_d1,dot_d2)
        use shtools, only: SHRead2
        implicit none
        integer, intent(out) :: exitstatus
        character*(*), intent(in) :: filename
        real*8, dimension(cilm_d0,cilm_d1,cilm_d2),intent(out) :: cilm
        integer, intent(out) :: lmax
        integer, intent(in) :: lmax_in
        real*8, intent(out) :: gm
        real*8, intent(out) :: r0_pot
        real*8, optional,dimension(dot_d0,dot_d1,dot_d2),intent(out) :: dot
        real*8, optional,intent(out) :: doystart
        real*8, optional,intent(out) :: doyend
        real*8, optional,intent(out) :: epoch
        integer, intent(in) :: cilm_d0
        integer, intent(in) :: cilm_d1
        integer, intent(in) :: cilm_d2
        integer, intent(in) :: dot_d0
        integer, intent(in) :: dot_d1
        integer, intent(in) :: dot_d2
        call SHRead2(filename,cilm,lmax,gm,r0_pot,dot=dot,doystart=doystart, &
                     doyend=doyend,epoch=epoch,exitstatus=exitstatus)
    end subroutine pySHRead2

    subroutine pySHRead2Error(exitstatus,filename,cilm,error,lmax,lmax_in,gm,&
                              r0_pot,dot,doystart,doyend,epoch,cilm_d0,cilm_d1,&
                              cilm_d2,error_d0,error_d1,error_d2,dot_d0,dot_d1,&
                              dot_d2)
        use shtools, only: SHRead2
        implicit none
        integer, intent(out) :: exitstatus
        character*(*), intent(in) :: filename
        real*8, dimension(cilm_d0,cilm_d1,cilm_d2),intent(out) :: cilm
        real*8, optional,dimension(error_d0,error_d1,error_d2),intent(out) ::&
                                                                          error
        integer, intent(out) :: lmax
        integer, intent(in) :: lmax_in
        real*8, intent(out) :: gm
        real*8, intent(out) :: r0_pot
        real*8, optional,dimension(dot_d0,dot_d1,dot_d2),intent(out) :: dot
        real*8, optional,intent(out) :: doystart
        real*8, optional,intent(out) :: doyend
        real*8, optional,intent(out) :: epoch
        integer, intent(in) :: cilm_d0
        integer, intent(in) :: cilm_d1
        integer, intent(in) :: cilm_d2
        integer, intent(in) :: error_d0
        integer, intent(in) :: error_d1
        integer, intent(in) :: error_d2
        integer, intent(in) :: dot_d0
        integer, intent(in) :: dot_d1
        integer, intent(in) :: dot_d2
        call SHRead2(filename,cilm,lmax,gm,r0_pot,error=error,dot=dot, &
                     doystart=doystart,doyend=doyend,epoch=epoch,&
                     exitstatus=exitstatus)
    end subroutine pySHRead2Error

    subroutine pySHReadJPL(exitstatus,filename,cilm,lmax,lmax_in,gm,&
                           formatstring,cilm_d0,cilm_d1,cilm_d2)
        use shtools, only: SHReadJPL
        implicit none
        integer, intent(out) :: exitstatus
        character*(*), intent(in) :: filename
        real*8, dimension(cilm_d0,cilm_d1,cilm_d2),intent(out) :: cilm
        integer, intent(out) :: lmax
        integer, intent(in) :: lmax_in
        real*8, optional,dimension(2),intent(out) :: gm
        character*6, optional,intent(in) :: formatstring
        integer, intent(in) :: cilm_d0
        integer, intent(in) :: cilm_d1
        integer, intent(in) :: cilm_d2
        call SHReadJPL(filename,cilm,lmax,gm=gm,formatstring=formatstring,&
                       exitstatus=exitstatus)
    end subroutine pySHReadJPL

    subroutine pySHReadJPLError(exitstatus,filename,cilm,error,lmax,lmax_in,&
                                gm,formatstring,cilm_d0,cilm_d1,cilm_d2,&
                                error_d0,error_d1,error_d2)
        use shtools, only: SHReadJPL
        implicit none
        integer, intent(out) :: exitstatus
        character*(*), intent(in) :: filename
        real*8, dimension(cilm_d0,cilm_d1,cilm_d2),intent(out) :: cilm
        real*8, optional,dimension(error_d0,error_d1,error_d2),intent(out) ::&
                                                                          error
        integer, intent(in) :: lmax
        integer, intent(in) :: lmax_in
        real*8, optional,dimension(2),intent(out) :: gm
        character*6, optional,intent(in) :: formatstring
        integer, intent(in) :: cilm_d0
        integer, intent(in) :: cilm_d1
        integer, intent(in) :: cilm_d2
        integer, intent(in) :: error_d0
        integer, intent(in) :: error_d1
        integer, intent(in) :: error_d2
        call SHReadJPL(filename,cilm,lmax,error=error,gm=gm, &
                       formatstring=formatstring,exitstatus=exitstatus)
    end subroutine pySHReadJPLError

    subroutine pySHCilmToVector(exitstatus,cilm,vector,lmax,vector_d0,cilm_d0,&
                                cilm_d1,cilm_d2)
        use shtools, only: SHCilmToVector
        implicit none
        integer, intent(out) :: exitstatus
        real*8, dimension(cilm_d0,cilm_d1,cilm_d2),intent(in) :: cilm
        real*8, dimension(vector_d0),intent(out) :: vector
        integer, intent(in) :: lmax
        integer, intent(in) :: vector_d0
        integer, intent(in) :: cilm_d0
        integer, intent(in) :: cilm_d1
        integer, intent(in) :: cilm_d2
        call SHCilmToVector(cilm,vector,lmax,exitstatus=exitstatus)
    end subroutine pySHCilmToVector

    subroutine pySHVectorToCilm(exitstatus,vector,cilm,lmax,vector_d0,cilm_d0,&
                                cilm_d1,cilm_d2)
        use shtools, only: SHVectorToCilm
        implicit none
        integer, intent(out) :: exitstatus
        real*8, dimension(vector_d0),intent(in) :: vector
        real*8, dimension(cilm_d0,cilm_d1,cilm_d2),intent(out) :: cilm
        integer, intent(in) :: lmax
        integer, intent(in) :: vector_d0
        integer, intent(in) :: cilm_d0
        integer, intent(in) :: cilm_d1
        integer, intent(in) :: cilm_d2
        call SHVectorToCilm(vector,cilm,lmax,exitstatus=exitstatus)
    end subroutine pySHVectorToCilm

    subroutine pySHCilmToCindex(exitstatus,cilm,cindex,degmax,cindex_d0,&
                                cindex_d1,cilm_d0,cilm_d1,cilm_d2)
        use shtools, only: SHCilmToCindex
        implicit none
        integer, intent(out) :: exitstatus
        real*8, dimension(cilm_d0,cilm_d1,cilm_d2),intent(in) :: cilm
        real*8, dimension(cindex_d0,cindex_d1),intent(out) :: cindex
        integer, optional,intent(in) :: degmax
        integer, intent(in) :: cindex_d0
        integer, intent(in) :: cindex_d1
        integer, intent(in) :: cilm_d0
        integer, intent(in) :: cilm_d1
        integer, intent(in) :: cilm_d2
        call SHCilmToCindex(cilm,cindex,degmax=degmax,exitstatus=exitstatus)
    end subroutine pySHCilmToCindex

    subroutine pySHCindexToCilm(exitstatus,cindex,cilm,degmax,cindex_d0,&
                                cindex_d1,cilm_d0,cilm_d1,cilm_d2)
        use shtools, only: SHCindexToCilm
        implicit none
        integer, intent(out) :: exitstatus
        real*8, dimension(cindex_d0,cindex_d1),intent(in) :: cindex
        real*8, dimension(cilm_d0,cilm_d1,cilm_d2),intent(out) :: cilm
        integer, optional,intent(in) :: degmax
        integer, intent(in) :: cindex_d0
        integer, intent(in) :: cindex_d1
        integer, intent(in) :: cilm_d0
        integer, intent(in) :: cilm_d1
        integer, intent(in) :: cilm_d2
        call SHCindexToCilm(cindex,cilm,degmax=degmax,exitstatus=exitstatus)
    end subroutine pySHCindexToCilm

    subroutine pySHrtoc(exitstatus,rcilm,ccilm,degmax,convention,switchcs,&
                        rcilm_d0,rcilm_d1,rcilm_d2,ccilm_d0,ccilm_d1,ccilm_d2)
        use shtools, only: SHrtoc
        implicit none
        integer, intent(out) :: exitstatus
        real*8, dimension(rcilm_d0,rcilm_d1,rcilm_d2),intent(in) :: rcilm
        real*8, dimension(ccilm_d0,ccilm_d1,ccilm_d2),intent(out) :: ccilm
        integer, optional,intent(in) :: degmax
        integer, optional,intent(in) :: convention
        integer, optional,intent(in) :: switchcs
        integer, intent(in) :: rcilm_d0
        integer, intent(in) :: rcilm_d1
        integer, intent(in) :: rcilm_d2
        integer, intent(in) :: ccilm_d0
        integer, intent(in) :: ccilm_d1
        integer, intent(in) :: ccilm_d2
        call SHrtoc(rcilm,ccilm,degmax=degmax,convention=convention, &
                    switchcs=switchcs,exitstatus=exitstatus)
    end subroutine pySHrtoc

    subroutine pySHctor(exitstatus,ccilm,rcilm,degmax,convention,switchcs,&
                        rcilm_d0,rcilm_d1,rcilm_d2,ccilm_d0,ccilm_d1,ccilm_d2)
        use shtools, only: SHctor
        implicit none
        integer, intent(out) :: exitstatus
        real*8, dimension(ccilm_d0,ccilm_d1,ccilm_d2),intent(in) :: ccilm
        real*8, dimension(rcilm_d0,rcilm_d1,rcilm_d2),intent(out) :: rcilm
        integer, optional,intent(in) :: degmax
        integer, optional,intent(in) :: convention
        integer, optional,intent(in) :: switchcs
        integer, intent(in) :: rcilm_d0
        integer, intent(in) :: rcilm_d1
        integer, intent(in) :: rcilm_d2
        integer, intent(in) :: ccilm_d0
        integer, intent(in) :: ccilm_d1
        integer, intent(in) :: ccilm_d2
        call SHctor(ccilm,rcilm,degmax=degmax,convention=convention, &
                    switchcs=switchcs,exitstatus=exitstatus)
    end subroutine pySHctor

    subroutine pydjpi2(exitstatus,dj,lmax,dj_d0,dj_d1,dj_d2)
        use shtools, only: djpi2
        implicit none
        integer, intent(out) :: exitstatus
        real*8, dimension(dj_d0,dj_d1,dj_d2),intent(out) :: dj
        integer, intent(in) :: lmax
        integer, intent(in) :: dj_d0
        integer, intent(in) :: dj_d1
        integer, intent(in) :: dj_d2
        call djpi2(dj,lmax,exitstatus=exitstatus)
    end subroutine pydjpi2

    subroutine pySHRotateCoef(exitstatus,x,cof,rcof,dj,lmax,rcof_d0,rcof_d1,&
                              dj_d0,dj_d1,dj_d2,cof_d0,cof_d1)
        use shtools, only: SHRotateCoef
        implicit none
        integer, intent(out) :: exitstatus
        real*8, dimension(3),intent(in) :: x
        real*8, dimension(cof_d0,cof_d1),intent(in) :: cof
        real*8, dimension(rcof_d0,rcof_d1),intent(out) :: rcof
        real*8, dimension(dj_d0,dj_d1,dj_d2),intent(in) :: dj
        integer, intent(in) :: lmax
        integer, intent(in) :: rcof_d0
        integer, intent(in) :: rcof_d1
        integer, intent(in) :: dj_d0
        integer, intent(in) :: dj_d1
        integer, intent(in) :: dj_d2
        integer, intent(in) :: cof_d0
        integer, intent(in) :: cof_d1
        call SHRotateCoef(x,cof,rcof,dj,lmax,exitstatus=exitstatus)
    end subroutine pySHRotateCoef

    subroutine pySHRotateRealCoef(exitstatus,cilmrot,cilm,lmax,x,dj,x_d0,dj_d0,&
                                  dj_d1,dj_d2,cilm_d0,cilm_d1,cilm_d2, &
                                  cilmrot_d0,cilmrot_d1,cilmrot_d2)
        use shtools, only: SHRotateRealCoef
        implicit none
        integer, intent(out) :: exitstatus
        real*8, dimension(cilmrot_d0,cilmrot_d1,cilmrot_d2),intent(out) ::&
                                                                        cilmrot
        real*8, dimension(cilm_d0,cilm_d1,cilm_d2),intent(in) :: cilm
        integer, intent(in) :: lmax
        real*8, dimension(x_d0),intent(in) :: x
        real*8, dimension(dj_d0,dj_d1,dj_d2),intent(in) :: dj
        integer, intent(in) :: x_d0
        integer, intent(in) :: dj_d0
        integer, intent(in) :: dj_d1
        integer, intent(in) :: dj_d2
        integer, intent(in) :: cilm_d0
        integer, intent(in) :: cilm_d1
        integer, intent(in) :: cilm_d2
        integer, intent(in) :: cilmrot_d0
        integer, intent(in) :: cilmrot_d1
        integer, intent(in) :: cilmrot_d2
        call SHRotateRealCoef(cilmrot,cilm,lmax,x,dj,exitstatus=exitstatus)
    end subroutine pySHRotateRealCoef

    subroutine pySHAdmitCorr(exitstatus,G,T,lmax,admit,admit_error,corr,G_d0,&
                             G_d1,G_d2,admit_d0,admit_error_d0,T_d0,T_d1,T_d2,&
                             corr_d0)
        use shtools, only: SHAdmitCorr
        implicit none
        integer, intent(out) :: exitstatus
        real*8, dimension(G_d0,G_d1,G_d2),intent(in) :: G
        real*8, dimension(T_d0,T_d1,T_d2),intent(in) :: T
        integer, intent(in) :: lmax
        real*8, dimension(admit_d0),intent(out) :: admit
        real*8, optional,dimension(admit_error_d0),intent(out) :: admit_error
        real*8, dimension(corr_d0),intent(out) :: corr
        integer, intent(in) :: G_d0
        integer, intent(in) :: G_d1
        integer, intent(in) :: G_d2
        integer, intent(in) :: admit_d0
        integer, intent(in) :: admit_error_d0
        integer, intent(in) :: T_d0
        integer, intent(in) :: T_d1
        integer, intent(in) :: T_d2
        integer, intent(in) :: corr_d0
        call SHAdmitCorr(G,T,lmax,admit,corr,admit_error=admit_error,&
                         exitstatus=exitstatus)
    end subroutine pySHAdmitCorr

    function pySHConfidence(l_conf,r)
        use shtools, only: SHConfidence
        implicit none
        integer, intent(in) :: l_conf
        real*8, intent(in) :: r
        real*8 :: pySHConfidence
        pySHConfidence=SHConfidence(l_conf,r)
    end function pySHConfidence

    subroutine pySHMultiTaperSE(exitstatus,mtse,sd,sh,lmax,tapers,taper_order,&
                                lmaxt,k,lat,lon,taper_wt,norm,csphase, &
                                taper_order_d0,taper_wt_d0,sh_d0,sh_d1,sh_d2, &
                                tapers_d0,tapers_d1,mtse_d0,sd_d0)
        use shtools, only: SHMultiTaperSE
        implicit none
        integer, intent(out) :: exitstatus
        real*8, dimension(mtse_d0),intent(out) :: mtse
        real*8, dimension(sd_d0),intent(out) :: sd
        real*8, dimension(sh_d0,sh_d1,sh_d2),intent(in) :: sh
        integer, intent(in) :: lmax
        real*8, dimension(tapers_d0,tapers_d1),intent(in) :: tapers
        integer, dimension(taper_order_d0),intent(in) :: taper_order
        integer, intent(in) :: lmaxt
        integer, intent(in) :: k
        real*8, optional,intent(in) :: lat
        real*8, optional,intent(in) :: lon
        real*8, optional,dimension(taper_wt_d0),intent(in) :: taper_wt
        integer, optional,intent(in) :: norm
        integer, optional,intent(in) :: csphase
        integer, intent(in) :: taper_order_d0
        integer, intent(in) :: taper_wt_d0
        integer, intent(in) :: sh_d0
        integer, intent(in) :: sh_d1
        integer, intent(in) :: sh_d2
        integer, intent(in) :: tapers_d0
        integer, intent(in) :: tapers_d1
        integer, intent(in) :: mtse_d0
        integer, intent(in) :: sd_d0
        if (taper_wt(1) < 0.d0) then
            call SHMultiTaperSE(mtse,sd,sh,lmax,tapers,taper_order, &
                                lmaxt,k,lat=lat,lon=lon,norm=norm, &
                                csphase=csphase,exitstatus=exitstatus)
        else
            call SHMultiTaperSE(mtse,sd,sh,lmax,tapers,taper_order,lmaxt, &
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
        implicit none
        integer, intent(out) :: exitstatus
        real*8, dimension(mtse_d0),intent(out) :: mtse
        real*8, dimension(sd_d0),intent(out) :: sd
        real*8, dimension(sh1_d0,sh1_d1,sh1_d2),intent(in) :: sh1
        integer, intent(in) :: lmax1
        real*8, dimension(sh2_d0,sh2_d1,sh2_d2),intent(in) :: sh2
        integer, intent(in) :: lmax2
        real*8, dimension(tapers_d0,tapers_d1),intent(in) :: tapers
        integer, dimension(taper_order_d0),intent(in) :: taper_order
        integer, intent(in) :: lmaxt
        integer, intent(in) :: k
        real*8, optional,intent(in) :: lat
        real*8, optional,intent(in) :: lon
        real*8, optional,dimension(taper_wt_d0),intent(in) :: taper_wt
        integer, optional,intent(in) :: norm
        integer, optional,intent(in) :: csphase
        integer, intent(in) :: sh1_d0
        integer, intent(in) :: sh1_d1
        integer, intent(in) :: sh1_d2
        integer, intent(in) :: sh2_d0
        integer, intent(in) :: sh2_d1
        integer, intent(in) :: sh2_d2
        integer, intent(in) :: taper_order_d0
        integer, intent(in) :: taper_wt_d0
        integer, intent(in) :: tapers_d0
        integer, intent(in) :: tapers_d1
        integer, intent(in) :: sd_d0
        integer, intent(in) :: mtse_d0
        if(taper_wt(1) < 0.d0) then
            call SHMultiTaperCSE(mtse,sd,sh1,lmax1,sh2,lmax2,tapers, &
                                 taper_order,lmaxt,k,lat=lat,lon=lon,&
                                 norm=norm,csphase=csphase, &
                                 exitstatus=exitstatus)
        else
            call SHMultiTaperCSE(mtse,sd,sh1,lmax1,sh2,lmax2,tapers, &
                                 taper_order,lmaxt,k,lat=lat,lon=lon,&
                                 taper_wt=taper_wt,norm=norm,csphase=csphase,&
                                 exitstatus=exitstatus)
        endif
    end subroutine pySHMultiTaperCSE

    subroutine pySHLocalizedAdmitCorr(exitstatus,g,t,tapers,taper_order,k,lat,&
                                      lon,lwin,lmax,admit,corr,admit_error,&
                                      corr_error,taper_wt,mtdef,k1linsig,&
                                      taper_order_d0,g_d0,g_d1,g_d2,&
                                      taper_wt_d0,corr_error_d0,admit_d0,&
                                      admit_error_d0,corr_d0,tapers_d0,&
                                      tapers_d1,t_d0,t_d1,t_d2)
        use shtools, only: SHLocalizedAdmitCorr
        implicit none
        integer, intent(out) :: exitstatus
        real*8, dimension(g_d0,g_d1,g_d2),intent(in) :: g
        real*8, dimension(t_d0,t_d1,t_d2),intent(in) :: t
        real*8, dimension(tapers_d0,tapers_d1),intent(in) :: tapers
        integer, dimension(taper_order_d0),intent(in) :: taper_order
        real*8, intent(in) :: lat
        real*8, intent(in) :: lon
        integer, intent(in) :: lwin
        integer, intent(in) :: lmax
        real*8, dimension(admit_d0),intent(out) :: admit
        real*8, dimension(corr_d0),intent(out) :: corr
        integer, intent(in) :: k
        real*8, optional,dimension(admit_error_d0),intent(out) :: admit_error
        real*8, optional,dimension(corr_error_d0),intent(out) :: corr_error
        real*8, optional,dimension(taper_wt_d0),intent(in) :: taper_wt
        integer, optional,intent(in) :: mtdef
        integer, optional,intent(in) :: k1linsig
        integer, intent(in) :: taper_order_d0
        integer, intent(in) :: g_d0
        integer, intent(in) :: g_d1
        integer, intent(in) :: g_d2
        integer, intent(in) :: taper_wt_d0
        integer, intent(in) :: corr_error_d0
        integer, intent(in) :: admit_d0
        integer, intent(in) :: admit_error_d0
        integer, intent(in) :: corr_d0
        integer, intent(in) :: tapers_d0
        integer, intent(in) :: tapers_d1
        integer, intent(in) :: t_d0
        integer, intent(in) :: t_d1
        integer, intent(in) :: t_d2
        if(taper_wt(1) < 0.d0) then
            if (k1linsig<0) then
                call SHLocalizedAdmitCorr(tapers,taper_order,lwin,lat,lon,g,t,&
                                          lmax,admit,corr,k,&
                                          admit_error=admit_error,&
                                          corr_error=corr_error,&
                                          mtdef=mtdef, exitstatus=exitstatus)
            else
                call SHLocalizedAdmitCorr(tapers,taper_order,lwin,lat,lon,g,t,&
                                          lmax,admit,corr,k,&
                                          admit_error=admit_error,&
                                          corr_error=corr_error, &
                                          mtdef=mtdef,k1linsig=k1linsig, &
                                          exitstatus=exitstatus)
            endif
        else
            if (k1linsig<0) then
                call SHLocalizedAdmitCorr(tapers,taper_order,lwin,lat,lon,g,t,&
                                          lmax,admit,corr,k,&
                                          admit_error=admit_error,&
                                          corr_error=corr_error,&
                                          taper_wt=taper_wt,mtdef=mtdef, &
                                          exitstatus=exitstatus)
            else
                call SHLocalizedAdmitCorr(tapers,taper_order,lwin,lat,lon,g,t,&
                                          lmax,admit,corr,k,&
                                          admit_error=admit_error,&
                                          corr_error=corr_error,&
                                          taper_wt=taper_wt,mtdef=mtdef,&
                                          k1linsig=k1linsig, &
                                          exitstatus=exitstatus)
            endif
        endif
    end subroutine pySHLocalizedAdmitCorr

    subroutine pySHReturnTapers(exitstatus,theta0,lmax,tapers,eigenvalues,&
                                taper_order,eigenvalues_d0,tapers_d0,&
                                tapers_d1,taper_order_d0)
        use shtools, only: SHReturnTapers
        implicit none
        integer, intent(out) :: exitstatus
        real*8, intent(in) :: theta0
        integer, intent(in) :: lmax
        real*8, dimension(tapers_d0,tapers_d1),intent(out) :: tapers
        real*8, dimension(eigenvalues_d0),intent(out) :: eigenvalues
        integer, dimension(taper_order_d0),intent(out) :: taper_order
        integer, intent(in) :: eigenvalues_d0
        integer, intent(in) :: tapers_d0
        integer, intent(in) :: tapers_d1
        integer, intent(in) :: taper_order_d0
        call SHReturnTapers(theta0,lmax,tapers,eigenvalues,taper_order,&
                            exitstatus=exitstatus)
    end subroutine pySHReturnTapers

    subroutine pySHReturnTapersM(exitstatus,theta0,lmax,m,tapers,eigenvalues, &
                                 tapers_d0,tapers_d1,eigenvalues_d0)
        use shtools, only: SHReturnTapersM
        implicit none
        integer, intent(out) :: exitstatus
        real*8, intent(in) :: theta0
        integer, intent(in) :: lmax
        integer, intent(in) :: m
        real*8, dimension(tapers_d0,tapers_d1),intent(out) :: tapers
        real*8, dimension(eigenvalues_d0),intent(out) :: eigenvalues
        integer, intent(in) :: tapers_d0
        integer, intent(in) :: tapers_d1
        integer, intent(in) :: eigenvalues_d0
        call SHReturnTapersM(theta0,lmax,m,tapers,eigenvalues,&
                             exitstatus=exitstatus)
    end subroutine pySHReturnTapersM

    subroutine pyComputeDm(exitstatus,dllm,lmax,m,theta0,dllm_d0,dllm_d1)
        use shtools, only: ComputeDm
        implicit none
        integer, intent(out) :: exitstatus
        real*8, dimension(dllm_d0,dllm_d1),intent(out) :: dllm
        integer, intent(in) :: lmax
        integer, intent(in) :: m
        real*8, intent(in) :: theta0
        integer, intent(in) :: dllm_d0
        integer, intent(in) :: dllm_d1
        call ComputeDm(dllm,lmax,m,theta0,exitstatus=exitstatus)
    end subroutine pyComputeDm

    subroutine pyComputeDG82(exitstatus,dG82,lmax,m,theta0,dG82_d0,dG82_d1)
        use shtools, only: ComputeDG82
        implicit none
        integer, intent(out) :: exitstatus
        real*8, dimension(dG82_d0,dG82_d1),intent(out) :: dG82
        integer, intent(in) :: lmax
        integer, intent(in) :: m
        real*8, intent(in) :: theta0
        integer, intent(in) :: dG82_d0
        integer, intent(in) :: dG82_d1
        call ComputeDG82(dG82,lmax,m,theta0,exitstatus=exitstatus)
    end subroutine pyComputeDG82

    function pySHFindLWin(theta0,m,alpha,taper_number)
        use shtools, only: SHFindLWin
        implicit none
        real*8, intent(in) :: theta0
        integer, intent(in) :: m
        real*8, intent(in) :: alpha
        integer, optional,intent(in) :: taper_number
        integer :: pySHFindLWin
        pySHFindLWin=SHFindLWin(theta0,m,alpha,taper_number=taper_number)
    end function pySHFindLWin

    subroutine pySHBiasK(exitstatus,tapers,lwin,k,incspectra,ldata,outcspectra, &
                         taper_wt,save_cg,taper_wt_d0,tapers_d0, &
                         tapers_d1,incspectra_d0,outcspectra_d0)
        use shtools, only: SHBiasK
        implicit none
        integer, intent(out) :: exitstatus
        real*8, dimension(tapers_d0,tapers_d1),intent(in) :: tapers
        integer, intent(in) :: lwin
        integer, intent(in) :: k
        real*8, dimension(incspectra_d0),intent(in) :: incspectra
        integer, intent(in) :: ldata
        real*8, dimension(outcspectra_d0),intent(out) :: outcspectra
        real*8, optional,dimension(taper_wt_d0),intent(in) :: taper_wt
        integer, optional,intent(in) :: save_cg
        integer, intent(in) :: taper_wt_d0
        integer, intent(in) :: tapers_d0
        integer, intent(in) :: tapers_d1
        integer, intent(in) :: incspectra_d0
        integer, intent(in) :: outcspectra_d0
        if (taper_wt(1) < 0.d0) then
            call SHBiasK(tapers,lwin,k,incspectra,ldata,outcspectra, &
                         save_cg=save_cg,exitstatus=exitstatus)
        else
            call SHBiasK(tapers,lwin,k,incspectra,ldata,outcspectra, &
                         taper_wt=taper_wt,save_cg=save_cg,&
                         exitstatus=exitstatus)
        endif
    end subroutine pySHBiasK

    subroutine pySHMTCouplingMatrix(exitstatus,Mmt, lmax, tapers_power, lwin, &
                                    k, taper_wt, Mmt_d0, Mmt_d1, &
                                    tapers_power_d0, tapers_power_d1, &
                                    taper_wt_d0)
        use shtools, only: SHMTCouplingMatrix
        implicit none
        integer, intent(out) :: exitstatus
        integer, intent(in) :: lmax, k, lwin
        real*8, intent(out) :: Mmt(Mmt_d0, Mmt_d1)
        real*8, intent(in) ::  tapers_power(tapers_power_d0, tapers_power_d1)
        real*8, optional,dimension(taper_wt_d0),intent(in) :: taper_wt
        integer, intent(in) :: Mmt_d0
        integer, intent(in) :: Mmt_d1
        integer, intent(in) :: tapers_power_d0
        integer, intent(in) :: tapers_power_d1
        integer, intent(in) :: taper_wt_d0
        if (taper_wt(1) < 0.d0) then
            call SHMTCouplingMatrix(Mmt, lmax, tapers_power, lwin, k, &
                                    exitstatus = exitstatus)
        else
            call SHMTCouplingMatrix(Mmt, lmax, tapers_power, lwin, k, &
                                    taper_wt = taper_wt, &
                                    exitstatus = exitstatus)
        endif
    end subroutine

    subroutine pySHBiasAdmitCorr(exitstatus,sgt,sgg,stt,lmax,tapers,lwin,k,&
                                 admit,corr,mtdef,taper_wt,taper_wt_d0,sgt_d0,&
                                 stt_d0,admit_d0,tapers_d0,tapers_d1,corr_d0,&
                                 sgg_d0)
        use shtools, only: SHBiasAdmitCorr
        implicit none
        integer, intent(out) :: exitstatus
        real*8, dimension(sgt_d0),intent(in) :: sgt
        real*8, dimension(sgg_d0),intent(in) :: sgg
        real*8, dimension(stt_d0),intent(in) :: stt
        integer, intent(in) :: lmax
        real*8, dimension(tapers_d0,tapers_d1),intent(in) :: tapers
        integer, intent(in) :: lwin
        integer, intent(in) :: k
        real*8, dimension(admit_d0),intent(out) :: admit
        real*8, dimension(corr_d0),intent(out) :: corr
        integer, optional,intent(in) :: mtdef
        real*8, optional,dimension(taper_wt_d0),intent(in) :: taper_wt
        integer, intent(in) :: taper_wt_d0
        integer, intent(in) :: sgt_d0
        integer, intent(in) :: stt_d0
        integer, intent(in) :: admit_d0
        integer, intent(in) :: tapers_d0
        integer, intent(in) :: tapers_d1
        integer, intent(in) :: corr_d0
        integer, intent(in) :: sgg_d0
        if (taper_wt(1) < 0.d0) then
            call SHBiasAdmitCorr(sgt,sgg,stt,lmax,tapers,lwin,k,admit,corr, &
                                 mtdef=mtdef,exitstatus=exitstatus)
        else
            call SHBiasAdmitCorr(sgt,sgg,stt,lmax,tapers,lwin,k,admit,corr, &
                                 mtdef=mtdef,taper_wt=taper_wt,&
                                 exitstatus=exitstatus)
        endif
    end subroutine pySHBiasAdmitCorr

    subroutine pySHMTDebias(exitstatus,mtdebias,mtspectra,lmax,tapers,lwin,k,&
                            nl,lmid,n,taper_wt,mtdebias_d0,mtdebias_d1,&
                            taper_wt_d0,mtspectra_d0,mtspectra_d1,tapers_d0,&
                            tapers_d1,lmid_d0)
        use shtools, only: SHMTDebias
        implicit none
        integer, intent(out) :: exitstatus
        real*8, dimension(mtdebias_d0,mtdebias_d1),intent(out) :: mtdebias
        real*8, dimension(mtspectra_d0,mtspectra_d1),intent(in) :: mtspectra
        integer, intent(in) :: lmax
        real*8, dimension(tapers_d0,tapers_d1),intent(in) :: tapers
        integer, intent(in) :: lwin
        integer, intent(in) :: k
        integer, intent(in) :: nl
        real*8, dimension(lmid_d0),intent(out) :: lmid
        integer, intent(out) :: n
        real*8, optional,dimension(taper_wt_d0),intent(in) :: taper_wt
        integer, intent(in) :: mtdebias_d0
        integer, intent(in) :: mtdebias_d1
        integer, intent(in) :: taper_wt_d0
        integer, intent(in) :: mtspectra_d0
        integer, intent(in) :: mtspectra_d1
        integer, intent(in) :: tapers_d0
        integer, intent(in) :: tapers_d1
        integer, intent(in) :: lmid_d0
        if (taper_wt(1) < 0.d0) then
            call SHMTDebias(mtdebias,mtspectra,lmax,tapers,lwin,k,nl,lmid,n,&
                            exitstatus = exitstatus)
        else
            call SHMTDebias(mtdebias,mtspectra,lmax,tapers,lwin,k,nl,lmid,n, &
                            taper_wt=taper_wt,exitstatus=exitstatus)
        endif
    end subroutine pySHMTDebias

    subroutine pySHMTVarOpt(exitstatus,l,tapers,taper_order,lwin,kmax,Sff,&
                            var_opt,var_unit,weight_opt,nocross, &
                            taper_order_d0,weight_opt_d0,weight_opt_d1, &
                            var_unit_d0,var_opt_d0,Sff_d0,tapers_d0,tapers_d1)
        use shtools, only: SHMTVarOpt
        implicit none
        integer, intent(out) :: exitstatus
        integer, intent(in) :: l
        real*8, dimension(tapers_d0,tapers_d1),intent(in) :: tapers
        integer, dimension(taper_order_d0),intent(in) :: taper_order
        integer, intent(in) :: lwin
        integer, intent(in) :: kmax
        real*8, dimension(Sff_d0),intent(in) :: Sff
        real*8, dimension(var_opt_d0),intent(out) :: var_opt
        real*8, dimension(var_unit_d0),intent(out) :: var_unit
        real*8, optional,dimension(weight_opt_d0,weight_opt_d1),intent(out) ::&
                                                                     weight_opt
        integer, optional,intent(in) :: nocross
        integer, intent(in) :: taper_order_d0
        integer, intent(in) :: weight_opt_d0
        integer, intent(in) :: weight_opt_d1
        integer, intent(in) :: var_unit_d0
        integer, intent(in) :: var_opt_d0
        integer, intent(in) :: Sff_d0
        integer, intent(in) :: tapers_d0
        integer, intent(in) :: tapers_d1
        call SHMTVarOpt(l,tapers,taper_order,lwin,kmax,Sff,var_opt,var_unit, &
                        weight_opt=weight_opt,nocross=nocross, &
                        exitstatus=exitstatus)
    end subroutine pySHMTVarOpt

    function pySHSjkPG(incspectra,l,m,mprime,hj_real,hk_real,mj,mk,lwin,hkcc, &
                       hk_real_d0,incspectra_d0,hj_real_d0)
        use shtools, only: SHSjkPG
        implicit none
        real*8, dimension(incspectra_d0),intent(in) :: incspectra
        integer, intent(in) :: l
        integer, intent(in) :: m
        integer, intent(in) :: mprime
        real*8, dimension(hj_real_d0),intent(in) :: hj_real
        real*8, dimension(hk_real_d0),intent(in) :: hk_real
        integer, intent(in) :: mj
        integer, intent(in) :: mk
        integer, intent(in) :: lwin
        integer, intent(in) :: hkcc
        integer, intent(in) :: hk_real_d0
        integer, intent(in) :: incspectra_d0
        integer, intent(in) :: hj_real_d0
        complex*16 :: pySHSjkPG
        pySHSjkPG=SHSjkPG(incspectra,l,m,mprime,hj_real,hk_real,mj,mk,lwin, &
                          hkcc)
    end function pySHSjkPG

    subroutine pySHReturnTapersMap(exitstatus,tapers,eigenvalues,dh_mask,n_dh,&
                                   lmax,sampling,ntapers,dh_mask_d0,dh_mask_d1,&
                                   tapers_d0,tapers_d1,eigenvalues_d0)
        use shtools, only: SHReturnTapersMap
        implicit none
        integer, intent(out) :: exitstatus
        real*8, dimension(tapers_d0,tapers_d1),intent(out) :: tapers
        real*8, dimension(eigenvalues_d0),intent(out) :: eigenvalues
        integer, dimension(dh_mask_d0,dh_mask_d1),intent(in) :: dh_mask
        integer, intent(in) :: n_dh
        integer, intent(in), optional :: sampling
        integer, intent(in) :: lmax
        integer, optional,intent(in) :: ntapers
        integer, intent(in) :: dh_mask_d0
        integer, intent(in) :: dh_mask_d1
        integer, intent(in) :: tapers_d0
        integer, intent(in) :: tapers_d1
        integer, intent(in) :: eigenvalues_d0
        call SHReturnTapersMap(tapers,eigenvalues,dh_mask,n_dh,lmax, &
                               sampling=sampling,ntapers=ntapers, &
                               exitstatus = exitstatus)
    end subroutine pySHReturnTapersMap

    subroutine pySHBiasKMask(exitstatus,tapers,lwin,k,incspectra,ldata,&
                             outcspectra,taper_wt,save_cg,taper_wt_d0,&
                             tapers_d0,tapers_d1,incspectra_d0,outcspectra_d0)
        use shtools, only: SHBiasKMask
        implicit none
        integer, intent(out) :: exitstatus
        real*8, dimension(tapers_d0,tapers_d1),intent(in) :: tapers
        integer, intent(in) :: lwin
        integer, intent(in) :: k
        real*8, dimension(incspectra_d0),intent(in) :: incspectra
        integer, intent(in) :: ldata
        real*8, dimension(outcspectra_d0),intent(out) :: outcspectra
        real*8, optional,dimension(taper_wt_d0),intent(in) :: taper_wt
        integer, optional,intent(in) :: save_cg
        integer, intent(in) :: taper_wt_d0
        integer, intent(in) :: tapers_d0
        integer, intent(in) :: tapers_d1
        integer, intent(in) :: incspectra_d0
        integer, intent(in) :: outcspectra_d0
        if (taper_wt(1) < 0.d0) then
            call SHBiasKMask(tapers,lwin,k,incspectra,ldata,outcspectra, &
                             save_cg=save_cg, exitstatus=exitstatus)
        else
            call SHBiasKMask(tapers,lwin,k,incspectra,ldata,outcspectra, &
                             taper_wt=taper_wt,save_cg=save_cg, &
                             exitstatus=exitstatus)
        endif
    end subroutine pySHBiasKMask

    subroutine pySHMultiTaperMaskSE(exitstatus,mtse,sd,sh,lmax,tapers,lmaxt,k,&
                                    taper_wt,norm,csphase, &
                                    taper_wt_d0,sh_d0,sh_d1,sh_d2, &
                                    tapers_d0,tapers_d1,mtse_d0,sd_d0)
        use shtools, only: SHMultiTaperMaskSE
        implicit none
        integer, intent(out) :: exitstatus
        real*8, dimension(mtse_d0),intent(out) :: mtse
        real*8, dimension(sd_d0),intent(out) :: sd
        real*8, dimension(sh_d0,sh_d1,sh_d2),intent(in) :: sh
        integer, intent(in) :: lmax
        real*8, dimension(tapers_d0,tapers_d1),intent(in) :: tapers
        integer, intent(in) :: lmaxt
        integer, intent(in) :: k
        real*8, optional,dimension(taper_wt_d0),intent(in) :: taper_wt
        integer, optional,intent(in) :: norm
        integer, optional,intent(in) :: csphase
        integer, intent(in) :: taper_wt_d0
        integer, intent(in) :: sh_d0
        integer, intent(in) :: sh_d1
        integer, intent(in) :: sh_d2
        integer, intent(in) :: tapers_d0
        integer, intent(in) :: tapers_d1
        integer, intent(in) :: mtse_d0
        integer, intent(in) :: sd_d0
        if (taper_wt(1) < 0.d0) then
            call SHMultiTaperMaskSE(mtse,sd,sh,lmax,tapers,lmaxt,k,norm=norm, &
                                    csphase=csphase, exitstatus=exitstatus)
        else
            call SHMultiTaperMaskSE(mtse,sd,sh,lmax,tapers,lmaxt,k, &
                                    taper_wt=taper_wt,norm=norm, &
                                    csphase=csphase, exitstatus=exitstatus)
        endif
    end subroutine pySHMultiTaperMaskSE

    subroutine pySHMultiTaperMaskCSE(exitstatus,mtse,sd,sh1,lmax1,sh2,lmax2,&
                                     tapers,lmaxt,k,taper_wt,norm,csphase,&
                                     sh1_d0,sh1_d1,sh1_d2,sh2_d0,sh2_d1,&
                                     sh2_d2,taper_wt_d0,tapers_d0,tapers_d1,&
                                     sd_d0,mtse_d0)
        use shtools, only: SHMultiTaperMaskCSE
        implicit none
        integer, intent(out) :: exitstatus
        real*8, dimension(mtse_d0),intent(out) :: mtse
        real*8, dimension(sd_d0),intent(out) :: sd
        real*8, dimension(sh1_d0,sh1_d1,sh1_d2),intent(in) :: sh1
        integer, intent(in) :: lmax1
        real*8, dimension(sh2_d0,sh2_d1,sh2_d2),intent(in) :: sh2
        integer, intent(in) :: lmax2
        real*8, dimension(tapers_d0,tapers_d1),intent(in) :: tapers
        integer, intent(in) :: lmaxt
        integer, intent(in) :: k
        real*8, optional,dimension(taper_wt_d0),intent(in) :: taper_wt
        integer, optional,intent(in) :: norm
        integer, optional,intent(in) :: csphase
        integer, intent(in) :: sh1_d0
        integer, intent(in) :: sh1_d1
        integer, intent(in) :: sh1_d2
        integer, intent(in) :: sh2_d0
        integer, intent(in) :: sh2_d1
        integer, intent(in) :: sh2_d2
        integer, intent(in) :: taper_wt_d0
        integer, intent(in) :: tapers_d0
        integer, intent(in) :: tapers_d1
        integer, intent(in) :: sd_d0
        integer, intent(in) :: mtse_d0
        if(taper_wt(1) < 0.d0) then
            call SHMultiTaperMaskCSE(mtse,sd,sh1,lmax1,sh2,lmax2,tapers, &
                                     lmaxt,k,norm=norm,csphase=csphase, &
                                     exitstatus=exitstatus)
        else
            call SHMultiTaperMaskCSE(mtse,sd,sh1,lmax1,sh2,lmax2,tapers, &
                                     lmaxt,k,taper_wt=taper_wt,norm=norm, &
                                     csphase=csphase,exitstatus=exitstatus)
        endif
    end subroutine pySHMultiTaperMaskCSE

    subroutine pyComputeDMap(exitstatus,Dij,dh_mask,n_dh,lmax,sampling,&
                             dh_mask_d0,dh_mask_d1,Dij_d0,Dij_d1)
        use shtools, only: ComputeDMap
        implicit none
        integer, intent(out) :: exitstatus
        real*8, dimension(Dij_d0,Dij_d1),intent(out) :: Dij
        integer, dimension(dh_mask_d0,dh_mask_d1),intent(in) :: dh_mask
        integer, intent(in) :: n_dh
        integer, intent(in), optional :: sampling
        integer, intent(in) :: lmax
        integer, intent(in) :: dh_mask_d0
        integer, intent(in) :: dh_mask_d1
        integer, intent(in) :: Dij_d0
        integer, intent(in) :: Dij_d1
        call ComputeDMap(Dij,dh_mask,n_dh,lmax,sampling=sampling,&
                         exitstatus=exitstatus)
    end subroutine pyComputeDMap

    subroutine pyCurve2Mask(exitstatus,dhgrid,n,sampling,profile,nprofile,NP, &
                            centralmeridian,profile_d0,profile_d1,dhgrid_d0, &
                            dhgrid_d1)
        use shtools, only: Curve2Mask
        implicit none
        integer, intent(out) :: exitstatus
        integer, dimension(dhgrid_d0,dhgrid_d1),intent(out) :: dhgrid
        integer, intent(in) :: n
        integer, intent(in) :: sampling
        real*8, dimension(profile_d0,profile_d1),intent(in) :: profile
        integer, intent(in) :: nprofile
        integer, intent(in) :: NP
        integer, optional, intent(in) :: centralmeridian
        integer, intent(in) :: profile_d0
        integer, intent(in) :: profile_d1
        integer, intent(in) :: dhgrid_d0
        integer, intent(in) :: dhgrid_d1
        call Curve2Mask(dhgrid,n,sampling,profile,nprofile,NP, &
                        centralmeridian=centralmeridian, &
                        exitstatus=exitstatus)
    end subroutine pyCurve2Mask

    subroutine pySHBias(exitstatus,Shh,lwin,incspectra,ldata,outcspectra,&
                        save_cg,Shh_d0,incspectra_d0,outcspectra_d0)
        use shtools, only: SHBias
        implicit none
        integer, intent(out) :: exitstatus
        real*8, dimension(Shh_d0),intent(in) :: Shh
        integer, intent(in) :: lwin
        real*8, dimension(incspectra_d0),intent(in) :: incspectra
        integer, intent(in) :: ldata
        real*8, dimension(outcspectra_d0),intent(out) :: outcspectra
        integer, optional,intent(in) :: save_cg
        integer, intent(in) :: Shh_d0
        integer, intent(in) :: incspectra_d0
        integer, intent(in) :: outcspectra_d0
        call SHBias(Shh,lwin,incspectra,ldata,outcspectra,save_cg=save_cg,&
                    exitstatus=exitstatus)
    end subroutine pySHBias

    subroutine pySphericalCapCoef(exitstatus,coef,theta,lmax,coef_d0)
        use shtools, only: SphericalCapCoef
        implicit none
        integer, intent(out) :: exitstatus
        real*8, dimension(coef_d0),intent(out) :: coef
        real*8, intent(in) :: theta
        integer, optional,intent(in) :: lmax
        integer, intent(in) :: coef_d0
        call SphericalCapCoef(coef,theta,lmax=lmax,exitstatus=exitstatus)
    end subroutine pySphericalCapCoef

    subroutine pyMakeGravGridDH(exitstatus,cilm,lmax,gm,r0,a,f,rad,theta,phi,&
                                total,pot,n,sampling,lmax_calc,omega,&
                                normal_gravity,phi_d0,phi_d1,total_d0,total_d1,&
                                rad_d0,rad_d1,cilm_d0,cilm_d1,cilm_d2, &
                                theta_d0,theta_d1,pot_d0,pot_d1)
        use shtools, only: MakeGravGridDH
        implicit none
        integer, intent(out) :: exitstatus
        real*8, dimension(cilm_d0,cilm_d1,cilm_d2),intent(in) :: cilm
        integer, intent(in) :: lmax
        real*8, intent(in) :: gm
        real*8, intent(in) :: r0
        real*8, intent(in) :: a
        real*8, intent(in) :: f
        real*8, dimension(rad_d0,rad_d1),intent(out) :: rad
        real*8, dimension(theta_d0,theta_d1),intent(out) :: theta
        real*8, dimension(phi_d0,phi_d1),intent(out) :: phi
        real*8, dimension(total_d0,total_d1),intent(out) :: total
        real*8, dimension(pot_d0,pot_d1),intent(out) :: pot
        integer, intent(out) :: n
        integer, optional,intent(in) :: sampling
        integer, optional,intent(in) :: lmax_calc
        real*8, optional,intent(in) :: omega
        integer, optional,intent(in) :: normal_gravity
        integer, intent(in) :: phi_d0
        integer, intent(in) :: phi_d1
        integer, intent(in) :: total_d0
        integer, intent(in) :: total_d1
        integer, intent(in) :: rad_d0
        integer, intent(in) :: rad_d1
        integer, intent(in) :: cilm_d0
        integer, intent(in) :: cilm_d1
        integer, intent(in) :: cilm_d2
        integer, intent(in) :: theta_d0
        integer, intent(in) :: theta_d1
        integer, intent(in) :: pot_d0
        integer, intent(in) :: pot_d1
        call MakeGravGridDH(cilm,lmax,gm,r0,a,f,rad,theta,phi,total,n, &
                            sampling=sampling,lmax_calc=lmax_calc,omega=omega,&
                            normal_gravity=normal_gravity,pot=pot,&
                            exitstatus=exitstatus)
    end subroutine pyMakeGravGridDH

    subroutine pyMakeGravGradGridDH(exitstatus,cilm,lmax,gm,r0,a,f,vxx,vyy,vzz,&
                                    vxy,vxz,vyz,n,sampling,lmax_calc,vyz_d0,&
                                    vyz_d1,vyy_d0,vyy_d1,cilm_d0,cilm_d1,&
                                    cilm_d2,vzz_d0,vzz_d1,vxy_d0,vxy_d1,vxx_d0,&
                                    vxx_d1,vxz_d0,vxz_d1)
        use shtools, only: MakeGravGradGridDH
        implicit none
        integer, intent(out) :: exitstatus
        real*8, dimension(cilm_d0,cilm_d1,cilm_d2),intent(in) :: cilm
        integer, intent(in) :: lmax
        real*8, intent(in) :: gm
        real*8, intent(in) :: r0
        real*8, intent(in) :: a
        real*8, intent(in) :: f
        real*8, dimension(vxx_d0,vxx_d1),intent(out) :: vxx
        real*8, dimension(vyy_d0,vyy_d1),intent(out) :: vyy
        real*8, dimension(vzz_d0,vzz_d1),intent(out) :: vzz
        real*8, dimension(vxy_d0,vxy_d1),intent(out) :: vxy
        real*8, dimension(vxz_d0,vxz_d1),intent(out) :: vxz
        real*8, dimension(vyz_d0,vyz_d1),intent(out) :: vyz
        integer, intent(out) :: n
        integer, optional,intent(in) :: sampling
        integer, optional,intent(in) :: lmax_calc
        integer, intent(in) :: vyz_d0
        integer, intent(in) :: vyz_d1
        integer, intent(in) :: vyy_d0
        integer, intent(in) :: vyy_d1
        integer, intent(in) :: cilm_d0
        integer, intent(in) :: cilm_d1
        integer, intent(in) :: cilm_d2
        integer, intent(in) :: vzz_d0
        integer, intent(in) :: vzz_d1
        integer, intent(in) :: vxy_d0
        integer, intent(in) :: vxy_d1
        integer, intent(in) :: vxx_d0
        integer, intent(in) :: vxx_d1
        integer, intent(in) :: vxz_d0
        integer, intent(in) :: vxz_d1
        call MakeGravGradGridDH(cilm,lmax,gm,r0,a,f,vxx,vyy,vzz,vxy,vxz,vyz,n,&
                                sampling=sampling,lmax_calc=lmax_calc,&
                                exitstatus=exitstatus)
    end subroutine pyMakeGravGradGridDH

    subroutine pyMakeGeoidGridDH(exitstatus,geoid,cilm,lmax,r0pot,GM,PotRef,&
                                 omega,r,sampling,order,nlat,nlong,lmax_calc,a,&
                                 f,cilm_d0,cilm_d1,cilm_d2,geoid_d0,geoid_d1)
        use shtools, only: MakeGeoidGrid
        implicit none
        integer, intent(out) :: exitstatus
        real*8, dimension(geoid_d0,geoid_d1),intent(out) :: geoid
        real*8, dimension(cilm_d0,cilm_d1,cilm_d2),intent(in) :: cilm
        integer, intent(in) :: lmax
        real*8, intent(in) :: r0pot
        real*8, intent(in) :: GM
        real*8, intent(in) :: PotRef
        real*8, intent(in) :: omega
        real*8, intent(in) :: r
        integer, intent(in) :: sampling
        integer, intent(in) :: order
        integer, intent(out) :: nlat
        integer, intent(out) :: nlong
        integer, optional,intent(in) :: lmax_calc
        real*8, optional,intent(in) :: a
        real*8, optional,intent(in) :: f
        integer, intent(in) :: cilm_d0
        integer, intent(in) :: cilm_d1
        integer, intent(in) :: cilm_d2
        integer, intent(in) :: geoid_d0
        integer, intent(in) :: geoid_d1
        if (sampling == 1) then
            call MakeGeoidGrid(geoid,cilm,lmax,r0pot,GM,PotRef,omega,r,2, &
                               order,nlat,nlong,lmax_calc=lmax_calc,a=a,f=f,&
                               exitstatus=exitstatus)
        elseif (sampling == 2) then
            call MakeGeoidGrid(geoid,cilm,lmax,r0pot,GM,PotRef,omega,r,3, &
                               order,nlat,nlong,lmax_calc=lmax_calc,a=a,f=f,&
                               exitstatus=exitstatus)
        endif
    end subroutine pyMakeGeoidGridDH

    subroutine pyCilmPlusDH(exitstatus,cilm,gridin,lmax,nmax,mass,d,rho,&
                            sampling,n,gridin_d0,gridin_d1,cilm_d0,cilm_d1,&
                            cilm_d2)
        use shtools, only: CilmPlus
        implicit none
        integer, intent(out) :: exitstatus
        real*8, dimension(cilm_d0,cilm_d1,cilm_d2),intent(out) :: cilm
        real*8, dimension(gridin_d0,gridin_d1),intent(in) :: gridin
        integer, intent(in) :: lmax
        integer, intent(in) :: nmax
        real*8, intent(in) :: mass
        real*8, intent(out) :: d
        real*8, intent(in) :: rho
        integer, intent(in) :: sampling
        integer, optional,intent(in) :: n
        integer, intent(in) :: gridin_d0
        integer, intent(in) :: gridin_d1
        integer, intent(in) :: cilm_d0
        integer, intent(in) :: cilm_d1
        integer, intent(in) :: cilm_d2
        if (sampling == 1) then
            call CilmPlus(cilm,gridin,lmax,nmax,mass,d,rho,2,n=n,&
                          exitstatus=exitstatus)
        else 
            call CilmPlus(cilm,gridin,lmax,nmax,mass,d,rho,3,n=n,&
                          exitstatus=exitstatus)
        endif
    end subroutine pyCilmPlusDH

    subroutine pyCilmMinusDH(exitstatus,cilm,gridin,lmax,nmax,mass,d,rho,&
                             sampling,n,gridin_d0,gridin_d1,cilm_d0,cilm_d1,&
                             cilm_d2)
        use shtools, only: CilmMinus
        implicit none
        integer, intent(out) :: exitstatus
        real*8, dimension(cilm_d0,cilm_d1,cilm_d2),intent(out) :: cilm
        real*8, dimension(gridin_d0,gridin_d1),intent(in) :: gridin
        integer, intent(in) :: lmax
        integer, intent(in) :: nmax
        real*8, intent(in) :: mass
        real*8, intent(out) :: d
        real*8, intent(in) :: rho
        integer, intent(in) :: sampling
        integer, optional,intent(in) :: n
        integer, intent(in) :: gridin_d0
        integer, intent(in) :: gridin_d1
        integer, intent(in) :: cilm_d0
        integer, intent(in) :: cilm_d1
        integer, intent(in) :: cilm_d2
        if (sampling == 1) then
            call CilmMinus(cilm,gridin,lmax,nmax,mass,d,rho,2,n=n, &
                           exitstatus=exitstatus)
        else 
            call CilmMinus(cilm,gridin,lmax,nmax,mass,d,rho,3,n=n, &
                           exitstatus=exitstatus)
        endif
    end subroutine pyCilmMinusDH

    subroutine pyCilmPlusRhoHDH(exitstatus,cilm,gridin,rho,lmax,nmax,mass,d,&
                                sampling,n,gridin_d0,gridin_d1,cilm_d0,cilm_d1,&
                                cilm_d2,rho_d0,rho_d1)
        use shtools, only: CilmPlusRhoH
        implicit none
        integer, intent(out) :: exitstatus
        real*8, dimension(cilm_d0,cilm_d1,cilm_d2),intent(out) :: cilm
        real*8, dimension(gridin_d0,gridin_d1),intent(in) :: gridin
        integer, intent(in) :: lmax
        integer, intent(in) :: nmax
        real*8, intent(in) :: mass
        real*8, intent(out) :: d
        real*8, dimension(rho_d0,rho_d1),intent(in) :: rho
        integer, intent(in) :: sampling
        integer, optional,intent(in) :: n
        integer, intent(in) :: gridin_d0
        integer, intent(in) :: gridin_d1
        integer, intent(in) :: cilm_d0
        integer, intent(in) :: cilm_d1
        integer, intent(in) :: cilm_d2
        integer, intent(in) :: rho_d0
        integer, intent(in) :: rho_d1
        if (sampling == 1) then
            call CilmPlusRhoH(cilm,gridin,lmax,nmax,mass,d,rho,2,n=n, &
                              exitstatus=exitstatus)
        else
            call CilmPlusRhoH(cilm,gridin,lmax,nmax,mass,d,rho,3,n=n,&
                              exitstatus=exitstatus)
        endif
    end subroutine pyCilmPlusRhoHDH

    subroutine pyCilmMinusRhoHDH(exitstatus,cilm,gridin,lmax,nmax,mass,d,rho,&
                                 sampling,n,gridin_d0,gridin_d1,cilm_d0,&
                                 cilm_d1,cilm_d2,rho_d0,rho_d1)
        use shtools, only: CilmMinusRhoH
        implicit none
        integer, intent(out) :: exitstatus
        real*8, dimension(cilm_d0,cilm_d1,cilm_d2),intent(out) :: cilm
        real*8, dimension(gridin_d0,gridin_d1),intent(in) :: gridin
        integer, intent(in) :: lmax
        integer, intent(in) :: nmax
        real*8, intent(in) :: mass
        real*8, intent(out) :: d
        real*8, dimension(rho_d0,rho_d1),intent(in) :: rho
        integer, intent(in) :: sampling
        integer, optional,intent(in) :: n
        integer, intent(in) :: gridin_d0
        integer, intent(in) :: gridin_d1
        integer, intent(in) :: cilm_d0
        integer, intent(in) :: cilm_d1
        integer, intent(in) :: cilm_d2
        integer, intent(in) :: rho_d0
        integer, intent(in) :: rho_d1
        if (sampling == 1) then
            call CilmMinusRhoH(cilm,gridin,lmax,nmax,mass,d,rho,2,n=n, &
                               exitstatus=exitstatus)
        else
            call CilmMinusRhoH(cilm,gridin,lmax,nmax,mass,d,rho,3,n=n, &
                               exitstatus=exitstatus)
        endif
    end subroutine pyCilmMinusRhoHDH

    subroutine pyBAtoHilmDH(exitstatus,cilm,ba,griddh,lmax,nmax,mass,r0,rho,&
                            sampling,filter_type,filter_deg,lmax_calc,ba_d0,&
                            ba_d1,ba_d2,griddh_d0,griddh_d1,cilm_d0,cilm_d1,&
                            cilm_d2)
        use shtools, only: BAtoHilm
        implicit none
        integer, intent(out) :: exitstatus
        real*8, dimension(cilm_d0,cilm_d1,cilm_d2),intent(out) :: cilm
        real*8, dimension(ba_d0,ba_d1,ba_d2),intent(in) :: ba
        real*8, dimension(griddh_d0,griddh_d1),intent(in) :: griddh
        integer, intent(in) :: lmax
        integer, intent(in) :: nmax
        real*8, intent(in) :: mass
        real*8, intent(in) :: r0
        real*8, intent(in) :: rho
        integer, intent(in) :: sampling
        integer, optional,intent(in) :: filter_type
        integer, optional,intent(in) :: filter_deg
        integer, optional,intent(in) :: lmax_calc
        integer, intent(in) :: ba_d0
        integer, intent(in) :: ba_d1
        integer, intent(in) :: ba_d2
        integer, intent(in) :: griddh_d0
        integer, intent(in) :: griddh_d1
        integer, intent(in) :: cilm_d0
        integer, intent(in) :: cilm_d1
        integer, intent(in) :: cilm_d2
        if (sampling == 1) then
            call BAtoHilm(cilm,ba,griddh,lmax,nmax,mass,r0,rho,2, &
                          filter_type=filter_type,filter_deg=filter_deg, &
                          lmax_calc=lmax_calc,exitstatus=exitstatus)
        else
            call BAtoHilm(cilm,ba,griddh,lmax,nmax,mass,r0,rho,3, &
                          filter_type=filter_type,filter_deg=filter_deg, &
                          lmax_calc=lmax_calc,exitstatus=exitstatus)
        endif
    end subroutine pyBAtoHilmDH
    
    subroutine pyBAtoHilmRhoHDH(exitstatus,cilm,ba,griddh,rho,lmax,nmax,mass,&
                                r0,sampling,filter_type,filter_deg,lmax_calc,&
                                ba_d0,ba_d1,ba_d2,griddh_d0,griddh_d1,cilm_d0,&
                                cilm_d1,cilm_d2,rho_d0,rho_d1)
        use shtools, only: BAtoHilmRhoH
        implicit none
        integer, intent(out) :: exitstatus
        real*8, dimension(cilm_d0,cilm_d1,cilm_d2),intent(out) :: cilm
        real*8, dimension(ba_d0,ba_d1,ba_d2),intent(in) :: ba
        real*8, dimension(griddh_d0,griddh_d1),intent(in) :: griddh
        integer, intent(in) :: lmax
        integer, intent(in) :: nmax
        real*8, intent(in) :: mass
        real*8, intent(in) :: r0
        real*8, dimension(rho_d0,rho_d1),intent(in) :: rho
        integer, optional,intent(in) :: sampling
        integer, optional,intent(in) :: filter_type
        integer, optional,intent(in) :: filter_deg
        integer, optional,intent(in) :: lmax_calc
        integer, intent(in) :: ba_d0
        integer, intent(in) :: ba_d1
        integer, intent(in) :: ba_d2
        integer, intent(in) :: griddh_d0
        integer, intent(in) :: griddh_d1
        integer, intent(in) :: cilm_d0
        integer, intent(in) :: cilm_d1
        integer, intent(in) :: cilm_d2
        integer, intent(in) :: rho_d0
        integer, intent(in) :: rho_d1
        if (sampling == 1) then
            call BAtoHilmRhoH(cilm,ba,griddh,lmax,nmax,mass,r0,rho,2,&
                              filter_type=filter_type,filter_deg=filter_deg, &
                              lmax_calc=lmax_calc,exitstatus=exitstatus)
        else
            call BAtoHilmRhoH(cilm,ba,griddh,lmax,nmax,mass,r0,rho,3,&
                              filter_type=filter_type,filter_deg=filter_deg, &
                              lmax_calc=lmax_calc,exitstatus=exitstatus)
        endif
    end subroutine pyBAtoHilmRhoHDH

    function pyDownContFilterMA(l,half,r,d)
        use shtools, only: DownContFilterMA
        implicit none
        integer, intent(in) :: l
        integer, intent(in) :: half
        real*8, intent(in) :: r
        real*8, intent(in) :: d
        real*8 :: pyDownContFilterMA
        pyDownContFilterMA=DownContFilterMA(l,half,r,d)
    end function pyDownContFilterMA

    function pyDownContFilterMC(l,half,r,d)
        use shtools, only: DownContFilterMC
        implicit none
        integer, intent(in) :: l
        integer, intent(in) :: half
        real*8, intent(in) :: r
        real*8, intent(in) :: d
        real*8 :: pyDownContFilterMC
        pyDownContFilterMC=DownContFilterMC(l,half,r,d)
    end function pyDownContFilterMC

    function pyNormalGravity(geocentric_lat,gm,omega,a,b)
        use shtools, only: NormalGravity
        implicit none
        real*8, intent(in) :: geocentric_lat
        real*8, intent(in) :: gm
        real*8, intent(in) :: omega
        real*8, intent(in) :: a
        real*8, intent(in) :: b
        real*8 :: pyNormalGravity
        pyNormalGravity=NormalGravity(geocentric_lat,gm,omega,a,b)
    end function pyNormalGravity

    subroutine pyMakeMagGridDH(exitstatus,cilm,lmax,r0,a,f,rad_grid,theta_grid,&
                               phi_grid,total_grid,n,sampling,lmax_calc,&
                               total_grid_d0,total_grid_d1,cilm_d0,cilm_d1,&
                               cilm_d2,rad_grid_d0,rad_grid_d1,theta_grid_d0, &
                               theta_grid_d1,phi_grid_d0,phi_grid_d1)
        use shtools, only: MakeMagGridDH
        implicit none
        integer, intent(out) :: exitstatus
        real*8, dimension(cilm_d0,cilm_d1,cilm_d2),intent(in) :: cilm
        integer, intent(in) :: lmax
        real*8, intent(in) :: r0
        real*8, intent(in) :: a
        real*8, intent(in) :: f
        real*8, dimension(rad_grid_d0,rad_grid_d1),intent(out) :: rad_grid
        real*8, dimension(theta_grid_d0,theta_grid_d1),intent(out) :: &
                                                                     theta_grid
        real*8, dimension(phi_grid_d0,phi_grid_d1),intent(out) :: phi_grid
        real*8, dimension(total_grid_d0,total_grid_d1),intent(out) :: &
                                                                     total_grid
        integer, intent(out) :: n
        integer, optional,intent(in) :: sampling
        integer, optional,intent(in) :: lmax_calc
        integer, intent(in) :: total_grid_d0
        integer, intent(in) :: total_grid_d1
        integer, intent(in) :: cilm_d0
        integer, intent(in) :: cilm_d1
        integer, intent(in) :: cilm_d2
        integer, intent(in) :: rad_grid_d0
        integer, intent(in) :: rad_grid_d1
        integer, intent(in) :: theta_grid_d0
        integer, intent(in) :: theta_grid_d1
        integer, intent(in) :: phi_grid_d0
        integer, intent(in) :: phi_grid_d1
        call MakeMagGridDH(cilm,lmax,r0,a,f,rad_grid,theta_grid,phi_grid, &
                           total_grid,n,sampling=sampling,lmax_calc=lmax_calc,&
                           exitstatus=exitstatus)
    end subroutine pyMakeMagGridDH

    subroutine pyMakeCircleCoord(exitstatus,coord,lat,lon,theta0,cinterval,cnum, &
                                 coord_d0,coord_d1)
        use shtools, only: MakeCircleCoord
        implicit none
        integer, intent(out) :: exitstatus
        real*8, dimension(coord_d0,coord_d1),intent(out) :: coord
        real*8, intent(in) :: lat
        real*8, intent(in) :: lon
        real*8, intent(in) :: theta0
        real*8, optional,intent(in) :: cinterval
        integer, optional,intent(out) :: cnum
        integer, intent(in) :: coord_d0
        integer, intent(in) :: coord_d1
        call MakeCircleCoord(coord,lat,lon,theta0,cinterval=cinterval,&
                             cnum=cnum, exitstatus=exitstatus)
    end subroutine pyMakeCircleCoord

    subroutine pyMakeEllipseCoord(exitstatus,coord,lat,lon,dec,A_theta,B_theta,cinterval,&
                                  cnum,coord_d0,coord_d1)
        use shtools, only: MakeEllipseCoord
        implicit none
        integer, intent(out) :: exitstatus
        real*8, dimension(coord_d0,coord_d1),intent(out) :: coord
        real*8, intent(in) :: lat
        real*8, intent(in) :: lon
        real*8, intent(in) :: dec
        real*8, intent(in) :: A_theta
        real*8, intent(in) :: B_theta
        real*8, optional,intent(in) :: cinterval
        integer, optional,intent(out) :: cnum
        integer, intent(in) :: coord_d0
        integer, intent(in) :: coord_d1
        call MakeEllipseCoord(coord,lat,lon,dec,A_theta,B_theta, &
                                cinterval=cinterval,cnum=cnum, exitstatus=exitstatus)
    end subroutine pyMakeEllipseCoord

    subroutine pyWigner3j(exitstatus,w3j,jmin,jmax,j2,j3,m1,m2,m3,w3j_d0)
        use shtools, only: Wigner3j
        implicit none
        integer, intent(out) :: exitstatus
        real*8, dimension(w3j_d0),intent(out) :: w3j
        integer, intent(out) :: jmin
        integer, intent(out) :: jmax
        integer, intent(in) :: j2
        integer, intent(in) :: j3
        integer, intent(in) :: m1
        integer, intent(in) :: m2
        integer, intent(in) :: m3
        integer, intent(in) :: w3j_d0
        call Wigner3j(w3j,jmin,jmax,j2,j3,m1,m2,m3, exitstatus=exitstatus)
    end subroutine pyWigner3j

    subroutine pyDHaj(exitstatus,n,aj,aj_d0)
        use shtools, only: DHaj
        implicit none
        integer, intent(out) :: exitstatus
        integer, intent(in) :: n
        real*8, dimension(aj_d0),intent(out) :: aj
        integer, intent(in) :: aj_d0
        call DHaj(n,aj,exitstatus=exitstatus)
    end subroutine pyDHaj
