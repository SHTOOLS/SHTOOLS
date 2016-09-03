    subroutine pyPlmBar(p,lmax,z,csphase,cnorm,p_d0)
        use shtools, only: PlmBar
        implicit none
        real*8, dimension(p_d0),intent(out) :: p
        integer, intent(in) :: lmax
        real*8, intent(in) :: z
        integer, optional,intent(in) :: csphase
        integer, optional,intent(in) :: cnorm
        integer, intent(in) :: p_d0
        integer :: exitstatus = 0
        external :: pystop
        call PlmBar(p,lmax,z,csphase=csphase,cnorm=cnorm,exitstatus=exitstatus)
        if (exitstatus /= 0) call pystop(exitstatus)
    end subroutine pyPlmBar

    subroutine pyPlmBar_d1(p,dp,lmax,z,csphase,cnorm,p_d0,dp_d0)
        use shtools, only: PlmBar_d1
        implicit none
        real*8, dimension(p_d0),intent(out) :: p
        real*8, dimension(dp_d0),intent(out) :: dp
        integer, intent(in) :: lmax
        real*8, intent(in) :: z
        integer, optional,intent(in) :: csphase
        integer, optional,intent(in) :: cnorm
        integer, intent(in) :: p_d0
        integer, intent(in) :: dp_d0
        call PlmBar_d1(p,dp,lmax,z,csphase=csphase,cnorm=cnorm)
    end subroutine pyPlmBar_d1

    subroutine pyPlBar(p,lmax,z,p_d0)
        use shtools, only: PlBar
        implicit none
        real*8, dimension(p_d0),intent(out) :: p
        integer, intent(in) :: lmax
        real*8, intent(in) :: z
        integer, intent(in) :: p_d0
        call PlBar(p,lmax,z)
    end subroutine pyPlBar

    subroutine pyPlBar_d1(p,dp,lmax,z,p_d0,dp_d0)
        use shtools, only: PlBar_d1
        implicit none
        real*8, dimension(p_d0),intent(out) :: p
        real*8, dimension(dp_d0),intent(out) :: dp
        integer, intent(in) :: lmax
        real*8, intent(in) :: z
        integer, intent(in) :: p_d0
        integer, intent(in) :: dp_d0
        call PlBar_d1(p,dp,lmax,z)
    end subroutine pyPlBar_d1

    subroutine pyPlmSchmidt(p,lmax,z,csphase,cnorm,p_d0)
        use shtools, only: PlmSchmidt
        implicit none
        real*8, dimension(p_d0),intent(out) :: p
        integer, intent(in) :: lmax
        real*8, intent(in) :: z
        integer, optional,intent(in) :: csphase
        integer, optional,intent(in) :: cnorm
        integer, intent(in) :: p_d0
        call PlmSchmidt(p,lmax,z,csphase=csphase,cnorm=cnorm)
    end subroutine pyPlmSchmidt

    subroutine pyPlSchmidt(p,lmax,z,p_d0)
        use shtools, only: PlSchmidt
        implicit none
        real*8, dimension(p_d0),intent(out) :: p
        integer, intent(in) :: lmax
        real*8, intent(in) :: z
        integer, intent(in) :: p_d0
        call PlSchmidt(p,lmax,z)
    end subroutine pyPlSchmidt

    subroutine pyPlmSchmidt_d1(p,dp,lmax,z,csphase,cnorm,p_d0,dp_d0)
        use shtools, only: PlmSchmidt_d1
        implicit none
        real*8, dimension(p_d0),intent(out) :: p
        real*8, dimension(dp_d0),intent(out) :: dp
        integer, intent(in) :: lmax
        real*8, intent(in) :: z
        integer, optional,intent(in) :: csphase
        integer, optional,intent(in) :: cnorm
        integer, intent(in) :: p_d0
        integer, intent(in) :: dp_d0
        call PlmSchmidt_d1(p,dp,lmax,z,csphase=csphase,cnorm=cnorm)
    end subroutine pyPlmSchmidt_d1

    subroutine pyPlSchmidt_d1(p,dp,lmax,z,p_d0,dp_d0)
        use shtools, only: PlSchmidt_d1
        implicit none
        real*8, dimension(p_d0),intent(out) :: p
        real*8, dimension(dp_d0),intent(out) :: dp
        integer, intent(in) :: lmax
        real*8, intent(in) :: z
        integer, intent(in) :: p_d0
        integer, intent(in) :: dp_d0
        call PlSchmidt_d1(p,dp,lmax,z)
    end subroutine pyPlSchmidt_d1

    subroutine pyPLegendre(p,lmax,z,p_d0)
        use shtools, only: PLegendre
        implicit none
        real*8, dimension(p_d0),intent(out) :: p
        integer, intent(in) :: lmax
        real*8, intent(in) :: z
        integer, intent(in) :: p_d0
        call PLegendre(p,lmax,z)
    end subroutine pyPLegendre

    subroutine pyPLegendreA(p,lmax,z,csphase,p_d0)
        use shtools, only: PLegendreA
        implicit none
        real*8, dimension(p_d0),intent(out) :: p
        integer, intent(in) :: lmax
        real*8, intent(in) :: z
        integer, optional,intent(in) :: csphase
        integer, intent(in) :: p_d0
        call PLegendreA(p,lmax,z,csphase=csphase)
    end subroutine pyPLegendreA

    subroutine pyPLegendre_d1(p,dp,lmax,z,p_d0,dp_d0)
        use shtools, only: PLegendre_d1
        implicit none
        real*8, dimension(p_d0),intent(out) :: p
        real*8, dimension(dp_d0),intent(out) :: dp
        integer, intent(in) :: lmax
        real*8, intent(in) :: z
        integer, intent(in) :: p_d0
        integer, intent(in) :: dp_d0
        call PLegendre_d1(p,dp,lmax,z)
    end subroutine pyPLegendre_d1

    subroutine pyPLegendreA_d1(p,dp,lmax,z,csphase,p_d0,dp_d0)
        use shtools, only: PLegendreA_d1
        implicit none
        real*8, dimension(p_d0),intent(out) :: p
        real*8, dimension(dp_d0),intent(out) :: dp
        integer, intent(in) :: lmax
        real*8, intent(in) :: z
        integer, optional,intent(in) :: csphase
        integer, intent(in) :: p_d0
        integer, intent(in) :: dp_d0
        call PLegendreA_d1(p,dp,lmax,z,csphase=csphase)
    end subroutine pyPLegendreA_d1

    subroutine pyCilmPlusDH(cilm,gridin,lmax,nmax,mass,d,rho,sampling,n, &
                            gridin_d0,gridin_d1,cilm_d0,cilm_d1,cilm_d2)
        use shtools, only: CilmPlus
        implicit none
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
            call CilmPlus(cilm,gridin,lmax,nmax,mass,d,rho,2,n=n)
        else 
            call CilmPlus(cilm,gridin,lmax,nmax,mass,d,rho,3,n=n)
        endif
    end subroutine pyCilmPlusDH

    subroutine pyCilmMinusDH(cilm,gridin,lmax,nmax,mass,d,rho,sampling,n, &
                             gridin_d0,gridin_d1,cilm_d0,cilm_d1,cilm_d2)
        use shtools, only: CilmMinus
        implicit none
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
            call CilmMinus(cilm,gridin,lmax,nmax,mass,d,rho,2,n=n)
        else 
            call CilmMinus(cilm,gridin,lmax,nmax,mass,d,rho,3,n=n)
        endif
    end subroutine pyCilmMinusDH

    subroutine pyCilmPlusRhoHDH(cilm,gridin,lmax,nmax,mass,d,rho,sampling,n, &
                                gridin_d0,gridin_d1,cilm_d0,cilm_d1,cilm_d2, &
                                rho_d0,rho_d1)
        use shtools, only: CilmPlusRhoH
        implicit none
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
            call CilmPlusRhoH(cilm,gridin,lmax,nmax,mass,d,rho,2,n=n)
        else
            call CilmPlusRhoH(cilm,gridin,lmax,nmax,mass,d,rho,3,n=n)
        endif
    end subroutine pyCilmPlusRhoHDH
    
    subroutine pyCilmMinusRhoHDH(cilm,gridin,lmax,nmax,mass,d,rho,sampling,n, &
                                 gridin_d0,gridin_d1,cilm_d0,cilm_d1,cilm_d2, &
                                 rho_d0,rho_d1)
        use shtools, only: CilmMinusRhoH
        implicit none
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
            call CilmMinusRhoH(cilm,gridin,lmax,nmax,mass,d,rho,2,n=n)
        else
            call CilmMinusRhoH(cilm,gridin,lmax,nmax,mass,d,rho,3,n=n)
        endif
    end subroutine pyCilmMinusRhoHDH
    
    subroutine pyBAtoHilmDH(cilm,ba,griddh,lmax,nmax,mass,r0,rho,sampling, &
                            filter_type,filter_deg,lmax_calc,ba_d0,ba_d1, &
                            ba_d2,griddh_d0,griddh_d1,cilm_d0,cilm_d1,cilm_d2)
        use shtools, only: BAtoHilm
        implicit none
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
                          lmax_calc=lmax_calc)
        else
            call BAtoHilm(cilm,ba,griddh,lmax,nmax,mass,r0,rho,3, &
                          filter_type=filter_type,filter_deg=filter_deg, &
                          lmax_calc=lmax_calc)
        endif   
    end subroutine pyBAtoHilmDH
    
    subroutine pyBAtoHilmRhoHDH(cilm,ba,griddh,lmax,nmax,mass,r0,rho,sampling,&
                                filter_type,filter_deg,lmax_calc,ba_d0,ba_d1, &
                                ba_d2,griddh_d0,griddh_d1,cilm_d0,cilm_d1, &
                                cilm_d2, rho_d0,rho_d1)
        use shtools, only: BAtoHilmRhoH
        implicit none
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
                              lmax_calc=lmax_calc)
        else
            call BAtoHilmRhoH(cilm,ba,griddh,lmax,nmax,mass,r0,rho,3,&
                              filter_type=filter_type,filter_deg=filter_deg, &
                              lmax_calc=lmax_calc)
        endif
    end subroutine pyBAtoHilmRhoHDH

    subroutine pyMakeGrid2d(grid,cilm,lmax,interval,nlat,nlong,norm,csphase,f,&
                            a,north,south,east,west,dealloc,cilm_d0,cilm_d1, &
                            cilm_d2,grid_d0,grid_d1)
        use shtools, only: MakeGrid2d
        implicit none
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
                            west=west,dealloc=dealloc)
        else
            call MakeGrid2d(grid,cilm,lmax,interval,nlat,nlong,norm=norm, &
                            csphase=csphase,f=f,a=a,north=north, &
                            south=south,east=east,west=west,dealloc=dealloc)
        endif
    end subroutine pyMakeGrid2d

    subroutine pyGLQGridCoord(latglq,longlq,lmax,nlat,nlong,latglq_d0, &
                             longlq_d0)
        use shtools, only: GLQGridCoord
        implicit none
        real*8, dimension(latglq_d0),intent(out) :: latglq
        real*8, dimension(longlq_d0),intent(out) :: longlq
        integer, intent(in) :: lmax
        integer, intent(out) :: nlat
        integer, intent(out) :: nlong
        integer, intent(in) :: latglq_d0
        integer, intent(in) :: longlq_d0
        call GLQGridCoord(latglq,longlq,lmax,nlat,nlong)
    end subroutine pyGLQGridCoord

    subroutine pyMakeGridGLQ(gridglq,cilm,lmax,zero,norm,csphase,lmax_calc, &
                             gridglq_d0,gridglq_d1,cilm_d0,cilm_d1,cilm_d2, &
                             zero_d0)
        use shtools, only: MakeGridGLQ
        implicit none
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
                         csphase=csphase,lmax_calc=lmax_calc)
    end subroutine pyMakeGridGLQ

    subroutine pySHExpandGLQ(cilm,lmax,gridglq,w,zero,norm,csphase,lmax_calc, &
                             cilm_d0,cilm_d1,cilm_d2,gridglq_d0,gridglq_d1, &
                             zero_d0,w_d0)
        use shtools, only: SHExpandGLQ
        implicit none
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
                         csphase=csphase,lmax_calc=lmax_calc)
    end subroutine pySHExpandGLQ

    subroutine pyshglq(lmax,zero,w,zero_d0,w_d0)
        use shtools, only: SHGLQ
        implicit none
        integer, intent(in) :: lmax
        real*8, dimension(zero_d0),intent(out) :: zero
        real*8, dimension(w_d0),intent(out) :: w
        integer, intent(in) :: zero_d0
        integer, intent(in) :: w_d0
        call SHGLQ(lmax,zero,w)
    end subroutine pySHGLQ

    subroutine pySHRead(filename,cilm,lmax,lmax_in,skip,cilm_d0,cilm_d1, &
                        cilm_d2)
        use shtools, only: SHRead
        implicit none
        character*(*), intent(in) :: filename
        real*8, dimension(cilm_d0,cilm_d1,cilm_d2),intent(out) :: cilm
        integer, intent(out) :: lmax
        integer, intent(in)  :: lmax_in
        integer, optional,intent(in) :: skip
        integer, intent(in) :: cilm_d0
        integer, intent(in) :: cilm_d1
        integer, intent(in) :: cilm_d2
        call SHRead(filename,cilm,lmax,skip=skip)
    end subroutine pySHRead
    
    subroutine pySHReadH(filename,cilm,lmax,lmax_in,nheader,header,skip, &
                         cilm_d0,cilm_d1,cilm_d2)
        use shtools, only: SHRead
        implicit none
        character*(*), intent(in) :: filename
        real*8, dimension(cilm_d0,cilm_d1,cilm_d2),intent(out) :: cilm
        integer, intent(out) :: lmax
        integer, intent(in)  :: lmax_in
        integer, intent(in) ::  nheader
        real*8, intent(out), dimension(nheader) :: header
        integer, optional,intent(in) :: skip
        integer, intent(in) :: cilm_d0
        integer, intent(in) :: cilm_d1
        integer, intent(in) :: cilm_d2
        call SHRead(filename,cilm,lmax,skip=skip,header=header)
    end subroutine pySHReadH
    
    subroutine pySHReadError(filename,cilm,error,lmax,lmax_in,skip,cilm_d0, &
                             cilm_d1,cilm_d2,error_d0,error_d1,error_d2)
        use shtools, only: SHRead
        implicit none
        character*(*), intent(in) :: filename
        real*8, dimension(cilm_d0,cilm_d1,cilm_d2),intent(out) :: cilm
        real*8, dimension(error_d0,error_d1,error_d2),intent(out) :: error
        integer, intent(out) :: lmax
        integer, intent(in)  :: lmax_in
        integer, optional,intent(in) :: skip
        integer, intent(in) :: cilm_d0
        integer, intent(in) :: cilm_d1
        integer, intent(in) :: cilm_d2
        integer, intent(in) :: error_d0
        integer, intent(in) :: error_d1
        integer, intent(in) :: error_d2
        call SHRead(filename,cilm,lmax,skip=skip,error=error)
    end subroutine pySHReadError
    
    subroutine pySHReadErrorH(filename,cilm,error,lmax,lmax_in,nheader, &
                              header,skip,cilm_d0,cilm_d1,cilm_d2,error_d0, &
                              error_d1,error_d2)
        use shtools, only: SHRead
        implicit none
        character*(*), intent(in) :: filename
        real*8, dimension(cilm_d0,cilm_d1,cilm_d2),intent(out) :: cilm
        real*8, dimension(error_d0,error_d1,error_d2),intent(out) :: error
        integer, intent(out) :: lmax
        integer, intent(in)  :: lmax_in
        integer, intent(in) :: nheader
        real*8, dimension(nheader),intent(out) :: header
        integer, optional,intent(in) :: skip
        integer, intent(in) :: cilm_d0
        integer, intent(in) :: cilm_d1
        integer, intent(in) :: cilm_d2
        integer, intent(in) :: error_d0
        integer, intent(in) :: error_d1
        integer, intent(in) :: error_d2
        call SHRead(filename,cilm,lmax,skip=skip,error=error,header=header)
    end subroutine pySHReadErrorH

    subroutine pyMakeMagGridDH(cilm,lmax,r0,a,f,rad_grid,theta_grid,phi_grid, &
                               total_grid,n,sampling,lmax_calc,total_grid_d0, &
                               total_grid_d1,cilm_d0,cilm_d1,cilm_d2, &
                               rad_grid_d0,rad_grid_d1,theta_grid_d0, &
                               theta_grid_d1,phi_grid_d0,phi_grid_d1)
        use shtools, only: MakeMagGridDH
        implicit none
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
                           total_grid,n,sampling=sampling,lmax_calc=lmax_calc)
    end subroutine pyMakeMagGridDH

    function pySHPowerL(c,l,c_d0,c_d1,c_d2)
        use shtools, only: SHPowerL
        implicit none
        real*8, dimension(c_d0,c_d1,c_d2),intent(in) :: c
        integer, intent(in) :: l
        integer, intent(in) :: c_d0
        integer, intent(in) :: c_d1
        integer, intent(in) :: c_d2
        real*8 :: pySHPowerL
        pySHPowerL=SHPowerL(c,l)
    end function pySHPowerL

    function pySHPowerDensityL(c,l,c_d0,c_d1,c_d2)
        use shtools, only: SHPowerDensityL
        implicit none
        real*8, dimension(c_d0,c_d1,c_d2),intent(in) :: c
        integer, intent(in) :: l
        integer, intent(in) :: c_d0
        integer, intent(in) :: c_d1
        integer, intent(in) :: c_d2
        real*8 :: pySHPowerDensityL
        pySHPowerDensityL=SHPowerDensityL(c,l)
    end function pySHPowerDensityL

    function pySHCrossPowerL(c1,c2,l,c2_d0,c2_d1,c2_d2,c1_d0,c1_d1,c1_d2)
        use shtools, only: SHCrossPowerL
        implicit none
        real*8, dimension(c1_d0,c1_d1,c1_d2),intent(in) :: c1
        real*8, dimension(c2_d0,c2_d1,c2_d2),intent(in) :: c2
        integer, intent(in) :: l
        integer, intent(in) :: c2_d0
        integer, intent(in) :: c2_d1
        integer, intent(in) :: c2_d2
        integer, intent(in) :: c1_d0
        integer, intent(in) :: c1_d1
        integer, intent(in) :: c1_d2
        real*8 :: pySHCrossPowerL
        pySHCrossPowerL=SHCrossPowerL(c1,c2,l)
    end function pySHCrossPowerL

    function pySHCrossPowerDensityL(c1,c2,l,c2_d0,c2_d1,c2_d2,c1_d0,c1_d1, &
                                    c1_d2)
        use shtools, only: SHCrossPowerDensityL
        implicit none
        real*8, dimension(c1_d0,c1_d1,c1_d2),intent(in) :: c1
        real*8, dimension(c2_d0,c2_d1,c2_d2),intent(in) :: c2
        integer, intent(in) :: l
        integer, intent(in) :: c2_d0
        integer, intent(in) :: c2_d1
        integer, intent(in) :: c2_d2
        integer, intent(in) :: c1_d0
        integer, intent(in) :: c1_d1
        integer, intent(in) :: c1_d2
        real*8 :: pySHCrossPowerDensityL
        pySHCrossPowerDensityL=SHCrossPowerDensityL(c1,c2,l)
    end function pySHCrossPowerDensityL

    subroutine pySHPowerSpectrum(c,lmax,spectra,c_d0,c_d1,c_d2,spectra_d0)
        use shtools, only: SHPowerSpectrum
        implicit none
        real*8, dimension(c_d0,c_d1,c_d2),intent(in) :: c
        integer, intent(in) :: lmax
        real*8, dimension(spectra_d0),intent(out) :: spectra
        integer, intent(in) :: c_d0
        integer, intent(in) :: c_d1
        integer, intent(in) :: c_d2
        integer, intent(in) :: spectra_d0
        call SHPowerSpectrum(c,lmax,spectra)
    end subroutine pySHPowerSpectrum

    subroutine pySHPowerSpectrumDensity(c,lmax,spectra,c_d0,c_d1,c_d2, &
                                        spectra_d0)
        use shtools, only: SHPowerSpectrumDensity
        implicit none
        real*8, dimension(c_d0,c_d1,c_d2),intent(in) :: c
        integer, intent(in) :: lmax
        real*8, dimension(spectra_d0),intent(out) :: spectra
        integer, intent(in) :: c_d0
        integer, intent(in) :: c_d1
        integer, intent(in) :: c_d2
        integer, intent(in) :: spectra_d0
        call SHPowerSpectrumDensity(c,lmax,spectra)
    end subroutine pySHPowerSpectrumDensity

    subroutine pySHCrossPowerSpectrum(c1,c2,lmax,cspectra,c2_d0,c2_d1,c2_d2, &
                                      c1_d0,c1_d1,c1_d2,cspectra_d0)
        use shtools, only: SHCrossPowerSpectrum
        implicit none
        real*8, dimension(c1_d0,c1_d1,c1_d2),intent(in) :: c1
        real*8, dimension(c2_d0,c2_d1,c2_d2),intent(in) :: c2
        integer, intent(in) :: lmax
        real*8, dimension(cspectra_d0),intent(out) :: cspectra
        integer, intent(in) :: c2_d0
        integer, intent(in) :: c2_d1
        integer, intent(in) :: c2_d2
        integer, intent(in) :: c1_d0
        integer, intent(in) :: c1_d1
        integer, intent(in) :: c1_d2
        integer, intent(in) :: cspectra_d0
        call SHCrossPowerSpectrum(c1,c2,lmax,cspectra)
    end subroutine pySHCrossPowerSpectrum

    subroutine pySHCrossPowerSpectrumDensity(c1,c2,lmax,cspectra,c2_d0,c2_d1, &
                                             c2_d2,c1_d0,c1_d1,c1_d2, &
                                             cspectra_d0)
        use shtools, only: SHCrossPowerSpectrumDensity
        implicit none
        real*8, dimension(c1_d0,c1_d1,c1_d2),intent(in) :: c1
        real*8, dimension(c2_d0,c2_d1,c2_d2),intent(in) :: c2
        integer, intent(in) :: lmax
        real*8, dimension(cspectra_d0),intent(out) :: cspectra
        integer, intent(in) :: c2_d0
        integer, intent(in) :: c2_d1
        integer, intent(in) :: c2_d2
        integer, intent(in) :: c1_d0
        integer, intent(in) :: c1_d1
        integer, intent(in) :: c1_d2
        integer, intent(in) :: cspectra_d0
        call SHCrossPowerSpectrumDensity(c1,c2,lmax,cspectra)
    end subroutine pySHCrossPowerSpectrumDensity

    subroutine pydjpi2(dj,lmax,dj_d0,dj_d1,dj_d2)
        use shtools, only: djpi2
        implicit none
        real*8, dimension(dj_d0,dj_d1,dj_d2),intent(out) :: dj
        integer, intent(in) :: lmax
        integer, intent(in) :: dj_d0
        integer, intent(in) :: dj_d1
        integer, intent(in) :: dj_d2
        call djpi2(dj,lmax)
    end subroutine pydjpi2

    subroutine pySHrtoc(rcilm,ccilm,degmax,convention,switchcs,rcilm_d0, &
                        rcilm_d1,rcilm_d2,ccilm_d0,ccilm_d1,ccilm_d2)
        use shtools, only: SHrtoc
        implicit none
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
                    switchcs=switchcs)
    end subroutine pySHrtoc

    subroutine pySHctor(ccilm,rcilm,degmax,convention,switchcs,rcilm_d0, &
                        rcilm_d1,rcilm_d2,ccilm_d0,ccilm_d1,ccilm_d2)
        use shtools, only: SHctor
        implicit none
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
                    switchcs=switchcs)
    end subroutine pySHctor

    subroutine pySHCilmToCindex(cilm,cindex,degmax,cindex_d0,cindex_d1, &
                                cilm_d0,cilm_d1,cilm_d2)
        use shtools, only: SHCilmToCindex
        implicit none
        real*8, dimension(cilm_d0,cilm_d1,cilm_d2),intent(in) :: cilm
        real*8, dimension(cindex_d0,cindex_d1),intent(out) :: cindex
        integer, optional,intent(in) :: degmax
        integer, intent(in) :: cindex_d0
        integer, intent(in) :: cindex_d1
        integer, intent(in) :: cilm_d0
        integer, intent(in) :: cilm_d1
        integer, intent(in) :: cilm_d2
        call SHCilmToCindex(cilm,cindex,degmax=degmax)
    end subroutine pySHCilmToCindex

    subroutine pySHCindexToCilm(cindex,cilm,degmax,cindex_d0,cindex_d1, &
                                cilm_d0,cilm_d1,cilm_d2)
        use shtools, only: SHCindexToCilm
        implicit none
        real*8, dimension(cindex_d0,cindex_d1),intent(in) :: cindex
        real*8, dimension(cilm_d0,cilm_d1,cilm_d2),intent(out) :: cilm
        integer, optional,intent(in) :: degmax
        integer, intent(in) :: cindex_d0
        integer, intent(in) :: cindex_d1
        integer, intent(in) :: cilm_d0
        integer, intent(in) :: cilm_d1
        integer, intent(in) :: cilm_d2
        call SHCindexToCilm(cindex,cilm,degmax=degmax)
    end subroutine pySHCindexToCilm

    subroutine pySHRotateCoef(x,cof,rcof,dj,lmax,rcof_d0,rcof_d1,dj_d0,dj_d1, &
                              dj_d2,cof_d0,cof_d1)
        use shtools, only: SHRotateCoef
        implicit none
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
        call SHRotateCoef(x,cof,rcof,dj,lmax)
    end subroutine pySHRotateCoef

    subroutine pySHRotateRealCoef(cilmrot,cilm,lmax,x,dj,x_d0,dj_d0,dj_d1, &
                                  dj_d2,cilm_d0,cilm_d1,cilm_d2, &
                                  cilmrot_d0,cilmrot_d1,cilmrot_d2)
        use shtools, only: SHRotateRealCoef
        implicit none
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
        call SHRotateRealCoef(cilmrot,cilm,lmax,x,dj)
    end subroutine pySHRotateRealCoef

    subroutine pySHExpandDH(grid,n,cilm,lmax,norm,sampling,csphase,lmax_calc, &
                            cilm_d0,cilm_d1,cilm_d2,grid_d0,grid_d1)
        use shtools, only: SHExpandDH
        implicit none
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
                        csphase=csphase,lmax_calc=lmax_calc)
    end subroutine pySHExpandDH

    subroutine pyMakeGridDH(griddh,n,cilm,lmax,norm,sampling,csphase, &
                            lmax_calc,cilm_d0,cilm_d1,cilm_d2,griddh_d0, &
                            griddh_d1)
        use shtools, only: MakeGridDH
        implicit none
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
                        csphase=csphase,lmax_calc=lmax_calc)
    end subroutine pyMakeGridDH

    function pyMakeGridPoint(cilm,lmax,lat,longitude,norm,csphase,dealloc, &
                             cilm_d0,cilm_d1,cilm_d2)
        use shtools, only: MakeGridPoint
        implicit none
        real*8, dimension(cilm_d0,cilm_d1,cilm_d2),intent(in) :: cilm
        integer, intent(in) :: lmax
        real*8, intent(in) :: lat
        real*8, intent(in) :: longitude
        integer, optional,intent(in) :: norm
        integer, optional,intent(in) :: csphase
        integer, optional,intent(in) :: dealloc
        integer, intent(in) :: cilm_d0
        integer, intent(in) :: cilm_d1
        integer, intent(in) :: cilm_d2
        real*8 :: pyMakeGridPoint
        pyMakeGridPoint=MakeGridPoint(cilm,lmax,lat,longitude,norm=norm, &
                                      csphase=csphase,dealloc=dealloc)
    end function pyMakeGridPoint

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

    subroutine pySHExpandLSQ(cilm,d,lat,lon,nmax,lmax,norm,chi2,csphase, &
                             d_d0,lon_d0,cilm_d0,cilm_d1,cilm_d2,lat_d0)
        use shtools, only: SHExpandLSQ
        implicit none
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
                         csphase=csphase)
    end subroutine pySHExpandLSQ

    subroutine pySHMultiply(shout,sh1,lmax1,sh2,lmax2,precomp,norm,csphase, &
                            sh1_d0,sh1_d1,sh1_d2,sh2_d0,sh2_d1,sh2_d2, &
                            shout_d0,shout_d1,shout_d2)
        use shtools, only: SHMultiply
        implicit none
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
                        norm=norm,csphase=csphase)
    end subroutine pySHMultiply

    subroutine pyComputeDm(dllm,lmax,m,theta0,dllm_d0,dllm_d1)
        use shtools, only: ComputeDm
        implicit none
        real*8, dimension(dllm_d0,dllm_d1),intent(out) :: dllm
        integer, intent(in) :: lmax
        integer, intent(in) :: m
        real*8, intent(in) :: theta0
        integer, intent(in) :: dllm_d0
        integer, intent(in) :: dllm_d1
        call ComputeDm(dllm,lmax,m,theta0)
    end subroutine pyComputeDm

    subroutine pySphericalCapCoef(coef,theta,lmax,coef_d0)
        use shtools, only: SphericalCapCoef
        implicit none
        real*8, dimension(coef_d0),intent(out) :: coef
        real*8, intent(in) :: theta
        integer, optional,intent(in) :: lmax
        integer, intent(in) :: coef_d0
        call SphericalCapCoef(coef,theta,lmax=lmax)
    end subroutine pySphericalCapCoef

    subroutine pySHReturnTapersM(theta0,lmax,m,tapers,eigenvalues, &
                                 tapers_d0,tapers_d1,eigenvalues_d0)
        use shtools, only: SHReturnTapersM
        implicit none
        real*8, intent(in) :: theta0
        integer, intent(in) :: lmax
        integer, intent(in) :: m
        real*8, dimension(tapers_d0,tapers_d1),intent(out) :: tapers
        real*8, dimension(eigenvalues_d0),intent(out) :: eigenvalues
        integer, intent(in) :: tapers_d0
        integer, intent(in) :: tapers_d1
        integer, intent(in) :: eigenvalues_d0
        call SHReturnTapersM(theta0,lmax,m,tapers,eigenvalues)
    end subroutine pySHReturnTapersM

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

    subroutine pySHAdmitCorr(G,T,lmax,admit,admit_error,corr,G_d0,G_d1,G_d2, &
                             admit_d0,admit_error_d0,T_d0,T_d1,T_d2,corr_d0)
        use shtools, only: SHAdmitCorr
        implicit none
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
        call SHAdmitCorr(G,T,lmax,admit,corr,admit_error=admit_error)
    end subroutine pySHAdmitCorr

    subroutine pySHLocalizedAdmitCorr(g,t,tapers,taper_order,k,lat,lon,lwin, &
                                      lmax,admit,corr,admit_error,corr_error, &
                                      taper_wt,mtdef,k1linsig,taper_order_d0, &
                                      g_d0,g_d1,g_d2,taper_wt_d0, &
                                      corr_error_d0,admit_d0,admit_error_d0, &
                                      corr_d0,tapers_d0,tapers_d1,t_d0,t_d1, &
                                      t_d2)
        use shtools, only: SHLocalizedAdmitCorr
        implicit none
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
                                          mtdef=mtdef)
            else
                call SHLocalizedAdmitCorr(tapers,taper_order,lwin,lat,lon,g,t,&
                                          lmax,admit,corr,k,&
                                          admit_error=admit_error,&
                                          corr_error=corr_error, &
                                          mtdef=mtdef,k1linsig=k1linsig)
            endif
        else
            if (k1linsig<0) then
                call SHLocalizedAdmitCorr(tapers,taper_order,lwin,lat,lon,g,t,&
                                          lmax,admit,corr,k,&
                                          admit_error=admit_error,&
                                          corr_error=corr_error,&
                                          taper_wt=taper_wt,mtdef=mtdef)
            else
                call SHLocalizedAdmitCorr(tapers,taper_order,lwin,lat,lon,g,t,&
                                          lmax,admit,corr,k,&
                                          admit_error=admit_error,&
                                          corr_error=corr_error,&
                                          taper_wt=taper_wt,mtdef=mtdef,&
                                          k1linsig=k1linsig)
            endif
        endif
    end subroutine pySHLocalizedAdmitCorr

    subroutine pyComputeDG82(dG82,lmax,m,theta0,dG82_d0,dG82_d1)
        use shtools, only: ComputeDG82
        implicit none
        real*8, dimension(dG82_d0,dG82_d1),intent(out) :: dG82
        integer, intent(in) :: lmax
        integer, intent(in) :: m
        real*8, intent(in) :: theta0
        integer, intent(in) :: dG82_d0
        integer, intent(in) :: dG82_d1
        call ComputeDG82(dG82,lmax,m,theta0)
    end subroutine pyComputeDG82

    subroutine pyWigner3j(w3j,jmin,jmax,j2,j3,m1,m2,m3,w3j_d0)
        use shtools, only: Wigner3j
        implicit none
        real*8, dimension(w3j_d0),intent(out) :: w3j
        integer, intent(out) :: jmin
        integer, intent(out) :: jmax
        integer, intent(in) :: j2
        integer, intent(in) :: j3
        integer, intent(in) :: m1
        integer, intent(in) :: m2
        integer, intent(in) :: m3
        integer, intent(in) :: w3j_d0
        call Wigner3j(w3j,jmin,jmax,j2,j3,m1,m2,m3)
    end subroutine pyWigner3j

    subroutine pySHBias(Shh,lwin,incspectra,ldata,outcspectra,save_cg,Shh_d0, &
                        incspectra_d0,outcspectra_d0)
        use shtools, only: SHBias
        implicit none
        real*8, dimension(Shh_d0),intent(in) :: Shh
        integer, intent(in) :: lwin
        real*8, dimension(incspectra_d0),intent(in) :: incspectra
        integer, intent(in) :: ldata
        real*8, dimension(outcspectra_d0),intent(out) :: outcspectra
        integer, optional,intent(in) :: save_cg
        integer, intent(in) :: Shh_d0
        integer, intent(in) :: incspectra_d0
        integer, intent(in) :: outcspectra_d0
        call SHBias(Shh,lwin,incspectra,ldata,outcspectra,save_cg=save_cg)
    end subroutine pySHBias

    subroutine pySHBiasK(tapers,lwin,k,incspectra,ldata,outcspectra, &
                         taper_wt,save_cg,taper_wt_d0,tapers_d0, &
                         tapers_d1,incspectra_d0,outcspectra_d0)
        use shtools, only: SHBiasK
        implicit none
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
                         save_cg=save_cg)
        else
            call SHBiasK(tapers,lwin,k,incspectra,ldata,outcspectra, &
                         taper_wt=taper_wt,save_cg=save_cg)
        endif
    end subroutine pySHBiasK

    subroutine pySHBiasKMask(tapers,lwin,k,incspectra,ldata,outcspectra, &
                         taper_wt,save_cg,taper_wt_d0,tapers_d0, &
                         tapers_d1,incspectra_d0,outcspectra_d0)
        use shtools, only: SHBiasKMask
        implicit none
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
                         save_cg=save_cg)
        else
            call SHBiasKMask(tapers,lwin,k,incspectra,ldata,outcspectra, &
                         taper_wt=taper_wt,save_cg=save_cg)
        endif
    end subroutine pySHBiasKMask

    subroutine pySHMultiTaperSE(mtse,sd,sh,lmax,tapers,taper_order,lmaxt,k, &
                                lat,lon,taper_wt,norm,csphase, &
                                taper_order_d0,taper_wt_d0,sh_d0,sh_d1,sh_d2, &
                                tapers_d0,tapers_d1,mtse_d0,sd_d0)
        use shtools, only: SHMultiTaperSE
        implicit none
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
                                csphase=csphase)
        else
            call SHMultiTaperSE(mtse,sd,sh,lmax,tapers,taper_order,lmaxt, &
                                k,lat=lat,lon=lon,taper_wt=taper_wt,norm=norm,&
                                csphase=csphase)
        endif
    end subroutine pySHMultiTaperSE

    subroutine pySHMultiTaperCSE(mtse,sd,sh1,lmax1,sh2,lmax2,tapers, &
                                 taper_order,lmaxt,k,lat,lon,taper_wt,norm, &
                                 csphase,sh1_d0,sh1_d1,sh1_d2,sh2_d0,sh2_d1, &
                                 sh2_d2,taper_order_d0,taper_wt_d0,tapers_d0,&
                                 tapers_d1,sd_d0,mtse_d0)
        use shtools, only: SHMultiTaperCSE
        implicit none
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
                                 norm=norm,csphase=csphase)
        else
            call SHMultiTaperCSE(mtse,sd,sh1,lmax1,sh2,lmax2,tapers, &
                                 taper_order,lmaxt,k,lat=lat,lon=lon,&
                                 taper_wt=taper_wt,norm=norm,csphase=csphase)
        endif
    end subroutine pySHMultiTaperCSE

    subroutine pySHMultiTaperMaskSE(mtse,sd,sh,lmax,tapers,lmaxt,k, &
                                    taper_wt,norm,csphase, &
                                    taper_wt_d0,sh_d0,sh_d1,sh_d2, &
                                    tapers_d0,tapers_d1,mtse_d0,sd_d0)
        use shtools, only: SHMultiTaperMaskSE
        implicit none
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
                                    csphase=csphase)
        else
            call SHMultiTaperMaskSE(mtse,sd,sh,lmax,tapers,lmaxt,k, &
                                    taper_wt=taper_wt,norm=norm, &
                                    csphase=csphase)
        endif
    end subroutine pySHMultiTaperMaskSE

    subroutine pySHMultiTaperMaskCSE(mtse,sd,sh1,lmax1,sh2,lmax2,tapers, &
                                     lmaxt,k,taper_wt,norm,csphase,sh1_d0, &
                                     sh1_d1,sh1_d2,sh2_d0,sh2_d1,sh2_d2, &
                                     taper_wt_d0,tapers_d0,tapers_d1,sd_d0, &
                                     mtse_d0)
        use shtools, only: SHMultiTaperMaskCSE
        implicit none
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
                                     lmaxt,k,norm=norm,csphase=csphase)
        else
            call SHMultiTaperMaskCSE(mtse,sd,sh1,lmax1,sh2,lmax2,tapers, &
                                     lmaxt,k,taper_wt=taper_wt,norm=norm, &
                                     csphase=csphase)
        endif
    end subroutine pySHMultiTaperMaskCSE

    subroutine pySHReadJPL(filename,cilm,lmax,lmax_in,gm,formatstring,cilm_d0,&
                           cilm_d1,cilm_d2)
        use shtools, only: SHReadJPL
        implicit none
        character*(*), intent(in) :: filename
        real*8, dimension(cilm_d0,cilm_d1,cilm_d2),intent(out) :: cilm
        integer, intent(out) :: lmax
        integer, intent(in) :: lmax_in
        real*8, optional,dimension(2),intent(out) :: gm
        character*6, optional,intent(in) :: formatstring
        integer, intent(in) :: cilm_d0
        integer, intent(in) :: cilm_d1
        integer, intent(in) :: cilm_d2
        call SHReadJPL(filename,cilm,lmax,gm=gm,formatstring=formatstring)
    end subroutine pySHReadJPL

    subroutine pySHReadJPLError(filename,cilm,error,lmax,lmax_in,gm, &
                                formatstring,cilm_d0,cilm_d1,cilm_d2,error_d0,&
                                error_d1,error_d2)
        use shtools, only: SHReadJPL
        implicit none
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
                       formatstring=formatstring)
    end subroutine pySHReadJPLError

    subroutine pySHRead2(filename,cilm,lmax,lmax_in,gm,r0_pot,dot,doystart, &
                         doyend,epoch,cilm_d0,cilm_d1,cilm_d2,dot_d0,dot_d1, &
                         dot_d2)
        use shtools, only: SHRead2
        implicit none
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
                     doyend=doyend,epoch=epoch)
    end subroutine pySHRead2

    subroutine pySHRead2Error(filename,cilm,error,lmax,lmax_in,gm,r0_pot,dot, &
                              doystart,doyend,epoch,cilm_d0,cilm_d1,cilm_d2, &
                              error_d0,error_d1,error_d2,dot_d0,dot_d1,dot_d2)
        use shtools, only: SHRead2
        implicit none
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
                     doystart=doystart,doyend=doyend,epoch=epoch)
    end subroutine pySHRead2Error
    
    subroutine pyMakeGeoidGridDH(geoid,cilm,lmax,r0pot,GM,PotRef,omega,r, &
                                 sampling,order,nlat,nlong,lmax_calc,a,f, &
                                 cilm_d0,cilm_d1,cilm_d2,geoid_d0,geoid_d1)
        use shtools, only: MakeGeoidGrid
        implicit none
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
                               order,nlat,nlong,lmax_calc=lmax_calc,a=a,f=f)
        elseif (sampling == 2) then
            call MakeGeoidGrid(geoid,cilm,lmax,r0pot,GM,PotRef,omega,r,3, &
                               order,nlat,nlong,lmax_calc=lmax_calc,a=a,f=f)
        endif
    end subroutine pyMakeGeoidGridDH

    subroutine pyMakeGeoidGrid2D(geoid,cilm,lmax,r0pot,GM,PotRef,omega,r, &
                                 order,nlat,nlong,interval,lmax_calc,a,f, &
                                 cilm_d0,cilm_d1,cilm_d2,geoid_d0,geoid_d1)
        use shtools, only: MakeGeoidGrid
        implicit none
        real*8, dimension(geoid_d0,geoid_d1),intent(out) :: geoid
        real*8, dimension(cilm_d0,cilm_d1,cilm_d2),intent(in) :: cilm
        integer, intent(in) :: lmax
        real*8, intent(in) :: r0pot
        real*8, intent(in) :: GM
        real*8, intent(in) :: PotRef
        real*8, intent(in) :: omega
        real*8, intent(in) :: r
        integer, intent(in) :: order
        integer, intent(out) :: nlat
        integer, intent(out) :: nlong
        real*8, optional,intent(in) :: interval
        integer, optional,intent(in) :: lmax_calc
        real*8, optional,intent(in) :: a
        real*8, optional,intent(in) :: f
        integer, intent(in) :: cilm_d0
        integer, intent(in) :: cilm_d1
        integer, intent(in) :: cilm_d2
        integer, intent(in) :: geoid_d0
        integer, intent(in) :: geoid_d1
        call MakeGeoidGrid(geoid,cilm,lmax,r0pot,GM,PotRef,omega,r,4,order, &
                           nlat,nlong,interval=interval,lmax_calc=lmax_calc, &
                           a=a,f=f)
    end subroutine pyMakeGeoidGrid2D

    subroutine pyPlmON(p,lmax,z,csphase,cnorm,p_d0)
        use shtools, only: PlmON
        implicit none
        real*8, dimension(p_d0),intent(out) :: p
        integer, intent(in) :: lmax
        real*8, intent(in) :: z
        integer, optional,intent(in) :: csphase
        integer, optional,intent(in) :: cnorm
        integer, intent(in) :: p_d0
        call PlmON(p,lmax,z,csphase=csphase,cnorm=cnorm)
    end subroutine pyPlmON

    subroutine pyPlON(p,lmax,z,p_d0)
        use shtools, only: PlON
        implicit none
        real*8, dimension(p_d0),intent(out) :: p
        integer, intent(in) :: lmax
        real*8, intent(in) :: z
        integer, intent(in) :: p_d0
        call PlON(p,lmax,z)
    end subroutine pyPlON

    subroutine pyPlmON_d1(p,dp,lmax,z,csphase,cnorm,p_d0,dp_d0)
        use shtools, only: PlmON_d1
        implicit none
        real*8, dimension(p_d0),intent(out) :: p
        real*8, dimension(dp_d0),intent(out) :: dp
        integer, intent(in) :: lmax
        real*8, intent(in) :: z
        integer, optional,intent(in) :: csphase
        integer, optional,intent(in) :: cnorm
        integer, intent(in) :: p_d0
        integer, intent(in) :: dp_d0
        call PlmON_d1(p,dp,lmax,z,csphase=csphase,cnorm=cnorm)
    end subroutine pyPlmON_d1

    subroutine pyPlON_d1(p,dp,lmax,z,p_d0,dp_d0)
        use shtools, only: PlON_d1
        implicit none
        real*8, dimension(p_d0),intent(out) :: p
        real*8, dimension(dp_d0),intent(out) :: dp
        integer, intent(in) :: lmax
        real*8, intent(in) :: z
        integer, intent(in) :: p_d0
        integer, intent(in) :: dp_d0
        call PlON_d1(p,dp,lmax,z)
    end subroutine pyPlON_d1

    subroutine pyMakeCircleCoord(coord,lat,lon,theta0,cinterval,cnum, &
                                 coord_d0,coord_d1)
        use shtools, only: MakeCircleCoord
        implicit none
        real*8, dimension(coord_d0,coord_d1),intent(out) :: coord
        real*8, intent(in) :: lat
        real*8, intent(in) :: lon
        real*8, intent(in) :: theta0
        real*8, optional,intent(in) :: cinterval
        integer, optional,intent(out) :: cnum
        integer, intent(in) :: coord_d0
        integer, intent(in) :: coord_d1
        call MakeCircleCoord(coord,lat,lon,theta0,cinterval=cinterval,&
                             cnum=cnum)
    end subroutine pyMakeCircleCoord

    subroutine pySHReturnTapers(theta0,lmax,tapers,eigenvalues,taper_order, &
                                eigenvalues_d0,tapers_d0,tapers_d1, &
                                taper_order_d0)
        use shtools, only: SHReturnTapers
        implicit none
        real*8, intent(in) :: theta0
        integer, intent(in) :: lmax
        real*8, dimension(tapers_d0,tapers_d1),intent(out) :: tapers
        real*8, dimension(eigenvalues_d0),intent(out) :: eigenvalues
        integer, dimension(taper_order_d0),intent(out) :: taper_order
        integer, intent(in) :: eigenvalues_d0
        integer, intent(in) :: tapers_d0
        integer, intent(in) :: tapers_d1
        integer, intent(in) :: taper_order_d0
        call SHReturnTapers(theta0,lmax,tapers,eigenvalues,taper_order)
    end subroutine pySHReturnTapers

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

    subroutine pySHMTVarOpt(l,tapers,taper_order,lwin,kmax,Sff,var_opt, &
                            var_unit,weight_opt,nocross, &
                            taper_order_d0,weight_opt_d0,weight_opt_d1, &
                            var_unit_d0,var_opt_d0,Sff_d0,tapers_d0,tapers_d1)
        use shtools, only: SHMTVarOpt
        implicit none
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
                        weight_opt=weight_opt,nocross=nocross)
    end subroutine pySHMTVarOpt

    subroutine pySHMTCouplingMatrix(Mmt, lmax, tapers_power, lwin, k, taper_wt,&
                                    Mmt_d0, Mmt_d1, tapers_power_d0, &
                                    tapers_power_d1, taper_wt_d0)
        use shtools, only: SHMTCouplingMatrix
        implicit none
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
            call SHMTCouplingMatrix(Mmt, lmax, tapers_power, lwin, k)
        else
            call SHMTCouplingMatrix(Mmt, lmax, tapers_power, lwin, k, &
                                    taper_wt=taper_wt)
        endif
    end subroutine

    subroutine pySHMTDebias(mtdebias,mtspectra,lmax,tapers,lwin,k,nl,lmid,n, &
                            taper_wt,mtdebias_d0,mtdebias_d1,taper_wt_d0, &
                            mtspectra_d0,mtspectra_d1,tapers_d0,tapers_d1, &
                            lmid_d0)
        use shtools, only: SHMTDebias
        implicit none
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
            call SHMTDebias(mtdebias,mtspectra,lmax,tapers,lwin,k,nl,lmid,n)
        else
            call SHMTDebias(mtdebias,mtspectra,lmax,tapers,lwin,k,nl,lmid,n, &
                            taper_wt=taper_wt)
        endif
    end subroutine pySHMTDebias

    subroutine pyMakeGravGridDH(cilm,lmax,gm,r0,a,f,rad,theta,phi,total,pot,&
                                n,sampling,lmax_calc,omega,normal_gravity, &
                                phi_d0,phi_d1,total_d0,total_d1, &
                                rad_d0,rad_d1,cilm_d0,cilm_d1,cilm_d2, &
                                theta_d0,theta_d1,pot_d0,pot_d1)
        use shtools, only: MakeGravGridDH
        implicit none
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
                            normal_gravity=normal_gravity,pot=pot)
    end subroutine pyMakeGravGridDH

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

    function pySHConfidence(l_conf,r)
        use shtools, only: SHConfidence
        implicit none
        integer, intent(in) :: l_conf
        real*8, intent(in) :: r
        real*8 :: pySHConfidence
        pySHConfidence=SHConfidence(l_conf,r)
    end function pySHConfidence

    function pySHMagPowerL(c,a,r,l,c_d0,c_d1,c_d2)
        use shtools, only: SHMagPowerL
        implicit none
        real*8, dimension(c_d0,c_d1,c_d2),intent(in) :: c
        real*8, intent(in) :: a
        real*8, intent(in) :: r
        integer, intent(in) :: l
        integer, intent(in) :: c_d0
        integer, intent(in) :: c_d1
        integer, intent(in) :: c_d2
        real*8 :: pySHMagPowerL
        pySHMagPowerL=SHMagPowerL(c,a,r,l)
    end function pySHMagPowerL

    subroutine pySHMagPowerSpectrum(c,a,r,lmax,spectra,c_d0,c_d1, &
                                    c_d2,spectra_d0)
        use shtools, only: SHMagPowerSpectrum
        implicit none
        real*8, dimension(c_d0,c_d1,c_d2),intent(in) :: c
        real*8, intent(in) :: a
        real*8, intent(in) :: r
        integer, intent(in) :: lmax
        real*8, dimension(spectra_d0),intent(out) :: spectra
        integer, intent(in) :: c_d0
        integer, intent(in) :: c_d1
        integer, intent(in) :: c_d2
        integer, intent(in) :: spectra_d0
        call SHMagPowerSpectrum(c,a,r,lmax,spectra)
    end subroutine pySHMagPowerSpectrum

    subroutine pySHExpandDHC(grid,n,cilm,lmax,norm,sampling,csphase,lmax_calc,&
                             cilm_d0,cilm_d1,cilm_d2,grid_d0,grid_d1)
        use shtools, only: SHExpandDHC
        implicit none
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
                         csphase=csphase,lmax_calc=lmax_calc)
    end subroutine pySHExpandDHC

    subroutine pyMakeGridDHC(griddh,n,cilm,lmax,norm,sampling,csphase, &
                             lmax_calc,cilm_d0,cilm_d1,cilm_d2,griddh_d0, &
                             griddh_d1)
        use shtools, only: MakeGridDHC
        implicit none
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
                         csphase=csphase,lmax_calc=lmax_calc)
    end subroutine pyMakeGridDHC

    subroutine pyMakeGridGLQC(gridglq,cilm,lmax,zero,norm,csphase,lmax_calc, &
                              gridglq_d0,gridglq_d1,cilm_d0,cilm_d1,cilm_d2, &
                              zero_d0)
        use shtools, only: MakeGridGLQC
        implicit none
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
                          csphase=csphase,lmax_calc=lmax_calc)
    end subroutine pyMakeGridGLQC

    subroutine pySHExpandGLQC(cilm,lmax,gridglq,w,zero,norm,csphase,lmax_calc,&
                              cilm_d0,cilm_d1,cilm_d2,gridglq_d0,gridglq_d1, &
                              zero_d0,w_d0)
        use shtools, only: SHExpandGLQC
        implicit none
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
                          csphase=csphase,lmax_calc=lmax_calc)
    end subroutine pySHExpandGLQC

    function pySHPowerLC(c,l,c_d0,c_d1,c_d2)
        use shtools, only: SHPowerLC
        implicit none
        complex*16, dimension(c_d0,c_d1,c_d2),intent(in) :: c
        integer, intent(in) :: l
        integer, intent(in) :: c_d0
        integer, intent(in) :: c_d1
        integer, intent(in) :: c_d2
        real*8 :: pySHPowerLC
        pySHPowerLC=SHPowerLC(c,l)
    end function pySHPowerLC

    function pySHPowerDensityLC(c,l,c_d0,c_d1,c_d2)
        use shtools, only: SHPowerDensityLC
        implicit none
        complex*16, dimension(c_d0,c_d1,c_d2),intent(in) :: c
        integer, intent(in) :: l
        integer, intent(in) :: c_d0
        integer, intent(in) :: c_d1
        integer, intent(in) :: c_d2
        real*8 :: pySHPowerDensityLC
        pySHPowerDensityLC=SHPowerDensityLC(c,l)
    end function pySHPowerDensityLC

    function pySHCrossPowerLC(c1,c2,l,c2_d0,c2_d1,c2_d2,c1_d0,c1_d1,c1_d2)
        use shtools, only: SHCrossPowerLC
        implicit none
        complex*16, dimension(c1_d0,c1_d1,c1_d2),intent(in) :: c1
        complex*16, dimension(c2_d0,c2_d1,c2_d2),intent(in) :: c2
        integer, intent(in) :: l
        integer, intent(in) :: c2_d0
        integer, intent(in) :: c2_d1
        integer, intent(in) :: c2_d2
        integer, intent(in) :: c1_d0
        integer, intent(in) :: c1_d1
        integer, intent(in) :: c1_d2
        complex*16 :: pySHCrossPowerLC
        pySHCrossPowerLC=SHCrossPowerLC(c1,c2,l)
    end function pySHCrossPowerLC

    function pySHCrossPowerDensityLC(c1,c2,l,c2_d0,c2_d1,c2_d2,c1_d0, &
                                     c1_d1,c1_d2)
        use shtools, only: SHCrossPowerDensityLC
        implicit none
        complex*16, dimension(c1_d0,c1_d1,c1_d2),intent(in) :: c1
        complex*16, dimension(c2_d0,c2_d1,c2_d2),intent(in) :: c2
        integer, intent(in) :: l
        integer, intent(in) :: c2_d0
        integer, intent(in) :: c2_d1
        integer, intent(in) :: c2_d2
        integer, intent(in) :: c1_d0
        integer, intent(in) :: c1_d1
        integer, intent(in) :: c1_d2
        complex*16 :: pySHCrossPowerDensityLC
        pySHCrossPowerDensityLC=SHCrossPowerDensityLC(c1,c2,l)
    end function pySHCrossPowerDensityLC

    subroutine pySHPowerSpectrumC(c,lmax,spectra,c_d0,c_d1,c_d2,spectra_d0)
        use shtools, only: SHPowerSpectrumC
        implicit none
        complex*16, dimension(c_d0,c_d1,c_d2),intent(in) :: c
        integer, intent(in) :: lmax
        real*8, dimension(spectra_d0),intent(out) :: spectra
        integer, intent(in) :: c_d0
        integer, intent(in) :: c_d1
        integer, intent(in) :: c_d2
        integer, intent(in) :: spectra_d0
        call SHPowerSpectrumC(c,lmax,spectra)
    end subroutine pySHPowerSpectrumC

    subroutine pySHPowerSpectrumDensityC(c,lmax,spectra,c_d0,c_d1,c_d2, &
                                         spectra_d0)
        use shtools, only: SHPowerSpectrumDensityC
        implicit none
        complex*16, dimension(c_d0,c_d1,c_d2),intent(in) :: c
        integer, intent(in) :: lmax
        real*8, dimension(spectra_d0),intent(out) :: spectra
        integer, intent(in) :: c_d0
        integer, intent(in) :: c_d1
        integer, intent(in) :: c_d2
        integer, intent(in) :: spectra_d0
        call SHPowerSpectrumDensityC(c,lmax,spectra)
    end subroutine pySHPowerSpectrumDensityC

    subroutine pySHCrossPowerSpectrumC(c1,c2,lmax,cspectra,c2_d0,c2_d1,c2_d2, &
                                       c1_d0,c1_d1,c1_d2,cspectra_d0)
        use shtools, only: SHCrossPowerSpectrumC
        implicit none
        complex*16, dimension(c1_d0,c1_d1,c1_d2),intent(in) :: c1
        complex*16, dimension(c2_d0,c2_d1,c2_d2),intent(in) :: c2
        integer, intent(in) :: lmax
        complex*16, dimension(cspectra_d0),intent(out) :: cspectra
        integer, intent(in) :: c2_d0
        integer, intent(in) :: c2_d1
        integer, intent(in) :: c2_d2
        integer, intent(in) :: c1_d0
        integer, intent(in) :: c1_d1
        integer, intent(in) :: c1_d2
        integer, intent(in) :: cspectra_d0
        call SHCrossPowerSpectrumC(c1,c2,lmax,cspectra)
    end subroutine pySHCrossPowerSpectrumC

    subroutine pySHCrossPowerSpectrumDensityC(c1,c2,lmax,cspectra,c2_d0,c2_d1,&
                                              c2_d2,c1_d0,c1_d1,c1_d2, &
                                              cspectra_d0)
        use shtools, only: SHCrossPowerSpectrumDensityC
        implicit none
        complex*16, dimension(c1_d0,c1_d1,c1_d2),intent(in) :: c1
        complex*16, dimension(c2_d0,c2_d1,c2_d2),intent(in) :: c2
        integer, intent(in) :: lmax
        complex*16, dimension(cspectra_d0),intent(out) :: cspectra
        integer, intent(in) :: c2_d0
        integer, intent(in) :: c2_d1
        integer, intent(in) :: c2_d2
        integer, intent(in) :: c1_d0
        integer, intent(in) :: c1_d1
        integer, intent(in) :: c1_d2
        integer, intent(in) :: cspectra_d0
        call SHCrossPowerSpectrumDensityC(c1,c2,lmax,cspectra)
    end subroutine pySHCrossPowerSpectrumDensityC

    subroutine pySHBiasAdmitCorr(sgt,sgg,stt,lmax,tapers,lwin,k,admit,corr, &
                                 mtdef,taper_wt,taper_wt_d0,sgt_d0,stt_d0, &
                                 admit_d0,tapers_d0,tapers_d1,corr_d0,sgg_d0)
        use shtools, only: SHBiasAdmitCorr
        implicit none
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
                                 mtdef=mtdef)
        else
            call SHBiasAdmitCorr(sgt,sgg,stt,lmax,tapers,lwin,k,admit,corr, &
                                 mtdef=mtdef,taper_wt=taper_wt)
        endif
    end subroutine pySHBiasAdmitCorr

    subroutine pySHCilmToVector(cilm,vector,lmax,vector_d0,cilm_d0,cilm_d1, &
                                cilm_d2)
        use shtools, only: SHCilmToVector
        implicit none
        real*8, dimension(cilm_d0,cilm_d1,cilm_d2),intent(in) :: cilm
        real*8, dimension(vector_d0),intent(out) :: vector
        integer, intent(in) :: lmax
        integer, intent(in) :: vector_d0
        integer, intent(in) :: cilm_d0
        integer, intent(in) :: cilm_d1
        integer, intent(in) :: cilm_d2
        call SHCilmToVector(cilm,vector,lmax)
    end subroutine pySHCilmToVector

    subroutine pySHVectorToCilm(vector,cilm,lmax,vector_d0,cilm_d0, &
                                cilm_d1,cilm_d2)
        use shtools, only: SHVectorToCilm
        implicit none
        real*8, dimension(vector_d0),intent(in) :: vector
        real*8, dimension(cilm_d0,cilm_d1,cilm_d2),intent(out) :: cilm
        integer, intent(in) :: lmax
        integer, intent(in) :: vector_d0
        integer, intent(in) :: cilm_d0
        integer, intent(in) :: cilm_d1
        integer, intent(in) :: cilm_d2
        call SHVectorToCilm(vector,cilm,lmax)
    end subroutine pySHVectorToCilm

    subroutine pyComputeDMap(Dij,dh_mask,n_dh,lmax,sampling,dh_mask_d0, &
                             dh_mask_d1,Dij_d0,Dij_d1)
        use shtools, only: ComputeDMap
        implicit none
        real*8, dimension(Dij_d0,Dij_d1),intent(out) :: Dij
        integer, dimension(dh_mask_d0,dh_mask_d1),intent(in) :: dh_mask
        integer, intent(in) :: n_dh
        integer, intent(in), optional :: sampling
        integer, intent(in) :: lmax
        integer, intent(in) :: dh_mask_d0
        integer, intent(in) :: dh_mask_d1
        integer, intent(in) :: Dij_d0
        integer, intent(in) :: Dij_d1
        call ComputeDMap(Dij,dh_mask,n_dh,lmax,sampling=sampling)
    end subroutine pyComputeDMap

    subroutine pySHReturnTapersMap(tapers,eigenvalues,dh_mask,n_dh,lmax, &
                                   sampling,ntapers,dh_mask_d0,dh_mask_d1, &
                                   tapers_d0,tapers_d1,eigenvalues_d0)
        use shtools, only: SHReturnTapersMap
        implicit none
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
                               sampling=sampling,ntapers=ntapers)
    end subroutine pySHReturnTapersMap

    subroutine pyCurve2Mask(dhgrid,n,sampling,profile,nprofile,NP, &
                            centralmeridian,profile_d0,profile_d1,dhgrid_d0, &
                            dhgrid_d1)
        use shtools, only: Curve2Mask
        implicit none
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
                        centralmeridian=centralmeridian)
    end subroutine pyCurve2Mask

    subroutine pyMakeEllipseCoord(coord,lat,lon,dec,A_theta,B_theta,cinterval,&
                                  cnum,coord_d0,coord_d1)
        use shtools, only: MakeEllipseCoord
        implicit none
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
                                cinterval=cinterval,cnum=cnum)
    end subroutine pyMakeEllipseCoord

    subroutine pyMakeGravGradGridDH(cilm,lmax,gm,r0,a,f,vxx,vyy,vzz,vxy,vxz, &
                                    vyz,n,sampling,lmax_calc,vyz_d0,vyz_d1, &
                                    vyy_d0,vyy_d1,cilm_d0,cilm_d1,cilm_d2, &
                                    vzz_d0,vzz_d1,vxy_d0,vxy_d1,vxx_d0, &
                                    vxx_d1,vxz_d0,vxz_d1)
        use shtools, only: MakeGravGradGridDH
        implicit none
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
                                sampling=sampling,lmax_calc=lmax_calc)
    end subroutine pyMakeGravGradGridDH
