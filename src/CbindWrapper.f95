    subroutine cbindPlmBar(p,p_d0,lmax,z,csphase,cnorm,exitstatus)  bind(c, name="cbind_plm_bar")
        use, intrinsic :: iso_c_binding
        use shtools, only: PlmBar
        implicit none
        integer(kind=c_int), intent(in) :: lmax
        integer(kind=c_int), intent(in) :: p_d0
        real(kind=c_double), dimension(p_d0),intent(out) :: p
        real(kind=c_double), intent(in) :: z
        integer(kind=c_int), optional,intent(in) :: csphase
        integer(kind=c_int), optional,intent(in) :: cnorm
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call PlmBar(p,lmax,z,csphase=csphase,cnorm=cnorm,exitstatus=exitstatus)
    end subroutine cbindPlmBar

    subroutine cbindPlmBar_d1(p,p_d0,dp1,dp1_d0,lmax,z,csphase,cnorm,exitstatus)  bind(c&
                               , name="cbind_plm_bar_d1")
        use, intrinsic :: iso_c_binding
        use shtools, only: PlmBar_d1
        implicit none
        integer(kind=c_int), intent(in) :: lmax
        integer(kind=c_int), intent(in) :: p_d0
        real(kind=c_double), dimension(p_d0),intent(out) :: p
        integer(kind=c_int), intent(in) :: dp1_d0
        real(kind=c_double), dimension(dp1_d0),intent(out) :: dp1
        real(kind=c_double), intent(in) :: z
        integer(kind=c_int), optional,intent(in) :: csphase
        integer(kind=c_int), optional,intent(in) :: cnorm
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call PlmBar_d1(p,dp1,lmax,z,csphase=csphase,cnorm=cnorm,exitstatus=exitstatus)
    end subroutine cbindPlmBar_d1

    subroutine cbindPlBar(p,p_d0,lmax,z,exitstatus)  bind(c, name="cbind_pl_bar")
        use, intrinsic :: iso_c_binding
        use shtools, only: PlBar
        implicit none
        integer(kind=c_int), intent(in) :: lmax
        integer(kind=c_int), intent(in) :: p_d0
        real(kind=c_double), dimension(p_d0),intent(out) :: p
        real(kind=c_double), intent(in) :: z
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call PlBar(p,lmax,z,exitstatus=exitstatus)
    end subroutine cbindPlBar

    subroutine cbindPlBar_d1(p,p_d0,dp1,dp1_d0,lmax,z,exitstatus)  bind(c, name="cbind_pl_bar_d1")
        use, intrinsic :: iso_c_binding
        use shtools, only: PlBar_d1
        implicit none
        integer(kind=c_int), intent(in) :: lmax
        integer(kind=c_int), intent(in) :: p_d0
        real(kind=c_double), dimension(p_d0),intent(out) :: p
        integer(kind=c_int), intent(in) :: dp1_d0
        real(kind=c_double), dimension(dp1_d0),intent(out) :: dp1
        real(kind=c_double), intent(in) :: z
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call PlBar_d1(p,dp1,lmax,z,exitstatus=exitstatus)
    end subroutine cbindPlBar_d1

    subroutine cbindPlmON(p,p_d0,lmax,z,csphase,cnorm,exitstatus)  bind(c, name="cbind_plm_on")
        use, intrinsic :: iso_c_binding
        use shtools, only: PlmON
        implicit none
        integer(kind=c_int), intent(in) :: lmax
        integer(kind=c_int), intent(in) :: p_d0
        real(kind=c_double), dimension(p_d0),intent(out) :: p
        real(kind=c_double), intent(in) :: z
        integer(kind=c_int), optional,intent(in) :: csphase
        integer(kind=c_int), optional,intent(in) :: cnorm
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call PlmON(p,lmax,z,csphase=csphase,cnorm=cnorm,exitstatus=exitstatus)
    end subroutine cbindPlmON

    subroutine cbindPlmON_d1(p,p_d0,dp1,dp1_d0,lmax,z,csphase,cnorm,exitstatus)  bind(c&
                              , name="cbind_plm_on_d1")
        use, intrinsic :: iso_c_binding
        use shtools, only: PlmON_d1
        implicit none
        integer(kind=c_int), intent(in) :: lmax
        integer(kind=c_int), intent(in) :: p_d0
        real(kind=c_double), dimension(p_d0),intent(out) :: p
        integer(kind=c_int), intent(in) :: dp1_d0
        real(kind=c_double), dimension(dp1_d0),intent(out) :: dp1
        real(kind=c_double), intent(in) :: z
        integer(kind=c_int), optional,intent(in) :: csphase
        integer(kind=c_int), optional,intent(in) :: cnorm
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call PlmON_d1(p,dp1,lmax,z,csphase=csphase,cnorm=cnorm,exitstatus=exitstatus)
    end subroutine cbindPlmON_d1

    subroutine cbindPlON(p,p_d0,lmax,z,exitstatus)  bind(c, name="cbind_pl_on")
        use, intrinsic :: iso_c_binding
        use shtools, only: PlON
        implicit none
        integer(kind=c_int), intent(in) :: lmax
        integer(kind=c_int), intent(in) :: p_d0
        real(kind=c_double), dimension(p_d0),intent(out) :: p
        real(kind=c_double), intent(in) :: z
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call PlON(p,lmax,z,exitstatus=exitstatus)
    end subroutine cbindPlON

    subroutine cbindPlON_d1(p,p_d0,dp1,dp1_d0,lmax,z,exitstatus)  bind(c, name="cbind_pl_on_d1")
        use, intrinsic :: iso_c_binding
        use shtools, only: PlON_d1
        implicit none
        integer(kind=c_int), intent(in) :: lmax
        integer(kind=c_int), intent(in) :: p_d0
        real(kind=c_double), dimension(p_d0),intent(out) :: p
        integer(kind=c_int), intent(in) :: dp1_d0
        real(kind=c_double), dimension(dp1_d0),intent(out) :: dp1
        real(kind=c_double), intent(in) :: z
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call PlON_d1(p,dp1,lmax,z,exitstatus=exitstatus)
    end subroutine cbindPlON_d1

    subroutine cbindPlmSchmidt(p,p_d0,lmax,z,csphase,cnorm,exitstatus)  bind(c, name="cbind_plm_schmidt")
        use, intrinsic :: iso_c_binding
        use shtools, only: PlmSchmidt
        implicit none
        integer(kind=c_int), intent(in) :: lmax
        integer(kind=c_int), intent(in) :: p_d0
        real(kind=c_double), dimension(p_d0),intent(out) :: p
        real(kind=c_double), intent(in) :: z
        integer(kind=c_int), optional,intent(in) :: csphase
        integer(kind=c_int), optional,intent(in) :: cnorm
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call PlmSchmidt(p,lmax,z,csphase=csphase,cnorm=cnorm,exitstatus=exitstatus)
    end subroutine cbindPlmSchmidt

    subroutine cbindPlmSchmidt_d1(p,p_d0,dp1,dp1_d0,lmax,z,csphase,cnorm,exitstatus)  bind(c&
                                   , name="cbind_plm_schmidt_d1")
        use, intrinsic :: iso_c_binding
        use shtools, only: PlmSchmidt_d1
        implicit none
        integer(kind=c_int), intent(in) :: lmax
        integer(kind=c_int), intent(in) :: p_d0
        real(kind=c_double), dimension(p_d0),intent(out) :: p
        integer(kind=c_int), intent(in) :: dp1_d0
        real(kind=c_double), dimension(dp1_d0),intent(out) :: dp1
        real(kind=c_double), intent(in) :: z
        integer(kind=c_int), optional,intent(in) :: csphase
        integer(kind=c_int), optional,intent(in) :: cnorm
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call PlmSchmidt_d1(p,dp1,lmax,z,csphase=csphase,cnorm=cnorm,exitstatus=exitstatus)
    end subroutine cbindPlmSchmidt_d1

    subroutine cbindPlSchmidt(p,p_d0,lmax,z,exitstatus)  bind(c, name="cbind_pl_schmidt")
        use, intrinsic :: iso_c_binding
        use shtools, only: PlSchmidt
        implicit none
        integer(kind=c_int), intent(in) :: lmax
        integer(kind=c_int), intent(in) :: p_d0
        real(kind=c_double), dimension(p_d0),intent(out) :: p
        real(kind=c_double), intent(in) :: z
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call PlSchmidt(p,lmax,z,exitstatus=exitstatus)
    end subroutine cbindPlSchmidt

    subroutine cbindPlSchmidt_d1(p,p_d0,dp1,dp1_d0,lmax,z,exitstatus)  bind(c, name="cbind_pl_schmidt_d1")
        use, intrinsic :: iso_c_binding
        use shtools, only: PlSchmidt_d1
        implicit none
        integer(kind=c_int), intent(in) :: lmax
        integer(kind=c_int), intent(in) :: p_d0
        real(kind=c_double), dimension(p_d0),intent(out) :: p
        integer(kind=c_int), intent(in) :: dp1_d0
        real(kind=c_double), dimension(dp1_d0),intent(out) :: dp1
        real(kind=c_double), intent(in) :: z
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call PlSchmidt_d1(p,dp1,lmax,z,exitstatus=exitstatus)
    end subroutine cbindPlSchmidt_d1

    subroutine cbindPLegendreA(p,p_d0,lmax,z,csphase,exitstatus)  bind(c, name="cbind_p_legendre_a")
        use, intrinsic :: iso_c_binding
        use shtools, only: PLegendreA
        implicit none
        integer(kind=c_int), intent(in) :: lmax
        integer(kind=c_int), intent(in) :: p_d0
        real(kind=c_double), dimension(p_d0),intent(out) :: p
        real(kind=c_double), intent(in) :: z
        integer(kind=c_int), optional,intent(in) :: csphase
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call PLegendreA(p,lmax,z,csphase=csphase,exitstatus=exitstatus)
    end subroutine cbindPLegendreA

    subroutine cbindPLegendreA_d1(p,p_d0,dp1,dp1_d0,lmax,z,csphase,exitstatus)  bind(c&
                                   , name="cbind_p_legendre_a_d1")
        use, intrinsic :: iso_c_binding
        use shtools, only: PLegendreA_d1
        implicit none
        integer(kind=c_int), intent(in) :: lmax
        integer(kind=c_int), intent(in) :: p_d0
        real(kind=c_double), dimension(p_d0),intent(out) :: p
        integer(kind=c_int), intent(in) :: dp1_d0
        real(kind=c_double), dimension(dp1_d0),intent(out) :: dp1
        real(kind=c_double), intent(in) :: z
        integer(kind=c_int), optional,intent(in) :: csphase
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call PLegendreA_d1(p,dp1,lmax,z,csphase=csphase,exitstatus=exitstatus)
    end subroutine cbindPLegendreA_d1

    subroutine cbindPLegendre(p,p_d0,lmax,z,exitstatus)  bind(c, name="cbind_p_legendre")
        use, intrinsic :: iso_c_binding
        use shtools, only: PLegendre
        implicit none
        integer(kind=c_int), intent(in) :: lmax
        integer(kind=c_int), intent(in) :: p_d0
        real(kind=c_double), dimension(p_d0),intent(out) :: p
        real(kind=c_double), intent(in) :: z
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call PLegendre(p,lmax,z,exitstatus=exitstatus)
    end subroutine cbindPLegendre

    subroutine cbindPLegendre_d1(p,p_d0,dp1,dp1_d0,lmax,z,exitstatus)  bind(c, name="cbind_p_legendre_d1")
        use, intrinsic :: iso_c_binding
        use shtools, only: PLegendre_d1
        implicit none
        integer(kind=c_int), intent(in) :: lmax
        integer(kind=c_int), intent(in) :: p_d0
        real(kind=c_double), dimension(p_d0),intent(out) :: p
        integer(kind=c_int), intent(in) :: dp1_d0
        real(kind=c_double), dimension(dp1_d0),intent(out) :: dp1
        real(kind=c_double), intent(in) :: z
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call PLegendre_d1(p,dp1,lmax,z,exitstatus=exitstatus)
    end subroutine cbindPLegendre_d1

    function cbindPlmIndex(l,m)  bind(c, name="cbind_plm_index")
        use, intrinsic :: iso_c_binding
        use shtools, only: PlmIndex
        implicit none
        integer(kind=c_int) :: cbindPlmIndex
        integer(kind=c_int), intent(in) :: l
        integer(kind=c_int), intent(in) :: m
        cbindPlmIndex=PlmIndex(l,m)
    end function cbindPlmIndex

    subroutine cbindSHExpandDH(grid,grid_d0,grid_d1,n,cilm,cilm_d,lmax,norm,sampling&
                                   ,csphase,lmax_calc,exitstatus)  bind(c, name="cbind_sh_expand_dh")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHExpandDH
        implicit none
        integer(kind=c_int), intent(in) :: grid_d0
        integer(kind=c_int), intent(in) :: grid_d1
        real(kind=c_double), dimension(grid_d0,grid_d1),intent(in) :: grid
        integer(kind=c_int), intent(in) :: cilm_d
        real(kind=c_double), dimension(2,cilm_d,cilm_d),intent(out) :: cilm
        integer(kind=c_int), intent(in) :: n
        integer(kind=c_int), intent(out) :: lmax
        integer(kind=c_int), optional,intent(in) :: norm
        integer(kind=c_int), optional,intent(in) :: sampling
        integer(kind=c_int), optional,intent(in) :: csphase
        integer(kind=c_int), optional,intent(in) :: lmax_calc
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call SHExpandDH(grid,n,cilm,lmax,norm=norm,sampling=sampling,csphase=csphase&
                            ,lmax_calc=lmax_calc,exitstatus=exitstatus)
    end subroutine cbindSHExpandDH

    subroutine cbindMakeGridDH(griddh,griddh_d0,griddh_d1,n,cilm,cilm_d,lmax,norm&
                                     ,sampling,csphase,lmax_calc,extend,exitstatus)  bind(c&
                                     , name="cbind_make_grid_dh")
        use, intrinsic :: iso_c_binding
        use shtools, only: MakeGridDH
        implicit none
        integer(kind=c_int), intent(in) :: cilm_d
        real(kind=c_double), dimension(2,cilm_d,cilm_d),intent(in) :: cilm
        integer(kind=c_int), intent(in) :: griddh_d0
        integer(kind=c_int), intent(in) :: griddh_d1
        real(kind=c_double), dimension(griddh_d0,griddh_d1),intent(out) :: griddh
        integer(kind=c_int), intent(in) :: lmax
        integer(kind=c_int), intent(out) :: n
        integer(kind=c_int), optional,intent(in) :: norm
        integer(kind=c_int), optional,intent(in) :: sampling
        integer(kind=c_int), optional,intent(in) :: csphase
        integer(kind=c_int), optional,intent(in) :: lmax_calc
        integer(kind=c_int), optional,intent(in) :: extend
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call MakeGridDH(griddh,n,cilm,lmax,norm=norm,sampling=sampling,csphase=csphase&
                              ,lmax_calc=lmax_calc,extend=extend,exitstatus=exitstatus)
    end subroutine cbindMakeGridDH

    subroutine cbindSHExpandDHC(grid,grid_d0,grid_d1,n,cilm,cilm_d,lmax,norm,sampling&
                                    ,csphase,lmax_calc,exitstatus)  bind(c, name="cbind_sh_expand_dhc")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHExpandDHC
        implicit none
        integer(kind=c_int), intent(in) :: grid_d0
        integer(kind=c_int), intent(in) :: grid_d1
        complex(kind=c_double_complex), dimension(grid_d0,grid_d1),intent(in) :: grid
        integer(kind=c_int), intent(in) :: cilm_d
        complex(kind=c_double_complex), dimension(2,cilm_d,cilm_d),intent(out) :: cilm
        integer(kind=c_int), intent(in) :: n
        integer(kind=c_int), intent(out) :: lmax
        integer(kind=c_int), optional,intent(in) :: norm
        integer(kind=c_int), optional,intent(in) :: sampling
        integer(kind=c_int), optional,intent(in) :: csphase
        integer(kind=c_int), optional,intent(in) :: lmax_calc
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call SHExpandDHC(grid,n,cilm,lmax,norm=norm,sampling=sampling,csphase=csphase&
                             ,lmax_calc=lmax_calc,exitstatus=exitstatus)
    end subroutine cbindSHExpandDHC

    subroutine cbindMakeGridDHC(griddh,griddh_d0,griddh_d1,n,cilm,cilm_d,lmax,norm&
                                      ,sampling,csphase,lmax_calc,extend,exitstatus)  bind(c&
                                      , name="cbind_make_grid_dhc")
        use, intrinsic :: iso_c_binding
        use shtools, only: MakeGridDHC
        implicit none
        integer(kind=c_int), intent(in) :: cilm_d
        complex(kind=c_double_complex), dimension(2,cilm_d,cilm_d),intent(in) :: cilm
        integer(kind=c_int), intent(in) :: griddh_d0
        integer(kind=c_int), intent(in) :: griddh_d1
        complex(kind=c_double_complex), dimension(griddh_d0,griddh_d1),intent(out) :: griddh
        integer(kind=c_int), intent(in) :: lmax
        integer(kind=c_int), intent(out) :: n
        integer(kind=c_int), optional,intent(in) :: norm
        integer(kind=c_int), optional,intent(in) :: sampling
        integer(kind=c_int), optional,intent(in) :: csphase
        integer(kind=c_int), optional,intent(in) :: lmax_calc
        integer(kind=c_int), optional,intent(in) :: extend
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call MakeGridDHC(griddh,n,cilm,lmax,norm=norm,sampling=sampling,csphase=csphase&
                               ,lmax_calc=lmax_calc,extend=extend,exitstatus=exitstatus)
    end subroutine cbindMakeGridDHC

    subroutine cbindSHGLQ(lmax,zero,zero_d0,w,w_d0,plx,plx_d0,plx_d1,norm,csphase&
                              ,cnorm,exitstatus)  bind(c, name="cbind_shglq")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHGLQ
        implicit none
        integer(kind=c_int), intent(in) :: lmax
        integer(kind=c_int), intent(in) :: zero_d0
        real(kind=c_double), dimension(zero_d0),intent(out) :: zero
        integer(kind=c_int), intent(in) :: w_d0
        real(kind=c_double), dimension(w_d0),intent(out) :: w
        integer(kind=c_int), intent(in) :: plx_d0
        integer(kind=c_int), intent(in) :: plx_d1
        real(kind=c_double), optional,dimension(plx_d0,plx_d1),intent(out) :: plx
        integer(kind=c_int), optional,intent(in) :: norm
        integer(kind=c_int), optional,intent(in) :: csphase
        integer(kind=c_int), optional,intent(in) :: cnorm
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call SHGLQ(lmax,zero,w,plx=plx,norm=norm,csphase=csphase,cnorm=cnorm,exitstatus=exitstatus)
    end subroutine cbindSHGLQ

    subroutine cbindSHExpandGLQ(cilm,cilm_d,lmax,gridglq,gridglq_d0,gridglq_d1,w,w_d0&
                                    ,plx,plx_d0,plx_d1,zero,zero_d0,norm,csphase,lmax_calc&
                                    ,exitstatus)  bind(c, name="cbind_sh_expand_glq")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHExpandGLQ
        implicit none
        integer(kind=c_int), intent(in) :: w_d0
        real(kind=c_double), dimension(w_d0),intent(in) :: w
        integer(kind=c_int), intent(in) :: gridglq_d0
        integer(kind=c_int), intent(in) :: gridglq_d1
        real(kind=c_double), dimension(gridglq_d0,gridglq_d1),intent(in) :: gridglq
        integer(kind=c_int), intent(in) :: plx_d0
        integer(kind=c_int), intent(in) :: plx_d1
        real(kind=c_double), optional,dimension(plx_d0,plx_d1),intent(in) :: plx
        integer(kind=c_int), intent(in) :: zero_d0
        real(kind=c_double), optional,dimension(zero_d0),intent(in) :: zero
        integer(kind=c_int), intent(in) :: cilm_d
        real(kind=c_double), dimension(2,cilm_d,cilm_d),intent(out) :: cilm
        integer(kind=c_int), intent(in) :: lmax
        integer(kind=c_int), optional,intent(in) :: norm
        integer(kind=c_int), optional,intent(in) :: csphase
        integer(kind=c_int), optional,intent(in) :: lmax_calc
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call SHExpandGLQ(cilm,lmax,gridglq,w,plx=plx,zero=zero,norm=norm,csphase=csphase&
                             ,lmax_calc=lmax_calc,exitstatus=exitstatus)
    end subroutine cbindSHExpandGLQ

    subroutine cbindMakeGridGLQ(gridglq,gridglq_d0,gridglq_d1,cilm,cilm_d,lmax,plx&
                                       ,plx_d0,plx_d1,zero,zero_d0,norm,csphase,lmax_calc&
                                       ,extend,exitstatus)  bind(c, name="cbind_make_grid_glq")
        use, intrinsic :: iso_c_binding
        use shtools, only: MakeGridGLQ
        implicit none
        integer(kind=c_int), intent(in) :: cilm_d
        real(kind=c_double), dimension(2,cilm_d,cilm_d),intent(in) :: cilm
        integer(kind=c_int), intent(in) :: plx_d0
        integer(kind=c_int), intent(in) :: plx_d1
        real(kind=c_double), optional,dimension(plx_d0,plx_d1),intent(in) :: plx
        integer(kind=c_int), intent(in) :: zero_d0
        real(kind=c_double), optional,dimension(zero_d0),intent(in) :: zero
        integer(kind=c_int), intent(in) :: gridglq_d0
        integer(kind=c_int), intent(in) :: gridglq_d1
        real(kind=c_double), dimension(gridglq_d0,gridglq_d1),intent(out) :: gridglq
        integer(kind=c_int), intent(in) :: lmax
        integer(kind=c_int), optional,intent(in) :: norm
        integer(kind=c_int), optional,intent(in) :: csphase
        integer(kind=c_int), optional,intent(in) :: lmax_calc
        integer(kind=c_int), optional,intent(in) :: extend
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call MakeGridGLQ(gridglq,cilm,lmax,plx=plx,zero=zero,norm=norm,csphase=csphase&
                                ,lmax_calc=lmax_calc,extend=extend,exitstatus=exitstatus)
    end subroutine cbindMakeGridGLQ

    subroutine cbindSHExpandGLQC(cilm,cilm_d,lmax,gridglq,gridglq_d0,gridglq_d1,w&
                                     ,w_d0,plx,plx_d0,plx_d1,zero,zero_d0,norm,csphase&
                                     ,lmax_calc,exitstatus)  bind(c, name="cbind_sh_expand_glqc")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHExpandGLQC
        implicit none
        integer(kind=c_int), intent(in) :: w_d0
        real(kind=c_double), dimension(w_d0),intent(in) :: w
        integer(kind=c_int), intent(in) :: gridglq_d0
        integer(kind=c_int), intent(in) :: gridglq_d1
        complex(kind=c_double_complex), dimension(gridglq_d0,gridglq_d1),intent(in) :: gridglq
        integer(kind=c_int), intent(in) :: plx_d0
        integer(kind=c_int), intent(in) :: plx_d1
        real(kind=c_double), optional,dimension(plx_d0,plx_d1),intent(in) :: plx
        integer(kind=c_int), intent(in) :: zero_d0
        real(kind=c_double), optional,dimension(zero_d0),intent(in) :: zero
        integer(kind=c_int), intent(in) :: cilm_d
        complex(kind=c_double_complex), dimension(2,cilm_d,cilm_d),intent(out) :: cilm
        integer(kind=c_int), intent(in) :: lmax
        integer(kind=c_int), optional,intent(in) :: norm
        integer(kind=c_int), optional,intent(in) :: csphase
        integer(kind=c_int), optional,intent(in) :: lmax_calc
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call SHExpandGLQC(cilm,lmax,gridglq,w,plx=plx,zero=zero,norm=norm,csphase=csphase&
                              ,lmax_calc=lmax_calc,exitstatus=exitstatus)
    end subroutine cbindSHExpandGLQC

    subroutine cbindMakeGridGLQC(gridglq,gridglq_d0,gridglq_d1,cilm,cilm_d,lmax,plx&
                                        ,plx_d0,plx_d1,zero,zero_d0,norm,csphase,lmax_calc&
                                        ,extend,exitstatus)  bind(c, name="cbind_make_grid_glqc")
        use, intrinsic :: iso_c_binding
        use shtools, only: MakeGridGLQC
        implicit none
        integer(kind=c_int), intent(in) :: cilm_d
        complex(kind=c_double_complex), dimension(2,cilm_d,cilm_d),intent(in) :: cilm
        integer(kind=c_int), intent(in) :: plx_d0
        integer(kind=c_int), intent(in) :: plx_d1
        real(kind=c_double), optional,dimension(plx_d0,plx_d1),intent(in) :: plx
        integer(kind=c_int), intent(in) :: zero_d0
        real(kind=c_double), optional,dimension(zero_d0),intent(in) :: zero
        integer(kind=c_int), intent(in) :: gridglq_d0
        integer(kind=c_int), intent(in) :: gridglq_d1
        complex(kind=c_double_complex), dimension(gridglq_d0,gridglq_d1),intent(out) :: gridglq
        integer(kind=c_int), intent(in) :: lmax
        integer(kind=c_int), optional,intent(in) :: norm
        integer(kind=c_int), optional,intent(in) :: csphase
        integer(kind=c_int), optional,intent(in) :: lmax_calc
        integer(kind=c_int), optional,intent(in) :: extend
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call MakeGridGLQC(gridglq,cilm,lmax,plx=plx,zero=zero,norm=norm,csphase=csphase&
                                 ,lmax_calc=lmax_calc,extend=extend,exitstatus=exitstatus)
    end subroutine cbindMakeGridGLQC

    subroutine cbindGLQGridCoord(latglq,latglq_d0,longlq,longlq_d0,lmax,nlat,nlong&
                                       ,extend,exitstatus)  bind(c, name="cbind_glq_grid_coord")
        use, intrinsic :: iso_c_binding
        use shtools, only: GLQGridCoord
        implicit none
        integer(kind=c_int), intent(in) :: lmax
        integer(kind=c_int), intent(out) :: nlat
        integer(kind=c_int), intent(out) :: nlong
        integer(kind=c_int), intent(in) :: latglq_d0
        real(kind=c_double), dimension(latglq_d0),intent(out) :: latglq
        integer(kind=c_int), intent(in) :: longlq_d0
        real(kind=c_double), dimension(longlq_d0),intent(out) :: longlq
        integer(kind=c_int), optional,intent(in) :: extend
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call GLQGridCoord(latglq,longlq,lmax,nlat,nlong,extend=extend,exitstatus=exitstatus)
    end subroutine cbindGLQGridCoord

    subroutine cbindSHExpandLSQ(cilm,cilm_d,d,d_d0,lat,lat_d0,lon,lon_d0,nmax,lmax&
                                    ,norm,chi2,csphase,weights,weights_d0,exitstatus)  bind(c&
                                    , name="cbind_sh_expand_lsq")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHExpandLSQ
        implicit none
        integer(kind=c_int), intent(in) :: d_d0
        real(kind=c_double), dimension(d_d0),intent(in) :: d
        integer(kind=c_int), intent(in) :: lat_d0
        real(kind=c_double), dimension(lat_d0),intent(in) :: lat
        integer(kind=c_int), intent(in) :: lon_d0
        real(kind=c_double), dimension(lon_d0),intent(in) :: lon
        integer(kind=c_int), intent(in) :: cilm_d
        real(kind=c_double), dimension(2,cilm_d,cilm_d),intent(out) :: cilm
        integer(kind=c_int), intent(in) :: nmax
        integer(kind=c_int), intent(in) :: lmax
        integer(kind=c_int), optional,intent(in) :: norm
        integer(kind=c_int), optional,intent(in) :: csphase
        real(kind=c_double), optional,intent(out) :: chi2
        integer(kind=c_int), intent(in) :: weights_d0
        real(kind=c_double), optional,dimension(weights_d0),intent(in) :: weights
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call SHExpandLSQ(cilm,d,lat,lon,nmax,lmax,norm=norm,chi2=chi2,csphase=csphase&
                             ,weights=weights,exitstatus=exitstatus)
    end subroutine cbindSHExpandLSQ

    subroutine cbindMakeGrid2d(grid,grid_d0,grid_d1,cilm,cilm_d,lmax,interval,nlat&
                                   ,nlong,norm,csphase,f,a,north,south,east,west,dealloc&
                                   ,exitstatus)  bind(c, name="cbind_make_grid2d")
        use, intrinsic :: iso_c_binding
        use shtools, only: MakeGrid2d
        implicit none
        integer(kind=c_int), intent(in) :: cilm_d
        real(kind=c_double), dimension(2,cilm_d,cilm_d),intent(in) :: cilm
        real(kind=c_double), intent(in) :: interval
        integer(kind=c_int), intent(in) :: grid_d0
        integer(kind=c_int), intent(in) :: grid_d1
        real(kind=c_double), dimension(grid_d0,grid_d1),intent(out) :: grid
        integer(kind=c_int), intent(in) :: lmax
        integer(kind=c_int), intent(out) :: nlat
        integer(kind=c_int), intent(out) :: nlong
        integer(kind=c_int), optional,intent(in) :: norm
        integer(kind=c_int), optional,intent(in) :: csphase
        integer(kind=c_int), optional,intent(in) :: dealloc
        real(kind=c_double), optional,intent(in) :: f
        real(kind=c_double), optional,intent(in) :: a
        real(kind=c_double), optional,intent(in) :: north
        real(kind=c_double), optional,intent(in) :: south
        real(kind=c_double), optional,intent(in) :: east
        real(kind=c_double), optional,intent(in) :: west
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call MakeGrid2d(grid,cilm,lmax,interval,nlat,nlong,norm=norm,csphase=csphase&
                            ,f=f,a=a,north=north,south=south,east=east,west=west,dealloc=dealloc&
                            ,exitstatus=exitstatus)
    end subroutine cbindMakeGrid2d

    function cbindMakeGridPoint(cilm,cilm_d,lmax,lat,lon,norm,csphase,dealloc)  bind(c&
                                    , name="cbind_make_grid_point")
        use, intrinsic :: iso_c_binding
        use shtools, only: MakeGridPoint
        implicit none
        real(kind=c_double) :: cbindMakeGridPoint
        integer(kind=c_int), intent(in) :: cilm_d
        real(kind=c_double), dimension(2,cilm_d,cilm_d),intent(in) :: cilm
        real(kind=c_double), intent(in) :: lat
        real(kind=c_double), intent(in) :: lon
        integer(kind=c_int), intent(in) :: lmax
        integer(kind=c_int), optional,intent(in) :: norm
        integer(kind=c_int), optional,intent(in) :: csphase
        integer(kind=c_int), optional,intent(in) :: dealloc
        cbindMakeGridPoint=MakeGridPoint(cilm,lmax,lat,lon,norm=norm,csphase=csphase&
                                             ,dealloc=dealloc)
    end function cbindMakeGridPoint

    function cbindMakeGridPointC(cilm,cilm_d,lmax,lat,lon,norm,csphase,dealloc)  bind(c&
                                     , name="cbind_make_grid_point_c")
        use, intrinsic :: iso_c_binding
        use shtools, only: MakeGridPointC
        implicit none
        complex(kind=c_double_complex) :: cbindMakeGridPointC
        integer(kind=c_int), intent(in) :: cilm_d
        complex(kind=c_double_complex), dimension(2,cilm_d,cilm_d),intent(in) :: cilm
        real(kind=c_double), intent(in) :: lat
        real(kind=c_double), intent(in) :: lon
        integer(kind=c_int), intent(in) :: lmax
        integer(kind=c_int), optional,intent(in) :: norm
        integer(kind=c_int), optional,intent(in) :: csphase
        integer(kind=c_int), optional,intent(in) :: dealloc
        cbindMakeGridPointC=MakeGridPointC(cilm,lmax,lat,lon,norm=norm,csphase=csphase&
                                               ,dealloc=dealloc)
    end function cbindMakeGridPointC

    subroutine cbindSHMultiply(shout,shout_d0,shout_d1,shout_d2,sh1,sh1_d0,sh1_d1&
                                    ,sh1_d2,lmax1,sh2,sh2_d0,sh2_d1,sh2_d2,lmax2,precomp&
                                    ,norm,csphase,exitstatus)  bind(c, name="cbind_sh_multiply")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHMultiply
        implicit none
        integer(kind=c_int), intent(in) :: shout_d0
        integer(kind=c_int), intent(in) :: shout_d1
        integer(kind=c_int), intent(in) :: shout_d2
        real(kind=c_double), dimension(shout_d0,shout_d1,shout_d2),intent(out) :: shout
        integer(kind=c_int), intent(in) :: sh1_d0
        integer(kind=c_int), intent(in) :: sh1_d1
        integer(kind=c_int), intent(in) :: sh1_d2
        real(kind=c_double), dimension(sh1_d0,sh1_d1,sh1_d2),intent(in) :: sh1
        integer(kind=c_int), intent(in) :: sh2_d0
        integer(kind=c_int), intent(in) :: sh2_d1
        integer(kind=c_int), intent(in) :: sh2_d2
        real(kind=c_double), dimension(sh2_d0,sh2_d1,sh2_d2),intent(in) :: sh2
        integer(kind=c_int), intent(in) :: lmax1
        integer(kind=c_int), intent(in) :: lmax2
        integer(kind=c_int), optional,intent(in) :: precomp
        integer(kind=c_int), optional,intent(in) :: norm
        integer(kind=c_int), optional,intent(in) :: csphase
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call SHMultiply(shout,sh1,lmax1,sh2,lmax2,precomp=precomp,norm=norm,csphase=csphase&
                             ,exitstatus=exitstatus)
    end subroutine cbindSHMultiply

    subroutine cbindSHRead(filename,filename_d1,cilm,cilm_d,lmax,skip,header,header_d0&
                                   ,error,error_d0,error_d1,error_d2,exitstatus)  bind(c&
                                   , name="cbind_sh_read")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHRead
        implicit none
        integer(kind=c_int), intent(in) :: filename_d1
        character(kind=c_char), dimension(filename_d1),intent(in) :: filename
        integer(kind=c_int), intent(out) :: lmax
        integer(kind=c_int), intent(in) :: cilm_d
        real(kind=c_double), dimension(2,cilm_d,cilm_d),intent(out) :: cilm
        integer(kind=c_int), intent(in) :: header_d0
        real(kind=c_double), optional,dimension(header_d0),intent(out) :: header
        integer(kind=c_int), intent(in) :: error_d0
        integer(kind=c_int), intent(in) :: error_d1
        integer(kind=c_int), intent(in) :: error_d2
        real(kind=c_double), optional,dimension(error_d0,error_d1,error_d2),intent(out) :: error
        integer(kind=c_int), optional,intent(in) :: skip
        integer(kind=c_int), optional,intent(out) :: exitstatus
        
        character(filename_d1) :: filename2
        filename2 = TRANSFER(filename,filename2)
        
        call SHRead(filename2,cilm,lmax,skip=skip,header=header,error=error,exitstatus=exitstatus)
    end subroutine cbindSHRead

    subroutine cbindSHRead2(filename,filename_d1,cilm,cilm_d,lmax,gm,r0_pot,error&
                                    ,error_d0,error_d1,error_d2,dot,dot_d0,dot_d1&
                                    ,dot_d2,doystart,doyend,epoch,exitstatus)  bind(c&
                                    , name="cbind_sh_read2")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHRead2
        implicit none
        integer(kind=c_int), intent(in) :: filename_d1
        character(kind=c_char), dimension(filename_d1),intent(in) :: filename
        integer(kind=c_int), intent(out) :: lmax
        integer(kind=c_int), intent(in) :: cilm_d
        real(kind=c_double), dimension(2,cilm_d,cilm_d),intent(out) :: cilm
        real(kind=c_double), intent(out) :: gm
        real(kind=c_double), intent(out) :: r0_pot
        integer(kind=c_int), intent(in) :: error_d0
        integer(kind=c_int), intent(in) :: error_d1
        integer(kind=c_int), intent(in) :: error_d2
        real(kind=c_double), optional,dimension(error_d0,error_d1,error_d2),intent(out) :: error
        integer(kind=c_int), intent(in) :: dot_d0
        integer(kind=c_int), intent(in) :: dot_d1
        integer(kind=c_int), intent(in) :: dot_d2
        real(kind=c_double), optional,dimension(dot_d0,dot_d1,dot_d2),intent(out) :: dot
        real(kind=c_double), optional,intent(out) :: doystart
        real(kind=c_double), optional,intent(out) :: doyend
        real(kind=c_double), optional,intent(out) :: epoch
        integer(kind=c_int), optional,intent(out) :: exitstatus
        character(filename_d1) :: filename2
        filename2 = TRANSFER(filename,filename2)
        
        call SHRead2(filename2,cilm,lmax,gm,r0_pot,error=error,dot=dot,doystart=doystart&
                             ,doyend=doyend,epoch=epoch,exitstatus=exitstatus)
    end subroutine cbindSHRead2

    subroutine cbindSHReadJPL(filename,filename_d1,cilm,cilm_d,lmax,error,error_d0&
                                      ,error_d1,error_d2,gm,formatstring,exitstatus)  bind(c&
                                      , name="cbind_sh_read_jpl")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHReadJPL
        implicit none
        integer(kind=c_int), intent(in) :: filename_d1
        character(kind=c_char), dimension(filename_d1),intent(in) :: filename
        integer(kind=c_int), intent(in) :: lmax
        integer(kind=c_int), intent(in) :: cilm_d
        real(kind=c_double), dimension(2,cilm_d,cilm_d),intent(out) :: cilm
        integer(kind=c_int), intent(in) :: error_d0
        integer(kind=c_int), intent(in) :: error_d1
        integer(kind=c_int), intent(in) :: error_d2
        real(kind=c_double), optional,dimension(error_d0,error_d1,error_d2),intent(out) :: error
        real(kind=c_double), optional,dimension(2),intent(out) :: gm
        character(kind=c_char), optional,intent(in) :: formatstring
        integer(kind=c_int), optional,intent(out) :: exitstatus
        
        character(filename_d1) :: filename2
        character(6) :: formatstring2
        
        filename2 = TRANSFER(filename,filename2)
        formatstring2 = TRANSFER(formatstring, formatstring2)
        
        call SHReadJPL(filename2,cilm,lmax,error=error,gm=gm,formatstring=formatstring2&
                               ,exitstatus=exitstatus)
    end subroutine cbindSHReadJPL

    subroutine cbindSHCilmToVector(cilm,cilm_d,vector,vector_d0,lmax,exitstatus)  bind(c&
                                       , name="cbind_sh_cilm_to_vector")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHCilmToVector
        implicit none
        integer(kind=c_int), intent(in) :: cilm_d
        real(kind=c_double), dimension(2,cilm_d,cilm_d),intent(in) :: cilm
        integer(kind=c_int), intent(in) :: vector_d0
        real(kind=c_double), dimension(vector_d0),intent(out) :: vector
        integer(kind=c_int), intent(in) :: lmax
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call SHCilmToVector(cilm,vector,lmax,exitstatus=exitstatus)
    end subroutine cbindSHCilmToVector

    subroutine cbindSHVectorToCilm(vector,vector_d0,cilm,cilm_d,lmax,exitstatus)  bind(c&
                                         , name="cbind_sh_vector_to_cilm")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHVectorToCilm
        implicit none
        integer(kind=c_int), intent(in) :: cilm_d
        real(kind=c_double), dimension(2,cilm_d,cilm_d),intent(out) :: cilm
        integer(kind=c_int), intent(in) :: vector_d0
        real(kind=c_double), dimension(vector_d0),intent(in) :: vector
        integer(kind=c_int), intent(in) :: lmax
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call SHVectorToCilm(vector,cilm,lmax,exitstatus=exitstatus)
    end subroutine cbindSHVectorToCilm

    subroutine cbindSHCilmToCindex(cilm,cilm_d,cindex,cindex_d0,cindex_d1,degmax,exitstatus)  bind(c&
                                       , name="cbind_sh_cilm_to_cindex")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHCilmToCindex
        implicit none
        integer(kind=c_int), intent(in) :: cilm_d
        real(kind=c_double), dimension(2,cilm_d,cilm_d),intent(in) :: cilm
        integer(kind=c_int), intent(in) :: cindex_d0
        integer(kind=c_int), intent(in) :: cindex_d1
        real(kind=c_double), dimension(cindex_d0,cindex_d1),intent(out) :: cindex
        integer(kind=c_int), optional,intent(in) :: degmax
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call SHCilmToCindex(cilm,cindex,degmax=degmax,exitstatus=exitstatus)
    end subroutine cbindSHCilmToCindex

    subroutine cbindSHCindexToCilm(cindex,cindex_d0,cindex_d1,cilm,cilm_d,degmax,exitstatus)  bind(c&
                                         , name="cbind_sh_cindex_to_cilm")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHCindexToCilm
        implicit none
        integer(kind=c_int), intent(in) :: cilm_d
        real(kind=c_double), dimension(2,cilm_d,cilm_d),intent(out) :: cilm
        integer(kind=c_int), intent(in) :: cindex_d0
        integer(kind=c_int), intent(in) :: cindex_d1
        real(kind=c_double), dimension(cindex_d0,cindex_d1),intent(in) :: cindex
        integer(kind=c_int), optional,intent(in) :: degmax
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call SHCindexToCilm(cindex,cilm,degmax=degmax,exitstatus=exitstatus)
    end subroutine cbindSHCindexToCilm

    subroutine cbindSHrtoc(rcilm,rcilm_d0,rcilm_d1,rcilm_d2,ccilm,ccilm_d0,ccilm_d1&
                                ,ccilm_d2,degmax,convention,switchcs,exitstatus)  bind(c&
                                , name="cbind_s_hrtoc")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHrtoc
        implicit none
        integer(kind=c_int), intent(in) :: rcilm_d0
        integer(kind=c_int), intent(in) :: rcilm_d1
        integer(kind=c_int), intent(in) :: rcilm_d2
        real(kind=c_double), dimension(rcilm_d0,rcilm_d1,rcilm_d2),intent(in) :: rcilm
        integer(kind=c_int), intent(in) :: ccilm_d0
        integer(kind=c_int), intent(in) :: ccilm_d1
        integer(kind=c_int), intent(in) :: ccilm_d2
        real(kind=c_double), dimension(ccilm_d0,ccilm_d1,ccilm_d2),intent(out) :: ccilm
        integer(kind=c_int), optional,intent(in) :: degmax
        integer(kind=c_int), optional,intent(in) :: convention
        integer(kind=c_int), optional,intent(in) :: switchcs
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call SHrtoc(rcilm,ccilm,degmax=degmax,convention=convention,switchcs=switchcs&
                         ,exitstatus=exitstatus)
    end subroutine cbindSHrtoc

    subroutine cbindSHctor(ccilm,ccilm_d0,ccilm_d1,ccilm_d2,rcilm,rcilm_d0,rcilm_d1&
                                ,rcilm_d2,degmax,convention,switchcs,exitstatus)  bind(c&
                                , name="cbind_s_hctor")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHctor
        implicit none
        integer(kind=c_int), intent(in) :: ccilm_d0
        integer(kind=c_int), intent(in) :: ccilm_d1
        integer(kind=c_int), intent(in) :: ccilm_d2
        real(kind=c_double), dimension(ccilm_d0,ccilm_d1,ccilm_d2),intent(in) :: ccilm
        integer(kind=c_int), intent(in) :: rcilm_d0
        integer(kind=c_int), intent(in) :: rcilm_d1
        integer(kind=c_int), intent(in) :: rcilm_d2
        real(kind=c_double), dimension(rcilm_d0,rcilm_d1,rcilm_d2),intent(out) :: rcilm
        integer(kind=c_int), optional,intent(in) :: degmax
        integer(kind=c_int), optional,intent(in) :: convention
        integer(kind=c_int), optional,intent(in) :: switchcs
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call SHctor(ccilm,rcilm,degmax=degmax,convention=convention,switchcs=switchcs&
                         ,exitstatus=exitstatus)
    end subroutine cbindSHctor

    subroutine cbinddjpi2(dj,dj_d0,dj_d1,dj_d2,lmax,exitstatus)  bind(c, name="cbinddjpi2")
        use, intrinsic :: iso_c_binding
        use shtools, only: djpi2
        implicit none
        integer(kind=c_int), intent(in) :: lmax
        integer(kind=c_int), intent(in) :: dj_d0
        integer(kind=c_int), intent(in) :: dj_d1
        integer(kind=c_int), intent(in) :: dj_d2
        real(kind=c_double), dimension(dj_d0,dj_d1,dj_d2),intent(out) :: dj
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call djpi2(dj,lmax,exitstatus=exitstatus)
    end subroutine cbinddjpi2

    subroutine cbindSHRotateCoef(x,cof,cof_d0,cof_d1,rcof,rcof_d0,rcof_d1,dj,dj_d0&
                                  ,dj_d1,dj_d2,lmax,exitstatus)  bind(c, name="cbind_sh_rotate_coef")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHRotateCoef
        implicit none
        integer(kind=c_int), intent(in) :: cof_d0
        integer(kind=c_int), intent(in) :: cof_d1
        real(kind=c_double), dimension(cof_d0,cof_d1),intent(in) :: cof
        integer(kind=c_int), intent(in) :: dj_d0
        integer(kind=c_int), intent(in) :: dj_d1
        integer(kind=c_int), intent(in) :: dj_d2
        real(kind=c_double), dimension(dj_d0,dj_d1,dj_d2),intent(in) :: dj
        real(kind=c_double), dimension(3),intent(in) :: x
        integer(kind=c_int), intent(in) :: rcof_d0
        integer(kind=c_int), intent(in) :: rcof_d1
        real(kind=c_double), dimension(rcof_d0,rcof_d1),intent(out) :: rcof
        integer(kind=c_int), intent(in) :: lmax
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call SHRotateCoef(x,cof,rcof,dj,lmax,exitstatus=exitstatus)
    end subroutine cbindSHRotateCoef

    subroutine cbindSHRotateRealCoef(cilmrot,cilmrot_d0,cilmrot_d1,cilmrot_d2,cilm&
                                            ,cilm_d,lmax,x,x_d0,dj,dj_d0,dj_d1,dj_d2&
                                            ,exitstatus)  bind(c, name="cbind_sh_rotate_real_coef")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHRotateRealCoef
        implicit none
        integer(kind=c_int), intent(in) :: cilm_d
        real(kind=c_double), dimension(2,cilm_d,cilm_d),intent(in) :: cilm
        integer(kind=c_int), intent(in) :: x_d0
        real(kind=c_double), dimension(x_d0),intent(in) :: x
        integer(kind=c_int), intent(in) :: dj_d0
        integer(kind=c_int), intent(in) :: dj_d1
        integer(kind=c_int), intent(in) :: dj_d2
        real(kind=c_double), dimension(dj_d0,dj_d1,dj_d2),intent(in) :: dj
        integer(kind=c_int), intent(in) :: cilmrot_d0
        integer(kind=c_int), intent(in) :: cilmrot_d1
        integer(kind=c_int), intent(in) :: cilmrot_d2
        real(kind=c_double), dimension(cilmrot_d0,cilmrot_d1,cilmrot_d2),intent(out) :: cilmrot
        integer(kind=c_int), intent(in) :: lmax
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call SHRotateRealCoef(cilmrot,cilm,lmax,x,dj,exitstatus=exitstatus)
    end subroutine cbindSHRotateRealCoef

    function cbindSHPowerL(c,c_d0,c_d1,c_d2,l)  bind(c, name="cbind_sh_power_l")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHPowerL
        implicit none
        real(kind=c_double) :: cbindSHPowerL
        integer(kind=c_int), intent(in) :: c_d0
        integer(kind=c_int), intent(in) :: c_d1
        integer(kind=c_int), intent(in) :: c_d2
        real(kind=c_double), dimension(c_d0,c_d1,c_d2),intent(in) :: c
        integer(kind=c_int), intent(in) :: l
        cbindSHPowerL=SHPowerL(c,l)
    end function cbindSHPowerL

    function cbindSHPowerDensityL(c,c_d0,c_d1,c_d2,l)  bind(c, name="cbind_sh_power_density_l")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHPowerDensityL
        implicit none
        real(kind=c_double) :: cbindSHPowerDensityL
        integer(kind=c_int), intent(in) :: c_d0
        integer(kind=c_int), intent(in) :: c_d1
        integer(kind=c_int), intent(in) :: c_d2
        real(kind=c_double), dimension(c_d0,c_d1,c_d2),intent(in) :: c
        integer(kind=c_int), intent(in) :: l
        cbindSHPowerDensityL=SHPowerDensityL(c,l)
    end function cbindSHPowerDensityL

    function cbindSHCrossPowerL(c1,c1_d0,c1_d1,c1_d2,c2,c2_d0,c2_d1,c2_d2,l)  bind(c&
                                  , name="cbind_sh_cross_power_l")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHCrossPowerL
        implicit none
        real(kind=c_double) :: cbindSHCrossPowerL
        integer(kind=c_int), intent(in) :: c1_d0
        integer(kind=c_int), intent(in) :: c1_d1
        integer(kind=c_int), intent(in) :: c1_d2
        real(kind=c_double), dimension(c1_d0,c1_d1,c1_d2),intent(in) :: c1
        integer(kind=c_int), intent(in) :: c2_d0
        integer(kind=c_int), intent(in) :: c2_d1
        integer(kind=c_int), intent(in) :: c2_d2
        real(kind=c_double), dimension(c2_d0,c2_d1,c2_d2),intent(in) :: c2
        integer(kind=c_int), intent(in) :: l
        cbindSHCrossPowerL=SHCrossPowerL(c1,c2,l)
    end function cbindSHCrossPowerL

    function cbindSHCrossPowerDensityL(c1,c1_d0,c1_d1,c1_d2,c2,c2_d0,c2_d1,c2_d2,l)  bind(c&
                                         , name="cbind_sh_cross_power_density_l")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHCrossPowerDensityL
        implicit none
        real(kind=c_double) :: cbindSHCrossPowerDensityL
        integer(kind=c_int), intent(in) :: c1_d0
        integer(kind=c_int), intent(in) :: c1_d1
        integer(kind=c_int), intent(in) :: c1_d2
        real(kind=c_double), dimension(c1_d0,c1_d1,c1_d2),intent(in) :: c1
        integer(kind=c_int), intent(in) :: c2_d0
        integer(kind=c_int), intent(in) :: c2_d1
        integer(kind=c_int), intent(in) :: c2_d2
        real(kind=c_double), dimension(c2_d0,c2_d1,c2_d2),intent(in) :: c2
        integer(kind=c_int), intent(in) :: l
        cbindSHCrossPowerDensityL=SHCrossPowerDensityL(c1,c2,l)
    end function cbindSHCrossPowerDensityL

    subroutine cbindSHPowerSpectrum(c,c_d0,c_d1,c_d2,lmax,spectra,spectra_d0,exitstatus)  bind(c&
                                     , name="cbind_sh_power_spectrum")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHPowerSpectrum
        implicit none
        integer(kind=c_int), intent(in) :: c_d0
        integer(kind=c_int), intent(in) :: c_d1
        integer(kind=c_int), intent(in) :: c_d2
        real(kind=c_double), dimension(c_d0,c_d1,c_d2),intent(in) :: c
        integer(kind=c_int), intent(in) :: lmax
        integer(kind=c_int), intent(in) :: spectra_d0
        real(kind=c_double), dimension(spectra_d0),intent(out) :: spectra
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call SHPowerSpectrum(c,lmax,spectra,exitstatus=exitstatus)
    end subroutine cbindSHPowerSpectrum

    subroutine cbindSHPowerSpectrumDensity(c,c_d0,c_d1,c_d2,lmax,spectra,spectra_d0&
                                            ,exitstatus)  bind(c, name="cbind_sh_power_spectrum_density")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHPowerSpectrumDensity
        implicit none
        integer(kind=c_int), intent(in) :: c_d0
        integer(kind=c_int), intent(in) :: c_d1
        integer(kind=c_int), intent(in) :: c_d2
        real(kind=c_double), dimension(c_d0,c_d1,c_d2),intent(in) :: c
        integer(kind=c_int), intent(in) :: lmax
        integer(kind=c_int), intent(in) :: spectra_d0
        real(kind=c_double), dimension(spectra_d0),intent(out) :: spectra
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call SHPowerSpectrumDensity(c,lmax,spectra,exitstatus=exitstatus)
    end subroutine cbindSHPowerSpectrumDensity

    subroutine cbindSHCrossPowerSpectrum(c1,c1_d0,c1_d1,c1_d2,c2,c2_d0,c2_d1,c2_d2&
                                           ,lmax,cspectra,cspectra_d0,exitstatus)  bind(c&
                                           , name="cbind_sh_cross_power_spectrum")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHCrossPowerSpectrum
        implicit none
        integer(kind=c_int), intent(in) :: c1_d0
        integer(kind=c_int), intent(in) :: c1_d1
        integer(kind=c_int), intent(in) :: c1_d2
        real(kind=c_double), dimension(c1_d0,c1_d1,c1_d2),intent(in) :: c1
        integer(kind=c_int), intent(in) :: c2_d0
        integer(kind=c_int), intent(in) :: c2_d1
        integer(kind=c_int), intent(in) :: c2_d2
        real(kind=c_double), dimension(c2_d0,c2_d1,c2_d2),intent(in) :: c2
        integer(kind=c_int), intent(in) :: lmax
        integer(kind=c_int), intent(in) :: cspectra_d0
        real(kind=c_double), dimension(cspectra_d0),intent(out) :: cspectra
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call SHCrossPowerSpectrum(c1,c2,lmax,cspectra,exitstatus=exitstatus)
    end subroutine cbindSHCrossPowerSpectrum

    subroutine cbindSHCrossPowerSpectrumDensity(c1,c1_d0,c1_d1,c1_d2,c2,c2_d0,c2_d1&
                                                  ,c2_d2,lmax,cspectra,cspectra_d0&
                                                  ,exitstatus)  bind(c, name="cbind_sh_cross_power_spectrum_density")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHCrossPowerSpectrumDensity
        implicit none
        integer(kind=c_int), intent(in) :: c1_d0
        integer(kind=c_int), intent(in) :: c1_d1
        integer(kind=c_int), intent(in) :: c1_d2
        real(kind=c_double), dimension(c1_d0,c1_d1,c1_d2),intent(in) :: c1
        integer(kind=c_int), intent(in) :: c2_d0
        integer(kind=c_int), intent(in) :: c2_d1
        integer(kind=c_int), intent(in) :: c2_d2
        real(kind=c_double), dimension(c2_d0,c2_d1,c2_d2),intent(in) :: c2
        integer(kind=c_int), intent(in) :: lmax
        integer(kind=c_int), intent(in) :: cspectra_d0
        real(kind=c_double), dimension(cspectra_d0),intent(out) :: cspectra
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call SHCrossPowerSpectrumDensity(c1,c2,lmax,cspectra,exitstatus=exitstatus)
    end subroutine cbindSHCrossPowerSpectrumDensity

    subroutine cbindSHAdmitCorr(G,G_d0,G_d1,G_d2,T,T_d0,T_d1,T_d2,lmax,admit,admit_d0&
                                 ,corr,corr_d0,admit_error,admit_error_d0,exitstatus)  bind(c&
                                 , name="cbind_sh_admit_corr")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHAdmitCorr
        implicit none
        integer(kind=c_int), intent(in) :: G_d0
        integer(kind=c_int), intent(in) :: G_d1
        integer(kind=c_int), intent(in) :: G_d2
        real(kind=c_double), dimension(G_d0,G_d1,G_d2),intent(in) :: G
        integer(kind=c_int), intent(in) :: T_d0
        integer(kind=c_int), intent(in) :: T_d1
        integer(kind=c_int), intent(in) :: T_d2
        real(kind=c_double), dimension(T_d0,T_d1,T_d2),intent(in) :: T
        integer(kind=c_int), intent(in) :: lmax
        integer(kind=c_int), intent(in) :: admit_d0
        real(kind=c_double), dimension(admit_d0),intent(out) :: admit
        integer(kind=c_int), intent(in) :: corr_d0
        real(kind=c_double), dimension(corr_d0),intent(out) :: corr
        integer(kind=c_int), intent(in) :: admit_error_d0
        real(kind=c_double), optional,dimension(admit_error_d0),intent(out) :: admit_error
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call SHAdmitCorr(G,T,lmax,admit,corr,admit_error=admit_error,exitstatus=exitstatus)
    end subroutine cbindSHAdmitCorr

    function cbindSHConfidence(l_conf,r)  bind(c, name="cbind_sh_confidence")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHConfidence
        implicit none
        real(kind=c_double) :: cbindSHConfidence
        real(kind=c_double), intent(in) :: r
        integer(kind=c_int), intent(in) :: l_conf
        cbindSHConfidence=SHConfidence(l_conf,r)
    end function cbindSHConfidence

    function cbindSHPowerLC(c,c_d0,c_d1,c_d2,l)  bind(c, name="cbind_sh_power_lc")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHPowerLC
        implicit none
        real(kind=c_double) :: cbindSHPowerLC
        integer(kind=c_int), intent(in) :: c_d0
        integer(kind=c_int), intent(in) :: c_d1
        integer(kind=c_int), intent(in) :: c_d2
        complex(kind=c_double_complex), dimension(c_d0,c_d1,c_d2),intent(in) :: c
        integer(kind=c_int), intent(in) :: l
        cbindSHPowerLC=SHPowerLC(c,l)
    end function cbindSHPowerLC

    function cbindSHPowerDensityLC(c,c_d0,c_d1,c_d2,l)  bind(c, name="cbind_sh_power_density_lc")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHPowerDensityLC
        implicit none
        real(kind=c_double) :: cbindSHPowerDensityLC
        integer(kind=c_int), intent(in) :: c_d0
        integer(kind=c_int), intent(in) :: c_d1
        integer(kind=c_int), intent(in) :: c_d2
        complex(kind=c_double_complex), dimension(c_d0,c_d1,c_d2),intent(in) :: c
        integer(kind=c_int), intent(in) :: l
        cbindSHPowerDensityLC=SHPowerDensityLC(c,l)
    end function cbindSHPowerDensityLC

    function cbindSHCrossPowerLC(c1,c1_d0,c1_d1,c1_d2,c2,c2_d0,c2_d1,c2_d2,l)  bind(c&
                                   , name="cbind_sh_cross_power_lc")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHCrossPowerLC
        implicit none
        complex(kind=c_double_complex) :: cbindSHCrossPowerLC
        integer(kind=c_int), intent(in) :: c1_d0
        integer(kind=c_int), intent(in) :: c1_d1
        integer(kind=c_int), intent(in) :: c1_d2
        complex(kind=c_double_complex), dimension(c1_d0,c1_d1,c1_d2),intent(in) :: c1
        integer(kind=c_int), intent(in) :: c2_d0
        integer(kind=c_int), intent(in) :: c2_d1
        integer(kind=c_int), intent(in) :: c2_d2
        complex(kind=c_double_complex), dimension(c2_d0,c2_d1,c2_d2),intent(in) :: c2
        integer(kind=c_int), intent(in) :: l
        cbindSHCrossPowerLC=SHCrossPowerLC(c1,c2,l)
    end function cbindSHCrossPowerLC

    function cbindSHCrossPowerDensityLC(c1,c1_d0,c1_d1,c1_d2,c2,c2_d0,c2_d1,c2_d2&
                                          ,l)  bind(c, name="cbind_sh_cross_power_density_lc")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHCrossPowerDensityLC
        implicit none
        complex(kind=c_double_complex) :: cbindSHCrossPowerDensityLC
        integer(kind=c_int), intent(in) :: c1_d0
        integer(kind=c_int), intent(in) :: c1_d1
        integer(kind=c_int), intent(in) :: c1_d2
        complex(kind=c_double_complex), dimension(c1_d0,c1_d1,c1_d2),intent(in) :: c1
        integer(kind=c_int), intent(in) :: c2_d0
        integer(kind=c_int), intent(in) :: c2_d1
        integer(kind=c_int), intent(in) :: c2_d2
        complex(kind=c_double_complex), dimension(c2_d0,c2_d1,c2_d2),intent(in) :: c2
        integer(kind=c_int), intent(in) :: l
        cbindSHCrossPowerDensityLC=SHCrossPowerDensityLC(c1,c2,l)
    end function cbindSHCrossPowerDensityLC

    subroutine cbindSHPowerSpectrumC(c,c_d0,c_d1,c_d2,lmax,spectra,spectra_d0,exitstatus)  bind(c&
                                      , name="cbind_sh_power_spectrum_c")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHPowerSpectrumC
        implicit none
        integer(kind=c_int), intent(in) :: c_d0
        integer(kind=c_int), intent(in) :: c_d1
        integer(kind=c_int), intent(in) :: c_d2
        complex(kind=c_double_complex), dimension(c_d0,c_d1,c_d2),intent(in) :: c
        integer(kind=c_int), intent(in) :: lmax
        integer(kind=c_int), intent(in) :: spectra_d0
        real(kind=c_double), dimension(spectra_d0),intent(out) :: spectra
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call SHPowerSpectrumC(c,lmax,spectra,exitstatus=exitstatus)
    end subroutine cbindSHPowerSpectrumC

    subroutine cbindSHPowerSpectrumDensityC(c,c_d0,c_d1,c_d2,lmax,spectra,spectra_d0&
                                             ,exitstatus)  bind(c, name="cbind_sh_power_spectrum_density_c")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHPowerSpectrumDensityC
        implicit none
        integer(kind=c_int), intent(in) :: c_d0
        integer(kind=c_int), intent(in) :: c_d1
        integer(kind=c_int), intent(in) :: c_d2
        complex(kind=c_double_complex), dimension(c_d0,c_d1,c_d2),intent(in) :: c
        integer(kind=c_int), intent(in) :: lmax
        integer(kind=c_int), intent(in) :: spectra_d0
        real(kind=c_double), dimension(spectra_d0),intent(out) :: spectra
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call SHPowerSpectrumDensityC(c,lmax,spectra,exitstatus=exitstatus)
    end subroutine cbindSHPowerSpectrumDensityC

    subroutine cbindSHCrossPowerSpectrumC(c1,c1_d0,c1_d1,c1_d2,c2,c2_d0,c2_d1,c2_d2&
                                            ,lmax,cspectra,cspectra_d0,exitstatus)  bind(c&
                                            , name="cbind_sh_cross_power_spectrum_c")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHCrossPowerSpectrumC
        implicit none
        integer(kind=c_int), intent(in) :: c1_d0
        integer(kind=c_int), intent(in) :: c1_d1
        integer(kind=c_int), intent(in) :: c1_d2
        complex(kind=c_double_complex), dimension(c1_d0,c1_d1,c1_d2),intent(in) :: c1
        integer(kind=c_int), intent(in) :: c2_d0
        integer(kind=c_int), intent(in) :: c2_d1
        integer(kind=c_int), intent(in) :: c2_d2
        complex(kind=c_double_complex), dimension(c2_d0,c2_d1,c2_d2),intent(in) :: c2
        integer(kind=c_int), intent(in) :: lmax
        integer(kind=c_int), intent(in) :: cspectra_d0
        complex(kind=c_double_complex), dimension(cspectra_d0),intent(out) :: cspectra
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call SHCrossPowerSpectrumC(c1,c2,lmax,cspectra,exitstatus=exitstatus)
    end subroutine cbindSHCrossPowerSpectrumC

    subroutine cbindSHCrossPowerSpectrumDensityC(c1,c1_d0,c1_d1,c1_d2,c2,c2_d0,c2_d1&
                                                   ,c2_d2,lmax,cspectra,cspectra_d0&
                                                   ,exitstatus)  bind(c, name="cbind_sh_cross_power_spectrum_density_c")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHCrossPowerSpectrumDensityC
        implicit none
        integer(kind=c_int), intent(in) :: c1_d0
        integer(kind=c_int), intent(in) :: c1_d1
        integer(kind=c_int), intent(in) :: c1_d2
        complex(kind=c_double_complex), dimension(c1_d0,c1_d1,c1_d2),intent(in) :: c1
        integer(kind=c_int), intent(in) :: c2_d0
        integer(kind=c_int), intent(in) :: c2_d1
        integer(kind=c_int), intent(in) :: c2_d2
        complex(kind=c_double_complex), dimension(c2_d0,c2_d1,c2_d2),intent(in) :: c2
        integer(kind=c_int), intent(in) :: lmax
        integer(kind=c_int), intent(in) :: cspectra_d0
        complex(kind=c_double_complex), dimension(cspectra_d0),intent(out) :: cspectra
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call SHCrossPowerSpectrumDensityC(c1,c2,lmax,cspectra,exitstatus=exitstatus)
    end subroutine cbindSHCrossPowerSpectrumDensityC

    subroutine cbindSHMultiTaperSE(mtse,mtse_d0,sd,sd_d0,sh,sh_d0,sh_d1,sh_d2,lmax&
                                       ,tapers,tapers_d0,tapers_d1,taper_order,taper_order_d0&
                                       ,lmaxt,k,alpha,alpha_d0,lat,lon,taper_wt,taper_wt_d0&
                                       ,norm,csphase,exitstatus)  bind(c, name="cbind_sh_multi_taper_se")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHMultiTaperSE
        implicit none
        integer(kind=c_int), intent(in) :: mtse_d0
        real(kind=c_double), dimension(mtse_d0),intent(out) :: mtse
        integer(kind=c_int), intent(in) :: sd_d0
        real(kind=c_double), dimension(sd_d0),intent(out) :: sd
        integer(kind=c_int), intent(in) :: sh_d0
        integer(kind=c_int), intent(in) :: sh_d1
        integer(kind=c_int), intent(in) :: sh_d2
        real(kind=c_double), dimension(sh_d0,sh_d1,sh_d2),intent(in) :: sh
        integer(kind=c_int), intent(in) :: tapers_d0
        integer(kind=c_int), intent(in) :: tapers_d1
        real(kind=c_double), dimension(tapers_d0,tapers_d1),intent(in) :: tapers
        integer(kind=c_int), intent(in) :: lmax
        integer(kind=c_int), intent(in) :: lmaxt
        integer(kind=c_int), intent(in) :: k
        integer(kind=c_int), intent(in) :: taper_order_d0
        integer(kind=c_int), dimension(taper_order_d0),intent(in) :: taper_order
        integer(kind=c_int), intent(in) :: alpha_d0
        real(kind=c_double), optional,dimension(alpha_d0),intent(in) :: alpha
        real(kind=c_double), optional,intent(in) :: lat
        real(kind=c_double), optional,intent(in) :: lon
        integer(kind=c_int), intent(in) :: taper_wt_d0
        real(kind=c_double), optional,dimension(taper_wt_d0),intent(in) :: taper_wt
        integer(kind=c_int), optional,intent(in) :: csphase
        integer(kind=c_int), optional,intent(in) :: norm
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call SHMultiTaperSE(mtse,sd,sh,lmax,tapers,taper_order,lmaxt,k,alpha=alpha&
                                ,lat=lat,lon=lon,taper_wt=taper_wt,norm=norm,csphase=csphase&
                                ,exitstatus=exitstatus)
    end subroutine cbindSHMultiTaperSE

    subroutine cbindSHMultiTaperCSE(mtse,mtse_d0,sd,sd_d0,sh1,sh1_d0,sh1_d1,sh1_d2&
                                        ,lmax1,sh2,sh2_d0,sh2_d1,sh2_d2,lmax2,tapers&
                                        ,tapers_d0,tapers_d1,taper_order,taper_order_d0&
                                        ,lmaxt,k,alpha,alpha_d0,lat,lon,taper_wt,taper_wt_d0&
                                        ,norm,csphase,exitstatus)  bind(c, name="cbind_sh_multi_taper_cse")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHMultiTaperCSE
        implicit none
        integer(kind=c_int), intent(in) :: mtse_d0
        real(kind=c_double), dimension(mtse_d0),intent(out) :: mtse
        integer(kind=c_int), intent(in) :: sd_d0
        real(kind=c_double), dimension(sd_d0),intent(out) :: sd
        integer(kind=c_int), intent(in) :: sh1_d0
        integer(kind=c_int), intent(in) :: sh1_d1
        integer(kind=c_int), intent(in) :: sh1_d2
        real(kind=c_double), dimension(sh1_d0,sh1_d1,sh1_d2),intent(in) :: sh1
        integer(kind=c_int), intent(in) :: sh2_d0
        integer(kind=c_int), intent(in) :: sh2_d1
        integer(kind=c_int), intent(in) :: sh2_d2
        real(kind=c_double), dimension(sh2_d0,sh2_d1,sh2_d2),intent(in) :: sh2
        integer(kind=c_int), intent(in) :: tapers_d0
        integer(kind=c_int), intent(in) :: tapers_d1
        real(kind=c_double), dimension(tapers_d0,tapers_d1),intent(in) :: tapers
        integer(kind=c_int), intent(in) :: lmax1
        integer(kind=c_int), intent(in) :: lmax2
        integer(kind=c_int), intent(in) :: lmaxt
        integer(kind=c_int), intent(in) :: k
        integer(kind=c_int), intent(in) :: taper_order_d0
        integer(kind=c_int), dimension(taper_order_d0),intent(in) :: taper_order
        integer(kind=c_int), intent(in) :: alpha_d0
        real(kind=c_double), optional,dimension(alpha_d0),intent(in) :: alpha
        real(kind=c_double), optional,intent(in) :: lat
        real(kind=c_double), optional,intent(in) :: lon
        integer(kind=c_int), intent(in) :: taper_wt_d0
        real(kind=c_double), optional,dimension(taper_wt_d0),intent(in) :: taper_wt
        integer(kind=c_int), optional,intent(in) :: csphase
        integer(kind=c_int), optional,intent(in) :: norm
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call SHMultiTaperCSE(mtse,sd,sh1,lmax1,sh2,lmax2,tapers,taper_order,lmaxt&
                                 ,k,alpha=alpha,lat=lat,lon=lon,taper_wt=taper_wt&
                                 ,norm=norm,csphase=csphase,exitstatus=exitstatus)
    end subroutine cbindSHMultiTaperCSE

    subroutine cbindSHLocalizedAdmitCorr(tapers,tapers_d0,tapers_d1,taper_order,taper_order_d0&
                                               ,lwin,lat,lon,g,g_d0,g_d1,g_d2,t,t_d0&
                                               ,t_d1,t_d2,lmax,admit,admit_d0,corr&
                                               ,corr_d0,k,admit_error,admit_error_d0&
                                               ,corr_error,corr_error_d0,taper_wt&
                                               ,taper_wt_d0,mtdef,k1linsig,exitstatus)  bind(c&
                                               , name="cbind_sh_localized_admit_corr")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHLocalizedAdmitCorr
        implicit none
        integer(kind=c_int), intent(in) :: tapers_d0
        integer(kind=c_int), intent(in) :: tapers_d1
        real(kind=c_double), dimension(tapers_d0,tapers_d1),intent(in) :: tapers
        real(kind=c_double), intent(in) :: lat
        real(kind=c_double), intent(in) :: lon
        integer(kind=c_int), intent(in) :: g_d0
        integer(kind=c_int), intent(in) :: g_d1
        integer(kind=c_int), intent(in) :: g_d2
        real(kind=c_double), dimension(g_d0,g_d1,g_d2),intent(in) :: g
        integer(kind=c_int), intent(in) :: t_d0
        integer(kind=c_int), intent(in) :: t_d1
        integer(kind=c_int), intent(in) :: t_d2
        real(kind=c_double), dimension(t_d0,t_d1,t_d2),intent(in) :: t
        integer(kind=c_int), intent(in) :: lwin
        integer(kind=c_int), intent(in) :: lmax
        integer(kind=c_int), intent(in) :: k
        integer(kind=c_int), intent(in) :: taper_order_d0
        integer(kind=c_int), dimension(taper_order_d0),intent(in) :: taper_order
        integer(kind=c_int), intent(in) :: admit_d0
        real(kind=c_double), dimension(admit_d0),intent(out) :: admit
        integer(kind=c_int), intent(in) :: corr_d0
        real(kind=c_double), dimension(corr_d0),intent(out) :: corr
        integer(kind=c_int), intent(in) :: admit_error_d0
        real(kind=c_double), optional,dimension(admit_error_d0),intent(out) :: admit_error
        integer(kind=c_int), intent(in) :: corr_error_d0
        real(kind=c_double), optional,dimension(corr_error_d0),intent(out) :: corr_error
        integer(kind=c_int), optional,intent(in) :: mtdef
        integer(kind=c_int), optional,intent(in) :: k1linsig
        integer(kind=c_int), intent(in) :: taper_wt_d0
        real(kind=c_double), optional,dimension(taper_wt_d0),intent(in) :: taper_wt
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call SHLocalizedAdmitCorr(tapers,taper_order,lwin,lat,lon,g,t,lmax,admit,corr&
                                        ,k,admit_error=admit_error,corr_error=corr_error&
                                        ,taper_wt=taper_wt,mtdef=mtdef,k1linsig=k1linsig&
                                        ,exitstatus=exitstatus)
    end subroutine cbindSHLocalizedAdmitCorr

    subroutine cbindSHReturnTapers(theta0,lmax,tapers,tapers_d0,tapers_d1,eigenvalues&
                                         ,eigenvalues_d0,taper_order,taper_order_d0&
                                         ,degrees,degrees_d0,exitstatus)  bind(c, name="cbind_sh_return_tapers")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHReturnTapers
        implicit none
        real(kind=c_double), intent(in) :: theta0
        integer(kind=c_int), intent(in) :: lmax
        integer(kind=c_int), intent(in) :: tapers_d0
        integer(kind=c_int), intent(in) :: tapers_d1
        real(kind=c_double), dimension(tapers_d0,tapers_d1),intent(out) :: tapers
        integer(kind=c_int), intent(in) :: eigenvalues_d0
        real(kind=c_double), dimension(eigenvalues_d0),intent(out) :: eigenvalues
        integer(kind=c_int), intent(in) :: taper_order_d0
        integer(kind=c_int), dimension(taper_order_d0),intent(out) :: taper_order
        integer(kind=c_int), intent(in) :: degrees_d0
        integer(kind=c_int), optional,dimension(degrees_d0),intent(in) :: degrees
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call SHReturnTapers(theta0,lmax,tapers,eigenvalues,taper_order,degrees=degrees&
                                  ,exitstatus=exitstatus)
    end subroutine cbindSHReturnTapers

    subroutine cbindSHReturnTapersM(theta0,lmax,m,tapers,tapers_d0,tapers_d1,eigenvalues&
                                          ,eigenvalues_d0,shannon,degrees,degrees_d0&
                                          ,ntapers,exitstatus)  bind(c, name="cbind_sh_return_tapers_m")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHReturnTapersM
        implicit none
        real(kind=c_double), intent(in) :: theta0
        integer(kind=c_int), intent(in) :: lmax
        integer(kind=c_int), intent(in) :: m
        integer(kind=c_int), intent(in) :: tapers_d0
        integer(kind=c_int), intent(in) :: tapers_d1
        real(kind=c_double), dimension(tapers_d0,tapers_d1),intent(out) :: tapers
        integer(kind=c_int), intent(in) :: eigenvalues_d0
        real(kind=c_double), dimension(eigenvalues_d0),intent(out) :: eigenvalues
        real(kind=c_double), optional,intent(out) :: shannon
        integer(kind=c_int), intent(in) :: degrees_d0
        integer(kind=c_int), optional,dimension(degrees_d0),intent(in) :: degrees
        integer(kind=c_int), optional,intent(out) :: ntapers
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call SHReturnTapersM(theta0,lmax,m,tapers,eigenvalues,shannon=shannon,degrees=degrees&
                                   ,ntapers=ntapers,exitstatus=exitstatus)
    end subroutine cbindSHReturnTapersM

    subroutine cbindComputeDm(dllm,dllm_d0,dllm_d1,lmax,m,theta0,degrees,degrees_d0&
                                  ,exitstatus)  bind(c, name="cbind_compute_dm")
        use, intrinsic :: iso_c_binding
        use shtools, only: ComputeDm
        implicit none
        integer(kind=c_int), intent(in) :: dllm_d0
        integer(kind=c_int), intent(in) :: dllm_d1
        real(kind=c_double), dimension(dllm_d0,dllm_d1),intent(out) :: dllm
        real(kind=c_double), intent(in) :: theta0
        integer(kind=c_int), intent(in) :: lmax
        integer(kind=c_int), intent(in) :: m
        integer(kind=c_int), intent(in) :: degrees_d0
        integer(kind=c_int), optional,dimension(degrees_d0),intent(in) :: degrees
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call ComputeDm(dllm,lmax,m,theta0,degrees=degrees,exitstatus=exitstatus)
    end subroutine cbindComputeDm

    subroutine cbindComputeDG82(dG82,dG82_d0,dG82_d1,lmax,m,theta0,exitstatus)  bind(c&
                                    , name="cbind_compute_dg82")
        use, intrinsic :: iso_c_binding
        use shtools, only: ComputeDG82
        implicit none
        integer(kind=c_int), intent(in) :: dG82_d0
        integer(kind=c_int), intent(in) :: dG82_d1
        real(kind=c_double), dimension(dG82_d0,dG82_d1),intent(out) :: dG82
        real(kind=c_double), intent(in) :: theta0
        integer(kind=c_int), intent(in) :: lmax
        integer(kind=c_int), intent(in) :: m
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call ComputeDG82(dG82,lmax,m,theta0,exitstatus=exitstatus)
    end subroutine cbindComputeDG82

    function cbindSHFindLWin(theta0,m,alpha,taper_number)  bind(c, name="cbind_sh_find_l_win")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHFindLWin
        implicit none
        integer(kind=c_int) :: cbindSHFindLWin
        real(kind=c_double), intent(in) :: theta0
        real(kind=c_double), intent(in) :: alpha
        integer(kind=c_int), intent(in) :: m
        integer(kind=c_int), optional,intent(in) :: taper_number
        cbindSHFindLWin=SHFindLWin(theta0,m,alpha,taper_number=taper_number)
    end function cbindSHFindLWin

    subroutine cbindSHBiasK(tapers,tapers_d0,tapers_d1,lwin,k,incspectra,incspectra_d0&
                                  ,ldata,outcspectra,outcspectra_d0,taper_wt,taper_wt_d0&
                                  ,save_cg,exitstatus)  bind(c, name="cbind_sh_bias_k")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHBiasK
        implicit none
        integer(kind=c_int), intent(in) :: tapers_d0
        integer(kind=c_int), intent(in) :: tapers_d1
        real(kind=c_double), dimension(tapers_d0,tapers_d1),intent(in) :: tapers
        integer(kind=c_int), intent(in) :: incspectra_d0
        real(kind=c_double), dimension(incspectra_d0),intent(in) :: incspectra
        integer(kind=c_int), intent(in) :: outcspectra_d0
        real(kind=c_double), dimension(outcspectra_d0),intent(out) :: outcspectra
        integer(kind=c_int), intent(in) :: lwin
        integer(kind=c_int), intent(in) :: ldata
        integer(kind=c_int), intent(in) :: k
        integer(kind=c_int), intent(in) :: taper_wt_d0
        real(kind=c_double), optional,dimension(taper_wt_d0),intent(in) :: taper_wt
        integer(kind=c_int), optional,intent(in) :: save_cg
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call SHBiasK(tapers,lwin,k,incspectra,ldata,outcspectra,taper_wt=taper_wt&
                           ,save_cg=save_cg,exitstatus=exitstatus)
    end subroutine cbindSHBiasK

    subroutine cbindSHMTCouplingMatrix(Mmt,Mmt_d0,Mmt_d1,lmax,tapers_power,tapers_power_d0&
                                          ,tapers_power_d1,lwin,k,taper_wt,taper_wt_d0&
                                          ,exitstatus)  bind(c, name="cbind_shmt_coupling_matrix")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHMTCouplingMatrix
        implicit none
        integer(kind=c_int), intent(in) :: Mmt_d0
        integer(kind=c_int), intent(in) :: Mmt_d1
        real(kind=c_double), dimension(Mmt_d0,Mmt_d1),intent(out) :: Mmt
        integer(kind=c_int), intent(in) :: tapers_power_d0
        integer(kind=c_int), intent(in) :: tapers_power_d1
        real(kind=c_double), dimension(tapers_power_d0,tapers_power_d1),intent(in) :: tapers_power
        integer(kind=c_int), intent(in) :: taper_wt_d0
        real(kind=c_double), optional,dimension(taper_wt_d0),intent(in) :: taper_wt
        integer(kind=c_int), intent(in) :: lmax
        integer(kind=c_int), intent(in) :: k
        integer(kind=c_int), intent(in) :: lwin
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call SHMTCouplingMatrix(Mmt,lmax,tapers_power,lwin,k,taper_wt=taper_wt,exitstatus=exitstatus)
    end subroutine cbindSHMTCouplingMatrix

    subroutine cbindSHBiasAdmitCorr(sgt,sgt_d0,sgg,sgg_d0,stt,stt_d0,lmax,tapers,tapers_d0&
                                       ,tapers_d1,lwin,k,admit,admit_d0,corr,corr_d0&
                                       ,mtdef,taper_wt,taper_wt_d0,exitstatus)  bind(c&
                                       , name="cbind_sh_bias_admit_corr")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHBiasAdmitCorr
        implicit none
        integer(kind=c_int), intent(in) :: sgt_d0
        real(kind=c_double), dimension(sgt_d0),intent(in) :: sgt
        integer(kind=c_int), intent(in) :: sgg_d0
        real(kind=c_double), dimension(sgg_d0),intent(in) :: sgg
        integer(kind=c_int), intent(in) :: stt_d0
        real(kind=c_double), dimension(stt_d0),intent(in) :: stt
        integer(kind=c_int), intent(in) :: tapers_d0
        integer(kind=c_int), intent(in) :: tapers_d1
        real(kind=c_double), dimension(tapers_d0,tapers_d1),intent(in) :: tapers
        integer(kind=c_int), intent(in) :: lmax
        integer(kind=c_int), intent(in) :: lwin
        integer(kind=c_int), intent(in) :: k
        integer(kind=c_int), intent(in) :: admit_d0
        real(kind=c_double), dimension(admit_d0),intent(out) :: admit
        integer(kind=c_int), intent(in) :: corr_d0
        real(kind=c_double), dimension(corr_d0),intent(out) :: corr
        integer(kind=c_int), optional,intent(in) :: mtdef
        integer(kind=c_int), intent(in) :: taper_wt_d0
        real(kind=c_double), optional,dimension(taper_wt_d0),intent(in) :: taper_wt
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call SHBiasAdmitCorr(sgt,sgg,stt,lmax,tapers,lwin,k,admit,corr,mtdef=mtdef&
                                ,taper_wt=taper_wt,exitstatus=exitstatus)
    end subroutine cbindSHBiasAdmitCorr

    subroutine cbindSHMTDebias(mtdebias,mtdebias_d0,mtdebias_d1,mtspectra,mtspectra_d0&
                                       ,mtspectra_d1,lmax,tapers,tapers_d0,tapers_d1&
                                       ,lwin,k,nl,lmid,lmid_d0,n,taper_wt,taper_wt_d0&
                                       ,exitstatus)  bind(c, name="cbind_shmt_debias")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHMTDebias
        implicit none
        integer(kind=c_int), intent(in) :: mtdebias_d0
        integer(kind=c_int), intent(in) :: mtdebias_d1
        real(kind=c_double), dimension(mtdebias_d0,mtdebias_d1),intent(out) :: mtdebias
        integer(kind=c_int), intent(in) :: lmid_d0
        real(kind=c_double), dimension(lmid_d0),intent(out) :: lmid
        integer(kind=c_int), intent(in) :: mtspectra_d0
        integer(kind=c_int), intent(in) :: mtspectra_d1
        real(kind=c_double), dimension(mtspectra_d0,mtspectra_d1),intent(in) :: mtspectra
        integer(kind=c_int), intent(in) :: tapers_d0
        integer(kind=c_int), intent(in) :: tapers_d1
        real(kind=c_double), dimension(tapers_d0,tapers_d1),intent(in) :: tapers
        integer(kind=c_int), intent(in) :: taper_wt_d0
        real(kind=c_double), optional,dimension(taper_wt_d0),intent(in) :: taper_wt
        integer(kind=c_int), intent(in) :: lmax
        integer(kind=c_int), intent(in) :: k
        integer(kind=c_int), intent(in) :: lwin
        integer(kind=c_int), intent(in) :: nl
        integer(kind=c_int), intent(out) :: n
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call SHMTDebias(mtdebias,mtspectra,lmax,tapers,lwin,k,nl,lmid,n,taper_wt=taper_wt&
                                ,exitstatus=exitstatus)
    end subroutine cbindSHMTDebias

    subroutine cbindSHMTVarOpt(l,tapers,tapers_d0,tapers_d1,taper_order,taper_order_d0&
                                ,lwin,kmax,Sff,Sff_d0,var_opt,var_opt_d0,var_unit&
                                ,var_unit_d0,weight_opt,weight_opt_d0,weight_opt_d1&
                                ,unweighted_covar,unweighted_covar_d0,unweighted_covar_d1&
                                ,nocross,exitstatus)  bind(c, name="cbind_shmt_var_opt")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHMTVarOpt
        implicit none
        integer(kind=c_int), intent(in) :: tapers_d0
        integer(kind=c_int), intent(in) :: tapers_d1
        real(kind=c_double), dimension(tapers_d0,tapers_d1),intent(in) :: tapers
        integer(kind=c_int), intent(in) :: Sff_d0
        real(kind=c_double), dimension(Sff_d0),intent(in) :: Sff
        integer(kind=c_int), intent(in) :: var_opt_d0
        real(kind=c_double), dimension(var_opt_d0),intent(out) :: var_opt
        integer(kind=c_int), intent(in) :: var_unit_d0
        real(kind=c_double), dimension(var_unit_d0),intent(out) :: var_unit
        integer(kind=c_int), intent(in) :: l
        integer(kind=c_int), intent(in) :: lwin
        integer(kind=c_int), intent(in) :: kmax
        integer(kind=c_int), intent(in) :: taper_order_d0
        integer(kind=c_int), dimension(taper_order_d0),intent(in) :: taper_order
        integer(kind=c_int), intent(in) :: weight_opt_d0
        integer(kind=c_int), intent(in) :: weight_opt_d1
        real(kind=c_double), optional,dimension(weight_opt_d0,weight_opt_d1),intent(out) :: weight_opt
        integer(kind=c_int), intent(in) :: unweighted_covar_d0
        integer(kind=c_int), intent(in) :: unweighted_covar_d1
        real(kind=c_double), optional,dimension(unweighted_covar_d0,unweighted_covar_d1)&
                           ,intent(out) :: unweighted_covar
        integer(kind=c_int), optional,intent(in) :: nocross
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call SHMTVarOpt(l,tapers,taper_order,lwin,kmax,Sff,var_opt,var_unit,weight_opt=weight_opt&
                         ,unweighted_covar=unweighted_covar,nocross=nocross,exitstatus=exitstatus)
    end subroutine cbindSHMTVarOpt

    subroutine cbindSHMTVar(l,tapers,tapers_d0,tapers_d1,taper_order,taper_order_d0&
                             ,lwin,kmax,Sff,Sff_d0,variance,taper_wt,taper_wt_d0,unweighted_covar&
                             ,unweighted_covar_d0,unweighted_covar_d1,nocross,exitstatus)  bind(c&
                             , name="cbind_shmt_var")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHMTVar
        implicit none
        integer(kind=c_int), intent(in) :: tapers_d0
        integer(kind=c_int), intent(in) :: tapers_d1
        real(kind=c_double), dimension(tapers_d0,tapers_d1),intent(in) :: tapers
        integer(kind=c_int), intent(in) :: Sff_d0
        real(kind=c_double), dimension(Sff_d0),intent(in) :: Sff
        real(kind=c_double), intent(out) :: variance
        integer(kind=c_int), intent(in) :: l
        integer(kind=c_int), intent(in) :: lwin
        integer(kind=c_int), intent(in) :: kmax
        integer(kind=c_int), intent(in) :: taper_order_d0
        integer(kind=c_int), dimension(taper_order_d0),intent(in) :: taper_order
        integer(kind=c_int), intent(in) :: taper_wt_d0
        real(kind=c_double), optional,dimension(taper_wt_d0),intent(in) :: taper_wt
        integer(kind=c_int), intent(in) :: unweighted_covar_d0
        integer(kind=c_int), intent(in) :: unweighted_covar_d1
        real(kind=c_double), optional,dimension(unweighted_covar_d0,unweighted_covar_d1)&
                           ,intent(out) :: unweighted_covar
        integer(kind=c_int), optional,intent(in) :: nocross
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call SHMTVar(l,tapers,taper_order,lwin,kmax,Sff,variance,taper_wt=taper_wt&
                      ,unweighted_covar=unweighted_covar,nocross=nocross,exitstatus=exitstatus)
    end subroutine cbindSHMTVar

    function cbindSHSjkPG(incspectra,incspectra_d0,l,m,mprime,hj_real,hj_real_d0,hk_real&
                                    ,hk_real_d0,mj,mk,lwin,hkcc)  bind(c, name="cbind_sh_sjk_pg")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHSjkPG
        implicit none
        complex(kind=c_double_complex) :: cbindSHSjkPG
        integer(kind=c_int), intent(in) :: incspectra_d0
        real(kind=c_double), dimension(incspectra_d0),intent(in) :: incspectra
        integer(kind=c_int), intent(in) :: hj_real_d0
        real(kind=c_double), dimension(hj_real_d0),intent(in) :: hj_real
        integer(kind=c_int), intent(in) :: hk_real_d0
        real(kind=c_double), dimension(hk_real_d0),intent(in) :: hk_real
        integer(kind=c_int), intent(in) :: lwin
        integer(kind=c_int), intent(in) :: l
        integer(kind=c_int), intent(in) :: m
        integer(kind=c_int), intent(in) :: mprime
        integer(kind=c_int), intent(in) :: mj
        integer(kind=c_int), intent(in) :: mk
        integer(kind=c_int), intent(in) :: hkcc
        cbindSHSjkPG=SHSjkPG(incspectra,l,m,mprime,hj_real,hk_real,mj,mk,lwin,hkcc)
    end function cbindSHSjkPG

    subroutine cbindSHReturnTapersMap(tapers,tapers_d0,tapers_d1,eigenvalues,eigenvalues_d0&
                                            ,dh_mask,dh_mask_d0,dh_mask_d1,n_dh,lmax&
                                            ,sampling,ntapers,degrees,degrees_d0,exitstatus)  bind(c&
                                            , name="cbind_sh_return_tapers_map")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHReturnTapersMap
        implicit none
        integer(kind=c_int), intent(in) :: tapers_d0
        integer(kind=c_int), intent(in) :: tapers_d1
        real(kind=c_double), dimension(tapers_d0,tapers_d1),intent(out) :: tapers
        integer(kind=c_int), intent(in) :: eigenvalues_d0
        real(kind=c_double), dimension(eigenvalues_d0),intent(out) :: eigenvalues
        integer(kind=c_int), intent(in) :: dh_mask_d0
        integer(kind=c_int), intent(in) :: dh_mask_d1
        integer(kind=c_int), dimension(dh_mask_d0,dh_mask_d1),intent(in) :: dh_mask
        integer(kind=c_int), intent(in) :: n_dh
        integer(kind=c_int), intent(in) :: lmax
        integer(kind=c_int), intent(in) :: sampling
        integer(kind=c_int), optional,intent(in) :: ntapers
        integer(kind=c_int), intent(in) :: degrees_d0
        integer(kind=c_int), optional,dimension(degrees_d0),intent(in) :: degrees
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call SHReturnTapersMap(tapers,eigenvalues,dh_mask,n_dh,lmax,sampling,ntapers=ntapers&
                                     ,degrees=degrees,exitstatus=exitstatus)
    end subroutine cbindSHReturnTapersMap

    subroutine cbindSHBiasKMask(tapers,tapers_d0,tapers_d1,lwin,k,incspectra,incspectra_d0&
                                      ,ldata,outcspectra,outcspectra_d0,taper_wt,taper_wt_d0&
                                      ,save_cg,exitstatus)  bind(c, name="cbind_sh_bias_k_mask")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHBiasKMask
        implicit none
        integer(kind=c_int), intent(in) :: tapers_d0
        integer(kind=c_int), intent(in) :: tapers_d1
        real(kind=c_double), dimension(tapers_d0,tapers_d1),intent(in) :: tapers
        integer(kind=c_int), intent(in) :: incspectra_d0
        real(kind=c_double), dimension(incspectra_d0),intent(in) :: incspectra
        integer(kind=c_int), intent(in) :: outcspectra_d0
        real(kind=c_double), dimension(outcspectra_d0),intent(out) :: outcspectra
        integer(kind=c_int), intent(in) :: lwin
        integer(kind=c_int), intent(in) :: ldata
        integer(kind=c_int), intent(in) :: k
        integer(kind=c_int), intent(in) :: taper_wt_d0
        real(kind=c_double), optional,dimension(taper_wt_d0),intent(in) :: taper_wt
        integer(kind=c_int), optional,intent(in) :: save_cg
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call SHBiasKMask(tapers,lwin,k,incspectra,ldata,outcspectra,taper_wt=taper_wt&
                               ,save_cg=save_cg,exitstatus=exitstatus)
    end subroutine cbindSHBiasKMask

    subroutine cbindSHMultiTaperMaskSE(mtse,mtse_d0,sd,sd_d0,sh,sh_d0,sh_d1,sh_d2&
                                           ,lmax,tapers,tapers_d0,tapers_d1,lmaxt&
                                           ,k,taper_wt,taper_wt_d0,norm,csphase,exitstatus)  bind(c&
                                           , name="cbind_sh_multi_taper_mask_se")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHMultiTaperMaskSE
        implicit none
        integer(kind=c_int), intent(in) :: mtse_d0
        real(kind=c_double), dimension(mtse_d0),intent(out) :: mtse
        integer(kind=c_int), intent(in) :: sd_d0
        real(kind=c_double), dimension(sd_d0),intent(out) :: sd
        integer(kind=c_int), intent(in) :: sh_d0
        integer(kind=c_int), intent(in) :: sh_d1
        integer(kind=c_int), intent(in) :: sh_d2
        real(kind=c_double), dimension(sh_d0,sh_d1,sh_d2),intent(in) :: sh
        integer(kind=c_int), intent(in) :: tapers_d0
        integer(kind=c_int), intent(in) :: tapers_d1
        real(kind=c_double), dimension(tapers_d0,tapers_d1),intent(in) :: tapers
        integer(kind=c_int), intent(in) :: lmax
        integer(kind=c_int), intent(in) :: lmaxt
        integer(kind=c_int), intent(in) :: k
        integer(kind=c_int), intent(in) :: taper_wt_d0
        real(kind=c_double), optional,dimension(taper_wt_d0),intent(in) :: taper_wt
        integer(kind=c_int), optional,intent(in) :: csphase
        integer(kind=c_int), optional,intent(in) :: norm
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call SHMultiTaperMaskSE(mtse,sd,sh,lmax,tapers,lmaxt,k,taper_wt=taper_wt,norm=norm&
                                    ,csphase=csphase,exitstatus=exitstatus)
    end subroutine cbindSHMultiTaperMaskSE

    subroutine cbindSHMultiTaperMaskCSE(mtse,mtse_d0,sd,sd_d0,sh1,sh1_d0,sh1_d1,sh1_d2&
                                            ,lmax1,sh2,sh2_d0,sh2_d1,sh2_d2,lmax2&
                                            ,tapers,tapers_d0,tapers_d1,lmaxt,k,taper_wt&
                                            ,taper_wt_d0,norm,csphase,exitstatus)  bind(c&
                                            , name="cbind_sh_multi_taper_mask_cse")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHMultiTaperMaskCSE
        implicit none
        integer(kind=c_int), intent(in) :: mtse_d0
        real(kind=c_double), dimension(mtse_d0),intent(out) :: mtse
        integer(kind=c_int), intent(in) :: sd_d0
        real(kind=c_double), dimension(sd_d0),intent(out) :: sd
        integer(kind=c_int), intent(in) :: sh1_d0
        integer(kind=c_int), intent(in) :: sh1_d1
        integer(kind=c_int), intent(in) :: sh1_d2
        real(kind=c_double), dimension(sh1_d0,sh1_d1,sh1_d2),intent(in) :: sh1
        integer(kind=c_int), intent(in) :: sh2_d0
        integer(kind=c_int), intent(in) :: sh2_d1
        integer(kind=c_int), intent(in) :: sh2_d2
        real(kind=c_double), dimension(sh2_d0,sh2_d1,sh2_d2),intent(in) :: sh2
        integer(kind=c_int), intent(in) :: tapers_d0
        integer(kind=c_int), intent(in) :: tapers_d1
        real(kind=c_double), dimension(tapers_d0,tapers_d1),intent(in) :: tapers
        integer(kind=c_int), intent(in) :: lmax1
        integer(kind=c_int), intent(in) :: lmax2
        integer(kind=c_int), intent(in) :: lmaxt
        integer(kind=c_int), intent(in) :: k
        integer(kind=c_int), intent(in) :: taper_wt_d0
        real(kind=c_double), optional,dimension(taper_wt_d0),intent(in) :: taper_wt
        integer(kind=c_int), optional,intent(in) :: csphase
        integer(kind=c_int), optional,intent(in) :: norm
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call SHMultiTaperMaskCSE(mtse,sd,sh1,lmax1,sh2,lmax2,tapers,lmaxt,k,taper_wt=taper_wt&
                                     ,norm=norm,csphase=csphase,exitstatus=exitstatus)
    end subroutine cbindSHMultiTaperMaskCSE

    subroutine cbindComputeDMap(Dij,Dij_d0,Dij_d1,dh_mask,dh_mask_d0,dh_mask_d1,n_dh&
                                   ,lmax,sampling,degrees,degrees_d0,exitstatus)  bind(c&
                                   , name="cbind_compute_d_map")
        use, intrinsic :: iso_c_binding
        use shtools, only: ComputeDMap
        implicit none
        integer(kind=c_int), intent(in) :: Dij_d0
        integer(kind=c_int), intent(in) :: Dij_d1
        real(kind=c_double), dimension(Dij_d0,Dij_d1),intent(out) :: Dij
        integer(kind=c_int), intent(in) :: dh_mask_d0
        integer(kind=c_int), intent(in) :: dh_mask_d1
        integer(kind=c_int), dimension(dh_mask_d0,dh_mask_d1),intent(in) :: dh_mask
        integer(kind=c_int), intent(in) :: n_dh
        integer(kind=c_int), intent(in) :: lmax
        integer(kind=c_int), optional,intent(in) :: sampling
        integer(kind=c_int), intent(in) :: degrees_d0
        integer(kind=c_int), optional,dimension(degrees_d0),intent(in) :: degrees
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call ComputeDMap(Dij,dh_mask,n_dh,lmax,sampling=sampling,degrees=degrees,exitstatus=exitstatus)
    end subroutine cbindComputeDMap

    subroutine cbindCurve2Mask(dhgrid,dhgrid_d0,dhgrid_d1,n,sampling,profile,profile_d0&
                                     ,profile_d1,nprofile,NP,extend,exitstatus)  bind(c&
                                     , name="cbind_curve2_mask")
        use, intrinsic :: iso_c_binding
        use shtools, only: Curve2Mask
        implicit none
        integer(kind=c_int), intent(in) :: dhgrid_d0
        integer(kind=c_int), intent(in) :: dhgrid_d1
        integer(kind=c_int), dimension(dhgrid_d0,dhgrid_d1),intent(out) :: dhgrid
        integer(kind=c_int), intent(in) :: profile_d0
        integer(kind=c_int), intent(in) :: profile_d1
        real(kind=c_double), dimension(profile_d0,profile_d1),intent(in) :: profile
        integer(kind=c_int), intent(in) :: n
        integer(kind=c_int), intent(in) :: sampling
        integer(kind=c_int), intent(in) :: nprofile
        integer(kind=c_int), optional,intent(in) :: extend
        integer(kind=c_int), optional,intent(out) :: exitstatus
        integer(kind=c_int) :: NP
        call Curve2Mask(dhgrid,n,sampling,profile,nprofile,NP,extend=extend,exitstatus=exitstatus)
    end subroutine cbindCurve2Mask

    subroutine cbindSHBias(Shh,Shh_d0,lwin,incspectra,incspectra_d0,ldata,outcspectra&
                              ,outcspectra_d0,save_cg,exitstatus)  bind(c, name="cbind_sh_bias")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHBias
        implicit none
        integer(kind=c_int), intent(in) :: Shh_d0
        real(kind=c_double), dimension(Shh_d0),intent(in) :: Shh
        integer(kind=c_int), intent(in) :: incspectra_d0
        real(kind=c_double), dimension(incspectra_d0),intent(in) :: incspectra
        integer(kind=c_int), intent(in) :: outcspectra_d0
        real(kind=c_double), dimension(outcspectra_d0),intent(out) :: outcspectra
        integer(kind=c_int), intent(in) :: lwin
        integer(kind=c_int), intent(in) :: ldata
        integer(kind=c_int), optional,intent(in) :: save_cg
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call SHBias(Shh,lwin,incspectra,ldata,outcspectra,save_cg=save_cg,exitstatus=exitstatus)
    end subroutine cbindSHBias

    subroutine cbindSphericalCapCoef(coef,coef_d0,theta,lmax,exitstatus)  bind(c, name="cbind_spherical_cap_coef")
        use, intrinsic :: iso_c_binding
        use shtools, only: SphericalCapCoef
        implicit none
        integer(kind=c_int), intent(in) :: coef_d0
        real(kind=c_double), dimension(coef_d0),intent(out) :: coef
        real(kind=c_double), intent(in) :: theta
        integer(kind=c_int), optional,intent(in) :: lmax
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call SphericalCapCoef(coef,theta,lmax=lmax,exitstatus=exitstatus)
    end subroutine cbindSphericalCapCoef

    subroutine cbindMakeGravGridDH(cilm,cilm_d,lmax,gm,r0,a,f,rad,rad_d0,rad_d1,theta&
                                       ,theta_d0,theta_d1,phi,phi_d0,phi_d1,total&
                                       ,total_d0,total_d1,n,sampling,lmax_calc,omega&
                                       ,normal_gravity,pot,pot_d0,pot_d1,extend,exitstatus)  bind(c&
                                       , name="cbind_make_grav_grid_dh")
        use, intrinsic :: iso_c_binding
        use shtools, only: MakeGravGridDH
        implicit none
        integer(kind=c_int), intent(in) :: cilm_d
        real(kind=c_double), dimension(2,cilm_d,cilm_d),intent(in) :: cilm
        real(kind=c_double), intent(in) :: gm
        real(kind=c_double), intent(in) :: r0
        real(kind=c_double), intent(in) :: a
        real(kind=c_double), intent(in) :: f
        integer(kind=c_int), intent(in) :: rad_d0
        integer(kind=c_int), intent(in) :: rad_d1
        real(kind=c_double), dimension(rad_d0,rad_d1),intent(out) :: rad
        integer(kind=c_int), intent(in) :: theta_d0
        integer(kind=c_int), intent(in) :: theta_d1
        real(kind=c_double), dimension(theta_d0,theta_d1),intent(out) :: theta
        integer(kind=c_int), intent(in) :: phi_d0
        integer(kind=c_int), intent(in) :: phi_d1
        real(kind=c_double), dimension(phi_d0,phi_d1),intent(out) :: phi
        integer(kind=c_int), intent(in) :: total_d0
        integer(kind=c_int), intent(in) :: total_d1
        real(kind=c_double), dimension(total_d0,total_d1),intent(out) :: total
        real(kind=c_double), optional,intent(in) :: omega
        integer(kind=c_int), intent(in) :: pot_d0
        integer(kind=c_int), intent(in) :: pot_d1
        real(kind=c_double), optional,dimension(pot_d0,pot_d1),intent(out) :: pot
        integer(kind=c_int), intent(in) :: lmax
        integer(kind=c_int), intent(out) :: n
        integer(kind=c_int), optional,intent(in) :: sampling
        integer(kind=c_int), optional,intent(in) :: lmax_calc
        integer(kind=c_int), optional,intent(in) :: normal_gravity
        integer(kind=c_int), optional,intent(in) :: extend
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call MakeGravGridDH(cilm,lmax,gm,r0,a,f,rad,theta,phi,total,n,sampling=sampling&
                                ,lmax_calc=lmax_calc,omega=omega,normal_gravity=normal_gravity&
                                ,pot=pot,extend=extend,exitstatus=exitstatus)
    end subroutine cbindMakeGravGridDH

    subroutine cbindMakeGravGradGridDH(cilm,cilm_d,lmax,gm,r0,a,f,vxx,vxx_d0,vxx_d1&
                                           ,vyy,vyy_d0,vyy_d1,vzz,vzz_d0,vzz_d1,vxy&
                                           ,vxy_d0,vxy_d1,vxz,vxz_d0,vxz_d1,vyz,vyz_d0&
                                           ,vyz_d1,n,sampling,lmax_calc,extend,exitstatus)  bind(c&
                                           , name="cbind_make_grav_grad_grid_dh")
        use, intrinsic :: iso_c_binding
        use shtools, only: MakeGravGradGridDH
        implicit none
        integer(kind=c_int), intent(in) :: cilm_d
        real(kind=c_double), dimension(2,cilm_d,cilm_d),intent(in) :: cilm
        real(kind=c_double), intent(in) :: gm
        real(kind=c_double), intent(in) :: r0
        real(kind=c_double), intent(in) :: a
        real(kind=c_double), intent(in) :: f
        integer(kind=c_int), intent(in) :: vxx_d0
        integer(kind=c_int), intent(in) :: vxx_d1
        real(kind=c_double), dimension(vxx_d0,vxx_d1),intent(out) :: vxx
        integer(kind=c_int), intent(in) :: vyy_d0
        integer(kind=c_int), intent(in) :: vyy_d1
        real(kind=c_double), dimension(vyy_d0,vyy_d1),intent(out) :: vyy
        integer(kind=c_int), intent(in) :: vzz_d0
        integer(kind=c_int), intent(in) :: vzz_d1
        real(kind=c_double), dimension(vzz_d0,vzz_d1),intent(out) :: vzz
        integer(kind=c_int), intent(in) :: vxy_d0
        integer(kind=c_int), intent(in) :: vxy_d1
        real(kind=c_double), dimension(vxy_d0,vxy_d1),intent(out) :: vxy
        integer(kind=c_int), intent(in) :: vxz_d0
        integer(kind=c_int), intent(in) :: vxz_d1
        real(kind=c_double), dimension(vxz_d0,vxz_d1),intent(out) :: vxz
        integer(kind=c_int), intent(in) :: vyz_d0
        integer(kind=c_int), intent(in) :: vyz_d1
        real(kind=c_double), dimension(vyz_d0,vyz_d1),intent(out) :: vyz
        integer(kind=c_int), intent(in) :: lmax
        integer(kind=c_int), intent(out) :: n
        integer(kind=c_int), optional,intent(in) :: sampling
        integer(kind=c_int), optional,intent(in) :: lmax_calc
        integer(kind=c_int), optional,intent(in) :: extend
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call MakeGravGradGridDH(cilm,lmax,gm,r0,a,f,vxx,vyy,vzz,vxy,vxz,vyz,n,sampling=sampling&
                                    ,lmax_calc=lmax_calc,extend=extend,exitstatus=exitstatus)
    end subroutine cbindMakeGravGradGridDH

    subroutine cbindMakeMagGradGridDH(cilm,cilm_d,lmax,r0,a,f,vxx,vxx_d0,vxx_d1,vyy&
                                          ,vyy_d0,vyy_d1,vzz,vzz_d0,vzz_d1,vxy,vxy_d0&
                                          ,vxy_d1,vxz,vxz_d0,vxz_d1,vyz,vyz_d0,vyz_d1&
                                          ,n,sampling,lmax_calc,extend,exitstatus)  bind(c&
                                          , name="cbind_make_mag_grad_grid_dh")
        use, intrinsic :: iso_c_binding
        use shtools, only: MakeMagGradGridDH
        implicit none
        integer(kind=c_int), intent(in) :: cilm_d
        real(kind=c_double), dimension(2,cilm_d,cilm_d),intent(in) :: cilm
        real(kind=c_double), intent(in) :: r0
        real(kind=c_double), intent(in) :: a
        real(kind=c_double), intent(in) :: f
        integer(kind=c_int), intent(in) :: vxx_d0
        integer(kind=c_int), intent(in) :: vxx_d1
        real(kind=c_double), dimension(vxx_d0,vxx_d1),intent(out) :: vxx
        integer(kind=c_int), intent(in) :: vyy_d0
        integer(kind=c_int), intent(in) :: vyy_d1
        real(kind=c_double), dimension(vyy_d0,vyy_d1),intent(out) :: vyy
        integer(kind=c_int), intent(in) :: vzz_d0
        integer(kind=c_int), intent(in) :: vzz_d1
        real(kind=c_double), dimension(vzz_d0,vzz_d1),intent(out) :: vzz
        integer(kind=c_int), intent(in) :: vxy_d0
        integer(kind=c_int), intent(in) :: vxy_d1
        real(kind=c_double), dimension(vxy_d0,vxy_d1),intent(out) :: vxy
        integer(kind=c_int), intent(in) :: vxz_d0
        integer(kind=c_int), intent(in) :: vxz_d1
        real(kind=c_double), dimension(vxz_d0,vxz_d1),intent(out) :: vxz
        integer(kind=c_int), intent(in) :: vyz_d0
        integer(kind=c_int), intent(in) :: vyz_d1
        real(kind=c_double), dimension(vyz_d0,vyz_d1),intent(out) :: vyz
        integer(kind=c_int), intent(in) :: lmax
        integer(kind=c_int), intent(out) :: n
        integer(kind=c_int), optional,intent(in) :: sampling
        integer(kind=c_int), optional,intent(in) :: lmax_calc
        integer(kind=c_int), optional,intent(in) :: extend
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call MakeMagGradGridDH(cilm,lmax,r0,a,f,vxx,vyy,vzz,vxy,vxz,vyz,n,sampling=sampling&
                                   ,lmax_calc=lmax_calc,extend=extend,exitstatus=exitstatus)
    end subroutine cbindMakeMagGradGridDH

    subroutine cbindMakeGeoidGrid(geoid,geoid_d0,geoid_d1,cilm,cilm_d,lmax,r0pot,GM&
                                       ,PotRef,omega,r,gridtype,order,nlat,nlong,interval&
                                       ,lmax_calc,a,f,extend,exitstatus)  bind(c, name="cbind_make_geoid_grid")
        use, intrinsic :: iso_c_binding
        use shtools, only: MakeGeoidGrid
        implicit none
        integer(kind=c_int), intent(in) :: geoid_d0
        integer(kind=c_int), intent(in) :: geoid_d1
        real(kind=c_double), dimension(geoid_d0,geoid_d1),intent(out) :: geoid
        integer(kind=c_int), intent(in) :: cilm_d
        real(kind=c_double), dimension(2,cilm_d,cilm_d),intent(in) :: cilm
        real(kind=c_double), intent(in) :: r0pot
        real(kind=c_double), intent(in) :: GM
        real(kind=c_double), intent(in) :: r
        real(kind=c_double), intent(in) :: PotRef
        real(kind=c_double), intent(in) :: omega
        integer(kind=c_int), intent(in) :: lmax
        integer(kind=c_int), intent(in) :: order
        integer(kind=c_int), intent(in) :: gridtype
        integer(kind=c_int), optional,intent(in) :: lmax_calc
        integer(kind=c_int), optional,intent(in) :: extend
        integer(kind=c_int), intent(out) :: nlat
        integer(kind=c_int), intent(out) :: nlong
        real(kind=c_double), optional,intent(in) :: interval
        real(kind=c_double), optional,intent(in) :: a
        real(kind=c_double), optional,intent(in) :: f
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call MakeGeoidGrid(geoid,cilm,lmax,r0pot,GM,PotRef,omega,r,gridtype,order&
                                ,nlat,nlong,interval=interval,lmax_calc=lmax_calc&
                                ,a=a,f=f,extend=extend,exitstatus=exitstatus)
    end subroutine cbindMakeGeoidGrid

    subroutine cbindCilmPlus(cilm,cilm_d,gridin,gridin_d0,gridin_d1,lmax,nmax,mass&
                                 ,d,rho,gridtype,w,w_d0,zero,zero_d0,plx,plx_d0,plx_d1&
                                 ,n,dref,exitstatus)  bind(c, name="cbind_cilm_plus")
        use, intrinsic :: iso_c_binding
        use shtools, only: CilmPlus
        implicit none
        integer(kind=c_int), intent(in) :: gridin_d0
        integer(kind=c_int), intent(in) :: gridin_d1
        real(kind=c_double), dimension(gridin_d0,gridin_d1),intent(in) :: gridin
        real(kind=c_double), intent(in) :: mass
        real(kind=c_double), intent(in) :: rho
        integer(kind=c_int), intent(in) :: w_d0
        real(kind=c_double), optional,dimension(w_d0),intent(in) :: w
        integer(kind=c_int), intent(in) :: zero_d0
        real(kind=c_double), optional,dimension(zero_d0),intent(in) :: zero
        integer(kind=c_int), intent(in) :: plx_d0
        integer(kind=c_int), intent(in) :: plx_d1
        real(kind=c_double), optional,dimension(plx_d0,plx_d1),intent(in) :: plx
        real(kind=c_double), optional,intent(in) :: dref
        integer(kind=c_int), intent(in) :: cilm_d
        real(kind=c_double), dimension(2,cilm_d,cilm_d),intent(out) :: cilm
        real(kind=c_double), intent(out) :: d
        integer(kind=c_int), intent(in) :: lmax
        integer(kind=c_int), intent(in) :: nmax
        integer(kind=c_int), intent(in) :: gridtype
        integer(kind=c_int), optional,intent(in) :: n
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call CilmPlus(cilm,gridin,lmax,nmax,mass,d,rho,gridtype,w=w,zero=zero,plx=plx&
                          ,n=n,dref=dref,exitstatus=exitstatus)
    end subroutine cbindCilmPlus

    subroutine cbindCilmMinus(cilm,cilm_d,gridin,gridin_d0,gridin_d1,lmax,nmax,mass&
                                  ,d,rho,gridtype,w,w_d0,zero,zero_d0,plx,plx_d0,plx_d1&
                                  ,n,dref,exitstatus)  bind(c, name="cbind_cilm_minus")
        use, intrinsic :: iso_c_binding
        use shtools, only: CilmMinus
        implicit none
        integer(kind=c_int), intent(in) :: gridin_d0
        integer(kind=c_int), intent(in) :: gridin_d1
        real(kind=c_double), dimension(gridin_d0,gridin_d1),intent(in) :: gridin
        real(kind=c_double), intent(in) :: mass
        real(kind=c_double), intent(in) :: rho
        integer(kind=c_int), intent(in) :: w_d0
        real(kind=c_double), optional,dimension(w_d0),intent(in) :: w
        integer(kind=c_int), intent(in) :: zero_d0
        real(kind=c_double), optional,dimension(zero_d0),intent(in) :: zero
        integer(kind=c_int), intent(in) :: plx_d0
        integer(kind=c_int), intent(in) :: plx_d1
        real(kind=c_double), optional,dimension(plx_d0,plx_d1),intent(in) :: plx
        real(kind=c_double), optional,intent(in) :: dref
        integer(kind=c_int), intent(in) :: cilm_d
        real(kind=c_double), dimension(2,cilm_d,cilm_d),intent(out) :: cilm
        real(kind=c_double), intent(out) :: d
        integer(kind=c_int), intent(in) :: lmax
        integer(kind=c_int), intent(in) :: nmax
        integer(kind=c_int), intent(in) :: gridtype
        integer(kind=c_int), optional,intent(in) :: n
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call CilmMinus(cilm,gridin,lmax,nmax,mass,d,rho,gridtype,w=w,zero=zero,plx=plx&
                           ,n=n,dref=dref,exitstatus=exitstatus)
    end subroutine cbindCilmMinus

    subroutine cbindCilmPlusRhoH(cilm,cilm_d,gridin,gridin_d0,gridin_d1,lmax,nmax&
                                     ,mass,d,rho,rho_d0,rho_d1,gridtype,w,w_d0,zero&
                                     ,zero_d0,plx,plx_d0,plx_d1,n,dref,exitstatus)  bind(c&
                                     , name="cbind_cilm_plus_rho_h")
        use, intrinsic :: iso_c_binding
        use shtools, only: CilmPlusRhoH
        implicit none
        integer(kind=c_int), intent(in) :: gridin_d0
        integer(kind=c_int), intent(in) :: gridin_d1
        real(kind=c_double), dimension(gridin_d0,gridin_d1),intent(in) :: gridin
        real(kind=c_double), intent(in) :: mass
        integer(kind=c_int), intent(in) :: rho_d0
        integer(kind=c_int), intent(in) :: rho_d1
        real(kind=c_double), dimension(rho_d0,rho_d1),intent(in) :: rho
        integer(kind=c_int), intent(in) :: w_d0
        real(kind=c_double), optional,dimension(w_d0),intent(in) :: w
        integer(kind=c_int), intent(in) :: zero_d0
        real(kind=c_double), optional,dimension(zero_d0),intent(in) :: zero
        integer(kind=c_int), intent(in) :: plx_d0
        integer(kind=c_int), intent(in) :: plx_d1
        real(kind=c_double), optional,dimension(plx_d0,plx_d1),intent(in) :: plx
        real(kind=c_double), optional,intent(in) :: dref
        integer(kind=c_int), intent(in) :: cilm_d
        real(kind=c_double), dimension(2,cilm_d,cilm_d),intent(out) :: cilm
        real(kind=c_double), intent(out) :: d
        integer(kind=c_int), intent(in) :: lmax
        integer(kind=c_int), intent(in) :: nmax
        integer(kind=c_int), intent(in) :: gridtype
        integer(kind=c_int), optional,intent(in) :: n
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call CilmPlusRhoH(cilm,gridin,lmax,nmax,mass,d,rho,gridtype,w=w,zero=zero&
                              ,plx=plx,n=n,dref=dref,exitstatus=exitstatus)
    end subroutine cbindCilmPlusRhoH

    subroutine cbindCilmMinusRhoH(cilm,cilm_d,gridin,gridin_d0,gridin_d1,lmax,nmax&
                                      ,mass,d,rho,rho_d0,rho_d1,gridtype,w,w_d0,zero&
                                      ,zero_d0,plx,plx_d0,plx_d1,n,dref,exitstatus)  bind(c&
                                      , name="cbind_cilm_minus_rho_h")
        use, intrinsic :: iso_c_binding
        use shtools, only: CilmMinusRhoH
        implicit none
        integer(kind=c_int), intent(in) :: gridin_d0
        integer(kind=c_int), intent(in) :: gridin_d1
        real(kind=c_double), dimension(gridin_d0,gridin_d1),intent(in) :: gridin
        real(kind=c_double), intent(in) :: mass
        integer(kind=c_int), intent(in) :: rho_d0
        integer(kind=c_int), intent(in) :: rho_d1
        real(kind=c_double), dimension(rho_d0,rho_d1),intent(in) :: rho
        integer(kind=c_int), intent(in) :: w_d0
        real(kind=c_double), optional,dimension(w_d0),intent(in) :: w
        integer(kind=c_int), intent(in) :: zero_d0
        real(kind=c_double), optional,dimension(zero_d0),intent(in) :: zero
        integer(kind=c_int), intent(in) :: plx_d0
        integer(kind=c_int), intent(in) :: plx_d1
        real(kind=c_double), optional,dimension(plx_d0,plx_d1),intent(in) :: plx
        real(kind=c_double), optional,intent(in) :: dref
        integer(kind=c_int), intent(in) :: cilm_d
        real(kind=c_double), dimension(2,cilm_d,cilm_d),intent(out) :: cilm
        real(kind=c_double), intent(out) :: d
        integer(kind=c_int), intent(in) :: lmax
        integer(kind=c_int), intent(in) :: nmax
        integer(kind=c_int), intent(in) :: gridtype
        integer(kind=c_int), optional,intent(in) :: n
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call CilmMinusRhoH(cilm,gridin,lmax,nmax,mass,d,rho,gridtype,w=w,zero=zero&
                               ,plx=plx,n=n,dref=dref,exitstatus=exitstatus)
    end subroutine cbindCilmMinusRhoH

    subroutine cbindBAtoHilm(cilm,cilm_d,ba,ba_d0,ba_d1,ba_d2,gridglq,gridglq_d0,gridglq_d1&
                                 ,lmax,nmax,mass,r0,rho,gridtype,w,w_d0,plx,plx_d0&
                                 ,plx_d1,zero,zero_d0,filter_type,filter_deg,lmax_calc&
                                 ,exitstatus)  bind(c, name="cbind_b_ato_hilm")
        use, intrinsic :: iso_c_binding
        use shtools, only: BAtoHilm
        implicit none
        integer(kind=c_int), intent(in) :: cilm_d
        real(kind=c_double), dimension(2,cilm_d,cilm_d),intent(out) :: cilm
        integer(kind=c_int), intent(in) :: ba_d0
        integer(kind=c_int), intent(in) :: ba_d1
        integer(kind=c_int), intent(in) :: ba_d2
        real(kind=c_double), dimension(ba_d0,ba_d1,ba_d2),intent(in) :: ba
        integer(kind=c_int), intent(in) :: gridglq_d0
        integer(kind=c_int), intent(in) :: gridglq_d1
        real(kind=c_double), dimension(gridglq_d0,gridglq_d1),intent(in) :: gridglq
        real(kind=c_double), intent(in) :: mass
        real(kind=c_double), intent(in) :: r0
        real(kind=c_double), intent(in) :: rho
        integer(kind=c_int), intent(in) :: plx_d0
        integer(kind=c_int), intent(in) :: plx_d1
        real(kind=c_double), optional,dimension(plx_d0,plx_d1),intent(in) :: plx
        integer(kind=c_int), intent(in) :: zero_d0
        real(kind=c_double), optional,dimension(zero_d0),intent(in) :: zero
        integer(kind=c_int), intent(in) :: w_d0
        real(kind=c_double), optional,dimension(w_d0),intent(in) :: w
        integer(kind=c_int), intent(in) :: lmax
        integer(kind=c_int), intent(in) :: nmax
        integer(kind=c_int), intent(in) :: gridtype
        integer(kind=c_int), optional,intent(in) :: filter_type
        integer(kind=c_int), optional,intent(in) :: filter_deg
        integer(kind=c_int), optional,intent(in) :: lmax_calc
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call BAtoHilm(cilm,ba,gridglq,lmax,nmax,mass,r0,rho,gridtype,w=w,plx=plx,zero=zero&
                          ,filter_type=filter_type,filter_deg=filter_deg,lmax_calc=lmax_calc&
                          ,exitstatus=exitstatus)
    end subroutine cbindBAtoHilm

    subroutine cbindBAtoHilmRhoH(cilm,cilm_d,ba,ba_d0,ba_d1,ba_d2,gridglq,gridglq_d0&
                                     ,gridglq_d1,lmax,nmax,mass,r0,rho,rho_d0,rho_d1&
                                     ,gridtype,w,w_d0,plx,plx_d0,plx_d1,zero,zero_d0&
                                     ,filter_type,filter_deg,lmax_calc,exitstatus)  bind(c&
                                     , name="cbind_b_ato_hilm_rho_h")
        use, intrinsic :: iso_c_binding
        use shtools, only: BAtoHilmRhoH
        implicit none
        integer(kind=c_int), intent(in) :: cilm_d
        real(kind=c_double), dimension(2,cilm_d,cilm_d),intent(out) :: cilm
        integer(kind=c_int), intent(in) :: ba_d0
        integer(kind=c_int), intent(in) :: ba_d1
        integer(kind=c_int), intent(in) :: ba_d2
        real(kind=c_double), dimension(ba_d0,ba_d1,ba_d2),intent(in) :: ba
        integer(kind=c_int), intent(in) :: gridglq_d0
        integer(kind=c_int), intent(in) :: gridglq_d1
        real(kind=c_double), dimension(gridglq_d0,gridglq_d1),intent(in) :: gridglq
        real(kind=c_double), intent(in) :: mass
        real(kind=c_double), intent(in) :: r0
        integer(kind=c_int), intent(in) :: rho_d0
        integer(kind=c_int), intent(in) :: rho_d1
        real(kind=c_double), dimension(rho_d0,rho_d1),intent(in) :: rho
        integer(kind=c_int), intent(in) :: plx_d0
        integer(kind=c_int), intent(in) :: plx_d1
        real(kind=c_double), optional,dimension(plx_d0,plx_d1),intent(in) :: plx
        integer(kind=c_int), intent(in) :: zero_d0
        real(kind=c_double), optional,dimension(zero_d0),intent(in) :: zero
        integer(kind=c_int), intent(in) :: w_d0
        real(kind=c_double), optional,dimension(w_d0),intent(in) :: w
        integer(kind=c_int), intent(in) :: lmax
        integer(kind=c_int), intent(in) :: nmax
        integer(kind=c_int), intent(in) :: gridtype
        integer(kind=c_int), optional,intent(in) :: filter_type
        integer(kind=c_int), optional,intent(in) :: filter_deg
        integer(kind=c_int), optional,intent(in) :: lmax_calc
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call BAtoHilmRhoH(cilm,ba,gridglq,lmax,nmax,mass,r0,rho,gridtype,w=w,plx=plx&
                              ,zero=zero,filter_type=filter_type,filter_deg=filter_deg&
                              ,lmax_calc=lmax_calc,exitstatus=exitstatus)
    end subroutine cbindBAtoHilmRhoH

    function cbindDownContFilterMA(l,half,r,d)  bind(c, name="cbind_down_cont_filter_ma")
        use, intrinsic :: iso_c_binding
        use shtools, only: DownContFilterMA
        implicit none
        real(kind=c_double) :: cbindDownContFilterMA
        integer(kind=c_int), intent(in) :: l
        integer(kind=c_int), intent(in) :: half
        real(kind=c_double), intent(in) :: r
        real(kind=c_double), intent(in) :: d
        cbindDownContFilterMA=DownContFilterMA(l,half,r,d)
    end function cbindDownContFilterMA

    function cbindDownContFilterMC(l,half,r,d)  bind(c, name="cbind_down_cont_filter_mc")
        use, intrinsic :: iso_c_binding
        use shtools, only: DownContFilterMC
        implicit none
        real(kind=c_double) :: cbindDownContFilterMC
        integer(kind=c_int), intent(in) :: l
        integer(kind=c_int), intent(in) :: half
        real(kind=c_double), intent(in) :: r
        real(kind=c_double), intent(in) :: d
        cbindDownContFilterMC=DownContFilterMC(l,half,r,d)
    end function cbindDownContFilterMC

    function cbindNormalGravity(geocentric_lat,gm,omega,a,b)  bind(c, name="cbind_normal_gravity")
        use, intrinsic :: iso_c_binding
        use shtools, only: NormalGravity
        implicit none
        real(kind=c_double) :: cbindNormalGravity
        real(kind=c_double), intent(in) :: geocentric_lat
        real(kind=c_double), intent(in) :: gm
        real(kind=c_double), intent(in) :: omega
        real(kind=c_double), intent(in) :: a
        real(kind=c_double), intent(in) :: b
        cbindNormalGravity=NormalGravity(geocentric_lat,gm,omega,a,b)
    end function cbindNormalGravity

    subroutine cbindMakeMagGridDH(cilm,cilm_d,lmax,r0,a,f,rad_grid,rad_grid_d0,rad_grid_d1&
                                      ,theta_grid,theta_grid_d0,theta_grid_d1,phi_grid&
                                      ,phi_grid_d0,phi_grid_d1,total_grid,total_grid_d0&
                                      ,total_grid_d1,n,sampling,lmax_calc,pot_grid&
                                      ,pot_grid_d0,pot_grid_d1,extend,exitstatus)  bind(c&
                                      , name="cbind_make_mag_grid_dh")
        use, intrinsic :: iso_c_binding
        use shtools, only: MakeMagGridDH
        implicit none
        integer(kind=c_int), intent(in) :: cilm_d
        real(kind=c_double), dimension(2,cilm_d,cilm_d),intent(in) :: cilm
        real(kind=c_double), intent(in) :: r0
        real(kind=c_double), intent(in) :: a
        real(kind=c_double), intent(in) :: f
        integer(kind=c_int), intent(in) :: rad_grid_d0
        integer(kind=c_int), intent(in) :: rad_grid_d1
        real(kind=c_double), dimension(rad_grid_d0,rad_grid_d1),intent(out) :: rad_grid
        integer(kind=c_int), intent(in) :: theta_grid_d0
        integer(kind=c_int), intent(in) :: theta_grid_d1
        real(kind=c_double), dimension(theta_grid_d0,theta_grid_d1),intent(out) :: theta_grid
        integer(kind=c_int), intent(in) :: phi_grid_d0
        integer(kind=c_int), intent(in) :: phi_grid_d1
        real(kind=c_double), dimension(phi_grid_d0,phi_grid_d1),intent(out) :: phi_grid
        integer(kind=c_int), intent(in) :: total_grid_d0
        integer(kind=c_int), intent(in) :: total_grid_d1
        real(kind=c_double), dimension(total_grid_d0,total_grid_d1),intent(out) :: total_grid
        integer(kind=c_int), intent(in) :: pot_grid_d0
        integer(kind=c_int), intent(in) :: pot_grid_d1
        real(kind=c_double), optional,dimension(pot_grid_d0,pot_grid_d1),intent(out) :: pot_grid
        integer(kind=c_int), intent(in) :: lmax
        integer(kind=c_int), intent(out) :: n
        integer(kind=c_int), optional,intent(in) :: sampling
        integer(kind=c_int), optional,intent(in) :: lmax_calc
        integer(kind=c_int), optional,intent(in) :: extend
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call MakeMagGridDH(cilm,lmax,r0,a,f,rad_grid,theta_grid,phi_grid,total_grid&
                               ,n,sampling=sampling,lmax_calc=lmax_calc,pot_grid=pot_grid&
                               ,extend=extend,exitstatus=exitstatus)
    end subroutine cbindMakeMagGridDH

    subroutine cbindSHMagPowerSpectrum(c,c_d0,c_d1,c_d2,a,r,lmax,spectra,spectra_d0&
                                        ,exitstatus)  bind(c, name="cbind_sh_mag_power_spectrum")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHMagPowerSpectrum
        implicit none
        integer(kind=c_int), intent(in) :: c_d0
        integer(kind=c_int), intent(in) :: c_d1
        integer(kind=c_int), intent(in) :: c_d2
        real(kind=c_double), dimension(c_d0,c_d1,c_d2),intent(in) :: c
        real(kind=c_double), intent(in) :: a
        real(kind=c_double), intent(in) :: r
        integer(kind=c_int), intent(in) :: lmax
        integer(kind=c_int), intent(in) :: spectra_d0
        real(kind=c_double), dimension(spectra_d0),intent(out) :: spectra
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call SHMagPowerSpectrum(c,a,r,lmax,spectra,exitstatus=exitstatus)
    end subroutine cbindSHMagPowerSpectrum

    function cbindSHMagPowerL(c,c_d0,c_d1,c_d2,a,r,l)  bind(c, name="cbind_sh_mag_power_l")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHMagPowerL
        implicit none
        real(kind=c_double) :: cbindSHMagPowerL
        integer(kind=c_int), intent(in) :: c_d0
        integer(kind=c_int), intent(in) :: c_d1
        integer(kind=c_int), intent(in) :: c_d2
        real(kind=c_double), dimension(c_d0,c_d1,c_d2),intent(in) :: c
        real(kind=c_double), intent(in) :: a
        real(kind=c_double), intent(in) :: r
        integer(kind=c_int), intent(in) :: l
        cbindSHMagPowerL=SHMagPowerL(c,a,r,l)
    end function cbindSHMagPowerL

    subroutine cbindMakeCircleCoord(coord,coord_d0,coord_d1,lat,lon,theta0,cinterval&
                                         ,cnum,exitstatus)  bind(c, name="cbind_make_circle_coord")
        use, intrinsic :: iso_c_binding
        use shtools, only: MakeCircleCoord
        implicit none
        real(kind=c_double), intent(in) :: lat
        real(kind=c_double), intent(in) :: lon
        real(kind=c_double), intent(in) :: theta0
        integer(kind=c_int), intent(in) :: coord_d0
        integer(kind=c_int), intent(in) :: coord_d1
        real(kind=c_double), dimension(coord_d0,coord_d1),intent(out) :: coord
        real(kind=c_double), optional,intent(in) :: cinterval
        integer(kind=c_int), optional,intent(out) :: cnum
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call MakeCircleCoord(coord,lat,lon,theta0,cinterval=cinterval,cnum=cnum,exitstatus=exitstatus)
    end subroutine cbindMakeCircleCoord

    subroutine cbindMakeEllipseCoord(coord,coord_d0,coord_d1,lat,lon,dec,A_theta,B_theta&
                                          ,cinterval,cnum,exitstatus)  bind(c, name="cbind_make_ellipse_coord")
        use, intrinsic :: iso_c_binding
        use shtools, only: MakeEllipseCoord
        implicit none
        real(kind=c_double), intent(in) :: lat
        real(kind=c_double), intent(in) :: lon
        real(kind=c_double), intent(in) :: A_theta
        real(kind=c_double), intent(in) :: B_theta
        real(kind=c_double), intent(in) :: dec
        integer(kind=c_int), intent(in) :: coord_d0
        integer(kind=c_int), intent(in) :: coord_d1
        real(kind=c_double), dimension(coord_d0,coord_d1),intent(out) :: coord
        real(kind=c_double), optional,intent(in) :: cinterval
        integer(kind=c_int), optional,intent(out) :: cnum
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call MakeEllipseCoord(coord,lat,lon,dec,A_theta,B_theta,cinterval=cinterval&
                                   ,cnum=cnum,exitstatus=exitstatus)
    end subroutine cbindMakeEllipseCoord

    subroutine cbindWigner3j(w3j,w3j_d0,jmin,jmax,j2,j3,m1,m2,m3,exitstatus)  bind(c&
                                , name="cbind_wigner3j")
        use, intrinsic :: iso_c_binding
        use shtools, only: Wigner3j
        implicit none
        integer(kind=c_int), intent(in) :: j2
        integer(kind=c_int), intent(in) :: j3
        integer(kind=c_int), intent(in) :: m1
        integer(kind=c_int), intent(in) :: m2
        integer(kind=c_int), intent(in) :: m3
        integer(kind=c_int), intent(out) :: jmin
        integer(kind=c_int), intent(out) :: jmax
        integer(kind=c_int), intent(in) :: w3j_d0
        real(kind=c_double), dimension(w3j_d0),intent(out) :: w3j
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call Wigner3j(w3j,jmin,jmax,j2,j3,m1,m2,m3,exitstatus=exitstatus)
    end subroutine cbindWigner3j

    function cbindRandomN(idum)  bind(c, name="cbind_random_n")
        use, intrinsic :: iso_c_binding
        use shtools, only: RandomN
        implicit none
        real(kind=c_double) :: cbindRandomN
        integer(kind=c_int), intent(inout) :: idum
        cbindRandomN=RandomN(idum)
    end function cbindRandomN

    function cbindRandomGaussian(idum)  bind(c, name="cbind_random_gaussian")
        use, intrinsic :: iso_c_binding
        use shtools, only: RandomGaussian
        implicit none
        real(kind=c_double) :: cbindRandomGaussian
        integer(kind=c_int), intent(inout) :: idum
        cbindRandomGaussian=RandomGaussian(idum)
    end function cbindRandomGaussian

    subroutine cbindPreGLQ(x1,x2,n,zero,zero_d0,w,w_d0,exitstatus)  bind(c, name="cbind_pre_glq")
        use, intrinsic :: iso_c_binding
        use shtools, only: PreGLQ
        implicit none
        real(kind=c_double), intent(in) :: x1
        real(kind=c_double), intent(in) :: x2
        integer(kind=c_int), intent(in) :: zero_d0
        real(kind=c_double), dimension(zero_d0),intent(out) :: zero
        integer(kind=c_int), intent(in) :: w_d0
        real(kind=c_double), dimension(w_d0),intent(out) :: w
        integer(kind=c_int), intent(in) :: n
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call PreGLQ(x1,x2,n,zero,w,exitstatus=exitstatus)
    end subroutine cbindPreGLQ

    function cbindNGLQ(degree)  bind(c, name="cbind_nglq")
        use, intrinsic :: iso_c_binding
        use shtools, only: NGLQ
        implicit none
        integer(kind=c_int) :: cbindNGLQ
        integer(kind=c_int), intent(in) :: degree
        cbindNGLQ=NGLQ(degree)
    end function cbindNGLQ

    function cbindNGLQSH(degree)  bind(c, name="cbind_nglqsh")
        use, intrinsic :: iso_c_binding
        use shtools, only: NGLQSH
        implicit none
        integer(kind=c_int) :: cbindNGLQSH
        integer(kind=c_int), intent(in) :: degree
        cbindNGLQSH=NGLQSH(degree)
    end function cbindNGLQSH

    function cbindNGLQSHN(degree,n)  bind(c, name="cbind_nglqshn")
        use, intrinsic :: iso_c_binding
        use shtools, only: NGLQSHN
        implicit none
        integer(kind=c_int) :: cbindNGLQSHN
        integer(kind=c_int), intent(in) :: degree
        integer(kind=c_int), intent(in) :: n
        cbindNGLQSHN=NGLQSHN(degree,n)
    end function cbindNGLQSHN

    subroutine cbindDHaj(n,aj,aj_d0,extend,exitstatus)  bind(c, name="cbind_d_haj")
        use, intrinsic :: iso_c_binding
        use shtools, only: DHaj
        implicit none
        integer(kind=c_int), intent(in) :: n
        integer(kind=c_int), intent(in) :: aj_d0
        real(kind=c_double), dimension(aj_d0),intent(out) :: aj
        integer(kind=c_int), optional,intent(in) :: extend
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call DHaj(n,aj,extend=extend,exitstatus=exitstatus)
    end subroutine cbindDHaj

    function cbindYilmIndexVector(i,l,m)  bind(c, name="cbind_yilm_index_vector")
        use, intrinsic :: iso_c_binding
        use shtools, only: YilmIndexVector
        implicit none
        integer(kind=c_int) :: cbindYilmIndexVector
        integer(kind=c_int), intent(in) :: i
        integer(kind=c_int), intent(in) :: l
        integer(kind=c_int), intent(in) :: m
        cbindYilmIndexVector=YilmIndexVector(i,l,m)
    end function cbindYilmIndexVector

    subroutine cbindEigValVecSym(ain,ain_d0,ain_d1,n,eig,eig_d0,evec,evec_d0,evec_d1&
                                    ,ul,K,exitstatus)  bind(c, name="cbind_eig_val_vec_sym")
        use, intrinsic :: iso_c_binding
        use shtools, only: EigValVecSym
        implicit none
        integer(kind=c_int), intent(in) :: ain_d0
        integer(kind=c_int), intent(in) :: ain_d1
        real(kind=c_double), dimension(ain_d0,ain_d1),intent(in) :: ain
        integer(kind=c_int), intent(in) :: n
        integer(kind=c_int), intent(in) :: eig_d0
        real(kind=c_double), dimension(eig_d0),intent(out) :: eig
        integer(kind=c_int), intent(in) :: evec_d0
        integer(kind=c_int), intent(in) :: evec_d1
        real(kind=c_double), dimension(evec_d0,evec_d1),intent(out) :: evec
        character(kind=c_char), optional,intent(in) :: ul
        integer(kind=c_int), optional,intent(in) :: K
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call EigValVecSym(ain,n,eig,evec,ul=ul,K=K,exitstatus=exitstatus)
    end subroutine cbindEigValVecSym

    subroutine cbindEigValVecSymTri(ain,ain_d0,ain_d1,n,eig,eig_d0,evec,evec_d0,evec_d1&
                                       ,ul,exitstatus)  bind(c, name="cbind_eig_val_vec_sym_tri")
        use, intrinsic :: iso_c_binding
        use shtools, only: EigValVecSymTri
        implicit none
        integer(kind=c_int), intent(in) :: ain_d0
        integer(kind=c_int), intent(in) :: ain_d1
        real(kind=c_double), dimension(ain_d0,ain_d1),intent(in) :: ain
        integer(kind=c_int), intent(in) :: n
        integer(kind=c_int), intent(in) :: eig_d0
        real(kind=c_double), dimension(eig_d0),intent(out) :: eig
        integer(kind=c_int), intent(in) :: evec_d0
        integer(kind=c_int), intent(in) :: evec_d1
        real(kind=c_double), dimension(evec_d0,evec_d1),intent(out) :: evec
        character(kind=c_char), optional,intent(in) :: ul
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call EigValVecSymTri(ain,n,eig,evec,ul=ul,exitstatus=exitstatus)
    end subroutine cbindEigValVecSymTri

    subroutine cbindEigValSym(ain,ain_d0,ain_d1,n,eval,eval_d0,ul)  bind(c, name="cbind_eig_val_sym")
        use, intrinsic :: iso_c_binding
        use shtools, only: EigValSym
        implicit none
        integer(kind=c_int), intent(in) :: ain_d0
        integer(kind=c_int), intent(in) :: ain_d1
        real(kind=c_double), dimension(ain_d0,ain_d1),intent(in) :: ain
        integer(kind=c_int), intent(in) :: n
        integer(kind=c_int), intent(in) :: eval_d0
        real(kind=c_double), dimension(eval_d0),intent(out) :: eval
        character(kind=c_char), optional,intent(in) :: ul
        call EigValSym(ain,n,eval,ul=ul)
    end subroutine cbindEigValSym

    subroutine cbindSHRotateTapers(tapersrot,tapersrot_d0,tapersrot_d1,tapers,tapers_d0&
                                            ,tapers_d1,taper_order,taper_order_d0&
                                            ,lmax,nrot,x,x_d0,dj,dj_d0,dj_d1,dj_d2&
                                            ,exitstatus)  bind(c, name="cbind_sh_rotate_tapers")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHRotateTapers
        implicit none
        integer(kind=c_int), intent(in) :: tapers_d0
        integer(kind=c_int), intent(in) :: tapers_d1
        real(kind=c_double), dimension(tapers_d0,tapers_d1),intent(in) :: tapers
        integer(kind=c_int), intent(in) :: x_d0
        real(kind=c_double), dimension(x_d0),intent(in) :: x
        integer(kind=c_int), intent(in) :: dj_d0
        integer(kind=c_int), intent(in) :: dj_d1
        integer(kind=c_int), intent(in) :: dj_d2
        real(kind=c_double), dimension(dj_d0,dj_d1,dj_d2),intent(in) :: dj
        integer(kind=c_int), intent(in) :: tapersrot_d0
        integer(kind=c_int), intent(in) :: tapersrot_d1
        real(kind=c_double), dimension(tapersrot_d0,tapersrot_d1),intent(out) :: tapersrot
        integer(kind=c_int), intent(in) :: taper_order_d0
        integer(kind=c_int), dimension(taper_order_d0),intent(in) :: taper_order
        integer(kind=c_int), intent(in) :: lmax
        integer(kind=c_int), intent(in) :: nrot
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call SHRotateTapers(tapersrot,tapers,taper_order,lmax,nrot,x,dj,exitstatus=exitstatus)
    end subroutine cbindSHRotateTapers

    subroutine cbindSlepianCoeffs(falpha,falpha_d0,galpha,galpha_d0,galpha_d1,flm&
                                        ,flm_d0,flm_d1,flm_d2,lmax,nmax,exitstatus)  bind(c&
                                        , name="cbind_slepian_coeffs")
        use, intrinsic :: iso_c_binding
        use shtools, only: SlepianCoeffs
        implicit none
        integer(kind=c_int), intent(in) :: falpha_d0
        real(kind=c_double), dimension(falpha_d0),intent(out) :: falpha
        integer(kind=c_int), intent(in) :: galpha_d0
        integer(kind=c_int), intent(in) :: galpha_d1
        real(kind=c_double), dimension(galpha_d0,galpha_d1),intent(in) :: galpha
        integer(kind=c_int), intent(in) :: flm_d0
        integer(kind=c_int), intent(in) :: flm_d1
        integer(kind=c_int), intent(in) :: flm_d2
        real(kind=c_double), dimension(flm_d0,flm_d1,flm_d2),intent(in) :: flm
        integer(kind=c_int), intent(in) :: lmax
        integer(kind=c_int), intent(in) :: nmax
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call SlepianCoeffs(falpha,galpha,flm,lmax,nmax,exitstatus=exitstatus)
    end subroutine cbindSlepianCoeffs

    subroutine cbindSlepianCoeffsToSH(flm,flm_d0,flm_d1,flm_d2,falpha,falpha_d0,galpha&
                                         ,galpha_d0,galpha_d1,lmax,nmax,exitstatus)  bind(c&
                                         , name="cbind_slepian_coeffs_to_sh")
        use, intrinsic :: iso_c_binding
        use shtools, only: SlepianCoeffsToSH
        implicit none
        integer(kind=c_int), intent(in) :: flm_d0
        integer(kind=c_int), intent(in) :: flm_d1
        integer(kind=c_int), intent(in) :: flm_d2
        real(kind=c_double), dimension(flm_d0,flm_d1,flm_d2),intent(out) :: flm
        integer(kind=c_int), intent(in) :: falpha_d0
        real(kind=c_double), dimension(falpha_d0),intent(in) :: falpha
        integer(kind=c_int), intent(in) :: galpha_d0
        integer(kind=c_int), intent(in) :: galpha_d1
        real(kind=c_double), dimension(galpha_d0,galpha_d1),intent(in) :: galpha
        integer(kind=c_int), intent(in) :: lmax
        integer(kind=c_int), intent(in) :: nmax
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call SlepianCoeffsToSH(flm,falpha,galpha,lmax,nmax,exitstatus=exitstatus)
    end subroutine cbindSlepianCoeffsToSH

    subroutine cbindSHSCouplingMatrix(kij,kij_d0,kij_d1,galpha,galpha_d0,galpha_d1&
                                         ,lmax,nmax,exitstatus)  bind(c, name="cbind_shs_coupling_matrix")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHSCouplingMatrix
        implicit none
        integer(kind=c_int), intent(in) :: kij_d0
        integer(kind=c_int), intent(in) :: kij_d1
        real(kind=c_double), dimension(kij_d0,kij_d1),intent(out) :: kij
        integer(kind=c_int), intent(in) :: galpha_d0
        integer(kind=c_int), intent(in) :: galpha_d1
        real(kind=c_double), dimension(galpha_d0,galpha_d1),intent(in) :: galpha
        integer(kind=c_int), intent(in) :: lmax
        integer(kind=c_int), intent(in) :: nmax
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call SHSCouplingMatrix(kij,galpha,lmax,nmax,exitstatus=exitstatus)
    end subroutine cbindSHSCouplingMatrix

    subroutine cbindSHSlepianVar(l,galpha,galpha_d0,galpha_d1,galpha_order,galpha_order_d0&
                                  ,lmax,kmax,Sff,Sff_d0,variance,exitstatus)  bind(c&
                                  , name="cbind_sh_slepian_var")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHSlepianVar
        implicit none
        integer(kind=c_int), intent(in) :: galpha_d0
        integer(kind=c_int), intent(in) :: galpha_d1
        real(kind=c_double), dimension(galpha_d0,galpha_d1),intent(in) :: galpha
        integer(kind=c_int), intent(in) :: Sff_d0
        real(kind=c_double), dimension(Sff_d0),intent(in) :: Sff
        real(kind=c_double), intent(out) :: variance
        integer(kind=c_int), intent(in) :: l
        integer(kind=c_int), intent(in) :: lmax
        integer(kind=c_int), intent(in) :: kmax
        integer(kind=c_int), intent(in) :: galpha_order_d0
        integer(kind=c_int), dimension(galpha_order_d0),intent(in) :: galpha_order
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call SHSlepianVar(l,galpha,galpha_order,lmax,kmax,Sff,variance,exitstatus=exitstatus)
    end subroutine cbindSHSlepianVar

    subroutine cbindSHSCouplingMatrixCap(kij,kij_d0,kij_d1,galpha,galpha_d0,galpha_d1&
                                            ,galpha_order,galpha_order_d0,lmax,nmax&
                                            ,exitstatus)  bind(c, name="cbind_shs_coupling_matrix_cap")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHSCouplingMatrixCap
        implicit none
        integer(kind=c_int), intent(in) :: kij_d0
        integer(kind=c_int), intent(in) :: kij_d1
        real(kind=c_double), dimension(kij_d0,kij_d1),intent(out) :: kij
        integer(kind=c_int), intent(in) :: galpha_d0
        integer(kind=c_int), intent(in) :: galpha_d1
        real(kind=c_double), dimension(galpha_d0,galpha_d1),intent(in) :: galpha
        integer(kind=c_int), intent(in) :: galpha_order_d0
        integer(kind=c_int), dimension(galpha_order_d0),intent(in) :: galpha_order
        integer(kind=c_int), intent(in) :: lmax
        integer(kind=c_int), intent(in) :: nmax
        integer(kind=c_int), optional,intent(out) :: exitstatus
        call SHSCouplingMatrixCap(kij,galpha,galpha_order,lmax,nmax,exitstatus=exitstatus)
    end subroutine cbindSHSCouplingMatrixCap

