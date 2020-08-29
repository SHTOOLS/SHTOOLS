    subroutine cbindPlmBar(p,lmax,z,csphase,cnorm,exitstatus,p_d0)  bind(c, name="cbind_plm_bar")
        use, intrinsic :: iso_c_binding
        use shtools, only: PlmBar
        implicit none
        integer(kind=c_int), intent(in) :: p_d0
        integer(kind=c_int), intent(in) :: lmax
        real(kind=c_double), dimension(p_d0),intent(out) :: p
        real(kind=c_double), intent(in) :: z
        integer(kind=c_int), intent(in) :: csphase
        integer(kind=c_int), intent(in) :: cnorm
        integer(kind=c_int), intent(out) :: exitstatus
        call PlmBar(p,lmax,z,csphase,cnorm,exitstatus)
    end subroutine cbindPlmBar

    subroutine cbindPlmBar_d1(p,dp1,lmax,z,csphase,cnorm,exitstatus,p_d0,dp1_d0)  bind(c&
                               , name="cbind_plm_bar_d1")
        use, intrinsic :: iso_c_binding
        use shtools, only: PlmBar_d1
        implicit none
        integer(kind=c_int), intent(in) :: dp1_d0
        integer(kind=c_int), intent(in) :: p_d0
        integer(kind=c_int), intent(in) :: lmax
        real(kind=c_double), dimension(p_d0),intent(out) :: p
        real(kind=c_double), dimension(dp1_d0),intent(out) :: dp1
        real(kind=c_double), intent(in) :: z
        integer(kind=c_int), intent(in) :: csphase
        integer(kind=c_int), intent(in) :: cnorm
        integer(kind=c_int), intent(out) :: exitstatus
        call PlmBar_d1(p,dp1,lmax,z,csphase,cnorm,exitstatus)
    end subroutine cbindPlmBar_d1

    subroutine cbindPlBar(p,lmax,z,exitstatus,p_d0)  bind(c, name="cbind_pl_bar")
        use, intrinsic :: iso_c_binding
        use shtools, only: PlBar
        implicit none
        integer(kind=c_int), intent(in) :: p_d0
        integer(kind=c_int), intent(in) :: lmax
        real(kind=c_double), dimension(p_d0),intent(out) :: p
        real(kind=c_double), intent(in) :: z
        integer(kind=c_int), intent(out) :: exitstatus
        call PlBar(p,lmax,z,exitstatus)
    end subroutine cbindPlBar

    subroutine cbindPlBar_d1(p,dp1,lmax,z,exitstatus,p_d0,dp1_d0)  bind(c, name="cbind_pl_bar_d1")
        use, intrinsic :: iso_c_binding
        use shtools, only: PlBar_d1
        implicit none
        integer(kind=c_int), intent(in) :: dp1_d0
        integer(kind=c_int), intent(in) :: p_d0
        integer(kind=c_int), intent(in) :: lmax
        real(kind=c_double), dimension(p_d0),intent(out) :: p
        real(kind=c_double), dimension(dp1_d0),intent(out) :: dp1
        real(kind=c_double), intent(in) :: z
        integer(kind=c_int), intent(out) :: exitstatus
        call PlBar_d1(p,dp1,lmax,z,exitstatus)
    end subroutine cbindPlBar_d1

    subroutine cbindPlmON(p,lmax,z,csphase,cnorm,exitstatus,p_d0)  bind(c, name="cbind_plm_on")
        use, intrinsic :: iso_c_binding
        use shtools, only: PlmON
        implicit none
        integer(kind=c_int), intent(in) :: p_d0
        integer(kind=c_int), intent(in) :: lmax
        real(kind=c_double), dimension(p_d0),intent(out) :: p
        real(kind=c_double), intent(in) :: z
        integer(kind=c_int), intent(in) :: csphase
        integer(kind=c_int), intent(in) :: cnorm
        integer(kind=c_int), intent(out) :: exitstatus
        call PlmON(p,lmax,z,csphase,cnorm,exitstatus)
    end subroutine cbindPlmON

    subroutine cbindPlmON_d1(p,dp1,lmax,z,csphase,cnorm,exitstatus,p_d0,dp1_d0)  bind(c&
                              , name="cbind_plm_on_d1")
        use, intrinsic :: iso_c_binding
        use shtools, only: PlmON_d1
        implicit none
        integer(kind=c_int), intent(in) :: dp1_d0
        integer(kind=c_int), intent(in) :: p_d0
        integer(kind=c_int), intent(in) :: lmax
        real(kind=c_double), dimension(p_d0),intent(out) :: p
        real(kind=c_double), dimension(dp1_d0),intent(out) :: dp1
        real(kind=c_double), intent(in) :: z
        integer(kind=c_int), intent(in) :: csphase
        integer(kind=c_int), intent(in) :: cnorm
        integer(kind=c_int), intent(out) :: exitstatus
        call PlmON_d1(p,dp1,lmax,z,csphase,cnorm,exitstatus)
    end subroutine cbindPlmON_d1

    subroutine cbindPlON(p,lmax,z,exitstatus,p_d0)  bind(c, name="cbind_pl_on")
        use, intrinsic :: iso_c_binding
        use shtools, only: PlON
        implicit none
        integer(kind=c_int), intent(in) :: p_d0
        integer(kind=c_int), intent(in) :: lmax
        real(kind=c_double), dimension(p_d0),intent(out) :: p
        real(kind=c_double), intent(in) :: z
        integer(kind=c_int), intent(out) :: exitstatus
        call PlON(p,lmax,z,exitstatus)
    end subroutine cbindPlON

    subroutine cbindPlON_d1(p,dp1,lmax,z,exitstatus,p_d0,dp1_d0)  bind(c, name="cbind_pl_on_d1")
        use, intrinsic :: iso_c_binding
        use shtools, only: PlON_d1
        implicit none
        integer(kind=c_int), intent(in) :: dp1_d0
        integer(kind=c_int), intent(in) :: p_d0
        integer(kind=c_int), intent(in) :: lmax
        real(kind=c_double), dimension(p_d0),intent(out) :: p
        real(kind=c_double), dimension(dp1_d0),intent(out) :: dp1
        real(kind=c_double), intent(in) :: z
        integer(kind=c_int), intent(out) :: exitstatus
        call PlON_d1(p,dp1,lmax,z,exitstatus)
    end subroutine cbindPlON_d1

    subroutine cbindPlmSchmidt(p,lmax,z,csphase,cnorm,exitstatus,p_d0)  bind(c, name="cbind_plm_schmidt")
        use, intrinsic :: iso_c_binding
        use shtools, only: PlmSchmidt
        implicit none
        integer(kind=c_int), intent(in) :: p_d0
        integer(kind=c_int), intent(in) :: lmax
        real(kind=c_double), dimension(p_d0),intent(out) :: p
        real(kind=c_double), intent(in) :: z
        integer(kind=c_int), intent(in) :: csphase
        integer(kind=c_int), intent(in) :: cnorm
        integer(kind=c_int), intent(out) :: exitstatus
        call PlmSchmidt(p,lmax,z,csphase,cnorm,exitstatus)
    end subroutine cbindPlmSchmidt

    subroutine cbindPlmSchmidt_d1(p,dp1,lmax,z,csphase,cnorm,exitstatus,p_d0,dp1_d0)  bind(c&
                                   , name="cbind_plm_schmidt_d1")
        use, intrinsic :: iso_c_binding
        use shtools, only: PlmSchmidt_d1
        implicit none
        integer(kind=c_int), intent(in) :: dp1_d0
        integer(kind=c_int), intent(in) :: p_d0
        integer(kind=c_int), intent(in) :: lmax
        real(kind=c_double), dimension(p_d0),intent(out) :: p
        real(kind=c_double), dimension(dp1_d0),intent(out) :: dp1
        real(kind=c_double), intent(in) :: z
        integer(kind=c_int), intent(in) :: csphase
        integer(kind=c_int), intent(in) :: cnorm
        integer(kind=c_int), intent(out) :: exitstatus
        call PlmSchmidt_d1(p,dp1,lmax,z,csphase,cnorm,exitstatus)
    end subroutine cbindPlmSchmidt_d1

    subroutine cbindPlSchmidt(p,lmax,z,exitstatus,p_d0)  bind(c, name="cbind_pl_schmidt")
        use, intrinsic :: iso_c_binding
        use shtools, only: PlSchmidt
        implicit none
        integer(kind=c_int), intent(in) :: p_d0
        integer(kind=c_int), intent(in) :: lmax
        real(kind=c_double), dimension(p_d0),intent(out) :: p
        real(kind=c_double), intent(in) :: z
        integer(kind=c_int), intent(out) :: exitstatus
        call PlSchmidt(p,lmax,z,exitstatus)
    end subroutine cbindPlSchmidt

    subroutine cbindPlSchmidt_d1(p,dp1,lmax,z,exitstatus,p_d0,dp1_d0)  bind(c, name="cbind_pl_schmidt_d1")
        use, intrinsic :: iso_c_binding
        use shtools, only: PlSchmidt_d1
        implicit none
        integer(kind=c_int), intent(in) :: dp1_d0
        integer(kind=c_int), intent(in) :: p_d0
        integer(kind=c_int), intent(in) :: lmax
        real(kind=c_double), dimension(p_d0),intent(out) :: p
        real(kind=c_double), dimension(dp1_d0),intent(out) :: dp1
        real(kind=c_double), intent(in) :: z
        integer(kind=c_int), intent(out) :: exitstatus
        call PlSchmidt_d1(p,dp1,lmax,z,exitstatus)
    end subroutine cbindPlSchmidt_d1

    subroutine cbindPLegendreA(p,lmax,z,csphase,exitstatus,p_d0)  bind(c, name="cbind_p_legendre_a")
        use, intrinsic :: iso_c_binding
        use shtools, only: PLegendreA
        implicit none
        integer(kind=c_int), intent(in) :: p_d0
        integer(kind=c_int), intent(in) :: lmax
        real(kind=c_double), dimension(p_d0),intent(out) :: p
        real(kind=c_double), intent(in) :: z
        integer(kind=c_int), intent(in) :: csphase
        integer(kind=c_int), intent(out) :: exitstatus
        call PLegendreA(p,lmax,z,csphase,exitstatus)
    end subroutine cbindPLegendreA

    subroutine cbindPLegendreA_d1(p,dp1,lmax,z,csphase,exitstatus,p_d0,dp1_d0)  bind(c&
                                   , name="cbind_p_legendre_a_d1")
        use, intrinsic :: iso_c_binding
        use shtools, only: PLegendreA_d1
        implicit none
        integer(kind=c_int), intent(in) :: dp1_d0
        integer(kind=c_int), intent(in) :: p_d0
        integer(kind=c_int), intent(in) :: lmax
        real(kind=c_double), dimension(p_d0),intent(out) :: p
        real(kind=c_double), dimension(dp1_d0),intent(out) :: dp1
        real(kind=c_double), intent(in) :: z
        integer(kind=c_int), intent(in) :: csphase
        integer(kind=c_int), intent(out) :: exitstatus
        call PLegendreA_d1(p,dp1,lmax,z,csphase,exitstatus)
    end subroutine cbindPLegendreA_d1

    subroutine cbindPLegendre(p,lmax,z,exitstatus,p_d0)  bind(c, name="cbind_p_legendre")
        use, intrinsic :: iso_c_binding
        use shtools, only: PLegendre
        implicit none
        integer(kind=c_int), intent(in) :: p_d0
        integer(kind=c_int), intent(in) :: lmax
        real(kind=c_double), dimension(p_d0),intent(out) :: p
        real(kind=c_double), intent(in) :: z
        integer(kind=c_int), intent(out) :: exitstatus
        call PLegendre(p,lmax,z,exitstatus)
    end subroutine cbindPLegendre

    subroutine cbindPLegendre_d1(p,dp1,lmax,z,exitstatus,p_d0,dp1_d0)  bind(c, name="cbind_p_legendre_d1")
        use, intrinsic :: iso_c_binding
        use shtools, only: PLegendre_d1
        implicit none
        integer(kind=c_int), intent(in) :: dp1_d0
        integer(kind=c_int), intent(in) :: p_d0
        integer(kind=c_int), intent(in) :: lmax
        real(kind=c_double), dimension(p_d0),intent(out) :: p
        real(kind=c_double), dimension(dp1_d0),intent(out) :: dp1
        real(kind=c_double), intent(in) :: z
        integer(kind=c_int), intent(out) :: exitstatus
        call PLegendre_d1(p,dp1,lmax,z,exitstatus)
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

    subroutine cbindSHExpandDH(grid,n,cilm,lmax,norm,sampling,csphase,lmax_calc,exitstatus&
                                   ,grid_d0,grid_d1,cilm_d)  bind(c, name="cbind_sh_expand_dh")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHExpandDH
        implicit none
        integer(kind=c_int), intent(in) :: cilm_d
        integer(kind=c_int), intent(in) :: grid_d0
        integer(kind=c_int), intent(in) :: grid_d1
        real(kind=c_double), dimension(grid_d0,grid_d1),intent(in) :: grid
        real(kind=c_double), dimension(2,cilm_d,cilm_d),intent(out) :: cilm
        integer(kind=c_int), intent(in) :: n
        integer(kind=c_int), intent(out) :: lmax
        integer(kind=c_int), intent(in) :: norm
        integer(kind=c_int), intent(in) :: sampling
        integer(kind=c_int), intent(in) :: csphase
        integer(kind=c_int), intent(in) :: lmax_calc
        integer(kind=c_int), intent(out) :: exitstatus
        call SHExpandDH(grid,n,cilm,lmax,norm,sampling,csphase,lmax_calc,exitstatus)
    end subroutine cbindSHExpandDH

    subroutine cbindMakeGridDH(griddh,n,cilm,lmax,norm,sampling,csphase,lmax_calc&
                                     ,extend,exitstatus,cilm_d,griddh_d0,griddh_d1)  bind(c&
                                     , name="cbind_make_grid_dh")
        use, intrinsic :: iso_c_binding
        use shtools, only: MakeGridDH
        implicit none
        integer(kind=c_int), intent(in) :: griddh_d0
        integer(kind=c_int), intent(in) :: griddh_d1
        integer(kind=c_int), intent(in) :: cilm_d
        real(kind=c_double), dimension(2,cilm_d,cilm_d),intent(in) :: cilm
        real(kind=c_double), dimension(griddh_d0,griddh_d1),intent(out) :: griddh
        integer(kind=c_int), intent(in) :: lmax
        integer(kind=c_int), intent(out) :: n
        integer(kind=c_int), intent(in) :: norm
        integer(kind=c_int), intent(in) :: sampling
        integer(kind=c_int), intent(in) :: csphase
        integer(kind=c_int), intent(in) :: lmax_calc
        integer(kind=c_int), intent(in) :: extend
        integer(kind=c_int), intent(out) :: exitstatus
        call MakeGridDH(griddh,n,cilm,lmax,norm,sampling,csphase,lmax_calc,extend&
                              ,exitstatus)
    end subroutine cbindMakeGridDH

    subroutine cbindSHExpandDHC(grid,n,cilm,lmax,norm,sampling,csphase,lmax_calc,exitstatus&
                                    ,grid_d0,grid_d1,cilm_d)  bind(c, name="cbind_sh_expand_dhc")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHExpandDHC
        implicit none
        integer(kind=c_int), intent(in) :: cilm_d
        integer(kind=c_int), intent(in) :: grid_d0
        integer(kind=c_int), intent(in) :: grid_d1
        complex(kind=c_double_complex), dimension(grid_d0,grid_d1),intent(in) :: grid
        complex(kind=c_double_complex), dimension(2,cilm_d,cilm_d),intent(out) :: cilm
        integer(kind=c_int), intent(in) :: n
        integer(kind=c_int), intent(out) :: lmax
        integer(kind=c_int), intent(in) :: norm
        integer(kind=c_int), intent(in) :: sampling
        integer(kind=c_int), intent(in) :: csphase
        integer(kind=c_int), intent(in) :: lmax_calc
        integer(kind=c_int), intent(out) :: exitstatus
        call SHExpandDHC(grid,n,cilm,lmax,norm,sampling,csphase,lmax_calc,exitstatus)
    end subroutine cbindSHExpandDHC

    subroutine cbindMakeGridDHC(griddh,n,cilm,lmax,norm,sampling,csphase,lmax_calc&
                                      ,extend,exitstatus,cilm_d,griddh_d0,griddh_d1)  bind(c&
                                      , name="cbind_make_grid_dhc")
        use, intrinsic :: iso_c_binding
        use shtools, only: MakeGridDHC
        implicit none
        integer(kind=c_int), intent(in) :: griddh_d0
        integer(kind=c_int), intent(in) :: griddh_d1
        integer(kind=c_int), intent(in) :: cilm_d
        complex(kind=c_double_complex), dimension(2,cilm_d,cilm_d),intent(in) :: cilm
        complex(kind=c_double_complex), dimension(griddh_d0,griddh_d1),intent(out) :: griddh
        integer(kind=c_int), intent(in) :: lmax
        integer(kind=c_int), intent(out) :: n
        integer(kind=c_int), intent(in) :: norm
        integer(kind=c_int), intent(in) :: sampling
        integer(kind=c_int), intent(in) :: csphase
        integer(kind=c_int), intent(in) :: lmax_calc
        integer(kind=c_int), intent(in) :: extend
        integer(kind=c_int), intent(out) :: exitstatus
        call MakeGridDHC(griddh,n,cilm,lmax,norm,sampling,csphase,lmax_calc,extend&
                               ,exitstatus)
    end subroutine cbindMakeGridDHC

    subroutine cbindSHGLQ(lmax,zero,w,plx,norm,csphase,cnorm,exitstatus,zero_d0,w_d0&
                              ,plx_d0,plx_d1)  bind(c, name="cbind_shglq")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHGLQ
        implicit none
        integer(kind=c_int), intent(in) :: plx_d0
        integer(kind=c_int), intent(in) :: plx_d1
        integer(kind=c_int), intent(in) :: w_d0
        integer(kind=c_int), intent(in) :: zero_d0
        integer(kind=c_int), intent(in) :: lmax
        real(kind=c_double), dimension(zero_d0),intent(out) :: zero
        real(kind=c_double), dimension(w_d0),intent(out) :: w
        real(kind=c_double), dimension(plx_d0,plx_d1),intent(out) :: plx
        integer(kind=c_int), intent(in) :: norm
        integer(kind=c_int), intent(in) :: csphase
        integer(kind=c_int), intent(in) :: cnorm
        integer(kind=c_int), intent(out) :: exitstatus
        call SHGLQ(lmax,zero,w,plx,norm,csphase,cnorm,exitstatus)
    end subroutine cbindSHGLQ

    subroutine cbindSHExpandGLQ(cilm,lmax,gridglq,w,plx,zero,norm,csphase,lmax_calc&
                                    ,exitstatus,w_d0,gridglq_d0,gridglq_d1,plx_d0&
                                    ,plx_d1,zero_d0,cilm_d)  bind(c, name="cbind_sh_expand_glq")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHExpandGLQ
        implicit none
        integer(kind=c_int), intent(in) :: cilm_d
        integer(kind=c_int), intent(in) :: zero_d0
        integer(kind=c_int), intent(in) :: plx_d0
        integer(kind=c_int), intent(in) :: plx_d1
        integer(kind=c_int), intent(in) :: gridglq_d0
        integer(kind=c_int), intent(in) :: gridglq_d1
        integer(kind=c_int), intent(in) :: w_d0
        real(kind=c_double), dimension(w_d0),intent(in) :: w
        real(kind=c_double), dimension(gridglq_d0,gridglq_d1),intent(in) :: gridglq
        real(kind=c_double), dimension(plx_d0,plx_d1),intent(in) :: plx
        real(kind=c_double), dimension(zero_d0),intent(in) :: zero
        real(kind=c_double), dimension(2,cilm_d,cilm_d),intent(out) :: cilm
        integer(kind=c_int), intent(in) :: lmax
        integer(kind=c_int), intent(in) :: norm
        integer(kind=c_int), intent(in) :: csphase
        integer(kind=c_int), intent(in) :: lmax_calc
        integer(kind=c_int), intent(out) :: exitstatus
        call SHExpandGLQ(cilm,lmax,gridglq,w,plx,zero,norm,csphase,lmax_calc,exitstatus)
    end subroutine cbindSHExpandGLQ

    subroutine cbindMakeGridGLQ(gridglq,cilm,lmax,plx,zero,norm,csphase,lmax_calc&
                                       ,extend,exitstatus,cilm_d,plx_d0,plx_d1,zero_d0&
                                       ,gridglq_d0,gridglq_d1)  bind(c, name="cbind_make_grid_glq")
        use, intrinsic :: iso_c_binding
        use shtools, only: MakeGridGLQ
        implicit none
        integer(kind=c_int), intent(in) :: gridglq_d0
        integer(kind=c_int), intent(in) :: gridglq_d1
        integer(kind=c_int), intent(in) :: zero_d0
        integer(kind=c_int), intent(in) :: plx_d0
        integer(kind=c_int), intent(in) :: plx_d1
        integer(kind=c_int), intent(in) :: cilm_d
        real(kind=c_double), dimension(2,cilm_d,cilm_d),intent(in) :: cilm
        real(kind=c_double), dimension(plx_d0,plx_d1),intent(in) :: plx
        real(kind=c_double), dimension(zero_d0),intent(in) :: zero
        real(kind=c_double), dimension(gridglq_d0,gridglq_d1),intent(out) :: gridglq
        integer(kind=c_int), intent(in) :: lmax
        integer(kind=c_int), intent(in) :: norm
        integer(kind=c_int), intent(in) :: csphase
        integer(kind=c_int), intent(in) :: lmax_calc
        integer(kind=c_int), intent(in) :: extend
        integer(kind=c_int), intent(out) :: exitstatus
        call MakeGridGLQ(gridglq,cilm,lmax,plx,zero,norm,csphase,lmax_calc,extend&
                                ,exitstatus)
    end subroutine cbindMakeGridGLQ

    subroutine cbindSHExpandGLQC(cilm,lmax,gridglq,w,plx,zero,norm,csphase,lmax_calc&
                                     ,exitstatus,w_d0,gridglq_d0,gridglq_d1,plx_d0&
                                     ,plx_d1,zero_d0,cilm_d)  bind(c, name="cbind_sh_expand_glqc")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHExpandGLQC
        implicit none
        integer(kind=c_int), intent(in) :: cilm_d
        integer(kind=c_int), intent(in) :: zero_d0
        integer(kind=c_int), intent(in) :: plx_d0
        integer(kind=c_int), intent(in) :: plx_d1
        integer(kind=c_int), intent(in) :: gridglq_d0
        integer(kind=c_int), intent(in) :: gridglq_d1
        integer(kind=c_int), intent(in) :: w_d0
        real(kind=c_double), dimension(w_d0),intent(in) :: w
        complex(kind=c_double_complex), dimension(gridglq_d0,gridglq_d1),intent(in) :: gridglq
        real(kind=c_double), dimension(plx_d0,plx_d1),intent(in) :: plx
        real(kind=c_double), dimension(zero_d0),intent(in) :: zero
        complex(kind=c_double_complex), dimension(2,cilm_d,cilm_d),intent(out) :: cilm
        integer(kind=c_int), intent(in) :: lmax
        integer(kind=c_int), intent(in) :: norm
        integer(kind=c_int), intent(in) :: csphase
        integer(kind=c_int), intent(in) :: lmax_calc
        integer(kind=c_int), intent(out) :: exitstatus
        call SHExpandGLQC(cilm,lmax,gridglq,w,plx,zero,norm,csphase,lmax_calc,exitstatus)
    end subroutine cbindSHExpandGLQC

    subroutine cbindMakeGridGLQC(gridglq,cilm,lmax,plx,zero,norm,csphase,lmax_calc&
                                        ,extend,exitstatus,cilm_d,plx_d0,plx_d1,zero_d0&
                                        ,gridglq_d0,gridglq_d1)  bind(c, name="cbind_make_grid_glqc")
        use, intrinsic :: iso_c_binding
        use shtools, only: MakeGridGLQC
        implicit none
        integer(kind=c_int), intent(in) :: gridglq_d0
        integer(kind=c_int), intent(in) :: gridglq_d1
        integer(kind=c_int), intent(in) :: zero_d0
        integer(kind=c_int), intent(in) :: plx_d0
        integer(kind=c_int), intent(in) :: plx_d1
        integer(kind=c_int), intent(in) :: cilm_d
        complex(kind=c_double_complex), dimension(2,cilm_d,cilm_d),intent(in) :: cilm
        real(kind=c_double), dimension(plx_d0,plx_d1),intent(in) :: plx
        real(kind=c_double), dimension(zero_d0),intent(in) :: zero
        complex(kind=c_double_complex), dimension(gridglq_d0,gridglq_d1),intent(out) :: gridglq
        integer(kind=c_int), intent(in) :: lmax
        integer(kind=c_int), intent(in) :: norm
        integer(kind=c_int), intent(in) :: csphase
        integer(kind=c_int), intent(in) :: lmax_calc
        integer(kind=c_int), intent(in) :: extend
        integer(kind=c_int), intent(out) :: exitstatus
        call MakeGridGLQC(gridglq,cilm,lmax,plx,zero,norm,csphase,lmax_calc,extend&
                                 ,exitstatus)
    end subroutine cbindMakeGridGLQC

    subroutine cbindGLQGridCoord(latglq,longlq,lmax,nlat,nlong,extend,exitstatus,latglq_d0&
                                       ,longlq_d0)  bind(c, name="cbind_glq_grid_coord")
        use, intrinsic :: iso_c_binding
        use shtools, only: GLQGridCoord
        implicit none
        integer(kind=c_int), intent(in) :: longlq_d0
        integer(kind=c_int), intent(in) :: latglq_d0
        integer(kind=c_int), intent(in) :: lmax
        integer(kind=c_int), intent(out) :: nlat
        integer(kind=c_int), intent(out) :: nlong
        real(kind=c_double), dimension(latglq_d0),intent(out) :: latglq
        real(kind=c_double), dimension(longlq_d0),intent(out) :: longlq
        integer(kind=c_int), intent(in) :: extend
        integer(kind=c_int), intent(out) :: exitstatus
        call GLQGridCoord(latglq,longlq,lmax,nlat,nlong,extend,exitstatus)
    end subroutine cbindGLQGridCoord

    subroutine cbindSHExpandLSQ(cilm,d,lat,lon,nmax,lmax,norm,chi2,csphase,weights&
                                    ,exitstatus,d_d0,lat_d0,lon_d0,cilm_d,weights_d0)  bind(c&
                                    , name="cbind_sh_expand_lsq")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHExpandLSQ
        implicit none
        integer(kind=c_int), intent(in) :: weights_d0
        integer(kind=c_int), intent(in) :: cilm_d
        integer(kind=c_int), intent(in) :: lon_d0
        integer(kind=c_int), intent(in) :: lat_d0
        integer(kind=c_int), intent(in) :: d_d0
        real(kind=c_double), dimension(d_d0),intent(in) :: d
        real(kind=c_double), dimension(lat_d0),intent(in) :: lat
        real(kind=c_double), dimension(lon_d0),intent(in) :: lon
        real(kind=c_double), dimension(2,cilm_d,cilm_d),intent(out) :: cilm
        integer(kind=c_int), intent(in) :: nmax
        integer(kind=c_int), intent(in) :: lmax
        integer(kind=c_int), intent(in) :: norm
        integer(kind=c_int), intent(in) :: csphase
        real(kind=c_double), intent(out) :: chi2
        real(kind=c_double), dimension(weights_d0),intent(in) :: weights
        integer(kind=c_int), intent(out) :: exitstatus
        call SHExpandLSQ(cilm,d,lat,lon,nmax,lmax,norm,chi2,csphase,weights,exitstatus)
    end subroutine cbindSHExpandLSQ

    subroutine cbindMakeGrid2d(grid,cilm,lmax,interval,nlat,nlong,norm,csphase,f,a&
                                   ,north,south,east,west,dealloc,exitstatus,cilm_d&
                                   ,grid_d0,grid_d1)  bind(c, name="cbind_make_grid2d")
        use, intrinsic :: iso_c_binding
        use shtools, only: MakeGrid2d
        implicit none
        integer(kind=c_int), intent(in) :: grid_d0
        integer(kind=c_int), intent(in) :: grid_d1
        integer(kind=c_int), intent(in) :: cilm_d
        real(kind=c_double), dimension(2,cilm_d,cilm_d),intent(in) :: cilm
        real(kind=c_double), intent(in) :: interval
        real(kind=c_double), dimension(grid_d0,grid_d1),intent(out) :: grid
        integer(kind=c_int), intent(in) :: lmax
        integer(kind=c_int), intent(out) :: nlat
        integer(kind=c_int), intent(out) :: nlong
        integer(kind=c_int), intent(in) :: norm
        integer(kind=c_int), intent(in) :: csphase
        integer(kind=c_int), intent(in) :: dealloc
        real(kind=c_double), intent(in) :: f
        real(kind=c_double), intent(in) :: a
        real(kind=c_double), intent(in) :: north
        real(kind=c_double), intent(in) :: south
        real(kind=c_double), intent(in) :: east
        real(kind=c_double), intent(in) :: west
        integer(kind=c_int), intent(out) :: exitstatus
        call MakeGrid2d(grid,cilm,lmax,interval,nlat,nlong,norm,csphase,f,a,north&
                            ,south,east,west,dealloc,exitstatus)
    end subroutine cbindMakeGrid2d

    function cbindMakeGridPoint(cilm,lmax,lat,lon,norm,csphase,dealloc,cilm_d)  bind(c&
                                    , name="cbind_make_grid_point")
        use, intrinsic :: iso_c_binding
        use shtools, only: MakeGridPoint
        implicit none
        integer(kind=c_int), intent(in) :: cilm_d
        real(kind=c_double) :: cbindMakeGridPoint
        real(kind=c_double), dimension(2,cilm_d,cilm_d),intent(in) :: cilm
        real(kind=c_double), intent(in) :: lat
        real(kind=c_double), intent(in) :: lon
        integer(kind=c_int), intent(in) :: lmax
        integer(kind=c_int), intent(in) :: norm
        integer(kind=c_int), intent(in) :: csphase
        integer(kind=c_int), intent(in) :: dealloc
        cbindMakeGridPoint=MakeGridPoint(cilm,lmax,lat,lon,norm,csphase,dealloc)
    end function cbindMakeGridPoint

    function cbindMakeGridPointC(cilm,lmax,lat,lon,norm,csphase,dealloc,cilm_d)  bind(c&
                                     , name="cbind_make_grid_point_c")
        use, intrinsic :: iso_c_binding
        use shtools, only: MakeGridPointC
        implicit none
        integer(kind=c_int), intent(in) :: cilm_d
        complex(kind=c_double_complex) :: cbindMakeGridPointC
        complex(kind=c_double_complex), dimension(2,cilm_d,cilm_d),intent(in) :: cilm
        real(kind=c_double), intent(in) :: lat
        real(kind=c_double), intent(in) :: lon
        integer(kind=c_int), intent(in) :: lmax
        integer(kind=c_int), intent(in) :: norm
        integer(kind=c_int), intent(in) :: csphase
        integer(kind=c_int), intent(in) :: dealloc
        cbindMakeGridPointC=MakeGridPointC(cilm,lmax,lat,lon,norm,csphase,dealloc)
    end function cbindMakeGridPointC

    subroutine cbindSHMultiply(shout,sh1,lmax1,sh2,lmax2,precomp,norm,csphase,exitstatus&
                                    ,shout_d0,shout_d1,shout_d2,sh1_d0,sh1_d1,sh1_d2&
                                    ,sh2_d0,sh2_d1,sh2_d2)  bind(c, name="cbind_sh_multiply")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHMultiply
        implicit none
        integer(kind=c_int), intent(in) :: sh2_d0
        integer(kind=c_int), intent(in) :: sh2_d1
        integer(kind=c_int), intent(in) :: sh2_d2
        integer(kind=c_int), intent(in) :: sh1_d0
        integer(kind=c_int), intent(in) :: sh1_d1
        integer(kind=c_int), intent(in) :: sh1_d2
        integer(kind=c_int), intent(in) :: shout_d0
        integer(kind=c_int), intent(in) :: shout_d1
        integer(kind=c_int), intent(in) :: shout_d2
        real(kind=c_double), dimension(shout_d0,shout_d1,shout_d2),intent(out) :: shout
        real(kind=c_double), dimension(sh1_d0,sh1_d1,sh1_d2),intent(in) :: sh1
        real(kind=c_double), dimension(sh2_d0,sh2_d1,sh2_d2),intent(in) :: sh2
        integer(kind=c_int), intent(in) :: lmax1
        integer(kind=c_int), intent(in) :: lmax2
        integer(kind=c_int), intent(in) :: precomp
        integer(kind=c_int), intent(in) :: norm
        integer(kind=c_int), intent(in) :: csphase
        integer(kind=c_int), intent(out) :: exitstatus
        call SHMultiply(shout,sh1,lmax1,sh2,lmax2,precomp,norm,csphase,exitstatus)
    end subroutine cbindSHMultiply

    subroutine cbindSHRead(filename,cilm,lmax,skip,header,error,exitstatus,cilm_d&
                                   ,header_d0,error_d0,error_d1,error_d2)  bind(c&
                                   , name="cbind_sh_read")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHRead
        implicit none
        integer(kind=c_int), intent(in) :: error_d0
        integer(kind=c_int), intent(in) :: error_d1
        integer(kind=c_int), intent(in) :: error_d2
        integer(kind=c_int), intent(in) :: header_d0
        integer(kind=c_int), intent(in) :: cilm_d
        character(kind=c_char), intent(in) :: filename
        integer(kind=c_int), intent(out) :: lmax
        real(kind=c_double), dimension(2,cilm_d,cilm_d),intent(out) :: cilm
        real(kind=c_double), dimension(header_d0),intent(out) :: header
        real(kind=c_double), dimension(error_d0,error_d1,error_d2),intent(out) :: error
        integer(kind=c_int), intent(in) :: skip
        integer(kind=c_int), intent(out) :: exitstatus
        call SHRead(filename,cilm,lmax,skip,header,error,exitstatus)
    end subroutine cbindSHRead

    subroutine cbindSHRead2(filename,cilm,lmax,gm,r0_pot,error,dot,doystart,doyend&
                                    ,epoch,exitstatus,cilm_d,error_d0,error_d1,error_d2&
                                    ,dot_d0,dot_d1,dot_d2)  bind(c, name="cbind_sh_read2")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHRead2
        implicit none
        integer(kind=c_int), intent(in) :: dot_d0
        integer(kind=c_int), intent(in) :: dot_d1
        integer(kind=c_int), intent(in) :: dot_d2
        integer(kind=c_int), intent(in) :: error_d0
        integer(kind=c_int), intent(in) :: error_d1
        integer(kind=c_int), intent(in) :: error_d2
        integer(kind=c_int), intent(in) :: cilm_d
        character(kind=c_char), intent(in) :: filename
        integer(kind=c_int), intent(out) :: lmax
        real(kind=c_double), dimension(2,cilm_d,cilm_d),intent(out) :: cilm
        real(kind=c_double), intent(out) :: gm
        real(kind=c_double), intent(out) :: r0_pot
        real(kind=c_double), dimension(error_d0,error_d1,error_d2),intent(out) :: error
        real(kind=c_double), dimension(dot_d0,dot_d1,dot_d2),intent(out) :: dot
        real(kind=c_double), intent(out) :: doystart
        real(kind=c_double), intent(out) :: doyend
        real(kind=c_double), intent(out) :: epoch
        integer(kind=c_int), intent(out) :: exitstatus
        call SHRead2(filename,cilm,lmax,gm,r0_pot,error,dot,doystart,doyend,epoch&
                             ,exitstatus)
    end subroutine cbindSHRead2

    subroutine cbindSHReadJPL(filename,cilm,lmax,error,gm,formatstring,exitstatus&
                                      ,cilm_d,error_d0,error_d1,error_d2)  bind(c&
                                      , name="cbind_sh_read_jpl")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHReadJPL
        implicit none
        integer(kind=c_int), intent(in) :: error_d0
        integer(kind=c_int), intent(in) :: error_d1
        integer(kind=c_int), intent(in) :: error_d2
        integer(kind=c_int), intent(in) :: cilm_d
        character(kind=c_char), intent(in) :: filename
        integer(kind=c_int), intent(in) :: lmax
        real(kind=c_double), dimension(2,cilm_d,cilm_d),intent(out) :: cilm
        real(kind=c_double), dimension(error_d0,error_d1,error_d2),intent(out) :: error
        real(kind=c_double), dimension(2),intent(out) :: gm
        character(kind=c_char), intent(in) :: formatstring
        integer(kind=c_int), intent(out) :: exitstatus
        call SHReadJPL(filename,cilm,lmax,error,gm,formatstring,exitstatus)
    end subroutine cbindSHReadJPL

    subroutine cbindSHCilmToVector(cilm,vector,lmax,exitstatus,cilm_d,vector_d0)  bind(c&
                                       , name="cbind_sh_cilm_to_vector")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHCilmToVector
        implicit none
        integer(kind=c_int), intent(in) :: vector_d0
        integer(kind=c_int), intent(in) :: cilm_d
        real(kind=c_double), dimension(2,cilm_d,cilm_d),intent(in) :: cilm
        real(kind=c_double), dimension(vector_d0),intent(out) :: vector
        integer(kind=c_int), intent(in) :: lmax
        integer(kind=c_int), intent(out) :: exitstatus
        call SHCilmToVector(cilm,vector,lmax,exitstatus)
    end subroutine cbindSHCilmToVector

    subroutine cbindSHVectorToCilm(vector,cilm,lmax,exitstatus,cilm_d,vector_d0)  bind(c&
                                         , name="cbind_sh_vector_to_cilm")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHVectorToCilm
        implicit none
        integer(kind=c_int), intent(in) :: vector_d0
        integer(kind=c_int), intent(in) :: cilm_d
        real(kind=c_double), dimension(2,cilm_d,cilm_d),intent(out) :: cilm
        real(kind=c_double), dimension(vector_d0),intent(in) :: vector
        integer(kind=c_int), intent(in) :: lmax
        integer(kind=c_int), intent(out) :: exitstatus
        call SHVectorToCilm(vector,cilm,lmax,exitstatus)
    end subroutine cbindSHVectorToCilm

    subroutine cbindSHCilmToCindex(cilm,cindex,degmax,exitstatus,cilm_d,cindex_d0&
                                       ,cindex_d1)  bind(c, name="cbind_sh_cilm_to_cindex")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHCilmToCindex
        implicit none
        integer(kind=c_int), intent(in) :: cindex_d0
        integer(kind=c_int), intent(in) :: cindex_d1
        integer(kind=c_int), intent(in) :: cilm_d
        real(kind=c_double), dimension(2,cilm_d,cilm_d),intent(in) :: cilm
        real(kind=c_double), dimension(cindex_d0,cindex_d1),intent(out) :: cindex
        integer(kind=c_int), intent(in) :: degmax
        integer(kind=c_int), intent(out) :: exitstatus
        call SHCilmToCindex(cilm,cindex,degmax,exitstatus)
    end subroutine cbindSHCilmToCindex

    subroutine cbindSHCindexToCilm(cindex,cilm,degmax,exitstatus,cilm_d,cindex_d0&
                                         ,cindex_d1)  bind(c, name="cbind_sh_cindex_to_cilm")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHCindexToCilm
        implicit none
        integer(kind=c_int), intent(in) :: cindex_d0
        integer(kind=c_int), intent(in) :: cindex_d1
        integer(kind=c_int), intent(in) :: cilm_d
        real(kind=c_double), dimension(2,cilm_d,cilm_d),intent(out) :: cilm
        real(kind=c_double), dimension(cindex_d0,cindex_d1),intent(in) :: cindex
        integer(kind=c_int), intent(in) :: degmax
        integer(kind=c_int), intent(out) :: exitstatus
        call SHCindexToCilm(cindex,cilm,degmax,exitstatus)
    end subroutine cbindSHCindexToCilm

    subroutine cbindSHrtoc(rcilm,ccilm,degmax,convention,switchcs,exitstatus,rcilm_d0&
                                ,rcilm_d1,rcilm_d2,ccilm_d0,ccilm_d1,ccilm_d2)  bind(c&
                                , name="cbind_s_hrtoc")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHrtoc
        implicit none
        integer(kind=c_int), intent(in) :: ccilm_d0
        integer(kind=c_int), intent(in) :: ccilm_d1
        integer(kind=c_int), intent(in) :: ccilm_d2
        integer(kind=c_int), intent(in) :: rcilm_d0
        integer(kind=c_int), intent(in) :: rcilm_d1
        integer(kind=c_int), intent(in) :: rcilm_d2
        real(kind=c_double), dimension(rcilm_d0,rcilm_d1,rcilm_d2),intent(in) :: rcilm
        real(kind=c_double), dimension(ccilm_d0,ccilm_d1,ccilm_d2),intent(out) :: ccilm
        integer(kind=c_int), intent(in) :: degmax
        integer(kind=c_int), intent(in) :: convention
        integer(kind=c_int), intent(in) :: switchcs
        integer(kind=c_int), intent(out) :: exitstatus
        call SHrtoc(rcilm,ccilm,degmax,convention,switchcs,exitstatus)
    end subroutine cbindSHrtoc

    subroutine cbindSHctor(ccilm,rcilm,degmax,convention,switchcs,exitstatus,ccilm_d0&
                                ,ccilm_d1,ccilm_d2,rcilm_d0,rcilm_d1,rcilm_d2)  bind(c&
                                , name="cbind_s_hctor")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHctor
        implicit none
        integer(kind=c_int), intent(in) :: rcilm_d0
        integer(kind=c_int), intent(in) :: rcilm_d1
        integer(kind=c_int), intent(in) :: rcilm_d2
        integer(kind=c_int), intent(in) :: ccilm_d0
        integer(kind=c_int), intent(in) :: ccilm_d1
        integer(kind=c_int), intent(in) :: ccilm_d2
        real(kind=c_double), dimension(ccilm_d0,ccilm_d1,ccilm_d2),intent(in) :: ccilm
        real(kind=c_double), dimension(rcilm_d0,rcilm_d1,rcilm_d2),intent(out) :: rcilm
        integer(kind=c_int), intent(in) :: degmax
        integer(kind=c_int), intent(in) :: convention
        integer(kind=c_int), intent(in) :: switchcs
        integer(kind=c_int), intent(out) :: exitstatus
        call SHctor(ccilm,rcilm,degmax,convention,switchcs,exitstatus)
    end subroutine cbindSHctor

    subroutine cbinddjpi2(dj,lmax,exitstatus,dj_d0,dj_d1,dj_d2)  bind(c, name="cbinddjpi2")
        use, intrinsic :: iso_c_binding
        use shtools, only: djpi2
        implicit none
        integer(kind=c_int), intent(in) :: dj_d0
        integer(kind=c_int), intent(in) :: dj_d1
        integer(kind=c_int), intent(in) :: dj_d2
        integer(kind=c_int), intent(in) :: lmax
        real(kind=c_double), dimension(dj_d0,dj_d1,dj_d2),intent(out) :: dj
        integer(kind=c_int), intent(out) :: exitstatus
        call djpi2(dj,lmax,exitstatus)
    end subroutine cbinddjpi2

    subroutine cbindSHRotateCoef(x,cof,rcof,dj,lmax,exitstatus,cof_d0,cof_d1,dj_d0&
                                  ,dj_d1,dj_d2,rcof_d0,rcof_d1)  bind(c, name="cbind_sh_rotate_coef")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHRotateCoef
        implicit none
        integer(kind=c_int), intent(in) :: rcof_d0
        integer(kind=c_int), intent(in) :: rcof_d1
        integer(kind=c_int), intent(in) :: dj_d0
        integer(kind=c_int), intent(in) :: dj_d1
        integer(kind=c_int), intent(in) :: dj_d2
        integer(kind=c_int), intent(in) :: cof_d0
        integer(kind=c_int), intent(in) :: cof_d1
        real(kind=c_double), dimension(cof_d0,cof_d1),intent(in) :: cof
        real(kind=c_double), dimension(dj_d0,dj_d1,dj_d2),intent(in) :: dj
        real(kind=c_double), dimension(3),intent(in) :: x
        real(kind=c_double), dimension(rcof_d0,rcof_d1),intent(out) :: rcof
        integer(kind=c_int), intent(in) :: lmax
        integer(kind=c_int), intent(out) :: exitstatus
        call SHRotateCoef(x,cof,rcof,dj,lmax,exitstatus)
    end subroutine cbindSHRotateCoef

    subroutine cbindSHRotateRealCoef(cilmrot,cilm,lmax,x,dj,exitstatus,cilm_d,x_d0&
                                            ,dj_d0,dj_d1,dj_d2,cilmrot_d0,cilmrot_d1&
                                            ,cilmrot_d2)  bind(c, name="cbind_sh_rotate_real_coef")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHRotateRealCoef
        implicit none
        integer(kind=c_int), intent(in) :: cilmrot_d0
        integer(kind=c_int), intent(in) :: cilmrot_d1
        integer(kind=c_int), intent(in) :: cilmrot_d2
        integer(kind=c_int), intent(in) :: dj_d0
        integer(kind=c_int), intent(in) :: dj_d1
        integer(kind=c_int), intent(in) :: dj_d2
        integer(kind=c_int), intent(in) :: x_d0
        integer(kind=c_int), intent(in) :: cilm_d
        real(kind=c_double), dimension(2,cilm_d,cilm_d),intent(in) :: cilm
        real(kind=c_double), dimension(x_d0),intent(in) :: x
        real(kind=c_double), dimension(dj_d0,dj_d1,dj_d2),intent(in) :: dj
        real(kind=c_double), dimension(cilmrot_d0,cilmrot_d1,cilmrot_d2),intent(out) :: cilmrot
        integer(kind=c_int), intent(in) :: lmax
        integer(kind=c_int), intent(out) :: exitstatus
        call SHRotateRealCoef(cilmrot,cilm,lmax,x,dj,exitstatus)
    end subroutine cbindSHRotateRealCoef

    function cbindSHPowerL(c,l,c_d0,c_d1,c_d2)  bind(c, name="cbind_sh_power_l")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHPowerL
        implicit none
        integer(kind=c_int), intent(in) :: c_d0
        integer(kind=c_int), intent(in) :: c_d1
        integer(kind=c_int), intent(in) :: c_d2
        real(kind=c_double) :: cbindSHPowerL
        real(kind=c_double), dimension(c_d0,c_d1,c_d2),intent(in) :: c
        integer(kind=c_int), intent(in) :: l
        cbindSHPowerL=SHPowerL(c,l)
    end function cbindSHPowerL

    function cbindSHPowerDensityL(c,l,c_d0,c_d1,c_d2)  bind(c, name="cbind_sh_power_density_l")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHPowerDensityL
        implicit none
        integer(kind=c_int), intent(in) :: c_d0
        integer(kind=c_int), intent(in) :: c_d1
        integer(kind=c_int), intent(in) :: c_d2
        real(kind=c_double) :: cbindSHPowerDensityL
        real(kind=c_double), dimension(c_d0,c_d1,c_d2),intent(in) :: c
        integer(kind=c_int), intent(in) :: l
        cbindSHPowerDensityL=SHPowerDensityL(c,l)
    end function cbindSHPowerDensityL

    function cbindSHCrossPowerL(c1,c2,l,c1_d0,c1_d1,c1_d2,c2_d0,c2_d1,c2_d2)  bind(c&
                                  , name="cbind_sh_cross_power_l")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHCrossPowerL
        implicit none
        integer(kind=c_int), intent(in) :: c2_d0
        integer(kind=c_int), intent(in) :: c2_d1
        integer(kind=c_int), intent(in) :: c2_d2
        integer(kind=c_int), intent(in) :: c1_d0
        integer(kind=c_int), intent(in) :: c1_d1
        integer(kind=c_int), intent(in) :: c1_d2
        real(kind=c_double) :: cbindSHCrossPowerL
        real(kind=c_double), dimension(c1_d0,c1_d1,c1_d2),intent(in) :: c1
        real(kind=c_double), dimension(c2_d0,c2_d1,c2_d2),intent(in) :: c2
        integer(kind=c_int), intent(in) :: l
        cbindSHCrossPowerL=SHCrossPowerL(c1,c2,l)
    end function cbindSHCrossPowerL

    function cbindSHCrossPowerDensityL(c1,c2,l,c1_d0,c1_d1,c1_d2,c2_d0,c2_d1,c2_d2)  bind(c&
                                         , name="cbind_sh_cross_power_density_l")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHCrossPowerDensityL
        implicit none
        integer(kind=c_int), intent(in) :: c2_d0
        integer(kind=c_int), intent(in) :: c2_d1
        integer(kind=c_int), intent(in) :: c2_d2
        integer(kind=c_int), intent(in) :: c1_d0
        integer(kind=c_int), intent(in) :: c1_d1
        integer(kind=c_int), intent(in) :: c1_d2
        real(kind=c_double) :: cbindSHCrossPowerDensityL
        real(kind=c_double), dimension(c1_d0,c1_d1,c1_d2),intent(in) :: c1
        real(kind=c_double), dimension(c2_d0,c2_d1,c2_d2),intent(in) :: c2
        integer(kind=c_int), intent(in) :: l
        cbindSHCrossPowerDensityL=SHCrossPowerDensityL(c1,c2,l)
    end function cbindSHCrossPowerDensityL

    subroutine cbindSHPowerSpectrum(c,lmax,spectra,exitstatus,c_d0,c_d1,c_d2,spectra_d0)  bind(c&
                                     , name="cbind_sh_power_spectrum")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHPowerSpectrum
        implicit none
        integer(kind=c_int), intent(in) :: spectra_d0
        integer(kind=c_int), intent(in) :: c_d0
        integer(kind=c_int), intent(in) :: c_d1
        integer(kind=c_int), intent(in) :: c_d2
        real(kind=c_double), dimension(c_d0,c_d1,c_d2),intent(in) :: c
        integer(kind=c_int), intent(in) :: lmax
        real(kind=c_double), dimension(spectra_d0),intent(out) :: spectra
        integer(kind=c_int), intent(out) :: exitstatus
        call SHPowerSpectrum(c,lmax,spectra,exitstatus)
    end subroutine cbindSHPowerSpectrum

    subroutine cbindSHPowerSpectrumDensity(c,lmax,spectra,exitstatus,c_d0,c_d1,c_d2&
                                            ,spectra_d0)  bind(c, name="cbind_sh_power_spectrum_density")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHPowerSpectrumDensity
        implicit none
        integer(kind=c_int), intent(in) :: spectra_d0
        integer(kind=c_int), intent(in) :: c_d0
        integer(kind=c_int), intent(in) :: c_d1
        integer(kind=c_int), intent(in) :: c_d2
        real(kind=c_double), dimension(c_d0,c_d1,c_d2),intent(in) :: c
        integer(kind=c_int), intent(in) :: lmax
        real(kind=c_double), dimension(spectra_d0),intent(out) :: spectra
        integer(kind=c_int), intent(out) :: exitstatus
        call SHPowerSpectrumDensity(c,lmax,spectra,exitstatus)
    end subroutine cbindSHPowerSpectrumDensity

    subroutine cbindSHCrossPowerSpectrum(c1,c2,lmax,cspectra,exitstatus,c1_d0,c1_d1&
                                           ,c1_d2,c2_d0,c2_d1,c2_d2,cspectra_d0)  bind(c&
                                           , name="cbind_sh_cross_power_spectrum")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHCrossPowerSpectrum
        implicit none
        integer(kind=c_int), intent(in) :: cspectra_d0
        integer(kind=c_int), intent(in) :: c2_d0
        integer(kind=c_int), intent(in) :: c2_d1
        integer(kind=c_int), intent(in) :: c2_d2
        integer(kind=c_int), intent(in) :: c1_d0
        integer(kind=c_int), intent(in) :: c1_d1
        integer(kind=c_int), intent(in) :: c1_d2
        real(kind=c_double), dimension(c1_d0,c1_d1,c1_d2),intent(in) :: c1
        real(kind=c_double), dimension(c2_d0,c2_d1,c2_d2),intent(in) :: c2
        integer(kind=c_int), intent(in) :: lmax
        real(kind=c_double), dimension(cspectra_d0),intent(out) :: cspectra
        integer(kind=c_int), intent(out) :: exitstatus
        call SHCrossPowerSpectrum(c1,c2,lmax,cspectra,exitstatus)
    end subroutine cbindSHCrossPowerSpectrum

    subroutine cbindSHCrossPowerSpectrumDensity(c1,c2,lmax,cspectra,exitstatus,c1_d0&
                                                  ,c1_d1,c1_d2,c2_d0,c2_d1,c2_d2,cspectra_d0)  bind(c&
                                                  , name="cbind_sh_cross_power_spectrum_density")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHCrossPowerSpectrumDensity
        implicit none
        integer(kind=c_int), intent(in) :: cspectra_d0
        integer(kind=c_int), intent(in) :: c2_d0
        integer(kind=c_int), intent(in) :: c2_d1
        integer(kind=c_int), intent(in) :: c2_d2
        integer(kind=c_int), intent(in) :: c1_d0
        integer(kind=c_int), intent(in) :: c1_d1
        integer(kind=c_int), intent(in) :: c1_d2
        real(kind=c_double), dimension(c1_d0,c1_d1,c1_d2),intent(in) :: c1
        real(kind=c_double), dimension(c2_d0,c2_d1,c2_d2),intent(in) :: c2
        integer(kind=c_int), intent(in) :: lmax
        real(kind=c_double), dimension(cspectra_d0),intent(out) :: cspectra
        integer(kind=c_int), intent(out) :: exitstatus
        call SHCrossPowerSpectrumDensity(c1,c2,lmax,cspectra,exitstatus)
    end subroutine cbindSHCrossPowerSpectrumDensity

    subroutine cbindSHAdmitCorr(G,T,lmax,admit,corr,admit_error,exitstatus,G_d0,G_d1&
                                 ,G_d2,T_d0,T_d1,T_d2,admit_d0,corr_d0,admit_error_d0)  bind(c&
                                 , name="cbind_sh_admit_corr")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHAdmitCorr
        implicit none
        integer(kind=c_int), intent(in) :: admit_error_d0
        integer(kind=c_int), intent(in) :: corr_d0
        integer(kind=c_int), intent(in) :: admit_d0
        integer(kind=c_int), intent(in) :: T_d0
        integer(kind=c_int), intent(in) :: T_d1
        integer(kind=c_int), intent(in) :: T_d2
        integer(kind=c_int), intent(in) :: G_d0
        integer(kind=c_int), intent(in) :: G_d1
        integer(kind=c_int), intent(in) :: G_d2
        real(kind=c_double), dimension(G_d0,G_d1,G_d2),intent(in) :: G
        real(kind=c_double), dimension(T_d0,T_d1,T_d2),intent(in) :: T
        integer(kind=c_int), intent(in) :: lmax
        real(kind=c_double), dimension(admit_d0),intent(out) :: admit
        real(kind=c_double), dimension(corr_d0),intent(out) :: corr
        real(kind=c_double), dimension(admit_error_d0),intent(out) :: admit_error
        integer(kind=c_int), intent(out) :: exitstatus
        call SHAdmitCorr(G,T,lmax,admit,corr,admit_error,exitstatus)
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

    function cbindSHPowerLC(c,l,c_d0,c_d1,c_d2)  bind(c, name="cbind_sh_power_lc")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHPowerLC
        implicit none
        integer(kind=c_int), intent(in) :: c_d0
        integer(kind=c_int), intent(in) :: c_d1
        integer(kind=c_int), intent(in) :: c_d2
        real(kind=c_double) :: cbindSHPowerLC
        complex(kind=c_double_complex), dimension(c_d0,c_d1,c_d2),intent(in) :: c
        integer(kind=c_int), intent(in) :: l
        cbindSHPowerLC=SHPowerLC(c,l)
    end function cbindSHPowerLC

    function cbindSHPowerDensityLC(c,l,c_d0,c_d1,c_d2)  bind(c, name="cbind_sh_power_density_lc")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHPowerDensityLC
        implicit none
        integer(kind=c_int), intent(in) :: c_d0
        integer(kind=c_int), intent(in) :: c_d1
        integer(kind=c_int), intent(in) :: c_d2
        real(kind=c_double) :: cbindSHPowerDensityLC
        complex(kind=c_double_complex), dimension(c_d0,c_d1,c_d2),intent(in) :: c
        integer(kind=c_int), intent(in) :: l
        cbindSHPowerDensityLC=SHPowerDensityLC(c,l)
    end function cbindSHPowerDensityLC

    function cbindSHCrossPowerLC(c1,c2,l,c1_d0,c1_d1,c1_d2,c2_d0,c2_d1,c2_d2)  bind(c&
                                   , name="cbind_sh_cross_power_lc")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHCrossPowerLC
        implicit none
        integer(kind=c_int), intent(in) :: c2_d0
        integer(kind=c_int), intent(in) :: c2_d1
        integer(kind=c_int), intent(in) :: c2_d2
        integer(kind=c_int), intent(in) :: c1_d0
        integer(kind=c_int), intent(in) :: c1_d1
        integer(kind=c_int), intent(in) :: c1_d2
        complex(kind=c_double_complex) :: cbindSHCrossPowerLC
        complex(kind=c_double_complex), dimension(c1_d0,c1_d1,c1_d2),intent(in) :: c1
        complex(kind=c_double_complex), dimension(c2_d0,c2_d1,c2_d2),intent(in) :: c2
        integer(kind=c_int), intent(in) :: l
        cbindSHCrossPowerLC=SHCrossPowerLC(c1,c2,l)
    end function cbindSHCrossPowerLC

    function cbindSHCrossPowerDensityLC(c1,c2,l,c1_d0,c1_d1,c1_d2,c2_d0,c2_d1,c2_d2)  bind(c&
                                          , name="cbind_sh_cross_power_density_lc")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHCrossPowerDensityLC
        implicit none
        integer(kind=c_int), intent(in) :: c2_d0
        integer(kind=c_int), intent(in) :: c2_d1
        integer(kind=c_int), intent(in) :: c2_d2
        integer(kind=c_int), intent(in) :: c1_d0
        integer(kind=c_int), intent(in) :: c1_d1
        integer(kind=c_int), intent(in) :: c1_d2
        complex(kind=c_double_complex) :: cbindSHCrossPowerDensityLC
        complex(kind=c_double_complex), dimension(c1_d0,c1_d1,c1_d2),intent(in) :: c1
        complex(kind=c_double_complex), dimension(c2_d0,c2_d1,c2_d2),intent(in) :: c2
        integer(kind=c_int), intent(in) :: l
        cbindSHCrossPowerDensityLC=SHCrossPowerDensityLC(c1,c2,l)
    end function cbindSHCrossPowerDensityLC

    subroutine cbindSHPowerSpectrumC(c,lmax,spectra,exitstatus,c_d0,c_d1,c_d2,spectra_d0)  bind(c&
                                      , name="cbind_sh_power_spectrum_c")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHPowerSpectrumC
        implicit none
        integer(kind=c_int), intent(in) :: spectra_d0
        integer(kind=c_int), intent(in) :: c_d0
        integer(kind=c_int), intent(in) :: c_d1
        integer(kind=c_int), intent(in) :: c_d2
        complex(kind=c_double_complex), dimension(c_d0,c_d1,c_d2),intent(in) :: c
        integer(kind=c_int), intent(in) :: lmax
        real(kind=c_double), dimension(spectra_d0),intent(out) :: spectra
        integer(kind=c_int), intent(out) :: exitstatus
        call SHPowerSpectrumC(c,lmax,spectra,exitstatus)
    end subroutine cbindSHPowerSpectrumC

    subroutine cbindSHPowerSpectrumDensityC(c,lmax,spectra,exitstatus,c_d0,c_d1,c_d2&
                                             ,spectra_d0)  bind(c, name="cbind_sh_power_spectrum_density_c")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHPowerSpectrumDensityC
        implicit none
        integer(kind=c_int), intent(in) :: spectra_d0
        integer(kind=c_int), intent(in) :: c_d0
        integer(kind=c_int), intent(in) :: c_d1
        integer(kind=c_int), intent(in) :: c_d2
        complex(kind=c_double_complex), dimension(c_d0,c_d1,c_d2),intent(in) :: c
        integer(kind=c_int), intent(in) :: lmax
        real(kind=c_double), dimension(spectra_d0),intent(out) :: spectra
        integer(kind=c_int), intent(out) :: exitstatus
        call SHPowerSpectrumDensityC(c,lmax,spectra,exitstatus)
    end subroutine cbindSHPowerSpectrumDensityC

    subroutine cbindSHCrossPowerSpectrumC(c1,c2,lmax,cspectra,exitstatus,c1_d0,c1_d1&
                                            ,c1_d2,c2_d0,c2_d1,c2_d2,cspectra_d0)  bind(c&
                                            , name="cbind_sh_cross_power_spectrum_c")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHCrossPowerSpectrumC
        implicit none
        integer(kind=c_int), intent(in) :: cspectra_d0
        integer(kind=c_int), intent(in) :: c2_d0
        integer(kind=c_int), intent(in) :: c2_d1
        integer(kind=c_int), intent(in) :: c2_d2
        integer(kind=c_int), intent(in) :: c1_d0
        integer(kind=c_int), intent(in) :: c1_d1
        integer(kind=c_int), intent(in) :: c1_d2
        complex(kind=c_double_complex), dimension(c1_d0,c1_d1,c1_d2),intent(in) :: c1
        complex(kind=c_double_complex), dimension(c2_d0,c2_d1,c2_d2),intent(in) :: c2
        integer(kind=c_int), intent(in) :: lmax
        complex(kind=c_double_complex), dimension(cspectra_d0),intent(out) :: cspectra
        integer(kind=c_int), intent(out) :: exitstatus
        call SHCrossPowerSpectrumC(c1,c2,lmax,cspectra,exitstatus)
    end subroutine cbindSHCrossPowerSpectrumC

    subroutine cbindSHCrossPowerSpectrumDensityC(c1,c2,lmax,cspectra,exitstatus,c1_d0&
                                                   ,c1_d1,c1_d2,c2_d0,c2_d1,c2_d2&
                                                   ,cspectra_d0)  bind(c, name="cbind_sh_cross_power_spectrum_density_c")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHCrossPowerSpectrumDensityC
        implicit none
        integer(kind=c_int), intent(in) :: cspectra_d0
        integer(kind=c_int), intent(in) :: c2_d0
        integer(kind=c_int), intent(in) :: c2_d1
        integer(kind=c_int), intent(in) :: c2_d2
        integer(kind=c_int), intent(in) :: c1_d0
        integer(kind=c_int), intent(in) :: c1_d1
        integer(kind=c_int), intent(in) :: c1_d2
        complex(kind=c_double_complex), dimension(c1_d0,c1_d1,c1_d2),intent(in) :: c1
        complex(kind=c_double_complex), dimension(c2_d0,c2_d1,c2_d2),intent(in) :: c2
        integer(kind=c_int), intent(in) :: lmax
        complex(kind=c_double_complex), dimension(cspectra_d0),intent(out) :: cspectra
        integer(kind=c_int), intent(out) :: exitstatus
        call SHCrossPowerSpectrumDensityC(c1,c2,lmax,cspectra,exitstatus)
    end subroutine cbindSHCrossPowerSpectrumDensityC

    subroutine cbindSHMultiTaperSE(mtse,sd,sh,lmax,tapers,taper_order,lmaxt,k,alpha&
                                       ,lat,lon,taper_wt,norm,csphase,exitstatus,mtse_d0&
                                       ,sd_d0,sh_d0,sh_d1,sh_d2,tapers_d0,tapers_d1&
                                       ,taper_order_d0,alpha_d0,taper_wt_d0)  bind(c&
                                       , name="cbind_sh_multi_taper_se")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHMultiTaperSE
        implicit none
        integer(kind=c_int), intent(in) :: taper_wt_d0
        integer(kind=c_int), intent(in) :: alpha_d0
        integer(kind=c_int), intent(in) :: taper_order_d0
        integer(kind=c_int), intent(in) :: tapers_d0
        integer(kind=c_int), intent(in) :: tapers_d1
        integer(kind=c_int), intent(in) :: sh_d0
        integer(kind=c_int), intent(in) :: sh_d1
        integer(kind=c_int), intent(in) :: sh_d2
        integer(kind=c_int), intent(in) :: sd_d0
        integer(kind=c_int), intent(in) :: mtse_d0
        real(kind=c_double), dimension(mtse_d0),intent(out) :: mtse
        real(kind=c_double), dimension(sd_d0),intent(out) :: sd
        real(kind=c_double), dimension(sh_d0,sh_d1,sh_d2),intent(in) :: sh
        real(kind=c_double), dimension(tapers_d0,tapers_d1),intent(in) :: tapers
        integer(kind=c_int), intent(in) :: lmax
        integer(kind=c_int), intent(in) :: lmaxt
        integer(kind=c_int), intent(in) :: k
        integer(kind=c_int), dimension(taper_order_d0),intent(in) :: taper_order
        real(kind=c_double), dimension(alpha_d0),intent(in) :: alpha
        real(kind=c_double), intent(in) :: lat
        real(kind=c_double), intent(in) :: lon
        real(kind=c_double), dimension(taper_wt_d0),intent(in) :: taper_wt
        integer(kind=c_int), intent(in) :: csphase
        integer(kind=c_int), intent(in) :: norm
        integer(kind=c_int), intent(out) :: exitstatus
        call SHMultiTaperSE(mtse,sd,sh,lmax,tapers,taper_order,lmaxt,k,alpha,lat,lon&
                                ,taper_wt,norm,csphase,exitstatus)
    end subroutine cbindSHMultiTaperSE

    subroutine cbindSHMultiTaperCSE(mtse,sd,sh1,lmax1,sh2,lmax2,tapers,taper_order&
                                        ,lmaxt,k,alpha,lat,lon,taper_wt,norm,csphase&
                                        ,exitstatus,mtse_d0,sd_d0,sh1_d0,sh1_d1,sh1_d2&
                                        ,sh2_d0,sh2_d1,sh2_d2,tapers_d0,tapers_d1&
                                        ,taper_order_d0,alpha_d0,taper_wt_d0)  bind(c&
                                        , name="cbind_sh_multi_taper_cse")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHMultiTaperCSE
        implicit none
        integer(kind=c_int), intent(in) :: taper_wt_d0
        integer(kind=c_int), intent(in) :: alpha_d0
        integer(kind=c_int), intent(in) :: taper_order_d0
        integer(kind=c_int), intent(in) :: tapers_d0
        integer(kind=c_int), intent(in) :: tapers_d1
        integer(kind=c_int), intent(in) :: sh2_d0
        integer(kind=c_int), intent(in) :: sh2_d1
        integer(kind=c_int), intent(in) :: sh2_d2
        integer(kind=c_int), intent(in) :: sh1_d0
        integer(kind=c_int), intent(in) :: sh1_d1
        integer(kind=c_int), intent(in) :: sh1_d2
        integer(kind=c_int), intent(in) :: sd_d0
        integer(kind=c_int), intent(in) :: mtse_d0
        real(kind=c_double), dimension(mtse_d0),intent(out) :: mtse
        real(kind=c_double), dimension(sd_d0),intent(out) :: sd
        real(kind=c_double), dimension(sh1_d0,sh1_d1,sh1_d2),intent(in) :: sh1
        real(kind=c_double), dimension(sh2_d0,sh2_d1,sh2_d2),intent(in) :: sh2
        real(kind=c_double), dimension(tapers_d0,tapers_d1),intent(in) :: tapers
        integer(kind=c_int), intent(in) :: lmax1
        integer(kind=c_int), intent(in) :: lmax2
        integer(kind=c_int), intent(in) :: lmaxt
        integer(kind=c_int), intent(in) :: k
        integer(kind=c_int), dimension(taper_order_d0),intent(in) :: taper_order
        real(kind=c_double), dimension(alpha_d0),intent(in) :: alpha
        real(kind=c_double), intent(in) :: lat
        real(kind=c_double), intent(in) :: lon
        real(kind=c_double), dimension(taper_wt_d0),intent(in) :: taper_wt
        integer(kind=c_int), intent(in) :: csphase
        integer(kind=c_int), intent(in) :: norm
        integer(kind=c_int), intent(out) :: exitstatus
        call SHMultiTaperCSE(mtse,sd,sh1,lmax1,sh2,lmax2,tapers,taper_order,lmaxt&
                                 ,k,alpha,lat,lon,taper_wt,norm,csphase,exitstatus)
    end subroutine cbindSHMultiTaperCSE

    subroutine cbindSHLocalizedAdmitCorr(tapers,taper_order,lwin,lat,lon,g,t,lmax&
                                               ,admit,corr,k,admit_error,corr_error&
                                               ,taper_wt,mtdef,k1linsig,exitstatus&
                                               ,tapers_d0,tapers_d1,g_d0,g_d1,g_d2&
                                               ,t_d0,t_d1,t_d2,taper_order_d0,admit_d0&
                                               ,corr_d0,admit_error_d0,corr_error_d0&
                                               ,taper_wt_d0)  bind(c, name="cbind_sh_localized_admit_corr")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHLocalizedAdmitCorr
        implicit none
        integer(kind=c_int), intent(in) :: taper_wt_d0
        integer(kind=c_int), intent(in) :: corr_error_d0
        integer(kind=c_int), intent(in) :: admit_error_d0
        integer(kind=c_int), intent(in) :: corr_d0
        integer(kind=c_int), intent(in) :: admit_d0
        integer(kind=c_int), intent(in) :: taper_order_d0
        integer(kind=c_int), intent(in) :: t_d0
        integer(kind=c_int), intent(in) :: t_d1
        integer(kind=c_int), intent(in) :: t_d2
        integer(kind=c_int), intent(in) :: g_d0
        integer(kind=c_int), intent(in) :: g_d1
        integer(kind=c_int), intent(in) :: g_d2
        integer(kind=c_int), intent(in) :: tapers_d0
        integer(kind=c_int), intent(in) :: tapers_d1
        real(kind=c_double), dimension(tapers_d0,tapers_d1),intent(in) :: tapers
        real(kind=c_double), intent(in) :: lat
        real(kind=c_double), intent(in) :: lon
        real(kind=c_double), dimension(g_d0,g_d1,g_d2),intent(in) :: g
        real(kind=c_double), dimension(t_d0,t_d1,t_d2),intent(in) :: t
        integer(kind=c_int), intent(in) :: lwin
        integer(kind=c_int), intent(in) :: lmax
        integer(kind=c_int), intent(in) :: k
        integer(kind=c_int), dimension(taper_order_d0),intent(in) :: taper_order
        real(kind=c_double), dimension(admit_d0),intent(out) :: admit
        real(kind=c_double), dimension(corr_d0),intent(out) :: corr
        real(kind=c_double), dimension(admit_error_d0),intent(out) :: admit_error
        real(kind=c_double), dimension(corr_error_d0),intent(out) :: corr_error
        integer(kind=c_int), intent(in) :: mtdef
        integer(kind=c_int), intent(in) :: k1linsig
        real(kind=c_double), dimension(taper_wt_d0),intent(in) :: taper_wt
        integer(kind=c_int), intent(out) :: exitstatus
        call SHLocalizedAdmitCorr(tapers,taper_order,lwin,lat,lon,g,t,lmax,admit,corr&
                                        ,k,admit_error,corr_error,taper_wt,mtdef,k1linsig&
                                        ,exitstatus)
    end subroutine cbindSHLocalizedAdmitCorr

    subroutine cbindSHReturnTapers(theta0,lmax,tapers,eigenvalues,taper_order,degrees&
                                         ,exitstatus,tapers_d0,tapers_d1,eigenvalues_d0&
                                         ,taper_order_d0,degrees_d0)  bind(c, name="cbind_sh_return_tapers")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHReturnTapers
        implicit none
        integer(kind=c_int), intent(in) :: degrees_d0
        integer(kind=c_int), intent(in) :: taper_order_d0
        integer(kind=c_int), intent(in) :: eigenvalues_d0
        integer(kind=c_int), intent(in) :: tapers_d0
        integer(kind=c_int), intent(in) :: tapers_d1
        real(kind=c_double), intent(in) :: theta0
        integer(kind=c_int), intent(in) :: lmax
        real(kind=c_double), dimension(tapers_d0,tapers_d1),intent(out) :: tapers
        real(kind=c_double), dimension(eigenvalues_d0),intent(out) :: eigenvalues
        integer(kind=c_int), dimension(taper_order_d0),intent(out) :: taper_order
        integer(kind=c_int), dimension(degrees_d0),intent(in) :: degrees
        integer(kind=c_int), intent(out) :: exitstatus
        call SHReturnTapers(theta0,lmax,tapers,eigenvalues,taper_order,degrees,exitstatus)
    end subroutine cbindSHReturnTapers

    subroutine cbindSHReturnTapersM(theta0,lmax,m,tapers,eigenvalues,shannon,degrees&
                                          ,ntapers,exitstatus,tapers_d0,tapers_d1&
                                          ,eigenvalues_d0,degrees_d0)  bind(c, name="cbind_sh_return_tapers_m")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHReturnTapersM
        implicit none
        integer(kind=c_int), intent(in) :: degrees_d0
        integer(kind=c_int), intent(in) :: eigenvalues_d0
        integer(kind=c_int), intent(in) :: tapers_d0
        integer(kind=c_int), intent(in) :: tapers_d1
        real(kind=c_double), intent(in) :: theta0
        integer(kind=c_int), intent(in) :: lmax
        integer(kind=c_int), intent(in) :: m
        real(kind=c_double), dimension(tapers_d0,tapers_d1),intent(out) :: tapers
        real(kind=c_double), dimension(eigenvalues_d0),intent(out) :: eigenvalues
        real(kind=c_double), intent(out) :: shannon
        integer(kind=c_int), dimension(degrees_d0),intent(in) :: degrees
        integer(kind=c_int), intent(out) :: ntapers
        integer(kind=c_int), intent(out) :: exitstatus
        call SHReturnTapersM(theta0,lmax,m,tapers,eigenvalues,shannon,degrees,ntapers&
                                   ,exitstatus)
    end subroutine cbindSHReturnTapersM

    subroutine cbindComputeDm(dllm,lmax,m,theta0,degrees,exitstatus,dllm_d0,dllm_d1&
                                  ,degrees_d0)  bind(c, name="cbind_compute_dm")
        use, intrinsic :: iso_c_binding
        use shtools, only: ComputeDm
        implicit none
        integer(kind=c_int), intent(in) :: degrees_d0
        integer(kind=c_int), intent(in) :: dllm_d0
        integer(kind=c_int), intent(in) :: dllm_d1
        real(kind=c_double), dimension(dllm_d0,dllm_d1),intent(out) :: dllm
        real(kind=c_double), intent(in) :: theta0
        integer(kind=c_int), intent(in) :: lmax
        integer(kind=c_int), intent(in) :: m
        integer(kind=c_int), dimension(degrees_d0),intent(in) :: degrees
        integer(kind=c_int), intent(out) :: exitstatus
        call ComputeDm(dllm,lmax,m,theta0,degrees,exitstatus)
    end subroutine cbindComputeDm

    subroutine cbindComputeDG82(dG82,lmax,m,theta0,exitstatus,dG82_d0,dG82_d1)  bind(c&
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
        integer(kind=c_int), intent(out) :: exitstatus
        call ComputeDG82(dG82,lmax,m,theta0,exitstatus)
    end subroutine cbindComputeDG82

    function cbindSHFindLWin(theta0,m,alpha,taper_number)  bind(c, name="cbind_sh_find_l_win")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHFindLWin
        implicit none
        integer(kind=c_int) :: cbindSHFindLWin
        real(kind=c_double), intent(in) :: theta0
        real(kind=c_double), intent(in) :: alpha
        integer(kind=c_int), intent(in) :: m
        integer(kind=c_int), intent(in) :: taper_number
        cbindSHFindLWin=SHFindLWin(theta0,m,alpha,taper_number)
    end function cbindSHFindLWin

    subroutine cbindSHBiasK(tapers,lwin,k,incspectra,ldata,outcspectra,taper_wt,save_cg&
                                  ,exitstatus,tapers_d0,tapers_d1,incspectra_d0,outcspectra_d0&
                                  ,taper_wt_d0)  bind(c, name="cbind_sh_bias_k")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHBiasK
        implicit none
        integer(kind=c_int), intent(in) :: taper_wt_d0
        integer(kind=c_int), intent(in) :: outcspectra_d0
        integer(kind=c_int), intent(in) :: incspectra_d0
        integer(kind=c_int), intent(in) :: tapers_d0
        integer(kind=c_int), intent(in) :: tapers_d1
        real(kind=c_double), dimension(tapers_d0,tapers_d1),intent(in) :: tapers
        real(kind=c_double), dimension(incspectra_d0),intent(in) :: incspectra
        real(kind=c_double), dimension(outcspectra_d0),intent(out) :: outcspectra
        integer(kind=c_int), intent(in) :: lwin
        integer(kind=c_int), intent(in) :: ldata
        integer(kind=c_int), intent(in) :: k
        real(kind=c_double), dimension(taper_wt_d0),intent(in) :: taper_wt
        integer(kind=c_int), intent(in) :: save_cg
        integer(kind=c_int), intent(out) :: exitstatus
        call SHBiasK(tapers,lwin,k,incspectra,ldata,outcspectra,taper_wt,save_cg,exitstatus)
    end subroutine cbindSHBiasK

    subroutine cbindSHMTCouplingMatrix(Mmt,lmax,tapers_power,lwin,k,taper_wt,exitstatus&
                                          ,Mmt_d0,Mmt_d1,tapers_power_d0,tapers_power_d1&
                                          ,taper_wt_d0)  bind(c, name="cbind_shmt_coupling_matrix")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHMTCouplingMatrix
        implicit none
        integer(kind=c_int), intent(in) :: taper_wt_d0
        integer(kind=c_int), intent(in) :: tapers_power_d0
        integer(kind=c_int), intent(in) :: tapers_power_d1
        integer(kind=c_int), intent(in) :: Mmt_d0
        integer(kind=c_int), intent(in) :: Mmt_d1
        real(kind=c_double), dimension(Mmt_d0,Mmt_d1),intent(out) :: Mmt
        real(kind=c_double), dimension(tapers_power_d0,tapers_power_d1),intent(in) :: tapers_power
        real(kind=c_double), dimension(taper_wt_d0),intent(in) :: taper_wt
        integer(kind=c_int), intent(in) :: lmax
        integer(kind=c_int), intent(in) :: k
        integer(kind=c_int), intent(in) :: lwin
        integer(kind=c_int), intent(out) :: exitstatus
        call SHMTCouplingMatrix(Mmt,lmax,tapers_power,lwin,k,taper_wt,exitstatus)
    end subroutine cbindSHMTCouplingMatrix

    subroutine cbindSHBiasAdmitCorr(sgt,sgg,stt,lmax,tapers,lwin,k,admit,corr,mtdef&
                                       ,taper_wt,exitstatus,sgt_d0,sgg_d0,stt_d0,tapers_d0&
                                       ,tapers_d1,admit_d0,corr_d0,taper_wt_d0)  bind(c&
                                       , name="cbind_sh_bias_admit_corr")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHBiasAdmitCorr
        implicit none
        integer(kind=c_int), intent(in) :: taper_wt_d0
        integer(kind=c_int), intent(in) :: corr_d0
        integer(kind=c_int), intent(in) :: admit_d0
        integer(kind=c_int), intent(in) :: tapers_d0
        integer(kind=c_int), intent(in) :: tapers_d1
        integer(kind=c_int), intent(in) :: stt_d0
        integer(kind=c_int), intent(in) :: sgg_d0
        integer(kind=c_int), intent(in) :: sgt_d0
        real(kind=c_double), dimension(sgt_d0),intent(in) :: sgt
        real(kind=c_double), dimension(sgg_d0),intent(in) :: sgg
        real(kind=c_double), dimension(stt_d0),intent(in) :: stt
        real(kind=c_double), dimension(tapers_d0,tapers_d1),intent(in) :: tapers
        integer(kind=c_int), intent(in) :: lmax
        integer(kind=c_int), intent(in) :: lwin
        integer(kind=c_int), intent(in) :: k
        real(kind=c_double), dimension(admit_d0),intent(out) :: admit
        real(kind=c_double), dimension(corr_d0),intent(out) :: corr
        integer(kind=c_int), intent(in) :: mtdef
        real(kind=c_double), dimension(taper_wt_d0),intent(in) :: taper_wt
        integer(kind=c_int), intent(out) :: exitstatus
        call SHBiasAdmitCorr(sgt,sgg,stt,lmax,tapers,lwin,k,admit,corr,mtdef,taper_wt&
                                ,exitstatus)
    end subroutine cbindSHBiasAdmitCorr

    subroutine cbindSHMTDebias(mtdebias,mtspectra,lmax,tapers,lwin,k,nl,lmid,n,taper_wt&
                                       ,exitstatus,mtdebias_d0,mtdebias_d1,lmid_d0&
                                       ,mtspectra_d0,mtspectra_d1,tapers_d0,tapers_d1&
                                       ,taper_wt_d0)  bind(c, name="cbind_shmt_debias")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHMTDebias
        implicit none
        integer(kind=c_int), intent(in) :: taper_wt_d0
        integer(kind=c_int), intent(in) :: tapers_d0
        integer(kind=c_int), intent(in) :: tapers_d1
        integer(kind=c_int), intent(in) :: mtspectra_d0
        integer(kind=c_int), intent(in) :: mtspectra_d1
        integer(kind=c_int), intent(in) :: lmid_d0
        integer(kind=c_int), intent(in) :: mtdebias_d0
        integer(kind=c_int), intent(in) :: mtdebias_d1
        real(kind=c_double), dimension(mtdebias_d0,mtdebias_d1),intent(out) :: mtdebias
        real(kind=c_double), dimension(lmid_d0),intent(out) :: lmid
        real(kind=c_double), dimension(mtspectra_d0,mtspectra_d1),intent(in) :: mtspectra
        real(kind=c_double), dimension(tapers_d0,tapers_d1),intent(in) :: tapers
        real(kind=c_double), dimension(taper_wt_d0),intent(in) :: taper_wt
        integer(kind=c_int), intent(in) :: lmax
        integer(kind=c_int), intent(in) :: k
        integer(kind=c_int), intent(in) :: lwin
        integer(kind=c_int), intent(in) :: nl
        integer(kind=c_int), intent(out) :: n
        integer(kind=c_int), intent(out) :: exitstatus
        call SHMTDebias(mtdebias,mtspectra,lmax,tapers,lwin,k,nl,lmid,n,taper_wt,exitstatus)
    end subroutine cbindSHMTDebias

    subroutine cbindSHMTVarOpt(l,tapers,taper_order,lwin,kmax,Sff,var_opt,var_unit&
                                ,weight_opt,unweighted_covar,nocross,exitstatus,tapers_d0&
                                ,tapers_d1,Sff_d0,var_opt_d0,var_unit_d0,taper_order_d0&
                                ,weight_opt_d0,weight_opt_d1,unweighted_covar_d0,unweighted_covar_d1)  bind(c&
                                , name="cbind_shmt_var_opt")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHMTVarOpt
        implicit none
        integer(kind=c_int), intent(in) :: unweighted_covar_d0
        integer(kind=c_int), intent(in) :: unweighted_covar_d1
        integer(kind=c_int), intent(in) :: weight_opt_d0
        integer(kind=c_int), intent(in) :: weight_opt_d1
        integer(kind=c_int), intent(in) :: taper_order_d0
        integer(kind=c_int), intent(in) :: var_unit_d0
        integer(kind=c_int), intent(in) :: var_opt_d0
        integer(kind=c_int), intent(in) :: Sff_d0
        integer(kind=c_int), intent(in) :: tapers_d0
        integer(kind=c_int), intent(in) :: tapers_d1
        real(kind=c_double), dimension(tapers_d0,tapers_d1),intent(in) :: tapers
        real(kind=c_double), dimension(Sff_d0),intent(in) :: Sff
        real(kind=c_double), dimension(var_opt_d0),intent(out) :: var_opt
        real(kind=c_double), dimension(var_unit_d0),intent(out) :: var_unit
        integer(kind=c_int), intent(in) :: l
        integer(kind=c_int), intent(in) :: lwin
        integer(kind=c_int), intent(in) :: kmax
        integer(kind=c_int), dimension(taper_order_d0),intent(in) :: taper_order
        real(kind=c_double), dimension(weight_opt_d0,weight_opt_d1),intent(out) :: weight_opt
        real(kind=c_double), dimension(unweighted_covar_d0,unweighted_covar_d1),intent(out) :: unweighted_covar
        integer(kind=c_int), intent(in) :: nocross
        integer(kind=c_int), intent(out) :: exitstatus
        call SHMTVarOpt(l,tapers,taper_order,lwin,kmax,Sff,var_opt,var_unit,weight_opt&
                         ,unweighted_covar,nocross,exitstatus)
    end subroutine cbindSHMTVarOpt

    subroutine cbindSHMTVar(l,tapers,taper_order,lwin,kmax,Sff,variance,taper_wt,unweighted_covar&
                             ,nocross,exitstatus,tapers_d0,tapers_d1,Sff_d0,taper_order_d0&
                             ,taper_wt_d0,unweighted_covar_d0,unweighted_covar_d1)  bind(c&
                             , name="cbind_shmt_var")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHMTVar
        implicit none
        integer(kind=c_int), intent(in) :: unweighted_covar_d0
        integer(kind=c_int), intent(in) :: unweighted_covar_d1
        integer(kind=c_int), intent(in) :: taper_wt_d0
        integer(kind=c_int), intent(in) :: taper_order_d0
        integer(kind=c_int), intent(in) :: Sff_d0
        integer(kind=c_int), intent(in) :: tapers_d0
        integer(kind=c_int), intent(in) :: tapers_d1
        real(kind=c_double), dimension(tapers_d0,tapers_d1),intent(in) :: tapers
        real(kind=c_double), dimension(Sff_d0),intent(in) :: Sff
        real(kind=c_double), intent(out) :: variance
        integer(kind=c_int), intent(in) :: l
        integer(kind=c_int), intent(in) :: lwin
        integer(kind=c_int), intent(in) :: kmax
        integer(kind=c_int), dimension(taper_order_d0),intent(in) :: taper_order
        real(kind=c_double), dimension(taper_wt_d0),intent(in) :: taper_wt
        real(kind=c_double), dimension(unweighted_covar_d0,unweighted_covar_d1),intent(out) :: unweighted_covar
        integer(kind=c_int), intent(in) :: nocross
        integer(kind=c_int), intent(out) :: exitstatus
        call SHMTVar(l,tapers,taper_order,lwin,kmax,Sff,variance,taper_wt,unweighted_covar&
                      ,nocross,exitstatus)
    end subroutine cbindSHMTVar

    function cbindSHSjkPG(incspectra,l,m,mprime,hj_real,hk_real,mj,mk,lwin,hkcc,incspectra_d0&
                                    ,hj_real_d0,hk_real_d0)  bind(c, name="cbind_sh_sjk_pg")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHSjkPG
        implicit none
        integer(kind=c_int), intent(in) :: hk_real_d0
        integer(kind=c_int), intent(in) :: hj_real_d0
        integer(kind=c_int), intent(in) :: incspectra_d0
        complex(kind=c_double_complex) :: cbindSHSjkPG
        real(kind=c_double), dimension(incspectra_d0),intent(in) :: incspectra
        real(kind=c_double), dimension(hj_real_d0),intent(in) :: hj_real
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

    subroutine cbindSHReturnTapersMap(tapers,eigenvalues,dh_mask,n_dh,lmax,sampling&
                                            ,ntapers,degrees,exitstatus,tapers_d0&
                                            ,tapers_d1,eigenvalues_d0,dh_mask_d0,dh_mask_d1&
                                            ,degrees_d0)  bind(c, name="cbind_sh_return_tapers_map")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHReturnTapersMap
        implicit none
        integer(kind=c_int), intent(in) :: degrees_d0
        integer(kind=c_int), intent(in) :: dh_mask_d0
        integer(kind=c_int), intent(in) :: dh_mask_d1
        integer(kind=c_int), intent(in) :: eigenvalues_d0
        integer(kind=c_int), intent(in) :: tapers_d0
        integer(kind=c_int), intent(in) :: tapers_d1
        real(kind=c_double), dimension(tapers_d0,tapers_d1),intent(out) :: tapers
        real(kind=c_double), dimension(eigenvalues_d0),intent(out) :: eigenvalues
        integer(kind=c_int), dimension(dh_mask_d0,dh_mask_d1),intent(in) :: dh_mask
        integer(kind=c_int), intent(in) :: n_dh
        integer(kind=c_int), intent(in) :: lmax
        integer(kind=c_int), intent(in) :: sampling
        integer(kind=c_int), intent(in) :: ntapers
        integer(kind=c_int), dimension(degrees_d0),intent(in) :: degrees
        integer(kind=c_int), intent(out) :: exitstatus
        call SHReturnTapersMap(tapers,eigenvalues,dh_mask,n_dh,lmax,sampling,ntapers&
                                     ,degrees,exitstatus)
    end subroutine cbindSHReturnTapersMap

    subroutine cbindSHBiasKMask(tapers,lwin,k,incspectra,ldata,outcspectra,taper_wt&
                                      ,save_cg,exitstatus,tapers_d0,tapers_d1,incspectra_d0&
                                      ,outcspectra_d0,taper_wt_d0)  bind(c, name="cbind_sh_bias_k_mask")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHBiasKMask
        implicit none
        integer(kind=c_int), intent(in) :: taper_wt_d0
        integer(kind=c_int), intent(in) :: outcspectra_d0
        integer(kind=c_int), intent(in) :: incspectra_d0
        integer(kind=c_int), intent(in) :: tapers_d0
        integer(kind=c_int), intent(in) :: tapers_d1
        real(kind=c_double), dimension(tapers_d0,tapers_d1),intent(in) :: tapers
        real(kind=c_double), dimension(incspectra_d0),intent(in) :: incspectra
        real(kind=c_double), dimension(outcspectra_d0),intent(out) :: outcspectra
        integer(kind=c_int), intent(in) :: lwin
        integer(kind=c_int), intent(in) :: ldata
        integer(kind=c_int), intent(in) :: k
        real(kind=c_double), dimension(taper_wt_d0),intent(in) :: taper_wt
        integer(kind=c_int), intent(in) :: save_cg
        integer(kind=c_int), intent(out) :: exitstatus
        call SHBiasKMask(tapers,lwin,k,incspectra,ldata,outcspectra,taper_wt,save_cg&
                               ,exitstatus)
    end subroutine cbindSHBiasKMask

    subroutine cbindSHMultiTaperMaskSE(mtse,sd,sh,lmax,tapers,lmaxt,k,taper_wt,norm&
                                           ,csphase,exitstatus,mtse_d0,sd_d0,sh_d0&
                                           ,sh_d1,sh_d2,tapers_d0,tapers_d1,taper_wt_d0)  bind(c&
                                           , name="cbind_sh_multi_taper_mask_se")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHMultiTaperMaskSE
        implicit none
        integer(kind=c_int), intent(in) :: taper_wt_d0
        integer(kind=c_int), intent(in) :: tapers_d0
        integer(kind=c_int), intent(in) :: tapers_d1
        integer(kind=c_int), intent(in) :: sh_d0
        integer(kind=c_int), intent(in) :: sh_d1
        integer(kind=c_int), intent(in) :: sh_d2
        integer(kind=c_int), intent(in) :: sd_d0
        integer(kind=c_int), intent(in) :: mtse_d0
        real(kind=c_double), dimension(mtse_d0),intent(out) :: mtse
        real(kind=c_double), dimension(sd_d0),intent(out) :: sd
        real(kind=c_double), dimension(sh_d0,sh_d1,sh_d2),intent(in) :: sh
        real(kind=c_double), dimension(tapers_d0,tapers_d1),intent(in) :: tapers
        integer(kind=c_int), intent(in) :: lmax
        integer(kind=c_int), intent(in) :: lmaxt
        integer(kind=c_int), intent(in) :: k
        real(kind=c_double), dimension(taper_wt_d0),intent(in) :: taper_wt
        integer(kind=c_int), intent(in) :: csphase
        integer(kind=c_int), intent(in) :: norm
        integer(kind=c_int), intent(out) :: exitstatus
        call SHMultiTaperMaskSE(mtse,sd,sh,lmax,tapers,lmaxt,k,taper_wt,norm,csphase&
                                    ,exitstatus)
    end subroutine cbindSHMultiTaperMaskSE

    subroutine cbindSHMultiTaperMaskCSE(mtse,sd,sh1,lmax1,sh2,lmax2,tapers,lmaxt,k&
                                            ,taper_wt,norm,csphase,exitstatus,mtse_d0&
                                            ,sd_d0,sh1_d0,sh1_d1,sh1_d2,sh2_d0,sh2_d1&
                                            ,sh2_d2,tapers_d0,tapers_d1,taper_wt_d0)  bind(c&
                                            , name="cbind_sh_multi_taper_mask_cse")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHMultiTaperMaskCSE
        implicit none
        integer(kind=c_int), intent(in) :: taper_wt_d0
        integer(kind=c_int), intent(in) :: tapers_d0
        integer(kind=c_int), intent(in) :: tapers_d1
        integer(kind=c_int), intent(in) :: sh2_d0
        integer(kind=c_int), intent(in) :: sh2_d1
        integer(kind=c_int), intent(in) :: sh2_d2
        integer(kind=c_int), intent(in) :: sh1_d0
        integer(kind=c_int), intent(in) :: sh1_d1
        integer(kind=c_int), intent(in) :: sh1_d2
        integer(kind=c_int), intent(in) :: sd_d0
        integer(kind=c_int), intent(in) :: mtse_d0
        real(kind=c_double), dimension(mtse_d0),intent(out) :: mtse
        real(kind=c_double), dimension(sd_d0),intent(out) :: sd
        real(kind=c_double), dimension(sh1_d0,sh1_d1,sh1_d2),intent(in) :: sh1
        real(kind=c_double), dimension(sh2_d0,sh2_d1,sh2_d2),intent(in) :: sh2
        real(kind=c_double), dimension(tapers_d0,tapers_d1),intent(in) :: tapers
        integer(kind=c_int), intent(in) :: lmax1
        integer(kind=c_int), intent(in) :: lmax2
        integer(kind=c_int), intent(in) :: lmaxt
        integer(kind=c_int), intent(in) :: k
        real(kind=c_double), dimension(taper_wt_d0),intent(in) :: taper_wt
        integer(kind=c_int), intent(in) :: csphase
        integer(kind=c_int), intent(in) :: norm
        integer(kind=c_int), intent(out) :: exitstatus
        call SHMultiTaperMaskCSE(mtse,sd,sh1,lmax1,sh2,lmax2,tapers,lmaxt,k,taper_wt&
                                     ,norm,csphase,exitstatus)
    end subroutine cbindSHMultiTaperMaskCSE

    subroutine cbindComputeDMap(Dij,dh_mask,n_dh,lmax,sampling,degrees,exitstatus&
                                   ,Dij_d0,Dij_d1,dh_mask_d0,dh_mask_d1,degrees_d0)  bind(c&
                                   , name="cbind_compute_d_map")
        use, intrinsic :: iso_c_binding
        use shtools, only: ComputeDMap
        implicit none
        integer(kind=c_int), intent(in) :: degrees_d0
        integer(kind=c_int), intent(in) :: dh_mask_d0
        integer(kind=c_int), intent(in) :: dh_mask_d1
        integer(kind=c_int), intent(in) :: Dij_d0
        integer(kind=c_int), intent(in) :: Dij_d1
        real(kind=c_double), dimension(Dij_d0,Dij_d1),intent(out) :: Dij
        integer(kind=c_int), dimension(dh_mask_d0,dh_mask_d1),intent(in) :: dh_mask
        integer(kind=c_int), intent(in) :: n_dh
        integer(kind=c_int), intent(in) :: lmax
        integer(kind=c_int), intent(in) :: sampling
        integer(kind=c_int), dimension(degrees_d0),intent(in) :: degrees
        integer(kind=c_int), intent(out) :: exitstatus
        call ComputeDMap(Dij,dh_mask,n_dh,lmax,sampling,degrees,exitstatus)
    end subroutine cbindComputeDMap

    subroutine cbindCurve2Mask(dhgrid,n,sampling,profile,nprofile,NP,extend,exitstatus&
                                     ,dhgrid_d0,dhgrid_d1,profile_d0,profile_d1)  bind(c&
                                     , name="cbind_curve2_mask")
        use, intrinsic :: iso_c_binding
        use shtools, only: Curve2Mask
        implicit none
        integer(kind=c_int), intent(in) :: profile_d0
        integer(kind=c_int), intent(in) :: profile_d1
        integer(kind=c_int), intent(in) :: dhgrid_d0
        integer(kind=c_int), intent(in) :: dhgrid_d1
        integer(kind=c_int), dimension(dhgrid_d0,dhgrid_d1),intent(out) :: dhgrid
        real(kind=c_double), dimension(profile_d0,profile_d1),intent(in) :: profile
        integer(kind=c_int), intent(in) :: n
        integer(kind=c_int), intent(in) :: sampling
        integer(kind=c_int), intent(in) :: nprofile
        integer(kind=c_int), intent(in) :: extend
        integer(kind=c_int), intent(out) :: exitstatus
        integer(kind=c_int) :: NP
        call Curve2Mask(dhgrid,n,sampling,profile,nprofile,NP,extend,exitstatus)
    end subroutine cbindCurve2Mask

    subroutine cbindSHBias(Shh,lwin,incspectra,ldata,outcspectra,save_cg,exitstatus&
                              ,Shh_d0,incspectra_d0,outcspectra_d0)  bind(c, name="cbind_sh_bias")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHBias
        implicit none
        integer(kind=c_int), intent(in) :: outcspectra_d0
        integer(kind=c_int), intent(in) :: incspectra_d0
        integer(kind=c_int), intent(in) :: Shh_d0
        real(kind=c_double), dimension(Shh_d0),intent(in) :: Shh
        real(kind=c_double), dimension(incspectra_d0),intent(in) :: incspectra
        real(kind=c_double), dimension(outcspectra_d0),intent(out) :: outcspectra
        integer(kind=c_int), intent(in) :: lwin
        integer(kind=c_int), intent(in) :: ldata
        integer(kind=c_int), intent(in) :: save_cg
        integer(kind=c_int), intent(out) :: exitstatus
        call SHBias(Shh,lwin,incspectra,ldata,outcspectra,save_cg,exitstatus)
    end subroutine cbindSHBias

    subroutine cbindSphericalCapCoef(coef,theta,lmax,exitstatus,coef_d0)  bind(c, name="cbind_spherical_cap_coef")
        use, intrinsic :: iso_c_binding
        use shtools, only: SphericalCapCoef
        implicit none
        integer(kind=c_int), intent(in) :: coef_d0
        real(kind=c_double), dimension(coef_d0),intent(out) :: coef
        real(kind=c_double), intent(in) :: theta
        integer(kind=c_int), intent(in) :: lmax
        integer(kind=c_int), intent(out) :: exitstatus
        call SphericalCapCoef(coef,theta,lmax,exitstatus)
    end subroutine cbindSphericalCapCoef

    subroutine cbindMakeGravGridDH(cilm,lmax,gm,r0,a,f,rad,theta,phi,total,n,sampling&
                                       ,lmax_calc,omega,normal_gravity,pot,extend&
                                       ,exitstatus,cilm_d,rad_d0,rad_d1,theta_d0,theta_d1&
                                       ,phi_d0,phi_d1,total_d0,total_d1,pot_d0,pot_d1)  bind(c&
                                       , name="cbind_make_grav_grid_dh")
        use, intrinsic :: iso_c_binding
        use shtools, only: MakeGravGridDH
        implicit none
        integer(kind=c_int), intent(in) :: pot_d0
        integer(kind=c_int), intent(in) :: pot_d1
        integer(kind=c_int), intent(in) :: total_d0
        integer(kind=c_int), intent(in) :: total_d1
        integer(kind=c_int), intent(in) :: phi_d0
        integer(kind=c_int), intent(in) :: phi_d1
        integer(kind=c_int), intent(in) :: theta_d0
        integer(kind=c_int), intent(in) :: theta_d1
        integer(kind=c_int), intent(in) :: rad_d0
        integer(kind=c_int), intent(in) :: rad_d1
        integer(kind=c_int), intent(in) :: cilm_d
        real(kind=c_double), dimension(2,cilm_d,cilm_d),intent(in) :: cilm
        real(kind=c_double), intent(in) :: gm
        real(kind=c_double), intent(in) :: r0
        real(kind=c_double), intent(in) :: a
        real(kind=c_double), intent(in) :: f
        real(kind=c_double), dimension(rad_d0,rad_d1),intent(out) :: rad
        real(kind=c_double), dimension(theta_d0,theta_d1),intent(out) :: theta
        real(kind=c_double), dimension(phi_d0,phi_d1),intent(out) :: phi
        real(kind=c_double), dimension(total_d0,total_d1),intent(out) :: total
        real(kind=c_double), intent(in) :: omega
        real(kind=c_double), dimension(pot_d0,pot_d1),intent(out) :: pot
        integer(kind=c_int), intent(in) :: lmax
        integer(kind=c_int), intent(out) :: n
        integer(kind=c_int), intent(in) :: sampling
        integer(kind=c_int), intent(in) :: lmax_calc
        integer(kind=c_int), intent(in) :: normal_gravity
        integer(kind=c_int), intent(in) :: extend
        integer(kind=c_int), intent(out) :: exitstatus
        call MakeGravGridDH(cilm,lmax,gm,r0,a,f,rad,theta,phi,total,n,sampling,lmax_calc&
                                ,omega,normal_gravity,pot,extend,exitstatus)
    end subroutine cbindMakeGravGridDH

    subroutine cbindMakeGravGradGridDH(cilm,lmax,gm,r0,a,f,vxx,vyy,vzz,vxy,vxz,vyz&
                                           ,n,sampling,lmax_calc,extend,exitstatus&
                                           ,cilm_d,vxx_d0,vxx_d1,vyy_d0,vyy_d1,vzz_d0&
                                           ,vzz_d1,vxy_d0,vxy_d1,vxz_d0,vxz_d1,vyz_d0&
                                           ,vyz_d1)  bind(c, name="cbind_make_grav_grad_grid_dh")
        use, intrinsic :: iso_c_binding
        use shtools, only: MakeGravGradGridDH
        implicit none
        integer(kind=c_int), intent(in) :: vyz_d0
        integer(kind=c_int), intent(in) :: vyz_d1
        integer(kind=c_int), intent(in) :: vxz_d0
        integer(kind=c_int), intent(in) :: vxz_d1
        integer(kind=c_int), intent(in) :: vxy_d0
        integer(kind=c_int), intent(in) :: vxy_d1
        integer(kind=c_int), intent(in) :: vzz_d0
        integer(kind=c_int), intent(in) :: vzz_d1
        integer(kind=c_int), intent(in) :: vyy_d0
        integer(kind=c_int), intent(in) :: vyy_d1
        integer(kind=c_int), intent(in) :: vxx_d0
        integer(kind=c_int), intent(in) :: vxx_d1
        integer(kind=c_int), intent(in) :: cilm_d
        real(kind=c_double), dimension(2,cilm_d,cilm_d),intent(in) :: cilm
        real(kind=c_double), intent(in) :: gm
        real(kind=c_double), intent(in) :: r0
        real(kind=c_double), intent(in) :: a
        real(kind=c_double), intent(in) :: f
        real(kind=c_double), dimension(vxx_d0,vxx_d1),intent(out) :: vxx
        real(kind=c_double), dimension(vyy_d0,vyy_d1),intent(out) :: vyy
        real(kind=c_double), dimension(vzz_d0,vzz_d1),intent(out) :: vzz
        real(kind=c_double), dimension(vxy_d0,vxy_d1),intent(out) :: vxy
        real(kind=c_double), dimension(vxz_d0,vxz_d1),intent(out) :: vxz
        real(kind=c_double), dimension(vyz_d0,vyz_d1),intent(out) :: vyz
        integer(kind=c_int), intent(in) :: lmax
        integer(kind=c_int), intent(out) :: n
        integer(kind=c_int), intent(in) :: sampling
        integer(kind=c_int), intent(in) :: lmax_calc
        integer(kind=c_int), intent(in) :: extend
        integer(kind=c_int), intent(out) :: exitstatus
        call MakeGravGradGridDH(cilm,lmax,gm,r0,a,f,vxx,vyy,vzz,vxy,vxz,vyz,n,sampling&
                                    ,lmax_calc,extend,exitstatus)
    end subroutine cbindMakeGravGradGridDH

    subroutine cbindMakeMagGradGridDH(cilm,lmax,r0,a,f,vxx,vyy,vzz,vxy,vxz,vyz,n,sampling&
                                          ,lmax_calc,extend,exitstatus,cilm_d,vxx_d0&
                                          ,vxx_d1,vyy_d0,vyy_d1,vzz_d0,vzz_d1,vxy_d0&
                                          ,vxy_d1,vxz_d0,vxz_d1,vyz_d0,vyz_d1)  bind(c&
                                          , name="cbind_make_mag_grad_grid_dh")
        use, intrinsic :: iso_c_binding
        use shtools, only: MakeMagGradGridDH
        implicit none
        integer(kind=c_int), intent(in) :: vyz_d0
        integer(kind=c_int), intent(in) :: vyz_d1
        integer(kind=c_int), intent(in) :: vxz_d0
        integer(kind=c_int), intent(in) :: vxz_d1
        integer(kind=c_int), intent(in) :: vxy_d0
        integer(kind=c_int), intent(in) :: vxy_d1
        integer(kind=c_int), intent(in) :: vzz_d0
        integer(kind=c_int), intent(in) :: vzz_d1
        integer(kind=c_int), intent(in) :: vyy_d0
        integer(kind=c_int), intent(in) :: vyy_d1
        integer(kind=c_int), intent(in) :: vxx_d0
        integer(kind=c_int), intent(in) :: vxx_d1
        integer(kind=c_int), intent(in) :: cilm_d
        real(kind=c_double), dimension(2,cilm_d,cilm_d),intent(in) :: cilm
        real(kind=c_double), intent(in) :: r0
        real(kind=c_double), intent(in) :: a
        real(kind=c_double), intent(in) :: f
        real(kind=c_double), dimension(vxx_d0,vxx_d1),intent(out) :: vxx
        real(kind=c_double), dimension(vyy_d0,vyy_d1),intent(out) :: vyy
        real(kind=c_double), dimension(vzz_d0,vzz_d1),intent(out) :: vzz
        real(kind=c_double), dimension(vxy_d0,vxy_d1),intent(out) :: vxy
        real(kind=c_double), dimension(vxz_d0,vxz_d1),intent(out) :: vxz
        real(kind=c_double), dimension(vyz_d0,vyz_d1),intent(out) :: vyz
        integer(kind=c_int), intent(in) :: lmax
        integer(kind=c_int), intent(out) :: n
        integer(kind=c_int), intent(in) :: sampling
        integer(kind=c_int), intent(in) :: lmax_calc
        integer(kind=c_int), intent(in) :: extend
        integer(kind=c_int), intent(out) :: exitstatus
        call MakeMagGradGridDH(cilm,lmax,r0,a,f,vxx,vyy,vzz,vxy,vxz,vyz,n,sampling&
                                   ,lmax_calc,extend,exitstatus)
    end subroutine cbindMakeMagGradGridDH

    subroutine cbindMakeGeoidGrid(geoid,cilm,lmax,r0pot,GM,PotRef,omega,r,gridtype&
                                       ,order,nlat,nlong,interval,lmax_calc,a,f,extend&
                                       ,exitstatus,geoid_d0,geoid_d1,cilm_d)  bind(c&
                                       , name="cbind_make_geoid_grid")
        use, intrinsic :: iso_c_binding
        use shtools, only: MakeGeoidGrid
        implicit none
        integer(kind=c_int), intent(in) :: cilm_d
        integer(kind=c_int), intent(in) :: geoid_d0
        integer(kind=c_int), intent(in) :: geoid_d1
        real(kind=c_double), dimension(geoid_d0,geoid_d1),intent(out) :: geoid
        real(kind=c_double), dimension(2,cilm_d,cilm_d),intent(in) :: cilm
        real(kind=c_double), intent(in) :: r0pot
        real(kind=c_double), intent(in) :: GM
        real(kind=c_double), intent(in) :: r
        real(kind=c_double), intent(in) :: PotRef
        real(kind=c_double), intent(in) :: omega
        integer(kind=c_int), intent(in) :: lmax
        integer(kind=c_int), intent(in) :: order
        integer(kind=c_int), intent(in) :: gridtype
        integer(kind=c_int), intent(in) :: lmax_calc
        integer(kind=c_int), intent(in) :: extend
        integer(kind=c_int), intent(out) :: nlat
        integer(kind=c_int), intent(out) :: nlong
        real(kind=c_double), intent(in) :: interval
        real(kind=c_double), intent(in) :: a
        real(kind=c_double), intent(in) :: f
        integer(kind=c_int), intent(out) :: exitstatus
        call MakeGeoidGrid(geoid,cilm,lmax,r0pot,GM,PotRef,omega,r,gridtype,order&
                                ,nlat,nlong,interval,lmax_calc,a,f,extend,exitstatus)
    end subroutine cbindMakeGeoidGrid

    subroutine cbindCilmPlus(cilm,gridin,lmax,nmax,mass,d,rho,gridtype,w,zero,plx&
                                 ,n,dref,exitstatus,gridin_d0,gridin_d1,w_d0,zero_d0&
                                 ,plx_d0,plx_d1,cilm_d)  bind(c, name="cbind_cilm_plus")
        use, intrinsic :: iso_c_binding
        use shtools, only: CilmPlus
        implicit none
        integer(kind=c_int), intent(in) :: cilm_d
        integer(kind=c_int), intent(in) :: plx_d0
        integer(kind=c_int), intent(in) :: plx_d1
        integer(kind=c_int), intent(in) :: zero_d0
        integer(kind=c_int), intent(in) :: w_d0
        integer(kind=c_int), intent(in) :: gridin_d0
        integer(kind=c_int), intent(in) :: gridin_d1
        real(kind=c_double), dimension(gridin_d0,gridin_d1),intent(in) :: gridin
        real(kind=c_double), intent(in) :: mass
        real(kind=c_double), intent(in) :: rho
        real(kind=c_double), dimension(w_d0),intent(in) :: w
        real(kind=c_double), dimension(zero_d0),intent(in) :: zero
        real(kind=c_double), dimension(plx_d0,plx_d1),intent(in) :: plx
        real(kind=c_double), intent(in) :: dref
        real(kind=c_double), dimension(2,cilm_d,cilm_d),intent(out) :: cilm
        real(kind=c_double), intent(out) :: d
        integer(kind=c_int), intent(in) :: lmax
        integer(kind=c_int), intent(in) :: nmax
        integer(kind=c_int), intent(in) :: gridtype
        integer(kind=c_int), intent(in) :: n
        integer(kind=c_int), intent(out) :: exitstatus
        call CilmPlus(cilm,gridin,lmax,nmax,mass,d,rho,gridtype,w,zero,plx,n,dref&
                          ,exitstatus)
    end subroutine cbindCilmPlus

    subroutine cbindCilmMinus(cilm,gridin,lmax,nmax,mass,d,rho,gridtype,w,zero,plx&
                                  ,n,dref,exitstatus,gridin_d0,gridin_d1,w_d0,zero_d0&
                                  ,plx_d0,plx_d1,cilm_d)  bind(c, name="cbind_cilm_minus")
        use, intrinsic :: iso_c_binding
        use shtools, only: CilmMinus
        implicit none
        integer(kind=c_int), intent(in) :: cilm_d
        integer(kind=c_int), intent(in) :: plx_d0
        integer(kind=c_int), intent(in) :: plx_d1
        integer(kind=c_int), intent(in) :: zero_d0
        integer(kind=c_int), intent(in) :: w_d0
        integer(kind=c_int), intent(in) :: gridin_d0
        integer(kind=c_int), intent(in) :: gridin_d1
        real(kind=c_double), dimension(gridin_d0,gridin_d1),intent(in) :: gridin
        real(kind=c_double), intent(in) :: mass
        real(kind=c_double), intent(in) :: rho
        real(kind=c_double), dimension(w_d0),intent(in) :: w
        real(kind=c_double), dimension(zero_d0),intent(in) :: zero
        real(kind=c_double), dimension(plx_d0,plx_d1),intent(in) :: plx
        real(kind=c_double), intent(in) :: dref
        real(kind=c_double), dimension(2,cilm_d,cilm_d),intent(out) :: cilm
        real(kind=c_double), intent(out) :: d
        integer(kind=c_int), intent(in) :: lmax
        integer(kind=c_int), intent(in) :: nmax
        integer(kind=c_int), intent(in) :: gridtype
        integer(kind=c_int), intent(in) :: n
        integer(kind=c_int), intent(out) :: exitstatus
        call CilmMinus(cilm,gridin,lmax,nmax,mass,d,rho,gridtype,w,zero,plx,n,dref&
                           ,exitstatus)
    end subroutine cbindCilmMinus

    subroutine cbindCilmPlusRhoH(cilm,gridin,lmax,nmax,mass,d,rho,gridtype,w,zero&
                                     ,plx,n,dref,exitstatus,gridin_d0,gridin_d1,rho_d0&
                                     ,rho_d1,w_d0,zero_d0,plx_d0,plx_d1,cilm_d)  bind(c&
                                     , name="cbind_cilm_plus_rho_h")
        use, intrinsic :: iso_c_binding
        use shtools, only: CilmPlusRhoH
        implicit none
        integer(kind=c_int), intent(in) :: cilm_d
        integer(kind=c_int), intent(in) :: plx_d0
        integer(kind=c_int), intent(in) :: plx_d1
        integer(kind=c_int), intent(in) :: zero_d0
        integer(kind=c_int), intent(in) :: w_d0
        integer(kind=c_int), intent(in) :: rho_d0
        integer(kind=c_int), intent(in) :: rho_d1
        integer(kind=c_int), intent(in) :: gridin_d0
        integer(kind=c_int), intent(in) :: gridin_d1
        real(kind=c_double), dimension(gridin_d0,gridin_d1),intent(in) :: gridin
        real(kind=c_double), intent(in) :: mass
        real(kind=c_double), dimension(rho_d0,rho_d1),intent(in) :: rho
        real(kind=c_double), dimension(w_d0),intent(in) :: w
        real(kind=c_double), dimension(zero_d0),intent(in) :: zero
        real(kind=c_double), dimension(plx_d0,plx_d1),intent(in) :: plx
        real(kind=c_double), intent(in) :: dref
        real(kind=c_double), dimension(2,cilm_d,cilm_d),intent(out) :: cilm
        real(kind=c_double), intent(out) :: d
        integer(kind=c_int), intent(in) :: lmax
        integer(kind=c_int), intent(in) :: nmax
        integer(kind=c_int), intent(in) :: gridtype
        integer(kind=c_int), intent(in) :: n
        integer(kind=c_int), intent(out) :: exitstatus
        call CilmPlusRhoH(cilm,gridin,lmax,nmax,mass,d,rho,gridtype,w,zero,plx,n,dref&
                              ,exitstatus)
    end subroutine cbindCilmPlusRhoH

    subroutine cbindCilmMinusRhoH(cilm,gridin,lmax,nmax,mass,d,rho,gridtype,w,zero&
                                      ,plx,n,dref,exitstatus,gridin_d0,gridin_d1,rho_d0&
                                      ,rho_d1,w_d0,zero_d0,plx_d0,plx_d1,cilm_d)  bind(c&
                                      , name="cbind_cilm_minus_rho_h")
        use, intrinsic :: iso_c_binding
        use shtools, only: CilmMinusRhoH
        implicit none
        integer(kind=c_int), intent(in) :: cilm_d
        integer(kind=c_int), intent(in) :: plx_d0
        integer(kind=c_int), intent(in) :: plx_d1
        integer(kind=c_int), intent(in) :: zero_d0
        integer(kind=c_int), intent(in) :: w_d0
        integer(kind=c_int), intent(in) :: rho_d0
        integer(kind=c_int), intent(in) :: rho_d1
        integer(kind=c_int), intent(in) :: gridin_d0
        integer(kind=c_int), intent(in) :: gridin_d1
        real(kind=c_double), dimension(gridin_d0,gridin_d1),intent(in) :: gridin
        real(kind=c_double), intent(in) :: mass
        real(kind=c_double), dimension(rho_d0,rho_d1),intent(in) :: rho
        real(kind=c_double), dimension(w_d0),intent(in) :: w
        real(kind=c_double), dimension(zero_d0),intent(in) :: zero
        real(kind=c_double), dimension(plx_d0,plx_d1),intent(in) :: plx
        real(kind=c_double), intent(in) :: dref
        real(kind=c_double), dimension(2,cilm_d,cilm_d),intent(out) :: cilm
        real(kind=c_double), intent(out) :: d
        integer(kind=c_int), intent(in) :: lmax
        integer(kind=c_int), intent(in) :: nmax
        integer(kind=c_int), intent(in) :: gridtype
        integer(kind=c_int), intent(in) :: n
        integer(kind=c_int), intent(out) :: exitstatus
        call CilmMinusRhoH(cilm,gridin,lmax,nmax,mass,d,rho,gridtype,w,zero,plx,n&
                               ,dref,exitstatus)
    end subroutine cbindCilmMinusRhoH

    subroutine cbindBAtoHilm(cilm,ba,gridglq,lmax,nmax,mass,r0,rho,gridtype,w,plx&
                                 ,zero,filter_type,filter_deg,lmax_calc,exitstatus&
                                 ,cilm_d,ba_d0,ba_d1,ba_d2,gridglq_d0,gridglq_d1,plx_d0&
                                 ,plx_d1,zero_d0,w_d0)  bind(c, name="cbind_b_ato_hilm")
        use, intrinsic :: iso_c_binding
        use shtools, only: BAtoHilm
        implicit none
        integer(kind=c_int), intent(in) :: w_d0
        integer(kind=c_int), intent(in) :: zero_d0
        integer(kind=c_int), intent(in) :: plx_d0
        integer(kind=c_int), intent(in) :: plx_d1
        integer(kind=c_int), intent(in) :: gridglq_d0
        integer(kind=c_int), intent(in) :: gridglq_d1
        integer(kind=c_int), intent(in) :: ba_d0
        integer(kind=c_int), intent(in) :: ba_d1
        integer(kind=c_int), intent(in) :: ba_d2
        integer(kind=c_int), intent(in) :: cilm_d
        real(kind=c_double), dimension(2,cilm_d,cilm_d),intent(out) :: cilm
        real(kind=c_double), dimension(ba_d0,ba_d1,ba_d2),intent(in) :: ba
        real(kind=c_double), dimension(gridglq_d0,gridglq_d1),intent(in) :: gridglq
        real(kind=c_double), intent(in) :: mass
        real(kind=c_double), intent(in) :: r0
        real(kind=c_double), intent(in) :: rho
        real(kind=c_double), dimension(plx_d0,plx_d1),intent(in) :: plx
        real(kind=c_double), dimension(zero_d0),intent(in) :: zero
        real(kind=c_double), dimension(w_d0),intent(in) :: w
        integer(kind=c_int), intent(in) :: lmax
        integer(kind=c_int), intent(in) :: nmax
        integer(kind=c_int), intent(in) :: gridtype
        integer(kind=c_int), intent(in) :: filter_type
        integer(kind=c_int), intent(in) :: filter_deg
        integer(kind=c_int), intent(in) :: lmax_calc
        integer(kind=c_int), intent(out) :: exitstatus
        call BAtoHilm(cilm,ba,gridglq,lmax,nmax,mass,r0,rho,gridtype,w,plx,zero,filter_type&
                          ,filter_deg,lmax_calc,exitstatus)
    end subroutine cbindBAtoHilm

    subroutine cbindBAtoHilmRhoH(cilm,ba,gridglq,lmax,nmax,mass,r0,rho,gridtype,w&
                                     ,plx,zero,filter_type,filter_deg,lmax_calc,exitstatus&
                                     ,cilm_d,ba_d0,ba_d1,ba_d2,gridglq_d0,gridglq_d1&
                                     ,rho_d0,rho_d1,plx_d0,plx_d1,zero_d0,w_d0)  bind(c&
                                     , name="cbind_b_ato_hilm_rho_h")
        use, intrinsic :: iso_c_binding
        use shtools, only: BAtoHilmRhoH
        implicit none
        integer(kind=c_int), intent(in) :: w_d0
        integer(kind=c_int), intent(in) :: zero_d0
        integer(kind=c_int), intent(in) :: plx_d0
        integer(kind=c_int), intent(in) :: plx_d1
        integer(kind=c_int), intent(in) :: rho_d0
        integer(kind=c_int), intent(in) :: rho_d1
        integer(kind=c_int), intent(in) :: gridglq_d0
        integer(kind=c_int), intent(in) :: gridglq_d1
        integer(kind=c_int), intent(in) :: ba_d0
        integer(kind=c_int), intent(in) :: ba_d1
        integer(kind=c_int), intent(in) :: ba_d2
        integer(kind=c_int), intent(in) :: cilm_d
        real(kind=c_double), dimension(2,cilm_d,cilm_d),intent(out) :: cilm
        real(kind=c_double), dimension(ba_d0,ba_d1,ba_d2),intent(in) :: ba
        real(kind=c_double), dimension(gridglq_d0,gridglq_d1),intent(in) :: gridglq
        real(kind=c_double), intent(in) :: mass
        real(kind=c_double), intent(in) :: r0
        real(kind=c_double), dimension(rho_d0,rho_d1),intent(in) :: rho
        real(kind=c_double), dimension(plx_d0,plx_d1),intent(in) :: plx
        real(kind=c_double), dimension(zero_d0),intent(in) :: zero
        real(kind=c_double), dimension(w_d0),intent(in) :: w
        integer(kind=c_int), intent(in) :: lmax
        integer(kind=c_int), intent(in) :: nmax
        integer(kind=c_int), intent(in) :: gridtype
        integer(kind=c_int), intent(in) :: filter_type
        integer(kind=c_int), intent(in) :: filter_deg
        integer(kind=c_int), intent(in) :: lmax_calc
        integer(kind=c_int), intent(out) :: exitstatus
        call BAtoHilmRhoH(cilm,ba,gridglq,lmax,nmax,mass,r0,rho,gridtype,w,plx,zero&
                              ,filter_type,filter_deg,lmax_calc,exitstatus)
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

    subroutine cbindMakeMagGridDH(cilm,lmax,r0,a,f,rad_grid,theta_grid,phi_grid,total_grid&
                                      ,n,sampling,lmax_calc,pot_grid,extend,exitstatus&
                                      ,cilm_d,rad_grid_d0,rad_grid_d1,theta_grid_d0&
                                      ,theta_grid_d1,phi_grid_d0,phi_grid_d1,total_grid_d0&
                                      ,total_grid_d1,pot_grid_d0,pot_grid_d1)  bind(c&
                                      , name="cbind_make_mag_grid_dh")
        use, intrinsic :: iso_c_binding
        use shtools, only: MakeMagGridDH
        implicit none
        integer(kind=c_int), intent(in) :: pot_grid_d0
        integer(kind=c_int), intent(in) :: pot_grid_d1
        integer(kind=c_int), intent(in) :: total_grid_d0
        integer(kind=c_int), intent(in) :: total_grid_d1
        integer(kind=c_int), intent(in) :: phi_grid_d0
        integer(kind=c_int), intent(in) :: phi_grid_d1
        integer(kind=c_int), intent(in) :: theta_grid_d0
        integer(kind=c_int), intent(in) :: theta_grid_d1
        integer(kind=c_int), intent(in) :: rad_grid_d0
        integer(kind=c_int), intent(in) :: rad_grid_d1
        integer(kind=c_int), intent(in) :: cilm_d
        real(kind=c_double), dimension(2,cilm_d,cilm_d),intent(in) :: cilm
        real(kind=c_double), intent(in) :: r0
        real(kind=c_double), intent(in) :: a
        real(kind=c_double), intent(in) :: f
        real(kind=c_double), dimension(rad_grid_d0,rad_grid_d1),intent(out) :: rad_grid
        real(kind=c_double), dimension(theta_grid_d0,theta_grid_d1),intent(out) :: theta_grid
        real(kind=c_double), dimension(phi_grid_d0,phi_grid_d1),intent(out) :: phi_grid
        real(kind=c_double), dimension(total_grid_d0,total_grid_d1),intent(out) :: total_grid
        real(kind=c_double), dimension(pot_grid_d0,pot_grid_d1),intent(out) :: pot_grid
        integer(kind=c_int), intent(in) :: lmax
        integer(kind=c_int), intent(out) :: n
        integer(kind=c_int), intent(in) :: sampling
        integer(kind=c_int), intent(in) :: lmax_calc
        integer(kind=c_int), intent(in) :: extend
        integer(kind=c_int), intent(out) :: exitstatus
        call MakeMagGridDH(cilm,lmax,r0,a,f,rad_grid,theta_grid,phi_grid,total_grid&
                               ,n,sampling,lmax_calc,pot_grid,extend,exitstatus)
    end subroutine cbindMakeMagGridDH

    subroutine cbindSHMagPowerSpectrum(c,a,r,lmax,spectra,exitstatus,c_d0,c_d1,c_d2&
                                        ,spectra_d0)  bind(c, name="cbind_sh_mag_power_spectrum")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHMagPowerSpectrum
        implicit none
        integer(kind=c_int), intent(in) :: spectra_d0
        integer(kind=c_int), intent(in) :: c_d0
        integer(kind=c_int), intent(in) :: c_d1
        integer(kind=c_int), intent(in) :: c_d2
        real(kind=c_double), dimension(c_d0,c_d1,c_d2),intent(in) :: c
        real(kind=c_double), intent(in) :: a
        real(kind=c_double), intent(in) :: r
        integer(kind=c_int), intent(in) :: lmax
        real(kind=c_double), dimension(spectra_d0),intent(out) :: spectra
        integer(kind=c_int), intent(out) :: exitstatus
        call SHMagPowerSpectrum(c,a,r,lmax,spectra,exitstatus)
    end subroutine cbindSHMagPowerSpectrum

    function cbindSHMagPowerL(c,a,r,l,c_d0,c_d1,c_d2)  bind(c, name="cbind_sh_mag_power_l")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHMagPowerL
        implicit none
        integer(kind=c_int), intent(in) :: c_d0
        integer(kind=c_int), intent(in) :: c_d1
        integer(kind=c_int), intent(in) :: c_d2
        real(kind=c_double) :: cbindSHMagPowerL
        real(kind=c_double), dimension(c_d0,c_d1,c_d2),intent(in) :: c
        real(kind=c_double), intent(in) :: a
        real(kind=c_double), intent(in) :: r
        integer(kind=c_int), intent(in) :: l
        cbindSHMagPowerL=SHMagPowerL(c,a,r,l)
    end function cbindSHMagPowerL

    subroutine cbindMakeCircleCoord(coord,lat,lon,theta0,cinterval,cnum,exitstatus&
                                         ,coord_d0,coord_d1)  bind(c, name="cbind_make_circle_coord")
        use, intrinsic :: iso_c_binding
        use shtools, only: MakeCircleCoord
        implicit none
        integer(kind=c_int), intent(in) :: coord_d0
        integer(kind=c_int), intent(in) :: coord_d1
        real(kind=c_double), intent(in) :: lat
        real(kind=c_double), intent(in) :: lon
        real(kind=c_double), intent(in) :: theta0
        real(kind=c_double), dimension(coord_d0,coord_d1),intent(out) :: coord
        real(kind=c_double), intent(in) :: cinterval
        integer(kind=c_int), intent(out) :: cnum
        integer(kind=c_int), intent(out) :: exitstatus
        call MakeCircleCoord(coord,lat,lon,theta0,cinterval,cnum,exitstatus)
    end subroutine cbindMakeCircleCoord

    subroutine cbindMakeEllipseCoord(coord,lat,lon,dec,A_theta,B_theta,cinterval,cnum&
                                          ,exitstatus,coord_d0,coord_d1)  bind(c, name="cbind_make_ellipse_coord")
        use, intrinsic :: iso_c_binding
        use shtools, only: MakeEllipseCoord
        implicit none
        integer(kind=c_int), intent(in) :: coord_d0
        integer(kind=c_int), intent(in) :: coord_d1
        real(kind=c_double), intent(in) :: lat
        real(kind=c_double), intent(in) :: lon
        real(kind=c_double), intent(in) :: A_theta
        real(kind=c_double), intent(in) :: B_theta
        real(kind=c_double), intent(in) :: dec
        real(kind=c_double), dimension(coord_d0,coord_d1),intent(out) :: coord
        real(kind=c_double), intent(in) :: cinterval
        integer(kind=c_int), intent(out) :: cnum
        integer(kind=c_int), intent(out) :: exitstatus
        call MakeEllipseCoord(coord,lat,lon,dec,A_theta,B_theta,cinterval,cnum,exitstatus)
    end subroutine cbindMakeEllipseCoord

    subroutine cbindWigner3j(w3j,jmin,jmax,j2,j3,m1,m2,m3,exitstatus,w3j_d0)  bind(c&
                                , name="cbind_wigner3j")
        use, intrinsic :: iso_c_binding
        use shtools, only: Wigner3j
        implicit none
        integer(kind=c_int), intent(in) :: w3j_d0
        integer(kind=c_int), intent(in) :: j2
        integer(kind=c_int), intent(in) :: j3
        integer(kind=c_int), intent(in) :: m1
        integer(kind=c_int), intent(in) :: m2
        integer(kind=c_int), intent(in) :: m3
        integer(kind=c_int), intent(out) :: jmin
        integer(kind=c_int), intent(out) :: jmax
        real(kind=c_double), dimension(w3j_d0),intent(out) :: w3j
        integer(kind=c_int), intent(out) :: exitstatus
        call Wigner3j(w3j,jmin,jmax,j2,j3,m1,m2,m3,exitstatus)
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

    subroutine cbindPreGLQ(x1,x2,n,zero,w,exitstatus,zero_d0,w_d0)  bind(c, name="cbind_pre_glq")
        use, intrinsic :: iso_c_binding
        use shtools, only: PreGLQ
        implicit none
        integer(kind=c_int), intent(in) :: w_d0
        integer(kind=c_int), intent(in) :: zero_d0
        real(kind=c_double), intent(in) :: x1
        real(kind=c_double), intent(in) :: x2
        real(kind=c_double), dimension(zero_d0),intent(out) :: zero
        real(kind=c_double), dimension(w_d0),intent(out) :: w
        integer(kind=c_int), intent(in) :: n
        integer(kind=c_int), intent(out) :: exitstatus
        call PreGLQ(x1,x2,n,zero,w,exitstatus)
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

    subroutine cbindDHaj(n,aj,extend,exitstatus,aj_d0)  bind(c, name="cbind_d_haj")
        use, intrinsic :: iso_c_binding
        use shtools, only: DHaj
        implicit none
        integer(kind=c_int), intent(in) :: aj_d0
        integer(kind=c_int), intent(in) :: n
        real(kind=c_double), dimension(aj_d0),intent(out) :: aj
        integer(kind=c_int), intent(in) :: extend
        integer(kind=c_int), intent(out) :: exitstatus
        call DHaj(n,aj,extend,exitstatus)
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

    subroutine cbindEigValVecSym(ain,n,eig,evec,ul,K,exitstatus,ain_d0,ain_d1,eig_d0&
                                    ,evec_d0,evec_d1)  bind(c, name="cbind_eig_val_vec_sym")
        use, intrinsic :: iso_c_binding
        use shtools, only: EigValVecSym
        implicit none
        integer(kind=c_int), intent(in) :: evec_d0
        integer(kind=c_int), intent(in) :: evec_d1
        integer(kind=c_int), intent(in) :: eig_d0
        integer(kind=c_int), intent(in) :: ain_d0
        integer(kind=c_int), intent(in) :: ain_d1
        real(kind=c_double), dimension(ain_d0,ain_d1),intent(in) :: ain
        integer(kind=c_int), intent(in) :: n
        real(kind=c_double), dimension(eig_d0),intent(out) :: eig
        real(kind=c_double), dimension(evec_d0,evec_d1),intent(out) :: evec
        character(kind=c_char), intent(in) :: ul
        integer(kind=c_int), intent(in) :: K
        integer(kind=c_int), intent(out) :: exitstatus
        call EigValVecSym(ain,n,eig,evec,ul,K,exitstatus)
    end subroutine cbindEigValVecSym

    subroutine cbindEigValVecSymTri(ain,n,eig,evec,ul,exitstatus,ain_d0,ain_d1,eig_d0&
                                       ,evec_d0,evec_d1)  bind(c, name="cbind_eig_val_vec_sym_tri")
        use, intrinsic :: iso_c_binding
        use shtools, only: EigValVecSymTri
        implicit none
        integer(kind=c_int), intent(in) :: evec_d0
        integer(kind=c_int), intent(in) :: evec_d1
        integer(kind=c_int), intent(in) :: eig_d0
        integer(kind=c_int), intent(in) :: ain_d0
        integer(kind=c_int), intent(in) :: ain_d1
        real(kind=c_double), dimension(ain_d0,ain_d1),intent(in) :: ain
        integer(kind=c_int), intent(in) :: n
        real(kind=c_double), dimension(eig_d0),intent(out) :: eig
        real(kind=c_double), dimension(evec_d0,evec_d1),intent(out) :: evec
        character(kind=c_char), intent(in) :: ul
        integer(kind=c_int), intent(out) :: exitstatus
        call EigValVecSymTri(ain,n,eig,evec,ul,exitstatus)
    end subroutine cbindEigValVecSymTri

    subroutine cbindEigValSym(ain,n,eval,ul,ain_d0,ain_d1,eval_d0)  bind(c, name="cbind_eig_val_sym")
        use, intrinsic :: iso_c_binding
        use shtools, only: EigValSym
        implicit none
        integer(kind=c_int), intent(in) :: eval_d0
        integer(kind=c_int), intent(in) :: ain_d0
        integer(kind=c_int), intent(in) :: ain_d1
        real(kind=c_double), dimension(ain_d0,ain_d1),intent(in) :: ain
        integer(kind=c_int), intent(in) :: n
        real(kind=c_double), dimension(eval_d0),intent(out) :: eval
        character(kind=c_char), intent(in) :: ul
        call EigValSym(ain,n,eval,ul)
    end subroutine cbindEigValSym

    subroutine cbindSHRotateTapers(tapersrot,tapers,taper_order,lmax,nrot,x,dj,exitstatus&
                                            ,tapers_d0,tapers_d1,x_d0,dj_d0,dj_d1&
                                            ,dj_d2,tapersrot_d0,tapersrot_d1,taper_order_d0)  bind(c&
                                            , name="cbind_sh_rotate_tapers")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHRotateTapers
        implicit none
        integer(kind=c_int), intent(in) :: taper_order_d0
        integer(kind=c_int), intent(in) :: tapersrot_d0
        integer(kind=c_int), intent(in) :: tapersrot_d1
        integer(kind=c_int), intent(in) :: dj_d0
        integer(kind=c_int), intent(in) :: dj_d1
        integer(kind=c_int), intent(in) :: dj_d2
        integer(kind=c_int), intent(in) :: x_d0
        integer(kind=c_int), intent(in) :: tapers_d0
        integer(kind=c_int), intent(in) :: tapers_d1
        real(kind=c_double), dimension(tapers_d0,tapers_d1),intent(in) :: tapers
        real(kind=c_double), dimension(x_d0),intent(in) :: x
        real(kind=c_double), dimension(dj_d0,dj_d1,dj_d2),intent(in) :: dj
        real(kind=c_double), dimension(tapersrot_d0,tapersrot_d1),intent(out) :: tapersrot
        integer(kind=c_int), dimension(taper_order_d0),intent(in) :: taper_order
        integer(kind=c_int), intent(in) :: lmax
        integer(kind=c_int), intent(in) :: nrot
        integer(kind=c_int), intent(out) :: exitstatus
        call SHRotateTapers(tapersrot,tapers,taper_order,lmax,nrot,x,dj,exitstatus)
    end subroutine cbindSHRotateTapers

    subroutine cbindSlepianCoeffs(falpha,galpha,flm,lmax,nmax,exitstatus,falpha_d0&
                                        ,galpha_d0,galpha_d1,flm_d0,flm_d1,flm_d2)  bind(c&
                                        , name="cbind_slepian_coeffs")
        use, intrinsic :: iso_c_binding
        use shtools, only: SlepianCoeffs
        implicit none
        integer(kind=c_int), intent(in) :: flm_d0
        integer(kind=c_int), intent(in) :: flm_d1
        integer(kind=c_int), intent(in) :: flm_d2
        integer(kind=c_int), intent(in) :: galpha_d0
        integer(kind=c_int), intent(in) :: galpha_d1
        integer(kind=c_int), intent(in) :: falpha_d0
        real(kind=c_double), dimension(falpha_d0),intent(out) :: falpha
        real(kind=c_double), dimension(galpha_d0,galpha_d1),intent(in) :: galpha
        real(kind=c_double), dimension(flm_d0,flm_d1,flm_d2),intent(in) :: flm
        integer(kind=c_int), intent(in) :: lmax
        integer(kind=c_int), intent(in) :: nmax
        integer(kind=c_int), intent(out) :: exitstatus
        call SlepianCoeffs(falpha,galpha,flm,lmax,nmax,exitstatus)
    end subroutine cbindSlepianCoeffs

    subroutine cbindSlepianCoeffsToSH(flm,falpha,galpha,lmax,nmax,exitstatus,flm_d0&
                                         ,flm_d1,flm_d2,falpha_d0,galpha_d0,galpha_d1)  bind(c&
                                         , name="cbind_slepian_coeffs_to_sh")
        use, intrinsic :: iso_c_binding
        use shtools, only: SlepianCoeffsToSH
        implicit none
        integer(kind=c_int), intent(in) :: galpha_d0
        integer(kind=c_int), intent(in) :: galpha_d1
        integer(kind=c_int), intent(in) :: falpha_d0
        integer(kind=c_int), intent(in) :: flm_d0
        integer(kind=c_int), intent(in) :: flm_d1
        integer(kind=c_int), intent(in) :: flm_d2
        real(kind=c_double), dimension(flm_d0,flm_d1,flm_d2),intent(out) :: flm
        real(kind=c_double), dimension(falpha_d0),intent(in) :: falpha
        real(kind=c_double), dimension(galpha_d0,galpha_d1),intent(in) :: galpha
        integer(kind=c_int), intent(in) :: lmax
        integer(kind=c_int), intent(in) :: nmax
        integer(kind=c_int), intent(out) :: exitstatus
        call SlepianCoeffsToSH(flm,falpha,galpha,lmax,nmax,exitstatus)
    end subroutine cbindSlepianCoeffsToSH

    subroutine cbindSHSCouplingMatrix(kij,galpha,lmax,nmax,exitstatus,kij_d0,kij_d1&
                                         ,galpha_d0,galpha_d1)  bind(c, name="cbind_shs_coupling_matrix")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHSCouplingMatrix
        implicit none
        integer(kind=c_int), intent(in) :: galpha_d0
        integer(kind=c_int), intent(in) :: galpha_d1
        integer(kind=c_int), intent(in) :: kij_d0
        integer(kind=c_int), intent(in) :: kij_d1
        real(kind=c_double), dimension(kij_d0,kij_d1),intent(out) :: kij
        real(kind=c_double), dimension(galpha_d0,galpha_d1),intent(in) :: galpha
        integer(kind=c_int), intent(in) :: lmax
        integer(kind=c_int), intent(in) :: nmax
        integer(kind=c_int), intent(out) :: exitstatus
        call SHSCouplingMatrix(kij,galpha,lmax,nmax,exitstatus)
    end subroutine cbindSHSCouplingMatrix

    subroutine cbindSHSlepianVar(l,galpha,galpha_order,lmax,kmax,Sff,variance,exitstatus&
                                  ,galpha_d0,galpha_d1,Sff_d0,galpha_order_d0)  bind(c&
                                  , name="cbind_sh_slepian_var")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHSlepianVar
        implicit none
        integer(kind=c_int), intent(in) :: galpha_order_d0
        integer(kind=c_int), intent(in) :: Sff_d0
        integer(kind=c_int), intent(in) :: galpha_d0
        integer(kind=c_int), intent(in) :: galpha_d1
        real(kind=c_double), dimension(galpha_d0,galpha_d1),intent(in) :: galpha
        real(kind=c_double), dimension(Sff_d0),intent(in) :: Sff
        real(kind=c_double), intent(out) :: variance
        integer(kind=c_int), intent(in) :: l
        integer(kind=c_int), intent(in) :: lmax
        integer(kind=c_int), intent(in) :: kmax
        integer(kind=c_int), dimension(galpha_order_d0),intent(in) :: galpha_order
        integer(kind=c_int), intent(out) :: exitstatus
        call SHSlepianVar(l,galpha,galpha_order,lmax,kmax,Sff,variance,exitstatus)
    end subroutine cbindSHSlepianVar

    subroutine cbindSHSCouplingMatrixCap(kij,galpha,galpha_order,lmax,nmax,exitstatus&
                                            ,kij_d0,kij_d1,galpha_d0,galpha_d1,galpha_order_d0)  bind(c&
                                            , name="cbind_shs_coupling_matrix_cap")
        use, intrinsic :: iso_c_binding
        use shtools, only: SHSCouplingMatrixCap
        implicit none
        integer(kind=c_int), intent(in) :: galpha_order_d0
        integer(kind=c_int), intent(in) :: galpha_d0
        integer(kind=c_int), intent(in) :: galpha_d1
        integer(kind=c_int), intent(in) :: kij_d0
        integer(kind=c_int), intent(in) :: kij_d1
        real(kind=c_double), dimension(kij_d0,kij_d1),intent(out) :: kij
        real(kind=c_double), dimension(galpha_d0,galpha_d1),intent(in) :: galpha
        integer(kind=c_int), dimension(galpha_order_d0),intent(in) :: galpha_order
        integer(kind=c_int), intent(in) :: lmax
        integer(kind=c_int), intent(in) :: nmax
        integer(kind=c_int), intent(out) :: exitstatus
        call SHSCouplingMatrixCap(kij,galpha,galpha_order,lmax,nmax,exitstatus)
    end subroutine cbindSHSCouplingMatrixCap

