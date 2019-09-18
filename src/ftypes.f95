module ftypes
!------------------------------------------------------------------------------
!
! Define data types for Fortran 95 reals and integers.
!
!------------------------------------------------------------------------------
    implicit none

    integer, parameter :: sp = selected_real_kind(p=6)
    integer, parameter :: dp = selected_real_kind(p=15)
    integer, parameter :: qp = selected_real_kind(p=33)
    integer, parameter :: int1 = selected_int_kind(2)
    integer, parameter :: int4 = selected_int_kind(9)
    integer, parameter :: int8 = selected_int_kind(18)

end module ftypes