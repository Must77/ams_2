! Fortran 绑定模块：对应 ams_c_api.h。
! 本求解器构建用 32 位 PetscInt（integer(c_int)）与双精度 PetscReal
! （real(c_double)）；若调用方整型为 64 位需先转换。
! 字符串参数必须以 c_null_char 结尾，如:
!   ierr = ams_petsc_init('-em_outer_rtol 1e-7'//c_null_char)
module ams_solver
    use iso_c_binding, only: c_int, c_double, c_char
    implicit none

    interface

        integer(c_int) function ams_petsc_init(options_string) &
            bind(C, name='ams_petsc_init')
            import :: c_int, c_char
            character(kind=c_char), dimension(*), intent(in) :: options_string
        end function ams_petsc_init

        integer(c_int) function ams_solve_upper( &
            n_edges, n_nodes, n_upper_row_ptr, n_upper_values, &
            edge_nodes_1based, node_coords, upper_row_ptr, upper_col_idx, &
            rhs_real, rhs_imag, upper_real, upper_imag, &
            n_result, out_real, out_imag) bind(C, name='ams_solve_upper')
            import :: c_int, c_double
            integer(c_int), value :: n_edges, n_nodes
            integer(c_int), value :: n_upper_row_ptr, n_upper_values
            integer(c_int), value :: n_result
            integer(c_int), dimension(*), intent(in) :: edge_nodes_1based
            real(c_double), dimension(*), intent(in) :: node_coords
            integer(c_int), dimension(*), intent(in) :: upper_row_ptr
            integer(c_int), dimension(*), intent(in) :: upper_col_idx
            real(c_double), dimension(*), intent(in) :: rhs_real, rhs_imag
            real(c_double), dimension(*), intent(in) :: upper_real, upper_imag
            real(c_double), dimension(*), intent(out) :: out_real, out_imag
        end function ams_solve_upper

        integer(c_int) function ams_petsc_finalize() &
            bind(C, name='ams_petsc_finalize')
            import :: c_int
        end function ams_petsc_finalize

    end interface
end module ams_solver
