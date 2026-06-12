! Fortran 绑定数值对照：4 边链小算例走完整 init/solve/finalize 路径。
! 系统与 adapter_self_check.cpp 同拓扑，但虚部数据与实部不成比例
! （成比例时解的虚部恒为 0，测不出虚部符号错误）。
! 期望解由 numpy 对 [Kr -Ki; -Ki -Kr](xr,xi)=(Sr,-Si) 直接求解得到
! （cond≈14.6），与 upper_example_check.cpp 的残差定义一致。
program fortran_example_check
    use iso_c_binding, only: c_int, c_double, c_null_char
    use ams_solver
    implicit none

    integer(c_int), parameter :: n_edges = 4, n_nodes = 5
    integer(c_int), parameter :: n_rp = 5, n_val = 7

    integer(c_int) :: edge_nodes(2*n_edges), rp(n_rp), ci(n_val)
    real(c_double) :: coords(3*n_nodes), vre(n_val), vim(n_val)
    real(c_double) :: sr(n_edges), si(n_edges)
    real(c_double) :: xr(n_edges), xi(n_edges)
    real(c_double) :: er(n_edges), ei(n_edges)
    real(c_double) :: diff2, ref2, rel
    integer(c_int) :: ierr
    integer :: k

    edge_nodes = [1, 2, 2, 3, 3, 4, 4, 5]
    coords = real([0,0,0, 1,0,0, 1,1,0, 0,1,0, 0,0,1], c_double)
    rp = [1, 3, 5, 7, 8]
    ci = [1, 3, 2, 4, 3, 4, 4]
    vre = real([10, 20, 30, 40, 50, 60, 70], c_double)
    vim = real([2, -1, 4, 1, -3, 2, 5], c_double)
    sr = real([100, 200, 300, 400], c_double)
    si = real([40, 10, -20, 30], c_double)

    er = [6.655811471009945e+00_c_double, 5.159694646004461e+00_c_double, &
          1.738893605180155e+00_c_double, 1.228142680379087e+00_c_double]
    ei = [1.102940138744213e+00_c_double, 1.296375269817600e+00_c_double, &
          8.698934637859068e-01_c_double, -1.268954483973123e+00_c_double]

    ierr = ams_petsc_init('' // c_null_char)
    if (ierr /= 0) then
        print *, 'ams_petsc_init failed with code', ierr
        stop 1
    end if

    ierr = ams_solve_upper(n_edges, n_nodes, n_rp, n_val, edge_nodes, &
                           coords, rp, ci, sr, si, vre, vim, &
                           n_edges, xr, xi)
    if (ierr /= 0) then
        print *, 'ams_solve_upper failed with code', ierr
        stop 1
    end if

    diff2 = 0
    ref2 = 0
    do k = 1, n_edges
        diff2 = diff2 + (xr(k) - er(k))**2 + (xi(k) - ei(k))**2
        ref2 = ref2 + er(k)**2 + ei(k)**2
    end do
    rel = sqrt(diff2 / ref2)
    print '(a,es12.4)', 'fortran_check_relative_l2 ', rel

    if (rel > 1.0e-6_c_double) then
        print *, 'FORTRAN BINDING CHECK FAILED'
        stop 1
    end if
    print *, 'FORTRAN BINDING CHECK PASSED'

    ierr = ams_petsc_finalize()
    if (ierr /= 0) stop 1

end program fortran_example_check
