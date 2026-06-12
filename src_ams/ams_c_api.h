#ifndef AMS_C_API_H
#define AMS_C_API_H

/* C 接口：供 Fortran/C 调用 PETSc AMS 求解器，不暴露任何 PETSc 类型。
   约定：所有数组由调用方分配；矩阵为 Fortran 1-based 上三角 CSR（对称，
   只存上三角含对角）；本接口假定 PetscInt 为 32 位 int、PetscReal 为
   double（编译期断言保证）。 */

#ifdef __cplusplus
extern "C" {
#endif

/* 初始化 PETSc 并注入求解器选项。options_string 传 NULL 或空串表示
   不加选项；大算例的推荐选项（含边长反缩放与 beta Poisson 修复）：
   "-fortran_upper_unscale_edge_len -ams_beta_mass_poisson "
   "-ams_beta_mass_shift 1e-6 -B_ksp_type gmres -B_ksp_max_it 15 "
   "-em_outer_max_it 60 -em_outer_rtol 1e-7 -A_ksp_gmres_restart 15"
   若 MPI/PETSc 已被调用方初始化则跳过初始化只注入选项。 */
int ams_petsc_init(const char* options_string);

/* 求解 (K_r + i K_i) x = (S_r + i S_i)。K 以 1-based 上三角 CSR 给出，
   edge_nodes_1based 为每条边两端节点号（1-based，长 2*n_edges），
   node_coords 为节点坐标（长 3*n_nodes）。out_real/out_imag 长 n_result
   （= n_edges）。语义与 solve_eg1_fortran_upper_1based 完全一致。
   返回 0 表示成功，非 0 为 PETSc 错误码。 */
int ams_solve_upper(int n_edges, int n_nodes, int n_upper_row_ptr,
                    int n_upper_values, const int* edge_nodes_1based,
                    const double* node_coords, const int* upper_row_ptr,
                    const int* upper_col_idx, const double* rhs_real,
                    const double* rhs_imag, const double* upper_real,
                    const double* upper_imag, int n_result, double* out_real,
                    double* out_imag);

/* 释放 PETSc。只有当 PETSc 是由 ams_petsc_init 初始化时才调用。 */
int ams_petsc_finalize(void);

#ifdef __cplusplus
}
#endif

#endif
