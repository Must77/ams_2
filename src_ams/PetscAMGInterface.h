#ifndef PETSCAMSINTERFACE
#define PETSCAMSINTERFACE
#include "em_ctx.h"
#include <iostream>
#include <fstream>
#include <vector>
using namespace std;



    PetscErrorCode test2Dim2(EMContext *ctx);

    PetscErrorCode testeg1(EMContext *ctx);

    PetscErrorCode expand_upper_triangle_1based_to_full_csr(
        PetscInt n_rows, PetscInt n_upper_row_ptr, PetscInt n_upper_values,
        const PetscInt *upper_row_ptr, const PetscInt *upper_col_idx,
        const PetscReal *upper_real, const PetscReal *upper_imag,
        std::vector<PetscInt> &row_ptr, std::vector<PetscInt> &col_idx,
        std::vector<PetscReal> &data_real, std::vector<PetscReal> &data_imag);

    PetscErrorCode prepare_fortran_upper_1based_inputs(
        PetscInt n_edges, PetscInt n_nodes, PetscInt n_upper_row_ptr,
        PetscInt n_upper_values, const PetscInt *edge_nodes_1based,
        const double *node_coords, const PetscInt *upper_row_ptr,
        const PetscInt *upper_col_idx, const PetscReal *rhs_real,
        const PetscReal *rhs_imag, const PetscReal *upper_real,
        const PetscReal *upper_imag, std::vector<double> &edgesN,
        std::vector<double> &nodes, std::vector<PetscInt> &row_ptr,
        std::vector<PetscInt> &col_idx, std::vector<PetscReal> &rhs_real_internal,
        std::vector<PetscReal> &rhs_imag_internal,
        std::vector<PetscReal> &data_real, std::vector<PetscReal> &data_imag);

    PetscErrorCode load_context_from_arrays(EMContext *ctx, PetscInt n_edges, PetscInt n_nodes,
                                            PetscInt n_row_ptr, PetscInt n_col_idx, PetscInt n_rhs,
                                            PetscInt n_values, const double *edgesN, const double *nodes,
                                            const PetscInt *row_ptr, const PetscInt *col_idx,
                                            const PetscReal *rhs_real, const PetscReal *rhs_imag,
                                            const PetscReal *data_real, const PetscReal *data_imag);

    PetscErrorCode setup_ams(EMContext *ctx);

    PetscErrorCode create_linear_system(EMContext *ctx);
    
    PetscErrorCode assemble_matrix(EMContext *ctx);
    
    PetscErrorCode assemble_rhs_csem(EMContext *ctx);

    PetscErrorCode matshell_mult_a(Mat A, Vec x, Vec y);

    PetscErrorCode matshell_createvecs_a(Mat A, Vec *right, Vec *left);
    
    PetscErrorCode pc_apply_b(PC pc, Vec b, Vec x);

    PetscErrorCode create_pc(EMContext *ctx);

    PetscErrorCode solve_linear_system(EMContext *ctx, const PETScBlockVector &s, PETScBlockVector &e, PetscInt max_it, PetscReal rtol);

    PetscErrorCode copy_result_to_arrays(EMContext *ctx, PetscInt n_result, PetscReal *out_real, PetscReal *out_imag);

    PetscErrorCode write_result_arrays(PetscInt n_result, const PetscReal *out_real, const PetscReal *out_imag);

    PetscErrorCode write_result(EMContext *ctx);

    PetscErrorCode destroy_pc(EMContext *ctx);

    PetscErrorCode destroy_linear_system(EMContext *ctx);

    //测试
    PetscErrorCode solve_eg1(PetscInt n_edges, PetscInt n_nodes,
                              PetscInt n_row_ptr, PetscInt n_col_idx,
                              PetscInt n_rhs, PetscInt n_values,
                              const double *edgesN, const double *nodes,
                              const PetscInt *row_ptr, const PetscInt *col_idx,
                              const PetscReal *rhs_real, const PetscReal *rhs_imag,
                              const PetscReal *data_real, const PetscReal *data_imag,
                              PetscInt n_result, PetscReal *out_real, PetscReal *out_imag);

    // rhs_imag follows the Fortran caller's physical imaginary part; this
    // adapter negates it for the internal real block system.
    PetscErrorCode solve_eg1_fortran_upper_1based(
        PetscInt n_edges, PetscInt n_nodes, PetscInt n_upper_row_ptr,
        PetscInt n_upper_values, const PetscInt *edge_nodes_1based,
        const double *node_coords, const PetscInt *upper_row_ptr,
        const PetscInt *upper_col_idx, const PetscReal *rhs_real,
        const PetscReal *rhs_imag, const PetscReal *upper_real,
        const PetscReal *upper_imag, PetscInt n_result, PetscReal *out_real,
        PetscReal *out_imag);

    // PetscErrorCode initAMS();

    // PetscErrorCode obainRCSMat();

    // PetscErrorCode obainRhs();

    // PetscErrorCode solve();

#endif
