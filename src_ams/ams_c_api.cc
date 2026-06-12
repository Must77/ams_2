#include "ams_c_api.h"

#include "PetscAMGInterface.h"

static_assert(sizeof(PetscInt) == sizeof(int),
              "ams_c_api 假定 32 位 PetscInt；64 位索引的 PETSc 需要改接口");
static_assert(sizeof(PetscReal) == sizeof(double),
              "ams_c_api 假定双精度 PetscReal");

namespace {
// 记录 PETSc 是否由本接口初始化，避免 finalize 调用方自己的 PETSc
bool petsc_owned_by_api = false;
}  // namespace

extern "C" {

int ams_petsc_init(const char* options_string) {
    PetscFunctionBegin;
    if (!PetscInitializeCalled) {
        PetscCall(PetscInitializeNoArguments());
        petsc_owned_by_api = true;
    }
    if (options_string && options_string[0]) {
        PetscCall(PetscOptionsInsertString(NULL, options_string));
    }
    PetscFunctionReturn(0);
}

int ams_solve_upper(int n_edges, int n_nodes, int n_upper_row_ptr,
                    int n_upper_values, const int* edge_nodes_1based,
                    const double* node_coords, const int* upper_row_ptr,
                    const int* upper_col_idx, const double* rhs_real,
                    const double* rhs_imag, const double* upper_real,
                    const double* upper_imag, int n_result, double* out_real,
                    double* out_imag) {
    PetscFunctionBegin;
    PetscCall(solve_eg1_fortran_upper_1based(
        n_edges, n_nodes, n_upper_row_ptr, n_upper_values, edge_nodes_1based,
        node_coords, upper_row_ptr, upper_col_idx, rhs_real, rhs_imag,
        upper_real, upper_imag, n_result, out_real, out_imag));
    PetscFunctionReturn(0);
}

int ams_petsc_finalize(void) {
    PetscFunctionBegin;
    if (petsc_owned_by_api && PetscInitializeCalled) {
        PetscCall(PetscFinalize());
        petsc_owned_by_api = false;
    }
    PetscFunctionReturn(0);
}

}  // extern "C"
