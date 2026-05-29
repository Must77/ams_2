#include "PetscAMGInterface.h"

#include <vector>

int main(int argc, char** argv) {
    PetscCall(PetscInitialize(&argc, &argv, NULL, NULL));

    const PetscInt n_rows = 4;
    const PetscInt n_upper_row_ptr = 5;
    const PetscInt n_upper_values = 7;

    const PetscInt upper_row_ptr[] = {1, 3, 5, 7, 8};
    const PetscInt upper_col_idx[] = {1, 3, 2, 4, 3, 4, 4};
    const PetscReal upper_real[] = {10, 20, 30, 40, 50, 60, 70};
    const PetscReal upper_imag[] = {1, 2, 3, 4, 5, 6, 7};

    std::vector<PetscInt> row_ptr;
    std::vector<PetscInt> col_idx;
    std::vector<PetscReal> data_real;
    std::vector<PetscReal> data_imag;

    PetscCall(expand_upper_triangle_1based_to_full_csr(
        n_rows, n_upper_row_ptr, n_upper_values, upper_row_ptr, upper_col_idx,
        upper_real, upper_imag, row_ptr, col_idx, data_real, data_imag));

    const std::vector<PetscInt> expected_row_ptr = {0, 2, 4, 7, 10};
    const std::vector<PetscInt> expected_col_idx = {0, 2, 1, 3, 0,
                                                    2, 3, 1, 2, 3};
    const std::vector<PetscReal> expected_real = {10, 20, 30, 40, 20,
                                                  50, 60, 40, 60, 70};
    const std::vector<PetscReal> expected_imag = {1, 2, 3, 4, 2,
                                                  5, 6, 4, 6, 7};

    PetscCheck(row_ptr == expected_row_ptr, PETSC_COMM_SELF, EM_ERR_USER,
               "Expanded row pointer does not match expected full CSR.");
    PetscCheck(col_idx == expected_col_idx, PETSC_COMM_SELF, EM_ERR_USER,
               "Expanded column index does not match expected full CSR.");
    PetscCheck(data_real == expected_real, PETSC_COMM_SELF, EM_ERR_USER,
               "Expanded real values do not match expected full CSR.");
    PetscCheck(data_imag == expected_imag, PETSC_COMM_SELF, EM_ERR_USER,
               "Expanded imaginary values do not match expected full CSR.");

    const PetscInt edge_nodes_1based[] = {1, 2, 2, 3, 3, 4, 4, 5};
    const double node_coords[] = {0, 0, 0, 1, 0, 0, 1, 1,
                                  0, 0, 1, 0, 0, 0, 1};
    const PetscReal rhs_real[] = {100, 200, 300, 400};
    const PetscReal rhs_imag[] = {10, 20, 30, 40};

    std::vector<double> edgesN;
    std::vector<double> nodes;
    std::vector<PetscReal> rhs_real_internal;
    std::vector<PetscReal> rhs_imag_internal;

    PetscCall(prepare_fortran_upper_1based_inputs(
        n_rows, 5, n_upper_row_ptr, n_upper_values, edge_nodes_1based,
        node_coords, upper_row_ptr, upper_col_idx, rhs_real, rhs_imag,
        upper_real, upper_imag, edgesN, nodes, row_ptr, col_idx,
        rhs_real_internal, rhs_imag_internal, data_real, data_imag));

    const std::vector<double> expected_edgesN = {0, 1, 1, 2, 2, 3, 3, 4};
    const std::vector<double> expected_nodes = {0, 0, 0, 1, 0, 0, 1, 1,
                                                0, 0, 1, 0, 0, 0, 1};
    const std::vector<PetscReal> expected_rhs_real = {100, 200, 300, 400};
    const std::vector<PetscReal> expected_rhs_imag = {-10, -20, -30, -40};

    PetscCheck(edgesN == expected_edgesN, PETSC_COMM_SELF, EM_ERR_USER,
               "Prepared edge-to-node array does not match zero-based input.");
    PetscCheck(nodes == expected_nodes, PETSC_COMM_SELF, EM_ERR_USER,
               "Prepared node coordinate array does not match expected input.");
    PetscCheck(row_ptr == expected_row_ptr, PETSC_COMM_SELF, EM_ERR_USER,
               "Prepared row pointer does not match expected full CSR.");
    PetscCheck(col_idx == expected_col_idx, PETSC_COMM_SELF, EM_ERR_USER,
               "Prepared column index does not match expected full CSR.");
    PetscCheck(rhs_real_internal == expected_rhs_real, PETSC_COMM_SELF,
               EM_ERR_USER, "Prepared RHS real array does not match input.");
    PetscCheck(rhs_imag_internal == expected_rhs_imag, PETSC_COMM_SELF,
               EM_ERR_USER, "Prepared RHS imaginary array was not negated.");
    PetscCheck(data_real == expected_real, PETSC_COMM_SELF, EM_ERR_USER,
               "Prepared real values do not match expected full CSR.");
    PetscCheck(data_imag == expected_imag, PETSC_COMM_SELF, EM_ERR_USER,
               "Prepared imaginary values do not match expected full CSR.");

    PetscCall(PetscFinalize());
    return 0;
}
