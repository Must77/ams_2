// #include "em_ctx.h"
#include <iostream>
#include <fstream>
#include <vector>
#include "PetscAMGInterface.h"

using namespace std;

static PetscErrorCode read_eg1_input_files(
    std::size_t& nE, std::size_t& nN, std::size_t& nrow,
    std::size_t& ncol, std::size_t& nrhs, std::size_t& nval,
    std::vector<double>& edgesN, std::vector<double>& nodes,
    std::vector<PetscInt>& row_ptr, std::vector<PetscInt>& col_idx,
    std::vector<PetscReal>& rhs_real, std::vector<PetscReal>& rhs_imag,
    std::vector<PetscReal>& data_real, std::vector<PetscReal>& data_imag) {
    /* ---------- 读 edgesN ---------- */
    std::ifstream fe("edgesN.txt");
    fe >> nE;  // 第一行：行数

    edgesN.reserve(nE * 2);
    double x{};
    while (fe >> x)
        edgesN.push_back(x);  // 逐个压扁
    fe.close();
    cout << "edgesN: " << edgesN.size() / 2 << endl;

    /* ---------- 读 nodes ---------- */
    std::ifstream fn("nodes.txt");
    fn >> nN;  // 第一行：行数

    nodes.reserve(nN * 3);
    while (fn >> x)
        nodes.push_back(x);  // 逐个压扁
    fn.close();
    cout << "nodes: " << nodes.size() / 3 << endl;

    /* ---------- 读 rowptr---------- */
    std::ifstream frow("rowPtr.txt");
    frow >> nrow;  // 第一行：行数
    row_ptr.reserve(nrow);
    PetscInt ix{};
    while (frow >> ix)
        row_ptr.push_back(ix);  // 逐个压扁
    frow.close();

    /* ---------- 读 colIdx---------- */
    std::ifstream fcol("colIdx.txt");
    fcol >> ncol;  // 第一行：行数
    col_idx.reserve(ncol);
    while (fcol >> ix)
        col_idx.push_back(ix);  // 逐个压扁
    fcol.close();

    /* ---------- 读 rhs---------- */
    std::ifstream frhs("rhs.txt");
    frhs >> nrhs;  // 第一行：行数
    rhs_real.resize(nrhs);
    rhs_imag.resize(nrhs);
    for (int i = 0; i < nrhs; ++i) {
        frhs >> rhs_real[i] >> rhs_imag[i];
    }
    frhs.close();

    /* ---------- 读 val_real_imag---------- */
    std::ifstream fval("val_real_imag.txt");
    fval >> nval;  // 第一行：行数
    data_real.resize(nval);
    data_imag.resize(nval);
    for (int i = 0; i < nval; ++i) {
        fval >> data_real[i] >> data_imag[i];
    }
    fval.close();

    return 0;
}

int main(int argc, char **argv){
    

    PetscCall(PetscInitialize(NULL, NULL, NULL, NULL));   //初始化

    std::size_t nE, nN, nrow, ncol, nrhs, nval;
    std::vector<double> edgesN;
    std::vector<double> nodes;
    std::vector<PetscInt> row_ptr;
    std::vector<PetscInt> col_idx;
    std::vector<PetscReal> rhs_real;
    std::vector<PetscReal> rhs_imag;
    std::vector<PetscReal> data_real;
    std::vector<PetscReal> data_imag;

    PetscCall(read_eg1_input_files(nE, nN, nrow, ncol, nrhs, nval, edgesN,
                                   nodes, row_ptr, col_idx, rhs_real, rhs_imag,
                                   data_real, data_imag));

    std::vector<PetscReal> out_real(nrhs);
    std::vector<PetscReal> out_imag(nrhs);

    PetscCall(solve_eg1(static_cast<PetscInt>(nE), static_cast<PetscInt>(nN),
                        static_cast<PetscInt>(nrow),
                        static_cast<PetscInt>(ncol),
                        static_cast<PetscInt>(nrhs),
                        static_cast<PetscInt>(nval), edgesN.data(),
                        nodes.data(), row_ptr.data(), col_idx.data(),
                        rhs_real.data(), rhs_imag.data(), data_real.data(),
                        data_imag.data(), static_cast<PetscInt>(nrhs),
                        out_real.data(), out_imag.data()));

    PetscCall(write_result_arrays(static_cast<PetscInt>(nrhs), out_real.data(),
                                  out_imag.data()));

    PetscCall(PetscFinalize());

    return 0;
}
