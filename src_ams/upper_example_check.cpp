#include "PetscAMGInterface.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace {

std::string join_path(const std::string& dir, const std::string& name) {
    if (!dir.empty() && dir[dir.size() - 1] == '/') {
        return dir + name;
    }
    return dir + "/" + name;
}

void require_open(const std::ifstream& stream, const std::string& path) {
    if (!stream) {
        throw std::runtime_error("failed to open " + path);
    }
}

PetscInt read_header_count(const std::string& line) {
    std::istringstream iss(line);
    std::string token;
    PetscInt value = 0;
    while (iss >> token) {
        std::istringstream maybe_int(token);
        if (maybe_int >> value) {
            return value;
        }
    }
    throw std::runtime_error("missing count in header: " + line);
}

std::vector<PetscInt> read_indexed_int_values(const std::string& path,
                                              PetscInt expected_count) {
    std::ifstream in(path.c_str());
    require_open(in, path);

    std::string header;
    std::getline(in, header);

    std::vector<PetscInt> values;
    values.reserve(expected_count);
    PetscInt id, value;
    while (in >> id >> value) {
        values.push_back(value);
    }

    if (static_cast<PetscInt>(values.size()) != expected_count) {
        throw std::runtime_error(path + ": unexpected integer value count");
    }
    return values;
}

std::vector<PetscReal> read_indexed_real_values(const std::string& path,
                                                PetscInt expected_count) {
    std::ifstream in(path.c_str());
    require_open(in, path);

    std::string header;
    std::getline(in, header);

    std::vector<PetscReal> values;
    values.reserve(expected_count);
    PetscInt id;
    PetscReal value;
    while (in >> id >> value) {
        values.push_back(value);
    }

    if (static_cast<PetscInt>(values.size()) != expected_count) {
        throw std::runtime_error(path + ": unexpected real value count");
    }
    return values;
}

void read_edges(const std::string& path, PetscInt expected_edges,
                std::vector<PetscInt>& edge_nodes_1based) {
    std::ifstream in(path.c_str());
    require_open(in, path);

    PetscInt n_edges;
    in >> n_edges;
    if (n_edges != expected_edges) {
        throw std::runtime_error(path + ": edge count does not match matrix");
    }

    edge_nodes_1based.reserve(expected_edges * 2);
    PetscInt id, n1, n2;
    for (PetscInt i = 0; i < expected_edges; ++i) {
        in >> id >> n1 >> n2;
        if (!in || id != i + 1) {
            throw std::runtime_error(path + ": invalid edge row");
        }
        edge_nodes_1based.push_back(n1);
        edge_nodes_1based.push_back(n2);
    }
}

PetscInt read_nodes(const std::string& path, std::vector<double>& nodes) {
    std::ifstream in(path.c_str());
    require_open(in, path);

    PetscInt n_nodes, dim, flag1, flag2;
    in >> n_nodes >> dim >> flag1 >> flag2;
    if (!in || dim != 3) {
        throw std::runtime_error(path + ": expected 3D node header");
    }

    nodes.reserve(n_nodes * 3);
    PetscInt id;
    double x, y, z;
    for (PetscInt i = 0; i < n_nodes; ++i) {
        in >> id >> x >> y >> z;
        if (!in || id != i + 1) {
            throw std::runtime_error(path + ": invalid node row");
        }
        nodes.push_back(x);
        nodes.push_back(y);
        nodes.push_back(z);
    }
    return n_nodes;
}

void read_rhs(const std::string& path, PetscInt expected_edges,
              std::vector<PetscReal>& rhs_real,
              std::vector<PetscReal>& rhs_imag) {
    std::ifstream in(path.c_str());
    require_open(in, path);

    std::string header;
    std::getline(in, header);

    rhs_real.reserve(expected_edges);
    rhs_imag.reserve(expected_edges);
    PetscInt id;
    PetscReal real, imag, unused1, unused2;
    for (PetscInt i = 0; i < expected_edges; ++i) {
        in >> id >> real >> imag >> unused1 >> unused2;
        if (!in || id != i + 1) {
            throw std::runtime_error(path + ": invalid RHS row");
        }
        rhs_real.push_back(real);
        rhs_imag.push_back(imag);
    }
}

void read_solution(const std::string& path, PetscInt expected_edges,
                   std::vector<PetscReal>& real,
                   std::vector<PetscReal>& imag) {
    std::ifstream in(path.c_str());
    require_open(in, path);

    real.reserve(expected_edges);
    imag.reserve(expected_edges);
    PetscReal r, i;
    while (in >> r >> i) {
        real.push_back(r);
        imag.push_back(i);
    }

    if (static_cast<PetscInt>(real.size()) != expected_edges) {
        throw std::runtime_error(path + ": unexpected solution length");
    }
}

void report_solution_error(const char* prefix, const std::vector<PetscReal>& real,
                           const std::vector<PetscReal>& imag,
                           const std::vector<PetscReal>& expected_real,
                           const std::vector<PetscReal>& expected_imag) {
    long double diff2 = 0;
    long double ref2 = 0;
    PetscReal max_abs = 0;
    PetscInt max_idx = 0;
    for (PetscInt i = 0; i < static_cast<PetscInt>(real.size()); ++i) {
        const long double dr = real[i] - expected_real[i];
        const long double di = imag[i] - expected_imag[i];
        const long double er = expected_real[i];
        const long double ei = expected_imag[i];
        const PetscReal abs_diff =
            static_cast<PetscReal>(std::sqrt(dr * dr + di * di));
        if (abs_diff > max_abs) {
            max_abs = abs_diff;
            max_idx = i;
        }
        diff2 += dr * dr + di * di;
        ref2 += er * er + ei * ei;
    }

    const PetscReal rel_l2 =
        ref2 > 0 ? static_cast<PetscReal>(std::sqrt(diff2 / ref2)) : 0;
    PetscCallAbort(PETSC_COMM_WORLD,
                   PetscPrintf(PETSC_COMM_WORLD, "%s_relative_l2 %.16e\n",
                               prefix, static_cast<double>(rel_l2)));
    PetscCallAbort(PETSC_COMM_WORLD,
                   PetscPrintf(PETSC_COMM_WORLD, "%s_max_abs %.16e at %d\n",
                               prefix, static_cast<double>(max_abs), max_idx));
}

}  // namespace

int main(int argc, char** argv) {
    const std::string data_dir = argc > 1 ? argv[1] : "../ams_example";

    PetscCall(PetscInitialize(&argc, &argv, NULL, NULL));

    try {
        const std::string jcol_path = join_path(data_dir, "jcol.txt");

        std::ifstream jcol_header_file(jcol_path.c_str());
        require_open(jcol_header_file, jcol_path);
        std::string jcol_header;
        std::getline(jcol_header_file, jcol_header);
        const PetscInt n_upper_row_ptr = read_header_count(jcol_header);
        const PetscInt n_edges = n_upper_row_ptr - 1;

        const std::string irw_path = join_path(data_dir, "irw.txt");
        std::ifstream irw_header_file(irw_path.c_str());
        require_open(irw_header_file, irw_path);
        std::string irw_header;
        std::getline(irw_header_file, irw_header);
        const PetscInt n_upper_values = read_header_count(irw_header);

        std::vector<PetscInt> upper_row_ptr =
            read_indexed_int_values(jcol_path, n_upper_row_ptr);
        std::vector<PetscInt> upper_col_idx =
            read_indexed_int_values(irw_path, n_upper_values);
        std::vector<PetscReal> upper_real = read_indexed_real_values(
            join_path(data_dir, "ssor_k1.txt"), n_upper_values);
        std::vector<PetscReal> upper_imag = read_indexed_real_values(
            join_path(data_dir, "ssor_k2.txt"), n_upper_values);

        std::vector<PetscInt> edge_nodes_1based;
        read_edges(join_path(data_dir, "edges.txt"), n_edges,
                   edge_nodes_1based);

        std::vector<double> node_coords;
        const PetscInt n_nodes = read_nodes(join_path(data_dir, "node.txt"),
                                            node_coords);

        std::vector<PetscReal> rhs_real;
        std::vector<PetscReal> rhs_imag;
        read_rhs(join_path(data_dir, "right-hand.txt"), n_edges, rhs_real,
                 rhs_imag);

        std::vector<PetscReal> expected_real;
        std::vector<PetscReal> expected_imag;
        read_solution(join_path(data_dir, "slove.txt"), n_edges, expected_real,
                      expected_imag);

        std::vector<PetscReal> ref_matvec_real(n_edges, 0);
        std::vector<PetscReal> ref_matvec_imag(n_edges, 0);
        for (PetscInt r = 0; r < n_edges; ++r) {
            const PetscInt begin = upper_row_ptr[r] - 1;
            const PetscInt end = upper_row_ptr[r + 1] - 1;
            for (PetscInt k = begin; k < end; ++k) {
                const PetscInt c = upper_col_idx[k] - 1;
                const PetscReal real = upper_real[k];
                const PetscReal imag = upper_imag[k];

                ref_matvec_real[r] += real * expected_real[c] -
                                      imag * expected_imag[c];
                ref_matvec_imag[r] += -imag * expected_real[c] -
                                      real * expected_imag[c];

                if (c != r) {
                    ref_matvec_real[c] += real * expected_real[r] -
                                          imag * expected_imag[r];
                    ref_matvec_imag[c] += -imag * expected_real[r] -
                                          real * expected_imag[r];
                }
            }
        }

        long double ref_residual2 = 0;
        long double rhs_norm2 = 0;
        PetscReal ref_max_residual = 0;
        for (PetscInt i = 0; i < n_edges; ++i) {
            const long double dr = ref_matvec_real[i] - rhs_real[i];
            const long double di = ref_matvec_imag[i] + rhs_imag[i];
            const PetscReal abs_residual =
                static_cast<PetscReal>(std::sqrt(dr * dr + di * di));
            ref_max_residual = std::max(ref_max_residual, abs_residual);
            ref_residual2 += dr * dr + di * di;
            rhs_norm2 += static_cast<long double>(rhs_real[i]) * rhs_real[i] +
                         static_cast<long double>(rhs_imag[i]) * rhs_imag[i];
        }
        const PetscReal ref_relative_residual =
            rhs_norm2 > 0
                ? static_cast<PetscReal>(std::sqrt(ref_residual2 / rhs_norm2))
                : 0;

        PetscCall(
            PetscPrintf(PETSC_COMM_WORLD, "running_fortran_upper_solve\n"));
        std::vector<PetscReal> out_real(n_edges);
        std::vector<PetscReal> out_imag(n_edges);
        PetscCall(solve_eg1_fortran_upper_1based(
            n_edges, n_nodes, n_upper_row_ptr, n_upper_values,
            edge_nodes_1based.data(), node_coords.data(), upper_row_ptr.data(),
            upper_col_idx.data(), rhs_real.data(), rhs_imag.data(),
            upper_real.data(), upper_imag.data(), n_edges, out_real.data(),
            out_imag.data()));

        PetscCall(PetscPrintf(PETSC_COMM_WORLD, "n_edges %d\n", n_edges));
        PetscCall(PetscPrintf(PETSC_COMM_WORLD, "n_nodes %d\n", n_nodes));
        PetscCall(PetscPrintf(PETSC_COMM_WORLD,
                              "reference_relative_residual %.16e\n",
                              static_cast<double>(ref_relative_residual)));
        PetscCall(PetscPrintf(PETSC_COMM_WORLD,
                              "reference_max_abs_residual %.16e\n",
                              static_cast<double>(ref_max_residual)));
        report_solution_error("solve", out_real, out_imag, expected_real,
                              expected_imag);
    } catch (const std::exception& e) {
        SETERRQ(PETSC_COMM_SELF, EM_ERR_USER, "%s", e.what());
    }

    PetscCall(PetscFinalize());
    return 0;
}
