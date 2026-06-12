// Dump the inner AMS system (B = C + M, discrete gradient G, vertex
// coordinates, rhs) built from collaborator upper-triangle data to PETSc
// binary files, so solver configurations can be tested without re-parsing
// ~1.7 GB of text each run. Honors -fortran_upper_unscale_edge_len.
#include "PetscAMGInterface.h"

#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

std::string join_path(const std::string& dir, const std::string& name) {
    if (dir.empty() || dir[dir.size() - 1] == '/') return dir + name;
    return dir + "/" + name;
}

void require_open(std::ifstream& in, const std::string& path) {
    if (!in.is_open()) throw std::runtime_error("cannot open " + path);
}

PetscInt read_header_count(const std::string& header) {
    std::istringstream iss(header);
    std::string tok;
    PetscInt last = -1;
    while (iss >> tok) {
        try { last = static_cast<PetscInt>(std::stoll(tok)); }
        catch (...) {}
    }
    if (last < 0) throw std::runtime_error("bad header: " + header);
    return last;
}

std::vector<PetscInt> read_indexed_int(const std::string& path, PetscInt n) {
    std::ifstream in(path.c_str());
    require_open(in, path);
    std::string header;
    std::getline(in, header);
    std::vector<PetscInt> out(n);
    long long idx, val;
    for (PetscInt i = 0; i < n; ++i) {
        if (!(in >> idx >> val)) throw std::runtime_error("short " + path);
        out[i] = static_cast<PetscInt>(val);
    }
    return out;
}

std::vector<PetscReal> read_indexed_real(const std::string& path, PetscInt n) {
    std::ifstream in(path.c_str());
    require_open(in, path);
    std::string header;
    std::getline(in, header);
    std::vector<PetscReal> out(n);
    long long idx;
    double val;
    for (PetscInt i = 0; i < n; ++i) {
        if (!(in >> idx >> val)) throw std::runtime_error("short " + path);
        out[i] = val;
    }
    return out;
}

}  // namespace

int main(int argc, char** argv) {
    const std::string data_dir = argc > 1 ? argv[1] : "../ams_example";
    const std::string out_dir = argc > 2 ? argv[2] : ".";

    PetscCall(PetscInitialize(&argc, &argv, NULL, NULL));
    try {
        const std::string jcol_path = join_path(data_dir, "jcol.txt");
        std::ifstream jh(jcol_path.c_str());
        require_open(jh, jcol_path);
        std::string header;
        std::getline(jh, header);
        const PetscInt n_row_ptr = read_header_count(header);
        const PetscInt n_edges = n_row_ptr - 1;

        const std::string irw_path = join_path(data_dir, "irw.txt");
        std::ifstream ih(irw_path.c_str());
        require_open(ih, irw_path);
        std::getline(ih, header);
        const PetscInt n_vals = read_header_count(header);

        std::vector<PetscInt> up_rp = read_indexed_int(jcol_path, n_row_ptr);
        std::vector<PetscInt> up_ci = read_indexed_int(irw_path, n_vals);
        std::vector<PetscReal> up_re =
            read_indexed_real(join_path(data_dir, "ssor_k1.txt"), n_vals);
        std::vector<PetscReal> up_im =
            read_indexed_real(join_path(data_dir, "ssor_k2.txt"), n_vals);

        std::vector<PetscInt> edge_nodes(n_edges * 2);
        {
            const std::string p = join_path(data_dir, "edges.txt");
            std::ifstream in(p.c_str());
            require_open(in, p);
            std::getline(in, header);
            long long id, a, b;
            for (PetscInt e = 0; e < n_edges; ++e) {
                if (!(in >> id >> a >> b)) throw std::runtime_error("short " + p);
                edge_nodes[e * 2 + 0] = static_cast<PetscInt>(a);
                edge_nodes[e * 2 + 1] = static_cast<PetscInt>(b);
            }
        }

        PetscInt n_nodes = 0;
        std::vector<double> node_coords;
        {
            const std::string p = join_path(data_dir, "node.txt");
            std::ifstream in(p.c_str());
            require_open(in, p);
            std::getline(in, header);
            std::istringstream iss(header);
            long long nn;
            iss >> nn;
            n_nodes = static_cast<PetscInt>(nn);
            node_coords.resize(n_nodes * 3);
            long long id;
            for (PetscInt v = 0; v < n_nodes; ++v) {
                if (!(in >> id >> node_coords[v * 3] >> node_coords[v * 3 + 1] >>
                      node_coords[v * 3 + 2]))
                    throw std::runtime_error("short " + p);
            }
        }

        std::vector<PetscReal> rhs_re(n_edges), rhs_im(n_edges);
        {
            const std::string p = join_path(data_dir, "right-hand.txt");
            std::ifstream in(p.c_str());
            require_open(in, p);
            std::getline(in, header);
            long long id;
            double re, im, c4, c5;
            for (PetscInt e = 0; e < n_edges; ++e) {
                if (!(in >> id >> re >> im >> c4 >> c5))
                    throw std::runtime_error("short " + p);
                rhs_re[e] = re;
                rhs_im[e] = im;
            }
        }

        std::vector<double> edgesN, nodes;
        std::vector<PetscInt> rp, ci;
        std::vector<PetscReal> bre, bim, dre, dim;
        PetscCall(prepare_fortran_upper_1based_inputs(
            n_edges, n_nodes, n_row_ptr, n_vals, edge_nodes.data(),
            node_coords.data(), up_rp.data(), up_ci.data(), rhs_re.data(),
            rhs_im.data(), up_re.data(), up_im.data(), edgesN, nodes, rp, ci,
            bre, bim, dre, dim));

        // B = C + M = data_real + data_imag, same CSR pattern
        std::vector<PetscScalar> bvals(dre.size());
        for (size_t k = 0; k < dre.size(); ++k) bvals[k] = dre[k] + dim[k];

        Mat B;
        PetscCall(MatCreateSeqAIJWithArrays(PETSC_COMM_SELF, n_edges, n_edges,
                                            rp.data(), ci.data(), bvals.data(),
                                            &B));

        // M = -Ki alone (mass / conductivity part), cancellation-free
        std::vector<PetscScalar> mvals(dim.begin(), dim.end());
        Mat Mmat;
        PetscCall(MatCreateSeqAIJWithArrays(PETSC_COMM_SELF, n_edges, n_edges,
                                            rp.data(), ci.data(), mvals.data(),
                                            &Mmat));

        // G: row e = -1 at first node, +1 at second (0-based ids in edgesN)
        std::vector<PetscInt> grp(n_edges + 1), gci(n_edges * 2);
        std::vector<PetscScalar> gv(n_edges * 2);
        for (PetscInt e = 0; e < n_edges; ++e) {
            grp[e] = e * 2;
            gci[e * 2 + 0] = (PetscInt)edgesN[e * 2 + 0];
            gci[e * 2 + 1] = (PetscInt)edgesN[e * 2 + 1];
            gv[e * 2 + 0] = -1.0;
            gv[e * 2 + 1] = 1.0;
        }
        grp[n_edges] = n_edges * 2;
        Mat G;
        PetscCall(MatCreateSeqAIJWithArrays(PETSC_COMM_SELF, n_edges, n_nodes,
                                            grp.data(), gci.data(), gv.data(),
                                            &G));

        Vec coords, brv, biv;
        PetscCall(VecCreateSeq(PETSC_COMM_SELF, n_nodes * 3, &coords));
        for (PetscInt i = 0; i < n_nodes * 3; ++i)
            PetscCall(VecSetValue(coords, i, nodes[i], INSERT_VALUES));
        PetscCall(VecAssemblyBegin(coords));
        PetscCall(VecAssemblyEnd(coords));
        PetscCall(VecCreateSeq(PETSC_COMM_SELF, n_edges, &brv));
        PetscCall(VecCreateSeq(PETSC_COMM_SELF, n_edges, &biv));
        for (PetscInt i = 0; i < n_edges; ++i) {
            PetscCall(VecSetValue(brv, i, bre[i], INSERT_VALUES));
            PetscCall(VecSetValue(biv, i, bim[i], INSERT_VALUES));
        }
        PetscCall(VecAssemblyBegin(brv));
        PetscCall(VecAssemblyEnd(brv));
        PetscCall(VecAssemblyBegin(biv));
        PetscCall(VecAssemblyEnd(biv));

        PetscViewer vw;
        PetscCall(PetscViewerBinaryOpen(PETSC_COMM_SELF,
                                        join_path(out_dir, "B.dat").c_str(),
                                        FILE_MODE_WRITE, &vw));
        PetscCall(MatView(B, vw));
        PetscCall(PetscViewerDestroy(&vw));
        PetscCall(PetscViewerBinaryOpen(PETSC_COMM_SELF,
                                        join_path(out_dir, "M.dat").c_str(),
                                        FILE_MODE_WRITE, &vw));
        PetscCall(MatView(Mmat, vw));
        PetscCall(PetscViewerDestroy(&vw));
        PetscCall(PetscViewerBinaryOpen(PETSC_COMM_SELF,
                                        join_path(out_dir, "G.dat").c_str(),
                                        FILE_MODE_WRITE, &vw));
        PetscCall(MatView(G, vw));
        PetscCall(PetscViewerDestroy(&vw));
        PetscCall(PetscViewerBinaryOpen(PETSC_COMM_SELF,
                                        join_path(out_dir, "aux.dat").c_str(),
                                        FILE_MODE_WRITE, &vw));
        PetscCall(VecView(coords, vw));
        PetscCall(VecView(brv, vw));
        PetscCall(VecView(biv, vw));
        PetscCall(PetscViewerDestroy(&vw));

        PetscCall(PetscPrintf(PETSC_COMM_SELF,
                              "dumped n_edges %d n_nodes %d nnz %d to %s\n",
                              n_edges, n_nodes, (PetscInt)dre.size(),
                              out_dir.c_str()));
    } catch (const std::exception& e) {
        SETERRQ(PETSC_COMM_SELF, EM_ERR_USER, "%s", e.what());
    }
    PetscCall(PetscFinalize());
    return 0;
}
