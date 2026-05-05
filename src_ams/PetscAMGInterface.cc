#include "PetscAMGInterface.h"

PetscErrorCode test2Dim2(EMContext* ctx) {
    PetscFunctionBegin;

    int nrow = 2;
    int nnz = 4;

    ctx->rptr.resize(nrow + 1);
    ctx->cidx.resize(nnz);
    ctx->data_real.resize(nnz);
    ctx->data_imag.resize(nnz);
    ctx->b_real.resize(nrow);
    ctx->b_imag.resize(nrow);

    ctx->rptr[0] = 0;
    ctx->rptr[1] = 2;
    ctx->rptr[2] = 4;

    ctx->cidx[0] = 0;
    ctx->cidx[1] = 1;
    ctx->cidx[2] = 0;
    ctx->cidx[3] = 1;

    ctx->data_real[0] = 1;
    ctx->data_real[1] = 2;
    ctx->data_real[2] = 2;
    ctx->data_real[3] = 2.5;

    ctx->data_imag[0] = 3;
    ctx->data_imag[1] = 1;
    ctx->data_imag[2] = 1;
    ctx->data_imag[3] = 2;

    ctx->b_real[0] = 1;
    ctx->b_real[1] = 3;

    ctx->b_imag[0] = 0;
    ctx->b_imag[1] = 0;

    // cout << "rptr: ";
    // for(auto&ix : ctx->rptr)
    //     cout << ix << "\t";
    // cout << endl;
    // cout << "cidx: ";
    // for(auto&ix : ctx->cidx)
    //     cout << ix << "\t";
    // cout << endl;
    // cout << "data_real: ";
    // for(auto&ix : ctx->data_real)
    //     cout << ix << "\t";
    // cout << endl;

    PetscFunctionReturn(0);
}

PetscErrorCode load_context_from_arrays(EMContext* ctx, PetscInt n_edges,
                                        PetscInt n_nodes, PetscInt n_row_ptr,
                                        PetscInt n_col_idx, PetscInt n_rhs,
                                        PetscInt n_values, const double* edgesN,
                                        const double* nodes,
                                        const PetscInt* row_ptr,
                                        const PetscInt* col_idx,
                                        const PetscReal* rhs_real,
                                        const PetscReal* rhs_imag,
                                        const PetscReal* data_real,
                                        const PetscReal* data_imag) {
    PetscFunctionBegin;

    ctx->edgesN.assign(edgesN, edgesN + n_edges * 2);
    ctx->nodes.assign(nodes, nodes + n_nodes * 3);

    ctx->rptr.assign(row_ptr, row_ptr + n_row_ptr);
    ctx->cidx.assign(col_idx, col_idx + n_col_idx);

    ctx->b_real.assign(rhs_real, rhs_real + n_rhs);
    ctx->b_imag.assign(rhs_imag, rhs_imag + n_rhs);

    ctx->data_real.assign(data_real, data_real + n_values);
    ctx->data_imag.assign(data_imag, data_imag + n_values);

    PetscFunctionReturn(0);
}

PetscErrorCode testeg1(EMContext* ctx) {
    PetscFunctionBegin;
    //
    /* ---------- 读 edgesN ---------- */
    std::ifstream fe("edgesN.txt");
    std::size_t nE;
    fe >> nE;  // 第一行：行数

    std::vector<double> edgesN;
    edgesN.reserve(nE * 2);
    double x{};
    while (fe >> x)
        edgesN.push_back(x);  // 逐个压扁
    fe.close();
    cout << "edgesN: " << edgesN.size() / 2 << endl;
    /* ---------- 读 nodes ---------- */
    std::ifstream fn("nodes.txt");
    std::size_t nN;
    fn >> nN;  // 第一行：行数

    std::vector<double> nodes;
    nodes.reserve(nN * 3);
    while (fn >> x)
        nodes.push_back(x);  // 逐个压扁
    fn.close();
    cout << "nodes: " << nodes.size() / 3 << endl;
    /* ---------- 读 rowptr---------- */
    std::ifstream frow("rowPtr.txt");
    std::size_t nrow;
    frow >> nrow;  // 第一行：行数
    std::vector<PetscInt> row_ptr;
    row_ptr.reserve(nrow);
    PetscInt ix{};
    while (frow >> ix)
        row_ptr.push_back(ix);  // 逐个压扁
    frow.close();

    /* ---------- 读 colIdx---------- */
    std::ifstream fcol("colIdx.txt");
    std::size_t ncol;
    fcol >> ncol;  // 第一行：行数
    std::vector<PetscInt> col_idx;
    col_idx.reserve(ncol);
    while (fcol >> ix)
        col_idx.push_back(ix);  // 逐个压扁
    fcol.close();

    /* ---------- 读 rhs---------- */
    std::ifstream frhs("rhs.txt");
    std::size_t nrhs;
    frhs >> nrhs;  // 第一行：行数
    std::vector<PetscReal> rhs_real(nrhs);
    std::vector<PetscReal> rhs_imag(nrhs);
    for (int i = 0; i < nrhs; ++i) {
        frhs >> rhs_real[i] >> rhs_imag[i];
    }
    frhs.close();
    /* ---------- 读 val_real_imag---------- */
    std::ifstream fval("val_real_imag.txt");
    std::size_t nval;
    fval >> nval;  // 第一行：行数
    std::vector<PetscReal> data_real(nval);
    std::vector<PetscReal> data_imag(nval);
    for (int i = 0; i < nval; ++i) {
        fval >> data_real[i] >> data_imag[i];
    }
    fval.close();

    PetscCall(load_context_from_arrays(ctx, static_cast<PetscInt>(nE),
                                       static_cast<PetscInt>(nN),
                                       static_cast<PetscInt>(nrow),
                                       static_cast<PetscInt>(ncol),
                                       static_cast<PetscInt>(nrhs),
                                       static_cast<PetscInt>(nval),
                                       edgesN.data(), nodes.data(),
                                       row_ptr.data(), col_idx.data(),
                                       rhs_real.data(), rhs_imag.data(),
                                       data_real.data(), data_imag.data()));

    PetscFunctionReturn(0);
}

PetscErrorCode setup_ams(EMContext* ctx) {
    // TetAccessor tet;
    std::vector<PetscReal> vals;
    std::vector<PetscInt> row_ptr, col_idx;
    PetscInt t, e, eidx, v, vidx, n_local_edges, n_edges, n_local_vertices,
        n_vertices, begin, end;

    PetscFunctionBegin;

    LogEventHelper leh(ctx->SetupAMS);

    n_edges = ctx->edgesN.size() / 2;

    n_vertices = ctx->nodes.size() / 3;

    row_ptr.resize(n_edges + 1);
    col_idx.resize(n_edges * 2);
    vals.resize(n_edges * 2);

    std::fill(row_ptr.begin(), row_ptr.end(), -1);

    for (t = 0; t < n_edges; ++t) {
        row_ptr[t] = t * 2;
        col_idx[t * 2 + 0] = ctx->edgesN[t * 2 + 0];
        col_idx[t * 2 + 1] = ctx->edgesN[t * 2 + 1];
        vals[t * 2 + 0] = -1;
        vals[t * 2 + 1] = 1;
    }
    row_ptr[n_edges] = n_edges * 2;

    PetscCall(MatCreateMPIAIJWithArrays(ctx->group_comm, n_edges, n_vertices,
                                        n_edges, n_vertices, &row_ptr[0],
                                        &col_idx[0], &vals[0], &ctx->G));
    PetscCall(PetscObjectSetName((PetscObject)(ctx->G), "G_"));

    ctx->v_coords.resize(n_vertices * 3);

    for (int i = 0; i < n_vertices; ++i) {
        ctx->v_coords[i * 3 + 0] = ctx->nodes[i * 3 + 0];
        ctx->v_coords[i * 3 + 1] = ctx->nodes[i * 3 + 1];
        ctx->v_coords[i * 3 + 2] = ctx->nodes[i * 3 + 2];
    }

    PetscFunctionReturn(0);
}

PetscErrorCode create_linear_system(EMContext* ctx) {
    PetscFunctionBegin;

    LogEventHelper leh(ctx->CreateLS);

    // rcs 格式 的基本信息
    int ndim = ctx->rptr.size() - 1;
    cout << "ndim: " << ndim << endl;
    cout << "nzeros: " << ctx->cidx.size() << endl;

    PetscCall(MatCreateMPIAIJWithArrays(
        ctx->group_comm, ndim, ndim, PETSC_DECIDE, PETSC_DECIDE, &ctx->rptr[0],
        &ctx->cidx[0], NULL, &ctx->C));
    PetscCall(MatDuplicate(ctx->C, MAT_SHARE_NONZERO_PATTERN, &ctx->M));

    PetscCall(VecCreateMPI(ctx->group_comm, ndim, PETSC_DECIDE, &ctx->s.re));
    PetscCall(VecDuplicate(ctx->s.re, &ctx->s.im));
    PetscCall(VecDuplicate(ctx->s.re, &ctx->dual_e.re));
    PetscCall(VecDuplicate(ctx->dual_e.re, &ctx->dual_e.im));

    PetscCall(VecDuplicate(ctx->s.re, &ctx->w));

    if ((InnerPCType)ctx->inner_pc_type == AMS) {
        ctx->use_ams = PETSC_TRUE;
    } else {
        ctx->use_ams = PETSC_FALSE;
    }

    PetscFunctionReturn(0);
}

PetscErrorCode assemble_matrix(EMContext* ctx) {
    PetscFunctionBegin;

    int ndim = ctx->rptr.size() - 1;
    int nzeros = ctx->cidx.size();
    PetscCall(MatCreateMPIAIJWithArrays(
        ctx->group_comm, ndim, ndim, PETSC_DECIDE, PETSC_DECIDE, &ctx->rptr[0],
        &ctx->cidx[0], &ctx->data_real[0], &ctx->C));
    PetscCall(MatCreateMPIAIJWithArrays(
        ctx->group_comm, ndim, ndim, PETSC_DECIDE, PETSC_DECIDE, &ctx->rptr[0],
        &ctx->cidx[0], &ctx->data_imag[0], &ctx->M));

    PetscCall(MatAssemblyBegin(ctx->C, MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(ctx->C, MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyBegin(ctx->M, MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(ctx->M, MAT_FINAL_ASSEMBLY));

    PetscFunctionReturn(0);
}
PetscErrorCode assemble_rhs_csem(EMContext* ctx) {
    int ndim = ctx->rptr.size() - 1;
    PetscCall(VecCreateMPIWithArray(ctx->group_comm, 1, ndim, ndim,
                                    &ctx->b_real[0], &ctx->s.re));
    PetscCall(VecCreateMPIWithArray(ctx->group_comm, 1, ndim, ndim,
                                    &ctx->b_imag[0], &ctx->s.im));

    PetscCall(VecAssemblyBegin(ctx->s.re));
    PetscCall(VecAssemblyEnd(ctx->s.re));
    PetscCall(VecAssemblyBegin(ctx->s.im));
    PetscCall(VecAssemblyEnd(ctx->s.im));

    PetscFunctionReturn(0);
}

PetscErrorCode matshell_mult_a(Mat A, Vec x, Vec y) {
    EMContext* ctx;
    Vec xr, xi, yr, yi;

    PetscFunctionBegin;

    PetscCall(VecNestGetSubVec(x, 0, &xr));
    PetscCall(VecNestGetSubVec(x, 1, &xi));
    PetscCall(VecNestGetSubVec(y, 0, &yr));
    PetscCall(VecNestGetSubVec(y, 1, &yi));

    PetscCall(MatShellGetContext(A, &ctx));

    PetscCall(MatMult(ctx->C, xr, yr));
    PetscCall(MatMult(ctx->M, xi, ctx->w));
    PetscCall(VecAXPY(yr, -1.0, ctx->w));
    PetscCall(MatMult(ctx->M, xr, yi));
    PetscCall(VecScale(yi, -1.0));
    PetscCall(MatMult(ctx->C, xi, ctx->w));
    PetscCall(VecAXPY(yi, -1.0, ctx->w));

    PetscFunctionReturn(0);
}

PetscErrorCode matshell_createvecs_a(Mat A, Vec* right, Vec* left) {
    Vec xx[2], x;
    EMContext* ctx;

    PetscFunctionBegin;

    PetscCall(MatShellGetContext(A, &ctx));

    xx[0] = ctx->w;
    xx[1] = ctx->w;
    PetscCall(VecCreateNest(ctx->group_comm, 2, NULL, xx, &x));

    if (right) {
        PetscCall(VecDuplicate(x, right));
    }
    if (left) {
        PetscCall(VecDuplicate(x, left));
    }

    PetscCall(VecDestroy(&x));

    PetscFunctionReturn(0);
}

PetscErrorCode pc_apply_b(PC pc, Vec b, Vec x) {
    EMContext* ctx;
    Vec br, bi, xr, xi;

    PetscFunctionBegin;

    PetscCall(PCShellGetContext(pc, (void**)&ctx));
    PetscCall(PetscViewerASCIIPushTab(ctx->LS_log));

    PetscCall(VecNestGetSubVec(b, 0, &br));
    PetscCall(VecNestGetSubVec(b, 1, &bi));
    PetscCall(VecNestGetSubVec(x, 0, &xr));
    PetscCall(VecNestGetSubVec(x, 1, &xi));

    PetscCall(KSPSolve(ctx->B_ksp, br, xr));   
    PetscCall(KSPSolve(ctx->B_ksp, bi, xi));

    PetscCall(PetscViewerASCIIPopTab(ctx->LS_log));

    PetscFunctionReturn(0);
}

PetscErrorCode create_pc(EMContext* ctx) {
    PetscInt nrow;
    PC A_pc, B_pc;
    PetscViewerAndFormat* vf;

    PetscFunctionBegin;

    LogEventHelper leh(ctx->CreatePC);

    PetscCall(MatDuplicate(ctx->C, MAT_SHARE_NONZERO_PATTERN, &ctx->B));
    PetscCall(MatCopy(ctx->C, ctx->B, SAME_NONZERO_PATTERN));
    PetscCall(MatAXPY(ctx->B, 1.0, ctx->M, SAME_NONZERO_PATTERN));  // B = C+M
    PetscCall(PetscObjectSetName((PetscObject)(ctx->B), "B"));
    PetscCall(MatSetOptionsPrefix(ctx->B, "B_"));
    PetscCall(MatSetOption(ctx->B, MAT_SPD, PETSC_TRUE));
    PetscCall(MatSetFromOptions(ctx->B));

    PetscCall(KSPCreate(ctx->group_comm, &ctx->B_ksp));
    PetscCall(KSPSetOptionsPrefix(ctx->B_ksp, "B_"));
    PetscCall(KSPSetOperators(ctx->B_ksp, ctx->B, ctx->B));

    if (ctx->use_ams) {
        PetscCall(KSPSetNormType(ctx->B_ksp, KSP_NORM_UNPRECONDITIONED));
    }

    PetscCall(KSPGetPC(ctx->B_ksp, &B_pc));
    if (ctx->use_ams) {
        PetscCall(setup_ams(ctx));
        PetscCall(KSPSetType(ctx->B_ksp, KSPFCG));
        PetscCall(KSPSetTolerances(ctx->B_ksp, 1E-2, PETSC_DEFAULT,
                                   PETSC_DEFAULT, 100));
        PetscCall(PCSetType(B_pc, PCHYPRE));
        PetscCall(PCHYPRESetType(B_pc, "ams"));
        PetscCall(PCHYPRESetDiscreteGradient(B_pc, ctx->G));
        PetscCall(PCSetCoordinates(B_pc, 3, (PetscInt)ctx->v_coords.size() / 3,
                                   &ctx->v_coords[0]));
    } else {
        PetscCall(KSPSetType(ctx->B_ksp, KSPPREONLY));
        // PetscCall(PCSetType(B_pc, PCCHOLESKY));
        PetscCall(PCSetType(B_pc, PCNONE));
        if (ctx->direct_solver_type == MUMPS) {
            PetscCall(PCFactorSetMatSolverType(B_pc, MATSOLVERMUMPS));
        } else if (ctx->direct_solver_type == SUPERLUDIST) {
            PetscCall(PCFactorSetMatSolverType(B_pc, MATSOLVERSUPERLU_DIST));
        }
    }
    PetscCall(PCSetFromOptions(B_pc));
    PetscCall(PCSetUp(B_pc));

    if (ctx->use_ams) {
        PetscCall(
            PetscViewerAndFormatCreate(ctx->LS_log, PETSC_VIEWER_DEFAULT, &vf));
        PetscCall(KSPMonitorSet(
            ctx->B_ksp,
            (PetscErrorCode (*)(KSP, PetscInt, PetscReal,
                                void*))KSPMonitorTrueResidual,
            vf, (PetscErrorCode (*)(void**))PetscViewerAndFormatDestroy));
    }
    PetscCall(KSPSetFromOptions(ctx->B_ksp));
    PetscCall(KSPSetUp(ctx->B_ksp));

    nrow = ctx->rptr.size() - 1;
    PetscCall(MatCreateShell(ctx->group_comm, nrow * 2, nrow * 2, PETSC_DECIDE,
                             PETSC_DECIDE, ctx, &ctx->A));
    PetscCall(MatShellSetOperation(ctx->A, MATOP_MULT,
                                   (void (*)(void))matshell_mult_a));
    PetscCall(MatShellSetOperation(ctx->A, MATOP_MULT_TRANSPOSE,
                                   (void (*)(void))matshell_mult_a));
    PetscCall(MatShellSetOperation(ctx->A, MATOP_CREATE_VECS,
                                   (void (*)(void))matshell_createvecs_a));
    PetscCall(PetscObjectSetName((PetscObject)(ctx->A), "A"));
    PetscCall(MatSetOptionsPrefix(ctx->A, "A_"));
    PetscCall(MatSetFromOptions(ctx->A));

    PetscCall(KSPCreate(ctx->group_comm, &ctx->A_ksp));
    PetscCall(KSPSetOptionsPrefix(ctx->A_ksp, "A_"));
    PetscCall(KSPSetOperators(ctx->A_ksp, ctx->A, ctx->A));
    PetscCall(KSPSetType(ctx->A_ksp, KSPFGMRES));
    PetscCall(KSPGetPC(ctx->A_ksp, &A_pc));
    PetscCall(PCSetType(A_pc, PCSHELL));
    PetscCall(PCShellSetContext(A_pc, ctx));
    PetscCall(PCShellSetApply(A_pc, pc_apply_b));

    PetscCall(
        PetscViewerAndFormatCreate(ctx->LS_log, PETSC_VIEWER_DEFAULT, &vf));
    PetscCall(KSPMonitorSet(
        ctx->A_ksp,
        (PetscErrorCode (*)(KSP, PetscInt, PetscReal,
                            void*))KSPMonitorTrueResidual,
        vf, (PetscErrorCode (*)(void**))PetscViewerAndFormatDestroy));

    PetscCall(KSPSetFromOptions(ctx->A_ksp));
    PetscCall(KSPSetUp(ctx->A_ksp));

    PetscFunctionReturn(0);
}

PetscErrorCode solve_linear_system(EMContext* ctx,
                                   const PETScBlockVector& s,
                                   PETScBlockVector& e,
                                   PetscInt max_it,
                                   PetscReal rtol) {
    cout << "maxiter: " << max_it << endl;
    cout << "rtol:" << rtol << endl;

    Vec xx[2], bb[2], x, b;

    PetscFunctionBegin;

    LogEventHelper leh(ctx->SolveLS);

    xx[0] = e.re;
    xx[1] = e.im;
    PetscCall(VecCreateNest(ctx->group_comm, 2, NULL, xx, &x));
    bb[0] = s.re;
    bb[1] = s.im;
    PetscCall(VecCreateNest(ctx->group_comm, 2, NULL, bb, &b));

    PetscCall(VecZeroEntries(x));

    PetscCall(KSPSetTolerances(ctx->A_ksp, rtol, PETSC_DEFAULT, PETSC_DEFAULT,
                               max_it));
    PetscCall(KSPSolve(ctx->A_ksp, b, x));

    PetscCall(VecDestroy(&x));
    PetscCall(VecDestroy(&b));

    PetscFunctionReturn(0);
}

PetscErrorCode destroy_pc(EMContext* ctx) {
    PetscFunctionBegin;

    if (ctx->use_ams) {
        PetscCall(MatDestroy(&ctx->G));
    }

    PetscCall(MatDestroy(&ctx->A));
    PetscCall(MatDestroy(&ctx->B));
    PetscCall(KSPDestroy(&ctx->A_ksp));
    PetscCall(KSPDestroy(&ctx->B_ksp));

    PetscFunctionReturn(0);
}

PetscErrorCode destroy_linear_system(EMContext* ctx) {
    PetscFunctionBegin;

    PetscCall(MatDestroy(&ctx->C));
    PetscCall(MatDestroy(&ctx->M));

    PetscCall(VecDestroy(&ctx->s.re));
    PetscCall(VecDestroy(&ctx->s.im));
    PetscCall(VecDestroy(&ctx->dual_e.re));
    PetscCall(VecDestroy(&ctx->dual_e.im));
    PetscCall(VecDestroy(&ctx->w));

    PetscFunctionReturn(0);
}

PetscErrorCode copy_result_to_arrays(EMContext* ctx, PetscInt n_result,
                                     PetscReal* out_real,
                                     PetscReal* out_imag) {
    Vec xr = ctx->dual_e.re;
    Vec xi = ctx->dual_e.im;

    PetscInt n;
    PetscFunctionBegin;
    PetscCall(VecGetLocalSize(xr, &n));  // 本地长度
    PetscCheck(n_result == n, PETSC_COMM_SELF, EM_ERR_USER,
               "Result output length does not match local solution size.");

    const PetscScalar *arr_r, *arr_i;
    PetscCall(VecGetArrayRead(xr, &arr_r));
    PetscCall(VecGetArrayRead(xi, &arr_i));

    for (PetscInt i = 0; i < n; ++i) {
        out_real[i] = arr_r[i];
        out_imag[i] = arr_i[i];
    }

    PetscCall(VecRestoreArrayRead(xr, &arr_r));
    PetscCall(VecRestoreArrayRead(xi, &arr_i));

    PetscFunctionReturn(0);
}

PetscErrorCode write_result_arrays(PetscInt n_result, const PetscReal* out_real,
                                   const PetscReal* out_imag) {
    PetscFunctionBegin;

    ofstream ofs("result.txt");
    for (PetscInt i = 0; i < n_result; ++i) {
        ofs << i << "\t" << out_real[i] << " " << out_imag[i] << endl;
    }
    ofs.close();

    PetscFunctionReturn(0);
}

PetscErrorCode write_result(EMContext* ctx) {
    PetscInt n;

    PetscFunctionBegin;
    PetscCall(VecGetLocalSize(ctx->dual_e.re, &n));  // 本地长度

    std::vector<PetscReal> s_real(n);
    std::vector<PetscReal> s_imag(n);
    PetscCall(copy_result_to_arrays(ctx, n, s_real.data(), s_imag.data()));

    PetscCall(write_result_arrays(n, s_real.data(), s_imag.data()));

    PetscFunctionReturn(0);
}

PetscErrorCode solve_eg1(PetscInt n_edges, PetscInt n_nodes,
                         PetscInt n_row_ptr, PetscInt n_col_idx,
                         PetscInt n_rhs, PetscInt n_values,
                         const double* edgesN, const double* nodes,
                         const PetscInt* row_ptr, const PetscInt* col_idx,
                         const PetscReal* rhs_real, const PetscReal* rhs_imag,
                         const PetscReal* data_real,
                         const PetscReal* data_imag, PetscInt n_result,
                         PetscReal* out_real, PetscReal* out_imag) {
    PetscFunctionBegin;  // petsc 开始

    EMContext ctx;

    PetscCall(process_options(&ctx));  // 一些 ams 的选项设置

    PetscCall(create_context(&ctx));  // 主要作用分配内存，空地址，并行没测试过

    // 需要数据
    // 标准 CSR 格式 的 jcol row val_real val_imag b_real b_imag
    // ams 所需要的 edgetonode[numedge.2] coods[numnodes,3]
    PetscCall(load_context_from_arrays(
        &ctx, n_edges, n_nodes, n_row_ptr, n_col_idx, n_rhs, n_values, edgesN,
        nodes, row_ptr, col_idx, rhs_real, rhs_imag, data_real, data_imag));

    PetscCall(PetscViewerASCIIPushTab(ctx.LS_log));  // petsc 检测信息
    //
    create_linear_system(&ctx);  // 创建线性系统，分配内存，创建系数csr 骨架

    PetscCall(assemble_matrix(&ctx));  // 把testeg1 厘米的数据写入pestc 专用 矩阵

    PetscCall(create_pc(&ctx));  // 创建预处理器，包括 创建 ams 的 G 矩阵

    PetscCall(assemble_rhs_csem(&ctx));  // 有端项目写入 pectc ，每次一个

    PetscCall(PetscViewerASCIIPrintf(ctx.LS_log, "Solving for dual mode:\n"));
    PetscCall(solve_linear_system(&ctx, ctx.s, ctx.dual_e, ctx.K_max_it,
                                  ctx.dual_rtol));  // 求解

    PetscCall(copy_result_to_arrays(&ctx, n_result, out_real, out_imag));

    // 释放内存
    PetscCall(destroy_pc(&ctx));

    PetscCall(destroy_linear_system(&ctx));

    PetscCall(PetscViewerASCIIPopTab(ctx.LS_log));

    PetscCall(destroy_context(&ctx));

    PetscFunctionReturn(0);
}
