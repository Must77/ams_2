#ifndef PETSCAMSINTERFACE
#define PETSCAMSINTERFACE
#include "em_ctx.h"
#include <iostream>
#include <fstream>
using namespace std;



    PetscErrorCode test2Dim2(EMContext *ctx);

    void testeg1(EMContext *ctx);

    PetscErrorCode setup_ams(EMContext *ctx);

    PetscErrorCode create_linear_system(EMContext *ctx);
    
    PetscErrorCode assemble_matrix(EMContext *ctx);
    
    PetscErrorCode assemble_rhs_csem(EMContext *ctx);

    PetscErrorCode matshell_mult_a(Mat A, Vec x, Vec y);

    PetscErrorCode matshell_createvecs_a(Mat A, Vec *right, Vec *left);
    
    PetscErrorCode pc_apply_b(PC pc, Vec b, Vec x);

    PetscErrorCode create_pc(EMContext *ctx);

    PetscErrorCode solve_linear_system(EMContext *ctx, const PETScBlockVector &s, PETScBlockVector &e, PetscInt max_it, PetscReal rtol);

    PetscErrorCode destroy_pc(EMContext *ctx);

    PetscErrorCode destroy_linear_system(EMContext *ctx);

    //测试
    PetscErrorCode solve_eg1();

    // PetscErrorCode initAMS();

    // PetscErrorCode obainRCSMat();

    // PetscErrorCode obainRhs();

    // PetscErrorCode solve();

#endif