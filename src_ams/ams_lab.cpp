// Fast AMS experiment driver: loads B/G/coords dumped by upper_dump_system
// and runs one inner-style solve with a fully option-controllable KSP/PC.
// Usage: AmsLab <dump_dir> [-lab_ksp_type fcg] [-lab_pc ams|none|jacobi]
//        plus any -lab_ksp_* / -lab_pc_* PETSc options.
#include <petscksp.h>

#include <string>

int main(int argc, char** argv) {
    const std::string dir = argc > 1 ? argv[1] : ".";

    PetscCall(PetscInitialize(&argc, &argv, NULL, NULL));

    Mat B, G;
    Vec coords, br, bi;
    PetscViewer vw;
    PetscCall(PetscViewerBinaryOpen(PETSC_COMM_SELF, (dir + "/B.dat").c_str(),
                                    FILE_MODE_READ, &vw));
    PetscCall(MatCreate(PETSC_COMM_SELF, &B));
    PetscCall(MatLoad(B, vw));
    PetscCall(PetscViewerDestroy(&vw));
    PetscCall(PetscViewerBinaryOpen(PETSC_COMM_SELF, (dir + "/G.dat").c_str(),
                                    FILE_MODE_READ, &vw));
    PetscCall(MatCreate(PETSC_COMM_SELF, &G));
    PetscCall(MatLoad(G, vw));
    PetscCall(PetscViewerDestroy(&vw));
    PetscCall(PetscViewerBinaryOpen(PETSC_COMM_SELF, (dir + "/aux.dat").c_str(),
                                    FILE_MODE_READ, &vw));
    PetscCall(VecCreate(PETSC_COMM_SELF, &coords));
    PetscCall(VecLoad(coords, vw));
    PetscCall(VecCreate(PETSC_COMM_SELF, &br));
    PetscCall(VecLoad(br, vw));
    PetscCall(VecCreate(PETSC_COMM_SELF, &bi));
    PetscCall(VecLoad(bi, vw));
    PetscCall(PetscViewerDestroy(&vw));

    PetscInt n, nc;
    PetscCall(MatGetSize(B, &n, NULL));
    PetscCall(VecGetSize(coords, &nc));
    PetscCall(PetscPrintf(PETSC_COMM_SELF, "loaded B %d x %d, coords %d\n", n,
                          n, nc / 3));

    char pcname[64] = "ams";
    PetscCall(PetscOptionsGetString(NULL, NULL, "-lab_pc", pcname,
                                    sizeof(pcname), NULL));

    KSP ksp;
    PC pc;
    PetscCall(KSPCreate(PETSC_COMM_SELF, &ksp));
    PetscCall(KSPSetOptionsPrefix(ksp, "lab_"));
    PetscCall(KSPSetOperators(ksp, B, B));
    PetscCall(KSPSetType(ksp, KSPFCG));
    PetscCall(KSPSetNormType(ksp, KSP_NORM_UNPRECONDITIONED));
    PetscCall(KSPSetTolerances(ksp, 1e-2, PETSC_DEFAULT, PETSC_DEFAULT, 100));
    PetscCall(KSPGetPC(ksp, &pc));
    if (std::string(pcname) == "ams") {
        PetscCall(PCSetType(pc, PCHYPRE));
        PetscCall(PCHYPRESetType(pc, "ams"));
        PetscCall(PCHYPRESetDiscreteGradient(pc, G));
        PetscBool beta_mass = PETSC_FALSE;
        PetscCall(PetscOptionsGetBool(NULL, NULL, "-lab_beta_mass",
                                      &beta_mass, NULL));
        if (beta_mass) {
            Mat Mmat, Abeta;
            PetscViewer mv;
            PetscCall(PetscViewerBinaryOpen(PETSC_COMM_SELF,
                                            (dir + "/M.dat").c_str(),
                                            FILE_MODE_READ, &mv));
            PetscCall(MatCreate(PETSC_COMM_SELF, &Mmat));
            PetscCall(MatLoad(Mmat, mv));
            PetscCall(PetscViewerDestroy(&mv));
            PetscReal mass_shift = 0.0;
            PetscCall(PetscOptionsGetReal(NULL, NULL, "-lab_mass_shift",
                                          &mass_shift, NULL));
            if (mass_shift > 0) {
                PetscReal mn;
                PetscCall(MatNorm(Mmat, NORM_INFINITY, &mn));
                PetscCall(MatShift(Mmat, mass_shift * mn));
                PetscCall(PetscPrintf(PETSC_COMM_SELF,
                                      "mass shifted by %g\n",
                                      (double)(mass_shift * mn)));
            }
            PetscCall(MatPtAP(Mmat, G, MAT_INITIAL_MATRIX, 2.0, &Abeta));
            PetscReal nrm;
            PetscCall(MatNorm(Abeta, NORM_INFINITY, &nrm));
            PetscCall(PetscPrintf(PETSC_COMM_SELF,
                                  "beta mass matrix norm %g\n", (double)nrm));
            PetscCall(PCHYPRESetBetaPoissonMatrix(pc, Abeta));
        }
        const PetscScalar* ca;
        PetscCall(VecGetArrayRead(coords, &ca));
        PetscCall(PCSetCoordinates(pc, 3, nc / 3, (PetscReal*)ca));
        PetscCall(VecRestoreArrayRead(coords, &ca));
    } else {
        PetscCall(PCSetType(pc, pcname));
    }
    PetscCall(KSPSetFromOptions(ksp));

    Vec x, b;
    PetscCall(MatCreateVecs(B, &x, &b));
    PetscBool use_random = PETSC_FALSE;
    PetscCall(PetscOptionsGetBool(NULL, NULL, "-lab_random_rhs", &use_random,
                                  NULL));
    if (use_random) {
        PetscRandom rnd;
        PetscCall(PetscRandomCreate(PETSC_COMM_SELF, &rnd));
        PetscCall(VecSetRandom(b, rnd));
        PetscCall(PetscRandomDestroy(&rnd));
    } else {
        PetscCall(VecCopy(br, b));
    }
    PetscReal bn;
    PetscCall(VecNorm(b, NORM_2, &bn));
    if (bn == 0) {
        PetscCall(PetscPrintf(PETSC_COMM_SELF, "rhs is zero, using random\n"));
        PetscRandom rnd;
        PetscCall(PetscRandomCreate(PETSC_COMM_SELF, &rnd));
        PetscCall(VecSetRandom(b, rnd));
        PetscCall(PetscRandomDestroy(&rnd));
    }

    PetscCall(KSPSolve(ksp, b, x));

    KSPConvergedReason reason;
    PetscInt its;
    PetscReal rnorm, xn;
    PetscCall(KSPGetConvergedReason(ksp, &reason));
    PetscCall(KSPGetIterationNumber(ksp, &its));
    PetscCall(KSPGetResidualNorm(ksp, &rnorm));
    PetscCall(VecNorm(x, NORM_2, &xn));
    PetscCall(PetscPrintf(
        PETSC_COMM_SELF, "reason %s its %d rnorm %g xnorm %g\n",
        KSPConvergedReasons[reason], its, (double)rnorm, (double)xn));

    PetscCall(PetscFinalize());
    return 0;
}
