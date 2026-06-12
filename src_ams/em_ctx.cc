#include "em_ctx.h"


#include <petsc.h>

PetscErrorCode create_context(EMContext *ctx) {
  PetscInt i;
  std::vector<PetscInt> np_per_group, proc_range;

  PetscFunctionBegin;

  PetscCall(MPI_Comm_dup(MPI_COMM_WORLD, &ctx->world_comm));
  PetscCall(MPI_Comm_size(ctx->world_comm, &ctx->world_size));
  PetscCall(MPI_Comm_rank(ctx->world_comm, &ctx->world_rank));

  if (ctx->n_groups > ctx->world_size) {
    SETERRQ(ctx->world_comm, EM_ERR_USER, "Number of groups should be less than the total number of MPI processes.");
  }

  proc_range.resize(ctx->n_groups + 1);
  np_per_group.resize(ctx->n_groups);

  proc_range[0] = 0;
  for (i = 0; i < ctx->n_groups; ++i) {
    np_per_group[i] = ctx->world_size / ctx->n_groups;
    if (i < (ctx->world_size % ctx->n_groups)) {
      np_per_group[i] += 1;
    }
    proc_range[i + 1] = proc_range[i] + np_per_group[i];
  }

  ctx->group_id = 0;
  for (i = 0; i < ctx->n_groups; ++i) {
    if (ctx->world_rank >= proc_range[i] && ctx->world_rank < proc_range[i + 1]) {
      ctx->group_id = i;
      break;
    }
  }
  PetscCall(MPI_Comm_split(ctx->world_comm, ctx->group_id, ctx->world_rank, &ctx->group_comm));

  PetscCall(MPI_Comm_size(ctx->group_comm, &ctx->group_size));
  PetscCall(MPI_Comm_rank(ctx->group_comm, &ctx->group_rank));

  ctx->C = NULL;
  ctx->M = NULL;
  ctx->A = NULL;
  ctx->B = NULL;
  ctx->A_ksp = NULL;
  ctx->B_ksp = NULL;

  ctx->G = NULL;

  ctx->s.re = NULL;
  ctx->s.im = NULL;
  ctx->dual_e.re = NULL;
  ctx->dual_e.im = NULL;
  ctx->w = NULL;

  PetscCall(PetscClassIdRegister("EMCTX", &ctx->EMCTX_ID));
  PetscCall(PetscLogEventRegister("CreateLS", ctx->EMCTX_ID, &ctx->CreateLS));
  PetscCall(PetscLogEventRegister("AssembleMat", ctx->EMCTX_ID, &ctx->AssembleMat));
  PetscCall(PetscLogEventRegister("AssembleRHS", ctx->EMCTX_ID, &ctx->AssembleRHS));
  PetscCall(PetscLogEventRegister("SetupAMS", ctx->EMCTX_ID, &ctx->SetupAMS));
  PetscCall(PetscLogEventRegister("CreatePC", ctx->EMCTX_ID, &ctx->CreatePC));
  PetscCall(PetscLogEventRegister("SolveLS", ctx->EMCTX_ID, &ctx->SolveLS));
  PetscCall(PetscLogEventRegister("EstimateError", ctx->EMCTX_ID, &ctx->EstimateError));
  PetscCall(PetscLogEventRegister("CalculateRSP", ctx->EMCTX_ID, &ctx->CalculateRSP));
  PetscCall(PetscLogDefaultBegin());

  PetscCall(PetscViewerASCIIOpen(ctx->group_comm, string_format("%s-group-%03d.log", ctx->oprefix, ctx->group_id).c_str(), &ctx->LS_log));

  PetscFunctionReturn(0);
}

PetscErrorCode destroy_context(EMContext *ctx) {
  PetscFunctionBegin;

  PetscCall(PetscViewerDestroy(&ctx->LS_log));

  PetscCall(MPI_Comm_free(&ctx->group_comm));
  PetscCall(MPI_Comm_free(&ctx->world_comm));

  PetscFunctionReturn(0);
}

PetscErrorCode process_options(EMContext *ctx) {
  PetscBool flg;
  const char *InnerPCType[] = {"mixed", "ams", "direct"};
  const char *DirectSolverType[] = {"mumps", "superlu_dist"};

  PetscFunctionBegin;

  ctx->n_groups = 1;

  // ctx->inner_pc_type = Direct 、 AMS
  ctx->inner_pc_type = AMS;

  // ctx->direct_solver_type = MUMPS;
  ctx->direct_solver_type = NoneP;

  ctx->pc_threshold = 500000;

  ctx->K_max_it = 100;
  ctx->dual_rtol = 1.0E-9;
  PetscCall(PetscOptionsGetInt(NULL, NULL, "-em_outer_max_it",
                               &ctx->K_max_it, NULL));
  PetscCall(PetscOptionsGetReal(NULL, NULL, "-em_outer_rtol",
                                &ctx->dual_rtol, NULL));

  ctx->max_dofs = 1000000;

  ctx->n_tx_divisions = 20;

  PetscFunctionReturn(0);
}
