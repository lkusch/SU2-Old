#include "../include/linear_solvers_structure_b.hpp"
#include "../include/linear_solvers_structure.hpp"
#include "../include/vector_structure.hpp"
#include "../include/matrix_structure.hpp"
#ifdef CODI_REVERSE_TYPE
void CSysSolve_b::Solve_b(AD::CheckpointHandler* data){
  /* --- Extract data from the checkpoint handler --- */

  int *LinSysRes_Indices;
  int *LinSysSol_Indices;

  data->getData(LinSysRes_Indices);
  data->getData(LinSysSol_Indices);

  unsigned long nBlk, nVar, nBlkDomain, size, i;

  data->getData(size);
  data->getData(nBlk);
  data->getData(nVar);
  data->getData(nBlkDomain);

  CSysMatrix* Jacobian;
  data->getData(Jacobian);

  CGeometry* geometry;
  data->getData(geometry);

  CConfig* config;
  data->getData(config);

  CSysVector LinSysRes_b(nBlk, nBlkDomain, nVar, 0.0);
  CSysVector LinSysSol_b(nBlk, nBlkDomain, nVar, 0.0);
  su2double Residual;

  unsigned long MaxIter = config->GetLinear_Solver_Iter();
  su2double SolverTol = config->GetLinear_Solver_Error();

  /* --- Initialize the right-hand side with the gradient of the solution of the primal linear system --- */

  for (i = 0; i < size; i ++){
    int& index = LinSysSol_Indices[i];
    LinSysRes_b[i] = AD::globalTape.getGradient(index);
    LinSysSol_b[i] = 0.0;
    AD::globalTape.gradient(index) = 0.0;
  }
  /* --- Set up preconditioner and matrix-vector product --- */

  CPreconditioner* precond  = NULL;

  switch(config->GetKind_DiscAdj_Linear_Prec()){
    case ILU:
      precond = new CILUPreconditioner(*Jacobian, geometry, config);
      break;
    case JACOBI:
      precond = new CJacobiPreconditioner(*Jacobian, geometry, config);
      break;
  }

  CMatrixVectorProduct* mat_vec = new CSysMatrixVectorProductTransposed(*Jacobian, geometry, config);

  CSysSolve *solver = new CSysSolve;

  /* --- Solve the system --- */

  switch(config->GetKind_DiscAdj_Linear_Solver()){
    case FGMRES:
      solver->FGMRES_LinSolver(LinSysRes_b, LinSysSol_b, *mat_vec, *precond, SolverTol , MaxIter, &Residual, false);
      break;
    case BCGSTAB:
      solver->BCGSTAB_LinSolver(LinSysRes_b, LinSysSol_b, *mat_vec, *precond, SolverTol , MaxIter, &Residual, false);
      break;
  }


  /* --- Update the gradients of the right-hand side of the primal linear system --- */

  for (i = 0; i < size; i ++){
    int& index = LinSysRes_Indices[i];
    AD::globalTape.gradient(index) += SU2_TYPE::GetPrimary(LinSysSol_b[i]);
  }

  delete mat_vec;
  delete precond;
  delete solver;
}


void CSysSolve_b::Delete_b(AD::CheckpointHandler* data){

  int *LinSysRes_Indices;
  int *LinSysSol_Indices;

  data->getData(LinSysRes_Indices);
  data->getData(LinSysSol_Indices);

  delete [] LinSysRes_Indices;
  delete [] LinSysSol_Indices;

  unsigned long nBlk, nVar, nBlkDomain, size;

  data->getData(size);
  data->getData(nBlk);
  data->getData(nVar);
  data->getData(nBlkDomain);

  CSysMatrix* Jacobian;
  data->getData(Jacobian);

  CGeometry* geometry;
  data->getData(geometry);

  CConfig* config;
  data->getData(config);
}
#endif