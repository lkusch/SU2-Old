/*!
 * \file COneShotSolver.cpp
 * \brief Main subroutines for solving the discrete adjoint problem.
 * \author T. Albring
 * \version 7.1.1 "Blackbird"
 *
 * SU2 Project Website: https://su2code.github.io
 *
 * The SU2 Project is maintained by the SU2 Foundation
 * (http://su2foundation.org)
 *
 * Copyright 2012-2020, SU2 Contributors (cf. AUTHORS.md)
 *
 * SU2 is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * SU2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with SU2. If not, see <http://www.gnu.org/licenses/>.
 */

#include "../../include/solvers/COneShotSolver.hpp"
#include "../../../Common/include/toolboxes/geometry_toolbox.hpp"

COneShotSolver::COneShotSolver(void) : CDiscAdjSolver () {

}

COneShotSolver::COneShotSolver(CGeometry *geometry, CConfig *config)  : CDiscAdjSolver(geometry, config) {

}

COneShotSolver::COneShotSolver(CGeometry *geometry, CConfig *config, CSolver *direct_solver, unsigned short Kind_Solver, unsigned short iMesh)  : CDiscAdjSolver(geometry, config, direct_solver, Kind_Solver, iMesh) {
 theta = 0.0;
 rho = 0.0;
 nConstr = config->GetnConstr();

 DConsVec = new su2double** [nConstr];
 for (unsigned short iConstr=0; iConstr<nConstr;iConstr++){
   DConsVec[iConstr] = new su2double* [nPoint];
   for (unsigned long iPoint = 0; iPoint < nPoint; iPoint++){
     DConsVec[iConstr][iPoint] = new su2double [nVar];
     for (unsigned short iVar = 0; iVar < nVar; iVar++){
       DConsVec[iConstr][iPoint][iVar]=0.0;
     }
   }
 }
}

COneShotSolver::~COneShotSolver(void) {
  for (unsigned short iConstr=0; iConstr<nConstr;iConstr++){
    for (unsigned long iPoint = 0; iPoint < nPoint; iPoint++){
      delete [] DConsVec[iConstr][iPoint];
    }
    delete [] DConsVec[iConstr];
  }
  delete [] DConsVec;
}

void COneShotSolver::SetRecording(CGeometry* geometry, CConfig *config){

  bool time_n1_needed = config->GetTime_Marching() == DT_STEPPING_2ND;
  bool time_n_needed = (config->GetTime_Marching() == DT_STEPPING_1ST) || time_n1_needed;

  unsigned long iPoint;
  unsigned short iVar;

  /*--- For the one-shot solver the solution is not reset in each iteration step to the initial solution ---*/

  if (time_n_needed) {
    for (iPoint = 0; iPoint < nPoint; iPoint++) {
      for (iVar = 0; iVar < nVar; iVar++) {
        AD::ResetInput(direct_solver->GetNodes()->GetSolution_time_n(iPoint)[iVar]);
      }
    }
  }
  if (time_n1_needed) {
    for (iPoint = 0; iPoint < nPoint; iPoint++) {
      for (iVar = 0; iVar < nVar; iVar++) {
        AD::ResetInput(direct_solver->GetNodes()->GetSolution_time_n1(iPoint)[iVar]);
      }
    }
  }

  /*--- Set the Jacobian to zero since this is not done inside the fluid iteration
   * when running the discrete adjoint solver. ---*/

  direct_solver->Jacobian.SetValZero();

  /*--- Set indices to zero ---*/

  RegisterVariables(geometry, config, true);

}

void COneShotSolver::SetGeometrySensitivityLagrangian(CGeometry *geometry){

    unsigned short iDim, nDim, nPoint;
    unsigned long iPoint;

    nPoint  = geometry->GetnPoint();
    nDim    = geometry->GetnDim();

    geometry->InitializeSensitivity();
    for (iPoint = 0; iPoint < nPoint; iPoint++) {
      for (iDim = 0; iDim < nDim; iDim++) {
        geometry->SetSensitivity(iPoint, iDim, nodes->GetSensitivity_AugmentedLagrangian(iPoint,iDim));
      }
    }
}

void COneShotSolver::SetGeometrySensitivityGradient(CGeometry *geometry){

    unsigned short iDim, nDim, nPoint;
    unsigned long iPoint;

    nPoint  = geometry->GetnPoint();
    nDim    = geometry->GetnDim();

    geometry->InitializeSensitivity();

    for (iPoint = 0; iPoint < nPoint; iPoint++) {
      for (iDim = 0; iDim < nDim; iDim++) {
        geometry->SetSensitivity(iPoint, iDim, nodes->GetSensitivity_ShiftedLagrangian(iPoint,iDim));
      }
    }
}

void COneShotSolver::SaveSensitivity(CGeometry *geometry){
    unsigned short iDim, nDim, nPoint;
    unsigned long iPoint;

    nPoint  = geometry->GetnPoint();
    nDim    = geometry->GetnDim();

    for (iPoint = 0; iPoint < nPoint; iPoint++) {
      for (iDim = 0; iDim < nDim; iDim++) {
        nodes->SetSensitivity_ShiftedLagrangian(iPoint, iDim, nodes->GetSensitivity(iPoint, iDim));
      }
    }
}

void COneShotSolver::ResetSensitivityLagrangian(CGeometry *geometry){
    unsigned short iDim, nDim, nPoint;
    unsigned long iPoint;

    nPoint  = geometry->GetnPoint();
    nDim    = geometry->GetnDim();

    for (iPoint = 0; iPoint < nPoint; iPoint++) {
      for (iDim = 0; iDim < nDim; iDim++) {
        nodes->SetSensitivity_AugmentedLagrangian(iPoint, iDim, 0.0);
      }
    }
}

void COneShotSolver::UpdateSensitivityLagrangian(CGeometry *geometry, su2double factor){
    unsigned short iDim, nDim, nPoint;
    unsigned long iPoint;

    nPoint  = geometry->GetnPoint();
    nDim    = geometry->GetnDim();
    for (iPoint = 0; iPoint < nPoint; iPoint++) {
      for (iDim = 0; iDim < nDim; iDim++) {
        nodes->SetSensitivity_AugmentedLagrangian(iPoint, iDim, nodes->GetSensitivity_AugmentedLagrangian(iPoint, iDim)+factor*nodes->GetSensitivity(iPoint, iDim));
      }
    }
}

void COneShotSolver::StoreMeshPoints(CConfig *config, CGeometry *geometry){
    unsigned long iVertex, jPoint;
    unsigned short iMarker;
    for (jPoint=0; jPoint<geometry->GetnPoint();jPoint++){
        geometry->nodes->SetCoord_Store(jPoint, geometry->nodes->GetCoord(jPoint));
    }
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
        for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          jPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          geometry->vertex[iMarker][iVertex]->SetNormal_Old(geometry->vertex[iMarker][iVertex]->GetNormal());
        }
    }
}

void COneShotSolver::LoadMeshPoints(CConfig *config, CGeometry *geometry){
    unsigned long iVertex, jPoint;
    unsigned short iMarker;
    for (jPoint=0; jPoint<geometry->GetnPoint();jPoint++){
        geometry->nodes->SetCoord(jPoint, geometry->nodes->GetCoord_Store(jPoint));
    }
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
        for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          jPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          geometry->vertex[iMarker][iVertex]->SetNormal(geometry->vertex[iMarker][iVertex]->GetNormal_Old());
        }
    }
}

void COneShotSolver::StoreSolution(){
  direct_solver->GetNodes()->Set_StoreSolution();
  nodes->Set_StoreSolution();
}

void COneShotSolver::StoreFormerSolution(){
  direct_solver->GetNodes()->Set_FormerSolution();
  nodes->Set_FormerSolution();
}

void COneShotSolver::LoadSolution(){
  unsigned long iPoint;
  for (iPoint = 0; iPoint < nPoint; iPoint++){
    direct_solver->GetNodes()->SetSolution(iPoint, direct_solver->GetNodes()->GetSolution_Store(iPoint));
    nodes->SetSolution(iPoint, nodes->GetSolution_Store(iPoint));
  }
}

void COneShotSolver::LoadSolutionStep(su2double stepsize){
  unsigned long iPoint, iVar;
  for (iPoint = 0; iPoint < nPoint; iPoint++){
    for (iVar = 0; iVar < nVar; iVar++){
      direct_solver->GetNodes()->SetSolution(iPoint, iVar, direct_solver->GetNodes()->GetSolution_Former(iPoint, iVar)+stepsize*direct_solver->GetNodes()->GetSolution_Delta_Store(iPoint, iVar));
      nodes->SetSolution(iPoint, iVar, nodes->GetSolution_Former(iPoint, iVar)+stepsize*nodes->GetSolution_Delta_Store(iPoint, iVar));
    }
  }
}

void COneShotSolver::ShiftFormerSolution(){
  unsigned long iPoint;
  for (iPoint = 0; iPoint < nPoint; iPoint++){
    direct_solver->GetNodes()->SetSolution_Store(iPoint, direct_solver->GetNodes()->GetSolution_Former(iPoint));
    nodes->SetSolution_Store(iPoint, nodes->GetSolution_Former(iPoint));
  }
}

void COneShotSolver::ShiftStoreSolution(){
  unsigned long iPoint;
  for (iPoint = 0; iPoint < nPoint; iPoint++){
    direct_solver->GetNodes()->SetSolution_Former(iPoint, direct_solver->GetNodes()->GetSolution_Store(iPoint));
    nodes->SetSolution_Former(iPoint, nodes->GetSolution_Store(iPoint));
  }
}

void COneShotSolver::StoreSaveSolution(){
  direct_solver->GetNodes()->Set_SaveSolution();
  nodes->Set_SaveSolution();
}

void COneShotSolver::LoadSaveSolution(){
  unsigned long iPoint;
  for (iPoint = 0; iPoint < nPoint; iPoint++){
    direct_solver->GetNodes()->SetSolution(iPoint, direct_solver->GetNodes()->GetSolution_Save(iPoint));
    nodes->SetSolution(iPoint, nodes->GetSolution_Save(iPoint));
  }
}

void COneShotSolver::CalculateAlphaBeta(){
  unsigned short iVar;
  unsigned long iPoint;
  su2double helper=0.0;
  su2double normDelta = 0.0;
  su2double normDeltaNew = 0.0;
  /* --- Estimate rho and theta values --- */
  for (iPoint = 0; iPoint < nPoint; iPoint++){
    for (iVar = 0; iVar < nVar; iVar++){
      normDelta += direct_solver->GetNodes()->GetSolution_Delta(iPoint,iVar)*direct_solver->GetNodes()->GetSolution_Delta(iPoint, iVar);
      normDeltaNew += (direct_solver->GetNodes()->GetSolution(iPoint, iVar)-direct_solver->GetNodes()->GetSolution_Store(iPoint, iVar))*(direct_solver->GetNodes()->GetSolution(iPoint, iVar)-direct_solver->GetNodes()->GetSolution_Store(iPoint, iVar));
      helper += direct_solver->GetNodes()->GetSolution_Delta(iPoint, iVar)*(nodes->GetSolution(iPoint, iVar)-nodes->GetSolution_Store(iPoint, iVar))-nodes->GetSolution_Delta(iPoint, iVar)*(direct_solver->GetNodes()->GetSolution(iPoint, iVar)-direct_solver->GetNodes()->GetSolution_Store(iPoint, iVar));
    }
  }
  if(sqrt(normDeltaNew)/sqrt(normDelta)>0.9*rho) rho = sqrt(normDeltaNew)/sqrt(normDelta);
  if(fabs(helper)/normDelta>0.9*theta) theta = fabs(helper)/normDelta;
  std::cout<<"Estimate: "<<sqrt(normDeltaNew)/sqrt(normDelta)<<" "<<fabs(helper)/normDelta<<std::endl;
  std::cout<<"Rho, Theta: "<<rho<<" "<<theta<<std::endl;
  std::cout<<"Alpha, Beta: "<<2*theta/((1-rho)*(1-rho))<<" "<<2./theta<<std::endl;
}

void COneShotSolver::SetAlphaBeta(CConfig *config){
  config->SetOneShotAlpha(2*theta/((1-rho)*(1-rho)));
  config->SetOneShotBeta(2./theta);
}

su2double COneShotSolver::CalculateLagrangianPart(CConfig *config, bool augmented){
  unsigned short iVar;
  unsigned long iPoint;
  su2double Lagrangian=0.0;
  su2double helper=0.0;

  for (iPoint = 0; iPoint < nPoint; iPoint++){
    for (iVar = 0; iVar < nVar; iVar++){
      direct_solver->GetNodes()->SetSolution_Delta(iPoint, iVar, direct_solver->GetNodes()->GetSolution(iPoint, iVar)-direct_solver->GetNodes()->GetSolution_Store(iPoint, iVar));
    }
    for (iVar = 0; iVar < nVar; iVar++){
      nodes->SetSolution_Delta(iPoint, iVar, nodes->GetSolution(iPoint, iVar)-nodes->GetSolution_Store(iPoint, iVar));
    }
  }

  /* --- Calculate augmented Lagrangian terms (alpha and beta) --- */
  if(augmented){
    for (iPoint = 0; iPoint < nPoint; iPoint++){
      for (iVar = 0; iVar < nVar; iVar++){
        helper+=direct_solver->GetNodes()->GetSolution_Delta(iPoint, iVar)*direct_solver->GetNodes()->GetSolution_Delta(iPoint, iVar);
      }
    }
    Lagrangian+=helper*(config->GetOneShotAlpha()/2);
    helper=0.0;
    for (iPoint = 0; iPoint < nPoint; iPoint++){
      for (iVar = 0; iVar < nVar; iVar++){
        helper+=nodes->GetSolution_Delta(iPoint, iVar)*nodes->GetSolution_Delta(iPoint, iVar);
      }
    }
    Lagrangian+=helper*(config->GetOneShotBeta()/2);
  }

  helper=0.0;
  for (iPoint = 0; iPoint < nPoint; iPoint++){
    for (iVar = 0; iVar < nVar; iVar++){
      helper+=direct_solver->GetNodes()->GetSolution_Delta(iPoint, iVar)*nodes->GetSolution_Store(iPoint, iVar);
    }
  }
  Lagrangian+=helper;
  return Lagrangian;
}

void COneShotSolver::SetAdjoint_OutputUpdate(CGeometry *geometry, CConfig *config) {

  unsigned long iPoint;

  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    direct_solver->GetNodes()->SetAdjointSolution(iPoint, direct_solver->GetNodes()->GetSolution_Delta(iPoint));
  }
}

void COneShotSolver::SetAdjoint_OutputZero(CGeometry *geometry, CConfig *config) {

  unsigned long iPoint;
  unsigned short iVar;
  su2double * ZeroSolution = new su2double[nVar];
  for (iVar = 0; iVar < nVar; iVar++){
      ZeroSolution[iVar] = 0.0;
  }

  for (iPoint = 0; iPoint < nPoint; iPoint++) {
    direct_solver->GetNodes()->SetAdjointSolution(iPoint, ZeroSolution);
  }

  delete [] ZeroSolution;
}

void COneShotSolver::ExtractAdjoint_Solution_Clean(CGeometry *geometry, CConfig *config){

  unsigned long iPoint;

  for (iPoint = 0; iPoint < nPoint; iPoint++) {

    /*--- Extract the adjoint solution ---*/

    direct_solver->GetNodes()->GetAdjointSolution(iPoint, Solution);

    /*--- Store the adjoint solution ---*/

    nodes->SetSolution(iPoint, Solution);
  }

}

void COneShotSolver::UpdateStateVariable(CConfig *config){
    unsigned long iPoint;
    unsigned short iVar;
    su2double fd_step=config->GetFDStep();
    for (iPoint = 0; iPoint < nPoint; iPoint++){
      for (iVar = 0; iVar < nVar; iVar++){
        Solution[iVar] = direct_solver->GetNodes()->GetSolution_Store(iPoint, iVar)+fd_step*nodes->GetSolution_Delta(iPoint, iVar);
      }
      direct_solver->GetNodes()->SetSolution(iPoint, Solution);
    }
}

void COneShotSolver::SetFiniteDifferenceSens(CGeometry *geometry, CConfig* config){
    unsigned short iDim, nDim, nPoint;
    unsigned long iPoint;

    nPoint  = geometry->GetnPoint();
    nDim    = geometry->GetnDim();

    for (iPoint = 0; iPoint < nPoint; iPoint++) {
      for (iDim = 0; iDim < nDim; iDim++) {
        nodes->SetSensitivity(iPoint, iDim, (nodes->GetSensitivity(iPoint, iDim)-nodes->GetSensitivity_ShiftedLagrangian(iPoint, iDim))*(1./config->GetFDStep()));
      }
    }
}

void COneShotSolver::StoreSolutionDelta(){
  unsigned short iVar;
  unsigned long iPoint;

  for (iPoint = 0; iPoint < nPoint; iPoint++){
    for (iVar = 0; iVar < nVar; iVar++){
      direct_solver->GetNodes()->SetSolution_Delta_Store(iPoint, iVar, direct_solver->GetNodes()->GetSolution_Delta(iPoint, iVar));
    }
    for (iVar = 0; iVar < nVar; iVar++){
      nodes->SetSolution_Delta_Store(iPoint, iVar,nodes->GetSolution_Delta(iPoint, iVar));
    }
  }
}

void COneShotSolver::SetConstrDerivative(unsigned short iConstr){
  unsigned short iVar;
  unsigned long iPoint;

  for (iPoint = 0; iPoint < nPoint; iPoint++){
    for (iVar = 0; iVar < nVar; iVar++){
      DConsVec[iConstr][iPoint][iVar]=nodes->GetSolution(iPoint, iVar);;
    }
  }

}

su2double COneShotSolver::MultiplyConstrDerivative(unsigned short iConstr, unsigned short jConstr){
  unsigned short iVar;
  unsigned long iPoint;
  su2double product = 0.0;
  for (iPoint = 0; iPoint < nPoint; iPoint++){
    for (iVar = 0; iVar < nVar; iVar++){
      product+= DConsVec[iConstr][iPoint][iVar]*DConsVec[jConstr][iPoint][iVar];
    }
  }
  return product;
}



