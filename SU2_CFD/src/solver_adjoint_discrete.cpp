/*!
 * \file solver_adjoint_discrete.cpp
 * \brief Main subroutines for solving the discrete adjoint problem.
 * \author T. Albring
 * \version 4.2.0 "Cardinal"
 *
 * SU2 Lead Developers: Dr. Francisco Palacios (Francisco.D.Palacios@boeing.com).
 *                      Dr. Thomas D. Economon (economon@stanford.edu).
 *
 * SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
 *                 Prof. Piero Colonna's group at Delft University of Technology.
 *                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *                 Prof. Rafael Palacios' group at Imperial College London.
 *
 * Copyright (C) 2012-2016 SU2, the open-source CFD code.
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

#include "../include/solver_structure.hpp"
#include <iostream>
#include <string>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <Eigen/Dense>

CDiscAdjSolver::CDiscAdjSolver(void) : CSolver (){

}

CDiscAdjSolver::CDiscAdjSolver(CGeometry *geometry, CConfig *config)  : CSolver(){

}

CDiscAdjSolver::CDiscAdjSolver(CGeometry *geometry, CConfig *config, CSolver *direct_solver, unsigned short Kind_Solver, unsigned short iMesh)  : CSolver(){

  unsigned short iVar, iMarker, iDim;
  unsigned long jVertex, iElem;

  bool restart = config->GetRestart();

  unsigned long iVertex, iPoint, index;
  string text_line, mesh_filename;
  ifstream restart_file;
  string filename, AdjExt;
  su2double dull_val;
  bool compressible = (config->GetKind_Regime() == COMPRESSIBLE);
  bool incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);

  turbulent = false, flow = false;

  switch(config->GetKind_Solver()){
    case DISC_ADJ_EULER: case DISC_ADJ_NAVIER_STOKES: flow = true; break;
    case DISC_ADJ_RANS: flow = true; turbulent = true; break;
  }

  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  nVar = direct_solver->GetnVar();
  nDim = geometry->GetnDim();

  /*--- Initialize arrays to NULL ---*/

  CSensitivity = NULL;

  Sens_Geo   = NULL;
  Sens_Mach  = NULL;
  Sens_AoA   = NULL;
  Sens_Press = NULL;
  Sens_Temp  = NULL;

  /*-- Store some information about direct solver ---*/
  this->KindDirect_Solver = Kind_Solver;
  this->direct_solver = direct_solver;


  nMarker      = config->GetnMarker_All();
  nPoint       = geometry->GetnPoint();
  nPointDomain = geometry->GetnPointDomain();

  /*--- Allocate the node variables ---*/

  node = new CVariable*[nPoint];

  /*--- Define some auxiliary vectors related to the residual ---*/

  Residual      = new su2double[nVar];         for (iVar = 0; iVar < nVar; iVar++) Residual[iVar]      = 1.0;
  Residual_RMS  = new su2double[nVar];         for (iVar = 0; iVar < nVar; iVar++) Residual_RMS[iVar]  = 1.0;
  Residual_Max  = new su2double[nVar];         for (iVar = 0; iVar < nVar; iVar++) Residual_Max[iVar]  = 1.0;

  /*--- Define some structures for locating max residuals ---*/

  Point_Max     = new unsigned long[nVar];  for (iVar = 0; iVar < nVar; iVar++) Point_Max[iVar]     = 0;
  Point_Max_Coord = new su2double*[nVar];
  for (iVar = 0; iVar < nVar; iVar++) {
    Point_Max_Coord[iVar] = new su2double[nDim];
    for (iDim = 0; iDim < nDim; iDim++) Point_Max_Coord[iVar][iDim] = 0.0;
  }

  /*--- Define some auxiliary vectors related to the solution ---*/

  Solution   = new su2double[nVar];

  for (iVar = 0; iVar < nVar; iVar++) Solution[iVar]   = 1e-16;

  /*--- Sensitivity definition and coefficient in all the markers ---*/


  machp=new su2double[4];
  points=new su2double[4];
  weights=new su2double[4];
  sigma=sqrt(0.0001);
  mu=0.8;
  points[0]=-1.65068012;
  points[1]=-0.52464762;
  points[2]=0.52464762;
  points[3]=1.65068012;
  weights[0]=0.08131284;
  weights[1]=0.80491409;
  weights[2]=0.80491409;
  weights[3]=0.08131284;
  for (unsigned short numQuad=0; numQuad<4;numQuad++){
      machp[numQuad]=sqrt(2)*sigma*points[numQuad]+mu;
  }

  //for quadrature
  CSensitivityQuad=new su2double**[100];
  LagrangeSensQuad=new su2double**[100];
  Lagrangian_Value_Quad=new su2double[100];
  Constraint_Save_Quad=new su2double*[100];
  Obj_Save_Quad=new su2double[100];
  for (unsigned short iCounter=0;iCounter<100;iCounter++){
      CSensitivityQuad[iCounter]=new su2double*[nMarker];
      LagrangeSensQuad[iCounter]=new su2double*[nMarker];
      Constraint_Save_Quad[iCounter]=new su2double[config->GetConstraintNum()];
      for (iMarker = 0; iMarker < nMarker; iMarker++) {
          CSensitivityQuad[iCounter][iMarker]= new su2double[geometry->nVertex[iMarker]];
          LagrangeSensQuad[iCounter][iMarker]= new su2double[geometry->nVertex[iMarker]];
      }
  }
  //end quadrature

  rkStore=new su2double* [config->GetLBFGSNum()];
  dukStore=new su2double* [config->GetLBFGSNum()];
  CSensitivity = new su2double* [nMarker];
  CSensitivityOld = new su2double* [nMarker];
  LagrangeSens=new su2double* [nMarker];
  ExpCSensitivityOld = new su2double* [nMarker];
  ExpLagrangeSens=new su2double* [nMarker];
  ProjectedSens=new su2double[38];
  ProjectedGradient=new su2double[38];
  ProjectedSensOld=new su2double[38];
  DesignVar=new su2double[38];
  DesignVarUpdate=new su2double[38];
  DesignVarUpdateSave=new su2double[38];
  DesignVarUpdateReal=new su2double[38];
  TotalIterations=0;
  for (iMarker=0; iMarker<38; iMarker++){
      DesignVar[iMarker]=0.0;
      DesignVarUpdate[iMarker]=0.0;
      DesignVarUpdateSave[iMarker]=0.0;
      DesignVarUpdateReal[iMarker]=0.0;
  }
  rho=0.01;
  normyold=1;
  BFGSCount=0;
  ewbasis=new su2double [3];
  evbasis= new su2double * [3];
  for (iMarker=0; iMarker<3; iMarker++){
      evbasis[iMarker]=new su2double[142];
  }
  lenNormal = new su2double[200];
  ConstraintFunc_Value=new su2double [config->GetConstraintNum()];
  Constraint_Save=new su2double [config->GetConstraintNum()];
  ExpConstraint_Save=new su2double [config->GetConstraintNum()];
  Constraint_Old=new su2double [config->GetConstraintNum()];
  cons_factor=new double [config->GetConstraintNum()];
  multiplier=new double [config->GetConstraintNum()];
  multiplierhelp=new double [config->GetConstraintNum()];
  multiplieroriginal=new double [config->GetConstraintNum()];
  for (iMarker=0; iMarker < config->GetConstraintNum();iMarker++){
      ConstraintFunc_Value[iMarker]=0.0;
      Constraint_Save[iMarker]=0.0;
      ExpConstraint_Save[iMarker]=0.0;
      Constraint_Old[iMarker]=0.0;
      multiplier[iMarker]=0.0;
      multiplierhelp[iMarker]=0.0;
      multiplieroriginal[iMarker]=0.0;
      cons_factor[iMarker]=0.0;
  }

  Hess=new su2double*[38];
  Bess=new su2double*[38];


  for (iMarker = 0; iMarker < nMarker; iMarker++) {
      if(flow){
        CSensitivity[iMarker]        = new su2double [geometry->nVertex[iMarker]];
        CSensitivityOld[iMarker]        = new su2double [geometry->nVertex[iMarker]];
        LagrangeSens[iMarker]        = new su2double [geometry->nVertex[iMarker]];
        ExpCSensitivityOld[iMarker]        = new su2double [geometry->nVertex[iMarker]];
        ExpLagrangeSens[iMarker]        = new su2double [geometry->nVertex[iMarker]];
      }else{
          CSensitivity[iMarker]        = new su2double [geometry->GetnElem()];
          CSensitivityOld[iMarker]        = new su2double [geometry->GetnElem()];
          LagrangeSens[iMarker]        = new su2double [geometry->GetnElem()];
          ExpCSensitivityOld[iMarker]        = new su2double [geometry->GetnElem()];
          ExpLagrangeSens[iMarker]        = new su2double [geometry->GetnElem()];
      }

  }
  for (iMarker = 0; iMarker < config->GetLBFGSNum(); iMarker++) {
      rkStore[iMarker]=new su2double[38];
      dukStore[iMarker]=new su2double[38];
  }

  for (iVertex = 0; iVertex < 38; iVertex++) {
     Hess[iVertex]= new su2double [38];
     Bess[iVertex]= new su2double [38];
  }

  if(flow){
      UpdateSens=new su2double[geometry->nVertex[0]];
      for (iVertex = 0; iVertex < geometry->nVertex[0]; iVertex++) {
          UpdateSens[iVertex]=0;
      }
  }else{
      UpdateSens=new su2double[geometry->GetnElem()];
      for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {
          UpdateSens[iElem]=0;
      }
  }

  for (iVertex = 0; iVertex < 38; iVertex++) {
      for (jVertex = 0; jVertex < 38; jVertex++) {
         Hess[iVertex][jVertex]= 0.0;
         Bess[iVertex][jVertex]= 0.0;
      }
      Hess[iVertex][iVertex]=1.0;
      Bess[iVertex][iVertex]=1.0;
      if(config->GetHInit()){
          Hess[iVertex][iVertex]=config->GetHScale();
          Bess[iVertex][iVertex]=1.0/config->GetHScale();
      }
  }
  Sens_Geo  = new su2double[nMarker];
  Sens_Mach = new su2double[nMarker];
  Sens_AoA  = new su2double[nMarker];
  Sens_Press = new su2double[nMarker];
  Sens_Temp  = new su2double[nMarker];
  for (iMarker = 0; iMarker < nMarker; iMarker++) {
      Sens_Geo[iMarker]  = 0.0;
      Sens_Mach[iMarker] = 0.0;
      Sens_AoA[iMarker]  = 0.0;
      Sens_Press[iMarker] = 0.0;
      Sens_Temp[iMarker]  = 0.0;
      if(flow){
          for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++){
              CSensitivity[iMarker][iVertex] = 0.0;
              CSensitivityOld[iMarker][iVertex] = 0.0;
              ExpCSensitivityOld[iMarker][iVertex] = 0.0;
              ExpLagrangeSens[iMarker][iVertex] = 0.0;
              LagrangeSens[iMarker][iVertex] = 0.0;
          }
      }else{
          for (iElem = 0; iElem < geometry->GetnElem(); iElem++){
              CSensitivity[iMarker][iElem] = 0.0;
              CSensitivityOld[iMarker][iElem] = 0.0;
              ExpCSensitivityOld[iMarker][iElem] = 0.0;
              ExpLagrangeSens[iMarker][iElem] = 0.0;
              LagrangeSens[iMarker][iElem] = 0.0;
          }
      }
  }

  /*--- Check for a restart and set up the variables at each node
   appropriately. Coarse multigrid levels will be intitially set to
   the farfield values bc the solver will immediately interpolate
   the solution from the finest mesh to the coarser levels. ---*/
  if (!restart || (iMesh != MESH_0)) {

    /*--- Restart the solution from zero ---*/
    for (iPoint = 0; iPoint < nPoint; iPoint++)
      node[iPoint] = new CDiscAdjVariable(Solution, nDim, nVar, config);

  }
  else {

    /*--- Restart the solution from file information ---*/
    mesh_filename = config->GetSolution_AdjFileName();
    filename = config->GetObjFunc_Extension(mesh_filename);

    restart_file.open(filename.data(), ios::in);

    /*--- In case there is no file ---*/
    if (restart_file.fail()) {
      if (rank == MASTER_NODE)
        cout << "There is no adjoint restart file!! " << filename.data() << "."<< endl;
      exit(EXIT_FAILURE);
    }

    /*--- In case this is a parallel simulation, we need to perform the
     Global2Local index transformation first. ---*/
    long *Global2Local;
    Global2Local = new long[geometry->GetGlobal_nPointDomain()];
    /*--- First, set all indices to a negative value by default ---*/
    for (iPoint = 0; iPoint < geometry->GetGlobal_nPointDomain(); iPoint++) {
      Global2Local[iPoint] = -1;
    }
    /*--- Now fill array with the transform values only for local points ---*/
    for (iPoint = 0; iPoint < nPointDomain; iPoint++) {
      Global2Local[geometry->node[iPoint]->GetGlobalIndex()] = iPoint;
    }

    /*--- Read all lines in the restart file ---*/
    long iPoint_Local; unsigned long iPoint_Global = 0;\

    /*--- Skip coordinates ---*/
    unsigned short skipVars = nDim;

    /*--- Skip flow adjoint variables ---*/
    if (Kind_Solver == RUNTIME_TURB_SYS){
      if (compressible){
        skipVars += nDim + 2;
      }
      if (incompressible){
        skipVars += nDim + 1;
      }
    }

    /*--- The first line is the header ---*/
    getline (restart_file, text_line);

    while (getline (restart_file, text_line)) {
      istringstream point_line(text_line);

      /*--- Retrieve local index. If this node from the restart file lives
       on a different processor, the value of iPoint_Local will be -1.
       Otherwise, the local index for this node on the current processor
       will be returned and used to instantiate the vars. ---*/
      iPoint_Local = Global2Local[iPoint_Global];
      if (iPoint_Local >= 0) {
        point_line >> index;
        for (iVar = 0; iVar < skipVars; iVar++){ point_line >> dull_val;}
        for (iVar = 0; iVar < nVar; iVar++){ point_line >> Solution[iVar];}
        node[iPoint_Local] = new CDiscAdjVariable(Solution, nDim, nVar, config);
      }
      iPoint_Global++;
    }

    /*--- Instantiate the variable class with an arbitrary solution
     at any halo/periodic nodes. The initial solution can be arbitrary,
     because a send/recv is performed immediately in the solver. ---*/
    for (iPoint = nPointDomain; iPoint < nPoint; iPoint++) {
      node[iPoint] = new CDiscAdjVariable(Solution, nDim, nVar, config);
    }

    /*--- Close the restart file ---*/
    restart_file.close();

    /*--- Free memory needed for the transformation ---*/
    delete [] Global2Local;
  }

  /*--- Store the direct solution ---*/

  for (iPoint = 0; iPoint < nPoint; iPoint++){
    node[iPoint]->SetSolution_Direct(direct_solver->node[iPoint]->GetSolution());
  }
}

CDiscAdjSolver::~CDiscAdjSolver(void){ 

  unsigned short iMarker;

  if (CSensitivity != NULL) {
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
      delete [] CSensitivity[iMarker];
    }
    delete [] CSensitivity;
  }

  if (Sens_Geo   != NULL) delete [] Sens_Geo;
  if (Sens_Mach  != NULL) delete [] Sens_Mach;
  if (Sens_AoA   != NULL) delete [] Sens_AoA;
  if (Sens_Press != NULL) delete [] Sens_Press;
  if (Sens_Temp  != NULL) delete [] Sens_Temp;

}

void CDiscAdjSolver::SetRecording(CGeometry* geometry, CConfig *config, unsigned short kind_recording){


  bool time_n_needed  = ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
      (config->GetUnsteady_Simulation() == DT_STEPPING_2ND)),
  time_n1_needed = config->GetUnsteady_Simulation() == DT_STEPPING_2ND;

  unsigned long iPoint;
  unsigned short iVar;

  /*--- Reset the solution to the initial (converged) solution ---*/

  /*for (iPoint = 0; iPoint < nPoint; iPoint++){
    direct_solver->node[iPoint]->SetSolution(node[iPoint]->GetSolution_Direct());
  }*/ //ResetSolution

  if (time_n_needed){
    for (iPoint = 0; iPoint < nPoint; iPoint++){
      for (iVar = 0; iVar < nVar; iVar++){
        AD::ResetInput(direct_solver->node[iPoint]->GetSolution_time_n()[iVar]);
      }
    }
  }
  if (time_n1_needed){
    for (iPoint = 0; iPoint < nPoint; iPoint++){
      for (iVar = 0; iVar < nVar; iVar++){
        AD::ResetInput(direct_solver->node[iPoint]->GetSolution_time_n1()[iVar]);
      }
    }
  }

  /*--- Set the Jacobian to zero since this is not done inside the meanflow iteration
   * when running the discrete adjoint solver. ---*/

  direct_solver->Jacobian.SetValZero();

  /*--- Set indices to zero ---*/

  RegisterVariables(geometry, config, true);

}

void CDiscAdjSolver::RegisterSolution(CGeometry *geometry, CConfig *config){
  unsigned long iPoint, nPoint = geometry->GetnPoint();

  bool time_n_needed  = ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
      (config->GetUnsteady_Simulation() == DT_STEPPING_2ND)),
  time_n1_needed = config->GetUnsteady_Simulation() == DT_STEPPING_2ND,
  input = true;

  /*--- Register solution at all necessary time instances and other variables on the tape ---*/

  for (iPoint = 0; iPoint < nPoint; iPoint++){
    direct_solver->node[iPoint]->RegisterSolution(input);
  }
  if (time_n_needed){
    for (iPoint = 0; iPoint < nPoint; iPoint++){
      direct_solver->node[iPoint]->RegisterSolution_time_n();
    }
  }
  if (time_n1_needed){
    for (iPoint = 0; iPoint < nPoint; iPoint++){
      direct_solver->node[iPoint]->RegisterSolution_time_n1();
    }
  }
}

void CDiscAdjSolver::RegisterVariables(CGeometry *geometry, CConfig *config, bool reset){

  /*--- Register farfield values as input ---*/

  if((config->GetKind_Regime() == COMPRESSIBLE) && (KindDirect_Solver == RUNTIME_FLOW_SYS)){

    su2double Velocity_Ref = config->GetVelocity_Ref();
    Alpha                  = config->GetAoA()*PI_NUMBER/180.0;
    Beta                   = config->GetAoS()*PI_NUMBER/180.0;
    Mach                   = config->GetMach();
    Pressure               = config->GetPressure_FreeStreamND();
    Temperature            = config->GetTemperature_FreeStreamND();

    su2double SoundSpeed = 0.0;
    
    if (nDim == 2) { SoundSpeed = config->GetVelocity_FreeStreamND()[0]*Velocity_Ref/(cos(Alpha)*Mach); }
    if (nDim == 3) { SoundSpeed = config->GetVelocity_FreeStreamND()[0]*Velocity_Ref/(cos(Alpha)*cos(Beta)*Mach); }

    if (!reset){
      AD::RegisterInput(Mach);
      AD::RegisterInput(Alpha);
      AD::RegisterInput(Temperature);
      AD::RegisterInput(Pressure);
    }

    /*--- Recompute the free stream velocity ---*/

    if (nDim == 2) {
      config->GetVelocity_FreeStreamND()[0] = cos(Alpha)*Mach*SoundSpeed/Velocity_Ref;
      config->GetVelocity_FreeStreamND()[1] = sin(Alpha)*Mach*SoundSpeed/Velocity_Ref;
    }
    if (nDim == 3) {
      config->GetVelocity_FreeStreamND()[0] = cos(Alpha)*cos(Beta)*Mach*SoundSpeed/Velocity_Ref;
      config->GetVelocity_FreeStreamND()[1] = sin(Beta)*Mach*SoundSpeed/Velocity_Ref;
      config->GetVelocity_FreeStreamND()[2] = sin(Alpha)*cos(Beta)*Mach*SoundSpeed/Velocity_Ref;
    }

    config->SetTemperature_FreeStreamND(Temperature);
    direct_solver->SetTemperature_Inf(Temperature);
    config->SetPressure_FreeStreamND(Pressure);
    direct_solver->SetPressure_Inf(Pressure);

  }


    /*--- Here it is possible to register other variables as input that influence the flow solution
     * and thereby also the objective function. The adjoint values (i.e. the derivatives) can be
     * extracted in the ExtractAdjointVariables routine. ---*/
}

void CDiscAdjSolver::RegisterOutput(CGeometry *geometry, CConfig *config){

  unsigned long iPoint, nPoint = geometry->GetnPoint();

  /*--- Register variables as output of the solver iteration ---*/

  bool input = false;

  /*--- Register output variables on the tape ---*/

  for (iPoint = 0; iPoint < nPoint; iPoint++){
    direct_solver->node[iPoint]->RegisterSolution(input);
  }
}



void CDiscAdjSolver::RegisterObj_Func(CConfig *config){

  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  /*--- Here we can add new (scalar) objective functions ---*/

  switch (config->GetKind_ObjFunc()){
  case DRAG_COEFFICIENT:
      ObjFunc_Value = direct_solver->GetTotal_CDrag();
      break;
  case LIFT_COEFFICIENT:
      ObjFunc_Value = direct_solver->GetTotal_CLift();
      break;
  case SIDEFORCE_COEFFICIENT:
      ObjFunc_Value = direct_solver->GetTotal_CSideForce();
      break;
  case EFFICIENCY:
      ObjFunc_Value = direct_solver->GetTotal_CEff();
      break;
  case MOMENT_X_COEFFICIENT:
      ObjFunc_Value = direct_solver->GetTotal_CMx();
      break;
  case MOMENT_Y_COEFFICIENT:
      ObjFunc_Value = direct_solver->GetTotal_CMy();
      break;
  case MOMENT_Z_COEFFICIENT:
      ObjFunc_Value = direct_solver->GetTotal_CMz();
      break;
  case EQUIVALENT_AREA:
      ObjFunc_Value = direct_solver->GetTotal_CEquivArea();
      break;
  case AVG_TOTAL_PRESSURE:
    ObjFunc_Value = direct_solver->GetOneD_TotalPress();
    break;
  case AVG_OUTLET_PRESSURE:
    ObjFunc_Value = direct_solver->GetOneD_FluxAvgPress();
    break;
  case MASS_FLOW_RATE:
    ObjFunc_Value = direct_solver->GetOneD_MassFlowRate();
    break;
  case MINIMUM_COMPLIANCE:
    ObjFunc_Value = direct_solver->GetMinimumCompliance();
    break;
 /*--- Template for new objective functions where TemplateObjFunction()
  *  is the routine that returns the obj. function value. The computation
  * must be done while the tape is active, i.e. between AD::StartRecording() and
  * AD::StopRecording() in DiscAdjMeanFlowIteration::Iterate(). The best place is somewhere
  * inside MeanFlowIteration::Iterate().
  *
  * case TEMPLATE_OBJECTIVE:
  *    ObjFunc_Value = TemplateObjFunction();
  *    break;
  * ---*/
  }
  ObjFunc_Value=config->GetScaleObj()*ObjFunc_Value;
  if (rank == MASTER_NODE){
    AD::RegisterOutput(ObjFunc_Value);
  }
}

void CDiscAdjSolver::RegisterConstraint_Func(CConfig *config, CGeometry *geometry){

  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  if(config->GetPosConstraint()){
    ConstraintFunc_Value[0] = config->GetTargetLift()-direct_solver->GetTotal_CLift();
  }else{
    ConstraintFunc_Value[0] = -config->GetTargetLift()+direct_solver->GetTotal_CLift(); //ATTENTION was set before
  }
  ConstraintFunc_Value[0]=config->GetScaleConstr()*ConstraintFunc_Value[0];
  if (rank == MASTER_NODE){
    AD::RegisterOutput(ConstraintFunc_Value[0]);
  }

  if(config->GetConstraintNum()>1){
      if(config->GetPosConstraint()){
        ConstraintFunc_Value[1] = -direct_solver->GetTotal_CMz();
      }else{
        ConstraintFunc_Value[1] = direct_solver->GetTotal_CMz();
      }
      ConstraintFunc_Value[1]=config->GetScaleConstr()*ConstraintFunc_Value[1];
      if (rank == MASTER_NODE){
        AD::RegisterOutput(ConstraintFunc_Value[1]);
      }
  }

  if(config->GetConstraintNum()>2){
      su2double *Plane_P0 = new su2double[3];
      su2double *Plane_Normal = new su2double[3];
      vector<su2double> *Xcoord_Airfoil = new vector<su2double>[1];
      vector<su2double> *Ycoord_Airfoil = new vector<su2double>[1];
      vector<su2double> *Zcoord_Airfoil = new vector<su2double>[1];
      vector<su2double> *Variable_Airfoil = new vector<su2double>[1];
      su2double MinXCoord, MaxXCoord;

     if (geometry->GetnDim() == 2) {
        MinXCoord = -1E6;
        MaxXCoord = 1E6;
        Plane_Normal[0] = 0.0;   Plane_P0[0] = 0.0;
        Plane_Normal[1] = 1.0;   Plane_P0[1] = 0.0;
        Plane_Normal[2] = 0.0;   Plane_P0[2] = 0.0;
      }

     geometry->ComputeAirfoil_Section(Plane_P0, Plane_Normal, MinXCoord, MaxXCoord, NULL,
                                         Xcoord_Airfoil[0], Ycoord_Airfoil[0], Zcoord_Airfoil[0], Variable_Airfoil[0], true, config);
      if(config->GetPosConstraint()){
        ConstraintFunc_Value[2] = -geometry->Compute_MaxThickness(Plane_P0, Plane_Normal, 0, config, Xcoord_Airfoil[0], Ycoord_Airfoil[0], Zcoord_Airfoil[0], true)+0.12;
      }else{
        ConstraintFunc_Value[2] = geometry->Compute_MaxThickness(Plane_P0, Plane_Normal, 0, config, Xcoord_Airfoil[0], Ycoord_Airfoil[0], Zcoord_Airfoil[0], true)-0.12;
      }

      ConstraintFunc_Value[2]=config->GetScaleConstr()*ConstraintFunc_Value[2];
      if (rank == MASTER_NODE){
        AD::RegisterOutput(ConstraintFunc_Value[2]);
      }

      delete [] Xcoord_Airfoil;
      delete [] Ycoord_Airfoil;
      delete [] Zcoord_Airfoil;
      delete [] Variable_Airfoil;
      delete [] Plane_P0;
      delete [] Plane_Normal;
  }
}

void CDiscAdjSolver::SetAdj_ObjFunc(CGeometry *geometry, CConfig *config, double initVal){
  int rank = MASTER_NODE;

  bool time_stepping = config->GetUnsteady_Simulation() != STEADY;
  unsigned long IterAvg_Obj = config->GetIter_Avg_Objective();
  unsigned long ExtIter = config->GetExtIter();
  su2double seeding = 1.0;

  if (time_stepping){
    if (ExtIter < IterAvg_Obj){
      seeding = 1.0/((su2double)IterAvg_Obj);
    }
    else{
      seeding = 0.0;
    }
  }

#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  if (rank == MASTER_NODE){
    SU2_TYPE::SetDerivative(ObjFunc_Value, initVal);
  } else {
    SU2_TYPE::SetDerivative(ObjFunc_Value, 0.0);
  }
}

void CDiscAdjSolver::SetAdj_ConstraintFuncAD(CGeometry *geometry, CConfig *config, su2double* initVal){
  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  for (unsigned short iValue=0; iValue<config->GetConstraintNum();iValue++){
      if (rank == MASTER_NODE){
        SU2_TYPE::SetDerivative(ConstraintFunc_Value[iValue], SU2_TYPE::GetValue(initVal[iValue]));
      } else {
        SU2_TYPE::SetDerivative(ConstraintFunc_Value[iValue], 0.0);
      }
  }
}

void CDiscAdjSolver::SetAdj_ConstraintFunc(CGeometry *geometry, CConfig *config, double* initVal){
  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  for (unsigned short iValue=0; iValue<config->GetConstraintNum();iValue++){
      if (rank == MASTER_NODE){
        SU2_TYPE::SetDerivative(ConstraintFunc_Value[iValue], initVal[iValue]);
      } else {
        SU2_TYPE::SetDerivative(ConstraintFunc_Value[iValue], 0.0);
      }
  }
}

su2double *CDiscAdjSolver::GetConstraintFunc_Value(){
    return ConstraintFunc_Value;
}

su2double CDiscAdjSolver::GetConstraint_Save(unsigned short indCons){
    return Constraint_Save[indCons];
}

su2double CDiscAdjSolver::GetObj_Save(){
    return Obj_Save;
}


void CDiscAdjSolver::StoreConstraint(CConfig *config){
   for (unsigned short iValue=0; iValue<config->GetConstraintNum();iValue++){
        Constraint_Old[iValue]=Constraint_Save[iValue];
        Constraint_Save[iValue]= ConstraintFunc_Value[iValue];
   }
}

double* CDiscAdjSolver::GetMultiplier(){
    return multiplier;
}

void CDiscAdjSolver::SetMultiplier(CConfig *config, double * value){
    unsigned short iHelp;
    for (iHelp=0; iHelp<config->GetConstraintNum();iHelp++){
        multiplier[iHelp]=value[iHelp];
        multiplierhelp[iHelp]=multiplier[iHelp];
        multiplieroriginal[iHelp]=multiplier[iHelp];
        cons_factor[iHelp]=SU2_TYPE::GetValue(config->GetConstraintFactor(iHelp));
    }

}

void CDiscAdjSolver::StoreMeshPoints(CConfig *config, CGeometry *geometry){
    unsigned long iVertex, jPoint;
    unsigned short iMarker;
    for (jPoint=0; jPoint<geometry->GetnPoint();jPoint++){
        geometry->node[jPoint]->SetCoord_Old(geometry->node[jPoint]->GetCoord());
    }
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
        for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          jPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          //geometry->node[jPoint]->SetCoord_Old(geometry->node[jPoint]->GetCoord());
          geometry->vertex[iMarker][iVertex]->SetNormal_Old(geometry->vertex[iMarker][iVertex]->GetNormal());
        }
    }
}

void CDiscAdjSolver::LoadMeshPoints(CConfig *config, CGeometry *geometry){
    unsigned long iVertex, jPoint;
    unsigned short iMarker;
    /*for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
        for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          jPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          geometry->node[jPoint]->SetCoord(geometry->node[jPoint]->GetCoord_Old());
          geometry->vertex[iMarker][iVertex]->SetNormal(geometry->vertex[iMarker][iVertex]->GetNormal_Old());
        }
    }*/
    for (jPoint=0; jPoint<geometry->GetnPoint();jPoint++){
        geometry->node[jPoint]->SetCoord(geometry->node[jPoint]->GetCoord_Old());
    }

}

void CDiscAdjSolver::DeformSurface(CConfig *config, CGeometry *geometry, su2double *rand, CVolumetricMovement *grid_movement){
    //surface_movement->SetExternal_Deformation(geometry, config, ZONE_0, 0);

      unsigned short iDim, nDim;
      su2double VarCoord[3];
      su2double Lref   = config->GetLength_Ref();
      unsigned long iVertex;
      unsigned short iMarker;

      nDim = geometry->GetnDim();
      unsigned long calcPoints=0;

      for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
        //if (config->GetMarker_All_Moving(iMarker) == YES) {
          if (config->GetMarker_All_DV(iMarker) == YES) {
          for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
            for (iDim = 0; iDim < nDim; iDim++){
              VarCoord[iDim]=0.0;
              if(calcPoints>=offPoints && calcPoints<numPoints-offPoints){
                    for (unsigned short ibase=0; ibase<3; ibase++){
                        VarCoord[iDim] += (sqrt(ewbasis[ibase])*evbasis[ibase][calcPoints-offPoints]*rand[ibase]*(geometry->vertex[iMarker][iVertex]->GetNormal_Old()[iDim]/lenNormal[calcPoints])*sqrt(0.01))/Lref; //sqrt(0.1)
                    }
              }
            }
            if (nDim == 2) VarCoord[nDim] = 0.0;
            calcPoints++;
            geometry->vertex[iMarker][iVertex]->SetVarCoord(VarCoord);

          }
        }
      }
      grid_movement->SetVolume_Deformation(geometry, config, true);
}

void CDiscAdjSolver::ComputeEVEW(CConfig *config, CGeometry *geometry){
    //compute the covariance matrix and its eigenvectors and eigenvalues
    su2double *Coords;// = new su2double[geometry->GetnDim()];
    su2double *Normal = new su2double[geometry->GetnDim()];
  //  su2double *NormalTest;
    unsigned long iMarker, iVertex, jPoint;
    numPoints=0;
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      if (config->GetMarker_All_DV(iMarker) == YES) {
        for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
            numPoints++;
        }
      }
    }
    su2double *x =new su2double[numPoints];
    su2double *y =new su2double[numPoints];
    numPoints=0;
    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      if (config->GetMarker_All_DV(iMarker) == YES) {
        for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {
          jPoint = geometry->vertex[iMarker][iVertex]->GetNode();
          Coords = geometry->node[jPoint]->GetCoord();

          geometry->vertex[iMarker][iVertex]->GetNormal(Normal);
     //     NormalTest=geometry->vertex[iMarker][iVertex]->GetNormal();
          lenNormal[numPoints]=0;
          for (unsigned short iDim = 0; iDim < geometry->GetnDim(); iDim++) {
              lenNormal[numPoints]+=(Normal)[iDim]*(Normal)[iDim];
          }
          for (unsigned short iDim = 0; iDim < geometry->GetnDim(); iDim++) {
              Normal[iDim]=Normal[iDim]/sqrt(lenNormal[numPoints]);

          }
          x[numPoints]=Coords[0];
          y[numPoints]=Coords[1];
          numPoints++;
        }
      }
    }
//    delete [] Coords;
    delete [] Normal;
    su2double covEntry;
    unsigned short icount, jcount;
    offPoints=29;
    Eigen::MatrixXd Cov(numPoints-2*offPoints,numPoints-2*offPoints);
    for (icount=offPoints; icount<numPoints-offPoints;icount++){
        for (jcount=offPoints; jcount<numPoints-offPoints; jcount++){
            covEntry=0.001*0.001*exp(-(1./(0.1*0.1))*((x[icount]-x[jcount])*(x[icount]-x[jcount])+(y[icount]-y[jcount])*(y[icount]-y[jcount])));
            Cov(icount-offPoints,jcount-offPoints)=SU2_TYPE::GetValue(covEntry);
        }
    }

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(Cov);

    double* evalues = new double[numPoints-2*offPoints];
    double** evectors = new double*[numPoints-2*offPoints];
    for (icount=0; icount<numPoints-2*offPoints;icount++){
        evalues[icount]=eigensolver.eigenvalues()[icount];
        evectors[icount]=new double[numPoints-2*offPoints];
        for (jcount=0; jcount<numPoints-2*offPoints;jcount++){
            evectors[icount][jcount]=eigensolver.eigenvectors().col(icount)(jcount);
        }
    }

    // find largest eigenvalues and corresponding evectors -> store in ewbasis and evbasis
    int numMaxEw=3;
    bool used;
    int* maxInd=new int[numMaxEw];
    for(unsigned short kcount=0;kcount<numMaxEw;kcount++){
        maxInd[kcount]=0;
        for (icount=0; icount<numPoints-2*offPoints;icount++){
            if(evalues[maxInd[kcount]]<=evalues[icount]){
                used = false;
                for (jcount=0; jcount<kcount;jcount++){
                    if(maxInd[jcount]==icount) used =true;
                }
                if (!used)  maxInd[kcount]=icount;
            }
        }
        ewbasis[kcount]=evalues[maxInd[kcount]];
        for (icount=0; icount<numPoints-2*offPoints;icount++){
            evbasis[kcount][icount]=evectors[maxInd[kcount]][icount];
        }
    }

    delete [] x;
    delete [] y;
    delete [] evalues;
    for (icount=0; icount<numPoints-2*offPoints;icount++){
        delete [] evectors[icount];
    }
    delete [] evectors;
    delete [] maxInd;


}

void CDiscAdjSolver::UpdateMultiplier(CConfig *config){
    //use this for +h and inequality
 /*   if(Constraint_Save>1E-3){ //neu: 6
            multiplier=multiplierhelp*(1+SU2_TYPE::GetPrimary(config->GetConstraintFactor())*SU2_TYPE::GetPrimary(Constraint_Save));
         //   multiplier=multiplierhelp-SU2_TYPE::GetPrimary(config->GetConstraintFactor())*SU2_TYPE::GetPrimary(Constraint_Save);
            multiplierhelp=multiplier;
    }else if(Constraint_Save<-1E-3 && config->GetEqualConstraint()){ //vorher 1E-3
            multiplier=-multiplierhelp*(1-SU2_TYPE::GetPrimary(config->GetConstraintFactor())*SU2_TYPE::GetPrimary(Constraint_Save)); //a
        //    multiplier=multiplierhelp-SU2_TYPE::GetPrimary(config->GetConstraintFactor())*SU2_TYPE::GetPrimary(Constraint_Save); //b
         //   multiplierhelp=multiplier; //a,b
            multiplierhelp=-multiplier; //c+a
    }else{
            multiplier=0;
            multiplierhelp=multiplieroriginal;
            std::cout<<"Multiplier set to 0"<<std::endl;
    }*/
    for (unsigned short iValue=0; iValue<config->GetConstraintNum();iValue++){
        if(config->GetPosConstraint()){
            if(Constraint_Save[iValue]<=0 && !(config->GetEqualConstraint(iValue))){
                    multiplier[iValue]=0;
                    multiplierhelp[iValue]=multiplieroriginal[iValue];
            }else{
                if(config->GetFactorIncrease()) cons_factor[iValue]=cons_factor[iValue]*(1+fabs(SU2_TYPE::GetValue(Constraint_Save[iValue]))/fabs(SU2_TYPE::GetValue(Constraint_Save[iValue])-SU2_TYPE::GetValue(Constraint_Old[iValue])));
                if(config->GetPosUpdate()){
                    multiplier[iValue]=multiplierhelp[iValue]+cons_factor[iValue]*SU2_TYPE::GetValue(Constraint_Save[iValue]);
                }else{
                    multiplier[iValue]=multiplierhelp[iValue]-cons_factor[iValue]*SU2_TYPE::GetValue(Constraint_Save[iValue]);
                }
                multiplierhelp[iValue]=multiplier[iValue];
            }
        }else{
            if(Constraint_Save[iValue]>=0 && !(config->GetEqualConstraint(iValue))){
                    multiplier[iValue]=0;
                    multiplierhelp[iValue]=multiplieroriginal[iValue];
            }else{
                if(config->GetFactorIncrease()) cons_factor[iValue]=cons_factor[iValue]*(1+fabs(SU2_TYPE::GetValue(Constraint_Save[iValue]))/fabs(SU2_TYPE::GetValue(Constraint_Save[iValue])-SU2_TYPE::GetValue(Constraint_Old[iValue])));
                if(config->GetPosUpdate()){
                    multiplier[iValue]=multiplierhelp[iValue]+cons_factor[iValue]*SU2_TYPE::GetValue(Constraint_Save[iValue]);
                }else{
                    multiplier[iValue]=multiplierhelp[iValue]-cons_factor[iValue]*SU2_TYPE::GetValue(Constraint_Save[iValue]);
                }
                multiplierhelp[iValue]=multiplier[iValue];
            }
        }
        std::cout<<"Update of Multiplier: "<<multiplier[iValue]<<" "<<cons_factor[iValue]<<std::endl;
        //cons_factor=cons_factor*1.1;
    }
}

void CDiscAdjSolver::ExtractAdjoint_Solution(CGeometry *geometry, CConfig *config){

  bool time_n_needed  = ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
      (config->GetUnsteady_Simulation() == DT_STEPPING_2ND));

  bool time_n1_needed = config->GetUnsteady_Simulation() == DT_STEPPING_2ND;

  unsigned short iVar;
  unsigned long iPoint;
  su2double residual;

  /*--- Set Residuals to zero ---*/

  for (iVar = 0; iVar < nVar; iVar++){
      SetRes_RMS(iVar,0.0);
      SetRes_Max(iVar,0.0,0);
  }

  for (iPoint = 0; iPoint < nPoint; iPoint++){

    /*--- Set the old solution ---*/

    node[iPoint]->Set_OldSolution();

    /*--- Extract the adjoint solution ---*/

    direct_solver->node[iPoint]->GetAdjointSolution(Solution);

    /*--- Store the adjoint solution ---*/

    node[iPoint]->SetSolution(Solution);
  }

  if (time_n_needed){
    for (iPoint = 0; iPoint < nPoint; iPoint++){

      /*--- Extract the adjoint solution at time n ---*/

      direct_solver->node[iPoint]->GetAdjointSolution_time_n(Solution);

      /*--- Store the adjoint solution at time n ---*/

      node[iPoint]->Set_Solution_time_n(Solution);
    }
  }
  if (time_n1_needed){
    for (iPoint = 0; iPoint < nPoint; iPoint++){

      /*--- Extract the adjoint solution at time n-1 ---*/

      direct_solver->node[iPoint]->GetAdjointSolution_time_n1(Solution);

      /*--- Store the adjoint solution at time n-1 ---*/

      node[iPoint]->Set_Solution_time_n1(Solution);
    }
  }

  /*--- Set the residuals ---*/

  for (iPoint = 0; iPoint < nPointDomain; iPoint++){
      for (iVar = 0; iVar < nVar; iVar++){
          residual = node[iPoint]->GetSolution(iVar) - node[iPoint]->GetSolution_Old(iVar);

          AddRes_RMS(iVar,residual*residual);
          AddRes_Max(iVar,fabs(residual),geometry->node[iPoint]->GetGlobalIndex(),geometry->node[iPoint]->GetCoord());
      }
  }

  SetResidual_RMS(geometry, config);
}

void CDiscAdjSolver::SetAdjointInputHelp(CGeometry *geometry, CConfig *config){ //do not overwrite oldSolution

  bool time_n_needed  = ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
      (config->GetUnsteady_Simulation() == DT_STEPPING_2ND));

  bool time_n1_needed = config->GetUnsteady_Simulation() == DT_STEPPING_2ND;


  unsigned short iVar;
  unsigned long iPoint;
  su2double residual;

  /*--- Set Residuals to zero --- */

  for (iVar = 0; iVar < nVar; iVar++){
      SetRes_RMS(iVar,0.0);
      SetRes_Max(iVar,0.0,0);
  }

  for (iPoint = 0; iPoint < nPoint; iPoint++){

    /*--- Extract the adjoint solution ---*/

    direct_solver->node[iPoint]->GetAdjointSolution(Solution);

    /*--- Store the adjoint solution ---*/

    node[iPoint]->SetSolution(Solution);
  }

  if (time_n_needed){
    for (iPoint = 0; iPoint < nPoint; iPoint++){

      /*--- Extract the adjoint solution at time n ---*/

      direct_solver->node[iPoint]->GetAdjointSolution_time_n(Solution);

      /*--- Store the adjoint solution at time n ---*/

      node[iPoint]->Set_Solution_time_n(Solution);
    }
  }
  if (time_n1_needed){
    for (iPoint = 0; iPoint < nPoint; iPoint++){

      /*--- Extract the adjoint solution at time n-1 ---*/

      direct_solver->node[iPoint]->GetAdjointSolution_time_n1(Solution);

      /*--- Store the adjoint solution at time n-1 ---*/

      node[iPoint]->Set_Solution_time_n1(Solution);
    }
  }

  /* --- Set the residuals --- */

  for (iPoint = 0; iPoint < nPointDomain; iPoint++){
      for (iVar = 0; iVar < nVar; iVar++){
          residual = node[iPoint]->GetSolution(iVar) - node[iPoint]->GetSolution_Old(iVar);

          AddRes_RMS(iVar,residual*residual);
          AddRes_Max(iVar,fabs(residual),geometry->node[iPoint]->GetGlobalIndex(),geometry->node[iPoint]->GetCoord());
      }
  }

  SetResidual_RMS(geometry, config);
}



void CDiscAdjSolver::StoreOldSolution(){
    unsigned long iPoint;
    for (iPoint = 0; iPoint < nPoint; iPoint++){
     // direct_solver->node[iPoint]->Set_OldSolution();
     // node[iPoint]->Set_OldSolution();
      direct_solver->node[iPoint]->Set_StoreSolution();
      node[iPoint]->Set_StoreSolution();
    }
}

void CDiscAdjSolver::LoadOldSolution(){
    TotalIterations=TotalIterations+1;
    unsigned long iPoint;
    for (iPoint = 0; iPoint < nPoint; iPoint++){
      //direct_solver->node[iPoint]->SetSolution(direct_solver->node[iPoint]->GetSolution_Old());
      //node[iPoint]->SetSolution(node[iPoint]->GetSolution_Old());
      direct_solver->node[iPoint]->SetSolution(direct_solver->node[iPoint]->GetSolution_Store());
      node[iPoint]->SetSolution(node[iPoint]->GetSolution_Store());
    }
}

void CDiscAdjSolver::StoreSolutionVec(unsigned short numQuad){
    unsigned long iPoint;
    for (iPoint = 0; iPoint < nPoint; iPoint++){
      direct_solver->node[iPoint]->SetSolutionVec(numQuad);
      node[iPoint]->SetSolutionVec(numQuad);
    }
}

void CDiscAdjSolver::LoadSolutionVec(unsigned short numQuad){
    unsigned long iPoint;
    for (iPoint = 0; iPoint < nPoint; iPoint++){
      direct_solver->node[iPoint]->SetSolution(direct_solver->node[iPoint]->GetSolutionVec(true,numQuad));
      node[iPoint]->SetSolution(node[iPoint]->GetSolutionVec(true,numQuad));
    }
}

void CDiscAdjSolver::StoreSolutionVecOld(unsigned short numQuad){
    unsigned long iPoint;
    for (iPoint = 0; iPoint < nPoint; iPoint++){
      direct_solver->node[iPoint]->SetSolutionVecOld(numQuad);
      node[iPoint]->SetSolutionVecOld(numQuad);
    }
}

void CDiscAdjSolver::LoadSolutionVecOld(unsigned short numQuad){
    unsigned long iPoint;
    for (iPoint = 0; iPoint < nPoint; iPoint++){
      direct_solver->node[iPoint]->SetSolution(direct_solver->node[iPoint]->GetSolutionVecOld(true,numQuad));
      node[iPoint]->SetSolution(node[iPoint]->GetSolutionVecOld(true,numQuad));
    }
}

void CDiscAdjSolver::LoadOldAdjoint(){
    unsigned long iPoint;
    unsigned short iVar;
    for (iPoint = 0; iPoint < nPoint; iPoint++){
            //node[iPoint]->SetSolution(node[iPoint]->GetSolution_Old());
        node[iPoint]->SetSolution(node[iPoint]->GetSolution_Store());
    }
}

void CDiscAdjSolver::StoreSaveSolution(){
    unsigned long iPoint;
    for (iPoint = 0; iPoint < nPoint; iPoint++){
      node[iPoint]->Set_SaveSolution();
      direct_solver->node[iPoint]->Set_SaveSolution();
    }
}

void CDiscAdjSolver::LoadSaveSolution(){
    unsigned long iPoint;
    for (iPoint = 0; iPoint < nPoint; iPoint++){
            node[iPoint]->SetSolution(node[iPoint]->GetSolution_Save());
            direct_solver->node[iPoint]->SetSolution(direct_solver->node[iPoint]->GetSolution_Save());
    }
}

void CDiscAdjSolver::OutputWritten(CGeometry *geometry){
    unsigned short iVar;
    unsigned long iPoint;
    unsigned long iVertex;
 /*   std::cout<<"fdot "<<SU2_TYPE::GetForwardDerivative(ObjFunc_Value)<<std::endl;
    std::cout<<"ydot "<<std::endl;
    for (iPoint = 0; iPoint < nPoint; iPoint++){
      for (iVar = 0; iVar < nVar; iVar++){
        std::cout<<SU2_TYPE::GetForwardDerivative(direct_solver->node[iPoint]->GetSolution(iVar))<<" ";
      }
    }*/
    //for (unsigned short numQuad=0; numQuad<4; numQuad++){
     //   std::cout<<"SolutionVec_Old["<<numQuad<<"]=";
    //    for (iPoint = 0; iPoint < nPoint; iPoint++){
    //      for (iVar=0; iVar<nVar; iVar++){
    //          std::cout<<direct_solver->node[iPoint]->GetSolutionVec(true,numQuad)[iVar]<<" ";
    //      }
          //node[iPoint]->SetSolution(node[iPoint]->GetSolutionVec(true,numQuad));
    //    }
   //     std::cout<<std::endl;
   // }
    /*---- DUMMY FUNCTION for writing arrays (call: solver_container[iZone][iMesh][ADJFLOW_SOL]->OutputWritten();) ----*/
    /*su2double norm=0;
    for (iPoint = 0; iPoint < nPoint; iPoint++){
      for (iVar = 0; iVar < nVar; iVar++){
        Solution[iVar] = 0.0;
        norm+=Solution[iVar]*Solution[iVar];
      }
      direct_solver->node[iPoint]->SetAdjointSolution(Solution);
    }
    std::cout<<"SetAdjointUpdate: "<<sqrt(norm)<<std::endl;*/
//    std::cout<<"CSensitivity: "<<std::endl;
    //for (iVertex = 0; iVertex < geometry->GetnVertex(0); iVertex++){
//    for (iVertex = 0; iVertex < geometry->GetnElem(); iVertex++){
 //     std::cout<<CSensitivity[0][iVertex]<<" ";
//    }
//    std::cout<<std::endl;
//    std::cout<<"Objective: "<<ObjFunc_Value<<std::endl;
    //std::cout<<"SolutionOld - direct_solver"<<std::endl;
    //for (iPoint = 0; iPoint < nPoint; iPoint++){
   //     for (iVar = 0; iVar < nVar; iVar++){
            //std::cout<<direct_solver->node[iPoint]->GetSolution(iVar)<<" ";
   //        std::cout<<direct_solver->node[0]->GetSolution(0)<<" ";

   //     }
   // }
   // std::cout<<std::endl;
    //std::cout<<"SolutionOld"<<std::endl;
//    for (iPoint = 0; iPoint < nPoint; iPoint++){
//        std::cout<<geometry->node[iPoint]->GetCoord(0)<<" "<<geometry->node[iPoint]->GetCoord(1)<<" "<<direct_solver->node[iPoint]->GetSolution(0)<<" "<<direct_solver->node[iPoint]->GetSolution(1)<<std::endl;
        //for (iVar = 0; iVar < nVar; iVar++){
            //std::cout<<node[iPoint]->GetSolution(iVar)<<" ";
    //std::cout<<node[0]->GetSolution(0)<<" ";

   //     }
//    }
 //   std::cout<<std::endl;
}

void CDiscAdjSolver::AssembleLagrangian(CConfig *config){
    unsigned short iVar;
    unsigned long iPoint;
    Lagrangian_Value=0.0;
    su2double helper=0.0;
    for (iPoint = 0; iPoint < nPoint; iPoint++){
      for (iVar = 0; iVar < nVar; iVar++){
        //helper+=(direct_solver->node[iPoint]->GetSolution(iVar)-direct_solver->node[iPoint]->GetSolution_Old(iVar))*(direct_solver->node[iPoint]->GetSolution(iVar)-direct_solver->node[iPoint]->GetSolution_Old(iVar));
        helper+=(direct_solver->node[iPoint]->GetSolution(iVar)-direct_solver->node[iPoint]->GetSolution_Store(iVar))*(direct_solver->node[iPoint]->GetSolution(iVar)-direct_solver->node[iPoint]->GetSolution_Store(iVar));
      }
    }
    if(config->GetOneShotConstraint()==true){
        for (unsigned short iValue=0; iValue<config->GetConstraintNum();iValue++){
            helper+=ConstraintFunc_Value[iValue]*ConstraintFunc_Value[iValue];
        }
        helper=sqrt(helper/(nPoint*nVar+config->GetConstraintNum()))*(config->GetOneShotAlpha()/2);
    }else{
        helper=sqrt(helper/(nPoint*nVar))*(config->GetOneShotAlpha()/2);
    }
    Lagrangian_Value+=helper;
    helper=0.0;
    for (iPoint = 0; iPoint < nPoint; iPoint++){
      for (iVar = 0; iVar < nVar; iVar++){
        //helper+=(node[iPoint]->GetSolution(iVar)-node[iPoint]->GetSolution_Old(iVar))*(node[iPoint]->GetSolution(iVar)-node[iPoint]->GetSolution_Old(iVar));
        helper+=(node[iPoint]->GetSolution(iVar)-node[iPoint]->GetSolution_Store(iVar))*(node[iPoint]->GetSolution(iVar)-node[iPoint]->GetSolution_Store(iVar));
      }
    }
    Lagrangian_Value+=sqrt(helper/(nPoint*nVar))*(config->GetOneShotBeta()/2);
    Lagrangian_Value+=ObjFunc_Value;
    Obj_Save=ObjFunc_Value;
    if(config->GetOneShotConstraint()==true){
        for (unsigned short iValue=0; iValue<config->GetConstraintNum();iValue++){
            Lagrangian_Value+=ConstraintFunc_Value[iValue]*multiplier[iValue];
        }
    }
    helper=0.0;
    for (iPoint = 0; iPoint < nPoint; iPoint++){
      for (iVar = 0; iVar < nVar; iVar++){
        //helper+=(direct_solver->node[iPoint]->GetSolution(iVar)-direct_solver->node[iPoint]->GetSolution_Old(iVar))*node[iPoint]->GetSolution_Old(iVar);
        helper+=(direct_solver->node[iPoint]->GetSolution(iVar)-direct_solver->node[iPoint]->GetSolution_Store(iVar))*node[iPoint]->GetSolution_Store(iVar);
      }
    }
    Lagrangian_Value+=helper;
}

void CDiscAdjSolver::UpdateStateVariable(CConfig *config){
    unsigned long iPoint;
    unsigned short iVar;
    su2double stepsize=config->GetFDStep();
    for (iPoint = 0; iPoint < nPoint; iPoint++){
      for (iVar = 0; iVar < nVar; iVar++){
        //Solution[iVar] = direct_solver->node[iPoint]->GetSolution_Old(iVar)+stepsize*(node[iPoint]->GetSolution(iVar)-node[iPoint]->GetSolution_Old(iVar));
        Solution[iVar] = direct_solver->node[iPoint]->GetSolution_Store(iVar)+stepsize*(node[iPoint]->GetSolution(iVar)-node[iPoint]->GetSolution_Store(iVar));
      }
      direct_solver->node[iPoint]->SetSolution(Solution);
    }
}

void CDiscAdjSolver::SetForwardDirection(CConfig *config){
    unsigned short iVar;
    unsigned long iPoint;
   // std::cout<<"ForwardDirection"<<std::endl;
    for (iPoint = 0; iPoint < nPoint; iPoint++){
      for (iVar = 0; iVar < nVar; iVar++){
        Solution[iVar] = (node[iPoint]->GetSolution_Save(iVar)-node[iPoint]->GetSolution_Store(iVar));
         // Solution[iVar] = 1.0;
       // std::cout<<Solution[iVar]<<" ";
      }
      direct_solver->node[iPoint]->SetForwardSolution(Solution);
    }
   // std::cout<<std::endl;
}

void CDiscAdjSolver::SetAdjointOutputUpdate(CGeometry *geometry, CConfig *config){
    unsigned short iVar;
    unsigned long iPoint;
    normy=0;
    for (iPoint = 0; iPoint < nPoint; iPoint++){
      for (iVar = 0; iVar < nVar; iVar++){
        //Solution[iVar] = (direct_solver->node[iPoint]->GetSolution(iVar)-direct_solver->node[iPoint]->GetSolution_Old(iVar));
        Solution[iVar] = (direct_solver->node[iPoint]->GetSolution(iVar)-direct_solver->node[iPoint]->GetSolution_Store(iVar));
        normy+=Solution[iVar]*Solution[iVar];
      }
      direct_solver->node[iPoint]->SetAdjointSolution(Solution);
    }
    normy=sqrt(normy/(nPoint*nVar));
    if((normy/normyold) >0.01*rho){
        rho=(normy/normyold);
    }else{
        rho=0.01*rho;
    }
    normyold=normy;
}

void CDiscAdjSolver::SetAdjointOutputZero(CGeometry *geometry, CConfig *config){
    unsigned short iVar;
    unsigned long iPoint;
    for (iPoint = 0; iPoint < nPoint; iPoint++){
      for (iVar = 0; iVar < nVar; iVar++){
        Solution[iVar] = 0.0;
      }
      direct_solver->node[iPoint]->SetAdjointSolution(Solution);
    }
}

void CDiscAdjSolver::ExtractAdjoint_Variables(CGeometry *geometry, CConfig *config){
  su2double Local_Sens_Press, Local_Sens_Temp, Local_Sens_AoA, Local_Sens_Mach;

  /*--- Extract the adjoint values of the farfield values ---*/

  if ((config->GetKind_Regime() == COMPRESSIBLE) && (KindDirect_Solver == RUNTIME_FLOW_SYS)){
    Local_Sens_Mach  = SU2_TYPE::GetDerivative(Mach);
    Local_Sens_AoA   = SU2_TYPE::GetDerivative(Alpha);
    Local_Sens_Temp  = SU2_TYPE::GetDerivative(Temperature);
    Local_Sens_Press = SU2_TYPE::GetDerivative(Pressure);

#ifdef HAVE_MPI
    SU2_MPI::Allreduce(&Local_Sens_Mach,  &Total_Sens_Mach,  1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(&Local_Sens_AoA,   &Total_Sens_AoA,   1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(&Local_Sens_Temp,  &Total_Sens_Temp,  1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(&Local_Sens_Press, &Total_Sens_Press, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
    Total_Sens_Mach  = Local_Sens_Mach;
    Total_Sens_AoA   = Local_Sens_AoA;
    Total_Sens_Temp  = Local_Sens_Temp;
    Total_Sens_Press = Local_Sens_Press;
#endif
  }

  /*--- Extract here the adjoint values of everything else that is registered as input in RegisterInput. ---*/

}

void CDiscAdjSolver::SetAdjoint_Output(CGeometry *geometry, CConfig *config){

  bool dual_time = (config->GetUnsteady_Simulation() == DT_STEPPING_1ST ||
      config->GetUnsteady_Simulation() == DT_STEPPING_2ND);

  unsigned short iVar;
  unsigned long iPoint;

  for (iPoint = 0; iPoint < nPoint; iPoint++){
    for (iVar = 0; iVar < nVar; iVar++){
      Solution[iVar] = node[iPoint]->GetSolution(iVar);
    }
    if (dual_time){
      for (iVar = 0; iVar < nVar; iVar++){
        Solution[iVar] += node[iPoint]->GetDual_Time_Derivative(iVar);
      }
    }
    direct_solver->node[iPoint]->SetAdjointSolution(Solution);
  }
}

void CDiscAdjSolver::ResetSensitivity(CGeometry *geometry){
    unsigned short iMarker=0;
    unsigned long iVertex;

    //for (iMarker = 0; iMarker < nMarker; iMarker++){
        for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++){
          LagrangeSens[iMarker][iVertex] = 0.0;
        }
    //}
}

void CDiscAdjSolver::UpdateLagrangeSensitivity(CGeometry *geometry, su2double factor){
    unsigned short iMarker=0;
    unsigned long iVertex;
    /*if(factor==200)     factor=2/((1-rho)*(1-rho));*/
    cout.precision(15);
    std::cout<<"factor: "<<factor<<std::endl;
   // for (iMarker = 0; iMarker < nMarker; iMarker++){
        for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++){
          LagrangeSens[iMarker][iVertex] += factor*CSensitivity[iMarker][iVertex];
          std::cout<<factor*CSensitivity[iMarker][iVertex]<<" ";
        }
   // }
    std::cout<<std::endl;
}

void CDiscAdjSolver::OverwriteSensitivityProjected(CGeometry *geometry){
    unsigned long iVertex;
    for (iVertex = 0; iVertex < geometry->GetnVertex(0); iVertex++){
          geometry->vertex[0][iVertex]->SetAuxVar(LagrangeSens[0][iVertex]);
    }
}

void CDiscAdjSolver::OverwriteGradientProjected(CGeometry *geometry){
    unsigned long iVertex;
    for (iVertex = 0; iVertex < geometry->GetnVertex(0); iVertex++){
          geometry->vertex[0][iVertex]->SetAuxVar(CSensitivityOld[0][iVertex]);
    }
}

void CDiscAdjSolver::SetProjectedSensitivity(unsigned long iDV, su2double value){
    ProjectedSens[iDV]=value;
}

void CDiscAdjSolver::SetProjectedGradient(unsigned long iDV, su2double value){
    ProjectedGradient[iDV]=value;
}

void CDiscAdjSolver::ApplyDesignVar(){
    unsigned long iVertex;
    for (iVertex = 0; iVertex < 38; iVertex++){
        DesignVar[iVertex]+=DesignVarUpdateReal[iVertex];
    }
}

void CDiscAdjSolver::SetDesignVariable(){
    su2double design [38]= {0.00430015707192, 0.00351786344767, 0.0032401667356, 0.00283807455377, 0.00206191727, 0.00108405876297, 0.000555539929847, 0.000786541534876, 0.00183052016302, 0.00338675330365, 0.00428013385902, 0.00345702317815, -0.000799563850502, -0.00492667722335, -0.00460574497224, -0.00466251584455, -0.00460957049057, -0.00465793791568, -0.000929145626517, -0.0031879407669, -0.00490755627308, -0.00484567484714, -0.0047489851551, -0.00449419269773, -0.00481480500091, -0.00162911460905, 0.00124351541834, 0.00201096724969, 0.0011449638727, -0.000714982158573, -0.00220847856125, -0.000765269648138, 0.00416351002763, 0.0048552422848, 0.00493010058661, 0.0049611775592, 0.00494926802939, 0.00196807919392};
    unsigned long iVertex;
    for (iVertex = 0; iVertex < 38; iVertex++){
        DesignVarUpdate[iVertex]=design[iVertex];
        DesignVar[iVertex]=design[iVertex];
    }
}

su2double CDiscAdjSolver::getDVValue(unsigned long iDV){
    return DesignVarUpdate[iDV];
}

void CDiscAdjSolver::SaveSurfaceSensitivity(CGeometry *geometry){
    unsigned short iMarker=0;
    unsigned long iVertex;

 //   for (iMarker = 0; iMarker < nMarker; iMarker++){
        for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++){
          CSensitivityOld[iMarker][iVertex] = CSensitivity[iMarker][iVertex];
        }
 //   }
}

void CDiscAdjSolver::StoreQuadValues(CGeometry *geometry, CConfig *config){
    unsigned short iMarker;
    unsigned long iVertex;

    for (iMarker = 0; iMarker < nMarker; iMarker++){
        for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++){
          CSensitivityQuad[countQuadrature][iMarker][iVertex] = CSensitivityOld[iMarker][iVertex];
          LagrangeSensQuad[countQuadrature][iMarker][iVertex] = LagrangeSens[iMarker][iVertex];
        }
    }
    for (unsigned short iValue=0; iValue<config->GetConstraintNum();iValue++){
        Constraint_Save_Quad[countQuadrature][iValue]=Constraint_Save[iValue];
    }
    Obj_Save_Quad[countQuadrature]=Obj_Save;
    Lagrangian_Value_Quad[countQuadrature]=Lagrangian_Value;
    countQuadrature=countQuadrature+1;

}

void CDiscAdjSolver::ResetStoredValues(CGeometry *geometry, CConfig *config){
    unsigned short iMarker;
    unsigned long iVertex;
    countQuadrature=0;
    for (unsigned short num=0; num<100;num++){
        for (iMarker = 0; iMarker < nMarker; iMarker++){
            for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++){
              CSensitivityQuad[num][iMarker][iVertex] = 0;
              LagrangeSensQuad[num][iMarker][iVertex] = 0;
            }
        }
        for (unsigned short iValue=0; iValue<config->GetConstraintNum();iValue++){
            Constraint_Save_Quad[num][iValue]=0;
        }
        Obj_Save_Quad[num]=0;
        Lagrangian_Value_Quad[num]=0;
    }
}

void CDiscAdjSolver::SetObjSave(CGeometry *geometry, CConfig *config, su2double obj_val){
    Obj_Save=obj_val;
}

void CDiscAdjSolver::SetLagrangian(CGeometry *geometry, CConfig *config, su2double lagrangian){
    Lagrangian_Value=lagrangian;
}

void CDiscAdjSolver::SetCSens(CGeometry *geometry, CConfig *config, su2double csens, unsigned short num){
    CSensitivity[0][num]=csens;
}

void CDiscAdjSolver::SetLagrangeSens(CGeometry *geometry, CConfig *config, su2double lsens, unsigned short num){
    LagrangeSens[0][num]=lsens;
}

void CDiscAdjSolver::SetConstraintSave(CGeometry *geometry, CConfig *config, su2double cons, unsigned short num){
    Constraint_Save[num]=cons;
}

su2double CDiscAdjSolver::GetStoredObjective(CGeometry *geometry, CConfig *config, unsigned short num){
    return Obj_Save_Quad[num];
}

su2double CDiscAdjSolver::GetStoredConstraints(CGeometry *geometry, CConfig *config, unsigned short num, unsigned short indx){
    return Constraint_Save_Quad[num][indx];
}

su2double CDiscAdjSolver::GetStoredLagragian(CGeometry *geometry, CConfig *config, unsigned short num){
    return Lagrangian_Value_Quad[num];
}

su2double CDiscAdjSolver::GetStoredCSens(CGeometry *geometry, CConfig *config, unsigned short num, unsigned short indx){
    return CSensitivityQuad[num][0][indx];
}

su2double CDiscAdjSolver::GetStoredLagrangeSens(CGeometry *geometry, CConfig *config, unsigned short num, unsigned short indx){
    return LagrangeSensQuad[num][0][indx];
}



void CDiscAdjSolver::ResetExpValues(CGeometry *geometry, CConfig *config){
    unsigned long iVertex;
        for (iVertex = 0; iVertex < geometry->GetnVertex(0); iVertex++){
          ExpCSensitivityOld[0][iVertex] = 0.0;
          ExpLagrangeSens[0][iVertex]=0.0;
        }
        for (unsigned short iValue=0; iValue<config->GetConstraintNum();iValue++){
            ExpConstraint_Save[iValue]=0.0;
        }
        ExpObjFunc_Value=0.0;
        ExpLagrangian_Value=0.0;
}

void CDiscAdjSolver::SumExpValues(CGeometry *geometry, CConfig *config,unsigned short numQuad){
    unsigned long iVertex;
        for (iVertex = 0; iVertex < geometry->GetnVertex(0); iVertex++){
          ExpCSensitivityOld[0][iVertex] += weights[numQuad]*CSensitivityOld[0][iVertex];
          ExpLagrangeSens[0][iVertex] += weights[numQuad]*LagrangeSens[0][iVertex];
        }
        for (unsigned short iValue=0; iValue<config->GetConstraintNum();iValue++){
            ExpConstraint_Save[iValue]+=weights[numQuad]*Constraint_Save[iValue];
        }
        ExpObjFunc_Value+=weights[numQuad]*Obj_Save;
        ExpLagrangian_Value+=weights[numQuad]*Lagrangian_Value;
}

su2double CDiscAdjSolver::GetMachP(unsigned short numQuad){
    return machp[numQuad];
}

void CDiscAdjSolver::DistributeExpValues(CGeometry *geometry, CConfig *config){
    unsigned long iVertex;
    su2double pi = 3.141592653589793;
        for (iVertex = 0; iVertex < geometry->GetnVertex(0); iVertex++){
          CSensitivityOld[0][iVertex] = 1./sqrt(pi)*ExpCSensitivityOld[0][iVertex];
          LagrangeSens[0][iVertex]=1./sqrt(pi)*ExpLagrangeSens[0][iVertex];
        }
        Obj_Save=1./sqrt(pi)*ExpObjFunc_Value;
        ObjFunc_Value=1./sqrt(pi)*ExpObjFunc_Value;
        for (unsigned short iValue=0; iValue<config->GetConstraintNum();iValue++){
            Constraint_Save[iValue]=1./sqrt(pi)*ExpConstraint_Save[iValue];
        }
        Lagrangian_Value=1./sqrt(pi)*ExpLagrangian_Value;
}

su2double CDiscAdjSolver::SensitivityNorm(CGeometry *geometry){
    unsigned long iVertex;
    su2double norm=0;
        for (iVertex = 0; iVertex < geometry->GetnVertex(0); iVertex++){
             norm+=CSensitivity[0][iVertex]*CSensitivity[0][iVertex];
        }
    norm=sqrt(norm);
    return norm;
}

void CDiscAdjSolver::SetSensitivity(CGeometry *geometry, CConfig *config){

  unsigned long iPoint;
  unsigned short iDim;
  su2double *Coord, Sensitivity, eps;

  bool time_stepping = (config->GetUnsteady_Simulation() != STEADY);

  for (iPoint = 0; iPoint < nPoint; iPoint++){
    Coord = geometry->node[iPoint]->GetCoord();

    for (iDim = 0; iDim < nDim; iDim++){

      Sensitivity = SU2_TYPE::GetDerivative(Coord[iDim]);

      /*--- Set the index manually to zero. ---*/

      AD::ResetInput(Coord[iDim]);

      /*--- If sharp edge, set the sensitivity to 0 on that region ---*/

      if (config->GetSens_Remove_Sharp()) {
        eps = config->GetLimiterCoeff()*config->GetRefElemLength();
        if ( geometry->node[iPoint]->GetSharpEdge_Distance() < config->GetSharpEdgesCoeff()*eps )
          Sensitivity = 0.0;
      }

      node[iPoint]->SetSensitivity(iDim, Sensitivity);
    }
  }
  SetSurface_Sensitivity(geometry, config);
}

void CDiscAdjSolver::SetMixedSensitivity(CGeometry *geometry, CConfig *config){

  unsigned long iPoint, iVar;
  unsigned short iDim;
  su2double *Coord, Sensitivity, eps;

  //std::cout<<"FORWARDDER"<<SU2_TYPE::GetMixedDerivative(ObjFunc_Value)<<std::endl;
  AD::ResetVectorPosition();

  for (iPoint = 0; iPoint < nPoint; iPoint++){
    Coord = geometry->node[iPoint]->GetCoord();

    for (iDim = 0; iDim < nDim; iDim++){
     // Sensitivity = Coord[iDim].getGradient().getGradient();
     // std::cout<<Sensitivity<<" ";
      Sensitivity = SU2_TYPE::GetMixedDerivative(Coord[iDim]);


      /*--- Set the index manually to zero. ---*/

     AD::ResetInput(Coord[iDim]);

      /*--- If sharp edge, set the sensitivity to 0 on that region ---*/

      if (config->GetSens_Remove_Sharp()) {
        eps = config->GetLimiterCoeff()*config->GetRefElemLength();
        if ( geometry->node[iPoint]->GetSharpEdge_Distance() < config->GetSharpEdgesCoeff()*eps )
          Sensitivity = 0.0;
      }
     // if (!time_stepping){
        node[iPoint]->SetSensitivity(iDim, Sensitivity);
     // } else {
     //   node[iPoint]->SetSensitivity(iDim, node[iPoint]->GetSensitivity(iDim) + Sensitivity);
     // }
    }
  }
/*  std::cout<<std::endl;
  std::cout<<"direct_solver->node"<<std::endl;
  for (iPoint = 0; iPoint < nPoint; iPoint++){

        for (iVar = 0; iVar < nVar; iVar++){
          std::cout<<SU2_TYPE::GetMixedDerivative(direct_solver->node[iPoint]->GetSolution(iVar))<<" ";
        }
  }*/
  SetSurface_Sensitivity(geometry, config);
}

void CDiscAdjSolver::SetSensDensity(CGeometry *geometry, CConfig *config){

    unsigned long iElem;
    su2double Density, Sensitivity;
    for (iElem = 0; iElem < geometry->GetnElem(); iElem++) {
      Density = geometry->elem[iElem]->GetDensity()[0];
      Sensitivity= SU2_TYPE::GetDerivative(Density);
      AD::ResetInput(Density);
      CSensitivity[0][iElem]=Sensitivity;
    }
}

void CDiscAdjSolver::WriteDesignVariable(){
    ofstream myfile;
    unsigned long iVertex;

    myfile.open ("designvar.txt");
    for (iVertex=0;iVertex<38;iVertex++)
    {
        myfile <<std::setprecision(16)<<DesignVarUpdate[iVertex]<<std::endl;

    }
    myfile.close();
}

void CDiscAdjSolver::DesignUpdate(CGeometry *geometry, CConfig *config){
    unsigned long iVertex,jVertex;
    for (iVertex = 0; iVertex < geometry->GetnVertex(0); iVertex++){
        UpdateSens[iVertex]=0.0;
        for (jVertex = 0; jVertex < geometry->GetnVertex(0); jVertex++){
            UpdateSens[iVertex]-=Hess[iVertex][jVertex]*CSensitivityOld[0][jVertex];
        }

    }
}

void CDiscAdjSolver::DesignStep(su2double values){
    unsigned long iVertex;
    for (iVertex = 0; iVertex < 38; iVertex++){
        DesignVarUpdate[iVertex]=values;
    }
}

void CDiscAdjSolver::DesignMinus(){
    unsigned long iVertex;
    for (iVertex = 0; iVertex < 38; iVertex++){
        DesignVarUpdate[iVertex]=DesignVarUpdateReal[iVertex]-DesignVarUpdateSave[iVertex];
    }
}

void CDiscAdjSolver::CalculatePhi(su2double steplen, su2double& Phi, su2double& dPhi){
    su2double helper=0.0;
/*    if(steplen==0.0){
        unsigned long iVertex;
        for (iVertex = 0; iVertex < 38; iVertex++){
            helper+=UpdateSens[iVertex]*ProjectedSensOld[iVertex];
        }
        dPhi=helper;
        Phi=Lagrangian_Value_Old;
    }else{*/
        unsigned long iVertex;
        for (iVertex = 0; iVertex < 38; iVertex++){
            helper+=UpdateSens[iVertex]*ProjectedSens[iVertex];
        }
        dPhi=helper;
        Phi=Lagrangian_Value;
//    }
}

su2double CDiscAdjSolver::QuadraticApproximation(su2double steplen){
    su2double helper=0.0;
    unsigned long iVertex;
    for (iVertex = 0; iVertex < 38; iVertex++){
        helper+=UpdateSens[iVertex]*ProjectedSens[iVertex];
    }
    GradPhiCubic=helper;
    PhiCubic=Lagrangian_Value_Old;
    PhiOld=Lagrangian_Value;
    StepOld=steplen;
    return ((-helper*steplen*steplen)/(2*(Lagrangian_Value-Lagrangian_Value_Old-helper*steplen)));
}



bool CDiscAdjSolver::CheckDescentDirection(su2double steplen){
    su2double helper=0.0; //Phi'(0)
    unsigned long iVertex;
    for (iVertex = 0; iVertex < 38; iVertex++){
        helper+=UpdateSens[iVertex]*ProjectedSens[iVertex];
    }
    if(helper<=0) return true;
    else{
        /*for (iVertex = 0; iVertex < 38; iVertex++){
            UpdateSens[iVertex]=-UpdateSens[iVertex];
        }*/
        return false;
    }
}

void CDiscAdjSolver::ChangeDirection(){
    unsigned long iVertex;
        for (iVertex = 0; iVertex < 38; iVertex++){
            UpdateSens[iVertex]=-UpdateSens[iVertex];
        }
}

su2double CDiscAdjSolver::CubicApproximation(su2double steplen){
    su2double avalue=StepOld*StepOld*(Lagrangian_Value-PhiCubic-GradPhiCubic*steplen)-steplen*steplen*(PhiOld-PhiCubic-GradPhiCubic*StepOld);
    su2double bvalue=-StepOld*StepOld*StepOld*(Lagrangian_Value-PhiCubic-GradPhiCubic*steplen)+steplen*steplen*steplen*(PhiOld-PhiCubic-GradPhiCubic*StepOld);
    avalue=1.0/(StepOld*StepOld*steplen*steplen*(steplen-StepOld))*avalue;
    bvalue=1.0/(StepOld*StepOld*steplen*steplen*(steplen-StepOld))*bvalue;
    su2double steplennew;
    steplennew=(-bvalue+sqrt(bvalue*bvalue-3*avalue*GradPhiCubic))/(3*avalue);
    PhiOld=Lagrangian_Value;
    StepOld=steplen;
    if(fabs(steplennew-steplen)<1E-15){
        std::cout<<"small correction: "<<steplennew<<std::endl;
        steplen=steplen*0.5;
    }else if(fabs(steplennew-steplen)>0.5 || steplennew==0){
        std::cout<<"big correction: "<<steplennew<<std::endl;
        steplen=steplen*0.5;
    }else{
        std::cout<<"CUBIC: "<<steplennew<<std::endl;
         steplen = steplennew;
    }
    return steplen;
}

void CDiscAdjSolver::DesignUpdateProjected(CGeometry *geometry, CConfig *config, unsigned short ExtIter, su2double steplen){
    unsigned long iVertex;
    su2double normsens=0;
    for (iVertex = 0; iVertex < 38; iVertex++){
        DesignVarUpdateSave[iVertex]=DesignVarUpdateReal[iVertex];
        DesignVarUpdate[iVertex]=0.0;
        normsens+=UpdateSens[iVertex]*UpdateSens[iVertex];
    }
    normsens=sqrt(normsens/(38*38));
    std::cout<<"Norm of Update: "<<normsens<<std::endl;

    for (iVertex = 0; iVertex < 38; iVertex++){
        DesignVarUpdate[iVertex]=UpdateSens[iVertex]*steplen;
        if(config->GetOSStepAdaptive()){
            if (normsens>10000){ //vorher normsens>10000000
                        DesignVarUpdate[iVertex]=DesignVarUpdate[iVertex]*1E-4; //was >100, 1E-2 //1E-3
            }
        /*    else if(normsens>1&&config->GetExtIter()>1416){ //vorher normsens>100000, 1E-1
                        DesignVarUpdate[iVertex]=DesignVarUpdate[iVertex]*1E-3; //was >1, 1E-1
            }*/
     /*       else if(normsens<=100000){
                        DesignVarUpdate[iVertex]=DesignVarUpdate[iVertex]*1.0; // was <1
            }*/
        }
        //set box constraints
        if((DesignVar[iVertex]+DesignVarUpdate[iVertex])>0.005)  DesignVarUpdate[iVertex]=0.005-DesignVar[iVertex];
        if((DesignVar[iVertex]+DesignVarUpdate[iVertex])<-0.005) DesignVarUpdate[iVertex]=-0.005-DesignVar[iVertex];

        DesignVarUpdateReal[iVertex]=DesignVarUpdate[iVertex];
    }
    std::cout<<std::endl;
}

su2double CDiscAdjSolver::DesignUpdateBounds(CGeometry *geometry, CConfig *config, unsigned short ExtIter, su2double steplen){
    unsigned long iVertex;
    su2double normsens=0;
    for (iVertex = 0; iVertex < 38; iVertex++){
        DesignVarUpdateSave[iVertex]=DesignVarUpdateReal[iVertex];
        DesignVarUpdate[iVertex]=0.0;
        normsens+=UpdateSens[iVertex]*UpdateSens[iVertex];
    }
    normsens=sqrt(normsens/(38*38));
    std::cout<<"Norm of Update: "<<normsens<<std::endl;

    for (iVertex = 0; iVertex < 38; iVertex++){
        while((DesignVar[iVertex]+UpdateSens[iVertex]*steplen)>0.005||(DesignVar[iVertex]+UpdateSens[iVertex]*steplen)<-0.005){
            steplen=steplen*0.5;
        }
    }
    for (iVertex = 0; iVertex < 38; iVertex++){
        DesignVarUpdate[iVertex]=UpdateSens[iVertex]*steplen;
        DesignVarUpdateReal[iVertex]=DesignVarUpdate[iVertex];
    }
    if(steplen<1E-30) std::cout<<"REACHED DESIGN VARIABLE BOUNDS"<<std::endl;
    return steplen;
}

bool CDiscAdjSolver::CheckFirstWolfe(su2double steplen){
    su2double helper=0.0;
    unsigned long iVertex;
    std::cout<<"LagrangeOld: "<<Lagrangian_Value_Old<<", LagrangeNew: "<<Lagrangian_Value<<", Stepsize: "<<steplen<<std::endl;
    for (iVertex = 0; iVertex < 38; iVertex++){
        helper+=DesignVarUpdateReal[iVertex]*ProjectedSens[iVertex];//UpdateSens[iVertex]*ProjectedSens[iVertex];
    }
    if (Lagrangian_Value<=(Lagrangian_Value_Old+1E-4*helper)){ //*steplen*helper)){
        return false;
    }
    else{
        std::cout<<"First Wolfe Condition not satisfied!"<<std::endl;
         return true;
       // return false;
    }
}

void CDiscAdjSolver::BFGSUpdateProjected(CGeometry *geometry, CConfig *config, unsigned short ExtIter){

    unsigned long iVertex, jVertex, kVertex, lVertex, mVertex;
    su2double *rk,*duk,*wone;
    rk=new su2double[38];
    duk=new su2double[38];
    wone=new su2double[38];
    su2double vk=0;
    su2double normrk=0;
    su2double normduk=0;

    //Output of Gradients and Information

    std::cout<<"Projected Gradient of Augmented Lagrangian "<<std::endl;
    for (iVertex=0;iVertex<38;iVertex++)
    {
        std::cout<<ProjectedSens[iVertex]<<" ";
    }
    std::cout<<std::endl;
    std::cout<<"iterationcount: "<<TotalIterations<<std::endl;
    std::cout<<"objfuncvalue: "<<Obj_Save<<std::endl;
    std::cout<<"constraintvalue: "<<Constraint_Save[0]<<std::endl;
    if(config->GetConstraintNum()>1) std::cout<<"momvalue: "<<Constraint_Save[1]<<std::endl;
    if(config->GetConstraintNum()>2) std::cout<<"thickness: "<<Constraint_Save[2]<<std::endl;
    std::cout<<"Projected Gradient N_u"<<std::endl;
    for (iVertex=0;iVertex<38;iVertex++)
    {
        std::cout<<ProjectedGradient[iVertex]<<" ";
    }
    std::cout<<std::endl;


    if(ExtIter>config->GetOneShotStart()){
        for (iVertex = 0; iVertex < 38; iVertex++){
            //DesignVar[iVertex]+=DesignVarUpdate[iVertex];
            rk[iVertex]=ProjectedSens[iVertex]-ProjectedSensOld[iVertex];
            duk[iVertex]=DesignVarUpdate[iVertex];
            vk+=rk[iVertex]*duk[iVertex];
            normrk+=rk[iVertex]*rk[iVertex];
            normduk+=duk[iVertex]*duk[iVertex];
        }
        std::cout<<std::endl;
        std::cout<<"vk "<<vk<<std::endl;
        std::cout<<"normduk "<<normduk<<", normrk "<<normrk<<", vk/normduk "<<vk/normduk<<std::endl;
        if(config->GetDampedBFGS()){
            //normalization
            for (iVertex=0;iVertex<38;iVertex++){
                rk[iVertex]=rk[iVertex]/sqrt(normrk);
                duk[iVertex]=duk[iVertex]/sqrt(normduk);
                if(config->NormalizeHB()){
                for (jVertex=0;jVertex<38;jVertex++){
                    Bess[iVertex][jVertex]=Bess[iVertex][jVertex]*sqrt(normduk)/sqrt(normrk);
                    Hess[iVertex][jVertex]=Hess[iVertex][jVertex]*sqrt(normrk)/sqrt(normduk);
                }
                }
            }
            //we store the Hessian approximation (Bess) and the inverse (Hess)
            su2double sBs=0.0;
            vk=0;
            su2double bmin=config->GetDampedMin();
            su2double bmax=config->GetDampedMax();
            su2double theta=1.0;
            su2double** MatA;
            MatA=new su2double*[38];
            for (iVertex=0;iVertex<38;iVertex++){
                MatA[iVertex]=new su2double[38];
                for (jVertex=0;jVertex<38;jVertex++){
                    MatA[iVertex][jVertex]=0.0;
                }
            }
            for (iVertex=0;iVertex<38;iVertex++){
                vk+=rk[iVertex]*duk[iVertex];
                for (jVertex=0;jVertex<38;jVertex++){
                    sBs+=duk[iVertex]*Bess[iVertex][jVertex]*duk[jVertex];
                }
            }
            if(config->GetDampedBFGSPow()){
               su2double gamma=config->GetDampedGamma();
               if(vk<gamma*sBs){
                   theta=0.9*((1-gamma)*sBs)/(sBs-vk);
               }
               std::cout<<"correction: "<<theta<<std::endl;
            }else{
            if(((vk/normduk)<bmin)||((vk/normduk)>bmax)){
                su2double useb;
                if((vk/normduk)<bmin){
                    useb=bmin;
                }else{
                    useb=bmax;
                }
                for (iVertex=0;iVertex<38;iVertex++){
                    for (jVertex=0;jVertex<38;jVertex++){
                        if(iVertex==jVertex){
                            theta+=duk[iVertex]*(Bess[iVertex][jVertex]-useb)*duk[jVertex];
                        }else{
                            theta+=duk[iVertex]*(Bess[iVertex][jVertex])*duk[jVertex];
                        }
                    }
                }
                theta=0.9*(theta/(sBs-vk));
                std::cout<<"correction: "<<theta<<std::endl;
            }
            }
            su2double normnewrk=0;
            for (iVertex=0;iVertex<38;iVertex++){
                rk[iVertex]=theta*rk[iVertex];
                for (jVertex=0;jVertex<38;jVertex++){
                    rk[iVertex]+=(1-theta)*Bess[iVertex][jVertex]*duk[jVertex];
                }
                normnewrk+=rk[iVertex]*rk[iVertex];
            }
            normnewrk=sqrt(normnewrk);
            if(config->NormalizeNewY()){
            for (iVertex=0;iVertex<38;iVertex++){
                rk[iVertex]=rk[iVertex]/normnewrk;
            }
            }
            vk=0;
            for (iVertex=0;iVertex<38;iVertex++){
                vk+=rk[iVertex]*duk[iVertex];
            }
            std::cout<<"vknew "<<vk<<std::endl;
            //normal Update
            //for Bess
            su2double* Bs=new su2double[38];
            for (iVertex=0;iVertex<38;iVertex++){
                Bs[iVertex]=0;
                for (jVertex=0;jVertex<38;jVertex++){
                    Bs[iVertex]+=Bess[iVertex][jVertex]*duk[jVertex];
                }
            }
            for (iVertex=0;iVertex<38;iVertex++){
                for (jVertex=0;jVertex<38;jVertex++){
                    MatA[iVertex][jVertex]=Bess[iVertex][jVertex]+(1.0/vk)*rk[iVertex]*rk[jVertex]-(1.0/sBs)*Bs[iVertex]*Bs[jVertex];
                        /*for (lVertex=0; lVertex<38; lVertex++){
                            for (mVertex=0; mVertex<38; mVertex++){
                                MatA[iVertex][jVertex]-=(1.0/sBs)*Bess[iVertex][lVertex]*duk[lVertex]*duk[mVertex]*Bess[mVertex][jVertex];
                            }
                        }*/
                }
            }
            delete [] Bs;
            for (iVertex=0;iVertex<38;iVertex++){
                for (jVertex=0;jVertex<38;jVertex++){
                    Bess[iVertex][jVertex]=MatA[iVertex][jVertex];
                }
            }
            //for Hess
            for (iVertex=0;iVertex<38;iVertex++){
                for (jVertex=0;jVertex<38;jVertex++){
                    MatA[iVertex][jVertex]=Hess[iVertex][jVertex]+(1.0/vk)*duk[iVertex]*duk[jVertex];
                    for (kVertex=0; kVertex<38; kVertex++){
                        MatA[iVertex][jVertex]+=-(1.0/vk)*duk[iVertex]*Hess[jVertex][kVertex]*rk[kVertex]-(1.0/vk)*duk[jVertex]*Hess[iVertex][kVertex]*rk[kVertex];
                        for (lVertex=0; lVertex<38; lVertex++){
                            MatA[iVertex][jVertex]+=(1.0/vk)*(1.0/vk)*duk[iVertex]*duk[jVertex]*rk[lVertex]*Hess[lVertex][kVertex]*rk[kVertex];
                        }
                    }
                }
            }
            for (iVertex=0;iVertex<38;iVertex++){
                for (jVertex=0;jVertex<38;jVertex++){
                    Hess[iVertex][jVertex]=MatA[iVertex][jVertex];
                }
            }
            for (iVertex=0;iVertex<38;iVertex++){
                rk[iVertex]=rk[iVertex]*sqrt(normrk);
                duk[iVertex]=duk[iVertex]*sqrt(normduk);
                if(config->NormalizeHB()){
                for (jVertex=0;jVertex<38;jVertex++){
                    Bess[iVertex][jVertex]=Bess[iVertex][jVertex]*sqrt(normrk)/sqrt(normduk);
                    Hess[iVertex][jVertex]=Hess[iVertex][jVertex]*sqrt(normduk)/sqrt(normrk);
                }
                }
            }
        }else{

        if (vk>0 && ((fabs(vk)>1E-3) || (config->GetCheckVk()==false))){
            //HESSOLD
     /*   for (iVertex = 0; iVertex < 38; iVertex++){
          wone[iVertex]=0.0;
          for (jVertex = 0; jVertex < 38; jVertex++){
            wone[iVertex]+=Hess[iVertex][jVertex]*rk[jVertex];
          }
        }
        for (iVertex = 0; iVertex < 38; iVertex++){
          wtwo+=rk[iVertex]*wone[iVertex];
        }
        for (iVertex = 0; iVertex < 38; iVertex++){
          for (jVertex = 0; jVertex < 38; jVertex++){
            Hess[iVertex][jVertex]=Hess[iVertex][jVertex]-(1.0/vk)*(wone[iVertex]*duk[jVertex]+wone[jVertex]*duk[iVertex])+(1.0/vk)*(1+wtwo/vk)*duk[iVertex]*duk[jVertex];
          }
        }*/
        su2double** MatA;//, **MatB, **MatC;
        MatA=new su2double*[38];
     //   MatB=new su2double*[38];
      //  MatC=new su2double*[38];
        for (iVertex=0;iVertex<38;iVertex++){
            MatA[iVertex]=new su2double[38];
        //    MatB[iVertex]=new su2double[38];
        //    MatC[iVertex]=new su2double[38];
            for (jVertex=0;jVertex<38;jVertex++){
                MatA[iVertex][jVertex]=0.0;
       //         MatB[iVertex][jVertex]=0.0;
       //         MatC[iVertex][jVertex]=0.0;
            }
        }
      /*  for (iVertex=0;iVertex<38;iVertex++){
            for (jVertex=0;jVertex<38;jVertex++){
                MatA[iVertex][jVertex]=duk[iVertex]*rk[jVertex];
            }
        }
        for (iVertex=0;iVertex<38;iVertex++){
            for (jVertex=0;jVertex<38;jVertex++){
                for (kVertex=0; kVertex<38; kVertex++){
                    if(iVertex==kVertex) MatB[iVertex][jVertex]+=(1.0-(1.0/vk)*MatA[iVertex][kVertex])*Hess[kVertex][jVertex];
                    else MatB[iVertex][jVertex]+=(-(1.0/vk)*MatA[iVertex][kVertex])*Hess[kVertex][jVertex];
                }
            }
        }
        for (iVertex=0;iVertex<38;iVertex++){
            for (jVertex=0;jVertex<38;jVertex++){
                for (kVertex=0; kVertex<38; kVertex++){
                    if(jVertex==kVertex) MatC[iVertex][jVertex]+=MatB[iVertex][kVertex]*(1.0-(1.0/vk)*MatA[jVertex][kVertex]);
                    else  MatC[iVertex][jVertex]+=MatB[iVertex][kVertex]*(-(1.0/vk)*MatA[jVertex][kVertex]);
                }
            }
        }
        for (iVertex=0;iVertex<38;iVertex++){
            for (jVertex=0;jVertex<38;jVertex++){
                Hess[iVertex][jVertex]=MatC[iVertex][jVertex]+(1.0/vk)*duk[iVertex]*duk[jVertex];
            }
        }*/
        //HESS2
        if(config->GetLBFGS()){
            unsigned long maxCount=config->GetLBFGSNum()-1;
            if(maxCount>BFGSCount) maxCount=BFGSCount;
            //std::cout<<maxCount<<std::endl;
            for(unsigned long kCount=0;kCount<maxCount;kCount++){
                for (iVertex=0;iVertex<38;iVertex++){
                    rkStore[kCount+1][iVertex]=rkStore[kCount][iVertex];
                    dukStore[kCount+1][iVertex]=dukStore[kCount][iVertex];
                }
            }
            if(config->GetHInit()) std::cout<<"Initialize H with "<<vk/normrk<<std::endl;
            for (iVertex=0;iVertex<38;iVertex++){
                rkStore[maxCount][iVertex]=rk[iVertex];
                dukStore[maxCount][iVertex]=duk[iVertex];
                for (jVertex=0;jVertex<38;jVertex++){
                    Hess[iVertex][jVertex]=0;
                    if(iVertex==jVertex) Hess[iVertex][jVertex]=1.0;
                    if(config->GetHInit()&& (iVertex==jVertex)){
                        Hess[iVertex][jVertex]=vk/normrk;
                    }
                }
            }
            for (unsigned long kCount=0;kCount<=maxCount;kCount++){
                vk=0;
                for (iVertex=0;iVertex<38;iVertex++){
                    vk+=rkStore[kCount][iVertex]*dukStore[kCount][iVertex];
                }
                for (iVertex=0;iVertex<38;iVertex++){
                    for (jVertex=0;jVertex<38;jVertex++){
                        MatA[iVertex][jVertex]=Hess[iVertex][jVertex]+(1.0/vk)*dukStore[kCount][iVertex]*dukStore[kCount][jVertex];
                        for (kVertex=0; kVertex<38; kVertex++){
                            MatA[iVertex][jVertex]+=-(1.0/vk)*dukStore[kCount][iVertex]*Hess[jVertex][kVertex]*rkStore[kCount][kVertex]-(1.0/vk)*dukStore[kCount][jVertex]*Hess[iVertex][kVertex]*rkStore[kCount][kVertex];
                            for (lVertex=0; lVertex<38; lVertex++){
                                MatA[iVertex][jVertex]+=(1.0/vk)*(1.0/vk)*dukStore[kCount][iVertex]*dukStore[kCount][jVertex]*rkStore[kCount][lVertex]*Hess[lVertex][kVertex]*rkStore[kCount][kVertex];
                            }
                        }
                    }
                }
                for (iVertex=0;iVertex<38;iVertex++){
                    for (jVertex=0;jVertex<38;jVertex++){
                        Hess[iVertex][jVertex]=MatA[iVertex][jVertex];
                    }
                }
            }
        }else{
            for (iVertex=0;iVertex<38;iVertex++){
                for (jVertex=0;jVertex<38;jVertex++){
                    MatA[iVertex][jVertex]=Hess[iVertex][jVertex]+(1.0/vk)*duk[iVertex]*duk[jVertex];
                    for (kVertex=0; kVertex<38; kVertex++){
                        MatA[iVertex][jVertex]+=-(1.0/vk)*duk[iVertex]*Hess[jVertex][kVertex]*rk[kVertex]-(1.0/vk)*duk[jVertex]*Hess[iVertex][kVertex]*rk[kVertex];
                        for (lVertex=0; lVertex<38; lVertex++){
                            MatA[iVertex][jVertex]+=(1.0/vk)*(1.0/vk)*duk[iVertex]*duk[jVertex]*rk[lVertex]*Hess[lVertex][kVertex]*rk[kVertex];
                        }
                    }
                }
            }
            for (iVertex=0;iVertex<38;iVertex++){
                for (jVertex=0;jVertex<38;jVertex++){
                    Hess[iVertex][jVertex]=MatA[iVertex][jVertex];
                }
            }
        }

        for (iVertex=0;iVertex<38;iVertex++){
            delete [] MatA[iVertex];
        //     delete [] MatB[iVertex];
        //     delete [] MatC[iVertex];
        }
        delete [] MatA;
       // delete [] MatB;
       // delete [] MatC;
        BFGSCount=BFGSCount+1;
        }else{
            std::cout<<"!!!!!!!!!!!!!!!!ATTENTION-HESSIAN NON-POSITIVE-DEFINITE!!!!!!!!!!!!!!!!!!!"<<std::endl;
            if(config->GetIdentityHessian()){
            for (iVertex = 0; iVertex < 38; iVertex++){
              for (jVertex = 0; jVertex < 38; jVertex++){
                Hess[iVertex][jVertex]=0.0;
                if(iVertex==jVertex) Hess[iVertex][jVertex]=1.0;
                if(config->GetHInit()&& (iVertex==jVertex)) Hess[iVertex][jVertex]=config->GetHScale();
              }
            }
            }
        }
        }
    }

    std::cout<<"Design Variable "<<std::endl;
    for (iVertex=0;iVertex<38;iVertex++)
    {
        std::cout<<DesignVar[iVertex]<<" ";
    }
    std::cout<<std::endl;

        for (iVertex = 0; iVertex < 38; iVertex++){
          ProjectedSensOld[iVertex] = ProjectedSens[iVertex];
       //   DesignVarOld[iVertex]=DesignVar[iVertex];
        }

        Lagrangian_Value_Old=Lagrangian_Value;

        for (iVertex = 0; iVertex < 38; iVertex++){

            UpdateSens[iVertex]=0.0;
            for (jVertex = 0; jVertex < 38; jVertex++){
                UpdateSens[iVertex]-=Hess[iVertex][jVertex]*ProjectedGradient[jVertex];

            }
        }
     delete [] rk;
     delete [] duk;
     delete [] wone;

}

void CDiscAdjSolver::SetSensitivityFD(CGeometry *geometry, CConfig *config){

  unsigned long iPoint;
  unsigned short iDim;
  su2double *Coord, Sensitivity;
  unsigned short iMarker=0;
  unsigned long iVertex;
  su2double *Normal, Area, Prod, Sens = 0.0, SensDim;
  su2double Total_Sens_Geo_local = 0.0;
  Total_Sens_Geo = 0.0;
  su2double stepsize=config->GetFDStep();
  su2double eps;

  for (iPoint = 0; iPoint < nPoint; iPoint++){
    Coord = geometry->node[iPoint]->GetCoord();

    for (iDim = 0; iDim < nDim; iDim++){
      Sensitivity = SU2_TYPE::GetDerivative(Coord[iDim]);

      AD::ResetInput(Coord[iDim]);

      if (config->GetSens_Remove_Sharp()) {
        eps = config->GetLimiterCoeff()*config->GetRefElemLength();
        if ( geometry->node[iPoint]->GetSharpEdge_Distance() < config->GetSharpEdgesCoeff()*eps )
          Sensitivity = 0.0;
      }

      node[iPoint]->SetSensitivity(iDim, Sensitivity);
    }
  }

  std::cout<<"FDStep: "<<stepsize<<std::endl;
//  for (iMarker = 0; iMarker < nMarker; iMarker++){
    Sens_Geo[iMarker] = 0.0;
    /* --- Loop over boundary markers to select those for Euler walls and NS walls --- */

    if(config->GetMarker_All_KindBC(iMarker) == EULER_WALL
       || config->GetMarker_All_KindBC(iMarker) == HEAT_FLUX
       || config->GetMarker_All_KindBC(iMarker) == ISOTHERMAL){
      for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++){
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
        Prod = 0.0;
        Area = 0.0;
        for (iDim = 0; iDim < nDim; iDim++){
          /* --- retrieve the gradient calculated with AD -- */
          SensDim = node[iPoint]->GetSensitivity(iDim);

          /* --- calculate scalar product for projection onto the normal vector ---*/
          Prod += Normal[iDim]*SensDim;
          Area += Normal[iDim]*Normal[iDim];
        }
        Area = sqrt(Area);

        /* --- projection of the gradient
         *     calculated with AD onto the normal
         *     vector of the surface --- */
        //Sens = Prod;
        Sens=Prod/Area;

        /*--- Compute sensitivity for each surface point ---*/
        CSensitivity[iMarker][iVertex] = -Sens;
        if (geometry->node[iPoint]->GetFlip_Orientation())
          CSensitivity[iMarker][iVertex] = -CSensitivity[iMarker][iVertex];
        CSensitivity[iMarker][iVertex] = (CSensitivity[iMarker][iVertex]-CSensitivityOld[iMarker][iVertex])/stepsize;

        if (geometry->node[iPoint]->GetDomain()){
          Sens_Geo[iMarker] += Sens*Sens;
        }
      }
      Total_Sens_Geo_local += sqrt(Sens_Geo[iMarker]);

    }
 // }

#ifdef HAVE_MPI
  SU2_MPI::Allreduce(&Total_Sens_Geo_local,&Total_Sens_Geo,1,MPI_DOUBLE,MPI_SUM, MPI_COMM_WORLD);
#else
  Total_Sens_Geo = Total_Sens_Geo_local;
#endif
}

void CDiscAdjSolver::SetSurface_Sensitivity(CGeometry *geometry, CConfig *config){
  unsigned short iMarker=0,iDim;
  unsigned long iVertex, iPoint;
  su2double *Normal, Prod, Sens = 0.0, SensDim, Area;
  su2double Total_Sens_Geo_local = 0.0;
  Total_Sens_Geo = 0.0;

 // for (iMarker = 0; iMarker < nMarker; iMarker++){
    Sens_Geo[iMarker] = 0.0;
    /*--- Loop over boundary markers to select those for Euler walls and NS walls ---*/

    if(config->GetMarker_All_KindBC(iMarker) == EULER_WALL
       || config->GetMarker_All_KindBC(iMarker) == HEAT_FLUX
       || config->GetMarker_All_KindBC(iMarker) == ISOTHERMAL){

      for (iVertex = 0; iVertex < geometry->GetnVertex(iMarker); iVertex++){
        iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
        Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
        Prod = 0.0;
        Area = 0.0;
        for (iDim = 0; iDim < nDim; iDim++){
          /*--- retrieve the gradient calculated with AD -- */
          SensDim = node[iPoint]->GetSensitivity(iDim);

          /*--- calculate scalar product for projection onto the normal vector ---*/
          Prod += Normal[iDim]*SensDim;

          Area += Normal[iDim]*Normal[iDim];
        }

        Area = sqrt(Area);

        /*--- projection of the gradient
         *     calculated with AD onto the normal
         *     vector of the surface ---*/
        Sens = Prod/Area;

        /*--- Compute sensitivity for each surface point ---*/
        CSensitivity[iMarker][iVertex] = -Sens;
        if (geometry->node[iPoint]->GetDomain()){
          Sens_Geo[iMarker] += Sens*Sens;
        }
      }
      Total_Sens_Geo_local += sqrt(Sens_Geo[iMarker]);

    }
 // }

#ifdef HAVE_MPI
  SU2_MPI::Allreduce(&Total_Sens_Geo_local,&Total_Sens_Geo,1,MPI_DOUBLE,MPI_SUM, MPI_COMM_WORLD);
#else
  Total_Sens_Geo = Total_Sens_Geo_local;
#endif
}

void CDiscAdjSolver::Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config_container, unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output){
  bool dual_time_1st = (config_container->GetUnsteady_Simulation() == DT_STEPPING_1ST);
  bool dual_time_2nd = (config_container->GetUnsteady_Simulation() == DT_STEPPING_2ND);
  bool dual_time = (dual_time_1st || dual_time_2nd);
  su2double *solution_n, *solution_n1;
  unsigned long iPoint;
  unsigned short iVar;
  if (dual_time){
      for (iPoint = 0; iPoint<geometry->GetnPoint(); iPoint++){
          solution_n = node[iPoint]->GetSolution_time_n();
          solution_n1 = node[iPoint]->GetSolution_time_n1();

          for (iVar=0; iVar < nVar; iVar++){
              node[iPoint]->SetDual_Time_Derivative(iVar, solution_n[iVar]+node[iPoint]->GetDual_Time_Derivative_n(iVar));
              node[iPoint]->SetDual_Time_Derivative_n(iVar, solution_n1[iVar]);

            }

        }

    }
}
