/*!
 * \file COneShotSingleZoneDriver.cpp
 * \brief The main subroutines for driving one-shot single-zone problems.
 * \author L. Kusch
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

#include "../../include/drivers/COneShotSinglezoneDriver.hpp"
#include "../../include/iteration/CIterationFactory.hpp"

COneShotSinglezoneDriver::COneShotSinglezoneDriver(char* confFile,
                                                 unsigned short val_nZone,
                                                 SU2_Comm MPICommunicator) : CDiscAdjSinglezoneDriver(confFile, val_nZone,
                                                                                    MPICommunicator) {
  unsigned short iDV, jDV, iDV_Value;

  /*---------- One-shot works on all design variables - find total number of design variables ---------*/
  nDV_Total = 0;
  for (iDV = 0; iDV  < config->GetnDV(); iDV++){
    for (iDV_Value = 0; iDV_Value < config->GetnDV_Value(iDV); iDV_Value++){
      nDV_Total++;
    }
  }

  Gradient = new su2double[nDV_Total];
  Gradient_Old = new su2double[nDV_Total];

  ShiftedLagrangianGradient = new su2double[nDV_Total];
  ShiftedLagrangianGradient_Old = new su2double[nDV_Total];

  AugmentedLagrangianGradient = new su2double[nDV_Total];
  AugmentedLagrangianGradient_Old = new su2double[nDV_Total];

  DesignVarUpdate = new su2double[nDV_Total];
  DesignVariable = new su2double[nDV_Total];
  SearchDirection = new su2double[nDV_Total];
  activeset = new bool[nDV_Total];

  BFGS_Inv = new su2double*[nDV_Total];

  ConstrFunc = new su2double[config->GetnConstr()];
  multiplier = new su2double[config->GetnConstr()];
  ConstrFunc_Store = new su2double[config->GetnConstr()];
  BCheck_Inv = new su2double*[config->GetnConstr()];

  nBFGSmax = config->GetLimitedMemoryIter();
  nBFGS = 0;
  ykvec = new su2double*[nBFGSmax];
  skvec = new su2double*[nBFGSmax];

  BFGS_Init = config->GetBFGSInitValue();

  for (iDV = 0; iDV  < nDV_Total; iDV++){
    Gradient[iDV] = 0.0;
    Gradient_Old[iDV] = 0.0;
    ShiftedLagrangianGradient[iDV] = 0.0;
    ShiftedLagrangianGradient_Old[iDV] = 0.0;
    AugmentedLagrangianGradient[iDV] = 0.0;
    AugmentedLagrangianGradient_Old[iDV] = 0.0;
    DesignVarUpdate[iDV] = 0.0;
    DesignVariable[iDV] = 0.0;
    SearchDirection[iDV] = 0.0;
    activeset[iDV]=false;
    BFGS_Inv[iDV] = new su2double[nDV_Total];
    for (jDV = 0; jDV < nDV_Total; jDV++){
      BFGS_Inv[iDV][jDV] = 0.0;
      if (iDV==jDV) BFGS_Inv[iDV][jDV] = BFGS_Init;
    }
  }

  for (unsigned short iConstr = 0; iConstr  < config->GetnConstr(); iConstr++){
    ConstrFunc[iConstr] = 0.0;
    ConstrFunc_Store[iConstr] = 0.0;
    multiplier[iConstr] = config->GetMultiplierStart(iConstr);
    BCheck_Inv[iConstr] = new su2double[config->GetnConstr()];
    for (unsigned short jConstr = 0; jConstr  < config->GetnConstr(); jConstr++){
      BCheck_Inv[iConstr][jConstr] = 0.0;
    }
    BCheck_Inv[iConstr][iConstr] = config->GetMultiplierFactor(iConstr);
  }

  for (unsigned short iBFGS = 0; iBFGS < nBFGSmax; iBFGS++){
    ykvec[iBFGS] = new su2double[nDV_Total];
    skvec[iBFGS] = new su2double[nDV_Total];
    for (jDV = 0; jDV < nDV_Total; jDV++){
      ykvec[iBFGS][jDV] = 0.0;
      skvec[iBFGS][jDV] = 0.0;
    }
  }

  /*----- calculate values for bound projection algorithm -------*/
  lb=-config->GetBound()/config->GetDesignScale();
  ub=config->GetBound()/config->GetDesignScale();
  epsilon=(ub-lb)/2.0;

  /*---- calculate line search parameter ----*/
  cwolfeone= 1E-4*config->GetDesignScale(); //Achtung: vorher nicht auskommentiert

  grid_movement[ZONE_0][INST_0] = new CVolumetricMovement(geometry_container[ZONE_0][INST_0][MESH_0], config_container[ZONE_0]);
  surface_movement[ZONE_0] = new CSurfaceMovement();

  config->SetnTime_Iter(config->GetnIter());
}

COneShotSinglezoneDriver::~COneShotSinglezoneDriver(void){

  /*----- free allocated memory -------*/
  unsigned short iDV;
  for (iDV = 0; iDV  < nDV_Total; iDV++){
    delete [] BFGS_Inv[iDV];
  }
  delete [] BFGS_Inv;
  delete [] Gradient;
  delete [] Gradient_Old;
  delete [] ShiftedLagrangianGradient;
  delete [] ShiftedLagrangianGradient_Old;

  delete [] AugmentedLagrangianGradient;
  delete [] AugmentedLagrangianGradient_Old;

  delete [] DesignVarUpdate;
  delete [] DesignVariable;
  delete [] SearchDirection;
  delete [] activeset;

  delete [] ConstrFunc;
  delete [] multiplier;
  delete [] ConstrFunc_Store;

  for (unsigned short iBFGS = 0; iBFGS < nBFGSmax; iBFGS++){
    delete [] ykvec[iBFGS];
    delete [] skvec[iBFGS];
  }
  delete [] ykvec;
  delete [] skvec;

}

void COneShotSinglezoneDriver::RunOneShot(){

  su2double stepsize = config->GetStepSize();
  unsigned short maxcounter = config->GetOneShotMaxCounter();
  unsigned short whilecounter = 0;
  bool testLagrange = config->GetOneShotLagrange(); //Lagrange function includes all updates and not only design update if testLagrange is set to false
  bool descent = true;
  bool partstep = config->GetOneShotPartStep();
  cout << "log10Adjoint[RMS Density]: " << log10(solver[ADJFLOW_SOL]->GetRes_RMS(0))<<std::endl;
  cout << "log10Primal[RMS Density]: " << log10(solver[FLOW_SOL]->GetRes_RMS(0))<<std::endl;

  /*--- Store the old solution and the old design for line search ---*/
  if(!testLagrange) solver[ADJFLOW_SOL]->StoreSolution();
  if(testLagrange) solver[ADJFLOW_SOL]->StoreFormerSolution();
  solver[ADJFLOW_SOL]->StoreMeshPoints(config, geometry);
  
  /*--- This is the line search loop that is only called once, if no update is performed ---*/
  do{
    if(TimeIter>config->GetOneShotStart()&&TimeIter<config->GetOneShotStop()){
      if(!CheckDescent()&&whilecounter==0&&config->GetCheckDescent()){
          //whilecounter=maxcounter-1;
          ComputeNegativeSearchDirection();
          std::cout<<"No Descent Direction!..."<<CheckDescent()<<std::endl;
          descent = false;
      }
      if(whilecounter>0){
        //Armijo line search (halfen step)
        stepsize=stepsize*0.5;
        if(whilecounter==maxcounter && partstep){
          stepsize = config->GetStepSize();
        }
        /*---Load the old design for line search---*/
        solver[ADJFLOW_SOL]->LoadMeshPoints(config, geometry);
      }else{
        if(!testLagrange && !partstep && config->GetnConstr()>0) UpdateMultiplier(1.0);
      }
      if(partstep) UpdateMultiplier(stepsize);
      if(!partstep||(whilecounter==0)){
        //Load the old solution for line search (either y_k or y_k-1)
        solver[ADJFLOW_SOL]->LoadSolution();
      }else{
        //Update the old solution with stepsize and load it
        solver[ADJFLOW_SOL]->LoadSolutionStep(stepsize);
        solver[ADJFLOW_SOL]->StoreSolution();
      }
      /*--- Do a design update based on the search direction (mesh deformation with stepsize) ---*/
      if (whilecounter!=maxcounter || (!config->GetZeroStep())) ComputeDesignVarUpdate(stepsize);
      else ComputeDesignVarUpdate(0.0);
      config->SetKind_SU2(SU2_COMPONENT::SU2_DEF); // set SU2_DEF as the solver
      SurfaceDeformation(surface_movement[ZONE_0], grid_movement[ZONE_0][INST_0]);
      config->SetKind_SU2(SU2_COMPONENT::SU2_CFD); // set SU2_CFD as the solver
    }
    /*--- Do a primal and adjoint update ---*/
    if(!testLagrange || (TimeIter>config->GetOneShotStart()&&TimeIter<config->GetOneShotStop())){
      PrimalDualStep();
      CalculateLagrangian(true);
    }
    whilecounter++;
  }
  while(TimeIter>config->GetOneShotStart()&&TimeIter<config->GetOneShotStop()&&(!CheckFirstWolfe())&&whilecounter<maxcounter+1);
  if (whilecounter==maxcounter+1){
    descent = false;
    if(config->GetZeroStep()) stepsize = 0.0;
  }
  
  std::cout<<"Line search information: "<<Lagrangian<<" "<<ObjFunc<<" "<<stepsize<<std::endl;
  if(testLagrange){

    if(TimeIter>config->GetOneShotStart()&&TimeIter<config->GetOneShotStop()){
      //Evaluate sensitivity of augmented Lagrangian for an update only in the design
      /*--- N_u ---*/
      solver[ADJFLOW_SOL]->ResetSensitivityLagrangian(geometry);
      solver[ADJFLOW_SOL]->UpdateSensitivityLagrangian(geometry,1.0);

      /*--- Alpha*Deltay^T*G_u ---*/
      ComputeAlphaTerm();
      solver[ADJFLOW_SOL]->UpdateSensitivityLagrangian(geometry,config->GetOneShotAlpha());

      /*--- Beta*DeltaBary^T*N_yu ---*/
      solver[ADJFLOW_SOL]->LoadSolution();
      solver[ADJFLOW_SOL]->UpdateStateVariable(config);

      ComputeBetaTerm();
      solver[ADJFLOW_SOL]->SetFiniteDifferenceSens(geometry, config);
      solver[ADJFLOW_SOL]->UpdateSensitivityLagrangian(geometry,config->GetOneShotBeta());

      /*--- Projection of the gradient L_u---*/
      solver[ADJFLOW_SOL]->SetGeometrySensitivityLagrangian(geometry); //Lagrangian

      ProjectMeshSensitivities();
      SetAugmentedLagrangianGradient();
    }

    solver[ADJFLOW_SOL]->ShiftFormerSolution();
    solver[ADJFLOW_SOL]->LoadSolution();
    
    UpdateMultiplier(1.0);
    PrimalDualStep();
    CalculateLagrangian(true);
  }
  if(partstep){
    solver[ADJFLOW_SOL]->ShiftStoreSolution();
    solver[ADJFLOW_SOL]->StoreSolutionDelta();
   
  }
  if(TimeIter>=config->GetOneShotStart()&&TimeIter<config->GetOneShotStop()){

    /*--- Update design variable ---*/
    UpdateDesignVariable();
    if(config->GetnConstr()>0) StoreConstrFunction();

    /*--- N_u ---*/
    solver[ADJFLOW_SOL]->SaveSensitivity(geometry);
    solver[ADJFLOW_SOL]->StoreSaveSolution();
    solver[ADJFLOW_SOL]->ResetSensitivityLagrangian(geometry);
    solver[ADJFLOW_SOL]->UpdateSensitivityLagrangian(geometry,1.0);   

    if(config->GetnConstr()>0 && !config->GetConstPrecond()) ComputePreconditioner();
    /*--- Alpha*Deltay^T*G_u ---*/
    ComputeAlphaTerm();
    solver[ADJFLOW_SOL]->UpdateSensitivityLagrangian(geometry,config->GetOneShotAlpha());

    /*--- Beta*DeltaBary^T*N_yu ---*/
    solver[ADJFLOW_SOL]->LoadSolution();
    solver[ADJFLOW_SOL]->UpdateStateVariable(config);

    ComputeBetaTerm();
    solver[ADJFLOW_SOL]->SetFiniteDifferenceSens(geometry, config);
    solver[ADJFLOW_SOL]->UpdateSensitivityLagrangian(geometry,config->GetOneShotBeta());
    if(!partstep) solver[ADJFLOW_SOL]->LoadSaveSolution();

    /*--- Projection of the gradient N_u---*/
    solver[ADJFLOW_SOL]->SetGeometrySensitivityGradient(geometry);
    ProjectMeshSensitivities();
    SetShiftedLagrangianGradient();

    /*--- Projection of the gradient L_u---*/
    solver[ADJFLOW_SOL]->SetGeometrySensitivityLagrangian(geometry); //Lagrangian

    ProjectMeshSensitivities();
    if (!testLagrange){ SetAugmentedLagrangianGradient(); }

    /*--- Use N_u to compute the active set (bound constraints) ---*/
    ComputeActiveSet();

    /*--- Do a BFGS update to approximate the inverse preconditioner ---*/
    if(TimeIter>config->GetOneShotStart() && !config->GetLimitedMemory()) BFGSUpdate();
    if(config->GetLimitedMemory()) LBFGSUpdate();

    if (testLagrange){ SetAugmentedLagrangianGradient(); }
    /*--- Compute the search direction for the line search procedure ---*/
    ComputeSearchDirection();

    StoreLagrangianInformation();
  }
}

void COneShotSinglezoneDriver::Run(){
 if (config->GetBoolPiggyBack()) RunPiggyBack();
 else if (config->GetBoolQuasiNewton()) RunBFGS();
 else RunOneShot();
}

void COneShotSinglezoneDriver::RunBFGS(){
  su2double stepsize = config->GetStepSize();
  unsigned short maxcounter = config->GetOneShotMaxCounter();
  unsigned short whilecounter = 0;
  unsigned short nConvIter = 300; //Achtung vorher 500 // dann 2000
  unsigned short iterCount;

  if(TimeIter>config->GetOneShotStart()){
      /*--- Store the old solution and the old design for line search ---*/
      solver[ADJFLOW_SOL]->StoreSolution();
      solver[ADJFLOW_SOL]->StoreMeshPoints(config, geometry);
  }

  do{
    if(whilecounter>0){
      /*---Load the old solution and the old design for line search---*/
      solver[ADJFLOW_SOL]->LoadSolution();
      solver[ADJFLOW_SOL]->LoadMeshPoints(config, geometry);
    }
    else{
      if(config->GetnConstr()>0) UpdateMultiplier(1.0);
    }
    if(TimeIter>config->GetOneShotStart()){
      /*--- Do a design update based on the search direction (mesh deformation with stepsize) ---*/
      ComputeDesignVarUpdate(stepsize);
      config->SetKind_SU2(SU2_COMPONENT::SU2_DEF); // set SU2_DEF as the solver
      SurfaceDeformation(surface_movement[ZONE_0], grid_movement[ZONE_0][INST_0]);
      config->SetKind_SU2(SU2_COMPONENT::SU2_CFD); // set SU2_CFD as the solver
    }
    for(iterCount=0;iterCount<nConvIter;iterCount++){
      PrimalDualStep();
    }
    CalculateLagrangian(false);
    cout << "log10Adjoint[RMS Density]: " << log10(solver[ADJFLOW_SOL]->GetRes_RMS(0))<<std::endl;
    cout << "log10Primal[RMS Density]: " << log10(solver[FLOW_SOL]->GetRes_RMS(0))<<std::endl;
    stepsize=stepsize*0.5;
    whilecounter++;
  }
  while(TimeIter>config->GetOneShotStart()&&(!CheckFirstWolfe())&&whilecounter<maxcounter+1);

  if(TimeIter>=config->GetOneShotStart()){
      /*--- Update design variable ---*/
      UpdateDesignVariable();
      if(config->GetnConstr()>0) StoreConstrFunction();

      /*---N_u---*/
      /*--- Projection of the gradient ---*/
      solver[ADJFLOW_SOL]->SaveSensitivity(geometry);
      solver[ADJFLOW_SOL]->SetGeometrySensitivityGradient(geometry);
      ProjectMeshSensitivities();

      SetShiftedLagrangianGradient();
      SetAugmentedLagrangianGradient();

      ComputeActiveSet();
      if(TimeIter>config->GetOneShotStart()) BFGSUpdate();

      /*--- Compute the search direction for the line search procedure ---*/
      ComputeSearchDirection();

      StoreLagrangianInformation();
  }
}

void COneShotSinglezoneDriver::PrimalDualStep(){

  /*--- Note: Unsteady cases not applicable to the one-shot method yet! ---*/

  SetRecording(RECORDING::SOLUTION_AND_MESH);

  /*--- Initialize the adjoint of the output variables of the iteration with the adjoint solution
   *    of the previous iteration. The values are passed to the AD tool. ---*/

  config->SetInnerIter(0);

  iteration->InitializeAdjoint(solver_container, geometry_container, config_container, ZONE_0, INST_0);

  /*--- Initialize the adjoint of the objective function with 1.0. ---*/

  SetAdj_ObjFunction();
  SetAdj_ConstrFunction(multiplier);

  /*--- Interpret the stored information by calling the corresponding routine of the AD tool. ---*/

  AD::ComputeAdjoint();

  /*--- Extract the computed adjoint values of the input variables and store them for the next iteration. ---*/
  iteration->Iterate(output_container[ZONE_0], integration_container, geometry_container,
                                          solver_container, numerics_container, config_container,
                                          surface_movement, grid_movement, FFDBox, ZONE_0, INST_0);

  /*--- Extract the computed sensitivity values. ---*/
  solver[ADJFLOW_SOL]->SetSensitivity(geometry,config);

  /*--- Clear the stored adjoint information to be ready for a new evaluation. ---*/

  AD::ClearAdjoints();
}

void COneShotSinglezoneDriver::RunPiggyBack() {

  PrimalDualStep();

  cout << "log10Adjoint[RMS Density]: " << log10(solver[ADJFLOW_SOL]->GetRes_RMS(0))<<std::endl;
  cout << "log10Primal[RMS Density]: " << log10(solver[FLOW_SOL]->GetRes_RMS(0))<<std::endl;

  if(TimeIter>config->GetOneShotStart()){

    /*---N_u---*/
    solver[ADJFLOW_SOL]->SaveSensitivity(geometry);
    solver[ADJFLOW_SOL]->SetGeometrySensitivityGradient(geometry);
    ProjectMeshSensitivities();
  }
}

void COneShotSinglezoneDriver::SetRecording(RECORDING kind_recording){

#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  AD::Reset();

  solver[ADJFLOW_SOL]->SetRecording(geometry, config);

  if (config->GetKind_Solver() == DISC_ADJ_RANS && !config->GetFrozen_Visc_Disc()) {
    solver[ADJTURB_SOL]->SetRecording(geometry, config);
  }
  if (config->GetKind_Solver() == ONE_SHOT_RANS && !config->GetFrozen_Visc_Disc()) {
    solver[ADJTURB_SOL]->SetRecording(geometry, config);
  }


  /*---Enable recording and register input of the flow iteration (conservative variables or node coordinates) --- */

  if (kind_recording != RECORDING::CLEAR_INDICES){

    AD::StartRecording();
    iteration->RegisterInput(solver_container, geometry_container, config_container, ZONE_0, INST_0, kind_recording);

  }

  iteration->SetDependencies(solver_container, geometry_container, numerics_container, config_container, ZONE_0, INST_0, kind_recording);

  /*--- Do one iteration of the direct flow solver ---*/

  DirectRun(kind_recording);

  RecordingState = kind_recording;

  iteration->RegisterOutput(solver_container, geometry_container, config_container, ZONE_0, INST_0);

  /*--- Extract the objective function and store it --- */

  SetObjFunction();
  SetConstrFunction();

  AD::StopRecording();

}

void COneShotSinglezoneDriver::SetProjection_AD(CSurfaceMovement *surface_movement, su2double* Gradient){

  su2double DV_Value, *VarCoord, Sensitivity, my_Gradient, localGradient;//, *Normal, Area = 0.0;
  unsigned short iDV_Value = 0, iMarker, nMarker, iDim, nDim, iDV, nDV, nDV_Value;
  unsigned long iVertex, nVertex, iPoint;
  unsigned long nDV_Count = 0;

  nMarker = config->GetnMarker_All();
  nDim    = geometry->GetnDim();
  nDV     = config->GetnDV();

  VarCoord = NULL;

  AD::Reset();

  /*--- Discrete adjoint gradient computation ---*/

  /*--- Start recording of operations ---*/

  AD::StartRecording();

  /*--- Register design variables as input and set them to zero
   * (since we want to have the derivative at alpha = 0, i.e. for the current design) ---*/

  for (iDV = 0; iDV < nDV; iDV++){

    nDV_Value =  config->GetnDV_Value(iDV);

    for (iDV_Value = 0; iDV_Value < nDV_Value; iDV_Value++){

      /*--- Initilization with su2double resets the index ---*/

      DV_Value = 0.0;

      AD::RegisterInput(DV_Value);

      config->SetDV_Value(iDV, iDV_Value, DV_Value);

    }
  }

  /*--- Call the surface deformation routine ---*/

  surface_movement->SetSurface_Deformation(geometry, config);

  /*--- Stop the recording --- */

  AD::StopRecording();

  /*--- Create a structure to identify points that have been already visited.
   * We need that to make sure to set the sensitivity of surface points only once
   *  (Markers share points, so we would visit them more than once in the loop over the markers below) ---*/

  bool* visited = new bool[geometry->GetnPoint()];
  for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++){
    visited[iPoint] = false;
  }

  /*--- Initialize the derivatives of the output of the surface deformation routine
   * with the discrete adjoints from the CFD solution ---*/

  for (iMarker = 0; iMarker < nMarker; iMarker++) {
    if (config->GetMarker_All_DV(iMarker) == YES) {
      nVertex = geometry->nVertex[iMarker];
      for (iVertex = 0; iVertex <nVertex; iVertex++) {
        iPoint      = geometry->vertex[iMarker][iVertex]->GetNode();
        if (!visited[iPoint]){
          VarCoord    = geometry->vertex[iMarker][iVertex]->GetVarCoord();
         /* Normal      = geometry->vertex[iMarker][iVertex]->GetNormal();

          Area = 0.0;
          for (iDim = 0; iDim < nDim; iDim++){
            Area += Normal[iDim]*Normal[iDim];
          }
          Area = sqrt(Area);*/

          for (iDim = 0; iDim < nDim; iDim++){
           // if (config->GetDiscrete_Adjoint()){
              Sensitivity = geometry->GetSensitivity(iPoint, iDim);
           /* } else {
              Sensitivity = -Normal[iDim]*geometry->vertex[iMarker][iVertex]->GetAuxVar()/Area;
            }*/
            SU2_TYPE::SetDerivative(VarCoord[iDim], SU2_TYPE::GetValue(Sensitivity));
          }
          visited[iPoint] = true;
        }
      }
    }
  }

  delete [] visited;

  /*--- Compute derivatives and extract gradient ---*/

  AD::ComputeAdjoint();

  for (iDV = 0; iDV  < nDV; iDV++){
    nDV_Value =  config->GetnDV_Value(iDV);

    for (iDV_Value = 0; iDV_Value < nDV_Value; iDV_Value++){
      DV_Value = config->GetDV_Value(iDV, iDV_Value);
      my_Gradient = SU2_TYPE::GetDerivative(DV_Value);
#ifdef HAVE_MPI
    SU2_MPI::Allreduce(&my_Gradient, &localGradient, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
      localGradient = my_Gradient;
#endif
      /*--- Angle of Attack design variable (this is different,tor<su2double>[1];
        vector<su2double> *Variable_Airfoil = new vector<su2double>[1];
       the value comes form the input file) ---*/

/*      if ((config->GetDesign_Variable(iDV) == ANGLE_OF_ATTACK) ||
          (config->GetDesign_Variable(iDV) == FFD_ANGLE_OF_ATTACK))  {
        Gradient[iDV][iDV_Value] = config->GetAoA_Sens();
      }*/

      Gradient[nDV_Count] = localGradient;
      nDV_Count++;
      std::cout<<std::setprecision(9)<<localGradient<<", ";
    }
  }
  std::cout<<std::endl;
  AD::Reset();

}

void COneShotSinglezoneDriver::SetProjection_FD(CSurfaceMovement *surface_movement, su2double* Gradient){
  unsigned short iDV, nDV, iFFDBox, nDV_Value, iMarker, iDim;
  unsigned long iVertex, iPoint;
  su2double delta_eps, my_Gradient, localGradient, *Normal, dS, *VarCoord, Sensitivity,
  dalpha[3], deps[3], dalpha_deps;
  bool *UpdatePoint, MoveSurface, Local_MoveSurface;
  CFreeFormDefBox **FFDBox;
    unsigned long nDV_Count = 0;

  int rank = SU2_MPI::GetRank();

  nDV = config->GetnDV();

  /*--- Boolean controlling points to be updated ---*/

  UpdatePoint = new bool[geometry->GetnPoint()];

  /*--- Definition of the FFD deformation class ---*/

  unsigned short nFFDBox = MAX_NUMBER_FFD;
  FFDBox = new CFreeFormDefBox*[nFFDBox];
  for (iFFDBox = 0; iFFDBox < MAX_NUMBER_FFD; iFFDBox++) FFDBox[iFFDBox] = NULL;

  for (iDV = 0; iDV  < nDV; iDV++){
    nDV_Value = config->GetnDV_Value(iDV);
    if (nDV_Value != 1){
      SU2_MPI::Error("The projection using finite differences currently only supports a fixed direction of movement for FFD points.", CURRENT_FUNCTION);
    }
  }

  /*--- Continuous adjoint gradient computation ---*/

  for (iDV = 0; iDV < nDV; iDV++) {
    config->SetDV_Value(iDV,0, 1E-4);
    MoveSurface = true;
    Local_MoveSurface = true;

    /*--- Free Form deformation based ---*/

    if ((config->GetDesign_Variable(iDV) == FFD_CONTROL_POINT_2D) ||
        (config->GetDesign_Variable(iDV) == FFD_CAMBER_2D) ||
        (config->GetDesign_Variable(iDV) == FFD_THICKNESS_2D) ||
        (config->GetDesign_Variable(iDV) == FFD_TWIST_2D) ||
        (config->GetDesign_Variable(iDV) == FFD_CONTROL_POINT) ||
        (config->GetDesign_Variable(iDV) == FFD_NACELLE) ||
        (config->GetDesign_Variable(iDV) == FFD_GULL) ||
        (config->GetDesign_Variable(iDV) == FFD_TWIST) ||
        (config->GetDesign_Variable(iDV) == FFD_ROTATION) ||
        (config->GetDesign_Variable(iDV) == FFD_CAMBER) ||
        (config->GetDesign_Variable(iDV) == FFD_THICKNESS) ||
        (config->GetDesign_Variable(iDV) == FFD_ANGLE_OF_ATTACK)) {

      /*--- Read the FFD information in the first iteration ---*/

      if (iDV == 0) {


        /*--- Read the FFD information from the grid file ---*/

        surface_movement->ReadFFDInfo(geometry, config, FFDBox, config->GetMesh_FileName());

        /*--- If the FFDBox was not defined in the input file ---*/
        if (!surface_movement->GetFFDBoxDefinition()) {
          SU2_MPI::Error("The input grid doesn't have the entire FFD information!", CURRENT_FUNCTION);
        }

        for (iFFDBox = 0; iFFDBox < surface_movement->GetnFFDBox(); iFFDBox++) {

          surface_movement->CheckFFDDimension(geometry, config, FFDBox[iFFDBox], iFFDBox);

          surface_movement->CheckFFDIntersections(geometry, config, FFDBox[iFFDBox], iFFDBox);

        }

      }

      /*--- Apply the control point change ---*/

      MoveSurface = false;

      for (iFFDBox = 0; iFFDBox < surface_movement->GetnFFDBox(); iFFDBox++) {

        /*--- Reset FFD box ---*/
        switch (config->GetDesign_Variable(iDV) ) {
          case FFD_CONTROL_POINT_2D : Local_MoveSurface = surface_movement->SetFFDCPChange_2D(geometry, config, FFDBox[iFFDBox], FFDBox, iDV, true); break;
          case FFD_CAMBER_2D :        Local_MoveSurface = surface_movement->SetFFDCamber_2D(geometry, config, FFDBox[iFFDBox], FFDBox, iDV, true); break;
          case FFD_THICKNESS_2D :     Local_MoveSurface = surface_movement->SetFFDThickness_2D(geometry, config, FFDBox[iFFDBox], FFDBox, iDV, true); break;
          case FFD_TWIST_2D :         Local_MoveSurface = surface_movement->SetFFDTwist_2D(geometry, config, FFDBox[iFFDBox], FFDBox, iDV, true); break;
          case FFD_CONTROL_POINT :    Local_MoveSurface = surface_movement->SetFFDCPChange(geometry, config, FFDBox[iFFDBox], FFDBox, iDV, true); break;
          case FFD_NACELLE :          Local_MoveSurface = surface_movement->SetFFDNacelle(geometry, config, FFDBox[iFFDBox], FFDBox, iDV, true); break;
          case FFD_GULL :             Local_MoveSurface = surface_movement->SetFFDGull(geometry, config, FFDBox[iFFDBox], FFDBox, iDV, true); break;
          case FFD_TWIST :            Local_MoveSurface = surface_movement->SetFFDTwist(geometry, config, FFDBox[iFFDBox], FFDBox, iDV, true); break;
          case FFD_ROTATION :         Local_MoveSurface = surface_movement->SetFFDRotation(geometry, config, FFDBox[iFFDBox], FFDBox, iDV, true); break;
          case FFD_CAMBER :           Local_MoveSurface = surface_movement->SetFFDCamber(geometry, config, FFDBox[iFFDBox], FFDBox, iDV, true); break;
          case FFD_THICKNESS :        Local_MoveSurface = surface_movement->SetFFDThickness(geometry, config, FFDBox[iFFDBox], FFDBox, iDV, true); break;
          case FFD_CONTROL_SURFACE :  Local_MoveSurface = surface_movement->SetFFDControl_Surface(geometry, config, FFDBox[iFFDBox], FFDBox, iDV, true); break;
       //   case FFD_ANGLE_OF_ATTACK :  Gradient[iDV][0] = config->GetAoA_Sens(); break;
        }

        /*--- Recompute cartesian coordinates using the new control points position ---*/

        if (Local_MoveSurface) {
          MoveSurface = true;
          surface_movement->SetCartesianCoord(geometry, config, FFDBox[iFFDBox], iFFDBox, true);
        }

      }

    }

    /*--- Hicks Henne design variable ---*/

    else if (config->GetDesign_Variable(iDV) == HICKS_HENNE) {
      surface_movement->SetHicksHenne(geometry, config, iDV, true);
    }

    /*--- Surface bump design variable ---*/

    else if (config->GetDesign_Variable(iDV) == SURFACE_BUMP) {
      surface_movement->SetSurface_Bump(geometry, config, iDV, true);
    }

    /*--- Kulfan (CST) design variable ---*/

    else if (config->GetDesign_Variable(iDV) == CST) {
      surface_movement->SetCST(geometry, config, iDV, true);
    }

    /*--- Displacement design variable ---*/

    else if (config->GetDesign_Variable(iDV) == TRANSLATION) {
      surface_movement->SetTranslation(geometry, config, iDV, true);
    }

    /*--- Angle of Attack design variable ---*/

 /*   else if (config->GetDesign_Variable(iDV) == ANGLE_OF_ATTACK) {
      Gradient[iDV][0] = config->GetAoA_Sens();
    }*/

    /*--- Scale design variable ---*/

    else if (config->GetDesign_Variable(iDV) == SCALE) {
      surface_movement->SetScale(geometry, config, iDV, true);
    }

    /*--- Rotation design variable ---*/

    else if (config->GetDesign_Variable(iDV) == ROTATION) {
      surface_movement->SetRotation(geometry, config, iDV, true);
    }

    /*--- NACA_4Digits design variable ---*/

    else if (config->GetDesign_Variable(iDV) == NACA_4DIGITS) {
      surface_movement->SetNACA_4Digits(geometry, config);
    }

    /*--- Parabolic design variable ---*/

    else if (config->GetDesign_Variable(iDV) == PARABOLIC) {
      surface_movement->SetParabolic(geometry, config);
    }

    /*--- Design variable not implement ---*/

    else {
      if (rank == MASTER_NODE)
        cout << "Design Variable not implement yet" << endl;
    }

    /*--- Load the delta change in the design variable (finite difference step). ---*/

    if ((config->GetDesign_Variable(iDV) != ANGLE_OF_ATTACK) &&
        (config->GetDesign_Variable(iDV) != FFD_ANGLE_OF_ATTACK)) {

      /*--- If the Angle of attack is not involved, reset the value of the gradient ---*/

      my_Gradient = 0.0; Gradient[nDV_Count] = 0.0;

      if (MoveSurface) {

        delta_eps = config->GetDV_Value(iDV);

        for (iPoint = 0; iPoint < geometry->GetnPoint(); iPoint++)
          UpdatePoint[iPoint] = true;

        for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
          if (config->GetMarker_All_DV(iMarker) == YES) {
            for (iVertex = 0; iVertex < geometry->nVertex[iMarker]; iVertex++) {

              iPoint = geometry->vertex[iMarker][iVertex]->GetNode();
              if ((iPoint < geometry->GetnPointDomain()) && UpdatePoint[iPoint]) {

                Normal = geometry->vertex[iMarker][iVertex]->GetNormal();
                su2double Prod = 0.0;
                VarCoord = geometry->vertex[iMarker][iVertex]->GetVarCoord();

                dS = 0.0;
                for (iDim = 0; iDim < geometry->GetnDim(); iDim++) {
                  Sensitivity = geometry->GetSensitivity(iPoint,iDim);
                  Prod+=Normal[iDim]*Sensitivity;
                  dS += Normal[iDim]*Normal[iDim];
                  deps[iDim] = VarCoord[iDim] / delta_eps;
                }
                dS = sqrt(dS);
                Sensitivity = -Prod/dS;

                dalpha_deps = 0.0;
                for (iDim = 0; iDim < geometry->GetnDim(); iDim++) {
                  dalpha[iDim] = Normal[iDim] / dS;
                  dalpha_deps -= dalpha[iDim]*deps[iDim];
                }

                my_Gradient += Sensitivity*dalpha_deps;
                UpdatePoint[iPoint] = false;
              }
            }
          }
        }

      }

#ifdef HAVE_MPI
    SU2_MPI::Allreduce(&my_Gradient, &localGradient, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
    localGradient = my_Gradient;
#endif
    Gradient[nDV_Count] = localGradient;
    nDV_Count++;
    std::cout<<std::setprecision(9)<<localGradient<<", ";
    }
  }
  std::cout<<std::endl;

  /*--- Delete memory for parameterization. ---*/

  if (FFDBox != NULL) {
    for (iFFDBox = 0; iFFDBox < MAX_NUMBER_FFD; iFFDBox++) {
      if (FFDBox[iFFDBox] != NULL) {
        delete FFDBox[iFFDBox];
      }
    }
    delete [] FFDBox;
  }

  delete [] UpdatePoint;
}

void COneShotSinglezoneDriver::SurfaceDeformation(CSurfaceMovement *surface_movement_cl, CVolumetricMovement *grid_movement_cl){

  unsigned short iMarker, iDV, iDV_Value, nDV_Value;
  bool allmoving=true;
  unsigned long nDV_Count = 0;

  for (iDV = 0; iDV < config->GetnDV(); iDV++) {
    nDV_Value =  config->GetnDV_Value(iDV);

    for (iDV_Value = 0; iDV_Value < nDV_Value; iDV_Value++){
      config->SetDV_Value(iDV,iDV_Value, config->GetDesignScale()*DesignVarUpdate[nDV_Count]);
      nDV_Count++;
    }
  }

  /*--- Surface grid deformation using design variables ---*/

  surface_movement_cl->SetSurface_Deformation(geometry, config);

  /*--- For scale, translation and rotation if all boundaries are moving they are set via volume method
   * Otherwise, the surface deformation has been set already in SetSurface_Deformation.  --- */
  /*--- Loop over markers, set flag to false if any are not moving ---*/
  for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++){
    if (config->GetMarker_All_DV(iMarker) == NO)
      allmoving = false;
  }

  /*--- Volumetric grid deformation/transformations ---*/

  if (config->GetDesign_Variable(0) == SCALE && allmoving) {

    grid_movement_cl->SetVolume_Scaling(geometry, config, false);

  } else if (config->GetDesign_Variable(0) == TRANSLATION && allmoving) {

    grid_movement_cl->SetVolume_Translation(geometry, config, false);

  } else if (config->GetDesign_Variable(0) == ROTATION && allmoving) {

    grid_movement_cl->SetVolume_Rotation(geometry, config, false);

  } else if (config->GetDesign_Variable(0) != FFD_SETTING) {

    grid_movement_cl->SetVolume_Deformation(geometry, config, false);

  }

}

void COneShotSinglezoneDriver::BFGSUpdate(){
  unsigned long iDV, jDV, kDV, lDV;

    su2double *yk, *sk;
    su2double vk=0;
    su2double normyk=0;
    su2double normsk=0;
    su2double helper1=0;
    su2double helper2=0;

    yk=new su2double[nDV_Total];
    sk=new su2double[nDV_Total];
    for (iDV = 0; iDV < nDV_Total; iDV++){
      yk[iDV]=ProjectionSet(iDV, AugmentedLagrangianGradient[iDV]-AugmentedLagrangianGradient_Old[iDV], false);
      sk[iDV]=ProjectionSet(iDV, DesignVarUpdate[iDV], false);
      helper1+=sk[iDV]*AugmentedLagrangianGradient[iDV];
      helper2+=sk[iDV]*AugmentedLagrangianGradient_Old[iDV];
      vk+=yk[iDV]*sk[iDV];
      normyk+=yk[iDV]*yk[iDV];
      normsk+=sk[iDV]*sk[iDV];
    }

    std::cout<<"yk*Gradient "<<helper1<<" yk*GradientOld "<<helper2<<std::endl;
    std::cout<<"vk "<<vk<<std::endl;
    std::cout<<"normsk "<<normsk<<", normyk "<<normyk<<", vk/normsk "<<vk/normsk<<std::endl;
    if (vk>0){
      su2double** MatA;
      MatA=new su2double*[nDV_Total];
      for (iDV=0;iDV<nDV_Total;iDV++){
        MatA[iDV]=new su2double[nDV_Total];
        for (jDV=0;jDV<nDV_Total;jDV++){
            MatA[iDV][jDV]=0.0;
        }
      }
      for (iDV=0;iDV<nDV_Total;iDV++){
        for (jDV=0;jDV<nDV_Total;jDV++){
          MatA[iDV][jDV]=ProjectionPAP(iDV,jDV,BFGS_Inv[iDV][jDV],false)+(1.0/vk)*sk[iDV]*sk[jDV];
          for (kDV=0; kDV<nDV_Total; kDV++){
            MatA[iDV][jDV]+=-(1.0/vk)*sk[iDV]*ProjectionPAP(kDV,jDV,BFGS_Inv[kDV][jDV],false)*yk[kDV]-(1.0/vk)*sk[jDV]*ProjectionPAP(iDV,kDV,BFGS_Inv[iDV][kDV],false)*yk[kDV];
            for (lDV=0; lDV<nDV_Total; lDV++){
              MatA[iDV][jDV]+=(1.0/vk)*(1.0/vk)*sk[iDV]*sk[jDV]*yk[lDV]*ProjectionPAP(lDV,kDV,BFGS_Inv[lDV][kDV],false)*yk[kDV];
            }
          }
        }
      }
      for (iDV=0;iDV<nDV_Total;iDV++){
        for (jDV=0;jDV<nDV_Total;jDV++){
          BFGS_Inv[iDV][jDV]=MatA[iDV][jDV];
        }
      }
      for (iDV=0;iDV<nDV_Total;iDV++){
        delete [] MatA[iDV];
      }
      delete [] MatA;
      if(config->GetBFGSInit()){
        for (iDV=0;iDV<nDV_Total;iDV++){
          BFGS_Init = vk/normyk;
        }
      }

    }else{
      std::cout<<"!Attention: Hessian not positive definite - Reset Preconditioner!"<<std::endl;
      if(config->GetOSHessianIdentity()){
        for (iDV = 0; iDV < nDV_Total; iDV++){
          for (jDV = 0; jDV < nDV_Total; jDV++){
            BFGS_Inv[iDV][jDV]=0.0;
            if(iDV==jDV){ BFGS_Inv[iDV][jDV]=ProjectionSet(iDV,BFGS_Init,false); }
          }
        }
      }
    }
    delete [] yk;
    delete [] sk;
}

void COneShotSinglezoneDriver::LBFGSUpdate(){
  unsigned long iDV;
  su2double vk=0.0;

  if (nBFGS == nBFGSmax) {
    nBFGS--;
    for (unsigned short iBFGS=0; iBFGS<nBFGSmax-1; iBFGS++){
      for (iDV = 0; iDV < nDV_Total; iDV++){
        ykvec[iBFGS][iDV]= ykvec[iBFGS+1][iDV];
        skvec[iBFGS][iDV]= skvec[iBFGS+1][iDV];
      }
    }
  }
  for (iDV = 0; iDV < nDV_Total; iDV++){
    ykvec[nBFGS][iDV]=ProjectionSet(iDV, AugmentedLagrangianGradient[iDV]-AugmentedLagrangianGradient_Old[iDV], false);
    skvec[nBFGS][iDV]=ProjectionSet(iDV, DesignVarUpdate[iDV], false);
    SearchDirection[iDV]=-ShiftedLagrangianGradient[iDV];
    vk+=ykvec[nBFGS][iDV]*skvec[nBFGS][iDV];
  }
  if(vk>0){
    nBFGS = nBFGS + 1;
    if(config->GetBFGSInit()){
      su2double helper = 0.0;
      for (iDV=0;iDV<nDV_Total;iDV++){
        helper+=ykvec[nBFGS-1][iDV]*ykvec[nBFGS-1][iDV];
      }
      for (iDV=0;iDV<nDV_Total;iDV++){
        for (unsigned short jDV=0;jDV<nDV_Total;jDV++){
          BFGS_Inv[iDV][jDV]= (ykvec[nBFGS-1][iDV]*skvec[nBFGS-1][jDV])/helper;
        }
      }
    }
  }else{
    nBFGS = 0;
  }

  LBFGSUpdateRecursive(nBFGS);
}

void COneShotSinglezoneDriver::LBFGSUpdateRecursive(unsigned short nCounter){
  unsigned long iDV, jDV;
  su2double *helper = new su2double [nDV_Total];
  su2double alpha = 0.0;
  su2double alphahelpone = 0.0;
  su2double alphahelptwo = 0.0;
  for (iDV = 0; iDV < nDV_Total; iDV++){
    SearchDirection[iDV]=ProjectionSet(iDV,SearchDirection[iDV],false);
  }
  if(nCounter == 0){
    for (iDV=0;iDV<nDV_Total;iDV++){
      helper[iDV]=0.0;
      for (jDV=0;jDV<nDV_Total;jDV++){
        helper[iDV]+=BFGS_Inv[iDV][jDV]*SearchDirection[jDV];
      }
    }
    for (iDV=0;iDV<nDV_Total;iDV++){
      SearchDirection[iDV] = helper[iDV];
    }
  }
  else{
    for (iDV=0;iDV<nDV_Total;iDV++){
      ykvec[nCounter-1][iDV] = ProjectionSet(iDV, ykvec[nCounter-1][iDV], false);
      skvec[nCounter-1][iDV] = ProjectionSet(iDV, skvec[nCounter-1][iDV], false);
      alphahelpone+=skvec[nCounter-1][iDV]*SearchDirection[iDV];
      alphahelptwo+=ykvec[nCounter-1][iDV]*skvec[nCounter-1][iDV];
    }
    alpha = alphahelpone/alphahelptwo;
    for (iDV=0;iDV<nDV_Total;iDV++){
      SearchDirection[iDV] -= alpha*ykvec[nCounter-1][iDV];
    }
    LBFGSUpdateRecursive(nCounter-1);
    alphahelpone = 0.0;
    for (iDV=0;iDV<nDV_Total;iDV++){
      alphahelpone+=ykvec[nCounter-1][iDV]*SearchDirection[iDV];
    }
    for (iDV=0;iDV<nDV_Total;iDV++){
      SearchDirection[iDV] += (alpha - alphahelpone/alphahelptwo)*skvec[nCounter-1][iDV];
      SearchDirection[iDV] = ProjectionSet(iDV,SearchDirection[iDV],false);
    }
  }
  delete [] helper;
}

bool COneShotSinglezoneDriver::CheckFirstWolfe(){
  unsigned short iDV;
  su2double admissible_step = 0.0;
  for (iDV=0;iDV<nDV_Total;iDV++){
    //ShiftedLagrangianGradient is the gradient at the old iterate (for One_Shot it is N_u and not L_u)
    admissible_step += DesignVarUpdate[iDV]*ShiftedLagrangianGradient[iDV];
    //admissible_step += (-1.0/alpha)*DesignVarUpdate[iDV]*DesignVarUpdate[iDV]; //Option was set before, now: first option again!
  }
  admissible_step *= cwolfeone;
  std::cout<<"Wolfe: "<<Lagrangian_Old<<" "<<Lagrangian<<" "<<admissible_step<<std::endl;
  return (Lagrangian<=Lagrangian_Old+admissible_step);
}

void COneShotSinglezoneDriver::ComputeDesignVarUpdate(su2double stepsize){
  unsigned short iDV;
  for (iDV=0;iDV<nDV_Total;iDV++){
    DesignVarUpdate[iDV]=BoundProjection(DesignVariable[iDV]+stepsize*SearchDirection[iDV])-DesignVariable[iDV];
  }
}

void COneShotSinglezoneDriver::ComputeSearchDirection(){
  unsigned short iDV, jDV;
  for (iDV=0;iDV<nDV_Total;iDV++){
    if(!config->GetLimitedMemory()){
      SearchDirection[iDV]=0.0;
      for (jDV=0;jDV<nDV_Total;jDV++){
        SearchDirection[iDV]+= BFGS_Inv[iDV][jDV]*ProjectionSet(jDV,-ShiftedLagrangianGradient[jDV],false); //ProjectionPAP(iDV, jDV, BFGS_Inv[iDV][jDV], false)*(-ShiftedLagrangianGradient[jDV]);//
      }
    }
    SearchDirection[iDV]=-ProjectionSet(iDV, ShiftedLagrangianGradient[iDV],true)+ProjectionSet(iDV, SearchDirection[iDV], false); // + SearchDirection[iDV];
  }
}

void COneShotSinglezoneDriver::ComputeNegativeSearchDirection(){
  unsigned short iDV, jDV;
  for (iDV=0;iDV<nDV_Total;iDV++){
    SearchDirection[iDV]=0.0;
    for (jDV=0;jDV<nDV_Total;jDV++){
      SearchDirection[iDV]+=BFGS_Inv[iDV][jDV]*ProjectionSet(jDV,ShiftedLagrangianGradient[jDV],false);
    }
    SearchDirection[iDV]=ProjectionSet(iDV, ShiftedLagrangianGradient[iDV],true)+ProjectionSet(iDV, SearchDirection[iDV], false);
  }
}

bool COneShotSinglezoneDriver::CheckDescent(){
  unsigned short iDV;
  su2double product = 0.0;
  for (iDV=0;iDV<nDV_Total;iDV++){
    product+=SearchDirection[iDV]*ProjectionSet(iDV, AugmentedLagrangianGradient_Old[iDV], false);
  }
  return (product<=0.0);
}

void COneShotSinglezoneDriver::StoreLagrangianInformation(){
  unsigned short iDV;
  for (iDV=0; iDV<nDV_Total; iDV++){
    AugmentedLagrangianGradient_Old[iDV] = AugmentedLagrangianGradient[iDV];
  }
  Lagrangian_Old = Lagrangian;
}

void COneShotSinglezoneDriver::UpdateDesignVariable(){
  unsigned short iDV;
  std::cout<<"Design Variable: ";
  for (iDV=0; iDV<nDV_Total; iDV++){
    DesignVariable[iDV] += DesignVarUpdate[iDV];
    std::cout<<DesignVariable[iDV]<<" ";
  }
  std::cout<<std::endl;
}

void COneShotSinglezoneDriver::CalculateLagrangian(bool augmented){
  std::cout<<"Objective function value: "<<ObjFunc<<std::endl;
  Lagrangian = 0.0;
  Lagrangian += ObjFunc; //TODO use for BFGS either only objective function or normal Lagrangian
  for (unsigned short iConstr = 0; iConstr < config->GetnConstr(); iConstr++){
    if(iConstr==0) std::cout<<"Constraint function value: ";
    Lagrangian += ConstrFunc[iConstr]*multiplier[iConstr];
    std::cout<<ConstrFunc[iConstr]<<" ";
    if (iConstr==config->GetnConstr()-1) std::cout<<std::endl;
  }
  if(augmented){
    su2double helper=0.0;
    for (unsigned short iConstr = 0; iConstr < config->GetnConstr(); iConstr++)
    {
      helper += ConstrFunc[iConstr]*ConstrFunc[iConstr];
    }
    Lagrangian += helper*(config->GetOneShotAlpha()/2);
    Lagrangian+=solver[ADJFLOW_SOL]->CalculateLagrangianPart(config, augmented);
  }
}

su2double COneShotSinglezoneDriver::ProjectionSet(unsigned short iDV, su2double value, bool active){
  if (active) {
      if(!activeset[iDV]) value = 0.0;
  } else {
      if(activeset[iDV]) value = 0.0;
  }
  return value;
}

su2double COneShotSinglezoneDriver::ProjectionPAP(unsigned short iDV, unsigned short jDV, su2double value, bool active){
  //returns for a Matrix entry a_iDV,jDV of a Matrix A the resulting entry of P*A*P (active or inactive set)
  if (active) {
      if(!activeset[iDV]||!activeset[jDV]) value = 0.0;
  } else {
      if(activeset[iDV]||activeset[jDV]) value = 0.0;
  }
  return value;
}

su2double COneShotSinglezoneDriver::BoundProjection(su2double value){
  if(value<=lb) value = lb;
  if(value>=ub) value = ub;
  return value;
}

void COneShotSinglezoneDriver::ComputeActiveSet(){
  //Compute ||x-P(x-gradx)||
  unsigned short iDV;
  su2double norm = 0.0;
  for (iDV=0; iDV<nDV_Total; iDV++) {
    norm+=(DesignVariable[iDV]-BoundProjection(DesignVariable[iDV]-config->GetStepSize()*ShiftedLagrangianGradient[iDV]))*(DesignVariable[iDV]-BoundProjection(DesignVariable[iDV]-config->GetStepSize()*ShiftedLagrangianGradient[iDV]));
  }
  norm=sqrt(norm);
  if(norm<(ub-lb)/2.0) {
    epsilon=norm;
  }
  for (iDV=0; iDV<nDV_Total; iDV++) {
    activeset[iDV]=false;
    if(ub-DesignVariable[iDV]<=epsilon) activeset[iDV]=true;
    if(DesignVariable[iDV]-lb<=epsilon) activeset[iDV]=true;
  }
  if (epsilon == 0) std::cout<<" All design variable bounds are reached "<<std::endl;
  else std::cout<<"Estimator for bounds: "<<epsilon<<std::endl;
}

void COneShotSinglezoneDriver::SetShiftedLagrangianGradient(){
  unsigned short iDV;
  for (iDV=0; iDV<nDV_Total; iDV++){
    ShiftedLagrangianGradient[iDV] = Gradient[iDV];
  }
}

void COneShotSinglezoneDriver::SetAugmentedLagrangianGradient(){
  unsigned short iDV;
  for (iDV=0; iDV<nDV_Total; iDV++){
    AugmentedLagrangianGradient[iDV] = Gradient[iDV];
  }
}

void COneShotSinglezoneDriver::ComputeAlphaTerm(){

    /*--- Note: Not applicable for unsteady code ---*/

    /*--- Initialize the adjoint of the output variables of the iteration with difference of the solution and the solution
     *    of the previous iteration. The values are passed to the AD tool. ---*/

    config->SetInnerIter(0);
    iteration->InitializeAdjoint_Update(solver_container, geometry_container, config_container, ZONE_0, INST_0);

    /*--- Initialize the adjoint of the objective function with 0.0. ---*/

    SetAdj_ObjFunction_Zero();
    SetAdj_ConstrFunction(ConstrFunc);

    /*--- Interpret the stored information by calling the corresponding routine of the AD tool. ---*/

    AD::ComputeAdjoint();

    /*--- Extract the computed adjoint values of the input variables and store them for the next iteration. ---*/
    iteration->Iterate_No_Residual(output_container[ZONE_0], integration_container, geometry_container,
                                            solver_container, numerics_container, config_container,
                                            surface_movement, grid_movement, FFDBox, ZONE_0, INST_0);

    /*--- Extract the computed sensitivity values. ---*/
    solver[ADJFLOW_SOL]->SetSensitivity(geometry,config);

    /*--- Clear the stored adjoint information to be ready for a new evaluation. ---*/

    AD::ClearAdjoints();

    AD::Reset();
}

void COneShotSinglezoneDriver::ComputePreconditioner(){

  unsigned short iConstr, jConstr;
  unsigned short nConstr = config->GetnConstr();

  su2double* seeding = new su2double[nConstr];
  for (iConstr = 0; iConstr < nConstr; iConstr++){
    seeding[iConstr] = 0.0;
  }
  su2double **BCheck = new su2double*[nConstr];
  for (iConstr = 0; iConstr  < nConstr; iConstr++){
    BCheck[iConstr] = new su2double[nConstr];
    for (jConstr = 0; jConstr  < nConstr; jConstr++){
      BCheck[iConstr][jConstr] = 0.0;
    }
  }

  for (unsigned short iConstr = 0; iConstr < nConstr; iConstr++){
    seeding[iConstr] = 1.0;

    config->SetInnerIter(0);
    iteration->InitializeAdjoint_Zero(solver_container, geometry_container, config_container, ZONE_0, INST_0);

    /*--- Initialize the adjoint of the objective function with 0.0. ---*/

    SetAdj_ObjFunction_Zero();
    SetAdj_ConstrFunction(seeding);

    /*--- Interpret the stored information by calling the corresponding routine of the AD tool. ---*/

    AD::ComputeAdjoint();

    /*--- Extract the computed adjoint values of the input variables and store them for the next iteration. ---*/
    iteration->Iterate_No_Residual(output_container[ZONE_0], integration_container, geometry_container,
                                          solver_container, numerics_container, config_container,
                                          surface_movement, grid_movement, FFDBox, ZONE_0, INST_0);

    solver[ADJFLOW_SOL]->SetConstrDerivative(iConstr);

    AD::ClearAdjoints();

    seeding[iConstr]=0.0;

  }

  su2double bcheck=0;
  for (iConstr = 0; iConstr  < nConstr; iConstr++){
    BCheck[iConstr][iConstr]=config->GetBCheckEpsilon();
    for (jConstr = 0; jConstr  < nConstr; jConstr++){
      BCheck[iConstr][jConstr]+=(1./config->GetMultiplierFactor(iConstr))*config->GetOneShotBeta()*solver[ADJFLOW_SOL]->MultiplyConstrDerivative(iConstr,jConstr);
    }
  }
  if (nConstr==1){
    BCheck_Inv[0][0]=1./BCheck[0][0];
  } else {
    bcheck=1./(BCheck[0][0]*BCheck[1][1]*BCheck[2][2]+BCheck[1][0]*BCheck[2][1]*BCheck[0][2]+BCheck[2][0]*BCheck[0][1]*BCheck[1][2]-BCheck[0][0]*BCheck[2][1]*BCheck[1][2]-BCheck[2][0]*BCheck[1][1]*BCheck[0][2]-BCheck[1][0]*BCheck[0][1]*BCheck[2][2]);
    BCheck_Inv[0][0]=bcheck*(BCheck[1][1]*BCheck[2][2]-BCheck[1][2]*BCheck[2][1]);
    BCheck_Inv[0][1]=bcheck*(BCheck[0][2]*BCheck[2][1]-BCheck[0][1]*BCheck[2][2]);
    BCheck_Inv[0][2]=bcheck*(BCheck[0][1]*BCheck[1][2]-BCheck[0][2]*BCheck[1][1]);
    BCheck_Inv[1][0]=bcheck*(BCheck[1][2]*BCheck[2][0]-BCheck[1][0]*BCheck[2][2]);
    BCheck_Inv[1][1]=bcheck*(BCheck[0][0]*BCheck[2][2]-BCheck[0][2]*BCheck[2][0]);
    BCheck_Inv[1][2]=bcheck*(BCheck[0][2]*BCheck[1][0]-BCheck[0][0]*BCheck[1][2]);
    BCheck_Inv[2][0]=bcheck*(BCheck[1][0]*BCheck[2][1]-BCheck[1][1]*BCheck[2][0]);
    BCheck_Inv[2][1]=bcheck*(BCheck[0][1]*BCheck[2][0]-BCheck[0][0]*BCheck[2][1]);
    BCheck_Inv[2][2]=bcheck*(BCheck[0][0]*BCheck[1][1]-BCheck[0][1]*BCheck[1][0]);
  }

  for (unsigned short iConstr = 0; iConstr  < config->GetnConstr(); iConstr++){
    delete [] BCheck[iConstr];
  }
  delete [] BCheck;
  delete [] seeding;
}

void COneShotSinglezoneDriver::ComputeBetaTerm(){

    /*--- Note: Not applicable for unsteady code ---*/

    /*--- For the one shot iteration we have to record for every steady state iteration. ---*/

    /*--- Store the computational graph of one direct iteration with the conservative variables and the mesh coordinates as input. ---*/

    SetRecording(RECORDING::SOLUTION_AND_MESH);

      /*--- Initialize the adjoint of the output variables of the iteration with the adjoint solution
     *    of the previous iteration. The values are passed to the AD tool. ---*/

    config->SetInnerIter(0);
    iteration->InitializeAdjoint(solver_container, geometry_container, config_container, ZONE_0, INST_0);

    /*--- Initialize the adjoint of the objective function with 1.0. ---*/

    SetAdj_ObjFunction();
    SetAdj_ConstrFunction(multiplier);

    /*--- Interpret the stored information by calling the corresponding routine of the AD tool. ---*/

    AD::ComputeAdjoint();

    /*--- Extract the computed adjoint values of the input variables and store them for the next iteration. ---*/
    iteration->Iterate_No_Residual(output_container[ZONE_0], integration_container, geometry_container,
                                            solver_container, numerics_container, config_container,
                                            surface_movement, grid_movement, FFDBox, ZONE_0, INST_0);

    /*--- Extract the computed sensitivity values. ---*/
    solver[ADJFLOW_SOL]->SetSensitivity(geometry,config);

    /*--- Clear the stored adjoint information to be ready for a new evaluation. ---*/

    AD::ClearAdjoints();

    AD::Reset();
}

void COneShotSinglezoneDriver::SetAdj_ObjFunction_Zero(){
  SU2_TYPE::SetDerivative(ObjFunc, 0.0);
}

void COneShotSinglezoneDriver::ProjectMeshSensitivities(){
  config->SetKind_SU2(SU2_COMPONENT::SU2_DOT); // set SU2_DOT as solver
  // get the dependency of the volumetric grid movement
  grid_movement[ZONE_0][INST_0]->SetVolume_Deformation(geometry, config, false, true);
  
  surface_movement[ZONE_0]->CopyBoundary(geometry, config);
  // project sensitivities (surface) on design variables
  if (config->GetProjectionAD()){
    SetProjection_AD(surface_movement[ZONE_0] , Gradient);
  }else{
    SetProjection_FD(surface_movement[ZONE_0] , Gradient);
  }

  config->SetKind_SU2(SU2_COMPONENT::SU2_CFD); // set SU2_CFD as solver
}

void COneShotSinglezoneDriver::SetAdj_ConstrFunction(su2double *seeding){

  for (unsigned short iConstr = 0; iConstr < config->GetnConstr(); iConstr++){
    if (rank == MASTER_NODE){
      SU2_TYPE::SetDerivative(ConstrFunc[iConstr], SU2_TYPE::GetValue(seeding[iConstr]));
    } else {
      SU2_TYPE::SetDerivative(ConstrFunc[iConstr], 0.0);
    }
  }

}

void COneShotSinglezoneDriver::SetConstrFunction(){

  bool heat         = (config->GetWeakly_Coupled_Heat());
  bool turbo        = (config->GetBoolTurbomachinery());

  unsigned short Kind_ConstrFunc;
  su2double FunctionValue = 0.0;

  for (unsigned short iConstr = 0; iConstr < config->GetnConstr(); iConstr++){
    if(iConstr==0) std::cout<<"Function Value: ";
    ConstrFunc[iConstr] = 0.0;

    /*--- Specific scalar constraint functions ---*/

    switch (config->GetKind_Solver()) {
    case DISC_ADJ_INC_EULER:       case DISC_ADJ_INC_NAVIER_STOKES:      case DISC_ADJ_INC_RANS:
    case DISC_ADJ_EULER:           case DISC_ADJ_NAVIER_STOKES:          case DISC_ADJ_RANS:
    case DISC_ADJ_FEM_EULER:       case DISC_ADJ_FEM_NS:                 case DISC_ADJ_FEM_RANS:
    case ONE_SHOT_EULER: 		 case ONE_SHOT_NAVIER_STOKES:	       case ONE_SHOT_RANS:

      solver[FLOW_SOL]->SetConstraintFunction(0.0, iConstr);

      /*--- Surface based obj. function ---*/

      Kind_ConstrFunc = config->GetKind_ConstrFunc(iConstr);
      if( Kind_ConstrFunc == AIRFOIL_AREA ) {
        unsigned short iPlane = 0;
        su2double *Plane_P0 = new su2double[3];
        su2double *Plane_Normal = new su2double[3];
        vector<su2double> *Xcoord_Airfoil = new vector<su2double>[1];
        vector<su2double> *Ycoord_Airfoil = new vector<su2double>[1];
        vector<su2double> *Zcoord_Airfoil = new vector<su2double>[1];
        vector<su2double> *Variable_Airfoil = new vector<su2double>[1];

        Plane_Normal[0] = 0.0;   Plane_P0[0] = 0.0;
        Plane_Normal[1] = 1.0;   Plane_P0[1] = 0.0;
        Plane_Normal[2] = 0.0;   Plane_P0[2] = 0.0;

        geometry->ComputeAirfoil_Section(Plane_P0, Plane_Normal, -1E6, 1E6, -1E6, 1E6, -1E6, 1E6, NULL,
                                                           Xcoord_Airfoil[iPlane], Ycoord_Airfoil[iPlane], Zcoord_Airfoil[iPlane],
                                                           Variable_Airfoil[iPlane], true, config);
        FunctionValue = geometry->Compute_Area(Plane_P0, Plane_Normal, config, Xcoord_Airfoil[iPlane], Ycoord_Airfoil[iPlane], Zcoord_Airfoil[iPlane]);

        delete [] Plane_P0;
        delete [] Plane_Normal;
        delete [] Xcoord_Airfoil;
        delete [] Ycoord_Airfoil;
        delete [] Zcoord_Airfoil;
        delete [] Variable_Airfoil;
      }
      if( Kind_ConstrFunc == MAX_THICKNESS ) {
        unsigned short iPlane = 0;
        su2double *Plane_P0 = new su2double[3];
        su2double *Plane_Normal = new su2double[3];
        vector<su2double> *Xcoord_Airfoil = new vector<su2double>[1];
        vector<su2double> *Ycoord_Airfoil = new vector<su2double>[1];
        vector<su2double> *Zcoord_Airfoil = new vector<su2double>[1];
        vector<su2double> *Variable_Airfoil = new vector<su2double>[1];

        Plane_Normal[0] = 0.0;   Plane_P0[0] = 0.0;
        Plane_Normal[1] = 1.0;   Plane_P0[1] = 0.0;
        Plane_Normal[2] = 0.0;   Plane_P0[2] = 0.0;

        geometry->ComputeAirfoil_Section(Plane_P0, Plane_Normal, -1E6, 1E6, -1E6, 1E6, -1E6, 1E6, NULL,
                                                           Xcoord_Airfoil[iPlane], Ycoord_Airfoil[iPlane], Zcoord_Airfoil[iPlane],
                                                           Variable_Airfoil[iPlane], true, config);
        FunctionValue = geometry->Compute_MaxThickness(Plane_P0, Plane_Normal, config, Xcoord_Airfoil[iPlane], Ycoord_Airfoil[iPlane], Zcoord_Airfoil[iPlane]);

        delete [] Plane_P0;
        delete [] Plane_Normal;
        delete [] Xcoord_Airfoil;
        delete [] Ycoord_Airfoil;
        delete [] Zcoord_Airfoil;
        delete [] Variable_Airfoil;

      }else{
        FunctionValue = solver[FLOW_SOL]->Evaluate_ConstrFunc(config, iConstr);
      }
    
      if (heat){
        if (config->GetKind_ObjFunc() == TOTAL_HEATFLUX) {
          FunctionValue += solver[HEAT_SOL]->GetTotal_HeatFlux();
        }
        else if (config->GetKind_ObjFunc() == AVG_TEMPERATURE) {
          FunctionValue += solver[HEAT_SOL]->GetTotal_AvgTemperature();
        }
      }
    }

    ConstrFunc[iConstr] = config->GetConstraintScale(iConstr)*(config->GetConstraintTarget(iConstr) - FunctionValue);
    //TODO-ONE-SHOT or use solver[FLOW_SOL]->GetConstraintFunction(iConstr);
    //TODO-ONE-SHOT turbo

    std::cout<<FunctionValue<<" ";
    if (rank == MASTER_NODE){
      AD::RegisterOutput(ConstrFunc[iConstr]);
    }
    if(iConstr==0) std::cout<<std::endl;
  }
}

void COneShotSinglezoneDriver::UpdateMultiplier(su2double stepsize){
  su2double helper;
  for(unsigned short iConstr = 0; iConstr < config->GetnConstr(); iConstr++){
    if(iConstr==0) std::cout<<"Multiplier Update: ";
    helper = 0.0;
    for(unsigned short jConstr = 0; jConstr < config->GetnConstr(); jConstr++){
       helper+= BCheck_Inv[iConstr][jConstr]*ConstrFunc_Store[jConstr];
    }
    multiplier[iConstr] = multiplier[iConstr]+stepsize*helper*config->GetMultiplierScale(iConstr);
    std::cout<<ConstrFunc_Store[iConstr]<<" "<<multiplier[iConstr]<<" ";
    if(iConstr==config->GetnConstr()-1) std::cout<<std::endl;
  }
}

void COneShotSinglezoneDriver::StoreConstrFunction(){
  for(unsigned short iConstr = 0; iConstr < config->GetnConstr(); iConstr++){
    ConstrFunc_Store[iConstr] = ConstrFunc[iConstr];
  }
}
