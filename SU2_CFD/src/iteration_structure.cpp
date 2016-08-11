/*!
 * \file iteration_structure.cpp
 * \brief Main subroutines used by SU2_CFD
 * \author F. Palacios, T. Economon
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

#include "../include/iteration_structure.hpp"

#include <cstdlib>
#include <stdio.h>

#include "../../SparseGrid/SG/oIndex.h"
//#include "../../SparseGrid/SG/oSolverControl.h"
#include "../../SparseGrid/SG/oFilesystem.h"
#include "../../SparseGrid/SG/util.h"
#include "../../SparseGrid/SG/interfaceOBJ.h"
#include "../../SparseGrid/SG/sparsegrid.h"

CIteration::CIteration(CConfig *config) { }
CIteration::~CIteration(void) { }

void CIteration::Preprocess(COutput *output,
                            CIntegration ***integration_container,
                            CGeometry ***geometry_container,
                            CSolver ****solver_container,
                            CNumerics *****numerics_container,
                            CConfig **config_container,
                            CSurfaceMovement **surface_movement,
                            CVolumetricMovement **grid_movement,
                            CFreeFormDefBox*** FFDBox,
                            unsigned short val_iZone) { }
void CIteration::Iterate(COutput *output,
                         CIntegration ***integration_container,
                         CGeometry ***geometry_container,
                         CSolver ****solver_container,
                         CNumerics *****numerics_container,
                         CConfig **config_container,
                         CSurfaceMovement **surface_movement,
                         CVolumetricMovement **grid_movement,
                         CFreeFormDefBox*** FFDBox,
                         unsigned short val_iZone) { }
void CIteration::Update(COutput *output,
                        CIntegration ***integration_container,
                        CGeometry ***geometry_container,
                        CSolver ****solver_container,
                        CNumerics *****numerics_container,
                        CConfig **config_container,
                        CSurfaceMovement **surface_movement,
                        CVolumetricMovement **grid_movement,
                        CFreeFormDefBox*** FFDBox,
                        unsigned short val_iZone)      { }
void CIteration::Monitor()     { }
void CIteration::Output()      { }
void CIteration::Postprocess() { }



CMeanFlowIteration::CMeanFlowIteration(CConfig *config) : CIteration(config) { }
CMeanFlowIteration::~CMeanFlowIteration(void) { }

void CMeanFlowIteration::Preprocess(COutput *output,
                                    CIntegration ***integration_container,
                                    CGeometry ***geometry_container,
                                    CSolver ****solver_container,
                                    CNumerics *****numerics_container,
                                    CConfig **config_container,
                                    CSurfaceMovement **surface_movement,
                                    CVolumetricMovement **grid_movement,
                                    CFreeFormDefBox*** FFDBox,
                                    unsigned short val_iZone) {
  
  unsigned long IntIter = 0; config_container[val_iZone]->SetIntIter(IntIter);
  unsigned long ExtIter = config_container[val_iZone]->GetExtIter();
  
  bool fsi = config_container[val_iZone]->GetFSI_Simulation();
  unsigned long FSIIter = config_container[val_iZone]->GetFSIIter();

  bool time_spectral = (config_container[val_iZone]->GetUnsteady_Simulation() == TIME_SPECTRAL);
  
  /*--- Set the initial condition ---*/
  /*--- For FSI problems with subiterations, this must only be done in the first subiteration ---*/
  if(!( (fsi) && (FSIIter > 0) ))
	 solver_container[val_iZone][MESH_0][FLOW_SOL]->SetInitialCondition(geometry_container[val_iZone], solver_container[val_iZone], config_container[val_iZone], ExtIter);
  
  /*--- Dynamic mesh update ---*/
  
  if ((config_container[val_iZone]->GetGrid_Movement()) && (!time_spectral)) {
    SetGrid_Movement(geometry_container[val_iZone], surface_movement[val_iZone], grid_movement[val_iZone], FFDBox[val_iZone], solver_container[val_iZone], config_container[val_iZone], val_iZone, IntIter, ExtIter);
  }
  
  /*--- Apply a Wind Gust ---*/
  
  if (config_container[val_iZone]->GetWind_Gust()) {
    SetWind_GustField(config_container[val_iZone], geometry_container[val_iZone], solver_container[val_iZone]);
  }
  
  
  /*--- Calculate and set Mixing Plane averaged quantities at interfaces ---*/
  
  if(config_container[val_iZone]->GetBoolMixingPlane())
    SetMixingPlane(geometry_container, solver_container, config_container, val_iZone);
  
  /*--- Compute turboperformance ---*/
  
  if(config_container[val_iZone]->GetBoolTurboPerf())
    SetTurboPerformance(geometry_container, solver_container, config_container, output, val_iZone);
}

void CMeanFlowIteration::Iterate(COutput *output,
                                 CIntegration ***integration_container,
                                 CGeometry ***geometry_container,
                                 CSolver ****solver_container,
                                 CNumerics *****numerics_container,
                                 CConfig **config_container,
                                 CSurfaceMovement **surface_movement,
                                 CVolumetricMovement **grid_movement,
                                 CFreeFormDefBox*** FFDBox,
                                 unsigned short val_iZone) {
  
  unsigned long IntIter = 0; config_container[val_iZone]->SetIntIter(IntIter);
  unsigned long ExtIter = config_container[val_iZone]->GetExtIter();
#ifdef HAVE_MPI
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  
  /*--- Set the value of the internal iteration ---*/
  
  IntIter = ExtIter;
  if ((config_container[val_iZone]->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
      (config_container[val_iZone]->GetUnsteady_Simulation() == DT_STEPPING_2ND)) IntIter = 0;
  
  /*--- Update global parameters ---*/
  
  if ((config_container[val_iZone]->GetKind_Solver() == EULER) ||
      (config_container[val_iZone]->GetKind_Solver() == DISC_ADJ_EULER)) {
    config_container[val_iZone]->SetGlobalParam(EULER, RUNTIME_FLOW_SYS, ExtIter);
  }
  if ((config_container[val_iZone]->GetKind_Solver() == NAVIER_STOKES) ||
      (config_container[val_iZone]->GetKind_Solver() == DISC_ADJ_NAVIER_STOKES)) {
    config_container[val_iZone]->SetGlobalParam(NAVIER_STOKES, RUNTIME_FLOW_SYS, ExtIter);
  }
  if ((config_container[val_iZone]->GetKind_Solver() == RANS) ||
      (config_container[val_iZone]->GetKind_Solver() == DISC_ADJ_RANS)) {
    config_container[val_iZone]->SetGlobalParam(RANS, RUNTIME_FLOW_SYS, ExtIter);
  }
  
  /*--- Solve the Euler, Navier-Stokes or Reynolds-averaged Navier-Stokes (RANS) equations (one iteration) ---*/
  
  integration_container[val_iZone][FLOW_SOL]->MultiGrid_Iteration(geometry_container, solver_container, numerics_container,
                                                                  config_container, RUNTIME_FLOW_SYS, IntIter, val_iZone);
  
  if ((config_container[val_iZone]->GetKind_Solver() == RANS) ||
      (config_container[val_iZone]->GetKind_Solver() == DISC_ADJ_RANS)) {
    
    /*--- Solve the turbulence model ---*/
    
    config_container[val_iZone]->SetGlobalParam(RANS, RUNTIME_TURB_SYS, ExtIter);
    integration_container[val_iZone][TURB_SOL]->SingleGrid_Iteration(geometry_container, solver_container, numerics_container,
                                                                     config_container, RUNTIME_TURB_SYS, IntIter, val_iZone);
    
    /*--- Solve transition model ---*/
    
    if (config_container[val_iZone]->GetKind_Trans_Model() == LM) {
      config_container[val_iZone]->SetGlobalParam(RANS, RUNTIME_TRANS_SYS, ExtIter);
      integration_container[val_iZone][TRANS_SOL]->SingleGrid_Iteration(geometry_container, solver_container, numerics_container,
                                                                        config_container, RUNTIME_TRANS_SYS, IntIter, val_iZone);
    }
    
  }
  
  /*--- Dual time stepping strategy ---*/

   if (((config_container[val_iZone]->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
        (config_container[val_iZone]->GetUnsteady_Simulation() == DT_STEPPING_2ND)) &&
	 !config_container[val_iZone]->GetDiscrete_Adjoint()) {

    for (IntIter = 1; IntIter < config_container[val_iZone]->GetUnst_nIntIter(); IntIter++) {
      
      /*--- Write the convergence history (only screen output) ---*/
      
      output->SetConvHistory_Body(NULL, geometry_container, solver_container, config_container, integration_container, true, 0.0, val_iZone);
      
      /*--- Set the value of the internal iteration ---*/
      
      config_container[val_iZone]->SetIntIter(IntIter);
      
      /*--- Pseudo-timestepping for the Euler, Navier-Stokes or Reynolds-averaged Navier-Stokes equations ---*/
      
      if ((config_container[val_iZone]->GetKind_Solver() == EULER) ||
          (config_container[val_iZone]->GetKind_Solver() == DISC_ADJ_EULER)) {
        config_container[val_iZone]->SetGlobalParam(EULER, RUNTIME_FLOW_SYS, ExtIter);
      }
      if ((config_container[val_iZone]->GetKind_Solver() == NAVIER_STOKES) ||
          (config_container[val_iZone]->GetKind_Solver() == DISC_ADJ_NAVIER_STOKES)) {
        config_container[val_iZone]->SetGlobalParam(NAVIER_STOKES, RUNTIME_FLOW_SYS, ExtIter);
      }
      if ((config_container[val_iZone]->GetKind_Solver() == RANS) ||
          (config_container[val_iZone]->GetKind_Solver() == DISC_ADJ_RANS)) {
        config_container[val_iZone]->SetGlobalParam(RANS, RUNTIME_FLOW_SYS, ExtIter);
      }
      
      /*--- Solve the Euler, Navier-Stokes or Reynolds-averaged Navier-Stokes (RANS) equations (one iteration) ---*/
      
      integration_container[val_iZone][FLOW_SOL]->MultiGrid_Iteration(geometry_container, solver_container, numerics_container,
                                                                      config_container, RUNTIME_FLOW_SYS, IntIter, val_iZone);
      
      /*--- Pseudo-timestepping the turbulence model ---*/
      
      if ((config_container[val_iZone]->GetKind_Solver() == RANS) ||
          (config_container[val_iZone]->GetKind_Solver() == DISC_ADJ_RANS)) {
        
        /*--- Solve the turbulence model ---*/
        
        config_container[val_iZone]->SetGlobalParam(RANS, RUNTIME_TURB_SYS, ExtIter);
        integration_container[val_iZone][TURB_SOL]->SingleGrid_Iteration(geometry_container, solver_container, numerics_container,
                                                                         config_container, RUNTIME_TURB_SYS, IntIter, val_iZone);
        
        /*--- Solve transition model ---*/
        
        if (config_container[val_iZone]->GetKind_Trans_Model() == LM) {
          config_container[val_iZone]->SetGlobalParam(RANS, RUNTIME_TRANS_SYS, ExtIter);
          integration_container[val_iZone][TRANS_SOL]->SingleGrid_Iteration(geometry_container, solver_container, numerics_container,
                                                                            config_container, RUNTIME_TRANS_SYS, IntIter, val_iZone);
        }
        
      }
      
      /*--- Call Dynamic mesh update if AEROELASTIC motion was specified ---*/
      if ((config_container[val_iZone]->GetGrid_Movement()) && (config_container[val_iZone]->GetAeroelastic_Simulation())) {
        SetGrid_Movement(geometry_container[val_iZone], surface_movement[val_iZone], grid_movement[val_iZone], FFDBox[val_iZone],
                         solver_container[val_iZone], config_container[val_iZone], val_iZone, IntIter, ExtIter);
        /*--- Apply a Wind Gust ---*/
        if (config_container[val_iZone]->GetWind_Gust()) {
          if (IntIter % config_container[val_iZone]->GetAeroelasticIter() ==0)
            SetWind_GustField(config_container[val_iZone], geometry_container[val_iZone], solver_container[val_iZone]);
        }
      }
      
      if (integration_container[val_iZone][FLOW_SOL]->GetConvergence()) break;
      
    }
    
  }
  
}

void CMeanFlowIteration::Update(COutput *output,
                                CIntegration ***integration_container,
                                CGeometry ***geometry_container,
                                CSolver ****solver_container,
                                CNumerics *****numerics_container,
                                CConfig **config_container,
                                CSurfaceMovement **surface_movement,
                                CVolumetricMovement **grid_movement,
                                CFreeFormDefBox*** FFDBox,
                                unsigned short val_iZone)      {
  
  unsigned short iMesh;
  su2double Physical_dt, Physical_t;
  unsigned long ExtIter = config_container[val_iZone]->GetExtIter();

  /*--- Dual time stepping strategy ---*/
  
  if ((config_container[val_iZone]->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
      (config_container[val_iZone]->GetUnsteady_Simulation() == DT_STEPPING_2ND)) {
    
    /*--- Update dual time solver on all mesh levels ---*/
    
    for (iMesh = 0; iMesh <= config_container[val_iZone]->GetnMGLevels(); iMesh++) {
      integration_container[val_iZone][FLOW_SOL]->SetDualTime_Solver(geometry_container[val_iZone][iMesh], solver_container[val_iZone][iMesh][FLOW_SOL], config_container[val_iZone], iMesh);
      integration_container[val_iZone][FLOW_SOL]->SetConvergence(false);
    }
    
    /*--- Update dual time solver for the turbulence model ---*/
    
    if ((config_container[val_iZone]->GetKind_Solver() == RANS) ||
        (config_container[val_iZone]->GetKind_Solver() == DISC_ADJ_RANS)) {
      integration_container[val_iZone][TURB_SOL]->SetDualTime_Solver(geometry_container[val_iZone][MESH_0], solver_container[val_iZone][MESH_0][TURB_SOL], config_container[val_iZone], MESH_0);
      integration_container[val_iZone][TURB_SOL]->SetConvergence(false);
    }
    
    /*--- Update dual time solver for the transition model ---*/
    
    if (config_container[val_iZone]->GetKind_Trans_Model() == LM) {
      integration_container[val_iZone][TRANS_SOL]->SetDualTime_Solver(geometry_container[val_iZone][MESH_0], solver_container[val_iZone][MESH_0][TRANS_SOL], config_container[val_iZone], MESH_0);
      integration_container[val_iZone][TRANS_SOL]->SetConvergence(false);
    }
    
    /*--- Verify convergence criteria (based on total time) ---*/
    
    Physical_dt = config_container[val_iZone]->GetDelta_UnstTime();
    Physical_t  = (ExtIter+1)*Physical_dt;
    if (Physical_t >=  config_container[val_iZone]->GetTotal_UnstTime())
      integration_container[val_iZone][FLOW_SOL]->SetConvergence(true);
    
  }
  
}

void CMeanFlowIteration::Monitor()     { }
void CMeanFlowIteration::Output()      { }
void CMeanFlowIteration::Postprocess() { }

void CMeanFlowIteration::SetWind_GustField(CConfig *config_container, CGeometry **geometry_container, CSolver ***solver_container) {
  // The gust is imposed on the flow field via the grid velocities. This method called the Field Velocity Method is described in the
  // NASA TM–2012-217771 - Development, Verification and Use of Gust Modeling in the NASA Computational Fluid Dynamics Code FUN3D
  // the desired gust is prescribed as the negative of the grid velocity.
  
  // If a source term is included to account for the gust field, the method is described by Jones et al. as the Split Velocity Method in
  // Simulation of Airfoil Gust Responses Using Prescribed Velocities.
  // In this routine the gust derivatives needed for the source term are calculated when applicable.
  // If the gust derivatives are zero the source term is also zero.
  // The source term itself is implemented in the class CSourceWindGust
  
  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  
  if (rank == MASTER_NODE)
    cout << endl << "Running simulation with a Wind Gust." << endl;
  unsigned short iDim, nDim = geometry_container[MESH_0]->GetnDim(); //We assume nDim = 2
  if (nDim != 2) {
    if (rank == MASTER_NODE) {
      cout << endl << "WARNING - Wind Gust capability is only verified for 2 dimensional simulations." << endl;
    }
  }
  
  /*--- Gust Parameters from config ---*/
  unsigned short Gust_Type = config_container->GetGust_Type();
  su2double xbegin = config_container->GetGust_Begin_Loc();    // Location at which the gust begins.
  su2double L = config_container->GetGust_WaveLength();        // Gust size
  su2double tbegin = config_container->GetGust_Begin_Time();   // Physical time at which the gust begins.
  su2double gust_amp = config_container->GetGust_Ampl();       // Gust amplitude
  su2double n = config_container->GetGust_Periods();           // Number of gust periods
  unsigned short GustDir = config_container->GetGust_Dir(); // Gust direction
  
  /*--- Variables needed to compute the gust ---*/
  unsigned short Kind_Grid_Movement = config_container->GetKind_GridMovement(ZONE_0);
  unsigned long iPoint;
  unsigned short iMGlevel, nMGlevel = config_container->GetnMGLevels();
  
  su2double x, y, x_gust, dgust_dx, dgust_dy, dgust_dt;
  su2double *Gust, *GridVel, *NewGridVel, *GustDer;
  
  su2double Physical_dt = config_container->GetDelta_UnstTime();
  unsigned long ExtIter = config_container->GetExtIter();
  su2double Physical_t = ExtIter*Physical_dt;
  
  su2double Uinf = solver_container[MESH_0][FLOW_SOL]->GetVelocity_Inf(0); // Assumption gust moves at infinity velocity
  
  Gust = new su2double [nDim];
  NewGridVel = new su2double [nDim];
  for (iDim = 0; iDim < nDim; iDim++) {
    Gust[iDim] = 0.0;
    NewGridVel[iDim] = 0.0;
  }
  
  GustDer = new su2double [3];
  for (unsigned short i = 0; i < 3; i++) {
    GustDer[i] = 0.0;
  }
  
  // Vortex variables
  unsigned long nVortex = 0;
  vector<su2double> x0, y0, vort_strenth, r_core; //vortex is positive in clockwise direction.
  if (Gust_Type == VORTEX) {
    InitializeVortexDistribution(nVortex, x0, y0, vort_strenth, r_core);
  }
  
  /*--- Check to make sure gust lenght is not zero or negative (vortex gust doesn't use this). ---*/
  if (L <= 0.0 && Gust_Type != VORTEX) {
    if (rank == MASTER_NODE) cout << "ERROR: The gust length needs to be positive" << endl;
#ifndef HAVE_MPI
    exit(EXIT_FAILURE);
#else
    MPI_Abort(MPI_COMM_WORLD,1);
    MPI_Finalize();
#endif
  }
  
  /*--- Loop over all multigrid levels ---*/
  
  for (iMGlevel = 0; iMGlevel <= nMGlevel; iMGlevel++) {
    
    /*--- Loop over each node in the volume mesh ---*/
    
    for (iPoint = 0; iPoint < geometry_container[iMGlevel]->GetnPoint(); iPoint++) {
      
      /*--- Reset the Grid Velocity to zero if there is no grid movement ---*/
      if (Kind_Grid_Movement == GUST) {
        for (iDim = 0; iDim < nDim; iDim++)
          geometry_container[iMGlevel]->node[iPoint]->SetGridVel(iDim, 0.0);
      }
      
      /*--- initialize the gust and derivatives to zero everywhere ---*/
      
      for (iDim = 0; iDim < nDim; iDim++) {Gust[iDim]=0.0;}
      dgust_dx = 0.0; dgust_dy = 0.0; dgust_dt = 0.0;
      
      /*--- Begin applying the gust ---*/
      
      if (Physical_t >= tbegin) {
        
        x = geometry_container[iMGlevel]->node[iPoint]->GetCoord()[0]; // x-location of the node.
        y = geometry_container[iMGlevel]->node[iPoint]->GetCoord()[1]; // y-location of the node.
        
        // Gust coordinate
        x_gust = (x - xbegin - Uinf*(Physical_t-tbegin))/L;
        
        /*--- Calculate the specified gust ---*/
        switch (Gust_Type) {
            
          case TOP_HAT:
            // Check if we are in the region where the gust is active
            if (x_gust > 0 && x_gust < n) {
              Gust[GustDir] = gust_amp;
              // Still need to put the gust derivatives. Think about this.
            }
            break;
            
          case SINE:
            // Check if we are in the region where the gust is active
            if (x_gust > 0 && x_gust < n) {
              Gust[GustDir] = gust_amp*(sin(2*PI_NUMBER*x_gust));
              
              // Gust derivatives
              //dgust_dx = gust_amp*2*PI_NUMBER*(cos(2*PI_NUMBER*x_gust))/L;
              //dgust_dy = 0;
              //dgust_dt = gust_amp*2*PI_NUMBER*(cos(2*PI_NUMBER*x_gust))*(-Uinf)/L;
            }
            break;
            
          case ONE_M_COSINE:
            // Check if we are in the region where the gust is active
            if (x_gust > 0 && x_gust < n) {
              Gust[GustDir] = gust_amp*(1-cos(2*PI_NUMBER*x_gust));
              
              // Gust derivatives
              //dgust_dx = gust_amp*2*PI_NUMBER*(sin(2*PI_NUMBER*x_gust))/L;
              //dgust_dy = 0;
              //dgust_dt = gust_amp*2*PI_NUMBER*(sin(2*PI_NUMBER*x_gust))*(-Uinf)/L;
            }
            break;
            
          case EOG:
            // Check if we are in the region where the gust is active
            if (x_gust > 0 && x_gust < n) {
              Gust[GustDir] = -0.37*gust_amp*sin(3*PI_NUMBER*x_gust)*(1-cos(2*PI_NUMBER*x_gust));
            }
            break;
            
          case VORTEX:
            
            /*--- Use vortex distribution ---*/
            // Algebraic vortex equation.
            for (unsigned long i=0; i<nVortex; i++) {
              su2double r2 = pow(x-(x0[i]+Uinf*(Physical_t-tbegin)), 2) + pow(y-y0[i], 2);
              su2double r = sqrt(r2);
              su2double v_theta = vort_strenth[i]/(2*PI_NUMBER) * r/(r2+pow(r_core[i],2));
              Gust[0] = Gust[0] + v_theta*(y-y0[i])/r;
              Gust[1] = Gust[1] - v_theta*(x-(x0[i]+Uinf*(Physical_t-tbegin)))/r;
            }
            break;
            
          case NONE: default:
            
            /*--- There is no wind gust specified. ---*/
            if (rank == MASTER_NODE) {
              cout << "No wind gust specified." << endl;
            }
            break;
            
        }
      }
      
      /*--- Set the Wind Gust, Wind Gust Derivatives and the Grid Velocities ---*/
      
      GustDer[0] = dgust_dx;
      GustDer[1] = dgust_dy;
      GustDer[2] = dgust_dt;
      
      solver_container[iMGlevel][FLOW_SOL]->node[iPoint]->SetWindGust(Gust);
      solver_container[iMGlevel][FLOW_SOL]->node[iPoint]->SetWindGustDer(GustDer);
      
      GridVel = geometry_container[iMGlevel]->node[iPoint]->GetGridVel();
      
      /*--- Store new grid velocity ---*/
      
      for (iDim = 0; iDim < nDim; iDim++) {
        NewGridVel[iDim] = GridVel[iDim] - Gust[iDim];
        geometry_container[iMGlevel]->node[iPoint]->SetGridVel(iDim, NewGridVel[iDim]);
      }
      
    }
  }
  
  delete [] Gust;
  delete [] GustDer;
  delete [] NewGridVel;
  
}

void CMeanFlowIteration::InitializeVortexDistribution(unsigned long &nVortex, vector<su2double>& x0, vector<su2double>& y0, vector<su2double>& vort_strength, vector<su2double>& r_core) {
  /*--- Read in Vortex Distribution ---*/
  std::string line;
  std::ifstream file;
  su2double x_temp, y_temp, vort_strength_temp, r_core_temp;
  file.open("vortex_distribution.txt");
  /*--- In case there is no vortex file ---*/
  if (file.fail()) {
    cout << "There is no vortex data file!!" << endl;
    cout << "Press any key to exit..." << endl;
    cin.get(); exit(EXIT_FAILURE);
  }
  
  // Ignore line containing the header
  getline(file, line);
  // Read in the information of the vortices (xloc, yloc, lambda(strength), eta(size, gradient))
  while (file.good())
  {
    getline(file, line);
    std::stringstream ss(line);
    if (line.size() != 0) { //ignore blank lines if they exist.
      ss >> x_temp;
      ss >> y_temp;
      ss >> vort_strength_temp;
      ss >> r_core_temp;
      x0.push_back(x_temp);
      y0.push_back(y_temp);
      vort_strength.push_back(vort_strength_temp);
      r_core.push_back(r_core_temp);
    }
  }
  file.close();
  // number of vortices
  nVortex = x0.size();
  
}

void CMeanFlowIteration::SetMixingPlane(CGeometry ***geometry_container, CSolver ****solver_container, CConfig **config_container, unsigned short iZone) {
  
  unsigned short jZone;
  unsigned short nZone = geometry_container[ZONE_0][MESH_0]->GetnZone();
  int intMarker, extMarker, intMarkerMix;
  string intMarker_Tag, extMarker_Tag;
  
  /*-- Loop on all the boundary to find MIXING_PLANE boundary --*/
  for (intMarker = 0; intMarker < config_container[iZone]->GetnMarker_All(); intMarker++) {
    for (intMarkerMix=0; intMarkerMix < config_container[iZone]->Get_nMarkerMixingPlane(); intMarkerMix++)
      if (config_container[iZone]->GetMarker_All_TagBound(intMarker) == config_container[iZone]->GetMarker_MixingPlane_Bound(intMarkerMix) ) {
        solver_container[iZone][MESH_0][FLOW_SOL]->Mixing_Process(geometry_container[iZone][MESH_0], solver_container[iZone][MESH_0], config_container[iZone], intMarker);
        extMarker_Tag = config_container[iZone]->GetMarker_MixingPlane_Donor(intMarkerMix);
        for (jZone = 0; jZone < nZone; jZone++){
          for (extMarker = 0; extMarker < config_container[jZone]->GetnMarker_All(); extMarker++)
            if (config_container[jZone]->GetMarker_All_TagBound(extMarker) == extMarker_Tag){
              solver_container[jZone][MESH_0][FLOW_SOL]->SetExtAveragedValue(solver_container[iZone][MESH_0][FLOW_SOL], intMarker, extMarker);
            }
        }
        
      }
    
  }
  
}

void CMeanFlowIteration::SetTurboPerformance(CGeometry ***geometry_container, CSolver ****solver_container, CConfig **config_container, COutput *output, unsigned short iZone) {
  
  unsigned short  jZone, inMarker, outMarker, inMarkerTP, Kind_TurboPerf;
  unsigned short nZone = geometry_container[iZone][MESH_0]->GetnZone();
  string inMarker_Tag, outMarker_Tag;
  
  
  /*-- Loop on all the boundary to find MIXING_PLANE boundary --*/
  for (inMarker = 0; inMarker < config_container[iZone]->GetnMarker_All(); inMarker++)
    for (inMarkerTP=0; inMarkerTP < config_container[iZone]->Get_nMarkerTurboPerf(); inMarkerTP++)
      if (config_container[iZone]->GetMarker_All_TagBound(inMarker) == config_container[iZone]->GetMarker_TurboPerf_BoundIn(inMarkerTP) ) {
        outMarker_Tag =	config_container[iZone]->GetMarker_TurboPerf_BoundOut(inMarkerTP);
        Kind_TurboPerf = config_container[iZone]->GetKind_TurboPerf(inMarkerTP);
        for (jZone = 0; jZone < nZone; jZone++)
          for (outMarker = 0; outMarker < config_container[jZone]->GetnMarker_All(); outMarker++)
            if (config_container[jZone]->GetMarker_All_TagBound(outMarker) == outMarker_Tag){
              solver_container[iZone][MESH_0][FLOW_SOL]->Mixing_Process(geometry_container[iZone][MESH_0], solver_container[iZone][MESH_0], config_container[iZone], inMarker);
              solver_container[jZone][MESH_0][FLOW_SOL]->Mixing_Process(geometry_container[jZone][MESH_0], solver_container[jZone][MESH_0], config_container[jZone], outMarker);
              solver_container[iZone][MESH_0][FLOW_SOL]->TurboPerformance(solver_container[jZone][MESH_0][FLOW_SOL], config_container[iZone], inMarker, outMarker, Kind_TurboPerf, inMarkerTP);
              solver_container[ZONE_0][MESH_0][FLOW_SOL]->StoreTurboPerformance(solver_container[iZone][MESH_0][FLOW_SOL], inMarkerTP);
            }
      }
}


CWaveIteration::CWaveIteration(CConfig *config) : CIteration(config) { }
CWaveIteration::~CWaveIteration(void) { }
void CWaveIteration::Preprocess(COutput *output,
                                CIntegration ***integration_container,
                                CGeometry ***geometry_container,
                                CSolver ****solver_container,
                                CNumerics *****numerics_container,
                                CConfig **config_container,
                                CSurfaceMovement **surface_movement,
                                CVolumetricMovement **grid_movement,
                                CFreeFormDefBox*** FFDBox,
                                unsigned short val_iZone) { }
void CWaveIteration::Iterate(COutput *output,
                             CIntegration ***integration_container,
                             CGeometry ***geometry_container,
                             CSolver ****solver_container,
                             CNumerics *****numerics_container,
                             CConfig **config_container,
                             CSurfaceMovement **surface_movement,
                             CVolumetricMovement **grid_movement,
                             CFreeFormDefBox*** FFDBox,
                             unsigned short val_iZone) {
  
  unsigned long IntIter = 0; config_container[ZONE_0]->SetIntIter(IntIter);
  unsigned long ExtIter = config_container[ZONE_0]->GetExtIter();
  
  /*--- Set the value of the internal iteration ---*/
  IntIter = ExtIter;
  if ((config_container[val_iZone]->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
      (config_container[val_iZone]->GetUnsteady_Simulation() == DT_STEPPING_2ND)) IntIter = 0;
  
  /*--- Wave equations ---*/
  config_container[val_iZone]->SetGlobalParam(WAVE_EQUATION, RUNTIME_WAVE_SYS, ExtIter);
  integration_container[val_iZone][WAVE_SOL]->SingleGrid_Iteration(geometry_container, solver_container, numerics_container,
                                                                   config_container, RUNTIME_WAVE_SYS, IntIter, val_iZone);
  
  /*--- Dual time stepping strategy ---*/
  if ((config_container[val_iZone]->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
      (config_container[val_iZone]->GetUnsteady_Simulation() == DT_STEPPING_2ND)) {
    
    for (IntIter = 1; IntIter < config_container[val_iZone]->GetUnst_nIntIter(); IntIter++) {
      output->SetConvHistory_Body(NULL, geometry_container, solver_container, config_container, integration_container, true, 0.0, val_iZone);
      config_container[val_iZone]->SetIntIter(IntIter);
      integration_container[val_iZone][WAVE_SOL]->SingleGrid_Iteration(geometry_container, solver_container, numerics_container,
                                                                       config_container, RUNTIME_WAVE_SYS, IntIter, val_iZone);
      if (integration_container[val_iZone][WAVE_SOL]->GetConvergence()) break;
    }
    
  }
  
}
void CWaveIteration::Update(COutput *output,
                            CIntegration ***integration_container,
                            CGeometry ***geometry_container,
                            CSolver ****solver_container,
                            CNumerics *****numerics_container,
                            CConfig **config_container,
                            CSurfaceMovement **surface_movement,
                            CVolumetricMovement **grid_movement,
                            CFreeFormDefBox*** FFDBox,
                            unsigned short val_iZone)      {
  
  unsigned short iMesh;
  su2double Physical_dt, Physical_t;
  unsigned long ExtIter = config_container[ZONE_0]->GetExtIter();
  
  /*--- Dual time stepping strategy ---*/
  if ((config_container[val_iZone]->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
      (config_container[val_iZone]->GetUnsteady_Simulation() == DT_STEPPING_2ND)) {
    
    /*--- Update dual time solver ---*/
    for (iMesh = 0; iMesh <= config_container[val_iZone]->GetnMGLevels(); iMesh++) {
      integration_container[val_iZone][WAVE_SOL]->SetDualTime_Solver(geometry_container[val_iZone][iMesh], solver_container[val_iZone][iMesh][WAVE_SOL], config_container[val_iZone], iMesh);
      integration_container[val_iZone][WAVE_SOL]->SetConvergence(false);
    }
    
    Physical_dt = config_container[val_iZone]->GetDelta_UnstTime(); Physical_t  = (ExtIter+1)*Physical_dt;
    if (Physical_t >=  config_container[val_iZone]->GetTotal_UnstTime()) integration_container[val_iZone][WAVE_SOL]->SetConvergence(true);
  }
}

void CWaveIteration::Monitor()     { }
void CWaveIteration::Output()      { }
void CWaveIteration::Postprocess() { }


CHeatIteration::CHeatIteration(CConfig *config) : CIteration(config) { }
CHeatIteration::~CHeatIteration(void) { }
void CHeatIteration::Preprocess(COutput *output,
                                CIntegration ***integration_container,
                                CGeometry ***geometry_container,
                                CSolver ****solver_container,
                                CNumerics *****numerics_container,
                                CConfig **config_container,
                                CSurfaceMovement **surface_movement,
                                CVolumetricMovement **grid_movement,
                                CFreeFormDefBox*** FFDBox,
                                unsigned short val_iZone) { }
void CHeatIteration::Iterate(COutput *output,
                             CIntegration ***integration_container,
                             CGeometry ***geometry_container,
                             CSolver ****solver_container,
                             CNumerics *****numerics_container,
                             CConfig **config_container,
                             CSurfaceMovement **surface_movement,
                             CVolumetricMovement **grid_movement,
                             CFreeFormDefBox*** FFDBox,
                             unsigned short val_iZone){
  
  unsigned long IntIter = 0; config_container[ZONE_0]->SetIntIter(IntIter);
  unsigned long ExtIter = config_container[ZONE_0]->GetExtIter();
  
  /*--- Set the value of the internal iteration ---*/
  IntIter = ExtIter;
  if ((config_container[val_iZone]->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
      (config_container[val_iZone]->GetUnsteady_Simulation() == DT_STEPPING_2ND)) IntIter = 0;
  
  /*--- Heat equation ---*/
  config_container[val_iZone]->SetGlobalParam(HEAT_EQUATION, RUNTIME_HEAT_SYS, ExtIter);
  integration_container[val_iZone][HEAT_SOL]->SingleGrid_Iteration(geometry_container, solver_container, numerics_container,
                                                                   config_container, RUNTIME_HEAT_SYS, IntIter, val_iZone);
  
  /*--- Dual time stepping strategy ---*/
  if ((config_container[val_iZone]->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
      (config_container[val_iZone]->GetUnsteady_Simulation() == DT_STEPPING_2ND)) {
    
    for (IntIter = 1; IntIter < config_container[val_iZone]->GetUnst_nIntIter(); IntIter++) {
      output->SetConvHistory_Body(NULL, geometry_container, solver_container, config_container, integration_container, true, 0.0, val_iZone);
      config_container[val_iZone]->SetIntIter(IntIter);
      integration_container[val_iZone][HEAT_SOL]->SingleGrid_Iteration(geometry_container, solver_container, numerics_container,
                                                                       config_container, RUNTIME_HEAT_SYS, IntIter, val_iZone);
      if (integration_container[val_iZone][HEAT_SOL]->GetConvergence()) break;
    }
  }
}

void CHeatIteration::Update(COutput *output,
                            CIntegration ***integration_container,
                            CGeometry ***geometry_container,
                            CSolver ****solver_container,
                            CNumerics *****numerics_container,
                            CConfig **config_container,
                            CSurfaceMovement **surface_movement,
                            CVolumetricMovement **grid_movement,
                            CFreeFormDefBox*** FFDBox,
                            unsigned short val_iZone)      {
  
  unsigned short iMesh;
  su2double Physical_dt, Physical_t;
  unsigned long ExtIter = config_container[ZONE_0]->GetExtIter();
  
  /*--- Dual time stepping strategy ---*/
  if ((config_container[val_iZone]->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
      (config_container[val_iZone]->GetUnsteady_Simulation() == DT_STEPPING_2ND)) {
    
    /*--- Update dual time solver ---*/
    for (iMesh = 0; iMesh <= config_container[val_iZone]->GetnMGLevels(); iMesh++) {
      integration_container[val_iZone][HEAT_SOL]->SetDualTime_Solver(geometry_container[val_iZone][iMesh], solver_container[val_iZone][iMesh][HEAT_SOL], config_container[val_iZone], iMesh);
      integration_container[val_iZone][HEAT_SOL]->SetConvergence(false);
    }
    
    Physical_dt = config_container[val_iZone]->GetDelta_UnstTime(); Physical_t  = (ExtIter+1)*Physical_dt;
    if (Physical_t >=  config_container[val_iZone]->GetTotal_UnstTime()) integration_container[val_iZone][HEAT_SOL]->SetConvergence(true);
  }
}
void CHeatIteration::Monitor()     { }
void CHeatIteration::Output()      { }
void CHeatIteration::Postprocess() { }


CPoissonIteration::CPoissonIteration(CConfig *config) : CIteration(config) { }
CPoissonIteration::~CPoissonIteration(void) { }
void CPoissonIteration::Preprocess(COutput *output,
                                   CIntegration ***integration_container,
                                   CGeometry ***geometry_container,
                                   CSolver ****solver_container,
                                   CNumerics *****numerics_container,
                                   CConfig **config_container,
                                   CSurfaceMovement **surface_movement,
                                   CVolumetricMovement **grid_movement,
                                   CFreeFormDefBox*** FFDBox,
                                   unsigned short val_iZone) { }
void CPoissonIteration::Iterate(COutput *output,
                                CIntegration ***integration_container,
                                CGeometry ***geometry_container,
                                CSolver ****solver_container,
                                CNumerics *****numerics_container,
                                CConfig **config_container,
                                CSurfaceMovement **surface_movement,
                                CVolumetricMovement **grid_movement,
                                CFreeFormDefBox*** FFDBox,
                                unsigned short val_iZone) {
  
  unsigned long IntIter = 0; config_container[ZONE_0]->SetIntIter(IntIter);
  unsigned long ExtIter = config_container[ZONE_0]->GetExtIter();
  
  /*--- Set the value of the internal iteration ---*/
  IntIter = ExtIter;
  if ((config_container[val_iZone]->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
      (config_container[val_iZone]->GetUnsteady_Simulation() == DT_STEPPING_2ND)) IntIter = 0;
  
  /*--- Poisson equation ---*/
  config_container[val_iZone]->SetGlobalParam(POISSON_EQUATION, RUNTIME_POISSON_SYS, ExtIter);
  integration_container[val_iZone][POISSON_SOL]->SingleGrid_Iteration(geometry_container, solver_container, numerics_container,
                                                                      config_container, RUNTIME_POISSON_SYS, IntIter, val_iZone);
  
  
}
void CPoissonIteration::Update(COutput *output,
                               CIntegration ***integration_container,
                               CGeometry ***geometry_container,
                               CSolver ****solver_container,
                               CNumerics *****numerics_container,
                               CConfig **config_container,
                               CSurfaceMovement **surface_movement,
                               CVolumetricMovement **grid_movement,
                               CFreeFormDefBox*** FFDBox,
                               unsigned short val_iZone)      { }
void CPoissonIteration::Monitor()     { }
void CPoissonIteration::Output()      { }
void CPoissonIteration::Postprocess() { }


CFEM_StructuralAnalysis::CFEM_StructuralAnalysis(CConfig *config) : CIteration(config) { }
CFEM_StructuralAnalysis::~CFEM_StructuralAnalysis(void) { }
void CFEM_StructuralAnalysis::Preprocess() { }
void CFEM_StructuralAnalysis::Iterate(COutput *output,
                         	 	  CIntegration ***integration_container,
                         	 	  CGeometry ***geometry_container,
                         	 	  CSolver ****solver_container,
                         	 	  CNumerics *****numerics_container,
                         	 	  CConfig **config_container,
                         	 	  CSurfaceMovement **surface_movement,
                         	 	  CVolumetricMovement **grid_movement,
                         	 	  CFreeFormDefBox*** FFDBox,
                                  unsigned short val_iZone
                         	 	  ) {

	int rank = MASTER_NODE;
#ifdef HAVE_MPI
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

	su2double loadIncrement;
	unsigned long IntIter = 0; config_container[val_iZone]->SetIntIter(IntIter);
	unsigned long ExtIter = config_container[val_iZone]->GetExtIter();

	bool fsi = config_container[val_iZone]->GetFSI_Simulation();

	unsigned long iIncrement;
	unsigned long nIncrements = config_container[val_iZone]->GetNumberIncrements();

	bool nonlinear = (config_container[val_iZone]->GetGeometricConditions() == LARGE_DEFORMATIONS);	// Geometrically non-linear problems
	bool linear = (config_container[val_iZone]->GetGeometricConditions() == SMALL_DEFORMATIONS);	// Geometrically non-linear problems

	bool initial_calc = config_container[val_iZone]->GetExtIter() == 0;				// Checks if it is the first calculation.
	bool first_iter = config_container[val_iZone]->GetIntIter() == 0;				// Checks if it is the first iteration
	bool restart = config_container[val_iZone]->GetRestart();												// Restart analysis
	bool initial_calc_restart = (SU2_TYPE::Int(config_container[val_iZone]->GetExtIter()) == config_container[val_iZone]->GetDyn_RestartIter()); // Initial calculation for restart

	su2double CurrentTime = config_container[val_iZone]->GetCurrent_DynTime();
	su2double Static_Time = config_container[val_iZone]->GetStatic_Time();

	bool statTime = (CurrentTime <= Static_Time);

	bool incremental_load = config_container[val_iZone]->GetIncrementalLoad();							// If an incremental load is applied

	/*--- This is to prevent problems when running a linear solver ---*/
	if (!nonlinear) incremental_load = false;

	/*--- Set the convergence monitor to false, to prevent the solver to stop in intermediate FSI subiterations ---*/
	integration_container[val_iZone][FEA_SOL]->SetConvergence(false);

	if (linear){

		/*--- Set the value of the internal iteration ---*/

		IntIter = ExtIter;

		/*--- FEA equations ---*/

		config_container[val_iZone]->SetGlobalParam(FEM_ELASTICITY, RUNTIME_FEA_SYS, ExtIter);

		/*--- Run the iteration ---*/

		integration_container[val_iZone][FEA_SOL]->Structural_Iteration(geometry_container, solver_container, numerics_container,
				config_container, RUNTIME_FEA_SYS, IntIter, val_iZone);

	}
	/*--- If the structure is held static and the solver is nonlinear, we don't need to solve for static time, but we need to compute Mass Matrix and Integration constants ---*/
	else if ((nonlinear) && ((!statTime) || (!fsi))){

		/*--- THIS IS THE DIRECT APPROACH (NO INCREMENTAL LOAD APPLIED) ---*/

		if (!incremental_load){

			/*--- Set the value of the internal iteration ---*/

			IntIter = 0;

			/*--- FEA equations ---*/

			config_container[val_iZone]->SetGlobalParam(FEM_ELASTICITY, RUNTIME_FEA_SYS, ExtIter);

			/*--- Run the iteration ---*/

			integration_container[val_iZone][FEA_SOL]->Structural_Iteration(geometry_container, solver_container, numerics_container,
					config_container, RUNTIME_FEA_SYS, IntIter, val_iZone);


			/*----------------- If the solver is non-linear, we need to subiterate using a Newton-Raphson approach ----------------------*/

			for (IntIter = 1; IntIter < config_container[val_iZone]->GetDyn_nIntIter(); IntIter++){

				/*--- Write the convergence history (only screen output) ---*/

				output->SetConvHistory_Body(NULL, geometry_container, solver_container, config_container, integration_container, true, 0.0, val_iZone);

				config_container[val_iZone]->SetIntIter(IntIter);

				integration_container[val_iZone][FEA_SOL]->Structural_Iteration(geometry_container, solver_container, numerics_container,
						config_container, RUNTIME_FEA_SYS, IntIter, val_iZone);

				if (integration_container[val_iZone][FEA_SOL]->GetConvergence()) break;

			}

		}
		/*--- The incremental load is only used in nonlinear cases ---*/
		else if (incremental_load){

			/*--- Set the initial condition: store the current solution as Solution_Old ---*/

			solver_container[val_iZone][MESH_0][FEA_SOL]->SetInitialCondition(geometry_container[val_iZone], solver_container[val_iZone], config_container[val_iZone], ExtIter);

			/*--- The load increment is 1.0 ---*/
			loadIncrement = 1.0;
			solver_container[val_iZone][MESH_0][FEA_SOL]->SetLoad_Increment(loadIncrement);

			/*--- Set the value of the internal iteration ---*/

			IntIter = 0;

			/*--- FEA equations ---*/

			config_container[val_iZone]->SetGlobalParam(FEM_ELASTICITY, RUNTIME_FEA_SYS, ExtIter);

			/*--- Run the first iteration ---*/

			integration_container[val_iZone][FEA_SOL]->Structural_Iteration(geometry_container, solver_container, numerics_container,
					config_container, RUNTIME_FEA_SYS, IntIter, val_iZone);


			/*--- Write the convergence history (only screen output) ---*/

			output->SetConvHistory_Body(NULL, geometry_container, solver_container, config_container, integration_container, true, 0.0, val_iZone);

			/*--- Run the second iteration ---*/

			IntIter = 1;

			config_container[val_iZone]->SetIntIter(IntIter);

			integration_container[val_iZone][FEA_SOL]->Structural_Iteration(geometry_container, solver_container, numerics_container,
					config_container, RUNTIME_FEA_SYS, IntIter, val_iZone);


			bool meetCriteria;
			su2double Residual_UTOL, Residual_RTOL, Residual_ETOL;
			su2double Criteria_UTOL, Criteria_RTOL, Criteria_ETOL;

			Criteria_UTOL = config_container[val_iZone]->GetIncLoad_Criteria(0);
			Criteria_RTOL = config_container[val_iZone]->GetIncLoad_Criteria(1);
			Criteria_ETOL = config_container[val_iZone]->GetIncLoad_Criteria(2);

			Residual_UTOL = log10(solver_container[val_iZone][MESH_0][FEA_SOL]->GetRes_FEM(0));
			Residual_RTOL = log10(solver_container[val_iZone][MESH_0][FEA_SOL]->GetRes_FEM(1));
			Residual_ETOL = log10(solver_container[val_iZone][MESH_0][FEA_SOL]->GetRes_FEM(2));

			meetCriteria = ( ( Residual_UTOL <  Criteria_UTOL ) &&
					( Residual_RTOL <  Criteria_RTOL ) &&
					( Residual_ETOL <  Criteria_ETOL ) );

			/*--- If the criteria is met and the load is not "too big", do the regular calculation ---*/
			if (meetCriteria){

				for (IntIter = 2; IntIter < config_container[val_iZone]->GetDyn_nIntIter(); IntIter++){

					/*--- Write the convergence history (only screen output) ---*/

					output->SetConvHistory_Body(NULL, geometry_container, solver_container, config_container, integration_container, true, 0.0, val_iZone);

					config_container[val_iZone]->SetIntIter(IntIter);

					integration_container[val_iZone][FEA_SOL]->Structural_Iteration(geometry_container, solver_container, numerics_container,
							config_container, RUNTIME_FEA_SYS, IntIter, val_iZone);

					if (integration_container[val_iZone][FEA_SOL]->GetConvergence()) break;

				}

			}

			/*--- If the criteria is not met, a whole set of subiterations for the different loads must be done ---*/

			else {

				/*--- Here we have to restart the solution to the original one of the iteration ---*/
				/*--- Retrieve the Solution_Old as the current solution before subiterating ---*/

				solver_container[val_iZone][MESH_0][FEA_SOL]->ResetInitialCondition(geometry_container[val_iZone], solver_container[val_iZone], config_container[val_iZone], ExtIter);

				/*--- For the number of increments ---*/
				for (iIncrement = 0; iIncrement < nIncrements; iIncrement++){

					loadIncrement = (iIncrement + 1.0) * (1.0 / nIncrements);

					/*--- Set the load increment and the initial condition, and output the parameters of UTOL, RTOL, ETOL for the previous iteration ---*/

					/*--- Set the convergence monitor to false, to force se solver to converge every subiteration ---*/
					integration_container[val_iZone][FEA_SOL]->SetConvergence(false);

					output->SetConvHistory_Body(NULL, geometry_container, solver_container, config_container, integration_container, true, 0.0, val_iZone);

					/*--- FEA equations ---*/

					config_container[val_iZone]->SetGlobalParam(FEM_ELASTICITY, RUNTIME_FEA_SYS, ExtIter);


					solver_container[val_iZone][MESH_0][FEA_SOL]->SetLoad_Increment(loadIncrement);

					if (rank == MASTER_NODE){
						cout << endl;
						cout << "-- Incremental load: increment " << iIncrement + 1 << " ------------------------------------------" << endl;
					}

					/*--- Set the value of the internal iteration ---*/
					IntIter = 0;
					config_container[val_iZone]->SetIntIter(IntIter);

					/*--- FEA equations ---*/

					config_container[val_iZone]->SetGlobalParam(FEM_ELASTICITY, RUNTIME_FEA_SYS, ExtIter);

					/*--- Run the iteration ---*/

					integration_container[val_iZone][FEA_SOL]->Structural_Iteration(geometry_container, solver_container, numerics_container,
							config_container, RUNTIME_FEA_SYS, IntIter, val_iZone);


					/*----------------- If the solver is non-linear, we need to subiterate using a Newton-Raphson approach ----------------------*/

					for (IntIter = 1; IntIter < config_container[val_iZone]->GetDyn_nIntIter(); IntIter++){

						/*--- Write the convergence history (only screen output) ---*/

						output->SetConvHistory_Body(NULL, geometry_container, solver_container, config_container, integration_container, true, 0.0, val_iZone);

						config_container[val_iZone]->SetIntIter(IntIter);

						integration_container[val_iZone][FEA_SOL]->Structural_Iteration(geometry_container, solver_container, numerics_container,
								config_container, RUNTIME_FEA_SYS, IntIter, val_iZone);

						if (integration_container[val_iZone][FEA_SOL]->GetConvergence()) break;

					}

				}

			}

		}


	}
	else if (
			(nonlinear && statTime) &&
			((first_iter && initial_calc) || (restart && initial_calc_restart))
	){

		/*--- We need to do the preprocessing to compute the Mass Matrix and integration constants ---*/
		solver_container[val_iZone][MESH_0][FEA_SOL]->Preprocessing(geometry_container[val_iZone][MESH_0], solver_container[val_iZone][MESH_0],
				config_container[val_iZone], numerics_container[val_iZone][MESH_0][FEA_SOL], MESH_0, 0, RUNTIME_FEA_SYS, false);

	}

}

void CFEM_StructuralAnalysis::Update(COutput *output,
	 	  CIntegration ***integration_container,
	 	  CGeometry ***geometry_container,
	 	  CSolver ****solver_container,
	 	  CNumerics *****numerics_container,
	 	  CConfig **config_container,
	 	  CSurfaceMovement **surface_movement,
	 	  CVolumetricMovement **grid_movement,
	 	  CFreeFormDefBox*** FFDBox,
	 	  unsigned short val_iZone){

	su2double Physical_dt, Physical_t;
  	unsigned long ExtIter = config_container[val_iZone]->GetExtIter();
	bool dynamic = (config_container[val_iZone]->GetDynamic_Analysis() == DYNAMIC);					// Dynamic problems

	/*----------------- Compute averaged nodal stress and reactions ------------------------*/

	solver_container[val_iZone][MESH_0][FEA_SOL]->Compute_NodalStress(geometry_container[val_iZone][MESH_0], solver_container[val_iZone][MESH_0], numerics_container[val_iZone][MESH_0][FEA_SOL], config_container[val_iZone]);

	/*----------------- Update structural solver ----------------------*/

	if (dynamic){
		integration_container[val_iZone][FEA_SOL]->SetFEM_StructuralSolver(geometry_container[val_iZone][MESH_0], solver_container[val_iZone][MESH_0], config_container[val_iZone], MESH_0);
		integration_container[val_iZone][FEA_SOL]->SetConvergence(false);

	    /*--- Verify convergence criteria (based on total time) ---*/

		Physical_dt = config_container[val_iZone]->GetDelta_DynTime();
		Physical_t  = (ExtIter+1)*Physical_dt;
		if (Physical_t >=  config_container[val_iZone]->GetTotal_DynTime())
			integration_container[val_iZone][FEA_SOL]->SetConvergence(true);
	}

}
void CFEM_StructuralAnalysis::Monitor()     { }
void CFEM_StructuralAnalysis::Output()      { }
void CFEM_StructuralAnalysis::Postprocess() { }


CAdjMeanFlowIteration::CAdjMeanFlowIteration(CConfig *config) : CIteration(config) { }
CAdjMeanFlowIteration::~CAdjMeanFlowIteration(void) { }
void CAdjMeanFlowIteration::Preprocess(COutput *output,
                                       CIntegration ***integration_container,
                                       CGeometry ***geometry_container,
                                       CSolver ****solver_container,
                                       CNumerics *****numerics_container,
                                       CConfig **config_container,
                                       CSurfaceMovement **surface_movement,
                                       CVolumetricMovement **grid_movement,
                                       CFreeFormDefBox*** FFDBox,
                                       unsigned short val_iZone) {
  
  unsigned short iMesh;
  bool time_spectral = (config_container[ZONE_0]->GetUnsteady_Simulation() == TIME_SPECTRAL);
  bool dynamic_mesh = config_container[ZONE_0]->GetGrid_Movement();
  unsigned long IntIter = 0; config_container[ZONE_0]->SetIntIter(IntIter);
  unsigned long ExtIter = config_container[ZONE_0]->GetExtIter();
  
  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  
  /*--- For the unsteady adjoint, load a new direct solution from a restart file. ---*/
  
  if (((dynamic_mesh && ExtIter == 0) || config_container[val_iZone]->GetUnsteady_Simulation()) && !time_spectral) {
    int Direct_Iter = SU2_TYPE::Int(config_container[val_iZone]->GetUnst_AdjointIter()) - SU2_TYPE::Int(ExtIter) - 1;
    if (rank == MASTER_NODE && val_iZone == ZONE_0 && config_container[val_iZone]->GetUnsteady_Simulation())
      cout << endl << " Loading flow solution from direct iteration " << Direct_Iter << "." << endl;
    solver_container[val_iZone][MESH_0][FLOW_SOL]->LoadRestart(geometry_container[val_iZone], solver_container[val_iZone], config_container[val_iZone], Direct_Iter);
  }
  
  /*--- Continuous adjoint Euler, Navier-Stokes or Reynolds-averaged Navier-Stokes (RANS) equations ---*/
  
  if ((ExtIter == 0) || config_container[val_iZone]->GetUnsteady_Simulation()) {
    
    if (config_container[val_iZone]->GetKind_Solver() == ADJ_EULER)
      config_container[val_iZone]->SetGlobalParam(ADJ_EULER, RUNTIME_FLOW_SYS, ExtIter);
    if (config_container[val_iZone]->GetKind_Solver() == ADJ_NAVIER_STOKES)
      config_container[val_iZone]->SetGlobalParam(ADJ_NAVIER_STOKES, RUNTIME_FLOW_SYS, ExtIter);
    if (config_container[val_iZone]->GetKind_Solver() == ADJ_RANS)
      config_container[val_iZone]->SetGlobalParam(ADJ_RANS, RUNTIME_FLOW_SYS, ExtIter);
    
    /*--- Solve the Euler, Navier-Stokes or Reynolds-averaged Navier-Stokes (RANS) equations (one iteration) ---*/
    
    if (rank == MASTER_NODE && val_iZone == ZONE_0)
      cout << "Begin direct solver to store flow data (single iteration)." << endl;
    
    if (rank == MASTER_NODE && val_iZone == ZONE_0)
      cout << "Compute residuals to check the convergence of the direct problem." << endl;
    
    integration_container[val_iZone][FLOW_SOL]->MultiGrid_Iteration(geometry_container, solver_container, numerics_container,
                                                                    config_container, RUNTIME_FLOW_SYS, 0, val_iZone);
    
    if (config_container[val_iZone]->GetKind_Solver() == ADJ_RANS) {
      
      /*--- Solve the turbulence model ---*/
      
      config_container[val_iZone]->SetGlobalParam(ADJ_RANS, RUNTIME_TURB_SYS, ExtIter);
      integration_container[val_iZone][TURB_SOL]->SingleGrid_Iteration(geometry_container, solver_container, numerics_container,
                                                                       config_container, RUNTIME_TURB_SYS, IntIter, val_iZone);
      
      /*--- Solve transition model ---*/
      
      if (config_container[val_iZone]->GetKind_Trans_Model() == LM) {
        config_container[val_iZone]->SetGlobalParam(RANS, RUNTIME_TRANS_SYS, ExtIter);
        integration_container[val_iZone][TRANS_SOL]->SingleGrid_Iteration(geometry_container, solver_container, numerics_container,
                                                                          config_container, RUNTIME_TRANS_SYS, IntIter, val_iZone);
      }
      
    }
    
    /*--- Output the residual (visualization purpouses to identify if
     the direct solution is converged)---*/
    if (rank == MASTER_NODE && val_iZone == ZONE_0)
      cout << "log10[Maximum residual]: " << log10(solver_container[val_iZone][MESH_0][FLOW_SOL]->GetRes_Max(0))
      <<", located at point "<< solver_container[val_iZone][MESH_0][FLOW_SOL]->GetPoint_Max(0) << "." << endl;
    
    /*--- Compute gradients of the flow variables, this is necessary for sensitivity computation,
     note that in the direct Euler problem we are not computing the gradients of the primitive variables ---*/
    
    if (config_container[val_iZone]->GetKind_Gradient_Method() == GREEN_GAUSS)
      solver_container[val_iZone][MESH_0][FLOW_SOL]->SetPrimitive_Gradient_GG(geometry_container[val_iZone][MESH_0], config_container[val_iZone]);
    if (config_container[val_iZone]->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES)
      solver_container[val_iZone][MESH_0][FLOW_SOL]->SetPrimitive_Gradient_LS(geometry_container[val_iZone][MESH_0], config_container[val_iZone]);
    
    /*--- Set contribution from cost function for boundary conditions ---*/
    
    for (iMesh = 0; iMesh <= config_container[val_iZone]->GetnMGLevels(); iMesh++) {
      
      /*--- Set the value of the non-dimensional coefficients in the coarse levels, using the fine level solution ---*/
      
      solver_container[val_iZone][iMesh][FLOW_SOL]->SetTotal_CDrag(solver_container[val_iZone][MESH_0][FLOW_SOL]->GetTotal_CDrag());
      solver_container[val_iZone][iMesh][FLOW_SOL]->SetTotal_CLift(solver_container[val_iZone][MESH_0][FLOW_SOL]->GetTotal_CLift());
      solver_container[val_iZone][iMesh][FLOW_SOL]->SetTotal_CT(solver_container[val_iZone][MESH_0][FLOW_SOL]->GetTotal_CT());
      solver_container[val_iZone][iMesh][FLOW_SOL]->SetTotal_CQ(solver_container[val_iZone][MESH_0][FLOW_SOL]->GetTotal_CQ());
      
      /*--- Compute the adjoint boundary condition on Euler walls ---*/
      
      solver_container[val_iZone][iMesh][ADJFLOW_SOL]->SetForceProj_Vector(geometry_container[val_iZone][iMesh], solver_container[val_iZone][iMesh], config_container[val_iZone]);
      
      /*--- Set the internal boundary condition on nearfield surfaces ---*/
      
      if ((config_container[val_iZone]->GetKind_ObjFunc() == EQUIVALENT_AREA) ||
          (config_container[val_iZone]->GetKind_ObjFunc() == NEARFIELD_PRESSURE))
        solver_container[val_iZone][iMesh][ADJFLOW_SOL]->SetIntBoundary_Jump(geometry_container[val_iZone][iMesh], solver_container[val_iZone][iMesh], config_container[val_iZone]);
      
    }
    
    if (rank == MASTER_NODE && val_iZone == ZONE_0)
      cout << "End direct solver, begin adjoint problem." << endl;
    
  }
  
}
void CAdjMeanFlowIteration::Iterate(COutput *output,
                                    CIntegration ***integration_container,
                                    CGeometry ***geometry_container,
                                    CSolver ****solver_container,
                                    CNumerics *****numerics_container,
                                    CConfig **config_container,
                                    CSurfaceMovement **surface_movement,
                                    CVolumetricMovement **grid_movement,
                                    CFreeFormDefBox*** FFDBox,
                                    unsigned short val_iZone) {
  
  unsigned long IntIter = 0; config_container[ZONE_0]->SetIntIter(IntIter);
  unsigned long ExtIter = config_container[ZONE_0]->GetExtIter();
  
  /*--- Set the value of the internal iteration ---*/
  
  IntIter = ExtIter;
  if ((config_container[val_iZone]->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
      (config_container[val_iZone]->GetUnsteady_Simulation() == DT_STEPPING_2ND)) {
    IntIter = 0;
  }
  
  if (config_container[val_iZone]->GetKind_Solver() == ADJ_EULER)
    config_container[val_iZone]->SetGlobalParam(ADJ_EULER, RUNTIME_ADJFLOW_SYS, ExtIter);
  if (config_container[val_iZone]->GetKind_Solver() == ADJ_NAVIER_STOKES)
    config_container[val_iZone]->SetGlobalParam(ADJ_NAVIER_STOKES, RUNTIME_ADJFLOW_SYS, ExtIter);
  if (config_container[val_iZone]->GetKind_Solver() == ADJ_RANS)
    config_container[val_iZone]->SetGlobalParam(ADJ_RANS, RUNTIME_ADJFLOW_SYS, ExtIter);
  
  /*--- Iteration of the flow adjoint problem ---*/
  
  integration_container[val_iZone][ADJFLOW_SOL]->MultiGrid_Iteration(geometry_container, solver_container, numerics_container,
                                                                     config_container, RUNTIME_ADJFLOW_SYS, IntIter, val_iZone);
  
  /*--- Iteration of the turbulence model adjoint ---*/
  
  if ((config_container[val_iZone]->GetKind_Solver() == ADJ_RANS) && (!config_container[val_iZone]->GetFrozen_Visc())) {
    
    /*--- Adjoint turbulence model solution ---*/
    
    config_container[val_iZone]->SetGlobalParam(ADJ_RANS, RUNTIME_ADJTURB_SYS, ExtIter);
    integration_container[val_iZone][ADJTURB_SOL]->SingleGrid_Iteration(geometry_container, solver_container, numerics_container,
                                                                        config_container, RUNTIME_ADJTURB_SYS, IntIter, val_iZone);
    
  }
  
  /*--- Dual time stepping strategy ---*/
  
  if ((config_container[val_iZone]->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
      (config_container[val_iZone]->GetUnsteady_Simulation() == DT_STEPPING_2ND)) {
    
    for (IntIter = 1; IntIter < config_container[val_iZone]->GetUnst_nIntIter(); IntIter++) {
      
      /*--- Write the convergence history (only screen output) ---*/
      
      output->SetConvHistory_Body(NULL, geometry_container, solver_container, config_container, integration_container, true, 0.0, val_iZone);
      
      /*--- Set the value of the internal iteration ---*/
      
      config_container[val_iZone]->SetIntIter(IntIter);
      
      /*--- All zones must be advanced and coupled with each pseudo timestep ---*/
      
      integration_container[val_iZone][ADJFLOW_SOL]->MultiGrid_Iteration(geometry_container, solver_container, numerics_container,
                                                                         config_container, RUNTIME_ADJFLOW_SYS, IntIter, val_iZone);
      
      /*--- Check to see if the convergence criteria has been met ---*/
      
      if (integration_container[val_iZone][ADJFLOW_SOL]->GetConvergence()) break;
    }
    
  }
  
}
void CAdjMeanFlowIteration::Update(COutput *output,
                                   CIntegration ***integration_container,
                                   CGeometry ***geometry_container,
                                   CSolver ****solver_container,
                                   CNumerics *****numerics_container,
                                   CConfig **config_container,
                                   CSurfaceMovement **surface_movement,
                                   CVolumetricMovement **grid_movement,
                                   CFreeFormDefBox*** FFDBox,
                                   unsigned short val_iZone)      {
  
  su2double Physical_dt, Physical_t;
  unsigned short iMesh;
  unsigned long ExtIter = config_container[ZONE_0]->GetExtIter();
  
  /*--- Dual time stepping strategy ---*/
  
  if ((config_container[val_iZone]->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
      (config_container[val_iZone]->GetUnsteady_Simulation() == DT_STEPPING_2ND)) {
    
    /*--- Update dual time solver ---*/
    
    for (iMesh = 0; iMesh <= config_container[val_iZone]->GetnMGLevels(); iMesh++) {
      integration_container[val_iZone][ADJFLOW_SOL]->SetDualTime_Solver(geometry_container[val_iZone][iMesh], solver_container[val_iZone][iMesh][ADJFLOW_SOL], config_container[val_iZone], iMesh);
      integration_container[val_iZone][ADJFLOW_SOL]->SetConvergence(false);
    }
    
    Physical_dt = config_container[val_iZone]->GetDelta_UnstTime(); Physical_t  = (ExtIter+1)*Physical_dt;
    if (Physical_t >=  config_container[val_iZone]->GetTotal_UnstTime()) integration_container[val_iZone][ADJFLOW_SOL]->SetConvergence(true);
    
  }
}

void CAdjMeanFlowIteration::Monitor()     { }
void CAdjMeanFlowIteration::Output()      { }
void CAdjMeanFlowIteration::Postprocess() { }

CDiscAdjMeanFlowIteration::CDiscAdjMeanFlowIteration(CConfig *config) : CIteration(config), CurrentRecording(NONE){
  
  meanflow_iteration = new CMeanFlowIteration(config);
  
  turbulent = config->GetKind_Solver() == DISC_ADJ_RANS;
  
}

CDiscAdjMeanFlowIteration::~CDiscAdjMeanFlowIteration(void) { }
void CDiscAdjMeanFlowIteration::Preprocess(COutput *output,
                                           CIntegration ***integration_container,
                                           CGeometry ***geometry_container,
                                           CSolver ****solver_container,
                                           CNumerics *****numerics_container,
                                           CConfig **config_container,
                                           CSurfaceMovement **surface_movement,
                                           CVolumetricMovement **grid_movement,
                                           CFreeFormDefBox*** FFDBox,
                                           unsigned short val_iZone) {
  
  unsigned long IntIter = 0, iPoint;
  config_container[ZONE_0]->SetIntIter(IntIter);
  unsigned short ExtIter = config_container[val_iZone]->GetExtIter();
  bool unsteady = config_container[val_iZone]->GetUnsteady_Simulation() != NONE;
  bool dual_time_1st = (config_container[val_iZone]->GetUnsteady_Simulation() == DT_STEPPING_1ST);
  bool dual_time_2nd = (config_container[val_iZone]->GetUnsteady_Simulation() == DT_STEPPING_2ND);
  bool dual_time = (dual_time_1st || dual_time_2nd);
  unsigned short iMesh;
  int Direct_Iter;

  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  /*--- For the unsteady adjoint, load direct solutions from restart files. ---*/

  if (config_container[val_iZone]->GetUnsteady_Simulation()) {

    Direct_Iter = SU2_TYPE::Int(config_container[val_iZone]->GetUnst_AdjointIter()) - SU2_TYPE::Int(ExtIter) - 2;

    /*--- For dual-time stepping we want to load the already converged solution at timestep n ---*/

    if (dual_time){
      Direct_Iter += 1;
    }

    if (dual_time_2nd){

      /*--- Load solution at timestep n-2 ---*/

      LoadUnsteady_Solution(geometry_container, solver_container,config_container, val_iZone, Direct_Iter-2);

      /*--- Push solution back to correct array ---*/

      for (iMesh=0; iMesh<=config_container[val_iZone]->GetnMGLevels();iMesh++){
        for(iPoint=0; iPoint<geometry_container[val_iZone][iMesh]->GetnPoint();iPoint++){
          solver_container[val_iZone][iMesh][FLOW_SOL]->node[iPoint]->Set_Solution_time_n();
          solver_container[val_iZone][iMesh][FLOW_SOL]->node[iPoint]->Set_Solution_time_n1();
          if (turbulent){
            solver_container[val_iZone][iMesh][TURB_SOL]->node[iPoint]->Set_Solution_time_n();
            solver_container[val_iZone][iMesh][TURB_SOL]->node[iPoint]->Set_Solution_time_n1();
          }
        }
      }
    }
    if (dual_time){

      /*--- Load solution at timestep n-1 ---*/

      LoadUnsteady_Solution(geometry_container, solver_container,config_container, val_iZone, Direct_Iter-1);

      /*--- Push solution back to correct array ---*/

      for (iMesh=0; iMesh<=config_container[val_iZone]->GetnMGLevels();iMesh++){
        for(iPoint=0; iPoint<geometry_container[val_iZone][iMesh]->GetnPoint();iPoint++){
          solver_container[val_iZone][iMesh][FLOW_SOL]->node[iPoint]->Set_Solution_time_n();
          if (turbulent){
            solver_container[val_iZone][iMesh][TURB_SOL]->node[iPoint]->Set_Solution_time_n();
          }
        }
      }
    }

    /*--- Load solution timestep n ---*/

    LoadUnsteady_Solution(geometry_container, solver_container,config_container, val_iZone, Direct_Iter);


    /*--- Store flow solution also in the adjoint solver in order to be able to reset it later ---*/

    for (iPoint = 0; iPoint < geometry_container[val_iZone][MESH_0]->GetnPoint(); iPoint++){
      solver_container[val_iZone][MESH_0][ADJFLOW_SOL]->node[iPoint]->SetSolution_Direct(solver_container[val_iZone][MESH_0][FLOW_SOL]->node[iPoint]->GetSolution());
    }
    if (turbulent){
      for (iPoint = 0; iPoint < geometry_container[val_iZone][MESH_0]->GetnPoint(); iPoint++){
        solver_container[val_iZone][MESH_0][ADJTURB_SOL]->node[iPoint]->SetSolution_Direct(solver_container[val_iZone][MESH_0][TURB_SOL]->node[iPoint]->GetSolution());
      }
    }
  }

  solver_container[val_iZone][MESH_0][ADJFLOW_SOL]->Preprocessing(geometry_container[val_iZone][MESH_0], solver_container[val_iZone][MESH_0],  config_container[val_iZone] , MESH_0, 0, RUNTIME_ADJFLOW_SYS, false);
  if (turbulent){
    solver_container[val_iZone][MESH_0][ADJTURB_SOL]->Preprocessing(geometry_container[val_iZone][MESH_0], solver_container[val_iZone][MESH_0],  config_container[val_iZone] , MESH_0, 0, RUNTIME_ADJTURB_SYS, false);
  }

//  if (CurrentRecording != FLOW_VARIABLES || unsteady){
    
    if (rank == MASTER_NODE){
      cout << "Direct iteration to store computational graph." << endl;
      cout << "Compute residuals to check the convergence of the direct problem." << endl;
    }
    
    /*--- Record one mean flow iteration with flow variables as input ---*/
    
    SetRecording(output, integration_container, geometry_container, solver_container, numerics_container,
                 config_container, surface_movement, grid_movement, FFDBox, val_iZone, ALL_VARIABLES);//FLOW_VARIABLES
    
    /*--- Print residuals in the first iteration ---*/
    
    if (rank == MASTER_NODE && ((ExtIter == 0) || unsteady )){
      cout << "log10[RMS Density]: "<< log10(solver_container[val_iZone][MESH_0][FLOW_SOL]->GetRes_RMS(0))
           <<", Drag: " <<solver_container[val_iZone][MESH_0][FLOW_SOL]->GetTotal_CDrag()
          <<", Lift: " << solver_container[val_iZone][MESH_0][FLOW_SOL]->GetTotal_CLift() << "." << endl;

      if (turbulent){
        cout << "log10[RMS k]: " << log10(solver_container[val_iZone][MESH_0][TURB_SOL]->GetRes_RMS(0)) << endl;
      }
    }
//  }
}



void CDiscAdjMeanFlowIteration::LoadUnsteady_Solution(CGeometry ***geometry_container,
                                           CSolver ****solver_container,
                                           CConfig **config_container,
                                           unsigned short val_iZone, int val_DirectIter) {
  unsigned short iMesh;

  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  if (val_DirectIter >= 0){
    if (rank == MASTER_NODE && val_iZone == ZONE_0)
      cout << " Loading flow solution from direct iteration " << val_DirectIter  << "." << endl;
    solver_container[val_iZone][MESH_0][FLOW_SOL]->LoadRestart(geometry_container[val_iZone], solver_container[val_iZone], config_container[val_iZone], val_DirectIter);
    solver_container[val_iZone][MESH_0][FLOW_SOL]->Preprocessing(geometry_container[val_iZone][MESH_0],solver_container[val_iZone][MESH_0], config_container[val_iZone], MESH_0, val_DirectIter, RUNTIME_FLOW_SYS, false);
    if (turbulent){
      solver_container[val_iZone][MESH_0][TURB_SOL]->LoadRestart(geometry_container[val_iZone], solver_container[val_iZone], config_container[val_iZone], val_DirectIter);
      solver_container[val_iZone][MESH_0][TURB_SOL]->Postprocessing(geometry_container[val_iZone][MESH_0],solver_container[val_iZone][MESH_0], config_container[val_iZone], MESH_0);
    }
  } else {
    /*--- If there is no solution file we set the freestream condition ---*/
    if (rank == MASTER_NODE && val_iZone == ZONE_0)
      cout << " Setting freestream conditions at direct iteration " << val_DirectIter << "." << endl;
    for (iMesh=0; iMesh<=config_container[val_iZone]->GetnMGLevels();iMesh++){
      solver_container[val_iZone][iMesh][FLOW_SOL]->SetFreeStream_Solution(config_container[val_iZone]);
      solver_container[val_iZone][iMesh][FLOW_SOL]->Preprocessing(geometry_container[val_iZone][iMesh],solver_container[val_iZone][iMesh], config_container[val_iZone], iMesh, val_DirectIter, RUNTIME_FLOW_SYS, false);
      if (turbulent){
        solver_container[val_iZone][iMesh][TURB_SOL]->SetFreeStream_Solution(config_container[val_iZone]);
        solver_container[val_iZone][iMesh][TURB_SOL]->Postprocessing(geometry_container[val_iZone][iMesh],solver_container[val_iZone][iMesh], config_container[val_iZone], iMesh);
      }
    }
  }
}


void CDiscAdjMeanFlowIteration::Iterate(COutput *output,
                                        CIntegration ***integration_container,
                                        CGeometry ***geometry_container,
                                        CSolver ****solver_container,
                                        CNumerics *****numerics_container,
                                        CConfig **config_container,
                                        CSurfaceMovement **surface_movement,
                                        CVolumetricMovement **volume_grid_movement,
                                        CFreeFormDefBox*** FFDBox,
                                        unsigned short val_iZone) {
  
  unsigned long ExtIter = config_container[ZONE_0]->GetExtIter();
  unsigned long IntIter=0, nIntIter = 1;
  bool dual_time_1st = (config_container[ZONE_0]->GetUnsteady_Simulation() == DT_STEPPING_1ST);
  bool dual_time_2nd = (config_container[ZONE_0]->GetUnsteady_Simulation() == DT_STEPPING_2ND);
  bool dual_time = (dual_time_1st || dual_time_2nd);


  config_container[val_iZone]->SetIntIter(IntIter);

  if(dual_time)
    nIntIter = config_container[val_iZone]->GetUnst_nIntIter();


  for(IntIter=0; IntIter< nIntIter; IntIter++){

    /*--- Set the internal iteration ---*/

    config_container[val_iZone]->SetIntIter(IntIter);

    /*--- Set the adjoint values of the flow and objective function ---*/

    InitializeAdjoint(solver_container, geometry_container, config_container, val_iZone);

    /*--- Run the adjoint computation ---*/

    AD::ComputeAdjoint();

    /*--- Extract the adjoints of the conservative input variables and store them for the next iteration ---*/

    solver_container[val_iZone][MESH_0][ADJFLOW_SOL]->ExtractAdjoint_Solution(geometry_container[val_iZone][MESH_0],
                                                                              config_container[val_iZone]);

    solver_container[val_iZone][MESH_0][ADJFLOW_SOL]->ExtractAdjoint_Variables(geometry_container[val_iZone][MESH_0],
                                                                               config_container[val_iZone]);

    if (config_container[ZONE_0]->GetKind_Solver() == DISC_ADJ_RANS) {
      solver_container[val_iZone][MESH_0][ADJTURB_SOL]->ExtractAdjoint_Solution(geometry_container[val_iZone][MESH_0],
                                                                                config_container[val_iZone]);
    }

    solver_container[val_iZone][MESH_0][ADJFLOW_SOL]->SetSensitivity(geometry_container[val_iZone][MESH_0],config_container[val_iZone]);

    /*--- Clear all adjoints to re-use the stored computational graph in the next iteration ---*/

    AD::ClearAdjoints();

    /*--- Set the convergence criteria (only residual possible) ---*/

    integration_container[val_iZone][ADJFLOW_SOL]->Convergence_Monitoring(geometry_container[val_iZone][MESH_0],config_container[val_iZone],
                                                                          IntIter,log10(solver_container[val_iZone][MESH_0][ADJFLOW_SOL]->GetRes_RMS(0)), MESH_0);

    if(integration_container[val_iZone][ADJFLOW_SOL]->GetConvergence()){
      break;
    }

    /*--- Write the convergence history (only screen output) ---*/

    if(dual_time && (IntIter != nIntIter-1))
      output->SetConvHistory_Body(NULL, geometry_container, solver_container, config_container, integration_container, true, 0.0, val_iZone);

  }


  if (dual_time){
    integration_container[val_iZone][ADJFLOW_SOL]->SetConvergence(false);
  }

  cout << "log10Ajdoint[RMS Density]: " << log10(solver_container[ZONE_0][MESH_0][ADJFLOW_SOL]->GetRes_RMS(0))<<std::endl;
  cout << "log10Primal[RMS Density]: " << log10(solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetRes_RMS(0))<<std::endl;
  
//  if ((ExtIter+1 >= config_container[val_iZone]->GetnExtIter()) || (integration_container[val_iZone][ADJFLOW_SOL]->GetConvergence()) ||
//      ((ExtIter % config_container[val_iZone]->GetWrt_Sol_Freq() == 0))){
    
//    /*--- Record one mean flow iteration with geometry variables as input ---*/
    
//    SetRecording(output, integration_container, geometry_container, solver_container, numerics_container,
//                 config_container, surface_movement, volume_grid_movement, FFDBox, val_iZone, GEOMETRY_VARIABLES);
    
//    /*--- Set the adjoint values of the flow and objective function ---*/
    
//    InitializeAdjoint(solver_container, geometry_container, config_container, val_iZone);
    
//    /*--- Run the adjoint computation ---*/
    
//    AD::ComputeAdjoint();
    
//    /*--- Extract the sensitivities (adjoint of node coordinates) ---*/
    
//    solver_container[val_iZone][MESH_0][ADJFLOW_SOL]->SetSensitivity(geometry_container[val_iZone][MESH_0],config_container[val_iZone]);
    
//  }
  
}

void CDiscAdjMeanFlowIteration::SetRecording(COutput *output,
                                             CIntegration ***integration_container,
                                             CGeometry ***geometry_container,
                                             CSolver ****solver_container,
                                             CNumerics *****numerics_container,
                                             CConfig **config_container,
                                             CSurfaceMovement **surface_movement,
                                             CVolumetricMovement **grid_movement,
                                             CFreeFormDefBox*** FFDBox,
                                             unsigned short val_iZone,
                                             unsigned short kind_recording)      {
  
  unsigned long IntIter = config_container[ZONE_0]->GetIntIter();
  unsigned long ExtIter = config_container[val_iZone]->GetExtIter(), DirectExtIter;
  bool unsteady = config_container[val_iZone]->GetUnsteady_Simulation() != NONE;
  unsigned short iMesh;

  DirectExtIter = 0;
  if (unsteady){
    DirectExtIter = SU2_TYPE::Int(config_container[val_iZone]->GetUnst_AdjointIter()) - SU2_TYPE::Int(ExtIter) - 1;
  }

  /*--- Reset the tape ---*/

  AD::Reset();

  /*--- We only need to reset the indices if the current recording is different from the recording we want to have ---*/

//  if (CurrentRecording != kind_recording && (CurrentRecording != NONE) ){

//    for (iMesh = 0; iMesh <= config_container[val_iZone]->GetnMGLevels(); iMesh++){
//      solver_container[val_iZone][iMesh][ADJFLOW_SOL]->SetRecording(geometry_container[val_iZone][MESH_0], config_container[val_iZone], kind_recording);
//    }

//    if (turbulent){
//      solver_container[val_iZone][MESH_0][ADJTURB_SOL]->SetRecording(geometry_container[val_iZone][MESH_0], config_container[val_iZone], kind_recording);
//    }

//    /*--- Clear indices of coupling variables ---*/

//    SetDependencies(solver_container, geometry_container, config_container, val_iZone, ALL_VARIABLES);
//
    /*--- Run one iteration while tape is passive - this clears all indices ---*/

//    meanflow_iteration->Iterate(output,integration_container,geometry_container,solver_container,numerics_container,
//                                config_container,surface_movement,grid_movement,FFDBox,val_iZone);

//  }
    /*--- Prepare for recording ---*/

    for (iMesh = 0; iMesh <= config_container[val_iZone]->GetnMGLevels(); iMesh++){
      solver_container[val_iZone][iMesh][ADJFLOW_SOL]->SetRecording(geometry_container[val_iZone][MESH_0], config_container[val_iZone], kind_recording);
    }

    if (turbulent){
      solver_container[val_iZone][MESH_0][ADJTURB_SOL]->SetRecording(geometry_container[val_iZone][MESH_0], config_container[val_iZone], kind_recording);
    }


  /*--- Start the recording of all operations ---*/
  
  AD::StartRecording();
  
  /*--- Register flow variables ---*/
  
  RegisterInput(solver_container, geometry_container, config_container, val_iZone, kind_recording);
  
  /*--- Compute coupling or update the geometry ---*/

  SetDependencies(solver_container, geometry_container, config_container, val_iZone, kind_recording);
  
  /*--- Set the correct direct iteration number ---*/

  if (unsteady){
    config_container[val_iZone]->SetExtIter(DirectExtIter);
  }

  /*--- Run the direct iteration ---*/

  meanflow_iteration->Iterate(output,integration_container,geometry_container,solver_container,numerics_container,
                              config_container,surface_movement,grid_movement,FFDBox, val_iZone);

  config_container[val_iZone]->SetExtIter(ExtIter);

  /*--- Register flow variables and objective function as output ---*/
  
  /*--- For flux-avg or area-avg objective functions the 1D values must be calculated first ---*/
  if (config_container[val_iZone]->GetKind_ObjFunc()==AVG_OUTLET_PRESSURE ||
      config_container[val_iZone]->GetKind_ObjFunc()==AVG_TOTAL_PRESSURE ||
      config_container[val_iZone]->GetKind_ObjFunc()==MASS_FLOW_RATE)
    output->OneDimensionalOutput(solver_container[val_iZone][MESH_0][FLOW_SOL],
                                 geometry_container[val_iZone][MESH_0], config_container[val_iZone]);
  
  RegisterOutput(solver_container, geometry_container, config_container, val_iZone);
  
  /*--- Stop the recording ---*/
  
  AD::StopRecording();
  
  /*--- Set the recording status ---*/
  
  CurrentRecording = kind_recording;

  /* --- Reset the number of the internal iterations---*/

  config_container[ZONE_0]->SetIntIter(IntIter);

}


void CDiscAdjMeanFlowIteration::RegisterInput(CSolver ****solver_container, CGeometry ***geometry_container, CConfig **config_container, unsigned short iZone, unsigned short kind_recording){
  
  
  if (kind_recording == FLOW_VARIABLES || kind_recording == ALL_VARIABLES){
    
    /*--- Register flow and turbulent variables as input ---*/
    
    solver_container[iZone][MESH_0][ADJFLOW_SOL]->RegisterSolution(geometry_container[iZone][MESH_0], config_container[iZone]);
    
    solver_container[iZone][MESH_0][ADJFLOW_SOL]->RegisterVariables(geometry_container[iZone][MESH_0], config_container[iZone]);
    
    if (turbulent){
      solver_container[iZone][MESH_0][ADJTURB_SOL]->RegisterSolution(geometry_container[iZone][MESH_0], config_container[iZone]);
    }
  }
  if ((kind_recording == GEOMETRY_VARIABLES) || (kind_recording == ALL_VARIABLES)){
    
    /*--- Register node coordinates as input ---*/
    
    geometry_container[iZone][MESH_0]->RegisterCoordinates(config_container[iZone]);
    
  }

}

void CDiscAdjMeanFlowIteration::SetDependencies(CSolver ****solver_container, CGeometry ***geometry_container, CConfig **config_container, unsigned short iZone, unsigned short kind_recording){


  if ((kind_recording == GEOMETRY_VARIABLES) || (kind_recording == ALL_VARIABLES)){

    /*--- Update geometry to get the influence on other geometry variables (normals, volume etc) ---*/
    
    geometry_container[iZone][MESH_0]->UpdateGeometry(geometry_container[iZone], config_container[iZone]);
    
  }

  /*--- Compute coupling between flow and turbulent equations ---*/

  if (turbulent){
    solver_container[iZone][MESH_0][FLOW_SOL]->Preprocessing(geometry_container[iZone][MESH_0],solver_container[iZone][MESH_0], config_container[iZone], MESH_0, NO_RK_ITER, RUNTIME_FLOW_SYS, true);
    solver_container[iZone][MESH_0][TURB_SOL]->Postprocessing(geometry_container[iZone][MESH_0],solver_container[iZone][MESH_0], config_container[iZone], MESH_0);
  }

}

void CDiscAdjMeanFlowIteration::RegisterOutput(CSolver ****solver_container, CGeometry ***geometry_container, CConfig **config_container, unsigned short iZone){
  
  /*--- Register objective function as output of the iteration ---*/
  
  solver_container[iZone][MESH_0][ADJFLOW_SOL]->RegisterObj_Func(config_container[iZone]);
  
  /*--- Register conservative variables as output of the iteration ---*/
  
  solver_container[iZone][MESH_0][ADJFLOW_SOL]->RegisterOutput(geometry_container[iZone][MESH_0],config_container[iZone]);
  
  if (turbulent){
    solver_container[iZone][MESH_0][ADJTURB_SOL]->RegisterOutput(geometry_container[iZone][MESH_0],
                                                                 config_container[iZone]);
  }
}

void CDiscAdjMeanFlowIteration::InitializeAdjoint(CSolver ****solver_container, CGeometry ***geometry_container, CConfig **config_container, unsigned short iZone){
  
  /*--- Initialize the adjoint of the objective function (typically with 1.0) ---*/
  
  solver_container[iZone][MESH_0][ADJFLOW_SOL]->SetAdj_ObjFunc(geometry_container[iZone][MESH_0], config_container[iZone],1.0);
  
  /*--- Initialize the adjoints the conservative variables ---*/
  
  solver_container[iZone][MESH_0][ADJFLOW_SOL]->SetAdjoint_Output(geometry_container[iZone][MESH_0],
                                                                  config_container[iZone]);
  
  if (turbulent){
    solver_container[iZone][MESH_0][ADJTURB_SOL]->SetAdjoint_Output(geometry_container[iZone][MESH_0],
                                                                    config_container[iZone]);
  }
}
void CDiscAdjMeanFlowIteration::Update(COutput *output,
                                       CIntegration ***integration_container,
                                       CGeometry ***geometry_container,
                                       CSolver ****solver_container,
                                       CNumerics *****numerics_container,
                                       CConfig **config_container,
                                       CSurfaceMovement **surface_movement,
                                       CVolumetricMovement **grid_movement,
                                       CFreeFormDefBox*** FFDBox,
                                       unsigned short val_iZone)      { }
void CDiscAdjMeanFlowIteration::Monitor()     { }
void CDiscAdjMeanFlowIteration::Output()      { }
void CDiscAdjMeanFlowIteration::Postprocess() { }

void FEM_StructuralIteration(COutput *output, CIntegration ***integration_container, CGeometry ***geometry_container,
                  	  	  	  	 CSolver ****solver_container, CNumerics *****numerics_container, CConfig **config_container,
                  	  	  	  	 CSurfaceMovement **surface_movement, CVolumetricMovement **grid_movement, CFreeFormDefBox*** FFDBox) {

	su2double Physical_dt, Physical_t;
	su2double loadIncrement;
	unsigned short iZone;
	unsigned short nZone = geometry_container[ZONE_0][MESH_0]->GetnZone();
	unsigned long IntIter = 0; config_container[ZONE_0]->SetIntIter(IntIter);
  	unsigned long ExtIter = config_container[ZONE_0]->GetExtIter();

  	unsigned long iIncrement;
  	unsigned long nIncrements = config_container[ZONE_0]->GetNumberIncrements();

	bool dynamic = (config_container[ZONE_0]->GetDynamic_Analysis() == DYNAMIC);					// Dynamic problems
	bool nonlinear = (config_container[ZONE_0]->GetGeometricConditions() == LARGE_DEFORMATIONS);	// Geometrically non-linear problems

	bool incremental_load = config_container[ZONE_0]->GetIncrementalLoad();							// If an incremental load is applied

	/*--- This is to prevent problems when running a linear solver ---*/
	if (!nonlinear) incremental_load = false;

	int rank = MASTER_NODE;
#ifdef HAVE_MPI
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif


	/*--- THIS IS THE DIRECT APPROACH (NO INCREMENTAL LOAD APPLIED) ---*/

	if (!incremental_load){

		/*--- Set the initial condition ---*/

//		for (iZone = 0; iZone < nZone; iZone++)
//			solver_container[iZone][MESH_0][FEA_SOL]->SetInitialCondition(geometry_container[iZone], solver_container[iZone], config_container[iZone], ExtIter);

		for (iZone = 0; iZone < nZone; iZone++) {

			/*--- Set the value of the internal iteration ---*/

			IntIter = ExtIter;
			if (nonlinear) IntIter = 0;

			/*--- FEA equations ---*/

			config_container[iZone]->SetGlobalParam(FEM_ELASTICITY, RUNTIME_FEA_SYS, ExtIter);

			/*--- Run the iteration ---*/

			integration_container[iZone][FEA_SOL]->Structural_Iteration(geometry_container, solver_container, numerics_container,
	                                                                		config_container, RUNTIME_FEA_SYS, IntIter, iZone);



		}

		/*----------------- If the solver is non-linear, we need to subiterate using a Newton-Raphson approach ----------------------*/

		if (nonlinear){
			for (IntIter = 1; IntIter < config_container[ZONE_0]->GetDyn_nIntIter(); IntIter++){

				for (iZone = 0; iZone < nZone; iZone++) {

					/*--- Write the convergence history (only screen output) ---*/

					output->SetConvHistory_Body(NULL, geometry_container, solver_container, config_container, integration_container, true, 0.0, ZONE_0);

					config_container[iZone]->SetIntIter(IntIter);

					integration_container[iZone][FEA_SOL]->Structural_Iteration(geometry_container, solver_container, numerics_container,
			                                                                		config_container, RUNTIME_FEA_SYS, IntIter, iZone);

				}

				if (integration_container[ZONE_0][FEA_SOL]->GetConvergence()) break;

			}

		}

	}
	/*--- The incremental load is only used in nonlinear cases ---*/
	else if (incremental_load){

		/*--- Set the initial condition: store the current solution as Solution_Old ---*/

		for (iZone = 0; iZone < nZone; iZone++)
			solver_container[iZone][MESH_0][FEA_SOL]->SetInitialCondition(geometry_container[iZone], solver_container[iZone], config_container[iZone], ExtIter);

		for (iZone = 0; iZone < nZone; iZone++) {

				/*--- The load increment is 1.0 ---*/
				loadIncrement = 1.0;
				solver_container[iZone][MESH_0][FEA_SOL]->SetLoad_Increment(loadIncrement);

				/*--- Set the value of the internal iteration ---*/

				IntIter = 0;

				/*--- FEA equations ---*/

				config_container[iZone]->SetGlobalParam(FEM_ELASTICITY, RUNTIME_FEA_SYS, ExtIter);

				/*--- Run the first iteration ---*/

				integration_container[iZone][FEA_SOL]->Structural_Iteration(geometry_container, solver_container, numerics_container,
		                                                                		config_container, RUNTIME_FEA_SYS, IntIter, iZone);


				/*--- Write the convergence history (only screen output) ---*/

				output->SetConvHistory_Body(NULL, geometry_container, solver_container, config_container, integration_container, true, 0.0, ZONE_0);

				/*--- Run the second iteration ---*/

				IntIter = 1;

				config_container[iZone]->SetIntIter(IntIter);

				integration_container[iZone][FEA_SOL]->Structural_Iteration(geometry_container, solver_container, numerics_container,
		                                                                		config_container, RUNTIME_FEA_SYS, IntIter, iZone);

		}

		bool meetCriteria;
		su2double Residual_UTOL, Residual_RTOL, Residual_ETOL;
		su2double Criteria_UTOL, Criteria_RTOL, Criteria_ETOL;

		Criteria_UTOL = config_container[ZONE_0]->GetIncLoad_Criteria(0);
		Criteria_RTOL = config_container[ZONE_0]->GetIncLoad_Criteria(1);
		Criteria_ETOL = config_container[ZONE_0]->GetIncLoad_Criteria(2);

		Residual_UTOL = log10(solver_container[ZONE_0][MESH_0][FEA_SOL]->GetRes_FEM(0));
		Residual_RTOL = log10(solver_container[ZONE_0][MESH_0][FEA_SOL]->GetRes_FEM(1));
		Residual_ETOL = log10(solver_container[ZONE_0][MESH_0][FEA_SOL]->GetRes_FEM(2));

		meetCriteria = ( ( Residual_UTOL <  Criteria_UTOL ) &&
				 	 	 ( Residual_RTOL <  Criteria_RTOL ) &&
						 ( Residual_ETOL <  Criteria_ETOL ) );

		/*--- If the criteria is met and the load is not "too big", do the regular calculation ---*/
		if (meetCriteria){

			for (IntIter = 2; IntIter < config_container[ZONE_0]->GetDyn_nIntIter(); IntIter++){

				for (iZone = 0; iZone < nZone; iZone++) {

				/*--- Write the convergence history (only screen output) ---*/

				output->SetConvHistory_Body(NULL, geometry_container, solver_container, config_container, integration_container, true, 0.0, ZONE_0);

				config_container[iZone]->SetIntIter(IntIter);

				integration_container[iZone][FEA_SOL]->Structural_Iteration(geometry_container, solver_container, numerics_container,
																				config_container, RUNTIME_FEA_SYS, IntIter, iZone);

				}

				if (integration_container[ZONE_0][FEA_SOL]->GetConvergence()) break;

			}

		}

		/*--- If the criteria is not met, a whole set of subiterations for the different loads must be done ---*/

		else {

			/*--- Here we have to restart the solution to the original one of the iteration ---*/
			/*--- Retrieve the Solution_Old as the current solution before subiterating ---*/

			for (iZone = 0; iZone < nZone; iZone++)
				solver_container[iZone][MESH_0][FEA_SOL]->ResetInitialCondition(geometry_container[iZone], solver_container[iZone], config_container[iZone], ExtIter);

			/*--- For the number of increments ---*/
			for (iIncrement = 0; iIncrement < nIncrements; iIncrement++){

				loadIncrement = (iIncrement + 1.0) * (1.0 / nIncrements);

				/*--- Set the load increment and the initial condition, and output the parameters of UTOL, RTOL, ETOL for the previous iteration ---*/

				for (iZone = 0; iZone < nZone; iZone++){

					/*--- Set the convergence monitor to false, to force se solver to converge every subiteration ---*/
					integration_container[iZone][FEA_SOL]->SetConvergence(false);

					output->SetConvHistory_Body(NULL, geometry_container, solver_container, config_container, integration_container, true, 0.0, ZONE_0);

					/*--- FEA equations ---*/

					config_container[iZone]->SetGlobalParam(FEM_ELASTICITY, RUNTIME_FEA_SYS, ExtIter);


					solver_container[iZone][MESH_0][FEA_SOL]->SetLoad_Increment(loadIncrement);
				}

				if (rank == MASTER_NODE){
					cout << endl;
					cout << "-- Incremental load: increment " << iIncrement + 1 << " ------------------------------------------" << endl;
				}

				for (iZone = 0; iZone < nZone; iZone++) {

					/*--- Set the value of the internal iteration ---*/
					IntIter = 0;
					config_container[iZone]->SetIntIter(IntIter);

					/*--- FEA equations ---*/

					config_container[iZone]->SetGlobalParam(FEM_ELASTICITY, RUNTIME_FEA_SYS, ExtIter);

					/*--- Run the iteration ---*/

					integration_container[iZone][FEA_SOL]->Structural_Iteration(geometry_container, solver_container, numerics_container,
			                                                                		config_container, RUNTIME_FEA_SYS, IntIter, iZone);



				}

				/*----------------- If the solver is non-linear, we need to subiterate using a Newton-Raphson approach ----------------------*/

				for (IntIter = 1; IntIter < config_container[ZONE_0]->GetDyn_nIntIter(); IntIter++){

					for (iZone = 0; iZone < nZone; iZone++) {

						/*--- Write the convergence history (only screen output) ---*/

						output->SetConvHistory_Body(NULL, geometry_container, solver_container, config_container, integration_container, true, 0.0, ZONE_0);

						config_container[iZone]->SetIntIter(IntIter);

						integration_container[iZone][FEA_SOL]->Structural_Iteration(geometry_container, solver_container, numerics_container,
					                                                                		config_container, RUNTIME_FEA_SYS, IntIter, iZone);

					}

					if (integration_container[ZONE_0][FEA_SOL]->GetConvergence()) break;

				}

			}

		}

	}



	/*----------------- Compute averaged nodal stress and reactions ------------------------*/

	for (iZone = 0; iZone < nZone; iZone++)
		solver_container[iZone][MESH_0][FEA_SOL]->Compute_NodalStress(geometry_container[iZone][MESH_0], solver_container[iZone][MESH_0], numerics_container[iZone][MESH_0][FEA_SOL], config_container[iZone]);

	/*----------------- Update structural solver ----------------------*/

	if (dynamic){
		for (iZone = 0; iZone < nZone; iZone++) {
			integration_container[iZone][FEA_SOL]->SetFEM_StructuralSolver(geometry_container[iZone][MESH_0], solver_container[iZone][MESH_0], config_container[iZone], MESH_0);
			integration_container[iZone][FEA_SOL]->SetConvergence(false);
		}

	    /*--- Verify convergence criteria (based on total time) ---*/

		Physical_dt = config_container[ZONE_0]->GetDelta_DynTime();
		Physical_t  = (ExtIter+1)*Physical_dt;
		if (Physical_t >=  config_container[ZONE_0]->GetTotal_DynTime())
			integration_container[ZONE_0][FEA_SOL]->SetConvergence(true);
	}


}

void SetGrid_Movement(CGeometry **geometry_container, CSurfaceMovement *surface_movement,
                      CVolumetricMovement *grid_movement, CFreeFormDefBox **FFDBox,
                      CSolver ***solver_container, CConfig *config_container,
                      unsigned short iZone, unsigned long IntIter, unsigned long ExtIter)   {
  
  unsigned short iDim, iMGlevel, nMGlevels = config_container->GetnMGLevels();
  unsigned short Kind_Grid_Movement = config_container->GetKind_GridMovement(iZone);
  unsigned long nIterMesh;
  unsigned long iPoint;
  bool stat_mesh = true;
  bool adjoint = config_container->GetContinuous_Adjoint();
  bool time_spectral = (config_container->GetUnsteady_Simulation() == TIME_SPECTRAL);
  
  /*--- For a time-spectral case, set "iteration number" to the zone number,
   so that the meshes are positioned correctly for each instance. ---*/
  if (time_spectral) {
    ExtIter = iZone;
    Kind_Grid_Movement = config_container->GetKind_GridMovement(ZONE_0);
  }
  
  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  
  /*--- Perform mesh movement depending on specified type ---*/
  switch (Kind_Grid_Movement) {
      
    case MOVING_WALL:
      
      /*--- Fixed wall velocities: set the grid velocities only one time
       before the first iteration flow solver. ---*/
      
      if (ExtIter == 0) {
        
        if (rank == MASTER_NODE)
          cout << endl << " Setting the moving wall velocities." << endl;
        
        surface_movement->Moving_Walls(geometry_container[MESH_0],
                                       config_container, iZone, ExtIter);
        
        /*--- Update the grid velocities on the coarser multigrid levels after
         setting the moving wall velocities for the finest mesh. ---*/
        
        grid_movement->UpdateMultiGrid(geometry_container, config_container);
        
      }
      
      break;
      
      
    case ROTATING_FRAME:
      
      /*--- Steadily rotating frame: set the grid velocities just once
       before the first iteration flow solver. ---*/
      
      if (ExtIter == 0) {
        
        if (rank == MASTER_NODE) {
          cout << endl << " Setting rotating frame grid velocities";
          cout << " for zone " << iZone << "." << endl;
        }
        
        /*--- Set the grid velocities on all multigrid levels for a steadily
         rotating reference frame. ---*/
        
        for (iMGlevel = 0; iMGlevel <= nMGlevels; iMGlevel++)
          geometry_container[iMGlevel]->SetRotationalVelocity(config_container, iZone);
        
      }
      
      break;
      
    case STEADY_TRANSLATION:
      
      /*--- Set the translational velocity and hold the grid fixed during
       the calculation (similar to rotating frame, but there is no extra
       source term for translation). ---*/
      
      if (ExtIter == 0) {
        
        if (rank == MASTER_NODE)
          cout << endl << " Setting translational grid velocities." << endl;
        
        /*--- Set the translational velocity on all grid levels. ---*/
        
        for (iMGlevel = 0; iMGlevel <= nMGlevels; iMGlevel++)
          geometry_container[iMGlevel]->SetTranslationalVelocity(config_container);
        
      }
      
      break;
      
    case RIGID_MOTION:
      
      if (rank == MASTER_NODE) {
        cout << endl << " Performing rigid mesh transformation." << endl;
      }
      
      /*--- Move each node in the volume mesh using the specified type
       of rigid mesh motion. These routines also compute analytic grid
       velocities for the fine mesh. ---*/
      
      grid_movement->Rigid_Translation(geometry_container[MESH_0],
                                       config_container, iZone, ExtIter);
      grid_movement->Rigid_Plunging(geometry_container[MESH_0],
                                    config_container, iZone, ExtIter);
      grid_movement->Rigid_Pitching(geometry_container[MESH_0],
                                    config_container, iZone, ExtIter);
      grid_movement->Rigid_Rotation(geometry_container[MESH_0],
                                    config_container, iZone, ExtIter);
      
      /*--- Update the multigrid structure after moving the finest grid,
       including computing the grid velocities on the coarser levels. ---*/
      
      grid_movement->UpdateMultiGrid(geometry_container, config_container);
      
      break;
      
    case DEFORMING:
      
      if (rank == MASTER_NODE)
        cout << endl << " Updating surface positions." << endl;
      
      /*--- Translating ---*/
      
      /*--- Compute the new node locations for moving markers ---*/
      
      surface_movement->Surface_Translating(geometry_container[MESH_0],
                                            config_container, ExtIter, iZone);
      /*--- Deform the volume grid around the new boundary locations ---*/
      
      if (rank == MASTER_NODE)
        cout << " Deforming the volume grid." << endl;
      grid_movement->SetVolume_Deformation(geometry_container[MESH_0],
                                           config_container, true);
      
      /*--- Plunging ---*/
      
      /*--- Compute the new node locations for moving markers ---*/
      
      surface_movement->Surface_Plunging(geometry_container[MESH_0],
                                         config_container, ExtIter, iZone);
      /*--- Deform the volume grid around the new boundary locations ---*/
      
      if (rank == MASTER_NODE)
        cout << " Deforming the volume grid." << endl;
      grid_movement->SetVolume_Deformation(geometry_container[MESH_0],
                                           config_container, true);
      
      /*--- Pitching ---*/
      
      /*--- Compute the new node locations for moving markers ---*/
      
      surface_movement->Surface_Pitching(geometry_container[MESH_0],
                                         config_container, ExtIter, iZone);
      /*--- Deform the volume grid around the new boundary locations ---*/
      
      if (rank == MASTER_NODE)
        cout << " Deforming the volume grid." << endl;
      grid_movement->SetVolume_Deformation(geometry_container[MESH_0],
                                           config_container, true);
      
      /*--- Rotating ---*/
      
      /*--- Compute the new node locations for moving markers ---*/
      
      surface_movement->Surface_Rotating(geometry_container[MESH_0],
                                         config_container, ExtIter, iZone);
      /*--- Deform the volume grid around the new boundary locations ---*/
      
      if (rank == MASTER_NODE)
        cout << " Deforming the volume grid." << endl;
      grid_movement->SetVolume_Deformation(geometry_container[MESH_0],
                                           config_container, true);
      
      /*--- Update the grid velocities on the fine mesh using finite
       differencing based on node coordinates at previous times. ---*/
      
      if (!adjoint) {
        if (rank == MASTER_NODE)
          cout << " Computing grid velocities by finite differencing." << endl;
        geometry_container[MESH_0]->SetGridVelocity(config_container, ExtIter);
      }
      
      /*--- Update the multigrid structure after moving the finest grid,
       including computing the grid velocities on the coarser levels. ---*/
      
      grid_movement->UpdateMultiGrid(geometry_container, config_container);
      
      break;
      
    case EXTERNAL: case EXTERNAL_ROTATION:
      
      /*--- Apply rigid rotation to entire grid first, if necessary ---*/
      
      if (Kind_Grid_Movement == EXTERNAL_ROTATION) {
        if (rank == MASTER_NODE)
          cout << " Updating node locations by rigid rotation." << endl;
        grid_movement->Rigid_Rotation(geometry_container[MESH_0],
                                      config_container, iZone, ExtIter);
      }
      
      /*--- Load new surface node locations from external files ---*/
      
      if (rank == MASTER_NODE)
        cout << " Updating surface locations from file." << endl;
      surface_movement->SetExternal_Deformation(geometry_container[MESH_0],
                                                config_container, iZone, ExtIter);
      
      /*--- Deform the volume grid around the new boundary locations ---*/
      
      if (rank == MASTER_NODE)
        cout << " Deforming the volume grid." << endl;
      grid_movement->SetVolume_Deformation(geometry_container[MESH_0],
                                           config_container, true);
      
      /*--- Update the grid velocities on the fine mesh using finite
       differencing based on node coordinates at previous times. ---*/
      
      if (!adjoint) {
        if (rank == MASTER_NODE)
          cout << " Computing grid velocities by finite differencing." << endl;
        geometry_container[MESH_0]->SetGridVelocity(config_container, ExtIter);
      }
      
      /*--- Update the multigrid structure after moving the finest grid,
       including computing the grid velocities on the coarser levels. ---*/
      
      grid_movement->UpdateMultiGrid(geometry_container, config_container);
      
      break;
      
    case AEROELASTIC: case AEROELASTIC_RIGID_MOTION:
      
      /*--- Apply rigid mesh transformation to entire grid first, if necessary ---*/
      if (IntIter == 0) {
        if (Kind_Grid_Movement == AEROELASTIC_RIGID_MOTION) {
          
          if (rank == MASTER_NODE) {
            cout << endl << " Performing rigid mesh transformation." << endl;
          }
          
          /*--- Move each node in the volume mesh using the specified type
           of rigid mesh motion. These routines also compute analytic grid
           velocities for the fine mesh. ---*/
          
          grid_movement->Rigid_Translation(geometry_container[MESH_0],
                                           config_container, iZone, ExtIter);
          grid_movement->Rigid_Plunging(geometry_container[MESH_0],
                                        config_container, iZone, ExtIter);
          grid_movement->Rigid_Pitching(geometry_container[MESH_0],
                                        config_container, iZone, ExtIter);
          grid_movement->Rigid_Rotation(geometry_container[MESH_0],
                                        config_container, iZone, ExtIter);
          
          /*--- Update the multigrid structure after moving the finest grid,
           including computing the grid velocities on the coarser levels. ---*/
          
          grid_movement->UpdateMultiGrid(geometry_container, config_container);
        }
        
      }
      
      /*--- Use the if statement to move the grid only at selected dual time step iterations. ---*/
      else if (IntIter % config_container->GetAeroelasticIter() ==0) {
        
        if (rank == MASTER_NODE)
          cout << endl << " Solving aeroelastic equations and updating surface positions." << endl;
        
        /*--- Solve the aeroelastic equations for the new node locations of the moving markers(surfaces) ---*/
        
        solver_container[MESH_0][FLOW_SOL]->Aeroelastic(surface_movement, geometry_container[MESH_0], config_container, ExtIter);
        
        /*--- Deform the volume grid around the new boundary locations ---*/
        
        if (rank == MASTER_NODE)
          cout << " Deforming the volume grid due to the aeroelastic movement." << endl;
        grid_movement->SetVolume_Deformation(geometry_container[MESH_0],
                                             config_container, true);
        
        /*--- Update the grid velocities on the fine mesh using finite
         differencing based on node coordinates at previous times. ---*/
        
        if (rank == MASTER_NODE)
          cout << " Computing grid velocities by finite differencing." << endl;
        geometry_container[MESH_0]->SetGridVelocity(config_container, ExtIter);
        
        /*--- Update the multigrid structure after moving the finest grid,
         including computing the grid velocities on the coarser levels. ---*/
        
        grid_movement->UpdateMultiGrid(geometry_container, config_container);
      }
      
      break;
      
    case ELASTICITY:
      
      if (ExtIter != 0) {
        
        if (rank == MASTER_NODE)
          cout << " Deforming the grid using the Linear Elasticity solution." << endl;
        
        /*--- Update the coordinates of the grid using the linear elasticity solution. ---*/
        for (iPoint = 0; iPoint < geometry_container[MESH_0]->GetnPoint(); iPoint++) {
          
          su2double *U_time_nM1 = solver_container[MESH_0][FEA_SOL]->node[iPoint]->GetSolution_time_n1();
          su2double *U_time_n   = solver_container[MESH_0][FEA_SOL]->node[iPoint]->GetSolution_time_n();
          
          for (iDim = 0; iDim < geometry_container[MESH_0]->GetnDim(); iDim++)
            geometry_container[MESH_0]->node[iPoint]->AddCoord(iDim, U_time_n[iDim] - U_time_nM1[iDim]);
          
        }
        
      }
      
      break;
      
    case FLUID_STRUCTURE:

      if (rank == MASTER_NODE)
        cout << endl << "Deforming the grid for Fluid-Structure Interaction applications." << endl;

      /*--- Deform the volume grid around the new boundary locations ---*/

      if (rank == MASTER_NODE)
        cout << "Deforming the volume grid." << endl;
      grid_movement->SetVolume_Deformation(geometry_container[MESH_0],
                                           config_container, true);

      nIterMesh = grid_movement->Get_nIterMesh();
      stat_mesh = (nIterMesh == 0);

      if (!adjoint && !stat_mesh) {
        if (rank == MASTER_NODE)
          cout << "Computing grid velocities by finite differencing." << endl;
        geometry_container[MESH_0]->SetGridVelocity(config_container, ExtIter);
      }
      else if (stat_mesh){
          if (rank == MASTER_NODE)
            cout << "The mesh is up-to-date. Using previously stored grid velocities." << endl;
      }

      /*--- Update the multigrid structure after moving the finest grid,
       including computing the grid velocities on the coarser levels. ---*/

      grid_movement->UpdateMultiGrid(geometry_container, config_container);

      break;

    case NO_MOVEMENT: case GUST: default:
      
      /*--- There is no mesh motion specified for this zone. ---*/
      if (rank == MASTER_NODE)
        cout << "No mesh motion specified." << endl;
      
      break;
  }
  
}

CDiscOneShotIteration::CDiscOneShotIteration(CConfig *config) : CIteration(config), CurrentRecording(NONE){

  meanflow_iteration = new CMeanFlowIteration(config);

  turbulent = config->GetKind_Solver() == DISC_ADJ_RANS;
  (SolverControl.RefsolutionN)=18;
  SolverControl.ObservationN=2;
  (SolverControl.maxorderQuad)=7;
  //SolverControl.maxorderQuad=20;
  SolverControl.OutputN=6;
  SolverControl.formulaquad=4;
  SolverControl.stdnoise=1.0;
  SolverControl.paradim=3;
  SolverControl.sizet=1;
  //SolverControl.ndgl=2;
  SolverControl.ndgl=1;
 // SolverControl.errortolquad=1E-04;
  SolverControl.errortolquad=1E-1;
  SolverControl.flagdimadap=0;

}

CDiscOneShotIteration::~CDiscOneShotIteration(void) { }
void CDiscOneShotIteration::Preprocess(COutput *output,
                                           CIntegration ***integration_container,
                                           CGeometry ***geometry_container,
                                           CSolver ****solver_container,
                                           CNumerics *****numerics_container,
                                           CConfig **config_container,
                                           CSurfaceMovement **surface_movement,
                                           CVolumetricMovement **grid_movement,
                                           CFreeFormDefBox*** FFDBox,
                                           unsigned short val_iZone) { }

void CDiscOneShotIteration::Iterate(COutput *output,
                                        CIntegration ***integration_container,
                                        CGeometry ***geometry_container,
                                        CSolver ****solver_container,
                                        CNumerics *****numerics_container,
                                        CConfig **config_container,
                                        CSurfaceMovement **surface_movement,
                                        CVolumetricMovement **volume_grid_movement,
                                        CFreeFormDefBox*** FFDBox,
                                        unsigned short val_iZone) {
    DiscAdjMeanFlowIterationConstraint(output, integration_container, geometry_container, solver_container, numerics_container, config_container,
                              surface_movement, volume_grid_movement, FFDBox);
}

void CDiscOneShotIteration::Update(COutput *output,
                                       CIntegration ***integration_container,
                                       CGeometry ***geometry_container,
                                       CSolver ****solver_container,
                                       CNumerics *****numerics_container,
                                       CConfig **config_container,
                                       CSurfaceMovement **surface_movement,
                                       CVolumetricMovement **grid_movement,
                                       CFreeFormDefBox*** FFDBox,
                                       unsigned short val_iZone)      { }
void CDiscOneShotIteration::Monitor()     { }
void CDiscOneShotIteration::Output()      { }
void CDiscOneShotIteration::Postprocess() { }

void CDiscOneShotIteration::DiscAdjMeanFlowIterationRun(COutput *output, CIntegration ***integration_container, CGeometry ***geometry_container,
                          CSolver ****solver_container, CNumerics *****numerics_container, CConfig **config_container,
                          CSurfaceMovement **surface_movement, CVolumetricMovement **volume_grid_movement, CFreeFormDefBox*** FFDBox) {

  unsigned short iZone, iMesh, ExtIter = config_container[ZONE_0]->GetExtIter();
  unsigned short nZone = config_container[ZONE_0]->GetnZone();

  bool turbulent = false, flow = false;

  switch(config_container[ZONE_0]->GetKind_Solver()){
    case DISC_ADJ_EULER: case DISC_ADJ_NAVIER_STOKES: flow = true; break;
    case DISC_ADJ_RANS: flow = true; turbulent = true; break;
  }

  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  /*--- store the old solution ---*/

  for (iZone = 0; iZone < nZone; iZone++){
          solver_container[iZone][MESH_0][ADJFLOW_SOL]->StoreOldSolution();
  }

  /*--- Deform the Grid according to Design Variable ---*/
  su2double steplen=config_container[0]->GetOSStepSize();
  unsigned short whilecounter=0;
  do{
  if ( whilecounter >0 ){
        #if defined CODI_REVERSE_TYPE
         if (config_container[ZONE_0]->GetOne_Shot()) AD::globalTape.reset();
        #endif

        /*  for (iZone = 0; iZone < nZone; iZone++){
            solver_container[iZone][MESH_0][ADJFLOW_SOL]->DesignMinus();
          }
          deformOneShot(output, geometry_container, config_container, surface_movement, volume_grid_movement, solver_container);
        */

        if(config_container[ZONE_0]->GetArmijo()){

            steplen=steplen*0.5;
        }else{
            //quadratic interpolation in first step, then cubic interpolation
            if(whilecounter==1){
                steplen=solver_container[ZONE_0][MESH_0][ADJFLOW_SOL]->QuadraticApproximation(steplen);
            }else{
                steplen=solver_container[ZONE_0][MESH_0][ADJFLOW_SOL]->CubicApproximation(steplen);
            }
        }
  }
  else{
      if(ExtIter>config_container[0]->GetOneShotStart()){
      if ((solver_container[ZONE_0][MESH_0][ADJFLOW_SOL]->CheckDescentDirection(steplen)==false)&&(config_container[ZONE_0]->GetCheckDescent()==true)){
          std::cout<<"No Descent Direction!"<<std::endl;
          steplen=-1E-10;
          //steplen=-1E-15;
          solver_container[ZONE_0][MESH_0][ADJFLOW_SOL]->ChangeDirection();
      }
      }
  }
  whilecounter=whilecounter+1;

  solver_container[ZONE_0][MESH_0][ADJFLOW_SOL]->LoadOldSolution();

  if(ExtIter>config_container[0]->GetOneShotStart()){
      for (iZone = 0; iZone < nZone; iZone++){
        if(config_container[0]->GetBoundsExact()) solver_container[iZone][MESH_0][ADJFLOW_SOL]->DesignUpdateProjected(geometry_container[iZone][MESH_0],config_container[iZone], ExtIter, steplen);
        else steplen=solver_container[iZone][MESH_0][ADJFLOW_SOL]->DesignUpdateBounds(geometry_container[iZone][MESH_0],config_container[iZone], ExtIter, steplen);
        if(whilecounter>1) solver_container[iZone][MESH_0][ADJFLOW_SOL]->DesignMinus();
        if(config_container[0]->GetOneShotConstraint()==true && whilecounter==1) solver_container[iZone][MESH_0][ADJFLOW_SOL]->UpdateMultiplier(config_container[0]);
      }
      deformOneShot(output, geometry_container, config_container, surface_movement, volume_grid_movement, solver_container);
  }

  if (ExtIter == 0 || config_container[ZONE_0]->GetOne_Shot()){
      AD::Reset();
      for (iMesh = 0; iMesh <= config_container[ZONE_0]->GetnMGLevels(); iMesh++){
        solver_container[ZONE_0][iMesh][ADJFLOW_SOL]->SetRecording(geometry_container[ZONE_0][MESH_0], config_container[ZONE_0], ALL_VARIABLES);
      }

    /*--- Start the recording of all operations ---*/

    AD::StartRecording();

    /*--- Register all necessary variables on the tape ---*/
    for (iZone = 0; iZone < nZone; iZone++){

      /*--- Register the node coordinates as input of the iteration ---*/

      geometry_container[iZone][MESH_0]->RegisterCoordinates(config_container[iZone]);

      if (flow){
          solver_container[iZone][MESH_0][ADJFLOW_SOL]->RegisterSolution(geometry_container[iZone][MESH_0],
                                                                     config_container[iZone]);
      }
    }

    /*--- Update the mesh structure to get the influence of node coordinates ---*/

    for (iZone = 0; iZone < nZone; iZone++){
      geometry_container[iZone][MESH_0]->UpdateGeometry(geometry_container[iZone],config_container[iZone]);
    }

 /*   if (rank == MASTER_NODE && !config_container[ZONE_0]->GetOne_Shot()){
      cout << "Begin direct solver to store computational graph (single iteration)." << endl;
    }*/

    /*--- One iteration of the flow solver ---*/

    if (flow){
      /*MeanFlowIteration(output, integration_container, geometry_container,
                      solver_container, numerics_container, config_container,
                      surface_movement, volume_grid_movement, FFDBox);*/
      for (iZone = 0; iZone < nZone; iZone++){
        meanflow_iteration->Iterate(output,integration_container,geometry_container,solver_container,numerics_container,
                                  config_container,surface_movement,volume_grid_movement,FFDBox, iZone);
      }
    }


    if (rank == MASTER_NODE){// && !config_container[ZONE_0]->GetOne_Shot()){
      if(flow){
          cout << "log10[RMS Density]: " << log10(solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetRes_RMS(0))<<std::endl;
          cout     <<"Drag: "<< solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetTotal_CDrag()
              <<", Lift: "<< solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetTotal_CLift() << "." << endl;
      }
    }

    for (iZone = 0; iZone < nZone; iZone++){
      /*--- Register objective function as output of the iteration ---*/
      if (flow){
        solver_container[iZone][MESH_0][ADJFLOW_SOL]->RegisterObj_Func(config_container[iZone]);
        if (config_container[0]->GetOneShotConstraint()==true) solver_container[iZone][MESH_0][ADJFLOW_SOL]->RegisterConstraint_Func(config_container[iZone], geometry_container[iZone][MESH_0]);
         /*--- Register conservative variables as output of the iteration ---*/
        solver_container[iZone][MESH_0][ADJFLOW_SOL]->RegisterOutput(geometry_container[iZone][MESH_0],
                                                                    config_container[iZone]);
      }
    }

    /*--- Stop the recording ---*/

    AD::StopRecording();

  }

  if(whilecounter==1 && config_container[0]->GetOneShotConstraint()==true){
  if(config_container[0]->GetMultiplierNorm()){

  su2double normgu, normfu;
  double* setConstraint=new double[3];
  setConstraint[0]=1.0;
  setConstraint[1]=0.0;
  setConstraint[2]=0.0;
  //Compute g_u
  for (iZone = 0; iZone < nZone; iZone++){
    solver_container[iZone][MESH_0][ADJFLOW_SOL]->LoadOldAdjoint();

  }
  for (iZone = 0; iZone < nZone; iZone++){
    if (flow){
      solver_container[iZone][MESH_0][ADJFLOW_SOL]->SetAdj_ObjFunc(geometry_container[iZone][MESH_0], config_container[iZone],0.0);
      solver_container[iZone][MESH_0][ADJFLOW_SOL]->SetAdj_ConstraintFunc(geometry_container[iZone][MESH_0], config_container[iZone],setConstraint);

      solver_container[iZone][MESH_0][ADJFLOW_SOL]->SetAdjointOutputZero(geometry_container[iZone][MESH_0],
                                                                      config_container[iZone]);
      }
  }


  AD::ComputeAdjoint();

  for (iZone = 0; iZone < nZone; iZone++){
    if (flow){
      solver_container[iZone][MESH_0][ADJFLOW_SOL]->SetSensitivity(geometry_container[iZone][MESH_0],config_container[iZone]);
      normgu=solver_container[iZone][MESH_0][ADJFLOW_SOL]->SensitivityNorm(geometry_container[iZone][MESH_0]);;
    }
  }

  AD::ClearAdjoints();
  setConstraint[0]=0.0;

  //compute f_u
  for (iZone = 0; iZone < nZone; iZone++){
      solver_container[iZone][MESH_0][ADJFLOW_SOL]->LoadOldAdjoint();

    }

    for (iZone = 0; iZone < nZone; iZone++){
      if (flow){
        solver_container[iZone][MESH_0][ADJFLOW_SOL]->SetAdj_ObjFunc(geometry_container[iZone][MESH_0], config_container[iZone],1.0);
        solver_container[iZone][MESH_0][ADJFLOW_SOL]->SetAdj_ConstraintFunc(geometry_container[iZone][MESH_0], config_container[iZone],setConstraint);
        solver_container[iZone][MESH_0][ADJFLOW_SOL]->SetAdjointOutputZero(geometry_container[iZone][MESH_0],
                                                                        config_container[iZone]);
      }
    }

    delete [] setConstraint;

    AD::ComputeAdjoint();

    for (iZone = 0; iZone < nZone; iZone++){
      if (flow){
        solver_container[iZone][MESH_0][ADJFLOW_SOL]->SetSensitivity(geometry_container[iZone][MESH_0],config_container[iZone]);
        normfu=solver_container[iZone][MESH_0][ADJFLOW_SOL]->SensitivityNorm(geometry_container[iZone][MESH_0]);
      }
    }

    AD::ClearAdjoints();

  std::cout<<"Norms: "<<normfu<<" "<<normgu<<std::endl;
  std::cout<<"Factor(Normfu/normgu): "<<normfu/normgu<<std::endl;

  for (iZone = 0; iZone < nZone; iZone++){
    solver_container[iZone][MESH_0][ADJFLOW_SOL]->LoadOldAdjoint();

  }
  double* multiplierVec=new double[3];
  multiplierVec[0]=SU2_TYPE::GetValue(normfu/normgu);
  multiplierVec[1]=0.0;
  multiplierVec[2]=0.0;
  solver_container[0][MESH_0][ADJFLOW_SOL]->SetMultiplier(config_container[0], multiplierVec);

  if(SU2_TYPE::GetValue(solver_container[0][MESH_0][ADJFLOW_SOL]->GetConstraintFunc_Value()[0])<0){
      if(config_container[0]->GetEqualConstraint(0)){
          solver_container[0][MESH_0][ADJFLOW_SOL]->SetMultiplier(config_container[0], multiplierVec);
      }else{
          multiplierVec[0]=0.0;
          solver_container[0][MESH_0][ADJFLOW_SOL]->SetMultiplier(config_container[0], multiplierVec);
      }
  }
  delete [] multiplierVec;
  }else{
      double* multiplierVec=new double[3];
      multiplierVec[0]=SU2_TYPE::GetValue(config_container[ZONE_0]->GetConstraintStart());
      multiplierVec[1]=0.0;
      multiplierVec[2]=0.0;
      if(ExtIter==0) solver_container[0][MESH_0][ADJFLOW_SOL]->SetMultiplier(config_container[0], multiplierVec);
      delete [] multiplierVec;
  }

  }

  //------------------COMPUTE N_u-------------------------
  for (iZone = 0; iZone < nZone; iZone++){

    /*--- Initialize the adjoint of the objective function (typically with 1.0) ---*/
    if (flow){
      solver_container[iZone][MESH_0][ADJFLOW_SOL]->SetAdj_ObjFunc(geometry_container[iZone][MESH_0], config_container[iZone],1.0);
      if(config_container[0]->GetOneShotConstraint()==true) solver_container[iZone][MESH_0][ADJFLOW_SOL]->SetAdj_ConstraintFunc(geometry_container[iZone][MESH_0], config_container[iZone],solver_container[iZone][MESH_0][ADJFLOW_SOL]->GetMultiplier());
      /*--- Initialize the adjoints the conservative variables ---*/
      solver_container[iZone][MESH_0][ADJFLOW_SOL]->SetAdjoint_Output(geometry_container[iZone][MESH_0],
                                                                      config_container[iZone]);
    }
  }

  /*--- Run the adjoint computation ---*/

  AD::ComputeAdjoint();

  for (iZone = 0; iZone < nZone; iZone++){

    /*--- Set the sensitivities by extracting the adjoints of the node coordinates ---*/
    /*--- N_u is only needed when a design update is wanted ---*/
    if (flow){
      solver_container[iZone][MESH_0][ADJFLOW_SOL]->ResetSensitivity(geometry_container[iZone][MESH_0]);
      solver_container[iZone][MESH_0][ADJFLOW_SOL]->SetSensitivity(geometry_container[iZone][MESH_0],config_container[iZone]); //contains N_u
      solver_container[iZone][MESH_0][ADJFLOW_SOL]->SaveSurfaceSensitivity(geometry_container[iZone][MESH_0]);
      solver_container[iZone][MESH_0][ADJFLOW_SOL]->UpdateLagrangeSensitivity(geometry_container[iZone][MESH_0],1.0);
      /*--- Extract the adjoints of the conservative input variables and store them for the next iteration ---*/
      solver_container[iZone][MESH_0][ADJFLOW_SOL]->ExtractAdjoint_Solution(geometry_container[iZone][MESH_0],
                                                                     config_container[iZone]);
    }
  }



  AD::ClearAdjoints();

  //Calculate Augmented Lagrangian L^a=alpha/2*||dy||^2+beta/2*||dybar||^2+f+ybar^T dy
  for (iZone = 0; iZone < nZone; iZone++){
      solver_container[iZone][MESH_0][ADJFLOW_SOL]->AssembleLagrangian(config_container[iZone]);
  }
  }
  while(ExtIter>config_container[0]->GetOneShotStart()&&solver_container[ZONE_0][MESH_0][ADJFLOW_SOL]->CheckFirstWolfe(steplen)&&(whilecounter<config_container[0]->GetSearchCounterMax())&&config_container[0]->GetLineSearch()&&steplen>1E-15);

  //-----BFGS Update (only when needed) -----------------------

  if(ExtIter>=config_container[0]->GetOneShotStart()){

      for (iZone = 0; iZone < nZone; iZone++){
        solver_container[iZone][MESH_0][ADJFLOW_SOL]->ApplyDesignVar();
        solver_container[iZone][MESH_0][ADJFLOW_SOL]->StoreSaveSolution();
        if(config_container[0]->GetOneShotConstraint()==true) solver_container[iZone][MESH_0][ADJFLOW_SOL]->StoreConstraint(config_container[0]);
      }

  //------------------COMPUTE alpha*Deltay^TG_u-------------------------

  //load old adjoint solution
  for (iZone = 0; iZone < nZone; iZone++){
    solver_container[iZone][MESH_0][ADJFLOW_SOL]->LoadOldAdjoint();

  }

  for (iZone = 0; iZone < nZone; iZone++){
    if (flow){
      solver_container[iZone][MESH_0][ADJFLOW_SOL]->SetAdj_ObjFunc(geometry_container[iZone][MESH_0], config_container[iZone],0.0);
      //if(SU2_TYPE::GetPrimary(solver_container[0][MESH_0][ADJFLOW_SOL]->GetConstraintFunc_Value()[0])>0){
      if(config_container[0]->GetOneShotConstraint()==true) solver_container[iZone][MESH_0][ADJFLOW_SOL]->SetAdj_ConstraintFuncAD(geometry_container[iZone][MESH_0], config_container[iZone],solver_container[iZone][MESH_0][ADJFLOW_SOL]->GetConstraintFunc_Value());
     /* }else{
        solver_container[iZone][MESH_0][ADJFLOW_SOL]->SetAdj_ConstraintFunc(geometry_container[iZone][MESH_0], config_container[iZone],0.0);
      }*/
    }
    solver_container[iZone][MESH_0][ADJFLOW_SOL]->SetAdjointOutputUpdate(geometry_container[iZone][MESH_0],config_container[iZone]);
  }

  AD::ComputeAdjoint();

  for (iZone = 0; iZone < nZone; iZone++){
    if (flow){
      solver_container[iZone][MESH_0][ADJFLOW_SOL]->SetSensitivity(geometry_container[iZone][MESH_0],config_container[iZone]); //contains Deltay^T G_u (Attention: overwrites Nu!)
      solver_container[iZone][MESH_0][ADJFLOW_SOL]->UpdateLagrangeSensitivity(geometry_container[iZone][MESH_0],config_container[ZONE_0]->GetOneShotAlpha());
    }
  }

  AD::ClearAdjoints();

  //Compute beta*Deltaybar^TNyu using Finite Differences

  for (iZone = 0; iZone < nZone; iZone++){
    solver_container[iZone][MESH_0][ADJFLOW_SOL]->LoadSaveSolution();
  }

  for (iZone = 0; iZone < nZone; iZone++){
      solver_container[iZone][MESH_0][ADJFLOW_SOL]->UpdateStateVariable(config_container[ZONE_0]);
  }
  cout << "log10Ajdoint[RMS Density]: " << log10(solver_container[ZONE_0][MESH_0][ADJFLOW_SOL]->GetRes_RMS(0))<<std::endl;
  cout << "log10Primal[RMS Density]: " << log10(solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetRes_RMS(0))<<std::endl;

  /*--- Clear all adjoints to re-use the stored computational graph in the next iteration ---*/


  /*--- reset the AD tape ---*/
#if defined CODI_REVERSE_TYPE
if (config_container[ZONE_0]->GetOne_Shot()) AD::globalTape.reset();
#endif

  /*--- Set the convergence criteria ---*/

  if (config_container[ZONE_0]->GetConvCriteria() == RESIDUAL){
    if (flow){
      integration_container[ZONE_0][ADJFLOW_SOL]->Convergence_Monitoring(geometry_container[ZONE_0][MESH_0],config_container[ZONE_0],
                                                                         ExtIter,log10(solver_container[ZONE_0][MESH_0][ADJFLOW_SOL]->GetRes_RMS(0)), MESH_0);
    }
  }
  else if (config_container[ZONE_0]->GetConvCriteria() == CAUCHY){
    if (flow){
      integration_container[ZONE_0][ADJFLOW_SOL]->Convergence_Monitoring(geometry_container[ZONE_0][MESH_0],config_container[ZONE_0],
                                                                         ExtIter,solver_container[ZONE_0][MESH_0][ADJFLOW_SOL]->GetTotal_Sens_Geo(), MESH_0);
    }
  }

  //Restore old adjoint Solution (Delta is stored in Solution)
  for (iZone = 0; iZone < nZone; iZone++){
    solver_container[iZone][MESH_0][ADJFLOW_SOL]->LoadOldAdjoint();
  }

  if (ExtIter == 0 || config_container[ZONE_0]->GetOne_Shot()){

    /*--- Start the recording of all operations ---*/
      AD::Reset();
      for (iMesh = 0; iMesh <= config_container[ZONE_0]->GetnMGLevels(); iMesh++){
        solver_container[ZONE_0][iMesh][ADJFLOW_SOL]->SetRecording(geometry_container[ZONE_0][MESH_0], config_container[ZONE_0], ALL_VARIABLES);
      }

    AD::StartRecording();

    /*--- Register all necessary variables on the tape ---*/
    for (iZone = 0; iZone < nZone; iZone++){

      /*--- Register the node coordinates as input of the iteration ---*/

      geometry_container[iZone][MESH_0]->RegisterCoordinates(config_container[iZone]);
      //getDVValue, SetDVValue

      if (flow){
          /*--- Register the conservative variables as input of the iteration ---*/
        solver_container[iZone][MESH_0][ADJFLOW_SOL]->RegisterSolution(geometry_container[iZone][MESH_0],
                                                                   config_container[iZone]);
      }

    }

    /*--- Update the mesh structure to get the influence of node coordinates ---*/

    for (iZone = 0; iZone < nZone; iZone++){
      geometry_container[iZone][MESH_0]->UpdateGeometry(geometry_container[iZone],config_container[iZone]);
    }

    if (rank == MASTER_NODE && !config_container[ZONE_0]->GetOne_Shot()){
      cout << "Begin direct solver to store computational graph (single iteration)." << endl;
    }

    /*--- One iteration of the flow solver ---*/

    if (flow){
        for (iZone = 0; iZone < nZone; iZone++){
          meanflow_iteration->Iterate(output,integration_container,geometry_container,solver_container,numerics_container,
                                    config_container,surface_movement,volume_grid_movement,FFDBox, iZone);
        }
    }

    if (rank == MASTER_NODE && !config_container[ZONE_0]->GetOne_Shot()){
      if(flow){
        cout << "log10[RMS Density]: " << log10(solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetRes_RMS(0))
             <<", Drag: "<< solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetTotal_CDrag()
            <<", Lift: "<< solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetTotal_CLift() << "." << endl;
      }
    }

    for (iZone = 0; iZone < nZone; iZone++){
      /*--- Register objective function as output of the iteration ---*/
      if (flow){
        solver_container[iZone][MESH_0][ADJFLOW_SOL]->RegisterObj_Func(config_container[iZone]);
        if(config_container[0]->GetOneShotConstraint()==true) solver_container[iZone][MESH_0][ADJFLOW_SOL]->RegisterConstraint_Func(config_container[iZone], geometry_container[iZone][MESH_0]);
        /*--- Register conservative variables as output of the iteration ---*/
        solver_container[iZone][MESH_0][ADJFLOW_SOL]->RegisterOutput(geometry_container[iZone][MESH_0],
                                                                    config_container[iZone]);
      }
    }

    /*--- Stop the recording ---*/

    AD::StopRecording();

  }

  //------------------COMPUTE new N_u-------------------------
  for (iZone = 0; iZone < nZone; iZone++){

    /*--- Initialize the adjoint of the objective function (typically with 1.0) ---*/
    if (flow){
      solver_container[iZone][MESH_0][ADJFLOW_SOL]->SetAdj_ObjFunc(geometry_container[iZone][MESH_0], config_container[iZone],1.0);
      if(config_container[0]->GetOneShotConstraint()==true) solver_container[iZone][MESH_0][ADJFLOW_SOL]->SetAdj_ConstraintFunc(geometry_container[iZone][MESH_0], config_container[iZone],solver_container[iZone][MESH_0][ADJFLOW_SOL]->GetMultiplier());
      /*--- Initialize the adjoints the conservative variables ---*/
      solver_container[iZone][MESH_0][ADJFLOW_SOL]->SetAdjoint_Output(geometry_container[iZone][MESH_0],
                                                                    config_container[iZone]);
    }
  }

  /*--- Run the adjoint computation ---*/

  AD::ComputeAdjoint();

  for (iZone = 0; iZone < nZone; iZone++){

    /*--- Set the sensitivities by extracting the adjoints of the node coordinates ---*/
    if (flow){
      solver_container[iZone][MESH_0][ADJFLOW_SOL]->SetSensitivityFD(geometry_container[iZone][MESH_0],config_container[iZone]); //contains N_u_new
      solver_container[iZone][MESH_0][ADJFLOW_SOL]->UpdateLagrangeSensitivity(geometry_container[iZone][MESH_0],config_container[ZONE_0]->GetOneShotBeta());
    }
  }


  for (iZone = 0; iZone < nZone; iZone++){
    /*--- Project Sensitivities ---*/
    solver_container[iZone][MESH_0][ADJFLOW_SOL]->OverwriteSensitivityProjected(geometry_container[iZone][MESH_0]);
    projectOneShot(geometry_container, config_container, surface_movement, volume_grid_movement, solver_container);

    solver_container[iZone][MESH_0][ADJFLOW_SOL]->OverwriteGradientProjected(geometry_container[iZone][MESH_0]);
    projectGradient(geometry_container, config_container, surface_movement, volume_grid_movement, solver_container);

    /*--- BFGS Update and Calculation of Descent Direction ---*/
    std::cout<<"steplen: "<<steplen<<std::endl;
    std::cout<<"searchsteps: "<<whilecounter<<std::endl;
    solver_container[iZone][MESH_0][ADJFLOW_SOL]->BFGSUpdateProjected(geometry_container[iZone][MESH_0],config_container[iZone],ExtIter);
  }

  /*--- Clear all adjoints to re-use the stored computational graph in the next iteration ---*/
  AD::ClearAdjoints();

  for (iZone = 0; iZone < nZone; iZone++){
    solver_container[iZone][MESH_0][ADJFLOW_SOL]->LoadSaveSolution();
  }
  }

#if defined CODI_REVERSE_TYPE
if (config_container[ZONE_0]->GetOne_Shot()) AD::globalTape.reset();
#endif




}

void CDiscOneShotIteration::DiscAdjMeanFlowIterationSparseGrid(COutput *output, CIntegration ***integration_container, CGeometry ***geometry_container,
                                                               CSolver ****solver_container, CNumerics *****numerics_container, CConfig **config_container,
                                                               CSurfaceMovement **surface_movement, CVolumetricMovement **volume_grid_movement, CFreeFormDefBox*** FFDBox){
    unsigned short iZone, ExtIter = config_container[ZONE_0]->GetExtIter();
    unsigned short nZone = config_container[ZONE_0]->GetnZone();

    su2double steplen=config_container[ZONE_0]->GetOSStepSize();

    // parameter initilization
    //all data related to the solver control
    SolverControl.errortolquad=100.0/(ExtIter);

    oIndex solref;
    oIndex erroriter;
    oIndex indexPoints;
    oIndex addPoints;
    oIndex indexintPoints;
    oIndex numind;
    oIndex indexset;
    oIndex indexwi[SolverControl.ndgl];
    oIndex quadwIf;
    oIndex quadwIAf;
    oIndex tspan(1,SolverControl.sizet);
    oIndex solSG;
    for (int k=0;k<SolverControl.sizet;k++){
        //tspan.set(0,k,k*1.0/((double)intPow(2,6)));
        tspan.set(0,k,k);
    }

    solver_container[ZONE_0][MESH_0][ADJFLOW_SOL]->StoreOldSolution();

    unsigned short whilecounter=0;
    do{
        if ( whilecounter >0 ){
              if(config_container[ZONE_0]->GetArmijo()){
                  steplen=steplen*0.5;
              /*}else{
                  //quadratic interpolation in first step, then cubic interpolation
                  if(whilecounter==1){
                      steplen=solver_container[ZONE_0][MESH_0][ADJFLOW_SOL]->QuadraticApproximation(steplen);
                  }else{
                      steplen=solver_container[ZONE_0][MESH_0][ADJFLOW_SOL]->CubicApproximation(steplen);
                  }
              */
              }
        }
        /*else{
            if(ExtIter>config_container[0]->GetOneShotStart()){
                if ((solver_container[ZONE_0][MESH_0][ADJFLOW_SOL]->CheckDescentDirection(steplen)==false)&&(config_container[ZONE_0]->GetCheckDescent()==true)){
                    std::cout<<"No Descent Direction!"<<std::endl;
                    //steplen=-1E-10;
                    //steplen=-1E-15;
                    solver_container[ZONE_0][MESH_0][ADJFLOW_SOL]->ChangeDirection();
                }
            }
        }*/
        whilecounter=whilecounter+1;

        solver_container[ZONE_0][MESH_0][ADJFLOW_SOL]->LoadOldSolution();

        if(ExtIter>=config_container[0]->GetOneShotStart()){
            if(ExtIter>config_container[0]->GetOneShotStart()){
                for (iZone = 0; iZone < nZone; iZone++){
                  solver_container[iZone][MESH_0][ADJFLOW_SOL]->DesignUpdateProjected(geometry_container[iZone][MESH_0],config_container[iZone], ExtIter, steplen);
                  if(whilecounter>1) solver_container[iZone][MESH_0][ADJFLOW_SOL]->DesignMinus();
                  if(config_container[0]->GetOneShotConstraint()==true && whilecounter==1) solver_container[iZone][MESH_0][ADJFLOW_SOL]->UpdateMultiplier(config_container[0]);
                }
                deformOneShot(output, geometry_container, config_container, surface_movement, volume_grid_movement, solver_container);
            }

            solver_container[ZONE_0][MESH_0][ADJFLOW_SOL]->ResetStoredValues(geometry_container[ZONE_0][MESH_0],config_container[0]);

            solver_container[ZONE_0][MESH_0][ADJFLOW_SOL]->ComputeEVEW(config_container[ZONE_0], geometry_container[ZONE_0][MESH_0]);
            solver_container[ZONE_0][MESH_0][ADJFLOW_SOL]->StoreMeshPoints(config_container[ZONE_0], geometry_container[ZONE_0][MESH_0]);
            SolverControl.ndgl=1;
            ASparseGridSingle( &indexset, &indexPoints,indexwi,&indexintPoints,&numind, &solref,&quadwIf,&quadwIAf,&erroriter,&SolverControl,
                               output, integration_container, geometry_container,
                               solver_container, numerics_container, config_container,
                                surface_movement, volume_grid_movement, FFDBox, whilecounter, &addPoints, this);
//new
            solver_container[ZONE_0][MESH_0][ADJFLOW_SOL]->LoadOldSolution();
            solver_container[ZONE_0][MESH_0][ADJFLOW_SOL]->LoadMeshPoints(config_container[ZONE_0], geometry_container[ZONE_0][MESH_0]);
//new
            OneShotStep(output, integration_container, geometry_container,
               solver_container, numerics_container, config_container,
                surface_movement, volume_grid_movement, FFDBox, whilecounter);

                         char sBuffer[2048];
                         sprintf(sBuffer,"./results/erroriter.txt" );
                         erroriter.mfprintf(sBuffer);
/*                         sprintf(sBuffer, "./results/indexset.txt");
                         indexset.mfprintf(sBuffer);
                         sprintf(sBuffer, "./results/indexset.bin");
                         indexset.mfprintb(sBuffer);
                         sprintf(sBuffer,"./results/indexPoints.txt" );
                         indexPoints.mfprintf(sBuffer);
                         sprintf(sBuffer,"./results/indexPoints.bin" );
                         indexPoints.mfprintb(sBuffer);
                         sprintf(sBuffer, "./results/indexintPoints.txt");
                         indexintPoints.mfprintf(sBuffer);
                         sprintf(sBuffer, "./results/indexintPoints.bin");
                         indexintPoints.mfprintb(sBuffer);
                         sprintf(sBuffer, "./results/addPoints.txt");
                         addPoints.mfprintf(sBuffer);
                         sprintf(sBuffer,"./results/erroriter.txt" );
                         erroriter.mfprintf(sBuffer);
                         sprintf(sBuffer,"./results/numind.txt" );
                         numind.mfprintf(sBuffer);
                         sprintf(sBuffer,"./results/numind.bin" );
                         numind.mfprintb(sBuffer);

                         for (int k=0;k<SolverControl.ndgl;k++){
                             sprintf(sBuffer, "./results/indexwi_%d.txt", k);
                             indexwi[k].mfprintf(sBuffer);
                             sprintf(sBuffer, "./results/indexwi_%d.bin", k);
                             indexwi[k].mfprintb(sBuffer);
                         }
                         sprintf(sBuffer,"./results/quadwI.txt");
                         quadwIf.mfprintf(sBuffer);
                         sprintf(sBuffer,"./results/quadwIA.txt" );
                         quadwIAf.mfprintf(sBuffer);*/

            //OBJ_SAVE
            SolverControl.ndgl=1;
            oIndex indexwo[SolverControl.ndgl];
            for (int k=0;k<SolverControl.ndgl;k++){
                 indexwo[k].zeros(SolverControl.sizet,indexintPoints.n);
            }
            for (int j=0;j<indexintPoints.n;j++){
                bool found=false;
                int k=0;
                while(k<addPoints.n&&found==false)
                {
                    if(indexintPoints.GetColumn(j).isequal(addPoints.GetColumn(k))){
                        found=true;
                        oIndex tmpw((SolverControl).ndgl,(SolverControl).sizet);
                        for (int l=0;l<(SolverControl).ndgl;l++){
                            for (int m=0;m<(SolverControl).sizet;m++){
                                tmpw.set(l,m,SU2_TYPE::GetValue(solver_container[ZONE_0][MESH_0][ADJFLOW_SOL]->GetStoredObjective(geometry_container[ZONE_0][MESH_0],config_container[0],k)));
                                indexwo[l].set(m,j,tmpw.get(l,m));
                            }
                        }
                    }
                    k++;
                }
            }
            //oIndex solSG;
            EvalSGnonnested(&indexset,indexset.n, &numind,&indexPoints, indexwo,SolverControl.formulaquad , &tspan,SolverControl.ndgl, &solSG,SolverControl.maxorderQuad);
            for (int k=0;k<SolverControl.ndgl;k++){
                su2double adder=0;
                std::cout<<"Set Obj "<<std::endl;
                for (int l=0;l<solSG.n;l++){
                    adder=adder+solSG.get(k,l);
                    std::cout<<solSG.get(k,l)<<" ";
                }
                std::cout<<"->"<<adder<<std::endl;
                solver_container[ZONE_0][MESH_0][ADJFLOW_SOL]->SetObjSave(geometry_container[ZONE_0][MESH_0],config_container[0],adder);
            }
            //LAGRANGIAN
            SolverControl.ndgl=1;
            oIndex indexwl[SolverControl.ndgl];
            for (int k=0;k<SolverControl.ndgl;k++){
                 indexwl[k].zeros(SolverControl.sizet,indexintPoints.n);
             }

            for (int j=0;j<indexintPoints.n;j++){
                bool found=false;
                int k=0;
                while(k<addPoints.n&&found==false)
                {
                    if(indexintPoints.GetColumn(j).isequal(addPoints.GetColumn(k))){
                        found=true;
                        oIndex tmpw((SolverControl).ndgl,(SolverControl).sizet);
                        for (int l=0;l<(SolverControl).ndgl;l++){
                            for (int m=0;m<(SolverControl).sizet;m++){
                                tmpw.set(l,m,SU2_TYPE::GetValue(solver_container[ZONE_0][MESH_0][ADJFLOW_SOL]->GetStoredLagragian(geometry_container[ZONE_0][MESH_0],config_container[0],k)));
                                indexwl[l].set(m,j,tmpw.get(l,m));
                            }
                        }
                    }
                    k++;
                }
            }
            EvalSGnonnested(&indexset,indexset.n, &numind,&indexPoints, indexwl,SolverControl.formulaquad , &tspan,SolverControl.ndgl, &solSG,SolverControl.maxorderQuad);

            for (int k=0;k<SolverControl.ndgl;k++){
                su2double adder=0;
                for (int l=0;l<solSG.n;l++){
                    adder=adder+solSG.get(k,l);
                }
                solver_container[ZONE_0][MESH_0][ADJFLOW_SOL]->SetLagrangian(geometry_container[ZONE_0][MESH_0],config_container[0],adder);
            }
            //CONSTRAINTS
            SolverControl.ndgl=config_container[0]->GetConstraintNum();
            oIndex indexwc[SolverControl.ndgl];
            for (int k=0;k<SolverControl.ndgl;k++){
                 indexwc[k].zeros(SolverControl.sizet,indexintPoints.n);
             }

            for (int j=0;j<indexintPoints.n;j++){
                bool found=false;
                int k=0;
                while(k<addPoints.n&&found==false)
                {
                    if(indexintPoints.GetColumn(j).isequal(addPoints.GetColumn(k))){
                        found=true;
                        oIndex tmpw((SolverControl).ndgl,(SolverControl).sizet);
                        for (int l=0;l<(SolverControl).ndgl;l++){
                            for (int m=0;m<(SolverControl).sizet;m++){
                                tmpw.set(l,m,SU2_TYPE::GetValue(solver_container[ZONE_0][MESH_0][ADJFLOW_SOL]->GetStoredConstraints(geometry_container[ZONE_0][MESH_0],config_container[0],k,l)));
                                indexwc[l].set(m,j,tmpw.get(l,m));
                            }
                        }
                    }
                    k++;
                }
            }
            EvalSGnonnested(&indexset,indexset.n, &numind,&indexPoints, indexwc,SolverControl.formulaquad , &tspan,SolverControl.ndgl, &solSG,SolverControl.maxorderQuad);

            for (int k=0;k<SolverControl.ndgl;k++){
                su2double adder=0;
                for (int l=0;l<solSG.n;l++){
                    adder=adder+solSG.get(k,l);
                }
                solver_container[ZONE_0][MESH_0][ADJFLOW_SOL]->SetConstraintSave(geometry_container[ZONE_0][MESH_0],config_container[0],adder,k);
            }

            //CSENS
            SolverControl.ndgl=geometry_container[ZONE_0][MESH_0]->GetnVertex(0);
            oIndex indexwcs[SolverControl.ndgl];
            for (int k=0;k<SolverControl.ndgl;k++){
                 indexwcs[k].zeros(SolverControl.sizet,indexintPoints.n);
             }

            for (int j=0;j<indexintPoints.n;j++){
                bool found=false;
                int k=0;
                while(k<addPoints.n&&found==false)
                {
                    if(indexintPoints.GetColumn(j).isequal(addPoints.GetColumn(k))){
                        found=true;
                        oIndex tmpw((SolverControl).ndgl,(SolverControl).sizet);
                        for (int l=0;l<(SolverControl).ndgl;l++){
                            for (int m=0;m<(SolverControl).sizet;m++){
                                tmpw.set(l,m,SU2_TYPE::GetValue(solver_container[ZONE_0][MESH_0][ADJFLOW_SOL]->GetStoredCSens(geometry_container[ZONE_0][MESH_0],config_container[0],k,l)));
                                indexwcs[l].set(m,j,tmpw.get(l,m));
                            }
                        }
                    }
                    k++;
                }
            }
            EvalSGnonnested(&indexset,indexset.n, &numind,&indexPoints, indexwcs,SolverControl.formulaquad , &tspan,SolverControl.ndgl, &solSG,SolverControl.maxorderQuad);

            for (int k=0;k<SolverControl.ndgl;k++){
                su2double adder=0;
                for (int l=0;l<solSG.n;l++){
                    adder=adder+solSG.get(k,l);
                }
                solver_container[ZONE_0][MESH_0][ADJFLOW_SOL]->SetCSens(geometry_container[ZONE_0][MESH_0],config_container[0],adder,k);
            }
            //LAGRANGESENS
            SolverControl.ndgl=geometry_container[ZONE_0][MESH_0]->GetnVertex(0);
            oIndex indexwls[SolverControl.ndgl];
            for (int k=0;k<SolverControl.ndgl;k++){
                 indexwls[k].zeros(SolverControl.sizet,indexintPoints.n);
             }

            for (int j=0;j<indexintPoints.n;j++){
                bool found=false;
                int k=0;
                while(k<addPoints.n&&found==false)
                {
                    if(indexintPoints.GetColumn(j).isequal(addPoints.GetColumn(k))){
                        found=true;
                        oIndex tmpw((SolverControl).ndgl,(SolverControl).sizet);
                        for (int l=0;l<(SolverControl).ndgl;l++){
                            for (int m=0;m<(SolverControl).sizet;m++){
                                tmpw.set(l,m,SU2_TYPE::GetValue(solver_container[ZONE_0][MESH_0][ADJFLOW_SOL]->GetStoredLagrangeSens(geometry_container[ZONE_0][MESH_0],config_container[0],k,l)));
                                indexwls[l].set(m,j,tmpw.get(l,m));
                            }
                        }
                    }
                    k++;
                }
            }
            EvalSGnonnested(&indexset,indexset.n, &numind,&indexPoints, indexwls,SolverControl.formulaquad , &tspan,SolverControl.ndgl, &solSG,SolverControl.maxorderQuad);

            for (int k=0;k<SolverControl.ndgl;k++){
                su2double adder=0;
                for (int l=0;l<solSG.n;l++){
                    adder=adder+solSG.get(k,l);
                }
                solver_container[ZONE_0][MESH_0][ADJFLOW_SOL]->SetLagrangeSens(geometry_container[ZONE_0][MESH_0],config_container[0],adder,k);
            }
        }
        else{
            //solver_container[ZONE_0][MESH_0][ADJFLOW_SOL]->LoadOldSolution();
  /*          su2double rand[3] = {0.0,0.0,1.0};
            solver_container[ZONE_0][MESH_0][ADJFLOW_SOL]->ComputeEVEW(config_container[ZONE_0], geometry_container[ZONE_0][MESH_0]);
            solver_container[ZONE_0][MESH_0][ADJFLOW_SOL]->StoreMeshPoints(config_container[ZONE_0], geometry_container[ZONE_0][MESH_0]);
            solver_container[ZONE_0][MESH_0][ADJFLOW_SOL]->DeformSurface(config_container[ZONE_0], geometry_container[ZONE_0][MESH_0],rand,volume_grid_movement[ZONE_0]);

            OneShotStep(output, integration_container, geometry_container,
               solver_container, numerics_container, config_container,
                surface_movement, volume_grid_movement, FFDBox, whilecounter);

            solver_container[ZONE_0][MESH_0][ADJFLOW_SOL]->LoadOldSolution();
            solver_container[ZONE_0][MESH_0][ADJFLOW_SOL]->LoadMeshPoints(config_container[ZONE_0], geometry_container[ZONE_0][MESH_0]);
            solver_container[ZONE_0][MESH_0][ADJFLOW_SOL]->DeformSurface(config_container[ZONE_0], geometry_container[ZONE_0][MESH_0],rand,volume_grid_movement[ZONE_0]);*/
            //geometry_container[ZONE_0][MESH_0]->UpdateGeometry(geometry_container[ZONE_0],config_container[ZONE_0]);

  /*          solver_container[ZONE_0][MESH_0][ADJFLOW_SOL]->ComputeEVEW(config_container[ZONE_0], geometry_container[ZONE_0][MESH_0]);
            solver_container[ZONE_0][MESH_0][ADJFLOW_SOL]->StoreMeshPoints(config_container[ZONE_0], geometry_container[ZONE_0][MESH_0]);

            OneShotStepTrial(output, integration_container, geometry_container,
               solver_container, numerics_container, config_container,
                surface_movement, volume_grid_movement, FFDBox, whilecounter);

            solver_container[ZONE_0][MESH_0][ADJFLOW_SOL]->LoadOldSolution();
            solver_container[ZONE_0][MESH_0][ADJFLOW_SOL]->LoadMeshPoints(config_container[ZONE_0], geometry_container[ZONE_0][MESH_0]);*/

            OneShotStep(output, integration_container, geometry_container,
               solver_container, numerics_container, config_container,
                surface_movement, volume_grid_movement, FFDBox, whilecounter);
        }

    }
    while(ExtIter>config_container[0]->GetOneShotStart()&&solver_container[ZONE_0][MESH_0][ADJFLOW_SOL]->CheckFirstWolfe(steplen)&&(whilecounter<config_container[0]->GetSearchCounterMax())&&config_container[0]->GetLineSearch()&&steplen>1E-15);


    if(ExtIter>=config_container[ZONE_0]->GetOneShotStart()){
        for (iZone = 0; iZone < nZone; iZone++){
          /*--- Project Sensitivities ---*/
          solver_container[iZone][MESH_0][ADJFLOW_SOL]->OverwriteSensitivityProjected(geometry_container[iZone][MESH_0]);
          projectOneShot(geometry_container, config_container, surface_movement, volume_grid_movement, solver_container);

          solver_container[iZone][MESH_0][ADJFLOW_SOL]->OverwriteGradientProjected(geometry_container[iZone][MESH_0]);
          projectGradient(geometry_container, config_container, surface_movement, volume_grid_movement, solver_container);

          /*--- BFGS Update and Calculation of Descent Direction ---*/
          solver_container[iZone][MESH_0][ADJFLOW_SOL]->ApplyDesignVar();
          solver_container[iZone][MESH_0][ADJFLOW_SOL]->BFGSUpdateProjected(geometry_container[iZone][MESH_0],config_container[iZone],ExtIter);
        }
    }
}

void CDiscOneShotIteration::DiscAdjMeanFlowIterationPiggyBack(COutput *output, CIntegration ***integration_container, CGeometry ***geometry_container,
                          CSolver ****solver_container, CNumerics *****numerics_container, CConfig **config_container,
                          CSurfaceMovement **surface_movement, CVolumetricMovement **volume_grid_movement, CFreeFormDefBox*** FFDBox) {

  unsigned short iZone, iMesh, ExtIter = config_container[ZONE_0]->GetExtIter();
  unsigned short nZone = config_container[ZONE_0]->GetnZone();

  bool turbulent = false, flow = false;

  switch(config_container[ZONE_0]->GetKind_Solver()){
    case DISC_ADJ_EULER: case DISC_ADJ_NAVIER_STOKES: flow = true; break;
    case DISC_ADJ_RANS: flow = true; turbulent = true; break;
  }

  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  AD::Reset();

  if (ExtIter == 0 || config_container[ZONE_0]->GetOne_Shot()){

      for (iMesh = 0; iMesh <= config_container[ZONE_0]->GetnMGLevels(); iMesh++){
        solver_container[ZONE_0][iMesh][ADJFLOW_SOL]->SetRecording(geometry_container[ZONE_0][MESH_0], config_container[ZONE_0], ALL_VARIABLES);
      }

    /*--- Start the recording of all operations ---*/

    AD::StartRecording();

    /*--- Register all necessary variables on the tape ---*/
    for (iZone = 0; iZone < nZone; iZone++){

      /*--- Register the node coordinates as input of the iteration ---*/

      geometry_container[iZone][MESH_0]->RegisterCoordinates(config_container[iZone]);

      for (iMesh = 0; iMesh <= config_container[iZone]->GetnMGLevels(); iMesh++){

        /*--- Register the conservative variables as input of the iteration ---*/
        if (flow){
          solver_container[iZone][iMesh][ADJFLOW_SOL]->RegisterSolution(geometry_container[iZone][iMesh],
                                                                     config_container[iZone]);
        }
        if (turbulent){
          solver_container[iZone][iMesh][ADJTURB_SOL]->RegisterSolution(geometry_container[iZone][iMesh],
                                                                     config_container[iZone]);
        }
      }
    }

    /*--- Update the mesh structure to get the influence of node coordinates ---*/

    for (iZone = 0; iZone < nZone; iZone++){
      geometry_container[iZone][MESH_0]->UpdateGeometry(geometry_container[iZone],config_container[iZone]);
    }

    if (rank == MASTER_NODE && !config_container[ZONE_0]->GetOne_Shot()){
      cout << "Begin direct solver to store computational graph (single iteration)." << endl;
    }


    /* --- Preprocessing of flow solver and postprocessing of turbulent solver.
     * We need this to get the dependency of the flow variables on the eddy viscosity. ---*/
    if (turbulent){
      for (iZone = 0; iZone < nZone; iZone++){
        for (iMesh = 0; iMesh <= config_container[iZone]->GetnMGLevels(); iMesh++){
          solver_container[iZone][iMesh][FLOW_SOL]->Preprocessing(geometry_container[iZone][iMesh], solver_container[iZone][iMesh], config_container[iZone], iMesh, NO_RK_ITER, RUNTIME_FLOW_SYS, false);
          solver_container[iZone][iMesh][TURB_SOL]->Postprocessing(geometry_container[iZone][iMesh],solver_container[iZone][iMesh], config_container[iZone], iMesh);
        }
      }
    }

    /*--- One iteration of the flow solver ---*/

    if (flow){
        for (iZone = 0; iZone < nZone; iZone++){
          meanflow_iteration->Iterate(output,integration_container,geometry_container,solver_container,numerics_container,
                                    config_container,surface_movement,volume_grid_movement,FFDBox, iZone);
        }
    }

    if (rank == MASTER_NODE){// && !config_container[ZONE_0]->GetOne_Shot()){
      if(flow){
        cout << "log10[RMS Density]: " << log10(solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetRes_RMS(0))<<std::endl;
        cout     <<"Drag: "<< solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetTotal_CDrag()
            <<", Lift: "<< solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetTotal_CLift() << "." << endl;
      }
      if (turbulent){
        cout << "log10[RMS k]: " << log10(solver_container[ZONE_0][MESH_0][TURB_SOL]->GetRes_RMS(0)) << std::endl;
      }
    }

    for (iZone = 0; iZone < nZone; iZone++){
      /*--- Register objective function as output of the iteration ---*/
      if (flow){
        solver_container[iZone][MESH_0][ADJFLOW_SOL]->RegisterObj_Func(config_container[iZone]);
      }
      for (iMesh = 0; iMesh <= config_container[iZone]->GetnMGLevels(); iMesh++){

        /*--- Register conservative variables as output of the iteration ---*/
        if (flow){
          solver_container[iZone][iMesh][ADJFLOW_SOL]->RegisterOutput(geometry_container[iZone][iMesh],
                                                                      config_container[iZone]);
        }
        if (turbulent){
          solver_container[iZone][iMesh][ADJTURB_SOL]->RegisterOutput(geometry_container[iZone][iMesh],
                                                                      config_container[iZone]);
        }
      }
    }

    /*--- Stop the recording ---*/

    AD::StopRecording();

  }
  for (iZone = 0; iZone < nZone; iZone++){

    /*--- Initialize the adjoint of the objective function (typically with 1.0) ---*/
    if (flow){
      solver_container[iZone][MESH_0][ADJFLOW_SOL]->SetAdj_ObjFunc(geometry_container[iZone][MESH_0], config_container[iZone], 1.0);
    }
    for (iMesh = 0; iMesh <= config_container[iZone]->GetnMGLevels(); iMesh++){

      /*--- Initialize the adjoints the conservative variables ---*/
      if (flow){
        solver_container[iZone][iMesh][ADJFLOW_SOL]->SetAdjoint_Output(geometry_container[iZone][iMesh],
                                                                      config_container[iZone]);
      }
      if (turbulent){
        solver_container[iZone][iMesh][ADJTURB_SOL]->SetAdjoint_Output(geometry_container[iZone][iMesh],
                                                                      config_container[iZone]);
      }
    }
  }

  /*--- Run the adjoint computation ---*/

  AD::ComputeAdjoint();

  for (iZone = 0; iZone < nZone; iZone++){

    /*--- Set the sensitivities by extracting the adjoints of the node coordinates ---*/
    if (flow){
      solver_container[iZone][MESH_0][ADJFLOW_SOL]->SetSensitivity(geometry_container[iZone][MESH_0],config_container[iZone]);
    }
    for (iMesh = 0; iMesh <= config_container[iZone]->GetnMGLevels(); iMesh++){

      /*--- Extract the adjoints of the conservative input variables and store them for the next iteration ---*/
      if (flow){
        solver_container[iZone][iMesh][ADJFLOW_SOL]->ExtractAdjoint_Solution(geometry_container[iZone][iMesh],
                                                                     config_container[iZone]);
      }
      if (turbulent){
        solver_container[iZone][iMesh][ADJTURB_SOL]->ExtractAdjoint_Solution(geometry_container[iZone][iMesh],
                                                                     config_container[iZone]);
      }
    }
  }

  /*--- Clear all adjoints to re-use the stored computational graph in the next iteration ---*/

  solver_container[ZONE_0][MESH_0][ADJFLOW_SOL]->OutputWritten(geometry_container[ZONE_0][MESH_0]);

  AD::ClearAdjoints();

  cout << "log10Ajdoint[RMS Density]: " << log10(solver_container[ZONE_0][MESH_0][ADJFLOW_SOL]->GetRes_RMS(0))<<std::endl;
  cout << "log10Primal[RMS Density]: " << log10(solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetRes_RMS(0))<<std::endl;

#if defined CODI_REVERSE_TYPE
if (config_container[ZONE_0]->GetOne_Shot()) AD::globalTape.reset();
#endif

  /*--- Set the convergence criteria ---*/

  if (config_container[ZONE_0]->GetConvCriteria() == RESIDUAL){
    if (flow){
      integration_container[ZONE_0][ADJFLOW_SOL]->Convergence_Monitoring(geometry_container[ZONE_0][MESH_0],config_container[ZONE_0],
                                                                         ExtIter,log10(solver_container[ZONE_0][MESH_0][ADJFLOW_SOL]->GetRes_RMS(0)), MESH_0);
    }
  }
  else if (config_container[ZONE_0]->GetConvCriteria() == CAUCHY){
    if (flow){
      integration_container[ZONE_0][ADJFLOW_SOL]->Convergence_Monitoring(geometry_container[ZONE_0][MESH_0],config_container[ZONE_0],
                                                                         ExtIter,solver_container[ZONE_0][MESH_0][ADJFLOW_SOL]->GetTotal_Sens_Geo(), MESH_0);
    }
  }
}

void CDiscOneShotIteration::DiscAdjMeanFlowIterationConstraint(COutput *output, CIntegration ***integration_container, CGeometry ***geometry_container,
                          CSolver ****solver_container, CNumerics *****numerics_container, CConfig **config_container,
                          CSurfaceMovement **surface_movement, CVolumetricMovement **volume_grid_movement, CFreeFormDefBox*** FFDBox) {

    unsigned short iZone, iMesh, ExtIter = config_container[ZONE_0]->GetExtIter();
    unsigned short nZone = config_container[ZONE_0]->GetnZone();

    su2double steplen=config_container[0]->GetOSStepSize();

 //   su2double alpha_next=steplen;

    //print residuals!
       cout << "log10Primal[RMS Density]: " << log10(solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetRes_RMS(0))<<std::endl;
       cout << "log10Ajdoint[RMS Density]: " << log10(solver_container[ZONE_0][MESH_0][ADJFLOW_SOL]->GetRes_RMS(0))<<std::endl;

    for (iZone = 0; iZone < nZone; iZone++){
            solver_container[iZone][MESH_0][ADJFLOW_SOL]->StoreOldSolution();
    }
    solver_container[ZONE_0][MESH_0][ADJFLOW_SOL]->StoreMeshPoints(config_container[ZONE_0], geometry_container[ZONE_0][MESH_0]);

    unsigned short whilecounter=0;
    do{

        if ( whilecounter >0 ){
              if(config_container[ZONE_0]->GetArmijo()){
                  steplen=steplen*0.5;
              }else{
                  //quadratic interpolation in first step, then cubic interpolation
                  if(whilecounter==1){
                      steplen=solver_container[ZONE_0][MESH_0][ADJFLOW_SOL]->QuadraticApproximation(steplen);
                  }else{
                      steplen=solver_container[ZONE_0][MESH_0][ADJFLOW_SOL]->CubicApproximation(steplen);
                  }
              }
        }
        else{
            if(ExtIter>config_container[0]->GetOneShotStart()){
                if ((solver_container[ZONE_0][MESH_0][ADJFLOW_SOL]->CheckDescentDirection(steplen)==false)&&(config_container[ZONE_0]->GetCheckDescent()==true)){
                    std::cout<<"No Descent Direction!"<<std::endl;
                    //steplen=-1E-10;
                    //steplen=-1E-15;
                    solver_container[ZONE_0][MESH_0][ADJFLOW_SOL]->ChangeDirection();
                }
            }
        }
        whilecounter=whilecounter+1;

      //  for (iZone = 0; iZone < nZone; iZone++){
                solver_container[ZONE_0][MESH_0][ADJFLOW_SOL]->LoadOldSolution();
                solver_container[ZONE_0][MESH_0][ADJFLOW_SOL]->LoadMeshPoints(config_container[ZONE_0], geometry_container[ZONE_0][MESH_0]);
     //   }

        if(ExtIter>config_container[0]->GetOneShotStart()){
            for (iZone = 0; iZone < nZone; iZone++){
              solver_container[iZone][MESH_0][ADJFLOW_SOL]->DesignUpdateProjected(geometry_container[iZone][MESH_0],config_container[iZone], ExtIter, steplen);
              //if(whilecounter>1) solver_container[iZone][MESH_0][ADJFLOW_SOL]->DesignMinus();
              if(config_container[0]->GetOneShotConstraint()==true && whilecounter==1) solver_container[iZone][MESH_0][ADJFLOW_SOL]->UpdateMultiplier(config_container[0]);
            }
            deformOneShot(output, geometry_container, config_container, surface_movement, volume_grid_movement, solver_container);
        }

        OneShotStep(output, integration_container, geometry_container,
               solver_container, numerics_container, config_container,
                surface_movement, volume_grid_movement, FFDBox, whilecounter);


    }
    while(ExtIter>config_container[0]->GetOneShotStart()&&solver_container[ZONE_0][MESH_0][ADJFLOW_SOL]->CheckFirstWolfe(steplen)&&(whilecounter<config_container[0]->GetSearchCounterMax())&&config_container[0]->GetLineSearch()&&steplen>1E-15);

/*    for (iZone = 0; iZone < nZone; iZone++){
            solver_container[iZone][MESH_0][ADJFLOW_SOL]->StoreOldSolution();
    }*/

    if(ExtIter>=config_container[ZONE_0]->GetOneShotStart()){
        for (iZone = 0; iZone < nZone; iZone++){
          /*--- Project Sensitivities ---*/
          solver_container[iZone][MESH_0][ADJFLOW_SOL]->OverwriteSensitivityProjected(geometry_container[iZone][MESH_0]);
          projectOneShot(geometry_container, config_container, surface_movement, volume_grid_movement, solver_container);

          solver_container[iZone][MESH_0][ADJFLOW_SOL]->OverwriteGradientProjected(geometry_container[iZone][MESH_0]);
          projectGradient(geometry_container, config_container, surface_movement, volume_grid_movement, solver_container);

          /*--- BFGS Update and Calculation of Descent Direction ---*/
          solver_container[iZone][MESH_0][ADJFLOW_SOL]->ApplyDesignVar();
          solver_container[iZone][MESH_0][ADJFLOW_SOL]->BFGSUpdateProjected(geometry_container[iZone][MESH_0],config_container[iZone],ExtIter);
        }
    }

    /*--- Deform the Grid according to Design Variable ---*/
/*    su2double steplen=config_container[0]->GetOSStepSize();
    unsigned short whilecounter=0;

    bool linesearch=false;

    su2double ftol=1E-4;
    su2double gtol=0.9;

    su2double stpmin = 0.0;
    su2double stpmax = 4;
    su2double xtol= 1E-15;

    su2double xtrapu = 4.0;
    su2double xtrapl = 1.1;
    //or xtrapf=4.0

    bool    brackt = false; //
    bool    stage = true;
    su2double    width = stpmax - stpmin;
    su2double    width1= width/0.5;

    su2double alpha0=0.0;
    su2double phi_alpha0, dPhi_alpha0; //finit und ginit
    solver_container[ZONE_0][MESH_0][ADJFLOW_SOL]->CalculatePhi(alpha0, phi_alpha0, dPhi_alpha0);

    if(ExtIter>config_container[0]->GetOneShotStart()){
        linesearch=true;
        if(dPhi_alpha0>0) linesearch=false;
    }

    su2double alpha_lo = alpha0; //current best (x)
    su2double alpha_hi= alpha0; //current endpoint (y)
    su2double phi_lo, phi_hi, dPhi_lo, dPhi_hi;
    phi_lo = phi_hi= phi_alpha0;
    dPhi_lo = dPhi_hi = dPhi_alpha0;

    su2double stmin, stmax;
    su2double alpha_next=steplen;
    bool infoc=true;


    do{
        if(brackt){
            stmin=std::min(alpha_lo,alpha_hi);
            stmax=std::max(alpha_lo, alpha_hi);
        }else{
            stmin=alpha_lo;
            stmax=alpha_next+xtrapu*(alpha_next-alpha_lo);
        }
        alpha_next = std::max(alpha_next, stpmin);
        alpha_next = std::min(alpha_next, stpmax);
        if ((brackt && (alpha_next <= stmin || alpha_next >= stmax)) || infoc == false || (brackt && stmax-stmin <= xtol*stmax)) {
            alpha_next=alpha_lo;
        }
        if ( whilecounter >0 ){
              #if defined CODI_REVERSE_TYPE
               if (config_container[ZONE_0]->GetOne_Shot()) AD::globalTape.reset();
              #endif
        }
        //call function with alpha_next
        whilecounter=whilecounter+1;

        solver_container[ZONE_0][MESH_0][ADJFLOW_SOL]->LoadOldSolution();

        if(ExtIter>config_container[0]->GetOneShotStart()){
            for (iZone = 0; iZone < nZone; iZone++){
              solver_container[iZone][MESH_0][ADJFLOW_SOL]->DesignUpdateProjected(geometry_container[iZone][MESH_0],config_container[iZone], ExtIter, alpha_next);
              if(whilecounter>1) solver_container[iZone][MESH_0][ADJFLOW_SOL]->DesignMinus();
              if(config_container[0]->GetOneShotConstraint()==true && whilecounter==1) solver_container[iZone][MESH_0][ADJFLOW_SOL]->UpdateMultiplier(config_container[0]);
            }
            deformOneShot(output, geometry_container, config_container, surface_movement, volume_grid_movement, solver_container);
        }

 /*
    while(ExtIter>config_container[0]->GetOneShotStart()&&solver_container[ZONE_0][MESH_0][ADJFLOW_SOL]->CheckFirstWolfe(steplen)&&config_container[0]->GetLineSearch()&&steplen>1E-15);
*/




}

void CDiscOneShotIteration::DiscAdjMeanFlowIterationRobust(COutput *output, CIntegration ***integration_container, CGeometry ***geometry_container,
                          CSolver ****solver_container, CNumerics *****numerics_container, CConfig **config_container,
                          CSurfaceMovement **surface_movement, CVolumetricMovement **volume_grid_movement, CFreeFormDefBox*** FFDBox) {

    unsigned short iZone, ExtIter = config_container[ZONE_0]->GetExtIter();
    unsigned short nZone = config_container[ZONE_0]->GetnZone(), numQuad;

    su2double steplen=config_container[ZONE_0]->GetOSStepSize();
    su2double machp;

    //config_container[ZONE_0]->SetMach(0.8);

    for (iZone = 0; iZone < nZone; iZone++){
            if(ExtIter==0){
                for (numQuad=0; numQuad<4;numQuad++){
                    machp=solver_container[ZONE_0][MESH_0][ADJFLOW_SOL]->GetMachP(numQuad);
                    solver_container[iZone][MESH_0][FLOW_SOL]->SetInflow_Mach(geometry_container[iZone][MESH_0], config_container[iZone], MESH_0, machp, true, numQuad);
                    solver_container[iZone][MESH_0][ADJFLOW_SOL]->StoreSolutionVec(numQuad);
                }
            }
            for (numQuad=0; numQuad<4;numQuad++){
                solver_container[iZone][MESH_0][ADJFLOW_SOL]->StoreSolutionVecOld(numQuad);
            }
    }

    solver_container[ZONE_0][MESH_0][ADJFLOW_SOL]->StoreMeshPoints(config_container[ZONE_0], geometry_container[ZONE_0][MESH_0]);


    unsigned short whilecounter=0;
    do{
        if ( whilecounter >0 ){
              if(config_container[ZONE_0]->GetArmijo()){
                  steplen=steplen*0.5;
             }
             /*else{
                  //quadratic interpolation in first step, then cubic interpolation
                  if(whilecounter==1){
                      steplen=solver_container[ZONE_0][MESH_0][ADJFLOW_SOL]->QuadraticApproximation(steplen);
                  }else{
                      steplen=solver_container[ZONE_0][MESH_0][ADJFLOW_SOL]->CubicApproximation(steplen);
                  }
              }*/
        }
       /* else{
            if(ExtIter>config_container[0]->GetOneShotStart()){
                if ((solver_container[ZONE_0][MESH_0][ADJFLOW_SOL]->CheckDescentDirection(steplen)==false)&&(config_container[ZONE_0]->GetCheckDescent()==true)){
                    std::cout<<"No Descent Direction!"<<std::endl;;
                    solver_container[ZONE_0][MESH_0][ADJFLOW_SOL]->ChangeDirection();
                }
            }
        }*/
        whilecounter=whilecounter+1;
        solver_container[ZONE_0][MESH_0][ADJFLOW_SOL]->LoadMeshPoints(config_container[ZONE_0], geometry_container[ZONE_0][MESH_0]);

        if(ExtIter>config_container[ZONE_0]->GetOneShotStart()){
            for (iZone = 0; iZone < nZone; iZone++){
              solver_container[iZone][MESH_0][ADJFLOW_SOL]->DesignUpdateProjected(geometry_container[iZone][MESH_0],config_container[iZone], ExtIter, steplen);
              //if(whilecounter>1) solver_container[iZone][MESH_0][ADJFLOW_SOL]->DesignMinus();
              if(config_container[0]->GetOneShotConstraint()==true && whilecounter==1) solver_container[iZone][MESH_0][ADJFLOW_SOL]->UpdateMultiplier(config_container[0]);
            }
            deformOneShot(output, geometry_container, config_container, surface_movement, volume_grid_movement, solver_container);
        }

        solver_container[ZONE_0][MESH_0][ADJFLOW_SOL]->ResetExpValues(geometry_container[ZONE_0][MESH_0],config_container[0]);

        for (numQuad=0; numQuad<4;numQuad++){
              solver_container[ZONE_0][MESH_0][ADJFLOW_SOL]->LoadSolutionVecOld(numQuad);
              solver_container[ZONE_0][MESH_0][ADJFLOW_SOL]->StoreOldSolution();
              machp=solver_container[ZONE_0][MESH_0][ADJFLOW_SOL]->GetMachP(numQuad);
              solver_container[ZONE_0][MESH_0][FLOW_SOL]->SetInflow_Mach(geometry_container[ZONE_0][MESH_0], config_container[ZONE_0], MESH_0, machp, false, numQuad);
               OneShotStep(output, integration_container, geometry_container,
                    solver_container, numerics_container, config_container,
                    surface_movement, volume_grid_movement, FFDBox, whilecounter);
               solver_container[ZONE_0][MESH_0][ADJFLOW_SOL]->StoreSolutionVec(numQuad);
               solver_container[ZONE_0][MESH_0][ADJFLOW_SOL]->SumExpValues(geometry_container[ZONE_0][MESH_0],config_container[0],numQuad);
        }

        solver_container[ZONE_0][MESH_0][ADJFLOW_SOL]->DistributeExpValues(geometry_container[ZONE_0][MESH_0],config_container[0]);

    }
    while(ExtIter>config_container[0]->GetOneShotStart()&&solver_container[ZONE_0][MESH_0][ADJFLOW_SOL]->CheckFirstWolfe(steplen)&&(whilecounter<config_container[0]->GetSearchCounterMax())&&config_container[0]->GetLineSearch()&&steplen>1E-15);

    if(ExtIter>=config_container[ZONE_0]->GetOneShotStart()){
        for (iZone = 0; iZone < nZone; iZone++){
          /*--- Project Sensitivities ---*/
          solver_container[iZone][MESH_0][ADJFLOW_SOL]->OverwriteSensitivityProjected(geometry_container[iZone][MESH_0]);
          projectOneShot(geometry_container, config_container, surface_movement, volume_grid_movement, solver_container);

          solver_container[iZone][MESH_0][ADJFLOW_SOL]->OverwriteGradientProjected(geometry_container[iZone][MESH_0]);
          projectGradient(geometry_container, config_container, surface_movement, volume_grid_movement, solver_container);

          /*--- BFGS Update and Calculation of Descent Direction ---*/
          solver_container[iZone][MESH_0][ADJFLOW_SOL]->ApplyDesignVar();
          solver_container[iZone][MESH_0][ADJFLOW_SOL]->BFGSUpdateProjected(geometry_container[iZone][MESH_0],config_container[iZone],ExtIter);
        }
    }

}

void CDiscOneShotIteration::OneShotStep(COutput *output, CIntegration ***integration_container, CGeometry ***geometry_container,
                 CSolver ****solver_container, CNumerics *****numerics_container, CConfig **config_container,
                 CSurfaceMovement **surface_movement, CVolumetricMovement **volume_grid_movement, CFreeFormDefBox*** FFDBox, unsigned short whilecounter) {

    unsigned short iZone, iMesh, ExtIter = config_container[ZONE_0]->GetExtIter();
    unsigned short nZone = config_container[ZONE_0]->GetnZone();

    bool turbulent = false, flow = false;

    switch(config_container[ZONE_0]->GetKind_Solver()){
      case DISC_ADJ_EULER: case DISC_ADJ_NAVIER_STOKES: flow = true; break;
      case DISC_ADJ_RANS: flow = true; turbulent = true; break;
    }

 /*   for (iZone = 0; iZone < nZone; iZone++){
            solver_container[iZone][MESH_0][ADJFLOW_SOL]->StoreOldSolution();
    }*/
    AD::Reset();
    for (iMesh = 0; iMesh <= config_container[ZONE_0]->GetnMGLevels(); iMesh++){
      solver_container[ZONE_0][iMesh][ADJFLOW_SOL]->SetRecording(geometry_container[ZONE_0][MESH_0], config_container[ZONE_0], ALL_VARIABLES);
    }

      AD::StartRecording();

      for (iZone = 0; iZone < nZone; iZone++){
        geometry_container[iZone][MESH_0]->RegisterCoordinates(config_container[iZone]);
        if (flow){
            solver_container[iZone][MESH_0][ADJFLOW_SOL]->RegisterSolution(geometry_container[iZone][MESH_0],
                                                                       config_container[iZone]);
        }
      }
      for (iZone = 0; iZone < nZone; iZone++){
        geometry_container[iZone][MESH_0]->UpdateGeometry(geometry_container[iZone],config_container[iZone]);
       // solver_container[iZone][MESH_0][FLOW_SOL]->Preprocessing(geometry_container[iZone][MESH_0], solver_container[iZone][MESH_0], config_container[iZone], MESH_0, NO_RK_ITER, RUNTIME_FLOW_SYS, false);

      }
      if (flow){
          for (iZone = 0; iZone < nZone; iZone++){
            meanflow_iteration->Iterate(output,integration_container,geometry_container,solver_container,numerics_container,
                                      config_container,surface_movement,volume_grid_movement,FFDBox, iZone);
          }
      }
        if(flow){
       //     cout << "log10Primal[RMS Density]: " << log10(solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetRes_RMS(0))<<std::endl;
            cout     <<"Drag: "<< solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetTotal_CDrag()
                <<", Lift: "<< solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetTotal_CLift() << "." << endl;
        }
      for (iZone = 0; iZone < nZone; iZone++){
        if (flow){
          solver_container[iZone][MESH_0][ADJFLOW_SOL]->RegisterObj_Func(config_container[iZone]);
          if (config_container[0]->GetOneShotConstraint()==true) solver_container[iZone][MESH_0][ADJFLOW_SOL]->RegisterConstraint_Func(config_container[iZone], geometry_container[iZone][MESH_0]);
          solver_container[iZone][MESH_0][ADJFLOW_SOL]->RegisterOutput(geometry_container[iZone][MESH_0],
                                                                      config_container[iZone]);
        }
      }
      AD::StopRecording();

      double* multiplierVec = new double[3];
      multiplierVec[0]=SU2_TYPE::GetValue(config_container[ZONE_0]->GetConstraintStart());
      multiplierVec[1]=0.0;
      multiplierVec[2]=0.0;
      if(ExtIter==0 && whilecounter==1 && config_container[0]->GetOneShotConstraint()==true) solver_container[0][MESH_0][ADJFLOW_SOL]->SetMultiplier(config_container[0],multiplierVec );
       delete [] multiplierVec;
    //------------------COMPUTE N_u-------------------------
    for (iZone = 0; iZone < nZone; iZone++){
      if (flow){
        solver_container[iZone][MESH_0][ADJFLOW_SOL]->SetAdj_ObjFunc(geometry_container[iZone][MESH_0], config_container[iZone],1.0);
        if(config_container[0]->GetOneShotConstraint()==true) solver_container[iZone][MESH_0][ADJFLOW_SOL]->SetAdj_ConstraintFunc(geometry_container[iZone][MESH_0], config_container[iZone],solver_container[iZone][MESH_0][ADJFLOW_SOL]->GetMultiplier());
        solver_container[iZone][MESH_0][ADJFLOW_SOL]->SetAdjoint_Output(geometry_container[iZone][MESH_0],
                                                                        config_container[iZone]);
      }
    }
    AD::ComputeAdjoint();
    for (iZone = 0; iZone < nZone; iZone++){
      if (flow){
        solver_container[iZone][MESH_0][ADJFLOW_SOL]->ResetSensitivity(geometry_container[iZone][MESH_0]);
        solver_container[iZone][MESH_0][ADJFLOW_SOL]->SetSensitivity(geometry_container[iZone][MESH_0],config_container[iZone]); //contains N_u
        solver_container[iZone][MESH_0][ADJFLOW_SOL]->SaveSurfaceSensitivity(geometry_container[iZone][MESH_0]);
        solver_container[iZone][MESH_0][ADJFLOW_SOL]->UpdateLagrangeSensitivity(geometry_container[iZone][MESH_0],1.0);
        solver_container[iZone][MESH_0][ADJFLOW_SOL]->ExtractAdjoint_Solution(geometry_container[iZone][MESH_0],
                                                                       config_container[iZone]);
      }
    }

    AD::ClearAdjoints();

   // cout << "log10Ajdoint[RMS Density]: " << log10(solver_container[ZONE_0][MESH_0][ADJFLOW_SOL]->GetRes_RMS(0))<<std::endl;

    //Calculate Augmented Lagrangian L^a=alpha/2*||dy||^2+beta/2*||dybar||^2+f+ybar^T dy
    for (iZone = 0; iZone < nZone; iZone++){
        solver_container[iZone][MESH_0][ADJFLOW_SOL]->AssembleLagrangian(config_container[iZone]);
    }

   if(ExtIter>=config_container[0]->GetOneShotStart()){

    for (iZone = 0; iZone < nZone; iZone++){
      solver_container[iZone][MESH_0][ADJFLOW_SOL]->StoreSaveSolution();
      if(config_container[0]->GetOneShotConstraint()==true) solver_container[iZone][MESH_0][ADJFLOW_SOL]->StoreConstraint(config_container[0]);
    }

    //------------------COMPUTE alpha*Deltay^TG_u-------------------------
    for (iZone = 0; iZone < nZone; iZone++){
      solver_container[iZone][MESH_0][ADJFLOW_SOL]->LoadOldAdjoint();
    }

    for (iZone = 0; iZone < nZone; iZone++){
      if (flow){
        solver_container[iZone][MESH_0][ADJFLOW_SOL]->SetAdj_ObjFunc(geometry_container[iZone][MESH_0], config_container[iZone],0.0);
        if(config_container[0]->GetOneShotConstraint()==true) solver_container[iZone][MESH_0][ADJFLOW_SOL]->SetAdj_ConstraintFuncAD(geometry_container[iZone][MESH_0], config_container[iZone],solver_container[iZone][MESH_0][ADJFLOW_SOL]->GetConstraintFunc_Value());
      }
      solver_container[iZone][MESH_0][ADJFLOW_SOL]->SetAdjointOutputUpdate(geometry_container[iZone][MESH_0],config_container[iZone]);
    }

    AD::ComputeAdjoint();

    for (iZone = 0; iZone < nZone; iZone++){
      if (flow){
        solver_container[iZone][MESH_0][ADJFLOW_SOL]->SetSensitivity(geometry_container[iZone][MESH_0],config_container[iZone]);
        solver_container[iZone][MESH_0][ADJFLOW_SOL]->UpdateLagrangeSensitivity(geometry_container[iZone][MESH_0],config_container[ZONE_0]->GetOneShotAlpha());
      }
    }

    AD::ClearAdjoints();

   //Compute beta*Deltaybar^TNyu using second order AD (Forward over Reverse)

    for (iZone = 0; iZone < nZone; iZone++){
          solver_container[iZone][MESH_0][ADJFLOW_SOL]->LoadOldSolution();
    }

    AD::Reset();

    for (iZone = 0; iZone < nZone; iZone++){
      solver_container[iZone][MESH_0][ADJFLOW_SOL]->SetForwardDirection(config_container[iZone]);  //set ydot=Deltaybar
    }

    for (iMesh = 0; iMesh <= config_container[ZONE_0]->GetnMGLevels(); iMesh++){
        solver_container[ZONE_0][iMesh][ADJFLOW_SOL]->SetRecording(geometry_container[ZONE_0][MESH_0], config_container[ZONE_0], ALL_VARIABLES);
    }

    AD::StartRecording();

    for (iZone = 0; iZone < nZone; iZone++){
      geometry_container[iZone][MESH_0]->RegisterCoordinates(config_container[iZone]);
      if (flow){
       solver_container[iZone][MESH_0][ADJFLOW_SOL]->RegisterSolution(geometry_container[iZone][MESH_0],
                                                                  config_container[iZone]);
      }
    }

    for (iZone = 0; iZone < nZone; iZone++){
      geometry_container[iZone][MESH_0]->UpdateGeometry(geometry_container[iZone],config_container[iZone]);
    }

    if (flow){
        for (iZone = 0; iZone < nZone; iZone++){
          meanflow_iteration->Iterate(output,integration_container,geometry_container,solver_container,numerics_container,
                                    config_container,surface_movement,volume_grid_movement,FFDBox, iZone);
        }
    }

    for (iZone = 0; iZone < nZone; iZone++){
      if (flow){
        solver_container[iZone][MESH_0][ADJFLOW_SOL]->RegisterObj_Func(config_container[iZone]);
        if(config_container[0]->GetOneShotConstraint()==true) solver_container[iZone][MESH_0][ADJFLOW_SOL]->RegisterConstraint_Func(config_container[iZone], geometry_container[iZone][MESH_0]);
        solver_container[iZone][MESH_0][ADJFLOW_SOL]->RegisterOutput(geometry_container[iZone][MESH_0],
                                                                      config_container[iZone]);
      }
    }

    AD::StopRecording();

    solver_container[ZONE_0][MESH_0][ADJFLOW_SOL]->OutputWritten(geometry_container[ZONE_0][MESH_0]);

    for (iZone = 0; iZone < nZone; iZone++){
      if (flow){
        solver_container[iZone][MESH_0][ADJFLOW_SOL]->SetAdj_ObjFunc(geometry_container[iZone][MESH_0], config_container[iZone],1.0);
        if(config_container[0]->GetOneShotConstraint()==true) solver_container[iZone][MESH_0][ADJFLOW_SOL]->SetAdj_ConstraintFunc(geometry_container[iZone][MESH_0], config_container[iZone],solver_container[iZone][MESH_0][ADJFLOW_SOL]->GetMultiplier());
      }
      solver_container[iZone][MESH_0][ADJFLOW_SOL]->SetAdjoint_Output(geometry_container[iZone][MESH_0], config_container[iZone]);
    }

    AD::ComputeAdjoint();

    for (iZone = 0; iZone < nZone; iZone++){
      if (flow){
        solver_container[iZone][MESH_0][ADJFLOW_SOL]->SetMixedSensitivity(geometry_container[iZone][MESH_0],config_container[iZone]); //read from ubardot (mixed derivative)
        solver_container[iZone][MESH_0][ADJFLOW_SOL]->UpdateLagrangeSensitivity(geometry_container[iZone][MESH_0],config_container[ZONE_0]->GetOneShotBeta());
      }
    }

    AD::ClearAdjoints();

   //Compute beta*Deltaybar^TNyu using Finite Differences

//    for (iZone = 0; iZone < nZone; iZone++){
//      solver_container[iZone][MESH_0][ADJFLOW_SOL]->LoadSaveSolution();
//    }

//    for (iZone = 0; iZone < nZone; iZone++){
//        solver_container[iZone][MESH_0][ADJFLOW_SOL]->UpdateStateVariable(config_container[ZONE_0]);
//    }
// //   cout << "log10Ajdoint[RMS Density]: " << log10(solver_container[ZONE_0][MESH_0][ADJFLOW_SOL]->GetRes_RMS(0))<<std::endl;

//  #if defined CODI_REVERSE_TYPE
//  if (config_container[ZONE_0]->GetOne_Shot()){
//      AD::globalTape.reset();
//   }
//  #endif

//  for (iMesh = 0; iMesh <= config_container[ZONE_0]->GetnMGLevels(); iMesh++){
//    solver_container[ZONE_0][iMesh][ADJFLOW_SOL]->SetRecording(geometry_container[ZONE_0][MESH_0], config_container[ZONE_0], ALL_VARIABLES);
//  }

//    for (iZone = 0; iZone < nZone; iZone++){
//      solver_container[iZone][MESH_0][ADJFLOW_SOL]->LoadOldAdjoint();
//    }

//      AD::StartRecording();

//      for (iZone = 0; iZone < nZone; iZone++){
//        geometry_container[iZone][MESH_0]->RegisterCoordinates(config_container[iZone]);
//        if (flow){
//          solver_container[iZone][MESH_0][ADJFLOW_SOL]->RegisterSolution(geometry_container[iZone][MESH_0],
//                                                                     config_container[iZone]);
//        }
//      }

//      for (iZone = 0; iZone < nZone; iZone++){
//        geometry_container[iZone][MESH_0]->UpdateGeometry(geometry_container[iZone],config_container[iZone]);
//      }

//      if (flow){
//          for (iZone = 0; iZone < nZone; iZone++){
//            meanflow_iteration->Iterate(output,integration_container,geometry_container,solver_container,numerics_container,
//                                      config_container,surface_movement,volume_grid_movement,FFDBox, iZone);
//          }
//      }

//      for (iZone = 0; iZone < nZone; iZone++){
//        if (flow){
//          solver_container[iZone][MESH_0][ADJFLOW_SOL]->RegisterObj_Func(config_container[iZone]);
//          if(config_container[0]->GetOneShotConstraint()==true) solver_container[iZone][MESH_0][ADJFLOW_SOL]->RegisterConstraint_Func(config_container[iZone], geometry_container[iZone][MESH_0]);
//          solver_container[iZone][MESH_0][ADJFLOW_SOL]->RegisterOutput(geometry_container[iZone][MESH_0],
//                                                                      config_container[iZone]);
//        }
//      }

//      AD::StopRecording();

//    //------------------COMPUTE new N_u-------------------------
//    for (iZone = 0; iZone < nZone; iZone++){
//      if (flow){
//        solver_container[iZone][MESH_0][ADJFLOW_SOL]->SetAdj_ObjFunc(geometry_container[iZone][MESH_0], config_container[iZone],1.0);
//        if(config_container[0]->GetOneShotConstraint()==true) solver_container[iZone][MESH_0][ADJFLOW_SOL]->SetAdj_ConstraintFunc(geometry_container[iZone][MESH_0], config_container[iZone],solver_container[iZone][MESH_0][ADJFLOW_SOL]->GetMultiplier());
//        solver_container[iZone][MESH_0][ADJFLOW_SOL]->SetAdjoint_Output(geometry_container[iZone][MESH_0],
//                                                                      config_container[iZone]);
//      }
//    }

//    AD::ComputeAdjoint();

//    for (iZone = 0; iZone < nZone; iZone++){
//      if (flow){
//        solver_container[iZone][MESH_0][ADJFLOW_SOL]->SetSensitivityFD(geometry_container[iZone][MESH_0],config_container[iZone]); //contains N_u_new
//        solver_container[iZone][MESH_0][ADJFLOW_SOL]->UpdateLagrangeSensitivity(geometry_container[iZone][MESH_0],config_container[ZONE_0]->GetOneShotBeta());
//      }
//    }

//    AD::ClearAdjoints();

    for (iZone = 0; iZone < nZone; iZone++){
      solver_container[iZone][MESH_0][ADJFLOW_SOL]->LoadSaveSolution();
    }
   }

   #if defined CODI_REVERSE_TYPE
   if (config_container[ZONE_0]->GetOne_Shot()){
       AD::globalTape.reset();
   }
   #endif
}

void CDiscOneShotIteration::OneShotStepTrial(COutput *output, CIntegration ***integration_container, CGeometry ***geometry_container,
                 CSolver ****solver_container, CNumerics *****numerics_container, CConfig **config_container,
                 CSurfaceMovement **surface_movement, CVolumetricMovement **volume_grid_movement, CFreeFormDefBox*** FFDBox, unsigned short whilecounter) {

    unsigned short iZone, iMesh, ExtIter = config_container[ZONE_0]->GetExtIter();
    unsigned short nZone = config_container[ZONE_0]->GetnZone();

    bool turbulent = false, flow = false;

    switch(config_container[ZONE_0]->GetKind_Solver()){
      case DISC_ADJ_EULER: case DISC_ADJ_NAVIER_STOKES: flow = true; break;
      case DISC_ADJ_RANS: flow = true; turbulent = true; break;
    }

 /*   for (iZone = 0; iZone < nZone; iZone++){
            solver_container[iZone][MESH_0][ADJFLOW_SOL]->StoreOldSolution();
    }*/
    AD::Reset();
    for (iMesh = 0; iMesh <= config_container[ZONE_0]->GetnMGLevels(); iMesh++){
      solver_container[ZONE_0][iMesh][ADJFLOW_SOL]->SetRecording(geometry_container[ZONE_0][MESH_0], config_container[ZONE_0], ALL_VARIABLES);
    }
      AD::StartRecording();

      for (iZone = 0; iZone < nZone; iZone++){
        geometry_container[iZone][MESH_0]->RegisterCoordinates(config_container[iZone]);
        if (flow){
            solver_container[iZone][MESH_0][ADJFLOW_SOL]->RegisterSolution(geometry_container[iZone][MESH_0],
                                                                       config_container[iZone]);
        }
      }
      for (iZone = 0; iZone < nZone; iZone++){
        geometry_container[iZone][MESH_0]->UpdateGeometry(geometry_container[iZone],config_container[iZone]);
       // solver_container[iZone][MESH_0][FLOW_SOL]->Preprocessing(geometry_container[iZone][MESH_0], solver_container[iZone][MESH_0], config_container[iZone], MESH_0, NO_RK_ITER, RUNTIME_FLOW_SYS, false);

      }
      solver_container[ZONE_0][MESH_0][ADJFLOW_SOL]->OutputWritten(geometry_container[ZONE_0][MESH_0]);
      if (flow){
          for (iZone = 0; iZone < nZone; iZone++){
            meanflow_iteration->Iterate(output,integration_container,geometry_container,solver_container,numerics_container,
                                      config_container,surface_movement,volume_grid_movement,FFDBox, iZone);
          }
      }
        if(flow){
       //     cout << "log10Primal[RMS Density]: " << log10(solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetRes_RMS(0))<<std::endl;
            cout     <<"Drag: "<< solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetTotal_CDrag()
                <<", Lift: "<< solver_container[ZONE_0][MESH_0][FLOW_SOL]->GetTotal_CLift() << "." << endl;
        }
      for (iZone = 0; iZone < nZone; iZone++){
        if (flow){
          solver_container[iZone][MESH_0][ADJFLOW_SOL]->RegisterObj_Func(config_container[iZone]);
          if (config_container[0]->GetOneShotConstraint()==true) solver_container[iZone][MESH_0][ADJFLOW_SOL]->RegisterConstraint_Func(config_container[iZone], geometry_container[iZone][MESH_0]);
          solver_container[iZone][MESH_0][ADJFLOW_SOL]->RegisterOutput(geometry_container[iZone][MESH_0],
                                                                      config_container[iZone]);
        }
      }
      AD::StopRecording();

      double* multiplierVec = new double[3];
      multiplierVec[0]=SU2_TYPE::GetValue(config_container[ZONE_0]->GetConstraintStart());
      multiplierVec[1]=0.0;
      multiplierVec[2]=0.0;
      if(ExtIter==0 && whilecounter==1 && config_container[0]->GetOneShotConstraint()==true) solver_container[0][MESH_0][ADJFLOW_SOL]->SetMultiplier(config_container[0],multiplierVec );
       delete [] multiplierVec;
    //------------------COMPUTE N_u-------------------------
    for (iZone = 0; iZone < nZone; iZone++){
      if (flow){
        solver_container[iZone][MESH_0][ADJFLOW_SOL]->SetAdj_ObjFunc(geometry_container[iZone][MESH_0], config_container[iZone],1.0);
        if(config_container[0]->GetOneShotConstraint()==true) solver_container[iZone][MESH_0][ADJFLOW_SOL]->SetAdj_ConstraintFunc(geometry_container[iZone][MESH_0], config_container[iZone],solver_container[iZone][MESH_0][ADJFLOW_SOL]->GetMultiplier());
        solver_container[iZone][MESH_0][ADJFLOW_SOL]->SetAdjoint_Output(geometry_container[iZone][MESH_0],
                                                                        config_container[iZone]);
      }
    }
    AD::ComputeAdjoint();
    for (iZone = 0; iZone < nZone; iZone++){
      if (flow){
        solver_container[iZone][MESH_0][ADJFLOW_SOL]->ResetSensitivity(geometry_container[iZone][MESH_0]);
        solver_container[iZone][MESH_0][ADJFLOW_SOL]->SetSensitivity(geometry_container[iZone][MESH_0],config_container[iZone]); //contains N_u
        solver_container[iZone][MESH_0][ADJFLOW_SOL]->SaveSurfaceSensitivity(geometry_container[iZone][MESH_0]);
        solver_container[iZone][MESH_0][ADJFLOW_SOL]->UpdateLagrangeSensitivity(geometry_container[iZone][MESH_0],1.0);
        solver_container[iZone][MESH_0][ADJFLOW_SOL]->ExtractAdjoint_Solution(geometry_container[iZone][MESH_0],
                                                                       config_container[iZone]);
      }
    }

    AD::ClearAdjoints();

   // cout << "log10Ajdoint[RMS Density]: " << log10(solver_container[ZONE_0][MESH_0][ADJFLOW_SOL]->GetRes_RMS(0))<<std::endl;

    //Calculate Augmented Lagrangian L^a=alpha/2*||dy||^2+beta/2*||dybar||^2+f+ybar^T dy
    for (iZone = 0; iZone < nZone; iZone++){
        solver_container[iZone][MESH_0][ADJFLOW_SOL]->AssembleLagrangian(config_container[iZone]);
    }

   if(ExtIter>=config_container[0]->GetOneShotStart()){

    for (iZone = 0; iZone < nZone; iZone++){
      solver_container[iZone][MESH_0][ADJFLOW_SOL]->StoreSaveSolution();
      if(config_container[0]->GetOneShotConstraint()==true) solver_container[iZone][MESH_0][ADJFLOW_SOL]->StoreConstraint(config_container[0]);
    }

    //------------------COMPUTE alpha*Deltay^TG_u-------------------------
    for (iZone = 0; iZone < nZone; iZone++){
      solver_container[iZone][MESH_0][ADJFLOW_SOL]->LoadOldAdjoint();

    }

    for (iZone = 0; iZone < nZone; iZone++){
      if (flow){
        solver_container[iZone][MESH_0][ADJFLOW_SOL]->SetAdj_ObjFunc(geometry_container[iZone][MESH_0], config_container[iZone],0.0);
        if(config_container[0]->GetOneShotConstraint()==true) solver_container[iZone][MESH_0][ADJFLOW_SOL]->SetAdj_ConstraintFuncAD(geometry_container[iZone][MESH_0], config_container[iZone],solver_container[iZone][MESH_0][ADJFLOW_SOL]->GetConstraintFunc_Value());
      }
      solver_container[iZone][MESH_0][ADJFLOW_SOL]->SetAdjointOutputUpdate(geometry_container[iZone][MESH_0],config_container[iZone]);
    }

    AD::ComputeAdjoint();

    for (iZone = 0; iZone < nZone; iZone++){
      if (flow){
        solver_container[iZone][MESH_0][ADJFLOW_SOL]->SetSensitivity(geometry_container[iZone][MESH_0],config_container[iZone]);
        solver_container[iZone][MESH_0][ADJFLOW_SOL]->UpdateLagrangeSensitivity(geometry_container[iZone][MESH_0],config_container[ZONE_0]->GetOneShotAlpha());
      }
    }

    AD::ClearAdjoints();

    //Compute beta*Deltaybar^TNyu using Finite Differences

    for (iZone = 0; iZone < nZone; iZone++){
      solver_container[iZone][MESH_0][ADJFLOW_SOL]->LoadSaveSolution();
    }

    for (iZone = 0; iZone < nZone; iZone++){
        solver_container[iZone][MESH_0][ADJFLOW_SOL]->UpdateStateVariable(config_container[ZONE_0]);
    }
 //   cout << "log10Ajdoint[RMS Density]: " << log10(solver_container[ZONE_0][MESH_0][ADJFLOW_SOL]->GetRes_RMS(0))<<std::endl;

/*  #if defined CODI_REVERSE_TYPE
  if (config_container[ZONE_0]->GetOne_Shot()){
      AD::globalTape.reset();
   }
  #endif*/

  AD::Reset();

  for (iMesh = 0; iMesh <= config_container[ZONE_0]->GetnMGLevels(); iMesh++){
    solver_container[ZONE_0][iMesh][ADJFLOW_SOL]->SetRecording(geometry_container[ZONE_0][MESH_0], config_container[ZONE_0], ALL_VARIABLES);
  }

    for (iZone = 0; iZone < nZone; iZone++){
      solver_container[iZone][MESH_0][ADJFLOW_SOL]->LoadOldAdjoint();
    }

      AD::StartRecording();

      for (iZone = 0; iZone < nZone; iZone++){
        geometry_container[iZone][MESH_0]->RegisterCoordinates(config_container[iZone]);
        if (flow){
          solver_container[iZone][MESH_0][ADJFLOW_SOL]->RegisterSolution(geometry_container[iZone][MESH_0],
                                                                     config_container[iZone]);
        }
      }

      for (iZone = 0; iZone < nZone; iZone++){
        geometry_container[iZone][MESH_0]->UpdateGeometry(geometry_container[iZone],config_container[iZone]);
      }
      solver_container[ZONE_0][MESH_0][ADJFLOW_SOL]->OutputWritten(geometry_container[ZONE_0][MESH_0]);
      if (flow){
          for (iZone = 0; iZone < nZone; iZone++){
            meanflow_iteration->Iterate(output,integration_container,geometry_container,solver_container,numerics_container,
                                      config_container,surface_movement,volume_grid_movement,FFDBox, iZone);
          }
      }

      for (iZone = 0; iZone < nZone; iZone++){
        if (flow){
          solver_container[iZone][MESH_0][ADJFLOW_SOL]->RegisterObj_Func(config_container[iZone]);
          if(config_container[0]->GetOneShotConstraint()==true) solver_container[iZone][MESH_0][ADJFLOW_SOL]->RegisterConstraint_Func(config_container[iZone], geometry_container[iZone][MESH_0]);
          solver_container[iZone][MESH_0][ADJFLOW_SOL]->RegisterOutput(geometry_container[iZone][MESH_0],
                                                                      config_container[iZone]);
        }
      }

      AD::StopRecording();

    //------------------COMPUTE new N_u-------------------------
    for (iZone = 0; iZone < nZone; iZone++){
      if (flow){
        solver_container[iZone][MESH_0][ADJFLOW_SOL]->SetAdj_ObjFunc(geometry_container[iZone][MESH_0], config_container[iZone],1.0);
        if(config_container[0]->GetOneShotConstraint()==true) solver_container[iZone][MESH_0][ADJFLOW_SOL]->SetAdj_ConstraintFunc(geometry_container[iZone][MESH_0], config_container[iZone],solver_container[iZone][MESH_0][ADJFLOW_SOL]->GetMultiplier());
        solver_container[iZone][MESH_0][ADJFLOW_SOL]->SetAdjoint_Output(geometry_container[iZone][MESH_0],
                                                                      config_container[iZone]);
      }
    }

    AD::ComputeAdjoint();

    for (iZone = 0; iZone < nZone; iZone++){
      if (flow){
        solver_container[iZone][MESH_0][ADJFLOW_SOL]->SetSensitivityFD(geometry_container[iZone][MESH_0],config_container[iZone]); //contains N_u_new
        solver_container[iZone][MESH_0][ADJFLOW_SOL]->UpdateLagrangeSensitivity(geometry_container[iZone][MESH_0],config_container[ZONE_0]->GetOneShotBeta());
      }
    }

    AD::ClearAdjoints();

    for (iZone = 0; iZone < nZone; iZone++){
      solver_container[iZone][MESH_0][ADJFLOW_SOL]->LoadSaveSolution();
    }
   }

   #if defined CODI_REVERSE_TYPE
   if (config_container[ZONE_0]->GetOne_Shot()){
       AD::globalTape.reset();
   }
#endif

}


void CDiscOneShotIteration::projectOneShot(CGeometry ***geometry_container, CConfig **config_container,
                       CSurfaceMovement **surface_movement, CVolumetricMovement **mesh_movement, CSolver ****solver_container) {
    unsigned short iMarker=0, iDim, iDV;
    unsigned long iVertex, iPoint;
    su2double delta_eps, my_Gradient, Gradient, *Normal, dS, *VarCoord, Sensitivity,
    dalpha[3], deps[3], dalpha_deps;

      /* --- No Output of Gradient! ---*/
      bool *UpdatePoint;

      int rank = MASTER_NODE;
    #ifdef HAVE_MPI
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    #endif

    /*--- Boolean controlling points to be updated ---*/

      UpdatePoint = new bool[geometry_container[ZONE_0][MESH_0]->GetnPoint()];

      if (rank == MASTER_NODE)
          cout << endl <<"---------- Start gradient evaluation using surface sensitivity ----------" << endl;

      for (iDV = 0; iDV < config_container[ZONE_0]->GetnDV(); iDV++) {

        config_container[ZONE_0]->SetDV_Value(iDV,0, 1E-4);

        /*--- Hicks Henne design variable ---*/

          surface_movement[ZONE_0]->SetHicksHenne(geometry_container[ZONE_0][MESH_0], config_container[ZONE_0], iDV, true);


        /*--- Continuous adjoint gradient computation ---*/
        if (rank == MASTER_NODE)

        /*--- Load the delta change in the design variable (finite difference step). ---*/

        delta_eps = config_container[ZONE_0]->GetDV_Value(iDV);
        my_Gradient = 0.0; Gradient = 0.0;

        /*--- Reset update points ---*/

        for (iPoint = 0; iPoint < geometry_container[ZONE_0][MESH_0]->GetnPoint(); iPoint++)
          UpdatePoint[iPoint] = true;

 //       for (iMarker = 0; iMarker < config_container[ZONE_0]->GetnMarker_All(); iMarker++) {
          if (config_container[ZONE_0]->GetMarker_All_DV(iMarker) == YES) {
            for (iVertex = 0; iVertex < geometry_container[ZONE_0][MESH_0]->nVertex[iMarker]; iVertex++) {

              iPoint = geometry_container[ZONE_0][MESH_0]->vertex[iMarker][iVertex]->GetNode();
              if ((iPoint < geometry_container[ZONE_0][MESH_0]->GetnPointDomain()) && UpdatePoint[iPoint]) {

                Normal = geometry_container[ZONE_0][MESH_0]->vertex[iMarker][iVertex]->GetNormal();
                VarCoord = geometry_container[ZONE_0][MESH_0]->vertex[iMarker][iVertex]->GetVarCoord();
                Sensitivity = geometry_container[ZONE_0][MESH_0]->vertex[iMarker][iVertex]->GetAuxVar();

                dS = 0.0;
                for (iDim = 0; iDim < geometry_container[ZONE_0][MESH_0]->GetnDim(); iDim++) {
                  dS += Normal[iDim]*Normal[iDim];
                  deps[iDim] = VarCoord[iDim] / delta_eps;
                }
                dS = sqrt(dS);

                dalpha_deps = 0.0;
                for (iDim = 0; iDim < geometry_container[ZONE_0][MESH_0]->GetnDim(); iDim++) {
                  dalpha[iDim] = Normal[iDim] / dS;
                  dalpha_deps -= dalpha[iDim]*deps[iDim];
                }

                /*--- Store the geometric sensitivity for this DV (rows) & this node (column) ---*/

                my_Gradient += Sensitivity*dalpha_deps;
                UpdatePoint[iPoint] = false;
              }
            }
          }
  //      }

  #ifdef HAVE_MPI
        SU2_MPI::Allreduce(&my_Gradient, &Gradient, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  #else
        Gradient = my_Gradient;
  #endif

        solver_container[ZONE_0][MESH_0][ADJFLOW_SOL]->SetProjectedSensitivity(iDV,Gradient);

      }
      delete [] UpdatePoint;
}

void CDiscOneShotIteration::deformOneShot(COutput *output, CGeometry ***geometry_container, CConfig **config_container,
                       CSurfaceMovement **surface_movement, CVolumetricMovement **grid_movement, CSolver ****solver_container) {
    int rank = MASTER_NODE;
  #ifdef HAVE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  #endif
    std::cout<<"Update of Design Variables for Deformation: "<<std::endl;
    unsigned short iDV;

    //Use updatecount to count number of non-zero design variables to check if (updatecount > user defined variable) (default: updatecount>0)
    unsigned short updatecount=0;

    /*--- Set the Design Variable (Stored in DesignVarUpdate)*/
    for (iDV = 0; iDV < config_container[ZONE_0]->GetnDV(); iDV++) {
          config_container[ZONE_0]->SetDV_Value(iDV,0, solver_container[ZONE_0][MESH_0][ADJFLOW_SOL]->getDVValue(iDV));
          std::cout<<solver_container[ZONE_0][MESH_0][ADJFLOW_SOL]->getDVValue(iDV)<<" ";
          if(solver_container[ZONE_0][MESH_0][ADJFLOW_SOL]->getDVValue(iDV)!=0)   updatecount=updatecount+1;
    }
    std::cout<<std::endl;

    if(updatecount>0){
    /*--- Surface grid deformation using design variables ---*/

    /*--- Surface grid deformation ---*/

    surface_movement[ZONE_0]->SetSurface_Deformation(geometry_container[ZONE_0][MESH_0], config_container[ZONE_0]);

    /*--- Volumetric grid deformation ---*/

    if (config_container[ZONE_0]->GetDesign_Variable(0) != FFD_SETTING) {

      grid_movement[ZONE_0]->SetVolume_Deformation(geometry_container[ZONE_0][MESH_0], config_container[ZONE_0], true);

    }

  }
    /*--- Computational grid preprocessing ---*/
    /*  ---> not writing deformed grid in this case ----*/
}

void CDiscOneShotIteration::projectGradient(CGeometry ***geometry_container, CConfig **config_container,
                       CSurfaceMovement **surface_movement, CVolumetricMovement **mesh_movement, CSolver ****solver_container) {
    unsigned short iMarker=0, iDim, iDV;
    unsigned long iVertex, iPoint;
    su2double delta_eps, my_Gradient, Gradient, *Normal, dS, *VarCoord, Sensitivity,
    dalpha[3], deps[3], dalpha_deps;
      bool *UpdatePoint;

      int rank = MASTER_NODE;
    #ifdef HAVE_MPI
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    #endif

      UpdatePoint = new bool[geometry_container[ZONE_0][MESH_0]->GetnPoint()];

      if (rank == MASTER_NODE)
          cout << endl <<"---------- Start gradient evaluation using surface sensitivity ----------" << endl;

      for (iDV = 0; iDV < config_container[ZONE_0]->GetnDV(); iDV++) {

        config_container[ZONE_0]->SetDV_Value(iDV,0, 1E-4);

        /*--- Hicks Henne design variable ---*/

          surface_movement[ZONE_0]->SetHicksHenne(geometry_container[ZONE_0][MESH_0], config_container[ZONE_0], iDV, true);


        /*--- Continuous adjoint gradient computation ---*/

        /*--- Load the delta change in the design variable (finite difference step). ---*/

        delta_eps = config_container[ZONE_0]->GetDV_Value(iDV);
        my_Gradient = 0.0; Gradient = 0.0;

        /*--- Reset update points ---*/

        for (iPoint = 0; iPoint < geometry_container[ZONE_0][MESH_0]->GetnPoint(); iPoint++)
          UpdatePoint[iPoint] = true;

   //     for (iMarker = 0; iMarker < config_container[ZONE_0]->GetnMarker_All(); iMarker++) {
          if (config_container[ZONE_0]->GetMarker_All_DV(iMarker) == YES) {
            for (iVertex = 0; iVertex < geometry_container[ZONE_0][MESH_0]->nVertex[iMarker]; iVertex++) {

              iPoint = geometry_container[ZONE_0][MESH_0]->vertex[iMarker][iVertex]->GetNode();
              if ((iPoint < geometry_container[ZONE_0][MESH_0]->GetnPointDomain()) && UpdatePoint[iPoint]) {

                Normal = geometry_container[ZONE_0][MESH_0]->vertex[iMarker][iVertex]->GetNormal();
                VarCoord = geometry_container[ZONE_0][MESH_0]->vertex[iMarker][iVertex]->GetVarCoord();
                Sensitivity = geometry_container[ZONE_0][MESH_0]->vertex[iMarker][iVertex]->GetAuxVar();
                //Coord= geometry_container[ZONE_0][MESH_0]->node[iPoint]->GetCoord()

                dS = 0.0;
                for (iDim = 0; iDim < geometry_container[ZONE_0][MESH_0]->GetnDim(); iDim++) {
                  dS += Normal[iDim]*Normal[iDim];
                  deps[iDim] = VarCoord[iDim] / delta_eps;
                }
                dS = sqrt(dS);

                dalpha_deps = 0.0;
                for (iDim = 0; iDim < geometry_container[ZONE_0][MESH_0]->GetnDim(); iDim++) {
                  dalpha[iDim] = Normal[iDim] / dS;
                  dalpha_deps -= dalpha[iDim]*deps[iDim];
                }

                /*--- Store the geometric sensitivity for this DV (rows) & this node (column) ---*/

                my_Gradient += Sensitivity*dalpha_deps;
                UpdatePoint[iPoint] = false;
              }
            }
          }
    //    }


  #ifdef HAVE_MPI
        SU2_MPI::Allreduce(&my_Gradient, &Gradient, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  #else
        Gradient = my_Gradient;
  #endif
        solver_container[ZONE_0][MESH_0][ADJFLOW_SOL]->SetProjectedGradient(iDV,Gradient);
    }
      delete [] UpdatePoint;

}

TopologyOptimization::TopologyOptimization(CConfig *config) : CIteration(config), CurrentRecording(NONE){

  mean_iteration = new CFEM_StructuralAnalysis(config);

//  turbulent = config->GetKind_Solver() == DISC_ADJ_RANS;

}

TopologyOptimization::~TopologyOptimization(void) { }
void TopologyOptimization::Preprocess(COutput *output,
                                           CIntegration ***integration_container,
                                           CGeometry ***geometry_container,
                                           CSolver ****solver_container,
                                           CNumerics *****numerics_container,
                                           CConfig **config_container,
                                           CSurfaceMovement **surface_movement,
                                           CVolumetricMovement **grid_movement,
                                           CFreeFormDefBox*** FFDBox,
                                           unsigned short val_iZone) { }

void TopologyOptimization::Iterate(COutput *output,
                                        CIntegration ***integration_container,
                                        CGeometry ***geometry_container,
                                        CSolver ****solver_container,
                                        CNumerics *****numerics_container,
                                        CConfig **config_container,
                                        CSurfaceMovement **surface_movement,
                                        CVolumetricMovement **volume_grid_movement,
                                        CFreeFormDefBox*** FFDBox,
                                        unsigned short val_iZone) {
    unsigned long IntIter = 0; config_container[val_iZone]->SetIntIter(IntIter);
    unsigned short iZone = val_iZone;

    unsigned short SolContainer_Position = config_container[val_iZone]->GetContainerPosition(RUNTIME_FEA_SYS);

    solver_container[iZone][MESH_0][SolContainer_Position]->Initialize_Density(geometry_container[iZone][MESH_0], config_container[iZone]);

    AD::Reset();

    solver_container[iZone][MESH_0][ADJFEA_SOL]->SetRecording(geometry_container[iZone][MESH_0], config_container[iZone], ALL_VARIABLES);

    AD::StartRecording();

    geometry_container[iZone][MESH_0]->RegisterDensity(config_container[iZone]);
    solver_container[iZone][MESH_0][ADJFEA_SOL]->RegisterSolution(geometry_container[iZone][MESH_0],
                                                                       config_container[iZone]);

    mean_iteration->Iterate(output,integration_container,geometry_container,solver_container,numerics_container,
                            config_container,surface_movement,volume_grid_movement,FFDBox, iZone);

    solver_container[iZone][MESH_0][ADJFEA_SOL]->RegisterObj_Func(config_container[iZone]);
    solver_container[iZone][MESH_0][ADJFEA_SOL]->RegisterOutput(geometry_container[iZone][MESH_0],
                                                                config_container[iZone]);

    AD::StopRecording();

    solver_container[iZone][MESH_0][ADJFEA_SOL]->SetAdj_ObjFunc(geometry_container[iZone][MESH_0], config_container[iZone],1.0);
    solver_container[iZone][MESH_0][ADJFEA_SOL]->SetAdjoint_Output(geometry_container[iZone][MESH_0],
                                                                    config_container[iZone]);

    AD::ComputeAdjoint();

    solver_container[iZone][MESH_0][ADJFEA_SOL]->SetSensDensity(geometry_container[iZone][MESH_0],config_container[iZone]);
    solver_container[iZone][MESH_0][ADJFEA_SOL]->ExtractAdjoint_Solution(geometry_container[iZone][MESH_0],
                                                                       config_container[iZone]);
    solver_container[iZone][MESH_0][ADJFEA_SOL]->OutputWritten(geometry_container[iZone][MESH_0]);

    /*solver_container[iZone][MESH_0][ADJFLOW_SOL]->SaveSurfaceSensitivity(geometry_container[iZone][MESH_0]);
    solver_container[iZone][MESH_0][ADJFLOW_SOL]->UpdateLagrangeSensitivity(geometry_container[iZone][MESH_0],1.0);

    */

    AD::ClearAdjoints();

  #if defined CODI_REVERSE_TYPE
  if (config_container[ZONE_0]->GetOne_Shot()){
      AD::globalTape.reset();
   }
  #endif
  //*/

}

void TopologyOptimization::Update(COutput *output,
                                       CIntegration ***integration_container,
                                       CGeometry ***geometry_container,
                                       CSolver ****solver_container,
                                       CNumerics *****numerics_container,
                                       CConfig **config_container,
                                       CSurfaceMovement **surface_movement,
                                       CVolumetricMovement **grid_movement,
                                       CFreeFormDefBox*** FFDBox,
                                       unsigned short val_iZone)      { }
void TopologyOptimization::Monitor()     { }
void TopologyOptimization::Output()      { }
void TopologyOptimization::Postprocess() { }



