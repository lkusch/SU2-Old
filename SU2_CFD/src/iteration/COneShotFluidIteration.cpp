/*!
 * \file COneShotFluidIteration.cpp
 * \brief Main subroutines used by SU2_CFD
 * \author L.Kusch
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

#include "../../include/iteration/COneShotFluidIteration.hpp"
#include "../../include/output/COutput.hpp"

void COneShotFluidIteration::RegisterInput(CSolver***** solver, CGeometry**** geometry, CConfig** config,
                                           unsigned short iZone, unsigned short iInst, unsigned short kind_recording) {

  /*--- For the one-shot strategy conservative variables as well as mesh coordinates are recorded. Furthermore, we need to record the mesh coordinates in every flow iteration,
   *  thus we make use of the COMBINED recording in each step ---*/

  bool frozen_visc = config[iZone]->GetFrozen_Visc_Disc();
  bool heat = config[iZone]->GetWeakly_Coupled_Heat();

  if (kind_recording == SOLUTION_VARIABLES || kind_recording == SOLUTION_AND_MESH) {
    /*--- Register flow and turbulent variables as input ---*/

    if (config[iZone]->GetFluidProblem()) {
      solver[iZone][iInst][MESH_0][ADJFLOW_SOL]->RegisterSolution(geometry[iZone][iInst][MESH_0], config[iZone]);

      solver[iZone][iInst][MESH_0][ADJFLOW_SOL]->RegisterVariables(geometry[iZone][iInst][MESH_0], config[iZone]);
    }

    if (turbulent && !frozen_visc) {
      solver[iZone][iInst][MESH_0][ADJTURB_SOL]->RegisterSolution(geometry[iZone][iInst][MESH_0], config[iZone]);
    }
    if (heat) {
      solver[iZone][iInst][MESH_0][ADJHEAT_SOL]->RegisterSolution(geometry[iZone][iInst][MESH_0], config[iZone]);
    }
    if (config[iZone]->AddRadiation()) {
      solver[iZone][iInst][MESH_0][ADJRAD_SOL]->RegisterSolution(geometry[iZone][iInst][MESH_0], config[iZone]);

      solver[iZone][iInst][MESH_0][ADJRAD_SOL]->RegisterVariables(geometry[iZone][iInst][MESH_0], config[iZone]);
    }
  }

  if (kind_recording == MESH_COORDS || kind_recording == SOLUTION_AND_MESH) {
    /*--- Register node coordinates as input ---*/

    geometry[iZone][iInst][MESH_0]->RegisterCoordinates(config[iZone]);
  }

  /*--- Register the variables of the mesh deformation ---*/
  //TODO-ONE-SHOT Do I need this?
  if (kind_recording == MESH_DEFORM) {
    /*--- Undeformed mesh coordinates ---*/
    solver[iZone][iInst][MESH_0][ADJMESH_SOL]->RegisterSolution(geometry[iZone][iInst][MESH_0], config[iZone]);

    /*--- Boundary displacements ---*/
    solver[iZone][iInst][MESH_0][ADJMESH_SOL]->RegisterVariables(geometry[iZone][iInst][MESH_0], config[iZone]);
  }
}

void COneShotFluidIteration::SetDependencies(CSolver***** solver, CGeometry**** geometry, CNumerics****** numerics,
                                             CConfig** config, unsigned short iZone, unsigned short iInst,
                                             unsigned short kind_recording) {

  /*--- For the one-shot strategy conservative variables as well as mesh coordinates are recorded. Furthermore, we need to record the mesh coordinates in every flow iteration,
   *  thus we make use of the COMBINED recording in each step ---*/

  bool frozen_visc = config[iZone]->GetFrozen_Visc_Disc();
  bool heat = config[iZone]->GetWeakly_Coupled_Heat();

  if ((kind_recording == MESH_COORDS) || (kind_recording == NONE) || (kind_recording == SOLUTION_AND_MESH)) {
    /*--- Update geometry to get the influence on other geometry variables (normals, volume etc) ---*/

    geometry[iZone][iInst][MESH_0]->UpdateGeometry(geometry[iZone][iInst], config[iZone]);
    //TODO-ONE-SHOT This was not here before...
    CGeometry::ComputeWallDistance(config, geometry); 
  }

  /*--- Compute coupling between flow and turbulent equations ---*/
  solver[iZone][iInst][MESH_0][FLOW_SOL]->Preprocessing(geometry[iZone][iInst][MESH_0], solver[iZone][iInst][MESH_0],
                                                        config[iZone], MESH_0, NO_RK_ITER, RUNTIME_FLOW_SYS, true);
  solver[iZone][iInst][MESH_0][FLOW_SOL]->InitiateComms(geometry[iZone][iInst][MESH_0], config[iZone], SOLUTION);
  solver[iZone][iInst][MESH_0][FLOW_SOL]->CompleteComms(geometry[iZone][iInst][MESH_0], config[iZone], SOLUTION);

  if (turbulent && !frozen_visc) {
    solver[iZone][iInst][MESH_0][TURB_SOL]->Postprocessing(geometry[iZone][iInst][MESH_0], solver[iZone][iInst][MESH_0],
                                                           config[iZone], MESH_0);
    solver[iZone][iInst][MESH_0][TURB_SOL]->InitiateComms(geometry[iZone][iInst][MESH_0], config[iZone], SOLUTION);
    solver[iZone][iInst][MESH_0][TURB_SOL]->CompleteComms(geometry[iZone][iInst][MESH_0], config[iZone], SOLUTION);
  }

  if (heat) {
    solver[iZone][iInst][MESH_0][HEAT_SOL]->Set_Heatflux_Areas(geometry[iZone][iInst][MESH_0], config[iZone]);
    solver[iZone][iInst][MESH_0][HEAT_SOL]->Preprocessing(geometry[iZone][iInst][MESH_0], solver[iZone][iInst][MESH_0],
                                                          config[iZone], MESH_0, NO_RK_ITER, RUNTIME_HEAT_SYS, true);
    solver[iZone][iInst][MESH_0][HEAT_SOL]->Postprocessing(geometry[iZone][iInst][MESH_0], solver[iZone][iInst][MESH_0],
                                                           config[iZone], MESH_0);
    solver[iZone][iInst][MESH_0][HEAT_SOL]->InitiateComms(geometry[iZone][iInst][MESH_0], config[iZone], SOLUTION);
    solver[iZone][iInst][MESH_0][HEAT_SOL]->CompleteComms(geometry[iZone][iInst][MESH_0], config[iZone], SOLUTION);
  }
  if (config[iZone]->AddRadiation()) {
    solver[iZone][iInst][MESH_0][RAD_SOL]->Postprocessing(geometry[iZone][iInst][MESH_0], solver[iZone][iInst][MESH_0],
                                                          config[iZone], MESH_0);
    solver[iZone][iInst][MESH_0][RAD_SOL]->InitiateComms(geometry[iZone][iInst][MESH_0], config[iZone], SOLUTION);
    solver[iZone][iInst][MESH_0][RAD_SOL]->CompleteComms(geometry[iZone][iInst][MESH_0], config[iZone], SOLUTION);
  }
}

void COneShotFluidIteration::InitializeAdjoint_Update(CSolver***** solver, CGeometry**** geometry, CConfig** config,
                                               unsigned short iZone, unsigned short iInst) {
  bool frozen_visc = config[iZone]->GetFrozen_Visc_Disc();
  bool heat = config[iZone]->GetWeakly_Coupled_Heat();

  /*--- Initialize the adjoints the conservative variables ---*/

  if (config[iZone]->GetFluidProblem()) {
    solver[iZone][iInst][MESH_0][ADJFLOW_SOL]->SetAdjoint_OutputUpdate(geometry[iZone][iInst][MESH_0], config[iZone]);
  }

  if (turbulent && !frozen_visc) {
    solver[iZone][iInst][MESH_0][ADJTURB_SOL]->SetAdjoint_OutputUpdate(geometry[iZone][iInst][MESH_0], config[iZone]);
  }

  if (heat) {
    solver[iZone][iInst][MESH_0][ADJHEAT_SOL]->SetAdjoint_OutputUpdate(geometry[iZone][iInst][MESH_0], config[iZone]);
  }

  if (config[iZone]->AddRadiation()) {
    solver[iZone][iInst][MESH_0][ADJRAD_SOL]->SetAdjoint_OutputUpdate(geometry[iZone][iInst][MESH_0], config[iZone]);
  }

  //TODO-ONE-SHOT
/*  if (config[iZone]->GetFluidProblem()) {
    solver[iZone][iInst][MESH_0][FLOW_SOL]->SetVertexTractionsAdjoint(geometry[iZone][iInst][MESH_0], config[iZone]);
  }
*/
}

void COneShotFluidIteration::InitializeAdjoint_Zero(CSolver***** solver, CGeometry**** geometry, CConfig** config,
                                               unsigned short iZone, unsigned short iInst) {
  bool frozen_visc = config[iZone]->GetFrozen_Visc_Disc();
  bool heat = config[iZone]->GetWeakly_Coupled_Heat();

  /*--- Initialize the adjoints the conservative variables ---*/

  if (config[iZone]->GetFluidProblem()) {
    solver[iZone][iInst][MESH_0][ADJFLOW_SOL]->SetAdjoint_OutputZero(geometry[iZone][iInst][MESH_0], config[iZone]);
  }

  if (turbulent && !frozen_visc) {
    solver[iZone][iInst][MESH_0][ADJTURB_SOL]->SetAdjoint_OutputZero(geometry[iZone][iInst][MESH_0], config[iZone]);
  }

  if (heat) {
    solver[iZone][iInst][MESH_0][ADJHEAT_SOL]->SetAdjoint_OutputZero(geometry[iZone][iInst][MESH_0], config[iZone]);
  }

  if (config[iZone]->AddRadiation()) {
    solver[iZone][iInst][MESH_0][ADJRAD_SOL]->SetAdjoint_OutputZero(geometry[iZone][iInst][MESH_0], config[iZone]);
  }

  //TODO-ONE-SHOT
/*  if (config[iZone]->GetFluidProblem()) {
    solver[iZone][iInst][MESH_0][FLOW_SOL]->SetVertexTractionsAdjointZero(geometry[iZone][iInst][MESH_0], config[iZone]);
  }
*/
}

void COneShotFluidIteration::Iterate_No_Residual(COutput* output, CIntegration**** integration, CGeometry**** geometry,
                                     CSolver***** solver, CNumerics****** numerics, CConfig** config,
                                     CSurfaceMovement** surface_movement, CVolumetricMovement*** volume_grid_movement,
                                     CFreeFormDefBox*** FFDBox, unsigned short iZone, unsigned short iInst) {
  bool frozen_visc = config[iZone]->GetFrozen_Visc_Disc();
  bool heat = config[iZone]->GetWeakly_Coupled_Heat();

  /*--- Extract the adjoints of the conservative input variables and store them for the next iteration ---*/

  if (config[iZone]->GetFluidProblem()) {
    solver[iZone][iInst][MESH_0][ADJFLOW_SOL]->ExtractAdjoint_Solution_Clean(geometry[iZone][iInst][MESH_0], config[iZone]);

    solver[iZone][iInst][MESH_0][ADJFLOW_SOL]->ExtractAdjoint_Variables(geometry[iZone][iInst][MESH_0], config[iZone]);
  }
  if (turbulent && !frozen_visc) {
    solver[iZone][iInst][MESH_0][ADJTURB_SOL]->ExtractAdjoint_Solution_Clean(geometry[iZone][iInst][MESH_0], config[iZone]);
  }
  if (heat) {
    solver[iZone][iInst][MESH_0][ADJHEAT_SOL]->ExtractAdjoint_Solution_Clean(geometry[iZone][iInst][MESH_0], config[iZone]);
  }
  if (config[iZone]->AddRadiation()) {
    solver[iZone][iInst][MESH_0][ADJRAD_SOL]->ExtractAdjoint_Solution_Clean(geometry[iZone][iInst][MESH_0], config[iZone]);

    solver[iZone][iInst][MESH_0][ADJRAD_SOL]->ExtractAdjoint_Variables(geometry[iZone][iInst][MESH_0], config[iZone]);
  }
}

