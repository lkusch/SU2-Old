/*!
 * \file COneShotFluidIteration.hpp
 * \brief Headers of the iteration classes used by SU2_CFD.
 *        Each CIteration class represents an available physics package.
 * \author F. Palacios, T. Economon, L.Kusch
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

#pragma once

#include "CDiscAdjFluidIteration.hpp"

class CFluidIteration;

/*!
 * \class COneShotFluidIteration
 * \brief Class for driving an iteration of the one-shot for fluid system.
 * \author L. Kusch
 */
class COneShotFluidIteration : public CDiscAdjFluidIteration {
 private:
  const bool turbulent;                      /*!< \brief Stores the turbulent flag. */

 public:
  /*!
   * \brief Constructor of the class.
   * \param[in] config - Definition of the particular problem.
   */
  explicit COneShotFluidIteration(const CConfig *config) : CDiscAdjFluidIteration(config),
    turbulent(config->GetKind_Solver() == ONE_SHOT_RANS) {}

  /*!
   * \brief Registers all output variables of the fluid iteration.
   * (Both, flow and geometry variables have to be registered in each iteration)
   * \param[in] solver - Container vector with all the solutions.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iZone - Index of the zone.
   * \param[in] iInst - Index of the instance.
   * \param[in] kind_recording - Kind of recording, either FLOW_VARIABLES, SOLUTION_AND_MESH or GEOMETRY_VARIABLES
   */
  void RegisterInput(CSolver***** solver, CGeometry**** geometry, CConfig** config, unsigned short iZone,
                     unsigned short iInst, unsigned short kind_recording) override;

  /*!
   * \brief Compute necessary variables that depend on the conservative variables AND the mesh node positions
   * (e.g. turbulent variables, normals, volumes). AND - difference to usual method
   * \param[in] solver - Container vector with all the solutions.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iZone - Index of the zone.
   * \param[in] iInst - Index of the zone.
   * \param[in] kind_recording - The kind of recording (geometry or flow).
   */
  void SetDependencies(CSolver***** solver, CGeometry**** geometry, CNumerics****** numerics, CConfig** config,
                       unsigned short iZone, unsigned short iInst, unsigned short kind_recording) override;
  /*!
   * \brief Registers all input variables of the fluid iteration. - The objective function is not set and
   * instead of an adjoint variable the difference in direct variables is set (needed for doubly augmented Lagrangian)
   * \param[in] solver - Container vector with all the solutions.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iZone - Index of the zone.
   * \param[in] iInst - Index of the instance.
   */
  void InitializeAdjoint_Update(CSolver***** solver, CGeometry**** geometry, CConfig** config, unsigned short iZone,
                         unsigned short iInst) override;

  /*!
   * \brief Registers all input variables of the fluid iteration. - Setting the adjoint output to zero
   * \param[in] solver - Container vector with all the solutions.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iZone - Index of the zone.
   * \param[in] iInst - Index of the instance.
   */
  void InitializeAdjoint_Zero(CSolver***** solver, CGeometry**** geometry, CConfig** config, unsigned short iZone,
                         unsigned short iInst) override;
  /*!
   * \brief Perform a single iteration of the adjoint fluid system without calculating residuals.
   * \param[in] output - Pointer to the COutput class.
   * \param[in] integration - Container vector with all the integration methods.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver - Container vector with all the solutions.
   * \param[in] numerics - Description of the numerical method (the way in which the equations are solved).
   * \param[in] config - Definition of the particular problem.
   * \param[in] surface_movement - Surface movement classes of the problem.
   * \param[in] grid_movement - Volume grid movement classes of the problem.
   * \param[in] FFDBox - FFD FFDBoxes of the problem.
   * \param[in] val_iZone - Index of the zone.
   * \param[in] val_iInst - Index of the instance
   */
  void Iterate_No_Residual(COutput* output, CIntegration**** integration, CGeometry**** geometry, CSolver***** solver,
               CNumerics****** numerics, CConfig** config, CSurfaceMovement** surface_movement,
               CVolumetricMovement*** grid_movement, CFreeFormDefBox*** FFDBox, unsigned short val_iZone,
               unsigned short val_iInst) override;
};
