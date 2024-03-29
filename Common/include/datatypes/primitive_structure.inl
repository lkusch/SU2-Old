/*!
 * \file primitive_structure.inl
 * \brief Inline subroutines for <i>datatype_structure.hpp<i>.
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
#pragma once
namespace SU2_TYPE{
  inline void SetValue(double& data, const double &val){data = val;}

  inline double GetValue(const double& data){ return data;}

  inline void SetSecondary(double& data, const double &val){}

  inline double GetDerivative(const double& data){return 0.0;}

  inline double GetSecondary(const double& data){return 0.0;}

  inline void SetDerivative(double &data, const double &val){}

  inline double GetForwardDerivative(const double& data){return 0.0;}
  inline void SetForwardDerivative(double &data, const double &val){}

  inline double GetMixedDerivative(const double& data){return 0.0;}
  inline void SetMixedDerivative(double &data, const double &val){}

}
