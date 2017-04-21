/*!
 * \file CKineticVariable.hpp
 * \brief Declaration of the kinetic variable class.
 * \author A. Mogavero
 * \version 5.0.0 "Raven"
 *
 * SU2 Lead Developers: Dr. Francisco Palacios (Francisco.D.Palacios@boeing.com).
 *                      Dr. Thomas D. Economon (economon@stanford.edu).
 *
 * SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
 *                 Prof. Piero Colonna's group at Delft University of Technology.
 *                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *                 Prof. Rafael Palacios' group at Imperial College London.
 *                 Prof. Edwin van der Weide's group at the University of Twente.
 *                 Prof. Vincent Terrapon's group at the University of Liege.
 *                 Dr. Marco Fossati's group at the University of Strathclyde.
 *
 * Copyright (C) 2012-2017 SU2, the open-source CFD code.
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

#include "variable_structure.hpp"

/*!
 * \class CKineticVariable
 * \brief Main class for defining the variables used in kinetic solver.
 * \ingroup Navier_Stokes_Equations
 * \author A. Mogavero
 * \version 5.0.0 "Raven"
 */
class CKineticVariable : public CNSVariable {
  su2double knudsenLocal; /*!< Local Knudsen number */

  /*!
   * \brief Calculates the magnitude of a vector.
   *
   * @param[in] v - Vector of size nDim
   * @return norm of order 2 of v
   */
  su2double CalcMagnitude(su2double* v)const;
public:
  /*!
   * \brief Constructor of the class.
   */
  CKineticVariable(void);

  /*!
   * \overload
   * \param[in] val_density - Value of the flow density (initialization value).
   * \param[in] val_velocity - Value of the flow velocity (initialization value).
   * \param[in] val_energy - Value of the flow energy (initialization value).
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nvar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CKineticVariable(su2double val_density, su2double *val_velocity,
              su2double val_energy, unsigned short val_nDim, unsigned short val_nvar, CConfig *config);

  /*!
   * \overload
   * \param[in] val_solution - Pointer to the flow value (initialization value).
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nvar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CKineticVariable(su2double *val_solution, unsigned short val_nDim, unsigned short val_nvar, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CKineticVariable(void);

  /*!
   * \brief Calculate local Knudsen number.
   * \details It uses the primitive values and gradients.
   * That therefore need to be calculated before calling this method.
   */
  void CalculateKnudsen();
};
