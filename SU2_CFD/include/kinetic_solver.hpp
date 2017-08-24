/*!
 * \file kinetic_solver.hpp
 * \brief Declaration of the kinetic solver class.
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

#include "solver_structure.hpp"

/*!
 * \class CKineticSolver
 * \brief Main class for defining the solver using kinetics schemes.
 * \ingroup Navier_Stokes_Equations
 * \author A. Mogavero
 * \version 5.0.0 "Raven"
 */
class CKineticSolver : public CNSSolver {
  /*!
   * \brief Impose the kinetic wall boundary condition.
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] conv_numerics - Description of the numerical method.
   * \param[in] visc_numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] val_marker - Surface marker where the boundary condition is applied.
   * \param[in] accom - accomodation coefficient 0=adiabatic 1=isothermal
   */
  void BC_Kinetic_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
      CNumerics *visc_numerics, CConfig *config, unsigned short val_marker, su2double accom);

public:
  /*!
   * \brief Constructor of the class.
   */
  CKineticSolver(void);

  /*!
   * \overload
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CKineticSolver(CGeometry *geometry, CConfig *config, unsigned short iMesh);

  /*!
   * \brief Destructor of the class.
   */
  ~CKineticSolver(void);

  unsigned long SetPrimitive_Variables(CSolver **solver_container, CConfig *config, bool Output);

  /*!
   * \brief Does nothing because the kinetic fluxes already contain the viscous fluxes
   * \param[in] geometry - Geometrical definition of the problem.
   * \param[in] solver_container - Container vector with all the solutions.
   * \param[in] numerics - Description of the numerical method.
   * \param[in] config - Definition of the particular problem.
   * \param[in] iMesh - Index of the mesh in multigrid computations.
   * \param[in] iRKStep - Current step of the Runge-Kutta iteration.
   */
  void Viscous_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                   CConfig *config, unsigned short iMesh, unsigned short iRKStep);

  void BC_Isothermal_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
      CNumerics *visc_numerics, CConfig *config, unsigned short val_marker);

  void Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config,
      unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output);
};
