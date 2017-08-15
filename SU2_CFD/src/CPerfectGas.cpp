/*!
 * CPerfectGas.cpp
 * \brief Source of the generic ideal gas model.
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

#include "../include/CPerfectGas.hpp"

CPerfectGas::CPerfectGas() :
	CGeneralIdealGas(),
	n_dof(0),
	Gamma(0){
}


CPerfectGas::CPerfectGas(su2double R, su2double g):
		CGeneralIdealGas(R),
	Gamma(g),
	n_dof(0){
	n_dof = 2/(Gamma - 1);
	Cv = 0.5*n_dof*Gas_Constant;
}


CPerfectGas::~CPerfectGas(void) {

}

su2double CPerfectGas::Energy(su2double T)const{
	return Cv*T;
}

su2double CPerfectGas::EnergyInv(su2double e)const{
	return e/Cv;
}

su2double CPerfectGas::SpecificHeatVol(su2double T)const{
	return Cv;
}

su2double CPerfectGas::EntropyTemp(su2double T)const{
	return Cv*log(T);
}

su2double CPerfectGas::EntropyTempInv(su2double s)const{
	return exp(s/Cv);
}

su2double CPerfectGas::EnthalpyInv(su2double h)const{
	return h/(Cv + Gas_Constant);
}


