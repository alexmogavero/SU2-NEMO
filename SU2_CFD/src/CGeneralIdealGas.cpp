/*!
 * CGeneralIdealGas.cpp
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

#include "../include/CGeneralIdealGas.hpp"

CGeneralIdealGas::CGeneralIdealGas() :
	CFluidModel(),
	Gas_Constant(0){
}


CGeneralIdealGas::CGeneralIdealGas(su2double R):
	CFluidModel(),
	Gas_Constant(R){
}


CGeneralIdealGas::~CGeneralIdealGas(void) {

}

void CGeneralIdealGas::SetTDState_rhoe (su2double rho, su2double e ) {
  
  Density = rho;
  StaticEnergy = e;
  Temperature = EnergyInv(e);
  Pressure = Density*Gas_Constant*Temperature;

  su2double Cv = SpecificHeatVol(Temperature);
  Cp = Cv + Gas_Constant;
  su2double Gamma = Cp/Cv;

  SoundSpeed2 = Gamma*Pressure/Density;

  Entropy = EntropyTemp(Temperature) + log(1.0/Density)*Gas_Constant;

  dPdrho_e = Gas_Constant*Temperature;
  dTdrho_e = 0.0;
  dTde_rho = 1/Cv;
  dPde_rho = Gas_Constant*Density*dTde_rho;

}

void CGeneralIdealGas::SetTDState_PT (su2double P, su2double T ) {
  su2double e = Energy(T);
  su2double rho = P/(T*Gas_Constant);
  SetTDState_rhoe(rho, e);
}

void CGeneralIdealGas::SetTDState_Prho (su2double P, su2double rho ) {
  su2double T = P/(Gas_Constant*rho);
  SetTDState_PT(P, T);

}

void CGeneralIdealGas::SetEnergy_Prho (su2double P, su2double rho ) {
	su2double T = P/(Gas_Constant*rho);
	StaticEnergy = Energy(T);

}

void CGeneralIdealGas::SetTDState_hs (su2double h, su2double s ) {
  su2double T = EnthalpyInv(h);
  su2double e = h - Gas_Constant*T;
  su2double v = exp((s - EntropyTemp(Temperature))/Gas_Constant);

  SetTDState_rhoe(1/v, e);
}

void CGeneralIdealGas::SetTDState_Ps (su2double P, su2double s ) {
	throw std::logic_error("Error: method SetTDState_Ps not implemented in class CGeneralIdealGas.");
}

void CGeneralIdealGas::SetTDState_rhoT (su2double rho, su2double T ) {
  su2double e = Energy(T);
  SetTDState_rhoe(rho, e);
}









