/*!
 * CMutationpp.cpp
 * \brief Source of mutation++ wrapper.
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

#include "../include/CMutationpp.hpp"
#include "../include/fluid_model.hpp"
#include <functional>
#include <boost/math/tools/roots.hpp>

CMutationpp::CMutationpp() :
	CFluidModel(),
	opt("air5"),
	mix(opt),
	comp(),
	Gas_Constant(0){
}


CMutationpp::CMutationpp(string optFile, vector<double> cmp):
	CFluidModel(),
	opt(optFile),
	mix(opt),
	comp(cmp),
	Gas_Constant(0){
  Gas_Constant = CalcGasConstant();
}

su2double CMutationpp::CalcGasConstant()const{
  su2double R = 0;
  for(size_t i=0; i<comp.size(); i++){
    R += comp[i]*mix.speciesMw(i);
  }

  return R;
}


CMutationpp::~CMutationpp(void) {

}

void CMutationpp::SetTDState_rhoe (su2double rho, su2double e ) {
  vector<double> rhoSpecie = SpecieDensity(rho);
  mix.setState(rhoSpecie.data(), &e, 0);

  UpdateState();
}

vector<double> CMutationpp::SpecieDensity(su2double rho)const{
  vector<double> rhoSpecie(comp.size(), 0);
  mix.convert<Mutation::Thermodynamics::CONC_TO_RHO>(comp.data(), rhoSpecie.data());
  double rhoFact = 0;
  for(size_t i=0; i<comp.size(); i++){
    rhoFact += rhoSpecie[i];
  }
  rhoFact = rho/rhoFact;
  for(size_t i=0; i<comp.size(); i++){
    rhoSpecie[i] *= rhoFact;
  }

  return rhoSpecie;
}

void CMutationpp::UpdateState(){
  Density = mix.density();
  StaticEnergy = mix.mixtureEnergyMass();
  Temperature = mix.T();
  Pressure = mix.P();

  Cp = mix.mixtureFrozenCpMass();

  SoundSpeed2 = pow(mix.frozenSoundSpeed(), 2);

  Entropy = mix.mixtureSMass();

  dPdrho_e = Pressure/Density;
  dTdrho_e = 0.0;
  dTde_rho = 1/mix.mixtureFrozenCvMass();
  dPde_rho = Pressure/Temperature*dTde_rho;
}

void CMutationpp::SetTDState_PT (su2double P, su2double T){
  su2double rho = P/(T*Gas_Constant);
  SetTDState_rhoT(rho, T);
}

void CMutationpp::SetTDState_Prho (su2double P, su2double rho ) {
  su2double T = P/(Gas_Constant*rho);
  SetTDState_rhoT(rho, T);

}

void CMutationpp::SetTDState_hs (su2double h, su2double s ) {
  throw std::logic_error("Error: method SetTDState_hs not implemented in class CMutationpp.");
}

void CMutationpp::SetTDState_Ps (su2double P, su2double s ) {
	throw std::logic_error("Error: method SetTDState_Ps not implemented in class CMutationpp.");
}

void CMutationpp::SetTDState_rhoT (su2double rho, su2double T ) {
  vector<double> rhoSpecie = SpecieDensity(rho);
  mix.setState(rhoSpecie.data(), &T, 1);

  UpdateState();
}

void CMutationpp::SetEnergy_Prho (su2double P, su2double rho ) {
  SetTDState_Prho(P, rho);
}



