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
#include "../include/CPerfectGas.hpp"
#include <functional>
#include <boost/math/tools/roots.hpp>

const su2double CMutationpp::Tmin = 55.0; //hard coded as it is hard coded in mutation++ it is 50 + 5 for tolerance

CMutationpp::CMutationpp() :
	CFluidModel(),
	opt("air5"),
	mix(opt),
	comp(),
	Gas_Constant(0),
	fallBackModel(),
	e_min(0){
}


CMutationpp::CMutationpp(string optFile, vector<double> cmp):
	CFluidModel(),
	opt(optFile),
	mix(opt),
	comp(cmp),
	Gas_Constant(0){
  Gas_Constant = CalcGasConstant();

  vector<double> rhoSpecie = SpecieDensity(1);
  mix.setState(rhoSpecie.data(), &Tmin, 1);
  su2double gamma = mix.mixtureFrozenCpMass()/mix.mixtureFrozenCvMass();

  fallBackModel = CPerfectGas(Gas_Constant, gamma);

  e_min = mix.mixtureEnergyMass();
  fallBackModel.SetTDState_rhoT(1, Tmin);
  form_e = e_min - fallBackModel.GetStaticEnergy();
  form_s = mix.mixtureSMass() - fallBackModel.GetEntropy();

  //Necessary to avoid error in mutation++ trying to do log of 0
  //TODO edit mutation++ to handle this
  for(size_t i=0; i<comp.size(); i++){
    if(comp[i]==0){
      comp[i] = 1e-20;
    }
  }
}

su2double CMutationpp::CalcGasConstant()const{
  su2double M = 0;
  for(size_t i=0; i<comp.size(); i++){
    M += comp[i]*mix.speciesMw(i);
  }

  return Mutation::RU/M;
}


CMutationpp::~CMutationpp(void) {

}

void CMutationpp::SetTDState_rhoe (su2double rho, su2double e ) {
  if(rho < 0 ){
    SetWrongState();
    return;
  }

  if(e < e_min){
    fallBackModel.SetTDState_rhoe(rho, e);
    UpdateFromFallBack();
    return;
  }

  vector<double> rhoSpecie = SpecieDensity(rho);
  su2double rhoe = rho*e;
  mix.setState(rhoSpecie.data(), &rhoe, 0);

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

void CMutationpp::SetWrongState(){
  Density = -1;
  StaticEnergy = 0; //Non non-physical but probably ok because the others are bad
  Temperature = -1;
  Pressure = -1;

  Cp = -1;

  SoundSpeed2 = -1;

  Entropy = 0; //Non non-physical but probably ok because the others are bad

  dPdrho_e = -1;
  dTdrho_e = 0.0;
  dTde_rho = -1;
  dPde_rho = -1;
}

void CMutationpp::UpdateFromFallBack(){
  Density = fallBackModel.GetDensity();
  StaticEnergy = fallBackModel.GetStaticEnergy() + form_e;
  Temperature = fallBackModel.GetTemperature();
  Pressure = fallBackModel.GetPressure();

  Cp = fallBackModel.GetCp();

  SoundSpeed2 = fallBackModel.GetSoundSpeed2();

  Entropy = fallBackModel.GetEntropy() + form_s;

  dPdrho_e = fallBackModel.GetdPdrho_e();
  dTdrho_e = fallBackModel.GetdTdrho_e();
  dTde_rho = fallBackModel.GetdTde_rho();
  dPde_rho = fallBackModel.GetdPde_rho();
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
  if(rho < 0 ){
    SetWrongState();
    return;
  }

  if(T < Tmin){
    fallBackModel.SetTDState_rhoT(rho, T);
    UpdateFromFallBack();
    return;
  }

  vector<double> rhoSpecie = SpecieDensity(rho);
  mix.setState(rhoSpecie.data(), &T, 1);

  UpdateState();
}

void CMutationpp::SetEnergy_Prho (su2double P, su2double rho ) {
  SetTDState_Prho(P, rho);
}

vector<su2double> CMutationpp::GetMaxwellMoment(unsigned short nDim)const{
  return mix.mixtureFrozenMaxwellMoment(3 - nDim);
}



