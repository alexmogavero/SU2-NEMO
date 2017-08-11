/*!
 * CElectronicExcitation
 * \brief Class that model a Maxwell-Boltzmann gas with electronic excitation model for specific heat.
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

#include "../include/CElectronicExcitation.hpp"
#include <functional>
#include <boost/numeric/odeint.hpp>

CElectronicExcitation::CElectronicExcitation() :
	CVibrationArmonics(){
}


CElectronicExcitation::CElectronicExcitation(su2double R, su2double g, su2double theta,
		std::vector<su2double> gEl, std::vector<su2double> thetaEl):
		CVibrationArmonics(R, g, theta),
		g(gEl),
		ThetaEl(thetaEl){
	if(gEl.size() != thetaEl.size()) throw std::logic_error("Error: gEl and thetaEl must be of the same size.");
}


CElectronicExcitation::~CElectronicExcitation(void) {

}

CElectronicExcitation::Int_Cv::Int_Cv(const CElectronicExcitation* f):
		obj(f){
}

void CElectronicExcitation::Int_Cv::operator() (const su2double &x , su2double &dxdt , const double T){
	dxdt = obj->Cv(T);
}

su2double CElectronicExcitation::Energy(su2double T)const{
	su2double e = 0;
	Int_Cv integ(this);
	std::size_t steps = boost::numeric::odeint::integrate<su2double, Int_Cv, su2double, su2double>(
			integ, e, 0.0, T, 0.1);

	return e;
}

su2double CElectronicExcitation::Cv(su2double T)const{
	su2double Z = 0;
	su2double A = 0;
	su2double B = 0;
	su2double C = 0;
	for(std::size_t i=0; i<g.size(); i++){
		Z += g[i]*exp(-ThetaEl[i]/T);
		A += g[i]*pow(ThetaEl[i]/T, 2.0)*exp(-ThetaEl[i]/T);
		B += g[i]*ThetaEl[i]*exp(-ThetaEl[i]/T);
		C += g[i]*(ThetaEl[i]/pow(T, 2.0))*exp(-ThetaEl[i]/T);
	}

	return Gas_Constant*(A/Z - B*C/pow(Z, 2.0));
}

su2double CElectronicExcitation::SpecificHeatVol(su2double T)const{
	return Cv(T) + CVibrationArmonics::SpecificHeatVol(T);
}

su2double CElectronicExcitation::EntropyTemp(su2double T)const{
	//TODO add implementation
}


