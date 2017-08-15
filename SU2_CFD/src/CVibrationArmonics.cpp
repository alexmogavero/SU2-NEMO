/*!
 * CVibrationArmonics.cpp
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

#include "../include/CVibrationArmonics.hpp"

CVibrationArmonics::CVibrationArmonics() :
	CPerfectGas(),
	ThetaVib(),
	weight(){
}


CVibrationArmonics::CVibrationArmonics(su2double R, su2double g,
		unsigned short n_modes, su2double* theta, su2double* w):
		CPerfectGas(R, g),
	ThetaVib(n_modes),
	weight(n_modes){
	for(unsigned short i=0; i<n_modes; i++){
		ThetaVib[i] = theta[i];
	}
	for(unsigned short i=0; i<n_modes; i++){
		weight[i] = w[i];
	}
}

CVibrationArmonics::CVibrationArmonics(su2double R, su2double g,
		std::vector<su2double> theta, std::vector<su2double> w):
		CPerfectGas(R, g),
	ThetaVib(theta),
	weight(w){
	if(theta.size() != weight.size()) throw std::logic_error("Error: theta and w must be of the same size.");
}


CVibrationArmonics::~CVibrationArmonics(void) {

}

su2double CVibrationArmonics::Energy(su2double T)const{
	su2double C = 0;
	for(std::size_t i=0; i<ThetaVib.size(); i++){
		C += weight[i]*ThetaVib[i]/(exp(ThetaVib[i]/T)-1);
	}
	return Gas_Constant*C + CPerfectGas::Energy(T);
}

su2double CVibrationArmonics::EnergyInv(su2double e)const{
	return CGeneralIdealGas::EnergyInv(e);
}

su2double CVibrationArmonics::SpecificHeatVol(su2double T)const{
	su2double cv = CPerfectGas::SpecificHeatVol(T);

	su2double C = 0;
	for(std::size_t i=0; i<ThetaVib.size(); i++){
		su2double theta_T = ThetaVib[i]/T;
		su2double exp_theta_T = exp(theta_T);
		C += weight[i]*pow(theta_T, 2)*exp_theta_T/pow(exp_theta_T - 1, 2);
	}
	return Gas_Constant*C + cv;
}

su2double CVibrationArmonics::EntropyTemp(su2double T)const{
	su2double C = 0;
	for(std::size_t i=0; i<ThetaVib.size(); i++){
		su2double e_T = exp(ThetaVib[i]/T);
		C += weight[i]*((ThetaVib[i] + ThetaVib[i]/(e_T - 1))*(1/T - 1) +
				log(abs(exp(ThetaVib[i]) - 1)) - log(abs(e_T - 1)));
	}

	return Gas_Constant*C +	CPerfectGas::EntropyTemp(T);
}

su2double CVibrationArmonics::EntropyTempInv(su2double s)const{
	return CGeneralIdealGas::EntropyTempInv(s);
}

su2double CVibrationArmonics::EnthalpyInv(su2double h)const{
	return CGeneralIdealGas::EnthalpyInv(h);
}

