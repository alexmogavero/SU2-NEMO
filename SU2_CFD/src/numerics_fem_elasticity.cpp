/*!
 * \file numerics_fem_linear_elasticity.cpp
 * \brief This file contains the routines for setting the tangent matrix and residual of a FEM linear elastic structural problem.
 * \author R. Sanchez
 * \version 4.0.0 "Cardinal"
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

#include "../include/numerics_structure.hpp"
#include <limits>

CFEM_Elasticity::CFEM_Elasticity(unsigned short val_nDim, unsigned short val_nVar,
                                   CConfig *config) : CNumerics(val_nDim, val_nVar, config) {

	E = config->GetElasticyMod();
	Nu = config->GetPoissonRatio();
	Rho_s = config->GetMaterialDensity();
	Mu = E / (2.0*(1.0 + Nu));
	Lambda = Nu*E/((1.0+Nu)*(1.0-2.0*Nu));
	Kappa = config->GetBulk_Modulus_Struct();

	unsigned short i;

	KAux_ab = new double* [nDim];
	for (i = 0; i < nDim; i++) {
		KAux_ab[i] = new double[nDim];
	}


	if (nDim == 2){
		Ba_Mat = new double* [3];
		Bb_Mat = new double* [3];
		D_Mat  = new double* [3];
		GradNi_Ref_Mat = new double* [4];	/*--- As of now, 4 is the maximum number of nodes for 2D problems ---*/
		GradNi_Curr_Mat = new double* [4];	/*--- As of now, 4 is the maximum number of nodes for 2D problems ---*/
		for (i = 0; i < 3; i++) {
			Ba_Mat[i]  		= new double[nDim];
			Bb_Mat[i]  		= new double[nDim];
			D_Mat[i]	   	= new double[3];
		}
		for (i = 0; i < 4; i++) {
			GradNi_Ref_Mat[i] 	= new double[nDim];
			GradNi_Curr_Mat[i] 	= new double[nDim];
		}
	}
	else if (nDim == 3){
		Ba_Mat = new double* [6];
		Bb_Mat = new double* [6];
		D_Mat  = new double* [6];
		GradNi_Ref_Mat = new double* [8];	/*--- As of now, 8 is the maximum number of nodes for 3D problems ---*/
		GradNi_Curr_Mat = new double* [8];	/*--- As of now, 8 is the maximum number of nodes for 3D problems ---*/
		for (i = 0; i < 6; i++) {
			Ba_Mat[i]  		= new double[nDim];
			Bb_Mat[i]  		= new double[nDim];
			D_Mat[i]      	= new double[6];
		}
		for (i = 0; i < 8; i++) {
			GradNi_Ref_Mat[i] 	= new double[nDim];
			GradNi_Curr_Mat[i] 	= new double[nDim];
		}
	}
}

CFEM_Elasticity::~CFEM_Elasticity(void) {

	unsigned short iVar, jVar;

	for (iVar = 0; iVar < nDim; iVar++){
		delete [] KAux_ab[iVar];
	}

	if (nDim == 2){
		for (iVar = 0; iVar < 3; iVar++){
			delete [] Ba_Mat[iVar];
			delete [] Bb_Mat[iVar];
			delete [] D_Mat[iVar];
		}
		for (iVar = 0; iVar < 4; iVar++){
			delete [] GradNi_Ref_Mat[iVar];
			delete [] GradNi_Curr_Mat[iVar];
		}
	}
	else if (nDim == 3){
		for (iVar = 0; iVar < 6; iVar++){
			delete [] Ba_Mat[iVar];
			delete [] Bb_Mat[iVar];
			delete [] D_Mat[iVar];
		}
		for (iVar = 0; iVar < 8; iVar++){
			delete [] GradNi_Ref_Mat[iVar];
			delete [] GradNi_Curr_Mat[iVar];
		}
	}

	delete [] KAux_ab;
	delete [] Ba_Mat;
	delete [] Bb_Mat;
	delete [] D_Mat;
	delete [] GradNi_Ref_Mat;
	delete [] GradNi_Curr_Mat;

}
