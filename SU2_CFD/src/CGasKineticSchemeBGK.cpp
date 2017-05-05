#include "../include/CGasKineticSchemeBGK.hpp"
#include "../include/CKineticVariable.hpp"

CGasKineticSchemeBGK::CGasKineticSchemeBGK(unsigned short val_nDim, unsigned short val_nVar, CConfig *config):
  CNumerics(val_nDim, val_nVar, config),
  FluidModel(NULL),
  node_i(NULL),
  node_j(NULL),
  node_I(NULL),
  config(config){

  Gamma = config->GetGamma();
  Gamma_Minus_One = Gamma - 1.0;

  implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  grid_movement = config->GetGrid_Movement();

  /*--- Artifical dissipation part ---*/
  Param_p = 0.3;
  Param_Kappa_0 = config->GetKappa_1st_Flow();

  /*--- Allocate some structures ---*/
  Diff_U = new su2double [nVar];
  Velocity_i = new su2double [nDim];
  Velocity_j = new su2double [nDim];
  MeanVelocity = new su2double [nDim];
  ProjFlux = new su2double [nVar];

  U_i = new su2double [nVar];
  U_j = new su2double [nVar];
  U_I = new su2double [nVar];
}

CGasKineticSchemeBGK::~CGasKineticSchemeBGK(void) {
  delete [] Diff_U;
  delete [] Velocity_i;
  delete [] Velocity_j;
  delete [] MeanVelocity;
  delete [] ProjFlux;

  delete [] U_i;
  delete [] U_j;

  if(node_I){
    delete node_I;
  }
}

void CGasKineticSchemeBGK::ComputeResidual(su2double *val_residual, CConfig *config){
  /*--- Pressure, density, enthalpy, energy, and velocity at points i and j ---*/
  Pressure_i = V_i[nDim+1];                       Pressure_j = V_j[nDim+1];
  Density_i = V_i[nDim+2];                        Density_j = V_j[nDim+2];
  Enthalpy_i = V_i[nDim+3];                       Enthalpy_j = V_j[nDim+3];
  Energy_i = Enthalpy_i - Pressure_i/Density_i;   Energy_j = Enthalpy_j - Pressure_j/Density_j;

  for (unsigned short iDim = 0; iDim < nDim; iDim++) {
    Velocity_i[iDim] = V_i[iDim+1];
    Velocity_j[iDim] = V_j[iDim+1];
  }

  /*--- Recompute conservative variables ---*/
  // TODO maybe this step can be bypassed because conserved quantities are calculated already by the solver
  U_i[0] = Density_i; U_j[0] = Density_j;
  for (unsigned short iDim = 0; iDim < nDim; iDim++) {
    U_i[iDim+1] = Density_i*Velocity_i[iDim];
    U_j[iDim+1] = Density_j*Velocity_j[iDim];
  }
  U_i[nDim+1] = Density_i*Energy_i;
  U_j[nDim+1] = Density_j*Energy_j;

  CalculateInterface();

  //theta calculation only valid for ideal gas
  //TODO move it to gas model class
  su2double theta_i = Density_i/(2*Pressure_i);
  su2double theta_j = Density_j/(2*Pressure_j);

  su2double Velocity2_I = 0;
  for(unsigned short iDim = 0; iDim<nDim; iDim++){
    Velocity2_I += pow(U_I[iDim+1]/U_I[0], 2);
  }
  su2double Energy_I = U_I[nVar-1] - 0.5*Velocity2_I;
  FluidModel->SetTDState_rhoe(U_I[0], Energy_I);
  su2double Pressure_I = FluidModel->GetPressure();
  su2double theta_I = U_I[0]/(2*Pressure_I);

  //Calculate the mean collision time
  //TODO check if it is ok to calculate it on the interface
  su2double tauColl = FluidModel->GetLaminarViscosity()/Pressure_I;

  su2double Flux_i[nVar], Flux_j[nVar], Flux_I[nVar];
  std::vector<std::vector<unsigned short> > exponents(nVar,std::vector<unsigned short>(nVar-1,0));
  exponents[0][0] = 1;
  for(unsigned short iDim=0; iDim<nDim; iDim++){
    exponents[iDim+1][iDim] = 1;
    exponents[iDim+1][0] += 1;
    exponents[nVar-1][iDim] = 2;
  }
  exponents[nVar-1][0] += 1;
  exponents[nVar-1][nVar-2] = 2;

  for(unsigned short iVar=0; iVar<nVar; iVar++){
    Flux_I[iVar] = MomentsMaxwellian(exponents[iVar], theta_I, ALL);
    Flux_i[iVar] = MomentsMaxwellian(exponents[iVar], theta_i, POSITIVE);
    Flux_j[iVar] = MomentsMaxwellian(exponents[iVar], theta_j, NEGATIVE);
  }
}

void CGasKineticSchemeBGK::CalculateInterface(){
  std::vector<su2double> U_L, U_R, U_I;

  U_L = PsiMaxwell(LEFT, POSITIVE);
  U_R = PsiMaxwell(LEFT, POSITIVE);

  for(std::size_t i=0; i<U_L.size(); i++){
    U_I[i] =  U_L[i] + U_R[i];
  }

  su2double vel[nDim];
  for(unsigned short iDim=0; iDim<nDim; iDim++){
    vel[iDim] = U_I[iDim+1]/U_I[0];
  }
  node_I = new CKineticVariable(U_I.data(), nDim, nVar, config);

  node_I->SetNon_Physical(false);

  bool RightSol = node_I->SetPrimVar(FluidModel);
  node_I->SetSecondaryVar(FluidModel);

  if (!RightSol) {
    node_I->SetNon_Physical(true);
  }
}

std::vector<su2double> CGasKineticSchemeBGK::PsiMaxwell(State state, IntLimits lim, bool uPsi)const{
  std::vector<su2double> out(nVar, 0);
  std::vector<unsigned short> exponents(nVar-1,0);

  if(uPsi) exponents[0]++;
  out[0] = MomentsMaxwellian(exponents, state, lim);

  for(unsigned short iDim=0; iDim<nDim; iDim++){
    exponents.assign(nVar-1, 0);
    exponents[iDim] = 1;
    if(uPsi) exponents[0]++;
    out[iDim+1] = MomentsMaxwellian(exponents, state, lim);
  }

  for(unsigned short iDim=0; iDim<nDim; iDim++){
    exponents.assign(nVar-1, 0);
    exponents[iDim] = 2;
    if(uPsi) exponents[0]++;
    out[nVar-1] += MomentsMaxwellian(exponents, state, lim);
  }
  exponents.assign(nVar-1, 0);
  exponents[nVar-2] = 2;
  if(uPsi) exponents[0]++;
  out[nVar-1] += MomentsMaxwellian(exponents, state, lim);
  out[nVar-1] /= 2;

  return out;
}

void CGasKineticSchemeBGK::ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j,
                                    CConfig *config) {

  //temporary
  unsigned short iDim, iVar, jVar;

  su2double U_i[5] = {0.0,0.0,0.0,0.0,0.0}, U_j[5] = {0.0,0.0,0.0,0.0,0.0};

  /*--- Pressure, density, enthalpy, energy, and velocity at points i and j ---*/

  Pressure_i = V_i[nDim+1];                       Pressure_j = V_j[nDim+1];
  Density_i = V_i[nDim+2];                        Density_j = V_j[nDim+2];
  Enthalpy_i = V_i[nDim+3];                       Enthalpy_j = V_j[nDim+3];
  SoundSpeed_i = V_i[nDim+4];                     SoundSpeed_j = V_j[nDim+4];
  Energy_i = Enthalpy_i - Pressure_i/Density_i;   Energy_j = Enthalpy_j - Pressure_j/Density_j;

  sq_vel_i = 0.0; sq_vel_j = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Velocity_i[iDim] = V_i[iDim+1];
    Velocity_j[iDim] = V_j[iDim+1];
    sq_vel_i += 0.5*Velocity_i[iDim]*Velocity_i[iDim];
    sq_vel_j += 0.5*Velocity_j[iDim]*Velocity_j[iDim];
  }

  /*--- Recompute conservative variables ---*/

  U_i[0] = Density_i; U_j[0] = Density_j;
  for (iDim = 0; iDim < nDim; iDim++) {
    U_i[iDim+1] = Density_i*Velocity_i[iDim]; U_j[iDim+1] = Density_j*Velocity_j[iDim];
  }
  U_i[nDim+1] = Density_i*Energy_i; U_j[nDim+1] = Density_j*Energy_j;

  /*--- Compute mean values of the variables ---*/

  MeanDensity = 0.5*(Density_i+Density_j);
  MeanPressure = 0.5*(Pressure_i+Pressure_j);
  MeanEnthalpy = 0.5*(Enthalpy_i+Enthalpy_j);
  for (iDim = 0; iDim < nDim; iDim++)
    MeanVelocity[iDim] =  0.5*(Velocity_i[iDim]+Velocity_j[iDim]);
  MeanEnergy = 0.5*(Energy_i+Energy_j);

  /*--- Get projected flux tensor ---*/

  GetInviscidProjFlux(&MeanDensity, MeanVelocity, &MeanPressure, &MeanEnthalpy, Normal, ProjFlux);

  /*--- Residual of the inviscid flux ---*/

  for (iVar = 0; iVar < nVar; iVar++)
    val_residual[iVar] = ProjFlux[iVar];

  /*--- Jacobians of the inviscid flux, scale = 0.5 because val_residual ~ 0.5*(fc_i+fc_j)*Normal ---*/

  if (implicit) {
    GetInviscidProjJac(MeanVelocity, &MeanEnergy, Normal, 0.5, val_Jacobian_i);
    for (iVar = 0; iVar < nVar; iVar++)
      for (jVar = 0; jVar < nVar; jVar++)
        val_Jacobian_j[iVar][jVar] = val_Jacobian_i[iVar][jVar];
  }

  /*--- Adjustment due to grid motion ---*/

  if (grid_movement) {
    ProjVelocity = 0.0;
    for (iDim = 0; iDim < nDim; iDim++)
      ProjVelocity += 0.5*(GridVel_i[iDim]+GridVel_j[iDim])*Normal[iDim];
    for (iVar = 0; iVar < nVar; iVar++) {
      val_residual[iVar] -= ProjVelocity * 0.5*(U_i[iVar]+U_j[iVar]);
      if (implicit) {
        val_Jacobian_i[iVar][iVar] -= 0.5*ProjVelocity;
        val_Jacobian_j[iVar][iVar] -= 0.5*ProjVelocity;
      }
    }
  }

  /*--- Computes differences btw. conservative variables,
   with a correction for the enthalpy ---*/

  for (iVar = 0; iVar < nDim+1; iVar++)
    Diff_U[iVar] = U_i[iVar]-U_j[iVar];
  Diff_U[nDim+1] = Density_i*Enthalpy_i-Density_j*Enthalpy_j;

  /*--- Compute the local spectral radius and the stretching factor ---*/

  ProjVelocity_i = 0.0; ProjVelocity_j = 0.0; Area = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    ProjVelocity_i += Velocity_i[iDim]*Normal[iDim];
    ProjVelocity_j += Velocity_j[iDim]*Normal[iDim];
    Area += Normal[iDim]*Normal[iDim];
  }
  Area = sqrt(Area);

  /*--- Adjustment due to grid motion ---*/
  if (grid_movement) {
    ProjGridVel = 0.0;
    for (iDim = 0; iDim < nDim; iDim++)
      ProjGridVel += 0.5*(GridVel_i[iDim]+GridVel_j[iDim])*Normal[iDim];
    ProjVelocity_i -= ProjGridVel;
    ProjVelocity_j -= ProjGridVel;
  }

  Local_Lambda_i = (fabs(ProjVelocity_i)+SoundSpeed_i*Area);
  Local_Lambda_j = (fabs(ProjVelocity_j)+SoundSpeed_j*Area);
  MeanLambda = 0.5*(Local_Lambda_i+Local_Lambda_j);

  Phi_i = pow(Lambda_i/(4.0*MeanLambda), Param_p);
  Phi_j = pow(Lambda_j/(4.0*MeanLambda), Param_p);
  StretchingFactor = 4.0*Phi_i*Phi_j/(Phi_i+Phi_j);

  sc0 = 3.0*(su2double(Neighbor_i)+su2double(Neighbor_j))/(su2double(Neighbor_i)*su2double(Neighbor_j));
  Epsilon_0 = Param_Kappa_0*sc0*su2double(nDim)/3.0;

  /*--- Compute viscous part of the residual ---*/

  for (iVar = 0; iVar < nVar; iVar++)
    val_residual[iVar] += Epsilon_0*Diff_U[iVar]*StretchingFactor*MeanLambda;

  /*--- Jacobian computation ---*/

  if (implicit) {
    cte = Epsilon_0*StretchingFactor*MeanLambda;
    for (iVar = 0; iVar < (nVar-1); iVar++) {
      val_Jacobian_i[iVar][iVar] += cte;
      val_Jacobian_j[iVar][iVar] -= cte;
    }

    /*--- Last row of Jacobian_i ---*/

    val_Jacobian_i[nVar-1][0] += cte*Gamma_Minus_One*sq_vel_i;
    for (iDim = 0; iDim < nDim; iDim++)
      val_Jacobian_i[nVar-1][iDim+1] -= cte*Gamma_Minus_One*Velocity_i[iDim];
    val_Jacobian_i[nVar-1][nVar-1] += cte*Gamma;

    /*--- Last row of Jacobian_j ---*/

    val_Jacobian_j[nVar-1][0] -= cte*Gamma_Minus_One*sq_vel_j;
    for (iDim = 0; iDim < nDim; iDim++)
      val_Jacobian_j[nVar-1][iDim+1] += cte*Gamma_Minus_One*Velocity_j[iDim];
    val_Jacobian_j[nVar-1][nVar-1] -= cte*Gamma;

  }

}
