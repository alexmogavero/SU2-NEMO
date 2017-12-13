#define _USE_MATH_DEFINES
#include <cmath>
#include "lapacke.h"
#include "../include/CGasKineticSchemeBGK.hpp"
#include "../include/CKineticVariable.hpp"

CGasKineticSchemeBGK::CGasKineticSchemeBGK(unsigned short val_nDim, unsigned short val_nVar, CConfig *config):
  CNumerics(val_nDim, val_nVar, config),
  node_I(NULL),
  node_iLoc(NULL),
  node_jLoc(NULL),
  node_iRot(NULL),
  node_jRot(NULL),
  rotMatrix(),
  config(config){
  su2double U[nVar];
  for(unsigned short i=0; i<nVar; i++){
    U[i] = 0;
  }
  node_I = new CKineticVariable(U, nDim, nVar, config);
}

CGasKineticSchemeBGK::~CGasKineticSchemeBGK(void) {
  if(node_I){
    delete node_I;
  }
  if(node_iLoc){
    delete node_iLoc;
  }
  if(node_jLoc){
    delete node_jLoc;
  }
  if(node_iRot){
    delete node_iRot;
  }
  if(node_jRot){
    delete node_jRot;
  }
}

void CGasKineticSchemeBGK::ComputeResidual(su2double *val_residual, CConfig *config){
  Clear();

  // Convert limited primitive variables to conservative
  su2double Pressure_i, Pressure_j, Density_i, Density_j, Enthalpy_i, Enthalpy_j, Energy_i, Energy_j;
  su2double Velocity_i[nDim], Velocity_j[nDim];
  unsigned short iDim, iVar;

  su2double U_i[5] = {0.0,0.0,0.0,0.0,0.0}, U_j[5] = {0.0,0.0,0.0,0.0,0.0};

  Pressure_i = V_i[nDim+1];                       Pressure_j = V_j[nDim+1];
  Density_i = V_i[nDim+2];                        Density_j = V_j[nDim+2];
  Enthalpy_i = V_i[nDim+3];                       Enthalpy_j = V_j[nDim+3];
  Energy_i = Enthalpy_i - Pressure_i/Density_i;   Energy_j = Enthalpy_j - Pressure_j/Density_j;

  for (iDim = 0; iDim < nDim; iDim++) {
    Velocity_i[iDim] = V_i[iDim+1];
    Velocity_j[iDim] = V_j[iDim+1];
  }

  U_i[0] = Density_i; U_j[0] = Density_j;
  for (iDim = 0; iDim < nDim; iDim++) {
    U_i[iDim+1] = Density_i*Velocity_i[iDim]; U_j[iDim+1] = Density_j*Velocity_j[iDim];
  }
  U_i[nDim+1] = Density_i*Energy_i; U_j[nDim+1] = Density_j*Energy_j;


  //Rotate Reference Frame
  node_iRot = node_i->duplicate();
  rotate(node_iRot);

  node_iLoc = node_i->duplicate();
  node_iLoc->SetSolution(U_i);

  node_iLoc->SetNon_Physical(false);
  bool RightSol = node_iLoc->SetPrimVar(FluidModel);
  if (!RightSol) {
    node_iLoc->SetNon_Physical(true);
  }
  rotate(node_iLoc);

  node_jRot = node_j->duplicate();
  rotate(node_jRot);

  node_jLoc = node_j->duplicate();
  node_jLoc->SetSolution(U_j);

  node_jLoc->SetNon_Physical(false);
  RightSol = node_jLoc->SetPrimVar(FluidModel);
  if (!RightSol) {
    node_jLoc->SetNon_Physical(true);
  }
  rotate(node_jLoc);

  CalculateInterface();

  ComputeFluxes( true, val_residual);

  if (abs(config->GetPrandtl_Lam()-1.0) > 1.e-12) {
    su2double val_residual_zero[nVar];
    ComputeFluxes( false, val_residual_zero);

    su2double q;
    su2double U = val_residual_zero[1]/val_residual_zero[0];
    q = val_residual[nVar-1] + 0.5 * pow(U,2) * val_residual[0]
    - U * val_residual[1] - U *  val_residual_zero[nVar -1]
    - 0.5 * pow(U,3) * val_residual_zero[0] + pow(U,2) * val_residual_zero[1];

    for(unsigned short i=1; i<nDim; i++){
      su2double V = val_residual_zero[i+1]/val_residual_zero[0];
      q += 0.5 * pow(V,2) * val_residual[0] - V * val_residual[i+1]
      -0.5 * U * pow(V,2) * val_residual_zero[0] + U*V*val_residual_zero[i+1];
    }

    val_residual[nVar-1] += (1/config->GetPrandtl_Lam() - 1) * q;
  }

  rotate(val_residual + 1, true);
}

void CGasKineticSchemeBGK::ComputeFluxes(bool order, su2double *val_residual){
  su2double Dt = min(node_i->GetDelta_Time(),node_j->GetDelta_Time());

  //Calculate the mean collision time
  su2double tauColl = node_I->GetLaminarViscosity()/node_I->GetPressure();
  double K = (5.0 - 3.0*Gamma) / (Gamma - 1.0) + (3.0 - nDim);
  double lL = (K+nDim) / (4.0*(node_iLoc->GetEnergy() - 0.5*node_iLoc->GetVelocity2()));
  double lR = (K+nDim) / (4.0*(node_jLoc->GetEnergy() - 0.5*node_jLoc->GetVelocity2()));
  double rho_lamL = node_iLoc->GetDensity()/lL;
  double rho_lamR = node_jLoc->GetDensity()/lR;
  tauColl += Dt*abs(rho_lamL - rho_lamR)/abs(rho_lamL + rho_lamR);

  std::vector<su2double> Flux_i, Flux_j, Flux_I;
  Flux_I = PsiMaxwell(INTERFACE, ALL, order);
  Flux_i = PsiMaxwell(LEFT, POSITIVE, order);
  Flux_j = PsiMaxwell(RIGHT, NEGATIVE, order);

  //calculate time integrals
  su2double e_dt = exp(-Dt/tauColl);
  su2double int_e_dt = tauColl*(1 - e_dt); //integral of exp(-Dt/tauColl)
  su2double int_1_e_dt = Dt - int_e_dt;
  su2double Dt_inv = 1/Dt;

  for(unsigned short iVar=0; iVar<nVar; iVar++){
    val_residual[iVar] = Dt_inv*(int_1_e_dt*Flux_I[iVar] + int_e_dt*(Flux_i[iVar] + Flux_j[iVar]))*Area;
  }

  if(config->GetSpatialOrder_Flow() == SECOND_ORDER ||
      config->GetSpatialOrder_Flow() == SECOND_ORDER_LIMITER){ // if second order
    std::vector<su2double> Flux_i, Flux_j, Flux_I, Flux_I_t;
    Flux_i = std::vector<su2double>(nVar, 0);
    Flux_j = std::vector<su2double>(nVar, 0);
    Flux_I = std::vector<su2double>(nVar, 0);

    // compute space derivatives at i and j
    std::vector<std::vector<su2double> > a_i(nDim, std::vector<su2double>(nVar, 0));
    std::vector<std::vector<su2double> > a_j(nDim, std::vector<su2double>(nVar, 0));
    std::vector<su2double> A_i(nVar, 0);
    std::vector<su2double> A_j(nVar, 0);

    if (order) {
      Derivatives(LEFT, a_i, A_i); // Dont need time derivatives
      Derivatives(RIGHT, a_j, A_j); // Dont need time derivatives
    } else {
      Derivatives(INTERFACE, a_i, A_i);
      a_j = a_i;
      A_j = A_i;
    }

    // compute space (and time) derivatives at the interface
    std::vector<su2double> ad_i(nVar, 0);
    std::vector<su2double> ad_j(nVar, 0);
    std::vector<std::vector<su2double> > ad(nDim-1, std::vector<su2double>(nVar, 0));
    std::vector<su2double> Ad(nVar, 0);

    std::vector<su2double> g(5, 0);
    g[0] = int_1_e_dt;
    g[1] = (e_dt - 1)/g[0];
    g[2] = (-Dt + 2*int_e_dt - Dt*e_dt)/g[0];
    g[3] = (1-e_dt)/g[0];
    g[4] = (Dt*e_dt - int_e_dt)/g[0];

    Interface_Derivatives(ad_i, ad_j, ad, Ad, a_i, a_j, g);

    // Assemble fluxes
    std::vector<unsigned short> exponents;
    exponents.assign(nVar-1, 0);
    exponents[0] += 1;
    if (order) exponents[0] += 1;
    Flux_I += ad_i*PsiPsiMaxwell(INTERFACE, POSITIVE, exponents);
    Flux_I += ad_j*PsiPsiMaxwell(INTERFACE, NEGATIVE, exponents);
    for(unsigned short i=1; i<nDim; i++){
      exponents.assign(nVar-1, 0);
      if (order) exponents[0]++;
      exponents[i]++;
      Flux_I += ad[i-1]*PsiPsiMaxwell(INTERFACE, ALL, exponents);
    }

    exponents.assign(nVar-1, 0);
    if (order) exponents[0]++;
    Flux_I_t = Ad*PsiPsiMaxwell(INTERFACE, ALL, exponents);

    for(unsigned short i=0; i<nDim; i++){
      exponents.assign(nVar-1, 0);
      if (order) exponents[0]++;
      exponents[i]++;
      Flux_i -= a_i[i]*PsiPsiMaxwell(LEFT, POSITIVE, exponents);
      Flux_j -= a_j[i]*PsiPsiMaxwell(RIGHT, NEGATIVE, exponents);
    }

    su2double int_t_e_dt = tauColl*(tauColl - (Dt + tauColl)*e_dt); //integral of t*exp(-t/tauColl)

    su2double int_A = int_t_e_dt - tauColl*int_1_e_dt; //integral of t*exp(-t/tauColl) - tauColl*(1-exp(-t/tauColl))
    su2double int_B = 0.5*pow(Dt,2) - tauColl*int_1_e_dt; //integral of t - tauColl*(1-exp(-t/tauColl))

    for(unsigned short iVar=0; iVar<nVar; iVar++){
      val_residual[iVar] += Dt_inv * int_A * Flux_I[iVar] * Area;
      val_residual[iVar] += Dt_inv * int_B * Flux_I_t[iVar] * Area;
      val_residual[iVar] += Dt_inv * int_t_e_dt *(Flux_i[iVar] + Flux_j[iVar]) * Area;
    }
  }

  if(config->GetViscous()){
    std::vector<std::vector<su2double> > a_i(nDim, std::vector<su2double>(nVar, 0)); //space derivatives
    std::vector<std::vector<su2double> > a_j(nDim, std::vector<su2double>(nVar, 0)); //space derivatives
    std::vector<su2double> A_i(nVar, 0); //Time derivatives
    std::vector<su2double> A_j(nVar, 0); //Time derivatives

    if (order) {
      Derivatives(LEFT, a_i, A_i);
      Derivatives(RIGHT, a_j, A_j);
    } else {
      Derivatives(INTERFACE, a_i, A_i);
      a_j = a_i;
      A_j = A_i;
    }

    Flux_i = std::vector<su2double>(nVar, 0);
    Flux_j = std::vector<su2double>(nVar, 0);
    std::vector<unsigned short> exponents;
    for(unsigned short i=0; i<nDim; i++){
      exponents.assign(nVar-1, 0);
      if (order) exponents[0]++;
      exponents[i]++;
      Flux_i += a_i[i]*PsiPsiMaxwell(LEFT, POSITIVE, exponents);
      Flux_j += a_j[i]*PsiPsiMaxwell(RIGHT, NEGATIVE, exponents);
    }
    exponents.assign(nVar-1, 0);
    if (order) exponents[0]++;
    Flux_i += A_i*PsiPsiMaxwell(LEFT, POSITIVE, exponents);
    Flux_j += A_j*PsiPsiMaxwell(RIGHT, NEGATIVE, exponents);

    for(unsigned short iVar=0; iVar<nVar; iVar++){
      val_residual[iVar] -= Dt_inv*tauColl*int_e_dt*(Flux_i[iVar] + Flux_j[iVar])*Area;
    }
  }
}

void CGasKineticSchemeBGK::CalculateInterface(){
  std::vector<su2double> U_L, U_R, U_I;

  U_L = PsiMaxwell(LEFT, POSITIVE);
  U_R = PsiMaxwell(RIGHT, NEGATIVE);

  U_I.resize(nVar);
  for(std::size_t i=0; i<U_L.size(); i++){
    U_I[i] =  U_L[i] + U_R[i];
  }

//  node_I = new CKineticVariable(U_I.data(), nDim, nVar, config);
  node_I->SetSolution(U_I.data());

  node_I->SetNon_Physical(false);

  bool RightSol = node_I->SetPrimVar(FluidModel);
  //node_I->SetSecondaryVar(FluidModel);

  if (!RightSol) {
    node_I->SetNon_Physical(true);
  }
}

std::vector<std::vector<su2double> > CGasKineticSchemeBGK::PsiPsiMaxwell(State state, IntLimits lim,
    std::vector<unsigned short> multipFactor){
  if(multipFactor.size()!=nVar-1) throw std::logic_error("Error: multipFactor must be of size nVar-1");

  std::vector<unsigned short> exponents(multipFactor);
  std::vector<std::vector<su2double> > out(nVar, std::vector<su2double>(nVar, 0));

  out[0] = PsiMaxwell(state, lim, exponents);

  for(unsigned short iDim=0; iDim<nDim; iDim++){
    exponents = multipFactor;
    exponents[iDim] += 1;
    out[iDim+1] = PsiMaxwell(state, lim, exponents);
  }

  for(unsigned short iDim=0; iDim<nDim; iDim++){
    exponents = multipFactor;
    exponents[iDim] += 2;
    out[nVar-1] += PsiMaxwell(state, lim, exponents);
  }
  exponents = multipFactor;
  exponents[nVar-2] += 2;
  out[nVar-1] += PsiMaxwell(state, lim, exponents);
  out[nVar-1] /= 2;

  return out;
}

std::vector<su2double> CGasKineticSchemeBGK::PsiMaxwell(State state, IntLimits lim,
    std::vector<unsigned short> multipFactor){
  std::vector<su2double> out(nVar, 0);

  if(multipFactor.size()!=nVar-1) throw std::logic_error("Error: multipFactor must be of size nVar-1");

  std::vector<unsigned short> exponents(multipFactor);

  out[0] = MomentsMaxwellian(exponents, state, lim);

  for(unsigned short iDim=0; iDim<nDim; iDim++){
    exponents = multipFactor;
    exponents[iDim] += 1;
    out[iDim+1] = MomentsMaxwellian(exponents, state, lim);
  }

  for(unsigned short iDim=0; iDim<nDim; iDim++){
    exponents = multipFactor;
    exponents[iDim] += 2;
    out[nVar-1] += MomentsMaxwellian(exponents, state, lim);
  }
  exponents = multipFactor;
  exponents[nVar-2] += 2;
  out[nVar-1] += MomentsMaxwellian(exponents, state, lim);
  out[nVar-1] /= 2;

  return out;
}

std::vector<su2double> CGasKineticSchemeBGK::PsiMaxwell(State state, IntLimits lim, bool uPsi){
  std::vector<unsigned short> mFactor(nVar-1, 0);
  if(uPsi) mFactor[0]++;

  return PsiMaxwell(state, lim, mFactor);
}

std::vector<su2double> CGasKineticSchemeBGK::PsiPsiMaxwell(State state, std::vector<unsigned short> exponents){
  return MatrixToVector(PsiPsiMaxwell(state, ALL, exponents));
}

std::vector<su2double> CGasKineticSchemeBGK::MatrixToVector(const std::vector<std::vector<su2double> >& mat){
  std::vector<su2double> out(mat.size()*mat[0].size(), 0);

  for(std::size_t i=0; i<mat.size(); i++){
    for(std::size_t j=0; j<mat[0].size(); j++){
      out[i*mat[0].size() + j] = mat[i][j];
    }
  }

  return out;
}

void CGasKineticSchemeBGK::Derivatives(State state, std::vector<std::vector<su2double> >& G, std::vector<su2double>& Ft){
  if(G.size() != nDim) throw std::logic_error("Error: G must be a matrix of size nDim x nVar.");
  if(G[0].size() != nVar) throw std::logic_error("Error: G must be a matrix of size nDim x nVar.");
  if(Ft.size() != nVar) throw std::logic_error("Error: Ft must be a vector of size nVar.");

  std::vector<su2double> M;
  std::vector<su2double> sysM;

  std::vector<unsigned short> exponents(nVar-1, 0);
  M = PsiPsiMaxwell(state, exponents);

  std::vector<su2double> out;

  CVariable* node;

  switch (state){
    case LEFT:
      node = node_iLoc;
      break;
    case RIGHT:
      node = node_jLoc;
      break;
    case INTERFACE:
      node = node_I;
      break;
  }

  // Space derivatives
  // std::vector<std::vector<su2double> > G(nDim, std::vector<su2double>(nVar, 0));
  for (unsigned int j=0; j<nDim; j++){
    for (unsigned int i=0; i<nVar; i++){
      G[j][i] = node->GetGradient(i, j);
    }

    int ipiv[nVar];
    sysM = M;
    LAPACKE_dsysv(LAPACK_COL_MAJOR, 'U', nVar, 1, sysM.data(), nVar, ipiv, G[j].data(), nVar);
  }

  //Time derivatives
  // std::vector<su2double> Ft(nVar, 0);
  for (unsigned int j=0; j<nDim; j++){
    exponents.assign(nVar-1,0);
    exponents[j] = 1;

    Ft -= G[j]*PsiPsiMaxwell(state, ALL, exponents);
  }

  int ipiv[nVar];
  sysM = M;
  LAPACKE_dsysv(LAPACK_COL_MAJOR, 'U', nVar, 1, sysM.data(), nVar, ipiv, Ft.data(), nVar);

  return;
}

void CGasKineticSchemeBGK::Interface_Derivatives(
  std::vector<su2double>& ad_i, std::vector<su2double>& ad_j,
  std::vector<std::vector<su2double> >& ad,
  std::vector<su2double>& Ad,
  std::vector<std::vector<su2double> >& a_i, std::vector<std::vector<su2double> >& a_j,
  std::vector<su2double>& g){

  std::vector<su2double> M;
  std::vector<su2double> sysM;

  std::vector<unsigned short> exponents(nVar-1, 0);
  M = PsiPsiMaxwell(INTERFACE, exponents);

  // Finite differences gradients MT C.19 (using interface at midpoint!!!)
  su2double dist_ij = 0.0;
  for (unsigned int i=0; i<nDim; i++)
    dist_ij += (Coord_j[i]-Coord_i[i])*(Coord_j[i]-Coord_i[i]);
  dist_ij = sqrt(dist_ij);

  for (unsigned int i=0; i<nVar; i++){
    //To calculate the finite differences the non reconstructed node has to be used
    ad_i[i] = (node_I->GetSolution(i)-node_iRot->GetSolution(i))/(dist_ij/2); //Finite differences gradient left
    ad_j[i] = (node_jRot->GetSolution(i)-node_I->GetSolution(i))/(dist_ij/2); //Finite differences gradient right
  }

  //\overline{a}_L MT C.22
  int ipiv[nVar];
  sysM = M;
  LAPACKE_dsysv(LAPACK_COL_MAJOR, 'U', nVar, 1, sysM.data(), nVar, ipiv, ad_i.data(), nVar);

  //\overline{a}_R MT C.23
  sysM = M;
  LAPACKE_dsysv(LAPACK_COL_MAJOR, 'U', nVar, 1, sysM.data(), nVar, ipiv, ad_j.data(), nVar);

  //\overline{b} \overline{c} MT C.25
  std::vector<std::vector<su2double> > psi_i = PsiPsiMaxwell(LEFT, POSITIVE, exponents);
  std::vector<std::vector<su2double> > psi_j = PsiPsiMaxwell(RIGHT, NEGATIVE, exponents);
  for (unsigned int i=1; i<nDim; i++){
    ad[i-1] = a_i[i] * psi_i +  a_j[i] * psi_j;

    sysM = M;
    LAPACKE_dsysv(LAPACK_COL_MAJOR, 'U', nVar, 1, sysM.data(), nVar, ipiv, ad[i-1].data(), nVar);
  }

  //\overline{A} VKI 4.37
  Ad = g[1]*PsiMaxwell(INTERFACE, ALL, exponents);

  exponents.assign(nVar-1, 0);
  exponents[0]++;
  Ad += g[2]*(ad_i*PsiPsiMaxwell(INTERFACE, POSITIVE, exponents) + ad_j*PsiPsiMaxwell(INTERFACE, NEGATIVE, exponents));

  for (unsigned int i=1; i<nDim; i++){
    exponents.assign(nVar-1, 0);
    exponents[i]++;
    Ad += g[2]*ad[i-1]*PsiPsiMaxwell(INTERFACE, ALL, exponents);
  }

  exponents.assign(nVar-1, 0);
  Ad += g[3]*(PsiMaxwell(LEFT, POSITIVE, exponents) + PsiMaxwell(RIGHT, NEGATIVE, exponents));

  for (unsigned int i=0; i<nDim; i++){
    exponents.assign(nVar-1, 0);
    exponents[i]++;
    Ad += g[4]*a_i[i]*PsiPsiMaxwell(LEFT, POSITIVE, exponents);
    Ad += g[4]*a_j[i]*PsiPsiMaxwell(RIGHT, NEGATIVE, exponents);
  }

  sysM = M;
  LAPACKE_dsysv(LAPACK_COL_MAJOR, 'U', nVar, 1, sysM.data(), nVar, ipiv, Ad.data(), nVar);

}

su2double CGasKineticSchemeBGK::MomentsMaxwellian(std::vector<unsigned short> exponents, State state, IntLimits lim){
  su2double mp, rho;

  moments_struct* moments;
  CVariable* node;

  switch (state){
    case LEFT:
      node = node_iLoc;
      moments = &moments_i;
    break;
    case RIGHT:
      node = node_jLoc;
      moments = &moments_j;
    break;
    case INTERFACE:
      node = node_I;
      moments = &moments_I;
    break;
  }

  if (moments->xi.empty()) CGasKineticSchemeBGK::ComputeMaxwellianMoments(node, moments);

  mp = 1.0;
  switch (lim) {
    case ALL:
      mp *= moments->A[0][exponents[0]];
      break;
    case NEGATIVE:
      mp *= moments->N[0][exponents[0]];
      break;
    case POSITIVE:
      mp *= moments->P[0][exponents[0]];
      break;
  }

  for (unsigned short i = 1; i<nDim; i++){
    mp *= moments->A[i][exponents[i]];
  }
  mp *= moments->xi[exponents[nDim]];

  rho = node->GetDensity();
  return mp*rho;
}

void CGasKineticSchemeBGK::ComputeMaxwellianMoments(CVariable* node, moments_struct*  moments){
  double K = (5.0 - 3.0*Gamma) / (Gamma - 1.0) + (3.0 - nDim);
  double l = 1 / (2.0*node->GetPressure()/node->GetDensity());

  moments->A.resize(nDim);
  moments->P.resize(nDim);
  moments->N.resize(nDim);

  for(unsigned short i = 0; i<nDim ; i++){
    moments->A[i].resize(9,1.0);
    moments->P[i].resize(9,1.0);
    moments->N[i].resize(9,1.0);
  }

  moments->xi.resize(7,1.0);

  for(unsigned short i = 0; i<nDim; i++){
    su2double U = node->GetVelocity(i);

    moments->A[i][0] = 1;
    moments->P[i][0] = 0.5*erfc(-sqrt(l)*U);
    moments->N[i][0] = 0.5*erfc(sqrt(l)*U);

    moments->A[i][1] = U;
    moments->P[i][1] = U*moments->P[i][0] + 0.5*exp(-l*U*U)/sqrt(M_PI*l);
    moments->N[i][1] = U*moments->N[i][0] - 0.5*exp(-l*U*U)/sqrt(M_PI*l);

    for(unsigned int n =0; n<7; n++){
      moments->A[i][n+2] = U*moments->A[i][n+1] + (n+1)/(2*l)*moments->A[i][n];
      moments->P[i][n+2] = U*moments->P[i][n+1] + (n+1)/(2*l)*moments->P[i][n];
      moments->N[i][n+2] = U*moments->N[i][n+1] + (n+1)/(2*l)*moments->N[i][n];
    }
  }

  FluidModel->SetTDState_PT(node->GetPressure(), node->GetTemperature());
  vector<su2double> out = FluidModel->GetMaxwellMoment(nDim);
  for(unsigned short n=2; n<7; n+=2){
  	moments->xi[n] = out[n];
  }

  //Add vibrational contribution
  if(config->GetKind_FluidModel()==HARMONIC_VIBR){
		su2double* theta = config->GetTheta_v();
		su2double* w = config->GetWeight_v();
		su2double T = node->GetTemperature();
		su2double tol = 1e-6;
		int nMax = 100;

		for(unsigned short s=0; s<config->GetnVibration_mode(); s++){
      su2double Z_v_inv = 1 - exp(-theta[s]/T);
      std::vector<bool> cont(3, true);

      for(int j=1; j<nMax; j++){
        bool finished = true;

        su2double xi_v2 = j*theta[s]/(l*T);
        su2double e_xi = exp(-l*xi_v2);

        for(unsigned short i=1; i<4; i++){
          if(cont[i]){
            su2double mom = w[s]*Z_v_inv*pow(xi_v2, i)*e_xi;
            if(mom/moments->xi[2*i] < tol){
              cont[i] = false;
            }
            moments->xi[2*i] += mom;
          }

          finished = finished && !cont[i];
        }

        if(finished) break;
      }
		}
  }
}

void CGasKineticSchemeBGK::SetNormal(su2double *val_normal){
  CNumerics::SetNormal(val_normal);

  rotMatrix.clear();

  Area = 0;
  for(unsigned short iDim=0; iDim<nDim; iDim++){
    Area += pow(val_normal[iDim], 2);
  }
  Area = sqrt(Area);

  std::vector<su2double> locX(nDim, 0);
  std::vector<su2double> locY(nDim, 0);
  std::vector<su2double> locZ(nDim, 0);

  for(unsigned short iDim=0; iDim<nDim; iDim++){
    locX[iDim] = val_normal[iDim]/Area;
  }

  if(nDim==2){
    locY[0] = -locX[1];
    locY[1] = locX[0];
  }else if(nDim==3){
    if(locX[1]==0.0 && locX[2]==0.0){
      locY[1] = 1.0;
    }else{
      //cross product between locX and (1,0,0)
      locY[1] = locX[2];
      locY[2] = -locX[1];
      su2double mag = sqrt(pow(locY[1],2) + pow(locY[2],2));
      locY[1] /=mag;
      locY[2] /=mag;
    }
    //locZ is cross prod between locX and locY
    locZ[0] = locX[1]*locY[2] - locX[2]*locY[1];
    locZ[1] = locX[2]*locY[0] - locX[0]*locY[2];
    locZ[2] = locX[0]*locY[1] - locX[1]*locY[0];
  }else{
    std::logic_error("Error: number of dimensions can be only 2 or 3.");
  }

  rotMatrix.push_back(locX);
  rotMatrix.push_back(locY);
  if(nDim==3){
    rotMatrix.push_back(locZ);
  }
}

void CGasKineticSchemeBGK::rotate(su2double* v, bool inverse)const{
  su2double vRot[nDim];
  for(unsigned short i=0; i<nDim; i++){
    vRot[i] = 0;
    for(unsigned short j=0; j<nDim; j++){
      if(inverse){
        vRot[i] += v[j]*rotMatrix[j][i]; // vRot = transp(R)*v
      }else{
        vRot[i] += v[j]*rotMatrix[i][j]; // vRot = R*v
      }
    }
  }

  for(unsigned short i=0; i<nDim; i++){
    v[i] = vRot[i];
  }
}

void CGasKineticSchemeBGK::rotate(su2double** t)const{
  su2double vRot[nDim][nDim]; //rotated and transposed vector
  for(unsigned short i=0; i<nDim; i++){
    for(unsigned short j=0; j<nDim; j++){
      vRot[i][j] = t[j][i];
    }
    rotate(vRot[i]);
  }

  for(unsigned short i=0; i<nDim; i++){
    for(unsigned short j=0; j<nDim; j++){
      t[i][j] = vRot[j][i];
    }
  }
}

void CGasKineticSchemeBGK::rotate(CVariable* node)const{
  su2double* v = node->GetSolution();
  rotate(++v);

  v = node->GetPrimitive();
  rotate(++v);

  su2double** t = node->GetGradient();
  for(unsigned short iVar=0; iVar<nVar; iVar++){
    rotate(t[iVar]);
  }
  rotate(++t);
}

void CGasKineticSchemeBGK::Clear(){
  if(node_iLoc) delete node_iLoc;
  if(node_jLoc) delete node_jLoc;
  node_iLoc = NULL;
  node_jLoc = NULL;
  if(node_iRot) delete node_iRot;
  if(node_jRot) delete node_jRot;
  node_iRot = NULL;
  node_jRot = NULL;

  moments_struct* mom[3] = {&moments_i, &moments_j, &moments_I};
  for(unsigned short i=0; i<3; i++){
    mom[i]->A.clear();
    mom[i]->P.clear();
    mom[i]->N.clear();
    mom[i]->xi.clear();
  }
}

std::vector<su2double> operator+(const std::vector<su2double>& a, const std::vector<su2double>& b){
  std::vector<su2double> out(a);
  out += b;
  return out;
}

std::vector<su2double> operator+=(std::vector<su2double>& a, const std::vector<su2double>& b){
  if(a.size() != b.size()) throw std::logic_error("Error: members of operation must be of the same size.");

  for(std::size_t i=0; i<a.size(); i++){
    a[i] += b[i];
  }
  return a;
}

std::vector<su2double> operator*=(std::vector<su2double>& a, const su2double& b){
  for(std::size_t i=0; i<a.size(); i++){
    a[i] *= b;
  }
  return a;
}

std::vector<su2double> operator/=(std::vector<su2double>& a, const su2double& b){
  su2double b_inv = 1/b;
  a *= b_inv;
  return a;
}

std::vector<su2double> operator*(const std::vector<su2double>& a, const std::vector<std::vector<su2double> >& b){
  if(a.size() != b.size()) throw std::logic_error("Error: members of operation must be of compatible sizes.");

  std::vector<su2double> out(b[0].size(), 0);

  for(std::size_t i=0; i<a.size(); i++){
    out += b[i]*a[i];
  }
  return out;
}

std::vector<su2double> operator*(const std::vector<su2double>& a, const su2double& b){
  std::vector<su2double> out(a);
  out *= b;
  return out;
}

std::vector<su2double> operator*(const su2double& a, const std::vector<su2double>& b){
  std::vector<su2double> out(b);
  out *= a;
  return out;
}

std::vector<su2double> operator-=(std::vector<su2double>& a, const std::vector<su2double>& b){
  if(a.size() != b.size()) throw std::logic_error("Error: members of operation must be of the same size.");

  for(std::size_t i=0; i<a.size(); i++){
    a[i] -= b[i];
  }
  return a;
}
