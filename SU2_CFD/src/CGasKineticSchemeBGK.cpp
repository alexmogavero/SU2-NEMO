#define _USE_MATH_DEFINES
#include <cmath>
#include "../include/CGasKineticSchemeBGK.hpp"
#include "../include/CKineticVariable.hpp"

CGasKineticSchemeBGK::CGasKineticSchemeBGK(unsigned short val_nDim, unsigned short val_nVar, CConfig *config):
  CNumerics(val_nDim, val_nVar, config),
  node_I(NULL),
  config(config){
}

CGasKineticSchemeBGK::~CGasKineticSchemeBGK(void) {
  if(node_I){
    delete node_I;
  }
}

void CGasKineticSchemeBGK::ComputeResidual(su2double *val_residual, CConfig *config){
  //TODO implement rotations in order to take into account the cases when the normal is not aligned with u
  CalculateInterface();

  //Calculate the mean collision time
  //TODO check if it is ok to calculate it on the interface
  su2double tauColl = node_I->GetLaminarViscosity()/node_I->GetPressure();

  std::vector<su2double> Flux_i, Flux_j, Flux_I;
  Flux_I = PsiMaxwell(INTERFACE, ALL, true);
  Flux_i = PsiMaxwell(LEFT, POSITIVE, true);
  Flux_j = PsiMaxwell(RIGHT, NEGATIVE, true);

  su2double Dt = 0.5*(node_i->GetDelta_Time() + node_j->GetDelta_Time());
  //calculate time integrals
  su2double int_ij = tauColl - tauColl*exp(-Dt/tauColl); //integral of exp(-Dt/tauColl)
  su2double int_I = Dt - int_ij;

  su2double Dt_inv = 1/Dt;
  for(unsigned short iVar; iVar<nVar; iVar++){
    val_residual[iVar] = Dt_inv*(int_I*Flux_I[iVar] + int_ij*(Flux_i[iVar] + Flux_j[iVar]))*Area;
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

  node_I = new CKineticVariable(U_I.data(), nDim, nVar, config);

  node_I->SetNon_Physical(false);

  bool RightSol = node_I->SetPrimVar(FluidModel);
  //node_I->SetSecondaryVar(FluidModel);

  if (!RightSol) {
    node_I->SetNon_Physical(true);
  }
}

std::vector<su2double> CGasKineticSchemeBGK::PsiMaxwell(State state, IntLimits lim, bool uPsi){
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

su2double CGasKineticSchemeBGK::MomentsMaxwellian(std::vector<unsigned short> exponents, State state, IntLimits lim){
  su2double mp, rho;
  
  switch (state){
    case LEFT:
      if (moments_i.xi.empty()) CGasKineticSchemeBGK::ComputeMaxwellianMoments(node_i, &moments_i);
      switch (lim) {
        case ALL:
          mp = 1.0;
          for (unsigned short i = 0; i<nDim; i++){
            mp *= moments_i.A[i][exponents[i]];
          }
          mp *= moments_i.xi[exponents[nDim]];
          break;
        case NEGATIVE:
          mp = 1.0;
          for (unsigned short i = 0; i<nDim; i++){
            mp *= moments_i.N[i][exponents[i]];
          }
          mp *= moments_i.xi[exponents[nDim]];
          break;
        case POSITIVE:
          mp = 1.0;
          for (unsigned short i = 0; i<nDim; i++){
            mp *= moments_i.P[i][exponents[i]];
          }
          mp *= moments_i.xi[exponents[nDim]];
          break;
      }
      rho = node_i->GetDensity();
    break;
    case RIGHT:
      if (moments_j.xi.empty()) CGasKineticSchemeBGK::ComputeMaxwellianMoments(node_j, &moments_j);
      switch (lim) {
        case ALL:
          mp = 1.0;
          for (unsigned short i = 0; i<nDim; i++){
            mp *= moments_j.A[i][exponents[i]];
          }
          mp *= moments_j.xi[exponents[nDim]];
          break;
        case NEGATIVE:
          mp = 1.0;
          for (unsigned short i = 0; i<nDim; i++){
            mp *= moments_j.N[i][exponents[i]];
          }
          mp *= moments_j.xi[exponents[nDim]];
          break;
        case POSITIVE:
          mp = 1.0;
          for (unsigned short i = 0; i<nDim; i++){
            mp *= moments_j.P[i][exponents[i]];
          }
          mp *= moments_j.xi[exponents[nDim]];
          break;
      }
      rho = node_j->GetDensity();
    break;
    case INTERFACE:
      if (moments_I.xi.empty()) CGasKineticSchemeBGK::ComputeMaxwellianMoments(node_I, &moments_I);
      switch (lim) {
        case ALL:
          mp = 1.0;
          for (unsigned short i = 0; i<nDim; i++){
            mp *= moments_I.A[i][exponents[i]];
          }
          mp *= moments_I.xi[exponents[nDim]];
          break;
        case NEGATIVE:
          mp = 1.0;
          for (unsigned short i = 0; i<nDim; i++){
            mp *= moments_I.N[i][exponents[i]];
          }
          mp *= moments_I.xi[exponents[nDim]];
          break;
        case POSITIVE:
          mp = 1.0;
          for (unsigned short i = 0; i<nDim; i++){
            mp *= moments_I.P[i][exponents[i]];
          }
          mp *= moments_I.xi[exponents[nDim]];
          break;
      }
      rho = node_I->GetDensity();
    break;
  }
  
  return mp*rho;
}

void CGasKineticSchemeBGK::ComputeMaxwellianMoments(CVariable* node, moments_struct*  moments){
  double K = (5.0 - 3.0*Gamma) / (Gamma - 1.0) + (3.0 - (nDim - 2.0));
  double l = (K+nDim) / (4.0*(node->GetEnergy() - 0.5*node->GetVelocity2()));
  
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
  
  moments->xi[2] = 0.5 * K / l;
  moments->xi[4] = 0.5 * moments->xi[2] * (K+2)/ l;
  moments->xi[6] = 0.5 * moments->xi[4] * (K+4)/ l;
}
