#define _USE_MATH_DEFINES
#include <cmath>
#include "../include/CGasKineticSchemeBGK.hpp"
#include "../include/CKineticVariable.hpp"

CGasKineticSchemeBGK::CGasKineticSchemeBGK(unsigned short val_nDim, unsigned short val_nVar, CConfig *config):
  CNumerics(val_nDim, val_nVar, config),
  FluidModel(NULL),
  node_i(NULL),
  node_j(NULL),
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
  su2double mp;
  
  switch (state){
    case LEFT:
      if (moments_i.xi.empty()) CGasKineticSchemeBGK::ComputeMaxwellianMoments(node_i, &moments_i);
      switch (lim) {
        case ALL:
          mp = moments_i.A[0][exponents[0]] * moments_i.A[1][exponents[1]] * moments_i.A[2][exponents[2]] *moments_i.xi[exponents[3]];
        case NEGATIVE:
          mp = moments_i.N[0][exponents[0]] * moments_i.N[1][exponents[1]] * moments_i.N[2][exponents[2]] *moments_i.xi[exponents[3]];
        case POSITIVE:
          mp = moments_i.P[0][exponents[0]] * moments_i.P[1][exponents[1]] * moments_i.P[2][exponents[2]] *moments_i.xi[exponents[3]];
      }
    case RIGHT:
      if (moments_j.xi.empty()) CGasKineticSchemeBGK::ComputeMaxwellianMoments(node_j, &moments_j);
      switch (lim) {
        case ALL:
          mp = moments_j.A[0][exponents[0]] * moments_j.A[1][exponents[1]] * moments_j.A[2][exponents[2]] *moments_j.xi[exponents[3]];
        case NEGATIVE:
          mp = moments_j.N[0][exponents[0]] * moments_j.N[1][exponents[1]] * moments_j.N[2][exponents[2]] *moments_j.xi[exponents[3]];
        case POSITIVE:
          mp = moments_j.P[0][exponents[0]] * moments_j.P[1][exponents[1]] * moments_j.P[2][exponents[2]] *moments_j.xi[exponents[3]];
      }
    case INTERFACE:
      if (moments_I.xi.empty()) CGasKineticSchemeBGK::ComputeMaxwellianMoments(node_I, &moments_I);
      switch (lim) {
        case ALL:
          mp = moments_I.A[0][exponents[0]] * moments_I.A[1][exponents[1]] * moments_I.A[2][exponents[2]] *moments_I.xi[exponents[3]];
        case NEGATIVE:
          mp = moments_I.N[0][exponents[0]] * moments_I.N[1][exponents[1]] * moments_I.N[2][exponents[2]] *moments_I.xi[exponents[3]];
        case POSITIVE:
          mp = moments_I.P[0][exponents[0]] * moments_I.P[1][exponents[1]] * moments_I.P[2][exponents[2]] *moments_I.xi[exponents[3]];
      }
  }
  
  return mp;
}

void CGasKineticSchemeBGK::ComputeMaxwellianMoments(CVariable* node, moments_struct*  moments){
  double K = (5.0 - 3.0*Gamma) / (Gamma - 1.0) + (3.0 - (nDim - 2.0));
  double l = (K+nDim) * node->GetDensity() / (4.0*((node->GetEnergy()-0.5*node->GetVelocity2()) - 0.5*node->GetDensity()*node->GetVelocity2()));
  
  
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
    
    moments->A[i][0] = U;
    moments->P[i][0] = U*moments->P[i][0] + 0.5*exp(-l*U*U)/sqrt(M_PI*l);
    moments->N[i][0] = U*moments->N[i][0] - 0.5*exp(-l*U*U)/sqrt(M_PI*l);
    
    for(unsigned int n =0; n<9; n++){
      moments->A[i][n+2] = U*moments->A[i][n+1] + (n+1)/(2*l)*moments->A[i][n];
      moments->P[i][n+2] = U*moments->P[i][n+1] + (n+1)/(2*l)*moments->P[i][n];
      moments->N[i][n+2] = U*moments->N[i][n+1] + (n+1)/(2*l)*moments->N[i][n];
    }
  }
  
  moments->xi[2] = 0.5 * K / l;
  moments->xi[4] = 0.5 * moments->xi[2] * (K+2)/ l;
  moments->xi[6] = 0.5 * moments->xi[4] * (K+4)/ l;
}
