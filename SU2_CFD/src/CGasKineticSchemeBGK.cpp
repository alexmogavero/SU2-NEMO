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
}

void CGasKineticSchemeBGK::ComputeResidual(su2double *val_residual, CConfig *config){
  Clear();

  //Rotate Reference Frame
  node_iLoc = node_i->duplicate();
  rotate(node_iLoc);

  node_jLoc = node_j->duplicate();
  rotate(node_jLoc);

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
  for(unsigned short iVar=0; iVar<nVar; iVar++){
    val_residual[iVar] = Dt_inv*(int_I*Flux_I[iVar] + int_ij*(Flux_i[iVar] + Flux_j[iVar]))*Area;
  }

  if(config->GetViscous()){
    std::vector<std::vector<su2double> > vFlux_i, vFlux_j; //u*u*Psi, u*v*Psi, u*w*Psi
    for(unsigned short i=0; i<nDim; i++){
      std::vector<unsigned short> exponents(nVar-1, 0);
      exponents[0]++;
      exponents[i]++;
      vFlux_i.push_back(PsiMaxwell(LEFT, POSITIVE, exponents));
      vFlux_j.push_back(PsiMaxwell(RIGHT, NEGATIVE, exponents));
    }

    std::vector<su2double> der_i(nDim+1, 1); //TODO calculate derivatives
    std::vector<su2double> der_j(nDim+1, 1); //TODO calculate derivatives
    
    std::vector<su2double> deritest;
    deritest = Derivatives(LEFT);

    for(unsigned short iVar=0; iVar<nVar; iVar++){
      val_residual[iVar] += Dt_inv*tauColl*int_ij*(der_i[0]*Flux_i[iVar] + der_j[0]*Flux_j[iVar])*Area;

      for(unsigned short iDim=0; iDim<nDim; iDim++){
        val_residual[iVar] += Dt_inv*tauColl*int_ij*(der_i[iDim+1]*vFlux_i[iDim][iVar] +
            der_j[iDim+1]*vFlux_j[iDim][iVar])*Area;
      }
    }
  }

  rotate(val_residual + 1, true);
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

std::vector<su2double> CGasKineticSchemeBGK::PsiPsiMaxwell(State state){
  std::vector<su2double> out(nVar*nVar, 0);
  std::vector<unsigned short> exponents(nVar-1, 0);
  
  std::vector<su2double> tmprow = PsiMaxwell(state, ALL, false); //1*psi
  for(unsigned short i=0; i<nVar; i++){
    out[i] = tmprow[i];
  }
  
  tmprow = PsiMaxwell(state, ALL, true); //u*psi
  for(unsigned short i=0; i<nVar; i++){
    out[nVar+i] = tmprow[i];
  }

  //v*psi, w*psi
  for(unsigned short i=1; i<nDim; i++){
    for(unsigned short j=1; j<nDim; j++){
      exponents.assign(nVar-1, 0);
      exponents[i]++;
      exponents[j]++;
      out[i+1 + nVar*(j+1)] = MomentsMaxwellian(exponents, state, ALL);
    }

    for(unsigned short iDim=0; iDim<nDim; iDim++){
      exponents.assign(nVar-1, 0);
      exponents[iDim] = 2;
      exponents[i]++;
      out[i+1 + nVar*(nVar-1)] += MomentsMaxwellian(exponents, state, ALL);
    }
    exponents.assign(nVar-1, 0);
    exponents[nVar-2] = 2;
    exponents[i]++;
    out[i+1 + nVar*(nVar-1)] += MomentsMaxwellian(exponents, state, ALL);
    out[i+1 + nVar*(nVar-1)] /= 2;
  }

  /*xi*psi that actually is only the component xi*xi
    the formula for 3D is:
      (1/4)*(u^4 + v^4 + w^4 + 2*u^2*v^2 + 2*u^2*w^2 + 2*v^2*w^2 +
        2*u^2*xi^2 + 2*v^2*xi^2 + 2*w^2*xi^2 + xi^4) */
  for(unsigned short i=0; i<nDim; i++){
    exponents.assign(nVar-1, 0);
    exponents[i] = 4;
    out[nVar-1 + nVar*(nVar-1)] += MomentsMaxwellian(exponents, state, ALL);

    for(unsigned short j=i+1; j<nDim; j++){
      exponents.assign(nVar-1, 0);
      exponents[i] = 2;
      exponents[j] = 2;
      out[nVar-1 + nVar*(nVar-1)] += 2*MomentsMaxwellian(exponents, state, ALL);
    }

    exponents.assign(nVar-1, 0);
    exponents[i] = 2;
    exponents[nVar-2] = 2;
    out[nVar-1 + nVar*(nVar-1)] += 2*MomentsMaxwellian(exponents, state, ALL);
  }

  exponents.assign(nVar-1, 0);
  exponents[nVar-2] = 4;
  out[nVar-1 + nVar*(nVar-1)] += MomentsMaxwellian(exponents, state, ALL);
  out[nVar-1 + nVar*(nVar-1)] /= 4;

  //Build the symmetrical part of the matrix
  for(unsigned short i=2; i<nVar; i++){
    for(unsigned short j=0; j<i; j++){
      out[j + nVar*i] = out[i + nVar*j];
    }
  }

  return out;
}

std::vector<su2double> CGasKineticSchemeBGK::Derivatives(State state){
  
  std::vector<su2double> M;
  M = PsiPsiMaxwell(state);
  
  CVariable* node;
  
  std::vector<su2double> b(nVar, 0);  //nVar?
  
  switch (state){
    case LEFT:
      node = node_i;
      break;
    case RIGHT:
      node = node_j;
      break;
    case INTERFACE:
      node = node_I;
      break;
  }  
  
  for (unsigned int i=0; i++; i<nVar){
    b[i] = node->GetGradient(i, 0);
  }
  
  int ipiv[nVar];
  
  LAPACKE_dsysv(LAPACK_COL_MAJOR, 'U', nVar, 1, M.data(), nVar, ipiv, b.data(), nVar);	
  
  return b;
}

su2double CGasKineticSchemeBGK::MomentsMaxwellian(std::vector<unsigned short> exponents, State state, IntLimits lim){
  su2double mp, rho;
  
  moments_struct* moments;
  CVariable* node;
  
  switch (state){
    case LEFT:
      node = node_i;
      moments = &moments_i;
    break;
    case RIGHT:
      node = node_j;
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
