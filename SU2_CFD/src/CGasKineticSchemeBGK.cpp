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

  //Rotate Reference Frame
  node_iRot = node_i->duplicate();
  rotate(node_iRot);

  node_jRot = node_j->duplicate();
  rotate(node_jRot);

  //Reconstruct
  if(config->GetSpatialOrder_Flow() == SECOND_ORDER){
    su2double dist_ij = 0.0;
    for (unsigned int i=0; i<nDim; i++)
      dist_ij += (Coord_j[i]-Coord_i[i])*(Coord_j[i]-Coord_i[i]);
    dist_ij = sqrt(dist_ij);

    node_iLoc = reconstruct(node_iRot, dist_ij*0.5);
    node_jLoc = reconstruct(node_jRot, -dist_ij*0.5);
  }else{
    node_iLoc = node_iRot->duplicate();
    node_jLoc = node_jRot->duplicate();
  }

  CalculateInterface();

  su2double Dt = 0.5*(node_i->GetDelta_Time() + node_j->GetDelta_Time());

  //Calculate the mean collision time
  su2double tauColl = node_I->GetLaminarViscosity()/node_I->GetPressure();
  double K = (5.0 - 3.0*Gamma) / (Gamma - 1.0) + (3.0 - nDim);
  double lL = (K+nDim) / (4.0*(node_iLoc->GetEnergy() - 0.5*node_iLoc->GetVelocity2()));
  double lR = (K+nDim) / (4.0*(node_jLoc->GetEnergy() - 0.5*node_jLoc->GetVelocity2()));
  double rho_lamL = node_iLoc->GetDensity()/lL;
  double rho_lamR = node_jLoc->GetDensity()/lR;
  tauColl += Dt*abs(rho_lamL - rho_lamR)/abs(rho_lamL + rho_lamR);

  std::vector<su2double> Flux_i, Flux_j, Flux_I;
  Flux_I = PsiMaxwell(INTERFACE, ALL, true);
  Flux_i = PsiMaxwell(LEFT, POSITIVE, true);
  Flux_j = PsiMaxwell(RIGHT, NEGATIVE, true);

  //calculate time integrals
  su2double int_ij = tauColl - tauColl*exp(-Dt/tauColl); //integral of exp(-Dt/tauColl)
  su2double int_I = Dt - int_ij;

  su2double Dt_inv = 1/Dt;
  for(unsigned short iVar=0; iVar<nVar; iVar++){
    val_residual[iVar] = Dt_inv*(int_I*Flux_I[iVar] + int_ij*(Flux_i[iVar] + Flux_j[iVar]))*Area;
  }

  if(config->GetSpatialOrder_Flow() == SECOND_ORDER){ // if second order
    std::vector<su2double> Flux_i, Flux_j, Flux_I, Flux_I_t;
    Flux_i = std::vector<su2double>(nVar, 0);
    Flux_j = std::vector<su2double>(nVar, 0);
    Flux_I = std::vector<su2double>(nVar, 0);

    // compute space derivatives at i and j
    std::vector<std::vector<su2double> > a_i(nDim, std::vector<su2double>(nVar, 0));
    std::vector<std::vector<su2double> > a_j(nDim, std::vector<su2double>(nVar, 0));
    std::vector<su2double> A_i(nVar, 0);
    std::vector<su2double> A_j(nVar, 0);
    Derivatives(LEFT, a_i, A_i); // Dont need time derivatives
    Derivatives(RIGHT, a_j, A_j); // Dont need time derivatives

    // compute space (and time) derivatives at the interface
    std::vector<su2double> ad_i(nVar, 0);
    std::vector<su2double> ad_j(nVar, 0);
    std::vector<std::vector<su2double> > ad(nDim-1, std::vector<su2double>(nVar, 0));
    std::vector<su2double> Ad(nVar, 0);

    std::vector<su2double> g(5, 0);
    g[0] = Dt - tauColl*(1-exp(-Dt/tauColl));
    g[1] = -(1-exp(-Dt/tauColl))/g[0];
    g[2] = (-Dt + 2*tauColl*(1-exp(-Dt/tauColl)) - Dt*exp(-Dt/tauColl))/g[0];
    g[3] = (1-exp(-Dt/tauColl))/g[0];
    g[4] = (Dt*exp(-Dt/tauColl) - tauColl*(1-exp(-Dt/tauColl)))/g[0];

    Interface_Derivatives(ad_i, ad_j, ad, Ad, a_i, a_j, g);

    // Assemble fluxes
    std::vector<unsigned short> exponents;
    exponents.assign(nVar-1, 0);
    exponents[0] += 2;
    Flux_I += ad_i*PsiPsiMaxwell(INTERFACE, POSITIVE, exponents);
    Flux_I += ad_j*PsiPsiMaxwell(INTERFACE, NEGATIVE, exponents);
    for(unsigned short i=1; i<nDim; i++){
      exponents.assign(nVar-1, 0);
      exponents[0]++;
      exponents[i]++;
      Flux_I += ad[i-1]*PsiPsiMaxwell(INTERFACE, ALL, exponents);
    }

    exponents.assign(nVar-1, 0);
    exponents[0]++;
    Flux_I_t = Ad*PsiPsiMaxwell(INTERFACE, ALL, exponents);

    for(unsigned short i=0; i<nDim; i++){
      exponents.assign(nVar-1, 0);
      exponents[0]++;
      exponents[i]++;
      Flux_i -= a_i[i]*PsiPsiMaxwell(LEFT, POSITIVE, exponents);
      Flux_j -= a_j[i]*PsiPsiMaxwell(RIGHT, NEGATIVE, exponents);
    }

    su2double int_t_exp = tauColl*(tauColl - (Dt + tauColl)*exp(-Dt/tauColl)); //integral of t*exp(-t/tauColl)

    for(unsigned short iVar=0; iVar<nVar; iVar++){
      val_residual[iVar] += Dt_inv * tauColl * (2*tauColl - Dt - exp(-Dt/tauColl) * (Dt+2*tauColl)) * Flux_I[iVar] * Area;
      val_residual[iVar] += Dt_inv * tauColl * (pow(Dt,2)/(2*tauColl) - Dt + int_ij) * Flux_I_t[iVar] * Area;
      val_residual[iVar] += Dt_inv * int_t_exp *(Flux_i[iVar] + Flux_j[iVar]) * Area;
    }
  }

  if(config->GetViscous()){
    std::vector<std::vector<su2double> > a_i(nDim, std::vector<su2double>(nVar, 0)); //space derivatives
    std::vector<std::vector<su2double> > a_j(nDim, std::vector<su2double>(nVar, 0)); //space derivatives
    std::vector<su2double> A_i(nVar, 0); //Time derivatives
    std::vector<su2double> A_j(nVar, 0); //Time derivatives
    Derivatives(LEFT, a_i, A_i);
    Derivatives(RIGHT, a_j, A_j);

    Flux_i = std::vector<su2double>(nVar, 0);
    Flux_j = std::vector<su2double>(nVar, 0);
    std::vector<unsigned short> exponents;
    for(unsigned short i=0; i<nDim; i++){
      exponents.assign(nVar-1, 0);
      exponents[0]++;
      exponents[i]++;
      Flux_i += a_i[i]*PsiPsiMaxwell(LEFT, POSITIVE, exponents);
      Flux_j += a_j[i]*PsiPsiMaxwell(RIGHT, NEGATIVE, exponents);
    }
    exponents.assign(nVar-1, 0);
    exponents[0]++;
    Flux_i += A_i*PsiPsiMaxwell(LEFT, POSITIVE, exponents);
    Flux_j += A_j*PsiPsiMaxwell(RIGHT, NEGATIVE, exponents);

    for(unsigned short iVar=0; iVar<nVar; iVar++){
      val_residual[iVar] -= Dt_inv*tauColl*int_ij*(Flux_i[iVar] + Flux_j[iVar])*Area;
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

CVariable* CGasKineticSchemeBGK::reconstruct(const CVariable* node, const su2double& d)const{
  CVariable* out = node->duplicate();

  su2double* v = out->GetSolution();
  for(unsigned int iVar=0; iVar<nVar; iVar++){
    v[iVar] += out->GetGradient(iVar, 0)*d;
  }

  out->SetNon_Physical(false);
  bool RightSol = out->SetPrimVar(FluidModel);
  if (!RightSol) {
    delete out;
    out = node->duplicate();
  }

  return out;
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
