
#include "../include/CKineticVariable.hpp"

CKineticVariable::CKineticVariable(void):
  CNSVariable(){
}

CKineticVariable::CKineticVariable(su2double val_density, su2double *val_velocity,
    su2double val_energy, unsigned short val_nDim, unsigned short val_nvar, CConfig *config):
      CNSVariable(val_density, val_velocity, val_energy, val_nDim, val_nvar, config),
      knudsenLocal(0.0){
}

CKineticVariable::CKineticVariable(su2double *val_solution, unsigned short val_nDim,
    unsigned short val_nvar, CConfig *config):
      CNSVariable(val_solution, val_nDim, val_nvar, config),
      knudsenLocal(0.0){
}

CKineticVariable::~CKineticVariable(void){
}

void CKineticVariable::CalculateKnudsen(){
  unsigned short iDens = nDim + 2;
  su2double knDens = CalcMagnitude(Gradient_Primitive[iDens])/Primitive[iDens];

  unsigned short iTemp = 0;
  su2double knTemp = CalcMagnitude(Gradient_Primitive[iTemp])/Primitive[iTemp];

  /*Calculates the gradient of the magnitude of U from the gradient of U
   * gradMagU = (1/magU)*(ux*grad(ux) + uy*grad(uy) + uz*grad(uz))
   */
  su2double gradMagU[nDim];
  for(unsigned short i=0; i<nDim; i++){
    gradMagU[i] = 0.0;
  }
  for(unsigned short i=0; i<nDim; i++){
    for(unsigned short j=0; j<nDim; j++){
      gradMagU[i] += Primitive[i+1]*Gradient_Primitive[i+1][j];
    }
  }

  su2double knU = CalcMagnitude(gradMagU)/Velocity2;

  knudsenLocal = max(knDens, max(knTemp, knU));
}

su2double CKineticVariable::CalcMagnitude(su2double* v)const{
  su2double out = 0;

  for(unsigned short i=0; i<nDim; i++){
    out += pow(v[i], 2.0);
  }
  return sqrt(out);
}
