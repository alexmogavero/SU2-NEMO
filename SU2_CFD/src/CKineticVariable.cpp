
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

CKineticVariable::CKineticVariable(const CKineticVariable& obj):
  CNSVariable(obj),
  knudsenLocal(obj.knudsenLocal){
}

CKineticVariable::~CKineticVariable(void){
}

inline CVariable* CKineticVariable::duplicate()const{
  return new CKineticVariable(*this);
}

void CKineticVariable::CalculateKnudsen(){
  knudsenLocal = CNSVariable::GetKnudsen();
}

su2double CKineticVariable::GetKnudsen()const{
	return knudsenLocal;
}

bool CKineticVariable::SetPrimVar(CFluidModel *FluidModel){
  return SetPrimVar(0, 0, FluidModel);
}
