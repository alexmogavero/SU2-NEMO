
#include "../include/kinetic_solver.hpp"
#include "../include/CKineticVariable.hpp"

CKineticSolver::CKineticSolver(void):
  CNSSolver(){
}

CKineticSolver::CKineticSolver(CGeometry *geometry, CConfig *config, unsigned short iMesh):
    CNSSolver(geometry, config, iMesh){
  for (unsigned long iPoint = 0; iPoint < nPoint; iPoint++){
    delete node[iPoint];
    node[iPoint] = new CKineticVariable(Density_Inf, Velocity_Inf, Energy_Inf, nDim, nVar, config);
  }

  /*--- Check that the initial solution is physical, report any non-physical nodes ---*/

  int counter_local = 0;

  for (unsigned long iPoint = 0; iPoint < nPoint; iPoint++) {

    su2double Density = node[iPoint]->GetSolution(0);

    su2double Velocity2 = 0.0;
    for (unsigned short iDim = 0; iDim < nDim; iDim++)
      Velocity2 += (node[iPoint]->GetSolution(iDim+1)/Density)*(node[iPoint]->GetSolution(iDim+1)/Density);

    su2double StaticEnergy= node[iPoint]->GetSolution(nDim+1)/Density - 0.5*Velocity2;

    FluidModel->SetTDState_rhoe(Density, StaticEnergy);
    su2double Pressure= FluidModel->GetPressure();
    su2double Temperature= FluidModel->GetTemperature();

    /*--- Use the values at the infinity ---*/

    if ((Pressure < 0.0) || (Density < 0.0) || (Temperature < 0.0)) {
      Solution[0] = Density_Inf;
      for (unsigned short iDim = 0; iDim < nDim; iDim++)
        Solution[iDim+1] = Velocity_Inf[iDim]*Density_Inf;
      Solution[nDim+1] = Energy_Inf*Density_Inf;
      node[iPoint]->SetSolution(Solution);
      node[iPoint]->SetSolution_Old(Solution);
      counter_local++;
    }

  }
}

CKineticSolver::~CKineticSolver(void){
}

unsigned long CKineticSolver::SetPrimitive_Variables(
    CSolver **solver_container, CConfig *config, bool Output){
  CNSSolver::SetPrimitive_Variables(solver_container, config, Output);

  for (unsigned long iPoint = 0; iPoint < nPoint; iPoint ++) {
    static_cast<CKineticVariable*>(node[iPoint])->CalculateKnudsen();
  }
}
