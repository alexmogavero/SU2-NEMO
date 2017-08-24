
#include "../include/kinetic_solver.hpp"
#include "../include/CKineticVariable.hpp"
#include "../include/CGasKineticSchemeBGK.hpp"

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

void CKineticSolver::Viscous_Residual(CGeometry *geometry, CSolver **solver_container, CNumerics *numerics,
                                   CConfig *config, unsigned short iMesh, unsigned short iRKStep){
}

void CKineticSolver::BC_Kinetic_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
    CNumerics *visc_numerics, CConfig *config, unsigned short val_marker, su2double accom) {

  unsigned short iDim;
  unsigned long iVertex, iPoint;

  su2double *Normal;
  su2double Twall;

  bool implicit = (config->GetKind_TimeIntScheme_Flow() == EULER_IMPLICIT);
  bool grid_movement  = config->GetGrid_Movement();

  /*--- Identify the boundary ---*/

  string Marker_Tag = config->GetMarker_All_TagBound(val_marker);

  /*--- Retrieve the specified wall temperature ---*/

  Twall = config->GetIsothermal_Temperature(Marker_Tag)/config->GetTemperature_Ref();

  CVariable* nodeB = NULL;

  /*--- Loop over boundary points ---*/

  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {

    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    if (geometry->node[iPoint]->GetDomain()) {

      /*--- Compute dual-grid area and boundary normal ---*/

      Normal = geometry->vertex[val_marker][iVertex]->GetNormal();

      conv_numerics->SetNormal(Normal);
      conv_numerics->SetNodes(node[iPoint], NULL);

      static_cast<CGasKineticSchemeBGK*>(conv_numerics)->GetInviscidProjFlux(Res_Conv, CGasKineticSchemeBGK::POSITIVE);

      LinSysRes.AddBlock(iPoint, Res_Conv);

      // TODO put theta calculation in CKineticVar
      double K = (5.0 - 3.0*Gamma) / (Gamma - 1.0) + (3.0 - nDim);
      double l_i = (K+nDim) / (4.0*(node[iPoint]->GetEnergy() - 0.5*node[iPoint]->GetVelocity2()));

      FluidModel->SetTDState_PT(node[iPoint]->GetPressure(), Twall);
      su2double E_w = FluidModel->GetStaticEnergy();
      double l_w = (K+nDim) / (4.0*E_w);

      su2double rho_w = 2*sqrt(M_PI*l_w)*Res_Conv[0];

      su2double rho_l_ref = accom*(rho_w/l_w) + 2*(Gamma-1)*(1-accom)*node[iPoint]->GetDensity()*node[iPoint]->GetEnergy();
      su2double rho_ref = 4*M_PI*Res_Conv[0]/rho_l_ref;
      su2double l_ref = rho_ref/rho_l_ref;
      su2double E_ref = (K+nDim) / (4.0*l_ref);

      //TODO check wether copying the gradients is ok
      nodeB = node[iPoint]->duplicate();
      nodeB->SetSolution(0, rho_ref);
      for(iDim=0; iDim<nDim; iDim++){
        nodeB->SetSolution(iDim+1, 0);
      }
      nodeB->SetSolution(nVar-1, rho_ref*E_ref);
      nodeB->SetNon_Physical(false);
      bool RightSol = nodeB->SetPrimVar(FluidModel);
      if (!RightSol) nodeB->SetNon_Physical(true);

      conv_numerics->SetNormal(Normal);
      conv_numerics->SetNodes(nodeB, NULL);

      static_cast<CGasKineticSchemeBGK*>(conv_numerics)->GetInviscidProjFlux(Res_Conv, CGasKineticSchemeBGK::NEGATIVE);

      LinSysRes.AddBlock(iPoint, Res_Conv);

      /*--- Calculate Jacobian for implicit time stepping ---*/

      if (implicit) {
        throw std::logic_error("Error: Implicit not implemented in kinetic solver.");
      }

      /*--- If the wall is moving, there are additional residual contributions
       due to pressure (p v_wall.n) and shear stress (tau.v_wall.n). ---*/

      if (grid_movement) {
        throw std::logic_error("Error: grid movement not implemented in kinetic solver.");
      }
    }
  }
}

void CKineticSolver::BC_Isothermal_Wall(CGeometry *geometry, CSolver **solver_container, CNumerics *conv_numerics,
    CNumerics *visc_numerics, CConfig *config, unsigned short val_marker) {
  BC_Kinetic_Wall(geometry, solver_container, conv_numerics,
    visc_numerics, config, val_marker, 1);
}

void CKineticSolver::Preprocessing(CGeometry *geometry, CSolver **solver_container, CConfig *config,
    unsigned short iMesh, unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output) {
  CNSSolver::Preprocessing(geometry, solver_container, config, iMesh, iRKStep, RunTime_EqSystem, Output);

  if (config->GetKind_Gradient_Method() == GREEN_GAUSS) SetSolution_Gradient_GG(geometry, config);
  if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) SetSolution_Gradient_LS(geometry, config);
}
