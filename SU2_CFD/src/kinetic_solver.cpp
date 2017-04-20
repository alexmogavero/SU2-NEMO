
#include "../include/kinetic_solver.hpp"

CKineticSolver::CKineticSolver(void):
  CNSSolver(){
}

CKineticSolver::CKineticSolver(CGeometry *geometry, CConfig *config, unsigned short iMesh):
    CNSSolver(geometry, config, iMesh){
}

CKineticSolver::~CKineticSolver(void){
}
