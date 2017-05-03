#pragma once

#include "numerics_structure.hpp"

/*!
 * \class CGasKineticSchemeBGK
 * \brief Class for computing the Kinetic Fluxes using BGK model.
 * \author A. Mogavero J. Herrera
 * \version 5.0.0 "Raven"
 */
class CGasKineticSchemeBGK : public CNumerics {
private:
  unsigned short iDim, iVar, jVar; /*!< \brief Iteration on dimension and variables. */
  su2double *Diff_U, /*!< \brief Difference of conservative variables. */
  *Velocity_i, *Velocity_j, /*!< \brief Velocity at node 0 and 1. */
  *MeanVelocity, ProjVelocity, ProjVelocity_i, ProjVelocity_j,  /*!< \brief Mean and projected velocities. */
  *ProjFlux,  /*!< \brief Projected inviscid flux tensor. */
  Density_i, Density_j, Energy_i, Energy_j,  /*!< \brief Mean Density and energies. */
  sq_vel_i, sq_vel_j,   /*!< \brief Modulus of the velocity and the normal vector. */
  MeanDensity, MeanPressure, MeanEnthalpy, MeanEnergy, /*!< \brief Mean values of primitive variables. */
  Param_p, Param_Kappa_0, /*!< \brief Artificial dissipation parameters. */
  Local_Lambda_i, Local_Lambda_j, MeanLambda, /*!< \brief Local eingenvalues. */
  Phi_i, Phi_j, sc0, StretchingFactor, /*!< \brief Streching parameters. */
  Epsilon_0, cte; /*!< \brief Artificial dissipation values. */
  bool implicit, /*!< \brief Implicit calculation. */
  grid_movement; /*!< \brief Modification for grid movement. */
  su2double ProjGridVel;

public:

  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimension of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CGasKineticSchemeBGK(unsigned short val_nDim, unsigned short val_nVar, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CGasKineticSchemeBGK(void);

  /*!
   * \brief Compute the flow residual using a Lax method.
   * \param[out] val_resconv - Pointer to the convective residual.
   * \param[out] val_resvisc - Pointer to the artificial viscosity residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j,
                       CConfig *config);
};
