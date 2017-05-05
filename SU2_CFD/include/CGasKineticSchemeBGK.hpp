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

protected:
  /*!
   * \brief define the integration limits
   */
  enum IntLimits{
    ALL,     //!< from \f$-\infty\f$ to \f$+\infty\f$
    NEGATIVE,//!< from \f$-\infty\f$ to 0
    POSITIVE //!< from 0 to \f$+\infty\f$
  };

  /*!
   * \brief identify the state
   */
  enum State{
    LEFT,     //!< left side of the edge i.e. i
    RIGHT,    //!< left side of the edge i.e. j
    INTERFACE //!< left side of the edge i.e. I
  };
  su2double* U_I; //!< vector of conserved quantities at the interface

  CFluidModel* FluidModel; //!< Thermodynamic model of the fluid

  //TODO define those in the base class
  CVariable* node_i; //!< Node that stores all the variables on the left size of the edge
  CVariable* node_j; //!< Node that stores all the variables on the right size of the edge

  CKineticVariable* node_I; //!< Node that stores all the variables at the interface

  /*!
   * \brief calculates the moments of the Maxwellian distribution
   * \details a generic moment of the Maxwellian is defined by means of the following:
   *  \f[
   *    \rho \langle\varphi\rangle = \iint_{-\infty}^{+\infty} \varphi f_0 d\mathbf{u}\mathbf{\xi}
   *  \f]
   *  where the moment can be always be:
   *  \f[
   *    \rho \langle\varphi\rangle = \rho\langle u^pv^qw^r\mathbf{\xi}^s\rangle =
   *    \rho\langle u^p\rangle \langle v^q\rangle\langle w^r\rangle \langle\mathbf{\xi}^s\rangle
   *  \f]
   *  the integration limits can be between \f$-\infty\f$ and \f$+\infty\f$, from 0
   *  to \f$+\infty\f$ or from \f$-\infty\f$ to 0
   * @param exponents - exponents that define the phi function (p,q,r,s)
   * @param theta - thermodynamic parameter dependent on temperature that define the Maxwellian
   * @param lim - flag that defines the integration limits
   * @return
   */
  su2double MomentsMaxwellian(std::vector<unsigned short> exponents, State state, IntLimits lim)const;

  /*!
   * \brief Calculates the conserved quantities at the interface.
   */
  void CalculateInterface()const;

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

  /*!
   * \brief Compute the flow residual using the GKS BGK method.
   * \param[out] val_residual - Pointer to the convective residual.
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, CConfig *config);
};
