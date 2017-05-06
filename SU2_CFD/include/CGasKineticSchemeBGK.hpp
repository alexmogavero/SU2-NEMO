#pragma once

#include "numerics_structure.hpp"

class CKineticVariable;

/*!
 * \class CGasKineticSchemeBGK
 * \brief Class for computing the Kinetic Fluxes using BGK model.
 * \author A. Mogavero J. Herrera
 * \version 5.0.0 "Raven"
 */
class CGasKineticSchemeBGK : public CNumerics {
private:
  CConfig* config; //!Configuration object of the whole CFD

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
   * \brief define the structure that holds the Maxwellian moments
   */
  struct moments_struct{
    std::vector<std::vector<double>> A;
    std::vector<std::vector<double>> P;
    std::vector<std::vector<double>> N;
    std::vector<double> xi;
  };
  
  moments_struct moments_i; //!< Structure that stores the Maxwellian moments on the left size of the edge
  moments_struct moments_j; //!< Structure that stores the Maxwellian moments on the right size of the edge
  moments_struct moments_I; //!< Structure that stores the Maxwellian moments at the interface

  /*!
   * \brief identify the state
   */
  enum State{
    LEFT,     //!< left side of the edge i.e. i
    RIGHT,    //!< left side of the edge i.e. j
    INTERFACE //!< left side of the edge i.e. I
  };

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
   * @param state - flag that defines the state
   * @param lim - flag that defines the integration limits
   * @return
   */
  su2double MomentsMaxwellian(std::vector<unsigned short> exponents, State state, IntLimits lim);
  
  /*!
   * \brief Actual function that computes the Maxwellian moments for a node 
   */
  void ComputeMaxwellianMoments(CVariable* node, moments_struct*  moments);

  /*!
   * \brief Calculate the moments of the Maxwellian distribution with function \f$ \varphi=\psi\f$ .
   * \details \f$\psi\f$ is defined so to calculate the conserved quantities.
   *  \f[
   *    \mathbf{w} = \iint_{-\infty}^{+\infty} \mathbf{\psi} f d\mathbf{u}\mathbf{\xi}
   *  \f]
   *  the fluxes can be calculated using the flag uPsi.
   * @param state identify the state for wich the momemnts will be calculated
   * @param lim defines the integration limits
   * @param if true moments are calculated for \f$ \varphi=u\psi\f$ instead
   * @return a vector with a moment for each conserved quantity.
   */
  std::vector<su2double> PsiMaxwell(State state, IntLimits lim, bool uPsi=false);

  /*!
   * \brief Calculates the state at the interface.
   * \details create the node `node_I` and calculate its primitive variables.
   */
  void CalculateInterface();

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
   * \brief Compute the flow residual using the GKS BGK method.
   * \param[out] val_residual - Pointer to the convective residual.
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, CConfig *config);
};
