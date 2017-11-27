#pragma once

#include "fluid_model.hpp"
#include "mutation++.h"

/*!
 * \class CMutationpp
 * \brief Wrapper gas to use Mutation++ gas model.
 * \author: A. Mogavero
 * \version 5.0.0 "Raven"
 */
class CMutationpp : public CFluidModel {

  su2double Gas_Constant; //!<Ideal gas constant of the mixture

  su2double CalcGasConstant()const; //!Calculates the ideal gas constant of the mixture.

  Mutation::MixtureOptions opt; //!<option object used to create mix
  Mutation::Mixture mix; //!<Mixture object

  vector<double> comp; //!<Mixture composition in molar fraction

  /*!
   * \brief Calculate the partial density of every specie
   * @param rho total density of mixture
   * @return a density for each specie
   */
  vector<double> SpecieDensity(su2double rho)const;

  /*!
   * \brief Set the state so that it is equal to the underlying mutation++ state
   */
  void UpdateState();

  /*!
   * \brief Set the state to all non-physical values.
   */
  void SetWrongState();

public:

     /*!
     * \brief Constructor of the class.
     */
    CMutationpp(void);

    /*!
     * \brief Constructor of the class.
     */
    CMutationpp(string R, vector<double>);


    /*!
     * \brief Destructor of the class.
     */
    virtual ~CMutationpp(void);

    /*!
     * \brief Set the Dimensionless State using Density and Internal Energy
     * \param[in] rho - first thermodynamic variable.
     * \param[in] e - second thermodynamic variable.
     */

    void SetTDState_rhoe (su2double rho, su2double e );

    /*!
     * \brief Set the Dimensionless State using Pressure  and Temperature
     * \param[in] P - first thermodynamic variable.
     * \param[in] T - second thermodynamic variable.
     */

    void SetTDState_PT (su2double P, su2double T );

    /*!
     * \brief Set the Dimensionless State using Pressure and Density
     * \param[in] P - first thermodynamic variable.
     * \param[in] rho - second thermodynamic variable.
     */

    void SetTDState_Prho (su2double P, su2double rho );

    /*!
     * \brief Set the Dimensionless Internal Energy using Pressure and Density
     * \param[in] P - first thermodynamic variable.
     * \param[in] rho - second thermodynamic variable.
     */

    void SetEnergy_Prho (su2double P, su2double rho );

    /*!
     * \brief Set the Dimensionless State using Enthalpy and Entropy
     * \param[in] th1 - first thermodynamic variable (h).
     * \param[in] th2 - second thermodynamic variable (s).
     *
     */
    void SetTDState_hs (su2double h, su2double s );


    /*!
     * \brief Set the Dimensionless State using Density and Temperature
     * \param[in] th1 - first thermodynamic variable (rho).
     * \param[in] th2 - second thermodynamic variable (T).
     *
     */
    void SetTDState_rhoT (su2double rho, su2double T );

    /*!
     * \brief Set the Dimensionless State using Pressure and Entropy
     * \param[in] th1 - first thermodynamic variable (P).
     * \param[in] th2 - second thermodynamic variable (s).
     */

    void SetTDState_Ps (su2double P, su2double s );

    /*!
     * \brief Get the moments of the Maxwellian distribution with respect to the energy invariant
     * \details It only consider the internal degrees of freedom, so the macroscopic translations are not accounted.
     * The parameter nDim is used to determine how many translational degrees of freedom have to be considered internal.
     * It calculates moments of order 2, 4 and 6.
     * @param nDim number of dimensions of the CFD analysis
     * @return The wanted moments for each order
     */
    vector<su2double> GetMaxwellMoment(unsigned short nDim)const;
};
