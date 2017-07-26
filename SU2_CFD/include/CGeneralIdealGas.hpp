#pragma once

#include "fluid_model.hpp"

/*!
 * \class CGeneralIdealGas
 * \brief Abstract class for defining a generic ideal gas model.
 * \author: A. Mogavero
 * \version 5.0.0 "Raven"
 */
class CGeneralIdealGas : public CFluidModel {

protected:
  su2double Gas_Constant;        /*!< \brief Gas Constant. */

  virtual su2double Energy(su2double T)const=0;
  virtual su2double EnergyInv(su2double e)const=0;
  virtual su2double EnthalpyInv(su2double h)const=0;
  virtual su2double SpecificHeatVol(su2double T)const=0;
  virtual su2double EntropyTemp(su2double T)const=0;
  virtual su2double EntropyTempInv(su2double s)const=0;

public:

     /*!
     * \brief Constructor of the class.
     */
    CGeneralIdealGas(void);

    /*!
     * \brief Constructor of the class.
     */
    CGeneralIdealGas(su2double R);


    /*!
     * \brief Destructor of the class.
     */
    virtual ~CGeneralIdealGas(void);

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
};
