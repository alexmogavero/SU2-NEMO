#pragma once

#include "CGeneralIdealGas.hpp"

/*!
 * \class CPerfectGas
 * \brief Class that model a perfect gas (i.e. constant Cp).
 * \author: A. Mogavero
 * \version 5.0.0 "Raven"
 */
class CPerfectGas : public CGeneralIdealGas {

protected:
  su2double n_dof;        /*!< \brief Number of degrees of freedom of the molecule. */
  su2double Gamma;        /*!< \brief Specifc heat ratio*/

  virtual su2double Energy(su2double T)const;
  virtual su2double EnergyInv(su2double e)const;
  virtual su2double EnthalpyInv(su2double h)const;
  virtual su2double SpecificHeatVol(su2double T)const;
  virtual su2double EntropyTemp(su2double T)const;
  virtual su2double EntropyTempInv(su2double s)const;

public:

     /*!
     * \brief Constructor of the class.
     */
    CPerfectGas(void);

    /*!
     * \brief Constructor of the class.
     */
    CPerfectGas(su2double R, su2double g);


    /*!
     * \brief Destructor of the class.
     */
    virtual ~CPerfectGas(void);


};
