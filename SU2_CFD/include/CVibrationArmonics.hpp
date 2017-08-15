#pragma once

#include "CPerfectGas.hpp"

/*!
 * \class CVibrationArmonics
 * \brief Class that model a Maxwell-Boltzmann gas that uses the infinite vibrational armonics model.
 * \author: A. Mogavero
 * \version 5.0.0 "Raven"
 */
class CVibrationArmonics : public CPerfectGas {
	std::vector<su2double> ThetaVib;        /*!< \brief Gas Constant. */
	std::vector<su2double> weight;       /*!< \brief Weight to apply to each vibrational contribution*/
protected:

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
    CVibrationArmonics(void);

    /*!
     * \brief Constructor of the class.
     */
    CVibrationArmonics(su2double R, su2double g,
    		unsigned short n_modes, su2double* theta, su2double* w);

    CVibrationArmonics(su2double R, su2double g,
    		std::vector<su2double> theta, std::vector<su2double> w);


    /*!
     * \brief Destructor of the class.
     */
    virtual ~CVibrationArmonics(void);


};
