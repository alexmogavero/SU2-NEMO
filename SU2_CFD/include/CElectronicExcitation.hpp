#pragma once

#include "CVibrationArmonics.hpp"

/*!
 * \class CElectronicExcitation
 * \brief Class that model a Maxwell-Boltzmann gas with electronic excitation model for specific heat.
 * \author: A. Mogavero
 * \version 5.0.0 "Raven"
 */
class CElectronicExcitation : public CVibrationArmonics {
	std::vector<su2double> g;
	std::vector<su2double> ThetaEl;

	su2double Cv(su2double T)const;

	/*! \brief Class to perform boost integration. */
	class Int_Cv {

	    const CElectronicExcitation* obj;

	public:
	    Int_Cv( const CElectronicExcitation* f );

	    void operator() ( const su2double &x , su2double &dxdt , const double T );
	};

protected:

  virtual su2double Energy(su2double T)const;
  virtual su2double SpecificHeatVol(su2double T)const;
  virtual su2double EntropyTemp(su2double T)const;

public:

     /*!
     * \brief Constructor of the class.
     */
    CElectronicExcitation(void);

    /*!
     * \brief Constructor of the class.
     */
    CElectronicExcitation(su2double R, su2double g,
    		unsigned short n_mode, su2double* theta, su2double* w,
    		std::vector<su2double> gEl, std::vector<su2double> thetaEl);


    /*!
     * \brief Destructor of the class.
     */
    virtual ~CElectronicExcitation(void);


};
