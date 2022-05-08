/********************************************************************
 *  Copyright (C) 2016 by Federico Marulli and Alfonso Veropalumbo  *
 *  federico.marulli3@unibo.it                                      *
 *                                                                  *
 *  This program is free software; you can redistribute it and/or   *
 *  modify it under the terms of the GNU General Public License as  *
 *  published by the Free Software Foundation; either version 2 of  * 
 *  the License, or (at your option) any later version.             *
 *                                                                  *
 *  This program is distributed in the hope that it will be useful, *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of  *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the   *
 *  GNU General Public License for more details.                    *
 *                                                                  *
 *  You should have received a copy of the GNU General Public       *
 *  License along with this program; if not, write to the Free      *
 *  Software Foundation, Inc.,                                      *
 *  59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.       *
 ********************************************************************/

/**
 *  @file Headers/Modelling_PowerSpectrum_Angular.h
 *
 *  @brief The class Modelling_PowerSpectrum_angular
 *
 *  This file defines the interface of the class
 *  Modelling_PowerSpectrum_angular, used to model the angular power spectrum
 *
 *  @author Federico Marulli, Massimiliano Romanello
 *
 *  @author federico.marulli3@unibo.it, massimilia.romanell2@unibo.it
 */

#ifndef __MODELLINGPOWSPECTRUMANG__
#define __MODELLINGPOWSPECTRUMANG__

#include "Modelling.h"
#include "PowerSpectrum_Angular.h"
#include "ModelFunction_PowerSpectrum_Angular.h"

// ===================================================================================================


namespace cbl {

  namespace modelling {

    namespace angularpk {
      
      /**
       *  @class Modelling_PowerSpectrum_angular
       *  Modelling_PowerSpectrum_Angular.h
       *  "Headers/Modelling_PowerSpectrum_Angular.h"
       *
       *  @brief The class Modelling_PowerSpectrum_angular
       *
       *  This file defines the interface of the base class
       *  Modelling_PowerSpectrum_angular, used for modelling
       *  the angular power spectrum 
       */
      class Modelling_PowerSpectrum_angular : public Modelling {
      protected:
	
	/// the container of parameters for angular power spectrum model computation
	std::shared_ptr<modelling::angularpk::STR_data_model> m_data_model;

      public:

	/**
	 *  @name Constructors/destructors
	 */
	///@{

	/**
	 *  @brief default constuctor
	 */
	Modelling_PowerSpectrum_angular () = default;

	/**
	 *  @brief constructor
	 *  
	 *  @param Pow the angular power spectrum to model
	 */
	Modelling_PowerSpectrum_angular (const std::shared_ptr<cbl::measure::angularpk::PowerSpectrum_angular> Pow){
  		m_data = Pow->dataset();
		};

	/**
	 *  @brief constructor
	 *  
	 *  @param dataset the dataset containing the angular power spectrum to model
	 */
	Modelling_PowerSpectrum_angular (const std::shared_ptr<data::Data> dataset){
  		m_data = dataset;
		};

	/**
	 *  @brief default destructor
	 *  
	 */
	virtual ~Modelling_PowerSpectrum_angular () = default;

	///@}
	

	/**
	 *  @name Member functions used to set the model parameters
	 */
	///@{



	/**
	 *  @brief set the model to fit the angular power spectrum
	 *
	 *  the model is the following:
	 *
	 *  \f[ C_l = \frac{b^2}{N^2}\int_{0}^{\infty} \frac{dN}{dz}\frac{dN}{dz}
	 *  P_{mat} \left(\frac{l+1/2}{r(z)}\right) \frac{H(z)}{c}\frac{1}{r^2(z)} dz\f]
	 *
	 *  the model has 1 parameters: 
	 *    - \f$ b \f$
	 *    
	 *  offset and slope of the dN/dz distribution are fixed by the user
	 *  the angular power spectrum is computed
	 *  using the input cosmological parameters
	 *
	 *  @param bias_prior prior for the parameter \f$ b \f$
	 *
	 *  
	 */
	void set_model_limber (const statistics::PriorDistribution bias_prior);

	/**
	 *  @brief set the model to fit the angular power spectrum
	 *
	 *  the model is the following:
	 *
	 *  \f[ C_l = \frac{b^2}{N^2}\int_0^\infty \frac{dN}{dz} \frac{dN}{dz}
	 *  P_{mat} \left(\frac{l+1/2}{r(z)}\right) \frac{H(z)}{c}\frac{1}{r^2(z)} dz\f]
	 *
	 *  the model has 3 parameters: 
	 *    - \f$ b \f$
	 *    - \f$ offset \f$
	 *    - \f$ slope \f$
	 *
	 *  offset and slope reproduce the dN/dz distribution.
	 *  The angular power spectrum is computed
	 *  using the input cosmological parameters
	 *
	 *  @param bias_prior prior for the parameter \f$ b \f$
	 *  @param offset_slope_prior prior for the corralated parameter \f$ offset \f$ 
	 *  and \f$ slope \f$
	 *  
	 */
	void set_model_limber (const statistics::PriorDistribution bias_prior, const statistics::PriorDistribution offset_slope_prior);

	/**
	 * @brief get the member \e m_data_model
	 * @return the container of parameters for two-point
	 * correlation function model computation
	 */
	std::shared_ptr<modelling::angularpk::STR_data_model> data_model () { return m_data_model; }

	/**
	 *  @brief Set the data used to construct models of
	 *  the angular power spectrum.
	 *
	 *  @param cosmology the cosmological model used to measure angular power spectrum
	 *
	 *  @param z_min the minimum redshift of the integral
	 *
	 *  @param z_max the maximum redshift of the integral
	 *
	 *  @param method_Pk method used to compute the power spectrum
	 *  (i.e. the Boltzmann solver); valid choices for method_Pk
	 *  are: CAMB [http://camb.info/], CLASS
	 *  [http://class-code.net/], MPTbreeze-v1
	 *  [http://arxiv.org/abs/1207.1465], EisensteinHu
	 *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
	 *
	 *  @param NL 0 \f$\rightarrow\f$ linear power spectrum; 1 \f$\rightarrow\f$
	 *  non-linear power spectrum
	 *
	 *  @param norm 0 \f$\rightarrow\f$ don't normalise the power
	 *  spectrum; 1 \f$\rightarrow\f$ normalise the power spectrum;
	 *  -1 \f$\rightarrow\f$ normalise only if sigma8 is set
	 *
	 *  @param k_min minimum wave vector module up to which the
	 *  power spectrum is computed in order to estimate the power
	 *  spectrum normalisation; this parameter is used only if
	 *  either norm=1, or norm=-1 and sigma8 is set
	 *
	 *  @param k_max maximum wave vector module up to which the
	 *  power spectrum is computed to estimate the power spectrum
	 *  normalisation; this parameter is used only if norm=1
	 *  @param fsky the fraction of the sky covered by the survey
	 *
	 *  @param ll the vector of the multipoles of the mixing matrix
	 *
	 *  @param mixing_matrix the mixing matrix
	 *
	 *  @param dN_par offset and slope of the dN/dz normalized distribution
	 *  
	 */
	void set_data_model (const cbl::cosmology::Cosmology cosmology, const double z_min=0., const double z_max=10., const std::string method_Pk="CAMB", const bool NL=false, const int norm=-1, const double k_min=0.001, const double k_max=100., const std::vector<double> dN_par={}, const double fsky=1., std::vector<double> ll={}, std::vector<std::vector<double>> mixing_matrix={});
	
	///@}

      };
    }
  }
}

#endif
