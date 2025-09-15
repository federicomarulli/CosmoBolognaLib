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
 *  @file Headers/Modelling_PowerSpectrum_multipoles.h
 *
 *  @brief The class Modelling_PowerSpectrum_multipoles
 *
 *  This file defines the interface of the class
 *  Modelling_PowerSpectrum_multipoles, used to model the multipoles
 *  of power spectrum
 */

#ifndef __MODELLINGPOWSPECMULTI__
#define __MODELLINGPOWSPECMULTI__


#include "Sigma.h"
#include "Measure.h"
#include "Modelling.h"
#include "ModelFunction_PowerSpectrum_multipoles.h"


// ===================================================================================================


namespace cbl {

  namespace modelling {
    
    /**
     *  @brief The namespace of the <B> 3D power spectrum modelling
     *  </B>
     *  
     *  The \e modelling::powspect namespace contains all the
     *  functions and classes to model the 3D power spectrum
     */
    namespace powspec {
      
      /**
       *  @class Modelling_PowerSpectrum_multipoles
       *  Modelling_PowerSpectrum_multipoles.h
       *  "Headers/Modelling_PowerSpectrum_multipoles.h"
       *
       *  @brief The class Modelling_PowerSpectrum_multipoles
       *
       *  This file defines the interface of the base class
       *  Modelling_PowerSpectrum_multipoles, used for modelling the
       *  multipoles of power spectrum
       */
      class Modelling_PowerSpectrum_multipoles : public Modelling {
      protected:
	
	/// the container of parameters for angular power spectrum model computation
	std::shared_ptr<modelling::powspec::STR_data_model> m_data_model;
      public:

	/**
	 *  @name Constructors/destructors
	 */
	///@{

	/**
	 *  @brief default constuctor
	 */
	Modelling_PowerSpectrum_multipoles () = default;

	/**
	 *  @brief constructor
	 *  
	 *  @param PowSpM_dataset the dataset containing the
	 *  multipoles of the power spectrum to model
	 */
	Modelling_PowerSpectrum_multipoles (const std::shared_ptr<data::Data> PowSpM_dataset)
	{
	  m_data = PowSpM_dataset;
	};

	/**
	 *  @brief default destructor
	 *  
	 */
	virtual ~Modelling_PowerSpectrum_multipoles () = default;

	///@}
	
	/**
	 *  @name Member functions used to set the model parameters
	 */
	///@{

         /**
	 *  @brief set the model to fit the multipole of the power spectrum
	 *
	 *  the model is the following:
	 *
	 *  \f[ P_l(k) = \frac{1}{q_\perp^2 q_\parallel} \frac{2l+1}{2} \int_{-1}^{1} \mathrm{d}\mu P(k,
	 *  \mu) L_l(\mu) \f]
	 *	 
	 *  the redshift-space power spectrum
         *  \f[ P(k, \mu) =\left(1+\frac{f}{b}\mu^2\right)^2b^2P^{lin}(k) \f]
	 *
	 *  the model has N cosmological parameters
	 *    
	 *
	 *  @param cosmo_param list of cosmological parameter for the model
	 *
	 *  @param cosmo_prior prior for the cosmological parameters
	 *
	 *  
	 */
	void set_model_multipoles_Kaiser (const std::vector<std::string> cosmo_param, const std::vector<statistics::PriorDistribution> cosmo_prior);
	
	/**
	 *  @brief Set the data used to construct models of
	 *  the multipoles of the power spectrum
	 *
	 *  @param cosmology the cosmological model used to compute the angular power spectrum
	 *
	 *  @param redshifts the vector of the redshift
	 *
	 *  @param ll the vector of the order of multipoles
	 *
	 *  @param bias the vector of the fiducial value of the bias
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
	 *  @param store_output if true the output files created by the
         *  Boltzmann solver are stored; if false the output files are
         *  removed
	 *
	 *  @param k_min minimum wave vector module up to which the
	 *  power spectrum is computed in order to estimate the power
	 *  spectrum normalisation; this parameter is used only if
	 *  either norm=1, or norm=-1 and sigma8 is set
	 *
	 *  @param k_max maximum wave vector module up to which the
	 *  power spectrum is computed to estimate the power spectrum
	 *  normalisation; this parameter is used only if norm=1
	 *
	 *  @param output_root the parameter of the parameter file used 
	 *  to compute the power spectrum; it can be any name, If this 
	 *  parameter is different from default value it will be used also
	 *  in the output directory name
	 *
	 *  @param file_par name of the parameter file; if a parameter file
	 *  is provided (i.e. file_par!=NULL), it will be used, ignoring the
	 *  cosmological parameters of the object
	 *
	 *  @param prec accuracy of the integration
	 *
	 *  @param unit true \f$\rightarrow\f$ force cosmological units
	 *
	 *
	 */
	void set_data_model (const std::shared_ptr<cbl::cosmology::Cosmology> cosmology, const std::vector<double> redshifts, const std::vector<int> ll, const std::vector<double> bias, const std::string method_Pk="CAMB", const bool NL=false, const int norm=0, const bool store_output=false, const double k_min=0.001, const double k_max=100., const std::string output_root="test", const std::string file_par=par::defaultString, const double prec=1.e-5, const bool unit=false);
	
	///@}
         

      };
    }
  }
}

#endif
