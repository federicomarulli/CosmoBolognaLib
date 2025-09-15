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
 *  @file Headers/Modelling_PowerSpectrum_Shear.h
 *
 *  @brief The class Modelling_PowerSpectrum_Shear
 *
 *  This file defines the interface of the class
 *  Modelling_PowerSpectrum_Shear, used to model the shear-shear lensing
 *
 */

#ifndef __MODELLINGSSLENS__
#define __MODELLINGSSLENS__

#include "Modelling.h"
#include "ModelFunction_PowerSpectrum_Shear.h"


// ===================================================================================================


namespace cbl {

  namespace modelling {

    namespace sslensing {
      
      /**
       *  @class Modelling_PowerSpectrum_Shear
       *  Modelling_PowerSpectrum_Shear.h
       *  "Headers/Modelling_PowerSpectrum_Shear.h"
       *
       *  @brief The class Modelling_PowerSpectrum_Shear
       *
       *  This file defines the interface of the base class
       *  Modelling_PowerSpectrum_Shear, used for modelling
       *  the shear-shear lensing
       */
      class Modelling_PowerSpectrum_shear : public Modelling {
      protected:
	
	/// the container of parameters for shear-shear lensing model computation
	std::shared_ptr<modelling::sslensing::STR_data_model> m_data_model;
	
      public:

	/**
	 *  @name Constructors/destructors
	 */
	///@{

	/**
	 *  @brief default constuctor
	 */
	Modelling_PowerSpectrum_shear () = default;

	/**
	 *  @brief constructor
	 *  
	 *  @param dataset the dataset containing shear-shear
	 *  lensing to model
	 */
	Modelling_PowerSpectrum_shear (const std::shared_ptr<data::Data> dataset)
	{
	  m_data = dataset;
	};
	
	/**
	 *  @brief default destructor
	 *  
	 */
	virtual ~Modelling_PowerSpectrum_shear () = default;

	///@}
	
	/**
	 *  @name Member functions used to set the model parameters
	 */
	///@{
	
	/**
	 *  @brief Set the data used to construct models of
	 *  the shear-shear angular power spectrum.
	 *
	 *  @param cosmology the cosmological model used to compute the shear-shear angular power spectrum
	 *
	 *  @param z_min the minimum redshift of the photometric redshift slice
	 *
	 *  @param z_max the maximum redshift of the photometric redshift slice
	 *
	 *  @param Nz_all the redshift distributions for the computation of the shear-shear angular power spectrum
	 *
	 *  @param pars_IA vector that contains the parameters for the function (\f$ \mathcal{A_{\rm IA}} \f$, \f$ \mathcal{C_{\rm IA}} \f$ , \f$ \eta_{\rm IA}  \f$ , \f$ \beta_{\rm IA} \f$ ) 
	 *
	 *  @param L_ratio the interpolated ratio between the redshift-dependent mean and the charateristic luminosity
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
	 */
         void set_data_model(const std::shared_ptr<cbl::cosmology::Cosmology> cosmology, const double z_min, const double z_max, const std::vector<std::shared_ptr<glob::FuncGrid>> Nz_all, const std::vector<double> pars_IA, const cbl::glob::FuncGrid L_ratio, const std::string method_Pk="CAMB", const bool NL=false, const int norm=0, const double k_min=1e-4, const double k_max=1250., const std::string output_root="test", const std::string file_par=par::defaultString, const double prec=1.e-5, const bool unit=false);
 
        /**
         *  @brief Set the model to fit the angular power spectrum
         *
         *  the model is the following:
         *
         *\f[ C_{ij}^{\rm LL} = C_{ij}^{\gamma \gamma}(\ell) +
            C_{ij}^{\rm I\gamma}(\ell) + C_{ij}^{\rm II}(\ell)\f]
         *
         *  \f[ C_{ij}^{\rm \gamma \gamma}(\ell) =
              \frac{(\ell+2)!}{(\ell-2)!} 
              \left(\frac{2}{2\ell+1}\right)^4
               c \int {\rm d}z
              \frac{W_{i}^{\rm \gamma}(z)W_{j}^{\rm \gamma}(z)}
              {H(z)r^2(z)}
              P_{\rm \delta \delta}
              \left[ \frac{\ell +1/2}{r(z)}, z \right] \f]
         *
         *  \f[  C_{ij}^{\rm I\gamma}(\ell) =
              \frac{(\ell+2)!}{(\ell-2)!} 
              \left(\frac{2}{2\ell+1}\right)^4
              c \int {\rm d}z
              \frac{W_{i}^{\rm \gamma}(z)W_{j}^{\rm I}(z) +
                  W_{i}^{\rm I}(z)W_{j}^{\rm \gamma}(z)}
              {H(z)r^2(z)} f_{IA}(z)
              P_{\rm \delta \delta}
              \left[ \frac{\ell +1/2}{r(z)}, z \right] \f]
         *
         *  \f[  C_{ij}^{\rm II}(\ell) =
              \frac{(\ell+2)!}{(\ell-2)!}
              \left(\frac{2}{2\ell+1}\right)^4 c
              \int {\rm d}z
              \frac{W_{i}^{\rm I}(z)W_{j}^{\rm I}(z)}
              {H(z)r^2(z)} f^2_{IA}(z)
              P_{\rm \delta \delta}
              \left[ \frac{\ell +1/2}{r(z)}, z \right] \f]
         *
         *  the model has n parameters:
         *    - N cosmological parameters
         *    
         *  the angular power spectrum is computed
         *  using the input cosmological parameters
         *
         *
         *  @param cosmo_param list of cosmological parameter for the model
         *
         *  @param cosmo_prior prior for the cosmological parameters
         *
         */
   	 void set_model_Cij_LL_eNLA (const std::vector<std::string> cosmo_param, const std::vector<statistics::PriorDistribution> cosmo_prior);
          
	///@}	
      };
    }
  }
}

#endif
