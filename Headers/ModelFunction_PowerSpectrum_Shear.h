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
 *  @file Headers/ModelFunction_PowerSpectrum_Shear.h
 *
 *  @brief Global functions to model the shear-shear angular power spectrum
 *
 *  This file contains all the prototypes of the functions used
 *  to model the shear-shear angular power spectrum
 *
 */

#ifndef __MODFUNCSSLENS__
#define __MODFUNCSSLENS__

#include "Modelling_MassObservableRelation.h"
#include "Modelling_PowerSpectrum_Angular.h"

// ============================================================================


namespace cbl {

  namespace modelling {

    /**
     *  @brief The namespace of the <B> shear-shear lensing
     *  modelling </B>
     *  
     *  The \e modelling::sslensing namespace contains all the
     *  functions and classes to model the shear-shear lensing
     */
    namespace sslensing {

      using ProxyFunction = std::function<double(double, double, std::shared_ptr<void>)>;
      
      struct STR_data_model {

	/// the cosmological model
	std::shared_ptr<cosmology::Cosmology> cosmology;

	/// cosmological parameters
	std::vector<std::string> Cpar;

	/// minimum redshift
	double z_min;

	/// maximum redshift
	double z_max;

	/// method used to compute the power spectrum (i.e. the Boltzmann solver)
	std::string method_Pk;

	/// linear or non-linear power spectrum
	bool NL;
	
	/// output root of the parameter file used to compute the dark matter power spectrum
	std::string output_root;

	/// norm 0 \f$\rightarrow\f$ don't normalise the power spectrum; 1 \f$\rightarrow\f$ normalise the power spectrum; -1 \f$\rightarrow\f$ normalise only if sigma8 is set
	int norm;

	/// minimum wave vector module up to which the power spectrum is computed
	double k_min;

	/// maximum wave vector module up to which the power spectrum is computed
	double k_max;

	/// name of the parameter file
	std::string file_par;
	
	/// accuracy of the integration 
	double prec;
	
	/// unit
	bool unit;

	/// pointer to a vector of FuncGrid class,
	std::vector<std::shared_ptr<glob::FuncGrid>> Nz_all;
	
	/// vector that contains the parameters for IA (\f$ \mathcal{A_{\rm IA}} \f$, \f$ \mathcal{C_{\rm IA}} \f$ , \f$ \eta_{\rm IA}  \f$ , \f$ \beta_{\rm IA} \f$ ) 
	std::vector<double> pars_IA;
	
	/// Luminosity ratio for the IA
	cbl::glob::FuncGrid L_ratio;

	
	STR_data_model () = default;
      };

      /**
       * @brief compute the clustering window function for the IA
       *
       *  \f[ W_i^{IA}(z) = \frac{H(z)}{c} n_{i}(z) \f]
       *
       * @param zz the redshift
       *
       * @param cosmo the cosmology
       *
       * @param nz the interpolated redshift distribution
       *
       * @return the clustering window function for the IA at given redshift
       *
       */ 
       double window_function_IA (const double zz, const std::shared_ptr<cbl::cosmology::Cosmology> cosmo, const std::shared_ptr<cbl::glob::FuncGrid> nz);
      
      /**
       * @brief compute the clustering window function for the shear lensing 
       *
       *  \f[ W_{i}^{\gamma}( z) =
            \frac{3}{2}\left ( \frac{H_0}{c}\right )^2
            \Omega_{{\rm m},0} (1 + z)r(z)
            \int_{z}^{z_{max}}{{\rm d}z^{\prime} n_{i}(z^{\prime})
            \frac{r(z^{\prime}) - r(z)}
            {r(z^{\prime})}} \f]
       *
       * @param zz the redshift
       *
       * @param z_max the z_max
       *
       * @param cosmo the cosmology
       *
       * @param nz the the interpolated redshift distribution
       *
       * @param DC the interpolated function of the comoving radial distance 
       *
       * @return the clustering window function for the shear lensing at a given redshift
       *
       */ 
       double window_function_shear (const double zz, const double z_max, const std::shared_ptr<cbl::cosmology::Cosmology> cosmo, const std::shared_ptr<cbl::glob::FuncGrid> nz, cbl::glob::FuncGrid DC );
      
      /**
       * @brief compute the term responsible for the Intrinsic Alignment 
       *
       *  \f[ f_{\rm IA}(z) = -\mathcal{A_{\rm IA}}\mathcal{C_{\rm IA}}\
            \frac{\Omega_{m,0}}{D(z)}\
            [(1 + z)/(1 + z_{\rm pivot})]^{\eta_{\rm IA}}\
            [\langle L \rangle(z) /L_{\star}(z)]^{\beta_{\rm IA}} \f]
       *
       * @param zz the redshift
       *
       * @param Om_0 Omega matter compute at z=0
       *
       * @param z_p the pivot redshift 
       *
       * @param pars vector that contains the parameters for the function (\f$ \mathcal{A_{\rm IA}} \f$, \f$ \mathcal{C_{\rm IA}} \f$ , \f$ \eta_{\rm IA}  \f$ , \f$ \beta_{\rm IA} \f$ ) 
       *
       * @param DN the interpolated function of the growth factor
       *
       * @param L_ratio the interpolated ratio between the redshift-dependent mean and the charateristic luminosity
       *
       * @return the term of the IA
       *
       */       
       double Intrinsic_Alignment (const double zz,const double Om_0,const double z_p, const std::vector<double> pars, cbl::glob::FuncGrid DN, cbl::glob::FuncGrid L_ratio);
      
      /**
       *  @brief the model for the angular power spectrum
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
       *  @param l the vector of multipoles at which the model is computed
       *
       *  @param inputs pointer to the structure that contains the power spectrum angular data model
       *
       *  @param parameter 1D vector containing the linear bias and the offset and the slope of the normalized dN/dz distribution
       *  
       *  @return the angular power spectrum 
       *
       */      
       std::vector<double> C_ij_LL_eNLA (const std::vector<double> l, const std::shared_ptr<void> inputs, std::vector<double> &parameter);

    }
  }
}

#endif
