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
 *  @file Headers/ModelFunction_PowerSpectrum_Angular.h
 *
 *  @brief Global functions to model the angular power spectrum
 *
 *  This file contains all the prototypes of the functions used
 *  to model the angular power spectrum
 *  
 *  @author Federico Marulli, Massimiliano Romanello
 *
 *  @author federico.marulli3@unibo.it, massimilia.romanell2@unibo.it
 */

#ifndef __MODFUNCPOWSPECTRUMANG__
#define __MODFUNCPOWSPECTRUMANG__

#include "Cosmology.h"

// ============================================================================


namespace cbl {

  namespace modelling {

    namespace angularpk {

      struct STR_data_model {
  
	std::shared_ptr<cosmology::Cosmology> cosmology;

	std::vector<cosmology::CosmologicalParameter> Cpar;

	double z_min;
	     	     
	double z_max;

	double z_min_bin2;
	     	     
	double z_max_bin2;
  
	std::string method_Pk;
        
	bool NL;
        
	int norm;

	double k_min;
	     	     
	double k_max;
	     
	double fsky;
	     	     	     	     
	std::vector<std::vector<double>> mixing_matrix;
	     
	std::vector<double> ll;

	bool interpolate_power_spectrum;
	     
	std::vector<double> dN_par;
	      
	std::vector<double> dN_par_bin2;

	STR_data_model () = default;
      };

      /**
       *  @brief the angular power spectrum convolved with the mixing matrix 
       *
       *  the function computes:
       *
       *  \f[C_l^{mixed}=\sum_{l'} R_{ll'}C_{l'}\f]
       *
       *  @param l_mixing the vector of multipoles of the mixing matrix
       *
       *  @param mixing_matrix the mixing matrix
       *
       *  @param l the vector of multipoles of the power spectrum
       *
       *  @param Cl the angular power spectrum
       *
       *  @param fsky the fraction of sky covered by the survey
       *
       *  @return the angular power spectrum convolved with the mixing matrix
       */
      std::vector<double> Cl_mixed(std::vector<double> l_mixing, std::vector<std::vector<double>> mixing_matrix, std::vector<double> l, std::vector<double> Cl, double fsky);
      
      /**
       *  @brief the integrand function in the Limber approximation for 
       *  the calculus of angular power spectrum 
       *
       *  the function computes:
       *
       *  \f[\frac{1}{N^2}\frac{dN}{dz}\frac{dN}{dz}
       *  P_{mat}\left(\frac{l+1/2}{r(z)}\right) \frac{H(z)}{c}\frac{1}{r^2(z)}\f]
       *
       *  @param redshift redshift
       *
       *  @param pp pointer to the structure that contains the
       *  power spectrum angular data model
       *
       *  @param par vector containing the input parameters
       *
       *  @return the integrand function of the angular power spectrum in the Limber approximation
       */
      double integrand_limber_exact (double redshift,  std::shared_ptr<void> pp, std::vector<double> par);  

      /**
       *  @brief the integrand function in the Limber approximation for 
       *  the calculus of angular power spectrum 
       *
       *  the function computes:
       *
       *  \f[C_l = \frac{1}{N^2}\int_{0}^{\infty} \frac{dN}{dz}\frac{dN}{dz}
       *  P_{mat} \left(\frac{l+1/2}{r(z)}\right) \frac{H(z)}{c}\frac{1}{r^2(z)} dz\f]
       *
       *  @param l the multipole at which angular power spectrum is computed
       *
       *  @param z_vector the redshift vector of the interpolation grid
       *
       *  @param kk the k vector of the interpolation grid
       *
       *  @param pk_interp the interpolation grid object
       *   
       *  @param parameter the list of parameters, containing the cosmological parameters
       *
       *  @param input pointer to the structure that contains the
       *  power spectrum angular data model
       *
       *  @return the integral of the angular power spectrum in the Limber approximation
       */
      double integral_limber_interp (double l, std::vector<double> z_vector, std::vector<double> kk, cbl::glob::FuncGrid2D pk_interp, std::vector<double> parameter, std::shared_ptr<void> input);

      /**
       *  @brief the integrand function in the Limber approximation for 
       *  the calculus of angular power spectrum 
       *
       *  the function computes:
       *
       *  \f[C_l = \frac{1}{N^2}\int_{0}^{\infty} \frac{dN}{dz}\frac{dN}{dz}
       *  P_{mat} \left(\frac{l+1/2}{r(z)}\right) \frac{H(z)}{c}\frac{1}{r^2(z)} dz\f]
       *
       *  @param l the multipole at which angular power spectrum is computed
       *
       *  @param parameter the list of parameters, containing the cosmological parameters and eventually the offset and slope of the normalized dN/dz distribution
       *
       *  @param input pointer to the structure that contains the
       *  power spectrum angular data model
       *
       *  @return the integral of the angular power spectrum in the Limber approximation
       */
      double integral_limber_exact (double l, std::vector<double> parameter, std::shared_ptr<void> input);


      /**
       *  @brief the model for the angular power spectrum
       *
       *  the model is the following:
       *
       *  \f[ C_l = \frac{b^2}{N^2}\int_{0}^{\infty} \frac{dN}{dz}\frac{dN}{dz}
       *  P_{mat} \left(\frac{l+1/2}{r(z)}\right) \frac{H(z)}{c}\frac{1}{r^2(z)} dz\f]
       *
       *  the model has n+3 parameters:
       *    - \f$ n cosmological parameters, including \Omega_m, \Omega_b and \sigma_8 \f$
       *    - \f$ b \f$
       *    - \f$ offset \f$
       *    - \f$ slope \f$
       *    
       *  the angular power spectrum is computed
       *  using the input cosmological parameters
       *
       *  @param l the vector of multipoles at which the model is computed
       *  @param inputs pointer to the structure that contains the power spectrum angular data model
       *
       *  @param parameter 1D vector containing the linear bias and the offset and the slope of the normalized dN/dz distribution
       *  
       *  @return the angular power spectrum in the Limber approximation
       */
      std::vector<double> Cl_limber (const std::vector<double> l, const std::shared_ptr<void> inputs, std::vector<double> &parameter);
       
    }
  }
}

#endif
