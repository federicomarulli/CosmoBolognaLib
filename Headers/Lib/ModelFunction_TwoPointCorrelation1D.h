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
 *  @file Headers/Lib/ModelFunction_TwoPointCorrelation1D.h
 *
 *  @brief Global functions to model 1D two-point correlation functions
 *  of any type
 *
 *  This file contains all the prototypes of the global functions used
 *  to model 1D two-point correlation functions of any type
 *  
 *  @author Federico Marulli, Alfonso Veropalumbo
 *
 *  @author federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#ifndef __MODFUNCTWOP1D__
#define __MODFUNCTWOP1D__

#include "ModelFunction_TwoPointCorrelation.h"


// ============================================================================


namespace cosmobl {

  namespace modelling {

    namespace twopt {
      
      /**
       *  @struct STR_data_HOD
       *  @brief the STR_data_HOD structure
       *
       *  This structure contains the data used for the HOD modelling
       *  of the two-point correlation function
       */
      struct STR_data_HOD {
	
	/// cosmology
	shared_ptr<cosmology::Cosmology> cosmology;

	/// redshift
	double redshift = 0.;
	
	/// author(s) of the mass function
	string author_MF = "Tinker";
	
	/// author(s) of the bias
	string author_bias = "Tinker";

	/// minimum halo mass
	double Mh_min = 0.;
	
        /// maximum halo mass
	double Mh_max = 1.e16;
	
	/// the upper limit of the line-of-sight integration
	double pi_max = 100.;

	/// the maximum separation used to count pairs; it is used to compute the upper limit of the line-of-sight integration in approximate projected clustering estimators
	double r_max_int = 100.;
	
	/// minimum separation up to which the correlation function is computed
	double r_min = 1.e-3;

	/// maximum separation up to which the correlation function is computed
	double r_max = 350.;
	
	/// minimum wave vector module up to which the power spectrum is computed
	double k_min = 0.;
	 
	/// maximum wave vector module up to which the power spectrum is computed
	double k_max = 100.;
	
	/// number of steps used to compute the binned dark matter correlation function
	int step = 200;
	
	/// method used to compute the power spectrum and &sigma;(mass);
	string method_Pk = "CAMB";

	/// false &rarr; linear power spectrum; true &rarr; non-linear power spectrum
	bool NL = true;
	
	/// output_root of the parameter file used to compute the power spectrum 
	string output_root = "test";
	
	///  &Delta;: the overdensity, defined as the mean interior density relative to the background
	double Delta = 200.;

	/// wave vector module
	double kk = 0.;
	
	///  method to interpolate the power spectrum
	string interpType = "Linear";
	
	/// norm 0 &rarr; don't normalize the power spectrum; 1 &rarr; normalise the power spectrum
	int norm = -1;
	
	/// accuracy of the GSL integration 
	double prec = 1.e-2;
	
	/// name of eiter the parameter file or the power spectrum file; if !=NULL, it will be used, ignoring the cosmological parameters of the object
	string input_file = par::defaultString;

	/// if true the input_file is the parameter file, if false is the power spectrum file
	bool is_parameter_file = true;
	
	/// the author(s) of the concentration-mass relation (see cosmobl::modelling::twopt::concentration)
	string author_cM = "Duffy";

	/// the density profile, see cosmobl::modelling::twopt::concentration
	string profile = "NFW";
	
        /// the halo definition, see cosmobl::modelling::twopt::concentration
	string halo_def = "vir";
	
	/// pointer to a function of func_grid_GSL class, used to interpolate the power spectrum
	shared_ptr<glob::FuncGrid> func_Pk;

	/// pointer to a function of func_grid_GSL class, used to interpolate \f$\sigma(M)\f$
	shared_ptr<glob::FuncGrid> func_sigma;
	
	/// pointer to a function of func_grid_GSL class, used to interpolate \f$d\ln\sigma(M)/d\ln M\f$
	shared_ptr<glob::FuncGrid> func_dlnsigma;

	/**
	 * @brief default constructor
	 * @return object of type STR_data_HOD
	 */
	STR_data_HOD () = default;
	
      };
      
    }
  }
}

#endif
