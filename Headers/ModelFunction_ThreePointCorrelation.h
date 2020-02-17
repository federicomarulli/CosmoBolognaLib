/********************************************************************
 *  Copyright (C) 2017 by Federico Marulli                          *
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
 *  @file Headers/ModelFunction_ThreePointCorrelation.h
 *
 *  @brief Functions to model the three-point correlation function
 *
 *  This file contains the prototypes of the functions used to model
 *  the three-point correlation function
 *  
 *  @author Federico Marulli, Michele Moresco
 *
 *  @author federico.marulli3@unbo.it, michele.moresco@unibo.it
 */

#ifndef __MODFUNCTHREEP__
#define __MODFUNCTHREEP__

#include "Cosmology.h"


// ============================================================================


namespace cbl {

  namespace modelling {

    namespace threept {

      /**
       *  @struct STR_data_model_threept
       *
       *  @brief the structure STR_data_model_threept
       *
       *  This structure contains the data used for statistical
       *  analyses of the two-point correlation function
       */
      struct STR_data_model_threept {

	/// Q dark matter
	std::vector<double> Q_DM;

	/// cosmology
	std::shared_ptr<cosmology::Cosmology> cosmology;

	/// 1st side of the triangle
	double r1;

	/// 2nd side of the triangle
	double r2;

	/// theta
	std::vector<double> theta;

	/// model for the 3PCF
	std::string model;

	/// k vector
	std::vector<double> kk;

	/// Dark matter power spectrum
	std::vector<double> Pk_DM;

	/// the redshift
	double redshift;

	/// method to compute the dark matter power spectrum
	std::string method_Pk;

	/// false &rarr; linear power spectrum; true &rarr; non-linear power spectrum
	bool NL;
	
	/// minimum wave vector module up to which the power spectrum is computed
	double k_min;
	
	/// maximum wave vector module up to which the power spectrum is computed
	double k_max;
	
	/// number of steps used to compute the binned dark matter power spectrum
	int step_k;
	
	/// minimum separation up to which the binned dark matter correlation function is computed
	double r_min;
	
	/// maximum separation up to which the binned dark matter correlation function is computed
	double r_max;
	
	/// number of steps used to compute the binned dark matter correlation function
	int step_r;
	
	/// vector of scales
	std::vector<double> rr;
	
	/// the output_dir directory where the output of external codes are written
	std::string output_dir;

	/// true \f$\rightarrow\f$ the output files created by the Boltmann solver are stored; false \f$\rightarrow\f$ the output files are removed
	bool store_output;
	
	/// output root of the parameter file used to compute the dark matter power spectrum
	std::string output_root;

	/// 0 &rarr; don't normalize the power spectrum; 1 &rarr; normalize the power spectrum
	int norm;
	
	/// accuracy of the GSL integration
	double prec;
	
	/// 1 &rarr; do the computation in real; 0 &rarr; do the computation in redshift space
	bool force_realSpace;
	
	/// maximum order of the expansion
	int max_ll; 
	
        /// 1 &rarr; use the \f$k_l\f$ in the model; 0 &rarr; don't use the \f$k_l\f$ in the model
	bool use_k;

	/// &sigma;<SUB>8</SUB> at redshift z
	double sigma8_z;
	
	/// the linear growth rate at redshift z
	double linear_growth_rate_z;

	/**
	 *  @brief default constructor
	 *  @return object of type STR_data_model_threept
	 */
	STR_data_model_threept () = default;
      };

    }
  }
}

#endif
