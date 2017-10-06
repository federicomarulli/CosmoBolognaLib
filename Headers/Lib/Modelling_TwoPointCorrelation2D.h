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
 *  @file Headers/Lib/Modelling_TwoPointCorrelation2D.h
 *
 *  @brief The class Modelling_TwoPointCorrelation2D
 *
 *  This file defines the interface of the class
 *  Modelling_TwoPointCorrelation2D, that contains all the common
 *  methods to model 2D two-point correlation functions
 *
 *  @author Federico Marulli, Alfonso Veropalumbo
 *
 *  @author federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#ifndef __MODELLINGTWOPOINT2D__
#define __MODELLINGTWOPOINT2D__


#include "Modelling_TwoPointCorrelation.h"


// ===================================================================================================


namespace cosmobl {

  namespace modelling {

    namespace twopt {

      /**
       *  @class Modelling_TwoPointCorrelation2D
       *  Modelling_TwoPointCorrelation2D.h
       *  "Headers/Lib/Modelling_TwoPointCorrelation2D.h"
       *
       *  @brief The class Modelling_TwoPointCorrelation2D
       *
       *  This file defines the interface of the base class
       *  Modelling_TwoPointCorrelation2D, used for modelling
       *  the 2D two-point correlation function in cartesian coordinates
       *
       */
      class Modelling_TwoPointCorrelation2D : public Modelling_TwoPointCorrelation {

      protected:
	
	/// \f$\alpha_\perp\f$
	shared_ptr<statistics::Parameter> m_alpha_perp;

	/// \f$\alpha_\parallel\f$
	shared_ptr<statistics::Parameter> m_alpha_par;

	/// \f$b\sigma_8\f$
	shared_ptr<statistics::Parameter> m_bsigma8;

	/// \f$f\sigma_8\f$
	shared_ptr<statistics::Parameter> m_fsigma8;

	/// \f$f\sigma_{12}\f$
	shared_ptr<statistics::Parameter> m_sigma12;

	
      public:

	/**
	 *  @name Constructors/destructors
	 */
	///@{

	/**
	 *  @brief default constuctor
	 *  @return object of class Modelling_TwoPointCorrelation2D
	 */
	Modelling_TwoPointCorrelation2D () = default;

	/**
	 *  @brief constructor
	 *  
	 *  @param twop the two-point correlation function to model
	 *
	 *  @return object of type Modelling_TwoPointCorrelation2D
	 */
	Modelling_TwoPointCorrelation2D (const shared_ptr<cosmobl::measure::twopt::TwoPointCorrelation> twop)
	  : Modelling_TwoPointCorrelation(twop) {}

	/**
	 *  @brief default destructor
	 *  @return none
	 */
	virtual ~Modelling_TwoPointCorrelation2D () = default;

	///@}


	/**
	 *  @name Member functions used to set the model parameters
	 */
	///@{
	
	/**
	 *  @brief set the scale range used for the fit
	 *
	 *  @param xmin the minimum x value
	 *  @param xmax the maximum x value
	 *  @param ymin the minimum y value
	 *  @param ymax the maximum y value
	 *  @return none
	 */
	void set_fit_range (const double xmin, const double xmax, const double ymin, const double ymax)
	{ m_data_fit = m_data->cut(xmin, xmax, ymin, ymax); m_fit_range = true; }

	/**
	 *  @brief set the parameters for the computation of the dark
	 *  matter two-point correlation function
	 *  
	 *  @param cosmology the cosmology used
	 *
	 *  @param redshift redshift
	 *
	 *  @param method_Pk method used to compute the power spectrum
	 *  and &sigma;(mass); valid choices for method_Pk are: CAMB
	 *  [http://camb.info/], classgal_v1 [http://class-code.net/],
	 *  MPTbreeze-v1 [http://arxiv.org/abs/1207.1465],
	 *  EisensteinHu
	 *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
	 *
	 *  @param sigmaNL damping of the wiggles in the linear power
	 *  spectrum
	 *
	 *  @param NL 0 &rarr; linear power spectrum; 1 &rarr;
	 *  non-linear power spectrum
	 *
	 *  @param FV 0 &rarr; exponential form for f(v); 1 &rarr;
	 *  Gaussian form for f(v); where f(v) is the velocity
	 *  distribution function
	 *
	 *  @param output_root output_root of the parameter file used
	 *  to compute the power spectrum and &sigma;(mass); it can be
	 *  any name
	 *
	 *  @param bias_nl 0 &rarr; linear bias; 1 &rarr; non-linear
	 *  bias
	 *
	 *  @param bA b<SUB>a</SUB> non-linear bias parameter
	 *
	 *  @param xiType 0 &rarr; standard; 1 &rarr; Chuang & Wang
	 *  model
	 *
	 *  @param k_star k<SUB>*</SUB> of the Chuang & Wang model
	 *
	 *  @param xiNL 0 &rarr; linear power spectrum; 1 &rarr;
	 *  non-linear power spectrum
	 *
	 *  @param v_min minimum velocity used in the convolution of
	 *  the correlation function
	 *
	 *  @param v_max maximum velocity used in the convolution of
	 *  the correlation function
	 *
	 *  @param step_v number of steps used in the convolution of
	 *  the correlation function
	 *
	 *  @param norm 0 &rarr; don't normalize the power spectrum; 1
	 *  &rarr; normalize the power spectrum
	 *
	 *  @param r_min minimum separation up to which the binned
	 *  dark matter correlation function is computed
	 *
	 *  @param r_max maximum separation up to which the binned
	 *  dark matter correlation function is computed
	 *
	 *  @param k_min minimum wave vector module up to which the
	 *  binned power spectrum is computed
	 *
	 *  @param k_max maximum wave vector module up to which the
	 *  binned power spectrum is computed
	 *
	 *  @param step number of steps used to compute the binned
	 *  dark matter correlation function
	 *
	 *  @param aa parameter \e a of Eq. 24 of Anderson et al. 2012
	 *
	 *  @param GSL 0 &rarr; the Numerical libraries are used; 1
	 *  &rarr; the GSL libraries are used
	 *
	 *  @param prec accuracy of the GSL integration
	 *
	 *  @param file_par name of the parameter file; if a parameter
	 *  file is provided (i.e. file_par!=NULL), it will be used,
	 *  ignoring the cosmological parameters of the object
	 *
	 *  @return none
	 */
	void set_data_model (const cosmology::Cosmology cosmology={}, const double redshift=0., const string method_Pk="CAMB", const double sigmaNL=0, const bool NL=true, const int FV=0, const string output_root="test", const bool bias_nl=false, const double bA=-1., const bool xiType=false, const double k_star=-1., const bool xiNL=false, const double v_min=-5000., const double v_max=5000., const int step_v=500, const int norm=-1, const double r_min=0.1, const double r_max=150., const double k_min=0., const double k_max=100., const int step=200, const double aa=0., const bool GSL=true, const double prec=1.e-2, const string file_par=par::defaultString);
	      
	///@}
	

	/**
	 *  @name Member functions used to write the outputs
	 */
	///@{

	/**
	 *  @brief compute and write the model 
	 *
	 *  @param output_dir the output directory
	 *
	 *  @param output_file the output file
	 *
	 *  @param xx vector of points at which the model is computed,
	 *  first axis
	 *
	 *  @param yy vector of points at which the model is computed,
	 *  second axis
	 *
	 *  @param parameter vector containing the input parameters
	 *  used to compute the model; if this vector is not provided,
	 *  the model will be computed using the best-fit parameters
	 *
	 *  @param start the starting position for each chain
	 *
	 *  @param thin the position step
	 *
	 *  @return none
	 */
	void write_model (const string output_dir, const string output_file, const vector<double> xx={}, const vector<double> yy={}, const vector<double> parameter={}, const int start=0, const int thin=10) override;
	
	/**
         *  @brief compute the best-fit model from MCMC chains
         *
         *  @param xx vector of points at which the model is computed,
         *  first axis
	 *
	 *  @param yy vector of points at which the model is computed,
	 *  second axis
	 *
	 *  @param start the starting position for each chain
	 *  @param thin the position step
         *  @param median_model the median model
         *  @param low_model the model at 16th percentile
         *  @param up_model the model at 84th percentile
         *  @return none
         */
        void compute_model_from_chains (const vector<double> xx, const vector<double> yy, const int start, const int thin, vector<vector<double>> &median_model, vector<vector<double>> &low_model, vector<vector<double>> &up_model) override;

	
	///@}

      };
    }
  }
}

#endif
