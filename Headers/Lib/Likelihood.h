/********************************************************************
 *  Copyright (C) 2014 by Federico Marulli and Alfonso Veropalumbo  *
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
 *  @file Headers/Lib/Likelihood.h
 *
 *  @brief The class Likelihood 
 *
 *  This file defines the interface of the class Likelihood, used for
 *  statistical analyses and Bayesian inference
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#ifndef __LIKEL__
#define __LIKEL__

#include "Chi2.h"


// ============================================================================================


namespace cosmobl {

  /**
   *  @brief The namespace of functions and classes used for statistical
   *  analysis
   *  
   * The \e statistic namespace contains all the functions and classes
   * used for statistical analyis
   */

  namespace statistics {

    typedef function<double(double , shared_ptr<void>)> likelihood_1par;
    typedef function<double(vector<double> , shared_ptr<void> )> likelihood_npar;

    double likelihood_gaussian_1D_model_1par (double model_parameters, const shared_ptr<void> fixed_parameters);
    double likelihood_gaussian_1D_error_1par (double model_parameters, const shared_ptr<void> fixed_parameters);
    double likelihood_gaussian_1D_covariance_1par (double model_parameters, const shared_ptr<void> fixed_parameters);
    double likelihood_gaussian_2D_error_1par (double model_parameters, const shared_ptr<void> fixed_parameters);

    double likelihood_gaussian_1D_model_npar (vector<double> model_parameters, const shared_ptr<void> fixed_parameters);
    double likelihood_gaussian_1D_error_npar (vector<double> model_parameters, const shared_ptr<void> fixed_parameters);
    double likelihood_gaussian_1D_covariance_npar (vector<double> model_parameters, const shared_ptr<void> fixed_parameters);
    double likelihood_gaussian_2D_error_npar (vector<double> model_parameters, const shared_ptr<void> fixed_parameters);

    class Likelihood
    {
      protected:
	shared_ptr<Data> m_data;
	shared_ptr<Model> m_model; 

	likelihood_1par m_likelihood_1par;
	likelihood_npar m_likelihood_npar;
	int m_nchains;
	int m_chain_size;
	bool m_npar;

      public:

	/**
	 *  @name Constructors/destructors
	 */
	///@{

	/**
	 *  @brief default constructor
	 *  @return object of class Likelihood
	 */
	Likelihood () {}

	/**
	 *  @brief constructor
	 *  @param data pointers to the data container
	 *  @param model pointers to the model 
	 *  @param lik_1par function of type likelihood_1par
	 *  @return object of class Likelihood
	 */
	Likelihood (const shared_ptr<Data> data, const shared_ptr<Model> model, const likelihood_1par lik_1par) :
	  m_data(data), m_model(model), m_likelihood_1par(lik_1par), m_npar(0) {}

	/**
	 *  @brief constructor
	 *  @param data pointers to the data container
	 *  @param model pointers to the model 
	 *  @param lik_npar function of type likelihood_npar
	 *  @return object of class Likelihood
	 */
	Likelihood (const shared_ptr<Data> data, const shared_ptr<Model> model, const likelihood_npar lik_npar) :
	  m_data(data), m_model(model), m_likelihood_npar(lik_npar), m_npar(1) {}

	/**
	 *  @brief default destructor
	 *  @return none
	 */
	~Likelihood () {}

	///@}


	/**
	 *  @brief function that maximize Likelihood, find best fit parameters and store them in model
	 *  @param parameter starting value of the parameter 
	 *  @param max_iter maximum number of iteration 
	 *  @param min minumum value for minima finding
	 *  @param max maximum value for minima finding
	 *  @return none
	 */
	void minimize (const double parameter, const unsigned int max_iter=100, const double min=-1.e30, const double max=1.e30);

	/**
	 *  @brief function that maximize Likelihood, find best fit
	 *  parameters and store them in model       
	 *  @param parameters vector containing parameters starting
	 *  values
	 *  @param max_iter maximum number of iteration 
	 *  @param tol the tolerance for minima finding convergence 
	 *  @return none
	 */
	void minimize (const vector<double> parameters, const unsigned int max_iter=100, const double tol=1.e-6); 

	/**
	 *  @brief function that samples likelihood, using stretch-move
	 *  algorithm on n-dimensional parameter space, and stores chain
	 *  parameters.
	 *  @param nchains number of chains to sample the parameter space 
	 *  @param chain_size number of step in each chain 
	 *  @param seed the seed for random number generator
	 *  @return averace acceptance ratio
	 */
	double sample (const int nchains, const int chain_size, const int seed);

	/**
	 *  @brief function that samples likelihood, using
	 *  Metropolis-Hastings algorithm on uni-dimensional parameter
	 *  space, and stores chain parameters.       
	 *  @param nchains number of chains to sample the parameter space 
	 *  @param chain_size number of step in each chain 
	 *  @param shift  percentage shift for proposed parameter 
	 *  @param nsubstep number of substeps in each iteration 
	 *  @return averace acceptance ratio
	 */
	double sample (const int nchains, const int chain_size, const double shift, const int nsubstep = 100); 

	/**
	 *  @brief function that write sampled parameters 
	 *  @param output_file file of output for parameters 
	 *  @param start starting position in the chains
	 *  @param stop final position in the chains
	 *  @param thin interval of parameter in output
	 *  @return none
	 */
	void write_chain (const string output_file, const double start=0.5, const double stop=1, const int thin=1);
    };
  }
}

#endif
