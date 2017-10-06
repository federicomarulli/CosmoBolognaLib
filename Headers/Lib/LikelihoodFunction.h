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
 *  @file Headers/Lib/LikelihoodFunction.h
 *
 *  @brief Likelihood function
 *
 *  This file defines the interface for likelihood function, used for
 *  statistical analyses and Bayesian inference
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#ifndef __LIKELFUN__
#define __LIKELFUN__

#include "LikelihoodParameters.h"


// ===================================================================================================


namespace cosmobl {

  namespace statistics {
    
    /**
     * @enum LikelihoodType
     * @brief the type of likelihood function
     */
    enum LikelihoodType {

      /// not set
      _Likelihood_NotSet_,

      /// Gaussian likelihood error
      _GaussianLikelihood_Error_,
      
      /// Gaussian likelihood covariance
      _GaussianLikelihood_Covariance_,
      
      /// Likelihood function defined by the user
      _UserDefinedLikelihood_

    };

    
    /**
     * @struct STR_likelihood_parameters
     * @brief the struct STR_likelihood_parameters
     *
     * This struct contains the data
     * and the model for the likelihood analysis
     */
    struct STR_likelihood_parameters
    {
      /// data containers
      shared_ptr<data::Data> data;

      /// model to test
      shared_ptr<Model> model;

      /// index(s) of the extra info vector containing the point(s) where to evaluate the model
      vector<int> x_index;

      /**
       *  @brief constructor
       *  @param _data pointers to the data container
       *  @param _model pointers to the model 
       *  @param _x_index vector contaning the x indeces
       *  @return object of type STR_likelihood_parameters
       */ 
      STR_likelihood_parameters (const shared_ptr<data::Data> _data, const shared_ptr<Model> _model, const vector<int> _x_index={})
      : data(_data), model(_model)
      {
	x_index.emplace_back(_x_index.size()>0 ? _x_index[0] : 0);
	x_index.emplace_back(_x_index.size()>1 ? _x_index[1] : 2);
      }
      
    };

    /** 
     *  @brief function to compute the gaussian loglikelihood 
     *  @param likelihood_parameters the parameters of the model
     *  @param inputs pointer to an object of type STR_params
     *  @return the value of the loglikelihood 
     */
    double LogLikelihood_Gaussian_1D_error (vector<double> &likelihood_parameters, const shared_ptr<void> inputs);

    /** 
     *  @brief function to compute the gaussian loglikelihood 
     *  @param likelihood_parameters the parameters of the model
     *  @param inputs pointer to an object of type STR_params
     *  @return the value of the loglikelihood 
     */
    double LogLikelihood_Gaussian_1D_covariance (vector<double> &likelihood_parameters, const shared_ptr<void> inputs);

    /** 
     *  @brief function to compute the gaussian loglikelihood 
     *  model with one parameter  &chi;&sup2; 
     *  @param likelihood_parameters the parameters of the model
     *  @param inputs pointer to an object of type STR_params
     *  @return the value of the loglikelihood 
     */
    double LogLikelihood_Gaussian_2D_error (vector<double> &likelihood_parameters, const shared_ptr<void> inputs);

    /**
     * @var typedef LogLikelihood_function
     * @brief definition of a function for computation of 
     * the LogLikelihood
     */
    typedef function<double (vector<double> &, const shared_ptr<void>)> LogLikelihood_function;
  }
}

#endif
