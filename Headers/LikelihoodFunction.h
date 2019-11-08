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
 *  @file Headers/LikelihoodFunction.h
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

#include "Model2D.h"


// ===================================================================================================


namespace cbl {

  namespace statistics {
    
    /**
     * @enum LikelihoodType
     * @brief the type of likelihood function
     */
    enum class LikelihoodType {

      /// not set
      _NotSet_,

      /// Gaussian likelihood error
      _Gaussian_Error_,
      
      /// Gaussian likelihood covariance
      _Gaussian_Covariance_,

      /// Poissonian likelihood
      _Poissonian_,
      
      /// Likelihood function defined by the user
      _UserDefined_

    };

    /**
     * @brief return a vector containing the
     * LikelihoodType names
     * @return a vector containing the
     * LikelihoodType names
     */
    inline std::vector<std::string> LikelihoodTypeNames () {return {"NotSet", "Gaussian_Error", "Gaussian_ColikelihoodTypeiance", "Poissonian", "UserDefined"}; }

    /**
     * @brief cast an enum of type LikelihoodType
     * from its index
     * @param likelihoodTypeIndex the likelihoodType index
     * @return object of class LikelihoodType
     */
    inline LikelihoodType LikelihoodTypeCast (const int likelihoodTypeIndex) {return castFromValue<LikelihoodType>(likelihoodTypeIndex);}

    /**
     * @brief cast an enum of type LikelihoodType
     * from its name
     * @param likelihoodTypeName the likelihoodType name
     * @return object of class LikelihoodType
     */
    inline LikelihoodType LikelihoodTypeCast (const std::string likelihoodTypeName) {return castFromName<LikelihoodType>(likelihoodTypeName, LikelihoodTypeNames());}

    /**
     * @brief cast an enum of type LikelihoodType
     * from indeces
     * @param likelihoodTypeIndeces the likelihoodType indeces
     * @return object of class LikelihoodType
     */
    inline std::vector<LikelihoodType> LikelihoodTypeCast (const std::vector<int> likelihoodTypeIndeces) {return castFromValues<LikelihoodType>(likelihoodTypeIndeces);} 

    /**
     * @brief cast enums of type LikelihoodType
     * from thier names
     * @param likelihoodTypeNames the likelihoodType names
     * @return vector of LikelihoodType enums
     */
    inline std::vector<LikelihoodType> LikelihoodTypeCast (const std::vector<std::string> likelihoodTypeNames) {return castFromNames<LikelihoodType>(likelihoodTypeNames, LikelihoodTypeNames());}

    /**
     * @struct STR_likelihood_inputs
     * @brief the struct STR_likelihood_inputs
     *
     * This struct contains the data and the model for the likelihood
     * analysis
     */
    struct STR_likelihood_inputs
    {
      /// data containers
      std::shared_ptr<data::Data> data;

      /// model to test
      std::shared_ptr<Model> model;

      /// x position where the model is computed
      std::vector<double> xx;

      /// y position where the model is computed
      std::vector<double> yy;

      /// weight for the bin - 1D
      std::vector<double> weights1D;

      /// weight for the bin - 2D
      std::vector<std::vector<double>> weights2D;

      /// interpolated likelihood function - 1D
      std::shared_ptr<glob::FuncGrid> interp_function1D;

      /// interpolated likelihood function - 2D
      std::shared_ptr<glob::FuncGrid2D> interp_function2D;

      /**
       *  @brief constructor
       *  
       *  @param input_data pointers to the data container
       *  @param input_model pointers to the model 
       *  @param input_x_index vector contaning the x indeces
       *  @param input_w_index weight index
       *
       *  @return object of type STR_likelihood_inputs
       */ 
      STR_likelihood_inputs (const std::shared_ptr<data::Data> input_data, const std::shared_ptr<Model> input_model, const std::vector<size_t> input_x_index={0, 2}, const int input_w_index=-1); 
      
    };

    /**
     *  @brief function to compute the loglikelihood on a grid
     *
     *  @param likelihood_parameter the parameters of the model
     *
     *  @param input pointer to an object of type STR_params
     *
     *  @return the value of the loglikelihood 
     */
    double LogLikelihood_1D_interpolated (std::vector<double> &likelihood_parameter, const std::shared_ptr<void> input);

    /**
     *  @brief function to compute the loglikelihood on a grid
     *
     *  @param likelihood_parameter the parameters of the model
     *
     *  @param input pointer to an object of type STR_params
     *
     *  @return the value of the loglikelihood 
     */
    double LogLikelihood_2D_interpolated (std::vector<double> &likelihood_parameter, const std::shared_ptr<void> input);

    /**
     *  @brief function to compute the gaussian loglikelihood 
     *
     *  @param likelihood_parameter the parameters of the model
     *
     *  @param input pointer to an object of type STR_params
     *
     *  @return the value of the loglikelihood 
     */
    double LogLikelihood_Gaussian_1D_error (std::vector<double> &likelihood_parameter, const std::shared_ptr<void> input);

    /**
     *  @brief function to compute the gaussian loglikelihood 
     *
     *  @param likelihood_parameter the parameters of the model
     *
     *  @param input pointer to an object of type STR_params
     *
     *  @return the value of the loglikelihood 
     */
    double LogLikelihood_Gaussian_1D_covariance (std::vector<double> &likelihood_parameter, const std::shared_ptr<void> input);

    /**
     *  @brief function to compute the gaussian loglikelihood 
     *
     *  @param likelihood_parameter the parameters of the model
     *
     *  @param input pointer to an object of type STR_params
     *
     *  @return the value of the loglikelihood 
     *
     *  @warning the Gaussian likelihood does not allow for null
     *  standard deviation
     */
    double LogLikelihood_Gaussian_1D_error (std::vector<double> &likelihood_parameter, const std::shared_ptr<void> input);

    /**
     *  @brief function to compute the gaussian loglikelihood 
     *
     *  @param likelihood_parameter the parameters of the model
     *
     *  @param input pointer to an object of type STR_params
     *
     *  @return the value of the loglikelihood
     */
    double LogLikelihood_Gaussian_1D_covariance (std::vector<double> &likelihood_parameter, const std::shared_ptr<void> input);

    /**
     *  @brief function to compute the gaussian loglikelihood model
     *  with one parameter \f$ \chi^2 \f$
     *
     *  @param likelihood_parameter the parameters of the model
     *
     *  @param input pointer to an object of type STR_params
     *
     *  @return the value of the loglikelihood 
     *
     *  @warning the Gaussian likelihood does not allow for null
     *  standard deviation
     */
    double LogLikelihood_Gaussian_2D_error (std::vector<double> &likelihood_parameter, const std::shared_ptr<void> input);

    /**
     *  @brief function to compute the poissonian loglikelihood 
     *
     *  @param likelihood_parameter the parameters of the model
     *
     *  @param input pointer to an object of type STR_params
     *
     *  @return the value of the loglikelihood 
     *
     *  @warning the poissian likelihood takes counts in input; it's
     *  care of the user to ensure that this is the case
     *
     */
    double LogLikelihood_Poissonian_1D_ (std::vector<double> &likelihood_parameter, const std::shared_ptr<void> input);

    /**
     *  @brief function to compute the poissonian loglikelihood
     *
     *  @param likelihood_parameter the parameters of the model
     *
     *  @param input pointer to an object of type STR_params
     *
     *  @return the value of the loglikelihood 
     *
     *  @warning the poissian likelihood takes counts in input; it's
     *  care of the user to ensure that this is the case
     */
    double LogLikelihood_Poissonian_2D_ (std::vector<double> &likelihood_parameter, const std::shared_ptr<void> input);

    /**
     * @var typedef LogLikelihood_function
     *
     * @brief definition of a function for computation of the
     * LogLikelihood
     */
    typedef std::function<double (std::vector<double> &, const std::shared_ptr<void>)> LogLikelihood_function;

    /**
     * @var typedef Likelihood_function
     *
     * @brief definition of a function for computation of the
     * Likelihood
     */
    typedef std::function<double (std::vector<double> &, const std::shared_ptr<void>)> Likelihood_function;
  }
}

#endif
