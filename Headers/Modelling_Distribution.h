/********************************************************************
 *  Copyright (C) 2021 by Federico Marulli and Giorgio Lesci        *
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
 *  @file Headers/Modelling_Distribution.h
 *
 *  @brief The class Modelling_Distribution
 *
 *  This file defines the interface of the class Modelling_Distribution, used to
 *  model any kind of distributions
 *
 *  @authors Giorgio Lesci (and Federico Marulli)
 *
 *  @authors giorgio.lesci2@unibo.it (and federico.marulli3@unibo.it)
 */

#ifndef __MODELLINGDISTR__
#define __MODELLINGDISTR__

#include "Cosmology.h"
#include "Measure.h"
#include "Modelling.h"
#include "Cluster.h"


// ===================================================================================================


namespace cbl {

  namespace modelling {

    /**
     *  @brief The namespace of the <B> distribution
     *  modelling </B>
     *  
     *  The \e modelling::distribution namespace contains all the functions
     *  to model distributions
     */
    namespace distribution {
    
      /**
       *  @struct STR_Distr_model
       *  @brief the structure STR_Distr_model
       *
       *  This structure contains the fixed parameters used for statistical
       *  analyses of generic distributions
       */
      struct STR_Distr_model {
      
        /// constant
        double k;

	/**
	 *  @brief default constructor
	 */
	STR_Distr_model () = default;
      };
    
      /**
       *  @class Modelling_Distribution
       *  Modelling_Distribution.h
       *  "Headers/Modelling_Distribution.h"
       *
       *  @brief The class Modelling_Distribution
       *
       *  This file defines the interface of the base class
       *  Modelling_Distribution, used for modelling
       *  any kind of statistical distribution.
       *
       */
      class Modelling_Distribution : public Modelling {
      
      protected:
      
        /// the container of fixed parameters
	STR_Distr_model m_data_model;

      public:

	/**
	 *  @name Constructors/destructors
	 */
	///@{

	/**
	 *  @brief default constuctor
	 *  
	 */
	Modelling_Distribution () = default;
	
	/**
	 *  @brief constuctor
	 *  @param dataset the dataset containing x, data and errors
	 */
	Modelling_Distribution (const std::shared_ptr<cbl::data::Data> dataset)
	{ m_data = dataset; }
	
	/**
	 *  @brief default destructor
	 *  
	 */
	virtual ~Modelling_Distribution () = default;

	///@}	
	
	/**
	 *  @name Member functions used to set the model parameters
	 */
	///@{

	/**
	 *  @brief set the parameters of a Gaussian PDF
	 *
	 *  @param mean_prior prior on the mean
	 *
	 *  @param std_prior prior on the standard deviation
	 *
	 *  @param mean_name string identifying the mean
	 *
	 *  @param std_name string identifying the standard deviation
	 *
	 */
	void set_model_Distribution (const statistics::PriorDistribution mean_prior, const statistics::PriorDistribution std_prior, const std::string mean_name, const std::string std_name);
	
	/**
	 *  @brief set the parameters of a Gaussian PDF, where
	 *  the standard deviation, \f$\sigma\f$, is expressed as
	 *  
	 *  \f$\sigma=\sigma_0k,\f$
	 *
	 *  where \f$k\f$ is fixed, \f$\sigma_0\f$ is a base parameter,
	 *  and \f$\sigma\f$ is a derived parameter.
	 *
	 *  @param k the value of \f$k\f$
	 *
	 *  @param mean_prior prior on the mean
	 *
	 *  @param std0_prior prior on \f$\sigma_0\f$
	 *
	 *  @param mean_name string identifying the mean
	 *
	 *  @param std0_name string identifying \f$\sigma_0\f$
	 *
	 *  @param std_name string identifying \f$\sigma\f$
	 *
	 */
	void set_model_Distribution (const double k, const statistics::PriorDistribution mean_prior, const statistics::PriorDistribution std0_prior, const std::string mean_name, const std::string std0_name, const std::string std_name);

	///@}
      };

      
      /**
       * @brief compute a Gaussian PDF
       *
       * @param x the points where the distribution is evaluated
       *
       * @param inputs model inputs
       *
       * @param parameter model parameters
       *
       * @return the Gaussian PDF values
       *
       */
      std::vector<double> model_gaussian (const std::vector<double> x, const std::shared_ptr<void> inputs, std::vector<double> &parameter);
      
      /**
       * @brief compute a Gaussian PDF, where
       * the standard deviation, \f$\sigma\f$, is expressed as
       *  
       * \f$\sigma=\sigma_0k,\f$
       *
       * where \f$k\f$ is fixed, \f$\sigma_0\f$ is a base parameter,
       * and \f$\sigma\f$ is a derived parameter.
       *
       * @param x the points where the distribution is evaluated
       *
       * @param inputs model inputs
       *
       * @param parameter model parameters
       *
       * @return the Gaussian PDF values
       *
       */
      std::vector<double> model_gaussian2 (const std::vector<double> x, const std::shared_ptr<void> inputs, std::vector<double> &parameter);
            
    }
  }
}

#endif
