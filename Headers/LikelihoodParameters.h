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
 *  @file Headers/LikelihoodParameters.h
 *
 *  @brief The class LikelihoodParameters
 *
 *  This file defines the interface of the class LikelihoodParameters
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unibo.it, alfonso.veropalumbo@unibo.it
 */

#ifndef __LPARAM__
#define __LPARAM__

#include "ModelParameters.h"

namespace cbl {

  namespace statistics {

    /**
     * @class LikelihoodParameters LikelihoodParameters.h
     * "Headers/LikelihoodParameters.h"
     *
     *  @brief The class LikelihoodParameters
     *
     *  This class is used to define the model parameters
     */
    class LikelihoodParameters : public ModelParameters {

      protected:

	/// false \f$\rightarrow\f$ free parameter; true \f$\rightarrow\f$ fixed parameters
	std::vector<bool> m_parameter_isFixed;

	/// the number of free parameters
	size_t m_nparameters_free = 0;

	/// the number of fixed parameters
	size_t m_nparameters_fixed = 0;

	/// the indexes of fixed parameters
	std::vector<unsigned int> m_fixed_parameter;

	/// the indexes of the free parameters
	std::vector<unsigned int> m_free_parameter;

	/// the model parameter fixed values
	std::vector<double> m_parameter_fixed_value;

	/// the best-fit parameter values, i.e. the maxima of the likelihood 
	std::vector<double> m_parameter_bestfit_value;

	/**
	 * @brief private member to set the parameter
	 * types
	 *
	 * @return none
	 */
	void m_set_parameter_type () override;

	
      public:

	/**
	 *  @name Constructors/destructors
	 */
	///@{

	/**
	 *  @brief default constructor
	 *
	 *  @return object of class LikelihoodParameters
	 */
	LikelihoodParameters () = default;

	/**
	 *  @brief constructor for LikelihoodParameters
	 *
	 *  @param nparameters the number of parameters
	 *
	 *  @param parameterTypes the parameter types
	 *
	 *  @param parameterNames the parameter names
	 *
	 *  @return object of class LikelihoodParameters
	 */
	LikelihoodParameters (const size_t nparameters, std::vector<ParameterType> parameterTypes, std::vector<std::string> parameterNames);

	/**
	 *  @brief default destructor
	 *
	 *  @return none
	 */
	~LikelihoodParameters () = default;

	///@}
	
	/**
	 * @brief reset the parameter vectors
	 *
	 * @return none
	 */
	void reset () override;

	/**
	 * @brief return the model parameter status
	 * 
	 * @param p the index of the parameter
	 *
	 * @return the parameter status
	 */
	std::string status (const int p) const;

	/**
	 * @brief return all the model parameter status
	 * 
	 * @return vector containing all the parameter statuss
	 */
	std::vector<std::string> status () const;

	/**
	 *  @name Member functions used to set private/protected of the ModelParameters
	 */
	///@{

	/**
	 * @brief return the number of free
	 * parameters
	 *
	 * @return the number of free parameters
	 */
	size_t nparameters_free () const override;

	/**
	 * @brief return the private member m_free_parameter
	 *
	 * @return the private member m_free_parameter
	 */
	std::vector<unsigned int> free_parameter () const { return m_free_parameter; }

	/**
	 * @brief return the number of fixed
	 * parameters
	 *
	 * @return the number of fixed parameters
	 */
	size_t nparameters_fixed () const override;

	/**
	 * @brief return the private member m_fixed_parameter
	 *
	 * @return the private member m_fixed_parameter
	 */
	std::vector<unsigned int> fixed_parameter () const { return m_fixed_parameter; }

	/**
	 * @brief return all the model parameter
	 * 
	 * @param parameter_value vector of free parameter
	 *
	 * @return all the parameter values
	 */
	std::vector<double> full_parameter (const std::vector<double> parameter_value) const override;

	/**
	 *  @brief set the parameter
	 *
	 *  @param nparameters the number of parameters
	 *
	 *  @param parameterTypes the parameter types
	 *
	 *  @param parameterNames the parameter names
	 *
	 *  @return none
	 */
	void set_parameters (const size_t nparameters, std::vector<ParameterType> parameterTypes, std::vector<std::string> parameterNames) override;

	///@}
	
	/**
	 *  @name Member functions to manage fixed/free parameters
	 */
	///@{
		
	/**
	 * @brief set m_fixed to false;
	 * 
	 * @param p the p-th parameter
	 *
	 * @return none
	 */
	void free (const int p) override;

	/**
	 * @brief fix the parameter at the input value;
	 * 
	 * @param p the p-th parameter
	 *
	 * @param value the input value
	 *
	 * @return none
	 */
	void fix (const int p, const double value) override;

	/**
	 * @brief fix the parameter at the bestfit value,
	 * contained in m_bestfit_value;
	 *
	 * @param p the p-th parameter
	 *
	 * @return none
	 */
	void fix_at_bestfit (const int p) override;

	///@}
	

	///@{
	
	/**
	 * @brief get the protected member m_value
	 *
	 *  @param p the p-th parameter
	 *
	 * @return the bestfit value of the parameter
	 */
	double bestfit_value (const int p) const override;
	
	/**
	 * @brief get the protected member m_value
	 *
	 * @return the parameter bestfit values
	 */
	std::vector<double> bestfit_value () const override;

	/**
	 *  @brief set the protected member m_bestfit_value
	 *
	 *  @param bestfit_value parameter bestfit values
	 *
	 *  @return none
	 */
	void set_bestfit_values (const std::vector<double> bestfit_value) override;

	/**
	 *  @brief write the best fit info
	 *
	 *  @return none
	 */
	void write_bestfit_info() override;

	///@}
    };
  }
}

#endif
