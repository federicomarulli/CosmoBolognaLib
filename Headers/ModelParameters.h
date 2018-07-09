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
 *  @file Headers/ModelParameters.h
 *
 *  @brief The class ModelParameters
 *
 *  This file defines the interface of the class ModelParameters
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#ifndef __MPARAM__
#define __MPARAM__

#include "Prior.h"
#include "PriorDistribution.h"
#include "PosteriorDistribution.h"

namespace cbl {

  namespace statistics {

    /**
     * @enum ParameterType
     * @brief the parameter type
     */
    enum class ParameterType {

      /// base parameter
      _Base_,

      /// derived parameter
      _Derived_,

    };

   /**
     * @brief return a vector containing the
     * ParameterType names
     * @return a vector containing the
     * ParameterType names
     */
    inline std::vector<std::string> ParameterTypeNames () { return {"Base", "Derived"}; }

    /**
     *
     * @brief cast an enum of type ParameterType
     * from its index
     * @param parameterTypeIndex the parameterType index
     * @return object of class ParameterType
     */
    inline ParameterType ParameterTypeCast (const int parameterTypeIndex) { return castFromValue<ParameterType>(parameterTypeIndex); }

    /**
     * @brief cast an enum of type ParameterType
     * from its name
     * @param parameterTypeName the parameterType name
     * @return object of class ParameterType
     */
    inline ParameterType ParameterTypeCast (const std::string parameterTypeName) { return castFromName<ParameterType>(parameterTypeName, ParameterTypeNames()); }

    /**
     * @brief cast an enum of type ParameterType
     * from indeces
     * @param parameterTypeIndeces the parameterType indeces
     * @return object of class ParameterType
     */
    inline std::vector<ParameterType> ParameterTypeCast (const std::vector<int> parameterTypeIndeces) { return castFromValues<ParameterType>(parameterTypeIndeces); } 

    /**
     * @brief cast enums of type ParameterType
     * from thier names
     * @param parameterTypeNames the parameterType names
     * @return vector of ParameterType enums
     */
    inline std::vector<ParameterType> ParameterTypeCast (const std::vector<std::string> parameterTypeNames) { return castFromNames<ParameterType>(parameterTypeNames, ParameterTypeNames()); }


    /**
     * @class ModelParameters ModelParameters.h
     * "Headers/ModelParameters.h"
     *
     *  @brief The class ModelParameters
     *
     *  This class is used to define the model parameters
     */
    class ModelParameters {

      protected:

	/// model parameters
	std::vector<ParameterType> m_parameter_type;

	/// model parameter types
	std::vector<std::string> m_parameter_name;

	/// number of parameters
	size_t m_nparameters = 0;

	/// number of base parameters
	size_t m_nparameters_base = 0;

	/// number of derived parameters
	size_t m_nparameters_derived = 0;

	/// indexes of base parameters
	std::vector<unsigned int> m_base_parameters;

	/// indexes of the derived parameters
	std::vector<unsigned int> m_derived_parameters;

	/**
	 * @brief private member to set the parameter
	 * types
	 *
	 * @return none
	 */
	virtual void m_set_parameter_type();

      public:

	/**
	 *  @name Constructors/destructors
	 */
	///@{

	/**
	 *  @brief default constructor
	 *
	 *  @return object of class ModelParameters
	 */
	ModelParameters () = default;

	/**
	 *  @brief constructor for ModelParameters
	 *
	 *  @param nparameters the number of parameters
	 *
	 *  @param parameterTypes the parameter types
	 *
	 *  @param parameterNames the parameter names
	 *
	 *  @return object of class ModelParameters
	 */
	ModelParameters (const size_t nparameters, std::vector<ParameterType> parameterTypes, std::vector<std::string> parameterNames);

	/**
	 *  @brief default destructor
	 *
	 *  @return none
	 */
	~ModelParameters () = default;

	///@}
	
	/**
	 * @brief reset the parameter vectors
	 *
	 * @return none
	 */
	virtual void reset ()
	{set_parameters(m_nparameters, m_parameter_type, m_parameter_name);} 
	
	/**
	 *  @name Member functions used to set private/protected of the ModelParameters
	 */
	///@{

	/**
	 * @brief return the total number of
	 * parameters
	 *
	 * @return the total number of parameters
	 */
	size_t nparameters () const;

	/**
	 * @brief return the number of base
	 * parameters
	 *
	 * @return the number of base parameters
	 */
	size_t nparameters_base() const;

	/**
	 * @brief return the number of derived
	 * parameters
	 *
	 * @return the number of derived parameters
	 */
	size_t nparameters_derived() const;

	/**
	 * @brief return the model parameter type
	 *
	 * @param p the index of the parameter
	 * 
	 * @return vector containing the parameter type
	 */
	ParameterType type (const int p) const {return m_parameter_type[p];}

	/**
	 * @brief return all the model parameter names
	 * 
	 * @return vector containing all the parameter names
	 */
	std::vector<ParameterType> type () const {return m_parameter_type;}

	/**
	 * @brief return the model parameter name
	 * 
	 * @param p the index of the parameter
	 *
	 * @return vector containing the parameter name
	 */
	std::string name (const int p) const {return m_parameter_name[p];}

	/**
	 * @brief return all the model parameter names
	 * 
	 * @return vector containing all the parameter names
	 */
	std::vector<std::string> name () const {return m_parameter_name;}

	/**
	 * @brief return all the model parameters
	 * 
	 * @param parameter_values vector of parameters
	 *
	 * @return all the parameter values
	 */
	virtual std::vector<double> full_parameters (const std::vector<double> parameter_values) const;

	/**
	 *  @brief set the parameter
	 *
	 *  @param nparameters the number of parameters
	 *
	 *  @param parameterTypes the parameter types
	 *
	 *  @param parameterNames the parameter names*
	 *
	 *  @return none
	 */
	virtual void set_parameters (const size_t nparameters, std::vector<ParameterType> parameterTypes, std::vector<std::string> parameterNames);

	/**
	 *  @brief set the parameter
	 *
	 *  @param nparameters the number of parameters
	 *  
	 *  @param priorDistributions the priorDistributions
	 *
	 *  @param parameterTypes the parameter types
	 *
	 *  @param parameterNames the parameter names*
	 *
	 *  @return none
	 */
	virtual void set_parameters (const size_t nparameters, const std::vector<std::shared_ptr<PriorDistribution>> priorDistributions, std::vector<ParameterType> parameterTypes, std::vector<std::string> parameterNames)
	{(void)nparameters; (void)priorDistributions; (void)parameterTypes, (void)parameterNames; ErrorCBL("Error in set_parameters() of ModelParameters.h!");}

	///@}
	
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
	virtual size_t nparameters_free () const 
	{ErrorCBL("Error in nparameters_free() of ModelParameters.h!"); return 1;}

	/**
	 * @brief return the number of fixed
	 * parameters
	 *
	 * @return the number of fixed parameters
	 */
	virtual size_t nparameters_fixed() const
	{ErrorCBL("Error in nparameters_fixed() of ModelParameters.h!"); return 1;}

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
	virtual void free (const int p)
	{(void)p; ErrorCBL("Error in free() of ModelParameters.h!");}

	/**
	 * @brief fix the parameter at the input value;
	 * 
	 * @param p the p-th parameter
	 *
	 * @param value the input value
	 *
	 * @return none
	 */
	virtual void fix (const int p, const double value)
	{(void)p; (void)value; ErrorCBL("Error in fix() of ModelParameters.h!");}

	/**
	 * @brief fix the parameter at the bestfit value,
	 * contained in m_bestfit_value;
	 *
	 * @param p the p-th parameter
	 *
	 * @return none
	 */
	virtual void fix_at_bestfit (const int p)
	{(void)p; ErrorCBL("Error in fix_at_bestfit() of ModelParameters.h!");}

	///@}

	///@{
	
	/**
	 * @brief get the protected member m_value
	 *
	 *  @param p the p-th parameter
	 *
	 * @return the bestfit value of the parameter
	 */
	virtual double bestfit_value (const int p) const
	{(void)p; ErrorCBL("Error in bestfit_value() of ModelParameters.h!"); return 0.;}
	
	/**
	 * @brief get the protected member m_value
	 *
	 * @return the parameter bestfit values
	 */
	virtual std::vector<double> bestfit_values () const
	{ErrorCBL("Error in bestfit_values() of ModelParameters.h!"); std::vector<double> vv; return vv;}

	/**
	 *  @brief set the protected member m_bestfit_value
	 *
	 *  @param bestfit_value parameter bestfit values
	 *
	 *  @return none
	 */
	virtual void set_bestfit_value (const std::vector<double> bestfit_value)
	{(void)bestfit_value; ErrorCBL("Error in set_bestfit_values() of ModelParameters.h!");}

	/**
	 *  @brief write the bestfit info
	 *
	 *  @return none
	 */
	virtual void write_bestfit_info()
	{ErrorCBL("Error in write_bestfit_info() of ModelParameters.h!");}

	///@}

	/**
	 *  @name Member functions used to interact with prior distribution
	 */
	///@{
	
	/**
	 * @brief set the prior distribution for the p-th 
	 * parameter
	 *
	 * @param p the p-th parameter
	 *
	 * @param priorDistribution the prior distribution
	 *
	 * @return none
	 */
	virtual void set_prior_distribution (const int p, const std::shared_ptr<PriorDistribution> priorDistribution)
	{(void)p; (void)priorDistribution; ErrorCBL("Error in set_prior_distribution() of ModelParameters.h!");}

	/**
	 * @brief set the prior distributions for the 
	 * parameters
	 *
	 * @param priorDistributions the prior distributions
	 *
	 * @return none
	 */
	virtual void set_prior_distribution (const std::vector<std::shared_ptr<PriorDistribution>> priorDistributions)
	{(void)priorDistributions; ErrorCBL("Error in set_prior_distribution() of ModelParameters.h!");}

	/**
	 * @brief set the prior distributions for the 
	 * parameters
	 *
	 * @param ran_generator the random generator
	 *
	 * @return none
	 */
	virtual void set_prior_distribution_seed (const std::shared_ptr<random::UniformRandomNumbers_Int> ran_generator)
	{(void)ran_generator; ErrorCBL("Error in set_prior_distribution_seed() of ModelParameters.h!");} 
	
	/**
	 * @brief get the prior distribution for the p-th 
	 * parameter
	 *
	 * @param p the p-th parameter
	 *
	 * @return pointer to the prior distribution of the p-th
	 * parameter
	 */
	virtual std::shared_ptr<PriorDistribution> prior_distribution (const int p) const
	{(void)p; ErrorCBL("Error in set_prior_distribution() of ModelParameters.h!"); return NULL;}

	/**
	 * @brief get the prior distribution for the p-th 
	 * parameter
	 *
	 * @return pointer to the prior distribution of the p-th
	 * parameter
	 */
	virtual std::vector<std::shared_ptr<PriorDistribution>> prior_distribution () const
	{ErrorCBL("Error in set_prior_distribution() of ModelParameters.h!"); std::vector<std::shared_ptr<PriorDistribution>> vv; return vv;}

	/**
	 * @brief get the prior function
	 *
	 * @return pointer to a class Prior
	 */
	virtual std::shared_ptr<Prior> prior () const
	{ErrorCBL("Error in prior() of ModelParameters.h!"); return NULL;}

	///@}
			
	/**
	 * @brief set the internal method
	 *  m_parameter_covariance
	 *
	 * @param start the starting position 
	 *
	 * @param thin number of jumped indexes in the chain
	 *
	 * @return none
	 */
	virtual void set_parameter_covariance(const int start=0, const int thin=1)
	{ (void)start; (void)thin; ErrorCBL("Error in set_parameter_covariance() of ModelParameters.h!");}

	/**
	 * @brief return the protected member m_parameter_covariance
	 *
	 * @param i index i
	 *
	 * @param j index j
	 *
	 * @return the value of the parameter covariance at i,j
	 */
	virtual double parameter_covariance (const int i, const int j) const
	{ (void)i; (void)j; ErrorCBL("Error in parameter_covariance() of ModelParameters.h!"); return 0.;}

	/**
	 * @brief return the protected member m_parameter_covariance
	 *
	 * @return vector containing the parameter covariance
	 * matrix
	 */
	virtual std::vector<std::vector<double>> parameter_covariance () const
	{ErrorCBL("Error in parameter_covariance() of ModelParameters.h!"); std::vector<std::vector<double>> vv; return vv;}

	
	/**
	 *  @name Member functions used to interact with posterior distribution
	 */
	///@{
			
	/**
	 * @brief set the prior distribution for the p-th 
	 * parameter
	 *
	 * @param start the starting position 
	 *
	 * @param thin number of jumped indexes in the chain
	 *
	 * @param nbins the number of bins
	 *
	 * @param seed seed for random extraction
	 *
	 * @return none
	 */
	virtual void set_posterior_distribution (const int start, const int thin, const int nbins, const int seed=34121)
	{(void)start; (void)thin; (void)nbins; (void)seed; ErrorCBL("Error in set_posterior_distribution() of ModelParameters.h!"); }

	/**
	 *  @brief get the posterior distribution for the 
	 *  chosen parameter
	 *
	 *  @param par the index of the parameter
	 *
	 *  @return the protected member m_parameter_posterior[param]
	 */
	virtual std::shared_ptr<PosteriorDistribution> posterior_distribution (const int par) const 
	{(void)par; ErrorCBL("Error in posterior_distribution() of ModelParameters.h!"); return NULL;}

	///@}

	/**
	 * @brief show the results on the standard output
	 *
	 * @param start the starting position 
	 *
	 * @param thin number of jumped indexes in the chain
	 *
	 * @param nbins the number of bins
	 *
	 * @param seed seed for random extraction
	 *
	 * @return none
	 */
	virtual void show_results (const int start, const int thin, const int nbins, const int seed=34121)
	{(void)start; (void)thin; (void)nbins; (void)seed;  ErrorCBL("Error in show_results() of ModelParameters.h!"); }

	/**
	 * @brief write results on files
	 *
	 * @param dir name of the output folder
	 *
	 * @param file name of the output file
	 *
	 * @param start the starting position 
	 *
	 * @param thin number of jumped indexes in the chain
	 *
	 * @param nbins the number of bins
	 *
	 * @param seed seed for random extraction
	 *
	 * @return none
	 */
	virtual void write_results (const std::string dir, const std::string file, const int start, const int thin, const int nbins, const int seed=34121)
	{(void)dir; (void)file; (void)start; (void)thin; (void)nbins; (void)seed;  ErrorCBL("Error in write_results() of ModelParameters.h!"); }

	/**
	 * @brief return the private member m_chain_size
	 *
	 * @return the chain size
	 */
	virtual size_t chain_size () const
	{ErrorCBL("Error in size() of ModelParameters.h!"); return 0;}

	/**
	 * @brief return the private member m_chain_nwalkers
	 *
	 * @return the chain size
	 */
	virtual size_t chain_nwalkers () const
	{ErrorCBL("Error in nwalkers() of ModelParameters.h!"); return 0; }

	/**
	 * @brief set the chain
	 *
	 * @param size the chain lenght 
         *
	 * @param nwalkers the number of parallel walkers
	 *
	 * @return none
	 */
	virtual void set_chain (const size_t size, const size_t nwalkers)
	{(void)size; (void)nwalkers; ErrorCBL("Error in set_chain() of ModelParameters.h!"); }

	/**
	 * @brief reset the chain using m_size and m_nwalkers
	 *
	 * @return none
	 */
	virtual void reset_chain ()
	{ErrorCBL("Error in reset_chain() of ModelParameters.h!"); }

	/**
	 * @brief expand the already existing chain
	 *
	 * @param append the lenght of the empty chunk of the chain 
	 *
	 * @return none
	 */
	virtual void expand_chain (const int append)
	{(void)append; ErrorCBL("Error in expand_chain() of ModelParameters.h!"); }

	/**
	 * @brief return the private member m_chain_values at the pos step
	 * for the ww-th walker, for the chosen parameter
	 *
	 * @param param the parameter index
	 * 
	 * @param pos the position in the chain
	 * 
	 * @param ww the walker index
	 *
	 * @return the chain value
	 */
	virtual double chain_value (const int param, const int pos, const int ww) const
	{(void)param; (void)pos; (void)ww; ErrorCBL("Error in chain_value() of ModelParameters.h!"); return 0.; }

	/**
	 * @brief return the private member m_values at the pp-th step
	 * for the ww-th step for all the parameters
	 * 
	 * @param pos the position in the chain
	 * 
	 * @param ww the walker index
	 *
	 * @return the chain value
	 */
	virtual std::vector<double> chain_value_parameters (const int pos, const int ww) const
	{(void)pos; (void)ww; ErrorCBL("Error in chain_value_parameters() of ModelParameters.h!"); std::vector<double> vv; return vv;}

	/**
	 * @brief return all the chain values for a parameter
	 * 
	 * @param par the parameter index
	 *
	 * @param start the starting position 
	 *
	 * @param thin number of jumped indexes in the chain
	 *
	 * @return the chain value
	 */
	virtual std::vector<double> parameter_chain_values (const int par, const int start=0, const int thin = 1) const
	{(void)par; (void)start; (void)thin; ErrorCBL("Error in parameter_chain_values() of ModelParameters.h!");  std::vector<double> vv; return vv;}

	/**
	 * @brief set the private member m_chain_values at the pp-th step
	 * for the ww-th step
	 *
	 * @param param the parameter index
	 * 
	 * @param pos the position in the chain
	 * 
	 * @param ww the walker index
	 *
	 * @param value the chain value
	 *
	 * @return none
	 */
	virtual void set_chain_value (const int param, const int pos, const int ww, const double value) 
	{(void)param; (void)pos; (void)ww; (void)value; ErrorCBL("Error in set_chain_value() of ModelParameters.h!"); }

	/**
	 *  @brief set the chain values
	 *
	 *  @param values the input chain values
	 *
	 *  @param nwalkers the number of parallel walkers
	 *
	 *  @return none
	 */
	virtual void set_chain_values (const std::vector<std::vector<double>> values, const int nwalkers)
	{(void)values; (void)nwalkers; ErrorCBL("Error in set_chain_values() of ModelParameters.h!"); }

	/**
	 *  @brief set the chain values
	 *
	 *  @param values the input chain values
	 *
	 *  @return none
	 */
	virtual void set_chain_values (const std::vector<std::vector<std::vector<double>>> values)
	{(void)values;  ErrorCBL("Error in set_chain_values() of ModelParameters.h!"); }

	/**
	 * @brief initialize the chain values
	 *
	 * @param values the starting values
	 *
	 * @return none
	 */
	virtual void initialize_chain (const std::vector<std::vector<double>> values)
	{(void)values;  ErrorCBL("Error in initialize_chain() of ModelParameters.h!"); }

	/**
	 * @brief initialize the chain values
	 * random sampling the parameter priors
	 *
	 * @return none
	 */
	virtual void initialize_chain_from_prior ()
	{ErrorCBL("Error in initialize_chain_from_prior() of ModelParameters.h!"); }

	/**
	 * @brief initialize the chain values
	 *
	 * @param center the ball center
	 *
	 * @param radius the ball radius
	 *
	 * @param seed the random number generator seed
	 *
	 * @return none
	 */
	virtual void initialize_chain_ball (const std::vector<double> center, const double radius, const double seed)
	{(void)center; (void)radius; (void)seed;  ErrorCBL("Error in initialize_chain_ball() of ModelParameters.h!"); }

	/**
	 * @brief initialize the chain values around bestfit 
	 * values
	 *
	 * @param radius the ball radius
	 *
	 * @param seed the random number generator seed
	 *
	 * @return none
	 */
	virtual void initialize_chain_ball_bestfit (const double radius, const double seed)
	{(void)radius; (void)seed;  ErrorCBL("Error in initialize_chain_ball() of ModelParameters.h!"); }

    };
  }
}

#endif
