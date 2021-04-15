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
 *  @authors federico.marulli3@unibo.it, alfonso.veropalumbo@unibo.it
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

       /// correlated parameters
        _Correlated_,
    };

    /**
     * @brief return a vector containing the ParameterType names
     *
     * @return a vector containing the ParameterType names
     */
    inline std::vector<std::string> ParameterTypeNames () { return {"Base", "Derived", "Correlated"}; }

    /**
     * @brief cast an enum of type ParameterType from its index
     *
     * @param parameterTypeIndex the parameterType index
     * @return object of class ParameterType
     */
    inline ParameterType ParameterTypeCast (const int parameterTypeIndex) { return castFromValue<ParameterType>(parameterTypeIndex); }

    /**
     * @brief cast an enum of type ParameterType from its name
     *
     * @param parameterTypeName the parameterType name
     * @return object of class ParameterType
     */
    inline ParameterType ParameterTypeCast (const std::string parameterTypeName) { return castFromName<ParameterType>(parameterTypeName, ParameterTypeNames()); }

    /**
     * @brief cast an enum of type ParameterType from indeces
     *
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
     *  @class ModelParameters ModelParameters.h
     *  "Headers/ModelParameters.h"
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

      /// number of correlated parameters
      size_t m_nparameters_correlated = 0;

      /// indexes of base parameters
      std::vector<unsigned int> m_base_parameter;

      /// indexes of the derived parameters
      std::vector<unsigned int> m_derived_parameter;


      /**
       * @brief private member to set the parameter
       * types
       *
       * @return none
       */
      virtual void m_set_parameter_type ();
	
	
    public:

      /**
       *  @name Constructors/destructors
       */
      ///@{

      /**
       *  @brief default constructor
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
       */
      ModelParameters (const size_t nparameters, std::vector<ParameterType> parameterTypes, std::vector<std::string> parameterNames);

      /**
       *  @brief constructor for ModelParameters
       *
       *  @param nparameters the number of parameters
       * 
       *  @param parameterNames the parameter names
       */
      ModelParameters (const size_t nparameters, std::vector<std::string> parameterNames);


      /**
       *  @brief default destructor
       */
      ~ModelParameters () = default;

      ///@}
	
      /**
       * @brief reset the parameter vectors
       *
       * @return none
       */
      virtual void reset ()
      { set_parameters(m_nparameters, m_parameter_type, m_parameter_name); } 

	
      /**
       *  @name Member functions used to set private/protected of the ModelParameters
       */
      ///@{
      
      /**
       * @brief return the total number of parameters
       *
       * @return the total number of parameters
       */
      size_t nparameters () const;

      /**
       * @brief return the number of base parameters
       *
       * @return the number of base parameters
       */
      size_t nparameters_base () const;

      /**
       * @brief return the private member m_base_parameters
       *
       * @return the private member m_base_parameters
       */
      std::vector<unsigned int> base_parameter () const { return m_base_parameter; }

      /**
       * @brief return the number of derived parameters
       *
       * @return the number of derived parameters
       */
      size_t nparameters_derived () const;

      /**
       * @brief return the private member m_derived_parameter
       *
       * @return the private member m_derived_parameter
       */
      std::vector<unsigned int> derived_parameter () const { return m_derived_parameter; }

      /**
       * @brief return the model parameter type
       *
       * @param p the index of the parameter
       * 
       * @return vector containing the parameter type
       */
      ParameterType type (const int p) const { return m_parameter_type[p]; }

      /**
       * @brief return all the model parameter names
       * 
       * @return vector containing all the parameter names
       */
      std::vector<ParameterType> type () const { return m_parameter_type; }

      /**
       * @brief return the model parameter name
       * 
       * @param p the index of the parameter
       *
       * @return the parameter name
       */
      std::string name (const int p) const { return m_parameter_name[p]; }

      /**
       * @brief return all the model parameter names
       * 
       * @return vector containing all the parameter names
       */
      std::vector<std::string> name () const { return m_parameter_name; }

      /**
       * @brief return the model parameter status
       * 
       * @param p the index of the parameter
       *
       * @return the parameter status
       */
      virtual std::string status (const int p) const
      { (void)p; ErrorCBL("", "status", "ModelParameters.h"); return par::defaultString;}

      /**
       * @brief return all the model parameter status
       * 
       * @return vector containing all the parameter statuss
       */
      virtual std::vector<std::string> status () const
      { ErrorCBL("", "status", "ModelParameters.h"); std::vector<std::string> vv; return vv;}

      /**
       * @brief return all the model parameters, for internal usage
       * 
       * @param parameter_value vector of parameters
       *
       * @return all the parameter values
       */
      virtual std::vector<double> full_parameter (const std::vector<double> parameter_value) const;

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
       *  @param parameterNames the parameter names
       *
       *  @return none
       */
      virtual void set_parameters (const size_t nparameters, std::vector<std::string> parameterNames);

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
      { (void)nparameters; (void)priorDistributions; (void)parameterTypes, (void)parameterNames; ErrorCBL("", "set_parameters", "ModelParameters.h"); }

      ///@}

	
      /**
       *  @name Member functions used to set private/protected of the ModelParameters
       */
      ///@{
      
      /**
       * @brief return the number of free parameters
       *
       * @return the number of free parameters
       */
      virtual size_t nparameters_free () const 
      { ErrorCBL("", "nparameters_free", "ModelParameters.h"); return 1; }

      /**
       * @brief return the private member m_free_parameter
       *
       * @return the private member m_free_parameter
       */
      virtual std::vector<unsigned int> free_parameter () const
      { ErrorCBL("", "free_parameter", "ModelParameters.h"); std::vector<unsigned int> vv; return vv; }

      /**
       * @brief return the number of fixed parameters
       *
       * @return the number of fixed parameters
       */
      virtual size_t nparameters_fixed () const
      { ErrorCBL("", "nparameters_fixed", "ModelParameters.h"); return 1; }

      /**
       * @brief return the private member m_fixed_parameters
       *
       * @return the private member m_fixed_parameters
       */
      virtual std::vector<unsigned int> fixed_parameter () const
      { ErrorCBL("", "fixed_parameters", "ModelParameters.h"); std::vector<unsigned int> vv; return vv; }

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
      { (void)p; ErrorCBL("", "free", "ModelParameters.h"); }

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
      { (void)p; (void)value; ErrorCBL("", "fix", "ModelParameters.h"); }

      /**
       * @brief fix the parameter at the bestfit value,
       * contained in m_bestfit_value
       *
       * @param p the p-th parameter
       *
       * @return none
       */
      virtual void fix_at_bestfit (const int p)
      { (void)p; ErrorCBL("", "fix_at_bestfit", "ModelParameters.h"); }

      ///@}

      ///@{
	
      /**
       * @brief get the protected member m_value
       *
       * @param p the p-th parameter
       *
       * @return the bestfit value of the parameter
       */
      virtual double bestfit_value (const int p) const
      { (void)p; ErrorCBL("", "bestfit_value", "ModelParameters.h"); return 0.; }
	
      /**
       * @brief get the protected member m_value
       *
       * @return the parameter bestfit values
       */
      virtual std::vector<double> bestfit_value () const
      { ErrorCBL("", "bestfit_value", "ModelParameters.h"); std::vector<double> vv; return vv; }

      /**
       *  @brief set the protected member m_bestfit_value
       *
       *  @param bestfit_value parameter bestfit values
       *
       *  @return none
       */
      virtual void set_bestfit_values (const std::vector<double> bestfit_value)
      { (void)bestfit_value; ErrorCBL("", "set_bestfit_values", "ModelParameters.h"); }
	
      /**
       *  @brief set the protected member m_bestfit_value
       *
       *  @param start the starting position 
       *
       *  @param thin number of jumped indexes in the chain
       *
       *  @param nbins the number of bins
       *
       *  @param seed seed for random extraction
       *
       *  @return none
       */
      virtual void set_bestfit_values (const int start, const int thin, const int nbins, const int seed) 
      { (void)start; (void)thin; (void)nbins; (void)seed; ErrorCBL("", "set_bestfit_values", "ModelParameters.h"); }
	
      /**
       *  @brief write the bestfit info
       *
       *  @return none
       */
      virtual void write_bestfit_info ()
      { ErrorCBL("", "write_bestfit_info", "ModelParameters.h"); }

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
      { (void)p; (void)priorDistribution; ErrorCBL("", "set_prior_distribution", "ModelParameters.h"); }

      /**
       * @brief set the prior distributions for the 
       * parameters
       *
       * @param priorDistributions the prior distributions
       *
       * @return none
       */
      virtual void set_prior_distribution (const std::vector<std::shared_ptr<PriorDistribution>> priorDistributions)
      { (void)priorDistributions; ErrorCBL("", "set_prior_distribution", "ModelParameters.h"); }

      /**
       * @brief set the prior distributions for the 
       * parameters
       *
       * @param ran_generator the random generator
       *
       * @return none
       */
      virtual void set_prior_distribution_seed (const std::shared_ptr<random::UniformRandomNumbers_Int> ran_generator)
      { (void)ran_generator; ErrorCBL("", "set_prior_distribution_seed", "ModelParameters.h"); } 
	
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
      { (void)p; ErrorCBL("", "prior_distribution", "ModelParameters.h"); return NULL; }

      /**
       * @brief get the prior distribution for the p-th parameter
       *
       * @return pointer to the prior distribution of the p-th
       * parameter
       */
      virtual std::vector<std::shared_ptr<PriorDistribution>> prior_distribution () const
      { ErrorCBL("", "prior_distribution", "ModelParameters.h"); std::vector<std::shared_ptr<PriorDistribution>> vv; return vv; }

      /**
       * @brief get the prior function
       *
       * @return pointer to a class Prior
       */
      virtual std::shared_ptr<Prior> prior () const
      { ErrorCBL("", "prior", "ModelParameters.h"); return NULL; }

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
      virtual void set_parameter_covariance (const int start=0, const int thin=1)
      { (void)start; (void)thin; ErrorCBL("", "set_parameter_covariance", "ModelParameters.h"); }

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
      { (void)i; (void)j; ErrorCBL("", "parameter_covariance", "ModelParameters.h"); return 0.; }

      /**
       * @brief return the protected member m_parameter_covariance
       *
       * @return vector containing the parameter covariance
       * matrix
       */
      virtual std::vector<std::vector<double>> parameter_covariance () const
      { ErrorCBL("", "parameter_covariance", "ModelParameters.h"); std::vector<std::vector<double>> vv; return vv; }

	
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
       * @param weight chain weight
       *
       * @return none
       */
      virtual void set_posterior_distribution (const int start, const int thin, const int nbins, const int seed=34121, const std::vector<double> weight={})
      { (void)start; (void)thin; (void)nbins; (void)seed; (void)weight; ErrorCBL("", "set_posterior_distribution", "ModelParameters.h"); }

      /**
       *  @brief get the posterior distribution for the 
       *  chosen parameter
       *
       *  @param par the index of the parameter
       *
       *  @return the protected member m_parameter_posterior[param]
       */
      virtual std::shared_ptr<PosteriorDistribution> posterior_distribution (const int par) const 
      { (void)par; ErrorCBL("", "posterior_distribution", "ModelParameters.h"); return NULL; }

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
       * @param show_mode true \f$\rightarrow\f$ show the posterior
       * mode; false \f$\rightarrow\f$ do not show the posterior
       * mode
       *
       * @param ns number of samples used to estimate the covariance
       * matrix
       *
       * @param nb number of data measurements, e.g. the bins of the
       * dataset
       *
       * @param weight chain weight
       *
       * @return none
       */
      virtual void show_results (const int start, const int thin, const int nbins, const int seed=34121, const bool show_mode=false, const int ns=-1, const int nb=-1, const std::vector<double> weight={})
      { (void)start; (void)thin; (void)nbins; (void)seed; (void)show_mode; (void)ns; (void)nb; (void)weight; ErrorCBL("", "show_results", "ModelParameters.h"); }

      /**
       * @brief store the results to file
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
       * @param compute_mode true \f$\rightarrow\f$ compute the
       * posterior mode; false \f$\rightarrow\f$ do not compute the
       * posterior mode
       *
       * @param ns number of samples used to estimate the covariance
       * matrix
       *
       * @param nb number of data measurements, e.g. the bins of the
       * dataset
       *
       * @param weight chain weight
       *
       * @return none
       */
      virtual void write_results (const std::string dir, const std::string file, const int start, const int thin, const int nbins, const int seed=34121, const bool compute_mode=false, const int ns=-1, const int nb=-1, const std::vector<double> weight={})
      { (void)dir; (void)file; (void)start; (void)thin; (void)nbins; (void)seed; (void)compute_mode; (void)ns; (void)nb; (void)weight; ErrorCBL("", "write_results", "ModelParameters.h"); }

      /**
       * @brief return the private member m_chain_size
       *
       * @return the chain size
       */
      virtual size_t chain_size () const
      { ErrorCBL("", "chain_size", "ModelParameters.h"); return 0; }

      /**
       * @brief return the private member m_chain_nwalkers
       *
       * @return the chain size
       */
      virtual size_t chain_nwalkers () const
      { ErrorCBL("", "chain_nwalkers", "ModelParameters.h"); return 0; }

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
      { (void)size; (void)nwalkers; ErrorCBL("", "set_chain", "ModelParameters.h"); }

      /**
       * @brief reset the chain using m_size and m_nwalkers
       *
       * @return none
       */
      virtual void reset_chain ()
      { ErrorCBL("", "reset_chain", "ModelParameters.h"); }

      /**
       * @brief expand the already existing chain
       *
       * @param append the lenght of the empty chunk of the chain 
       *
       * @return none
       */
      virtual void expand_chain (const int append)
      { (void)append; ErrorCBL("", "expand_chain", "ModelParameters.h"); }

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
      { (void)param; (void)pos; (void)ww; ErrorCBL("", "chain_value", "ModelParameters.h"); return 0.; }

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
      { (void)pos; (void)ww; ErrorCBL("", "chain_value_parameters", "ModelParameters.h"); std::vector<double> vv; return vv; }

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
      { (void)par; (void)start; (void)thin; ErrorCBL("", "parameter_chain_values", "ModelParameters.h");  std::vector<double> vv; return vv; }

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
      { (void)param; (void)pos; (void)ww; (void)value; ErrorCBL("", "set_chain_value", "ModelParameters.h"); }

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
      { (void)values; (void)nwalkers; ErrorCBL("", "set_chain_values", "ModelParameters.h"); }

      /**
       *  @brief set the chain values
       *
       *  @param values the input chain values
       *
       *  @return none
       */
      virtual void set_chain_values (const std::vector<std::vector<std::vector<double>>> values)
      { (void)values;  ErrorCBL("", "set_chain_values", "ModelParameters.h"); }

      /**
       * @brief initialize the chain values
       *
       * @param values the starting values
       *
       * @return none
       */
      virtual void initialize_chain (const std::vector<std::vector<double>> values)
      { (void)values;  ErrorCBL("", "initialize_chain", "ModelParameters.h"); }

      /**
       * @brief initialize the chain values
       * random sampling the parameter priors
       *
       * @return none
       */
      virtual void initialize_chain_from_prior ()
      { ErrorCBL("", "initialize_chain_from_prior", "ModelParameters.h"); }

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
      { (void)center; (void)radius; (void)seed; ErrorCBL("", "initialize_chain_ball", "ModelParameters.h"); }

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
      { (void)radius; (void)seed;  ErrorCBL("", "initialize_chain_ball_bestfit", "ModelParameters.h"); }

    };
  }
}

#endif
