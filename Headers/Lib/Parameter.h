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
 *  @file Headers/Lib/Parameter.h
 *
 *  @brief The class Parameter
 *
 *  This file defines the interface of the class Parameter
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#ifndef __PARAM__
#define __PARAM__

#include "Chain.h"
#include "Prior.h"
#include "Posterior.h"


// ===================================================================================================


namespace cosmobl {

  namespace statistics {

    /**
     * @enum ParameterType
     * @brief the parameter type
     */
    enum ParameterType {

      /// base parameter
      _BaseParameter_,

      /// derived parameter
      _DerivedParameter_,
      
    };


    /**
     *  @class Parameter Parameter.h "Headers/Lib/Parameter.h"
     *
     *  @brief The class Parameter
     *
     *  This class is used to define the parameters of models
     */
    class Parameter {

    protected:

      /// the parameter name
      string m_name;

      /// the parameter type: it can be either _BaseParameter_ or _DerivedParameter_
      ParameterType m_pType;

      /// 1 &rarr; bestfit is set, 0 &rarr; bestfit is not set
      bool m_bestfit_set = false;

      /// the best-fit parameter value
      double m_bestfit_value;

      /// the parameter chain
      shared_ptr<Chain> m_chain;

      /// the parameter posterior
      shared_ptr<statistics::Posterior> m_posterior;

      /**
       * @brief set the prior 
       * 
       * @return none
       */
      virtual void m_set_prior ()
      { ErrorCBL("Error in m_set_prior of Parameter.h"); }

    public:

      /**
       *  @name Constructors/destructors
       */
      ///@{

      /**
       *  @brief default constructor
       *
       *  @return object of class Parameter
       */
      Parameter () = default;

      /**
       *  @brief base constructor
       *
       *  @param pType the parameter type: it can be either
       *  _BaseParameter_ or _DerivedParameter_
       *
       *  @param name parameter name
       *
       *  @return object of class Parameter
       */
      Parameter (const ParameterType pType, const string name) : m_name(name), m_pType(pType) {}

      /**
       *  @brief default destructor
       *
       *  @return none
       */
      virtual ~Parameter () = default;

      ///@}
      

      /**
       *  @name Member functions used to set private/protected  
       *  of the the Parameter
       */
      ///@{

      /**
       * @brief get the protected member m_name
       *
       * @return the name of the parameter
       */
      string name () const;

      /**
       * @brief get the protected member m_pType
       *
       * @return the type of parameter
       */
      ParameterType parameterType () const;

      /**
       * @brief get the protected member m_bestfit_value
       *
       * @return the best fit value of the parameter
       */
      double bestfit_value () const;

      /**
       * @brief get the defined posterior
       * 
       * @return a shared pointer to an object of class statistics::Posterior
       */
      shared_ptr<statistics::Posterior> posterior () const;

      /**
       *  @brief set the protected member m_name
       *
       *  @param name parameter name
       *
       *  @return none
       */
      void set_name (const string name);

      /**
       *  @brief set the protected member m_pType
       *
       *  @param pType the parameter type
       *
       *  @return none
       */
      void set_pType (const ParameterType pType);

      /**
       *  @brief set the protected member m_bestfit_value
       *
       *  @param bestfit_value parameter value
       *
       *  @return none
       */
      void set_bestfit_value (const double bestfit_value);

      ///@}
	

      /**
       *  @name Member functions used to interact with Chain
       */
      ///@{
	
      /**
       * @brief return the chain
       *
       * @return chain pointer to the a Chain object
       */
      shared_ptr<Chain> chain() const;

      /**
       * @brief return the chain size
       *
       * @return the chain size
       */
      int chain_size() const;

      /**
       * @brief return the number
       * of parallel walkers
       *
       * @return the number of parallel walkers
       */
      int nwalkers() const;
	
      /**
       * @brief return the chain value
       *
       * @param pp the position in the chain
       *
       * @param ww the walker
       *
       * @return chain value
       */
      double chain_value (const int pp, const int ww) const;
		
      /**
       * @brief return the chain values
       *
       * @param start the starting position 
       *
       * @param thin number of jumped indexes in the chain
       *
       * @return vector containing the flattened chain values
       */
      vector<double> chain_values (const int start=0, const int thin=1) const;
	
      /**
       * @brief set the chain
       *
       * @param chain pointer to the a Chain object
       *
       * @return none
       */
      void set_chain (shared_ptr<Chain> chain);

      /**
       * @brief set the chain
       *
       * @param chain_size the chain size
       *
       * @param nwalkers the number of parallel walkers
       *
       * @return none
       */
      void set_chain (const int chain_size, const int nwalkers);

      /**
       * @brief set the chain value
       *
       * @param pp index of position in the chain
       *
       * @param ww index of the walker
       *
       * @param chain_value the chain value
       *
       * @return none
       */
      void set_chain_value (const int pp, const int ww, const double chain_value);

      /**
       * @brief set the chain
       *
       * @param chain_values the chain values
       *
       * @param nwalkers the number of parallel walkers
       *
       * @return none
       */
      void set_chain (const vector<double> chain_values, const int nwalkers);

      /**
       * @brief set the chain
       *
       * @param chain_values the chain values
       *
       * @return none
       */
      void set_chain (const vector<vector<double>> chain_values);
	
      /**
       * @brief expand the chain
       *
       * @param append the empty chunk size
       *
       * @return none
       */
      void expand_chain (const int append);

      /**
       * @brief set the posterior distribution
       * from the chain
       *
       * @param start the chain values
       *
       * @param thin the number of parallel walkers
       *
       * @param seed the distribution seed
       *
       * @return none
       */
      virtual void set_posterior (const int start, const int thin, const int seed)
      { (void)start; (void)thin; (void)seed; ErrorCBL("Error in set_posterior of Parameter.h"); }

      /**
       * @brief set the posterior distribution
       * from the chain
       *
       * @param chain_values the chain values
       *
       * @param nwalkers the number of parallel walkers
       *
       * @param seed the distribution seed
       *
       * @return none
       */
      virtual void set_posterior (const vector<double> chain_values, const int nwalkers, const int seed)
      { (void)chain_values; (void)nwalkers; (void)seed; ErrorCBL("Error in set_posterior of Parameter.h"); }

      /**
       * @brief set the posterior distribution
       * from the chain
       *
       * @param chain_values the chain values
       *
       * @param seed the distribution seed
       *
       * @return none
       */
      virtual void set_posterior (const vector<vector<double>> chain_values, const int seed)
      { (void)chain_values; (void)seed; ErrorCBL("Error in set_posterior of Parameter.h"); }

      ///@}

      /**
       *  @name Member functions used to interact with
       *  posterior distribution
       */
      ///@{

      /**
       * @brief value of the posterior at a proposed value
       *
       * @param value proporsed parameter value
       *
       * @return posterior value
       */ 
      double PosteriorProbability (const double value) const;

      /**
       * @brief get the posterior distribution mean
       *
       * @return the mean of the parameter
       * posterior distribution
       */
      double posterior_mean () const;

      /**
       * @brief get the posterior distribution median
       *
       * @return the median value of the parameter
       * posterior
       */   
      double posterior_median () const;

      /**
       * @brief get the posterior distribution standard deviation
       *
       * @return the standard deviation value of the 
       * parameter posterior
       */
      double posterior_std () const;

      /**
       * @brief get the posterior percentile
       *
       * @param i the i-th percentile
       *
       * @return the posterior i-th percentile
       */
      double posterior_percentile (const unsigned int i) const;

      /**
       * @brief get the posterior mode
       *
       * @return the posterior mode
       */
      double posterior_mode () const;

      /**
       *  @brief extract a parameter value from the posterior
       *
       *  @return value random-sampled from posterior distribution
       */
      double posterior_sample ();

      /**
       *  @brief extract values from the posterior 	 
       *  
       *  @param sample_size the size  of the extracted sample
       *
       *  @return a vector of values random-sampled from 
       *  posterior distribution
       */
      vector<double> posterior_sample (const int sample_size);

      ///@}

      /**
       *  @name Member functions used to interact with
       *  prior distribution
       */
      ///@{
	
      /**
       * @brief return private member m_fixed;
       *
       * @return false &rarr; the parameter is free; true &rarr; the
       * parameter is fixed
       */
      virtual bool fixed () const
      { ErrorCBL("Error in fixed of Parameter.h"); return 0; }

      /**
       * @brief set the parameter as fixed
       *
       * @param fix 1 &rarr; the parameter is fixed; 0
       * &rarr; the parameter is free
       *
       * @return none
       */
      virtual void set_fixed (const bool fix)
      { (void)fix; ErrorCBL("Error in set_fixed of Parameter.h"); }

      /**
       * @brief free the parameter
       *
       * @return none
       */
      virtual void free ()
      { ErrorCBL("Error in free of Parameter.h"); }

      /**
       * @brief fix the parameter 
       * at the m_value;
       *
       * @return none
       */
      virtual void fix ()
      { ErrorCBL("Error in fix of Parameter.h"); }

      /**
       * @brief fix the parameter at the 
       * input value;
       * 
       * @param value the input value
       *
       * @return none
       */
      virtual void fix (const double value)
      { (void)value; ErrorCBL("Error in fix of Parameter.h"); }

      /**
       * @brief fix the parameter at the bestfit value,
       * contained in m_bestfit_value;
       *
       * @return none
       */
      virtual void fix_at_bestfit()
      { ErrorCBL("Error in fix_at_bestfit of Parameter.h"); }

      /**
       * @brief get the protected member m_value
       *
       * @return the value of the parameter
       */
      virtual double value () const
      { ErrorCBL("Error in bestfit_value of Parameter.h"); return 0; }

      /**
       *  @brief set the protected member m_bestfit_value
       *
       *  @param value parameter value
       *
       *  @return none
       */
      virtual void set_value (const double value)
      { (void)value; ErrorCBL("Error in value of Parameter.h"); }

      /**
       * @brief get the defined prior
       * 
       * @return a shared pointer to an object of class glob::Distribution
       */
      virtual shared_ptr<statistics::Prior> prior () const
      {; ErrorCBL("Error in prior of Parameter.h");  return NULL; }
	
      /**
       *  @brief set user defined prior
       *
       *  @param prior user defined prior
       *
       *  @return none
       */
      virtual void set_prior (const statistics::Prior prior)
      { (void)prior; ErrorCBL("Error in set_prior of Parameter.h"); }

      /**
       *  @brief set an uniform prior
       *
       *  @param xmin lower limit of the distribution
       *
       *  @param xmax upper limit of the distribution
       *
       *  @param seed the distribution seed for random sampling
       *
       *  @return none
       */
      virtual void set_prior (const double xmin, const double xmax, const int seed)
      { (void)xmin; (void)xmax; (void)seed; ErrorCBL("Error in set_prior of Parameter.h"); }

      /**
       *  @brief set Poisson or Gaussiand distribution prior: 
       *  see the Documentation for the class Distribution
       *
       *  @param priorType the type of distribution to be created
       *
       *  @param prior_params parameters of the distribution function or discrete
       *  list of values for discrete distribution
       *
       *  @param xmin lower limit of the distribution
       *
       *  @param xmax upper limit of the distribution
       *
       *  @param seed the distribution seed for random sampling
       *
       *  @return none
       */
      virtual void set_prior (const glob::DistributionType priorType, const vector<double> prior_params, const double xmin, const double xmax, const int seed)
      { (void)priorType; (void)prior_params; (void)xmin; (void)xmax; (void)seed;  ErrorCBL("Error in set_prior of Parameter.h"); }

      /**
       * @brief value of the prior at a proposed value
       *
       * @param value proporsed parameter value
       *
       * @return prior value
       */ 
      virtual double PriorProbability (const double value) const
      { (void)value; ErrorCBL("Error in PriorProbability of Parameter.h"); return 0; }

      /**
       * @brief value of the prior at a proposed value
       *
       * @param value proporsed parameter value
       *
       * @return prior value
       */ 
      virtual double LogPriorProbability (const double value) const
      { (void)value; ErrorCBL("Error in LogPriorProbability of Parameter.h"); return 0; }

      /**
       *  @brief set prior seed
       *
       *  @param seed the prior seed
       *
       *  @return none
       */
      virtual void set_prior_seed (const int seed)
      { (void)seed; ErrorCBL("Error in set_prior_seed of Parameter.h"); }

      /**
       *  @brief extract a parameter value from the prior
       *
       *  @return value extracted from the prior distribution
       */
      virtual double prior_sample ()
      { ErrorCBL("Error in prior_sample of Parameter.h"); return 0; }

      /**
       *  @brief extract values from the prior 	 
       *
       *  @param sample_size the size of the extracted sample
       *
       *  @return a vector of values extracted from the 
       *  prior distribution
       */
      virtual vector<double> prior_sample (const int sample_size)
      { (void)sample_size; ErrorCBL("Error in prior_sample of Parameter.h"); vector<double> vv; return vv; }

      /**
       *  @brief return the size of 
       *  the prior range times epsilon 	 
       *
       *  @param epsilon the relative fraction of the interval size
       *
       *  @return the prior range
       */
      virtual double prior_range (const double epsilon=1.)
      { (void)epsilon; ErrorCBL("Error in prior_range of Parameter.h"); return 0; }

      ///@}

    };
  }
}

#endif
