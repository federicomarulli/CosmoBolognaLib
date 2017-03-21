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

#include "Prior.h"

namespace cosmobl {

  namespace statistics {

    /**
     * @enum ParameterType
     * @brief the parameter type
     */
    enum ParameterType {
      
      /// fixed parameter
      _fixed_,
      
      /// free parameter
      _free_
      
    };

  
    /**
     * @class Parameter Parameter.h "Headers/Lib/Parameter.h"
     *
     *  @brief The class Parameter
     *
     *  This class is used to define the parameters of models
     */
    class Parameter {
      
    protected:

      /// the parameter value used for internal computations
      double m_value;

      /// the best-fit parameter value, mode of the posterior
      double m_best_value;

      /// the parameter name
      string m_name;

      /// the parameter type: it can be either _free_ or _fixed_
      ParameterType m_pType;

      /// the parameter prior
      shared_ptr<Prior> m_prior;

      /// the number of chains 
      int m_nchains;
	
      /// the chain size
      int m_chain_size;
	
      /// the vectors containing the parameter chains
      vector< shared_ptr<Chain> > m_chains;

      /// the mean value of the posterior
      double m_mean;
	
      /// the parameter standard deviation 
      double m_std;
	
      /// the median value of the posterior
      double m_median;

      /// the binned parameter range
      vector<double> m_var;
	
      /// the parameter distribution
      vector<double> m_dist;

      
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
       *  @brief constructor for parameter with uniform non limited
       *  prior
       *
       *  @param value parameter value
       *  @param pType the parameter type: it can be either _free_ or
       *  _fixed_
       *  @param name parameter name
       *
       *  @return object of class Parameter
       */
      Parameter (const double value, const ParameterType pType, const string name);

      /**
       *  @brief constructor for parameter with uniform prior, with
       *  limits
       *
       *  @param value parameter value
       *  @param xmin lower value for the parameter
       *  @param xmax upper value for the parameter
       *  @param pType the parameter type: it can be either _free_ or
       *  _fixed_
       *  @param name parameter name
       *
       *  @return object of class Parameter
       */
      Parameter (const double value, const double xmin, const double xmax, const ParameterType pType=_free_, const string name="parameter");

      /**
       *  @brief constructor for parameter with a Gaussian/Poisson
       *  prior, with limits
       *
       *  @param value parameter value
       *  @param priorType the type of prior to be created
       *  @param prior_params parameters of the prior function	
       *  value for the parameter
       *  @param xmin lower value for the parameter
       *  @param xmax upper value for the parameter
       *  @param pType the parameter type: it can be either _free_ or
       *  _fixed_
       *  @param name parameter name
       *
       *  @return object of class Parameter
       */
      Parameter (const double value, const PriorType priorType, const vector<double> prior_params, const double xmin, const double xmax, const ParameterType pType=_free_, const string name="parameter");

      /**
       *  @brief constructor for parameter with discrete distribution
       *  prior
       *
       *  @param value parameter value
       *  @param priorType the type of prior to be created
       *  @param discrete_values discrete values for the parameter
       *  @param weights weights for discrete values
       *  @param pType the parameter type: it can be either _free_ or
       *  _fixed_
       *  @param name parameter name
       *
       *  @return object of class Parameter
       */
      Parameter (const double value, const PriorType priorType, const vector<double> discrete_values, const vector<double> weights, const ParameterType pType=_free_, const string name="parameter");

      /**
       *  @brief default constructor
       * 
       *  @param value parameter value
       *  @param prior prior for the parameter 
       *  @param pType the parameter type: it can be either _free_ or
       *  _fixed_
       *  @param name parameter name
       *
       *  @return object of class Parameter
       */
      Parameter (const double value, const Prior prior, const ParameterType pType=_free_, string name="parameter");

      /**
       *  @brief default constructor for fixed parameters
       * 
       *  @param value parameter value
       *
       *  @return object of class Parameter
       */
      Parameter (const double value) : Parameter(value, _fixed_, "parameter") {}

      /**
       *  @brief default destructor
       *
       *  @return none
       */
      ~Parameter () = default;

      ///@}

      
      /**
       *  @name Member functions used to get private/protected parameters
       */
      ///@{

      /**
       * @brief get the protected member m_value
       *
       * @return the value of the parameter
       */
      double value () const { return m_value; }

      /**
       * @brief get the protected member m_best_value
       *
       * @return the best-fit value of the parameter
       */  
      double best_value () const { return m_best_value; }

      /**
       * @brief get the protected member m_name
       *
       * @return return the parameter name
       */
      string name () const { return m_name; }

      /**
       * @brief get the i-th chain
       *
       * @param i the i-th chain
       *
       * @return the i-th chain 
       */
      shared_ptr<Chain> chain (const int i) const { return m_chains[i];}

      /**
       * @brief get all the chains
       *
       * @return a vector containing all the chains
       */
      vector<shared_ptr<Chain> > chains () const { return m_chains; }
      
      /**
       * @brief get the protected member m_mean
       *
       * @return the mean value of the parameter
       */
      double mean () const { return m_mean; }

      /**
       * @brief get the protected member m_median
       *
       * @return the median value of the parameter
       */   
      double median () const { return m_median; }

      /**
       * @brief get the protected member m_std
       *
       * @return the standard deviation value of the parameter
       */
      double std () const { return m_std; }

      /**
       * @brief get the defined prior
       * 
       * @return a shared pointer to an object of class Prior
       */
      shared_ptr<Prior> prior () const { return m_prior; }

      /**
       * @brief get the prior interval lenght
       * 
       * @return the prior interval lenght
       */
      double interval_size () const { return (m_prior->xmax()-m_prior->xmin()); }
      
      /**
       * @brief get the value of the prior at the parameter value
       *
       * @return prior value
       */
      double PriorProbability () const { return (isFixed()) ? 1. : m_prior->operator()(m_value); }

      /**
       * @brief value of the prior at a proposed value
       *
       * @param value proporsed parameter value
       *
       * @return prior value
       */ 
      double PriorProbability (const double value) const { return m_prior->operator()(value); }
      
      ///@}

      
      /**
       *  @name Member functions used to set private/protected parameters
       */
      ///@{
      
      /**
       *  @brief set prior seed
       *  @param seed the prior seed
       *  @return none
       */
      void set_prior_seed (const int seed) { if (m_pType==_free_) m_prior->set_seed(seed); }

      /**
       *  @brief set the protected member m_value
       *  @param value parameter value
       *  @return none
       */
      void set_value (const double value); 

      /**
       *  @brief set the protected member m_best_value
       *  @param value parameter best value
       *  @return none
       */
      void set_best_value (const double value) { m_best_value = value; }

      /**
       *  @brief set the protected member m_name
       *  @param name parameter name
       *  @return none
       */
      void set_name (const string name) { m_name = name; }

      /**
       *  @brief set the protected member m_pType
       *  @param pType the parameter type
       *  @return none
       */
      void set_pType (const ParameterType pType) { m_pType = pType; }

      /**
       *  @brief set user defined prior
       *  @param prior user defined prior
       *  @return none
       */
      void set_prior (const shared_ptr<Prior> prior) { m_prior = prior; }

      /**
       *  @brief set the chains
       *
       *  @param nchains the number of chains
       *  @param chain_size the chains size
       *
       *  @return none
       */
      void set_chains (const int nchains, const int chain_size);

      /**
       *  @brief set the chain values
       *
       *  @param chain the n-th chain
       *  @param position the i-th position in the chain
       *  @param value the value for the chain
       *
       *  @return none
       */
      void set_chain_value (const int chain, const int position, const double value);

      /**
       *  @brief set the chains values
       *
       *  @param position the i-th position in the chain
       *  @param values the value for the chain
       *
       *  @return none
       */
      void set_chains_values (const int position, const vector<double> values);

      /**
       *  @brief set the chains values from prior
       *
       *  @param position the i-th position in the chain
       *
       *  @return none
       */
      void set_chains_values_from_prior (const int position);
   
      /**
       *  @brief set the chains values from prior
       *
       *  @param position the i-th position in the chain
       *
       *  @param radius the radius

       *  @param seed the prior seed
       * 
       *  @return none
       */
      void set_chains_values_sphere (const int position, const double radius, const int seed);   

      ///@}
      

      /**
       *  @name Member functions used to operate on parameters
       */
      ///@{

      /**
       *  @brief get the protected member m_pType
       *  @return the parameter type
       */
      bool isFixed () const { return (m_pType==_fixed_) ? true : false; }

      /**
       *  @brief extract a parameter value from the prior
       *  @return a parameter value
       */
      double sample_from_prior () { return (isFixed()) ? m_value : m_prior->sample(); }

      
      /**
       *  @brief extract values from the prior 	 
       *  @param sample_size the size of the extracted sample
       *  @return a vector of values
       */
      vector<double> sample_from_prior (const int sample_size);
      
      /**
       *  @brief extract a parameter value from the prior
       *  @return a parameter value
       */
      vector<double> sample_sphere (const int sample_size, const double radius, const int seed);

      /**
       *  @brief merge the chains
       *
       *  @param max maximum step
       *  @param min minimum step
       *  @param thin pick values every thin times
       *
       *  @return the merged chain
       */
      shared_ptr<Chain> merge_chains (const int max = -1, const int min=-1, const int thin = 1);

      /**
       *  @brief compute the convergence 
       *
       *  @param max maximum step
       *  @param min minimum step
       *  @param thin pick values every thin times
       *
       *  @return the convergence
       */
      double chains_convergence (const int max = -1, const int min=-1, const int thin = 1);

      /**
       *  @brief generate a random parameter value in the prior range
       *  @return the parameter random value
       */
      double random_value () const;

      ///@}
      
    };
  }
}

#endif
