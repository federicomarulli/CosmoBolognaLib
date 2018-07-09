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
 *  @file Headers/Sampler.h
 *
 *  @brief The class Sampler
 *
 *  This file defines the interface of the class Sampler, used for
 *  statistical analyses and Bayesian inference
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#ifndef __SAMP__
#define __SAMP__

#include "Distribution.h"


// ============================================================================================


namespace cbl {

  namespace statistics {

    /**
     * @enum SamplerType
     * @brief the parameter type
     */
    enum class SamplerType {

      /// Metropolis-Hastings sampler
      _MetropolisHastings_,

      /// stretch-move sampler (Goodman & Weare 2010, Foreman-Mackey et al. 2012) 
      _StretchMove_,
      
    };

   /**
     * @brief return a vector containing the
     * SamplerType names
     * @return a vector containing the
     * SamplerType names
     */
    inline std::vector<std::string> SamplerTypeNames () { return {"MetropolisHastings", "StretchMove"}; }

    /**
     *
     * @brief cast an enum of type SamplerType
     * from its index
     * @param samplerTypeIndex the samplerType index
     * @return object of class SamplerType
     */
    inline SamplerType SamplerTypeCast (const int samplerTypeIndex) { return castFromValue<SamplerType>(samplerTypeIndex); }

    /**
     * @brief cast an enum of type SamplerType
     * from its name
     * @param samplerTypeName the samplerType name
     * @return object of class SamplerType
     */
    inline SamplerType SamplerTypeCast (const std::string samplerTypeName) { return castFromName<SamplerType>(samplerTypeName, SamplerTypeNames()); }

    /**
     * @brief cast an enum of type SamplerType
     * from indeces
     * @param samplerTypeIndeces the samplerType indeces
     * @return object of class SamplerType
     */
    inline std::vector<SamplerType> SamplerTypeCast (const std::vector<int> samplerTypeIndeces) { return castFromValues<SamplerType>(samplerTypeIndeces); } 

    /**
     * @brief cast enums of type SamplerType
     * from their names
     * @param samplerTypeNames the samplerType names
     * @return vector of SamplerType enums
     */
    inline std::vector<SamplerType> SamplerTypeCast (const std::vector<std::string> samplerTypeNames) { return castFromNames<SamplerType>(samplerTypeNames, SamplerTypeNames()); }
    
    /**
     *  @class Sampler Sampler.h "Headers/Sampler.h"
     *
     *  @brief The class Sampler
     *
     *  This class is used to handle objects of type sample. It samples
     *  generic functions
     */
    class Sampler
    {
      
    protected:

      /// number of chains
      int m_nwalkers;

      /// size of the chains
      int m_chain_size;

      /// number of parameters
      int m_npar;

      /// number of free parameters
      int m_npar_free;

      /// the function to be sampled
      std::function<double(std::vector<double> &)> m_function;

      /// chain acceptance ratio
      std::vector<double> m_acceptance;

      /// use python-defined function
      bool m_use_python;

      /// value of the function at sampled points
      std::vector<std::vector<double>> m_function_chain;

      /// the chains
      std::vector<std::vector<std::vector<double>>> m_chains;

      /**
       * @brief return the random generator for the 
       * stretch-move 
       *
       * @param seed the random generator seed
       * @param aa the stretch-move distribution parameter
       * 
       * @return pointer to the random number generator
       */
      std::shared_ptr<random::DistributionRandomNumbers>m_set_gz (const int seed, const double aa=2);

      /**
       * @brief initialize chains, generating random
       * points in a sphere around starting position
       *
       * @param start vector containing the starting position
       * for the parameters
       * 
       * @return none
       */
      void m_initialize_chains (const std::vector<std::vector<double>> start);

      /**
       *  @brief samples lthe input function, using stretch-move
       *  algorithm
       *  on n-dimensional parameter space. Parallel version
       *  
       *  @param chain_size number of step in each chain 
       *  @param nwalkers number of parallel walkers 
       *  @param start vector containing the starting position
       *  for the parameters
       *  @param seed the seed for random number generator
       *  @param aa the stretch-move distribution parameter
       *
       *  @return none
       */
      void m_sample_stretch_move_parallel_cpp (const int chain_size, const int nwalkers, const std::vector<std::vector<double>> start, const int seed=4241, const double aa=2);

      /**
       *  @brief samples lthe input function, using stretch-move
       *  algorithm on n-dimensional parameter space. Parallel version
       *  Special function for function written in python
       *
       *  @param chain_size number of step in each chain 
       *  @param nwalkers number of parallel walkers   
       *  @param start vector containing the starting position
       *  for the parameters
       *  @param seed the seed for random number generator
       *  @param aa the stretch-move distribution parameter
       *
       *  @return none
       */
      void m_sample_stretch_move_parallel_py (const int chain_size, const int nwalkers, const std::vector<std::vector<double>> start, const int seed=4241, const double aa=2);

    public:

      /**
       *  @name Constructors/destructors
       */
      ///@{

      /**
       *  @brief default constructor
       *  @return object of class Likelihood
       */
      Sampler () {}

      /**
       *  @brief constructor
       *  @param npar the number of parameters
       *  @param function the function to sample
       *  @return object of class Sampler
       */
      Sampler (const int npar, const std::function<double(std::vector<double> &)> function) : m_npar(npar), m_npar_free(npar) {set_function(function);}

      /**
       *  @brief constructor
       *  @param npar the number of parameters
       *  @param npar_free the number of free parameters
       *  @param function the function to sample
       *  @return object of class Sampler
       */
      Sampler (const int npar, const int npar_free, const std::function<double(std::vector<double> &)> function) : m_npar(npar), m_npar_free(npar_free) {set_function(function);}

      /**
       *  @brief default destructor
       *  @return none
       */
      ~Sampler () = default;

      ///@}

      /**
       * @brief evaluate the function and 
       * return its value and derived parameters
       *
       * @param pp the function parameters
       *
       * @return none
       */
      double operator () (std::vector<double> &pp) 
      {	
	return m_function(pp);
      }

      /**
       * @brief return the chain value
       *
       * @param par the parameter
       *
       * @param chain the chain number
       *
       * @param step the step in the chain
       *
       * @return none
       */
      double get_chain (const int par, const int chain, const int step) {return m_chains[step][chain][par];}

      /**
       * @brief return the chains
       *
       * @return vector containing the chains
       */
      std::vector<std::vector<std::vector<double>>> get_chain () {return m_chains;}

      /**
       * @brief return the function value
       *
       * @param chain the chain number
       *
       * @param step the step in the chain
       *
       * @return none
       */
      double get_function (const int chain, const int step) {return m_function_chain[step][chain];}

      /**
       *  @brief get the chain values and the function
       *  
       *  @param[out] chains vector containing the chains
       *  @param[out] function the function computed at each step of the chain 
       *  @param[out] acceptance the acceptance rate
       *  @param start the starting position for each chain
       *  @param thin the position step
       *
       *  @return none
       */
      void get_chain_function_acceptance(std::vector<std::vector<double>> &chains, std::vector<double> &function, std::vector<double> &acceptance, const int start=0, const int thin=1);

      /**
       * @brief function to set the chains
       *
       * @param npar the number of parameters
       * @param npar_free the number of free parameters
       * @param chain_size number of step in each chain 
       * @param nwalkers number of parallel walkers 
       *
       * @return none
       */
      void set_chain (const int npar, const int npar_free, const int chain_size, const int nwalkers);

      /**
       * @brief set the function 
       *
       * @param function the function to be sampled
       *
       * @return none
       */
      void set_function (const std::function<double(std::vector<double> &)> function);

      /**
       *  @brief samples lthe input function, using stretch-move
       *  algorithm on n-dimensional parameter space 
       *
       *  @param chain_size number of step in each chain 
       *  @param nwalkers number of parallel walkers  
       *  @param start vector containing the starting position
       *  for the parameters
       *  @param seed the seed for random number generator
       *  @param aa the stretch-move distribution parameter
       *  
       *  @return averace acceptance ratio
       */
      void sample_stretch_move (const int chain_size, const int nwalkers, const std::vector<std::vector<double>> start, const int seed=4241, const double aa=2);

      /**
       *  @brief samples lthe input function, using stretch-move
       *  algorithm on n-dimensional parameter space. Parallel version
       *
       *  @param chain_size number of step in each chain 
       *  @param nwalkers number of parallel walkers  *
       *  @param start vector containing the starting position
       *  for the parameters
       *  @param seed the seed for random number generator
       *  @param aa the stretch-move distribution parameter
       *
       *  @return none
       */
      void sample_stretch_move_parallel (const int chain_size, const int nwalkers, const std::vector<std::vector<double>> start, const int seed=4241, const double aa=2);

      /**
       *  @brief write the chains in an output file
       *  
       *  @param dir_output the output directory
       *  @param file the output file
       *  @param start the starting position for each chain
       *  @param thin the position step
       *
       *  @return none
       */
      void write_chain(const std::string dir_output, const std::string file, const int start, const int thin);
      
    };
  }
}

#endif
