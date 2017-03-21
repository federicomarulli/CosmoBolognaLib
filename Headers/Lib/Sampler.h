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
 *  @file Headers/Lib/Sampler.h
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

#include "Chi2.h"


// ============================================================================================


namespace cosmobl {

  namespace statistics {
    
    /**
     *  @class Likelihood Likelihood.h "Headers/Lib/Likelihood.h"
     *
     *  @brief The class Likelihood
     *
     *  This class is used to handle objects of type likelihood. It is
     *  used for all kind of likelihood analyses, sample and minimization
     */
    class Sampler
    {
      
    protected:

      /// number of chains
      int m_nchains;

      /// size of the chains
      int m_chain_size;

      /// number of parameters
      int m_npar;

      /// the function to be sampled
      function<double(vector<double>)> m_function;

      /// chain acceptance ratio
      vector<double> m_acceptance;

      /// value of the function at sampled points
      vector<vector<double>> m_function_chain;

      /// the chains
      vector<vector<vector<double>>> m_chains;

      /**
       * @brief return the random generator for the 
       * stretch-move
       *
       * @param seed the random generator seed
       * @param aa the stretch-move distribution parameter
       * 
       * @return pointer to the random number generator
       */
      shared_ptr<random::DistributionRandomNumbers>m_set_gz(const int seed, const double aa=2);

      /**
       * @brief initialize chains, generating random
       * points in a sphere around starting position
       *
       * @param seed the random generator seed
       * @param start vector containing the starting position
       * for the parameters
       * @param radius the sphere radius
       * 
       * @return none
       */
      void m_initialize_chains(const int seed, const vector<double> start, const double radius=1.e-3);

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
       *  @param data pointers to the data container
       *  @param model pointers to the model 
       *  @param likelihood_type type of likelihood
       *  @return object of class Likelihood
       */
      Sampler (const int npar, const function<double(vector<double>)> function) : m_npar(npar), m_function(function) {}

      /**
       *  @brief default destructor
       *  @return none
       */
      ~Sampler () = default;

      ///@}

      /**
       * @brief evaluate the function
       *
       * @param pp the function parameters
       *
       * @return none
       */
      double operator () (vector<double> pp) 
      {	
	return m_function(pp);
      }

      /**
       * @brief function to set the chains
       *
       * @param npar the number of parameters
       * @param nchains the number of chains
       * @param chain_size the size of the chain
       *
       * @return none
       */
      void set_chain(const int npar, const int nchains, const int chain_size);

      /**
       * @brief set the function 
       *
       * @param function the function to be sampled
       *
       * @return none
       */
      void set_function(const function<double(vector<double>)> function);

      /**
       *  @brief samples lthe input function, using stretch-move
       *  algorithm on n-dimensional parameter space 
       *    
       *  @param nchains number of chains to sample the parameter space 
       *  @param chain_size number of step in each chain 
       *  @param start vector containing the starting position
       * for the parameters
       *  @param radius the sphere radius
       *  @param seed the seed for random number generator
       *  @param aa the stretch-move distribution parameter
       *  
       *  @return averace acceptance ratio
       */
      void sample_stretch_move (const int nchains, const int chain_size, const vector<double> start, double radius=1.e-3, const int seed=4241, const double aa=2);

      /**
       *  @brief samples lthe input function, using stretch-move
       *  algorithm on n-dimensional parameter space. Parallel version
       *  
       *  @param nchains number of chains to sample the parameter space 
       *  @param chain_size number of step in each chain 
       *  @param start vector containing the starting position
       *  for the parameters
       *  @param radius the sphere radius
       *  @param seed the seed for random number generator
       *  @param aa the stretch-move distribution parameter
       *
       *  @return none
       */
      void sample_stretch_move_parallel (const int nchains, const int chain_size, const vector<double> start, double radius=1.e-3, const int seed=4241, const double aa=2);

      /**
       *  @brief write the chains in an output file
       *  
       *  @param dir_output the output directory
       *  @param file the output file
       *  @param start the starting position for each chain
       *  @param stop the ending position for each chain
       *  @param thin the position step
       *
       *  @return none
       */
      void write_chain(const string dir_output, const string file, const int start, const int stop, const int thin);
      
    };
  }
}

#endif
