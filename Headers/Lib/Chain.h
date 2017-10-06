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
 *  @file Headers/Lib/Chain.h
 *
 *  @brief The class Chain
 *
 *  This file defines the interface of the class Chain
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#ifndef __CHAIN__
#define __CHAIN__

#include "Model2D.h"
#include "Posterior.h"

namespace cosmobl {

  namespace statistics {

    /**
     *  @class Chain Chain.h "Headers/Lib/Chain.h"
     *
     *  @brief The class Chain
     *
     *  This class is used to define the chains, output of the montecarlo
     *  process
     */
    class Chain {

      protected:

	/// the lenght of the chain
	int m_size;

	/// the number of parallel walkers
	int m_nwalkers;

	/// content of the chain
	vector<double> m_values;

	/**
	 *  @brief get the position in the vector
	 *  m_values from position index and walker index
	 *
	 *  @param pp the position in the chain
	 *
	 *  @param ww the walker
	 *
	 *  @return object of class Chain.  
	 */
	int m_inds_to_index(const int pp, const int ww) const { return pp*m_nwalkers+ww;}


      public:

	/**
	 *  @name Constructors/destructors
	 */
	///@{

	/**
	 *  @brief default constructor
	 *
	 *  @return object of class Chain.  
	 */
	Chain () = default;

	/**
	 *  @brief constructor
	 *
	 *  @param size size of the chain 
	 *
	 *  @param nwalkers the number of parallel walkers
	 *
	 *  @return object of class Chain.
	 */
	Chain (const int size, const int nwalkers);

	/**
	 *  @brief constructor
	 *
	 *  @param values the input chain values
	 *
	 *  @param nwalkers the number of parallel walkers
	 *
	 *  @return object of class Chain.
	 */
	Chain (const vector<double> values, const int nwalkers);

	/**
	 *  @brief constructor
	 *
	 *  @param values the input chain values
	 *
	 *  @return object of class Chain.
	 */
	Chain (vector<vector<double>> values);

	/**
	 *  @brief default destructor
	 *  @return none
	 */
	~Chain () = default;

	///@}

	/**
	 * @brief return the private member m_size
	 *
	 * @return the chain size
	 */
	int size () const { return m_size; }

	/**
	 * @brief return the private member m_nwalkers
	 *
	 * @return the chain size
	 */
	int nwalkers () const { return m_nwalkers; }

	/**
	 * @brief set the chain
	 *
	 * @param size the chain lenght 
         *
	 * @param nwalkers the number of parallel walkers
	 *
	 * @return none
	 */
	void set_chain (const int size, const int nwalkers);

	/**
	 * @brief reset the chain using m_size and m_nwalkers
	 *
	 * @return none
	 */
	void reset();

	/**
	 * @brief expand the already existing chain
	 *
	 * @param append the lenght of the empty chunk of the chain 
	 *
	 * @return none
	 */
	void expand(const int append);

	/**
	 * @brief return the private member m_values at the pp-th step
	 * for the ww-th step
	 * 
	 * @param pp the chain step
	 * 
	 * @param ww the walker index
	 *
	 * @return the chain value
	 */
	double value (const int pp, const int ww) const { return m_values[m_inds_to_index(pp, ww)]; }

	/**
	 * @brief set the private member m_values at the pp-th step
	 * for the ww-th step
	 * 
	 * @param pp the chain step
	 * 
	 * @param ww the walker index
	 *
	 * @param value the chain value
	 *
	 * @return none
	 */
	void set_value (const int pp, const int ww, const double value) { m_values[m_inds_to_index(pp, ww)] = value; }

	/**
	 *  @brief set the chain values
	 *
	 *  @param values the input chain values
	 *
	 *  @param nwalkers the number of parallel walkers
	 *
	 *  @return none
	 */
	void set_values(const vector<double> values, const int nwalkers);

	/**
	 *  @brief set the chain values
	 *
	 *  @param values the input chain values
	 *
	 *  @return none
	 */
	void set_values(const vector<vector<double>> values);

	/**
	 * @brief return the chain
	 * 
	 * @param start the starting point 
	 * 
	 * @param thin number of jumped indexes in the chain
	 *
	 * @return the  chain values
	 */
	vector<double> values(const int start, const int thin = 1) const;

	/**
	 * @brief get the posterior distribution from the chain values
	 * 
	 * @param start the starting point 
	 * 
	 * @param thin number of jumped indexes in the chain
	 *
	 * @param seed the random number generator seed 
	 *
	 * @return the posterior distribution
	*/
	shared_ptr<statistics::Posterior> PosteriorDistribution(const int start, const int thin, const int seed=43241) const;

	/**
	 * @brief get the posterior distribution from the chain values
	 * 
	 * @param seed the random number generator seed 
	 *
	 * @return the posterior distribution
	*/
	shared_ptr<statistics::Posterior> PosteriorDistribution(const int seed=43241) const;

    };
  }
}

#endif
