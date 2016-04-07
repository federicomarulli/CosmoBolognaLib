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

#include "Data2D.h"

namespace cosmobl {

  /**
   *  @brief The namespace of functions and classes used for statistical
   *  analysis
   *  
   * The \e statistic namespace contains all the functions and classes
   * used for statistical analyis
   */
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
	int m_chain_size;

	/// content of the chain
	vector<double> m_values;

	/// the chain mean value 
	double m_mean;

	/// the standard deviation of chain values 
	double m_std;

	/// the chain median value
	double m_median;

	/// var 
	vector<double> m_var;

	/// dist
	vector<double> m_dist;

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
	 *  @param chain_size size of the chain 
	 *
	 *  @return object of class Chain.
	 */
	Chain (const int chain_size) : m_chain_size(chain_size) { m_values.resize(m_chain_size, 0); }

	/**
	 *  @brief default destructor
	 *  @return none
	 */
	~Chain () = default;

	///@}


	/**
	 *  @brief compute statistics (mean, std, median) of the chain
	 *  
	 *  @param max maximum step of the chain to use
	 *  @param min minumium step of the chain to use
	 *
	 *  @return none
	 */
	void Statistics (const int max=-1, const int min=-1);

	/**
	 *  @brief compute chain distribution 
	 *  
	 *  @param nbin numbers of bin
	 *
	 *  @return none
	 */
	void ComputeDistribution(const int nbin);

	/**
	 *  @brief return the private member m_var
	 *  
	 *  @return the range of the binned chain values 
	 */
	vector<double> var () const { return m_var; }

	/**
	 *  @brief return the private member m_dist
	 *  
	 *  @return the distribution of chain values 
	 */
	vector<double> dist () const { return m_dist; }

	/**
	 *  @brief return the private member m_mean 
	 *  
	 *  @return the chain mean value
	 */
	double mean () const { return m_mean; }

	/**
	 *  @brief return the private member m_median
	 *  
	 *  @return the chain median value
	 */
	double median () const { return m_median;}

	/**
	 *  @brief return the private member m_std 
	 *  
	 *  @return the chain standard deviation value
	 */
	double std () const { return m_std; }

	/**
	 *  @brief set the chain size
	 *  
	 *  @param chain_size the chain size
	 *
	 *  @return none
	 */
	void set_chain_size (const int chain_size) { m_chain_size=chain_size; m_values.resize(m_chain_size, 0); }

	/**
	 *  @brief set the i-th chain value
	 *
	 *  @param i the i-th chain step
	 *  @param value the value at the i-th step
	 *
	 *  @return none
	 */
	void set_chain_value (const int i, const double value) { m_values[i] = value; }

	/**
	 * @brief return the private member m_chain_size
	 *
	 * @return the chain size
	 */
	int chain_size () const { return m_chain_size; }

	/**
	 * @brief return the private member m_values at the i-th step
	 * 
	 * @param i the i-th step
	 *
	 * @return the chain value at i-th step
	 */
	double chain_value (const int i) const { return m_values[i]; } 

    };
  }
}

#endif
