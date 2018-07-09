/********************************************************************
 *  Copyright (C) 2010 by Federico Marulli                          *
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
 *  @file Headers/NumberCounts2D.h
 *
 *  @brief The class NumberCounts2D
 *
 *  This file defines the interface of the class NumberCounts2D,
 *  used to measure the number counts
 *
 *  @author Federico Marulli
 *
 *  @author federico.marulli3@unbo.it
 */

#ifndef __NCOUNTS2D__
#define __NCOUNTS2D__


#include "NumberCounts.h"


// ===================================================================================================


namespace cbl {

  namespace measure {

    /**
     *  @brief The namespace of the <B> number counts
     *  </B>
     *  
     *  The \e measure::numbercounts namespace contains all the functions and
     *  classes to measure the number counts
     */
    namespace numbercounts {

      /**
       *  @class NumberCounts2D NumberCounts2D.h
       *  "Headers/NumberCounts2D.h"
       *
       *  @brief The class NumberCounts2D
       *
       *  This is the base class used to measure the 
       *  number counts of one variable
       *
       */
      class NumberCounts2D : public NumberCounts {

	protected:

	  /// the first catalogue variable to bin
	  catalogue::Var m_Var1;

	  /// the second catalogue variable to bin
	  catalogue::Var m_Var2;	 
	 
	  /**
	   *  @name Protected member functions to measure the number counts
	   */
	  ///@{

	  /**
	   *  @brief measure the number counts with Poisson 
	   *  errors
	   *
	   *  @return none
	   */
	  std::shared_ptr<data::Data> m_measurePoisson () override;

	  /**
	   *  @brief measure the number counts with Jackknife 
	   *  covariance matrix
	   *
	   *  @param dir_output_resample output directory of the
	   *  resampled correlation function
	   *
	   *  @return none
	   */
	  std::shared_ptr<data::Data> m_measureJackknife (const std::string dir_output_resample=par::defaultString) override;

	  /**
	   *  @brief measure the number counts with
	   *  Bootstrap covariance matrix
	   *
	   *  @param dir_output_resample output directory of the
	   *  resampled correlation function
	   *
	   *  @param nResamplings number of resampling used for bootstrap
	   *
	   *  @param seed the seed for random number generation
	   *
	   *  @return none
	   */
	  std::shared_ptr<data::Data> m_measureBootstrap (const std::string dir_output_resample=par::defaultString, const int nResamplings=0, const int seed=3213) override;

	  ///@}

	public:

	  /**
	   *  @name Constructors/destructors
	   */
	  ///@{

	  /**
	   *  @brief default constructor
	   *
	   *  @return object of class NumberCounts2D
	   */
	  NumberCounts2D () {}

	  /**
	   *  @brief default destructor
	   *  @return none
	   */
	  virtual ~NumberCounts2D () = default;

	  /**
	   * @brief constructor
	   *
	   * @param var1 the first variable type 
	   *
	   * @param bin_type1 the first bin type
	   *
	   * @param var2 the second variable type 
	   *
	   * @param bin_type2 the second bin type
	   * 
	   * @param data object of class Catalogue 
	   *
	   * @param nbins1 the number of bins for the first variable
	   *
	   * @param nbins2 the number of bins for the second variable
	   *
	   * @param minVar1 minimum range  for the first variable
	   * 
	   * @param maxVar1 maximmum range for the first variable
	   * 
	   * @param minVar2 minimum range for the second variable 
	   * 
	   * @param maxVar2 maximmum range for the second variable
	   * 
	   * @param shift1 bin shift for the first variable
	   * 
	   * @param shift2 bin shift for the second variable
	   * 
	   * @param hist_type the type of histogram
	   *
	   * @param fact factor used to normalized the distribution
	   *
	   * @return object of class Histogram1D
	   */
	  NumberCounts2D (const catalogue::Var var1, const BinType bin_type1, const catalogue::Var var2, const BinType bin_type2, const catalogue::Catalogue data, const size_t nbins1, const size_t nbins2, const double minVar1=par::defaultDouble, const double maxVar1=par::defaultDouble, const double minVar2=par::defaultDouble, const double maxVar2=par::defaultDouble, const double shift1=0.5, const double shift2=0.5, const glob::HistogramType hist_type=glob::HistogramType::_N_V_, const double fact = 1.); 

	  ///@}


	  /**
	   *  @name Member functions to measure the number counts
	   */
	  ///@{

	  /**
	   *  @brief measure the number counts
	   *
	   *  @param errorType type of error
	   *
	   *  @param dir_output_resample output directory of the
	   *  resampled correlation function
	   *
	   *  @param nResamplings number of resampling used for bootstrap
	   *
	   *  @param seed the seed for random number generation
	   *
	   *  @return none
	   */
	  void measure (const ErrorType errorType=ErrorType::_Poisson_, const std::string dir_output_resample=par::defaultString, const int nResamplings=0, const int seed=3213) override;


	  /**
	   *  @brief compute the covariance matrix
	   *  @param histo vector containing the measures
	   *  used to compute the covariance matrix
	   *  @param JK true &rarr; compute the jackknife covariance
	   *  matrix; false compute the standard covariance matrix
	   *  @return none
	   */
	  void compute_covariance (const std::vector<std::shared_ptr<glob::Histogram>> histo, const bool JK) override;

	  ///@}

	  /**
	   *  @name Input/Output member functions (customized in all the derived classes)
	   */
	  ///@{

	  /**
	   *  @brief write the measured number counts
	   *  @param dir output directory
	   *  @param file output file
	   *  @param rank cpu index (for MPI usage)
	   *  @return none
	   */
	  void write (const std::string dir=par::defaultString, const std::string file=par::defaultString, const int rank=0) const override;

	  /**
	   *  @brief write measured covariance matrix
	   *  @param dir output directory
	   *  @param file output file
	   *  @return none
	   */
	  void write_covariance (const std::string dir, const std::string file) const override;

	  ///@}
      }; 
    }
  }
}

#endif
