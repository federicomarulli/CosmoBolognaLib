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
 *  @file Headers/NumberCounts1D.h
 *
 *  @brief The class NumberCounts1D
 *
 *  This file defines the interface of the class NumberCounts1D,
 *  used to measure the number counts
 *
 *  @author Federico Marulli
 *
 *  @author federico.marulli3@unbo.it
 */

#ifndef __NCOUNTS1D__
#define __NCOUNTS1D__


#include "Data1D_extra.h"
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
       *  @class NumberCounts1D NumberCounts1D.h
       *  "Headers/NumberCounts1D.h"
       *
       *  @brief The class NumberCounts1D
       *
       *  This is the base class used to measure the 
       *  number counts of one variable
       *
       */
      class NumberCounts1D : public NumberCounts {

	protected:

	  /// the catalogue variable to bin
	  catalogue::Var m_Var;

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
	   *  @return object of class NumberCounts1D
	   */
	  NumberCounts1D () {}

	  /**
	   *  @brief default destructor
	   *  @return none
	   */
	  virtual ~NumberCounts1D () = default;

	  /**
	   *  @brief constructor
	   *
	   *  @param var the variable type 
	   *
	   *  @param bin_type the bin type
	   *
	   *  @param data object of class Catalogue 
	   *
	   *  @param nbins the number of bins
	   *
	   *  @param minVar minimum range 
	   *  
	   *  @param maxVar maximmum range
	   *  
	   *  @param shift the shift of the bin
	   *  
	   *  @param hist_type the type of histogram
	   *
	   *  @param fact factor used to normalized the distribution
	   *
	   *  @return object of class NumberCounts1D
	   */
	  NumberCounts1D (const catalogue::Var var, const BinType bin_type, const catalogue::Catalogue data, const size_t nbins, const double minVar=par::defaultDouble, const double maxVar=par::defaultDouble, const double shift = 0.5, const glob::HistogramType hist_type=glob::HistogramType::_N_V_, const double fact = 1.); 

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
	   *  @param hist vector containing the measures
	   *  used to compute the covariance matrix
	   *  @param JK true &rarr; compute the jackknife covariance
	   *  matrix; false compute the standard covariance matrix
	   *  @return none
	   */
	  void compute_covariance (const std::vector<std::shared_ptr<glob::Histogram>> hist, const bool JK) override;

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
	   *  @brief write the measured covariance
	   *  matrix
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
