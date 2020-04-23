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
 *  @file Headers/NumberCounts.h
 *
 *  @brief The class NumberCounts
 *
 *  This file defines the interface of the class NumberCounts,
 *  used to measure the number counts
 *
 *  @author Federico Marulli
 *
 *  @author federico.marulli3@unibo.it
 */

#ifndef __NCOUNTS__
#define __NCOUNTS__


#include "Catalogue.h"
#include "Measure.h"
#include "Histogram.h"


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
       *  @class NumberCounts NumberCounts.h
       *  "Headers/NumberCounts.h"
       *
       *  @brief The class NumberCounts
       *
       *  This is the base class used to measure the 
       *  number counts
       *
       */
      class NumberCounts : public Measure {

	protected :

	  /**
	   *  @name Input catalogue
	   */
	  ///@{

	  /// input data catalogue
	  std::shared_ptr<catalogue::Catalogue> m_data;

	  ///@}

	  /**
	   *  @name Binned data
	   */
	  ///@{
	  
	  /// the histogram type
	  glob::HistogramType m_HistogramType;

	  /// the normalization factor
	  double m_fact;

	  /// number counts type
	  std::shared_ptr<glob::Histogram> m_histogram;

	  ///@}
	  
	  ///@}

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
	  virtual std::shared_ptr<data::Data> m_measurePoisson () 
	  { ErrorCBL("", "m_measurePoisson", "NumberCounts.h"); return NULL; }

	  /**
	   *  @brief measure the number counts with Jackknife 
	   *  covariance matrix
	   *
	   *  @param dir_output_resample output directory of the
	   *  resampling correlation functions
	   *
	   *  @return none
	   */
	  virtual std::shared_ptr<data::Data> m_measureJackknife (const std::string dir_output_resample=par::defaultString)
	  { (void)dir_output_resample; ErrorCBL("", "m_measureJackknife", "NumberCounts.h"); return NULL; }

	  /**
	   *  @brief measure the number counts with
	   *  Bootstrap covariance matrix
	   *
	   *  @param dir_output_resample output directory of the
	   *  resampling correlation functions
	   *
	   *  @param nResamplings number of resampling used for bootstrap
	   *
	   *  @param seed the seed for random number generation
	   *
	   *  @return none
	   */
	  virtual std::shared_ptr<data::Data> m_measureBootstrap (const std::string dir_output_resample=par::defaultString, const int nResamplings=0, const int seed=3213)
	  { (void)dir_output_resample; (void)nResamplings; (void)seed; ErrorCBL("", "m_measureBootstrap", "NumberCounts.h");  return NULL; }

	  ///@}

	public:

	  /**
	   *  @name Constructors/destructors
	   */
	  ///@{

	  /**
	   *  @brief default constructor
	   *
	   *  @return object of class NumberCounts
	   */
	  NumberCounts () {}

	  /**
	   *  @brief default destructor
	   *  @return none
	   */
	  virtual ~NumberCounts () = default;

	  ///@}

	  /**
	   *  @name Member functions to set protected members
	   */
	  ///@{

	  /**
	   *  @brief add a data catalogue
	   *  @param data object of class Catalogue 
	   *  @return none
	   */
	  void set_data (const catalogue::Catalogue data)
	  { m_data = std::make_shared<catalogue::Catalogue>(catalogue::Catalogue(std::move(data))); }

	  ///@}

	  /**
	   *  @name Member functions to get protected members
	   */
	  ///@{

	  /**
	   *  @brief function to get the protected member m_data
	   *  @return return the protected member m_data
	   */
	  catalogue::Catalogue catalogue () { return m_data; }

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
	   *  resampling correlation functions
	   *
	   *  @param nResamplings number of resampling used for bootstrap
	   *
	   *  @param seed the seed for random number generation
	   *
	   *  @param conv true &rarr; compute the Gaussian convolvolution of
	   *  the distribution; false &rarr; do not convolve
	   *
	   *  @param sigma &sigma; of the Gaussian kernel
	   *
	   *  @return none
	   */
	  virtual void measure (const ErrorType errorType=ErrorType::_Poisson_, const std::string dir_output_resample=par::defaultString, const int nResamplings=0, const int seed=3213, const bool conv=false, const double sigma=0.)
	  { (void)errorType; (void)dir_output_resample; (void)nResamplings; (void)seed; (void)conv; (void)sigma; ErrorCBL("", "measure", "NumberCounts.h"); }

	  ///@}

	  /**
	   *  @name input/output member functions (customized in all the derived classes)
	   */
	  ///@{

	  /**
	   *  @brief write the measured number counts
	   *  @param dir output directory
	   *  @param file output file
	   *  @param rank cpu index (for mpi usage)
	   *  @return none
	   */
	  virtual void write (const std::string dir=par::defaultString, const std::string file=par::defaultString, const int rank=0) const
	  { (void)dir; (void)file; (void)rank; ErrorCBL("", "write", "NumberCounts.h"); }


	  ///@}

	  /**
	   *  @name Member functions to estimate the errors and covariance matrices
	   */
	  ///@{ 

	  /**
	   *  @brief write the measured covariance
	   *  matrix
	   *  @param dir output directory
	   *  @param file output file
	   *  @return none
	   */
	  virtual void write_covariance (const std::string dir, const std::string file) const
	  { (void)dir; (void)file; ErrorCBL("", "write_covariance", "NumberCounts.h"); }

	  /**
	   *  @brief compute the covariance matrix
	   *  @param histo vector containing the measures 
	   *  used to compute the covariance matrix
	   *  @param JK true &rarr; compute the jackknife covariance
	   *  matrix; false compute the standard covariance matrix
	   *  @return none
	   */
	  virtual void compute_covariance (const std::vector<std::shared_ptr<glob::Histogram>> histo, const bool JK)
	  { (void)histo; (void)JK; ErrorCBL("", "compute_covariance", "NumberCounts.h"); }
	  	  
	  /**
	   *  @brief apply a Gaussian filter to the distribution
	   *  @param &sigma; of the Gaussian kernel
	   *  @return none
	   */
	  virtual std::shared_ptr<data::Data> Gaussian_smoothing (const double sigma)
	  { (void)sigma; ErrorCBL("", "Gaussian_smoothing", "NumberCounts.h"); return NULL; }

	  ///@}

	  /**
	   *  @name Functions to get the private members of the class
	   */
	  ///@{
	  	  
	  /**
	   * @brief return the binned counts
	   *
	   * @return pointer to an object of class Histogram
	   */
	  std::shared_ptr<glob::Histogram> histogram () { return m_histogram; }

	  /**
	   * @brief return the type of histogram
	   * normalization
	   *
	   * @return the type of histogram normalization
	   */
	  glob::HistogramType HistogramType () { return m_HistogramType; }

	  /**
	   * @brief return the normalization factor
	   *
	   * @return the normalization factor
	   */
	  double fact () { return m_fact; }

	  ///@}


      };
    }
  }
}

#endif
