/********************************************************************
 *  Copyright (C) 2010 by Federico Marulli and Alfonso Veropalumbo  *
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
 *  @file Headers/TwoPointCorrelation1D_filtered.h
 *
 *  @brief The class TwoPointCorrelation1D_filtered
 *
 *  This file defines the interface of the class
 *  TwoPointCorrelation1D_filtered, used to measure the filtered
 *  monopole of the two-point correlation function
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unibo.it, alfonso.veropalumbo@unibo.it
 */

#ifndef __TWOPOINTFIL__
#define __TWOPOINTFIL__


#include "TwoPointCorrelation1D_monopole.h"


// ===================================================================================================


namespace cbl {
  
  namespace measure {

    namespace twopt {
  
      /**
       *  @class TwoPointCorrelation1D_filtered
       *  TwoPointCorrelation1D_filtered.h
       *  "Headers/TwoPointCorrelation1D_filtered.h"
       *
       *  @brief The class TwoPointCorrelation1D_filtered
       *
       *  This class is used to handle objects of type <EM>
       *  TwoPointCorrelation1D_filtered </EM>. It is used to measure
       *  the filtered monopole of the two-point correlation function,
       *  \f$w(r_c)=2\pi \int dr \xi(r) W(r,r_c) r^2\f$, where \f$W(x) =
       *  (2x)^2(1-x)(0.5-x)r_c^{-3}\f$, and \f$\xi(r)\f$ is the
       *  monopole of the two-point correlation function
       *   
       */
      class TwoPointCorrelation1D_filtered : public TwoPointCorrelation1D_monopole
      {

      protected:
      
	/// scales at which the filtered correlation function is computed
	std::vector<double> m_rc;

	/**
	 *  @brief measure the filtered two-point correlation function,
	 *  \f$w(r_c)=2\pi \int dr \xi(r) W(r,r_c) r^2\f$, where \f$W(x)
	 *  = (2x)^2(1-x)(0.5-x)r_c^{-3}\f$
	 *  
	 *  @param data pointer to an object of type Pair containing the
	 *  data-data pairs
	 *
	 *  @return pointer to an object of type Data
	 */
	std::shared_ptr<data::Data> Filtered (const std::shared_ptr<data::Data> data) override;

	/**
	 *  @brief measure the filtered two-point correlation function
	 *  with Poisson errors
	 *
	 *  @param dir_output_pairs output directory used to store the
	 *  number of pairs
	 *
	 *  @param dir_input_pairs vector of input directories used to
	 *  store the number of pairs (if the pairs are read from files)
	 *
	 *  @param count_dd true &rarr; count the number of data-data
	 *  pairs; false &rarr; read the number of data-data pairs from
	 *  file
	 *
	 *  @param count_rr true &rarr; count the number of
	 *  random-random pairs; false &rarr; read the number of
	 *  random-random pairs from file
	 *
	 *  @param count_dr true &rarr; count the number of data-random
	 *  pairs; false &rarr; read the number of data-random pairs
	 *
	 *  @param tcount true &rarr; activate the time counter; false
	 *  &rarr; no time counter
	 *
	 *  @param estimator the estimator used to measure the two-point
	 *  correlation function
	 *
	 *  @return none
	 */
	void measurePoisson (const std::string dir_output_pairs=par::defaultString, const std::vector<std::string> dir_input_pairs={}, const bool count_dd=true, const bool count_rr=true, const bool count_dr=true, const bool tcount=true, const Estimator estimator=Estimator::_LandySzalay_) override;

	/**
	 *  @brief measure the filtered two-point correlation function
	 *  estimating the covariance with Jackknife resampling
	 *
	 *  @param dir_output_pairs output directory used to store the
	 *  number of pairs
	 *
	 *  @param dir_input_pairs vector of input directories used to
	 *  store the number of pairs (if the pairs are read from files)
	 *
	 *  @param dir_output_resample output directory used to store
	 *  the Jackknife resampling correlation functions, with
	 *  Poisson errors; if an empty string (i.e. "" or "NULL") is
	 *  provided, no output will be stored
	 *
	 *  @param count_dd true &rarr; count the number of data-data
	 *  pairs; false &rarr; read the number of data-data pairs from
	 *  file
	 *
	 *  @param count_rr true &rarr; count the number of
	 *  random-random pairs; false &rarr; read the number of
	 *  random-random pairs from file
	 *
	 *  @param count_dr true &rarr; count the number of data-random
	 *  pairs; false &rarr; read the number of data-random pairs
	 *
	 *  @param tcount true &rarr; activate the time counter; false
	 *  &rarr; no time counter
	 *
	 *  @param estimator the estimator used to measure the two-point
	 *  correlation function
	 *
	 *  @return none
	 */
	void measureJackknife (const std::string dir_output_pairs = par::defaultString, const std::vector<std::string> dir_input_pairs={}, const std::string dir_output_resample=par::defaultString, const bool count_dd=true, const bool count_rr=true, const bool count_dr=true, const bool tcount=true, const Estimator estimator=Estimator::_LandySzalay_) override;

	/**
	 *  @brief measure the filtered two-point correlation function
	 *  estimating the covariance with Bootstrap resampling
	 *
	 *  @param nMocks number of mocks to be generated with bootstrap
	 *  resampling
	 *
	 *  @param dir_output_pairs output directory used to store the
	 *  number of pairs
	 *
	 *  @param dir_input_pairs vector of input directories used to
	 *  store the number of pairs (if the pairs are read from files)
	 *
	 *  @param dir_output_resample output directory used to store
	 *  the Bootstrap resampling correlation functions, with
	 *  Poisson errors; if an empty string (i.e. "" or "NULL") is
	 *  provided, no output will be stored
	 *
	 *  @param count_dd true &rarr; count the number of data-data
	 *  pairs; false &rarr; read the number of data-data pairs from
	 *  file
	 *
	 *  @param count_rr true &rarr; count the number of
	 *  random-random pairs; false &rarr; read the number of
	 *  random-random pairs from file
	 *
	 *  @param count_dr true &rarr; count the number of data-random
	 *  pairs; false &rarr; read the number of data-random pairs
	 *
	 *  @param tcount true &rarr; activate the time counter; false
	 *  &rarr; no time counter
	 *
	 *  @param estimator the estimator used to measure the two-point
	 *  correlation function
	 *
	 *  @param seed the seed for random number generation
	 *
	 *  @return none
	 */
	void measureBootstrap (const int nMocks, const std::string dir_output_pairs = par::defaultString, const std::vector<std::string> dir_input_pairs={}, const std::string dir_output_resample=par::defaultString, const bool count_dd=true, const bool count_rr=true, const bool count_dr=true, const bool tcount=true, const Estimator estimator=Estimator::_LandySzalay_, const int seed=3213) override;

      public:
    
	/**
	 *  @name Constructors/destructors
	 */
	///@{

	/**
	 *  @brief default constructor
	 *  @return object of class TwoPointCorrelation1D_monopole
	 */
	TwoPointCorrelation1D_filtered () { m_twoPType = TwoPType::_filtered_; }

	/**
	 *  @brief constructor for the filtered two-point correlation
	 *  function
	 *
	 *  @param data object of class Catalogue containing the input
	 *  catalogue
	 *
	 *  @param random of class Catalogue containing the random data
	 *  catalogue
	 *
	 *  @param binType_D1 binning type in the first dimension
	 *
	 *  @param Min_D1 minimum separation in the first dimensionused
	 *  to count the pairs
	 *
	 *  @param Max_D1 maximum separation in the first dimension used
	 *  to count the pairs
	 *
	 *  @param nbins_D1 number of bins in the first dimension
	 *
	 *  @param shift_D1 shift parameter in the first dimension,
	 *  i.e. the radial shift is binSize*shift
	 *
	 *  @param Min_D2 minimum separation in the second dimensionused
	 *  to count the pairs
	 *
	 *  @param Max_D2 maximum separation in the second dimension used
	 *  to count the pairs
	 *
	 *  @param nbins_D2 number of bins in the second dimension
	 *
	 *  @param shift_D2 shift parameter in the second dimension,
	 *  i.e. the radial shift is binSize*shift
	 *
	 *  @param angularUnits angular units
	 *  @param angularWeight angular weight function
	 *
	 *  @param compute_extra_info true &rarr; compute extra
	 *  information related to the pairs, such as the mean pair
	 *  separation and redshift
	 *  @param random_dilution_fraction fraction between the number
	 *  of objects in the diluted and original random samples, used
	 *  to improve performances in random-random pair counts
	 *
	 *  @return a pointer to an object of class TwoPointCorrelation of
	 *  a given type
	 */
	TwoPointCorrelation1D_filtered (const catalogue::Catalogue data, const catalogue::Catalogue random, const BinType binType_D1, const double Min_D1, const double Max_D1, const int nbins_D1, const double shift_D1, const double Min_D2, const double Max_D2, const int nbins_D2, const double shift_D2, const CoordinateUnits angularUnits=CoordinateUnits::_radians_, std::function<double(double)> angularWeight=nullptr, const bool compute_extra_info=false, const double random_dilution_fraction=1.)
	  : TwoPointCorrelation1D_monopole(data, random, BinType::_linear_, Min_D2, Max_D2, nbins_D2, shift_D2, angularUnits, angularWeight, compute_extra_info, random_dilution_fraction)
	  { m_twoPType = TwoPType::_filtered_;  set_parameters(binType_D1, Min_D1, Max_D1, nbins_D1, shift_D1); }

	/**
	 *  @brief constructor for the filtered two-point correlation
	 *  function
	 *  
	 *  @param data object of class Catalogue containing the input
	 *  catalogue
	 *
	 *  @param random of class Catalogue containing the random data
	 *  catalogue
	 *
	 *  @param binType_D1 binning type in the first dimension
	 *
	 *  @param Min_D1 minimum separation in the first dimension used
	 *  to count the pairs
	 *
	 *  @param Max_D1 maximum separation in the first dimension used
	 *  to count the pairs
	 *
	 *  @param binSize_D1 the bin size in the first dimension
	 *
	 *  @param shift_D1 shift parameter in the first dimension,
	 *  i.e. the radial shift is binSize*shift
	 *
	 *  @param Min_D2 minimum separation in the second dimension used
	 *  to count the pairs
	 *
	 *  @param Max_D2 maximum separation in the second dimension used
	 *  to count the pairs
	 *
	 *  @param binSize_D2 the bin size in the second dimension
	 *
	 *  @param shift_D2 shift parameter in the second dimension,
	 *  i.e. the radial shift is binSize*shift
	 *
	 *  @param angularUnits angular units
	 *  @param angularWeight angular weight function
	 *
	 *  @param compute_extra_info true &rarr; compute extra
	 *  information related to the pairs, such as the mean pair
	 *  separation and redshift
	 *  @param random_dilution_fraction fraction between the number
	 *  of objects in the diluted and original random samples, used
	 *  to improve performances in random-random pair counts
	 *
	 *  @return a pointer to an object of class TwoPointCorrelation of
	 *  a given type
	 */
	TwoPointCorrelation1D_filtered (const catalogue::Catalogue data, const catalogue::Catalogue random, const BinType binType_D1, const double Min_D1, const double Max_D1, const double binSize_D1, const double shift_D1, const double Min_D2, const double Max_D2, const double binSize_D2, const double shift_D2, const CoordinateUnits angularUnits=CoordinateUnits::_radians_, std::function<double(double)> angularWeight=nullptr, const bool compute_extra_info=false, const double random_dilution_fraction=1.)
	  : TwoPointCorrelation1D_monopole(data, random, BinType::_linear_, Min_D2, Max_D2, binSize_D2, shift_D2, angularUnits, angularWeight, compute_extra_info, random_dilution_fraction)
	  { m_twoPType = TwoPType::_filtered_; set_parameters(binType_D1, Min_D1, Max_D1, binSize_D1, shift_D1); }     

	/**
	 *  @brief default destructor
	 *  @return none
	 */
	~TwoPointCorrelation1D_filtered () = default;

	///@}

      
	/**
	 *  @name Member functions to set the binning parameters
	 */
	///@{
      
	/**
	 *  @brief set the binning parameters
	 *  @param binType binning type
	 *  @param rMin minimum separation used to count the pairs
	 *  @param rMax maximum separation used to count the pairs
	 *  @param nbins number of bins
	 *  @param shift shift parameter, i.e. the radial shift is
	 *  binSize*shift
	 *  @return none
	 */
	void set_parameters (const BinType binType, const double rMin, const double rMax, const int nbins, const double shift);

	/**
	 *  @brief set the binning parameters
	 *  @param binType binning type
	 *  @param rMin minimum separation used to count the pairs
	 *  @param rMax maximum separation used to count the pairs
	 *  @param binSize size of the bin
	 *  @param shift shift parameter, i.e. the radial shift is
	 *  binSize*shift
	 *  @return none
	 */
	void set_parameters (const BinType binType, const double rMin, const double rMax, const double binSize, const double shift);

	///@}


	/**
	 *  @name Member functions to count the number of pairs and measure the two-point correlation function
	 */
	///@{

	/**
	 *  @brief measure the filtered of the two-point correlation
	 *  function
	 *
	 *  @param errorType type of error
	 *  
	 *  @param dir_output_pairs output directory used to store the
	 *  number of pairs
	 *
	 *  @param dir_input_pairs vector of input directories used to
	 *  store the number of pairs (if the pairs are read from files)
	 *
	 *  @param dir_output_resample output directory of the
	 *  resampling correlation functions; if an empty string (i.e. ""
	 *  or "NULL") is provided, no output will be stored
	 *
	 *  @param nMocks number of resampling used for bootstrap
	 *
	 *  @param count_dd true &rarr; count the number of data-data
	 *  pairs; false &rarr; read the number of data-random pairs from
	 *  file
	 *
	 *  @param count_rr true &rarr; count the number of random-random
	 *  pairs; false &rarr; read the number of random-random pairs
	 *
	 *  @param count_dr true &rarr; count the number of data-random
	 *  pairs; false &rarr; read the number of data-random pairs
	 *
	 *  @param tcount true &rarr; activate the time counter; false
	 *  &rarr; no time counter
	 *
	 *  @param estimator the estimator used to measure the two-point
	 *  correlation function
	 *
	 *  @param seed the seed for random number generation
	 *
	 *  @return none
	 */
	void measure (const ErrorType errorType = ErrorType::_Poisson_, const std::string dir_output_pairs=par::defaultString, const std::vector<std::string> dir_input_pairs={}, const std::string dir_output_resample=par::defaultString, const int nMocks=0, const bool count_dd=true, const bool count_rr=true, const bool count_dr=true, const bool tcount=true, const Estimator estimator=Estimator::_LandySzalay_, const int seed=3213) override;
    
	///@}

    
	/**
	 *  @name Input/Output methods
	 */
	///@{

	/**
	 *  @brief read the filtered two-point correlation function
	 *  @param dir input directory
	 *  @param file input file
	 *  @return none
	 */
	void read (const std::string dir, const std::string file) override;

	/**
	 *  @brief write the filtered two-point correlation function
	 *  @param dir output directory
	 *  @param file output file
	 *  @param rank cpu index (for MPI usage)
	 *  @return none
	 */
	void write (const std::string dir=par::defaultString, const std::string file=par::defaultString, const int rank=0) const override;
      
	///@}
      
      };
    }
  }
}

#endif
