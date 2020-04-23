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
 *  @file Headers/TwoPointCorrelation2D_polar.h
 *
 *  @brief The class TwoPointCorrelation2D_polar
 *
 *  This file defines the interface of the class
 *  TwoPointCorrelation2D_polar, used to measure the 2D two-point
 *  correlation function in polar coordinates
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unibo.it, alfonso.veropalumbo@unibo.it
 */

#ifndef __TWOPOINT2DPOL__
#define __TWOPOINT2DPOL__


#include "TwoPointCorrelation2D.h"


// ===================================================================================================


namespace cbl {
  
  namespace measure {
    
    namespace twopt {
    
      /**
       *  @class TwoPointCorrelation2D_polar TwoPointCorrelation2D_polar.h
       *  "Headers/TwoPointCorrelation2D_polar.h"
       *
       *  @brief The class TwoPointCorrelation2D_polar
       *
       *  This class is used to handle objects of type <EM>
       *  TwoPointCorrelation2D_polar </EM>. It is used to measure the
       *  2D two-point correlation function in polar coordinates,
       *  \f$\xi(r,\mu)\f$, that is as a function of absolute
       *  separation, \f$r=\sqrt{r_p^2+\pi^2}\f$, and the cosine of the
       *  angle between the separation vector and the line of sight,
       *  \f$\mu\equiv\cos\theta=s_\parallel/s\f$.
       */
      class TwoPointCorrelation2D_polar : public TwoPointCorrelation2D {

      protected:

	/**
	 *  @brief measure the 2D two-point correlation function in
	 *  polar coordinates, with Poisson errors
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
	 *  @brief measure the 2D two-point correlation function in
	 *  polar coordinates, estimating the covariance with Jackknife
	 *  resampling
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
	void measureJackknife (const std::string dir_output_pairs, const std::vector<std::string> dir_input_pairs={}, const std::string dir_output_resample = "NULL", const bool count_dd=true, const bool count_rr=true, const bool count_dr=true, const bool tcount=true, const Estimator estimator=Estimator::_LandySzalay_) override;

	/**
	 *  @brief measure the 2D two-point correlation function in
	 *  polar coordinates, estimating the covariance with Bootstrap
	 *  resampling
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
	void measureBootstrap (const int nMocks, const std::string dir_output_pairs, const std::vector<std::string> dir_input_pairs={}, const std::string dir_output_resample="NULL", const bool count_dd=true, const bool count_rr=true, const bool count_dr=true, const bool tcount=true, const Estimator estimator=Estimator::_LandySzalay_, const int seed=3213) override;

      
      public:

	/**
	 *  @name Constructors/destructors
	 */
	///@{

	/**
	 *  @brief default constructor
	 *  @return object of class TwoPointCorrelation2D_polar
	 */
	TwoPointCorrelation2D_polar () { m_twoPType = TwoPType::_2D_polar_; }

	/**
	 *  @brief constructor
	 *  @param data object of class Catalogue containing the input
	 *  catalogue
	 *  @param random of class Catalogue containing the random data
	 *  catalogue
	 *  @param binType_rad binning type in absolute separations
	 *  @param rMin minimum absolute separation used to count
	 *  the pairs
	 *  @param rMax maximum absolute separation used to count
	 *  the pairs
	 *  @param nbins_rad number of bins in the absolute
	 *  separation
	 *  @param shift_rad shift parameter in the absolute
	 *  separation, i.e. the radial shift is binSize*shift
	 *  @param binType_mu binning type in angular separations
	 *  @param muMin minimum angular used to count the pairs
	 *  @param muMax maximum angular used to count the pairs
	 *  @param nbins_mu number of bins in the angular
	 *  separation
	 *  @param shift_mu shift parameter in the angular
	 *  separation, i.e. the radial shift is binSize*shift
	 *  @param angularUnits angular units
	 *  @param angularWeight angular weight function
	 *  @param compute_extra_info true &rarr; compute extra
	 *  information related to the pairs, such as the mean pair
	 *  separation and redshift
	 *  @param random_dilution_fraction fraction between the number
	 *  of objects in the diluted and original random samples, used
	 *  to improve performances in random-random pair counts
	 *  @return object of class TwoPointCorrelation2D_polar
	 */
	TwoPointCorrelation2D_polar (catalogue::Catalogue data, catalogue::Catalogue random, const BinType binType_rad, const double rMin, const double rMax, const int nbins_rad, const double shift_rad, const BinType binType_mu, const double muMin, const double muMax, const int nbins_mu, const double shift_mu, const CoordinateUnits angularUnits=CoordinateUnits::_radians_, std::function<double(double)> angularWeight=nullptr, const bool compute_extra_info=false, const double random_dilution_fraction=1.)
	  : TwoPointCorrelation2D(data, random, compute_extra_info, random_dilution_fraction) { m_twoPType = TwoPType::_2D_polar_; set_parameters(binType_rad, rMin, rMax, nbins_rad, shift_rad, binType_mu, muMin, muMax, nbins_mu, shift_mu, angularUnits, angularWeight, compute_extra_info); }

	/**
	 *  @brief constructor
	 *  @param data object of class Catalogue containing the input
	 *  catalogue
	 *  @param random of class Catalogue containing the random data
	 *  catalogue
	 *  @param binType_rad binning type in absolute separations
	 *  @param rMin minimum absolute separation used to count
	 *  the pairs
	 *  @param rMax maximum absolute separation used to count
	 *  the pairs
	 *  @param binSize_rad bin size in the absolute separation
	 *  @param shift_rad shift parameter in the absolute
	 *  separation, i.e. the radial shift is binSize*shift
	 *  @param binType_mu binning type in angular separations
	 *  @param muMin minimum angular separation used to count
	 *  the pairs
	 *  @param muMax maximum angular separation used to count
	 *  the pairs
	 *  @param binSize_mu bin size in the angular separation
	 *  @param shift_mu shift parameter in the angular
	 *  separation, i.e. the radial shift is binSize*shift
	 *  @param angularUnits angular units
	 *  @param angularWeight angular weight function
	 *  @param compute_extra_info true &rarr; compute extra
	 *  information related to the pairs, such as the mean pair
	 *  separation and redshift
	 *  @param random_dilution_fraction fraction between the number
	 *  of objects in the diluted and original random samples, used
	 *  to improve performances in random-random pair counts
	 *  @return object of class TwoPointCorrelation2D_polar
	 */
	TwoPointCorrelation2D_polar (catalogue::Catalogue data, catalogue::Catalogue random, const BinType binType_rad, const double rMin, const double rMax, const double binSize_rad, const double shift_rad, const BinType binType_mu, const double muMin, const double muMax, const double binSize_mu, const double shift_mu, const CoordinateUnits angularUnits=CoordinateUnits::_radians_, std::function<double(double)> angularWeight=nullptr, const bool compute_extra_info=false, const double random_dilution_fraction=1.)
	  : TwoPointCorrelation2D(data, random, compute_extra_info, random_dilution_fraction) { m_twoPType = TwoPType::_2D_polar_; set_parameters(binType_rad, rMin, rMax, binSize_rad, shift_rad, binType_mu, muMin, muMax, binSize_mu, shift_mu, angularUnits, angularWeight, compute_extra_info); }

	/**
	 *  @brief default destructor
	 *  @return none
	 */
	~TwoPointCorrelation2D_polar () = default;

	///@}


	/**
	 *  @name Member functions to set the binning parameters
	 */
	///@{

	/**
	 *  @brief set the binning parameters
	 *  @param binType_rad binning type in absolute separations
	 *  @param rMin minimum absolute separation used to count
	 *  the pairs
	 *  @param rMax maximum absolute separation used to count
	 *  the pairs
	 *  @param nbins_rad number of bins in the absolute
	 *  separation
	 *  @param shift_rad shift parameter in the absolute
	 *  separation, i.e. the radial shift is binSize*shift
	 *  @param binType_mu binning type in angular separations
	 *  @param muMin minimum angular used to count the pairs
	 *  @param muMax maximum angular used to count the pairs
	 *  @param nbins_mu number of bins in the angular
	 *  separation
	 *  @param shift_mu shift parameter in the angular
	 *  separation, i.e. the radial shift is binSize*shift
	 *  @param angularUnits angular units
	 *  @param angularWeight angular weight function
	 *  @param compute_extra_info true &rarr; compute extra
	 *  information related to the pairs, such as the mean pair
	 *  separation and redshift
	 *  @return none
	 */
	void set_parameters (const BinType binType_rad, const double rMin, const double rMax, const int nbins_rad, const double shift_rad, const BinType binType_mu, const double muMin, const double muMax, const int nbins_mu, const double shift_mu, const CoordinateUnits angularUnits=CoordinateUnits::_radians_, std::function<double(double)> angularWeight=nullptr, const bool compute_extra_info=false);

	/**
	 *  @brief set the binning parameters
	 *  @param binType_rad binning type in absolute separations
	 *  @param rMin minimum absolute separation used to count
	 *  the pairs
	 *  @param rMax maximum absolute separation used to count
	 *  the pairs
	 *  @param binSize_rad bin size in the absolute separation
	 *  @param shift_rad shift parameter in the absolute
	 *  separation, i.e. the radial shift is binSize*shift
	 *  @param binType_mu binning type in angular separations
	 *  @param muMin minimum angular separation used to count
	 *  the pairs
	 *  @param muMax maximum angular separation used to count
	 *  the pairs
	 *  @param binSize_mu bin size in the angular separation
	 *  @param shift_mu shift parameter in the angular
	 *  separation, i.e. the radial shift is binSize*shift
	 *  @param angularUnits angular units
	 *  @param angularWeight angular weight function
	 *  @param compute_extra_info true &rarr; compute extra
	 *  information related to the pairs, such as the mean pair
	 *  separation and redshift
	 *  @return none
	 */
	void set_parameters (const BinType binType_rad, const double rMin, const double rMax, const double binSize_rad, const double shift_rad, const BinType binType_mu, const double muMin, const double muMax, const double binSize_mu, const double shift_mu, const CoordinateUnits angularUnits=CoordinateUnits::_radians_, std::function<double(double)> angularWeight=nullptr, const bool compute_extra_info=false);

	///@}


	/**
	 *  @name Methods to set the binning parameters
	 */
	///@{

	/**
	 *  @brief measure the 2D two-point correlation function in
	 *  polar coordinates
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
	 *  resampling correlation functions; if an empty string
	 *  (i.e. "" or "NULL") is provided, no output will be stored
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
	void measure (const ErrorType errorType=ErrorType::_Poisson_, const std::string dir_output_pairs=par::defaultString, const std::vector<std::string> dir_input_pairs={}, const std::string dir_output_resample=par::defaultString, const int nMocks=0, const bool count_dd=true, const bool count_rr=true, const bool count_dr=true, const bool tcount=true, const Estimator estimator=Estimator::_LandySzalay_, const int seed=3213) override;

	///@}


	/**
	 *  @name Input/Output methods
	 */
	///@{

	/**
	 *  @brief read the 2D two-point correlation function
	 *  @param dir input directory
	 *  @param file input file
	 *  @return none
	 */
	void read (const std::string dir, const std::string file) override;     

	/**
	 *  @brief write the 2D two-point correlation function
	 *  @param dir output directory
	 *  @param file output file
	 *  @param full false &rarr; simply store the data; true &rarr;
	 *  duplicate the data in the other three quadrands (usefull
	 *  e.g. when storing the 2D correlation function)
	 *  @param rank cpu index (for MPI usage)
	 *  @return none
	 */
	void write (const std::string dir, const std::string file, const bool full, const int rank=0) const override;

	/**
	 *  @brief write the 2D two-point correlation function
	 *  @param dir output directory
	 *  @param file output file
	 *  @param rank cpu index (for MPI usage)
	 *  @return none
	 */
	void write (const std::string dir=par::defaultString, const std::string file=par::defaultString, const int rank=0) const override
	{ write(dir, file, true, rank); }
      
	///@}

      
	/**
	 *  @name Member functions to compute, read and write the covariance matrix (customised in all the derived classes)
	 */
	///@{ 

	/**
	 *  @brief read the measured covariance matrix
	 *  @param dir input directory
	 *  @param file input file
	 *  @return none
	 */
	void read_covariance (const std::string dir, const std::string file) override
	{ (void)dir; (void)file; ErrorCBL("", "read_covariance", "TwoPointCorrelation2D_polar.h", glob::ExitCode::_workInProgress_); }

	/**
	 *  @brief write the measured two-point correlation
	 *  @param dir output directory
	 *  @param file output file
	 *  @return none
	 */
	void write_covariance (const std::string dir, const std::string file) const override
	{ (void)dir; (void)file; ErrorCBL("", "write_covariance", "TwoPointCorrelation2D_polar.h", glob::ExitCode::_workInProgress_); }
      
	/**
	 *  @brief compute the covariance matrix
	 *  @param xi vector containing the measure correlation
	 *  functions used to compute the covariance matrix
	 *  @param JK true &rarr; compute the jackknife covariance
	 *  matrix; false compute the standard covariance matrix
	 *  @return none
	 */
	void compute_covariance (const std::vector<std::shared_ptr<data::Data>> xi, const bool JK) override
	{ (void)xi; (void)JK; ErrorCBL("", "compute_covariance", "TwoPointCorrelation2D_polar.h", glob::ExitCode::_workInProgress_); }
      
	/**
	 *  @brief compute the covariance matrix
	 *  @param file vector containing the input files with the
	 *  measured correlation functions used to compute the
	 *  covariance matrix
	 *  @param JK true &rarr; compute the jackknife covariance
	 *  matrix; false compute the standard covariance matrix
	 *  @return none
	 */
	void compute_covariance (const std::vector<std::string> file, const bool JK) override
	{ (void)file; (void)JK; ErrorCBL("", "compute_covariance", "TwoPointCorrelation2D_cartesian.h", glob::ExitCode::_workInProgress_); }
      
	///@} 

      
      };
    }
  }
}

#endif
