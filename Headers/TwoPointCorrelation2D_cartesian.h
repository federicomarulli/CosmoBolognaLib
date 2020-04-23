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
 *  @file Headers/TwoPointCorrelation2D_cartesian.h
 *
 *  @brief The class TwoPointCorrelation2D_cartesian
 *
 *  This file defines the interface of the class
 *  TwoPointCorrelation2D_cartesian, used to measure the 2D two-point
 *  correlation function in cartesian coordinates
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unibo.it, alfonso.veropalumbo@unibo.it
 */

#ifndef __TWOPOINT2DCART__
#define __TWOPOINT2DCART__


#include "TwoPointCorrelation2D.h"


// ===================================================================================================


namespace cbl {
  
  namespace measure {

    namespace twopt {
    
      /**
       *  @class TwoPointCorrelation2D_cartesian
       *  TwoPointCorrelation2D_cartesian.h
       *  "Headers/TwoPointCorrelation2D_cartesian.h"
       *
       *  @brief The class TwoPointCorrelation2D_cartesian
       *
       *  This class is used to handle objects of type <EM>
       *  TwoPointCorrelation2D_cartesian </EM>. It is used to measure
       *  the 2D two-point correlation function in cartesian
       *  coordinates, \f$\xi(r_p,\pi)\f$, that is as a function of
       *  perpendicular, \f$r_p\f$, and parallel, \f$\pi\f$,
       *  line-of-sight separations.
       */
      class TwoPointCorrelation2D_cartesian : public TwoPointCorrelation2D {

      protected:
      
	/**
	 *  @brief measure the 2D two-point correlation function in
	 *  Cartesian coordinates, with Poisson errors
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
	 *  Cartesian coordinates, estimating the covariance with
	 *  Jackknife resampling
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
	void measureJackknife (const std::string dir_output_pairs=par::defaultString, const std::vector<std::string> dir_input_pairs={}, const std::string dir_output_resample=par::defaultString, const bool count_dd=true, const bool count_rr=true, const bool count_dr=true, const bool tcount=true, const Estimator estimator=Estimator::_LandySzalay_) override;

	/**
	 *  @brief measure the 2D two-point correlation function in
	 *  Cartesian coordinates, estimating the covariance with
	 *  Bootstrap resampling
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
	void measureBootstrap (const int nMocks, const std::string dir_output_pairs=par::defaultString, const std::vector<std::string> dir_input_pairs={}, const std::string dir_output_resample=par::defaultString, const bool count_dd=true, const bool count_rr=true, const bool count_dr=true, const bool tcount=true, const Estimator estimator=Estimator::_LandySzalay_, const int seed=3213) override;

      public:

	/**
	 *  @name Constructors/destructors
	 */
	///@{

	/**
	 *  @brief default constructor
	 *  @return object of class TwoPointCorrelation2D_cartesian
	 */
	TwoPointCorrelation2D_cartesian () { m_twoPType = TwoPType::_2D_Cartesian_; }

	/**
	 *  @brief constructor
	 *  @param data object of class Catalogue containing the input
	 *  catalogue
	 *  @param random of class Catalogue containing the random data
	 *  catalogue
	 *  @param binType_rp binning type in perpendicular separations
	 *  @param rpMin minimum perpendicular separation used to count
	 *  the pairs
	 *  @param rpMax maximum perpendicular separation used to count
	 *  the pairs
	 *  @param nbins_rp number of bins in the perpendicular
	 *  separation
	 *  @param shift_rp shift parameter in the perpendicular
	 *  separation, i.e. the radial shift is binSize*shift
	 *  @param binType_pi binning type in parallel separations
	 *  @param piMin minimum parallel separation used to count
	 *  the pairs
	 *  @param piMax maximum parallel separation used to count
	 *  the pairs
	 *  @param nbins_pi number of bins in the parallel
	 *  separation
	 *  @param shift_pi shift parameter in the parallel
	 *  separation, i.e. the radial shift is binSize*shift
	 *  @param angularUnits angular units
	 *  @param angularWeight angular weight function
	 *  @param compute_extra_info true &rarr; compute extra
	 *  information related to the pairs, such as the mean pair
	 *  separation and redshift
	 *  @param random_dilution_fraction fraction between the number
	 *  of objects in the diluted and original random samples, used
	 *  to improve performances in random-random pair counts
	 *  @return object of class TwoPointCorrelation2D_cartesian
	 */
	TwoPointCorrelation2D_cartesian (catalogue::Catalogue data, catalogue::Catalogue random, const BinType binType_rp, const double rpMin, const double rpMax, const int nbins_rp, const double shift_rp, const BinType binType_pi, const double piMin, const double piMax, const int nbins_pi, const double shift_pi, const CoordinateUnits angularUnits=CoordinateUnits::_radians_, std::function<double(double)> angularWeight=nullptr, const bool compute_extra_info=false, const double random_dilution_fraction=1.)
	  : TwoPointCorrelation2D(data, random, compute_extra_info, random_dilution_fraction)
	  { m_twoPType = TwoPType::_2D_Cartesian_; set_parameters(binType_rp, rpMin, rpMax, nbins_rp, shift_rp, binType_pi, piMin, piMax, nbins_pi, shift_pi, angularUnits, angularWeight, compute_extra_info); }

	/**
	 *  @brief constructor
	 *  @param data object of class Catalogue containing the input
	 *  catalogue
	 *  @param random of class Catalogue containing the random data
	 *  catalogue
	 *  @param binType_rp binning type in perpendicular separations
	 *  @param rpMin minimum perpendicular separation used to count
	 *  the pairs
	 *  @param rpMax maximum perpendicular separation used to count
	 *  the pairs
	 *  @param binSize_rp bin size in the perpendicular separation
	 *  @param shift_rp shift parameter in the perpendicular
	 *  separation, i.e. the radial shift is binSize*shift
	 *  @param binType_pi binning type in parallel separations
	 *  @param piMin minimum parallel separation used to count
	 *  the pairs
	 *  @param piMax maximum parallel separation used to count
	 *  the pairs
	 *  @param binSize_pi bin size in the parallel separation
	 *  @param shift_pi shift parameter in the parallel
	 *  separation, i.e. the radial shift is binSize*shift
	 *  @param angularUnits angular units
	 *  @param angularWeight angular weight function
	 *  @param compute_extra_info true &rarr; compute extra
	 *  information related to the pairs, such as the mean pair
	 *  separation and redshift
	 *  @param random_dilution_fraction fraction between the number
	 *  of objects in the diluted and original random samples, used
	 *  to improve performances in random-random pair counts
	 *  @return object of class TwoPointCorrelation2D_cartesian
	 */
	TwoPointCorrelation2D_cartesian (catalogue::Catalogue data, catalogue::Catalogue random, const BinType binType_rp, const double rpMin, const double rpMax, const double binSize_rp, const double shift_rp, const BinType binType_pi, const double piMin, const double piMax, const double binSize_pi, const double shift_pi, const CoordinateUnits angularUnits=CoordinateUnits::_radians_, std::function<double(double)> angularWeight=nullptr, const bool compute_extra_info=false, const double random_dilution_fraction=1.)
	  : TwoPointCorrelation2D(data, random, compute_extra_info, random_dilution_fraction) { m_twoPType = TwoPType::_2D_Cartesian_; set_parameters(binType_rp, rpMin, rpMax, binSize_rp, shift_rp, binType_pi, piMin, piMax, binSize_pi, shift_pi, angularUnits, angularWeight, compute_extra_info); }

	/**
	 *  @brief default destructor
	 *  @return none
	 */
	~TwoPointCorrelation2D_cartesian () = default;

	///@}


	/**
	 *  @name Member functions to set the binning parameters
	 */
	///@{

	/**
	 *  @brief set the binning parameters
	 *  @param binType_rp binning type in perpendicular separations
	 *  @param rpMin minimum perpendicular separation used to count
	 *  the pairs
	 *  @param rpMax maximum perpendicular separation used to count
	 *  the pairs
	 *  @param nbins_rp number of bins in the perpendicular
	 *  separation
	 *  @param shift_rp shift parameter in the perpendicular
	 *  separation, i.e. the radial shift is binSize*shift
	 *  @param binType_pi binning type in parallel separations
	 *  @param piMin minimum parallel separation used to count
	 *  the pairs
	 *  @param piMax maximum parallel separation used to count
	 *  the pairs
	 *  @param nbins_pi number of bins in the parallel
	 *  separation
	 *  @param shift_pi shift parameter in the parallel
	 *  separation, i.e. the radial shift is binSize*shift
	 *  @param angularUnits angular units
	 *  @param angularWeight angular weight function
	 *  @param compute_extra_info true &rarr; compute extra
	 *  information related to the pairs, such as the mean pair
	 *  separation and redshift
	 *  @return object of class TwoPointCorrelation2D_cartesian
	 */
	void set_parameters (const BinType binType_rp, const double rpMin, const double rpMax, const int nbins_rp, const double shift_rp, const BinType binType_pi, const double piMin, const double piMax, const int nbins_pi, const double shift_pi, const CoordinateUnits angularUnits=CoordinateUnits::_radians_, std::function<double(double)> angularWeight=nullptr, const bool compute_extra_info=false);

	/**
	 *  @brief set the binning parameters
	 *  @param binType_rp binning type in perpendicular separations
	 *  @param rpMin minimum perpendicular separation used to count
	 *  the pairs
	 *  @param rpMax maximum perpendicular separation used to count
	 *  the pairs
	 *  @param binSize_rp bin size in the perpendicular separation
	 *  @param shift_rp shift parameter in the perpendicular
	 *  separation, i.e. the radial shift is binSize*shift
	 *  @param binType_pi binning type in parallel separations
	 *  @param piMin minimum parallel separation used to count
	 *  the pairs
	 *  @param piMax maximum parallel separation used to count
	 *  the pairs
	 *  @param binSize_pi bin size in the parallel separation
	 *  @param shift_pi shift parameter in the parallel
	 *  separation, i.e. the radial shift is binSize*shift
	 *  @param angularUnits angular units
	 *  @param angularWeight angular weight function
	 *  @param compute_extra_info true &rarr; compute extra
	 *  information related to the pairs, such as the mean pair
	 *  separation and redshift
	 *  @return none
	 */
	void set_parameters (const BinType binType_rp, const double rpMin, const double rpMax, const double binSize_rp, const double shift_rp, const BinType binType_pi, const double piMin, const double piMax, const double binSize_pi, const double shift_pi, const CoordinateUnits angularUnits=CoordinateUnits::_radians_, std::function<double(double)> angularWeight=nullptr, const bool compute_extra_info=false);

	///@}


	/**
	 *  @name Methods to set the binning parameters
	 */
	///@{

	/**
	 *  @brief measure the 2D two-point correlation function in
	 *  Cartesian coordinates
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
	{ (void)dir; (void)file; ErrorCBL("", "read_covariance", "TwoPointCorrelation2D_cartesian.h", glob::ExitCode::_workInProgress_); }

	/**
	 *  @brief write the measured two-point correlation
	 *  @param dir output directory
	 *  @param file output file
	 *  @return none
	 */
	void write_covariance (const std::string dir, const std::string file) const override
	{ (void)dir; (void)file; ErrorCBL("", "write_covariance", "TwoPointCorrelation2D_cartesian.h", glob::ExitCode::_workInProgress_); }
      
	/**
	 *  @brief compute the covariance matrix
	 *  @param xi vector containing the measure correlation
	 *  functions used to compute the covariance matrix
	 *  @param JK true &rarr; compute the jackknife covariance
	 *  matrix; false compute the standard covariance matrix
	 *  @return none
	 */
	void compute_covariance (const std::vector<std::shared_ptr<data::Data>> xi, const bool JK) override
	{ (void)xi; (void)JK; ErrorCBL("", "compute_covariance", "TwoPointCorrelation2D_cartesian.h", glob::ExitCode::_workInProgress_); }
      
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
