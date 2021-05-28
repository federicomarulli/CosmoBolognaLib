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
 *  @file Headers/TwoPointCorrelation_projected.h
 *
 *  @brief The class TwoPointCorrelation_projected
 *
 *  This file defines the interface of the class
 *  TwoPointCorrelation_projected, used to measure the projected
 *  two-point correlation function
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unibo.it, alfonso.veropalumbo@unibo.it
 */

#ifndef __TWOPOINTPROJ__
#define __TWOPOINTPROJ__


#include "TwoPointCorrelation2D_cartesian.h"


// ===================================================================================================


namespace cbl {
  
  namespace measure {

    namespace twopt {
  
      /**
       *  @class TwoPointCorrelation_projected
       *  TwoPointCorrelation_projected.h
       *  "Headers/TwoPointCorrelation_projected.h"
       *
       *  @brief The class TwoPointCorrelation_projected
       *
       *  This class is used to handle objects of type <EM>
       *  TwoPointCorrelation_projected </EM>. It is used to measure the
       *  projected two-point correlation function,
       *  \f$w(r_p)=2\int^{\pi_{\rm max}}_{0} d\pi' \, \xi(r_p,
       *  \pi')\f$.
       */
      class TwoPointCorrelation_projected : public TwoPointCorrelation2D_cartesian {

      protected:

	/// the upper limit of the integral
	double m_piMax_integral;

	/**
	 *  @brief return a data object with extra info
	 *  
	 *  @param rp vector containing the binned perpendicular
	 *  separations
	 *  @param ww vector containing the binned projected correlation
	 *  function
	 *  @param error vector containing the errors
	 *
	 *  @return pointer to an object of type Data
	 */
	std::shared_ptr<data::Data> data_with_extra_info (const std::vector<double> rp, const std::vector<double> ww, const std::vector<double> error) const;
      
	/**
	 *  @brief measure projected correlation function
	 *
	 *  @param rp projected separation
	 *
	 *  @param pi line of sight separation
	 *
	 *  @param xi 2d cartesian 2pcf
	 *
	 *  @param error_xi error on the 2d cartesian 2pcf
	 *
	 *  @return pointer to an object of type Data
	 */
	std::shared_ptr<data::Data> Projected (const std::vector<double> rp, const std::vector<double> pi, const std::vector<std::vector<double>> xi, const std::vector<std::vector<double>> error_xi) override;

	/**
	 *  @brief measure the projected two-point correlation function
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
	 *  @param fact factor used to compute the cell size of the
	 *  chain mesh: it is multiplied by the maximum distance
	 *  considered for the couples and can be setted by the user
	 *  to optimize the count of the couples
	 *
	 *  
	 */
	void measurePoisson (const std::string dir_output_pairs= par::defaultString, const std::vector<std::string> dir_input_pairs={}, const bool count_dd=true, const bool count_rr=true, const bool count_dr=true, const bool tcount=true, const Estimator estimator=Estimator::_LandySzalay_, const double fact=0.1) override;

	/**
	 *  @brief measure the projected two-point correlation function
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
	 *  @param fact factor used to compute the cell size of the
	 *  chain mesh: it is multiplied by the maximum distance
	 *  considered for the couples and can be setted by the user
	 *  to optimize the count of the couples
	 *
	 *  
	 */
	void measureJackknife (const std::string dir_output_pairs= par::defaultString, const std::vector<std::string> dir_input_pairs={}, const std::string dir_output_resample= par::defaultString, const bool count_dd=true, const bool count_rr=true, const bool count_dr=true, const bool tcount=true, const Estimator estimator=Estimator::_LandySzalay_, const double fact=0.1) override;

	/**
	 *  @brief measure the projected two-point correlation function
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
	 *  @param fact factor used to compute the cell size of the
	 *  chain mesh: it is multiplied by the maximum distance
	 *  considered for the couples and can be setted by the user
	 *  to optimize the count of the couples
	 *
	 *  @param seed the seed for random number generation
	 *
	 *  
	 */
	void measureBootstrap (const int nMocks, const std::string dir_output_pairs, const std::vector<std::string> dir_input_pairs={}, const std::string dir_output_resample = par::defaultString, const bool count_dd=true, const bool count_rr=true, const bool count_dr=true, const bool tcount=true, const Estimator estimator=Estimator::_LandySzalay_, const double fact=0.1, const int seed=3213) override;

	/**
	 *  @brief measure the jackknife resampling of the two-point correlation
	 *  function, &xi;(r) 
	 *
	 *  @param dd vector of data-data pairs, divided per regions
	 *
	 *  @param rr vector of random-random pairs, divided per regions
	 *
	 *  @return vector of pointers to objects of type Data
	 */
	std::vector<std::shared_ptr<data::Data>> XiJackknife (const std::vector<std::shared_ptr<pairs::Pair>> dd, const std::vector<std::shared_ptr<pairs::Pair>> rr) override;

	/**
	 *  @brief measure the jackknife resampling of the two-point correlation
	 *  function, &xi;(r)         
	 *
	 *  @param dd vector of data-data pairs, divided per regions
	 *
	 *  @param rr vector of random-random pairs, divided per
	 *  regions
	 *
	 *  @param dr vector of data-random pairs, divided per regions
	 * 
	 *  @return a vector of pointers to objects of type Data
	 */
	std::vector<std::shared_ptr<data::Data>> XiJackknife (const std::vector<std::shared_ptr<pairs::Pair>> dd, const std::vector<std::shared_ptr<pairs::Pair>> rr, const std::vector<std::shared_ptr<pairs::Pair>> dr) override;

	/**
	 *  @brief measure the bootstrap resampling of the two-point correlation
	 *  function, &xi;(r)  
	 *
	 *  @param nMocks number of bootstrap resampling
	 *
	 *  @param dd vector of data-data pairs, divided per regions
	 *
	 *  @param rr vector of random-random pairs, divided per
	 *  regions
	 *
	 *  @param seed the seed for random number generation
	 *
	 *  @return a vector of pointers to objects of type Data
	 */
	std::vector<std::shared_ptr<data::Data>> XiBootstrap (const int nMocks, const std::vector<std::shared_ptr<pairs::Pair>> dd, const std::vector<std::shared_ptr<pairs::Pair>> rr, const int seed=3213) override;

	/**
	 *  @brief measure the bootstrap resampling of the two-point correlation
	 *  function, &xi;(r)  
	 *
	 *  @param nMocks number of bootstrap resampling
	 *
	 *  @param dd vector of data-data pairs, divided per regions
	 *
	 *  @param rr vector of random-random pairs, divided per
	 *  regions
	 *
	 *  @param dr vector of data-random pairs, divided per regions
	 *
	 *  @param seed the seed for random number generation
	 *
	 *  @return a vector of pointers to objects of type Data
	 */
	std::vector<std::shared_ptr<data::Data>> XiBootstrap (const int nMocks, const std::vector<std::shared_ptr<pairs::Pair>> dd, const std::vector<std::shared_ptr<pairs::Pair>> rr, const std::vector<std::shared_ptr<pairs::Pair>> dr, const int seed=3213) override;

      
      public:
    
	/**
	 *  @name Constructors/destructors
	 */
	///@{

	/**
	 *  @brief default constructor
	 *  _projected
	 */
	TwoPointCorrelation_projected () { m_twoPType = TwoPType::_projected_; m_piMax_integral=50.; }

	/**
	 *  @brief constructor
	 *  @param data object of class Catalogue containing the input
	 *  catalogue
	 *  @param random of class Catalogue containing the random data
	 *  catalogue
	 *  @param binType_rp binning type
	 *  @param rpMin minimum perpendicular separation used to count
	 *  the pairs
	 *  @param rpMax maximum perpendicular separation used to count
	 *  the pairs
	 *  @param nbins_rp number of bins in the perpendicular
	 *  separation
	 *  @param shift_rp shift parameter in the perpendicular
	 *  separation, i.e. the radial shift is binSize*shift
	 *  @param piMin minimum parallel separation used to count
	 *  the pairs
	 *  @param piMax maximum parallel separation used to count
	 *  the pairs
	 *  @param nbins_pi number of bins in the parallel
	 *  separation
	 *  @param shift_pi shift parameter in the parallel
	 *  separation, i.e. the radial shift is binSize*shift
	 *  @param piMax_integral upper limits of the integral
	 *  @param angularUnits angular units
	 *  @param angularWeight angular weight function
	 *  @param compute_extra_info true &rarr; compute extra
	 *  information related to the pairs, such as the mean pair
	 *  separation and redshift
	 *  @param random_dilution_fraction fraction between the number
	 *  of objects in the diluted and original random samples, used
	 *  to improve performances in random-random pair counts
	 *  _projected
	 */
	TwoPointCorrelation_projected (catalogue::Catalogue data, catalogue::Catalogue random, const BinType binType_rp, const double rpMin, const double rpMax, const int nbins_rp, const double shift_rp, const double piMin, const double piMax, const int nbins_pi, const double shift_pi, const double piMax_integral, const CoordinateUnits angularUnits=CoordinateUnits::_radians_, std::function<double(double)> angularWeight=nullptr, const bool compute_extra_info=false, const double random_dilution_fraction=1.)
	  : TwoPointCorrelation2D_cartesian(data, random, binType_rp, rpMin, rpMax, nbins_rp, shift_rp, BinType::_linear_, piMin, piMax, nbins_pi, shift_pi, angularUnits, angularWeight, compute_extra_info, random_dilution_fraction)
	  { m_twoPType = TwoPType::_projected_; m_piMax_integral = piMax_integral; }
      
	/**
	 *  @brief constructor
	 *  @param data object of class Catalogue containing the input
	 *  catalogue
	 *  @param random of class Catalogue containing the random data
	 *  catalogue
	 *  @param binType_rp binning type
	 *  @param rpMin minimum perpendicular separation used to count
	 *  the pairs
	 *  @param rpMax maximum perpendicular separation used to count
	 *  the pairs
	 *  @param binSize_rp bin size in the perpendicular separation
	 *  @param shift_rp shift parameter in the perpendicular
	 *  separation, i.e. the radial shift is binSize*shift
	 *  @param piMin minimum parallel separation used to count
	 *  the pairs
	 *  @param piMax maximum parallel separation used to count
	 *  the pairs
	 *  @param binSize_pi bin size in the parallel separation
	 *  @param shift_pi shift parameter in the parallel
	 *  separation, i.e. the radial shift is binSize*shift
	 *  @param piMax_integral upper limits of the integral
	 *  @param angularUnits angular units
	 *  @param angularWeight angular weight function
	 *  @param compute_extra_info true &rarr; compute extra
	 *  information related to the pairs, such as the mean pair
	 *  separation and redshift
	 *  @param random_dilution_fraction fraction between the number
	 *  of objects in the diluted and original random samples, used
	 *  to improve performances in random-random pair counts
	 *  2D_projected
	 */
	TwoPointCorrelation_projected (catalogue::Catalogue data, catalogue::Catalogue random, const BinType binType_rp, const double rpMin, const double rpMax, const double binSize_rp, const double shift_rp, const double piMin, const double piMax, const double binSize_pi, const double shift_pi, const double piMax_integral, const CoordinateUnits angularUnits=CoordinateUnits::_radians_, std::function<double(double)> angularWeight=nullptr, const bool compute_extra_info=false, const double random_dilution_fraction=1.)
	  : TwoPointCorrelation2D_cartesian(data, random, binType_rp, rpMin, rpMax, binSize_rp, shift_rp, BinType::_linear_, piMin, piMax, binSize_pi, shift_pi, angularUnits, angularWeight, compute_extra_info, random_dilution_fraction)
	  { m_twoPType = TwoPType::_projected_; m_piMax_integral = piMax_integral; }
      
	/**
	 *  @brief default destructor
	 *  
	 */
	~TwoPointCorrelation_projected () = default;

	///@}

	/**
	 *  @brief get the y coordinates
	 *  @return the y coordinates
	 */
	std::vector<double> yy () const 
	  { cbl::ErrorCBL("", "yy", "TwoPointCorrelation_projected.h"); std::vector<double> vv; return vv; }

	/**
	 *  @brief get the the binned correlation function 
	 *  @return the binned correlation function 
	 */
	std::vector<double> xi1D () const
	  { std::vector<double> vv; m_dataset->get_data(vv); return vv; } 

	/**
	 *  @brief get the error on the binned correlation function
	 *  function
	 *  @return the error on the binned correlation function
	 *  function
	 */
	virtual std::vector<double> error1D () const
	{ std::vector<double> vv; m_dataset->get_error(vv); return vv; }

	/**
	 *  @brief get the the binned correlation function 
	 *  @return the binned correlation function 
	 */
	std::vector<std::vector<double>> xi2D () const 
	  { cbl::ErrorCBL("", "xi2D", "TwoPointCorrelation_projected.h"); std::vector<std::vector<double>> vv; return vv; }

	/**
	 *  @brief get the error on the binned correlation function
	 *  function
	 *  @return the error on the binned correlation function
	 *  function
	 */
	std::vector<std::vector<double>> error2D () const 
	  { cbl::ErrorCBL("", "error2D", "TwoPointCorrelation_projected.h"); std::vector<std::vector<double>> vv; return vv; }
      

	/**
	 *  @name Member functions to count the number of pairs and measure the two-point correlation function
	 */
	///@{

	/**
	 *  @brief measure the projected two-point correlation function
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
	 *  @param fact factor used to compute the cell size of the
	 *  chain mesh: it is multiplied by the maximum distance
	 *  considered for the couples and can be setted by the user
	 *  to optimize the count of the couples
	 *
	 *  @param seed the seed for random number generation
	 *
	 *  
	 */
	void measure (const ErrorType errorType=ErrorType::_Poisson_, const std::string dir_output_pairs=par::defaultString, const std::vector<std::string> dir_input_pairs={}, const std::string dir_output_resample=par::defaultString, const int nMocks=0, const bool count_dd=true, const bool count_rr=true, const bool count_dr=true, const bool tcount=true, const Estimator estimator=Estimator::_LandySzalay_, const double fact=0.1, const int seed=3213) override;

	///@}

    
	/**
	 *  @name Input/Output methods
	 */
	///@{

	/**
	 *  @brief read the projected two-point correlation function
	 *  @param dir input directory
	 *  @param file input file
	 *  
	 */
	void read (const std::string dir, const std::string file) override;

	/**
	 *  @brief write the projected two-point correlation function
	 *  @param dir output directory
	 *  @param file output file
	 *  @param rank cpu index (for MPI usage)
	 *  
	 */
	void write (const std::string dir=par::defaultString, const std::string file=par::defaultString, const int rank=0) const override;
    
	///@}


	/**
	 *  @name Member functions to compute, read and write the covariance matrix
	 */
	///@{ 

	/**
	 *  @brief read the measured covariance matrix
	 *  @param dir input directory
	 *  @param file input file
	 *  
	 */
	void read_covariance (const std::string dir, const std::string file) override;

	/**
	 *  @brief write the measured two-point correlation
	 *  @param dir output directory
	 *  @param file output file
	 *  
	 */
	void write_covariance (const std::string dir, const std::string file) const override;

	/**
	 *  @brief compute the covariance matrix
	 *  @param xi vector containing the measure correlation
	 *  functions used to compute the covariance matrix
	 *  @param JK true &rarr; compute the jackknife covariance
	 *  matrix; false compute the standard covariance matrix
	 *  
	 */
	void compute_covariance (const std::vector<std::shared_ptr<data::Data>> xi, const bool JK) override;

	/**
	 *  @brief compute the covariance matrix
	 *  @param file vector containing the input files with the
	 *  measured correlation functions used to compute the
	 *  covariance matrix
	 *  @param JK true &rarr; compute the jackknife covariance
	 *  matrix; false compute the standard covariance matrix
	 *  
	 */
	void compute_covariance (const std::vector<std::string> file, const bool JK) override;

	///@}
      
      };
    }
  }
}

#endif
