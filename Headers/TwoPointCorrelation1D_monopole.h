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
 *  @file Headers/TwoPointCorrelation1D_monopole.h
 *
 *  @brief The class TwoPointCorrelation1D_monopole
 *
 *  This file defines the interface of the class
 *  TwoPointCorrelation1D_monopole, used to measure the monopole of the
 *  two-point correlation function
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unibo.it, alfonso.veropalumbo@unibo.it
 */

#ifndef __TWOPOINTMON__
#define __TWOPOINTMON__


#include "TwoPointCorrelation1D.h"


// ===================================================================================================


namespace cbl {

  namespace measure {
  
    namespace twopt {
  
      /**
       *  @class TwoPointCorrelation1D_monopole
       *  TwoPointCorrelation1D_monopole.h
       *  "Headers/TwoPointCorrelation1D_monopole.h"
       *
       *  @brief The class TwoPointCorrelation1D_monopole
       *
       *  This class is used to handle objects of type <EM>
       *  TwoPointCorrelation1D_monopole </EM>. It is used to measure
       *  the monopole of the two-point correlation function,
       *  \f$\xi(r)\f$, defined as \f$dP_{12}=n^2[1+\xi(r)]dV_1dV_2\f$,
       *  where \f$n\f$ is the average number density, and \f$dP_{12}\f$
       *  is the probability of finding a pair with one object in the
       *  volume \f$dV_1\f$ and the other one in the volume \f$dV_2\f$,
       *  separated by a comoving distance r.
       * 
       *  This class estimates the monopole with the so-called <EM>
       *  direct </EM> method (e.g. http://arxiv.org/abs/1105.2037 ,
       *  appendix E). This method might be preferable to the <EM>
       *  integrated </EM> method, since it is less affected by
       *  numerical issues (related to the binning and the to
       *  integration algorithm). However, it provides an unbiased
       *  estimate of the monopole only if the random-random (RR) pairs
       *  do not depend on \f$\mu\f$. Otherwise, an
       *  \f$RR(\mu,s)/RR(s)\f$ weight has to be added to the pair
       *  count.
       *
       *  The first three multipole moments -- monopole, quadrupole and
       *  hexadecapole -- estimated with the <EM> integrated </EM>
       *  method can be obtained with the class
       *  TwoPointCorrelation_multipoles.
       */
      class TwoPointCorrelation1D_monopole : public virtual TwoPointCorrelation1D {

      protected:

	/**
	 *  @brief measure the monopole of the two-point correlation
	 *  function with Poisson errors
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
	 *  
	 */
	void measurePoisson (const std::string dir_output_pairs = par::defaultString, const std::vector<std::string> dir_input_pairs={}, const bool count_dd=true, const bool count_rr=true, const bool count_dr=true, const bool tcount=true, const Estimator estimator=Estimator::_LandySzalay_) override;

	/**
	 *  @brief measure the monopole of the two-point correlation
	 *  function estimating the covariance with Jackknife
	 *  resampling
	 *
	 *  @param dir_output_pairs output directory used to store the
	 *  number of pairs
	 *
	 *  @param dir_input_pairs vector of input directories used to
	 *  store the number of pairs (if the pairs are read from
	 *  files)
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
	 *  
	 */
	void measureJackknife (const std::string dir_output_pairs = par::defaultString, const std::vector<std::string> dir_input_pairs={}, const std::string dir_output_resample = par::defaultString, const bool count_dd=true, const bool count_rr=true, const bool count_dr=true, const bool tcount=true, const Estimator estimator=Estimator::_LandySzalay_) override;

	/**
	 *  @brief measure the monopole of the two-point correlation
	 *  function estimating the covariance with Jackknife resampling
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
	 *  
	 */
	void measureJackknifeTest (const std::string dir_output_pairs = par::defaultString, const std::vector<std::string> dir_input_pairs={}, const std::string dir_output_resample = par::defaultString, const bool count_dd=true, const bool count_rr=true, const bool count_dr=true, const bool tcount=true, const Estimator estimator=Estimator::_LandySzalay_) override;

	/**
	 *  @brief measure the monopole of the two-point correlation
	 *  function estimating the covariance with Bootstrap resampling
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
	 *  
	 */
	void measureBootstrap (const int nMocks, const std::string dir_output_pairs = par::defaultString, const std::vector<std::string> dir_input_pairs={}, const std::string dir_output_resample = par::defaultString, const bool count_dd=true, const bool count_rr=true, const bool count_dr=true, const bool tcount=true, const Estimator estimator=Estimator::_LandySzalay_, const int seed=3213) override;

      public:
    
	/**
	 *  @name Constructors/destructors
	 */
	///@{

	/**
	 *  @brief default constructor
	 *  1D_monopole
	 */
	TwoPointCorrelation1D_monopole () { m_twoPType = TwoPType::_monopole_; }

	/**
	 *  @brief constructor
	 *  @param data object of class Catalogue containing the input
	 *  catalogue
	 *  @param random of class Catalogue containing the random data
	 *  catalogue
	 *  @param binType binning type
	 *  @param rMin minimum separation used to count the pairs
	 *  @param rMax maximum separation used to count the pairs
	 *  @param nbins number of bins
	 *  @param shift shift parameter, i.e. the radial shift is
	 *  binSize*shift
	 *  @param angularUnits angular units
	 *  @param angularWeight angular weight function
	 *  @param compute_extra_info true &rarr; compute extra
	 *  information related to the pairs, such as the mean pair
	 *  separation and redshift
	 *  @param random_dilution_fraction fraction between the number
	 *  of objects in the diluted and original random samples, used
	 *  to improve performances in random-random pair counts
	 *  1D_monopole
	 */
	TwoPointCorrelation1D_monopole (const catalogue::Catalogue data, const catalogue::Catalogue random, const BinType binType, const double rMin, const double rMax, const int nbins, const double shift, const CoordinateUnits angularUnits=CoordinateUnits::_radians_, std::function<double(double)> angularWeight=nullptr, const bool compute_extra_info=false, const double random_dilution_fraction=1.)
	  : TwoPointCorrelation(data, random, compute_extra_info, random_dilution_fraction), TwoPointCorrelation1D(data, random, compute_extra_info, random_dilution_fraction)
	  { m_twoPType = TwoPType::_monopole_; set_parameters(binType, rMin, rMax, nbins, shift, angularUnits, angularWeight, compute_extra_info); }
      
	/**
	 *  @brief constructor
	 *  @param data object of class Catalogue containing the input
	 *  catalogue
	 *  @param random of class Catalogue containing the random data
	 *  catalogue
	 *  @param binType binning type
	 *  @param rMin minimum separation used to count the pairs
	 *  @param rMax maximum separation used to count the pairs
	 *  @param binSize the bin size
	 *  @param shift shift parameter, i.e. the radial shift is
	 *  binSize*shift
	 *  @param angularUnits angular units
	 *  @param angularWeight angular weight function
	 *  @param compute_extra_info true &rarr; compute extra
	 *  information related to the pairs, such as the mean pair
	 *  separation and redshift
	 *  @param random_dilution_fraction fraction between the number
	 *  of objects in the diluted and original random samples, used
	 *  to improve performances in random-random pair counts
	 *  1D_monopole
	 */
	TwoPointCorrelation1D_monopole (const catalogue::Catalogue data, const catalogue::Catalogue random, const BinType binType, const double rMin, const double rMax, const double binSize, const double shift, const CoordinateUnits angularUnits=CoordinateUnits::_radians_, std::function<double(double)> angularWeight=nullptr, const bool compute_extra_info=false, const double random_dilution_fraction=1.)
	  : TwoPointCorrelation(data, random, compute_extra_info, random_dilution_fraction), TwoPointCorrelation1D(data, random, compute_extra_info, random_dilution_fraction)
	  { m_twoPType = TwoPType::_monopole_; set_parameters(binType, rMin, rMax, binSize, shift, angularUnits, angularWeight, compute_extra_info); }

	/**
	 *  @brief default destructor
	 *  
	 */
	~TwoPointCorrelation1D_monopole () = default;

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
	 *  @param angularUnits angular units
	 *  @param angularWeight angular weight function
	 *  @param compute_extra_info true &rarr; compute extra
	 *  information related to the pairs, such as the mean pair
	 *  separation and redshift
	 *  
	 */
	void set_parameters (const BinType binType, const double rMin, const double rMax, const int nbins, const double shift, const CoordinateUnits angularUnits=CoordinateUnits::_radians_, std::function<double(double)> angularWeight=nullptr, const bool compute_extra_info=false);

	/**
	 *  @brief set the binning parameters
	 *  @param binType binning type
	 *  @param rMin minimum separation used to count the pairs
	 *  @param rMax maximum separation used to count the pairs
	 *  @param binSize the bin size
	 *  @param shift shift parameter, i.e. the radial shift is
	 *  binSize*shift
	 *  @param angularUnits angular units
	 *  @param angularWeight angular weight function
	 *  @param compute_extra_info true &rarr; compute extra
	 *  information related to the pairs, such as the mean pair
	 *  separation and redshift
	 *  
	 */
	void set_parameters (const BinType binType, const double rMin, const double rMax, const double binSize, const double shift, const CoordinateUnits angularUnits=CoordinateUnits::_radians_, std::function<double(double)> angularWeight=nullptr, const bool compute_extra_info=false);

	///@}


	/**
	 *  @name Member functions to count the number of pairs and measure the two-point correlation function
	 */
	///@{

	/**
	 *  @brief measure the monopole of the two-point correlation
	 *  function (with the direct estimator)
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
	 *  
	 */
	void measure (const ErrorType errorType=ErrorType::_Poisson_, const std::string dir_output_pairs=par::defaultString, const std::vector<std::string> dir_input_pairs={},  const std::string dir_output_resample=par::defaultString, const int nMocks=0, const bool count_dd=true, const bool count_rr=true, const bool count_dr=true, const bool tcount=true, const Estimator estimator=Estimator::_LandySzalay_, const int seed=3213) override;

	///@}

    
	/**
	 *  @name Input/Output methods
	 */
	///@{

	/**
	 *  @brief read the monopole of the two-point correlation
	 *  @param dir input directory
	 *  @param file input file
	 *  
	 */
	void read (const std::string dir, const std::string file) override;

	/**
	 *  @brief write the monopole of the two-point correlation
	 *  function
	 *  @param dir output directory
	 *  @param file output file
	 *  @param rank cpu index (for MPI usage)
	 *  
	 */
	void write (const std::string dir=par::defaultString, const std::string file=par::defaultString, const int rank=0) const override;

	///@}

      };
    }
  }
}

#endif
