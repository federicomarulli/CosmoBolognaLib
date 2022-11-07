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
 *  @file Headers/TwoPointCorrelation_wedges.h
 *
 *  @brief The class TwoPointCorrelation_wedges
 *
 *  This file defines the interface of the class
 *  TwoPointCorrelation_wedges, used to measure the wedges
 *  of the two-point correlation function
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unibo.it, alfonso.veropalumbo@unibo.it
 */

#ifndef __TWOPOINTWED__
#define __TWOPOINTWED__


#include "TwoPointCorrelation2D_polar.h"


// ===================================================================================================


namespace cbl {

  namespace measure {
  
    namespace twopt {

      /**
       *  @class TwoPointCorrelation_wedges TwoPointCorrelation_wedges.h
       *  "Headers/TwoPointCorrelation_wedges.h"
       *
       *  @brief The class TwoPointCorrelation_wedges
       *
       *  This class is used to handle objects of type <EM>
       *  TwoPointCorrelation_wedges </EM>. It is used to measure the
       *  wedges of the two-point correlation function,
       *  \f$\xi_{\Delta\mu} = \frac{1}{\Delta\mu}
       *  \int_{\mu_0}^{\mu_1} \xi(s,\mu)d\mu\f$. The default case is
       *  the one with two wedges, perpendicular - \f$\xi_\perp\f$ and
       *  parallel - \f$\xi_\parallel\f$, where \f$\mu_0=0\f$ and
       *  \f$\mu_1=0.5\f$ for \f$\xi_\perp\f$, and \f$\mu_0=0.5\f$ and
       *  \f$\mu_1=1\f$ for \f$\xi_\parallel\f$. In the most general
       *  case, we have \f$N\f$ wedges \f$\xi_{w,i}\f$, with \f$\mu_0
       *  = (i-1)/N\f$ and \f$\mu_1=i/N\f$, for \f$1\le i \le N\f$.
       */
      class TwoPointCorrelation_wedges : public TwoPointCorrelation2D_polar {

      protected:

	/// number of measured wedges
	int m_nWedges;

	/// integral limits used to measure the wedges
	std::vector<std::vector<double>> m_mu_integral_limits;
	
	/**
	 *  @brief pointer to the data containing the wedges of the
	 *  two-point correlation function
	 *  
	 *  @param rr input vector contaning the absolute separations
	 *  of the 2D two-point correlation function in polar
	 *  coordinates used to compute the wedges; these absolute
	 *  separations are also the separations at which the wedges
	 *  will be measured
	 *
	 *  @param mu input vector contaning the angular separations
	 *  of the 2D two-point correlation function in polar
	 *  coordinates used to compute the wedges
	 *
	 *  @param xi input matrix contaning the values of the 2D
	 *  two-point correlation function in polar coordinates used
	 *  to compute the wedges
	 *
	 *  @param error_xi input matrix contaning the values of the
	 *  errors of the 2D two-point correlation function in polar
	 *  coordinates used to compute the wedges
	 *
	 *  @return pointer to an object of type Data
	 */
	std::shared_ptr<data::Data> Wedges (const std::vector<double> rr, const std::vector<double> mu, const std::vector<std::vector<double> > xi, const std::vector<std::vector<double> > error_xi) override;

	/**
	 *  @brief pointer to the data object containing the values of
	 *  the wedges with the extra info retrieved from the 2D
	 *  two-point correlation function in polar coordinates
	 *  
	 *  @param rad vector containing the absolute separations
	 *
	 *  @param wedges vector containing the values of the wedges
	 *
	 *  @param error vector containing the errors
	 *
	 *  @return pointer to an object of type Data
	 */
	std::shared_ptr<data::Data> data_with_extra_info (const std::vector<double> rad, const std::vector<double> wedges, const std::vector<double> error) const;
      
	/**
	 *  @brief measure the wedges of the two-point correlation
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
	 *  @param fact factor used to compute the cell size of the
	 *  chain mesh: it is multiplied by the maximum distance
	 *  considered for the couples and can be setted by the user
	 *  to optimize the count of the couples
	 *
	 *  
	 */
	void measurePoisson (const std::string dir_output_pairs = par::defaultString, const std::vector<std::string> dir_input_pairs={}, const bool count_dd=true, const bool count_rr=true, const bool count_dr=true, const bool tcount=true, const Estimator estimator=Estimator::_LandySzalay_, const double fact=0.1) override;

	/**
	 *  @brief measure the wedges of the two-point correlation
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
	 *  @param fact factor used to compute the cell size of the
	 *  chain mesh: it is multiplied by the maximum distance
	 *  considered for the couples and can be setted by the user
	 *  to optimize the count of the couples
	 *
	 *  
	 */
	void measureJackknife (const std::string dir_output_pairs = par::defaultString, const std::vector<std::string> dir_input_pairs={}, const std::string dir_output_resample = par::defaultString, const bool count_dd=true, const bool count_rr=true, const bool count_dr=true, const bool tcount=true, const Estimator estimator=Estimator::_LandySzalay_, const double fact=0.1) override;

	/**
	 *  @brief measure the wedges of the two-point correlation
	 *  function estimating the covariance with Bootstrap resampling
	 *
	 *  @param nMocks number of mocks to be generated with
	 *  bootstrap resampling
	 *
	 *  @param dir_output_pairs output directory used to store the
	 *  number of pairs
	 *
	 *  @param dir_input_pairs vector of input directories used to
	 *  store the number of pairs (if the pairs are read from
	 *  files)
	 *
	 *  @param dir_output_resample output directory used to store
	 *  the Bootstrap resampling correlation functions, with
	 *  Poisson errors; if an empty string (i.e. "" or "NULL") is
	 *  provided, no output will be stored
	 *
	 *  @param count_dd true &rarr; count the number of data-data
	 *  pairs; false &rarr; read the number of data-data pairs
	 *  from file
	 *
	 *  @param count_rr true &rarr; count the number of
	 *  random-random pairs; false &rarr; read the number of
	 *  random-random pairs from file
	 *
	 *  @param count_dr true &rarr; count the number of
	 *  data-random pairs; false &rarr; read the number of
	 *  data-random pairs
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
	void measureBootstrap (const int nMocks, const std::string dir_output_pairs = par::defaultString, const std::vector<std::string> dir_input_pairs={}, const std::string dir_output_resample = par::defaultString, const bool count_dd=true, const bool count_rr=true, const bool count_dr=true, const bool tcount=true, const Estimator estimator=Estimator::_LandySzalay_, const double fact=0.1, const int seed=3213) override;

	/**
	 *  @brief measure the jackknife resampling of the wedges of the
	 *  two-point correlation function
	 *
	 *  @param dd vector of data-data pairs, divided per regions
	 *
	 *  @param rr vector of random-random pairs, divided per regions
	 *
	 *  @return a vector of pointers to objects of type Data
	 */
	std::vector<std::shared_ptr<data::Data> > XiJackknife (const std::vector<std::shared_ptr<pairs::Pair> > dd, const std::vector<std::shared_ptr<pairs::Pair> > rr) override;

	/**
	 *  @brief measure the jackknife resampling of the wedges of
	 *  two-point correlation function
	 *
	 *  @param dd vector of data-data pairs, divided per regions
	 *
	 *  @param rr vector of random-random pairs, divided per regions
	 *
	 *  @param dr vector of data-random pairs, divided per regions  
	 *
	 *  @return a vector of pointers to objects of type Data
	 */
	std::vector<std::shared_ptr<data::Data> > XiJackknife (const std::vector<std::shared_ptr<pairs::Pair> > dd, const std::vector<std::shared_ptr<pairs::Pair> > rr, const std::vector<std::shared_ptr<pairs::Pair> > dr) override;

	/**
	 *  @brief measure the bootstrap resampling of the wedges of the
	 *  two-point correlation function
	 * 
	 *  @param nMocks number of bootstrap resamplings
	 *
	 *  @param dd vector of data-data pairs, divided per regions
	 *
	 *  @param rr vector of random-random pairs, divided per regions
	 *
	 *  @param seed the seed for random number generation
	 *
	 *  @return a vector of pointers to objects of type Data
	 */
	std::vector<std::shared_ptr<data::Data> > XiBootstrap (const int nMocks, const std::vector<std::shared_ptr<pairs::Pair> > dd, const std::vector<std::shared_ptr<pairs::Pair> > rr, const int seed=3213) override;

	/**
	 *  @brief measure the bootstrap resampling of the wedges of
	 *  the two-point correlation function
	 *
	 *  @param nMocks number of bootstrap resamplings
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
	std::vector<std::shared_ptr<data::Data> > XiBootstrap (const int nMocks, const std::vector<std::shared_ptr<pairs::Pair> > dd, const std::vector<std::shared_ptr<pairs::Pair> > rr, const std::vector<std::shared_ptr<pairs::Pair> > dr, const int seed=3213) override;

      public:

	/**
	 *  @name Constructors/destructors
	 */
	///@{

	/**
	 *  @brief default constructor
	 *  _wedges
	 */
	TwoPointCorrelation_wedges () { m_twoPType = TwoPType::_wedges_; }
	
	/**
	 *  @brief constructor, using the number of bins
	 *
	 *  this constructor is used to used to measure the wedges of
	 *  the two-point correlation function, \f$\xi_{\Delta\mu} =
	 *  \frac{1}{\Delta\mu} \int_{\mu_0}^{\mu_1}
	 *  \xi(s,\mu)d\mu\f$. The default case is the one with two
	 *  wedges, perpendicular - \f$\xi_\perp\f$ and parallel -
	 *  \f$\xi_\parallel\f$, where \f$\mu_0=0\f$ and
	 *  \f$\mu_1=0.5\f$ for \f$\xi_\perp\f$, and \f$\mu_0=0.5\f$
	 *  and \f$\mu_1=1\f$ for \f$\xi_\parallel\f$. In the most
	 *  general case, we have \f$N\f$ wedges \f$\xi_{w,i}\f$, with
	 *  \f$\mu_0 = (i-1)/N\f$ and \f$\mu_1=i/N\f$, for \f$1\le i
	 *  \le N\f$.
	 *
	 *  @param data object of class Catalogue containing the input
	 *  catalogue
	 *
	 *  @param random of class Catalogue containing the random data
	 *  catalogue
	 *
	 *  @param binType_rad binning type in absolute separations
	 *
	 *  @param rMin minimum absolute separation used to count the
	 *  pairs
	 *
	 *  @param rMax maximum absolute separation used to count the
	 *  pairs
	 *
	 *  @param nbins_rad number of bins in the absolute
	 *  separation
	 *
	 *  @param shift_rad shift parameter in the absolute
	 *  separation, i.e. the radial shift is binSize*shift
	 *
	 *  @param nWedges number of wedges to be measured; the
	 *  default is two wedges, \f$\xi_\perp\f$ and
	 *  \f$\xi_\parallel\f$

	 *  @param nbins_mu number of bins in the angular separation
	 *  used to integrate
	 *
	 *  @param shift_mu shift parameter in the angular separation,
	 *  i.e. the radial shift is binSize*shift
	 *  
	 *  @param mu_integral_limits the \f$\mu\f$ integral limits
	 *  used to measure the wedges
	 *
	 *  @param angularUnits angular units
	 *
	 *  @param angularWeight angular weight function
	 *
	 *  @param compute_extra_info true &rarr; compute extra
	 *  information related to the pairs, such as the mean pair
	 *  separation and redshift
	 *
	 *  @param random_dilution_fraction fraction between the
	 *  number of objects in the diluted and original random
	 *  samples, used to improve performances in random-random
	 *  pair counts
	 *
	 *  _wedges
	 */
	TwoPointCorrelation_wedges (catalogue::Catalogue data, catalogue::Catalogue random, const BinType binType_rad, const double rMin, const double rMax, const int nbins_rad, const double shift_rad, const int nWedges=2, const int nbins_mu=50, const double shift_mu=0.5, const std::vector<std::vector<double>> mu_integral_limits={{0., 0.5}, {0.5, 1.}}, const CoordinateUnits angularUnits=CoordinateUnits::_radians_, std::function<double(double)> angularWeight=nullptr, const bool compute_extra_info=false, const double random_dilution_fraction=1.)
	  : TwoPointCorrelation2D_polar(data, random, binType_rad, rMin, rMax, nbins_rad, shift_rad, BinType::_linear_, 0., 1., nbins_mu, shift_mu, angularUnits, angularWeight, compute_extra_info, random_dilution_fraction)
	  {
	    m_twoPType = TwoPType::_wedges_;
	    
	    if (nWedges<2)
	      cbl::ErrorCBL("nWedges must be larger than 2!", "TwoPointCorrelation_wedges", "TwoPointCorrelation_wedges.h");

	    if (int(mu_integral_limits.size())!=nWedges)
	      cbl::ErrorCBL("mu_integral_limits.size()="+cbl::conv(mu_integral_limits.size(), cbl::par::fINT)+" must be equal to nWedges="+cbl::conv(nWedges, cbl::par::fINT), "TwoPointCorrelation_wedges", "TwoPointCorrelation_wedges.h");
	    
	    m_nWedges = nWedges;
	    m_mu_integral_limits = mu_integral_limits;
	  }

	/**
	 *  @brief constructor, using the bin size
	 *
	 *  this constructor is used to used to measure the wedges of
	 *  the two-point correlation function, \f$\xi_{\Delta\mu} =
	 *  \frac{1}{\Delta\mu} \int_{\mu_0}^{\mu_1}
	 *  \xi(s,\mu)d\mu\f$. The default case is the one with two
	 *  wedges, perpendicular - \f$\xi_\perp\f$ and parallel -
	 *  \f$\xi_\parallel\f$, where \f$\mu_0=0\f$ and
	 *  \f$\mu_1=0.5\f$ for \f$\xi_\perp\f$, and \f$\mu_0=0.5\f$
	 *  and \f$\mu_1=1\f$ for \f$\xi_\parallel\f$. In the most
	 *  general case, we have \f$N\f$ wedges \f$\xi_{w,i}\f$, with
	 *  \f$\mu_0 = (i-1)/N\f$ and \f$\mu_1=i/N\f$, for \f$1\le i
	 *  \le N\f$.
	 *
	 *  @param data object of class Catalogue containing the input
	 *  catalogue
	 *
	 *  @param random of class Catalogue containing the random
	 *  data catalogue
	 *
	 *  @param binType_rad binning type in absolute separations
	 *
	 *  @param rMin minimum absolute separation used to count the
	 *  pairs
	 *
	 *  @param rMax maximum absolute separation used to count the
	 *  pairs
	 *
	 *  @param binSize_rad bin size in the absolute separation
	 *
	 *  @param shift_rad shift parameter in the absolute
	 *  separation, i.e. the radial shift is binSize*shift
	 *
	 *  @param nWedges number of wedges to be measured; the
	 *  default is two wedges, \f$\xi_\perp\f$ and
	 *  \f$\xi_\parallel\f$
	 *
	 *  @param binSize_mu bin size in the angular separation
	 *
	 *  @param shift_mu shift parameter in the angular separation,
	 *  i.e. the radial shift is binSize*shift
	 *
	 *  @param mu_integral_limits the \f$\mu\f$ integral limits
	 *  used to measure the wedges
	 *
	 *  @param angularUnits angular units
	 *
	 *  @param angularWeight angular weight function
	 *
	 *  @param compute_extra_info true &rarr; compute extra
	 *  information related to the pairs, such as the mean pair
	 *  separation and redshift
	 *
	 *  @param random_dilution_fraction fraction between the
	 *  number of objects in the diluted and original random
	 *  samples, used to improve performances in random-random
	 *  pair counts
	 *
	 *  _wedges
	 */
	TwoPointCorrelation_wedges (catalogue::Catalogue data, catalogue::Catalogue random, const BinType binType_rad, const double rMin, const double rMax, const double binSize_rad, const double shift_rad, const int nWedges=2, const double binSize_mu=1., const double shift_mu=0.5, const std::vector<std::vector<double>> mu_integral_limits={{0., 0.5}, {0.5, 1.}}, const CoordinateUnits angularUnits=CoordinateUnits::_radians_, std::function<double(double)> angularWeight=nullptr, const bool compute_extra_info=false, const double random_dilution_fraction=1.)
	  : TwoPointCorrelation2D_polar(data, random, binType_rad, rMin, rMax, binSize_rad, shift_rad, BinType::_linear_, 0., 1., binSize_mu, shift_mu, angularUnits, angularWeight, compute_extra_info, random_dilution_fraction)
	  {
	    m_twoPType = TwoPType::_wedges_;

	    if (nWedges<2)
	      cbl::ErrorCBL("nWedges must be larger than 2!", "TwoPointCorrelation_wedges", "TwoPointCorrelation_wedges.h");

	    if (int(mu_integral_limits.size())!=nWedges)
	      cbl::ErrorCBL("mu_integral_limits.size()="+cbl::conv(mu_integral_limits.size(), cbl::par::fINT)+" must be equal to nWedges="+cbl::conv(nWedges, cbl::par::fINT), "TwoPointCorrelation_wedges", "TwoPointCorrelation_wedges.h");
	    
	    m_nWedges = nWedges;
	    m_mu_integral_limits = mu_integral_limits;
	  }

	/**
	 *  @brief default destructor
	 *  
	 */
	~TwoPointCorrelation_wedges () = default;

	///@}

      	/**
	 *  @name Member functions to get the private/protected
	 *  variables of the class
	 */
	///@{
	
	/**
	 *  @brief get the number of measured wedges
	 *  @return the number of measured wedges
	 */
	int nWedges () const
	{ return m_nWedges; }
	
	/**
	 *  @brief get the integral limits used to measure the wedges
	 *  @return the integral limits used to measure the wedges
	 */
	std::vector<std::vector<double>> mu_integral_limits () const
	{ return m_mu_integral_limits; }
	
	/**
	 *  @brief get the x coordinates
	 *  @return the x coordinates
	 */
	std::vector<double> xx () const override;

	/**
	 *  @brief get the perpendicular wedge
	 *  @return the perpendicular wedge
	 */
	std::vector<double> xiPerpendicular () const override;

	/**
	 *  @brief get the errors on the perpendicular wedge
	 *  @return the errors on the perpendicular wedge
	 */
	std::vector<double> errorPerpendicular () const override;

	/**
	 *  @brief get the parallel wedge 
	 *  @return the parallel wedge 
	 */
	std::vector<double> xiParallel () const override;

	/**
	 *  @brief get the errors on the parallel wedge
	 *  @return the error on the parallel wedge
	 */
	std::vector<double> errorParallel () const override;

	/**
	 *  @brief get the y coordinates
	 *  @return the y coordinates
	 */
	std::vector<double> yy () const 
	  { cbl::ErrorCBL("", "yy", "TwoPointCorrelation_wedges.h"); std::vector<double> vv; return vv; }

	/**
	 *  @brief get the the binned correlation function 
	 *  @return the binned correlation function 
	 */
	std::vector<double> xi1D () const
	  { cbl::ErrorCBL("", "xi1D", "TwoPointCorrelation_wedges.h"); std::vector<double> vv; return vv; }

	/**
	 *  @brief get the error on the binned correlation function
	 *  function
	 *  @return the error on the binned correlation function
	 *  function
	 */
	std::vector<double> error1D () const
	  { cbl::ErrorCBL("", "error1D", "TwoPointCorrelation_wedges.h"); std::vector<double> vv; return vv; }

	/**
	 *  @brief get the the binned 2D correlation function 
	 *  @return the binned 2D correlation function 
	 */
	std::vector<std::vector<double> > xi2D () const 
	  { cbl::ErrorCBL("", "xi2D", "TwoPointCorrelation_wedges.h"); std::vector<std::vector<double> > vv; return vv; }

	/**
	 *  @brief get the errors on the binned 2D correlation function
	 *  function
	 *  @return the errors on the binned 2D correlation function
	 *  function
	 */
	std::vector<std::vector<double> > error2D () const 
	  { cbl::ErrorCBL("", "error2D", "TwoPointCorrelation_wedges.h"); std::vector<std::vector<double> > vv; return vv; }

	///@}

	
	/**
	 *  @name Member functions to count the number of pairs and measure the two-point correlation function
	 */
	///@{

	/**
	 *  @brief measure the wedges of the two-point correlation
	 *  function
	 *
	 *  @param errorType type of error
	 *  
	 *  @param dir_output_pairs output directory used to store the
	 *  number of pairs
	 *
	 *  @param dir_input_pairs vector of input directories used to
	 *  store the number of pairs (if the pairs are read from
	 *  files)
	 *
	 *  @param dir_output_resample output directory of the
	 *  resampling correlation functions; if an empty string
	 *  (i.e. "" or "NULL") is provided, no output will be stored
	 *
	 *  @param nMocks number of resampling used for bootstrap
	 *
	 *  @param count_dd true &rarr; count the number of data-data
	 *  pairs; false &rarr; read the number of data-random pairs
	 *  from file
	 *
	 *  @param count_rr true &rarr; count the number of
	 *  random-random pairs; false &rarr; read the number of
	 *  random-random pairs
	 *
	 *  @param count_dr true &rarr; count the number of
	 *  data-random pairs; false &rarr; read the number of
	 *  data-random pairs
	 *
	 *  @param tcount true &rarr; activate the time counter; false
	 *  &rarr; no time counter
	 *
	 *  @param estimator the estimator used to measure the
	 *  two-point correlation function
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
	void measure (const ErrorType errorType=ErrorType::_Poisson_, const std::string dir_output_pairs=par::defaultString, const std::vector<std::string> dir_input_pairs={},  const std::string dir_output_resample=par::defaultString, const int nMocks=0, const bool count_dd=true, const bool count_rr=true, const bool count_dr=true, const bool tcount=true, const Estimator estimator=Estimator::_LandySzalay_, const double fact=0.1, const int seed=3213) override;

	///@}


	/**
	 *  @name Member functions for input/output 
	 */
	///@{

	/**
	 *  @brief read the wedges of the two-point correlation function
	 *  @param dir input directory
	 *  @param file input file
	 *  
	 */
	void read (const std::string dir, const std::string file) override
	{ (void)dir; (void)file; ErrorCBL("", "read", "TwoPointCorrelation_wedges.h", glob::ExitCode::_workInProgress_); }  

	/**
	 *  @brief write the wedges of the two-point correlation
	 *  function
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
