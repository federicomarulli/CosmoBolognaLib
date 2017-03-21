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
 *  @file Headers/Lib/TwoPointCorrelation_wedges.h
 *
 *  @brief The class TwoPointCorrelation_wedges
 *
 *  This file defines the interface of the class
 *  TwoPointCorrelation_wedges, used to measure the wedges
 *  of the two-point correlation function
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#ifndef __TWOPOINTWED__
#define __TWOPOINTWED__


#include "TwoPointCorrelation2D_polar.h"


// ===================================================================================================


namespace cosmobl {

  namespace twopt {


    /**
     *  @class TwoPointCorrelation_wedges TwoPointCorrelation_wedges.h
     *  "Headers/Lib/TwoPointCorrelation_wedges.h"
     *
     *  @brief The class TwoPointCorrelation_wedges
     *
     *  This class is used to handle objects of type <EM>
     *  TwoPointCorrelation_wedges </EM>. It is used to measure the
     *  wedges \f$\xi_\perp, \xi_\parallel\f$ as \f$\xi_{\Delta\mu} =
     *  \frac{1}{2} \int_{\mu_0}^{\mu_1} \xi(s,\mu)d\mu\f$ with
     *  \f$\mu_0=0\f$ and \f$\mu_1=0.5\f$ for \f$\xi_\perp\f$ and \f$\mu_0=0.5\f$
     *  and \f$\mu_1=1\f$ for \f$\xi_\parallel\f$
     */
    class TwoPointCorrelation_wedges : public TwoPointCorrelation2D_polar {

    protected:

      /**
       *  @brief measure the wedges of the two-point correlation
       *  function
       *  
       *  @param rr absolute separation 
       *
       *  @param mu angular separation
       *
       *  @param xi 2d cartesian 2pcf
       *
       *  @param error_xi error on the 2d polar 2pcf
       *
       *  @return pointer to an object of type Data
       */
      shared_ptr<data::Data> Wedges (const vector<double> rr, const vector<double> mu, const vector<vector<double> > xi, const vector<vector<double> > error_xi) override;

      /**
       *  @brief return a data object with extra info
       *  
       *  @param rad vector containing the binned separations
       *
       *  @param wedges vector containing the binned wedges of the
       *  correlation function
       *
       *  @param error vector containing the errors
       *
       *  @return pointer to an object of type Data
       */
      shared_ptr<data::Data> data_with_extra_info (const vector<double> rad, const vector<double> wedges, const vector<double> error) const;
      
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
       *  @return none
       */
      void measurePoisson (const string dir_output_pairs = par::defaultString, const vector<string> dir_input_pairs={}, const bool count_dd=true, const bool count_rr=true, const bool count_dr=true, const bool tcount=true, const Estimator estimator=_LandySzalay_) override;

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
       *  the Jackknife resampling correlation functions, with Poisson
       *  errors
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
      void measureJackknife (const string dir_output_pairs = par::defaultString, const vector<string> dir_input_pairs={}, const string dir_output_resample = par::defaultString, const bool count_dd=true, const bool count_rr=true, const bool count_dr=true, const bool tcount=true, const Estimator estimator=_LandySzalay_) override;

      /**
       *  @brief measure the wedges of the two-point correlation
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
       *  the Bootstrap resampling correlation function, with Poisson
       *  errors
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
      void measureBootstrap (const int nMocks, const string dir_output_pairs = par::defaultString, const vector<string> dir_input_pairs={}, const string dir_output_resample = par::defaultString, const bool count_dd=true, const bool count_rr=true, const bool count_dr=true, const bool tcount=true, const Estimator estimator=_LandySzalay_) override;

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
      vector<shared_ptr<data::Data> > XiJackknife (const vector<shared_ptr<pairs::Pair> > dd, const vector<shared_ptr<pairs::Pair> > rr) override;

      /**
       *  @brief measure the jackknife resampling of the wedges of
       *  two-point correlation function
       *
       *  @param dd vector of data-data pairs, divided per regions
       *
       *  @param rr vector of random-random pairs, divided per regions
       *
       *  @param dr vector of random-random pairs, divided per regions  
       *
       *  @return a vector of pointers to objects of type Data
       */
      vector<shared_ptr<data::Data> > XiJackknife (const vector<shared_ptr<pairs::Pair> > dd, const vector<shared_ptr<pairs::Pair> > rr, const vector<shared_ptr<pairs::Pair> > dr) override;

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
       *  @return a vector of pointers to objects of type Data
       */
      vector<shared_ptr<data::Data> > XiBootstrap (const int nMocks, const vector<shared_ptr<pairs::Pair> > dd, const vector<shared_ptr<pairs::Pair> > rr) override;

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
       *  @param dr vector of random-random pairs, divided per regions  
       *
       *  @return a vector of pointers to objects of type Data
       */
      vector<shared_ptr<data::Data> > XiBootstrap (const int nMocks, const vector<shared_ptr<pairs::Pair> > dd, const vector<shared_ptr<pairs::Pair> > rr, const vector<shared_ptr<pairs::Pair> > dr) override;

    public:

      /**
       *  @name Constructors/destructors
       */
      ///@{

      /**
       *  @brief default constructor
       *  @return object of class TwoPointCorrelation_wedges
       */
      TwoPointCorrelation_wedges () { m_twoPType = _1D_wedges_; }

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
       *  @return object of class TwoPointCorrelation_wedges
       */
      TwoPointCorrelation_wedges (catalogue::Catalogue data, catalogue::Catalogue random, const binType binType_rad, const double rMin, const double rMax, const int nbins_rad, const double shift_rad, const double muMin, const double muMax, const int nbins_mu, const double shift_mu, const CoordUnits angularUnits=_radians_, function<double(double)> angularWeight=nullptr, const bool compute_extra_info=false, const double random_dilution_fraction=1.)
	: TwoPointCorrelation2D_polar(data, random, binType_rad, rMin, rMax, nbins_rad, shift_rad, _linear_, muMin, muMax, nbins_mu, shift_mu, angularUnits, angularWeight, compute_extra_info, random_dilution_fraction)
	{ m_twoPType = _1D_wedges_; }

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
       *  @return object of class TwoPointCorrelation_wedges
       */
      TwoPointCorrelation_wedges (catalogue::Catalogue data, catalogue::Catalogue random, const binType binType_rad, const double rMin, const double rMax, const double binSize_rad, const double shift_rad, const double muMin, const double muMax, const double binSize_mu, const double shift_mu, const CoordUnits angularUnits=_radians_, function<double(double)> angularWeight=nullptr, const bool compute_extra_info=false, const double random_dilution_fraction=1.)
	: TwoPointCorrelation2D_polar(data, random, binType_rad, rMin, rMax, binSize_rad, shift_rad, _linear_, muMin, muMax, binSize_mu, shift_mu, angularUnits, angularWeight, compute_extra_info, random_dilution_fraction)
	{ m_twoPType = _1D_wedges_; }

      /**
       *  @brief default destructor
       *  @return none
       */
      ~TwoPointCorrelation_wedges () = default;

      ///@}

      
      /**
       *  @brief get the x coordinates
       *  @return the x coordinates
       */
      vector<double> xx () const  override;

      /**
       *  @brief get the perpendicular wedge
       *  @return the perpendicular wedge
       */
      vector<double> xiPerpendicular () const override;

      /**
       *  @brief get the errors on the perpendicular wedge
       *  @return the errors on the perpendicular wedge
       */
      vector<double> errorPerpendicular () const override;

      /**
       *  @brief get the parallel wedge 
       *  @return the parallel wedge 
       */
      vector<double> xiParallel () const override;

      /**
       *  @brief get the errors on the parallel wedge
       *  @return the error on the parallel wedge
       */
      vector<double> errorParallel () const override;

      /**
       *  @brief get the y coordinates
       *  @return the y coordinates
       */
      vector<double> yy () const 
	{ cosmobl::ErrorCBL("Error in yy() of TwoPointCorrelation_wedges.h!"); vector<double> vv; return vv; }

      /**
       *  @brief get the the binned correlation function 
       *  @return the binned correlation function 
       */
      vector<double> xi1D () const
	{ cosmobl::ErrorCBL("Error in xi1D() of TwoPointCorrelation_wedges.h!"); vector<double> vv; return vv; }

      /**
       *  @brief get the error on the binned correlation function
       *  function
       *  @return the error on the binned correlation function
       *  function
       */
      vector<double> error1D () const
	{ cosmobl::ErrorCBL("Error in error1D() of TwoPointCorrelation_wedges.h!"); vector<double> vv; return vv; }

      /**
       *  @brief get the the binned 2D correlation function 
       *  @return the binned 2D correlation function 
       */
      vector<vector<double> > xi2D () const 
	{ cosmobl::ErrorCBL("Error in xi2D() of TwoPointCorrelation_wedges.h!"); vector<vector<double> > vv; return vv; }

      /**
       *  @brief get the errors on the binned 2D correlation function
       *  function
       *  @return the errors on the binned 2D correlation function
       *  function
       */
      vector<vector<double> > error2D () const 
	{ cosmobl::ErrorCBL("Error in error2D() of TwoPointCorrelation_multipoles.h!"); vector<vector<double> > vv; return vv; }

      
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
       *  store the number of pairs (if the pairs are read from files)
       *
       *  @param dir_output_resample output directory of the
       *  resampled correlation function
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
       *  @return none
       */
      void measure (const ErrorType errorType=ErrorType::_Poisson_, const string dir_output_pairs=par::defaultString, const vector<string> dir_input_pairs={},  const string dir_output_resample=par::defaultString, const int nMocks=0, const bool count_dd=true, const bool count_rr=true, const bool count_dr=true, const bool tcount=true, const Estimator estimator=_LandySzalay_) override;

      ///@}


      /**
       *  @name Member functions for input/output 
       */
      ///@{

      /**
       *  @brief read the wedges of the two-point correlation function
       *  @param dir input directory
       *  @param file input file
       *  @return none
       */
      void read (const string dir, const string file) override
      { (void)dir; (void)file; ErrorCBL("Error in TwoPointCorrelation_wedges::read of TwoPointCorrelation_wedges.h: work in progress!"); }  

      /**
       *  @brief write the wedges of the two-point correlation
       *  function
       *  @param dir output directory
       *  @param file output file
       *  @param rank cpu index (for MPI usage)
       *  @return none
       */
      void write (const string dir=par::defaultString, const string file=par::defaultString, const int rank=0) const override;

      ///@}

      
      /**
       *  @name Member functions to compute, read and write the covariance matrix
       */
      ///@{ 

      /**
       *  @brief read the measured covariance matrix
       *  @param dir input directory
       *  @param file input file
       *  @return none
       */
      void read_covariance (const string dir, const string file) override;

      /**
       *  @brief write the measured two-point correlation
       *  @param dir output directory
       *  @param file output file
       *  @return none
       */
      void write_covariance (const string dir, const string file) const override;

      /**
       *  @brief compute the covariance matrix
       *  @param xi vector containing the measure correlation
       *  functions used to compute the covariance matrix
       *  @param JK true &rarr; compute the jackknife covariance
       *  matrix; false compute the standard covariance matrix
       *  @return none
       */
      void compute_covariance (const vector<shared_ptr<data::Data>> xi, const bool JK) override;

      /**
       *  @brief compute the covariance matrix
       *  @param file vector containing the input files with the
       *  measured correlation functions used to compute the
       *  covariance matrix
       *  @param JK true &rarr; compute the jackknife covariance
       *  matrix; false compute the standard covariance matrix
       *  @return none
       */
      void compute_covariance (const vector<string> file, const bool JK) override;

      ///@} 

    };
  }
}

#endif
