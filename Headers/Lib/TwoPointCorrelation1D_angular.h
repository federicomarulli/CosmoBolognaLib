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
 *  @file Headers/Lib/TwoPointCorrelation1D_angular.h
 *
 *  @brief The class TwoPointCorrelation1D_angular
 *
 *  This file defines the interface of the class
 *  TwoPointCorrelation1D_angular, used to measure the monopole of the
 *  two-point correlation function
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#ifndef __TWOPOINTANG__
#define __TWOPOINTANG__


#include "TwoPointCorrelation1D.h"


// ===================================================================================================


namespace cosmobl {

  namespace twopt {
  
    /**
     *  @class TwoPointCorrelation1D_angular TwoPointCorrelation1D_angular.h
     *  "Headers/Lib/TwoPointCorrelation1D_angular.h"
     *
     *  @brief The class TwoPointCorrelation1D_angular
     *
     *  This class is used to handle objects of type <EM>
     *  TwoPointCorrelation1D_angular </EM>. It is used to measure the
     *  angular two-point correlation function, \f$w(\theta)\f$,
     *  defined as \f$dP_{12}=n^2[1+w(\theta)]dV_1dV_2\f$, where
     *  \f$n\f$ is the average number density, and \f$dP_{12}\f$ is
     *  the probability of finding a pair with one object in the
     *  volume \f$dV_1\f$ and the other one in the volume \f$dV_2\f$,
     *  separated by an angle \f$\theta\f$.
     */
    class TwoPointCorrelation1D_angular : public TwoPointCorrelation1D {

      /**
       *  @brief measure the angular two-point correlation
       *  function, &xi;(r) with Poisson error
       *  @param dir_output_pairs output directory used to store the
       *  number of pairs
       *  @param dir_input_pairs vector of input directories used to
       *  store the number of pairs (if the pairs are read from files)
       *  @param count_dd 1 &rarr; count the number of data-data
       *  opairs; 0 &rarr; read the number of data-data pairs from
       *  file
       *  @param count_rr 1 &rarr; count the number of random-random
       *  opairs; 0 &rarr; read the number of random-random pairs from
       *  file
       *  @param count_dr 1 &rarr; count the number of data-random
       *  opairs; 0 &rarr; read the number of data-random pairs from
       *  file
       *  @param tcount 1 &rarr; activate the time counter; 0 &rarr;
       *  don't activate the time counter; 
       *  @return none
       */
      void measurePoisson (const string dir_output_pairs=par::defaultString, const vector<string> dir_input_pairs={}, const int count_dd=1, const int count_rr=1, const int count_dr=1, const bool tcount=1) override;

      /**
       *  @brief measure the angular two-point correlation
       *  function, &xi;(r), estimate the covariance with Jackknife resampling
       *
       *  @param dir_output_pairs output directory used to store the
       *  number of pairs
       *
       *  @param dir_input_pairs vector of input directories used to store the
       *  number of pairs (if the pairs are read from files)
       *
       *  @param dir_output_ResampleXi output directory used to store the
       *  Jackknife resampling, with Poisson error
       *
       *  @param count_dd 1 &rarr; count the number of data-data
       *  opairs; 0 &rarr; read the number of data-data pairs from
       *  file
       *
       *  @param count_rr 1 &rarr; count the number of random-random
       *  opairs; 0 &rarr; read the number of random-random pairs from
       *  file
       *
       *  @param count_dd 1 &rarr; count the number of data-random
       *  opairs; 0 &rarr; read the number of data-random pairs from
       *  file
       *
       *  @param count_rr 1 &rarr; count the number of random-random
       *  pairs; 0 &rarr; read the number of random-random pairs
       *
       *  @param count_dr 1 &rarr; count the number of data-random
       *  pairs; 0 &rarr; read the number of data-random pairs
       *
       *  @param tcount 1 &rarr; activate the time counter; 0 &rarr;
       *  don't activate the time counter; 
       *
       *  @return none
       */
      void measureJackknife (const string dir_output_pairs=par::defaultString, const vector<string> dir_input_pairs={}, const string dir_output_ResampleXi = par::defaultString, const int count_dd=1, const int count_rr=1, const int count_dr=1, const bool tcount=1) override;

      /**
       *  @brief measure the angular two-point correlation
       *  function, &xi;(r), estimate the covariance with Bootstrap resampling
       *
       *  @param nMocks number of mocks to be generated with bootstrap resampling
       *
       *  @param dir_output_pairs output directory used to store the
       *  number of pairs
       *
       *  @param dir_input_pairs vector of input directories used to store the
       *  number of pairs (if the pairs are read from files)
       *
       *  @param dir_output_ResampleXi output directory used to store the
       *  bootstrap resampling,  with Poisson error
       *
       *  @param count_dd 1 &rarr; count the number of data-data
       *  opairs; 0 &rarr; read the number of data-data pairs from
       *  file
       *
       *  @param count_rr 1 &rarr; count the number of random-random
       *  opairs; 0 &rarr; read the number of random-random pairs from
       *  file
       *
       *  @param count_dd 1 &rarr; count the number of data-random
       *  opairs; 0 &rarr; read the number of data-random pairs from
       *  file
       *
       *  @param count_rr 1 &rarr; count the number of random-random
       *  pairs; 0 &rarr; read the number of random-random pairs
       *
       *  @param count_dr 1 &rarr; count the number of data-random
       *  pairs; 0 &rarr; read the number of data-random pairs
       *
       *  @param tcount 1 &rarr; activate the time counter; 0 &rarr;
       *  don't activate the time counter; 
       *
       *  @return none
       */
      void measureBootstrap (const int nMocks, const string dir_output_pairs = par::defaultString, const vector<string> dir_input_pairs={}, const string dir_output_ResampleXi = par::defaultString, const int count_dd=1, const int count_rr=1, const int count_dr=1, const bool tcount=1) override;

    public:
    
      /**
       *  @name Constructors/destructors
       */
      ///@{

      /**
       *  @brief default constructor
       *  @return object of class TwoPointCorrelation1D_angular
       */
      TwoPointCorrelation1D_angular () { m_twoPType = _1D_angular_; }

      /**
       *  @brief constructor
       *  @param data object of class Catalogue containing the input
       *  catalogue
       *  @param random of class Catalogue containing the random data
       *  catalogue
       *  @param binType binning type: 0 &rarr; linear; 1 &rarr;
       *  logarithmic
       *  @param thetaMin minimum angular separation used to count the pairs
       *  @param thetaMax maximum angular separation used to count the pairs
       *  @param nbins number of bins
       *  @param shift shift parameter, i.e. the radial shift is
       *  binSize*shift
       *  @param angularUnits angular units
       *  @param angularWeight angular weight function
       *  @return object of class TwoPointCorrelation1D_angular
       */
      TwoPointCorrelation1D_angular (const catalogue::Catalogue data, const catalogue::Catalogue random, const binType binType, const double thetaMin, const double thetaMax, const int nbins, const double shift, const CoordUnits angularUnits=_radians_, function<double(double)> angularWeight=nullptr)
	: TwoPointCorrelation1D(data, random) { m_twoPType = _1D_angular_; set_parameters(binType, thetaMin, thetaMax, nbins, shift, angularUnits, angularWeight); }

      /**
       *  @brief constructor
       *  @param data object of class Catalogue containing the input
       *  catalogue
       *  @param random of class Catalogue containing the random data
       *  catalogue
       *  @param binType binning type: 0 &rarr; linear; 1 &rarr;
       *  logarithmic
       *  @param thetaMin minimum angular separation used to count the pairs
       *  @param thetaMax maximum angular separation used to count the pairs
       *  @param binSize the bin size
       *  @param shift shift parameter, i.e. the radial shift is
       *  binSize*shift
       *  @param angularUnits angular units
       *  @param angularWeight angular weight function
       *  @return object of class TwoPointCorrelation1D_angular
       */
      TwoPointCorrelation1D_angular (const catalogue::Catalogue data, const catalogue::Catalogue random, const binType binType, const double thetaMin, const double thetaMax, const double binSize, const double shift, const CoordUnits angularUnits=_radians_, function<double(double)> angularWeight=nullptr)
	: TwoPointCorrelation1D(data, random) { m_twoPType = _1D_angular_; set_parameters(binType, thetaMin, thetaMax, binSize, shift, angularUnits, angularWeight); }

      /**
       *  @brief default destructor
       *  @return none
       */
      ~TwoPointCorrelation1D_angular () = default;

      ///@}
      

      /**
       *  @name Member functions to set the binning parameters
       */
      ///@{
      
      /**
       *  @brief set the binning parameters
       *  @param binType binning type: 0 &rarr; linear; 1 &rarr;
       *  logarithmic
       *  @param thetaMin minimum angular separation used to count the pairs
       *  @param thetaMax maximum angular separation used to count the pairs
       *  @param nbins number of bins
       *  @param shift shift parameter, i.e. the radial shift is
       *  binSize*shift
       *  @param angularUnits angular units
       *  @param angularWeight angular weight function
       *  @return none
       */
      void set_parameters (const binType binType, const double thetaMin, const double thetaMax, const int nbins, const double shift, const CoordUnits angularUnits=_radians_, function<double(double)> angularWeight=nullptr);

      /**
       *  @brief set the binning parameters
       *  @param binType binning type: 0 &rarr; linear; 1 &rarr;
       *  logarithmic
       *  @param thetaMin minimum angular separation used to count the pairs
       *  @param thetaMax maximum angular separation used to count the pairs
       *  @param binSize the bin size
       *  @param shift shift parameter, i.e. the radial shift is
       *  binSize*shift
       *  @param angularUnits angular units
       *  @param angularWeight angular weight function
       *  @return none
       */
      void set_parameters (const binType binType, const double thetaMin, const double thetaMax, const double binSize, const double shift, const CoordUnits angularUnits=_radians_, function<double(double)> angularWeight=nullptr);

      ///@}


      /**
       *  @name Member functions to count the number of pairs and measure the two-point correlation function
       */
      ///@{

      /**
       *  @brief measure the angular the two-point correlation
       *  function, w(theta)
       *  @param dir_output_pairs output directory used to store the
       *  number of pairs
       *  @param dir_input_pairs vector of input directories used to
       *  @param errorType the type of error to be computed
       *  store the number of pairs (if the pairs are read from files)
       *  @param dir_output_ResampleXi output directory used to store xi from
       *  resampled catalogues
       *  @param nMocks number of resampling to be generate with bootsrap*
       *  @param count_dd 1 &rarr; count the number of data-data
       *  opairs; 0 &rarr; read the number of data-data pairs from
       *  file
       *  @param count_rr 1 &rarr; count the number of random-random
       *  opairs; 0 &rarr; read the number of random-random pairs from
       *  file
       *  @param count_dr 1 &rarr; count the number of data-random
       *  opairs; 0 &rarr; read the number of data-random pairs from
       *  file
       *  @param tcount 1 &rarr; activate the time counter; 0 &rarr;
       *  don't activate the time counter; 
       *  @return none
       */
      void measure (const ErrorType errorType=ErrorType::_Poisson_, const string dir_output_pairs=par::defaultString, const vector<string> dir_input_pairs={},  const string dir_output_ResampleXi=par::defaultString, const int nMocks=0, const int count_dd=1, const int count_rr=1, const int count_dr=1, const bool tcount=1) override;
    
      ///@}

    
      /**
       *  @name Input/Output methods
       */
      ///@{

      /**
       *  @brief read the angular two-point correlation function
       *  @param dir input directory
       *  @param file input file
       *  @return none
       */
      void read (const string dir, const string file) override;

      /**
       *  @brief write the angular two-point correlation function
       *  @param dir output directory
       *  @param file output file
       *  @param rank cpu index (for MPI usage)
       *  @return none
       */
      void write (const string dir=par::defaultString, const string file=par::defaultString, const int rank=0) const override;
    
      ///@}

    };
  }
}

#endif
