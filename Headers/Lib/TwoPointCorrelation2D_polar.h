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
 *  @file Headers/Lib/TwoPointCorrelation2D_polar.h
 *
 *  @brief The class TwoPointCorrelation2D_polar
 *
 *  This file defines the interface of the class
 *  TwoPointCorrelation2D_polar, used to measure the 2D two-point
 *  correlation function in polar coordinates
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#ifndef __TWOPOINT2DPOL__
#define __TWOPOINT2DPOL__


#include "TwoPointCorrelation2D.h"


// ===================================================================================================


namespace cosmobl {

  namespace twopt {
    
    /**
     *  @class TwoPointCorrelation2D_polar TwoPointCorrelation2D_polar.h
     *  "Headers/Lib/TwoPointCorrelation2D_polar.h"
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
       *  @brief measure the 2d polar two-point correlation
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
       *  @brief measure the 2d polar two-point correlation
       *  function, &xi;(r), estimate the covariance with Jackknife resampling
       *
       *  @param dir_output_pairs output directory used to store the
       *  number of pairs
       *
       *  @param dir_input_pairs vector of input directories used to store the
       *  number of pairs (if the pairs are read from files)
       *
       *  @param dir_output_JackknifeXi output directory used to store the
       *  Jackknife resampling Xi, with Poisson error
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
      void measureJackknife (const string dir_output_pairs, const vector<string> dir_input_pairs={}, const string dir_output_JackknifeXi = "NULL", const int count_dd=1, const int count_rr=1, const int count_dr=1, const bool tcount=1) override;

      /**
       *  @brief measure the 2d polar two-point correlation
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
       *  @param dir_output_BootstrapXi output directory used to store the
       *  Jackknife resampling Xi, with Poisson error
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
      void measureBootstrap (const int nMocks, const string dir_output_pairs, const vector<string> dir_input_pairs={}, const string dir_output_BootstrapXi = "NULL", const int count_dd=1, const int count_rr=1, const int count_dr=1, const bool tcount=1) override;

    public:

      /**
       *  @name Constructors/destructors
       */
      ///@{

      /**
       *  @brief default constructor
       *  @return object of class TwoPointCorrelation2D_polar
       */
      TwoPointCorrelation2D_polar () { m_twoPType = _2D_polar_; }

      /**
       *  @brief constructor
       *  @param data object of class Catalogue containing the input
       *  catalogue
       *  @param random of class Catalogue containing the random data
       *  catalogue
       *  @param binType_rad binning type in absolute separations:
       *  0 &rarr; linear; 1 &rarr; logarithmic
       *  @param rMin minimum absolute separation used to count
       *  the pairs
       *  @param rMax maximum absolute separation used to count
       *  the pairs
       *  @param nbins_rad number of bins in the absolute
       *  separation
       *  @param shift_rad shift parameter in the absolute
       *  separation, i.e. the radial shift is binSize*shift
       *  @param binType_mu binning type in angular separations:
       *  0 &rarr; linear; 1 &rarr; logarithmic
       *  @param muMin minimum angular used to count the pairs
       *  @param muMax maximum angular used to count the pairs
       *  @param nbins_mu number of bins in the angular
       *  separation
       *  @param shift_mu shift parameter in the angular
       *  separation, i.e. the radial shift is binSize*shift
       *  @param angularUnits angular units
       *  @param angularWeight angular weight function
       *  @return object of class TwoPointCorrelation2D_polar
       */
      TwoPointCorrelation2D_polar (catalogue::Catalogue data, catalogue::Catalogue random, const binType binType_rad, const double rMin, const double rMax, const int nbins_rad, const double shift_rad, const binType binType_mu, const double muMin, const double muMax, const int nbins_mu, const double shift_mu, const CoordUnits angularUnits=_radians_, function<double(double)> angularWeight=nullptr)
	: TwoPointCorrelation2D(data, random) { m_twoPType = _2D_polar_; set_parameters(binType_rad, rMin, rMax, nbins_rad, shift_rad, binType_mu, muMin, muMax, nbins_mu, shift_mu, angularUnits, angularWeight); }

      /**
       *  @brief constructor
       *  @param data object of class Catalogue containing the input
       *  catalogue
       *  @param random of class Catalogue containing the random data
       *  catalogue
       *  @param binType_rad binning type in absolute separations:
       *  0 &rarr; linear; 1 &rarr; logarithmic
       *  @param rMin minimum absolute separation used to count
       *  the pairs
       *  @param rMax maximum absolute separation used to count
       *  the pairs
       *  @param binSize_rad bin size in the absolute separation
       *  @param shift_rad shift parameter in the absolute
       *  separation, i.e. the radial shift is binSize*shift
       *  @param binType_mu binning type in angular separations:
       *  0 &rarr; linear; 1 &rarr; logarithmic
       *  @param muMin minimum angular separation used to count
       *  the pairs
       *  @param muMax maximum angular separation used to count
       *  the pairs
       *  @param binSize_mu bin size in the angular separation
       *  @param shift_mu shift parameter in the angular
       *  separation, i.e. the radial shift is binSize*shift
       *  @param angularUnits angular units
       *  @param angularWeight angular weight function
       *  @return object of class TwoPointCorrelation2D_polar
       */
      TwoPointCorrelation2D_polar (catalogue::Catalogue data, catalogue::Catalogue random, const binType binType_rad, const double rMin, const double rMax, const double binSize_rad, const double shift_rad, const binType binType_mu, const double muMin, const double muMax, const double binSize_mu, const double shift_mu, const CoordUnits angularUnits=_radians_, function<double(double)> angularWeight=nullptr)
	: TwoPointCorrelation2D(data, random) { m_twoPType = _2D_polar_; set_parameters(binType_rad, rMin, rMax, binSize_rad, shift_rad, binType_mu, muMin, muMax, binSize_mu, shift_mu, angularUnits, angularWeight); }

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
       *  @param binType_rad binning type in absolute separations:
       *  0 &rarr; linear; 1 &rarr; logarithmic
       *  @param rMin minimum absolute separation used to count
       *  the pairs
       *  @param rMax maximum absolute separation used to count
       *  the pairs
       *  @param nbins_rad number of bins in the absolute
       *  separation
       *  @param shift_rad shift parameter in the absolute
       *  separation, i.e. the radial shift is binSize*shift
       *  @param binType_mu binning type in angular separations:
       *  0 &rarr; linear; 1 &rarr; logarithmic
       *  @param muMin minimum angular used to count the pairs
       *  @param muMax maximum angular used to count the pairs
       *  @param nbins_mu number of bins in the angular
       *  separation
       *  @param shift_mu shift parameter in the angular
       *  separation, i.e. the radial shift is binSize*shift
       *  @param angularUnits angular units
       *  @param angularWeight angular weight function
       *  @return none
       */
      void set_parameters (const binType binType_rad, const double rMin, const double rMax, const int nbins_rad, const double shift_rad, const binType binType_mu, const double muMin, const double muMax, const int nbins_mu, const double shift_mu, const CoordUnits angularUnits=_radians_, function<double(double)> angularWeight=nullptr);

      /**
       *  @brief set the binning parameters
       *  @param binType_rad binning type in absolute separations:
       *  0 &rarr; linear; 1 &rarr; logarithmic
       *  @param rMin minimum absolute separation used to count
       *  the pairs
       *  @param rMax maximum absolute separation used to count
       *  the pairs
       *  @param binSize_rad bin size in the absolute separation
       *  @param shift_rad shift parameter in the absolute
       *  separation, i.e. the radial shift is binSize*shift
       *  @param binType_mu binning type in angular separations:
       *  0 &rarr; linear; 1 &rarr; logarithmic
       *  @param muMin minimum angular separation used to count
       *  the pairs
       *  @param muMax maximum angular separation used to count
       *  the pairs
       *  @param binSize_mu bin size in the angular separation
       *  @param shift_mu shift parameter in the angular
       *  separation, i.e. the radial shift is binSize*shift
       *  @param angularUnits angular units
       *  @param angularWeight angular weight function
       *  @return none
       */
      void set_parameters (const binType binType_rad, const double rMin, const double rMax, const double binSize_rad, const double shift_rad, const binType binType_mu, const double muMin, const double muMax, const double binSize_mu, const double shift_mu, const CoordUnits angularUnits=_radians_, function<double(double)> angularWeight=nullptr);

      ///@}


      /**
       *  @name Methods to set the binning parameters
       */
      ///@{

      /**
       *  @brief measure the 2d polar two-point correlation
       *  function, &xi;(r)
       *
       *  @param errType type of &xi;(r) error
       *  
       *  @param dir_output_pairs output directory used to store the
       *  number of pairs
       *
       *  @param dir_input_pairs vector of input directories used to store the
       *  number of pairs (if the pairs are read from files)
       *
       *  @param dir_output_ResampleXi output directory of the resampled &xi;(r)
       *
       *  @param nMocks number of resampling for bootstrap
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
      void measure (const ErrorType errType=ErrorType::_Poisson_, const string dir_output_pairs=par::defaultString, const vector<string> dir_input_pairs={}, const string dir_output_ResampleXi=par::defaultString, const int nMocks=0, const int count_dd=1, const int count_rr=1, const int count_dr=1, const bool tcount=1) override;

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
      void read (const string dir, const string file) override;     

      /**
       *  @brief write the 2D two-point correlation function
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
