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
 *  @file Headers/Lib/TwoPointCorrelation2D_cartesian.h
 *
 *  @brief The class TwoPointCorrelation2D_cartesian
 *
 *  This file defines the interface of the class
 *  TwoPointCorrelation2D_cartesian, used to measure the 2D two-point
 *  correlation function in cartesian coordinates
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#ifndef __TWOPOINT2DCART__
#define __TWOPOINT2DCART__


#include "TwoPointCorrelation2D.h"


// ===================================================================================================


namespace cosmobl {

  namespace twopt {
    
    /**
     *  @class TwoPointCorrelation2D_cartesian TwoPointCorrelation2D_cartesian.h
     *  "Headers/Lib/TwoPointCorrelation2D_cartesian.h"
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
	 *  @brief measure the 2d cartesian two-point correlation
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
	 *  @brief measure the 2d cartesian two-point correlation
	 *  function, &xi;(r), estimate the covariance with Jackknife resampling
	 *
	 *  @param dir_output_pairs output directory used to store the
	 *  number of pairs
	 *
	 *  @param dir_input_pairs vector of input directories used to store the
	 *  number of pairs (if the pairs are read from files)
	 *
	 *  @param dir_output_ResampleXi output directory used to store the
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
	void measureJackknife (const string dir_output_pairs=par::defaultString, const vector<string> dir_input_pairs={}, const string dir_output_ResampleXi=par::defaultString, const int count_dd=1, const int count_rr=1, const int count_dr=1, const bool tcount=1) override;

	/**
	 *  @brief measure the 2d cartesian two-point correlation
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
	void measureBootstrap (const int nMocks, const string dir_output_pairs=par::defaultString, const vector<string> dir_input_pairs={}, const string dir_output_ResampleXi=par::defaultString, const int count_dd=1, const int count_rr=1, const int count_dr=1, const bool tcount=1) override;

      public:

	/**
	 *  @name Constructors/destructors
	 */
	///@{

	/**
	 *  @brief default constructor
	 *  @return object of class TwoPointCorrelation2D_cartesian
	 */
	TwoPointCorrelation2D_cartesian () { m_twoPType = _2D_Cartesian_; }

	/**
	 *  @brief constructor
	 *  @param data object of class Catalogue containing the input
	 *  catalogue
	 *  @param random of class Catalogue containing the random data
	 *  catalogue
	 *  @param binType_rp binning type in perpendicular separations:
	 *  0 &rarr; linear; 1 &rarr; logarithmic
	 *  @param rpMin minimum perpendicular separation used to count
	 *  the pairs
	 *  @param rpMax maximum perpendicular separation used to count
	 *  the pairs
	 *  @param nbins_rp number of bins in the perpendicular
	 *  separation
	 *  @param shift_rp shift parameter in the perpendicular
	 *  separation, i.e. the radial shift is binSize*shift
	 *  @param binType_pi binning type in parallel separations:
	 *  0 &rarr; linear; 1 &rarr; logarithmic
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
	 *  @return object of class TwoPointCorrelation2D_cartesian
	 */
	TwoPointCorrelation2D_cartesian (catalogue::Catalogue data, catalogue::Catalogue random, const binType binType_rp, const double rpMin, const double rpMax, const int nbins_rp, const double shift_rp, const binType binType_pi, const double piMin, const double piMax, const int nbins_pi, const double shift_pi, const CoordUnits angularUnits=_radians_, function<double(double)> angularWeight=nullptr)
	  : TwoPointCorrelation2D(data, random) { m_twoPType = _2D_Cartesian_; set_parameters(binType_rp, rpMin, rpMax, nbins_rp, shift_rp, binType_pi, piMin, piMax, nbins_pi, shift_pi, angularUnits, angularWeight); }

	/**
	 *  @brief constructor
	 *  @param data object of class Catalogue containing the input
	 *  catalogue
	 *  @param random of class Catalogue containing the random data
	 *  catalogue
	 *  @param binType_rp binning type in perpendicular separations:
	 *  0 &rarr; linear; 1 &rarr; logarithmic
	 *  @param rpMin minimum perpendicular separation used to count
	 *  the pairs
	 *  @param rpMax maximum perpendicular separation used to count
	 *  the pairs
	 *  @param binSize_rp bin size in the perpendicular separation
	 *  @param shift_rp shift parameter in the perpendicular
	 *  separation, i.e. the radial shift is binSize*shift
	 *  @param binType_pi binning type in parallel separations:
	 *  0 &rarr; linear; 1 &rarr; logarithmic
	 *  @param piMin minimum parallel separation used to count
	 *  the pairs
	 *  @param piMax maximum parallel separation used to count
	 *  the pairs
	 *  @param binSize_pi bin size in the parallel separation
	 *  @param shift_pi shift parameter in the parallel
	 *  separation, i.e. the radial shift is binSize*shift
	 *  @param angularUnits angular units
       *  @param angularWeight angular weight function
	 *  @return object of class TwoPointCorrelation2D_cartesian
	 */
	TwoPointCorrelation2D_cartesian (catalogue::Catalogue data, catalogue::Catalogue random, const binType binType_rp, const double rpMin, const double rpMax, const double binSize_rp, const double shift_rp, const binType binType_pi, const double piMin, const double piMax, const double binSize_pi, const double shift_pi, const CoordUnits angularUnits=_radians_, function<double(double)> angularWeight=nullptr)
	  : TwoPointCorrelation2D(data, random) { m_twoPType = _2D_Cartesian_; set_parameters(binType_rp, rpMin, rpMax, binSize_rp, shift_rp, binType_pi, piMin, piMax, binSize_pi, shift_pi, angularUnits, angularWeight); }

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
	 *  @param binType_rp binning type in perpendicular separations:
	 *  0 &rarr; linear; 1 &rarr; logarithmic
	 *  @param rpMin minimum perpendicular separation used to count
	 *  the pairs
	 *  @param rpMax maximum perpendicular separation used to count
	 *  the pairs
	 *  @param nbins_rp number of bins in the perpendicular
	 *  separation
	 *  @param shift_rp shift parameter in the perpendicular
	 *  separation, i.e. the radial shift is binSize*shift
	 *  @param binType_pi binning type in parallel separations:
	 *  0 &rarr; linear; 1 &rarr; logarithmic
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
	 *  @return object of class TwoPointCorrelation2D_cartesian
	 */
	void set_parameters (const binType binType_rp, const double rpMin, const double rpMax, const int nbins_rp, const double shift_rp, const binType binType_pi, const double piMin, const double piMax, const int nbins_pi, const double shift_pi, const CoordUnits angularUnits=_radians_, function<double(double)> angularWeight=nullptr);

	/**
	 *  @brief set the binning parameters
	 *  @param binType_rp binning type in perpendicular separations:
	 *  0 &rarr; linear; 1 &rarr; logarithmic
	 *  @param rpMin minimum perpendicular separation used to count
	 *  the pairs
	 *  @param rpMax maximum perpendicular separation used to count
	 *  the pairs
	 *  @param binSize_rp bin size in the perpendicular separation
	 *  @param shift_rp shift parameter in the perpendicular
	 *  separation, i.e. the radial shift is binSize*shift
	 *  @param binType_pi binning type in parallel separations:
	 *  0 &rarr; linear; 1 &rarr; logarithmic
	 *  @param piMin minimum parallel separation used to count
	 *  the pairs
	 *  @param piMax maximum parallel separation used to count
	 *  the pairs
	 *  @param binSize_pi bin size in the parallel separation
	 *  @param shift_pi shift parameter in the parallel
	 *  separation, i.e. the radial shift is binSize*shift
	 *  @param angularUnits angular units
       *  @param angularWeight angular weight function
	 *  @return none
	 */
	void set_parameters (const binType binType_rp, const double rpMin, const double rpMax, const double binSize_rp, const double shift_rp, const binType binType_pi, const double piMin, const double piMax, const double binSize_pi, const double shift_pi, const CoordUnits angularUnits=_radians_, function<double(double)> angularWeight=nullptr);

	///@}


	/**
	 *  @name Methods to set the binning parameters
	 */
	///@{

	/**
	 *  @brief measure the monopole of the two-point correlation
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
