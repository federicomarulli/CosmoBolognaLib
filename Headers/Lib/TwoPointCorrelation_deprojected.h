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
 *  @file Headers/Lib/TwoPointCorrelation_deprojected.h
 *
 *  @brief The class TwoPointCorrelation_deprojected
 *
 *  This file defines the interface of the class
 *  TwoPointCorrelation_deprojected, used to measure the deprojected
 *  two-point correlation function
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#ifndef __TWOPOINTDEPROJ__
#define __TWOPOINTDEPROJ__


#include "TwoPointCorrelation_projected.h"


// ===================================================================================================


namespace cosmobl {

  namespace twopt {
  
    /**
     *  @class TwoPointCorrelation_deprojected TwoPointCorrelation_deprojected.h
     *  "Headers/Lib/TwoPointCorrelation_deprojected.h"
     *
     *  @brief The class TwoPointCorrelation_deprojected
     *
     *  This class is used to handle objects of type <EM>
     *  TwoPointCorrelation_deprojected </EM>. It is used to measure
     *  the deprojected two-point correlation function,
     *  \f$\xi(r)=-\frac{1}{\pi}\int^{r_{\rm max}}_r dr_p'
     *  \frac{dw_p(r_p')/dr_p}{\sqrt{r_p'^2-r^2}}\f$.
     */
    class TwoPointCorrelation_deprojected : public TwoPointCorrelation_projected {
    
    protected:
      /**
       *  @brief measure deprojected correlation function
       *
       *  @param rp projected separation
       *
       *  @param xi 2d cartesian 2pcf
       *
       *  @param error_xi error on the 2d cartesian 2pcf
       *
       *  @return pointer to an object of type Data
       */
      shared_ptr<Data> DeProjectedTwoP(const vector<double> rp, const vector<double> xi, const vector<double> error_xi) override;

      /**
       *  @brief measure the deprojected two-point correlation
       *  function, &xi;(r) with Poisson error
       *
       *  @param dir_output_pairs output directory used to store the
       *  number of pairs
       *
       *  @param dir_input_pairs vector of input directories used to store the
       *  number of pairs (if the pairs are read from files)
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
      void measurePoisson (const string dir_output_pairs= par::defaultString, const vector<string> dir_input_pairs={}, const int count_dd=1, const int count_rr=1, const int count_dr=1, const bool tcount=1) override;

      /**
       *  @brief measure the deprojected two-point correlation
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
      void measureJackknife (const string dir_output_pairs= par::defaultString, const vector<string> dir_input_pairs={}, const string dir_output_ResampleXi= par::defaultString, const int count_dd=1, const int count_rr=1, const int count_dr=1, const bool tcount=1) override;


      /**
       *  @brief measure the deprojected two-point correlation
       *  function, &xi;(r), estimate the covariance with Bootstrap
       *  resampling
       *
       *  @param nMocks number of mocks to be generated with
       *  bootstrap resampling
       *
       *  @param dir_output_pairs output directory used to store the
       *  number of pairs
       *
       *  @param dir_input_pairs vector of input directories used to store the
       *  number of pairs (if the pairs are read from files)
       *
       *  @param dir_output_ResampleXi output directory used to store the
       *  Bootstrap resampling Xi, with Poisson error
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
      void measureBootstrap (const int nMocks, const string dir_output_pairs, const vector<string> dir_input_pairs={}, const string dir_output_ResampleXi = par::defaultString, const int count_dd=1, const int count_rr=1, const int count_dr=1, const bool tcount=1) override;

      /**
       *  @brief measure the jackknife resampling of the two-point
       *  correlation function, &xi;(r)
       *
       *  @param dd vector of data-data pairs, divider per regions
       *
       *  @param rr vector of random-random pairs, divider per regions
       *
       *  @return none
       */
      vector<shared_ptr<Data> > XiJackknife (const vector<shared_ptr<pairs::Pair> > dd, const vector<shared_ptr<pairs::Pair> > rr) override;

      /**
       *  @brief measure the jackknife resampling of the two-point
       *  correlation function, &xi;(r)
       *
       *  @param dd vector of data-data pairs, divider per regions
       *
       *  @param rr vector of random-random pairs, divider per regions
       *
       *  @param dr vector of random-random pairs, divider per regions   *
       *
       *  @return none
       */
      vector<shared_ptr<Data> > XiJackknife (const vector<shared_ptr<pairs::Pair> > dd, const vector<shared_ptr<pairs::Pair> > rr, const vector<shared_ptr<pairs::Pair> > dr) override;

      /**
       *  @brief measure the bootstrap resampling of the two-point
       *  correlation function, &xi;(r)
       *
       *  @param nMocks number of bootstrap resampling
       *
       *  @param dd vector of data-data pairs, divider per regions
       *
       *  @param rr vector of random-random pairs, divider per regions      
       *
       *  @return none
       */
       vector<shared_ptr<Data> > XiBootstrap (const int nMocks, const vector<shared_ptr<pairs::Pair> > dd, const vector<shared_ptr<pairs::Pair> > rr) override;

      /**
       *  @brief measure the bootstrap resampling of the two-point
       *  correlation function, &xi;(r)
       *
       *  @param nMocks number of bootstrap resampling
       *
       *  @param dd vector of data-data pairs, divider per regions
       *
       *  @param rr vector of random-random pairs, divider per regions 
       *
       *  @param dr vector of random-random pairs, divider per regions  
       *
       *  @return none
       */
      vector<shared_ptr<Data> > XiBootstrap (const int nMocks, const vector<shared_ptr<pairs::Pair> > dd, const vector<shared_ptr<pairs::Pair> > rr, const vector<shared_ptr<pairs::Pair> > dr) override;

    public:
    
      /**
       *  @name Constructors/destructors
       */
      ///@{

      /**
       *  @brief default constructor
       *  @return object of class TwoPointCorrelation_deprojected
       */
      TwoPointCorrelation_deprojected () { m_twoPType = _1D_deprojected_; }

      /**
       *  @brief constructor
       *  @param data object of class Catalogue containing the input
       *  catalogue
       *  @param random of class Catalogue containing the random data
       *  catalogue
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
       *  @return object of class TwoPointCorrelation_deprojected
       */
      TwoPointCorrelation_deprojected (catalogue::Catalogue data, catalogue::Catalogue random, const double rpMin, const double rpMax, const int nbins_rp, const double shift_rp, const double piMin, const double piMax, const int nbins_pi, const double shift_pi, const double piMax_integral, const CoordUnits angularUnits=_radians_, function<double(double)> angularWeight=nullptr)
	: TwoPointCorrelation_projected(data, random, _logarithmic_, rpMin, rpMax, nbins_rp, shift_rp, piMin, piMax, nbins_pi, shift_pi, piMax_integral, angularUnits, angularWeight)
	{ m_twoPType = _1D_deprojected_; }
      
      /**
       *  @brief constructor
       *  @param data object of class Catalogue containing the input
       *  catalogue
       *  @param random of class Catalogue containing the random data
       *  catalogue
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
       *  @return object of class TwoPointCorrelation_2D_deprojected
       */
      TwoPointCorrelation_deprojected (catalogue::Catalogue data, catalogue::Catalogue random, const double rpMin, const double rpMax, const double binSize_rp, const double shift_rp, const double piMin, const double piMax, const double binSize_pi, const double shift_pi, const double piMax_integral, const CoordUnits angularUnits=_radians_, function<double(double)> angularWeight=nullptr)
	: TwoPointCorrelation_projected(data, random, _logarithmic_, rpMin, rpMax, binSize_rp, shift_rp, piMin, piMax, binSize_pi, shift_pi, piMax_integral, angularUnits, angularWeight)
	{ m_twoPType = _1D_deprojected_; }
      
      /**
       *  @brief default destructor
       *  @return none
       */
      ~TwoPointCorrelation_deprojected () = default;

      ///@}

      /**
       *  @brief get the y coordinates
       *  @return the y coordinates
       */
      vector<double> yy () const 
      { cosmobl::ErrorMsg("Error in yy() of TwoPointCorrelation_deprojected.h!"); vector<double> vv; return vv; }

      /**
       *  @brief get the the binned correlation function 
       *  @return the binned correlation function 
       */
      vector<double> xi1D () const { return m_dataset->fx(); }

      /**
       *  @brief get the error on the binned correlation function
       *  function
       *  @return the error on the binned correlation function
       *  function
       */
      virtual vector<double> error1D () const { return m_dataset->error_fx(); }

      /**
       *  @brief get the the binned correlation function 
       *  @return the binned correlation function 
       */
      vector<vector<double> > xi2D () const 
      { cosmobl::ErrorMsg("Error in xi2D() of TwoPointCorrelation_deprojected.h!"); vector<vector<double> > vv; return vv; }

      /**
       *  @brief get the error on the binned correlation function
       *  function
       *  @return the error on the binned correlation function
       *  function
       */
      vector<vector<double> > error2D () const 
      { cosmobl::ErrorMsg("Error in error2D() of TwoPointCorrelation_deprojected.h!"); vector<vector<double> > vv; return vv; }
      
      /**
       *  @name Member functions to count the number of pairs and measure the two-point correlation function
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
       void measure (const ErrorType errType=ErrorType::_Poisson_, const string dir_output_pairs=par::defaultString, const vector<string> dir_input_pairs={}, const string dir_output_ResampleXi=par::defaultString, const int nMocks = 0., const int count_dd=1, const int count_rr=1, const int count_dr=1, const bool tcount=1) override;
              
      ///@}

    
      /**
       *  @name Input/Output methods
       */
      ///@{

      /**
       *  @brief read the deprojected two-point correlation function
       *  @param dir input directory
       *  @param file input file
       *  @return none
       */
      void read (const string dir, const string file) override;
      
      /**
       *  @brief write the deprojected two-point correlation function
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
