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
 *  @file Headers/Lib/TwoPointCorrelation1D_filtered.h
 *
 *  @brief The class TwoPointCorrelation1D_filtered
 *
 *  This file defines the interface of the class
 *  TwoPointCorrelation1D_filtered, used to measure the filtered 
 *  monopole of the two-point correlation function
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#ifndef __TWOPOINTFIL__
#define __TWOPOINTFIL__


#include "TwoPointCorrelation1D_monopole.h"


// ===================================================================================================


namespace cosmobl {

  namespace twopt {
  
    /**
     *  @class TwoPointCorrelation1D_filtered TwoPointCorrelation1D_filtered.h
     *  "Headers/Lib/TwoPointCorrelation1D_filtered.h"
     *
     *  @brief The class TwoPointCorrelation1D_filtered
     *
     *  This class is used to handle objects of type <EM>
     *  TwoPointCorrelation1D_filtered </EM>. It is used to measure
     *  the filtered monopole of the two-point correlation function,
     *  \f$w(r_c)=2\pi \int dr \xi(r) W(r,r_c) r^2\f$, where \f$W(x) =
     *  (2x)^2(1-x)(0.5-x)r_c^{-3}\f$, and \f$\xi(r)\f$ is the
     *  monopole of the two-point correlation function
     *   
     */
    class TwoPointCorrelation1D_filtered : public TwoPointCorrelation1D_monopole {

    public:
    
      /**
       *  @name Constructors/destructors
       */
      ///@{

      /**
       *  @brief default constructor
       *  @return object of class TwoPointCorrelation1D_monopole
       */
      TwoPointCorrelation1D_filtered () { m_twoPType = _1D_filtered_; }

      /**
       *  @brief constructor
       *  @param data object of class Catalogue containing the input
       *  catalogue
       *  @param random of class Catalogue containing the random data
       *  catalogue
       *  @param bintype binning type: 0 &rarr; linear; 1 &rarr;
       *  logarithmic
       *  @param rMin minimum separation used to count the pairs
       *  @param rMax maximum separation used to count the pairs
       *  @param nbins number of bins
       *  @param monopole_rMin minimum separation used to count the pairs
       *  @param monopole_rMax maximum separation used to count the pairs
       *  @param monopole_nbins number of bins
       *  @param monopole_shift shift parameter, i.e. the radial shift is
       *  binSize*shift*  @return object of class TwoPointCorrelation1D_monopole
       */
      TwoPointCorrelation1D_filtered (const catalogue::Catalogue data, const catalogue::Catalogue random, const binType bintype, const double rMin, const double rMax, const int nbins, const double monopole_rMin, const double monopole_rMax, const int monopole_nbins, const double monopole_shift)
	: TwoPointCorrelation1D_monopole(data, random, _linear_, monopole_rMin, monopole_rMax, monopole_nbins, monopole_shift)
      { m_twoPType = _1D_filtered_; set_parameters(bintype, rMin, rMax, nbins); }

       /**
       *  @brief constructor
       *  @param data object of class Catalogue containing the input
       *  catalogue
       *  @param random of class Catalogue containing the random data
       *  catalogue
       *  @param bintype binning type: 0 &rarr; linear; 1 &rarr;
       *  logarithmic
       *  @param rMin minimum separation used to filter the xi
       *  @param rMax maximum separation used to filter the xi
       *  @param nbins the number of bins
       *  @param monopole_rMin minimum separation used to count the pairs
       *  @param monopole_rMax maximum separation used to count the pairs
       *  @param monopole_binSize the size of the bin for the monopole
       *  @param monopole_shift shift parameter, i.e. the radial shift is
       *  binSize*shift
       *  @return object of class TwoPointCorrelation1D_filtered
       */
      TwoPointCorrelation1D_filtered (const catalogue::Catalogue data, const catalogue::Catalogue random, const binType bintype, const double rMin, const double rMax, const int nbins, const double monopole_rMin, const double monopole_rMax, const double monopole_binSize, const double monopole_shift)
	: TwoPointCorrelation1D_monopole(data, random, _linear_ , monopole_rMin, monopole_rMax, monopole_binSize, monopole_shift)
      { m_twoPType = _1D_filtered_; set_parameters(bintype, rMin, rMax, nbins); }     

      /**
       *  @brief constructor
       *  @param data object of class Catalogue containing the input
       *  catalogue
       *  @param random of class Catalogue containing the random data
       *  catalogue
       *  @param bintype binning type: 0 &rarr; linear; 1 &rarr;
       *  logarithmic
       *  @param rMin minimum separation used to count the pairs
       *  @param rMax maximum separation used to count the pairs
       *  @param binSize the bin size
       *  @param monopole_rMin minimum separation used to count the pairs
       *  @param monopole_rMax maximum separation used to count the pairs
       *  @param monopole_nbins number of bins
       *  @param monopole_shift shift parameter, i.e. the radial shift is
       *  binSize*shift*  @return object of class TwoPointCorrelation1D_monopole
       */
      TwoPointCorrelation1D_filtered (const catalogue::Catalogue data, const catalogue::Catalogue random, const binType bintype, const double rMin, const double rMax, const double binSize, const double monopole_rMin, const double monopole_rMax, const int monopole_nbins, const double monopole_shift)
	: TwoPointCorrelation1D_monopole(data, random,_linear_, monopole_rMin, monopole_rMax, monopole_nbins, monopole_shift)
      { m_twoPType = _1D_filtered_; set_parameters(bintype, rMin, rMax, binSize); }

       /**
       *  @brief constructor
       *  @param data object of class Catalogue containing the input
       *  catalogue
       *  @param random of class Catalogue containing the random data
       *  catalogue
       *  @param bintype binning type: 0 &rarr; linear; 1 &rarr;
       *  logarithmic
       *  @param rMin minimum separation used to filter the xi
       *  @param rMax maximum separation used to filter the xi
       *  @param binSize the number of bins
       *  logarithmic
       *  @param monopole_rMin minimum separation used to count the pairs
       *  @param monopole_rMax maximum separation used to count the pairs
       *  @param monopole_binSize the size of the bin for the monopole
       *  @param monopole_shift shift parameter, i.e. the radial shift is
       *  binSize*shift
       *  @return object of class TwoPointCorrelation1D_filtered
       */
      TwoPointCorrelation1D_filtered (const catalogue::Catalogue data, const catalogue::Catalogue random, const binType bintype, const double rMin, const double rMax, const double binSize, const double monopole_rMin, const double monopole_rMax, const double monopole_binSize, const double monopole_shift)
	: TwoPointCorrelation1D_monopole(data, random, _linear_, monopole_rMin, monopole_rMax, monopole_binSize, monopole_shift)
      { m_twoPType = _1D_filtered_; set_parameters(bintype, rMin, rMax, binSize); }     

      /**
       *  @brief default destructor
       *  @return none
       */
      ~TwoPointCorrelation1D_filtered () = default;

      ///@}

      
      /**
       *  @name Member functions to set the binning parameters
       */
      ///@{
      
      /**
       *  @brief set the binning parameters
       *  @param binType binning type: 0 &rarr; linear; 1 &rarr;
       *  logarithmic
       *  @param rMin minimum separation used to count the pairs
       *  @param rMax maximum separation used to count the pairs
       *  @param nbins number of bins
       *  @return none
       */
      void set_parameters (const binType binType, const double rMin, const double rMax, const int nbins);

      /**
       *  @brief set the binning parameters
       *  @param binType binning type: 0 &rarr; linear; 1 &rarr;
       *  logarithmic
       *  @param rMin minimum separation used to count the pairs
       *  @param rMax maximum separation used to count the pairs
       *  @param binSize the bin size
       *  @return none
       */
      void set_parameters (const binType binType, const double rMin, const double rMax, const double binSize);

      ///@}


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
      void measure (const ErrorType errType = ErrorType::_Poisson_, const string dir_output_pairs=par::defaultString, const vector<string> dir_input_pairs={}, const string dir_output_ResampleXi=par::defaultString, int nMocks = 0, const int count_dd=1, const int count_rr=1, const int count_dr=1, const bool tcount=1) override;
    
      ///@}

    
      /**
       *  @name Input/Output methods
       */
      ///@{

      /**
       *  @brief read the filtered two-point correlation function
       *  @param dir input directory
       *  @param file input file
       *  @return none
       */
      void read (const string dir, const string file) override;

      /**
       *  @brief write the filtered two-point correlation function
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
