/********************************************************************
 *  Copyright (C) 2017 by Federico Marulli and Carlo Giocoli        *
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
 *  @file Headers/TwoPointCorrelationCross1D_monopole.h
 *
 *  @brief The class TwoPointCorrelationCross1D_monopole
 *
 *  This file defines the interface of the class
 *  TwoPointCorrelationCross1D_monopole, used to measure the monopole
 *  of the cross two-point correlation function
 *
 *  @authors Federico Marulli, Carlo Giocoli
 *
 *  @authors federico.marulli3@unbo.it, carlo.giocoli@unibo.it
 */

#ifndef __TWOPOINTCROSSMON__
#define __TWOPOINTCROSSMON__

#include "TwoPointCorrelationCross1D.h"
#include "TwoPointCorrelation1D_monopole.h"


// ===================================================================================================


namespace cbl {

  namespace measure {
  
    namespace twopt {
  
      /**
       *  @class TwoPointCorrelationCross1D_monopole
       *  TwoPointCorrelationCross1D_monopole.h
       *  "Headers/TwoPointCorrelationCross1D_monopole.h"
       *
       *  @brief The class TwoPointCorrelationCross1D_monopole
       *
       *  This class is used to handle objects of type <EM>
       *  TwoPointCorrelationCross1D_monopole </EM>. It is used to
       *  measure the monopole of the cross two-point correlation
       *  function.
       * 
       */
      class TwoPointCorrelationCross1D_monopole : public virtual TwoPointCorrelationCross1D, public virtual TwoPointCorrelation1D_monopole {

      protected:
	
	/**
	 *  @brief measure the monopole of the cross two-point
	 *  correlation function with Poisson errors
	 *
	 *  @param dir_output_pairs output directory used to store the
	 *  number of pairs
	 *
	 *  @param dir_input_pairs vector of input directories used to
	 *  store the number of pairs (if the pairs are read from files)
	 *
	 *  @param count_d1d2 true &rarr; count the number of
	 *  data1-data2 pairs; false &rarr; read the number of
	 *  data1-data2 pairs from file
	 *
	 *  @param count_rr true &rarr; count the number of
	 *  random-random pairs; false &rarr; read the number of
	 *  random-random pairs from file
	 *
	 *  @param count_d1r true &rarr; count the number of
	 *  data1-random pairs; false &rarr; read the number of
	 *  data1-random pairs
	 *
	 *  @param count_d2r true &rarr; count the number of
	 *  data2-random pairs; false &rarr; read the number of
	 *  data2-random pairs
	 *
	 *  @param tcount true &rarr; activate the time counter; false
	 *  &rarr; no time counter
	 *
	 *  @param estimator the estimator used to measure the cross
	 *  two-point correlation function
	 *
	 *  @return none
	 */
	void measurePoisson (const std::string dir_output_pairs = par::defaultString, const std::vector<std::string> dir_input_pairs={}, const bool count_d1d2=true, const bool count_rr=true, const bool count_d1r=true, const bool count_d2r=true, const bool tcount=true, const Estimator estimator=Estimator::_SzapudiSzalay_) override;

      public:
    
	/**
	 *  @name Constructors/destructors
	 */
	///@{

	/**
	 *  @brief default constructor
	 *  @return object of class TwoPointCorrelationCross1D_monopole
	 */
	TwoPointCorrelationCross1D_monopole () { m_twoPType = TwoPType::_monopole_; }

	/**
	 *  @brief constructor
	 *
	 *  @param data1 object of class Catalogue containing the
	 *  first input catalogue
	 *
	 *  @param data2 object of class Catalogue containing the
	 *  second input catalogue
	 *
	 *  @param random of class Catalogue containing the random data
	 *  catalogue
	 *
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
	 *  @return object of class TwoPointCorrelationCross1D_monopole
	 */
	TwoPointCorrelationCross1D_monopole (const catalogue::Catalogue data1, const catalogue::Catalogue data2, const catalogue::Catalogue random, const BinType binType, const double rMin, const double rMax, const int nbins, const double shift, const CoordinateUnits angularUnits=CoordinateUnits::_radians_, std::function<double(double)> angularWeight=nullptr, const bool compute_extra_info=false, const double random_dilution_fraction=1.)
	  : TwoPointCorrelation(data1, random, compute_extra_info, random_dilution_fraction), TwoPointCorrelationCross(data1, data2, random, compute_extra_info, random_dilution_fraction), TwoPointCorrelation1D(data1, random, compute_extra_info, random_dilution_fraction), TwoPointCorrelationCross1D(data1, data2, random, compute_extra_info, random_dilution_fraction)
	  { m_twoPType = TwoPType::_monopole_; set_parameters(binType, rMin, rMax, nbins, shift, angularUnits, angularWeight, compute_extra_info); }
      
	/**
	 *  @brief constructor
	 *
	 *  @param data1 object of class Catalogue containing the
	 *  first input catalogue
	 *
	 *  @param data2 object of class Catalogue containing the
	 *  second input catalogue
	 *
	 *  @param random of class Catalogue containing the random data
	 *  catalogue
	 *
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
	 *  @return object of class TwoPointCorrelationCross1D_monopole
	 */
	TwoPointCorrelationCross1D_monopole (const catalogue::Catalogue data1, const catalogue::Catalogue data2, const catalogue::Catalogue random, const BinType binType, const double rMin, const double rMax, const double binSize, const double shift, const CoordinateUnits angularUnits=CoordinateUnits::_radians_, std::function<double(double)> angularWeight=nullptr, const bool compute_extra_info=false, const double random_dilution_fraction=1.)
	  : TwoPointCorrelation(data1, random, compute_extra_info, random_dilution_fraction), TwoPointCorrelationCross(data1, data2, random, compute_extra_info, random_dilution_fraction), TwoPointCorrelation1D(data1, random, compute_extra_info, random_dilution_fraction), TwoPointCorrelationCross1D(data1, data2, random, compute_extra_info, random_dilution_fraction)
	  { m_twoPType = TwoPType::_monopole_; set_parameters(binType, rMin, rMax, binSize, shift, angularUnits, angularWeight, compute_extra_info); }

	/**
	 *  @brief default destructor
	 *  @return none
	 */
	~TwoPointCorrelationCross1D_monopole () = default;

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
	 *  @return none
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
	 *  @return none
	 */
	void set_parameters (const BinType binType, const double rMin, const double rMax, const double binSize, const double shift, const CoordinateUnits angularUnits=CoordinateUnits::_radians_, std::function<double(double)> angularWeight=nullptr, const bool compute_extra_info=false);

	///@}

	
	/**
	 *  @name Member functions to count the number of pairs and measure the cross two-point correlation function
	 */
	///@{

	/**
	 *  @brief measure the monopole of the cross two-point
	 *  correlation function (with the direct estimator)
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
	 *  resampling correlation functions, with Poisson errors; if
	 *  an empty string (i.e. "" or "NULL") is provided, no output
	 *  will be stored
	 *
	 *  @param nMocks number of resampling used for bootstrap
	 *
	 *  @param count_d1d2 true &rarr; count the number of
	 *  data1-data2 pairs; false &rarr; read the number of
	 *  data1-data2 pairs from file
	 *
	 *  @param count_rr true &rarr; count the number of
	 *  random-random pairs; false &rarr; read the number of
	 *  random-random pairs
	 *
	 *  @param count_d1r true &rarr; count the number of
	 *  data1-random pairs; false &rarr; read the number of
	 *  data1-random pairs
	 *
	 *  @param count_d2r true &rarr; count the number of
	 *  data2-random pairs; false &rarr; read the number of
	 *  data2-random pairs
	 *
	 *  @param tcount true &rarr; activate the time counter; false
	 *  &rarr; no time counter
	 *
	 *  @param estimator the estimator used to measure monopole of
	 *  the cross two-point correlation function
	 *
	 *  @return none
	 */
	void measure (const ErrorType errorType=ErrorType::_Poisson_, const std::string dir_output_pairs=par::defaultString, const std::vector<std::string> dir_input_pairs={},  const std::string dir_output_resample=par::defaultString, const int nMocks=0, const bool count_d1d2=true, const bool count_rr=true, const bool count_d1r=true, const bool count_d2r=true, const bool tcount=true, const Estimator estimator=Estimator::_SzapudiSzalay_) override;
	
	///@}
		/**
	 *  @name Input/Output methods
	 */

	
	///@{

	/**
	 *  @brief read the monopole of the cross two-point correlation
	 *  @param dir input directory
	 *  @param file input file
	 *  @return none
	 */
	void read (const std::string dir, const std::string file) override;

	/**
	 *  @brief write the monopole of the cross two-point correlation
	 *  function
	 *  @param dir output directory
	 *  @param file output file
	 *  @param rank cpu index (for MPI usage)
	 *  @return none
	 */
	void write (const std::string dir=par::defaultString, const std::string file=par::defaultString, const int rank=0) const override;

	///@}
      };
    }
  }
}

#endif
	
