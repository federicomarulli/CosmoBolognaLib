/********************************************************************
 *  Copyright (C) 2010 by Federico Marulli and Carlo Giocoli        *
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
 *  @file Headers/TwoPointCorrelationCross1D.h
 *
 *  @brief The class TwoPointCorrelationCross1D
 *
 *  This file defines the interface of the class
 *  TwoPointCorrelationCross1D, used to measure cross 1D two-point
 *  correlation functions
 *
 *  @authors Federico Marulli, Carlo Giocoli
 *
 *  @authors federico.marulli3@unbo.it, carlo.giocoli@unibo.it
 */

#ifndef __TWOPOINTCROSS1D__
#define __TWOPOINTCROSS1D__


#include "TwoPointCorrelationCross.h"
#include "TwoPointCorrelation1D.h"


// ===================================================================================================


namespace cbl {

  namespace measure {
  
    namespace twopt {
  
      /**
       *  @class TwoPointCorrelationCross1D TwoPointCorrelationCross1D.h
       *  "Headers/TwoPointCorrelationCross1D.h"
       *
       *  @brief The class TwoPointCorrelationCross1D
       *
       *  This class is used to handle objects of type <EM>
       *  TwoPointCorrelationCross1D </EM>. It is used to measure the
       *  1D two-point correlation function.
       */
      class TwoPointCorrelationCross1D : public virtual TwoPointCorrelationCross, public virtual TwoPointCorrelation1D {
	
      protected :

	/**
	 *  @name Member functions to measure the two-point cross correlation function
	 */
	///@{
	
	/**
	 *  @brief get a dataset containing the two-point cross
	 *  correlation function measured with the Szapudi-Szalay
	 *  estimator, and its Poisson errors
	 *  
	 *  @param d1d2 pointer to an object of type Pair containing
	 *  the data1-data2 pairs
	 *
	 *  @param rr pointer to an object of type Pair containing the
	 *  random-random pairs
	 *
	 *  @param d1r pointer to an object of type Pair containing the
	 *  data1-random pairs
	 *
	 *  @param d2r pointer to an object of type Pair containing the
	 *  data2-random pairs
	 *
	 *  @param nData1 number of objects in the first data catalogue
	 *
	 *  @param nData2 number of objects in the second data catalogue
	 *
	 *  @param nData1_weighted weighted number of objects in the
	 *  first data catalogue
	 *
	 *  @param nData2_weighted weighted number of objects in the
	 *  second data catalogue
	 *	 
	 *  @param nRandom number of objects in the random catalogue
	 *
	 *  @param nRandom_weighted weighted number of objects in the
	 *  random catalogue
	 *
	 *  @return pointer to an object of type Data
	 */
	std::shared_ptr<data::Data> correlation_SzapudiSzalayEstimator (const std::shared_ptr<pairs::Pair> d1d2, const std::shared_ptr<pairs::Pair> rr, const std::shared_ptr<pairs::Pair> d1r, const std::shared_ptr<pairs::Pair> d2r, const int nData1=0, const double nData1_weighted=0., const int nData2=0, const double nData2_weighted=0., const int nRandom=0, const double nRandom_weighted=0.) override;

 	///@}


      public:
    
	/**
	 *  @name Constructors/destructors
	 */
	///@{

	/**
	 *  @brief default constructor
	 *  @return object of class TwoPointCorrelationCross1D
	 */
	TwoPointCorrelationCross1D () { m_dataset = data::Data::Create(data::DataType::_1D_); }

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
	 *  @param compute_extra_info true &rarr; compute extra
	 *  information related to the pairs, such as the mean pair
	 *  separation and redshift
	 *  @param random_dilution_fraction fraction between the number
	 *  of objects in the diluted and original random samples, used
	 *  to improve performances in random-random pair counts
	 *  @return object of class TwoPointCorrelationCross1D
	 */
	TwoPointCorrelationCross1D (const catalogue::Catalogue data1, const catalogue::Catalogue data2, const catalogue::Catalogue random, const bool compute_extra_info=false, const double random_dilution_fraction=1.) 
	  : TwoPointCorrelationCross(data1, data2, random, compute_extra_info, random_dilution_fraction)
	  { m_dataset = (!compute_extra_info) ? data::Data::Create(data::DataType::_1D_) : data::Data::Create(data::DataType::_1D_extra_); }

	/**
	 *  @brief default destructor
	 *  @return none
	 */
	~TwoPointCorrelationCross1D () = default;
      
	///@}	

	/**
	 *  @name Member functions to count measure the cross two-point correlation function
	 */
	///@{

	/**
	 *  @brief measure the cross two-point correlation function
	 *
	 *  @param errorType type of &xi;(r) error
	 *  
	 *  @param dir_output_pairs output directory used to store the
	 *  number of pairs
	 *
	 *  @param dir_input_pairs vector of input directories used to
	 *  store the number of pairs (if the pairs are read from files)
	 *
	 *  @param dir_output_resample output directory of the
	 *  resampled &xi;(r)
	 *
	 *  @param nMocks number of resampling for bootstrap
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
	 *  @param estimator the estimator used to measure the two-point
	 *  correlation function
	 *
	 *  @return none
	 */
	virtual void measure (const ErrorType errorType=ErrorType::_Poisson_, const std::string dir_output_pairs=par::defaultString, const std::vector<std::string> dir_input_pairs={}, const std::string dir_output_resample=par::defaultString, const int nMocks=0, const bool count_d1d2=true, const bool count_rr=true, const bool count_d1r=true, const bool count_d2r=true, const bool tcount=true, const Estimator estimator=Estimator::_SzapudiSzalay_) = 0;
      
	///@}
	
      };
    }
  }
}

#endif
