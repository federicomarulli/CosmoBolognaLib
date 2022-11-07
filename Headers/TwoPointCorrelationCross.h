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
 *  @file Headers/TwoPointCorrelationCross.h
 *
 *  @brief The class TwoPointCorrelationCross
 *
 *  This file defines the interface of the class
 *  TwoPointCorrelationCross, used to measure the cross two-point
 *  correlation function
 *
 *  @author Federico Marulli, Carlo Giocoli
 *
 *  @author federico.marulli3@unibo.it, carlo.giocoli@unibo.it
 */

#ifndef __TWOPOINTCROSS__
#define __TWOPOINTCROSS__


#include "TwoPointCorrelation.h"


// ===================================================================================================


namespace cbl {

  namespace measure {
  
    namespace twopt {

      /**
       *  @class TwoPointCorrelationCross TwoPointCorrelationCross.h
       *  "Headers/TwoPointCorrelationCross.h"
       *
       *  @brief The class TwoPointCorrelationCross
       *
       *  This is the base class used to measure the two-point
       *  correlation function
       *
       */
      class TwoPointCorrelationCross : public virtual TwoPointCorrelation {

      protected :

	/**
	 *  @name Input second catalogue
	 */
	///@{
    
	/// input second data catalogue
	std::shared_ptr<catalogue::Catalogue> m_data2;
    
	///@}

	
	/**
	 *  @name Object pairs
	 */
	///@{
    
	/// number of data1-data2 pairs
	std::shared_ptr<pairs::Pair> m_d1d2;

	/// number of data1-random pairs
	std::shared_ptr<pairs::Pair> m_d1r;

	/// number of data2-random pairs
	std::shared_ptr<pairs::Pair> m_d2r;

	///@}

	
	/**
	 *  @name Member functions to get the private/protected members
	 */
	///@{
      
	/**
	 *  @brief get the protected member m_data2
	 *  @return the input second data catalogue
	 */
	std::shared_ptr<catalogue::Catalogue> data2 () const { return m_data2; }
	
	/**
	 *  @brief get the protected member m_d1d2
	 *  @return the number of data1-data2 pairs
	 */
	std::shared_ptr<pairs::Pair> d1d2 () const { return m_d1d2; }

	/**
	 *  @brief get the protected member m_d1r
	 *  @return the number of data1-random pairs
	 */
	std::shared_ptr<pairs::Pair> d1r () const { return m_d1r; }

	/**
	 *  @brief get the protected member m_d2r
	 *  @return the number of data2-random pairs
	 */
	std::shared_ptr<pairs::Pair> d2r () const { return m_d2r; }
	
	///@}

	
	/**
	 *  @name Member functions to count the number of pairs 
	 */
	///@{

	/**
	 *  @brief count the data-data, random-random and data-random
	 *  pairs, used to construct the estimator of the two-point
	 *  correlation function
	 *
	 *  @param type the type of two-point correlation function; it
	 *  can be: \_monopole\_, \_projected\_,
	 *  \_deprojected\_, \_multipoles\_, \_angular\_,
	 *  \_2D_Cartesian\_, \_2D_polar\_
	 *
	 *  @param dir_output_pairs output directory used to store the
	 *  number of pairs
	 *
	 *  @param dir_input_pairs vector of input directories used to
	 *  store the number of pairs (if the pairs are read from files)
	 *
	 *  @param count_d1d2 true &rarr; count the number of
	 *  data1-data2 pairs; false &rarr; read the number of
	 *  data-data pairs from file
	 *
	 *  @param count_rr true &rarr; count the number of
	 *  random-random pairs; false &rarr; read the number of
	 *  random-random pairs from file
	 *
	 *  @param count_d1r true &rarr; count the number of
	 *  data1-random pairs; false &rarr; read the number of
	 *  data1-random pairs from file
	 *
	 *  @param count_d2r true &rarr; count the number of
	 *  data2-random pairs; false &rarr; read the number of
	 *  data2-random pairs from file
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
	 *  
	 */
	void count_allPairs (const TwoPType type, const std::string dir_output_pairs=par::defaultString, const std::vector<std::string> dir_input_pairs={}, const bool count_d1d2=true, const bool count_rr=true, const bool count_d1r=true, const bool count_d2r=true, const bool tcount=true, const Estimator estimator=Estimator::_SzapudiSzalay_, const double fact=0.1);

	///@}

	
	/**
	 *  @name Member functions to compute the cross two-point correlation function
	 */
	///@{
	
	/**
	 *  @brief get a dataset containing the cross two-point
	 *  correlation function measured with the Szapudi-Szalay
	 *  estimator, and its Poisson errors
	 *  
	 *  @param d1d2 pointer to an object of type Pair containing the
	 *  data1-data2 pairs
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
	 *  @param nData1_weighted weighted number of objects in the
	 *  first data catalogue
	 *
	 *  @param nData2 number of objects in the second data catalogue
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
	virtual std::shared_ptr<data::Data> correlation_SzapudiSzalayEstimator (const std::shared_ptr<pairs::Pair> d1d2, const std::shared_ptr<pairs::Pair> rr, const std::shared_ptr<pairs::Pair> d1r, const std::shared_ptr<pairs::Pair> d2r, const int nData1, const double nData1_weighted, const int nData2, const double nData2_weighted, const int nRandom, const double nRandom_weighted) = 0;
	
	/**
	 *  @brief measure the two-point correlation function with
	 *  Poisson errors
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
	 *  @param estimator the estimator used to measure the two-point
	 *  correlation function
	 *
	 *  @param fact factor used to compute the cell size of the
	 *  chain mesh: it is multiplied by the maximum distance
	 *  considered for the couples and can be setted by the user
	 *  to optimize the count of the couples
	 */
	virtual void measurePoisson (const std::string dir_output_pairs = par::defaultString, const std::vector<std::string> dir_input_pairs={}, const bool count_d1d2=true, const bool count_rr=true, const bool count_d1r=true, const bool count_d2r=true, const bool tcount=true, const Estimator estimator=Estimator::_LandySzalay_, double fact=0.1)
	{ (void)dir_output_pairs; (void)dir_input_pairs; (void)count_d1d2; (void)count_rr; (void)count_d1r; (void)count_d2r; (void)tcount; (void)estimator; (void)fact; cbl::ErrorCBL("", "measurePoisson", "TwoPointCorrelation.h"); }
	
	///@}

	
      public:
    
	/**
	 *  @name Constructors/destructors
	 */
	///@{

	/**
	 *  @brief default constructor
	 *  Cross
	 */
	TwoPointCorrelationCross () = default;

	/**
	 *  @brief constructor
	 *
	 *  @param data1 object of class Catalogue containing the
	 *  first input catalogue
	 *
	 *  @param data2 object of class Catalogue containing the
	 *  second input catalogue
	 *
	 *  @param random object of class Catalogue containing the
	 *  random data catalogue
	 *
	 *  @param compute_extra_info true &rarr; compute extra
	 *  information related to the pairs, such as the mean pair
	 *  separation and redshift
	 *
	 *  @param random_dilution_fraction fraction between the number
	 *  of objects in the diluted and original random samples, used
	 *  to improve performances in random-random pair counts
	 *
	 *  
	 */
	TwoPointCorrelationCross (const catalogue::Catalogue data1, const catalogue::Catalogue data2, const catalogue::Catalogue random, const bool compute_extra_info=false, const double random_dilution_fraction=1.) 
	  : TwoPointCorrelation(data1, random, compute_extra_info, random_dilution_fraction), m_data2(std::make_shared<catalogue::Catalogue>(catalogue::Catalogue(std::move(data2)))) {}

	/**
	 *  @brief default destructor
	 *  
	 */
	virtual ~TwoPointCorrelationCross () = default;
	
	///@}


	/**
	 *  @name Member functions to estimate the errors and covariance matrices
	 */
	///@{

	/**
	 *  @brief the Poisson errors 
	 *
	 *  This function computes the Poisson errors associated to
	 *  the Szapudi&Szalay estimators of the cross two-point
	 *  correlation function.
	 *
	 *  @param estimator the estimator used to measure the cross
	 *  two-point correlation function
	 *
	 *  @param d1d2 number of data1-data2 pairs
	 *  @param rr number of random-random pairs
	 *  @param d1r number of data1-random pairs
	 *  @param d2r number of data2-random pairs
	 *  @param nData1 number of data points in the first catalogue
	 *  @param nData2 number of data points in the second
	 *  catalogue
	 *  @param nRandom number of random points
	 *
	 *  @return the Poisson error
	 *  
	 *  @warning This function currently works with only the natural
	 *  and Landy&Szalay estimators
	 */
	double PoissonError (const Estimator estimator, const double d1d2, const double rr, const double d1r, const double d2r, const int nData1, const int nData2, const int nRandom) const;
    
	///@}
	
      };
    }
  }
}

#endif
