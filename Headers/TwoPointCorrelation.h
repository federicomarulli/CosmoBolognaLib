/********************************************************************
 *  Copyright (C) 2010 by Federico Marulli                          *
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
 *  @file Headers/TwoPointCorrelation.h
 *
 *  @brief The class TwoPointCorrelation
 *
 *  This file defines the interface of the class TwoPointCorrelation,
 *  used to measure the two-point correlation function
 *
 *  @author Federico Marulli
 *
 *  @author federico.marulli3@unibo.it
 */

#ifndef __TWOPOINT__
#define __TWOPOINT__


#include "Measure.h"
#include "ChainMesh_Catalogue.h"
#include "Pair1D_extra.h"
#include "Pair2D_extra.h"


// ===================================================================================================


namespace cbl {

  namespace measure {
  
    /**
     *  @brief The namespace of the <B> two-point correlation function
     *  </B>
     *  
     *  The \e measure::twopt namespace contains all the functions and
     *  classes to measure the two-point correlation function
     */
    namespace twopt {

      /**
       * @enum TwoPType
       * @brief the two-point correlation function type
       */
      enum class TwoPType { 

	/// the angle-averaged two-point correlation function, i.e. the monopole, &xi;(r)
	_monopole_,

	/// the projected two-point correlation function, w(r<SUB>p</SUB>)
	_projected_,
    
	/// the deprojected two-point correlation function, &xi;(r)
	_deprojected_,
    
	/// the multipoles of the two-point correlation function, &xi;<SUB>i</SUB>(r), computed with the integrated estimator
	_multipoles_integrated_,
    
	/// the multipoles of the two-point correlation function, &xi;<SUB>i</SUB>(r), computed with the direct estimator
	_multipoles_direct_,

	/// the wedges of the two-point correlation function, &xi;<SUB>i</SUB>(r)
	_wedges_,
      
	/// filtered two-point correlation function
	_filtered_,

	/// angular two-point correlation function
	_angular_,
    
	/// 2D two-point correlation function in Cartesian coordinates, &xi;(r<SUB>p</SUB>,&pi;)
	_2D_Cartesian_,

	/// 2D two-point correlation function in polar coordinates, &xi;(r,&mu;)
	_2D_polar_
	  
      };

      /**
       * @brief return a vector containing the
       * TwoPType names
       * @return a vector containing the
       * TwoPType names
       */
      inline std::vector<std::string> TwoPTypeNames () {return {"monopole", "projected", "deprojected", "multipoles_integrated", "multipoles_direct", "wedges", "filtered", "angular", "2D_Cartesian", "2D_polar"}; }

      /**
       * @brief cast an enum of type TwoPType
       * from its index
       * @param twoPTypeIndex the twoPType index
       * @return object of class TwoPType
       */
      inline TwoPType TwoPTypeCast (const int twoPTypeIndex) {return castFromValue<TwoPType>(twoPTypeIndex);}

      /**
       * @brief cast an enum of type TwoPType
       * from its name
       * @param twoPTypeName the twoPType name
       * @return object of class TwoPType
       */
      inline TwoPType TwoPTypeCast (const std::string twoPTypeName)
      { return castFromName<TwoPType>(twoPTypeName, TwoPTypeNames()); }

      /**
       * @brief cast an enum of type TwoPType
       * from indeces
       * @param twoPTypeIndeces the twoPType indeces
       * @return vector of objects of class TwoPType
       */
      inline std::vector<TwoPType> TwoPTypeCast (const std::vector<int> twoPTypeIndeces)
      { return castFromValues<TwoPType>(twoPTypeIndeces); } 

      /**
       * @brief cast an enum of type TwoPType
       * from thier names
       * @param twoPTypeNames the twoPType names
       * @return vector of objects of class TwoPType
       */
      inline std::vector<TwoPType> TwoPTypeCast (const std::vector<std::string> twoPTypeNames)
      { return castFromNames<TwoPType>(twoPTypeNames, TwoPTypeNames()); }

      /**
       * @enum Estimator
       * @brief the two-point correlation estimator
       */
      enum class Estimator { 

	/// natural estimator
	_natural_,

	/// Landy&Szalay estimator
	_LandySzalay_,

	/// Szapudi&Szalay estimator
	_SzapudiSzalay_

      };

      /**
       * @brief return a vector containing the
       * Estimator names
       * @return a vector containing the
       * Estimator names
       */
      inline std::vector<std::string> EstimatorNames () {return {"natural", "LandySzalay", "SzapudiSzalay"}; }

      /**
       * @brief cast an enum of type Estimator
       * from its index
       * @param estimatorIndex the estimator index
       * @return object of class Estimator
       */
      inline Estimator EstimatorCast (const int estimatorIndex) {return castFromValue<Estimator>(estimatorIndex);}

      /**
       * @brief cast an enum of type Estimator
       * from its name
       * @param estimatorName the estimator name
       * @return object of class Estimator
       */
      inline Estimator EstimatorCast (const std::string estimatorName) {return castFromName<Estimator>(estimatorName, EstimatorNames());}

      /**
       * @brief cast an enum of type Estimator
       * from indeces
       * @param estimatorIndeces the estimator indeces
       * @return object of class Estimator
       */
      inline std::vector<Estimator> EstimatorCast (const std::vector<int> estimatorIndeces) {return castFromValues<Estimator>(estimatorIndeces);} 

      /**
       * @brief cast an enum of type Estimator
       * from thier names
       * @param estimatorNames the estimator names
       * @return vector of Estimator enums
       */
      inline std::vector<Estimator> EstimatorCast (const std::vector<std::string> estimatorNames) {return castFromNames<Estimator>(estimatorNames, EstimatorNames());}

      /**
       *  @class TwoPointCorrelation TwoPointCorrelation.h
       *  "Headers/TwoPointCorrelation.h"
       *
       *  @brief The class TwoPointCorrelation
       *
       *  This is the base class used to measure the two-point
       *  correlation function
       *
       */
      class TwoPointCorrelation : public Measure {

      protected :

	/**
	 *  @name Two-point correlation function data
	 */
	///@{
      
	/// two-point correlation function type
	TwoPType m_twoPType;

	///@}

      
	/**
	 *  @name Input and random catalogues
	 */
	///@{
    
	/// input data catalogue
	std::shared_ptr<catalogue::Catalogue> m_data;

	/// output data catalogue
	std::shared_ptr<catalogue::Catalogue> m_random;
    
	///@}


	/**
	 *  @name Object pairs
	 */
	///@{
    
	/// number of data-data pairs
	std::shared_ptr<pairs::Pair> m_dd;

	/// number of random-random pairs
	std::shared_ptr<pairs::Pair> m_rr;

	/// number of data-random pairs
	std::shared_ptr<pairs::Pair> m_dr;

	///@}

	/**
	 *  @name Object pairs for resampling
	 */
	///@{
    
	/// number of data-data pairs
	std::vector<std::shared_ptr<pairs::Pair>> m_dd_res;

	/// number of random-random pairs
	std::vector<std::shared_ptr<pairs::Pair>> m_rr_res;

	/// number of data-random pairs
	std::vector<std::shared_ptr<pairs::Pair>> m_dr_res;

	///@}
      
	/**
	 *  @name Other parameters
	 */
	///@{
	 
	/// true &rarr; compute extra information related to the pairs, such as the mean pair separation and redshift
	bool m_compute_extra_info;

	/// fraction between the number of random objects in the diluted and original samples, used to improve performances in random-random pair counts
	double m_random_dilution_fraction;
      
	///@}

      
	/**
	 *  @name Internal input/output member functions (customized in all the derived classes)
	 */
	///@{
    
	/**
	 *  @brief write the number of pairs
	 *  @param PP pointer to an object of class Pair
	 *  @param dir output directory
	 *  @param file output file
	 *  @return none, or an error message if the derived object does
	 *  not have this member
	 */
	virtual void write_pairs (const std::shared_ptr<pairs::Pair> PP, const std::string dir, const std::string file) const = 0;
 
	/**
	 *  @brief read the number of pairs
	 *  @param [out] PP pointer to an object of class Pair
	 *  @param [in] dir input directory
	 *  @param [in] file input file
	 *  @return none, or an error message if the derived object does
	 *  not have this member
	 */
	virtual void read_pairs (std::shared_ptr<pairs::Pair> PP, const std::vector<std::string> dir, const std::string file) const = 0;

	/**
	 *  @brief write the number of pairs
	 *  @param PP pointer to a vector of objects of class Pair
	 *  @param dir output directory
	 *  @param file output file
	 *  @return none, or an error message if the derived object does
	 *  not have this member
	 */
	virtual void write_pairs (const std::vector<std::shared_ptr<pairs::Pair> > PP, const std::string dir, const std::string file) const = 0;
 
	/**
	 *  @brief read the number of pairs
	 *  @param [out] PP pointer to a vector of objects of class Pair
	 *  @param [in] dir input directory
	 *  @param [in] file input file
	 *  @return none, or an error message if the derived object does
	 *  not have this member
	 */
	virtual void read_pairs (std::vector<std::shared_ptr<pairs::Pair> > PP, const std::vector<std::string> dir, const std::string file) const = 0;

	///@}
	
	/**
	 *  @name Member functions to reset the number of pairs 
	 */
	///@{
	
	/**
	 * @brief reset the pair counts
	 *
	 * @details set the pair counts to zero,
	 * by calling the cbl::Pairs::reset() method
	 * for the internal variables m_dd, m_rr, m_dr
	 */
	void resets ();
	
	///@}


	/**
	 *  @name Member functions to count the number of pairs 
	 */
	///@{

	/**
	 *  @brief count the number of pairs
	 *
	 *  @param cat1 pointer to an object of class Catalogue,
	 *  containing the first catalogue
	 *
	 *  @param ChM object of class ChainMesh_Catalogue, used to
	 *  construct the chain-mesh
	 *
	 *  @param pp pointer to an object of class Pair
	 *
	 *  @param cross true &rarr; count the number of cross pairs
	 *  (e.g. data-random pairs); false &rarr; count the number of
	 *  pairs of objects of the same type (e.g. data-data and
	 *  random-random pairs)
	 *
	 *  @param tcount true &rarr; activate the time counter; false
	 *  &rarr; no time counter
	 */
	void count_pairs (const std::shared_ptr<catalogue::Catalogue> cat1, const chainmesh::ChainMesh_Catalogue &ChM, std::shared_ptr<pairs::Pair> pp, const bool cross=true, const bool tcount=false);

	/**
	 *  @brief count the number of pairs, used for
	 *  Jackknife/Bootstrap methods
	 *
	 *  @param cat1 pointer to an object of class Catalogue,
	 *  containing the first catalogue
	 *
	 *  @param ChM object of class ChainMesh_Catalogue, used to
	 *  construct the chain-mesh
	 *
	 *  @param pp pointer to an object of class Pair
	 *
	 *  @param pp_regions vector containing pointers to object of
	 *  class Pair
	 *
	 *  @param cross true &rarr; count the number of cross pairs
	 *  (e.g. data-random pairs); false &rarr; count the number of
	 *  pairs of objects of the same type (e.g. data-data and
	 *  random-random pairs)
	 *
	 *  @param tcount true &rarr; activate the time counter; false
	 *  &rarr; no time counter
	 */
	void count_pairs_region (const std::shared_ptr<catalogue::Catalogue> cat1, const chainmesh::ChainMesh_Catalogue &ChM, std::shared_ptr<pairs::Pair> pp, std::vector< std::shared_ptr<pairs::Pair> > pp_regions, const bool cross=true, const bool tcount=false);
      
	/**
	 *  @brief count the number of pairs, used for
	 *  Jackknife/Bootstrap methods
	 *
	 *  @param cat1 pointer to an object of class Catalogue,
	 *  containing the first catalogue
	 *
	 *  @param ChM object of class ChainMesh_Catalogue, used to
	 *  construct the chain-mesh
	 *
	 *  @param pp pointer to an object of class Pair
	 *
	 *  @param pp_res vector containing pointers to object of
	 *  class Pair
	 *
	 *  @param weight weights of the region
	 *
	 *  @param cross true &rarr; count the number of cross pairs
	 *  (e.g. data-random pairs); false &rarr; count the number of
	 *  pairs of objects of the same type (e.g. data-data and
	 *  random-random pairs)
	 *
	 *  @param tcount true &rarr; activate the time counter; false
	 *  &rarr; no time counter
	 */
	void count_pairs_region_test (const std::shared_ptr<catalogue::Catalogue> cat1, const chainmesh::ChainMesh_Catalogue &ChM, std::shared_ptr<pairs::Pair> pp, std::vector< std::shared_ptr<pairs::Pair> > pp_res, const std::vector<double> weight, const bool cross=true, const bool tcount=false);
      
	/**
	 *  @brief count the number of pairs, used for
	 *  Jackknife/Bootstrap methods, 1D pairs
	 *
	 *  @param cat1 pointer to an object of class Catalogue,
	 *  containing the first catalogue
	 *
	 *  @param ChM object of class ChainMesh_Catalogue, used to
	 *  construct the chain-mesh
	 *
	 *  @param pp pointer to an object of class Pair
	 *
	 *  @param pp_res vector containing pointers to object of
	 *  class Pair
	 *
	 *  @param weight weights of the region
	 *
	 *  @param cross true &rarr; count the number of cross pairs
	 *  (e.g. data-random pairs); false &rarr; count the number of
	 *  pairs of objects of the same type (e.g. data-data and
	 *  random-random pairs)
	 *
	 *  @param tcount true &rarr; activate the time counter; false
	 *  &rarr; no time counter
	 */
	void count_pairs_region_test_1D (const std::shared_ptr<catalogue::Catalogue> cat1, const chainmesh::ChainMesh_Catalogue &ChM, std::shared_ptr<pairs::Pair> pp, std::vector< std::shared_ptr<pairs::Pair> > pp_res, const std::vector<double> weight, const bool cross=true, const bool tcount=false);
      
	/**
	 *  @brief count the number of pairs, used for
	 *  Jackknife/Bootstrap methods, 2D pairs
	 *
	 *  @param cat1 pointer to an object of class Catalogue,
	 *  containing the first catalogue
	 *
	 *  @param ChM object of class ChainMesh_Catalogue, used to
	 *  construct the chain-mesh
	 *
	 *  @param pp pointer to an object of class Pair
	 *
	 *  @param pp_res vector containing pointers to object of
	 *  class Pair
	 *
	 *  @param weight weights of the region
	 *
	 *  @param cross true &rarr; count the number of cross pairs
	 *  (e.g. data-random pairs); false &rarr; count the number of
	 *  pairs of objects of the same type (e.g. data-data and
	 *  random-random pairs)
	 *
	 *  @param tcount true &rarr; activate the time counter; false
	 *  &rarr; no time counter
	 */
	void count_pairs_region_test_2D (const std::shared_ptr<catalogue::Catalogue> cat1, const chainmesh::ChainMesh_Catalogue &ChM, std::shared_ptr<pairs::Pair> pp, std::vector< std::shared_ptr<pairs::Pair> > pp_res, const std::vector<double> weight, const bool cross=true, const bool tcount=false);

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
	 *  @param count_dd true &rarr; count the number of data-data
	 *  pairs; false &rarr; read the number of data-data pairs from
	 *  file
	 *
	 *  @param count_rr true &rarr; count the number of
	 *  random-random pairs; false &rarr; read the number of
	 *  random-random pairs from file
	 *
	 *  @param count_dr true &rarr; count the number of
	 *  data-random pairs; false &rarr; read the number of
	 *  data-random pairs from file
	 *
	 *  @param tcount true &rarr; activate the time counter; false
	 *  &rarr; no time counter
	 *
	 *  @param estimator the estimator used to measure the two-point
	 *  correlation function
	 */
	void count_allPairs (const TwoPType type, const std::string dir_output_pairs=par::defaultString, const std::vector<std::string> dir_input_pairs={}, const bool count_dd=true, const bool count_rr=true, const bool count_dr=true, const bool tcount=true, const Estimator estimator=Estimator::_LandySzalay_);

	/**
	 *  @brief count the data-data, random-random and data-random
	 *  pairs, used to construct the estimator of the two-point
	 *  correlation function
	 *
	 *  @param dd_regions data-data pairs in the sub-regions
	 *
	 *  @param rr_regions random-random pairs in the sub-regions
	 *
	 *  @param dr_regions data-random pairs in the sub-regions
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
	 *  @param count_dd true &rarr; count the number of data-data
	 *  pairs; false &rarr; read the number of data-data pairs from
	 *  file
	 *
	 *  @param count_rr true &rarr; count the number of
	 *  random-random pairs; false &rarr; read the number of
	 *  random-random pairs from file
	 *
	 *  @param count_dr true &rarr; count the number of
	 *  data-random pairs; false &rarr; read the number of
	 *  data-random pairs from file
	 *
	 *  @param tcount true &rarr; activate the time counter; false
	 *  &rarr; no time counter
	 *
	 *  @param estimator the estimator used to measure the two-point
	 *  correlation function
	 */
	void count_allPairs_region (std::vector<std::shared_ptr<pairs::Pair> > &dd_regions, std::vector<std::shared_ptr<pairs::Pair> > &rr_regions, std::vector<std::shared_ptr<pairs::Pair> > &dr_regions, const TwoPType type, const std::string dir_output_pairs=par::defaultString, const std::vector<std::string> dir_input_pairs={}, const bool count_dd=true, const bool count_rr=true, const bool count_dr=true, const bool tcount=true, const Estimator estimator=Estimator::_LandySzalay_);

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
	 *  @param weight region weights
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
	 *  @param count_dr true &rarr; count the number of
	 *  data-random pairs; false &rarr; read the number of
	 *  data-random pairs from file
	 *
	 *  @param tcount true &rarr; activate the time counter; false
	 *  &rarr; no time counter
	 *
	 *  @param estimator the estimator used to measure the two-point
	 *  correlation function
	 *
	 */
	void count_allPairs_region_test (const TwoPType type, const std::vector<double> weight, const std::string dir_output_pairs=par::defaultString, const std::vector<std::string> dir_input_pairs={}, const bool count_dd=true, const bool count_rr=true, const bool count_dr=true, const bool tcount=true, const Estimator estimator=Estimator::_LandySzalay_);

	///@}


	/**
	 *  @name Member functions to compute the two-point correlation function
	 */
	///@{

	/**
	 *  @brief get a dataset containing the two-point correlation
	 *  function measured with the natural estimator, and its
	 *  Poisson errors
	 *  
	 *  @param dd pointer to an object of type Pair containing the
	 *  data-data pairs
	 *
	 *  @param rr pointer to an object of type Pair containing the
	 *  random-random pairs
	 *
	 *  @param nData number of objects in the data catalogue
	 *
	 *  @param nData_weighted weighted number of objects in the
	 *  data catalogue
	 *
	 *  @param nRandom number of objects in the random catalogue
	 *
	 *  @param nRandom_weighted weighted number of objects in the
	 *  random catalogue
	 *
	 *  @return pointer to an object of type Data
	 */
	virtual std::shared_ptr<data::Data> correlation_NaturalEstimator (const std::shared_ptr<pairs::Pair> dd, const std::shared_ptr<pairs::Pair> rr, const int nData=0, const double nData_weighted=0., const int nRandom=0, const double nRandom_weighted=0.) = 0;

	/**
	 *  @brief get a dataset containing the two-point correlation
	 *  function measured with the Landy-Szalay estimator, and its
	 *  Poisson errors
	 *  
	 *  @param dd pointer to an object of type Pair containing the
	 *  data-data pairs
	 *
	 *  @param rr pointer to an object of type Pair containing the
	 *  random-random pairs
	 *
	 *  @param dr pointer to an object of type Pair containing the
	 *  data-random pairs
	 *
	 *  @param nData number of objects in the data catalogue
	 *
	 *  @param nData_weighted weighted number of objects in the data
	 *  catalogue
	 *
	 *  @param nRandom number of objects in the random catalogue
	 *
	 *  @param nRandom_weighted weighted number of objects in the
	 *  random catalogue
	 *
	 *  @return pointer to an object of type Data
	 */
	virtual std::shared_ptr<data::Data> correlation_LandySzalayEstimator (const std::shared_ptr<pairs::Pair> dd, const std::shared_ptr<pairs::Pair> rr, const std::shared_ptr<pairs::Pair> dr, const int nData, const double nData_weighted, const int nRandom, const double nRandom_weighted) = 0;

	/**
	 *  @brief measure the filtered two-point correlation function,
	 *  \f$w(r_c)=2\pi \int dr \xi(r) W(r,r_c) r^2\f$, where \f$W(x)
	 *  = (2x)^2(1-x)(0.5-x)r_c^{-3}\f$
	 *  
	 *  @param data pointer to an object of type Pair containing the
	 *  data-data pairs
	 *
	 *  @return pointer to an object of type Data
	 */
	virtual std::shared_ptr<data::Data> Filtered (const std::shared_ptr<data::Data> data)
	{ (void)data; ErrorCBL("", "Filtered", "TwoPointCorrelation.h"); std::shared_ptr<data::Data> dd; return dd; }
      
	/**
	 *  @brief measure the projected two-point correlation function
	 *
	 *  @param rp projected separation
	 *
	 *  @param pi line of sight separation
	 *
	 *  @param xi the 2D two-point correlation function
	 *
	 *  @param error_xi errors on the 2D two-point correlation
	 *  function
	 *
	 *  @return pointer to an object of type Data
	 */
	virtual std::shared_ptr<data::Data> Projected (const std::vector<double> rp, const std::vector<double> pi, const std::vector<std::vector<double> > xi, const std::vector<std::vector<double> > error_xi)
	{ (void)rp; (void)pi; (void)xi; (void)error_xi; ErrorCBL("", "Projected", "TwoPointCorrelation.h"); std::shared_ptr<data::Data> data; return data; }

	/**
	 *  @brief measure the deprojected two-point correlation
	 *  function
	 *  
	 *  @param rp perpendicular separation
	 *
	 *  @param ww projected two-point correlation function
	 *
	 *  @param error_ww error on the projected two-point correlation
	 *  function
	 *
	 *  @return pointer to an object of type Data
	 */
	virtual std::shared_ptr<data::Data> Deprojected (const std::vector<double> rp, const std::vector<double> ww, const std::vector<double> error_ww)
	{ (void)rp; (void)ww; (void)error_ww; ErrorCBL("", "Deprojected", "TwoPointCorrelation.h"); std::shared_ptr<data::Data> data; return data; }

	/**
	 *  @brief measure the multipoles of the two-point correlation
	 *  function
	 *  
	 *  @param rr absolute separation 
	 *
	 *  @param mu angular separation
	 *
	 *  @param xi the 2D two-point correlation function
	 *
	 *  @param error_xi errors on the 2D two-point correlation
	 *  function
	 *
	 *  @return pointer to an object of type Data
	 */
	virtual std::shared_ptr<data::Data> Multipoles (const std::vector<double> rr, const std::vector<double> mu, const std::vector<std::vector<double>> xi, const std::vector<std::vector<double>> error_xi)
	{ (void)rr; (void)mu; (void)xi; (void)error_xi; ErrorCBL("", "Multipoles", "TwoPointCorrelation.h"); std::shared_ptr<data::Data> data; return data; }

	/**
	 *  @brief measure the wedges of the two-poinr correlation
	 *  function
	 *  
	 *  @param rr absolute separation 
	 *
	 *  @param mu angular separation
	 *
	 *  @param xi the 2D two-point correlation function
	 *
	 *  @param error_xi errors on the 2D two-point correlation
	 *  function
	 *
	 *  @return pointer to an object of type Data
	 */
	virtual std::shared_ptr<data::Data> Wedges (const std::vector<double> rr, const std::vector<double> mu, const std::vector<std::vector<double>> xi, const std::vector<std::vector<double> > error_xi)
	{ (void)rr; (void)mu; (void)xi; (void)error_xi; ErrorCBL("", "Wedges", "TwoPointCorrelation.h"); std::shared_ptr<data::Data> data; return data; }
 

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
	 *  @return none, or an error message if the derived object does
	 *  not have this member
	 */
	virtual void measurePoisson (const std::string dir_output_pairs = par::defaultString, const std::vector<std::string> dir_input_pairs={}, const bool count_dd=true, const bool count_rr=true, const bool count_dr=true, const bool tcount=true, const Estimator estimator=Estimator::_LandySzalay_)
	{ (void)dir_output_pairs; (void)dir_input_pairs; (void)count_dd; (void)count_rr; (void)count_dr; (void)tcount; (void)estimator; cbl::ErrorCBL("", "measurePoisson", "TwoPointCorrelation.h"); }

	/**
	 *  @brief measure the two-point correlation function estimating
	 *  the covariance with Jackknife resampling
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
	 *  @return none, or an error message if the derived object does
	 *  not have this member
	 */
	virtual void measureJackknife (const std::string dir_output_pairs = par::defaultString, const std::vector<std::string> dir_input_pairs={}, const std::string dir_output_resample=par::defaultString, const bool count_dd=true, const bool count_rr=true, const bool count_dr=true, const bool tcount=true, const Estimator estimator=Estimator::_LandySzalay_)
	{ (void)dir_output_pairs; (void)dir_input_pairs; (void)dir_output_resample; (void)count_dd; (void)count_rr; (void)count_dr; (void)tcount; (void)estimator; cbl::ErrorCBL("", "measureJackknife", "TwoPointCorrelation.h"); }

	/**
	 *  @brief measure the two-point correlation function estimating
	 *  the covariance with Jackknife resampling, test
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
	 *  @return none, or an error message if the derived object does
	 *  not have this member
	 */
	virtual void measureJackknifeTest (const std::string dir_output_pairs = par::defaultString, const std::vector<std::string> dir_input_pairs={}, const std::string dir_output_resample=par::defaultString, const bool count_dd=true, const bool count_rr=true, const bool count_dr=true, const bool tcount=true, const Estimator estimator=Estimator::_LandySzalay_)
	{ (void)dir_output_pairs; (void)dir_input_pairs; (void)dir_output_resample; (void)count_dd; (void)count_rr; (void)count_dr; (void)tcount; (void)estimator; cbl::ErrorCBL("", "measureJackknifeTest", "TwoPointCorrelation.h"); }

	/**
	 *  @brief measure the two-point correlation function estimating
	 *  the covariance with Bootstrap resampling
	 *
	 *  @param nMocks number of mocks to be generated with Bootstrap
	 *  resampling
	 *
	 *  @param dir_output_pairs output directory used to store the
	 *  number of pairs
	 *
	 *  @param dir_input_pairs vector of input directories used to
	 *  store the number of pairs (if the pairs are read from files)
	 *
	 *  @param dir_output_resample output directory used to store
	 *  the Bootstrap resampling correlation functions, with
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
	 *  @param seed the seed for random number generation
	 *
	 *  @return none, or an error message if the derived object does
	 *  not have this member
	 */
	virtual void measureBootstrap (const int nMocks, const std::string dir_output_pairs=par::defaultString, const std::vector<std::string> dir_input_pairs={}, const std::string dir_output_resample=par::defaultString, const bool count_dd=true, const bool count_rr=true, const bool count_dr=true, const bool tcount=true, const Estimator estimator=Estimator::_LandySzalay_, const int seed=3213)
	{ (void)nMocks; (void)dir_output_pairs; (void)dir_input_pairs; (void)dir_output_resample; (void)count_dd; (void)count_rr; (void)count_dr; (void)tcount; (void)estimator; (void)seed; cbl::ErrorCBL("", "measureBootstrap", "TwoPointCorrelation.h"); }

	/**
	 *  @brief measure the Jackknife resampling of the two-point
	 *  correlation function, &xi;(r)
	 *
	 *  @param dd vector of data-data pairs, divided per regions
	 *
	 *  @param rr vector of random-random pairs, divided per regions
	 *
	 *  @return vector of pointers to objects of type Data 
	 */
	virtual std::vector<std::shared_ptr<data::Data> > XiJackknife (const std::vector<std::shared_ptr<pairs::Pair> > dd, const std::vector<std::shared_ptr<pairs::Pair> > rr)
	{ (void)dd; (void)rr; cbl::ErrorCBL("", "XiJackknife", "TwoPointCorrelation.h"); std::vector<std::shared_ptr<data::Data> > data; return data; }

	/**
	 *  @brief measure the Jackknife resampling of the two-point
	 *  correlation function, &xi;(r)
	 *
	 *  @param dd vector of data-data pairs, divided per regions
	 *
	 *  @param rr vector of random-random pairs, divided per
	 *  regions
	 *
	 *  @param dr vector of data-random pairs, divided per
	 *  regions
	 *
	 *  @return vector of pointers to objects of type Data 
	 */
	virtual std::vector<std::shared_ptr<data::Data> > XiJackknife (const std::vector<std::shared_ptr<pairs::Pair> > dd, const std::vector<std::shared_ptr<pairs::Pair> > rr, const std::vector<std::shared_ptr<pairs::Pair> > dr)
	{ (void)dd; (void)rr; (void)dr; cbl::ErrorCBL("", "XiJackknife", "TwoPointCorrelation.h"); std::vector<std::shared_ptr<data::Data> > data; return data; }

	/**
	 *  @brief measure the Jackknife resampling of the two-point
	 *  correlation function, &xi;(r)
	 *
	 *  @param dd vector of data-data pairs, divided per regions
	 *
	 *  @param rr vector of random-random pairs, divided per regions
	 *
	 *  @return vector of pointers to objects of type Data 
	 */
	virtual std::vector<std::shared_ptr<data::Data> > XiJackknifeTest (const std::vector<std::shared_ptr<pairs::Pair> > dd, const std::vector<std::shared_ptr<pairs::Pair> > rr)
	{ (void)dd; (void)rr; cbl::ErrorCBL("", "XiJackknifeTest", "TwoPointCorrelation.h"); std::vector<std::shared_ptr<data::Data> > data; return data; }

	/**
	 *  @brief measure the Jackknife resampling of the two-point
	 *  correlation function, &xi;(r)
	 *
	 *  @param dd vector of data-data pairs, divided per regions
	 *
	 *  @param rr vector of random-random pairs, divided per
	 *  regions
	 *
	 *  @param dr vector of data-random pairs, divided per regions
	 *
	 *  @return vector of pointers to objects of type Data 
	 */
	virtual std::vector<std::shared_ptr<data::Data> > XiJackknifeTest (const std::vector<std::shared_ptr<pairs::Pair> > dd, const std::vector<std::shared_ptr<pairs::Pair> > rr, const std::vector<std::shared_ptr<pairs::Pair> > dr)
	{ (void)dd; (void)rr; (void)dr; cbl::ErrorCBL("", "XiJackknifeTest", "TwoPointCorrelation.h"); std::vector<std::shared_ptr<data::Data> > data; return data; }

	/**
	 *  @brief measure the Bootstrap resampling of the two-point
	 *  correlation function, &xi;(r)
	 *
	 *  @param nMocks number of Bootstrap resampling
	 *
	 *  @param dd vector of data-data pairs, divided per regions
	 *
	 *  @param rr vector of random-random pairs, divided per
	 *  regions
	 *
	 *  @param seed the seed for random number generation
	 *
	 *  @return none, or an error message if the derived object does
	 *  not have this member
	 */
	virtual std::vector<std::shared_ptr<data::Data> > XiBootstrap (const int nMocks, const std::vector<std::shared_ptr<pairs::Pair> > dd, const std::vector<std::shared_ptr<pairs::Pair> > rr, const int seed=3213)
	{ (void)nMocks; (void)dd; (void)rr; (void)seed; cbl::ErrorCBL("", "XiBootstrap", "TwoPointCorrelation.h"); std::vector<std::shared_ptr<data::Data> > data; return data; }

	/**
	 *  @brief measure the Bootstrap resampling of the two-point
	 *  correlation function, &xi;(r)
	 *
	 *  @param nMocks number of Bootstrap resampling
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
	 *  @return none, or an error message if the derived object does
	 *  not have this member
	 */
	virtual std::vector<std::shared_ptr<data::Data> > XiBootstrap (const int nMocks, const std::vector<std::shared_ptr<pairs::Pair> > dd, const std::vector<std::shared_ptr<pairs::Pair> > rr, const std::vector<std::shared_ptr<pairs::Pair> > dr, const int seed=3213)
	{ (void)nMocks; (void)dd; (void)rr; (void)dr; (void)seed; cbl::ErrorCBL("", "XiBootstrap", "TwoPointCorrelation.h"); std::vector<std::shared_ptr<data::Data> > data; return data; }

	///@}

	
      public:
    
	/**
	 *  @name Constructors/destructors
	 */
	///@{

	/**
	 *  @brief default constructor
	 *  
	 */
	TwoPointCorrelation () = default;

	/**
	 *  @brief constructor
	 *
	 *  @param data object of class Catalogue containing the input
	 *  catalogue
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
	TwoPointCorrelation (const catalogue::Catalogue data, const catalogue::Catalogue random, const bool compute_extra_info=false, const double random_dilution_fraction=1.) 
	  : m_data(std::make_shared<catalogue::Catalogue>(catalogue::Catalogue(std::move(data)))), m_random(std::make_shared<catalogue::Catalogue>(catalogue::Catalogue(std::move(random)))), m_compute_extra_info(compute_extra_info), m_random_dilution_fraction(random_dilution_fraction) {}
    
	/**
	 *  @brief default destructor
	 *  
	 */
	virtual ~TwoPointCorrelation () = default;

	/**
	 *  @brief static factory used to construct two-point
	 *  correlation functions of any type
	 *
	 *  @param type the type of two-point correlation function; it
	 *  can be: \_monopole\_, \_angular\_
	 *
	 *  @param data object of class Catalogue containing the input
	 *  catalogue
	 *
	 *  @param random of class Catalogue containing the random data
	 *  catalogue
	 *
	 *  @param binType binning type: 0 &rarr; linear; 1 &rarr;
	 *  logarithmic
	 *
	 *  @param Min minimum separation used to count the pairs
	 *
	 *  @param Max maximum separation used to count the pairs
	 *
	 *  @param nbins number of bins
	 *
	 *  @param shift shift parameter, i.e. the radial shift is
	 *  binSize*shift
	 *
	 *  @param angularUnits angular units
	 *
	 *  @param angularWeight angular weight function
	 *
	 *  @param compute_extra_info true &rarr; compute extra
	 *  information related to the pairs, such as the mean pair
	 *  separation and redshift
	 *
	 *  @param random_dilution_fraction fraction between the number
	 *  of objects in the diluted and original random samples, used
	 *  to improve performances in random-random pair counts
	 *
	 *  @return a pointer to an object of class
	 *  TwoPointCorrelation of a given type
	 */
	static std::shared_ptr<TwoPointCorrelation> Create (const TwoPType type, const catalogue::Catalogue data, const catalogue::Catalogue random, const BinType binType, const double Min, const double Max, const int nbins, const double shift, const CoordinateUnits angularUnits=CoordinateUnits::_radians_, std::function<double(double)> angularWeight=nullptr, const bool compute_extra_info=false, const double random_dilution_fraction=1.);

	/**
	 *  @brief static factory used to construct two-point
	 *  correlation functions of any type
	 *
	 *  @param type the type of two-point correlation function; it
	 *  can be: \_monopole\_, \_angular\_
	 *
	 *  @param data object of class Catalogue containing the input
	 *  catalogue
	 *
	 *  @param random of class Catalogue containing the random data
	 *  catalogue
	 *
	 *  @param binType binning type: 0 &rarr; linear; 1 &rarr;
	 *  logarithmic
	 *
	 *  @param Min minimum separation used to count the pairs
	 *
	 *  @param Max maximum separation used to count the pairs
	 *
	 *  @param binSize the bin size
	 *
	 *  @param shift shift parameter, i.e. the radial shift is
	 *  binSize*shift
	 *
	 *  @param angularUnits angular units
	 *
	 *  @param angularWeight angular weight function
	 *
	 *  @param compute_extra_info true &rarr; compute extra
	 *  information related to the pairs, such as the mean pair
	 *  separation and redshift
	 *
	 *  @param random_dilution_fraction fraction between the number
	 *  of objects in the diluted and original random samples, used
	 *  to improve performances in random-random pair counts
	 *
	 *  @return a pointer to an object of class
	 *  TwoPointCorrelation of a given type
	 */
	static std::shared_ptr<TwoPointCorrelation> Create (const TwoPType type, const catalogue::Catalogue data, const catalogue::Catalogue random, const BinType binType, const double Min, const double Max, const double binSize, const double shift, const CoordinateUnits angularUnits=CoordinateUnits::_radians_, std::function<double(double)> angularWeight=nullptr, const bool compute_extra_info=false, const double random_dilution_fraction=1.);

	/**
	 *  @brief static factory used to construct two-point
	 *  correlation functions of any type
	 *
	 *  @param type the type of two-point correlation function; it
	 *  can be: \_2D_Cartesian\_, \_2D_polar\_
	 *
	 *  @param data object of class Catalogue containing the input
	 *  catalogue
	 *
	 *  @param random of class Catalogue containing the random data
	 *  catalogue
	 *
	 *  @param binType_D1 binning type in the first dimension: 0
	 *  &rarr; linear; 1 &rarr; logarithmic
	 *
	 *  @param Min_D1 minimum separation in the first dimensionused
	 *  to count the pairs
	 *
	 *  @param Max_D1 maximum separation in the first dimension used
	 *  to count the pairs
	 *
	 *  @param nbins_D1 number of bins in the first dimension
	 *
	 *  @param shift_D1 shift parameter in the first dimension,
	 *  i.e. the radial shift is binSize*shift
	 *
	 *  @param binType_D2 binning type in the second dimension: 0
	 *  &rarr; linear; 1 &rarr; logarithmic
	 *
	 *  @param Min_D2 minimum separation in the second dimensionused
	 *  to count the pairs
	 *
	 *  @param Max_D2 maximum separation in the second dimension
	 *  used to count the pairs
	 *
	 *  @param nbins_D2 number of bins in the second dimension
	 *
	 *  @param shift_D2 shift parameter in the second dimension,
	 *  i.e. the radial shift is binSize*shift
	 * 
	 *  @param angularUnits angular units
	 *  @param angularWeight angular weight function
	 *
	 *  @param compute_extra_info true &rarr; compute extra
	 *  information related to the pairs, such as the mean pair
	 *  separation and redshift
	 *
	 *  @param random_dilution_fraction fraction between the number
	 *  of objects in the diluted and original random samples, used
	 *  to improve performances in random-random pair counts
	 *
	 *  @return a pointer to an object of class
	 *  TwoPointCorrelation of a given type
	 */
	static std::shared_ptr<TwoPointCorrelation> Create (const TwoPType type, const catalogue::Catalogue data, const catalogue::Catalogue random, const BinType binType_D1, const double Min_D1, const double Max_D1, const int nbins_D1, const double shift_D1, const BinType binType_D2, const double Min_D2, const double Max_D2, const int nbins_D2, const double shift_D2, const CoordinateUnits angularUnits=CoordinateUnits::_radians_, std::function<double(double)> angularWeight=nullptr, const bool compute_extra_info=false, const double random_dilution_fraction=1.);

	/**
	 *  @brief static factory used to construct two-point
	 *  correlation functions of any type
	 *
	 *  @param type the type of two-point correlation function; it
	 *  can be: \_2D_Cartesian\_, \_2D_polar\_
	 *
	 *  @param data object of class Catalogue containing the input
	 *  catalogue
	 *
	 *  @param random of class Catalogue containing the random data
	 *  catalogue
	 *
	 *  @param binType_D1 binning type in the first dimension: 0
	 *  &rarr; linear; 1 &rarr; logarithmic
	 *
	 *  @param Min_D1 minimum separation in the first dimension used
	 *  to count the pairs
	 *
	 *  @param Max_D1 maximum separation in the first dimension used
	 *  to count the pairs
	 *
	 *  @param binSize_D1 the bin size in the first dimension
	 *
	 *  @param shift_D1 shift parameter in the first dimension,
	 *  i.e. the radial shift is binSize*shift
	 *
	 *  @param binType_D2 binning type in the second dimension: 0
	 *  &rarr; linear; 1 &rarr; logarithmic
	 *
	 *  @param Min_D2 minimum separation in the second dimension
	 *  used to count the pairs
	 *
	 *  @param Max_D2 maximum separation in the second dimension
	 *  used to count the pairs
	 *
	 *  @param binSize_D2 the bin size in the second dimension
	 *
	 *  @param shift_D2 shift parameter in the second dimension,
	 *  i.e. the radial shift is binSize*shift
	 *
	 *  @param angularUnits angular units
	 *  @param angularWeight angular weight function
	 *
	 *  @param compute_extra_info true &rarr; compute extra
	 *  information related to the pairs, such as the mean pair
	 *  separation and redshift
	 *
	 *  @param random_dilution_fraction fraction between the number
	 *  of objects in the diluted and original random samples, used
	 *  to improve performances in random-random pair counts
	 *
	 *  @return a pointer to an object of class
	 *  TwoPointCorrelation of a given type
	 */
	static std::shared_ptr<TwoPointCorrelation> Create (const TwoPType type, const catalogue::Catalogue data, const catalogue::Catalogue random, const BinType binType_D1, const double Min_D1, const double Max_D1, const double binSize_D1, const double shift_D1, const BinType binType_D2, const double Min_D2, const double Max_D2, const double binSize_D2, const double shift_D2, const CoordinateUnits angularUnits=CoordinateUnits::_radians_, std::function<double(double)> angularWeight=nullptr, const bool compute_extra_info=false, const double random_dilution_fraction=1.);
      
	/**
	 *  @brief static factory used to construct two-point
	 *  correlation functions of any type
	 *
	 *  @param type the type of two-point correlation function; it
	 *  can be: _multipoles_integrated_, _filtered_
	 *
	 *  @param data object of class Catalogue containing the input
	 *  catalogue
	 *
	 *  @param random of class Catalogue containing the random data
	 *  catalogue
	 *
	 *  @param binType_D1 binning type in the first dimension: 0
	 *  &rarr; linear; 1 &rarr; logarithmic
	 *
	 *  @param Min_D1 minimum separation in the first dimensionused
	 *  to count the pairs
	 *
	 *  @param Max_D1 maximum separation in the first dimension used
	 *  to count the pairs
	 *
	 *  @param nbins_D1 number of bins in the first dimension
	 *
	 *  @param shift_D1 shift parameter in the first dimension,
	 *  i.e. the radial shift is binSize*shift
	 *
	 *  @param Min_D2 minimum separation in the second dimensionused
	 *  to count the pairs
	 *
	 *  @param Max_D2 maximum separation in the second dimension
	 *  used to count the pairs
	 *
	 *  @param nbins_D2 number of bins in the second dimension
	 *
	 *  @param shift_D2 shift parameter in the second dimension,
	 *  i.e. the radial shift is binSize*shift
	 *
	 *  @param angularUnits angular units
	 *  @param angularWeight angular weight function
	 *
	 *  @param compute_extra_info true &rarr; compute extra
	 *  information related to the pairs, such as the mean pair
	 *  separation and redshift
	 *
	 *  @param random_dilution_fraction fraction between the number
	 *  of objects in the diluted and original random samples, used
	 *  to improve performances in random-random pair counts
	 *
	 *  @return a pointer to an object of class
	 *  TwoPointCorrelation of a given type
	 */
	static std::shared_ptr<TwoPointCorrelation> Create (const TwoPType type, const catalogue::Catalogue data, const catalogue::Catalogue random, const BinType binType_D1, const double Min_D1, const double Max_D1, const int nbins_D1, const double shift_D1, const double Min_D2, const double Max_D2, const int nbins_D2, const double shift_D2, const CoordinateUnits angularUnits=CoordinateUnits::_radians_, std::function<double(double)> angularWeight=nullptr, const bool compute_extra_info=false, const double random_dilution_fraction=1.);

	/**
	 *  @brief static factory used to construct two-point
	 *  correlation functions of any type
	 *
	 *  @param type the type of two-point correlation function; it
	 *  can be: _multipoles_integrated_, _filtered_
	 *
	 *  @param data object of class Catalogue containing the input
	 *  catalogue
	 *
	 *  @param random of class Catalogue containing the random data
	 *  catalogue
	 *
	 *  @param binType_D1 binning type in the first dimension: 0
	 *  &rarr; linear; 1 &rarr; logarithmic
	 *
	 *  @param Min_D1 minimum separation in the first dimension used
	 *  to count the pairs
	 *
	 *  @param Max_D1 maximum separation in the first dimension used
	 *  to count the pairs
	 *
	 *  @param binSize_D1 the bin size in the first dimension
	 *
	 *  @param shift_D1 shift parameter in the first dimension,
	 *  i.e. the radial shift is binSize*shift
	 *
	 *  @param Min_D2 minimum separation in the second dimension used
	 *  to count the pairs
	 *
	 *  @param Max_D2 maximum separation in the second dimension used
	 *  to count the pairs
	 *
	 *  @param binSize_D2 the bin size in the second dimension
	 *
	 *  @param shift_D2 shift parameter in the second dimension,
	 *  i.e. the radial shift is binSize*shift
	 *
	 *  @param angularUnits angular units
	 *  @param angularWeight angular weight function
	 *
	 *  @param compute_extra_info true &rarr; compute extra
	 *  information related to the pairs, such as the mean pair
	 *  separation and redshift
	 *
	 *  @param random_dilution_fraction fraction between the number
	 *  of objects in the diluted and original random samples, used
	 *  to improve performances in random-random pair counts
	 *
	 *  @return a pointer to an object of class
	 *  TwoPointCorrelation of a given type
	 */
	static std::shared_ptr<TwoPointCorrelation> Create (const TwoPType type, const catalogue::Catalogue data, const catalogue::Catalogue random, const BinType binType_D1, const double Min_D1, const double Max_D1, const double binSize_D1, const double shift_D1, const double Min_D2, const double Max_D2, const double binSize_D2, const double shift_D2, const CoordinateUnits angularUnits=CoordinateUnits::_radians_, std::function<double(double)> angularWeight=nullptr, const bool compute_extra_info=false, const double random_dilution_fraction=1.);

	/**
	 *  @brief static factory used to construct two-point
	 *  correlation functions of any type
	 *
	 *  @param type the type of two-point correlation function; it
	 *  can be: \_projected\_, \_deprojected\_
	 *
	 *  @param data object of class Catalogue containing the input
	 *  catalogue
	 *
	 *  @param random of class Catalogue containing the random data
	 *  catalogue
	 *
	 *  @param binType_D1 binning type in the first dimension: 0
	 *  &rarr; linear; 1 &rarr; logarithmic
	 *
	 *  @param Min_D1 minimum separation in the first dimensionused
	 *  to count the pairs
	 *
	 *  @param Max_D1 maximum separation in the first dimension used
	 *  to count the pairs
	 *
	 *  @param nbins_D1 number of bins in the first dimension
	 *
	 *  @param shift_D1 shift parameter in the first dimension,
	 *  i.e. the radial shift is binSize*shift
	 *
	 *  @param Min_D2 minimum separation in the second dimensionused
	 *  to count the pairs
	 *
	 *  @param Max_D2 maximum separation in the second dimension
	 *  used to count the pairs
	 *
	 *  @param nbins_D2 number of bins in the second dimension
	 *
	 *  @param shift_D2 shift parameter in the second dimension,
	 *  i.e. the radial shift is binSize*shift
	 *
	 *  @param piMax_integral upper limits of the integral
	 *
	 *  @param angularUnits angular units
	 *  @param angularWeight angular weight function
	 *
	 *  @param compute_extra_info true &rarr; compute extra
	 *  information related to the pairs, such as the mean pair
	 *  separation and redshift
	 *
	 *  @param random_dilution_fraction fraction between the number
	 *  of objects in the diluted and original random samples, used
	 *  to improve performances in random-random pair counts
	 *
	 *  @return a pointer to an object of class
	 *  TwoPointCorrelation of a given type
	 */
	static std::shared_ptr<TwoPointCorrelation> Create (const TwoPType type, const catalogue::Catalogue data, const catalogue::Catalogue random, const BinType binType_D1, const double Min_D1, const double Max_D1, const int nbins_D1, const double shift_D1, const double Min_D2, const double Max_D2, const int nbins_D2, const double shift_D2, const double piMax_integral, const CoordinateUnits angularUnits=CoordinateUnits::_radians_, std::function<double(double)> angularWeight=nullptr, const bool compute_extra_info=false, const double random_dilution_fraction=1.);

	/**
	 *  @brief static factory used to construct two-point
	 *  correlation functions of any type
	 *
	 *  @param type the type of two-point correlation function; it
	 *  can be: \_projected\_, \_deprojected\_
	 *
	 *  @param data object of class Catalogue containing the input
	 *  catalogue
	 *
	 *  @param random of class Catalogue containing the random data
	 *  catalogue
	 *
	 *  @param binType_D1 binning type in the first dimension: 0
	 *  &rarr; linear; 1 &rarr; logarithmic
	 *
	 *  @param Min_D1 minimum separation in the first dimension used
	 *  to count the pairs
	 *
	 *  @param Max_D1 maximum separation in the first dimension used
	 *  to count the pairs
	 *
	 *  @param binSize_D1 the bin size in the first dimension
	 *
	 *  @param shift_D1 shift parameter in the first dimension,
	 *  i.e. the radial shift is binSize*shift
	 *
	 *  @param Min_D2 minimum separation in the second dimension used
	 *  to count the pairs
	 *
	 *  @param Max_D2 maximum separation in the second dimension used
	 *  to count the pairs
	 *
	 *  @param binSize_D2 the bin size in the second dimension
	 *
	 *  @param shift_D2 shift parameter in the second dimension,
	 *  i.e. the radial shift is binSize*shift
	 *
	 *  @param piMax_integral upper limits of the integral
	 *
	 *  @param angularUnits angular units
	 *  @param angularWeight angular weight function
	 *
	 *  @param compute_extra_info true &rarr; compute extra
	 *  information related to the pairs, such as the mean pair
	 *  separation and redshift
	 *
	 *  @param random_dilution_fraction fraction between the number
	 *  of objects in the diluted and original random samples, used
	 *  to improve performances in random-random pair counts
	 *
	 *  @return a pointer to an object of class
	 *  TwoPointCorrelation of a given type
	 */
	static std::shared_ptr<TwoPointCorrelation> Create (const TwoPType type, const catalogue::Catalogue data, const catalogue::Catalogue random, const BinType binType_D1, const double Min_D1, const double Max_D1, const double binSize_D1, const double shift_D1, const double Min_D2, const double Max_D2, const double binSize_D2, const double shift_D2, const double piMax_integral, const CoordinateUnits angularUnits=CoordinateUnits::_radians_, std::function<double(double)> angularWeight=nullptr, const bool compute_extra_info=false, const double random_dilution_fraction=1.);

	/**
	 *  @brief static factory used to construct two-point
	 *  correlation functions of any type
	 *
	 *  @param type the type of two-point correlation function; it
	 *  can be: _wedges_
	 *
	 *  @param data object of class Catalogue containing the input
	 *  catalogue
	 *
	 *  @param random of class Catalogue containing the random data
	 *  catalogue
	 *
	 *  @param binType_D1 binning type in the first dimension: 0
	 *  &rarr; linear; 1 &rarr; logarithmic
	 *
	 *  @param Min_D1 minimum separation in the first dimensionused
	 *  to count the pairs
	 *
	 *  @param Max_D1 maximum separation in the first dimension used
	 *  to count the pairs
	 *
	 *  @param nbins_D1 number of bins in the first dimension
	 *
	 *  @param shift_D1 shift parameter in the first dimension,
	 *  i.e. the radial shift is binSize*shift
	 *
	 *  @param nWedges number of wedges to be measured; the
	 *  default is two wedges, \f$\xi_\perp\f$ and
	 *  \f$\xi_\parallel\f$
	 *
	 *  @param nbins_D2 number of bins in the second dimension
	 *
	 *  @param shift_D2 shift parameter in the second dimension,
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
	 *  @param random_dilution_fraction fraction between the number
	 *  of objects in the diluted and original random samples, used
	 *  to improve performances in random-random pair counts
	 *
	 *  @return a pointer to an object of class
	 *  TwoPointCorrelation of a given type
	 */
	static std::shared_ptr<TwoPointCorrelation> Create (const TwoPType type, const catalogue::Catalogue data, const catalogue::Catalogue random, const BinType binType_D1, const double Min_D1, const double Max_D1, const int nbins_D1, const double shift_D1, const int nWedges, const int nbins_D2, const double shift_D2, const std::vector<std::vector<double>> mu_integral_limits={{0., 0.5}, {0.5, 1}}, const CoordinateUnits angularUnits=CoordinateUnits::_radians_, std::function<double(double)> angularWeight=nullptr, const bool compute_extra_info=false, const double random_dilution_fraction=1.);

	/**
	 *  @brief static factory used to construct two-point
	 *  correlation functions of any type
	 *
	 *  @param type the type of two-point correlation function; it
	 *  can be: _wedges_,
	 *
	 *  @param data object of class Catalogue containing the input
	 *  catalogue
	 *
	 *  @param random of class Catalogue containing the random data
	 *  catalogue
	 *
	 *  @param binType_D1 binning type in the first dimension: 0
	 *  &rarr; linear; 1 &rarr; logarithmic
	 *
	 *  @param Min_D1 minimum separation in the first dimension used
	 *  to count the pairs
	 *
	 *  @param Max_D1 maximum separation in the first dimension used
	 *  to count the pairs
	 *
	 *  @param binSize_D1 the bin size in the first dimension
	 *
	 *  @param shift_D1 shift parameter in the first dimension,
	 *  i.e. the radial shift is binSize*shift
	 *
	 *  @param binSize_D2 the bin size in the second dimension
	 *
	 *  @param nWedges number of wedges to be measured; the
	 *  default is two wedges, \f$\xi_\perp\f$ and
	 *  \f$\xi_\parallel\f$
	 *
	 *  @param shift_D2 shift parameter in the second dimension,
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
	 *  @return a pointer to an object of class
	 *  TwoPointCorrelation of a given type
	 */
	static std::shared_ptr<TwoPointCorrelation> Create (const TwoPType type, const catalogue::Catalogue data, const catalogue::Catalogue random, const BinType binType_D1, const double Min_D1, const double Max_D1, const double binSize_D1, const double shift_D1, const int nWedges, const double binSize_D2, const double shift_D2, const std::vector<std::vector<double>> mu_integral_limits={{0., 0.5}, {0.5, 1}}, const CoordinateUnits angularUnits=CoordinateUnits::_radians_, std::function<double(double)> angularWeight=nullptr, const bool compute_extra_info=false, const double random_dilution_fraction=1.);

	///@}

    
	/**
	 *  @name Member functions to get the private/protected members
	 */
	///@{
      
	/**
	 *  @brief get the protected member m_twoPType
	 *  @return the two-point correlation function type
	 */
	TwoPType twoPType () const { return m_twoPType; }
      
	/**
	 *  @brief get the protected member m_data
	 *  @return the input data catalogue
	 */
	std::shared_ptr<catalogue::Catalogue> data () const { return m_data; }
   
	/**
	 *  @brief get the protected member m_random
	 *  @return the input random catalogue
	 */
	std::shared_ptr<catalogue::Catalogue> random () const { return m_random; }

	/**
	 *  @brief get the protected member m_dd
	 *  @return the number of data-data pairs
	 */
	std::shared_ptr<pairs::Pair> dd () const { return m_dd; }

	/**
	 *  @brief get the protected member m_rr
	 *  @return the number of random-random pairs
	 */
	std::shared_ptr<pairs::Pair> rr () const { return m_rr; }

	/**
	 *  @brief get the protected member m_dr
	 *  @return the number of data-random pairs
	 */
	std::shared_ptr<pairs::Pair> dr () const { return m_dr; }

	/**
	 *  @brief get the protected member m_compute_extra_info
	 *  @return true &rarr; compute extra information related to the
	 *  pairs, such as the mean pair separation and redshift
	 */
	bool compute_extra_info () const { return m_compute_extra_info; }

	/**
	 *  @brief get the protected member m_random_dilution_fraction
	 *  @return the fraction between the number of objects in the
	 *  diluted and original random samples, used to improve
	 *  performances in random-random pair counts
	 */
	bool random_dilution_fraction () const { return m_random_dilution_fraction; }
      
	/**
	 *  @brief get the x coordinates
	 *  @return the x coordinates
	 */
	virtual std::vector<double> xx () const = 0;

	/**
	 *  @brief get the y coordinates
	 *  @return the y coordinates
	 */
	virtual std::vector<double> yy () const { cbl::ErrorCBL("", "yy", "TwoPointCorrelation.h"); std::vector<double> vv; return vv; }

	/**
	 *  @brief get the the binned correlation function 
	 *  @return the binned correlation function 
	 */
	virtual std::vector<double> xi1D () const { cbl::ErrorCBL("", "xi1D", "TwoPointCorrelation.h"); std::vector<double> vv; return vv; }

	/**
	 *  @brief get the error on the binned correlation function
	 *  function
	 *  @return the error on the binned correlation function
	 *  function
	 */
	virtual std::vector<double> error1D () const { cbl::ErrorCBL("", "error1D", "TwoPointCorrelation.h"); std::vector<double> vv; return vv; }
      
	/**
	 *  @brief get the the binned correlation function 
	 *  @return the binned correlation function 
	 */
	virtual std::vector<std::vector<double> > xi2D () const { cbl::ErrorCBL("", "xi2D", "TwoPointCorrelation.h"); std::vector<std::vector<double> > vv; return vv; }

	/**
	 *  @brief get the error on the binned correlation function
	 *  @return the error on the binned correlation function
	 */
	virtual std::vector<std::vector<double> > error2D () const { cbl::ErrorCBL("", "error2D", "TwoPointCorrelation.h"); std::vector<std::vector<double> > vv; return vv; }

	/**
	 *  @brief get the monopole of the polar xi
	 *  @return the xiMonopole
	 */
	virtual std::vector<double> xiMonopole () const 
	{ cbl::ErrorCBL("", "xiMonopole", "TwoPointCorrelation.h"); std::vector<double> vv; return vv; }

	/**
	 *  @brief get the error on the monopole of the polar xi
	 *  @return the xiMonopole
	 */
	virtual std::vector<double> errorMonopole () const 
	{ cbl::ErrorCBL("", "errorMonopole", "TwoPointCorrelation.h"); std::vector<double> vv; return vv; }

	/**
	 *  @brief get the quadrupole of the polar xi
	 *  @return the xiQuadrupole
	 */
	virtual std::vector<double> xiQuadrupole () const 
	{ cbl::ErrorCBL("", "xiQuadrupole", "TwoPointCorrelation.h"); std::vector<double> vv; return vv; }

	/**
	 *  @brief get the error on the quadrupole of the polar xi
	 *  @return the xiQuadrupole
	 */
	virtual std::vector<double> errorQuadrupole () const 
	{ cbl::ErrorCBL("", "errorQuadrupole", "TwoPointCorrelation.h"); std::vector<double> vv; return vv; }
	
	/**
	 *  @brief get the hexadecapole of the polar xi
	 *  @return the xiHexadecapole
	 */
	virtual std::vector<double> xiHexadecapole () const 
	{ cbl::ErrorCBL("", "xiHexadecapole", "TwoPointCorrelation.h"); std::vector<double> vv; return vv; }

	/**
	 *  @brief get the error on the hexadecapole of the polar xi
	 *  @return the error on xiOctupole
	 */
	virtual std::vector<double> errorHexadecapole () const 
	{ cbl::ErrorCBL("", "errorHexadecapole", "TwoPointCorrelation.h"); std::vector<double> vv; return vv; }

	/**
	 *  @brief get the perpendicular wedge of the polar xi
	 *  @return the perpendicular wedge of the polar xi
	 */
	virtual std::vector<double> xiPerpendicular () const 
	{ cbl::ErrorCBL("", "xiPerpendicular", "TwoPointCorrelation.h"); std::vector<double> vv; return vv; }

	/**
	 *  @brief get the perpendicular wedge of the polar xi
	 *  @return the error on the perpendicular wedge of the polar xi
	 */
	virtual std::vector<double> errorPerpendicular () const 
	{ cbl::ErrorCBL("", "errorPerpendicular", "TwoPointCorrelation.h"); std::vector<double> vv; return vv; }

	/**
	 *  @brief get the parallel wedge of the polar xi
	 *  @return the xiPar coordinates
	 */
	virtual std::vector<double> xiParallel () const 
	{ cbl::ErrorCBL("", "xiParallel", "TwoPointCorrelation.h"); std::vector<double> vv; return vv; }

	/**
	 *  @brief get the error on the parallel wedge of the polar xi
	 *  @return the error on the parallel wedge of the polar xi
	 */
	virtual std::vector<double> errorParallel () const 
	{ cbl::ErrorCBL("", "errorParallel", "TwoPointCorrelation.h"); std::vector<double> vv; return vv; }

	///@}
      

	/**
	 *  @name Member functions to set protected members
	 */
	///@{
    
	/**
	 *  @brief add a data catalogue
	 *  @param data object of class Catalogue 
	 */
	void set_data (const catalogue::Catalogue data) { m_data = std::make_shared<catalogue::Catalogue>(catalogue::Catalogue(std::move(data))); }

	/**
	 *  @brief add a random catalogue
	 *  @param random object of class Catalogue 
	 *  
	 */
	void set_random (const catalogue::Catalogue random) { m_random = std::make_shared<catalogue::Catalogue>(catalogue::Catalogue(std::move(random))); }

	///@}

    
	/**
	 *  @name Member functions to measure the two-point correlation function
	 */
	///@{

	/**
	 *  @brief measure the two-point correlation function
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
	 *  resampling correlation functions; if an empty string
	 *  (i.e. "" or "NULL") is provided, no output will be stored
	 *
	 *  @param nMocks number of resampling used for Bootstrap
	 *
	 *  @param count_dd true &rarr; count the number of data-data
	 *  pairs; false &rarr; read the number of data-random pairs
	 *  from file
	 *
	 *  @param count_rr true &rarr; count the number of
	 *  random-random pairs; false &rarr; read the number of
	 *  random-random pairs
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
	 *  @param seed the seed for random number generation
	 *
	 *  @return none, or an error message if the derived object does
	 *  not have this member
	 */
	virtual void measure (const ErrorType errorType=ErrorType::_Poisson_, const std::string dir_output_pairs=par::defaultString, const std::vector<std::string> dir_input_pairs={}, const std::string dir_output_resample=par::defaultString, const int nMocks=0, const bool count_dd=true, const bool count_rr=true, const bool count_dr=true, const bool tcount=true, const Estimator estimator=Estimator::_LandySzalay_, const int seed=3213) = 0;

	///@}

  
	/**
	 *  @name Input/Output member functions (customized in all the derived classes)
	 */
	///@{

	/**
	 *  @brief read the measured two-point correlation
	 *  @param dir input directory
	 *  @param file input file
	 *  @return none, or an error message if the derived object does
	 *  not have this member
	 */
	virtual void read (const std::string dir, const std::string file) = 0;

	/**
	 *  @brief write the measured two-point correlation
	 *  @param dir output directory
	 *  @param file output file
	 *  @param rank cpu index (for MPI usage)
	 *  @return none, or an error message if the derived object does
	 *  not have this member
	 */
	virtual void write (const std::string dir=par::defaultString, const std::string file=par::defaultString, const int rank=0) const = 0;

	/**
	 *  @brief write the measured two-point correlation
	 *  @param dir output directory
	 *  @param file output file
	 *  @param full false &rarr; simply store the data; true &rarr;
	 *  duplicate the data in the other three quadrands (usefull
	 *  e.g. when storing the 2D correlation function)
	 *  @param rank cpu index (for MPI usage)
	 *  @return none, or an error message if the derived object does
	 *  not have this member
	 */
	virtual void write (const std::string dir, const std::string file, const bool full, const int rank=0) const
	{ (void)dir; (void)file; (void)full; (void)rank; cbl::ErrorCBL("", "write", "TwoPointCorrelation.h"); }
      
	///@}
      

	/**
	 *  @name Member functions to estimate the errors and covariance matrices
	 */
	///@{ 

	/**
	 *  @brief read the measured covariance matrix
	 *  @param dir input directory
	 *  @param file input file
	 *  @return none, or an error message if the derived object
	 *  does not have this member
	 */
	virtual void read_covariance (const std::string dir, const std::string file) = 0;

	/**
	 *  @brief write the measured two-point correlation
	 *  @param dir output directory
	 *  @param file output file
	 *  @return none, or an error message if the derived object
	 *  does not have this member
	 */
	virtual void write_covariance (const std::string dir, const std::string file) const = 0;

	/**
	 *  @brief compute the covariance matrix
	 *  @param xi vector containing the measure correlation
	 *  functions used to compute the covariance matrix
	 *  @param JK true &rarr; compute the Jackknife covariance
	 *  matrix; false compute the standard covariance matrix
	 *  @return none, or an error message if the derived object
	 *  does not have this member
	 */
	virtual void compute_covariance (const std::vector<std::shared_ptr<data::Data>> xi, const bool JK) = 0;

	/**
	 *  @brief compute the covariance matrix
	 *  @param file vector containing the input files with the
	 *  measured correlation functions used to compute the
	 *  covariance matrix
	 *  @param JK true &rarr; compute the Jackknife covariance
	 *  matrix; false compute the standard covariance matrix
	 *  @return none, or an error message if the derived object
	 *  does not have this member
	 */
	virtual void compute_covariance (const std::vector<std::string> file, const bool JK) = 0;
    
	/**
	 *  @brief the Poisson errors 
	 *
	 *  This function computes the Poisson errors associated to the
	 *  natural or Landy&Szalay estimators of the two-point
	 *  correlation function. 
	 *  
	 *  The formula implemented is the following:
	 *   
	 *  \f[ \delta\xi = \sqrt{\left(N_1\frac{\sqrt{DD}}{RR}\right)^2
	 *  + \left(N_2\frac{\sqrt{DR}}{RR}\right)^2 +
	 *  \left(\frac{N_1DD-N_2DR}{RR^{1.5}}\right)^2} \f]
	 *
	 *  where 
	 *
	 *  \f[ N_1 = \frac{f_Rn_R(f_Rn_R-1)}{n_D(n_D-1)} \f]
	 *
	 *  and
	 *       
	 *  \f[ N_2 = \frac{f_Rn_R(f_Rn_R-1)}{n_Rn_D} \f]
	 *
	 *  DD, RR and DR are the un-normalised (weighted) numbers of
	 *  data-data, random-random and data-random pairs,
	 *  respectively. \f$n_D\f$ and \f$n_R\f$ are the total
	 *  (weighted) number of data and random objects,
	 *  respectively. \f$f_R=n_{R,dil}/n_R\f$ is the fraction
	 *  between the total (weighted) number of random objects and
	 *  the diluted one (used to improve performances when counting
	 *  the random-random pairs).
	 *
	 *  @param estimator the estimator used to measure the two-point
	 *  correlation function
	 *
	 *  @param dd number of data-data pairs
	 *  @param rr number of random-random pairs
	 *  @param dr number of data-random pairs
	 *  @param nData number of data points
	 *  @param nRandom number of random points
	 *
	 *  @return the Poisson error
	 *  
	 *  @warning This function currently works with only the natural
	 *  and Landy&Szalay estimators
	 */
	double PoissonError (const Estimator estimator, const double dd, const double rr, const double dr, const int nData, const int nRandom) const;
    
	///@}

      };
    }
  }
}

#endif
