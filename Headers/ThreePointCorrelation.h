/********************************************************************
 *  Copyright (C) 2010 by Federico Marulli, Michele Moresco         *
 *  and Alfonso Veropalumbo                                         *
 *                                                                  *
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
 *  @file Headers/ThreePointCorrelation.h
 *
 *  @brief The class ThreePointCorrelation
 *
 *  This file defines the interface of the class ThreePointCorrelation,
 *  used to measure the three-point correlation function
 *
 *  @authors Federico Marulli, Michele Moresco, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unibo.it, michele.moresco@unibo.it,
 *  alfonso.veropalumbo@unibo.it
 */

#ifndef __THREEPOINT__
#define __THREEPOINT__ 


#include "Measure.h"
#include "ChainMesh_Catalogue.h"
#include "Triplet1D.h"
#include "Triplet2D.h"


// ===================================================================================================


namespace cbl {
  
  namespace measure {

    /**
     *  @brief The namespace of the <B> three-point correlation function
     *  </B>
     *  
     * The \e measure::threept namespace contains all the functions and
     * classes to measure the three-point correlation function
     */
    namespace threept {

      /**
       * @enum ThreePType
       * @brief the three-point correlation function type
       */
      enum class ThreePType { 

	/// the connected three-point correlation function in angular coordinates
	_angular_connected_,

	/// the reduced three-point correlation function in angular coordinates
	_angular_reduced_,
      
	/// the connected three-point correlation function in comoving coordinates
	_comoving_connected_,

	/// the reduced three-point correlation function in comoving coordinates
	_comoving_reduced_
      
      };

      /**
       * @brief return a vector containing the
       * ThreePType names
       * @return a vector containing the
       * ThreePType names
       */
      inline std::vector<std::string> ThreePTypeNames () {return {"angular_connected", "angular_reduced", "comoving_connected", "comoving_reduced"}; }

      /**
       * @brief cast an enum of type ThreePType
       * from its index
       * @param threePTypeIndex the threePType index
       * @return object of class ThreePType
       */
      inline ThreePType ThreePTypeCast (const int threePTypeIndex) {return castFromValue<ThreePType>(threePTypeIndex);}

      /**
       * @brief cast an enum of type ThreePType
       * from its name
       * @param threePTypeName the threePType name
       * @return object of class ThreePType
       */
      inline ThreePType ThreePTypeCast (const std::string threePTypeName) {return castFromName<ThreePType>(threePTypeName, ThreePTypeNames());}

      /**
       * @brief cast an enum of type ThreePType
       * from indeces
       * @param threePTypeIndeces the threePType indeces
       * @return object of class ThreePType
       */
      inline std::vector<ThreePType> ThreePTypeCast (const std::vector<int> threePTypeIndeces) {return castFromValues<ThreePType>(threePTypeIndeces);} 

      /**
       * @brief cast an enum of type ThreePType
       * from thier names
       * @param threePTypeNames the threePType names
       * @return vector of ThreePType enums
       */
      inline std::vector<ThreePType> ThreePTypeCast (const std::vector<std::string> threePTypeNames) {return castFromNames<ThreePType>(threePTypeNames, ThreePTypeNames());}

      /**
       *  @class ThreePointCorrelation ThreePointCorrelation.h
       * "Headers/ThreePointCorrelation.h"
       *
       *  @brief The class ThreePointCorrelation
       *
       *  This is the base class used to measure the three-point
       *  correlation function
       */
      class ThreePointCorrelation : public Measure {
	
      protected :

	/**
	 *  @name Three-point correlation function data
	 */
	///@{
      
	/// three-point correlation function type
	ThreePType m_threePType;

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
	 *  @name Matrix to store the three-point correlation function for each Jackknife resampling
	 */
	///@{

	/// jackknife resampling matrix
	std::vector< std::vector<double> > m_resampling_threept;

	/**
	 *  @name Object triplets
	 */
	///@{

	/// number of data-data-data triplets
	std::shared_ptr<triplets::Triplet> m_ddd;

	/// number of random-random-random triplets
	std::shared_ptr<triplets::Triplet> m_rrr;

	/// number of data-data-random triplets
	std::shared_ptr<triplets::Triplet> m_ddr;
    
	/// number of data-random-random triplets
	std::shared_ptr<triplets::Triplet> m_drr;

 
	/// number of data-data-data triplets
	std::vector<std::shared_ptr<triplets::Triplet>> m_ddd_regions;

	/// number of random-random-random triplets
	std::vector<std::shared_ptr<triplets::Triplet>> m_rrr_regions;

	/// number of data-data-random triplets
	std::vector<std::shared_ptr<triplets::Triplet>> m_ddr_regions;
    
	/// number of data-random-random triplets
	std::vector<std::shared_ptr<triplets::Triplet>> m_drr_regions;

	///@}

 
	/**
	 *  @name Member functions to count the number of triplets and for I/O 
	 */
	///@{
      
	/**
	 * @brief method to count the number of triplets
	 *
	 * @param cat1 object of class \e Catalogue: the input catalogue
	 * to analyse
	 *
	 * @param ChainMesh_rMAX1 object of class \e ChainMesh_Catalogue:
	 * the chain mesh used to count the pairs relative to the first
	 * side of the triangle (1-2)
	 *
	 * @param ChainMesh_rMAX2 object of class \e ChainMesh_Catalogue:
	 * the chain mesh used to count the pairs relative to the second
	 * side of the triangle (1-3)
	 *
	 * @param tt pointer to an object of class \e Triplet
	 *
	 * @param tcount 1 &rarr; activate the CPU time counter; 0 &rarr; no time counter
	 * @return none
	 *
	 * @warning the angular three-point correlation function is not
	 * implemented yet
	 */
	void count_triplets (const std::shared_ptr<catalogue::Catalogue> cat1, const chainmesh::ChainMesh_Catalogue &ChainMesh_rMAX1, const chainmesh::ChainMesh_Catalogue &ChainMesh_rMAX2, std::shared_ptr<triplets::Triplet> tt, const bool tcount=false);

	/**
	 *  @brief count the data-data-data, random-random-random,
	 *  data-data-random and data-random-random triplets, used to
	 *  construct the estimator of the three-point correlation
	 *  function
	 *  
	 *  @param dir_output_triplets name of the output directory used to
	 *  store the number of triplets
	 * 
	 *  @param dir_input_triplets name of the input directories
	 *  containing the number of triplets
	 *
	 *  @param count_ddd 1 &rarr; count the data-data-data triplets;
	 *  0 &rarr; read the data-data-data triplets from a file
	 *
	 *  @param count_rrr 1 &rarr; count the random-random-random
	 *  triplets; 0 &rarr; read the random-random-random triplets
	 *  from a file
	 *
	 *  @param count_ddr 1 &rarr; count the data-data-random
	 *  triplets; 0 &rarr; read the data-data-random triplets from a
	 *  file
	 *
	 *  @param count_drr 1 &rarr; count the data-random-random
	 *  triplets; 0 &rarr; read the data-random-random triplets from
	 *  a file
	 *
	 *  @param tcount 1 &rarr; activate the CPU time counter; 0
	 *  &rarr; no time counter
	 *
	 *  @return none
	 */
	void count_allTriplets (const std::string dir_output_triplets=par::defaultString, const std::vector<std::string> dir_input_triplets={}, const bool count_ddd=true, const bool count_rrr=true, const bool count_ddr=true, const bool count_drr=true, const bool tcount=false);
            
	/**
	 * @brief method to count the number of triplets
	 *
	 * @param cat1 object of class \e Catalogue: the input catalogue
	 * to analyse
	 *
	 * @param ChainMesh_rMAX1 object of class \e
	 * ChainMesh_Catalogue: the chain mesh used to count the pairs
	 * relative to the first side of the triangle (1-2)
	 *
	 * @param ChainMesh_rMAX2 object of class \e
	 * ChainMesh_Catalogue: the chain mesh used to count the pairs
	 * relative to the second side of the triangle (1-3)
	 *
	 * @param tt pointer to an object of class \e Triplet
	 *
	 * @param tt_regions pointer vector to objects of class \e
	 * Triplet
	 *
	 * @param weight region weights
	 *
	 * @param tcount true &rarr; activate the CPU time counter;
	 * false &rarr; no time counter @return none
	 *
	 * @warning the angular three-point correlation function is
	 * not implemented yet
	 */
	void count_triplets_region (const std::shared_ptr<catalogue::Catalogue> cat1, const chainmesh::ChainMesh_Catalogue &ChainMesh_rMAX1, const chainmesh::ChainMesh_Catalogue &ChainMesh_rMAX2, std::shared_ptr<triplets::Triplet> tt, std::vector<std::shared_ptr<triplets::Triplet>> tt_regions, const std::vector<std::vector<double>> weight, const bool tcount=false);

	/**
	 *  @brief count the data-data-data, random-random-random,
	 *  data-data-random and data-random-random triplets, used to
	 *  construct the estimator of the three-point correlation
	 *  function
	 *  
	 *  @param weight region weights
	 *
	 *  @param dir_output_triplets name of the output directory used to
	 *  store the number of triplets
	 * 
	 *  @param dir_input_triplets name of the input directories
	 *  containing the number of triplets
	 *
	 *  @param count_ddd 1 &rarr; count the data-data-data triplets;
	 *  0 &rarr; read the data-data-data triplets from a file
	 *
	 *  @param count_rrr 1 &rarr; count the random-random-random
	 *  triplets; 0 &rarr; read the random-random-random triplets
	 *  from a file
	 *
	 *  @param count_ddr 1 &rarr; count the data-data-random
	 *  triplets; 0 &rarr; read the data-data-random triplets from a
	 *  file
	 *
	 *  @param count_drr 1 &rarr; count the data-random-random
	 *  triplets; 0 &rarr; read the data-random-random triplets from
	 *  a file
	 *
	 *  @param tcount 1 &rarr; activate the CPU time counter; 0
	 *  &rarr; no time counter
	 *
	 *  @return none
	 */
	void count_allTriplets_region (const std::vector<std::vector<double>> weight, const std::string dir_output_triplets=par::defaultString, const std::vector<std::string> dir_input_triplets={}, const bool count_ddd=true, const bool count_rrr=true, const bool count_ddr=true, const bool count_drr=true, const bool tcount=false);

	/**
	 *  @name Internal input/output member functions (customized in all the derived classes)
	 */
	///@{
    
	/**
	 *  @brief write the number of triplets
	 *  @param TT pointer to an object of class Triplet
	 *  @param dir output directory
	 *  @param file output file
	 *  @return none
	 */
	void write_triplets (const std::shared_ptr<triplets::Triplet> TT, const std::string dir, const std::string file) const;
 
	/**
	 *  @brief read the number of triplets
	 *  @param [out] TT pointer to an object of class Triplet
	 *  @param [in] dir input directory
	 *  @param [in] file input file
	 *  @return none
	 */
	void read_triplets (std::shared_ptr<triplets::Triplet> TT, const std::vector<std::string> dir, const std::string file);
    
	///@}


      
      public:

	/**
	 *  @name Constructors/destructors
	 */
	///@{

	/**
	 * @brief default constructor
	 * @return object of class ThreePointCorrelation
	 */
	ThreePointCorrelation () = default;

	/**
	 *  @brief constructor
	 *  @param data object of class Catalogue containing the input
	 *  catalogue
	 *  @param random of class Catalogue containing the random data
	 *  catalogue
	 *  @return object of class ThreePointCorrelation
	 */
	ThreePointCorrelation (const catalogue::Catalogue data, const catalogue::Catalogue random) 
	  : m_data(std::make_shared<catalogue::Catalogue>(catalogue::Catalogue(std::move(data)))), m_random(std::make_shared<catalogue::Catalogue>(catalogue::Catalogue(std::move(random)))) {}

	/**
	 * @brief default destructor
	 * @return none
	 */
	virtual ~ThreePointCorrelation () = default;
      
	/**
	 *  @brief static factory used to construct three-point correlation
	 *  functions of any type
	 *
	 *  @param type the type of three-point correlation function; it
	 *  can be: ThreePType::_angular_connected_, ThreePType::_angular_reduced_,
	 *  ThreePType::_comoving_connected_, ThreePType::_comoving_reduced_
	 *
	 *  @param data object of class Catalogue containing the input
	 *  catalogue
	 *  @param random of class Catalogue containing the random data
	 *  catalogue
	 *  @param tripletType the triplet type; it can be:
	 *  TripletType::_comoving_theta_, TripletType::_comoving_side_
	 *  @param side_s the size of r<SUB>12</SUB>
	 *  @param side_u the ratio r<SUB>13</SUB>/r<SUB>12</SUB>
	 *  @param perc_increase the ratio
	 *  @param nbins the number of bins
	 *  @return a pointer to an object of class
	 *  ThreePointCorrelation of a given type
	 */
	static std::shared_ptr<ThreePointCorrelation> Create (const ThreePType type, const catalogue::Catalogue data, const catalogue::Catalogue random, const triplets::TripletType tripletType, const double side_s, const double side_u, const double perc_increase, const int nbins);
      
	/**
	 *  @brief static factory used to construct three-point correlation
	 *  functions of any type
	 *
	 *  @param type the type of three-point correlation function; it
	 *  can be: ThreePType::_angular_connected_, ThreePType::_angular_reduced_,
	 *  ThreePType::_comoving_connected_, ThreePType::_comoving_reduced_
	 *
	 *  @param data object of class Catalogue containing the input
	 *  catalogue
	 *  @param random of class Catalogue containing the random data
	 *  catalogue
	 *  @param tripletType the triplet type; it can be:
	 *  TripletType::_comoving_theta_, TripletType::_comoving_side_
	 *  @param r12 the size of r<SUB>12</SUB>
	 *  @param r12_binSize the size of r<SUB>12</SUB> bin
	 *  @param r13 the size of r<SUB>13</SUB>
	 *  @param r13_binSize the size of r<SUB>13</SUB> bin
	 *  @param nbins the number of bins
	 *  @return a pointer to an object of class
	 *  ThreePointCorrelation of a given type
	 */
	static std::shared_ptr<ThreePointCorrelation> Create (const ThreePType type, const catalogue::Catalogue data, const catalogue::Catalogue random, const triplets::TripletType tripletType, const double r12, const double r12_binSize, const double r13, const double r13_binSize, const int nbins);
      
	///@}
      

	/**
	 *  @name Member functions to get the private/protected members
	 */
	///@{

	/**
	 *  @brief get the protected member m_threePType
	 *  @return the three-point correlation function type
	 */
	ThreePType threePType () const { return m_threePType; }
	
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
	 *  @brief get the protected member m_ddd
	 *  @return the number of data-data-data triplets
	 */
	std::shared_ptr<triplets::Triplet> ddd () const { return m_ddd; }

	/**
	 *  @brief get the protected member m_rrr
	 *  @return the number of random-random-random triplets
	 */
	std::shared_ptr<triplets::Triplet> rrr () const { return m_rrr; }

	/**
	 *  @brief get the protected member m_ddr
	 *  @return the number of data-data-random triplets
	 */
	std::shared_ptr<triplets::Triplet> ddr () const { return m_ddr; }

	/**
	 *  @brief get the protected member m_drr
	 *  @return the number of data-random-random triplets
	 */
	std::shared_ptr<triplets::Triplet> drr () const { return m_drr; }

	/**
	 *  @brief get the protected member m_scale
	 *  @return the scale bins
	 */
	virtual std::vector<double> scale () const
	{ cbl::ErrorCBL("", "scale", "ThreePointCorrelation.h"); std::vector<double> vv; return vv; }
      
	/**
	 *  @brief get the protected member m_zeta
	 *  @return the binned connected three-point correlation
	 *  function
	 */
	virtual std::vector<double> zeta () const  
	{ cbl::ErrorCBL("", "zeta", "ThreePointCorrelation.h"); std::vector<double> vv; return vv; }

	/**
	 *  @brief get the protected member m_QQ
	 *  @return the binned reduced three-point correlation function
	 */
	virtual std::vector<double> QQ () const  
	{ cbl::ErrorCBL("", "QQ", "ThreePointCorrelation.h"); std::vector<double> vv; return vv; }

	/**
	 *  @brief get the protected member
	 *  ThreePointCorrelation_comoving_connected::m_error
	 *  @return the error on the connected three-point
	 *  correlation function
	 */
	virtual std::vector<double> error () const  
	{ cbl::ErrorCBL("", "error", "ThreePointCorrelation.h"); std::vector<double> vv; return vv; }

      
	///@}

      
	/**
	 *  @name Member functions to set protected members
	 */
	///@{
    
	/**
	 *  @brief add a data catalogue
	 *  @param data object of class Catalogue 
	 *  @return none
	 */
	void set_data (const catalogue::Catalogue data) { m_data = std::make_shared<catalogue::Catalogue>(catalogue::Catalogue(std::move(data))); }

	/**
	 *  @brief add a random catalogue
	 *  @param random object of class Catalogue 
	 *  @return none
	 */
	void set_random (const catalogue::Catalogue random) { m_random = std::make_shared<catalogue::Catalogue>(catalogue::Catalogue(std::move(random))); }

	///@}

      
	/**
	 *  @name Member functions to measure the threep-point correlation function
	 */
	///@{

	/**
	 * @brief method to measure the three-point correlation function
	 *
	 * @param dir_output_triplets name of the output directory used to
	 * store the number of triplets
	 * 
	 * @param dir_input_triplets name of the input directories
	 * containing the number of triplets
	 *
	 * @param count_ddd 1 &rarr; count the data-data-data
	 * triplets; 0 &rarr; read the data-data-data triplets
	 * from a file
	 *
	 * @param count_rrr 1 &rarr; count the random-random-random
	 * triplets; 0 &rarr; read the random-random-random triplets
	 * from a file
	 *
	 * @param count_ddr 1 &rarr; count the data-data-random
	 * triplets; 0 &rarr; read the data-data-random triplets
	 * from a file
	 *
	 * @param count_drr 1 &rarr; count the data-random-random
	 * triplets; 0 &rarr; read the data-random-random triplets
	 * from a file
	 *
	 * @param tcount 1 &rarr; activate the CPU time counter; 0
	 * &rarr; no time counter
	 *
	 * @param seed the seed for random number generation
	 *
	 * @return none
	 */
	virtual void measure (const std::string dir_output_triplets, const std::vector<std::string> dir_input_triplets={}, const bool count_ddd=true, const bool count_rrr=true, const bool count_ddr=true, const bool count_drr=true, const bool tcount=false, const int seed=3213)
	{ (void)dir_output_triplets; (void)dir_input_triplets; (void)count_ddd; (void)count_rrr; (void)count_ddr; (void)count_drr; (void)tcount; (void)seed; cbl::ErrorCBL("", "measure", "ThreePointCorrelation.h"); }

	/**
	 * @brief method to measure the three-point correlation function
	 *
	 * @param weight region weights
	 *
	 * @param doJK &rarr; normalize to 1/(n-1); 1 &rarr; normalize
	 *  to n-1/n (for Jackknife)
	 *
	 * @param dir_output_triplets name of the output directory used to
	 * store the number of triplets
	 * 
	 * @param dir_input_triplets name of the input directories
	 * containing the number of triplets
	 *
	 * @param count_ddd 1 &rarr; count the data-data-data
	 * triplets; 0 &rarr; read the data-data-data triplets
	 * from a file
	 *
	 * @param count_rrr 1 &rarr; count the random-random-random
	 * triplets; 0 &rarr; read the random-random-random triplets
	 * from a file
	 *
	 * @param count_ddr 1 &rarr; count the data-data-random
	 * triplets; 0 &rarr; read the data-data-random triplets
	 * from a file
	 *
	 * @param count_drr 1 &rarr; count the data-random-random
	 * triplets; 0 &rarr; read the data-random-random triplets
	 * from a file
	 *
	 * @param tcount 1 &rarr; activate the CPU time counter; 0
	 * &rarr; no time counter
	 *
	 * @param seed the seed for random number generation
	 *
	 * @return none
	 */
	virtual void measure (const std::vector<std::vector<double>> weight, const bool doJK, const std::string dir_output_triplets=par::defaultString, const std::vector<std::string> dir_input_triplets={}, const bool count_ddd=true, const bool count_rrr=true, const bool count_ddr=true, const bool count_drr=true, const bool tcount=false, const int seed=3213)
	{ (void)weight; (void)doJK; (void)dir_output_triplets; (void)dir_input_triplets; (void)count_ddd; (void)count_rrr; (void)count_ddr; (void)count_drr; (void)tcount; (void)seed; cbl::ErrorCBL("", "measure", "ThreePointCorrelation.h"); }

	/**
	 * @brief method to measure the three-point correlation function
	 *
	 * @param errorType type of error 
	 *
	 * @param dir_output_triplets name of the output directory used to
	 * store the number of triplets
	 * 
	 * @param dir_input_triplets name of the input directories
	 * containing the number of triplets
	 *
	 * @param nResamplings number of resamplings
	 *
	 * @param count_ddd 1 &rarr; count the data-data-data
	 * triplets; 0 &rarr; read the data-data-data triplets
	 * from a file
	 *
	 * @param count_rrr 1 &rarr; count the random-random-random
	 * triplets; 0 &rarr; read the random-random-random triplets
	 * from a file
	 *
	 * @param count_ddr 1 &rarr; count the data-data-random
	 * triplets; 0 &rarr; read the data-data-random triplets
	 * from a file
	 *
	 * @param count_drr 1 &rarr; count the data-random-random
	 * triplets; 0 &rarr; read the data-random-random triplets
	 * from a file
	 *
	 * @param tcount 1 &rarr; activate the CPU time counter; 0
	 * &rarr; no time counter
	 *
	 * @param seed the seed for random number generation
	 *
	 * @return none
	 */
	virtual void measure (const ErrorType errorType, const std::string dir_output_triplets=par::defaultString, const std::vector<std::string> dir_input_triplets={}, const int nResamplings=100, const bool count_ddd=true, const bool count_rrr=true, const bool count_ddr=true, const bool count_drr=true, const bool tcount=false, const int seed=3213)
	{ (void)errorType; (void)dir_output_triplets; (void)dir_input_triplets; (void)nResamplings; (void)count_ddd; (void)count_rrr; (void)count_ddr; (void)count_drr; (void)tcount; (void)seed; cbl::ErrorCBL("", "measure", "ThreePointCorrelation.h"); }

	/**
	 * @brief method to measure the three-point correlation function
	 *
	 * @param dir_output_triplets name of the output directory used to
	 * store the number of triplets
	 * 
	 * @param dir_output_2pt name of the output directory used to
	 * store the two-point correlation functions
	 *
	 * @param dir_input_triplets name of the input directories
	 * containing the number of triplets
	 *
	 * @param count_ddd 1 &rarr; count the data-data-data
	 * triplets; 0 &rarr; read the data-data-data triplets
	 * from a file
	 *
	 * @param count_rrr 1 &rarr; count the random-random-random
	 * triplets; 0 &rarr; read the random-random-random triplets
	 * from a file
	 *
	 * @param count_ddr 1 &rarr; count the data-data-random
	 * triplets; 0 &rarr; read the data-data-random triplets
	 * from a file
	 *
	 * @param count_drr 1 &rarr; count the data-random-random
	 * triplets; 0 &rarr; read the data-random-random triplets
	 * from a file
	 *
	 * @param tcount 1 &rarr; activate the CPU time counter; 0
	 * &rarr; no time counter
	 *
	 * @param seed the seed for random number generation
	 *
	 * @return none
	 */
	virtual void measure (const std::string dir_output_triplets, const std::string dir_output_2pt, const std::vector<std::string> dir_input_triplets={}, const bool count_ddd=true, const bool count_rrr=true, const bool count_ddr=true, const bool count_drr=true, const bool tcount=false, const int seed=3213)
	{ (void)dir_output_triplets; (void)dir_output_2pt; (void)dir_input_triplets; (void)count_ddd; (void)count_rrr; (void)count_ddr; (void)count_drr; (void)tcount; (void)seed; cbl::ErrorCBL("", "measure", "ThreePointCorrelation.h!"); }
                
	/**
	 * @brief method to measure the three-point correlation function
	 *
	 * @param weight region weights
	 *
	 * @param doJK &rarr; normalize to 1/(n-1); 1 &rarr; normalize
	 *  to n-1/n (for Jackknife)
	 *
	 * @param dir_output_triplets name of the output directory
	 * used to store the number of triplets
	 * 
	 * @param dir_output_2pt name of the output directory used to
	 * store the two-point correlation functions
	 *
	 * @param dir_input_triplets name of the input directories
	 * containing the number of triplets
	 *
	 * @param count_ddd 1 &rarr; count the data-data-data
	 * triplets; 0 &rarr; read the data-data-data triplets
	 * from a file
	 *
	 * @param count_rrr 1 &rarr; count the random-random-random
	 * triplets; 0 &rarr; read the random-random-random triplets
	 * from a file
	 *
	 * @param count_ddr 1 &rarr; count the data-data-random
	 * triplets; 0 &rarr; read the data-data-random triplets
	 * from a file
	 *
	 * @param count_drr 1 &rarr; count the data-random-random
	 * triplets; 0 &rarr; read the data-random-random triplets
	 * from a file
	 *
	 * @param tcount 1 &rarr; activate the CPU time counter; 0
	 * &rarr; no time counter
	 *
	 * @param seed the seed for random number generation
	 *
	 * @return none
	 */
	virtual void measure (const std::vector<std::vector<double>> weight, const bool doJK, const std::string dir_output_triplets, const std::string dir_output_2pt, const std::vector<std::string> dir_input_triplets={}, const bool count_ddd=true, const bool count_rrr=true, const bool count_ddr=true, const bool count_drr=true, const bool tcount=false, const int seed=3213)
	{ (void)weight; (void)doJK; (void)dir_output_triplets; (void)dir_output_2pt; (void)dir_input_triplets; (void)count_ddd; (void)count_rrr; (void)count_ddr; (void)count_drr; (void)tcount; (void)seed; cbl::ErrorCBL("", "measure", "ThreePointCorrelation.h"); }

	/**
	 * @brief method to measure the three-point correlation function
	 *
	 * @param errorType type of error 
	 *
	 * @param dir_output_triplets name of the output directory used to
	 * store the number of triplets
	 * 
	 * @param dir_output_2pt name of the output directory used to
	 * store the two-point correlation functions
	 *
	 * @param dir_input_triplets name of the input directories
	 * containing the number of triplets
	 *
	 * @param nResamplings number of resamplings
	 *
	 * @param count_ddd 1 &rarr; count the data-data-data
	 * triplets; 0 &rarr; read the data-data-data triplets
	 * from a file
	 *
	 * @param count_rrr 1 &rarr; count the random-random-random
	 * triplets; 0 &rarr; read the random-random-random triplets
	 * from a file
	 *
	 * @param count_ddr 1 &rarr; count the data-data-random
	 * triplets; 0 &rarr; read the data-data-random triplets
	 * from a file
	 *
	 * @param count_drr 1 &rarr; count the data-random-random
	 * triplets; 0 &rarr; read the data-random-random triplets
	 * from a file
	 *
	 * @param tcount 1 &rarr; activate the CPU time counter; 0
	 * &rarr; no time counter
	 *
	 * @param seed the seed for random number generation
	 *
	 * @return none
	 */
	virtual void measure (const ErrorType errorType, const std::string dir_output_triplets, const std::string dir_output_2pt, const std::vector<std::string> dir_input_triplets={}, const int nResamplings=100, const bool count_ddd=true, const bool count_rrr=true, const bool count_ddr=true, const bool count_drr=true, const bool tcount=false, const int seed=3213)
	{ (void)errorType; (void)dir_output_triplets; (void)dir_output_2pt; (void)dir_input_triplets; (void)nResamplings; (void)count_ddd; (void)count_rrr; (void)count_ddr; (void)count_drr; (void)tcount; (void)seed; cbl::ErrorCBL("", "measure", "ThreePointCorrelation.h"); }

	///@}
      

	/**
	 *  @name Input/Output member functions (customized in all the derived classes)
	 */
	///@{
    
	/**
	 *  @brief write the measured three-point correlation
	 *  @param dir output directory
	 *  @param file output file
	 *  @return none
	 */
	virtual void write (const std::string dir, const std::string file) const
	{ (void)dir; (void)file; cbl::ErrorCBL("", "write", "ThreePointCorrelation.h"); }
      
	/**
	 *  @brief write the measured three-point correlation
	 *  @param dir output directory
	 *  @param file output file
	 *  @param connected 0 &rarr; write the reducted 3pt correlation
	 *  function; 1 &rarr; write both the reduced and connected 3pt
	 *  correlation function
	 *  @return none
	 */
	virtual void write (const std::string dir, const std::string file, const bool connected) const
	{ (void)dir; (void)file; (void)connected; cbl::ErrorCBL("", "write", "ThreePointCorrelation.h"); }

	virtual void write_jackknife_resampling (const std::string dir) const
	{ (void)dir;}

	/**
	 *  @brief write the measured three-point correlation covariance
	 *  @param dir output directory
	 *  @param file output file
	 *  @return none
	 */
	virtual void write_covariance (const std::string dir, const std::string file) const = 0;

	///@}
      
      };
    }
  }
}

#endif
