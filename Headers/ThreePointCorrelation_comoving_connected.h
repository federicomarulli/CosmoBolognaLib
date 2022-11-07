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
 *  @file Headers/ThreePointCorrelation_comoving_connected.h
 *
 *  @brief The class ThreePointCorrelation_comoving_connected
 *
 *  This file defines the interface of the class
 *  ThreePointCorrelation_comoving_connected, used to measure the
 *  connected three-point correlation function in comoving coordinates
 *
 *  @authors Federico Marulli, Michele Moresco, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unibo.it, michele.moresco@unibo.it,
 *  alfonso.veropalumbo@unibo.it
 */

#ifndef __THREEPOINTCOMCON__
#define __THREEPOINTCOMCON__ 


#include "ThreePointCorrelation.h"


// ===================================================================================================


namespace cbl {

  namespace measure {

    namespace threept {
    
      /**
       *  @class ThreePointCorrelation_comoving_connected
       * ThreePointCorrelation_comoving_connected.h
       * "Headers/ThreePointCorrelation_comoving_connected.h"
       *
       *  @brief The class ThreePointCorrelation_comoving_connected
       *
       *  This is the base class used to measure the connected three-point
       *  correlation function in comoving coordinates
       */
      class ThreePointCorrelation_comoving_connected : public ThreePointCorrelation {

      protected :
    
	/**
	 *  @name Three-point correlation function data
	 */
	///@{
    
	/// scale bins
	std::vector<double> m_scale;

	/// binned connected three-point correlation function
	std::vector<double> m_zeta;

	/// error on the binned connected three-point correlation function
	std::vector<double> m_error;
    
	///@}

    
      public:
    
	/**
	 *  @name Constructors/destructors
	 */
	///@{
	
	/**
	 *  @brief default constructor
	 */
	ThreePointCorrelation_comoving_connected () = default;

	/**
	 *  @brief constructor
	 *  @param data object of class Catalogue containing the input
	 *  catalogue
	 *  @param random of class Catalogue containing the random data
	 *  catalogue
	 *  @param tripletType the triplet type; it can be:
	 *  TripletType::_comoving_theta_, TripletType::_comoving_side_
	 *  @param side_s the size of r<SUB>12</SUB>
	 *  @param side_u the ratio r<SUB>13</SUB>/r<SUB>12</SUB>
	 *  @param perc_increase the ratio
	 *  &Delta;r<SUB>12</SUB>/r<SUB>12</SUB>=&Delta;r<SUB>13</SUB>/r<SUB>13</SUB>
	 *  @param nbins number of bins
	 */
	ThreePointCorrelation_comoving_connected (const catalogue::Catalogue data, const catalogue::Catalogue random, const triplets::TripletType tripletType, const double side_s, const double side_u, const double perc_increase, const int nbins)
	  : ThreePointCorrelation(data, random) { set_parameters(tripletType, side_s, side_u, perc_increase, nbins); }

	/**
	 *  @brief constructor
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
	 *  @param nbins number of bins
	 */
	ThreePointCorrelation_comoving_connected (const catalogue::Catalogue data, const catalogue::Catalogue random, const triplets::TripletType tripletType, const double r12, const double r12_binSize, const double r13, const double r13_binSize, const int nbins)
	  : ThreePointCorrelation(data, random) { set_parameters(tripletType, r12, r12_binSize, r13, r13_binSize, nbins); }

	/**
	 *  @brief default destructor
	 */
	~ThreePointCorrelation_comoving_connected () = default;

	///@}

      
	/**
	 *  @name Member functions to set the binning parameters
	 */
	///@{
      
	/**
	 *  @brief set the binning parameters
	 *  @param tripletType the triplet type; it can be:
	 *  TripletType::_comoving_theta_, TripletType::_comoving_side_
	 *  @param side_s the size of r<SUB>12</SUB>
	 *  @param side_u the ratio r<SUB>13</SUB>/r<SUB>12</SUB>
	 *  @param perc_increase the ratio
	 *  &Delta;r<SUB>12</SUB>/r<SUB>12</SUB>=&Delta;r<SUB>13</SUB>/r<SUB>13</SUB>
	 *  @param nbins number of bins
	 */
	void set_parameters (const triplets::TripletType tripletType, const double side_s, const double side_u, const double perc_increase, const int nbins);

	/**
	 *  @brief set the binning parameters
	 *  @param tripletType the triplet type; it can be:
	 *  TripletType::_comoving_theta_, TripletType::_comoving_side_
	 *  @param r12 the size of r<SUB>12</SUB>
	 *  @param r12_binSize the size of r<SUB>12</SUB> bin
	 *  @param r13 the size of r<SUB>13</SUB>
	 *  @param r13_binSize the size of r<SUB>13</SUB> bin
	 *  @param nbins number of bins
	 */
	void set_parameters (const triplets::TripletType tripletType, const double r12, const double r12_binSize, const double r13, const double r13_binSize, const int nbins);

	///@}

      
	/**
	 *  @name Member functions to get protected parameters
	 */
	///@{
	
	/**
	 *  @brief get the protected member
	 *  ThreePointCorrelation_comoving_connected::m_scale
	 *  @return the scale bins
	 */
	std::vector<double> scale () const override { return m_scale; }

	/**
	 *  @brief get the protected member
	 *  ThreePointCorrelation_comoving_connected::m_zeta
	 *  @return the binned connected three-point correlation
	 *  function
	 */
	std::vector<double> zeta () const override { return m_zeta; }

	/**
	 *  @brief get the protected member
	 *  ThreePointCorrelation_comoving_connected::m_error
	 *  @return the error on the connected three-point
	 *  correlation function
	 */
	std::vector<double> error () const override { return m_error; }

	///@}
	
	/**
	 *  @name Member functions to measure the three-point correlation function
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
	 * @param fact factor used to compute the cell size of the
	 * chain mesh: it is multiplied by the maximum distance
	 * considered for the couples and can be setted by the user
	 * to optimize the count of the couples
	 *
	 * @param seed the seed for random number generation
	 */
	void measure (const std::string dir_output_triplets, const std::vector<std::string> dir_input_triplets={}, const bool count_ddd=true, const bool count_rrr=true, const bool count_ddr=true, const bool count_drr=true, const bool tcount=true, const double fact=0.1, const int seed=3213) override;

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
	 * @param fact factor used to compute the cell size of the
	 * chain mesh: it is multiplied by the maximum distance
	 * considered for the couples and can be setted by the user
	 * to optimize the count of the couples
	 *
	 * @param seed the seed for random number generation
	 */
	void measure (const std::vector<std::vector<double>> weight, const bool doJK, const std::string dir_output_triplets=par::defaultString, const std::vector<std::string> dir_input_triplets={}, const bool count_ddd=true, const bool count_rrr=true, const bool count_ddr=true, const bool count_drr=true, const bool tcount=true, const double fact=0.1, const int seed=3213) override;

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
	 * @param fact factor used to compute the cell size of the
	 * chain mesh: it is multiplied by the maximum distance
	 * considered for the couples and can be setted by the user
	 * to optimize the count of the couples
	 *
	 * @param seed the seed for random number generation
	 */
	void measure (const ErrorType errorType, const std::string dir_output_triplets, const std::vector<std::string> dir_input_triplets={}, const int nResamplings = 100, const bool count_ddd=true, const bool count_rrr=true, const bool count_ddr=true, const bool count_drr=true, const bool tcount=true, const double fact=0.1, const int seed=3213) override;
	
	///@}
    
	/**
	 *  @name Input/Output methods
	 */
	///@{
    
	/**
	 *  @brief write the monopole of the two-point correlation
	 *  function
	 *  @param dir output directory
	 *  @param file output file
	 */
	void write (const std::string dir, const std::string file) const override;

        /**
         *  @brief write the measured three-point correlation covariance
         *  @param dir output directory
         *  @param file output file
         */
        void write_covariance (const std::string dir, const std::string file) const override;
    
	///@}
      
      };
    }
  }
}

#endif
