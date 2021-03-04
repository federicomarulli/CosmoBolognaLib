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
 *  @file Headers/ThreePointCorrelation_angular_connected.h
 *
 *  @brief The class ThreePointCorrelation_angular_connected
 *
 *  This file defines the interface of the class
 *  ThreePointCorrelation_angular_connected, used to measure the
 *  connected three-point correlation function in angular coordinates
 *
 *  @authors Federico Marulli, Michele Moresco, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unibo.it, michele.moresco@unibo.it,
 *  alfonso.veropalumbo@unibo.it
 */

#ifndef __THREEPOINTANGCON__
#define __THREEPOINTANGCON__ 


#include "ThreePointCorrelation.h"


// ===================================================================================================


namespace cbl {
  
  namespace measure {

    namespace threept {
    
      /**
       *  @class ThreePointCorrelation_angular_connected
       *  ThreePointCorrelation_angular_connected.h
       *  "Headers/ThreePointCorrelation_angular_connected.h"
       *
       *  @brief The class ThreePointCorrelation_angular_connected
       *
       *  This is the base class used to measure the connected three-point
       *  correlation function in angular coordinates
       *
       *  @warning This class has not been fully implemented yet
       */
      class ThreePointCorrelation_angular_connected : public ThreePointCorrelation {

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
	ThreePointCorrelation_angular_connected () = default;
	
	/**
	 *  @brief constructor
	 *  @param data object of class Catalogue containing the input
	 *  catalogue
	 *  @param random of class Catalogue containing the random data
	 *  catalogue
	 *  @param side_s the size of r<SUB>12</SUB>
	 *  @param side_u the ratio r<SUB>13</SUB>/r<SUB>12</SUB>
	 *  @param perc_increase the ratio
	 *  &Delta;r<SUB>12</SUB>/r<SUB>12</SUB>=&Delta;r<SUB>13</SUB>/r<SUB>13</SUB>
	 *  @param nbins number of bins
	 */
	ThreePointCorrelation_angular_connected (const catalogue::Catalogue data, const catalogue::Catalogue random, const double side_s, const double side_u, const double perc_increase, const int nbins)
	  : ThreePointCorrelation(data, random) { ErrorCBL("", "ThreePointCorrelation_angular_connected", "ThreePointCorrelation_angular_connected.h", cbl::glob::ExitCode::_workInProgress_); set_parameters(side_s, side_u, perc_increase, nbins); }
	
	/**
	 *  @brief constructor
	 *  @param data object of class Catalogue containing the input
	 *  catalogue
	 *  @param random of class Catalogue containing the random data
	 *  catalogue
	 *  @param r12 the size of r<SUB>12</SUB>
	 *  @param r12_binSize the size of r<SUB>12</SUB> bin
	 *  @param r13 the size of r<SUB>13</SUB>
	 *  @param r13_binSize the size of r<SUB>13</SUB> bin
	 *  @param nbins number of bins
	 */
	ThreePointCorrelation_angular_connected (const catalogue::Catalogue data, const catalogue::Catalogue random, const double r12, const double r12_binSize, const double r13, const double r13_binSize, const int nbins)
	  : ThreePointCorrelation(data, random) { ErrorCBL("", "ThreePointCorrelation_angular_connected", "ThreePointCorrelation_angular_connected.h", cbl::glob::ExitCode::_workInProgress_); set_parameters(r12, r12_binSize, r13, r13_binSize, nbins); }

	/**
	 *  @brief default destructor
	 */
	~ThreePointCorrelation_angular_connected () = default;

	///@}

      
	/**
	 *  @name Member functions to set the binning parameters
	 */
	///@{
      
	/**
	 *  @brief set the binning parameters
	 *  @param side_s the size of r<SUB>12</SUB>
	 *  @param side_u the ratio r<SUB>13</SUB>/r<SUB>12</SUB>
	 *  @param perc_increase the ratio
	 *  &Delta;r<SUB>12</SUB>/r<SUB>12</SUB>=&Delta;r<SUB>13</SUB>/r<SUB>13</SUB>
	 *  @param nbins number of bins
	 *
	 *  @warning This method has not been implemented yet
	 */
	void set_parameters (const double side_s, const double side_u, const double perc_increase, const int nbins);

	/**
	 *  @brief set the binning parameters
	 *  @param r12 the size of r<SUB>12</SUB>
	 *  @param r12_binSize the size of r<SUB>12</SUB> bin 
	 *  @param r13 the size of r<SUB>13</SUB>
	 *  @param r13_binSize the size of r<SUB>13</SUB> bin 
	 *  @param nbins number of bins
	 *  
	 *  @warning This method has not been implemented yet
	 */
	void set_parameters (const double r12, const double r12_binSize, const double r13, const double r13_binSize, const int nbins);

	///@}

      
	/**
	 *  @name Member functions to get protected parameters
	 */
	///@{
	
	/**
	 *  @brief get the protected member
	 *  ThreePointCorrelation_angular_connected::m_scale
	 *  @return the scale bins
	 */
	std::vector<double> scale () const override { return m_scale; }

	/**
	 *  @brief get the protected member
	 *  ThreePointCorrelation_angular_connected::m_zeta
	 *  @return the binned connected three-point correlation
	 *  function
	 */
	std::vector<double> zeta () const override { return m_zeta; }

	/**
	 *  @brief get the protected member
	 *  ThreePointCorrelation_angular_connected::m_error
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
	 * @param seed the seed for random number generation
	 *
	 * @warning This method has not been implemented yet
	 */
	void measure (const std::string dir_output_triplets, const std::vector<std::string> dir_input_triplets={}, const bool count_ddd=true, const bool count_rrr=true, const bool count_ddr=true, const bool count_drr=true, const bool tcount=false, const int seed=3213) override;
    
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
	 *  
	 *  @warning This method has not been implemented yet
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
