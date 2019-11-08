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
 *  @file Headers/ThreePointCorrelation_comoving_multipoles_single.h
 *
 *  @brief The class ThreePointCorrelation_comoving_multipoles_single
 *
 *  This file defines the interface of the class
 *  ThreePointCorrelation_comoving_multipoles_single,
 *  used to measure the legendre coefficients of the
 *  three-point correlation function for a single triangle configuration
 *
 *  @authors Federico Marulli, Michele Moresco, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, michele.moresco@unibo.it,
 *  alfonso.veropalumbo@unibo.it
 */

#ifndef __THREEPOINTMULTS__
#define __THREEPOINTMULTS__ 


#include "ThreePointCorrelation_comoving_multipoles.h"


// ===================================================================================================


namespace cbl {

  namespace measure {

    namespace threept {

      /**
       *  @class ThreePointCorrelation_comoving_multipoles_single 
       *  ThreePointCorrelation_comoving_multipoles_single.h
       * "Headers/ThreePointCorrelation_comoving_multipoles_single.h"
       *
       *  @brief The class ThreePointCorrelation_comoving_multipoles_single
       *
       *  This is the class used to measure the three-point
       *  correlation function for single triangle configuration
       */
      class ThreePointCorrelation_comoving_multipoles_single : public ThreePointCorrelation_comoving_multipoles {

	protected :

	  /**
	   *  @name variables for triangles
	   */
	  ///@{

	  /// First side minimum scale
	  double m_r12Min;

	  /// First side maximum scale
	  double m_r12Max;

	  /// Second side minimum scale
	  double m_r13Min;

	  /// Second side maximum scale
	  double m_r13Max;

	  ///@}


	  /**
	   * @brief compute the triples using
	   * Slepian, Eisenstein 2015 approach
	   *
	   * This function computes the multipoles expansion
	   * of the binned triplet counts from a the joined data+random 
	   * catalogue. The random catalogue objects are weighted by
	   * \f[ \hat{w}_i = -\frac{N_D}{N_R} w_i; \f]
	   *
	   * with \f$ N_D, N_R \f$ the weighted number of data and random objects
	   * and \f$ w_i \f$ the original weight of the random object.
	   *
	   * The multipoles expansion of triplets in the joined catalogue corresponds
	   * to the numerator of the Szapudi-Szalay estimator of the three-point
	   * correlation function.
	   *
	   * By looking at negative weights we can compute the multipoles
	   * expansion 
	   *
	   * @param NNN the multipoles expansion of the joined catalogue
	   * triplets
	   *
	   * @param RRR the multipoles expansion of the random catalogue
	   * triplets
	   *
	   * @param r12_min the minimum separation of the first shell
	   *
	   * @param r12_max the maximum separation of the first shell
	   * 
	   * @param r13_min the minimum separation of the second shell
	   * 
	   * @param r13_max the maximum separation of the second shell
	   * 
	   * @param norders the number of multipoles, \f$ l_{max}+1 \f$
	   * 
	   * @param catalogue the catalogue
	   *
	   * @return the triplets
	   */
	  void  m_count_triplets (std::vector<double> &NNN, std::vector<double> &RRR, const double r12_min, const double r12_max, const double r13_min, const double r13_max, const int norders, const catalogue::Catalogue& catalogue) const;

	  /**
	   *  @name Internal input/output member functions (customized in all the derived classes)
	   */
	  ///@{

	  /**
	   *  @brief write the number of triplets
	   *  @param TL pointer to an object of class Triplet
	   *  @param dir output directory
	   *  @param file output file
	   *  @return none
	   */
	  void m_write_triplets (const std::vector<double> TL, const std::string dir, const std::string file) const;

	  /**
	   *  @brief read the number of triplets
	   *  @param [out] TL pointer to an object of class Triplet
	   *  @param [in] dir input directory
	   *  @param [in] file input file
	   *  @return none
	   */
	  void m_read_triplets (std::vector<double> &TL, const std::vector<std::string> dir, const std::string file);

	  ///@}


	public:

	  /**
	   * @brief default constructor
	   *
	   * @return object of type ThreePointCorrelation_comoving_multipoles
	   */
	  ThreePointCorrelation_comoving_multipoles_single () {}

	  /**
	   * @brief constructor of ThreePointCorrelation_comoving_multipoles
	   *
	   * @details constructor of ThreePointCorrelation_comoving_multipoles
	   * to compute single configuration three-point correlation function multipoles
	   *
	   * @param catalogue the data catalogue
	   * @param random_catalogue the random catalogue
	   * @param r12Min the minimum triangle first side
	   * @param r12Max the maximum triangle first side
	   * @param r13Min the minimum triangle second side
	   * @param r13Max the maximum triangle second side
	   * @param nOrders the number of Legendre multipoles
	   * @param split factor to split the random sample. 
	   * 	it must be a multiple m_data.nobjects()
	   * @param seed seed to shuffle the random sample
	   *
	   * @return object of type ThreePointCorrelation_comoving_multipoles
	   * @warning this function will raise an error if m_random.nObjects() < split*m_data.nObjects
	   * if m_random.nObjects() > split*m_data.nObjects, only random points up to split*m_data.nObjects
	   * Negative values of the split factor allow to use the whole random sample.
	   */
	  ThreePointCorrelation_comoving_multipoles_single (cbl::catalogue::Catalogue catalogue, cbl::catalogue::Catalogue random_catalogue, const double r12Min, const double r12Max, const double r13Min, const double r13Max, const size_t nOrders, const double split=-1, const size_t seed=234);

	  /**
	   * @brief default destructor
	   *
	   * @return None
	   */
	  ~ThreePointCorrelation_comoving_multipoles_single () {}

	  /**
	   * @brief set parameters for single configuration
	   * three-point correlation function multipoles
	   *
	   * @param r12Min the minimum triangle first side
	   * @param r12Max the maximum triangle first side
	   * @param r13Min the minimum triangle second side
	   * @param r13Max the maximum triangle second side
	   * @param nOrders the number of Legendre multipoles
	   *
	   * @return None
	   */
	  void set_parameters (const double r12Min, const double r12Max, const double r13Min, const double r13Max, const size_t nOrders);

	  /**
	   * @brief measure the three-point correlation function multipoles
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
	   * @param count_triplets 1 &rarr; count the triplets
	   * triplets; 0 &rarr; read the triplets from a file
	   *
	   * @param tcount 1 &rarr; activate the CPU time counter; 0
	   * &rarr; no time counter
	   *
	   * @param seed the seed for random number generation
	   *
	   * @return none
	   *
	   * @warning no error have been implemented so far, any choice will
	   * be ignored.
	   */
	  void measure (const ErrorType errorType, const std::string dir_output_triplets=par::defaultString, const std::vector<std::string> dir_input_triplets={}, const int nResamplings=100, const bool count_triplets=true, const bool tcount=false, const int seed=3213);

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

	  void write (const std::string dir, const std::string file) const;

	  ///@}

	  /**
	   *  @name members function to resum the triplet counts
	   */
	  ///@{

	  /**
	   *  @brief resum the three-point correlation function, write 
	   *  output in file
	   *  @param dir output directory
	   *  @param file output file
	   *  @param tripletType the triplet type
	   *  @param nBins the number of bins
	   *  @param bin true \f$\rightarrow\f$ average legendre polynomials,
	   *  	false \f$\rightarrow\f$ compute legendre polynomial at the bin center
	   *  @return none
	   */

	  void resum (const std::string dir, const std::string file, const cbl::triplets::TripletType tripletType, const int nBins, const bool bin=true) const;

	  ///@}

      };
    }
  }
}

#endif 
