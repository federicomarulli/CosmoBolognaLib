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
 *  @file Headers/ThreePointCorrelation_comoving_multipoles_all.h
 *
 *  @brief The class ThreePointCorrelation_comoving_multipoles_all
 *
 *  This file defines the interface of the class 
 *  ThreePointCorrelation_comoving_multipoles,
 *  used to measure the legendre coefficients of the
 *  three-point correlation function for all triangle configuration
 *
 *  @authors Federico Marulli, Michele Moresco, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unibo.it, michele.moresco@unibo.it,
 *  alfonso.veropalumbo@unibo.it
 */

#ifndef __THREEPOINTMULTA__
#define __THREEPOINTMULTA__ 


#include "ThreePointCorrelation_comoving_multipoles.h"


// ===================================================================================================


namespace cbl {

  namespace measure {

    namespace threept {

      /**
       *  @class ThreePointCorrelation_comoving_multipoles_all
       *  ThreePointCorrelation_comoving_multipoles_all.h
       *  "Headers/ThreePointCorrelation_comoving_multipoles_all.h"
       *
       *  @brief The class
       *  ThreePointCorrelation_comoving_multipoles_all
       *
       *  This is the class used to measure the three-point
       *  correlation function multipoles for all triangle
       *  configurations
       */
      class ThreePointCorrelation_comoving_multipoles_all : public ThreePointCorrelation_comoving_multipoles {

      protected :

	/**
	 *  @name variables for triangles
	 */
	///@{
	
	/// minimum scale
	double m_rMin;

	/// minimum scale
	double m_rMax;

	/// shell size
	double m_binSize;

	/// number of bins
	size_t m_nBins;

	///@}
	
	/**
	 * @param NNN the multipoles expansion of the triplets for the
	 * separation bins
	 *
	 * @param RRR the multipoles expansion of the random triplets
	 * for the separation bins
	 *
	 * @param rmin the minimum separation
	 *
	 * @param rmax the maximum separation
	 * 
	 * @param nbins the number of separation bins 
	 * 
	 * @param norders the number of multipoles, \f$ l_{max}+1 \f$
	 * 
	 * @param catalogue the catalogue
	 */
	void m_count_triplets (std::vector<double> &NNN, std::vector<double> &RRR, const double rmin, const double rmax, const int nbins, const int norders, const catalogue::Catalogue& catalogue) const;

	/**
	 *  @name Internal input/output member functions (customized in all the derived classes)
	 */
	///@{
	
	/**
	 *  @brief write the number of triplets
	 *  @param TL pointer to an object of class Triplet
	 *  @param dir output directory
	 *  @param file output file
	 */
	void m_write_triplets (const std::vector<double> TL, const std::string dir, const std::string file) const;

	/**
	 *  @brief read the number of triplets
	 *  @param [out] TL pointer to an object of class Triplet
	 *  @param [in] dir input directory
	 *  @param [in] file input file
	 */
	void m_read_triplets (std::vector<double> &TL, const std::vector<std::string> dir, const std::string file);

	///@}

      public:

	/**
	 * @brief default constructor
	 */
	ThreePointCorrelation_comoving_multipoles_all () = default;

	/**
	 * @brief constructor
	 *
	 * @details constructor of
	 * ThreePointCorrelation_comoving_multipoles to compute all
	 * configurations three-point correlation function multipoles
	 *
	 * @param catalogue the data catalogue
	 * @param random_catalogue the random catalogue
	 * @param rMin the minimum triangle side
	 * @param rMax the maximum triangle side
	 * @param binSize the triangle side width
	 * @param nOrders the number of Legendre multipoles
	 * @param split factor to split the random sample. 
	 * 	it must be a multiple m_data.nobjects()
	 * @param seed seed to shuffle the random sample
	 *
	 * @warning this function will raise an error if
	 * m_random.nObjects() < split*m_data.nObjects if
	 * m_random.nObjects() > split*m_data.nObjects, only random
	 * points up to split*m_data.nObjects Negative values of the
	 * split factor allow to use the whole random sample.
	 */
	ThreePointCorrelation_comoving_multipoles_all (const cbl::catalogue::Catalogue catalogue, const cbl::catalogue::Catalogue random_catalogue, const double rMin, const double rMax, const double binSize, const size_t nOrders, const double split=-1, const size_t seed=234);

	/**
	 * @brief constructor
	 *
	 * @details constructor of
	 * ThreePointCorrelation_comoving_multipoles_all to compute
	 * all configurations three-point correlation function
	 * multipoles, rebinning previously measured triplets
	 *
	 * @param threept object of type
	 * ThreePointCorrelation_comoving_multipoles_all
	 *
	 * @param newBinSize the triangle side width
	 *
	 * @warning this function rebin Legendre triplet coefficients
	 * previously computed. The new bin size could disagree with
	 * the original triangle sides.  In that case the minimum
	 * separation will be changed accordingly.
	 */
	ThreePointCorrelation_comoving_multipoles_all (const ThreePointCorrelation_comoving_multipoles_all &threept, const double newBinSize);

	/**
	 * @brief default destructor
	 */
	~ThreePointCorrelation_comoving_multipoles_all () = default;

	/**
	 * @brief set parameters for all configurations 
	 * three-point correlation function multipoles
	 *
	 * @param rMin the minimum triangle side
	 * @param rMax the maximum triangle side
	 * @param binSize the triangle side width
	 * @param nOrders the number of Legendre multipoles
	 */
	void set_parameters (const double rMin, const double rMax, const double binSize, const size_t nOrders);

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
	 * @warning no error have been implemented so far, any choice
	 * will be ignored
	 */
	void measure (const ErrorType errorType, const std::string dir_output_triplets, const std::vector<std::string> dir_input_triplets, const int nResamplings, const bool count_triplets, const bool tcount, const int seed) override;

	/**
	 *  @name Input/Output member functions (customized in all the derived classes)
	 */
	///@{
	
	/**
	 *  @brief write the measured three-point correlation
	 *  @param dir output directory
	 *  @param file output file
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
	 *
	 *  @param dir output directory
	 *  @param file output file
	 *  @param tripletType the triplet type
	 *  @param nBins the number of bins
	 *
	 *  @param bin true \f$\rightarrow\f$ average legendre
	 *  polynomials, false \f$\rightarrow\f$ compute legendre
	 *  polynomial at the bin center
	 */
	void resum (const std::string dir, const std::string file, const cbl::triplets::TripletType tripletType, const int nBins, const bool bin=true) const;

	///@}

      };
    }
  }
}

#endif 
