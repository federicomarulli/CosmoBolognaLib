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
 *  @file Headers/ThreePointCorrelation_comoving_multipoles.h
 *
 *  @brief The class ThreePointCorrelation_comoving_multipoles
 *
 *  This file defines the interface of the class 
 *  ThreePointCorrelation_comoving_multipoles,
 *  used to measure the legendre coefficients of the
 *  three-point correlation function
 *
 *  @authors Federico Marulli, Michele Moresco, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unibo.it, michele.moresco@unibo.it,
 *  alfonso.veropalumbo@unibo.it
 */

#ifndef __THREEPOINTMULT__
#define __THREEPOINTMULT__ 


#include "ThreePointCorrelation.h"


// ===================================================================================================


namespace cbl {

  namespace measure {

    namespace threept {

      /**
       *  @class ThreePointCorrelation_comoving_multipoles 
       *  ThreePointCorrelation_comoving_multipoles.h
       * "Headers/ThreePointCorrelation_comoving_multipoles.h"
       *
       *  @brief The class ThreePointCorrelation_comoving_multipoles
       *
       *  This is the base class used to measure the three-point
       *  correlation function multipoles
       */
      class ThreePointCorrelation_comoving_multipoles : public ThreePointCorrelation {
	
      protected :
	
	/**
	 *  @name variables for triangles
	 */
	///@{

	/// number of legendre polynomial in output
	size_t m_nOrders;

	///@}
	
	/**
	 *  @name Input and random catalogues
	 */
	///@{

	/// input data catalogue
	std::shared_ptr<catalogue::Catalogue> m_data;

	/// input_random catalogue
	std::shared_ptr<catalogue::Catalogue> m_random;

	/// data-random catalogue
	std::shared_ptr<catalogue::Catalogue> m_joined;

	/// factor to split the random sample. It must be a multiple m_data.nObjects()
	double m_splitFactor;

	/// m_random.nObjects/(m_data.nObjects()*m_splitFactor)
	size_t m_nSplit;

	///@}

	/**
	 *  @name triplet vectors
	 */
	///@{
	
	/// number of (data-random) triplets
	std::vector<double> m_nnn;

	/// number of random triplets
	std::vector<double> m_rrr;

	///@}
	  
	/**
	 *  @name three-point correlation vector
	 */
	///@{

	/// three-point correlation legendre coefficients
	std::vector<double> m_zeta;

	///@}

	/**
	 * @brief join data and random catalogues
	 *
	 * @details join data and random catalogues, weighting the
	 * random objects \f$w = - N_{data}/N_{random} \f$, with \f$
	 * N_{data}\f$ and \f$N_{random}\f$ are the weighted number
	 * of objects respectively in the data and random sample
	 *
	 * @param data the data catalogue
	 *
	 * @param random the random catalogue
	 *
	 * @param startPos the first object in the random catalogue
	 * to be used
	 *
	 * @param endPos the last object in the random catalogue to be used
	 *
	 * @return the joined catalogue
	 *
	 * @warning an error is raised if startPos or endPos are
	 * greater than m_random.nObjects()
	 */
	std::shared_ptr<cbl::catalogue::Catalogue> m_join_catalogues (const cbl::catalogue::Catalogue& data, const cbl::catalogue::Catalogue& random, const size_t startPos, const size_t endPos) const;

	/**
	 *  @brief compute edge-correction and return legendre
	 *  coefficients of the 3PCF
	 *
	 *  @param nnn the Data-Random triplet legendre coefficients
	 *
	 *  @param rrr the random triplet legendre coefficients
	 *
	 *  @param normalization the normalization factor: \f$ n_R^3/n_G^3\f$
	 *
	 *  @return the legendre coefficients of the 3PCF
	 */
	std::vector<double> m_SzapudiSzalay_multipoles (const std::vector<double> nnn, const std::vector<double> rrr, const double normalization=1.) const;

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
	virtual void m_write_triplets (const std::vector<double> TL, const std::string dir, const std::string file) const = 0;

	/**
	 *  @brief read the number of triplets
	 *  @param [out] TL pointer to an object of class Triplet
	 *  @param [in] dir input directory
	 *  @param [in] file input file
	 */
	virtual void m_read_triplets (std::vector<double> &TL, const std::vector<std::string> dir, const std::string file) = 0;

	///@}

      public:

	/**
	 * @brief default constructor
	 */
	ThreePointCorrelation_comoving_multipoles () = default;

	/**
	 * @brief default destructor
	 */
	~ThreePointCorrelation_comoving_multipoles () = default;

	/**
	 * @brief set the catalogues for the analysis
	 *
	 * @details set the internal objects m_data, m_random and
	 * compute the joined catalogue
	 *
	 * @param catalogue the data catalogue
	 *
	 * @param random_catalogue the random catalogue
	 *
	 * @param split factor to split the random sample. It must
	 * be a multiple m_data.nObjects()
	 *
	 * @param seed seed to shuffle the random sample
	 *
	 * @warning this function will raise an error if
	 * m_random.nObjects() < split*m_data.nObjects if
	 * m_random.nObjects() > split*m_data.nObjects, only random
	 * points up to split*m_data.nObjects. Negative values of
	 * the split factor allow to use the whole random sample.
	 */
	void set_catalogues (cbl::catalogue::Catalogue catalogue, cbl::catalogue::Catalogue random_catalogue, const double split=-1, const int seed=234);

	/**
	 * @brief set parameters for single configuration
	 * three-point correlation function multipoles
	 *
	 * @param r12Min the minimum triangle first side
	 * @param r12Max the maximum triangle first side
	 * @param r13Min the minimum triangle second side
	 * @param r13Max the maximum triangle second side
	 * @param nOrders the number of Legendre multipoles
	 */
	virtual void set_parameters (const double r12Min, const double r12Max, const double r13Min, const double r13Max, const size_t nOrders)
	{
	  (void)r12Min; (void)r12Max; (void)r13Min; (void)r13Max; (void)nOrders;
	  ErrorCBL("", "set_parameters", "ThreePointCorrelation_comoving_multipoles.cpp");
	}

	/**
	 * @brief set parameters for all configurations 
	 * three-point correlation function multipoles
	 *
	 * @param rMin the minimum triangle side
	 * @param rMax the maximum triangle side
	 * @param binSize the triangle side width
	 * @param nOrders the number of Legendre multipoles
	 */
	void set_parameters (const double rMin, const double rMax, const double binSize, const size_t nOrders)
	{
	  (void)rMin; (void)rMax; (void)binSize; (void)nOrders;
	  ErrorCBL("", "set_parameters", "ThreePointCorrelation_comoving_multipoles.h");
	}
	
	
	/**
	 *  @name Input/Output member functions (customized in all the
	 *  derived classes)
	 */
	///@{
	
	/**
	 *  @brief write the measured three-point correlation
	 *  @param dir output directory
	 *  @param file output file
	 */
	virtual void write (const std::string dir, const std::string file) const = 0;

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
	 *
	 *  @param file output file
	 *
	 *  @param tripletType the triplet type
	 *
	 *  @param nBins the number of bins
	 *
	 *  @param bin true \f$\rightarrow\f$ average legendre
	 *  polynomials, false \f$\rightarrow\f$ compute legendre
	 *  polynomial at the bin center
	 */
	virtual void resum (const std::string dir, const std::string file, const cbl::triplets::TripletType tripletType, const int nBins, const bool bin=true) const = 0;

	///@}

      };
    }
  }
}

#endif 
