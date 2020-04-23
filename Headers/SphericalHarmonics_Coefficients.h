/*******************************************************************
 *  Copyright (C) 2010 by Federico Marulli                         *
 *  federico.marulli3@unibo.it                                     *
 *                                                                 *
 *  This program is free software; you can redistribute it and/or  *
 *  modify it under the terms of the GNU General Public License as *
 *  published by the Free Software Foundation; either version 2 of *
 *  the License, or (at your option) any later version.            *
 *                                                                 *
 *  This program is distributed in the hope that it will be useful,*
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the  *
 *  GNU General Public License for more details.                   *
 *                                                                 *
 *  You should have received a copy of the GNU General Public      *
 *  License along with this program; if not, write to the Free     *
 *  Software Foundation, Inc.,                                     *
 *  59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.      *
 *******************************************************************/

/**
 *  @file Headers/SphericalHarmonics_Coefficients.h
 *
 *  @brief Generic functions that use one or more classes of the
 *  CosmoBolognaLib
 *
 *  This file contains the prototypes of the class SphericalHarmonics_Coefficients
 *
 *  @author Federico Marulli 
 *
 *  @author federico.marulli3@unibo.it
 */

#ifndef __SHCOEFF__
#define __SHCOEFF__

#include "Func.h"

namespace cbl {

 namespace glob {

    /**
     *  @class SphericalHarmonics_Coefficients SphericalHarmonics_Coefficients.h "Headers/SphericalHarmonics_Coefficients.h"
     *
     *  @brief The class SphericalHarmonics_Coefficients
     *
     *  This class is used to handle objects of type 
     *  <EM> SphericalHarmonics_Coefficients </EM>. 
     *  It contains all methods to compute coefficients of
     *  spherical harmonics expansion \f$a_{lm}\f$ in any position of the
     *  unit sphere.
     *
     *  Coefficients can be binned according to the magnitude of the separation
     *  vector and accumulated. 
     */
    class SphericalHarmonics_Coefficients {

      protected:

	/// the number of separation bins
	int m_nbins;

	/// the number of multipoles, \f$ l_{max}+1 \f$
	int m_norder;

	/// the maximum multipole \f$ l_{max} \f$
	int m_lmax;

	/// the total number of spherical harmonics
	int m_n_sph;

	/// the number of spherical harmonics for a given choice of \f$l\f$
	std::vector<int> m_n_sph_l;

	/// the spherical harmonics expansion coefficients in separation bins
	std::vector<std::vector<std::complex<double>>> m_alm;

	/// the normalization
	std::vector<double> m_normalization;

	/// vector for temporary computation of associated legendre polynomials
 	std::vector<double> m_Plm;
  	
	/// vector for temporary computation of spherical harmonics
	std::vector<std::complex<double>> m_sph;

      public:
	/**
	 *  @name Constructors/destructors
	 */
	///@{

	/**
	 * @brief default constructor
	 *
	 * @return object of type SphericalHarmonics_Coefficients
	 */
	SphericalHarmonics_Coefficients () {}

	/**
	 * @brief default constructor
	 *
	 * @param norder the number of multipoles, \f$ l_{max}+1 \f$
	 *
	 * @param nbins the number of separation bins
	 *
	 * @return object of type SphericalHarmonics_Coefficients
	 */
	SphericalHarmonics_Coefficients (const int norder, const int nbins=1) { initialize(norder, nbins); }

	/**
	 * @brief default descructor
	 *
	 * @return none
	 */
	~SphericalHarmonics_Coefficients () {}

	///@}

	/**
	 * @brief return the real part of the n-th coefficient
	 * of the expansion for a given separation bin
	 *
	 * @param n the n-th coefficient of the spherical harmonics
	 * expansion
	 *
	 * @param bin the separation bin
	 *
	 * @return  the real part of the n-th coefficient
	 * of the expansion for a given separation bin
	 */
	double real (const int n, const int bin=0) { return m_alm[bin][n].real();} 

	/**
	 * @brief return the imaginary part of the n-th coefficient
	 * of the expansion for a given separation bin
	 *
	 * @param n the n-th coefficient of the spherical harmonics
	 * expansion
	 *
	 * @param bin the separation bin
	 *
	 * @return  the imaginary part of the n-th coefficient
	 * of the expansion for a given separation bin
	 */
	double imag (const int n, const int bin=0) { return m_alm[bin][n].imag();} 

	/**
	 * @brief initialize the internal quantities
	 *
	 * @param norder the number of multipoles, \f$ l_{max}+1 \f$
	 *
	 * @param nbins the number of separation bins 
	 *
	 * @return none
	 */
	void initialize (const int norder, const int nbins=1);

	/**
	 * @brief reset the internal quantities
	 *
	 * @return none
	 */
	void reset ();

	/**
	 * @brief compute the \f$ a_{lm}\f$ for the normalized
	 * coordinates \f$ \lbrace x, y, z \rbrace \f$
	 *
	 * @param xx the x coordinate
	 *
	 * @param yy the y coordinate
	 *
	 * @param zz the z coordinate
	 *
	 * @return vector containing the \f$ a_{lm}\f$ for the 
	 * normalized coordinates \f$ \lbrace x, y, z \rbrace \f$
	 */
	std::vector<std::complex<double>> alm(const double xx, const double yy, const double zz);

	/**
	 * @brief add the \f$ a_{lm}\f$ to a specific separation
	 * bin with a weight
	 *
	 * @param alm the spherical harmonics expansion coefficients
	 *
	 * @param ww the weight
	 *
	 * @param bin the separation bin
	 *
	 * @return none
	 */
	void add (const std::vector<std::complex<double>> alm, const double ww, const int bin=0);

	/**
	 * @brief compute add the \f$ a_{lm}\f$ for the normalized
	 * coordinates \f$ \lbrace x, y, z \rbrace \f$ to a specific 
	 * separation bin with a weight
	 *
	 * @param xx the x coordinate
	 *
	 * @param yy the y coordinate
	 *
	 * @param zz the z coordinate
	 *
	 * @param ww the weight
	 *
	 * @param bin the separation bin
	 *
	 * @return none
	 */
	void add (const double xx, const double yy, const double zz, const double ww, const int bin=0);

	/**
	 * @brief compute the product of the 
	 * \f$ a_{lm} \f$ in two separation bin.
	 *
	 * This function computes the product of the 
	 * \f$ a_{lm} \f$ in two separation bin:
	 *
	 * \f[
	 *   \zeta_l(r_1, r_2) = \sum_{m=-l}^l a_{lm}(r_1) a^*_{lm} (r_2)
	 * \f]
	 *
	 * @param l the coefficient order
	 *
	 * @param bin1 the first separation bin
	 *
	 * @param bin2 the second separation bin
	 *
	 * @return none
	 */
	double power (const int l, const int bin1, const int bin2);
    };

 }
}

#endif
