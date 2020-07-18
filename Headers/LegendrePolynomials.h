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
 *  @file Headers/LegendrePolynomials.h
 *
 *  @brief Class to manage Legendre polymials computation
 *
 *  This file contains the prototypes of the class LegendrePolynomials
 *
 *  @author Alfonso Veropalumbo
 *
 *  @author alfonso.veropalumbo@uniroma3.it
 */

#ifndef __LEGPOLS__
#define __LEGPOLS__

#include "Func.h"

namespace cbl {

 namespace glob {

    /**
     *  @class LegendrePolynomials LegendrePolynomials.h "Headers/LegendrePolynomials.h"
     *
     *  @brief The class LegendrePolynomials
     *
     *  This class is used to handle Legendre polynomials. 
     *  It contains all methods to compute Legendre polynomials as well
     *  as their integrals.
     *
     */
    class LegendrePolynomials {

      protected:

	/// Number of Legendre polynomial
	size_t m_nOrders;

	/// Coefficients of the Legendre polynomials
	Eigen::MatrixXd m_coefficients;

	/**
	 * @brief set internal attribute m_coefficients
	 *
	 * Set the coefficients of the Legendre Polynomials according
	 * to the following recursive relation:
	 *
	 * \f[ P_{n}(x)=2^{n} \sum_{k=0}^{n}
	 *	x^{k}\left(\begin{array}{l} n \\ k
	 *	\end{array}\right)\left(\begin{array}{c}
	 *	\frac{n+k-1}{2} \\ n \end{array}\right) \f]
	 *
	 * @param lMax maximum order of Legendre polynomials
	 *
	 * @return None
	 */
	void m_set_coefficients (const int lMax);
	
      public:

	/**
	 *  @name Constructors/destructors
	 */
	///@{

	/**
	 * @brief Default constructor of LegendrePolynomials
	 *
	 * @return Object of type LegendrePolynomials
	 */
	LegendrePolynomials () { set(0); }

	/**
	 * @brief Constructor of LegendrePolynomials
	 *
	 * @param lMax maximum order of Legendre polynomials
	 *
	 * @param safe True \f$\rightarrow \f$ check input range is
	 * among -1 a 1 False \f$ \rightarrow \f$ do not check input
	 *
	 * @return Object of type LegendrePolynomials
	 */
	LegendrePolynomials (const int lMax, const bool safe=false);

	/**
	 * @brief Default destructor
	 *
	 * @return None
	 */
	~LegendrePolynomials () {}

	///@}
	
	/**
	 * @brief set maximum order of expansion
	 *
	 * @param lMax maximum order of Legendre polynomials
	 *
	 * @return None
	 */
	void set (const int lMax);

	/**
	 * @brief evaluate the Legendre polynomial of order ell
	 * at x
	 *
	 * @param x the point to evaluate Legendre polynomial
	 * @param ell the order of the Legendre polynomial
	 *
	 * @return the Legendre polynomial of order ell
	 * at x
	 */
	double operator () (const double x, const int ell);

	/**
	 * @brief evaluate the Legendre polynomial up to lMax
	 * at x
	 *
	 * @param x the point to evaluate Legendre polynomials
	 *
	 * @return vector containing the Legendre polynomials
	 * at x
	 */
	std::vector<double> operator () (const double x);

	/**
	 * @brief evaluate the bin-averaged Legendre polynomial of order ell
	 *
	 * @param x_min the lower bin edge
	 * @param x_max the up bin edge
	 * @param ell the order of the Legendre polynomial
	 *
	 * @return the bin-averaged Legendre polynomial of order ell
	 */
	double integral (const double x_min, const double x_max, const int ell);

	/**
	 * @brief evaluate the bin-averaged Legendre polynomials
	 *
	 * @param x_min the lower bin edge
	 * @param x_max the up bin edge
	 *
	 * @return the bin-averaged Legendre polynomials
	 */
	std::vector<double> integral (const double x_min, const double x_max);

	/**
	 * @brief evaluate the Legendre polynomials for
	 * triangle angles.
	 *
	 * @param r12 first triangle side
	 * @param r13 second triangle side
	 * @param r23 third triangle side
	 *
	 * @return the Legendre polynomials averaged over a triangle
	 */
	std::vector<std::vector<double>> triangle (const double r12, const double r13, const double r23);

	/**
	 * @brief evaluate the bin-averaged Legendre polynomials
	 * over a triangle. Triangle side can vary from a minimum to
	 * a maximum value
	 *
	 * @param r12_min the first side lower edge
	 * @param r12_max the first side upper edge
	 * @param r13_min the second side lower edge
	 * @param r13_max the second side upper edge
	 * @param r23_min the third side lower edge
	 * @param r23_max the third side upper edge
	 * @param rel_err the relative integration error
	 * @param nevals maximum number of function evaluation
	 *
	 * @return the Legendre polynomials averaged over a triangle
	 */
	std::vector<double> triangle_integral (const double r12_min, const double r12_max, const double r13_min, const double r13_max, const double r23_min, const double r23_max, const double rel_err=1.e-4, const int nevals=10000);
    };
 }
}

#endif
