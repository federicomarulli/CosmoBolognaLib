/*******************************************************************
 *  Copyright (C) 2010 by Federico Marulli and Alfonso Veropalumbo *
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
 *  @file Headers/FuncGrid_Bspline.h
 *
 *  @brief Class used to handle functions interpolated using a 
 *  basis spline http://mathworld.wolfram.com/B-Spline.html
 *
 *  @author Federico Marulli, Alfonso Veropalumbo
 *
 *  @author federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#ifndef __FUNCGRIDBSPL__
#define __FUNCGRIDBSPL__ 

#include "GSLwrapper.h"

// =====================================================================================


namespace cbl {

  namespace glob {
    
    /**
     *  @class FuncGrid FuncGrid_Bspline.h "Headers/FuncGrid_Bspline.h"
     *
     *  @brief The class FuncGrid_Bspline
     *
     *  This class is used to handle functions interpolated using a 
     *  basis spline (http://mathworld.wolfram.com/B-Spline.html)
     *  In particular this class wraps the GSL Bspline implementation.
     *  For further details, please refer to
     *  https://www.gnu.org/software/gsl/doc/html/bspline.html
     */
    class FuncGrid_Bspline
    {
      private:

        /// x values
	std::vector<double> m_x;

        /// fx values
	std::vector<double> m_fx;

	/// number of breakpoints
	int m_nbreakpoints;

	/// number of coefficients;
	int m_ncoefficients;

	/// basis spline order
	int m_order;

	/// pointer to GSL bspline workspace 
  	std::shared_ptr<gsl_bspline_workspace> m_bspline;

	/// pointer to GSL vector containing the basis spline coefficients
  	std::shared_ptr<gsl_vector> m_Bcoeff; 
  	
	/// pointer to GSL vector containing the basis spline weights
	std::shared_ptr<gsl_vector> m_weights; 
	
	/// pointer to GSL matrix containing the linear fit covariance
	std::shared_ptr<gsl_matrix> m_covariance;

	/// Integral of the function over the m_x range
	double m_integral;

      /**
       *  @name Functions to set internal variables
       */
      ///@{
	/**
	 * @brief set internal variables 
	 * for gsl bspline wrapping
	 *
	 * Function that sets internal variables 
	 * for gsl bspline wrapping
	 *
	 * @param x vector containing the x values
	 * @param fx vector containing the f(x) values
	 * @param nbreakpoints number of breakpoints
	 * @param order the basis spline order.
	 *
	 * @return None
	 */
	void m_set_bspline(const std::vector<double> x, const std::vector<double> fx, const int nbreakpoints, const int order);

	/**
	 * @brief set the b-spline knots
	 *
	 * set the knots uniformely spaced on a grid
	 *
	 * @param xmin minimum breakpoint
	 * @param xmax maximum breakpoint
	 *
	 * @return None
	 */
	void m_set_knots(const double xmin=cbl::par::defaultDouble, const double xmax=cbl::par::defaultDouble);

	/**
	 * @brief set the b-spline knots
	 *
	 * set the knots from a vector of breakpoints
	 * knots are the 
	 *
	 * @param breakpoints the breakpoints
	 *
	 * @return None
	 */
	void m_set_knots(const std::vector<double> breakpoints);

	/**
	 * @brief compute basis spline coefficients
	 * via a linear fit
	 *
	 * The fit uses an error: frac*f(x)
	 * with frac given in input. Default value is 0.1
	 *
	 * @param frac fraction of the fx to use as fit error.
	 * Default value is 0.1
	 *
	 * @return None
	 */
	void m_linear_fit(const double frac=0.1);

	/**
	 * @brief compute the integral of the function
	 * over the range [Min(m_x)-Max(m_x)] and set internal
	 * variable m_integral
	 * This is used to renormalize the output
	 *
	 * @return None
	 */
	void m_compute_func_integral();
	
      ///@}

      public:
      /**
       *  @name Constructors/destructors
       */
      ///@{

	/**
	 * @brief default constructor
	 * @return Object of type FuncGrid_Bspline
	 */
	FuncGrid_Bspline () {}

	/**
	 * @brief constructor
	 *
	 * @param x vector containing the x values
	 * @param fx vector containing the f(x) values
	 * @param nbreakpoints number of breakpoints
	 * @param order the basis spline order.
	 * Default is 4, leading to cubic basis spline
	 * @param frac fraction of the fx to use as fit error.
	 * Default is 0.1
	 * @param xmin minimum breakpoint
	 * @param xmax maximum breakpoint
	 *
	 * @return Object of type FuncGrid_Bspline
	 */
	FuncGrid_Bspline (const std::vector<double> x, const std::vector<double> fx, const int nbreakpoints, const int order=4, const double frac=0.1, const double xmin=cbl::par::defaultDouble, const double xmax=cbl::par::defaultDouble);

	/**
	 * @brief constructor
	 *
	 * @param x vector containing the x values
	 * @param fx vector containing the f(x) values
	 * @param breakpoints vector containing the breakpoints
	 * @param order the basis spline order.
	 * Default is 4, leading to cubic basis spline
	 * @param frac fraction of the fx to use as fit error.
	 * Default is 0.1
	 *
	 * @return Object of type FuncGrid_Bspline
	 */
	FuncGrid_Bspline (const std::vector<double> x, const std::vector<double> fx, const std::vector<double> breakpoints, const int order=4, const double frac=0.1);

	/**
	 *  @brief default destructor
	 *  @return none
	 */
	~FuncGrid_Bspline() = default;

      ///@}

      /**
       *  @name Functions to set the basis spline
       */
      ///@{

	/**
	 * @brief set internal members
	 *
	 * @param x vector containing the x values
	 * @param fx vector containing the f(x) values
	 * @param nbreakpoints number of breakpoints
	 * @param order the basis spline order
	 * Default is 4, leading to cubic basis spline
	 * @param frac fraction of the fx to use as fit error.
	 * Default is 0.1
	 * @param xmin minimum breakpoint
	 * @param xmax maximum breakpoint
	 *
	 * @return None
	 */
	void set (const std::vector<double> x, const std::vector<double> fx, const int nbreakpoints, const int order=4, const double frac=0.1, const double xmin=cbl::par::defaultDouble, const double xmax=cbl::par::defaultDouble);

        /**
	 * @brief set internal members
	 *
	 * @param x vector containing the x values
	 * @param fx vector containing the f(x) values
	 * @param breakpoints vector containing the breakpoints
	 * @param order the basis spline order
	 * Default is 4, leading to cubic basis spline
	 * @param frac fraction of the fx to use as fit error.
	 * Default is 0.1
	 *
	 * @return None
	 */
	void set (const std::vector<double> x, const std::vector<double> fx, const std::vector<double> breakpoints, const int order=4, const double frac=0.1);
      	
      ///@}

      /**
       *  @name Functions to get interpolated value
       */
      ///@{

      	/**
       	 *  @brief overloading of the () operator
         *  @param xx the value at which the function will be
         *  evaluated
	 *  @param integral output value of the function integral
         *  @return the function evaluated at xx
         */   
      	double operator () (const double xx, const double integral=cbl::par::defaultDouble) const;

        /**
         *  @brief evaluate the function at the xx points
         *  @param xx the values at which the function will be
         *  evaluated
	 *  @param integral output value of the function integral
         *  @return the function evaluated at the xx points
         */   
        std::vector<double> eval_func (const std::vector<double> xx, const double integral=cbl::par::defaultDouble) const;
      	
      ///@}
    };

  }
}
#endif
