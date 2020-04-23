/********************************************************************
 *  Copyright (C) 2014 by Federico Marulli and Alfonso Veropalumbo  *
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
 *  @file Headers/TaperedCovarianceMatrix.h
 *
 *  @brief The class TaperedCovarianceMatrix
 *
 *  This file defines the interface of the class TaperedCovarianceMatrix
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unibo.it, alfonso.veropalumbo@unibo.it
 */

#ifndef __TCOVMAT__
#define __TCOVMAT__

#include "CovarianceMatrix.h"

namespace cbl {

  /**
   *  @brief The namespace of the functions and classes used to handle
   *  <B> data </B> of any kind
   *  
   *  The \e data namespace contains all the main functions and
   *  classes to handle data of any kind
   */
  namespace data {

     
    /**
     *  @class TaperedCovarianceMatrix TaperedCovarianceMatrix.h
     *  "Headers/TaperedCovarianceMatrix.h"
     *
     *  @brief The class TaperedCovarianceMatrix
     *
     *  This is the base class used to manage 
     *  tapered covariance matrices
     */
    class TaperedCovarianceMatrix : public CovarianceMatrix
    {
      protected:

	/// tapering factor
	double m_tapering_factor;

	/// tapering function
	Eigen::MatrixXd m_tapering_function;

	/**
	 * @brief privat member that sets the tapering matrix
	 *
	 * @param tapering_factor the tapering factor
	 *
	 * @return None
	 */
	void m_set_tapering(const double tapering_factor);

	/**
	 * @brief set internal attributes
	 *
	 * @param matrix the covariance matrix
	 *
	 * @param nmeasures number of measures
	 *
	 * @param prec the precision required in the inversion of the
	 * covariance matrix
	 *
	 * @return None
	 */
 	void m_set (const std::vector<double> matrix, const double nmeasures=-1, const double prec=1.e-10);

      public:

      /**
       *  @name Constructors/destructors
       */
      ///@{

      /**
       *  @brief default constructor
       *
       *  @return an object of class TaperedCovarianceMatrix
       */
      TaperedCovarianceMatrix () { m_set_default(); }

      /**
       *  @brief constructor which sets the
       *  covariance matrix
       *
       *  @param tapering_factor the tapering factor, in bin units
       *
       *  @param covariance object of type cbl::data::CovarianceMatrix
       *
       *  @return an object of class TaperedCovarianceMatrix
       */
      TaperedCovarianceMatrix (const double tapering_factor, const CovarianceMatrix covariance) 
      { set(tapering_factor, covariance); }

      /**
       *  @brief default destructor
       *  @return none
       */
      virtual ~TaperedCovarianceMatrix () = default;

      ///@}

      /**
       *  @name Member functions to set the private/protected members
       */
      ///@{ 


      /**
       * @brief Apply covariance tapering
       *
       * We smooth the covariance matrix using the tapering technique.
       * We follow procedure descrived in Paz et al. 2015
       * [https://arxiv.org/abs/1508.03162]
       *
       * This is done by doing the Hadamard product (i.e. the entry- wise product)
       * of the covariance matric \f$C\f$ and the tapering matrix \f$T\f$.
       *
       * The tapering matrix is defined as an isotropic covariance matrix 
       * by means of a taper function \f$K\f$:
       *
       * Following Kauffam et al. 2008, the tapering function is:
       * \f[
       *    K(x)=\left\{\begin{array}{ll}
       *    \left(1-\frac{x}{T_{\mathrm{p}}}\right)^{4}\left(4 \frac{x}{T_{\mathrm{p}}}+1\right) & \text { if } x<T_{\mathrm{p}} \\
       *    0 & \text { if } x \geq T_{\mathrm{p}}
       *    \end{array}\right.
       * \f]
       *
       * \f$T_p\f$ is the tapering parameter
       *
       * @param tapering_factor the tapering parameter
       *
       * @param covariance object of type cbl::data::CovarianceMatrix
       *
       * @return None
       */
      void set(const double tapering_factor, const CovarianceMatrix covariance);

      ///@}
    };
  }
}

#endif
