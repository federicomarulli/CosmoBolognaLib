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
 *  @file Data/TaperedCovarianceMatrix.cpp
 *
 *  @brief Methods of the class TaperedCovarianceMatrix
 *
 *  This file contains the implementation of the methods of the class
 *  TaperedCovarianceMatrix
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unibo.it, alfonso.veropalumbo@unibo.it
 */

#include "Data.h"
#include "EigenWrapper.h"
#include "CovarianceMatrix.h"
#include "TaperedCovarianceMatrix.h"

using namespace std;

using namespace cbl;
using namespace data;

// ======================================================================================


void cbl::data::TaperedCovarianceMatrix::m_set_tapering(const double tapering_factor)
{
  m_tapering_factor = tapering_factor;

  m_tapering_function.resize(m_order, m_order);

  for (int i=0; i<static_cast<int>(m_order); i++)
    for (int j=0; j<static_cast<int>(m_order); j++) {
      double xx = abs(i-j);
      if (xx<m_tapering_factor)
	m_tapering_function(i, j) = pow(1-xx/m_tapering_factor, 4)*(4*xx/m_tapering_factor+1);
      else
	m_tapering_function(i, j) = 0;
    }
}


// ======================================================================================


void cbl::data::TaperedCovarianceMatrix::m_set (const vector<double> covariance, const double nmeasures, const double prec)
{
  (void)prec;

  m_hartlap_factor = hartlap_factor(m_order, nmeasures);

  m_matrix = cbl::wrapper::eigen::SquareMatrixToEigen(covariance).cwiseProduct(m_tapering_function);

  m_precision = m_matrix.inverse().cwiseProduct(m_tapering_function);
  
  m_variance = m_matrix.diagonal();
  
  m_std = m_variance.cwiseSqrt();

  m_correlation = m_matrix.cwiseQuotient((m_std * m_std.transpose()));
}


// ======================================================================================


void cbl::data::TaperedCovarianceMatrix::set(const double tapering_factor, const cbl::data::CovarianceMatrix covariance)
{
  m_order = covariance.order();
  m_set_tapering(tapering_factor);
  m_set(flatten(covariance()));
}
