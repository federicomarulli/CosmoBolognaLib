/********************************************************************
 *  Copyright (C) 2017 by Federico Marulli and Carlo Giocoli        *
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
 *  @file
 *  CatalogueAnalysis/TwoPointCorrelation/TwoPointCorrelationCross1D_monopole.cpp
 *
 *  @brief Methods of the class TwoPointCorrelationCross1D_monopole
 *  used to measure the monopole of the cross two-point correlation
 *  function
 *
 *  This file contains the implementation of the methods of the class
 *  TwoPointCorrelationCross1D_monopole used to measure the monopole
 *  of the cross two-point correlation function
 *
 *  @authors Federico Marulli, Carlo Giocoli
 *
 *  @authors federico.marulli3@unbo.it, carlo.giocoli@unibo.it
 */


#include "TwoPointCorrelationCross1D_monopole.h"

using namespace std;

using namespace cbl;
using namespace catalogue;
using namespace chainmesh;
using namespace data;
using namespace pairs;
using namespace measure;
using namespace twopt;


// ============================================================================================


void cbl::measure::twopt::TwoPointCorrelationCross1D_monopole::measure (const ErrorType errorType, const string dir_output_pairs, const vector<string> dir_input_pairs, const string dir_output_resample, const int nMocks, const bool count_d1d2, const bool count_rr, const bool count_d1r, const bool count_d2r, const bool tcount, const Estimator estimator)
{
  (void)dir_output_resample, (void)nMocks;
  
  switch (errorType) {
  case (ErrorType::_Poisson_) :
    measurePoisson(dir_output_pairs, dir_input_pairs, count_d1d2, count_rr, count_d1r, count_d2r, tcount, estimator);
    break;
  default:
    ErrorCBL("Error in measure() of TwoPointCorrelationCross1D_monopole.cpp, unknown type of error");
  }
}


// ============================================================================================


void cbl::measure::twopt::TwoPointCorrelationCross1D_monopole::measurePoisson (const string dir_output_pairs, const vector<string> dir_input_pairs, const bool count_d1d2, const bool count_rr, const bool count_d1r, const bool count_d2r, const bool tcount, const Estimator estimator)
{
  // ----------- count the data1-data2, random-random, data1-random and data2-random pairs, or read them from file ----------- 
  
  count_allPairs(m_twoPType, dir_output_pairs, dir_input_pairs, count_d1d2, count_rr, count_d1r, count_d2r, tcount, estimator);
  
  
  // ----------- compute the monopole of the two-point cross correlation function ----------- 

  if (estimator==Estimator::_SzapudiSzalay_)
    m_dataset = correlation_SzapudiSzalayEstimator(m_d1d2, m_rr, m_d1r, m_d2r);
  else
    ErrorCBL("Error in measurePoisson() of TwoPointCorrelationCross1D_monopole.cpp: the chosen estimator is not implemented!");
  
}
