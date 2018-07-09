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
 *  CatalogueAnalysis/TwoPointCorrelation/TwoPointCorrelationCross1D.cpp
 *
 *  @brief Methods of the class TwoPointCorrelationCross1D used to
 *  measure the 1D cross of the two-point correlation function
 *
 *  This file contains the implementation of the methods of the class
 *  TwoPointCorrelationCross1D used to measure the 1D cross two-point
 *  correlation function
 *
 *  @authors Federico Marulli, Carlo Giocoli
 *
 *  @authors federico.marulli3@unbo.it, carlo.giocoli@unibo.it
 */


#include "TwoPointCorrelationCross1D.h"
#include "Data1D_extra.h"

using namespace std;

using namespace cbl;
using namespace catalogue;
using namespace chainmesh;
using namespace pairs;
using namespace measure;
using namespace twopt;


// ============================================================================


shared_ptr<data::Data> cbl::measure::twopt::TwoPointCorrelationCross1D::correlation_SzapudiSzalayEstimator (const shared_ptr<pairs::Pair> d1d2, const shared_ptr<pairs::Pair> rr, const shared_ptr<pairs::Pair> d1r, const shared_ptr<pairs::Pair> d2r, const int nData1, const double nData1_weighted, const int nData2, const double nData2_weighted, const int nRandom, const double nRandom_weighted)
{
  vector<double> rad(m_d1d2->nbins()), xi(m_d1d2->nbins(), -1.), error(m_d1d2->nbins(), 1000.);

  // number of objects in the first data catalogue
  int nD1 = (nData1>0) ? nData1 : m_data->nObjects();

  // number of objects in the second data catalogue
  int nD2 = (nData2>0) ? nData2 : m_data2->nObjects();
  
  // weighted number of objects in the first data catalogue
  double nD1w = (nData1_weighted>0) ? nData1_weighted : m_data->weightedN();

  // weighted number of objects in the second data catalogue
  double nD2w = (nData2_weighted>0) ? nData2_weighted : m_data2->weightedN();
  
  // number of objects in the random catalogue
  int nR = (nRandom>0) ? nRandom : m_random->nObjects();

  // weighted number of objects in the random catalogue
  double nRw = (nRandom_weighted>0) ? nRandom_weighted : m_random->weightedN();

  // inverse of the total number of data1-data2 pairs
  double nD1D2i = 1./(nD1w*nD2w);
  
  // inverse of the total number of random-random pairs
  double nRRi = 1./(nRw*m_random_dilution_fraction*(nRw*m_random_dilution_fraction-1.)*0.5);

  // inverse of the total number of data1-random pairs
  double nD1Ri = 1./(nD1w*nRw);

  // inverse of the total number of data2-random pairs
  double nD2Ri = 1./(nD2w*nRw);
  
  for (int i=0; i<d1d2->nbins(); i++) {
  
    rad[i] = d1d2->scale(i);

    if (d1d2->PP1D_weighted(i)>0) {

      if (rr->PP1D_weighted(i)<1.e-30) 
	ErrorCBL("Error in correlation_SzapudiSzalayEstimator() of TwoPointCorrelation1D.cpp: there are no random objects in the bin "+conv(i, par::fINT)+"; please, either increase the total number of random objects or enlarge the bin size!");

      // normalised number of data1-data2 weighted pairs
      double D1D2_norm = d1d2->PP1D_weighted(i)*nD1D2i;

      // normalised number of random-random weighted pairs
      double RR_norm = rr->PP1D_weighted(i)*nRRi;

      // normalised number of data1-random weighted pairs
      double D1R_norm = d1r->PP1D_weighted(i)*nD1Ri;

      // normalised number of data2-random weighted pairs
      double D2R_norm = d2r->PP1D_weighted(i)*nD2Ri;
      
      // Szapudi & Szalay estimator
      xi[i] = max(-1., (D1D2_norm-D1R_norm-D2R_norm+RR_norm)/RR_norm);

      // Poisson error
      error[i] = PoissonError(Estimator::_SzapudiSzalay_, d1d2->PP1D(i), rr->PP1D(i), d1r->PP1D(i), d2r->PP1D(i), nD1, nD2, nR);
      
    }
  }
  
  return (!m_compute_extra_info) ? move(unique_ptr<data::Data1D>(new data::Data1D(rad, xi, error))) : data_with_extra_info(d1d2, rad, xi, error);
}
