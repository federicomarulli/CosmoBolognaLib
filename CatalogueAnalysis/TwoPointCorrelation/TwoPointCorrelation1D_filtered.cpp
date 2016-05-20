/********************************************************************
 *  Copyright (C) 2015 by Federico Marulli and Alfonso Veropalumbo  *
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
 *  @file CatalogueAnalysis/TwoPointCorrelation/TwoPointCorrelation1D_filtered.cpp
 *
 *  @brief Methods of the class TwoPointCorrelation1D_filtered used to
 *  measure the filtered monopole of the two-point correlation function
 *
 *  This file contains the implementation of the methods of the class
 *  TwoPointCorrelation1D_filtered used to measure the filtered monopole of the
 *  two-point correlation function
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */


#include "TwoPointCorrelation1D_filtered.h"

using namespace cosmobl;
using namespace catalogue;
using namespace pairs;
using namespace twopt;


// ============================================================================================


void cosmobl::twopt::TwoPointCorrelation1D_filtered::set_parameters (const binType binType, const double rMin, const double rMax, const int nbins) 
{
  vector<double> rc, wc, error;
 
  rc = linear_bin_vector(nbins, rMin, rMax);
  wc.resize(nbins, 0);
  error.resize(nbins, 0);

  m_dataset->set_xx(rc);
  m_dataset->set_fx(wc);
  m_dataset->set_error_fx(error);
}


// ============================================================================================


void cosmobl::twopt::TwoPointCorrelation1D_filtered::set_parameters (const binType binType, const double rMin, const double rMax, const double binSize)
{
  int nbins = (rMax-rMin)/binSize;

  vector<double> rc, wc, error;

  rc = linear_bin_vector(nbins, rMin, rMax);
  wc.resize(nbins, 0);
  error.resize(nbins, 0);

  m_dataset->set_xx(rc);
  m_dataset->set_fx(wc);
  m_dataset->set_error_fx(error);
}


// ============================================================================================


void cosmobl::twopt::TwoPointCorrelation1D_filtered::measure (const ErrorType errType, const string dir_output_pairs, const vector<string> dir_input_pairs, const string dir_output_ResampleXi, const int nMocks, const int count_dd, const int count_rr, const int count_dr, const bool tcount)
{
  TwoPointCorrelation1D_monopole::measurePoisson(dir_output_pairs, dir_input_pairs, count_dd, count_rr, count_dr, tcount);

  double binSize  = 1./m_dd->binSize_inv_D1();

  vector<double> rad = TwoPointCorrelation1D_monopole::m_dataset->xx();
  vector<double> xi = TwoPointCorrelation1D_monopole::m_dataset->fx();
  vector<double> error_xi = TwoPointCorrelation1D_monopole::m_dataset->error_fx();

  vector<double> rc = m_dataset->xx();
  vector<double> wc = m_dataset->fx();
  vector<double> error_wc = m_dataset->error_fx();

  for (size_t i=0; i<rc.size(); i++) 
    for (size_t j=0; j<rad.size(); j++) {
      wc[i] += rad[j]*rad[j]*binSize*xi[j]*Filter(rad[j],rc[i]);
      error_wc[i] += pow(rad[j]*rad[j]*binSize*error_xi[j]*Filter(rad[j],rc[i]),2);
    }

  for_each( error_wc.begin(), error_wc.end(), [] (double &vv) { vv = sqrt(vv);} );

  m_dataset->set_xx(rc);
  m_dataset->set_fx(wc);
  m_dataset->set_error_fx(error_wc);
}


// ============================================================================================


void cosmobl::twopt::TwoPointCorrelation1D_filtered::read (const string dir, const string file) 
{
  m_dataset->read(dir+file);
}


// ============================================================================================


void cosmobl::twopt::TwoPointCorrelation1D_filtered::write (const string dir, const string file, const int rank) const 
{
  m_dataset->write(dir, file, "rc", "wc", rank);
}

