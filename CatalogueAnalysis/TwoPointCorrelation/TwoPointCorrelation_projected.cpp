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
 *  @file CatalogueAnalysis/TwoPointCorrelation/TwoPointCorrelation_projected.cpp
 *
 *  @brief Methods of the class TwoPointCorrelation_projected used to
 *  measure the projected two-point correlation function
 *
 *  This file contains the implementation of the methods of the class
 *  TwoPointCorrelation_projected used to measure the projected
 *  two-point correlation function
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */


#include "TwoPointCorrelation_projected.h"
#include "Data1D.h"

using namespace cosmobl;
using namespace catalogue;
using namespace pairs;
using namespace twopt;


// ============================================================================================


shared_ptr<Data> cosmobl::twopt::TwoPointCorrelation_projected::ProjectedTwoP (const vector<double> rp, const vector<double> pi, const vector<vector<double> > xi, const vector<vector<double> > error_xi)
{
  vector<double> ww, error;

  ww.resize(0); ww.resize(rp.size(), 0.);
  error.resize(0); error.resize(rp.size(), 0.);
 
  double binSize = pi[1]-pi[0];

  int pim = nint((m_piMax_integral-Min(pi))/binSize); // to convert from Mpc/h into the vector index

  for (size_t i=0; i<rp.size(); i++) {

    ww[i] = 0.;
    error[i] = 0.;
    
    for (int j=0; j<pim; j++) {  
      if (xi[i][j]>-1.e29) { // check!!!
	ww[i] = ww[i]+2.*binSize*xi[i][j];
	error[i] += pow(2.*binSize*error_xi[i][j], 2); // check!!!!
      }
    }
  }

  for_each( error.begin(), error.end(), [] (double &vv) { vv = sqrt(vv);} );

  return move(unique_ptr<Data1D>(new Data1D(rp, ww, error)));

}


// ============================================================================================


void cosmobl::twopt::TwoPointCorrelation_projected::read(const string dir, const string file) 
{
  m_dataset->read(dir+file);
}


// ============================================================================================


void cosmobl::twopt::TwoPointCorrelation_projected::write (const string dir, const string file, const int rank) const 
{
  checkDim(m_dataset->xx(), m_dd->nbins_D1(), "rp");
  
  m_dataset->write(dir, file, "rp", "xi_projected", rank);
}


// ============================================================================================


void cosmobl::twopt::TwoPointCorrelation_projected::measure (const ErrorType errType, const string dir_output_pairs, const vector<string> dir_input_pairs, const string dir_output_ResampleXi, const int nMocks, const int count_dd, const int count_rr, const int count_dr, const bool tcount)
{
  switch (errType) {
  case (ErrorType::_Poisson_) :
    measurePoisson(dir_output_pairs, dir_input_pairs, count_dd, count_rr, count_dr, tcount);
    break;
  case (ErrorType::_Jackknife_) :
    measureJackknife(dir_output_pairs, dir_input_pairs, dir_output_ResampleXi, count_dd, count_rr, count_dr, tcount);
    break;
  case (ErrorType::_Bootstrap_) :
    measureBootstrap(nMocks, dir_output_pairs, dir_input_pairs, dir_output_ResampleXi, count_dd, count_rr, count_dr, tcount);
    break;
  }
}


// ============================================================================================


void cosmobl::twopt::TwoPointCorrelation_projected::measurePoisson (const string dir_output_pairs, const vector<string> dir_input_pairs, const int count_dd, const int count_rr, const int count_dr, const bool tcount)
{
  // ----------- measure the 2D two-point correlation function, xi(rp,pi) ----------- 

  TwoPointCorrelation2D_cartesian::measurePoisson(dir_output_pairs, dir_input_pairs, count_dd, count_rr, count_dr, tcount);

  
  // ----------- integrate the 2D two-point correlation function along the parallel direction ----------- 
  
  m_dataset = ProjectedTwoP(TwoPointCorrelation2D_cartesian::xx(), TwoPointCorrelation2D_cartesian::yy(), TwoPointCorrelation2D_cartesian::xi2D(), TwoPointCorrelation2D_cartesian::error2D());
}


// ============================================================================================


void cosmobl::twopt::TwoPointCorrelation_projected::measureJackknife (const string dir_output_pairs, const vector<string> dir_input_pairs, const string dir_output_ResampleXi, const int count_dd, const int count_rr, const int count_dr, const bool tcount)
{
  if (dir_output_ResampleXi != par::defaultString && dir_output_ResampleXi != "") {
    string mkdir = "mkdir -p "+dir_output_ResampleXi;
    if (system(mkdir.c_str())) {}
  }

  vector<shared_ptr<Data> > data;
  vector<shared_ptr<pairs::Pair> > dd_regions, rr_regions, dr_regions;
  count_allPairs_region (dd_regions, rr_regions, dr_regions, TwoPType::_2D_Cartesian_, dir_output_pairs,dir_input_pairs, count_dd, count_rr, count_dr,  tcount);

  auto data_cart = (count_dr>-1) ? LandySzalayEstimatorTwoP(m_dd, m_rr, m_dr, m_data->weightedN(), m_random->weightedN()) : NaturalEstimatorTwoP(m_dd, m_rr, m_data->weightedN(), m_random->weightedN());

  if (count_dr>-1) {
    data = XiJackknife(dd_regions, rr_regions, dr_regions);
  }
  else {
    data = XiJackknife(dd_regions, rr_regions);
  }

  vector<vector<double> > ww, covariance;
  for (size_t i=0; i<data.size(); i++) {
    ww.push_back(data[i]->fx());

    if (dir_output_ResampleXi != par::defaultString && dir_output_ResampleXi != "") {
      string file = "xi_projected_Jackkknife_"+conv(i, par::fINT)+".dat";
      data[i]->write(dir_output_ResampleXi, file, "rp", "xi_projected", 0);
    }
  }

  covariance_matrix(ww,covariance, 1);

  m_dataset = ProjectedTwoP(data_cart->xx(), data_cart->yy(), data_cart->fxy(), data_cart->error_fxy());
  m_dataset->set_covariance_fx(covariance);

}


// ============================================================================================


void cosmobl::twopt::TwoPointCorrelation_projected::measureBootstrap (const int nMocks, const string dir_output_pairs, const vector<string> dir_input_pairs, const string dir_output_ResampleXi, const int count_dd, const int count_rr, const int count_dr, const bool tcount)
{
  if (dir_output_ResampleXi != par::defaultString && dir_output_ResampleXi != "") {
    string mkdir = "mkdir -p "+dir_output_ResampleXi;
    if(system(mkdir.c_str())){}
  }
  
  vector<shared_ptr<Data> > data;
  vector<shared_ptr<pairs::Pair> > dd_regions, rr_regions, dr_regions;
  count_allPairs_region (dd_regions, rr_regions, dr_regions, TwoPType::_2D_Cartesian_, dir_output_pairs,dir_input_pairs, count_dd, count_rr, count_dr,  tcount);

  auto data_cart = (count_dr>-1) ? LandySzalayEstimatorTwoP(m_dd, m_rr, m_dr, m_data->weightedN(), m_random->weightedN()) : NaturalEstimatorTwoP(m_dd, m_rr, m_data->weightedN(), m_random->weightedN());

  if (count_dr>-1) 
    data = XiBootstrap(nMocks, dd_regions, rr_regions, dr_regions);
  else
    data = XiBootstrap(nMocks, dd_regions, rr_regions);

  vector<vector<double> > ww, covariance;
  for (size_t i=0; i<data.size(); i++) {
    ww.push_back(data[i]->fx());
    if (dir_output_ResampleXi != par::defaultString){
      string file = "xi_projected_Bootstrap_"+conv(i, par::fINT)+".dat";
      data[i]->write(dir_output_ResampleXi, file, "rp", "xi_projected", 0);
    }
  }
  covariance_matrix(ww, covariance, 0);

  m_dataset = ProjectedTwoP(data_cart->xx(), data_cart->yy(), data_cart->fxy(), data_cart->error_fxy());
  m_dataset->set_covariance_fx(covariance);
}


// ============================================================================================


vector<shared_ptr<Data> > cosmobl::twopt::TwoPointCorrelation_projected::XiJackknife (const vector<shared_ptr<pairs::Pair> > dd, const vector<shared_ptr<pairs::Pair> > rr)
{
  vector<shared_ptr<Data> > data;
  
  auto data2d = TwoPointCorrelation2D_cartesian::XiJackknife(dd, rr);

  for (size_t i=0; i<data2d.size(); i++)
    data.push_back(move(ProjectedTwoP(data2d[i]->xx(), data2d[i]->yy(), data2d[i]->fxy(), data2d[i]->error_fxy())));

  return data;
}


// ============================================================================================


vector<shared_ptr<Data> > cosmobl::twopt::TwoPointCorrelation_projected::XiJackknife (const vector<shared_ptr<pairs::Pair> > dd, const vector<shared_ptr<pairs::Pair> > rr, const vector<shared_ptr<pairs::Pair> > dr)
{
  vector<shared_ptr<Data> > data;

  auto data2d = TwoPointCorrelation2D_cartesian::XiJackknife(dd, rr, dr);

  for (size_t i=0; i<data2d.size(); i++)
    data.push_back(move(ProjectedTwoP(data2d[i]->xx(), data2d[i]->yy(), data2d[i]->fxy(), data2d[i]->error_fxy())));
  
  return data;
}


// ============================================================================================


vector<shared_ptr<Data> > cosmobl::twopt::TwoPointCorrelation_projected::XiBootstrap (const int nMocks, const vector<shared_ptr<pairs::Pair> > dd, const vector<shared_ptr<pairs::Pair> > rr)
{
  vector<shared_ptr<Data> > data;

  auto data2d = TwoPointCorrelation2D_cartesian::XiBootstrap(nMocks, dd, rr);

  for (size_t i=0; i<data2d.size(); i++)
    data.push_back(move(ProjectedTwoP(data2d[i]->xx(), data2d[i]->yy(), data2d[i]->fxy(), data2d[i]->error_fxy())));

  return data;
}


// ============================================================================================


vector<shared_ptr<Data> > cosmobl::twopt::TwoPointCorrelation_projected::XiBootstrap (const int nMocks, const vector<shared_ptr<pairs::Pair> > dd, const vector<shared_ptr<pairs::Pair> > rr, const vector<shared_ptr<pairs::Pair> > dr)
{
  vector<shared_ptr<Data> > data;

  auto data2d = TwoPointCorrelation2D_cartesian::XiBootstrap(nMocks, dd, rr, dr);
  
  for (size_t i=0; i<data2d.size(); i++)
    data.push_back(move(ProjectedTwoP(data2d[i]->xx(), data2d[i]->yy(), data2d[i]->fxy(), data2d[i]->error_fxy())));
  
  return data;
}
