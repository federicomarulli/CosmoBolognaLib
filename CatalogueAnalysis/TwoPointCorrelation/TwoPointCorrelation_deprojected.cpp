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
 *  @file
 *  CatalogueAnalysis/TwoPointCorrelation/TwoPointCorrelation_deprojected.cpp
 *
 *  @brief Methods of the class TwoPointCorrelation_deprojected used
 *  to measure the projected two-point correlation function
 *
 *  This file contains the implementation of the methods of the class
 *  TwoPointCorrelation_deprojected used to measure the projected
 *  two-point correlation function
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */


#include "TwoPointCorrelation_deprojected.h"
#include "Data1D_extra.h"

using namespace cosmobl;
using namespace catalogue;
using namespace chainmesh;
using namespace pairs;
using namespace twopt;


// ============================================================================================


shared_ptr<data::Data> cosmobl::twopt::TwoPointCorrelation_deprojected::Deprojected (const vector<double> rp, const vector<double> ww, const vector<double> error_ww)
{

  vector<double> rad(rp.size(),0), xi(rp.size(),0), error(rp.size(),0);

  for (size_t i=0; i<rp.size(); i++) {
    
    double ri = rp[i];
    rad[i] = ri;

    for (size_t j=i; j<rp.size()-1; j++) {
      double rj = rp[j];
      double rj1 = rp[j+1];
      double fact = 1./(rj1-rj)*log((rj1+sqrt(max(0., rj1*rj1-ri*ri)))/(rj+sqrt(max(0., rj*rj-ri*ri))));

      xi[i] -= (ww[j+1]-ww[j])*fact;
      error[i] += pow(error_ww[j+1]*fact, 2)+pow(error_ww[j]*fact, 2);
    }
    
  }

  for_each( xi.begin(), xi.end(), [] (double &vv) { vv /= par::pi;} );
  for_each( error.begin(), error.end(), [] (double &vv) { vv = sqrt(vv/par::pi);} );

  return (!m_compute_extra_info) ? move(unique_ptr<data::Data1D>(new data::Data1D(rad, xi, error))) : data_with_extra_info(rad, xi, error);

}


// ============================================================================================


void cosmobl::twopt::TwoPointCorrelation_deprojected::measure (const ErrorType errorType, const string dir_output_pairs, const vector<string> dir_input_pairs, const string dir_output_resample, const int nMocks, const bool count_dd, const bool count_rr, const bool count_dr, const bool tcount, const Estimator estimator)
{
  switch (errorType) {
  case (ErrorType::_Poisson_) :
    measurePoisson(dir_output_pairs, dir_input_pairs, count_dd, count_rr, count_dr, tcount, estimator);
    break;
  case (ErrorType::_Jackknife_) :
    measureJackknife(dir_output_pairs, dir_input_pairs, dir_output_resample, count_dd, count_rr, count_dr, tcount, estimator);
    break;
  case (ErrorType::_Bootstrap_) :
    measureBootstrap(nMocks, dir_output_pairs, dir_input_pairs, dir_output_resample, count_dd, count_rr, count_dr, tcount, estimator);
    break;
  }
}


// ============================================================================================


void cosmobl::twopt::TwoPointCorrelation_deprojected::measurePoisson (const string dir_output_pairs, const vector<string> dir_input_pairs, const bool count_dd, const bool count_rr, const bool count_dr, const bool tcount, const Estimator estimator)
{
  // ----------- measure the 2D two-point correlation function, xi(rp,pi) ----------- 

  TwoPointCorrelation_projected::measurePoisson(dir_output_pairs, dir_input_pairs, count_dd, count_rr, count_dr, tcount, estimator);

  
  // ----------- integrate the 2D two-point correlation function along the parallel direction ----------- 
  
  m_dataset = Deprojected(TwoPointCorrelation_projected::dataset()->xx(), TwoPointCorrelation_projected::dataset()->fx(), TwoPointCorrelation_projected::dataset()->error_fx());
}


// ============================================================================================


void cosmobl::twopt::TwoPointCorrelation_deprojected::measureJackknife (const string dir_output_pairs, const vector<string> dir_input_pairs, const string dir_output_resample, const bool count_dd, const bool count_rr, const bool count_dr, const bool tcount, const Estimator estimator)
{
  if (dir_output_resample != par::defaultString && dir_output_resample != "") {
    string mkdir = "mkdir -p "+dir_output_resample;
    if(system(mkdir.c_str())){}
  }

  vector<shared_ptr<data::Data>> data;
  vector<shared_ptr<pairs::Pair>> dd_regions, rr_regions, dr_regions;
  count_allPairs_region(dd_regions, rr_regions, dr_regions, TwoPType::_2D_Cartesian_, dir_output_pairs, dir_input_pairs, count_dd, count_rr, count_dr, tcount, estimator);

  auto data_cart = (estimator==_natural_) ? NaturalEstimator(m_dd, m_rr, m_data->weightedN(), m_random->weightedN()) : LandySzalayEstimator(m_dd, m_rr, m_dr, m_data->weightedN(), m_random->weightedN());

  auto data_proj = TwoPointCorrelation_projected::Projected(data_cart->xx(), data_cart->yy(), data_cart->fxy(), data_cart->error_fxy());

  if (estimator==_natural_) 
    data = XiJackknife(dd_regions, rr_regions);
  else if (estimator==_LandySzalay_)
    data = XiJackknife(dd_regions, rr_regions, dr_regions);
  else
    ErrorCBL("Error in measureJackknife() of TwoPointCorrelation_deprojected.cpp: the chosen estimator is not implemented!");
  
  vector<vector<double>> ww, covariance;
  for (size_t i=0; i<data.size(); i++) {
    ww.push_back(data[i]->fx());

    if (dir_output_resample != par::defaultString && dir_output_resample != "") {
      string file = "xi_deprojected_Jackknife_"+conv(i, par::fINT)+".dat";
      data[i]->write(dir_output_resample, file, "[1] separation at the bin centre # [2] deprojected two-point correlation function # [3] error", 0);
    }
  }

  covariance_matrix(ww, covariance, true);

  m_dataset = Deprojected(data_proj->xx(), data_proj->fx(), data_proj->error_fx());
  m_dataset->set_covariance(covariance);

}


// ============================================================================================


void cosmobl::twopt::TwoPointCorrelation_deprojected::measureBootstrap (const int nMocks, const string dir_output_pairs, const vector<string> dir_input_pairs, const string dir_output_resample, const bool count_dd, const bool count_rr, const bool count_dr, const bool tcount, const Estimator estimator)
{
  if (dir_output_resample != par::defaultString && dir_output_resample != "") {
    string mkdir = "mkdir -p "+dir_output_resample;
    if (system(mkdir.c_str())) {}
  }
  
  vector<shared_ptr<data::Data>> data;
  vector<shared_ptr<pairs::Pair>> dd_regions, rr_regions, dr_regions;
  count_allPairs_region(dd_regions, rr_regions, dr_regions, TwoPType::_2D_Cartesian_, dir_output_pairs, dir_input_pairs, count_dd, count_rr, count_dr, tcount, estimator);

  auto data_cart = (estimator==_natural_) ? NaturalEstimator(m_dd, m_rr, m_data->weightedN(), m_random->weightedN()) : LandySzalayEstimator(m_dd, m_rr, m_dr, m_data->weightedN(), m_random->weightedN());

  auto data_proj = TwoPointCorrelation_projected::Projected(data_cart->xx(), data_cart->yy(), data_cart->fxy(), data_cart->error_fxy());

  if (estimator==_natural_) 
    data = XiBootstrap(nMocks, dd_regions, rr_regions);
  else if (estimator==_LandySzalay_)
    data = XiBootstrap(nMocks, dd_regions, rr_regions, dr_regions);
  else
    ErrorCBL("Error in measureBootstrap() of TwoPointCorrelation_deprojected.cpp: the chosen estimator is not implemented!");
  
  vector<vector<double>> ww, covariance;
  for (size_t i=0; i<data.size(); i++) {
    ww.push_back(data[i]->fx());

    if (dir_output_resample != par::defaultString && dir_output_resample != "") {
      string file = "xi_deprojected_Bootstrap_"+conv(i, par::fINT)+".dat";
      data[i]->write(dir_output_resample, file, "[1] separation at the bin centre # [2] deprojected two-point correlation function # [3] error", 0);
    }
  }
  
  covariance_matrix(ww, covariance, false);

  m_dataset = Deprojected(data_proj->xx(), data_proj->fx(), data_proj->error_fx());
  m_dataset->set_covariance(covariance);
}


// ============================================================================================


vector<shared_ptr<data::Data>> cosmobl::twopt::TwoPointCorrelation_deprojected::XiJackknife (const vector<shared_ptr<pairs::Pair>> dd, const vector<shared_ptr<pairs::Pair>> rr)
{
  vector<shared_ptr<data::Data>> data;
  
  auto data_proj = TwoPointCorrelation_projected::XiJackknife(dd, rr);

  for (size_t i=0; i<data_proj.size(); i++)
    data.push_back(move(Deprojected(data_proj[i]->xx(), data_proj[i]->fx(), data_proj[i]->error_fx())));
  
  return data;
}


// ============================================================================================


vector<shared_ptr<data::Data>> cosmobl::twopt::TwoPointCorrelation_deprojected::XiJackknife (const vector<shared_ptr<pairs::Pair>> dd, const vector<shared_ptr<pairs::Pair>> rr, const vector<shared_ptr<pairs::Pair>> dr)
{
  vector<shared_ptr<data::Data>> data;
  
  auto data_proj = TwoPointCorrelation_projected::XiJackknife(dd, rr, dr);
  
  for (size_t i=0; i<data_proj.size(); i++)
    data.push_back(move(Deprojected(data_proj[i]->xx(), data_proj[i]->fx(), data_proj[i]->error_fx())));
  
  return data;
}


// ============================================================================================


vector<shared_ptr<data::Data>> cosmobl::twopt::TwoPointCorrelation_deprojected::XiBootstrap (const int nMocks, const vector<shared_ptr<pairs::Pair>> dd, const vector<shared_ptr<pairs::Pair>> rr)
{
  vector<shared_ptr<data::Data>> data;

  auto data_proj = TwoPointCorrelation_projected::XiBootstrap(nMocks, dd, rr);

  for (size_t i=0; i<data_proj.size(); i++)
    data.push_back(move(Deprojected(data_proj[i]->xx(), data_proj[i]->fx(), data_proj[i]->error_fx())));
  
  return data;
}


// ============================================================================================


vector<shared_ptr<data::Data>> cosmobl::twopt::TwoPointCorrelation_deprojected::XiBootstrap (const int nMocks, const vector<shared_ptr<pairs::Pair>> dd, const vector<shared_ptr<pairs::Pair>> rr, const vector<shared_ptr<pairs::Pair>> dr)
{
  vector<shared_ptr<data::Data>> data;
  
  auto data_proj = TwoPointCorrelation_projected::XiBootstrap(nMocks, dd, rr, dr);
  
  for (size_t i=0; i<data_proj.size(); i++)
    data.push_back(move(Deprojected(data_proj[i]->xx(), data_proj[i]->fx(), data_proj[i]->error_fx())));
  
  return data;
}


// ============================================================================================


void cosmobl::twopt::TwoPointCorrelation_deprojected::read (const string dir, const string file) 
{
  m_dataset->read(dir+file);
}


// ============================================================================================


void cosmobl::twopt::TwoPointCorrelation_deprojected::write (const string dir, const string file, const int rank) const 
{
  checkDim(m_dataset->xx(), m_dd->nbins_D1(), "rad");

  string header = "[1] separation at the bin centre # [2] deprojected two-point correlation function # [3] error";
  if (m_compute_extra_info) header += " # [4] mean separation # [5] standard deviation of the separation distribution # [6] mean redshift # [7] standard deviation of the redshift distribution";
  
  m_dataset->write(dir, file, header, rank);
}


// ============================================================================


void cosmobl::twopt::TwoPointCorrelation_deprojected::read_covariance_matrix (const string dir, const string file)
{
  m_dataset->set_covariance(dir+file);
}


// ============================================================================


void cosmobl::twopt::TwoPointCorrelation_deprojected::write_covariance_matrix (const string dir, const string file) const
{
  m_dataset->write_covariance(dir, file, "r");
}


// ============================================================================


void cosmobl::twopt::TwoPointCorrelation_deprojected::compute_covariance_matrix (const vector<shared_ptr<data::Data>> xi_collection, const bool doJK)
{
  vector<vector<double>> xi;

  for(size_t i=0;i<xi_collection.size();i++)
    xi.push_back(xi_collection[i]->fx());

  vector<vector<double>> cov_mat;
  cosmobl::covariance_matrix(xi, cov_mat, doJK);
  
  m_dataset->set_covariance(cov_mat);
}


// ============================================================================


void cosmobl::twopt::TwoPointCorrelation_deprojected::compute_covariance_matrix (const vector<string> file_xi, const bool doJK)
{
  vector<double> rad, mean;
  vector<vector<double>> cov_mat;

  cosmobl::covariance_matrix (file_xi, rad, mean, cov_mat, doJK); 
  m_dataset->set_covariance(cov_mat);
}
