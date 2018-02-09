/********************************************************************
 *  Copyright (C) 2015 by Federico Marulli                          *
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
 *  @file NumberCounts/NumberCounts.cpp
 *
 *  @brief Methods of the class NumberCounts 
 *
 *  This file contains the implementation of the methods of the class
 *  NumberCounts, used to handle catalogues of astronomical
 *  sources
 *
 *  @author Federico Marulli
 *
 *  @author federico.marulli3@unbo.it
 */

#include "NumberCounts.h"
#include "Data1D.h"

using namespace cosmobl;
using namespace catalogue;
using namespace measure::numbercounts;

// ============================================================================


shared_ptr<cosmobl::data::Data> cosmobl::measure::numbercounts::NumberCounts::m_measurePoisson (const int nbin, const cosmobl::binType bintype, const double Volume) const
{
  vector<double> variable = m_data->var(m_variable);
  vector<double> weight = m_data->var(catalogue::Var::_Weight_);

  if (variable.size()<=0 || nbin<=0) ErrorCBL("Error in m_measurePoisson of NumberCounts.cpp!");

  vector<double> xx(nbin, 0), fx(nbin, 0);
  vector<vector<double>> covariance(nbin, vector<double>(nbin, 0));
  
  double minVar = (m_MinVar>cosmobl::par::defaultDouble) ? m_MinVar : Min(variable)*0.9999;
  double maxVar = (m_MaxVar>cosmobl::par::defaultDouble) ? m_MaxVar : Max(variable)*1.0001;

  const bool linear = (bintype == binType::_linear_) ? true : false;

  const bool dndv = (m_ncType == NCType::_dn_dV_) ? true : false;
  
  // using GSL to create the histogram 

  gsl_histogram *nc = gsl_histogram_alloc(nbin);

  if (linear) gsl_histogram_set_ranges_uniform(nc, minVar, maxVar);

  else {
    vector<double> vv = logarithmic_bin_vector(nbin+1, minVar, maxVar);
    double *vvv = new double[nbin+1]; for (int i=0; i<nbin+1; i++) vvv[i] = vv[i];
    gsl_histogram_set_ranges(nc, vvv, nbin+1);
  }
  
  /// fill the histogram
  for (size_t i=0; i<variable.size(); i++)
    gsl_histogram_accumulate(nc, variable[i], weight[i]);
  
  double x1, x2;

  for (int i=0; i<nbin; i++) {

    gsl_histogram_get_range(nc, i, &x1, &x2);
    double val = gsl_histogram_get(nc, i);
    
    if (linear) xx[i] = 0.5*(x1+x2);
    else xx[i] = pow(10., 0.5*(log10(x1)+log10(x2)));

    if (dndv) {
      fx[i] = val/((x2-x1)*Volume);
      covariance[i][i] = pow(sqrt(val)/((x2-x1)*Volume), 2);
    }
    else {
      fx[i] = val/((log10(x2)-log10(x1))*Volume);
      covariance[i][i] = pow(sqrt(val)/((log10(x2)-log10(x1))*Volume), 2);
    }
    
  }

  gsl_histogram_free(nc);

  auto data = make_shared<data::Data1D>(data::Data1D(xx, fx, covariance));
  return data;
}


// ============================================================================


shared_ptr<cosmobl::data::Data>cosmobl::measure::numbercounts::NumberCounts::m_measureJackknife (const int nbin, const cosmobl::binType bintype, const double Volume, const string dir_output_resample) const
{
  vector<double> variable = m_data->var(m_variable);
  vector<double> weight = m_data->var(catalogue::Var::_Weight_);
  vector<long> regions = m_data->region_list();

  const int nRegions = (int)m_data->nRegions();

  if (variable.size()<=0 || nbin<=0 || nRegions<2) ErrorCBL("Error in m_measurePoisson of NumberCounts.cpp!");

  bool write_JK = false;
  if(dir_output_resample!=par::defaultString){
    string mkdir = "mkdir -p "+dir_output_resample;
    if (system(mkdir.c_str())) {}
    write_JK = true;
  }
  (void)write_JK;

  vector<double> xx(nbin, 0), fx(nbin, 0);
  
  double minVar = (m_MinVar>cosmobl::par::defaultDouble) ? m_MinVar : Min(variable)*0.9999;
  double maxVar = (m_MaxVar>cosmobl::par::defaultDouble) ? m_MaxVar : Max(variable)*1.0001;

  const bool linear = (bintype == binType::_linear_) ? true : false;

  const bool dndv = (m_ncType == NCType::_dn_dV_) ? true : false;
  
  // using GSL to create the histogram 

  vector<gsl_histogram *> nc_region;
  for(int i=0; i<nRegions; i++){
    gsl_histogram *nc = gsl_histogram_alloc(nbin);
    nc_region.push_back(nc);
  }

  if (linear) {
    for (int i=0; i<nRegions; i++)
      gsl_histogram_set_ranges_uniform(nc_region[i], minVar, maxVar);
  }

  else {
    vector<double> vv = logarithmic_bin_vector(nbin+1, minVar, maxVar);
    for (int i=0; i<nRegions; i++)
      gsl_histogram_set_ranges(nc_region[i], vv.data(), nbin+1);
  }

  
  /// Fill the histogram 
  for (size_t i=0; i<variable.size(); i++)
    gsl_histogram_accumulate(nc_region[regions[i]], variable[i], weight[i]);
  
  vector<vector<double>> fx_JK(nRegions, vector<double>(nbin, 0));

  /// set the measure
  for (int i=0; i<nbin; i++) {

    double x1, x2;
    double val = 0;

    for(int j=0; j<nRegions; j++)
      val += gsl_histogram_get(nc_region[j], i);
    
    gsl_histogram_get_range(nc_region[0], i, &x1, &x2);

    if (linear) xx[i] = 0.5*(x1+x2);
    else xx[i] = pow(10., 0.5*(log10(x1)+log10(x2)));

    if (dndv) fx[i] = val/((x2-x1)*Volume);
    
    else fx[i] = val/((log10(x2)-log10(x1))*Volume);

    for(int j=0; j<nRegions; j++){
      if (dndv) fx_JK[j][i]  = fx[i]-gsl_histogram_get(nc_region[j], i)/((x2-x1)*Volume);

      else fx_JK[j][i] =  fx[i]-gsl_histogram_get(nc_region[j], i)/((log10(x2)-log10(x1))*Volume);
    }

  }

  vector<vector<double>> covariance;
  covariance_matrix( fx_JK, covariance, true);

  if(write_JK)
    for(int i=0; i<nRegions; i++){
      auto data = make_shared<data::Data1D>(data::Data1D(xx, fx_JK[i], covariance));
      data->write(dir_output_resample, "number_counts_Jackknife_region="+conv(i, par::fINT)+".dat", "# var number_counts", 10);
    }

  auto data = make_shared<data::Data1D>(data::Data1D(xx, fx, covariance));
  return data;
}


// ============================================================================


shared_ptr<cosmobl::data::Data> cosmobl::measure::numbercounts::NumberCounts::m_measureBootstrap (const int nbin, const cosmobl::binType bintype, const double Volume, const string dir_output_resample, const int nMocks, const int seed) const 
{
  ErrorCBL("Error in m_measureBootstrap, work in progress!");
  (void)nbin; (void)bintype; (void)Volume;
  (void)dir_output_resample; (void)nMocks; (void)seed;

  auto data = make_shared<data::Data1D>(data::Data1D());
  return data;
}


// ============================================================================


shared_ptr<cosmobl::data::Data>cosmobl::measure::numbercounts::NumberCounts::m_measureJackknifeObjects (const int nbin, const cosmobl::binType bintype, const double Volume, const string dir_output_resample) const
{
  ErrorCBL("Error in m_measureJackknifeObjects, work in progress!");
  (void)nbin; (void)bintype; (void)Volume;
  (void)dir_output_resample;

  auto data = make_shared<data::Data1D>(data::Data1D());
  return data;
}


// ============================================================================


shared_ptr<cosmobl::data::Data> cosmobl::measure::numbercounts::NumberCounts::m_measureBootstrapObjects (const int nbin, const cosmobl::binType bintype, const double Volume, const string dir_output_resample, const int nMocks, const int seed) const 
{
  ErrorCBL("Error in m_measureBootstrapObjects, work in progress!");
  (void)nbin; (void)bintype; (void)Volume;
  (void)dir_output_resample; (void)nMocks; (void)seed;

  auto data = make_shared<data::Data1D>(data::Data1D());
  return data;
}

// ============================================================================


void cosmobl::measure::numbercounts::NumberCounts::measure (const cosmobl::measure::ErrorType errorType, const int nbin, const cosmobl::binType bintype, const double Volume, const string dir_output_resample, const int nMocks, const int seed)
{
  switch (errorType) {
    case (ErrorType::_Poisson_) :
      m_dataset = m_measurePoisson(nbin, bintype, Volume);
      break;
    case (ErrorType::_Jackknife_) :
      m_dataset = m_measureJackknife(nbin, bintype, Volume, dir_output_resample);
      break;
    case (ErrorType::_Bootstrap_) :
      m_dataset = m_measureBootstrap(nbin, bintype, Volume, dir_output_resample, nMocks, seed);
      break;

    default:
      ErrorCBL("Error in measure() of NumberCounts.cpp, unknown type of error");
  }
}


// ============================================================================


void cosmobl::measure::numbercounts::NumberCounts::read (const string dir, const string file)
{
  m_dataset->read(dir+file);
}


// ============================================================================

void cosmobl::measure::numbercounts::NumberCounts::write (const string dir, const string file, const int rank)
{
  vector<double> xx; m_dataset->xx(xx);

  string header = "[1] variable # [2] number counts # [3] error";
  
  m_dataset->write(dir, file, header, 10, rank);
}


// ============================================================================


void cosmobl::measure::numbercounts::NumberCounts::read_covariance (const string dir, const string file)
{
  m_dataset->set_covariance(dir+file);
}


// ============================================================================


void cosmobl::measure::numbercounts::NumberCounts::write_covariance (const string dir, const string file) const
{
  m_dataset->write_covariance(dir, file);
}


// ============================================================================


void cosmobl::measure::numbercounts::NumberCounts::compute_covariance (const vector<shared_ptr<data::Data>> number_counts, const bool JK)
{
  vector<vector<double>> nc;

  for (size_t i=0; i<number_counts.size(); i++) {
    vector<double> vv;
    number_counts[i]->data(vv);
    nc.push_back(vv);
  }

  vector<vector<double>> cov_mat;
  cosmobl::covariance_matrix(nc, cov_mat, JK);
  
  m_dataset->set_covariance(cov_mat);
}


// ============================================================================


void cosmobl::measure::numbercounts::NumberCounts::compute_covariance (const vector<string> file, const bool JK)
{
  vector<double> xx, fx_mean;
  vector<vector<double>> cov_mat;

  cosmobl::covariance_matrix(file, xx, fx_mean, cov_mat, JK);
  
  m_dataset->set_covariance(cov_mat);
}
