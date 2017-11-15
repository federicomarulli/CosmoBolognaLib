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
 *  CatalogueAnalysis/TwoPointCorrelation/TwoPointCorrelation1D_monopole.cpp
 *
 *  @brief Methods of the class TwoPointCorrelation1D_monopole used to
 *  measure the monopole of the two-point correlation function
 *
 *  This file contains the implementation of the methods of the class
 *  TwoPointCorrelation1D_monopole used to measure the monopole of the
 *  two-point correlation function
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */


#include "TwoPointCorrelation1D_monopole.h"

using namespace cosmobl;
using namespace catalogue;
using namespace chainmesh;
using namespace data;
using namespace pairs;
using namespace measure;
using namespace twopt;


// ============================================================================================


void cosmobl::measure::twopt::TwoPointCorrelation1D_monopole::set_parameters (const binType binType, const double rMin, const double rMax, const int nbins, const double shift, const CoordUnits angularUnits, function<double(double)> angularWeight, const bool compute_extra_info) 
{
  if (!compute_extra_info) 
    m_dd = (binType==_logarithmic_) ? move(Pair::Create(_comoving_log_, _standard_, rMin, rMax, nbins, shift, angularUnits, angularWeight))
      : move(Pair::Create(_comoving_lin_, _standard_, rMin, rMax, nbins, shift, angularUnits, angularWeight));
  else 
    m_dd = (binType==_logarithmic_) ? move(Pair::Create(_comoving_log_, _extra_, rMin, rMax, nbins, shift, angularUnits, angularWeight))
      : move(Pair::Create(_comoving_lin_, _extra_, rMin, rMax, nbins, shift, angularUnits, angularWeight));
    
  m_rr = (binType==_logarithmic_) ? move(Pair::Create(_comoving_log_, _standard_, rMin, rMax, nbins, shift, angularUnits))
    : move(Pair::Create(_comoving_lin_, _standard_, rMin, rMax, nbins, shift, angularUnits));
  
  m_dr = (binType==_logarithmic_) ? move(Pair::Create(_comoving_log_, _standard_, rMin, rMax, nbins, shift, angularUnits))
    : move(Pair::Create(_comoving_lin_, _standard_, rMin, rMax, nbins, shift, angularUnits));
}


// ============================================================================================


void cosmobl::measure::twopt::TwoPointCorrelation1D_monopole::set_parameters (const binType binType, const double rMin, const double rMax, const double binSize, const double shift, const CoordUnits angularUnits, function<double(double)> angularWeight, const bool compute_extra_info)
{
  if (!compute_extra_info) 
    m_dd = (binType==_logarithmic_) ? move(Pair::Create(_comoving_log_, _standard_, rMin, rMax, binSize, shift, angularUnits, angularWeight))
      : move(Pair::Create(_comoving_lin_, _standard_, rMin, rMax, binSize, shift, angularUnits, angularWeight));
  else 
    m_dd = (binType==_logarithmic_) ? move(Pair::Create(_comoving_log_, _extra_, rMin, rMax, binSize, shift, angularUnits, angularWeight))
      : move(Pair::Create(_comoving_lin_, _extra_, rMin, rMax, binSize, shift, angularUnits, angularWeight));
    
  m_rr = (binType==_logarithmic_) ? move(Pair::Create(_comoving_log_, _standard_, rMin, rMax, binSize, shift, angularUnits))
    : move(Pair::Create(_comoving_lin_, _standard_, rMin, rMax, binSize, shift, angularUnits));
  
  m_dr = (binType==_logarithmic_) ? move(Pair::Create(_comoving_log_, _standard_, rMin, rMax, binSize, shift, angularUnits))
    : move(Pair::Create(_comoving_lin_, _standard_, rMin, rMax, binSize, shift, angularUnits));
}


// ============================================================================================


void cosmobl::measure::twopt::TwoPointCorrelation1D_monopole::read (const string dir, const string file) 
{
  m_dataset->read(dir+file);
}


// ============================================================================================


void cosmobl::measure::twopt::TwoPointCorrelation1D_monopole::write (const string dir, const string file, const int rank) const 
{
  vector<double> xx; m_dataset->xx(xx);

  checkDim(xx, m_dd->nbins(), "rad");

  string header = "[1] separation at the bin centre # [2] spherically averagerded two-point correlation function # [3] error";
  if (m_compute_extra_info) header += " # [4] mean separation # [5] standard deviation of the separation distribution # [6] mean redshift # [7] standard deviation of the redshift distribution";
  
  m_dataset->write(dir, file, header, 5, rank);
}


// ============================================================================================


void cosmobl::measure::twopt::TwoPointCorrelation1D_monopole::measure (const ErrorType errorType, const string dir_output_pairs, const vector<string> dir_input_pairs, const string dir_output_resample, const int nMocks, const bool count_dd, const bool count_rr, const bool count_dr, const bool tcount, const Estimator estimator, const int seed)
{
  switch (errorType) {
    case (ErrorType::_Poisson_) :
      measurePoisson(dir_output_pairs, dir_input_pairs, count_dd, count_rr, count_dr, tcount, estimator);
      break;
    case (ErrorType::_Jackknife_) :
      measureJackknife(dir_output_pairs, dir_input_pairs, dir_output_resample, count_dd, count_rr, count_dr, tcount, estimator);
      break;
    case (ErrorType::_Bootstrap_) :
      measureBootstrap(nMocks, dir_output_pairs, dir_input_pairs, dir_output_resample, count_dd, count_rr, count_dr, tcount, estimator, seed);
      break;
    case (ErrorType::_JackknifeTest_) :
      measureJackknifeTest(dir_output_pairs, dir_input_pairs, dir_output_resample, count_dd, count_rr, count_dr, tcount, estimator);
      break;

    default:
      ErrorCBL("Error in measure() of TwoPointCorrelation1D_monopole.cpp, unknown type of error");
  }
}


// ============================================================================================


void cosmobl::measure::twopt::TwoPointCorrelation1D_monopole::measurePoisson (const string dir_output_pairs, const vector<string> dir_input_pairs, const bool count_dd, const bool count_rr, const bool count_dr, const bool tcount, const Estimator estimator)
{
  // ----------- count the data-data, random-random and data-random pairs, or read them from file ----------- 
  
  count_allPairs(m_twoPType, dir_output_pairs, dir_input_pairs, count_dd, count_rr, count_dr, tcount, estimator);
  
  
  // ----------- compute the monopole of the two-point correlation function ----------- 

  if (estimator==_natural_)
    m_dataset = correlation_NaturalEstimator(m_dd, m_rr);
  else if (estimator==_LandySzalay_)
    m_dataset = correlation_LandySzalayEstimator(m_dd, m_rr, m_dr);
  else
    ErrorCBL("Error in measurePoisson() of TwoPointCorrelation1D_monopole.cpp: the chosen estimator is not implemented!");
  
}


// ============================================================================================


void cosmobl::measure::twopt::TwoPointCorrelation1D_monopole::measureJackknife (const string dir_output_pairs, const vector<string> dir_input_pairs, const string dir_output_JackknifeXi, const bool count_dd, const bool count_rr, const bool count_dr, const bool tcount, const Estimator estimator)
{
  if (dir_output_JackknifeXi!=par::defaultString && dir_output_JackknifeXi!="") {
    string mkdir = "mkdir -p "+dir_output_JackknifeXi;
    if (system(mkdir.c_str())) {}
  }
  
  vector<long> region_list = m_data->region_list();
  size_t nRegions = region_list.size();
  
  vector<vector<double> > xi_SubSamples, covariance;

  vector<shared_ptr<Pair> > dd_regions, rr_regions, dr_regions;

  count_allPairs_region(dd_regions, rr_regions, dr_regions, m_twoPType, dir_output_pairs, dir_input_pairs, count_dd, count_rr, count_dr, tcount, estimator);
  
  vector<shared_ptr<Data> > data_SS = (estimator==_natural_) ? XiJackknife(dd_regions, rr_regions) : XiJackknife(dd_regions, rr_regions, dr_regions);
  
  for (size_t i=0; i<nRegions; i++) {

    if (dir_output_JackknifeXi!=par::defaultString && dir_output_JackknifeXi!="") {
      string file = "xi_Jackknife_"+conv(i, par::fINT)+".dat";
      string header = "[1] separation at the bin centre # [2] spherically averagerded two-point correlation function # [3] error";
      if (m_compute_extra_info) header += " # [4] mean separation # [5] standard deviation of the separation distribution # [6] mean redshift # [7] standard deviation of the redshift distribution";
      data_SS[i]->write(dir_output_JackknifeXi, file, header, 10, 0);
    }

    vector<double> dd; data_SS[i]->data(dd);

    xi_SubSamples.push_back(dd);
  }

  covariance_matrix(xi_SubSamples, covariance, 1);
  
  if (estimator==_natural_)
    m_dataset = correlation_NaturalEstimator(m_dd, m_rr);
  else if (estimator==_LandySzalay_)
    m_dataset = correlation_LandySzalayEstimator(m_dd, m_rr, m_dr);
  else
    ErrorCBL("Error in measureJackknife() of TwoPointCorrelation1D_monopole.cpp: the chosen estimator is not implemented!");
  
  m_dataset->set_covariance(covariance);
  
}


// ============================================================================================


void cosmobl::measure::twopt::TwoPointCorrelation1D_monopole::measureJackknifeTest (const string dir_output_pairs, const vector<string> dir_input_pairs, const string dir_output_JackknifeXi, const bool count_dd, const bool count_rr, const bool count_dr, const bool tcount, const Estimator estimator)
{
  if (dir_output_JackknifeXi!=par::defaultString && dir_output_JackknifeXi!="") {
    string mkdir = "mkdir -p "+dir_output_JackknifeXi;
    if (system(mkdir.c_str())) {}
  }
  
  vector<long> region_list = m_data->region_list();
  size_t nRegions = region_list.size();
  vector<double> weights(nRegions, 1.);
  
  vector<vector<double> > xi_SubSamples, covariance;

  count_allPairs_region_test(m_twoPType, weights, dir_output_pairs, dir_input_pairs, count_dd, count_rr, count_dr, tcount, estimator);
  
  vector<shared_ptr<Data> > data_SS = (estimator==_natural_) ? XiJackknifeTest(m_dd_res, m_rr_res) : XiJackknifeTest(m_dd_res, m_rr_res, m_dr_res);
  
  for (size_t i=0; i<nRegions; i++) {

    if (dir_output_JackknifeXi!=par::defaultString && dir_output_JackknifeXi!="") {
      string file = "xi_Jackknife_"+conv(i, par::fINT)+".dat";
      string header = "[1] separation at the bin centre # [2] spherically averagerded two-point correlation function # [3] error";
      if (m_compute_extra_info) header += " # [4] mean separation # [5] standard deviation of the separation distribution # [6] mean redshift # [7] standard deviation of the redshift distribution";
      data_SS[i]->write(dir_output_JackknifeXi, file, header, 10, 0);
    }

    vector<double> dd; data_SS[i]->data(dd);

    xi_SubSamples.push_back(dd);
  }

  covariance_matrix(xi_SubSamples, covariance, 1);
  
  if (estimator==_natural_)
    m_dataset = correlation_NaturalEstimator(m_dd, m_rr);
  else if (estimator==_LandySzalay_)
    m_dataset = correlation_LandySzalayEstimator(m_dd, m_rr, m_dr);
  else
    ErrorCBL("Error in measureJackknife() of TwoPointCorrelation1D_monopole.cpp: the chosen estimator is not implemented!");
  
  m_dataset->set_covariance(covariance);
  
}

// ============================================================================================


void cosmobl::measure::twopt::TwoPointCorrelation1D_monopole::measureBootstrap (const int nMocks, const string dir_output_pairs, const vector<string> dir_input_pairs, const string dir_output_BootstrapXi, const bool count_dd, const bool count_rr, const bool count_dr, const bool tcount, const Estimator estimator, const int seed)
{
  if (nMocks <=0)
    ErrorCBL("Error in measureBootstrap() of TwoPointCorrelation1D_monopole.cpp, number of mocks must be >0");

  if (dir_output_BootstrapXi!=par::defaultString && dir_output_BootstrapXi!="") {
    string mkdir = "mkdir -p "+dir_output_BootstrapXi;
    if (system(mkdir.c_str())) {}
  }

  vector<vector<double> > xi_SubSamples, covariance;

  vector<shared_ptr<Pair> > dd_regions, rr_regions, dr_regions;
  
  count_allPairs_region(dd_regions, rr_regions, dr_regions, m_twoPType, dir_output_pairs, dir_input_pairs, count_dd, count_rr, count_dr, tcount, estimator);

  vector<shared_ptr<Data> > data_SS = (estimator==_natural_) ? XiBootstrap(nMocks, dd_regions, rr_regions, seed) : XiBootstrap(nMocks, dd_regions, rr_regions, dr_regions, seed);

  for (int i=0; i<nMocks; i++) {

     if (dir_output_BootstrapXi!=par::defaultString && dir_output_BootstrapXi!="") {
      string file = "xi_Bootstrap_"+conv(i, par::fINT)+".dat";
      string header = "[1] separation at the bin centre # [2] spherically averagerded two-point correlation function # [3] error";
      if (m_compute_extra_info) header += " # [4] mean separation # [5] standard deviation of the separation distribution # [6] mean redshift # [7] standard deviation of the redshift distribution";
      data_SS[i]->write(dir_output_BootstrapXi, file, header, 10, 0);
    }

    vector<double> dd; data_SS[i]->data(dd);

    xi_SubSamples.push_back(dd);
  }

  covariance_matrix(xi_SubSamples, covariance, 0);
  
  if (estimator==_natural_)
    m_dataset = correlation_NaturalEstimator(m_dd, m_rr);
  else if (estimator==_LandySzalay_)
    m_dataset = correlation_LandySzalayEstimator(m_dd, m_rr, m_dr);
  else
    ErrorCBL("Error in measureBootstrap() of TwoPointCorrelation1D_monopole.cpp: the chosen estimator is not implemented!");

  m_dataset->set_covariance(covariance);

}
