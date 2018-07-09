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
 *  CatalogueAnalysis/TwoPointCorrelation/TwoPointCorrelation1D_angular.cpp
 *
 *  @brief Methods of the class TwoPointCorrelation1D_angular used to
 *  measure the angular two-point correlation function
 *
 *  This file contains the implementation of the methods of the class
 *  TwoPointCorrelation1D_angular used to measure the angular
 *  two-point correlation function
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */


#include "TwoPointCorrelation1D_angular.h"

using namespace std;

using namespace cbl;
using namespace catalogue;
using namespace chainmesh;
using namespace data;
using namespace pairs;
using namespace measure;
using namespace twopt;


// ============================================================================================


void cbl::measure::twopt::TwoPointCorrelation1D_angular::set_parameters (const BinType binType, const double thetaMin, const double thetaMax, const int nbins, const double shift, const CoordinateUnits angularUnits, function<double(double)> angularWeight, const bool compute_extra_info) 
{
  if (!compute_extra_info) 
    m_dd = (binType==BinType::_logarithmic_) ? move(Pair::Create(PairType::_angular_log_, PairInfo::_standard_, thetaMin, thetaMax, nbins, shift, angularUnits, angularWeight))
      : move(Pair::Create(PairType::_angular_lin_, PairInfo::_standard_, thetaMin, thetaMax, nbins, shift, angularUnits, angularWeight));
  else 
    m_dd = (binType==BinType::_logarithmic_) ? move(Pair::Create(PairType::_angular_log_, PairInfo::_extra_, thetaMin, thetaMax, nbins, shift, angularUnits, angularWeight))
      : move(Pair::Create(PairType::_angular_lin_, PairInfo::_extra_, thetaMin, thetaMax, nbins, shift, angularUnits, angularWeight));
  
  m_rr = (binType==BinType::_logarithmic_) ? move(Pair::Create(PairType::_angular_log_, PairInfo::_standard_, thetaMin, thetaMax, nbins, shift, angularUnits))
    : move(Pair::Create(PairType::_angular_lin_, PairInfo::_extra_, thetaMin, thetaMax, nbins, shift, angularUnits));
  
  m_dr = (binType==BinType::_logarithmic_) ? move(Pair::Create(PairType::_angular_log_, PairInfo::_standard_, thetaMin, thetaMax, nbins, shift, angularUnits))
    : move(Pair::Create(PairType::_angular_lin_, PairInfo::_extra_, thetaMin, thetaMax, nbins, shift, angularUnits));
}


// ============================================================================================


void cbl::measure::twopt::TwoPointCorrelation1D_angular::set_parameters (const BinType binType, const double thetaMin, const double thetaMax, const double binSize, const double shift, const CoordinateUnits angularUnits, function<double(double)> angularWeight, const bool compute_extra_info)
{
  if (!compute_extra_info) 
    m_dd = (binType==BinType::_logarithmic_) ? move(Pair::Create(PairType::_angular_log_, PairInfo::_standard_, thetaMin, thetaMax, binSize, shift, angularUnits, angularWeight))
      : move(Pair::Create(PairType::_angular_lin_, PairInfo::_standard_, thetaMin, thetaMax, binSize, shift, angularUnits, angularWeight));
  else 
    m_dd = (binType==BinType::_logarithmic_) ? move(Pair::Create(PairType::_angular_log_, PairInfo::_extra_, thetaMin, thetaMax, binSize, shift, angularUnits, angularWeight))
      : move(Pair::Create(PairType::_angular_lin_, PairInfo::_extra_, thetaMin, thetaMax, binSize, shift, angularUnits, angularWeight));
  
  m_rr = (binType==BinType::_logarithmic_) ? move(Pair::Create(PairType::_angular_log_, PairInfo::_standard_, thetaMin, thetaMax, binSize, shift, angularUnits))
    : move(Pair::Create(PairType::_angular_lin_, PairInfo::_standard_, thetaMin, thetaMax, binSize, shift, angularUnits));
  
  m_dr = (binType==BinType::_logarithmic_) ? move(Pair::Create(PairType::_angular_log_, PairInfo::_standard_, thetaMin, thetaMax, binSize, shift, angularUnits))
    : move(Pair::Create(PairType::_angular_lin_, PairInfo::_standard_, thetaMin, thetaMax, binSize, shift, angularUnits));
}


// ============================================================================================


void cbl::measure::twopt::TwoPointCorrelation1D_angular::read (const string dir, const string file) 
{
  m_dataset->read(dir+file);
}


// ============================================================================================


void cbl::measure::twopt::TwoPointCorrelation1D_angular::write (const string dir, const string file, const int rank) const 
{
  vector<double> xx; m_dataset->xx(xx);

  checkDim(xx, m_dd->nbins(), "theta");

  string header = "[1] angular separation at the bin centre # [2] angular two-point correlation function # [3] error";
  if (m_compute_extra_info) header += " # [4] mean angular separation # [5] standard deviation of the distribution of angular separations # [6] mean redshift # [7] standard deviation of the redshift distribution";
  
  m_dataset->write(dir, file, header, 5, rank);
}


// ============================================================================================


void cbl::measure::twopt::TwoPointCorrelation1D_angular::measure (const ErrorType errorType, const string dir_output_pairs, const vector<string> dir_input_pairs, const string dir_output_resample, int nMocks, bool count_dd, const bool count_rr, const bool count_dr, const bool tcount, const Estimator estimator, const int seed)
{
  switch (errorType) {
  case(ErrorType::_Poisson_):
    measurePoisson(dir_output_pairs, dir_input_pairs, count_dd, count_rr, count_dr, tcount, estimator);
    break;
  case(ErrorType::_Jackknife_):
    measureJackknife(dir_output_pairs, dir_input_pairs, dir_output_resample, count_dd, count_rr, count_dr, tcount, estimator);
    break;
  case(ErrorType::_Bootstrap_):
    measureBootstrap(nMocks, dir_output_pairs, dir_input_pairs, dir_output_resample, count_dd, count_rr, count_dr, tcount, estimator, seed);
    break;
  default:
    ErrorCBL("Error in measure() of TwoPointCorrelation1D_angular.cpp, unknown type of error");
  }
}


// ============================================================================================


void cbl::measure::twopt::TwoPointCorrelation1D_angular::measurePoisson (const string dir_output_pairs, const vector<string> dir_input_pairs, const bool count_dd, const bool count_rr, const bool count_dr, const bool tcount, const Estimator estimator)
{ 
  // ----------- count the data-data, random-random and data-random pairs, or read them from file ----------- 
  
  count_allPairs(m_twoPType, dir_output_pairs, dir_input_pairs, count_dd, count_rr, count_dr, tcount, estimator);

  
  // ----------- compute the angular two-point correlation function ----------- 

  if (estimator==Estimator::_natural_)
    m_dataset = correlation_NaturalEstimator(m_dd, m_rr);
  else if (estimator==Estimator::_LandySzalay_)
    m_dataset = correlation_LandySzalayEstimator(m_dd, m_rr, m_dr);
  else
    ErrorCBL("Error in measurePoisson() of TwoPointCorrelation1D_angular.cpp: the chosen estimator is not implemented!");
}


// ============================================================================================


void cbl::measure::twopt::TwoPointCorrelation1D_angular::measureJackknife (const string dir_output_pairs, const vector<string> dir_input_pairs, const string dir_output_JackknifeXi, const bool count_dd, const bool count_rr, const bool count_dr, const bool tcount, const Estimator estimator)
{
  if (dir_output_JackknifeXi!=par::defaultString) {
    string mkdir = "mkdir -p "+dir_output_JackknifeXi;
    if (system(mkdir.c_str())) {}
  }
  
  vector<long> region_list = m_data->region_list();
  size_t nRegions = region_list.size();

  vector<vector<double> > xi_SubSamples, covariance;

  vector<shared_ptr<Pair> > dd_regions, rr_regions, dr_regions;
  
  count_allPairs_region(dd_regions, rr_regions, dr_regions, m_twoPType, dir_output_pairs, dir_input_pairs, count_dd, count_rr, count_dr, tcount, estimator);

  vector<shared_ptr<Data> > data_SS = (estimator==Estimator::_natural_) ? XiJackknife(dd_regions, rr_regions) : XiJackknife(dd_regions, rr_regions, dr_regions);

  for (size_t i=0; i<nRegions; i++) {

    if (dir_output_JackknifeXi!=par::defaultString && dir_output_JackknifeXi!="") {
      string file = "xi_Jackknife_"+conv(i, par::fINT)+".dat";
      string header = "[1] angular separation at the bin centre # [2] angular two-point correlation function # [3] error";
      if (m_compute_extra_info) header += " # [4] mean separation # [5] standard deviation of the separation distribution # [6] mean redshift # [7] standard deviation of the redshift distribution";
      data_SS[i]->write(dir_output_JackknifeXi, file, header, 10, 0);
    }

    vector<double> dd; data_SS[i]->data(dd);

    xi_SubSamples.push_back(dd);
  }

  covariance_matrix(xi_SubSamples, covariance, 1);

  if (estimator==Estimator::_natural_)
    m_dataset = correlation_NaturalEstimator(m_dd, m_rr);
  else if (estimator==Estimator::_LandySzalay_)
    m_dataset = correlation_LandySzalayEstimator(m_dd, m_rr, m_dr);
  else
    ErrorCBL("Error in measureJackknife() of TwoPointCorrelation1D_angular.cpp: the chosen estimator is not implemented!");

  m_dataset->set_covariance(covariance);

}


// ============================================================================================


void cbl::measure::twopt::TwoPointCorrelation1D_angular::measureBootstrap (const int nMocks, const string dir_output_pairs, const vector<string> dir_input_pairs, const string dir_output_BootstrapXi, const bool count_dd, const bool count_rr, const bool count_dr, const bool tcount, const Estimator estimator, const int seed)
{
  if (nMocks <=0)
    ErrorCBL("Error in measureBootstrap() of TwoPointCorrelation1D_angular.cpp, number of mocks must be >0");

  if (dir_output_BootstrapXi!=par::defaultString) {
    string mkdir = "mkdir -p "+dir_output_BootstrapXi;
    if (system(mkdir.c_str())) {}
  }

  vector<vector<double> > xi_SubSamples, covariance;

  vector<shared_ptr<Pair> > dd_regions, rr_regions, dr_regions;
  
  count_allPairs_region(dd_regions, rr_regions, dr_regions, m_twoPType, dir_output_pairs, dir_input_pairs, count_dd, count_rr, count_dr, tcount, estimator);

  vector<shared_ptr<Data> > data_SS = (estimator==Estimator::_natural_) ? XiBootstrap(nMocks, dd_regions, rr_regions, seed) : XiBootstrap(nMocks, dd_regions, rr_regions, dr_regions, seed);

  for (int i=0; i<nMocks; i++) {

     if (dir_output_BootstrapXi!=par::defaultString && dir_output_BootstrapXi!="") {
      string file = "xi_Bootstrap_"+conv(i, par::fINT)+".dat";
      string header = "[1] angular separation at the bin centre # [2] angular two-point correlation function # [3] error";
      if (m_compute_extra_info) header += " # [4] mean separation # [5] standard deviation of the separation distribution # [6] mean redshift # [7] standard deviation of the redshift distribution";
      data_SS[i]->write(dir_output_BootstrapXi, file, header, 10, 0);
    }

    vector<double> dd; data_SS[i]->data(dd);

    xi_SubSamples.push_back(dd);
  }

  covariance_matrix(xi_SubSamples, covariance, 0);

  if (estimator==Estimator::_natural_)
    m_dataset = correlation_NaturalEstimator(m_dd, m_rr);
  else if (estimator==Estimator::_LandySzalay_)
    m_dataset = correlation_LandySzalayEstimator(m_dd, m_rr, m_dr);
  else
    ErrorCBL("Error in measureBootstrap() of TwoPointCorrelation1D_angular.cpp: the chosen estimator is not implemented!");

  m_dataset->set_covariance(covariance);

}
