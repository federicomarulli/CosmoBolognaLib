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
 *  @authors federico.marulli3@unibo.it, alfonso.veropalumbo@unibo.it
 */


#include "TwoPointCorrelation1D_monopole.h"

using namespace std;

using namespace cbl;
using namespace catalogue;
using namespace chainmesh;
using namespace data;
using namespace pairs;
using namespace measure;
using namespace twopt;


// ============================================================================================


void cbl::measure::twopt::TwoPointCorrelation1D_monopole::set_parameters (const BinType binType, const double rMin, const double rMax, const int nbins, const double shift, const CoordinateUnits angularUnits, std::function<double(double)> angularWeight, const bool compute_extra_info) 
{
  if (!compute_extra_info) 
    m_dd = (binType==BinType::_logarithmic_) ? move(Pair::Create(PairType::_comoving_log_, PairInfo::_standard_, rMin, rMax, nbins, shift, angularUnits, angularWeight))
      : move(Pair::Create(PairType::_comoving_lin_, PairInfo::_standard_, rMin, rMax, nbins, shift, angularUnits, angularWeight));
  else 
    m_dd = (binType==BinType::_logarithmic_) ? move(Pair::Create(PairType::_comoving_log_, PairInfo::_extra_, rMin, rMax, nbins, shift, angularUnits, angularWeight))
      : move(Pair::Create(PairType::_comoving_lin_, PairInfo::_extra_, rMin, rMax, nbins, shift, angularUnits, angularWeight));
    
  m_rr = (binType==BinType::_logarithmic_) ? move(Pair::Create(PairType::_comoving_log_, PairInfo::_standard_, rMin, rMax, nbins, shift, angularUnits))
    : move(Pair::Create(PairType::_comoving_lin_, PairInfo::_standard_, rMin, rMax, nbins, shift, angularUnits));
  
  m_dr = (binType==BinType::_logarithmic_) ? move(Pair::Create(PairType::_comoving_log_, PairInfo::_standard_, rMin, rMax, nbins, shift, angularUnits))
    : move(Pair::Create(PairType::_comoving_lin_, PairInfo::_standard_, rMin, rMax, nbins, shift, angularUnits));
}


// ============================================================================================


void cbl::measure::twopt::TwoPointCorrelation1D_monopole::set_parameters (const BinType binType, const double rMin, const double rMax, const double binSize, const double shift, const CoordinateUnits angularUnits, std::function<double(double)> angularWeight, const bool compute_extra_info)
{
  if (!compute_extra_info) 
    m_dd = (binType==BinType::_logarithmic_) ? move(Pair::Create(PairType::_comoving_log_, PairInfo::_standard_, rMin, rMax, binSize, shift, angularUnits, angularWeight))
      : move(Pair::Create(PairType::_comoving_lin_, PairInfo::_standard_, rMin, rMax, binSize, shift, angularUnits, angularWeight));
  else 
    m_dd = (binType==BinType::_logarithmic_) ? move(Pair::Create(PairType::_comoving_log_, PairInfo::_extra_, rMin, rMax, binSize, shift, angularUnits, angularWeight))
      : move(Pair::Create(PairType::_comoving_lin_, PairInfo::_extra_, rMin, rMax, binSize, shift, angularUnits, angularWeight));
    
  m_rr = (binType==BinType::_logarithmic_) ? move(Pair::Create(PairType::_comoving_log_, PairInfo::_standard_, rMin, rMax, binSize, shift, angularUnits))
    : move(Pair::Create(PairType::_comoving_lin_, PairInfo::_standard_, rMin, rMax, binSize, shift, angularUnits));
  
  m_dr = (binType==BinType::_logarithmic_) ? move(Pair::Create(PairType::_comoving_log_, PairInfo::_standard_, rMin, rMax, binSize, shift, angularUnits))
    : move(Pair::Create(PairType::_comoving_lin_, PairInfo::_standard_, rMin, rMax, binSize, shift, angularUnits));
}


// ============================================================================================


void cbl::measure::twopt::TwoPointCorrelation1D_monopole::read (const std::string dir, const std::string file) 
{
  m_dataset->read(dir+file);
}


// ============================================================================================


void cbl::measure::twopt::TwoPointCorrelation1D_monopole::write (const std::string dir, const std::string file, const int rank) const 
{
  vector<double> xx = m_dataset->xx();

  checkDim(xx, m_dd->nbins(), "rad");

  string header = "[1] separation at the bin centre # [2] spherically averagerded two-point correlation function # [3] error";
  if (m_compute_extra_info) header += " # [4] mean separation # [5] standard deviation of the separation distribution # [6] mean redshift # [7] standard deviation of the redshift distribution";
  
  m_dataset->write(dir, file, header, 5, rank);
}


// ============================================================================================


void cbl::measure::twopt::TwoPointCorrelation1D_monopole::measure (const ErrorType errorType, const std::string dir_output_pairs, const std::vector<std::string> dir_input_pairs, const std::string dir_output_resample, const int nMocks, const bool count_dd, const bool count_rr, const bool count_dr, const bool tcount, const Estimator estimator, const double fact, const int seed)
{
  switch (errorType) {
    case (ErrorType::_Poisson_) :
      measurePoisson(dir_output_pairs, dir_input_pairs, count_dd, count_rr, count_dr, tcount, estimator, fact);
      break;
    case (ErrorType::_Jackknife_) :
      measureJackknife(dir_output_pairs, dir_input_pairs, dir_output_resample, count_dd, count_rr, count_dr, tcount, estimator, fact);
      break;
    case (ErrorType::_Bootstrap_) :
      measureBootstrap(nMocks, dir_output_pairs, dir_input_pairs, dir_output_resample, count_dd, count_rr, count_dr, tcount, estimator, fact, seed);
      break;
    case (ErrorType::_JackknifeTest_) :
      measureJackknifeTest(dir_output_pairs, dir_input_pairs, dir_output_resample, count_dd, count_rr, count_dr, tcount, estimator, fact);
      break;

    default:
      ErrorCBL("unknown type of error", "measure", "TwoPointCorrelation1D_monopole.cpp");
  }
}


// ============================================================================================


void cbl::measure::twopt::TwoPointCorrelation1D_monopole::measurePoisson (const std::string dir_output_pairs, const std::vector<std::string> dir_input_pairs, const bool count_dd, const bool count_rr, const bool count_dr, const bool tcount, const Estimator estimator, const double fact)
{
  // ----------- count the data-data, random-random and data-random pairs, or read them from file ----------- 
  
  count_allPairs(m_twoPType, dir_output_pairs, dir_input_pairs, count_dd, count_rr, count_dr, tcount, estimator, fact);
  
  
  // ----------- compute the monopole of the two-point correlation function ----------- 

  if (estimator==Estimator::_natural_)
    m_dataset = correlation_NaturalEstimator(m_dd, m_rr);
  else if (estimator==Estimator::_LandySzalay_)
    m_dataset = correlation_LandySzalayEstimator(m_dd, m_rr, m_dr);
  else
    ErrorCBL("the chosen estimator is not implemented!", "measurePoisson", "TwoPointCorrelation1D_monopole.cpp");
  
}


// ============================================================================================


void cbl::measure::twopt::TwoPointCorrelation1D_monopole::measureJackknife (const std::string dir_output_pairs, const std::vector<std::string> dir_input_pairs, const std::string dir_output_resample, const bool count_dd, const bool count_rr, const bool count_dr, const bool tcount, const Estimator estimator, const double fact)
{
  if (dir_output_resample!=par::defaultString && dir_output_resample!="") {
    string mkdir = "mkdir -p "+dir_output_resample;
    if (system(mkdir.c_str())) {}
  }
  
  vector<long> region_list = m_data->region_list();
  size_t nRegions = region_list.size();
  
  vector<vector<double> > xi_SubSamples, covariance;

  vector<shared_ptr<Pair> > dd_regions, rr_regions, dr_regions;

  count_allPairs_region(dd_regions, rr_regions, dr_regions, m_twoPType, dir_output_pairs, dir_input_pairs, count_dd, count_rr, count_dr, tcount, estimator, fact);
  
  vector<shared_ptr<Data> > data_SS = (estimator==Estimator::_natural_) ? XiJackknife(dd_regions, rr_regions) : XiJackknife(dd_regions, rr_regions, dr_regions);
  
  for (size_t i=0; i<nRegions; i++) {

    if (dir_output_resample!=par::defaultString && dir_output_resample!="") {
      string file = "xi_Jackknife_"+conv(i, par::fINT)+".dat";
      string header = "[1] separation at the bin centre # [2] spherically averagerded two-point correlation function # [3] error";
      if (m_compute_extra_info) header += " # [4] mean separation # [5] standard deviation of the separation distribution # [6] mean redshift # [7] standard deviation of the redshift distribution";
      data_SS[i]->write(dir_output_resample, file, header, 10, 0);
    }

    vector<double> dd; data_SS[i]->get_data(dd);

    xi_SubSamples.push_back(dd);
  }

  covariance_matrix(xi_SubSamples, covariance, 1);
  
  if (estimator==Estimator::_natural_)
    m_dataset = correlation_NaturalEstimator(m_dd, m_rr);
  else if (estimator==Estimator::_LandySzalay_)
    m_dataset = correlation_LandySzalayEstimator(m_dd, m_rr, m_dr);
  else
    ErrorCBL("the chosen estimator is not implemented!", "measureJackknife", "TwoPointCorrelation1D_monopole.cpp");
  
  m_dataset->set_covariance(covariance);
  
}


// ============================================================================================


void cbl::measure::twopt::TwoPointCorrelation1D_monopole::measureJackknifeTest (const std::string dir_output_pairs, const std::vector<std::string> dir_input_pairs, const std::string dir_output_resample, const bool count_dd, const bool count_rr, const bool count_dr, const bool tcount, const Estimator estimator, const double fact)
{
  if (dir_output_resample!=par::defaultString && dir_output_resample!="") {
    string mkdir = "mkdir -p "+dir_output_resample;
    if (system(mkdir.c_str())) {}
  }
  
  vector<long> region_list = m_data->region_list();
  size_t nRegions = region_list.size();
  vector<double> weights(nRegions, 1.);
  
  vector<vector<double> > xi_SubSamples, covariance;

  count_allPairs_region_test(m_twoPType, weights, dir_output_pairs, dir_input_pairs, count_dd, count_rr, count_dr, tcount, estimator, fact);
  
  vector<shared_ptr<Data> > data_SS = (estimator==Estimator::_natural_) ? XiJackknifeTest(m_dd_res, m_rr_res) : XiJackknifeTest(m_dd_res, m_rr_res, m_dr_res);
  
  for (size_t i=0; i<nRegions; i++) {

    if (dir_output_resample!=par::defaultString && dir_output_resample!="") {
      string file = "xi_Jackknife_"+conv(i, par::fINT)+".dat";
      string header = "[1] separation at the bin centre # [2] spherically averagerded two-point correlation function # [3] error";
      if (m_compute_extra_info) header += " # [4] mean separation # [5] standard deviation of the separation distribution # [6] mean redshift # [7] standard deviation of the redshift distribution";
      data_SS[i]->write(dir_output_resample, file, header, 10, 0);
    }

    vector<double> dd; data_SS[i]->get_data(dd);

    xi_SubSamples.push_back(dd);
  }

  covariance_matrix(xi_SubSamples, covariance, 1);
  
  if (estimator==Estimator::_natural_)
    m_dataset = correlation_NaturalEstimator(m_dd, m_rr);
  else if (estimator==Estimator::_LandySzalay_)
    m_dataset = correlation_LandySzalayEstimator(m_dd, m_rr, m_dr);
  else
    ErrorCBL("the chosen estimator is not implemented!", "measureJackknifeTest", "TwoPointCorrelation1D_monopole.cpp");
  
  m_dataset->set_covariance(covariance);
  
}

// ============================================================================================


void cbl::measure::twopt::TwoPointCorrelation1D_monopole::measureBootstrap (const int nMocks, const std::string dir_output_pairs, const std::vector<std::string> dir_input_pairs, const std::string dir_output_resample, const bool count_dd, const bool count_rr, const bool count_dr, const bool tcount, const Estimator estimator, const double fact, const int seed)
{
  if (nMocks<=0)
    ErrorCBL("number of mocks must be >0", "measureBootstrap", "TwoPointCorrelation1D_monopole.cpp");

  if (dir_output_resample!=par::defaultString && dir_output_resample!="") {
    string mkdir = "mkdir -p "+dir_output_resample;
    if (system(mkdir.c_str())) {}
  }

  vector<vector<double> > xi_SubSamples, covariance;

  vector<shared_ptr<Pair> > dd_regions, rr_regions, dr_regions;
  
  count_allPairs_region(dd_regions, rr_regions, dr_regions, m_twoPType, dir_output_pairs, dir_input_pairs, count_dd, count_rr, count_dr, tcount, estimator, fact);

  vector<shared_ptr<Data> > data_SS = (estimator==Estimator::_natural_) ? XiBootstrap(nMocks, dd_regions, rr_regions, seed) : XiBootstrap(nMocks, dd_regions, rr_regions, dr_regions, seed);

  for (int i=0; i<nMocks; i++) {

     if (dir_output_resample!=par::defaultString && dir_output_resample!="") {
      string file = "xi_Bootstrap_"+conv(i, par::fINT)+".dat";
      string header = "[1] separation at the bin centre # [2] spherically averagerded two-point correlation function # [3] error";
      if (m_compute_extra_info) header += " # [4] mean separation # [5] standard deviation of the separation distribution # [6] mean redshift # [7] standard deviation of the redshift distribution";
      data_SS[i]->write(dir_output_resample, file, header, 10, 0);
    }

    vector<double> dd; data_SS[i]->get_data(dd);

    xi_SubSamples.push_back(dd);
  }

  covariance_matrix(xi_SubSamples, covariance, 0);
  
  if (estimator==Estimator::_natural_)
    m_dataset = correlation_NaturalEstimator(m_dd, m_rr);
  else if (estimator==Estimator::_LandySzalay_)
    m_dataset = correlation_LandySzalayEstimator(m_dd, m_rr, m_dr);
  else
    ErrorCBL("the chosen estimator is not implemented!", "measureBootstrap", "TwoPointCorrelation1D_monopole.cpp");

  m_dataset->set_covariance(covariance);

}
