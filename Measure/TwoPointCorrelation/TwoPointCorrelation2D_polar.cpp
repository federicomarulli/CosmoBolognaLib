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
 *  CatalogueAnalysis/TwoPointCorrelation/TwoPointCorrelation2D_polar.cpp
 *
 *  @brief Methods of the class TwoPointCorrelation2D_polar used to
 *  measure the monopole of the two-point correlation function
 *
 *  This file contains the implementation of the methods of the class
 *  TwoPointCorrelation2D_polar used to measure the monopole of the
 *  two-point correlation function
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */


#include "TwoPointCorrelation2D_polar.h"

using namespace std;

using namespace cbl;
using namespace catalogue;
using namespace chainmesh;
using namespace data;
using namespace pairs;
using namespace measure;
using namespace twopt;


// ============================================================================================


void cbl::measure::twopt::TwoPointCorrelation2D_polar::set_parameters (const BinType binType_rad, const double rMin, const double rMax, const int nbins_rad, const double shift_rad, const BinType binType_mu, const double muMin, const double muMax, const int nbins_mu, const double shift_mu, const CoordinateUnits angularUnits, std::function<double(double)> angularWeight, const bool compute_extra_info) 
{
  if (muMin<0.) ErrorCBL("Error in cbl::measure::twopt::TwoPointCorrelation2D_polar::set_parameter() of TwoPointCorrelation2D_polar.cpp: mMun must be >0 !");
  if (muMin>1.) ErrorCBL("Error in cbl::measure::twopt::TwoPointCorrelation2D_polar::set_parameter() of TwoPointCorrelation2D_polar.cpp: mMun must be <1 !");
  if (rMin<0.) ErrorCBL("Error in cbl::measure::twopt::TwoPointCorrelation2D_polar::set_parameter() of TwoPointCorrelation2D_polar.cpp: rMun must be >0 !");
  
  if (!compute_extra_info) {
    if (binType_rad==BinType::_logarithmic_) {
      if (binType_mu==BinType::_logarithmic_) {
	m_dd = move(Pair::Create(PairType::_comovingPolar_loglog_, PairInfo::_standard_, rMin, rMax, nbins_rad, shift_rad, muMin, muMax, nbins_mu, shift_mu, angularUnits, angularWeight));
      }
      else {
	m_dd = move(Pair::Create(PairType::_comovingPolar_loglin_, PairInfo::_standard_, rMin, rMax, nbins_rad, shift_rad, muMin, muMax, nbins_mu, shift_mu, angularUnits, angularWeight));
      }
    }
    else {
      if (binType_mu==BinType::_logarithmic_) {
	m_dd = move(Pair::Create(PairType::_comovingPolar_linlog_, PairInfo::_standard_, rMin, rMax, nbins_rad, shift_rad, muMin, muMax, nbins_mu, shift_mu, angularUnits, angularWeight));
      }
      else {
	m_dd = move(Pair::Create(PairType::_comovingPolar_linlin_, PairInfo::_standard_, rMin, rMax, nbins_rad, shift_rad, muMin, muMax, nbins_mu, shift_mu, angularUnits, angularWeight));
      }
    }
  }
  else {
    if (binType_rad==BinType::_logarithmic_) {
      if (binType_mu==BinType::_logarithmic_) {
	m_dd = move(Pair::Create(PairType::_comovingPolar_loglog_,PairInfo::_extra_, rMin, rMax, nbins_rad, shift_rad, muMin, muMax, nbins_mu, shift_mu, angularUnits, angularWeight));
      }
      else {
	m_dd = move(Pair::Create(PairType::_comovingPolar_loglin_, PairInfo::_extra_, rMin, rMax, nbins_rad, shift_rad, muMin, muMax, nbins_mu, shift_mu, angularUnits, angularWeight));
      }
    }
    else {
      if (binType_mu==BinType::_logarithmic_) {
	m_dd = move(Pair::Create(PairType::_comovingPolar_linlog_, PairInfo::_extra_, rMin, rMax, nbins_rad, shift_rad, muMin, muMax, nbins_mu, shift_mu, angularUnits, angularWeight));
      }
      else {
	m_dd = move(Pair::Create(PairType::_comovingPolar_linlin_, PairInfo::_extra_, rMin, rMax, nbins_rad, shift_rad, muMin, muMax, nbins_mu, shift_mu, angularUnits, angularWeight));
      }
    }
  }

  
  if (binType_rad==BinType::_logarithmic_) {
    if (binType_mu==BinType::_logarithmic_) {
      m_rr = move(Pair::Create(PairType::_comovingPolar_loglog_, PairInfo::_standard_, rMin, rMax, nbins_rad, shift_rad, muMin, muMax, nbins_mu, shift_mu, angularUnits));
      m_dr = move(Pair::Create(PairType::_comovingPolar_loglog_, PairInfo::_standard_, rMin, rMax, nbins_rad, shift_rad, muMin, muMax, nbins_mu, shift_mu, angularUnits));
    }
    else {
      m_rr = move(Pair::Create(PairType::_comovingPolar_loglin_, PairInfo::_standard_, rMin, rMax, nbins_rad, shift_rad, muMin, muMax, nbins_mu, shift_mu, angularUnits));
      m_dr = move(Pair::Create(PairType::_comovingPolar_loglin_, PairInfo::_standard_, rMin, rMax, nbins_rad, shift_rad, muMin, muMax, nbins_mu, shift_mu, angularUnits));
    }
  }
  else {
    if (binType_mu==BinType::_logarithmic_) {
      m_rr = move(Pair::Create(PairType::_comovingPolar_linlog_, PairInfo::_standard_, rMin, rMax, nbins_rad, shift_rad, muMin, muMax, nbins_mu, shift_mu, angularUnits));
      m_dr = move(Pair::Create(PairType::_comovingPolar_linlog_, PairInfo::_standard_, rMin, rMax, nbins_rad, shift_rad, muMin, muMax, nbins_mu, shift_mu, angularUnits));
    }
    else {
      m_rr = move(Pair::Create(PairType::_comovingPolar_linlin_, PairInfo::_standard_, rMin, rMax, nbins_rad, shift_rad, muMin, muMax, nbins_mu, shift_mu, angularUnits));
      m_dr = move(Pair::Create(PairType::_comovingPolar_linlin_, PairInfo::_standard_, rMin, rMax, nbins_rad, shift_rad, muMin, muMax, nbins_mu, shift_mu, angularUnits));
    }
  }
 
}


// ============================================================================================


void cbl::measure::twopt::TwoPointCorrelation2D_polar::set_parameters (const BinType binType_rad, const double rMin, const double rMax, const double binSize_rad, const double shift_rad, const BinType binType_mu, const double muMin, const double muMax, const double binSize_mu, const double shift_mu, const CoordinateUnits angularUnits, std::function<double(double)> angularWeight, const bool compute_extra_info)
{
  if (muMin<0.) ErrorCBL("Error in cbl::measure::twopt::TwoPointCorrelation2D_polar::set_parameter() of TwoPointCorrelation2D_polar.cpp: mMun must be >0 !");
  if (muMin>1.) ErrorCBL("Error in cbl::measure::twopt::TwoPointCorrelation2D_polar::set_parameter() of TwoPointCorrelation2D_polar.cpp: mMun must be <1 !");
  if (rMin<0.) ErrorCBL("Error in cbl::measure::twopt::TwoPointCorrelation2D_polar::set_parameter() of TwoPointCorrelation2D_polar.cpp: rMun must be >0 !");
  
  if (!compute_extra_info) {
    if (binType_rad==BinType::_logarithmic_) {
      if (binType_mu==BinType::_logarithmic_) {
	m_dd = move(Pair::Create(PairType::_comovingPolar_loglog_, PairInfo::_standard_, rMin, rMax, binSize_rad, shift_rad, muMin, muMax, binSize_mu, shift_mu, angularUnits, angularWeight));
      }
      else {
	m_dd = move(Pair::Create(PairType::_comovingPolar_loglin_, PairInfo::_standard_, rMin, rMax, binSize_rad, shift_rad, muMin, muMax, binSize_mu, shift_mu, angularUnits, angularWeight));
      }
    }
    else {
      if (binType_mu==BinType::_logarithmic_) {
	m_dd = move(Pair::Create(PairType::_comovingPolar_linlog_, PairInfo::_standard_, rMin, rMax, binSize_rad, shift_rad, muMin, muMax, binSize_mu, shift_mu, angularUnits, angularWeight));
      }
      else {
	m_dd = move(Pair::Create(PairType::_comovingPolar_linlin_, PairInfo::_standard_, rMin, rMax, binSize_rad, shift_rad, muMin, muMax, binSize_mu, shift_mu, angularUnits, angularWeight));
      }
    }
  }
  else {
    if (binType_rad==BinType::_logarithmic_) {
      if (binType_mu==BinType::_logarithmic_) {
	m_dd = move(Pair::Create(PairType::_comovingPolar_loglog_,PairInfo::_extra_, rMin, rMax, binSize_rad, shift_rad, muMin, muMax, binSize_mu, shift_mu, angularUnits, angularWeight));
      }
      else {
	m_dd = move(Pair::Create(PairType::_comovingPolar_loglin_, PairInfo::_extra_, rMin, rMax, binSize_rad, shift_rad, muMin, muMax, binSize_mu, shift_mu, angularUnits, angularWeight));
      }
    }
    else {
      if (binType_mu==BinType::_logarithmic_) {
	m_dd = move(Pair::Create(PairType::_comovingPolar_linlog_, PairInfo::_extra_, rMin, rMax, binSize_rad, shift_rad, muMin, muMax, binSize_mu, shift_mu, angularUnits, angularWeight));
      }
      else {
	m_dd = move(Pair::Create(PairType::_comovingPolar_linlin_, PairInfo::_extra_, rMin, rMax, binSize_rad, shift_rad, muMin, muMax, binSize_mu, shift_mu, angularUnits, angularWeight));
      }
    }
  }


  if (binType_rad==BinType::_logarithmic_) {
    if (binType_mu==BinType::_logarithmic_) {
      m_rr = move(Pair::Create(PairType::_comovingPolar_loglog_, PairInfo::_standard_, rMin, rMax, binSize_rad, shift_rad, muMin, muMax, binSize_mu, shift_mu, angularUnits));
      m_dr = move(Pair::Create(PairType::_comovingPolar_loglog_, PairInfo::_standard_, rMin, rMax, binSize_rad, shift_rad, muMin, muMax, binSize_mu, shift_mu, angularUnits));
    }
    else {
      m_rr = move(Pair::Create(PairType::_comovingPolar_loglin_, PairInfo::_standard_, rMin, rMax, binSize_rad, shift_rad, muMin, muMax, binSize_mu, shift_mu, angularUnits));
      m_dr = move(Pair::Create(PairType::_comovingPolar_loglin_, PairInfo::_standard_, rMin, rMax, binSize_rad, shift_rad, muMin, muMax, binSize_mu, shift_mu, angularUnits));
    }
  }
  else {
    if (binType_mu==BinType::_logarithmic_) {
      m_rr = move(Pair::Create(PairType::_comovingPolar_linlog_, PairInfo::_standard_, rMin, rMax, binSize_rad, shift_rad, muMin, muMax, binSize_mu, shift_mu, angularUnits));
      m_dr = move(Pair::Create(PairType::_comovingPolar_linlog_, PairInfo::_standard_, rMin, rMax, binSize_rad, shift_rad, muMin, muMax, binSize_mu, shift_mu, angularUnits));
    }
    else {
      m_rr = move(Pair::Create(PairType::_comovingPolar_linlin_, PairInfo::_standard_, rMin, rMax, binSize_rad, shift_rad, muMin, muMax, binSize_mu, shift_mu, angularUnits));
      m_dr = move(Pair::Create(PairType::_comovingPolar_linlin_, PairInfo::_standard_, rMin, rMax, binSize_rad, shift_rad, muMin, muMax, binSize_mu, shift_mu, angularUnits));
    }
  }
}


// ============================================================================================


void cbl::measure::twopt::TwoPointCorrelation2D_polar::read (const std::string dir, const std::string file) 
{
  m_dataset->read(dir+file);
}


// ============================================================================================


void cbl::measure::twopt::TwoPointCorrelation2D_polar::write (const std::string dir, const std::string file, const bool full, const int rank) const 
{
  vector<double> xx = m_dataset->xx(), yy = m_dataset->yy();

  checkDim(xx, m_dd->nbins_D1(), "rp"); checkDim(yy, m_dd->nbins_D2(), "pi");

  string header = "[1] separation at the bin centre # [2] angular separation at the bin centre # [3] 2D two-point correlation function # [4] error";
  if (m_compute_extra_info) header += " # [5] mean separation # [6] standard deviation of the separation distribution # [7] mean angular separation # [8] standard deviation of the angular separation distribution # [9] mean redshift # [10] standard deviation of the redshift distribution";

  m_dataset->write(dir, file, header, full, 5, rank);
}


// ============================================================================================


void cbl::measure::twopt::TwoPointCorrelation2D_polar::measure (const ErrorType errorType, const std::string dir_output_pairs, const std::vector<std::string> dir_input_pairs, const std::string dir_output_resample, int nMocks, bool count_dd, const bool count_rr, const bool count_dr, const bool tcount, const Estimator estimator, const int seed)
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
  default:
    ErrorCBL("Error in measure() of TwoPointCorrelation2D_polar.cpp, unknown type of error");
  }
}


// ============================================================================================


void cbl::measure::twopt::TwoPointCorrelation2D_polar::measurePoisson (const std::string dir_output_pairs, const std::vector<std::string> dir_input_pairs, const bool count_dd, const bool count_rr, const bool count_dr, const bool tcount, const Estimator estimator)
{ 
  // ----------- count the data-data, random-random and data-random pairs, or read them from file ----------- 
  
  count_allPairs(m_twoPType, dir_output_pairs, dir_input_pairs, count_dd, count_rr, count_dr, tcount, estimator);

  
  // ----------- compute the monopole of the two-point correlation function ----------- 

  if (estimator==Estimator::_natural_)
    m_dataset = correlation_NaturalEstimator(m_dd, m_rr);
  else if (estimator==Estimator::_LandySzalay_)
    m_dataset = correlation_LandySzalayEstimator(m_dd, m_rr, m_dr);
  else
    ErrorCBL("Error in measurePoisson() of TwoPointCorrelation2D_polar.cpp: the chosen estimator is not implemented!");

}


// ============================================================================================


void cbl::measure::twopt::TwoPointCorrelation2D_polar::measureJackknife (const std::string dir_output_pairs, const std::vector<std::string> dir_input_pairs, const std::string dir_output_JackknifeXi, const bool count_dd, const bool count_rr, const bool count_dr, const bool tcount, const Estimator estimator)
{
  if (dir_output_JackknifeXi!=par::defaultString) { 
    string mkdir = "mkdir -p "+dir_output_JackknifeXi;
    if (system(mkdir.c_str())) {}
  }
  
  vector<long> region_list = m_data->region_list();
  size_t nRegions = region_list.size();

  vector< vector<double > > xi_SubSamples, covariance;

  vector<shared_ptr<Pair> > dd_regions, rr_regions, dr_regions;
  count_allPairs_region(dd_regions, rr_regions, dr_regions, m_twoPType, dir_output_pairs, dir_input_pairs, count_dd, count_rr, count_dr, tcount, estimator);

  vector<shared_ptr<Data> > data_SS = (estimator==Estimator::_natural_) ? XiJackknife(dd_regions, rr_regions) : XiJackknife(dd_regions, rr_regions, dr_regions);

  for (size_t i=0; i<nRegions; i++) {
    
    if (dir_output_JackknifeXi!=par::defaultString && dir_output_JackknifeXi!="") {
      string file = "xi_Jackknife_"+conv(i, par::fINT)+".dat";
      string header = "[1] separation at the bin centre # [2] angular separation at the bin centre # [3] 2D two-point correlation function # [4] error";
      if (m_compute_extra_info) header += " # [5] mean separation # [6] standard deviation of the separation distribution # [7] mean angular separation # [8] standard deviation of the angular separation distribution # [9] mean redshift # [10] standard deviation of the redshift distribution";
      data_SS[i]->write(dir_output_JackknifeXi, file, header, false, 10, 0);
    }

    xi_SubSamples.push_back(data_SS[i]->data());
  }

  covariance_matrix(xi_SubSamples, covariance, 1);

  if (estimator==Estimator::_natural_)
    m_dataset = correlation_NaturalEstimator(m_dd, m_rr);
  else if (estimator==Estimator::_LandySzalay_)
    m_dataset = correlation_LandySzalayEstimator(m_dd, m_rr, m_dr);
  else
    ErrorCBL("Error in measurePoisson() of TwoPointCorrelation2D_polar.cpp: the chosen estimator is not implemented!");

  m_dataset->set_covariance(covariance);
}


// ============================================================================================


void cbl::measure::twopt::TwoPointCorrelation2D_polar::measureBootstrap (const int nMocks, const std::string dir_output_pairs, const std::vector<std::string> dir_input_pairs, const std::string dir_output_BootstrapXi, const bool count_dd, const bool count_rr, const bool count_dr, const bool tcount, const Estimator estimator, const int seed)
{
  if (nMocks<=0)
    ErrorCBL("Error in measureBootstrap() of TwoPointCorrelation2D_polar.cpp, number of mocks must be >0");

  if (dir_output_BootstrapXi!=par::defaultString) {
    string mkdir = "mkdir -p "+dir_output_BootstrapXi;
    if (system(mkdir.c_str())) {}
  }

  vector< vector<double > > xi_SubSamples, covariance;

  vector<shared_ptr<Pair> > dd_regions, rr_regions, dr_regions;
  count_allPairs_region(dd_regions, rr_regions, dr_regions, m_twoPType, dir_output_pairs, dir_input_pairs, count_dd, count_rr, count_dr, tcount, estimator);

  vector<shared_ptr<Data> > data_SS = (estimator==Estimator::_natural_) ? XiBootstrap(nMocks, dd_regions, rr_regions, seed) : XiBootstrap(nMocks, dd_regions, rr_regions, dr_regions, seed);

  for (int i=0; i<nMocks; i++) {

    if (dir_output_BootstrapXi!=par::defaultString && dir_output_BootstrapXi!="") {
      string file = "xi_Bootstrap_"+conv(i, par::fINT)+".dat";
      string header = "[1] separation at the bin centre # [2] angular separation at the bin centre # [3] 2D two-point correlation function # [4] error";
      if (m_compute_extra_info) header += " # [5] mean separation # [6] standard deviation of the separation distribution # [7] mean angular separation # [8] standard deviation of the angular separation distribution # [9] mean redshift # [10] standard deviation of the redshift distribution";
      data_SS[i]->write(dir_output_BootstrapXi, file, header, false, 10, 0);
    }

    xi_SubSamples.push_back(data_SS[i]->data());
  }

  covariance_matrix(xi_SubSamples, covariance, 0);
  

  if (estimator==Estimator::_natural_)
    m_dataset = correlation_NaturalEstimator(m_dd, m_rr);
  else if (estimator==Estimator::_LandySzalay_)
    m_dataset = correlation_LandySzalayEstimator(m_dd, m_rr, m_dr);
  else
    ErrorCBL("Error in measurePoisson() of TwoPointCorrelation2D_polar.cpp: the chosen estimator is not implemented!");

  m_dataset->set_covariance(covariance);
}
