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
 *  CatalogueAnalysis/TwoPointCorrelation/TwoPointCorrelation_projected.cpp
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
 *  @authors federico.marulli3@unibo.it, alfonso.veropalumbo@unibo.it
 */


#include "TwoPointCorrelation_projected.h"
#include "Data1D_extra.h"

using namespace std;

using namespace cbl;
using namespace catalogue;
using namespace chainmesh;
using namespace pairs;
using namespace measure;
using namespace twopt;


// ============================================================================


std::shared_ptr<data::Data> cbl::measure::twopt::TwoPointCorrelation_projected::data_with_extra_info (const std::vector<double> rp, const std::vector<double> ww, const std::vector<double> error) const
{
  auto dd2D = cbl::measure::twopt::TwoPointCorrelation2D_cartesian::m_dd;
  
  vector<double> weightTOT(dd2D->nbins_D1(), 0.), scale_mean(dd2D->nbins_D1(), 0.), scale_sigma(dd2D->nbins_D1(), 0.), z_mean(dd2D->nbins_D1(), 0.), z_sigma(dd2D->nbins_D1(), 0.);

  double fact_err, fact_scale, fact_z;
  
  for (int i=0; i<dd2D->nbins_D1(); ++i) {

    for (int j=0; j<dd2D->nbins_D2(); ++j)
      weightTOT[i] += dd2D->PP2D_weighted(i, j);
    
    for (int j=0; j<dd2D->nbins_D2(); ++j) {
      scale_mean[i] += dd2D->scale_D1_mean(i, j)*dd2D->PP2D_weighted(i, j)/weightTOT[i];
      z_mean[i] += dd2D->z_mean(i, j)*dd2D->PP2D_weighted(i, j)/weightTOT[i];
    }

    scale_sigma[i] = pow(dd2D->scale_D1_sigma(i, 0), 2)*dd2D->PP2D_weighted(i, 0);
    z_sigma[i] = pow(dd2D->z_sigma(i, 0), 2)*dd2D->PP2D_weighted(i, 0);
    
    for (int j=1; j<dd2D->nbins_D2(); ++j) {
      if (dd2D->PP2D_weighted(i, j)>0) {
	fact_err = dd2D->PP2D_weighted(i, j)*dd2D->PP2D_weighted(i, j-1)/(dd2D->PP2D_weighted(i, j)+dd2D->PP2D_weighted(i, j-1));
	fact_scale = pow(dd2D->scale_D1_mean(i, j)-dd2D->scale_D1_mean(i, j-1), 2)*fact_err;
	fact_z = pow(dd2D->z_mean(i, j)-dd2D->z_mean(i, j-1), 2)*fact_err;
	scale_sigma[i] += pow(dd2D->scale_D1_sigma(i, j), 2)*dd2D->PP2D_weighted(i, j)+fact_scale;
	z_sigma[i] += pow(dd2D->z_sigma(i, j),2)*weightTOT[i]+fact_z;
      }
    }
  }
  
  vector<vector<double>> extra(4);
  
  for (int i=0; i<dd2D->nbins_D1(); ++i) {
    extra[0].push_back(scale_mean[i]);
    extra[1].push_back(sqrt(scale_sigma[i]/weightTOT[i]));
    extra[2].push_back(z_mean[i]);
    extra[3].push_back(sqrt(z_sigma[i]/weightTOT[i]));
  }
    
  return move(unique_ptr<data::Data1D_extra>(new data::Data1D_extra(rp, ww, error, extra)));
}


// ============================================================================================


std::shared_ptr<data::Data> cbl::measure::twopt::TwoPointCorrelation_projected::Projected (const std::vector<double> rp, const std::vector<double> pi, const std::vector<std::vector<double>> xi, const std::vector<std::vector<double>> error_xi)
{
  vector<double> ww, error;

  ww.resize(0); ww.resize(rp.size(), 0.);
  error.resize(0); error.resize(rp.size(), 0.);
 
  const double binSize = 1./m_dd->binSize_inv_D2();

  const int pim = nint((m_piMax_integral-Min(pi))/binSize); // to convert from Mpc/h into the vector index
  
  for (size_t i=0; i<rp.size(); i++) {

    ww[i] = 0.;
    error[i] = 0.;
    
    for (int j=0; j<pim; j++) {  
      ww[i] = ww[i]+2.*binSize*xi[i][j];
      if (ww[i]>-1.) error[i] += pow(2.*binSize*error_xi[i][j], 2); // check!!!!
    }
  }

  for_each( error.begin(), error.end(), [] (double &vv) { vv = sqrt(vv);} );

  return (!m_compute_extra_info) ? move(unique_ptr<data::Data1D>(new data::Data1D(rp, ww, error))) : data_with_extra_info(rp, ww, error);

}


// ============================================================================================


void cbl::measure::twopt::TwoPointCorrelation_projected::read (const std::string dir, const std::string file) 
{
  m_dataset->read(dir+file);
}


// ============================================================================================


void cbl::measure::twopt::TwoPointCorrelation_projected::write (const std::string dir, const std::string file, const int rank) const 
{
  vector<double> xx = m_dataset->xx();

  checkDim(xx, m_dd->nbins_D1(), "rp");

  string header = "[1] perpendicular separation at the bin centre # [2] projected two-point correlation function # [3] error";
  if (m_compute_extra_info) header += " # [4] mean perpendicular separation # [5] standard deviation of the distribution of perpendicular separations # [6] mean redshift # [7] standard deviation of the redshift distribution";

  m_dataset->write(dir, file, header, 5, rank);
}


// ============================================================================================


void cbl::measure::twopt::TwoPointCorrelation_projected::measure (const ErrorType errorType, const std::string dir_output_pairs, const std::vector<std::string> dir_input_pairs, const std::string dir_output_resample, const int nMocks, const bool count_dd, const bool count_rr, const bool count_dr, const bool tcount, const Estimator estimator, const double fact, const int seed)
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
    default:
      ErrorCBL("unknown type of error!", "measure", "TwoPointCorrelation_projected.cpp");
  }
}


// ============================================================================================


void cbl::measure::twopt::TwoPointCorrelation_projected::measurePoisson (const std::string dir_output_pairs, const std::vector<std::string> dir_input_pairs, const bool count_dd, const bool count_rr, const bool count_dr, const bool tcount, const Estimator estimator, const double fact)
{
  // ----------- measure the 2D two-point correlation function, xi(rp,pi) ----------- 

  TwoPointCorrelation2D_cartesian::measurePoisson(dir_output_pairs, dir_input_pairs, count_dd, count_rr, count_dr, tcount, estimator, fact);

  
  // ----------- integrate the 2D two-point correlation function along the parallel direction ----------- 
  
  m_dataset = Projected(TwoPointCorrelation2D_cartesian::xx(), TwoPointCorrelation2D_cartesian::yy(), TwoPointCorrelation2D_cartesian::xi2D(), TwoPointCorrelation2D_cartesian::error2D());
}


// ============================================================================================


void cbl::measure::twopt::TwoPointCorrelation_projected::measureJackknife (const std::string dir_output_pairs, const std::vector<std::string> dir_input_pairs, const std::string dir_output_resample, const bool count_dd, const bool count_rr, const bool count_dr, const bool tcount, const Estimator estimator, const double fact)
{
  if (dir_output_resample!=par::defaultString && dir_output_resample!="") {
    string mkdir = "mkdir -p "+dir_output_resample;
    if (system(mkdir.c_str())) {}
  }

  vector<shared_ptr<data::Data>> data;
  vector<shared_ptr<pairs::Pair>> dd_regions, rr_regions, dr_regions;
  count_allPairs_region(dd_regions, rr_regions, dr_regions, TwoPType::_2D_Cartesian_, dir_output_pairs,dir_input_pairs, count_dd, count_rr, count_dr, tcount, estimator, fact);

  auto data_cart = (estimator==Estimator::_natural_) ? correlation_NaturalEstimator(m_dd, m_rr) : correlation_LandySzalayEstimator(m_dd, m_rr, m_dr);
  
  if (estimator==Estimator::_natural_) 
    data = XiJackknife(dd_regions, rr_regions);
  else if (estimator==Estimator::_LandySzalay_)
    data = XiJackknife(dd_regions, rr_regions, dr_regions);
  else
    ErrorCBL("the chosen estimator is not implemented!", "measureJackknife", "TwoPointCorrelation_projected.cpp");
  
  vector<vector<double>> ww, covariance;
  for (size_t i=0; i<data.size(); i++) {
    ww.push_back(data[i]->data());
    if (dir_output_resample != par::defaultString && dir_output_resample != "") {
      string file = "xi_projected_Jackkknife_"+conv(i, par::fINT)+".dat";
      string header = "[1] perpendicular separation at the bin centre # [2] projected two-point correlation function # [3] error";
      if (m_compute_extra_info) header += " # [4] mean perpendicular separation # [5] standard deviation of the distribution of perpendicular separations # [6] mean redshift # [7] standard deviation of the redshift distribution";
      data[i]->write(dir_output_resample, file, header, 10, 0);
    }
  }
  
  covariance_matrix(ww, covariance, true);

  vector<double> xx_cart = data_cart->xx(), yy_cart = data_cart->yy();
  vector<vector<double> > dd_cart, error_cart;
  data_cart->get_data(dd_cart); data_cart->get_error(error_cart);

  m_dataset = Projected(xx_cart, yy_cart, dd_cart, error_cart);
  m_dataset->set_covariance(covariance);

}


// ============================================================================================


void cbl::measure::twopt::TwoPointCorrelation_projected::measureBootstrap (const int nMocks, const std::string dir_output_pairs, const std::vector<std::string> dir_input_pairs, const std::string dir_output_resample, const bool count_dd, const bool count_rr, const bool count_dr, const bool tcount, const Estimator estimator, const double fact, const int seed)
{
  if (nMocks<=0)
    ErrorCBL("the number of mocks must be >0!", "measureBootstrap", "TwoPointCorrelation1D_monopole.cpp");

  if (dir_output_resample!=par::defaultString && dir_output_resample!="") {
    string mkdir = "mkdir -p "+dir_output_resample;
    if (system(mkdir.c_str())) {}
  }
  
  vector<shared_ptr<data::Data>> data;
  vector<shared_ptr<pairs::Pair>> dd_regions, rr_regions, dr_regions;
  count_allPairs_region(dd_regions, rr_regions, dr_regions, TwoPType::_2D_Cartesian_, dir_output_pairs,dir_input_pairs, count_dd, count_rr, count_dr, tcount, estimator, fact);

  auto data_cart = (estimator==Estimator::_natural_) ? correlation_NaturalEstimator(m_dd, m_rr) : correlation_LandySzalayEstimator(m_dd, m_rr, m_dr);

  if (estimator==Estimator::_natural_)
    data = XiBootstrap(nMocks, dd_regions, rr_regions, seed);
  else if (estimator==Estimator::_LandySzalay_)
    data = XiBootstrap(nMocks, dd_regions, rr_regions, dr_regions, seed);
  else
    ErrorCBL("the chosen estimator is not implemented!", "measureBootstrap", "TwoPointCorrelation_projected.cpp");
  
  vector<vector<double>> ww, covariance;
  for (size_t i=0; i<data.size(); i++) {
    ww.push_back(data[i]->data());
    if (dir_output_resample!=par::defaultString && dir_output_resample!="") {
      string file = "xi_projected_Bootstrap_"+conv(i, par::fINT)+".dat";
      string header = "[1] perpendicular separation at the bin centre # [2] projected two-point correlation function # [3] error";
      if (m_compute_extra_info) header += " # [4] mean perpendicular separation # [5] standard deviation of the distribution of perpendicular separations # [6] mean redshift # [7] standard deviation of the redshift distribution";
      data[i]->write(dir_output_resample, file, header, 10, 0);
    }
  }
  
  covariance_matrix(ww, covariance, false);

  vector<double> xx_cart = data_cart->xx(), yy_cart = data_cart->yy();
  vector<vector<double> > dd_cart, error_cart;
  data_cart->get_data(dd_cart); data_cart->get_error(error_cart);

  m_dataset = Projected(xx_cart, yy_cart, dd_cart, error_cart);
  m_dataset->set_covariance(covariance);
}


// ============================================================================================


std::vector<std::shared_ptr<data::Data>> cbl::measure::twopt::TwoPointCorrelation_projected::XiJackknife (const std::vector<std::shared_ptr<pairs::Pair>> dd, const std::vector<std::shared_ptr<pairs::Pair>> rr)
{
  vector<shared_ptr<data::Data>> data;
  
  auto data2d = TwoPointCorrelation2D_cartesian::XiJackknife(dd, rr);

  for (size_t i=0; i<data2d.size(); i++) {
    vector<double> xx_cart = data2d[i]->xx(), yy_cart = data2d[i]->yy();
    vector<vector<double>> dd_cart, error_cart;
    data2d[i]->get_data(dd_cart); data2d[i]->get_error(error_cart);
    data.push_back(move(Projected(xx_cart, yy_cart, dd_cart, error_cart)));
  }

  return data;
}


// ============================================================================================


std::vector<std::shared_ptr<data::Data>> cbl::measure::twopt::TwoPointCorrelation_projected::XiJackknife (const std::vector<std::shared_ptr<pairs::Pair>> dd, const std::vector<std::shared_ptr<pairs::Pair>> rr, const std::vector<std::shared_ptr<pairs::Pair>> dr)
{
  vector<shared_ptr<data::Data>> data;
 
  auto data2d = TwoPointCorrelation2D_cartesian::XiJackknife(dd, rr, dr);
  
  for (size_t i=0; i<data2d.size(); i++) {
    vector<double> xx_cart = data2d[i]->xx(), yy_cart = data2d[i]->yy();
    vector<vector<double>> dd_cart, error_cart;
    data2d[i]->get_data(dd_cart); data2d[i]->get_error(error_cart);
    data.push_back(move(Projected(xx_cart, yy_cart, dd_cart, error_cart)));
  }

  return data;
}


// ============================================================================================


std::vector<std::shared_ptr<data::Data>> cbl::measure::twopt::TwoPointCorrelation_projected::XiBootstrap (const int nMocks, const std::vector<std::shared_ptr<pairs::Pair>> dd, const std::vector<std::shared_ptr<pairs::Pair>> rr, const int seed)
{
  vector<shared_ptr<data::Data>> data;

  auto data2d = TwoPointCorrelation2D_cartesian::XiBootstrap(nMocks, dd, rr, seed);

  for (size_t i=0; i<data2d.size(); i++) {
    vector<double> xx_cart = data2d[i]->xx(), yy_cart = data2d[i]->yy();
    vector<vector<double>> dd_cart, error_cart;
    data2d[i]->get_data(dd_cart); data2d[i]->get_error(error_cart);
    data.push_back(move(Projected(xx_cart, yy_cart, dd_cart, error_cart)));
  }

  return data;
}


// ============================================================================================


std::vector<std::shared_ptr<data::Data>> cbl::measure::twopt::TwoPointCorrelation_projected::XiBootstrap (const int nMocks, const std::vector<std::shared_ptr<pairs::Pair>> dd, const std::vector<std::shared_ptr<pairs::Pair>> rr, const std::vector<std::shared_ptr<pairs::Pair>> dr, const int seed)
{
  vector<shared_ptr<data::Data>> data;

  auto data2d = TwoPointCorrelation2D_cartesian::XiBootstrap(nMocks, dd, rr, dr, seed);
  
  for (size_t i=0; i<data2d.size(); i++) {
    vector<double> xx_cart = data2d[i]->xx(), yy_cart = data2d[i]->yy();
    vector<vector<double>> dd_cart, error_cart;
    data2d[i]->get_data(dd_cart); data2d[i]->get_error(error_cart);
    data.push_back(move(Projected(xx_cart, yy_cart, dd_cart, error_cart)));
  }
  
  return data;
}


// ============================================================================


void cbl::measure::twopt::TwoPointCorrelation_projected::read_covariance (const std::string dir, const std::string file)
{
  m_dataset->set_covariance(dir+file);
}


// ============================================================================


void cbl::measure::twopt::TwoPointCorrelation_projected::write_covariance (const std::string dir, const std::string file) const
{
  m_dataset->write_covariance(dir, file);
}


// ============================================================================


void cbl::measure::twopt::TwoPointCorrelation_projected::compute_covariance (const std::vector<std::shared_ptr<data::Data>> xi, const bool JK)
{
  vector<vector<double>> Xi;

  for (size_t i=0; i<xi.size(); i++)
    Xi.push_back(xi[i]->data());

  vector<vector<double>> cov_mat;
  cbl::covariance_matrix(Xi, cov_mat, JK);
  
  m_dataset->set_covariance(cov_mat);
}


// ============================================================================


void cbl::measure::twopt::TwoPointCorrelation_projected::compute_covariance (const std::vector<std::string> file, const bool JK)
{
  vector<double> rad, mean;
  vector<vector<double>> cov_mat;

  cbl::covariance_matrix(file, rad, mean, cov_mat, JK); 
  m_dataset->set_covariance(cov_mat);
}
