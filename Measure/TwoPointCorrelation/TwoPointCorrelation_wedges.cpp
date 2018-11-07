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
 *  CatalogueAnalysis/TwoPointCorrelation/TwoPointCorrelation_wedges.cpp
 *
 *  @brief Methods of the class TwoPointCorrelation_multipoles used to
 *  measure the wedges of the two-point correlation function
 *
 *  This file contains the implementation of the methods of the class
 *  TwoPointCorrelation_wedges used to measure the wedges of the
 *  two-point correlation function
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */


#include "TwoPointCorrelation_wedges.h"
#include "Data1D_extra.h"

using namespace std;

using namespace cbl;
using namespace catalogue;
using namespace chainmesh;
using namespace pairs;
using namespace measure;
using namespace twopt;


// ============================================================================


std::shared_ptr<data::Data> cbl::measure::twopt::TwoPointCorrelation_wedges::data_with_extra_info (const std::vector<double> rad, const std::vector<double> wedges, const std::vector<double> error) const
{
  auto dd2D = cbl::measure::twopt::TwoPointCorrelation2D_polar::m_dd;
  
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
  
  return move(unique_ptr<data::Data1D_extra>(new data::Data1D_extra(rad, wedges, error, extra)));
}


// ============================================================================================


std::vector<double> cbl::measure::twopt::TwoPointCorrelation_wedges::xx () const
{
  vector<double> rad, xx; m_dataset->xx(rad);

  for (size_t i=0; i<xx.size()/2; i++)
    rad.push_back(xx[i]);

  return rad;
}


// ============================================================================================


std::vector<double> cbl::measure::twopt::TwoPointCorrelation_wedges::xiPerpendicular () const
{
  vector<double> vv; 
  m_dataset->data(vv);

  size_t sz = vv.size();

  vector<double> xiPerp;
  
  for (size_t i=0; i<sz/2; i++)
    xiPerp.push_back(vv[i]);

  return xiPerp;
}


// ============================================================================================


std::vector<double> cbl::measure::twopt::TwoPointCorrelation_wedges::errorPerpendicular () const
{
vector<double> vv; 
  m_dataset->error(vv);

  size_t sz = vv.size();

  vector<double> error_xiPerp;
  
  for (size_t i=0; i<sz/2; i++)
    error_xiPerp.push_back(vv[i]);

  return error_xiPerp;
}


// ============================================================================================


std::vector<double> cbl::measure::twopt::TwoPointCorrelation_wedges::xiParallel () const 
{
vector<double> vv; 
  m_dataset->data(vv);

  size_t sz = vv.size();
  
  vector<double> xiParallel;

  for (size_t i=sz/2; i<sz; i++)
    xiParallel.push_back(vv[i]);

  return xiParallel;
}


// ============================================================================================


std::vector<double> cbl::measure::twopt::TwoPointCorrelation_wedges::errorParallel () const 
{
  vector<double> vv; 
  m_dataset->error(vv);

  size_t sz = vv.size();

  vector<double> error_xiParallel;

  for (size_t i=sz/2; i<sz; i++)
    error_xiParallel.push_back(vv[i]);

  return error_xiParallel;
}


// ============================================================================================


void cbl::measure::twopt::TwoPointCorrelation_wedges::write (const std::string dir, const std::string file, const int rank) const 
{
  (void)rank;
  
  vector<double> rad; m_dataset->xx(rad);
  vector<double> xiw = m_dataset->data();
  vector<double> error = m_dataset->error();

  checkDim(rad, m_dd->nbins_D1()*2, "rad");
  
  string file_out = dir+file;
  ofstream fout(file_out.c_str()); checkIO(fout, file_out);

  string header = "[1] separation at the bin centre # [2] perpendicular wedge # [3] error on the perpendicular wedge # [4] parallel wedge # [5] error on the parallel wedge";
  if (m_compute_extra_info)
    header += " # [6] mean separation # [7] standard deviation of the separation distribution # [8] mean redshift # [9] standard deviation of the redshift distribution";
  
  fout << "### " << header << " ###" <<endl;

  for (int i=0; i<m_dd->nbins_D1(); i++) {
    fout << setiosflags(ios::fixed) << setprecision(5) << setw(10) << right << rad[i]
	 << "  " << setiosflags(ios::fixed) << setprecision(5) << setw(10) << right  << xiw[i]
	 << "  " << setiosflags(ios::fixed) << setprecision(5) << setw(10) << right  << error[i]
	 << "  " << setiosflags(ios::fixed) << setprecision(5) << setw(10) << right  << xiw[i+m_dd->nbins_D1()]
	 << "  " << setiosflags(ios::fixed) << setprecision(5) << setw(10) << right  << error[i+m_dd->nbins_D1()];
    if (m_compute_extra_info)
      for (size_t ex=0; ex<m_dataset->extra_info().size(); ++ex)
	fout << "  " << setiosflags(ios::fixed) << setprecision(5) << setw(10) << right  << m_dataset->extra_info(ex, i);
    fout << endl;
  }
  
  fout.close(); cout << endl; coutCBL << "I wrote the file: " << file_out << endl << endl;
}


// ============================================================================================


std::shared_ptr<data::Data> cbl::measure::twopt::TwoPointCorrelation_wedges::Wedges (const std::vector<double> rr, const std::vector<double> mu, const std::vector<std::vector<double> > xi, const std::vector<std::vector<double> > error_xi)
{
  double binSize = mu[1]-mu[0];
  int muSize = mu.size();

  int half = max(0, min(int(0.5/binSize), muSize));
  half = (mu[half]<0.5) ? half+1 : half;
  
  vector<double> rad(2*rr.size(),0), wedges(2*rr.size(),0), error(2*rr.size(),0);

  for (size_t i=0; i<rr.size(); i++) {
    rad[i] = rr[i];
    rad[i+rr.size()] = rr[i];

    for (int j=0; j<half; j++) {
      wedges[i] += 2.*xi[i][j]*binSize;   	                        // xi_perp
      if (wedges[i]>-1.) error[i] += 2.*pow(error_xi[i][j]*binSize, 2); // error[xi_perp]
    }

    for (size_t j=half; j<mu.size(); j++) {
      wedges[i+rr.size()] += 2.*xi[i][j]*binSize;                                           // xi_par
      if (wedges[i+rr.size()]>-1.) error[i+rr.size()] += 2.*pow(error_xi[i][j]*binSize, 2); // error[xi_par]
    }  
  }

  for_each( error.begin(), error.end(), [] (double &vv) { vv = sqrt(vv);} );

  return (!m_compute_extra_info) ? unique_ptr<data::Data1D>(new data::Data1D(rad, wedges, error)) : data_with_extra_info(rad, wedges, error);
}


// ============================================================================================


void cbl::measure::twopt::TwoPointCorrelation_wedges::measure(const ErrorType errorType, const std::string dir_output_pairs, const std::vector<std::string> dir_input_pairs, const std::string dir_output_resample, const int nMocks, const bool count_dd, const bool count_rr, const bool count_dr, const bool tcount, const Estimator estimator, const int seed)
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
      ErrorCBL("Error in measure() of TwoPointCorrelation_multipoles.cpp, unknown type of error");
  }
}


// ============================================================================================


void cbl::measure::twopt::TwoPointCorrelation_wedges::measurePoisson (const std::string dir_output_pairs, const std::vector<std::string> dir_input_pairs, const bool count_dd, const bool count_rr, const bool count_dr, const bool tcount, const Estimator estimator)
{
  // ----------- measure the 2D two-point correlation function, xi(rp,pi) ----------- 

  TwoPointCorrelation2D_polar::measurePoisson(dir_output_pairs, dir_input_pairs, count_dd, count_rr, count_dr, tcount, estimator);

  
  // ----------- integrate the 2D two-point correlation function along the parallel direction ----------- 

  m_dataset = Wedges(TwoPointCorrelation2D_polar::xx(), TwoPointCorrelation2D_polar::yy(), TwoPointCorrelation2D_polar::xi2D(), TwoPointCorrelation2D_polar::error2D());
}


// ============================================================================================


void cbl::measure::twopt::TwoPointCorrelation_wedges::measureJackknife (const std::string dir_output_pairs, const std::vector<std::string> dir_input_pairs, const std::string dir_output_resample, const bool count_dd, const bool count_rr, const bool count_dr, const bool tcount, const Estimator estimator)
{
  if (dir_output_resample != par::defaultString) {
    string mkdir = "mkdir -p "+dir_output_resample;
    if (system(mkdir.c_str())) {}
  }

  vector<shared_ptr<data::Data> > data;
  vector<shared_ptr<pairs::Pair> > dd_regions, rr_regions, dr_regions;
  count_allPairs_region (dd_regions, rr_regions, dr_regions, TwoPType::_2D_polar_, dir_output_pairs,dir_input_pairs, count_dd, count_rr, count_dr,  tcount, estimator);

  auto data_polar = (estimator==Estimator::_natural_) ? correlation_NaturalEstimator(m_dd, m_rr) : correlation_LandySzalayEstimator(m_dd, m_rr, m_dr);

  if (estimator==Estimator::_natural_) 
    data = XiJackknife(dd_regions, rr_regions);
  else if (estimator==Estimator::_LandySzalay_)
    data = XiJackknife(dd_regions, rr_regions, dr_regions);
  else
    ErrorCBL("Error in measureJackknife() of TwoPointCorrelation_wedges.cpp: the chosen estimator is not implemented!");
  
  vector<vector<double> > ww, covariance;
  
  for (size_t i=0; i<data.size(); i++) {
    ww.push_back(data[i]->data());
    
    if (dir_output_resample != par::defaultString) {
      string file_out = dir_output_resample+"xi_wedges_Jackknife_"+conv(i, par::fINT)+".dat";

      vector<double> rad; data[i]->xx(rad);
      vector<double> xiw = data[i]->data();
      vector<double> error = data[i]->error();

      checkDim(rad, m_dd->nbins_D1()*2, "rad");
  
      ofstream fout(file_out.c_str()); checkIO(fout, file_out);

      fout << "### [1] separation at the bin centre # [2] perpendicular wedge # [3] error on the perpendicular wedge # [4] parallel wedge # [5] error on the parallel wedge ###" << endl;

      for (int i=0; i<m_dd->nbins_D1(); i++) 
	fout << setiosflags(ios::fixed) << setprecision(5) << setw(10) << right << rad[i]
	     << "  " << setiosflags(ios::fixed) << setprecision(5) << setw(10) << right << xiw[i]
	     << "  " << setiosflags(ios::fixed) << setprecision(5) << setw(10) << right << error[i]
	     << "  " << setiosflags(ios::fixed) << setprecision(5) << setw(10) << right << xiw[i+m_dd->nbins_D1()]
	     << "  " << setiosflags(ios::fixed) << setprecision(5) << setw(10) << right << error[i+m_dd->nbins_D1()] << endl;
      fout.close(); cout << endl; coutCBL << "I wrote the file: " << file_out << endl << endl;
    }
  }

  covariance_matrix(ww, covariance, true);

  vector<double> xx_polar, yy_polar;
  vector<vector<double>> dd_polar, error_polar;
  data_polar->xx(xx_polar); data_polar->yy(yy_polar);
  data_polar->data(dd_polar); data_polar->error(error_polar);

  m_dataset = Wedges(xx_polar, yy_polar, dd_polar, error_polar);
  m_dataset->set_covariance(covariance);

}


// ============================================================================================


void cbl::measure::twopt::TwoPointCorrelation_wedges::measureBootstrap (const int nMocks, const std::string dir_output_pairs, const std::vector<std::string> dir_input_pairs, const std::string dir_output_resample, const bool count_dd, const bool count_rr, const bool count_dr, const bool tcount, const Estimator estimator, const int seed)
{
  if (dir_output_resample != par::defaultString) {
    string mkdir = "mkdir -p "+dir_output_resample;
    if (system(mkdir.c_str())) {}
  }
  
  vector<shared_ptr<data::Data> > data;
  vector<shared_ptr<pairs::Pair> > dd_regions, rr_regions, dr_regions;
  count_allPairs_region (dd_regions, rr_regions, dr_regions, TwoPType::_2D_polar_, dir_output_pairs,dir_input_pairs, count_dd, count_rr, count_dr,  tcount, estimator);

  auto data_polar = (estimator==Estimator::_natural_) ? correlation_NaturalEstimator(m_dd, m_rr) : correlation_LandySzalayEstimator(m_dd, m_rr, m_dr);

  if (estimator==Estimator::_natural_)
    data = XiBootstrap(nMocks, dd_regions, rr_regions, seed);
  else if (estimator==Estimator::_LandySzalay_)
    data = XiBootstrap(nMocks, dd_regions, rr_regions, dr_regions, seed);
  else
    ErrorCBL("Error in measureBootstrap() of TwoPointCorrelation_wedges.cpp: the chosen estimator is not implemented!");
  
  vector<vector<double> > ww, covariance;
  
  for (size_t i=0; i<data.size(); i++) {
    ww.push_back(data[i]->data());
    
    if (dir_output_resample != par::defaultString) {
      string file_out = dir_output_resample+"xi_wedges_Bootstrap_"+conv(i, par::fINT)+".dat";

      vector<double> rad; data[i]->xx(rad);
      vector<double> xiw = data[i]->data();
      vector<double> error = data[i]->error();

      checkDim(rad, m_dd->nbins_D1()*2, "rad");
  
      ofstream fout(file_out.c_str()); checkIO(fout, file_out);

      fout << "### [1] separation at the bin centre # [2] perpendicular wedge # [3] error on the perpendicular wedge # [4] parallel wedge # [5] error on the parallel wedge ###" << endl;

      for (int i=0; i<m_dd->nbins_D1(); i++) 
	fout << setiosflags(ios::fixed) << setprecision(5) << setw(10) << right << rad[i]
	     << "  " << setiosflags(ios::fixed) << setprecision(5) << setw(10) << right << xiw[i]
	     << "  " << setiosflags(ios::fixed) << setprecision(5) << setw(10) << right << error[i]
	     << "  " << setiosflags(ios::fixed) << setprecision(5) << setw(10) << right << xiw[i+m_dd->nbins_D1()]
	     << "  " << setiosflags(ios::fixed) << setprecision(5) << setw(10) << right << error[i+m_dd->nbins_D1()] << endl;
   
      fout.close(); cout << endl; coutCBL << "I wrote the file: " << file_out << endl << endl;
    }
  }
  
  covariance_matrix(ww, covariance, false);

  vector<double> xx_polar, yy_polar;
  vector<vector<double>> dd_polar, error_polar;
  data_polar->xx(xx_polar); data_polar->yy(yy_polar);
  data_polar->data(dd_polar); data_polar->error(error_polar);

  m_dataset = Wedges(xx_polar, yy_polar, dd_polar, error_polar);
  m_dataset->set_covariance(covariance);
}


// ============================================================================================


std::vector<std::shared_ptr<data::Data> > cbl::measure::twopt::TwoPointCorrelation_wedges::XiJackknife (const std::vector<std::shared_ptr<pairs::Pair> > dd, const std::vector<std::shared_ptr<pairs::Pair> > rr)
{
  vector<shared_ptr<data::Data> > data;

  auto data2d = TwoPointCorrelation2D_polar::XiJackknife(dd, rr);

  for (size_t i=0; i<data2d.size(); i++){
    vector<double> xx_polar, yy_polar;
    vector<vector<double>> dd_polar, error_polar;
    data2d[i]->xx(xx_polar); data2d[i]->yy(yy_polar);
    data2d[i]->data(dd_polar); data2d[i]->error(error_polar);
    data.push_back(move(Wedges(xx_polar, yy_polar, dd_polar, error_polar)));
  }
  
  return data;
}


// ============================================================================================


std::vector<std::shared_ptr<data::Data> > cbl::measure::twopt::TwoPointCorrelation_wedges::XiJackknife (const std::vector<std::shared_ptr<pairs::Pair> > dd, const std::vector<std::shared_ptr<pairs::Pair> > rr, const std::vector<std::shared_ptr<pairs::Pair> > dr)
{
  vector<shared_ptr<data::Data> > data;
  
  auto data2d = TwoPointCorrelation2D_polar::XiJackknife(dd, rr, dr);

  for (size_t i=0; i<data2d.size(); i++){
    vector<double> xx_polar, yy_polar;
    vector<vector<double>> dd_polar, error_polar;
    data2d[i]->xx(xx_polar); data2d[i]->yy(yy_polar);
    data2d[i]->data(dd_polar); data2d[i]->error(error_polar);
    data.push_back(move(Wedges(xx_polar, yy_polar, dd_polar, error_polar)));
  }
  
  return data;
}


// ============================================================================================


std::vector<std::shared_ptr<data::Data> > cbl::measure::twopt::TwoPointCorrelation_wedges::XiBootstrap (const int nMocks, const std::vector<std::shared_ptr<pairs::Pair> > dd, const std::vector<std::shared_ptr<pairs::Pair> > rr, const int seed)
{
  vector<shared_ptr<data::Data> > data;
  
  auto data2d = TwoPointCorrelation2D_polar::XiBootstrap(nMocks, dd, rr, seed);

  for (size_t i=0; i<data2d.size(); i++){
    vector<double> xx_polar, yy_polar;
    vector<vector<double>> dd_polar, error_polar;
    data2d[i]->xx(xx_polar); data2d[i]->yy(yy_polar);
    data2d[i]->data(dd_polar); data2d[i]->error(error_polar);
    data.push_back(move(Wedges(xx_polar, yy_polar, dd_polar, error_polar)));
  }
  
  return data;
}


// ============================================================================================


std::vector<std::shared_ptr<data::Data> > cbl::measure::twopt::TwoPointCorrelation_wedges::XiBootstrap (const int nMocks, const std::vector<std::shared_ptr<pairs::Pair> > dd, const std::vector<std::shared_ptr<pairs::Pair> > rr, const std::vector<std::shared_ptr<pairs::Pair> > dr, const int seed)
{
  vector<shared_ptr<data::Data> > data;
  
  auto data2d = TwoPointCorrelation2D_polar::XiBootstrap(nMocks, dd, rr, dr, seed);

  for (size_t i=0; i<data2d.size(); i++){
    vector<double> xx_polar, yy_polar;
    vector<vector<double>> dd_polar, error_polar;
    data2d[i]->xx(xx_polar); data2d[i]->yy(yy_polar);
    data2d[i]->data(dd_polar); data2d[i]->error(error_polar);
    data.push_back(move(Wedges(xx_polar, yy_polar, dd_polar, error_polar)));
  }
  
  return data;
}


// ============================================================================


void cbl::measure::twopt::TwoPointCorrelation_wedges::read_covariance (const std::string dir, const std::string file)
{
  m_dataset->set_covariance(dir+file);
}


// ============================================================================


void cbl::measure::twopt::TwoPointCorrelation_wedges::write_covariance (const std::string dir, const std::string file) const
{
  m_dataset->write_covariance(dir, file);
}


// ============================================================================


void cbl::measure::twopt::TwoPointCorrelation_wedges::compute_covariance (const std::vector<std::shared_ptr<data::Data>> xi, const bool JK)
{
  vector<vector<double>> Xi;

  for (size_t i=0; i<xi.size(); i++)
    Xi.push_back(xi[i]->data());

  vector<vector<double>> cov_mat;
  cbl::covariance_matrix(Xi, cov_mat, JK);
  
  m_dataset->set_covariance(cov_mat);
}


// ============================================================================


void cbl::measure::twopt::TwoPointCorrelation_wedges::compute_covariance (const std::vector<std::string> file, const bool JK)
{
  vector<double> rad, mean;
  vector<vector<double>> cov_mat;

  cbl::covariance_matrix(file, rad, mean, cov_mat, JK); 
  m_dataset->set_covariance(cov_mat);
}
