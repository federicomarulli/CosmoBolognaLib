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
 *  CatalogueAnalysis/ThreePointCorrelation/ThreePointCorrelation_comoving_reduced.cpp
 *
 *  @brief Methods of the class ThreePointCorrelation_comoving_reduced
 *  used to measure the monopole of the three-point correlation
 *  function
 *
 *  This file contains the implementation of the methods of the class
 *  ThreePointCorrelation_comoving_reduced used to measure the
 *  monopole of the three-point correlation function
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */


#include "ThreePointCorrelation_comoving_reduced.h"
#include "TwoPointCorrelation1D_monopole.h"
#include "Data1D.h"

using namespace cosmobl;
using namespace catalogue;
using namespace measure;
using namespace triplets;
using namespace threept;


// ============================================================================


void cosmobl::measure::threept::ThreePointCorrelation_comoving_reduced::measure (const string dir_output_triplets, const string dir_output_2pt, const vector<string> dir_input_triplets, const bool count_ddd, const bool count_rrr, const bool count_ddr, const bool count_drr, const bool tcount, const int seed) 
{   
  (void)seed;
  
  // ----------- compute the connected three-point correlation function -----------
  
  ThreePointCorrelation_comoving_connected::measure(dir_output_triplets, dir_input_triplets, count_ddd, count_rrr, count_ddr, count_drr, tcount);
  m_scale = ThreePointCorrelation_comoving_connected::m_scale;
  
  
  // ----------- compute the two-point correlation function ----------- 

  Catalogue data = *m_data;
  Catalogue random = *m_random;
  double rMin = m_ddd->r12();
  double rMax = m_ddd->r12()+m_ddd->r13()+m_ddd->r13_binSize();
  double binSize = 0.05;
  double shift = 0.5;
  twopt::TwoPointCorrelation1D_monopole TwoP {data, random, _logarithmic_, rMin, rMax, binSize, shift};
  TwoP.measure(ErrorType::_Poisson_,dir_output_triplets, {}, par::defaultString,0, 1, 1, 1, tcount);

  vector<double> log_r(TwoP.dd()->nbins()), log_xi(TwoP.dd()->nbins());
  for (int i=0; i<TwoP.dd()->nbins(); i++) {
    log_r[i] = log10(TwoP.xx()[i]);
    log_xi[i] = log10(1.+TwoP.xi1D()[i]);
  }
  
  vector<double> values_interp(m_ddd->nbins()+2, 0.), xi_real_lin(m_ddd->nbins()+2, 0.);
  
  values_interp[0] = log10(m_ddd->r12());
  values_interp[1] = log10(m_ddd->r13());
  
  vector<double> theta(m_ddd->nbins(),0.);

  for (int i=0; i<m_ddd->nbins(); i++) {
    theta[i]=(i+0.5)*m_ddd->binSize();
    double tmp_value = sqrt(m_ddd->r12()*m_ddd->r12()+m_ddd->r13()*m_ddd->r13()-2.*m_ddd->r12()*m_ddd->r13()*cos(theta[i]));
    //double tmp_value = (m_ddd->side_s()+((i+0.5)*m_ddd->binSize()));
    values_interp[i+2] = log10(tmp_value);
  }

  string file_2pt = dir_output_2pt+"2ptCorrelation_3pt.dat";
  ofstream fout(file_2pt.c_str()); checkIO(fout, file_2pt);
 
  for (size_t i=0; i<values_interp.size(); i++) {
    xi_real_lin[i] = pow(10., interpolated(values_interp[i], log_r, log_xi, "Linear"))-1.;
    //coutCBL << pow(10.,values_interp[i]) << " --- " << xi_real_lin[i] << endl;
    fout << pow(10.,values_interp[i]) << "     " << xi_real_lin[i] << endl;
  }
  
  fout.clear(); fout.close(); coutCBL <<"I wrote the file "<<file_2pt<<endl;

  
  // ----------- compute the reduced three-point correlation function -----------

  m_QQ.resize(m_ddd->nbins()); m_error.resize(m_ddd->nbins());
  
  for (int i=0; i<m_ddd->nbins(); i++) {
    m_QQ[i] = m_zeta[i]/((xi_real_lin[0]*xi_real_lin[1])+(xi_real_lin[0]*xi_real_lin[i+2])+(xi_real_lin[1]*xi_real_lin[i+2]));
    m_error[i] = 0.001; // work in progress...
  }

  m_dataset = move(unique_ptr<data::Data1D>(new data::Data1D(theta, m_QQ, m_error)));
}


// ============================================================================


void cosmobl::measure::threept::ThreePointCorrelation_comoving_reduced::measure (const vector<vector<double>> weight, const bool doJK, const string dir_output_triplets, const string dir_output_2pt, const vector<string> dir_input_triplets, const bool count_ddd, const bool count_rrr, const bool count_ddr, const bool count_drr, const bool tcount, const int seed) 
{
  (void)seed;

  // ----------- compute the connected three-point correlation function -----------
  
  ThreePointCorrelation_comoving_connected::measure(weight, doJK, dir_output_triplets, dir_input_triplets, count_ddd, count_rrr, count_ddr, count_drr, tcount);
  m_scale = ThreePointCorrelation_comoving_connected::m_scale;
  
  
  // ----------- compute the two-point correlation function ----------- 

  Catalogue data = *m_data;
  Catalogue random = *m_random;
  double rMin = m_ddd->r12();
  double rMax = m_ddd->r12()+m_ddd->r13()+m_ddd->r13_binSize();
  double binSize = 0.05;
  double shift = 0.5;
  twopt::TwoPointCorrelation1D_monopole TwoP {data, random, _logarithmic_, rMin, rMax, binSize, shift};
  TwoP.measure(ErrorType::_Poisson_,dir_output_triplets, {}, par::defaultString,0, 1, 1, 1, tcount);

  vector<double> log_r(TwoP.dd()->nbins()), log_xi(TwoP.dd()->nbins());
  for (int i=0; i<TwoP.dd()->nbins(); i++) {
    log_r[i] = log10(TwoP.xx()[i]);
    log_xi[i] = log10(1.+TwoP.xi1D()[i]);
  }
  
  vector<double> values_interp(m_ddd->nbins()+2, 0.), xi_real_lin(m_ddd->nbins()+2, 0.);

  values_interp[0] = log10(m_ddd->r12());
  values_interp[1] = log10(m_ddd->r13());
  
  vector<double> theta(m_ddd->nbins(),0.);

  for (int i=0; i<m_ddd->nbins(); i++) {
    theta[i]=(i+0.5)*m_ddd->binSize();
    double tmp_value = sqrt(m_ddd->r12()*m_ddd->r12()+m_ddd->r13()*m_ddd->r13()-2.*m_ddd->r12()*m_ddd->r13()*cos(theta[i]));
    //double tmp_value = (m_ddd->side_s()+((i+0.5)*m_ddd->binSize()));
    values_interp[i+2] = log10(tmp_value);
  }

  string file_2pt = dir_output_2pt+"2ptCorrelation_3pt.dat";
  ofstream fout(file_2pt.c_str()); checkIO(fout, file_2pt);
 
  for (size_t i=0; i<values_interp.size(); i++) {
    xi_real_lin[i] = pow(10., interpolated(values_interp[i], log_r, log_xi, "Linear"))-1.;
    //coutCBL << pow(10.,values_interp[i]) << " --- " << xi_real_lin[i] << endl;
    fout << pow(10.,values_interp[i]) << "     " << xi_real_lin[i] << endl;
  }
  
  fout.clear(); fout.close(); coutCBL <<"I wrote the file "<<file_2pt<<endl;

  
  // ----------- compute the reduced three-point correlation function -----------

  m_QQ.resize(m_ddd->nbins()); m_error.resize(m_ddd->nbins());
  
  for (int i=0; i<m_ddd->nbins(); i++) {
    m_QQ[i] = m_zeta[i]/((xi_real_lin[0]*xi_real_lin[1])+(xi_real_lin[0]*xi_real_lin[i+2])+(xi_real_lin[1]*xi_real_lin[i+2]));
    m_error[i] = 0.001; // work in progress...
  }

  /// Compute resamplings and covariance matrix
  
  vector<vector<double>> resampling_threept(weight.size(), vector<double>(m_ddd->nbins(), 0));
  vector<long> region_list = m_data->region_list();

  vector<double> nData_reg_weighted, nRandom_reg_weighted;

  const int nRegions = m_data->nRegions(); 

  for (int i=0; i<nRegions; i++) {
    nData_reg_weighted.push_back(m_data->weightedN_condition(Var::_Region_, region_list[i], region_list[i]+1, 0));
    nRandom_reg_weighted.push_back(m_random->weightedN_condition(Var::_Region_, region_list[i], region_list[i]+1, 0));
  }

  for (size_t i=0; i<weight.size(); i++) {

    double nData = 0;
    double nRan = 0;

    for (size_t j=0; j<weight[i].size(); j++) {
      nData += weight[i][j]*nData_reg_weighted[j];
      nRan += weight[i][j]*nRandom_reg_weighted[j];
    }

    double norm1 = (double(nData)*double(nData-1)*double(nData-2))/6.;
    double norm2 = (double(nData)*double(nData-1)*double(nRan))*0.5;
    double norm3 = (double(nData)*double(nRan)*double(nRan-1))*0.5;
    double norm4 = (double(nRan)*double(nRan-1)*double(nRan-2))/6.;

    for (int j=0; j<m_ddd->nbins(); j++) 
      if (m_ddd_regions[i]->TT1D(j)>0 && m_rrr_regions[i]->TT1D(j)>0) {
	resampling_threept[i][j] = ((m_ddd_regions[i]->TT1D(j)/norm1)/(m_rrr_regions[i]->TT1D(j)/norm4))-3.*((m_ddr_regions[i]->TT1D(j)/norm2)/(m_rrr_regions[i]->TT1D(j)/norm4))+3.*(((m_drr_regions[i]->TT1D(j)/norm3)/(m_rrr_regions[i]->TT1D(j)/norm4)))-1.;
	resampling_threept[i][j] /= ((xi_real_lin[0]*xi_real_lin[1])+(xi_real_lin[0]*xi_real_lin[j+2])+(xi_real_lin[1]*xi_real_lin[j+2]));
      }
  }

  vector<vector<double>> cov_mat;
  cosmobl::covariance_matrix(resampling_threept, cov_mat, doJK);
  m_dataset = move(unique_ptr<data::Data1D>(new data::Data1D(m_scale, m_zeta, cov_mat)));

  m_dataset->error(m_error);
}


// ============================================================================================



void cosmobl::measure::threept::ThreePointCorrelation_comoving_reduced::measure (const cosmobl::measure::ErrorType errorType, const string dir_output_triplets, const string dir_output_2pt, const vector<string> dir_input_triplets, const int nResamplings, const bool count_ddd, const bool count_rrr, const bool count_ddr, const bool count_drr, const bool tcount, const int seed) 
{  

  switch (errorType) {
    
    case cosmobl::measure::ErrorType::_None_:
      {
	measure(dir_output_triplets, dir_output_2pt, dir_input_triplets, count_ddd, count_rrr, count_ddr, count_drr, tcount);
	break;
      }
    
    case cosmobl::measure::ErrorType::_Jackknife_:
      {
	const int nRegions = m_data->nRegions();

	vector<vector<double>> weight(nRegions, vector<double>(nRegions, 1));
	for (int i=0; i<nRegions; i++)
	  weight[i][i] = 0;

	measure(weight, true, dir_output_triplets, dir_output_2pt, dir_input_triplets, count_ddd, count_rrr, count_ddr, count_drr, tcount);
	break;
      }

    case cosmobl::measure::ErrorType::_Bootstrap_:
      {
	const int nRegions = m_data->nRegions();

	random::UniformRandomNumbers ran(0., nRegions-1, seed);
	
	int val = 3; // see Norberg et al. 2009

	vector<vector<double>> weight(nResamplings, vector<double>(nRegions, 0));
	for (int i=0; i<nResamplings; i++)
	  for (int j=0; j<val*nRegions; j++)
	    weight[i][ran()] ++;

	measure(weight, false, dir_output_triplets, dir_output_2pt, dir_input_triplets, count_ddd, count_rrr, count_ddr, count_drr, tcount);
	break;
      }

    default:
      ErrorCBL("Error in measure() of ThreePointCorrelation_comoving_reduced, no such kind of error type!");
  }

}


// ============================================================================


void cosmobl::measure::threept::ThreePointCorrelation_comoving_reduced::write (const string dir, const string file, bool connected) const
{      
  checkDim(m_scale, m_ddd->nbins(), "scale");
  
  string file_out = dir+file;
  ofstream fout(file_out.c_str()); checkIO(fout, file_out);

  if (!connected) {  
    fout << "# scale  Q  error(work in progress)" << endl;
    for (size_t i=0; i<m_scale.size(); i++) 
      fout << setiosflags(ios::fixed) << setprecision(4) << setw(8) << m_scale[i] << "  " << setw(8) << m_QQ[i] << "  " << setw(8) << m_error[i] << endl;
  }
  
  else {
    fout << "# scale  z  error(work in progrss)  Q  error(work in progress)" << endl;
    for (size_t i=0; i<m_scale.size(); i++) 
      fout << setiosflags(ios::fixed) << setprecision(4) << setw(8) << m_scale[i] << "  " << setw(8) << m_zeta[i] << "  " << setw(8) << ThreePointCorrelation_comoving_connected::m_error[i] << setw(8) << m_QQ[i] << "  " << setw(8) << m_error[i] << endl;
  }
  
  fout.close(); coutCBL << endl << "I wrote the file: " << file_out << endl << endl;
}  


// ============================================================================


void cosmobl::measure::threept::ThreePointCorrelation_comoving_reduced::write_covariance (const string dir, const string file) const
{
  m_dataset->write_covariance(dir, file);
}
