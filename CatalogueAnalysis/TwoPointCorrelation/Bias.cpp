/*******************************************************************
 *  Copyright (C) 2010 by Federico Marulli                         *
 *  federico.marulli3@unibo.it                                     *
 *                                                                 *
 *  This program is free software; you can redistribute it and/or  *
 *  modify it under the terms of the GNU General Public License as *
 *  published by the Free Software Foundation; either version 2 of *
 *  the License, or (at your option) any later version.            *
 *                                                                 *
 *  This program is distributed in the hope that it will be useful,*
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the  *
 *  GNU General Public License for more details.                   *
 *                                                                 *
 *  You should have received a copy of the GNU General Public      *
 *  License along with this program; if not, write to the Free     *
 *  Software Foundation, Inc.,                                     *
 *  59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.      *
 *******************************************************************/

/**
 *  @file CatalogueAnalysis/TwoPointCorrelation/Bias.cpp
 *
 *  @brief Methods of the class TwoPointCorrelation used to measure
 *  the bias
 *
 *  This file contains the implementation of the methods of the class
 *  TwoPointCorrelation used to estimate the bias function
 *
 *  @author Federico Marulli
 *
 *  @author federico.marulli3@unbo.it
 */

#include "TwoPointCorrelation.h"
using namespace cosmobl;


// ============================================================================


double cosmobl::TwoPointCorrelation::sigmaR_obj (double RR, int corrType) 
{
  
  vector<double> rrr, xxx;

  if (corrType==1) { // using the measured correlation function
    for (unsigned int i=0; i<m_rad_lin.size(); i++)
      if (m_rad_lin[i]>0 && m_xi_lin[i]>0) {
	rrr.push_back(m_rad_lin[i]);
	xxx.push_back(m_xi_lin[i]);
      }
  }
  
  if (corrType==2) { // using the projected correlation function
    for (unsigned int i=0; i<m_rp_proj.size(); i++)
      if (m_rp_proj[i]>0 && m_xi_proj[i]>0) {
	rrr.push_back(m_rp_proj[i]);
	xxx.push_back(m_xi_proj[i]);
      }
  }

  if (corrType==3) { // using the deprojected correlation function
    for (unsigned int i=0; i<m_xi_real_lin.size(); i++)
      if (m_rad_lin[i]>0 && m_xi_real_lin[i]>0) {
	rrr.push_back(m_rad_lin[i]);
	xxx.push_back(m_xi_real_lin[i]);
      }
    corrType = 1;
  }
  
  if (rrr.size()==0) ErrorMsg("Error in cosmobl::TwoPointCorrelation::sigmaR_DM of Bias.cpp!");
  
  return sigmaR (RR, corrType, rrr, xxx);
}


// ============================================================================


void cosmobl::TwoPointCorrelation::measure_bias (Cosmology &cosm, double &redshift, string &aut, string &Model, string file_bias, bool proj) 
{
  if (proj==0 && m_rad_log.size()==0) 
    ErrorMsg("Error in cosmobl::TwoPointCorrelation::measure_bias of Bias.cpp! rad_log.size=0");
  if (proj==1 && m_rp_proj.size()==0)
    ErrorMsg("Error in cosmobl::TwoPointCorrelation::measure_bias of Bias.cpp! rp_proj.size=0");
 
  string author = "CAMB";
  double s8z = cosm.sigma8_Pk(author, redshift, Model);

  ofstream fout;  
  if (file_bias!="NULL") {fout.open(file_bias.c_str()); checkIO(file_bias,0);}
 
  if (proj) {
    m_rr_bias_lin_wp.erase (m_rr_bias_lin_wp.begin(), m_rr_bias_lin_wp.end()); 
    m_bias_lin_wp.erase (m_bias_lin_wp.begin(), m_bias_lin_wp.end()); 
    m_error_bias_lin_wp.erase (m_error_bias_lin_wp.begin(), m_error_bias_lin_wp.end()); 
    m_rr_bias_nl_wp.erase (m_rr_bias_nl_wp.begin(), m_rr_bias_nl_wp.end()); 
    m_bias_nl_wp.erase (m_bias_nl_wp.begin(), m_bias_nl_wp.end()); 
    m_error_bias_nl_wp.erase (m_error_bias_nl_wp.begin(), m_error_bias_nl_wp.end());

    for (unsigned int i=0; i<m_rp_proj.size(); i++) {

      if (m_rp_proj[i]>0) {
	double wpDM_lin = cosm.wp_DM(m_rp_proj[i],aut,redshift,Model,0);
	double wpDM_nl = cosm.wp_DM(m_rp_proj[i],aut,redshift,Model);
	
	if (m_xi_proj[i]/wpDM_lin>0) {
	  cout <<"rp = "<<m_rp_proj[i]<<": b_lin = "<<sqrt(m_xi_proj[i]/wpDM_lin)<<endl;
	  m_rr_bias_lin_wp.push_back(m_rp_proj[i]);
	  m_bias_lin_wp.push_back(sqrt(m_xi_proj[i]/wpDM_lin));
	  m_error_bias_lin_wp.push_back(0.5*sqrt(wpDM_lin/m_xi_proj[i])*m_error_xi_proj[i]/wpDM_lin);
	}

	if (m_xi_proj[i]/wpDM_nl>0) {
	  m_rr_bias_nl_wp.push_back(m_rp_proj[i]);
	  m_bias_nl_wp.push_back(sqrt(m_xi_proj[i]/wpDM_nl));
	  m_error_bias_nl_wp.push_back(0.5*sqrt(wpDM_nl/m_xi_proj[i])*m_error_xi_proj[i]/wpDM_nl);
	}
      }

      else {
	m_rr_bias_lin_wp.push_back(-1.);
	m_bias_lin_wp.push_back(-1.);
	m_error_bias_lin_wp.push_back(-1.);
	m_rr_bias_nl_wp.push_back(-1.);
	m_bias_nl_wp.push_back(-1.);
	m_error_bias_nl_wp.push_back(-1.);
      }

    }

    if (file_bias!="NULL") {

      if (m_rr_bias_lin_wp.size()!=m_rr_bias_nl_wp.size()) {
	string Err = "Error in cosmobl::TwoPointCorrelation::measure_bias of Bias.cpp: " + conv(m_rr_bias_lin_wp.size(),par::fINT) + " != " + conv(m_rr_bias_nl_wp.size(),par::fINT);	
	ErrorMsg(Err);
      }

      for (unsigned int bb=0; bb<m_rr_bias_lin_wp.size(); bb++) 
	if (m_rr_bias_lin_wp[bb]>0)
	  fout <<m_rr_bias_lin_wp[bb]<<"   "<<m_bias_lin_wp[bb]<<"   "<<m_error_bias_lin_wp[bb]<<"   "<<m_bias_nl_wp[bb]<<"   "<<m_error_bias_nl_wp[bb]<<"   "<<m_bias_lin_wp[bb]*s8z<<"   "<<m_error_bias_lin_wp[bb]*s8z<<"   "<<m_bias_nl_wp[bb]*s8z<<"   "<<m_error_bias_nl_wp[bb]*s8z<<endl;
    }
  }
  
  else {
    m_rr_bias_lin_xi.erase (m_rr_bias_lin_xi.begin(), m_rr_bias_lin_xi.end()); 
    m_bias_lin_xi.erase (m_bias_lin_xi.begin(), m_bias_lin_xi.end()); 
    m_error_bias_lin_xi.erase (m_error_bias_lin_xi.begin(), m_error_bias_lin_xi.end()); 
    m_rr_bias_nl_xi.erase (m_rr_bias_nl_xi.begin(), m_rr_bias_nl_xi.end()); 
    m_bias_nl_xi.erase (m_bias_nl_xi.begin(), m_bias_nl_xi.end()); 
    m_error_bias_nl_xi.erase (m_error_bias_nl_xi.begin(), m_error_bias_nl_xi.end());

    for (unsigned int i=0; i<m_rad_log.size(); i++) {
      if (m_rad_log[i]>0) {
	double xiDM_lin = cosm.xi_DM(m_rad_log[i],aut,redshift,Model,0);
	double xiDM_nl = cosm.xi_DM(m_rad_log[i],aut,redshift,Model);

	if (m_xi_log[i]/xiDM_lin>0) {
	  m_rr_bias_lin_xi.push_back(m_rad_log[i]);
	  m_bias_lin_xi.push_back(sqrt(m_xi_log[i]/xiDM_lin));
	  m_error_bias_lin_xi.push_back(0.5*sqrt(xiDM_lin/m_xi_log[i])*m_error_xi_log[i]/xiDM_lin);
	}
	if (xiDM_nl/m_xi_log[i]>0) {
	  m_rr_bias_nl_xi.push_back(m_rad_log[i]);
	  m_bias_nl_xi.push_back(sqrt(m_xi_log[i]/xiDM_nl));
	  m_error_bias_nl_xi.push_back(0.5*sqrt(xiDM_nl/m_xi_log[i])*m_error_xi_log[i]/xiDM_nl);
	}
      }

      else {
	m_rr_bias_lin_xi.push_back(-1.);
	m_bias_lin_xi.push_back(-1.);
	m_error_bias_lin_xi.push_back(-1.);
	m_rr_bias_nl_xi.push_back(-1.);
	m_bias_nl_xi.push_back(-1.);
	m_error_bias_nl_xi.push_back(-1.);
      }

    }

    if (file_bias!="NULL") {
      
      if (m_rr_bias_lin_xi.size()!=m_rr_bias_nl_xi.size()) {
	string Err = "Error in cosmobl::TwoPointCorrelation::measure_bias of Bias.cpp: " + conv(m_rr_bias_lin_xi.size(),par::fINT) + " != " + conv(m_rr_bias_nl_xi.size(),par::fINT);	
	ErrorMsg(Err);
      }
      
      for (unsigned int bb=0; bb<m_rr_bias_lin_xi.size(); bb++) 
	if (m_rr_bias_lin_xi[bb]>0)
	  fout <<m_rr_bias_lin_xi[bb]<<"   "<<m_bias_lin_xi[bb]<<"   "<<m_error_bias_lin_xi[bb]<<"   "<<m_bias_nl_xi[bb]<<"   "<<m_error_bias_nl_xi[bb]<<"   "<<m_bias_lin_xi[bb]*s8z<<"   "<<m_error_bias_lin_xi[bb]*s8z<<"   "<<m_bias_nl_xi[bb]*s8z<<"   "<<m_error_bias_nl_xi[bb]*s8z<<endl;
    }
    
  }

  if (file_bias!="NULL") {fout.clear(); fout.close(); cout <<"I wrote the file: "<<file_bias<<endl;}
}


// ============================================================================


double cosmobl::TwoPointCorrelation::mean_bias (double &rMin, double &rMax, bool proj, bool NL) 
{
  if (proj==0 && m_bias_lin_xi.size()==0) 
    ErrorMsg("Error in cosmobl::TwoPointCorrelation::mean_bias of Bias.cpp! bias_lin_xi.size=0");
  if (proj==1 && m_bias_lin_wp.size()==0)
    ErrorMsg("Error in cosmobl::TwoPointCorrelation::mean_bias of Bias.cpp! bias_lin_wp.size=0");

  vector<double> biasv, inv_error;

  if (proj) {
    int dim = (NL) ? m_bias_nl_wp.size() : m_bias_lin_wp.size();
    for (int i=0; i<dim; i++) {
      double rr = (NL) ? m_rr_bias_nl_wp[i] : m_rr_bias_lin_wp[i];
      if (rMin<rr && rr<rMax) {	
	double bb = (NL) ? m_bias_nl_wp[i] : m_bias_lin_wp[i];
	double err = (NL) ? m_error_bias_nl_wp[i] : m_error_bias_lin_wp[i];
	biasv.push_back(bb);
	inv_error.push_back(1./(err*err));
      }
    }
  }
  else {
    int dim = (NL) ? m_bias_nl_xi.size() : m_bias_lin_xi.size();
    for (int i=0; i<dim; i++) {
      double rr = (NL) ? m_rr_bias_nl_xi[i] : m_rr_bias_lin_xi[i];
      if (rMin<rr && rr<rMax) {	
	double bb = (NL) ? m_bias_nl_xi[i] : m_bias_lin_xi[i];
	double err = (NL) ? m_error_bias_nl_xi[i] : m_error_bias_lin_xi[i];
	biasv.push_back(bb);
	inv_error.push_back(1./(err*err));
      }
    }
  } 
  
  m_bias = Average(biasv,inv_error);

  return m_bias;
}


// ============================================================================


double cosmobl::TwoPointCorrelation::error_mean_bias (double &rMin, double &rMax, bool proj, bool NL) 
{
  if (proj==0 && m_bias_lin_xi.size()==0) 
    ErrorMsg("Error in cosmobl::TwoPointCorrelation::error_mean_bias of Bias.cpp! rad_log.size=0");
  if (proj==1 && m_bias_lin_wp.size()==0)
    ErrorMsg("Error in cosmobl::TwoPointCorrelation::error_mean_bias of Bias.cpp! rp_proj.size=0");

  vector<double> error_bias;

  if (proj) {
    for (unsigned int i=0; i<m_rp_proj.size(); i++) {
      double rr = (NL) ? m_rr_bias_nl_wp[i] : m_rr_bias_lin_wp[i];
      if (rMin<rr && rr<rMax) {	
	double err = (NL) ? m_error_bias_nl_wp[i] : m_error_bias_lin_wp[i];
	error_bias.push_back(err);
      }
    }
  }
  else {
    for (unsigned int i=0; i<m_rad_log.size(); i++) {
      double rr = (NL) ? m_rr_bias_nl_xi[i] : m_rr_bias_lin_xi[i];
      if (rMin<rr && rr<rMax) {	
	double err = (NL) ? m_error_bias_nl_xi[i] : m_error_bias_lin_xi[i];
	error_bias.push_back(err);
      }
    }
  } 

  m_bias_min = m_bias-Average(error_bias); // check!!!
  m_bias_max = m_bias+Average(error_bias); // check!!!

  return Average(error_bias);
}


// ============================================================================


void cosmobl::TwoPointCorrelation::fit_bias (double &rMin, double &rMax, bool proj, bool NL) 
{
  if (proj==0 && m_bias_lin_xi.size()==0) 
    ErrorMsg("Error in cosmobl::TwoPointCorrelation::fit_bias of Bias.cpp! bias_lin_xi.size=0");
  if (proj==1 && m_bias_lin_wp.size()==0)
    ErrorMsg("Error in cosmobl::TwoPointCorrelation::fit_bias of Bias.cpp! bias_lin_wp.size=0");
  
  vector<double> rad, biasv, weight;

  if (proj) {
    int dim = (NL) ? m_bias_nl_wp.size() : m_bias_lin_wp.size();
    for (int i=0; i<dim; i++) {
      double rr = (NL) ? m_rr_bias_nl_wp[i] : m_rr_bias_lin_wp[i];
      if (rMin<rr && rr<rMax) {	
	double bb = (NL) ? m_bias_nl_wp[i] : m_bias_lin_wp[i];
	double err = (NL) ? m_error_bias_nl_wp[i] : m_error_bias_lin_wp[i];
	rad.push_back(rr);
	biasv.push_back(bb);
	weight.push_back(1./(err*err));
      }
    }
  }
  else {
    int dim = (NL) ? m_bias_nl_xi.size() : m_bias_lin_xi.size();
    for (int i=0; i<dim; i++) {
      double rr = (NL) ? m_rr_bias_nl_xi[i] : m_rr_bias_lin_xi[i];
      if (rMin<rr && rr<rMax) {	
	double bb = (NL) ? m_bias_nl_xi[i] : m_bias_lin_xi[i];
	double err = (NL) ? m_error_bias_nl_xi[i] : m_error_bias_lin_xi[i];
	rad.push_back(rr);
	biasv.push_back(bb);
	weight.push_back(1./(err*err));
      }
    }
  } 

  double c0, c1, cov00, cov01, cov11, chisq;

  double *radG = new double[rad.size()];
  double *weightG = new double[rad.size()];
  double *biasG = new double[rad.size()];
 
  for (unsigned int i=0; i<rad.size(); i++) {
    radG[i] = rad[i]; 
    weightG[i] = weight[i]; 
    biasG[i] = biasv[i];
  }

  gsl_fit_wlinear (radG, 1, weightG, 1, biasG, 1, rad.size(), &c0, &c1, &cov00, &cov01, &cov11, &chisq);
  
  /*
  cout <<"best fit: bias(r) = ("<<c0<<") + ("<<c1<<") * r"<<endl;
  cout <<"covariance matrix:"<<endl;
  cout <<"[ "<<cov00<<", "<<cov01<<endl<<"#   "<<cov01<<", "<<cov11<<"]"<<endl; 
  cout <<"chisq = "<<chisq<<endl;
  */

  // estimate the mean bias from the best-fit model
 
  int step = 100;
  double r_min = radG[0];
  double r_max = radG[rad.size()-1];
  if (r_min>r_max*0.999) {
    string Err = "Error in cosmobl::TwoPointCorrelation::fit_bias of Bias.cpp: r_mim = " + conv(r_min,par::fDP3) + ", r_max = " + conv(r_max,par::fDP3) + "!";
    ErrorMsg(Err);
  }

  vector<double> radF = linear_bin_vector(step, r_min, r_max);


  double yf, yf_err;
  vector<double> bm, bm1, bm2;

  for (unsigned int i=0; i<radF.size(); i++) {
    gsl_fit_linear_est (radF[i], c0, c1, cov00, cov01, cov11, &yf, &yf_err);
    bm.push_back(yf);
    bm1.push_back(yf-yf_err);
    bm2.push_back(yf+yf_err);
  }

  m_bias = Average(bm);
  m_bias_min = Average(bm1);
  m_bias_max = Average(bm2);
}

