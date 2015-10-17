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
 *  @file CatalogueAnalysis/ModelTwoPointCorrelation/DispersionModelXiMeasured.cpp
 *
 *  @brief Methods of the class ModelTwoPointCorrelation used to model
 *  the two-point correlation function with the dispersion model
 *
 *  This file contains the implementation of the methods of the class
 *  ModelTwoPointCorrelation used to model the two-point correlation
 *  function with the dispersion model. The dark matter correlation
 *  function is estimated from the measured two-point correlation
 *  function.
 *
 *  @author Federico Marulli
 *
 *  @author federico.marulli3@unbo.it
 */

#include "ModelTwoPointCorrelation.h"
using namespace cosmobl;


// ============================================================================


void cosmobl::ModelTwoPointCorrelation::compute_xi2D_modelXiMeasured (vector< vector<double> > &Xi, double &beta, double &sigma12, Cosmology &cosm, double &Redshift, bool &NL, int &XiReal, double &pimax, string &dir_realxi, double &rApp, int &FV, bool bias_nl, double bA, double v_min, double v_max, int step_v)
{
  if (Xi.size()==0) 
    ErrorMsg("Error in cosmobl::ModelTwoPointCorrelation::compute_xi2D_modelXiMeasured of DispersionModel_XiMeasured.cpp: Xi.size()=0!");
  

  // ----- set default limits ----- 

  if (m_lim_index_fit.size()==0) {
    m_lim_index_fit.push_back(0);
    m_lim_index_fit.push_back(m_TwoP->sizeof_xi_lin()-1);
    m_lim_index_fit.push_back(0);
    m_lim_index_fit.push_back(m_TwoP->sizeof_xi_lin()-1);
  }
 

  // ----- compute the real space xi(r) ----- 

  m_TwoP->derive_real_xi(XiReal, pimax, dir_realxi, rApp);

  vector<double> rad_real_lin, xi_real_lin_extr;
  for (int i=0; i<m_TwoP->sizeof_xi_real_lin_extr(); i++) 
    if (m_TwoP->rad_real_lin(i)>0) {
      rad_real_lin.push_back(m_TwoP->rad_real_lin(i));
      xi_real_lin_extr.push_back(m_TwoP->xi_real_lin_extr(i));
    }


  // ----- compute the barred functions: xi_(r) and xi__(r) -----
   
  double r0 = m_TwoP->r0_linextr();
  double gamma = m_TwoP->gamma_linextr();

  vector<double> xi_, xi__;
  for (unsigned int i=0; i<rad_real_lin.size(); i++) {
    xi_.push_back(barred_xi_direct(rad_real_lin[i], rad_real_lin, xi_real_lin_extr, rApp, r0, gamma));
    xi__.push_back(barred_xi__direct(rad_real_lin[i], rad_real_lin, xi_real_lin_extr, rApp, r0, gamma));
  }
 

  // ----- compute xi(rp,pi) ----- 

  double var = (1.+Redshift)/cosm.HH(Redshift);
  double bias = 1.;

  for (int i=m_lim_index_fit[0]; i<m_lim_index_fit[1]; i++) {
    double rp = m_TwoP->rad_lin(i); 
    for (int j=m_lim_index_fit[2]; j<m_lim_index_fit[3]; j++) {
      double pi = m_TwoP->rad_lin(j);
      if (!NL) 
	Xi[i][j] = xi2D_lin_model(rp, pi, beta, bias, rad_real_lin, xi_real_lin_extr, xi_, xi__, -1, bias_nl, bA);
      else 
	Xi[i][j] = xi2D_model(rp, pi, beta, bias, sigma12, rad_real_lin, xi_real_lin_extr, xi_, xi__, var, FV, -1, bias_nl, bA, v_min, v_max, step_v);
    }
  }

}


// ============================================================================


void cosmobl::ModelTwoPointCorrelation::compute_multipoles_modelXiMeasured (vector<double> &xi0_lin, vector<double> &xi2_lin, vector<double> &xi4_lin, vector<double> &xi0_log, vector<double> &xi2_log, vector<double> &xi4_log, double &beta, double &sigma12, Cosmology &cosm, double &Redshift, bool &NL, int &XiReal, double &pimax, string &dir_realxi, double &rApp, int &FV, bool bias_nl, double bA, double v_min, double v_max, int step_v, int step_cos)
{
  xi0_lin.erase(xi0_lin.begin(), xi0_lin.end());
  xi2_lin.erase(xi2_lin.begin(), xi2_lin.end());
  xi4_lin.erase(xi4_lin.begin(), xi4_lin.end());
  xi0_log.erase(xi0_log.begin(), xi0_log.end());
  xi2_log.erase(xi2_log.begin(), xi2_log.end());
  xi4_log.erase(xi4_log.begin(), xi4_log.end());
 
  double val = -1.;


  // -------------------------------------------------------------------------------
  // ---- compute xi(rp,pi) with the dispersion model, using the measured xi(r) ----
  // -------------------------------------------------------------------------------

  int Nl = m_TwoP->sizeof_xi_lin()-1;
  vector< vector<double> > Xi (Nl, vector<double> (Nl,-1.e30));

  compute_xi2D_modelXiMeasured(Xi, beta, sigma12, cosm, Redshift, NL, XiReal, pimax, dir_realxi, rApp, FV, bias_nl, bA, v_min, v_max, step_v);

       
  // -----------------------------
  // ---- eliminate null bins ----
  // -----------------------------
  
  vector< vector<double> > XiR = Xi;

  vector<double> rad_linR1, rad_linR2;
  for (int i=0; i<Nl; i++) {
    rad_linR1.push_back(m_TwoP->rad_lin(i)); 
    rad_linR2.push_back(m_TwoP->rad_lin(i)); 
  }

  SubMatrix (rad_linR1, rad_linR2, XiR, -1.); 


  // -------------------------------------------------
  // ---- interpolate xi(rp,pi) in the r,mu plane ----
  // -------------------------------------------------
  
  double cos_min = 0.;
  double cos_max = 1.;
  vector<double> cos_lin = linear_bin_vector(step_cos, cos_min, cos_max);
  

  // ---- logarithimic binning ----

  vector<double> rad_log;
  for (int i=0; i<m_TwoP->sizeof_xi_log()-1; i++)
    rad_log.push_back(m_TwoP->rad_log(i));
 
  vector< vector<double> > xi_coslog_interp(rad_log.size(), vector<double>(cos_lin.size(),0.));
  
  for (unsigned int i=0; i<rad_log.size()-1; i++) {
    for (unsigned int j=0; j<cos_lin.size(); j++) {
      double rp = rad_log[i]*sqrt(1.-cos_lin[j]*cos_lin[j]);
      double pi = rad_log[i]*cos_lin[j];
      
      interpolation_extrapolation_2D (rp, pi, rad_linR1, rad_linR2, XiR, "Linear", 400, &val);
      xi_coslog_interp[i][j] = val;
    }
  }


  // ---- linear binning ----

  vector< vector<double> > xi_coslin_interp(rad_linR1.size(), vector<double>(cos_lin.size(),0.));

  for (unsigned int i=0; i<rad_linR1.size(); i++) 
    for (unsigned int j=0; j<cos_lin.size(); j++) {
      double rp = rad_linR1[i]*sqrt(1.-cos_lin[j]*cos_lin[j]);
      double pi = rad_linR1[i]*cos_lin[j];
      
      interpolation_extrapolation_2D (rp, pi, rad_linR1, rad_linR2, XiR, "Linear", -1, &val);
      xi_coslin_interp[i][j] = val;
    }
  
  
  // ---------------------------------------------------
  // ---- measure the multipoles of the correlation ----
  // ---------------------------------------------------  

  for (unsigned int i=0; i<xi_coslog_interp.size(); i++) {
    xi0_log.push_back(multipole_xi0(i,cos_lin,xi_coslog_interp));
    xi2_log.push_back(multipole_xi2(i,cos_lin,xi_coslog_interp));
    xi4_log.push_back(multipole_xi4(i,cos_lin,xi_coslog_interp));
  }
  
  for (unsigned int i=0; i<xi_coslin_interp.size(); i++) {
    xi0_lin.push_back(multipole_xi0(i,cos_lin,xi_coslin_interp));
    xi2_lin.push_back(multipole_xi2(i,cos_lin,xi_coslin_interp));
    xi4_lin.push_back(multipole_xi4(i,cos_lin,xi_coslin_interp));
  }
}


// ============================================================================


void cosmobl::ModelTwoPointCorrelation::compute_effective_multipoles_modelXiMeasured (vector<double> &xi0_lin, vector<double> &xi2_lin, vector<double> &xi4_lin, vector<double> &xi0_log, vector<double> &xi2_log, vector<double> &xi4_log, double &beta, double &sigma12, Cosmology &cosm, double &Redshift, bool &NL, int &XiReal, double &pimax, string &dir_realxi, double &rApp, int &FV, bool bias_nl, double bA, double v_min, double v_max, int step_v)
{
  xi0_lin.erase(xi0_lin.begin(), xi0_lin.end());
  xi2_lin.erase(xi2_lin.begin(), xi2_lin.end());
  xi4_lin.erase(xi4_lin.begin(), xi4_lin.end());
  xi0_log.erase(xi0_log.begin(), xi0_log.end());
  xi2_log.erase(xi2_log.begin(), xi2_log.end());
  xi4_log.erase(xi4_log.begin(), xi4_log.end());


  // -------------------------------------------------------------------------------
  // ---- compute xi(rp,pi) with the dispersion model, using the measured xi(r) ----
  // -------------------------------------------------------------------------------

  int Nl = m_TwoP->sizeof_xi_lin()-1;
  vector< vector<double> > Xi (Nl, vector<double> (Nl,-1.e30));
  
  compute_xi2D_modelXiMeasured(Xi, beta, sigma12, cosm, Redshift, NL, XiReal, pimax, dir_realxi, rApp, FV, bias_nl, bA, v_min, v_max, step_v);


  // ---------------------------------------------------
  // ---- measure the multipoles of the correlation ----
  // ---------------------------------------------------
  
  vector<double> rad_lin;
  for (int i=0; i<Nl; i++) rad_lin.push_back(m_TwoP->rad_lin(i)); 

  double delta_s = m_TwoP->linbinsz();

  for (int i=0; i<m_TwoP->nlinbins(); i++) {
    double ss = m_TwoP->rad_lin(i);
    xi0_lin.push_back(multipole_xi0(ss,rad_lin,rad_lin,Xi,delta_s));
    xi2_lin.push_back(multipole_xi2(ss,rad_lin,rad_lin,Xi,delta_s));
    xi4_lin.push_back(multipole_xi4(ss,rad_lin,rad_lin,Xi,delta_s));
  }

  for (int i=0; i<m_TwoP->nlogbins(); i++) {
    double ss = m_TwoP->rad_log(i);
    xi0_log.push_back(multipole_xi0(ss,rad_lin,rad_lin,Xi,delta_s));
    xi2_log.push_back(multipole_xi2(ss,rad_lin,rad_lin,Xi,delta_s));
    xi4_log.push_back(multipole_xi4(ss,rad_lin,rad_lin,Xi,delta_s));
  }

}


// ============================================================================


void cosmobl::ModelTwoPointCorrelation::measure_beta_KaiserLimit_XiMeasured (int &XiReal, double &pimax, string &dir_realxi, double &rApp, double beta_guess)
{

  // ----- compute the real space xi(r) ----- 

  m_TwoP->derive_real_xi(XiReal, pimax, dir_realxi, rApp);

  vector<double> rad_real_lin, xi_real_lin, error_xi_real_lin;
  for (int i=0; i<m_TwoP->sizeof_xi_real_lin(); i++) 
    if (m_TwoP->rad_real_lin(i)>0) {
      rad_real_lin.push_back(m_TwoP->rad_real_lin(i));
      xi_real_lin.push_back(m_TwoP->xi_real_lin(i));
      error_xi_real_lin.push_back(m_TwoP->error_xi_real_lin(i));
    }
 

  // ----- get the measured values of xi(s) ----- 

  vector<double> rad_lin, xi_lin, error_xi_lin;

  for (int i=0; i<m_TwoP->sizeof_xi_lin(); i++) 
    if (m_TwoP->rad_lin(i)>0 && m_TwoP->xi_lin(i)>0) {
      rad_lin.push_back(m_TwoP->rad_lin(i));
      xi_lin.push_back(m_TwoP->xi_lin(i));
      error_xi_lin.push_back(m_TwoP->error_xi_lin(i));
    }


  // ----- estimate xi(s)/xi(r) ----- 

  vector<double> ratio, error_ratio;
  double XiR, ErrXiR, Err; 
  
  for (unsigned int i=0; i<rad_lin.size(); i++) {
    interpolation_extrapolation(rad_lin[i],rad_real_lin,xi_real_lin,"Poly",4,&XiR,&Err);
    interpolation_extrapolation(rad_lin[i],rad_real_lin,error_xi_real_lin,"Poly",4,&ErrXiR,&Err);
   
    double RR = xi_lin[i]/XiR;
    double ER = sqrt(pow(error_xi_lin[i]/xi_lin[i],2)+pow(ErrXiR/XiR,2))*RR;

    ratio.push_back(RR);
    error_ratio.push_back(ER);
  }


  // ----- linear Kaiser model for xi(s)/xi(r) ----- 

  double (*p_ratio) (double, void *, vector<double>);
  p_ratio = xi_ratio;


  // ----- diagonal chi^2 -----

  vector<double> pp;

  Chi2 chi2 (rad_lin, ratio, error_ratio, p_ratio, &pp); 


  // ----- set the limits for the fit -----

  if (m_lim_fit.size()==2)
    chi2.set_limits (m_lim_fit[0], m_lim_fit[1]);
  

  // ----- set the priors -----

  vector< vector<double> > priors(1,vector<double>(2)); priors[0][0] = 0.; priors[0][1] = 2.;

  chi2.set_par_limits(priors); 
  

  // ----- estimate the best-fit value of beta -----

  vector<double> pstart; pstart.push_back(beta_guess);
  
  chi2.get_bestfit(pstart); 
  m_beta_best = chi2.bestfit(0); 

  vector<double> bb(1); bb[0] = m_beta_best; 
  m_E_min = chi2.get_chi2(bb);
}


// ============================================================================


void cosmobl::ModelTwoPointCorrelation::measure_beta_XiMeasured (int &CHI, int &XiReal, double &pimax, string &dir_realxi, double &rApp, double beta_guess)
{ 
  cout <<endl<<"I start to model the RSD..."<<endl;


  // ----- set default limits ----- 

  if (m_lim_index_fit.size()==0) {
    m_lim_index_fit.push_back(0);
    m_lim_index_fit.push_back(m_TwoP->sizeof_xi_lin()-1);
    m_lim_index_fit.push_back(0);
    m_lim_index_fit.push_back(m_TwoP->sizeof_xi_lin()-1);
  }


  // ----- compute the real space xi(r) ----- 

  m_TwoP->derive_real_xi(XiReal, pimax, dir_realxi, rApp);
  

  // ----- get the values of real space xi(r) ----- 

  vector<double> rad_real_lin, xi_real_lin_extr, xi_real_lin;
  for (int i=0; i<m_TwoP->sizeof_xi_real_lin_extr(); i++) 
    if (m_TwoP->rad_real_lin(i)>0) {
      rad_real_lin.push_back(m_TwoP->rad_real_lin(i));
      xi_real_lin_extr.push_back(m_TwoP->xi_real_lin_extr(i));
    }
  

  // ----- interpolate the correlation functions at the scales used for the fit ----- 

  double r0 = m_TwoP->r0_linextr();
  double gamma = m_TwoP->gamma_linextr();

  vector<double> Rp, Pi, xi_real, xi_, xi__, P2, P4;
  vector< vector<double> > xi(m_lim_index_fit[1]-m_lim_index_fit[0]), error_xi(m_lim_index_fit[1]-m_lim_index_fit[0]);
  
  for (int i=m_lim_index_fit[0]; i<m_lim_index_fit[1]; i++) {
    double rp = m_TwoP->rad_lin(i);
    Rp.push_back(rp);

    for (int j=m_lim_index_fit[2]; j<m_lim_index_fit[3]; j++) {
      double pi = m_TwoP->rad_lin(j);      
      if (i==m_lim_index_fit[0]) Pi.push_back(pi); 

      double rr = sqrt(rp*rp+pi*pi);
      double cos_i = pi/rr;
      
      xi[i-m_lim_index_fit[0]].push_back(m_TwoP->xi_2d_lin(i,j));
      error_xi[i-m_lim_index_fit[0]].push_back(m_TwoP->error_xi_2d_lin(i,j));

      double Xi, Err; interpolation_extrapolation(rr,rad_real_lin,xi_real_lin_extr,"Poly",4,&Xi,&Err);
      xi_real.push_back(Xi);

      xi_.push_back(barred_xi_direct(rr,rad_real_lin,xi_real_lin_extr,rApp,r0,gamma));
      xi__.push_back(barred_xi__direct(rr,rad_real_lin,xi_real_lin_extr,rApp,r0,gamma));
      
      //xi_.push_back(barred_xi_(rr,rad_real_lin,xi_real_lin_extr,rApp,r0,gamma));
      //xi__.push_back(barred_xi__(rr,rad_real_lin,xi_real_lin_extr,rApp,r0,gamma));

      P2.push_back(P_2(cos_i));
      P4.push_back(P_4(cos_i));
    }
  }
 

  
  // ----- linear Kaiser model for xi(rp,pi) -----  

  cosmobl::glob::STR_xi2D_model pp;
  pp.xi_real = xi_real;
  pp.xi_ = xi_;
  pp.xi__ = xi__;
  pp.P2 = P2;
  pp.P4 = P4;
  pp.bias_nl = 0;

  double (*p_xi2D) (double, double, void *, vector<double>);
  p_xi2D = xi2D_lin_model;  
  
  
  // ----- diagonal chi^2: 0 -> chi^2, 1 -> log(chi^2) ----- 
  
  Chi2 chi2 (Rp, Pi, xi, error_xi, p_xi2D, &pp, 1, CHI); // the first two inputs are useless here


  // ----- set the priors -----
  
  vector< vector<double> > priors(1,vector<double>(2)); priors[0][0] = 0.; priors[0][1] = 2.;
    
  chi2.set_par_limits(priors);
    

  // ----- estimate the best-fit value of beta -----
  
  vector<double> pstart(1,beta_guess);
  chi2.get_bestfit(pstart); 
  m_beta_best = chi2.bestfit(0); 
  
  vector<double> bb(1); bb[0] = m_beta_best; 
  m_E_min = chi2.get_chi2(bb);
  
}


// ============================================================================


void cosmobl::ModelTwoPointCorrelation::measure_beta_sigma12_DispersionModel_XiMeasured (int &CHI, int &FV, Cosmology &cosm, double &Redshift, int &XiReal, double &pimax, string &dir_realxi, double &rApp, string fileMCMC, double v_min, double v_max, int step_v)
{ 
  cout <<endl<<"I start modelling the RSD..."<<endl;
 

  // ----- set default limits ----- 

  if (m_lim_index_fit.size()==0) {
    m_lim_index_fit.push_back(0);
    m_lim_index_fit.push_back(m_TwoP->sizeof_xi_lin()-1);
    m_lim_index_fit.push_back(0);
    m_lim_index_fit.push_back(m_TwoP->sizeof_xi_lin()-1);
  }


  // ----- compute the real space xi(r) ----- 

  m_TwoP->derive_real_xi(XiReal, pimax, dir_realxi, rApp);


  // ----- get the values of real space xi(r) ----- 

  vector<double> rad_real_lin, xi_real_lin_extr, xi_real_lin;
  for (int i=0; i<m_TwoP->sizeof_xi_real_lin_extr(); i++) 
    if (m_TwoP->rad_real_lin(i)>0) {
      rad_real_lin.push_back(m_TwoP->rad_real_lin(i));
      xi_real_lin_extr.push_back(m_TwoP->xi_real_lin_extr(i));
    }
 

  // ----- interpolate the correlation functions at the scales used in the fit ----- 

  double r0 = m_TwoP->r0_linextr();
  double gamma = m_TwoP->gamma_linextr();

  vector<double> Rp, Pi, Rp_model, Pi_model, xi_real, xi_, xi__, P2, P4, vel;
  vector< vector<double> > xi(m_lim_index_fit[1]-m_lim_index_fit[0]), error_xi(m_lim_index_fit[1]-m_lim_index_fit[0]);

  double delta_v = (v_max-v_min)/step_v;
  double var = (1.+Redshift)/cosm.HH(Redshift);
  
  for (int i=m_lim_index_fit[0]; i<m_lim_index_fit[1]; i++) {
    double rp = m_TwoP->rad_lin(i);
    Rp.push_back(rp);

    for (int j=m_lim_index_fit[2]; j<m_lim_index_fit[3]; j++) {
      double pi = m_TwoP->rad_lin(j);
      if (i==m_lim_index_fit[0]) Pi.push_back(pi); 

      xi[i-m_lim_index_fit[0]].push_back(m_TwoP->xi_2d_lin(i,j));
      error_xi[i-m_lim_index_fit[0]].push_back(m_TwoP->error_xi_2d_lin(i,j));


      // xi(rp,pi) convolution with the velocity distribution function f(v)		
      
      double Vel = v_min;

      for (int k=0; k<step_v; k++) {
	vel.push_back(Vel);

	double pi_new = pi-Vel*var;
	
	Rp_model.push_back(rp);
	Pi_model.push_back(pi_new); 

	double rr = sqrt(rp*rp+pi_new*pi_new);
	double cos_i = pi_new/rr;
	
	double Xi, Err; interpolation_extrapolation(rr,rad_real_lin,xi_real_lin_extr,"Poly",4,&Xi,&Err);
	xi_real.push_back(Xi);
	
	xi_.push_back(barred_xi_direct(rr,rad_real_lin,xi_real_lin_extr,rApp,r0,gamma));
	xi__.push_back(barred_xi__direct(rr,rad_real_lin,xi_real_lin_extr,rApp,r0,gamma));
	
	//xi_.push_back(barred_xi_(rr,rad_real_lin,xi_real_lin_extr,rApp,r0,gamma));
	//xi__.push_back(barred_xi__(rr,rad_real_lin,xi_real_lin_extr,rApp,r0,gamma));
	
	P2.push_back(P_2(cos_i));
	P4.push_back(P_4(cos_i));

	Vel += delta_v;
      }
    }
  }


  // ----- dispersion model for xi(rp,pi) ----- 

  cosmobl::glob::STR_xi2D_model pp;
  pp.rp = Rp_model;
  pp.pi = Pi_model;
  pp.xi_real = xi_real;
  pp.xi_ = xi_;
  pp.xi__ = xi__;
  pp.P2 = P2;
  pp.P4 = P4;
  pp.vel = vel;
  pp.step_v = step_v;
  pp.delta_v = delta_v;
  pp.FV = FV;
  pp.bias_nl = 0;

  double (*p_xi2D) (double, double, void *, vector<double>);
  p_xi2D = xi2D_model;  
  

  // ----- diagonal chi^2: 0 -> chi^2, 1 -> log(chi^2) ----- 
  
  Chi2 chi2 (Rp, Pi, xi, error_xi, p_xi2D, &pp, 1, CHI); // the first two inputs are useless here
  

  // ----- set the priors -----
  
  vector< vector<double> > priors(2,vector<double>(2)); 
  priors[0][0] = 0.; priors[0][1] = 2.;
  priors[1][0] = 50.; priors[1][1] = 5000.; 
  // the minimum prior, priors[1][0], is used to minimize numerical problems in the sigma_v convolution (check!!!)

  chi2.set_par_limits(priors);


  // ----- estimate the best-fit value of beta -----
  
  vector<double> pstart;
  pstart.push_back(0.5);
  pstart.push_back(300.);

  chi2.get_bestfit(pstart); 
  
  m_beta_best = chi2.bestfit(0); 
  m_sigma12_best = chi2.bestfit(1); 

  vector<double> bb(2); bb[0] = m_beta_best, bb[1] = m_sigma12_best; 
  m_E_min = chi2.get_chi2(bb);
 
  
  // ----- performe an MCMC analysis to improve the best-fit values and estimate marginalized errors ----- 

  if (fileMCMC!="NULL") { 

    ErrorMsg("WORK IN PROGRESS... (see DispersionModel_XiMeasured.cpp");

    /*
    cout <<"I'm doing the MCMC analysis..."<<endl;

    MyState ss (XX);
 
    Plog_chi2_xirppi plog (xi, error_xi, xi_real, xi_, xi__, P2, P4, vel, delta_v, step_v, FV);
  
    Proposal_chi2_xirppi propose (10102,0.01);
  
    Doub accept;
  
    ofstream fout (fileMCMC.c_str()); checkIO (fileMCMC,0);

    vector<double> beta_MCMC, sigma12_MCMC;

    // without this "burn-in", the MCMC parameter distribution is not gaussian: it's important to check the convergence!!!
    //for (int i=0; i<1000; i++) accept = mcmcstep_chi2_xirppi(1,ss,plog,propose); 
  
    int intS = 0;

    for (int i=0; i<300; i++) {

      accept = mcmcstep_chi2_xirppi(10,ss,plog,propose);

      beta_MCMC.push_back(ss.Par[0]);
      sigma12_MCMC.push_back(ss.Par[1]);

      intS ++;
      if (intS==10) {
	intS = 0;
	cout <<"MCMC step: "<<i<<" --> beta = "<<ss.Par[0]<<", sigma12 = "<<ss.Par[1]<<", plog = "<<ss.plog<<endl;
      }

      fout <<i<<"   "<<ss.Par[0]<<"   "<<ss.Par[1]<<"   "<<ss.plog<<"   "<<accept<<endl;
    }
    fout.clear(); fout.close(); cout <<"I wrote the file: "<<fileMCMC<<endl<<endl;
   

    cout <<"MCMC ----> beta = "<<ss.Par[0]<<", sigma12 = "<<ss.Par[1]<<endl<<endl;

    cout <<"mean beta from the MCMC points: "<<Average(beta_MCMC)<<" +o- "<<Sigma(beta_MCMC)<<endl;
    cout <<"mean sigma12 from the MCMC points: "<<Average(sigma12_MCMC)<<" +o- "<<Sigma(sigma12_MCMC)<<endl;

    //m_beta_best = XX[0];
    //m_sigma12_best = XX[1];
   
    m_beta_best = Average(beta_MCMC);
    error_beta = Sigma(beta_MCMC);
    m_sigma12_best = Average(sigma12_MCMC);
    error_sigma12 = Sigma(sigma12_MCMC);
    */
 
  }

}


// ============================================================================


void cosmobl::ModelTwoPointCorrelation::measure_beta_sigma12_DispersionModel_XiMeasured_multipoles (int &CHI, int &FV, Cosmology &cosm, double &Redshift, int &XiReal, double &pimax, string &dir_realxi, double &rApp, double v_min, double v_max, int step_v)
{   
  cout <<endl<<"I start modelling the RSD..."<<endl;

 
  // ----- set default limits ----- 

  if (m_lim_index_fit.size()==0) {
    m_lim_index_fit.push_back(0);
    m_lim_index_fit.push_back(m_TwoP->sizeof_xi_lin()-1);
  }


  // ----- compute the real space xi(r) ----- 

  m_TwoP->derive_real_xi(XiReal, pimax, dir_realxi, rApp);
  

  // ----- get the values of real space xi(r) ----- 

  vector<double> rp, pi, rad_real_lin, xi_real_lin_extr, xi_real_lin;
  for (int i=0; i<m_TwoP->sizeof_xi_real_lin_extr(); i++) 
    if (m_TwoP->rad_real_lin(i)>0) {
      rad_real_lin.push_back(m_TwoP->rad_real_lin(i));
      xi_real_lin_extr.push_back(m_TwoP->xi_real_lin_extr(i));
    }
  

  // ----- store the (first) measured multipoles in a single vector

  vector<double> rad, xi_mult, error_xi_mult;
  vector<int> type;

  for (int i=m_lim_index_fit[0]; i<m_lim_index_fit[1]; i++) {
    rad.push_back(m_TwoP->rad_log(i));
    xi_mult.push_back(m_TwoP->xi0_log(i));
    error_xi_mult.push_back(m_TwoP->error_xi0_log(i));
    type.push_back(1);
  }
  
  for (int i=m_lim_index_fit[0]; i<m_lim_index_fit[1]; i++) {
    rad.push_back(m_TwoP->rad_log(i));
    xi_mult.push_back(m_TwoP->xi2_log(i));
    error_xi_mult.push_back(m_TwoP->error_xi2_log(i));
    type.push_back(2);
  }


  // ----- interpolate the correlation functions at the scales used in the fit ----- 

  double r0 = m_TwoP->r0_linextr();
  double gamma = m_TwoP->gamma_linextr();

  vector<double> xi_real, xi_, xi__, P2, P4, vel;

  double delta_v = (v_max-v_min)/step_v;
  double var = (1.+Redshift)/cosm.HH(Redshift);

  for (int i=0; i<m_TwoP->sizeof_xi_lin()-1; i++) {
    double Rp = m_TwoP->rad_lin(i);

    rp.push_back(m_TwoP->rad_lin(i));
    pi.push_back(m_TwoP->rad_lin(i));

    for (int j=0; j<m_TwoP->sizeof_xi_lin()-1; j++) {
      double Pi = m_TwoP->rad_lin(j);


      // xi(rp,pi) convolution with the velocity distribution function f(v)		
      
      double Vel = v_min;

      for (int k=0; k<step_v; k++) {
	vel.push_back(Vel);

	double pi_new = Pi-Vel*var;
	
	double rr = sqrt(Rp*Rp+pi_new*pi_new);
	double cos_i = pi_new/rr;

	double Xi, Err; interpolation_extrapolation(rr,rad_real_lin,xi_real_lin_extr,"Poly",4,&Xi,&Err);
	xi_real.push_back(Xi);
	
	xi_.push_back(barred_xi_direct(rr,rad_real_lin,xi_real_lin_extr,rApp,r0,gamma));
	xi__.push_back(barred_xi__direct(rr,rad_real_lin,xi_real_lin_extr,rApp,r0,gamma));
	
	//xi_.push_back(barred_xi_(rr,rad_real_lin,xi_real_lin_extr,rApp,r0,gamma));
	//xi__.push_back(barred_xi__(rr,rad_real_lin,xi_real_lin_extr,rApp,r0,gamma));
	
	P2.push_back(P_2(cos_i));
	P4.push_back(P_4(cos_i));

	Vel += delta_v;
      }
    }
  }

  cosmobl::glob::STR_xi2D_model pp;
  pp.rp = rp;
  pp.pi = pi;
  pp.xi_real = xi_real;
  pp.xi_ = xi_;
  pp.xi__ = xi__;
  pp.P2 = P2;
  pp.P4 = P4;
  pp.vel = vel;
  pp.step_v = step_v;
  pp.delta_v = delta_v;
  pp.FV = FV;
  pp.lim_index_fit = m_lim_index_fit;
  pp.dim = m_TwoP->sizeof_xi_lin()-1;
  pp.type = type;
  
 

  // ----- dispersion model for the multipoles ----- 
 
  double (*p_mult) (double, void *, vector<double>);
  p_mult = multipoles;  
 

  // diagonal chi^2: CHI=0 -> chi^2, CHI=1 -> log(chi^2)

  Chi2 chi2 (rad, xi_mult, error_xi_mult, p_mult, &pp, 1, CHI);


  // ----- set the priors -----
  
  vector< vector<double> > priors(2,vector<double>(2)); 
  priors[0][0] = 0.; priors[0][1] = 2.;
  priors[1][0] = 0.; priors[1][1] = 500.;

  chi2.set_par_limits(priors);


  // ----- estimate the best-fit value of beta -----
  
  vector<double> pstart;
  pstart.push_back(0.3);
  pstart.push_back(150.);
  
  chi2.get_bestfit(pstart); 
 
  m_beta_best = chi2.bestfit(0); 
  m_sigma12_best = chi2.bestfit(1); 

  vector<double> bb(2); bb[0] = m_beta_best, bb[1] = m_sigma12_best; 
  m_E_min = chi2.get_chi2(bb);
 
}


// ============================================================================


double cosmobl::ModelTwoPointCorrelation::relative_error_beta_catalogue (double &Volume, double &rMin, double &rMax, bool proj, bool NL) 
{ // from Eq. 20 of Bianchi et al. 2012

  double density = m_TwoP->nObjects()/Volume; 

  double Bias = m_TwoP->mean_bias(rMin,rMax,proj,NL);

  cout <<"bias = "<<Bias<<", density = "<<density<<", Volume = "<<Volume<<endl;
  
  return relative_error_beta (Bias, Volume, density);
}


// ============================================================================


double cosmobl::ModelTwoPointCorrelation::relative_error_beta_catalogue_box (double &rMin, double &rMax, bool proj, bool NL) 
{ // from Eq. 22 of Bianchi et al. 2012
  
  double Volume = (m_TwoP->Xmax()-m_TwoP->Xmin())*(m_TwoP->Ymax()-m_TwoP->Ymin())*(m_TwoP->Zmax()-m_TwoP->Zmin());

  double density = m_TwoP->nObjects()/Volume; 

  double Bias = m_TwoP->mean_bias(rMin, rMax, proj, NL);

  cout <<"bias = "<<Bias<<", density = "<<density<<", Volume = "<<Volume<<endl;
  
  return relative_error_beta (Bias, Volume, density);
}


