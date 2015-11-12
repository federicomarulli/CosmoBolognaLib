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
 *  @file CatalogueAnalysis/ModelTwoPointCorrelation/DispersionModel.cpp
 *
 *  @brief Methods of the class ModelTwoPointCorrelation used to model
 *  the two-point correlation function with the dispersion model
 *
 *  This file contains the implementation of the methods of the class
 *  ModelTwoPointCorrelation used to model the two-point correlation
 *  function with the dispersion model. The dark matter correlation
 *  function is estimated from theory (e.g. with CAMB, CLASS or
 *  MPTbreeze)
 *
 *  @author Federico Marulli
 *
 *  @author federico.marulli3@unbo.it
 */

#include "ModelTwoPointCorrelation.h"
using namespace cosmobl;


// ============================================================================


void cosmobl::ModelTwoPointCorrelation::compute_xi2D_model (vector< vector<double> > &Xi, double &f_sigma8, double &bias_sigma8, double &sigma12, Cosmology &cosm, double &redshift, bool &NL, string &author, string &Model, int &FV, bool bias_nl, double bA, bool xiType, bool xiNL, double v_min, double v_max, int step_v, int norm, double r_min, double r_max, double k_min, double k_max, double aa, bool GSL, double prec, string file_par)
{
  // ----- set default limits ----- 

  if (m_lim_index_fit.size()==0) {
    m_lim_index_fit.push_back(0);
    m_lim_index_fit.push_back(m_TwoP->sizeof_xi_lin()-1);
    m_lim_index_fit.push_back(0);
    m_lim_index_fit.push_back(m_TwoP->sizeof_xi_lin()-1);
  }
 

  // ----- compute the real space xi(r) ----- 

  double k_star = cosm.k_star (author, redshift, Model, k_max, file_par); 

  vector<double> rad_real_lin, xi_real_lin, xi_, xi__;

  cosm.get_xi (rad_real_lin, xi_real_lin, author, redshift, Model, xiType, k_star, xiNL, norm, r_min, r_max, k_min, k_max, aa, GSL, prec, file_par);

  cosm.get_barred_xi (rad_real_lin, xi_real_lin, xi_, xi__, author, redshift, xiType, k_star, xiNL, norm, r_min, r_max, k_min, k_max, aa, GSL, prec, file_par);
 

  // ----- compute xi(rp,pi) ----- 

  double var = (1.+redshift)/cosm.HH(redshift);

  double beta = f_sigma8/bias_sigma8;
  double bias = bias_sigma8/cosm.sigma8_Pk(author,redshift,Model);

  for (int i=m_lim_index_fit[0]; i<m_lim_index_fit[1]; i++) {
    double rp = m_TwoP->rad_lin(i); 
    for (int j=m_lim_index_fit[2]; j<m_lim_index_fit[3]; j++) {
      double pi = m_TwoP->rad_lin(j);
      if (!NL) 
	Xi[i][j] = xi2D_lin_model (rp, pi, beta, bias, rad_real_lin, xi_real_lin, xi_, xi__, -1, bias_nl, bA);
      else 
	Xi[i][j] = xi2D_model (rp, pi, beta, bias, sigma12, rad_real_lin, xi_real_lin, xi_, xi__, var, FV, -1, bias_nl, bA, v_min, v_max, step_v);
    }
  }

}


// ============================================================================


void cosmobl::ModelTwoPointCorrelation::measure_fsigma8_KaiserLimit (Cosmology &cosm, double &bias_sigma8, string &author, double &redshift, string &Model, int &CHI, double fsigma8_guess, bool xiType, double k_star, bool xiNL, int norm, double r_min, double r_max, double k_min, double k_max, double aa, bool GSL, double prec, string file_par, __attribute__((unused)) int rank)
{
  // ----- get the measured values of xi(s) ----- 
 
  vector<double> rad_log, xi_log, error_xi_log;
  
  for (int i=0; i<m_TwoP->sizeof_xi_log(); i++) 
    if (m_TwoP->rad_log(i)>0 && m_TwoP->xi_log(i)>0) {
      rad_log.push_back(m_TwoP->rad_log(i));
      xi_log.push_back(m_TwoP->xi_log(i));
      error_xi_log.push_back(m_TwoP->error_xi_log(i));
    }

 
  // ----- compute the real space DM xi(r) ----- 

  vector<double> rr, XiDM;
  
  cosm.get_xi (rr, XiDM, author, redshift, Model, xiType, k_star, xiNL, norm, r_min, r_max, k_min, k_max, aa, GSL, prec, file_par);

  
  // ----- interpolate the real space DM xi(r) at the scales used for the fit ----- 

  vector<double> lgrr, lgXi;
  for (unsigned int i=0; i<rr.size(); i++) 
    if (rr[i]>0 && XiDM[i]>0) {
      lgrr.push_back(log10(rr[i]));
      lgXi.push_back(log10(XiDM[i]));
    }

  vector<double> xi_DM;

  for (unsigned int i=0; i<rad_log.size(); i++) 
    xi_DM.push_back(pow(10.,interpolated(log10(rad_log[i]), lgrr, lgXi, "Linear", -1)));
 

  // ----- linear Kaiser model for the monopole xi(s) ----- 

  cosmobl::glob::STR_xi0_model pp;
  pp.bias_sigma8 = bias_sigma8;
  pp.sigma8z = cosm.sigma8_Pk(author,redshift,Model);
  pp.xi_DM = xi_DM;

  double (*p_Kaiser) (double, void *, vector<double>);
  p_Kaiser = multipole_xi0_model;
 

  // ----- diagonal chi^2: CHI=0 -> chi^2, CHI=1 -> log(chi^2) -----

  Chi2 chi2 (rad_log, xi_log, error_xi_log, p_Kaiser, &pp, 0, CHI); 

 
  // ----- set the limits for the fit -----
  
  if (m_lim_fit.size()==2) chi2.set_limits (m_lim_fit[0], m_lim_fit[1]);
  else chi2.set_limits (0, m_TwoP->sizeof_xi_log()-1);
  

  // ----- set the priors -----

  vector< vector<double> > priors(1,vector<double>(2)); 
  priors[0][0] = 0.; priors[0][1] = 2.;
  chi2.set_par_limits(priors); 


  // ----- estimate the best-fit value of beta -----

  vector<double> ppp; ppp.push_back(fsigma8_guess); 

  chi2.get_bestfit(ppp); 
 
  m_f_sigma8_best = chi2.bestfit(0);
  m_beta_best = chi2.bestfit(0)/bias_sigma8; 

  vector<double> bb(1); bb[0] = m_f_sigma8_best;
  m_E_min = chi2.get_chi2(bb);
  
}


// =====================================================================================


void cosmobl::ModelTwoPointCorrelation::measure_fsigma8_bsigma8 (Cosmology &cosm, string &author, double &redshift, string &Model, int &CHI, int &FV, bool &NL, vector< vector<double> > Priors, vector<double> Pstart, bool bias_nl, bool xiType, double k_star, bool xiNL, double v_min, double v_max, int step_v, int norm, double r_min, double r_max, double k_min, double k_max, double aa, bool GSL, double prec, string file_par, __attribute__((unused)) int rank)
{
  // ----- get the measured values of xi(rp,pi) ----- 

  vector<double> rad_lin(m_TwoP->sizeof_xi_lin()-1);
  vector< vector<double> > xi_2d_lin(m_TwoP->sizeof_xi_lin()-1), error_xi_2d_lin(m_TwoP->sizeof_xi_lin()-1);
  
  for (unsigned int i=0; i<rad_lin.size(); i++) {
    rad_lin[i] = m_TwoP->rad_lin(i);
    for (unsigned int j=0; j<rad_lin.size(); j++) {
      xi_2d_lin[i].push_back(m_TwoP->xi_2d_lin(i,j));
      error_xi_2d_lin[i].push_back(m_TwoP->error_xi_2d_lin(i,j));
    }
  }

  
  // ----- compute the real space xi(r) ----- 

  vector<double> rr, Xi, Xi_, Xi__;
  
  cosm.get_xi(rr, Xi, author, redshift, Model, xiType, k_star, xiNL, norm, r_min, r_max, k_min, k_max, aa, GSL, prec, file_par);

  cosm.get_barred_xi(rr, Xi, Xi_, Xi__, author, redshift, xiType, k_star, xiNL, norm, r_min, r_max, k_min, k_max, aa, GSL, prec, file_par);
   
  vector<double> lgrr, lgXi, lgXi_, lgXi__;
  for (unsigned int i=0; i<rr.size(); i++) 
    if (rr[i]>0 && Xi[i]>0 && Xi_[i]>0 && Xi__[i]>0) {
      lgrr.push_back(log10(rr[i]));
      lgXi.push_back(log10(Xi[i]));
      lgXi_.push_back(log10(Xi_[i]));
      lgXi__.push_back(log10(Xi__[i]));
    }
  
  
  // ----- set default limits ----- 

  if (m_lim_index_fit.size()==0) {
    m_lim_index_fit.push_back(0);
    m_lim_index_fit.push_back(m_TwoP->sizeof_xi_lin()-1);
    m_lim_index_fit.push_back(0);
    m_lim_index_fit.push_back(m_TwoP->sizeof_xi_lin()-1);
  }


  // ----- interpolate the correlation functions at the scales used for the fit ----- 

  vector<double> Rp, Pi, Rp_model, Pi_model, rad_real, xi_real, xi_, xi__, vel, P2, P4;
  vector< vector<double> > xi(m_lim_index_fit[1]-m_lim_index_fit[0]), error_xi(m_lim_index_fit[1]-m_lim_index_fit[0]);

  double delta_v = (v_max-v_min)/step_v;
  double var = (1.+redshift)/cosm.HH(redshift);

  for (int i=m_lim_index_fit[0]; i<m_lim_index_fit[1]; i++) {
    double rp = m_TwoP->rad_lin(i);
    Rp.push_back(rp);

    for (int j=m_lim_index_fit[2]; j<m_lim_index_fit[3]; j++) {
      double pi = m_TwoP->rad_lin(j);
      if (i==m_lim_index_fit[0]) Pi.push_back(pi); 

      xi[i-m_lim_index_fit[0]].push_back(m_TwoP->xi_2d_lin(i,j));
      error_xi[i-m_lim_index_fit[0]].push_back(m_TwoP->error_xi_2d_lin(i,j));
 
      double lgXXi, lgXXi_, lgXXi__; 

      if (!NL) {

	double rr = sqrt(rp*rp+pi*pi);
	double cos_i = pi/rr;

	lgXXi = interpolated(log10(rr), lgrr, lgXi, "Linear", -1);
	lgXXi_ = interpolated(log10(rr), lgrr, lgXi_, "Linear", -1);
	lgXXi__ = interpolated(log10(rr), lgrr, lgXi__, "Linear", -1);

	xi_real.push_back(pow(10.,lgXXi));
	xi_.push_back(pow(10.,lgXXi_));
	xi__.push_back(pow(10.,lgXXi__));

	P2.push_back(P_2(cos_i));
	P4.push_back(P_4(cos_i));
      }
      
      else {

        // xi(rp,pi) convolution with the velocity distribution function f(v)		

	double Vel = v_min;
	
	for (int k=0; k<step_v; k++) {
	  vel.push_back(Vel);

	  double pi_new = pi-Vel*var;

	  Rp_model.push_back(rp);
	  Pi_model.push_back(pi_new); 

	  double rr = sqrt(rp*rp+pi_new*pi_new);
	  double cos_i = pi_new/rr;

	  lgXXi = interpolated(log10(rr), lgrr, lgXi, "Linear", -1);
	  lgXXi_ = interpolated(log10(rr), lgrr, lgXi_, "Linear", -1);
	  lgXXi__ = interpolated(log10(rr), lgrr, lgXi__, "Linear", -1);

	  xi_real.push_back(pow(10.,lgXXi));
	  xi_.push_back(pow(10.,lgXXi_));
	  xi__.push_back(pow(10.,lgXXi__));

	  P2.push_back(P_2(cos_i));
	  P4.push_back(P_4(cos_i));

	  Vel += delta_v;
	}
	
      }
    }
  }
 
  
  cosmobl::glob::STR_xi2D_model pp;
  pp.rp = Rp_model;
  pp.pi = Pi_model;
  pp.xi_real = xi_real;
  pp.xi_ = xi_;
  pp.xi__ = xi__;
  pp.P2 = P2;
  pp.P4 = P4;
  pp.bias_nl = bias_nl;
  if (NL) {
    pp.vel = vel;
    pp.step_v = step_v;
    pp.delta_v = delta_v;
    pp.FV = FV;
  }


  // ----- dispersion model for xi(rp,pi) ----- 
 
  double (*p_xi2D) (double, double, void *, vector<double>);
  if (NL) p_xi2D = xi2D_model; 
  else p_xi2D = xi2D_lin_model;  

  Chi2 chi2 (Rp, Pi, xi, error_xi, p_xi2D, &pp, 1, CHI);
  

  // ----- set the priors -----

  vector<double> pstart;

  if (!NL) {
    if (Priors.size()>0 && Priors.size()!=2) ErrorMsg("Error in ModelTwoPointCorrelation::measure_fsigma8_bsigma8 of DispersionModel.cpp!");
    vector< vector<double> > priors(2,vector<double>(2)); 
    priors[0][0] = (Priors.size()>0) ? Priors[0][0] : 0.; priors[0][1] = (Priors.size()>0) ? Priors[0][1] : 10.;
    priors[1][0] = (Priors.size()>0) ? Priors[1][0] : 0.; priors[1][1] = (Priors.size()>0) ? Priors[1][1] : 10.;
    chi2.set_par_limits(priors);
    if (Pstart.size()>0) {
      pstart.push_back(Pstart[0]);
      pstart.push_back(Pstart[1]);
    } else {
      pstart.push_back(1.);
      pstart.push_back(1.);
    }
  }
  else {
    if (!bias_nl) {
      if (Priors.size()>0 && Priors.size()!=3) ErrorMsg("Error in ModelTwoPointCorrelation::measure_fsigma8_bsigma8 of DispersionModel.cpp!");
      vector< vector<double> > priors(3,vector<double>(3)); 
      priors[0][0] = (Priors.size()>0) ? Priors[0][0] : 0.; priors[0][1] = (Priors.size()>0) ? Priors[0][1] : 10.;
      priors[1][0] = (Priors.size()>0) ? Priors[1][0] : 1.e-20; priors[1][1] = (Priors.size()>0) ? Priors[1][1] : 1000.;
      priors[2][0] = (Priors.size()>0) ? Priors[2][0] : 0.; priors[2][1] = (Priors.size()>0) ? Priors[2][1] : 10.;
      chi2.set_par_limits(priors);
      if (Pstart.size()>0) {
	pstart.push_back(Pstart[0]);
	pstart.push_back(Pstart[1]);
	pstart.push_back(Pstart[2]);
      } else {
	pstart.push_back(1.);
	pstart.push_back(200.);
	pstart.push_back(1.);
      }
    }
    else {
      if (Priors.size()>0 && Priors.size()!=4) ErrorMsg("Error in ModelTwoPointCorrelation::measure_fsigma8_bsigma8 of DispersionModel.cpp!");
      vector< vector<double> > priors(4,vector<double>(4)); 
      priors[0][0] = (Priors.size()>0) ? Priors[0][0] : 0.; priors[0][1] = (Priors.size()>0) ? Priors[0][1] : 10.;
      priors[1][0] = (Priors.size()>0) ? Priors[1][0] : 1.e-20; priors[1][1] = (Priors.size()>0) ? Priors[1][1] : 1000.;
      priors[2][0] = (Priors.size()>0) ? Priors[2][0] : 0.; priors[2][1] = (Priors.size()>0) ? Priors[2][1] : 10.;
      priors[3][0] = (Priors.size()>0) ? Priors[3][0] : -10.; priors[3][1] = (Priors.size()>0) ? Priors[3][1] : 10.;
      chi2.set_par_limits(priors);
      if (Pstart.size()>0) {
	pstart.push_back(Pstart[0]);
	pstart.push_back(Pstart[1]);
	pstart.push_back(Pstart[2]);
	pstart.push_back(Pstart[3]);
      } else {
	pstart.push_back(1.);
	pstart.push_back(200.);
	pstart.push_back(1.);
	pstart.push_back(0.);
      }
    }
  }

 
  // ----- estimate the best-fit values of f*sigma8, bias*sigma8 (and sigma12) -----

  cout <<"Please wait, I'm minimizing the chi^2 function..."<<endl;
  chi2.get_bestfit(pstart); 
  
  m_beta_best = chi2.bestfit(0);

  if (!NL) m_bias_sigma8_best = chi2.bestfit(1)*cosm.sigma8_Pk(author, redshift, Model);
  else {
    m_bias_sigma8_best = chi2.bestfit(2)*cosm.sigma8_Pk(author, redshift, Model);
    m_sigma12_best = chi2.bestfit(1);
    if (bias_nl) m_bA_best = chi2.bestfit(3);
  }

  m_f_sigma8_best = m_beta_best*m_bias_sigma8_best;
  

  vector<double> bb; 
  bb.push_back(m_beta_best); 
  if (NL) bb.push_back(m_sigma12_best);
  bb.push_back(m_bias_sigma8_best/cosm.sigma8_Pk(author, redshift, Model));
  if (bias_nl) bb.push_back(m_bA_best);

  m_E_min = chi2.get_chi2(bb);

}
