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
 *  @file CatalogueAnalysis/TwoPointCorrelation/Multipoles.cpp
 *
 *  @brief Methods of the class TwoPointCorrelation used to measure
 *  the multipoles of the two-point correlation function
 *
 *  This file contains the implementation of the methods of the class
 *  TwoPointCorrelation used to measure the multipoles of the
 *  two-point correlation function
 *
 *  @author Federico Marulli
 *
 *  @author federico.marulli3@unbo.it
 */

#include "TwoPointCorrelation.h"
using namespace cosmobl;


// ============================================================================


void cosmobl::TwoPointCorrelation::measure_multipoles_xirmu (double rApprox)
{
  cout <<"I'm measuring the multipoles..."<<endl;

  if (m_xi_coslog.size()==0) ErrorMsg("Error in cosmobl::TwoPointCorrelation::measure_multipoles_xirmu of Multipoles.cpp!");

  erase_multipoles ();

  
  // -----------------------------
  // ---- eliminate null bins ----
  // -----------------------------
  
  vector<double> rad_logR, rad_linR, cos_linR1, cos_linR2;
  rad_logR = m_rad_log;
  rad_linR = m_rad_lin;
  cos_linR1 = m_cos_lin;
  cos_linR2 = m_cos_lin;

  vector< vector<double> > xi_coslogR, xi_coslinR; 
  xi_coslogR = m_xi_coslog;
  xi_coslinR = m_xi_coslin;

  SubMatrix (rad_logR, cos_linR1, xi_coslogR, -1.); 
  SubMatrix (rad_linR, cos_linR2, xi_coslinR, -1.); 


  // ------------------------------------------
  // ---- extrapolate at small separations ----
  // ------------------------------------------
  
  vector<double> rad_logE, rad_linE;
  rad_logE = rad_logR;
  rad_linE = rad_linR;

  vector< vector<double> > xi_coslogE, xi_coslinE; 
  xi_coslogE = xi_coslogR;
  xi_coslinE = xi_coslinR;

  vector<int> ll;
  int indM = -1; 
  while (rad_logE[++indM]<rApprox) ll.push_back(indM);
  
  Erase(rad_logE, ll);
  Erase_lines(xi_coslogE, ll);
  
  vector<int> lll;
  indM = -1; 
  while (rad_linE[++indM]<rApprox) lll.push_back(indM);

  Erase(rad_linE, lll);
  Erase_lines(xi_coslinE, lll);

  
  for (unsigned int i=0; i<xi_coslogR.size(); i++) 
    if (rad_logR[i]<rApprox)
      for (unsigned int j=0; j<xi_coslogR[i].size(); j++) 
	xi_coslogR[i][j] = interpolated_2D(rad_logR[i], cos_linR1[j], rad_logE, cos_linR1, xi_coslogE, "Poly", 4);
	
  for (unsigned int i=0; i<xi_coslinR.size(); i++) 
    if (rad_linR[i]<rApprox)
      for (unsigned int j=0; j<xi_coslinR[i].size(); j++) 
	xi_coslinR[i][j] = interpolated_2D(rad_linR[i], cos_linR2[j], rad_linE, cos_linR2, xi_coslinE, "Poly", 4);
     

  // ---------------------------------------------------
  // ---- measure the multipoles of the correlation ----
  // ---------------------------------------------------

  for (unsigned int i=0; i<xi_coslogR.size(); i++) {
    m_xi0_log.push_back(multipole_xi0(i, m_cos_lin, xi_coslogR));
    m_xi2_log.push_back(multipole_xi2(i, m_cos_lin, xi_coslogR));
    m_xi4_log.push_back(multipole_xi4(i, m_cos_lin, xi_coslogR));
    m_error_xi0_log.push_back(error_multipole_xi0(i, m_cos_lin, xi_coslogR));
    m_error_xi2_log.push_back(error_multipole_xi2(i, m_cos_lin, xi_coslogR));
    m_error_xi4_log.push_back(error_multipole_xi4(i, m_cos_lin, xi_coslogR));
  }

  for (unsigned int i=0; i<xi_coslinR.size(); i++) {
    m_xi0_lin.push_back(multipole_xi0(i, m_cos_lin, xi_coslinR));
    m_xi2_lin.push_back(multipole_xi2(i, m_cos_lin, xi_coslinR));
    m_xi4_lin.push_back(multipole_xi4(i, m_cos_lin, xi_coslinR));
    m_error_xi0_lin.push_back(error_multipole_xi0(i, m_cos_lin, xi_coslinR));
    m_error_xi2_lin.push_back(error_multipole_xi2(i, m_cos_lin, xi_coslinR));
    m_error_xi4_lin.push_back(error_multipole_xi4(i, m_cos_lin, xi_coslinR));
  }
  
}


// ============================================================================

 
void cosmobl::TwoPointCorrelation::measure_multipoles_xirppi (int step_cos)
{

  cout <<"I'm measuring the multipoles..."<<endl;

  if (m_xi_2d_lin.size()==0) ErrorMsg("Error in cosmobl::TwoPointCorrelation::measure_multipoles_xirppi of Multipoles.cpp!");

  erase_multipoles ();

  
  // -----------------------------
  // ---- eliminate null bins ----
  // -----------------------------

  vector<double> rad_linR1 = m_rad_lin, rad_linR2 = m_rad_lin, rad_linR3 = m_rad_lin, rad_linR4 = m_rad_lin;
  vector< vector<double> > xi_2d_linR = m_xi_2d_lin, error_xi_2d_linR = m_error_xi_2d_lin; 

  SubMatrix (rad_linR1, rad_linR2, xi_2d_linR, -1.); 
  SubMatrix (rad_linR3, rad_linR4, error_xi_2d_linR, -1.); 


  // -------------------------------------------------
  // ---- interpolate xi(rp,pi) in the r,mu plane ----
  // -------------------------------------------------
 
  double cos_min = 0.;
  double cos_max = 1.;
  vector<double> cos_lin = linear_bin_vector(step_cos, cos_min, cos_max);


  // ---- logarithimic binning ----

  vector< vector<double> > m_xi_coslog_interp(m_rad_log.size(), vector<double>(cos_lin.size(), 0.));
  vector< vector<double> > error_m_xi_coslog_interp(m_rad_log.size(), vector<double>(cos_lin.size(), 0.));

  for (unsigned int i=0; i<m_rad_log.size()-1; i++) 
    for (unsigned int j=0; j<cos_lin.size(); j++) {
      double rp = m_rad_log[i]*sqrt(1.-cos_lin[j]*cos_lin[j]);
      double pi = m_rad_log[i]*cos_lin[j];
      
      m_xi_coslog_interp[i][j] = interpolated_2D(rp, pi, rad_linR1, rad_linR2, xi_2d_linR, "Poly", 4);

      error_m_xi_coslog_interp[i][j] = interpolated_2D(rp, pi, rad_linR3, rad_linR4, error_xi_2d_linR, "Poly", 4);
    }
  
  
  // ---- linear binning ----
  
  vector< vector<double> > xi_coslin_interp(rad_linR1.size(), vector<double>(cos_lin.size(), 0.));
  vector< vector<double> > error_xi_coslin_interp(rad_linR1.size(), vector<double>(cos_lin.size(), 0.));

  for (unsigned int i=0; i<rad_linR1.size(); i++) 
    for (unsigned int j=0; j<cos_lin.size(); j++) {
      double rp = m_rad_lin[i]*sqrt(1.-cos_lin[j]*cos_lin[j]);
      double pi = m_rad_lin[i]*cos_lin[j];

      xi_coslin_interp[i][j] = interpolated_2D(rp, pi, rad_linR1, rad_linR2, xi_2d_linR, "Poly", 4);

      error_xi_coslin_interp[i][j] = interpolated_2D(rp, pi, rad_linR3, rad_linR4, error_xi_2d_linR, "Poly", 4);
    }


  // ---------------------------------------------------
  // ---- measure the multipoles of the correlation ----
  // ---------------------------------------------------

  for (unsigned int i=0; i<m_xi_coslog_interp.size(); i++) 
    if (0<m_rad_log[i] && m_rad_log[i]<m_rMAX/sqrt(2.)) {
      m_xi0_log.push_back(multipole_xi0(i, cos_lin, m_xi_coslog_interp));
      m_xi2_log.push_back(multipole_xi2(i, cos_lin, m_xi_coslog_interp));
      m_xi4_log.push_back(multipole_xi4(i, cos_lin, m_xi_coslog_interp));
      m_error_xi0_log.push_back(error_multipole_xi0(i, cos_lin, m_xi_coslog_interp));
      m_error_xi2_log.push_back(error_multipole_xi2(i, cos_lin, m_xi_coslog_interp));
      m_error_xi4_log.push_back(error_multipole_xi4(i, cos_lin, m_xi_coslog_interp));
    }

  for (unsigned int i=0; i<xi_coslin_interp.size(); i++) 
    if (0<m_rad_lin[i] && m_rad_lin[i]<m_rMAX/sqrt(2.)) {
      m_xi0_lin.push_back(multipole_xi0(i, cos_lin, xi_coslin_interp));
      m_xi2_lin.push_back(multipole_xi2(i, cos_lin, xi_coslin_interp));
      m_xi4_lin.push_back(multipole_xi4(i, cos_lin, xi_coslin_interp));
      m_error_xi0_lin.push_back(error_multipole_xi0(i, cos_lin, xi_coslin_interp));
      m_error_xi2_lin.push_back(error_multipole_xi2(i, cos_lin, xi_coslin_interp));
      m_error_xi4_lin.push_back(error_multipole_xi4(i, cos_lin, xi_coslin_interp));
    }
}


// ============================================================================


void cosmobl::TwoPointCorrelation::measure_effective_multipoles_xirppi () // check the propagated errors!
{
  cout <<"I'm measuring the multipoles..."<<endl;

  if (m_xi_2d_lin.size()==0) ErrorMsg("Error in cosmobl::TwoPointCorrelation::measure_multipoles_xirppi of Multipoles.cpp!");

  erase_multipoles ();

  double delta_s, rad;


  // -----------------------------
  // ---- eliminate null bins ----
  // -----------------------------

  vector<double> rad_linR1 = m_rad_lin, rad_linR2 = m_rad_lin, rad_linR3 = m_rad_lin, rad_linR4 = m_rad_lin;
  vector< vector<double> > xi_2d_linR = m_xi_2d_lin, error_xi_2d_linR = m_error_xi_2d_lin; 

  SubMatrix (rad_linR1, rad_linR2, xi_2d_linR, -1.); 
  SubMatrix (rad_linR3, rad_linR4, error_xi_2d_linR, -1.); 


  // ---------------------------------------------------
  // ---- measure the multipoles of the correlation ----
  // ---------------------------------------------------

  delta_s = m_linbinsz;

  for (int i=0; i<m_nlinbins; i++) {
    rad = m_rad_lin[i];
    if (0<rad && rad<m_rMAX_eff/sqrt(2.)) {
      m_xi0_lin.push_back(multipole_xi0(rad, rad_linR1, rad_linR2, xi_2d_linR, delta_s));
      m_xi2_lin.push_back(multipole_xi2(rad, rad_linR1, rad_linR2, xi_2d_linR, delta_s));
      m_xi4_lin.push_back(multipole_xi4(rad, rad_linR1, rad_linR2, xi_2d_linR, delta_s));
      m_error_xi0_lin.push_back(error_multipole_xi0(rad, rad_linR1, rad_linR2, xi_2d_linR, delta_s));
      m_error_xi2_lin.push_back(error_multipole_xi2(rad, rad_linR1, rad_linR2, xi_2d_linR, delta_s));
      m_error_xi4_lin.push_back(error_multipole_xi4(rad, rad_linR1, rad_linR2, xi_2d_linR, delta_s));
    }
  }

  for (int i=0; i<m_nlogbins; i++) {
    rad = m_rad_log[i];
    if (0<rad && rad<m_rMAX_eff/sqrt(2.)) {
      m_xi0_log.push_back(multipole_xi0(rad, rad_linR1, rad_linR2, xi_2d_linR, delta_s));
      m_xi2_log.push_back(multipole_xi2(rad, rad_linR1, rad_linR2, xi_2d_linR, delta_s));
      m_xi4_log.push_back(multipole_xi4(rad, rad_linR1, rad_linR2, xi_2d_linR, delta_s));
      m_error_xi0_log.push_back(error_multipole_xi0(rad, rad_linR1, rad_linR2, xi_2d_linR, delta_s));
      m_error_xi2_log.push_back(error_multipole_xi2(rad, rad_linR1, rad_linR2, xi_2d_linR, delta_s));
      m_error_xi4_log.push_back(error_multipole_xi4(rad, rad_linR1, rad_linR2, xi_2d_linR, delta_s));
    }
  }
  
}


// ============================================================================


void cosmobl::TwoPointCorrelation::measure_normalized_quadrupole (double rApprox)
{
  measure_multipoles_xirmu (rApprox);

  m_quad.erase(m_quad.begin(), m_quad.end());

  double Int = 0., error_Int = 0., rr, num, den, error_num, error_den, error_tot; 

  for (int i=0; i<m_nlinbins; i++) {
    rr = m_rad_lin[i];  

    num = m_xi2_lin[i];
    den = m_xi0_lin[i]-3./pow(rr,3.)*Int;

    if (den!=0) {
      m_quad.push_back(num/den);

      error_num = m_error_xi2_lin[i];
      error_den = sqrt(pow(m_error_xi0_lin[i],2.)+pow(3./pow(rr,3.)*sqrt(error_Int),2.));
      error_tot = num/den*sqrt(pow(error_num/num,2.)+pow(error_den/den,2.));
      m_error_quad.push_back(error_tot);
    }

    Int += m_xi0_lin[i]*pow(rr,2)*m_linbinsz;
    error_Int += pow(m_error_xi0_lin[i]*pow(rr,2)*m_linbinsz,2);
  }
}
