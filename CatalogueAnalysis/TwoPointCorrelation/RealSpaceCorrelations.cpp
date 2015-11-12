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
 *  @file CatalogueAnalysis/TwoPointCorrelation/RealSpaceCorrelations.cpp
 *        
 *  @brief Methods of the class TwoPointCorrelation used to estimate
 *  the real-space correlation function
 *
 *  This file contains the implementation of the methods of the class
 *  TwoPointCorrelation used to estimate the real-space correlation
 *  function, either by modelling the projected correlation function
 *  or with the deprojection method
 *
 *  @author Federico Marulli
 *
 *  @author federico.marulli3@unbo.it
 */

#include "TwoPointCorrelation.h"
using namespace cosmobl;


// ============================================================================


void cosmobl::TwoPointCorrelation::compute_real_space_xi_deprojected (double &pimax, string &dir) // following Saunders 1992
{
  cout <<"I'm computing the real space xi..."<<endl;

  //if (xi_proj.size()==0)
  measure_projected_xi (pimax,dir);

  int nlogMax = m_nlogbins;
  int nlinMax = m_nlinbins;
  
  m_rad_real_log.erase(m_rad_real_log.begin(), m_rad_real_log.end());
  m_rad_real_lin.erase(m_rad_real_lin.begin(), m_rad_real_lin.end());
  m_xi_real_log.erase(m_xi_real_log.begin(), m_xi_real_log.end());
  m_xi_real_lin.erase(m_xi_real_lin.begin(), m_xi_real_lin.end());

  for (int i=0; i<nlogMax; i++) {m_rad_real_log.push_back(0.); m_xi_real_log.push_back(0.);} 
  for (int i=0; i<nlinMax; i++) {m_rad_real_lin.push_back(0.); m_xi_real_lin.push_back(0.);} 


  // ---- deprojected xi in logarithmic bins ---- 

  vector<double> log_rr_log, log_xi_real_log;
  
  for (int i=0; i<nlogMax; i++) {
    double ri = pow(10,i*m_logbinsz+m_shift_log+log10(m_rMIN));
    m_rad_real_log[i] = ri;
    for (int j=i; j<nlogMax-1; j++) {
      double rj = pow(10,j*m_logbinsz+m_shift_log+log10(m_rMIN));
      double rj1 = pow(10,(j+1)*m_logbinsz+m_shift_log+log10(m_rMIN));
      m_xi_real_log[i] -= (m_xi_proj[j+1]-m_xi_proj[j])/(rj1-rj)*log((rj1+pow(max(0.,rj1*rj1-ri*ri),0.5))/(rj+pow(max(0.,rj*rj-ri*ri),0.5)));
    }
    m_xi_real_log[i] /= par::pi;

    if (m_xi_real_log[i]>0) { // check!!!
      log_rr_log.push_back(log10(ri));
      log_xi_real_log.push_back(log10(m_xi_real_log[i]));
    }
  }

  string file_deproj = dir+"xi_deprojected_log"; 
  ofstream fout_deproj (file_deproj.c_str()); checkIO(file_deproj,0);
  for (int i=0; i<nlogMax; i++)
    fout_deproj <<m_rad_real_log[i]<<"   "<<m_xi_real_log[i]<<"   -1"<<endl; // check!!!
  fout_deproj.clear(); fout_deproj.close();
  cout <<"I wrote the file: "<<file_deproj<<endl;
  

  // ---- deprojected xi in linear bins ---- 
  
  if (log_rr_log.size()==0) 
    ErrorMsg("Error in cosmobl::TwoPointCorrelation::compute_real_space_xi_deprojected of RealSpaceCorrelations.cpp: log_rr_log.size()=0!");

  string file_deproj2 = dir+"xi_deprojected_lin";
  ofstream fout_deproj2 (file_deproj2.c_str()); checkIO(file_deproj2,0);
  
  double log_xi_real_lin;

  for (int i=0; i<nlinMax; i++) {
    double log_rr = log10(m_rad_lin[i]);
    m_rad_real_lin[i] = m_rad_lin[i];

    log_xi_real_lin = interpolated(log_rr, log_rr_log, log_xi_real_log, "Linear", -1);
    m_xi_real_lin[i] = pow(10.,log_xi_real_lin);

    fout_deproj2 <<pow(10.,log_rr)<<"   "<<m_xi_real_lin[i]<<"   -1"<<endl; // check!!!
  }
  fout_deproj2.clear(); fout_deproj2.close();
  cout <<"I wrote the file: "<<file_deproj2<<endl;
}


// ============================================================================


void cosmobl::TwoPointCorrelation::derive_real_xi (int &XiReal, double &pimax, string &dir_realxi, double &rApp)
{
  m_rad_real_lin.erase(m_rad_real_lin.begin(), m_rad_real_lin.end());
  m_xi_real_lin.erase(m_xi_real_lin.begin(), m_xi_real_lin.end());
  m_error_xi_real_lin.erase(m_error_xi_real_lin.begin(), m_error_xi_real_lin.end());

  
  if (XiReal==0) { // deprojected correlation function
    string MK = "mkdir -p "+dir_realxi; if (system (MK.c_str())) {};
    compute_real_space_xi_deprojected (pimax, dir_realxi);
  }
  
  else if (XiReal==1) { // correlation function measured directly in real space
    string file_real = dir_realxi+"csi_lin"; cout <<"file with the xi_real data: "<<file_real<<endl;
    ifstream fin_real (file_real.c_str()); checkIO(file_real,1);
    
    double RRR, XXX, ERR;
    while (fin_real >>RRR>>XXX>>ERR) {
      m_rad_real_lin.push_back(RRR);
      m_xi_real_lin.push_back(XXX);
      m_error_xi_real_lin.push_back(ERR);
    }
    fin_real.clear(); fin_real.close();
  }

  else ErrorMsg("Error in cosmobl::TwoPointCorrelation::derive_real_xi of RealSpaceCorrelations.cpp -> XiReal must be <2 !");
  

  // ---------- extrapolate (linearly) the real space correlation function at small scales ---------- 

  m_xi_real_lin_extr.erase(m_xi_real_lin_extr.begin(), m_xi_real_lin_extr.end());

  // for the first bins we use a power-low approximation, here we define the maximum bin considered
  int kk = int((rApp-m_shift_lin)/m_linbinsz);

  // minimum bin used for the fit 
  int k_min = int((3.-m_shift_lin)/m_linbinsz);

  // maximum bin used for the fit 
  int k_max = int((10.-m_shift_lin)/m_linbinsz); 

  // compute the best power-low for the first scales 
  vector<double> lgXX, lgYY; 
  for (int i=k_min; i<k_max; i++) {
    if (m_xi_real_lin[i]>0) {
      lgXX.push_back(log10(m_rad_lin[i]));
      lgYY.push_back(log10(m_xi_real_lin[i]));
    }
  }
  VecDoub LGXX(lgXX.size()), LGYY(lgXX.size()), ERR(lgXX.size());
  for (unsigned int i=0; i<lgXX.size(); i++) {
    LGXX[i] = lgXX[i];
    LGYY[i] = lgYY[i];
    ERR[i] = 1.; // check!!!
  }
  
  Fitlin fitlin (LGXX, LGYY, ERR, linearfit); 
  fitlin.fit();  
  m_r0_linextr = pow(10.,-fitlin.a[0]/fitlin.a[1]);                                                  
  m_gamma_linextr = -fitlin.a[1];  
  
  /*
  Fitmed fitmed (LGXX,LGYY);
  m_r0_linextr = pow(10.,-fitmed.a/fitmed.b);                                                  
  m_gamma_linextr = -fitmed.b;  
  */

  for (int i=0; i<int(m_xi_real_lin.size()); i++) {
   
    if (i<kk) {
      double RR = m_rad_lin[i];
      m_xi_real_lin_extr.push_back(pow(RR/m_r0_linextr,-m_gamma_linextr));
    }
    else m_xi_real_lin_extr.push_back(m_xi_real_lin[i]);

  }

  // check!!!
  m_rad_real_lin.push_back(0);
  m_xi_real_lin.push_back(0);
  m_error_xi_real_lin.push_back(0);
  m_xi_real_lin_extr.push_back(0);
}

