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
 *  @file CatalogueAnalysis/TwoPointCorrelation/Init.cpp
 *
 *  @brief Methods of the class TwoPointCorrelation used to initialize
 *  the private members
 *
 *  This file contains the implementation of the methods of the class
 *  TwoPointCorrelation used to initialize the private members
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#include "TwoPointCorrelation.h"
using namespace cosmobl;

	
// ============================================================================


cosmobl::TwoPointCorrelation::TwoPointCorrelation (shared_ptr<Catalogue> data, shared_ptr<Catalogue> random, bool do_all) 
  : m_data(move(data)), m_random(move(random)), m_read_only(0), m_do_all(do_all) 
{
  m_nGal = m_data->weightedN(); 
  m_nRan = m_random->weightedN();
}


// ============================================================================


void cosmobl::TwoPointCorrelation::setParameters (double &rMIN, double &rMAX, double &logbinSize, double &binSize, double &cosSize, bool ANG)
{ 
  if (rMIN<1.e-30) ErrorMsg("Error in cosmobl::TwoPointCorrelation::setParameters of Init.cpp: rMIN must be >0!");


  if (ANG) {
    
    m_thetaMIN = rMIN*par::pi/180.;
    m_thetaMAX = rMAX*par::pi/180.; 

   
    // logarithmic binning

    m_logthetabinsz = logbinSize;                                	           // size of the logarithmic bins
    m_nlogthetabins = nint((log10(m_thetaMAX)-log10(m_thetaMIN))/m_logthetabinsz); // number of logarithmic bins
    m_shift_theta_log = m_logthetabinsz*0.5;                                       // position in the bin where the w(theta) is stored

    
    // linear binning
    
    m_linthetabinsz = binSize*par::pi/180.;        	             // size of the linear bins                  
    m_nlinthetabins = nint((m_thetaMAX-m_thetaMIN)/m_linthetabinsz); // number of bins
    m_shift_theta_lin = m_linthetabinsz*0.5;                         // position in the bin where the w(theta) is stored

    
    // effictive maximum separation to be used to count the pairs 
    
    m_thetaMAX_eff = max(pow(10.,(m_nlogthetabins-0.5)*m_logthetabinsz+m_shift_theta_log+log10(m_thetaMIN)),(m_nlinthetabins-0.5)*m_linthetabinsz+m_shift_theta_lin);
    m_thetaMAX_eff = max(m_thetaMAX_eff,m_thetaMAX);
    
  }
   
  else {
  
    m_rMIN = rMIN;
    m_rMAX = rMAX; 
    m_logbinSize = logbinSize;
    m_binSize = binSize;
    m_cosSize = cosSize;


    // logarithmic binning

    m_logbinsz = logbinSize;                                 // size of the logarithmic bins
    m_nlogbins = nint((log10(rMAX)-log10(rMIN))/m_logbinsz); // number of logarithmic bins
    m_shift_log = m_logbinsz*0.5;                            // position in the bin where the \xi is stored
   

    // linear binning

    m_linbinsz = binSize;                      // size of the linear bins                  
    m_nlinbins = nint((rMAX-rMIN)/m_linbinsz); // number of bins
    m_shift_lin = binSize*0.5;                 // position in the bin where the \xi is stored


    // cos binning

    m_cosbinsz = cosSize;
    m_ncosbins = nint(1./cosSize);
    m_shift_cos = m_cosbinsz*0.5;     


    // effictive maximum separation to be used to count the pairs 
    
    m_rMAX_eff = max(pow(10.,(m_nlogbins-0.5)*m_logbinsz+m_shift_log+log10(rMIN)),(m_nlinbins-0.5)*m_linbinsz+m_shift_lin);
    m_rMAX_eff = max(m_rMAX_eff,m_rMAX);

  }

}


// ============================================================================


void cosmobl::TwoPointCorrelation::allocate_vectors_xi (bool &doGR)
{
  m_gg_log.erase(m_gg_log.begin(),m_gg_log.end());
  m_rr_log.erase(m_rr_log.begin(),m_rr_log.end());
  if (doGR) m_gr_log.erase(m_gr_log.begin(),m_gr_log.end());
  m_rad_log.erase(m_rad_log.begin(),m_rad_log.end());
  m_xi_log.erase(m_xi_log.begin(),m_xi_log.end());
  m_error_xi_log.erase(m_error_xi_log.begin(),m_error_xi_log.end());
  m_gg_lin.erase(m_gg_lin.begin(),m_gg_lin.end());
  m_rr_lin.erase(m_rr_lin.begin(),m_rr_lin.end());
  if (doGR) m_gr_lin.erase(m_gr_lin.begin(),m_gr_lin.end());
  m_rad_lin.erase(m_rad_lin.begin(),m_rad_lin.end());
  m_xi_lin.erase(m_xi_lin.begin(),m_xi_lin.end());
  m_error_xi_lin.erase(m_error_xi_lin.begin(),m_error_xi_lin.end());
  m_gg_2d.erase(m_gg_2d.begin(),m_gg_2d.end());
  m_rr_2d.erase(m_rr_2d.begin(),m_rr_2d.end());
  if (doGR) m_gr_2d.erase(m_gr_2d.begin(),m_gr_2d.end());
  m_xi_2d_lin.erase(m_xi_2d_lin.begin(),m_xi_2d_lin.end());
  m_error_xi_2d_lin.erase(m_error_xi_2d_lin.begin(),m_error_xi_2d_lin.end());
  m_gg_slog.erase(m_gg_slog.begin(),m_gg_slog.end());
  m_rr_slog.erase(m_rr_slog.begin(),m_rr_slog.end());
  if (doGR) m_gr_slog.erase(m_gr_slog.begin(),m_gr_slog.end());
  m_xi_2d_loglin.erase(m_xi_2d_loglin.begin(),m_xi_2d_loglin.end());
  m_error_xi_2d_loglin.erase(m_error_xi_2d_loglin.begin(),m_error_xi_2d_loglin.end());
  m_cos_lin.erase(m_cos_lin.begin(),m_cos_lin.end());
  m_gg_coslog.erase(m_gg_coslog.begin(),m_gg_coslog.end());
  m_rr_coslog.erase(m_rr_coslog.begin(),m_rr_coslog.end());
  if (doGR) m_gr_coslog.erase(m_gr_coslog.begin(),m_gr_coslog.end());
  m_gg_coslin.erase(m_gg_coslin.begin(),m_gg_coslin.end());
  m_rr_coslin.erase(m_rr_coslin.begin(),m_rr_coslin.end());
  if (doGR) m_gr_coslin.erase(m_gr_coslin.begin(),m_gr_coslin.end());
  m_xi_coslin.erase(m_xi_coslin.begin(),m_xi_coslin.end());
  m_error_xi_coslin.erase(m_error_xi_coslin.begin(),m_error_xi_coslin.end());
  m_xi_coslog.erase(m_xi_coslog.begin(),m_xi_coslog.end());
  m_error_xi_coslog.erase(m_error_xi_coslog.begin(),m_error_xi_coslog.end());
  
  for (int i=0; i<m_nlogbins+1; i++) { // "+1" is to avoid some "if" (i.e. it's to improve the performance)
    m_gg_log.push_back(0.);
    m_rr_log.push_back(0.);
    if (doGR) m_gr_log.push_back(0.);
    m_rad_log.push_back(-1.e30);
    m_xi_log.push_back(-1.e30);
    m_error_xi_log.push_back(-1.e30);
  }
  for (int i=0; i<m_nlinbins+1; i++) {
    m_gg_lin.push_back(0.);
    m_rr_lin.push_back(0.);
    if (doGR) m_gr_lin.push_back(0.);
    m_rad_lin.push_back(-1.e30);
    m_xi_lin.push_back(-1.e30);
    m_error_xi_lin.push_back(-1.e30);
  }
  for (int i=0; i<m_ncosbins+1; i++) 
    m_cos_lin.push_back(-1.e30);
 
  vector<double> vv1 (m_nlinbins+1,0.);
  vector<double> vv2 (m_ncosbins+1,0.);
  vector<double> vvv1 (m_nlinbins+1,-1.e30);
  vector<double> vvv2 (m_ncosbins+1,-1.e30);
  for (int i=0; i<m_nlinbins+1; i++) {
    m_gg_2d.push_back(vv1);
    m_rr_2d.push_back(vv1);
    if (doGR) m_gr_2d.push_back(vv1);
    m_xi_2d_lin.push_back(vvv1);
    m_error_xi_2d_lin.push_back(vvv1);
    m_gg_coslin.push_back(vv2);
    m_rr_coslin.push_back(vv2);
    if (doGR) m_gr_coslin.push_back(vv2);
    m_xi_coslin.push_back(vvv2);
    m_error_xi_coslin.push_back(vvv2);
  }
  for (int i=0; i<m_nlogbins+1; i++) {
    m_gg_slog.push_back(vv1);
    m_rr_slog.push_back(vv1);
    if (doGR) m_gr_slog.push_back(vv1);
    m_xi_2d_loglin.push_back(vvv1);
    m_error_xi_2d_loglin.push_back(vvv1);
    m_gg_coslog.push_back(vv2);
    m_rr_coslog.push_back(vv2);
    if (doGR) m_gr_coslog.push_back(vv2);
    m_xi_coslog.push_back(vvv2);
    m_error_xi_coslog.push_back(vvv2);
  }
}


// ============================================================================


void cosmobl::TwoPointCorrelation::allocate_vectors_ACF (bool &doGR)
{
  m_gg_theta_log.erase(m_gg_theta_log.begin(),m_gg_theta_log.end());
  m_rr_theta_log.erase(m_rr_theta_log.begin(),m_rr_theta_log.end());
  if (doGR) m_gr_theta_log.erase(m_gr_theta_log.begin(),m_gr_theta_log.end());
  m_theta_log.erase(m_theta_log.begin(),m_theta_log.end());
  m_wtheta_log.erase(m_wtheta_log.begin(),m_wtheta_log.end());
  m_error_wtheta_log.erase(m_error_wtheta_log.begin(),m_error_wtheta_log.end());

  m_gg_theta_lin.erase(m_gg_theta_lin.begin(),m_gg_theta_lin.end());
  m_rr_theta_lin.erase(m_rr_theta_lin.begin(),m_rr_theta_lin.end());
  if (doGR) m_gr_theta_lin.erase(m_gr_theta_lin.begin(),m_gr_theta_lin.end());
  m_theta_lin.erase(m_theta_lin.begin(),m_theta_lin.end());
  m_wtheta_lin.erase(m_wtheta_lin.begin(),m_wtheta_lin.end());
  m_error_wtheta_lin.erase(m_error_wtheta_lin.begin(),m_error_wtheta_lin.end());
  
  for (int i=0; i<m_nlogthetabins+1; i++) { // "+1" is to avoid some "if" (i.e. it's to improve the performance)
    m_gg_theta_log.push_back(0.);
    m_rr_theta_log.push_back(0.);
    if (doGR) m_gr_theta_log.push_back(0.);
    m_theta_log.push_back(-1.e30);
    m_wtheta_log.push_back(-1.e30);
    m_error_wtheta_log.push_back(-1.e30);
  }

  for (int i=0; i<m_nlinthetabins+1; i++) {
    m_gg_theta_lin.push_back(0.);
    m_rr_theta_lin.push_back(0.);
    if (doGR) m_gr_theta_lin.push_back(0.);
    m_theta_lin.push_back(-1.e30);
    m_wtheta_lin.push_back(-1.e30);
    m_error_wtheta_lin.push_back(-1.e30);
  }
  
}


// ============================================================================


void cosmobl::TwoPointCorrelation::erase_multipoles ()
{
  m_xi0_log.erase(m_xi0_log.begin(),m_xi0_log.end());
  m_xi2_log.erase(m_xi2_log.begin(),m_xi2_log.end());
  m_xi4_log.erase(m_xi4_log.begin(),m_xi4_log.end());
  m_xi0_lin.erase(m_xi0_lin.begin(),m_xi0_lin.end());
  m_xi2_lin.erase(m_xi2_lin.begin(),m_xi2_lin.end());
  m_xi4_lin.erase(m_xi4_lin.begin(),m_xi4_lin.end());
  m_error_xi0_log.erase(m_error_xi0_log.begin(),m_error_xi0_log.end());
  m_error_xi2_log.erase(m_error_xi2_log.begin(),m_error_xi2_log.end());
  m_error_xi4_log.erase(m_error_xi4_log.begin(),m_error_xi4_log.end());
  m_error_xi0_lin.erase(m_error_xi0_lin.begin(),m_error_xi0_lin.end());
  m_error_xi2_lin.erase(m_error_xi2_lin.begin(),m_error_xi2_lin.end());
  m_error_xi4_lin.erase(m_error_xi4_lin.begin(),m_error_xi4_lin.end());
}


// ============================================================================


void cosmobl::TwoPointCorrelation::set_rad_log (vector<double> vect) 
{
  m_rad_log.erase(m_rad_log.begin(),m_rad_log.end());
  for (unsigned int i=0; i<vect.size(); i++) m_rad_log.push_back(vect[i]);
}


// ============================================================================


void cosmobl::TwoPointCorrelation::set_xi_log (vector<double> vect) 
{
  m_xi_log.erase(m_xi_log.begin(),m_xi_log.end());
  for (unsigned int i=0; i<vect.size(); i++) m_xi_log.push_back(vect[i]);
}


// ============================================================================


void cosmobl::TwoPointCorrelation::set_error_xi_log (vector<double> vect) 
{
  m_error_xi_log.erase(m_error_xi_log.begin(),m_error_xi_log.end());
  for (unsigned int i=0; i<vect.size(); i++) m_error_xi_log.push_back(vect[i]);
}
  

// ============================================================================


void cosmobl::TwoPointCorrelation::set_rad_lin (vector<double> vect) 
{
  m_rad_lin.erase(m_rad_lin.begin(),m_rad_lin.end());
  for (unsigned int i=0; i<vect.size(); i++) m_rad_lin.push_back(vect[i]);
}


// ============================================================================


void cosmobl::TwoPointCorrelation::set_xi_lin (vector<double> vect) 
{
  m_xi_lin.erase(m_xi_lin.begin(),m_xi_lin.end());
  for (unsigned int i=0; i<vect.size(); i++) m_xi_lin.push_back(vect[i]);
}


// ============================================================================


void cosmobl::TwoPointCorrelation::set_error_xi_lin (vector<double> vect) 
{
  m_error_xi_lin.erase(m_error_xi_lin.begin(),m_error_xi_lin.end());
  for (unsigned int i=0; i<vect.size(); i++) m_error_xi_lin.push_back(vect[i]);
}


// ============================================================================


void cosmobl::TwoPointCorrelation::set_error_xi_2d_lin (vector< vector<double> > vect) 
{
  m_error_xi_2d_lin.erase(m_error_xi_2d_lin.begin(),m_error_xi_2d_lin.end());
  vector<double> vvv (m_nlinbins+1,-1.e30);
  for (int i=0; i<m_nlinbins+1; i++) m_error_xi_2d_lin.push_back(vvv);
 
  for (int i=0; i<m_nlinbins; i++)
    for (int j=0; j<m_nlinbins; j++) 
      m_error_xi_2d_lin[i][j] = vect[i][j];
}


// ============================================================================


void cosmobl::TwoPointCorrelation::set_error_xi_2d_loglin (vector< vector<double> > vect) 
{
  m_error_xi_2d_loglin.erase(m_error_xi_2d_loglin.begin(),m_error_xi_2d_loglin.end());
  vector<double> vvv (m_nlinbins+1,-1.e30);
  for (int i=0; i<m_nlogbins+1; i++) m_error_xi_2d_loglin.push_back(vvv);
  
  for (int i=0; i<m_nlogbins; i++)
    for (int j=0; j<m_nlinbins; j++) 
      m_error_xi_2d_loglin[i][j] = vect[i][j];
}

void cosmobl::TwoPointCorrelation::set_error_xi_coslog (vector< vector<double> > vect) 
{
  m_error_xi_coslog.erase(m_error_xi_coslog.begin(),m_error_xi_coslog.end());
  vector<double> vvv (m_ncosbins+1,-1.e30);
  for (int i=0; i<m_nlogbins+1; i++) m_error_xi_coslog.push_back(vvv);

  for (int i=0; i<m_nlogbins; i++)
    for (int j=0; j<m_ncosbins; j++) 
      m_error_xi_coslog[i][j] = vect[i][j];
}


// ============================================================================


void cosmobl::TwoPointCorrelation::set_error_xi_coslin (vector< vector<double> > vect) 
{
  m_error_xi_coslin.erase(m_error_xi_coslin.begin(),m_error_xi_coslin.end());
  vector<double> vvv (m_ncosbins+1,-1.e30);
  for (int i=0; i<m_nlogbins+1; i++) m_error_xi_coslin.push_back(vvv);

  for (int i=0; i<m_nlogbins; i++)
    for (int j=0; j<m_ncosbins; j++) 
      m_error_xi_coslin[i][j] = vect[i][j];
}


// ============================================================================


void cosmobl::TwoPointCorrelation::set_error_proj (vector<double> vect) 
{
  m_error_xi_proj.erase(m_error_xi_proj.begin(),m_error_xi_proj.end());
  for (unsigned int i=0; i<vect.size(); i++) m_error_xi_proj.push_back(vect[i]);
}


// ============================================================================


void cosmobl::TwoPointCorrelation::set_correlations (vector<double> xi_log, vector<double> xi_lin, vector< vector<double> > xi_2d_lin, vector< vector<double> > xi_2d_loglin, vector< vector<double> > xi_coslog, vector< vector<double> > xi_coslin)
{
  if (m_xi_log.size()!=xi_log.size() || m_xi_lin.size()!=xi_lin.size() || m_xi_2d_lin.size()!=xi_2d_lin.size() || m_xi_2d_loglin.size()!=xi_2d_loglin.size() || m_xi_coslog.size()!=xi_coslog.size() || m_xi_coslin.size()!=xi_coslin.size())
    ErrorMsg("Error in cosmobl::TwoPointCorrelation::set_correlations of Init.cpp");
    
  for (unsigned int i=0; i<m_xi_2d_lin.size(); i++) 
    for (unsigned int j=0; j<m_xi_2d_lin[i].size(); j++) 
      if (m_xi_2d_lin[i].size()!=m_xi_2d_lin[i].size()) 
	ErrorMsg("Error in cosmobl::TwoPointCorrelation::set_correlations of Init.cpp");
   
  for (unsigned int i=0; i<m_xi_2d_loglin.size(); i++) 
    for (unsigned int j=0; j<m_xi_2d_loglin[i].size(); j++) 
      if (m_xi_2d_loglin[i].size()!=m_xi_2d_loglin[i].size()) 
	ErrorMsg("Error in cosmobl::TwoPointCorrelation::set_correlations of Init.cpp");

  for (unsigned int i=0; i<m_xi_coslog.size(); i++) 
    for (unsigned int j=0; j<m_xi_coslog[i].size(); j++) 
      if (m_xi_coslog[i].size()!=m_xi_coslog[i].size()) 
	ErrorMsg("Error in cosmobl::TwoPointCorrelation::set_correlations of Init.cpp");

  for (unsigned int i=0; i<m_xi_coslin.size(); i++) 
    for (unsigned int j=0; j<m_xi_coslin[i].size(); j++) 
      if (m_xi_coslin[i].size()!=m_xi_coslin[i].size()) 
	ErrorMsg("Error in cosmobl::TwoPointCorrelation::set_correlations of Init.cpp");

  m_xi_log = xi_log;
  m_xi_lin = xi_lin;
  m_xi_2d_lin = xi_2d_lin;
  m_xi_2d_loglin = xi_2d_loglin;  
  m_xi_coslin = xi_coslin;
  m_xi_coslog = xi_coslog;
 
}


// ============================================================================


void cosmobl::TwoPointCorrelation::set_xi_real_lin (vector<double> vect)
{
  m_xi_real_lin.erase(m_xi_real_lin.begin(),m_xi_real_lin.end());
  for (unsigned int i=0; i<vect.size(); i++) m_xi_real_lin.push_back(vect[i]);
}  


// ============================================================================


void cosmobl::TwoPointCorrelation::set_xi_real_log (vector<double> vect)
{
  m_xi_real_log.erase(m_xi_real_log.begin(),m_xi_real_log.end());
  for (unsigned int i=0; i<vect.size(); i++) m_xi_real_log.push_back(vect[i]);
}  


// ============================================================================


void cosmobl::TwoPointCorrelation::set_xi_real_lin_interp (vector<double> vect)
{
  m_xi_real_lin_interp.erase(m_xi_real_lin_interp.begin(),m_xi_real_lin_interp.end());
  for (unsigned int i=0; i<vect.size(); i++) m_xi_real_lin_interp.push_back(vect[i]);
}  


// ============================================================================


void cosmobl::TwoPointCorrelation::set_xi_real_lin_extr (vector<double> vect)
{
  m_xi_real_lin_extr.erase(m_xi_real_lin_extr.begin(),m_xi_real_lin_extr.end());
  for (unsigned int i=0; i<vect.size(); i++) m_xi_real_lin_extr.push_back(vect[i]);
}  


// ============================================================================


void cosmobl::TwoPointCorrelation::set_xi_proj (vector<double> rp_proj, vector<double> xi_proj, vector<double> error_xi_proj)
{
  if (rp_proj.size()!=xi_proj.size() || rp_proj.size()!=error_xi_proj.size())
    ErrorMsg("Error in cosmobl::TwoPointCorrelation::set_xi_proj of Init.cpp");

  m_rp_proj.erase(m_rp_proj.begin(),m_rp_proj.end());
  m_xi_proj.erase(m_xi_proj.begin(),m_xi_proj.end());
  m_error_xi_proj.erase(m_error_xi_proj.begin(),m_error_xi_proj.end());
  for (unsigned int i=0; i<rp_proj.size(); i++) {
    m_rp_proj.push_back(rp_proj[i]);
    m_xi_proj.push_back(xi_proj[i]);
    m_error_xi_proj.push_back(error_xi_proj[i]);
  }
}


// ============================================================================


void cosmobl::TwoPointCorrelation::set_gg_pairs (vector<double> gg_log, vector<double> gg_lin, vector< vector<double> > gg_2d, vector< vector<double> > gg_slog, vector< vector<double> > gg_coslog, vector< vector<double> > gg_coslin)
{
  if (m_gg_log.size()!=gg_log.size() || m_gg_lin.size()!=gg_lin.size() || m_gg_2d.size()!=gg_2d.size() || m_gg_slog.size()!=gg_slog.size() || m_gg_coslog.size()!=gg_coslog.size() || m_gg_coslin.size()!=gg_coslin.size()){
    cout << m_gg_log.size() << " " << gg_log.size() << endl;
    cout << m_gg_lin.size() << " " << gg_lin.size() << endl;
    cout << m_gg_2d.size() << " " << gg_2d.size() << endl;
    cout << m_gg_slog.size() << " " << gg_slog.size() << endl;
    cout << m_gg_coslog.size() << " " << gg_coslog.size() << endl;
    cout << m_gg_coslin.size() << " " << gg_coslin.size() << endl;
    ErrorMsg("Error in cosmobl::TwoPointCorrelation::set_gg_pairs of Init.cpp 1");
  }
    
  
  for (unsigned int i=0; i<m_gg_2d.size(); i++) 
    for (unsigned int j=0; j<m_gg_2d[i].size(); j++) 
      if (m_gg_2d[i].size()!=m_gg_2d[i].size()) 
	ErrorMsg("Error in cosmobl::TwoPointCorrelation::set_gg_pairs of Init.cpp 2");
   
  for (unsigned int i=0; i<m_gg_slog.size(); i++) 
    for (unsigned int j=0; j<m_gg_slog[i].size(); j++) 
      if (m_gg_slog[i].size()!=m_gg_slog[i].size()) 
	ErrorMsg("Error in cosmobl::TwoPointCorrelation::set_gg_pairs of Init.cpp 3");

  for (unsigned int i=0; i<m_gg_coslog.size(); i++) 
    for (unsigned int j=0; j<m_gg_coslog[i].size(); j++) 
      if (m_gg_coslog[i].size()!=m_gg_coslog[i].size()) 
	ErrorMsg("Error in cosmobl::TwoPointCorrelation::set_gg_pairs of Init.cpp 4");

  for (unsigned int i=0; i<m_gg_coslin.size(); i++) 
    for (unsigned int j=0; j<m_gg_coslin[i].size(); j++) 
      if (m_gg_coslin[i].size()!=m_gg_coslin[i].size()) 
	ErrorMsg("Error in cosmobl::TwoPointCorrelation::set_gg_pairs of Init.cpp 5");

  m_gg_log = gg_log;
  m_gg_lin = gg_lin;
  m_gg_2d = gg_2d;
  m_gg_slog = gg_slog;  
  m_gg_coslin = gg_coslin;
  m_gg_coslog = gg_coslog;
 
}


// ============================================================================


void cosmobl::TwoPointCorrelation::set_rr_pairs (vector<double> rr_log, vector<double> rr_lin, vector< vector<double> > rr_2d, vector< vector<double> > rr_slog, vector< vector<double> > rr_coslog, vector< vector<double> > rr_coslin)
{
  if (m_rr_log.size()!=rr_log.size() || m_rr_lin.size()!=rr_lin.size() || m_rr_2d.size()!=rr_2d.size() || m_rr_slog.size()!=rr_slog.size() || m_rr_coslog.size()!=rr_coslog.size() || m_rr_coslin.size()!=rr_coslin.size())
    ErrorMsg("Error in cosmobl::TwoPointCorrelation::set_rr_pairs of Init.cpp");
    
  
  for (unsigned int i=0; i<m_rr_2d.size(); i++) 
    for (unsigned int j=0; j<m_rr_2d[i].size(); j++) 
      if (m_rr_2d[i].size()!=m_rr_2d[i].size()) 
	ErrorMsg("Error in cosmobl::TwoPointCorrelation::set_rr_pairs of Init.cpp");
   
  for (unsigned int i=0; i<m_rr_slog.size(); i++) 
    for (unsigned int j=0; j<m_rr_slog[i].size(); j++) 
      if (m_rr_slog[i].size()!=m_rr_slog[i].size()) 
	ErrorMsg("Error in cosmobl::TwoPointCorrelation::set_rr_pairs of Init.cpp");

  for (unsigned int i=0; i<m_rr_coslog.size(); i++) 
    for (unsigned int j=0; j<m_rr_coslog[i].size(); j++) 
      if (m_rr_coslog[i].size()!=m_rr_coslog[i].size()) 
	ErrorMsg("Error in cosmobl::TwoPointCorrelation::set_rr_pairs of Init.cpp");

  for (unsigned int i=0; i<m_rr_coslin.size(); i++) 
    for (unsigned int j=0; j<m_rr_coslin[i].size(); j++) 
      if (m_rr_coslin[i].size()!=m_rr_coslin[i].size()) 
	ErrorMsg("Error in cosmobl::TwoPointCorrelation::set_rr_pairs of Init.cpp");

  m_rr_log = rr_log;
  m_rr_lin = rr_lin;
  m_rr_2d = rr_2d;
  m_rr_slog = rr_slog;  
  m_rr_coslin = rr_coslin;
  m_rr_coslog = rr_coslog;
 
}


// ============================================================================


void cosmobl::TwoPointCorrelation::set_gr_pairs (vector<double> gr_log, vector<double> gr_lin, vector< vector<double> > gr_2d, vector< vector<double> > gr_slog, vector< vector<double> > gr_coslog, vector< vector<double> > gr_coslin)
{
  if (m_gr_log.size()!=gr_log.size() || m_gr_lin.size()!=gr_lin.size() || m_gr_2d.size()!=gr_2d.size() || m_gr_slog.size()!=gr_slog.size() || m_gr_coslog.size()!=gr_coslog.size() || m_gr_coslin.size()!=gr_coslin.size())
    ErrorMsg("Error in cosmobl::TwoPointCorrelation::set_gr_pairs of Init.cpp");
  
  for (unsigned int i=0; i<m_gr_2d.size(); i++) 
    for (unsigned int j=0; j<m_gr_2d[i].size(); j++) 
      if (m_gr_2d[i].size()!=m_gr_2d[i].size()) 
	ErrorMsg("Egror in cosmobl::TwoPointCorrelation::set_gr_pairs of Init.cpp");
   
  for (unsigned int i=0; i<m_gr_slog.size(); i++) 
    for (unsigned int j=0; j<m_gr_slog[i].size(); j++) 
      if (m_gr_slog[i].size()!=m_gr_slog[i].size()) 
	ErrorMsg("Egror in cosmobl::TwoPointCorrelation::set_gr_pairs of Init.cpp");

  for (unsigned int i=0; i<m_gr_coslog.size(); i++) 
    for (unsigned int j=0; j<m_gr_coslog[i].size(); j++) 
      if (m_gr_coslog[i].size()!=m_gr_coslog[i].size()) 
	ErrorMsg("Egror in cosmobl::TwoPointCorrelation::set_gr_pairs of Init.cpp");

  for (unsigned int i=0; i<m_gr_coslin.size(); i++) 
    for (unsigned int j=0; j<m_gr_coslin[i].size(); j++) 
      if (m_gr_coslin[i].size()!=m_gr_coslin[i].size()) 
	ErrorMsg("Egror in cosmobl::TwoPointCorrelation::set_gr_pairs of Init.cpp");

  m_gr_log = gr_log;
  m_gr_lin = gr_lin;
  m_gr_2d = gr_2d;
  m_gr_slog = gr_slog;  
  m_gr_coslin = gr_coslin;
  m_gr_coslog = gr_coslog;
 
}
