/********************************************************************
 *  Copyright (C) 2010 by Federico Marulli                          *
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
 *  @file Func/FuncMultipoles.cpp
 *
 *  @brief Functions used to analyse the multipoles of the two-point
 *  correlation function
 *
 *  This file contains the implementation of the functions used to
 *  model the multipoles of the two-point correlation function
 *
 *  @author Federico Marulli 
 *
 *  @author federico.marulli3@unbo.it
 */

#include "Func.h"
using namespace cosmobl;


// =====================================================================================


double cosmobl::multipole_xi0 (int indexR, vector<double> mu, vector< vector<double> > xi) 
{
  double bin = mu[1]-mu[0];
  double xi0 = 0.;

  for (unsigned int j=0; j<xi[indexR].size(); j++) 
    xi0 += xi[indexR][j]*bin;

  return xi0;
}


// ============================================================================


double cosmobl::multipole_xi2 (int indexR, vector<double> mu, vector< vector<double> > xi) 
{
  double bin = mu[1]-mu[0];
  double xi2 = 0.;

  for (unsigned int j=0; j<xi[indexR].size(); j++) 
    xi2 += xi[indexR][j]*P_2(mu[j])*bin;

  return 5.*xi2;
}


// ============================================================================


double cosmobl::multipole_xi4 (int indexR, vector<double> mu, vector< vector<double> > xi) 
{
  double bin = mu[1]-mu[0];
  double xi4 = 0.;

  for (unsigned int j=0; j<xi[indexR].size(); j++) 
    xi4 += xi[indexR][j]*P_4(mu[j])*bin;

  return 9.*xi4;
}


// ============================================================================


double cosmobl::error_multipole_xi0 (int indexR, vector<double> mu, vector<vector<double>> error) 
{
  double bin = mu[1]-mu[0];
  double err = 0.;

  for (unsigned int j=0; j<error[indexR].size(); j++) 
    err += pow(error[indexR][j],2.);
 
  return sqrt(err)*bin;
}


// ============================================================================


double cosmobl::error_multipole_xi2 (int indexR, vector<double> mu, vector<vector<double>> error) 
{
  double bin = mu[1]-mu[0];
  double err = 0.;

  for (unsigned int j=0; j<error[indexR].size(); j++) 
    err += pow(error[indexR][j]*P_2(mu[j]),2.);
  
  return 5.*sqrt(err)*bin;
}


// ============================================================================


double cosmobl::error_multipole_xi4 (int indexR, vector<double> mu, vector<vector<double>> error) 
{
  double bin = mu[1]-mu[0];
  double err = 0.;

  for (unsigned int j=0; j<error[indexR].size(); j++) 
    err += pow(error[indexR][j]*P_4(mu[j]),2.); 
    
  return 9.*sqrt(err)*bin;
}


// ============================================================================


double cosmobl::multipole_xi0 (double &ss, vector<double> rp, vector<double> pi, vector<vector<double>> xi, double &delta_s) 
{
  double xi0 = 0.;
  int Nbin = 0;

  for (unsigned int i=0; i<rp.size(); i++)
    for (unsigned int j=0; j<pi.size(); j++) {

      double sss = sqrt(rp[i]*rp[i]+pi[j]*pi[j]);
      double mu = pi[j]/sss;

      if (ss-delta_s*0.5<sss && sss<ss+delta_s*0.5) {
	Nbin ++;
	xi0 += xi[i][j]*sqrt(1.-mu*mu);
      }
    }

  return (Nbin>0) ? 0.5*par::pi*xi0/Nbin : -1000.;
}


// ============================================================================


double cosmobl::multipole_xi2 (double &ss, vector<double> rp, vector<double> pi, vector<vector<double>> xi, double &delta_s) 
{
  double xi2 = 0.;
  int Nbin = 0;

  for (unsigned int i=0; i<rp.size(); i++)
    for (unsigned int j=0; j<pi.size(); j++) {
      
      double sss = sqrt(rp[i]*rp[i]+pi[j]*pi[j]);
      double mu = pi[j]/sss;
      
      if (ss-delta_s*0.5<sss && sss<ss+delta_s*0.5) {
	Nbin ++;
	xi2 += xi[i][j]*P_2(mu)*sqrt(1.-mu*mu);
      }
    }  
    
  return (Nbin>0) ? 0.5*par::pi*5.*xi2/Nbin : -1000.;
}


// ============================================================================


double cosmobl::multipole_xi4 (double &ss, vector<double> rp, vector<double> pi, vector<vector<double>> xi, double &delta_s) 
{
  double xi4 = 0.;
  int Nbin = 0;

  for (unsigned int i=0; i<rp.size(); i++)
    for (unsigned int j=0; j<pi.size(); j++) {
      
      double sss = sqrt(rp[i]*rp[i]+pi[j]*pi[j]);
      double mu = pi[j]/sss;
      
      if (ss-delta_s*0.5<sss && sss<ss+delta_s*0.5) {
	Nbin ++;
	xi4 += xi[i][j]*P_4(mu)*sqrt(1.-mu*mu);
      }
    }
  
  return (Nbin>0) ? 0.5*par::pi*9.*xi4/Nbin : -1000.;
}


// ============================================================================


double cosmobl::error_multipole_xi0 (double &ss, vector<double> rp, vector<double> pi, vector<vector<double>> error, double &delta_s) 
{
  double err = 0.;
  int Nbin = 0;

  for (unsigned int i=0; i<rp.size(); i++)
    for (unsigned int j=0; j<pi.size(); j++) {
      
      double sss = sqrt(rp[i]*rp[i]+pi[j]*pi[j]);
      double mu = pi[j]/sss;
      
      if (ss-delta_s*0.5<sss && sss<ss+delta_s*0.5) {
	Nbin ++;
	err += pow(error[i][j]*sqrt(1.-mu*mu),2.);
      }
    }
  
  return (Nbin>0) ? 0.5*par::pi*sqrt(err)/Nbin : -1000.;
}


// ============================================================================


double cosmobl::error_multipole_xi2 (double &ss, vector<double> rp, vector<double> pi, vector<vector<double>> error, double &delta_s) 
{
  double err = 0.;
  int Nbin = 0;

  for (unsigned int i=0; i<rp.size(); i++)
    for (unsigned int j=0; j<pi.size(); j++) {
      
      double sss = sqrt(rp[i]*rp[i]+pi[j]*pi[j]);
      double mu = pi[j]/sss;
      
      if (ss-delta_s*0.5<sss && sss<ss+delta_s*0.5) {
	Nbin ++;
	err += pow(error[i][j]*P_2(mu)*sqrt(1.-mu*mu),2.);
      }
    }
  
  return (Nbin>0) ? 0.5*par::pi*5.*sqrt(err)/Nbin : -1000.;
}


// ============================================================================


double cosmobl::error_multipole_xi4 (double &ss, vector<double> rp, vector<double> pi, vector<vector<double>> error, double &delta_s) 
{
  double err = 0.;
  int Nbin = 0;

  for (unsigned int i=0; i<rp.size(); i++)
    for (unsigned int j=0; j<pi.size(); j++) {
      
      double sss = sqrt(rp[i]*rp[i]+pi[j]*pi[j]);
      double mu = pi[j]/sss;
      
      if (ss-delta_s*0.5<sss && sss<ss+delta_s*0.5) {
	Nbin ++;
	err += pow(error[i][j]*P_4(mu)*sqrt(1.-mu*mu),2.);
      }
    }
  
  return (Nbin>0) ? 0.5*par::pi*9.*sqrt(err)/Nbin : -1000.;
}


// ============================================================================

/// @cond glob

double cosmobl::multipoles (double rr, void *pp, vector<double> par)
{ 
  int index = par[par.size()-1];
  
 
  // -----------------------------------------------------
  // ---- compute xi(rp,pi) with the dispersion model ----
  // -----------------------------------------------------

  cosmobl::glob::STR_xi2D_model &vec = *static_cast<cosmobl::glob::STR_xi2D_model *>(pp);
  
  vector< vector<double> > Xi (vec.dim, vector<double> (vec.dim,-1.e30));
 
  int index2 = 0;
  
  for (int i=0; i<vec.dim; i++) 
    for (int j=0; j<vec.dim; j++) {
      par[par.size()-1] = index2++;
      Xi[i][j] = xi2D_model(vec.rp[i],vec.pi[j],pp,par);
    }

  /*
  // -----------------------------
  // ---- eliminate null bins ----
  // -----------------------------
 
  vector< vector<double> > XiR = Xi;
  vector<double> rpp, pii;

  for (unsigned int i=0; i<vec.rp.size(); i++) {
    rpp.push_back(vec.rp[i]); 
    pii.push_back(vec.pi[i]); 
  }

  SubMatrix (rpp, pii, XiR, -1); 
  */

  
  // -------------------------------------------------
  // ---- interpolate xi(rp,pi) in the r,mu plane ----
  // -------------------------------------------------
  
  double cos_min = 0.;
  double cos_max = 1.;
  int step_cos = 3000;
  vector<double> cos_lin = linear_bin_vector(step_cos, cos_min, cos_max);

  double rp, pi;
  vector< vector<double> > xi_cos(1);

  for (unsigned int i=0; i<cos_lin.size(); i++) {
    rp = rr*sqrt(1.-cos_lin[i]*cos_lin[i]);
    pi = rr*cos_lin[i];
    xi_cos[0].push_back(interpolated_2D(rp, pi, vec.rp, vec.pi, Xi/*XiR*/, "Linear", -1));
  }


  // --------------------------------------------------
  // ---- measure the multipole of the correlation ----
  // --------------------------------------------------
  
  if (vec.type[index]==1) return multipole_xi0(0,cos_lin,xi_cos); 
  else if (vec.type[index]==2) return multipole_xi2(0,cos_lin,xi_cos); 
  else { ErrorMsg ("Error in the function multipoles of FuncMultipoles.cpp!"); return 0; } 
}

/// @endcond

// ============================================================================


double cosmobl::multipole_xi0_model (double beta, double xi_real) 
{ 
  return xi_ratio(beta)*xi_real;
}


// ============================================================================


double cosmobl::multipole_xi0_model (double f_sigma8, double bias_sigma8, double sigma8z, double xi_DM) 
{ 
  return xi_ratio(f_sigma8,bias_sigma8)*xi_DM*pow(bias_sigma8/sigma8z,2);
}


// ============================================================================

/// @cond glob

double cosmobl::multipole_xi0_model (__attribute__((unused)) double xx, void *pp, vector<double> par) 
{ 
  // xx is not used!

  cosmobl::glob::STR_xi0_model &vec = *static_cast<cosmobl::glob::STR_xi0_model *>(pp);

  if (par.size()==2) return multipole_xi0_model (par[0], vec.bias_sigma8, vec.sigma8z, vec.xi_DM[par[par.size()-1]]); 

  else { ErrorMsg("Error in multipole_xi0_model of FuncMultipoles.cpp!"); return 0; }
}

/// @endcond

// ============================================================================


double cosmobl::multipole_xi2_model (double beta, double xi_real, double xi_) 
{ 
  return (4./3.*beta+4./7.*beta*beta)*(xi_real-xi_);
}


// ============================================================================


double cosmobl::multipole_xi4_model (double beta, double xi_real, double xi_, double xi__) 
{ 
  return (8./35.*beta*beta)*(xi_real+2.5*xi_-3.5*xi__);
}
