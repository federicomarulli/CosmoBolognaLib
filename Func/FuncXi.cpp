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
 *  @file Func/FuncXi.cpp
 *
 *  @brief Functions used to model the two-point correlation function
 *
 *  This file contains the implementation of the functions used to
 *  model the two-point correlation function
 *
 *  @author Federico Marulli 
 *
 *  @author federico.marulli3@unbo.it
 */

#include "Func.h"
using namespace cosmobl;


// =====================================================================================

/// @cond glob

double cosmobl::glob::func_xi_GSL (double kk, void *params)
{
  struct cosmobl::glob::STR_xi *pp = (struct cosmobl::glob::STR_xi *) params;
  
  double lgk = log10(kk);
  
  double lgPkK = interpolated(lgk, pp->lgkk, pp->lgPk, "Linear");
    
  double Int = pow(10.,lgPkK)*sin(kk*pp->rr)*kk/pp->rr;

  return Int * exp(-kk*kk*pp->aa*pp->aa); // eq. 24 of Anderson et al. 2012  
}


// =====================================================================================


double cosmobl::glob::func_SSM_GSL (double kk, void *params)
{
  struct cosmobl::glob::STR_SSM *pp = (struct cosmobl::glob::STR_SSM *) params;

  double fact = (pp->unit) ? 1. : pp->hh;
  double lgk = log10(kk/fact);
  
  double lgPkK = interpolated(lgk, pp->lgkk, pp->lgPk, "Linear");
  double rr = Radius(pp->mass, pp->rho);

  return pow(10.,lgPkK)*pow(TopHat_WF(kk*rr)*kk,2)/pow(fact,pp->n_spec); 
}

/// @endcond

// =====================================================================================


double cosmobl::xi_from_Pk (const double rr, const vector<double> lgkk, const vector<double> lgPk, const double k_min, const double k_max, const double aa, const bool GSL, const double prec) 
{ 
  double Int = -1.;
  
  if (GSL) { 
    int limit_size = 1000;
    gsl_integration_workspace *ww = gsl_integration_workspace_alloc (limit_size);
  
    //gsl_integration_workspace *cycle_workspace = gsl_integration_workspace_alloc (limit_size);
    //gsl_integration_qawo_table *wf = gsl_integration_qawo_table_alloc (rr, 0., GSL_INTEG_SINE, limit_size);            

    cosmobl::glob::STR_xi str;
    str.rr = rr;
    str.aa = aa;
    str.lgkk = lgkk;
    str.lgPk = lgPk;

    gsl_function Func;
    Func.function = &glob::func_xi_GSL;
    Func.params = &str;

    double error = -1.;
    gsl_integration_qag (&Func, k_min, k_max, 0., prec, limit_size, 6, ww, &Int, &error); 
    //gsl_integration_qagiu (&Func, k_min, 0., prec, limit_size, ww, &Int, &error); 
    //gsl_integration_qawf (&Func, k_min, prec, limit_size, ww, cycle_workspace,  wf, &Int, &error); 

    gsl_integration_workspace_free (ww);
  }

  else { // using Numerical libraries
    cosmobl::classfunc::func_xi func (lgkk, lgPk, rr, aa);
    Midpnt<cosmobl::classfunc::func_xi> q1(func,k_min,k_max); 
    Int = qromo(q1);
  }

  return 1./(2.*pow(par::pi,2))*Int;
} 


// =====================================================================================


double cosmobl::xi_from_Pk (const double rr, const string file, const int c1, const int c2, const double k_min, const double k_max, const double aa, const bool GSL, const double prec) 
{
  int C1 = c1-1, C2 = c2-1;

  ifstream fin (file.c_str()); checkIO (file,1); 
  
  double KK, PK, AA;
  vector<double> lgkk, lgPk;
  string line; 

  while (getline(fin,line)) {
    stringstream ss(line);
    vector<double> cc;
    while (ss>>AA) cc.push_back(AA);
    if (C1<int(cc.size()) && C2<int(cc.size())) {
      KK = cc[C1];
      PK = cc[C2];
      if (KK>0 && PK>0) {
	lgkk.push_back(log10(KK));
	lgPk.push_back(log10(PK));
      }
    }
  }
  fin.clear(); fin.close();

  return xi_from_Pk(rr, lgkk, lgPk, k_min, k_max, aa, GSL, prec); 
} 


// ============================================================================


double cosmobl::Pk_from_xi (const double kk, const vector<double> lgrr, const vector<double> lgxi, const double r_min, const double r_max) 
{
  cosmobl::classfunc::func_Pk func (lgrr, lgxi, kk);
    
  Midpnt<cosmobl::classfunc::func_Pk> q1(func, r_min, r_max); 
  double Int = qromo(q1);
  
  return 4.*par::pi*Int;
} 


// ============================================================================


double cosmobl::Pk_from_xi (const double kk, const string file, const  int c1, const int c2, const double r_min, const double r_max) 
{
  int C1 = c1-1, C2 = c2-1;

  ifstream fin (file.c_str()); checkIO (file,1); 
  
  double RR, XI, aa;
  vector<double> lgrr, lgxi;
  string line; 

  while (getline(fin,line)) {
    stringstream ss(line);
    vector<double> cc;
    while (ss>>aa) cc.push_back(aa);
    if (C1<int(cc.size()) && C2<int(cc.size())) {
      RR = cc[C1];
      XI = cc[C2];
      if (RR>0 && XI>0) {
	lgrr.push_back(log10(RR));
	lgxi.push_back(log10(XI));
      }
    }
  }
  fin.clear(); fin.close();
  
  return Pk_from_xi (kk, lgrr, lgxi, r_min, r_max) ;
} 


// ============================================================================


double cosmobl::wp (const double rp, const vector<double> rr, const vector<double> xi, const double r_max) 
{
  cosmobl::classfunc::func_wp func (rr, xi, rp);
  
  Midsql<cosmobl::classfunc::func_wp> qq(func, rp, r_max);

  return 2.*qromo(qq);
}


// ============================================================================


double cosmobl::wp (const double rp, const string file, const double r_max) 
{
  ifstream fin (file.c_str()); checkIO (file,1); 
  
  double RR, XI;
  vector<double> rr, xi;

  while (fin>>RR>>XI) {
    rr.push_back(RR);
    xi.push_back(XI);
  }
  fin.clear(); fin.close();

  return wp(rp, rr, xi, r_max);
}


// ============================================================================


double cosmobl::sigmaR (const double RR, const int corrType, const vector<double> rr, const vector<double> corr) 
{
  double sigmaR = -1;

  if (corrType==1) { // using the spherically averaged correlation function
    cosmobl::classfunc::func_sigma_xi func (rr, corr, RR); 
    Midpnt<cosmobl::classfunc::func_sigma_xi> qq(func,0.,2.*RR); 
    if (1./pow(RR,3)*qromo(qq)<0) ErrorMsg("Error in sigmaR of Func.cpp!");
    sigmaR = sqrt(1./pow(RR,3)*qromo(qq));
  }

  else if (corrType==2) { // using the projected correlation function
    cosmobl::classfunc::func_sigma_wp func (rr, corr, RR); 
    Midpnt<cosmobl::classfunc::func_sigma_wp> q1(func,0.,1.);
    Midinf<cosmobl::classfunc::func_sigma_wp> q2(func,1.,100.);
    if (1./pow(RR,3)*(qromo(q1)+qromo(q2))<0) ErrorMsg("Error in sigmaR of Func.cpp!");
    sigmaR = sqrt(1./pow(RR,3)*(qromo(q1)+qromo(q2)));
  }
  
  else ErrorMsg("Error in sigmaR of Func.cpp!");
  
  return sigmaR;
}


// ============================================================================


double cosmobl::xi_projected_powerlaw (const double rp, const double r0, const double gamma) 
{
  return rp*pow(r0/rp,gamma)*exp(gammln(0.5))*exp(gammln((gamma-1.)*0.5))/exp(gammln(gamma*0.5));
}


// ============================================================================


double cosmobl::xi_ratio (const double beta)
{ 
  return 1.+2./3.*beta+0.2*beta*beta;
}


// ============================================================================


double cosmobl::xi_ratio (const double f_sigma8, const double bias_sigma8)
{ 
  return (bias_sigma8!=0) ? 1.+2./3.*f_sigma8/bias_sigma8+0.2*pow(f_sigma8/bias_sigma8,2) : -1.e30;
}


// ============================================================================

/// @cond glob

double cosmobl::xi_ratio (double xx, shared_ptr<void> pp, vector<double> par) 
{ 
  if (par.size()==2) return xi_ratio(par[0]);
  
  else if (par.size()==3) return xi_ratio(par[0], par[1]);
  
  else { ErrorMsg("Error in xi_ratio of FuncXi.cpp!"); return 0; }
}

/// @endcond

// ============================================================================


double cosmobl::error_xi_ratio (const double beta, const double error_beta) 
{ 
  return (2./3.+0.4*beta)*error_beta;
}


// ============================================================================


double cosmobl::barred_xi_direct (const double RR, const vector<double> rr, const vector<double> xi, const double rAPP, const double r0, const double gamma) 
{ 
  vector<double> xi_;

  // power-law approximation for r<rAPP
            
  double kk = 0;
  while (rr[kk]<rAPP)                                          
    xi_.push_back(3.*pow(rr[kk++]/r0,-gamma)/(-gamma+3.));                            
  kk --;

  
  // estimate the function at r>rApp
     
  double Int = (rAPP>0) ? xi_[kk]*pow(rr[kk],3)/3. : 0.;
 
  double bin = rr[1]-rr[0];

  for (unsigned int i=kk+1; i<rr.size(); i++) {
    xi_.push_back(3./pow(rr[i],3)*(Int+xi[i]*pow(rr[i],2)*bin));
    Int += xi[i]*pow(rr[i],2)*bin;
  }

  return interpolated(RR, rr, xi_, "Linear");
}


// ============================================================================


double cosmobl::barred_xi__direct (const double RR, const vector<double> rr, const vector<double> xi, const double rAPP, const double r0, const double gamma) 
{ 
  vector<double> xi__;


  // power-low approximation for r<rAPP
            
  double kk = 0;
  while (rr[kk]<rAPP)                                          
    xi__.push_back(5.*pow(rr[kk++]/r0,-gamma)/(-gamma+5.));                            
  kk --;

  
  // estimate the function at r>rApp
 
  double Int = (rAPP>0) ? xi__[kk]*pow(rr[kk],5)/5. : 0.;
 
  double bin = rr[1]-rr[0];

  for (unsigned int i=kk+1; i<rr.size(); i++) {
    xi__.push_back(5./pow(rr[i],5)*(Int+xi[i]*pow(rr[i],4)*bin));
    Int += xi[i]*pow(rr[i],4)*bin;
  }

  return interpolated(RR, rr, xi__, "Linear");
}


// ============================================================================


double cosmobl::barred_xi_ (const double RR, const vector<double> rr, const vector<double> xi, const double rApp, const double r0, const double gamma) 
{   
  vector<double> log_r, log_xi;

  for (unsigned int i=0; i<xi.size(); i++) 
    if (xi[i]>0) {
      log_r.push_back(log10(rr[i]));
      log_xi.push_back(log10(xi[i]));
    }
  
  cosmobl::classfunc::func_xi_ func_ (log_r, log_xi, 0);

  double int_an = (RR<rApp) ?  1./((3.-gamma)*pow(r0,-gamma))*pow(RR,3.-gamma) : 1./((3.-gamma)*pow(r0,-gamma))*pow(rApp,3.-gamma);

  double int_num = (RR>rApp) ? qromb(func_,rApp,RR) : 0.;

  return 3./pow(RR,3.)*(int_an+int_num);
}


// ============================================================================


double cosmobl::barred_xi__ (const double RR, const vector<double> rr, const vector<double> xi, const double rApp, const double r0, const double gamma) 
{   
  vector<double> log_r, log_xi;

  for (unsigned int i=0; i<xi.size(); i++) 
    if (xi[i]>0) {
      log_r.push_back(log10(rr[i]));
      log_xi.push_back(log10(xi[i]));
    }

  cosmobl::classfunc::func_xi_ func__ (log_r, log_xi, 1);

  double int_an = (RR<rApp) ? 1./((5.-gamma)*pow(r0,-gamma))*pow(RR,5.-gamma) : 1./((5.-gamma)*pow(r0,-gamma))*pow(rApp,5.-gamma);

  double int_num = (RR>rApp) ? qromb(func__,rApp,RR) : 0.;
    
  return 5./pow(RR,5.)*(int_an+int_num);
}


// ============================================================================

/// @cond glob

double cosmobl::xi2D_lin_model (double rp, double pi, shared_ptr<void> pp, vector<double> par)
{ 
  if (par.size()!=2 && par.size()!=3 && par.size()!=4) {
    string Err = "Error in xi2D_lin_model! par.size() = " + conv(par.size(),par::fINT) + "!";
    ErrorMsg(Err);
  }

  double beta = par[0];  
  double bias = (par.size()==3) ? par[1] : 1;
  int index = par[par.size()-1];

  shared_ptr<cosmobl::glob::STR_xi2D_model> vec = static_pointer_cast<cosmobl::glob::STR_xi2D_model>(pp);

  
  if (vec->bias_nl) {
    if (par.size()!=4) {
      string Err = "Error in xi2D_lin_model! par.size() = " + conv(par.size(),par::fINT) + "!";
      ErrorMsg(Err);
    }
    double bA = par[3];
    double rr = sqrt(pow(rp,2)+pow(pi,2));
    bias *= b_nl(rr,bA);
  }

  double b2 = bias*bias;
  
  
  double xi_real = vec->xi_real[index]*b2;
  double xi_ = vec->xi_[index]*b2;
  double xi__ = vec->xi__[index]*b2;
 
  double xi_0 = multipole_xi0_model(beta, xi_real);
  double xi_2 = multipole_xi2_model(beta, xi_real, xi_);
  double xi_4 = multipole_xi4_model(beta, xi_real, xi_, xi__);

  return xi_0+xi_2*vec->P2[index]+xi_4*vec->P4[index]; 
}

/// @endcond

// ============================================================================


double cosmobl::xi2D_lin_model (const double beta, const double bias, const double xi_real, const double xi_, const double xi__, const double P_2, const double P_4)
{
  double bias2 = bias*bias;
  double Xi_real = xi_real * bias2;
  double Xi_ = xi_ * bias2;
  double Xi__ = xi__ * bias2;

  double xi_0 = multipole_xi0_model(beta, Xi_real);
  double xi_2 = multipole_xi2_model(beta, Xi_real, Xi_);
  double xi_4 = multipole_xi4_model(beta, Xi_real, Xi_, Xi__);

  return xi_0+xi_2*P_2+xi_4*P_4; 
}


// ============================================================================


double cosmobl::xi2D_lin_model (const double rp, const double pi, const double beta, const double bias, const vector<double> rad_real, const vector<double> xi_real, const vector<double> xi_, const vector<double> xi__, const int index, const bool bias_nl, const double bA)
{ 
  double rr = sqrt(rp*rp+pi*pi);
  double cos = pi/rr;

  double xiR = (index>-1) ? xi_real[index] : -1.;
  double xiR_ = (index>-1) ? xi_[index] : -1.;
  double xiR__ = (index>-1) ? xi__[index] : -1.;

  if (xiR<0) {
    xiR = interpolated(rr, rad_real, xi_real, "Linear");
    xiR_ = interpolated(rr, rad_real, xi_, "Linear");
    xiR__ = interpolated(rr, rad_real, xi__, "Linear");
  }

  double Bias = bias;
  if (bias_nl) Bias *= b_nl(rr, bA);
    
  double bias2 = Bias*Bias;
  xiR *= bias2;
  xiR_ *= bias2;
  xiR__ *= bias2;

  double xi_0 = multipole_xi0_model (beta, xiR);
  double xi_2 = multipole_xi2_model (beta, xiR, xiR_);
  double xi_4 = multipole_xi4_model (beta, xiR, xiR_, xiR__);

  return xi_0+xi_2*P_2(cos)+xi_4*P_4(cos);
}


// ============================================================================

/// @cond glob

double cosmobl::xi2D_model (double rp, double pi, shared_ptr<void> pp, vector<double> par)
{ 
  if (par.size()<3) {
    string Err = "Error in xi2D_model! par.size() = " + conv(par.size(),par::fINT) + "!";
    ErrorMsg(Err);
  }

  shared_ptr<cosmobl::glob::STR_xi2D_model> vec = static_pointer_cast<cosmobl::glob::STR_xi2D_model>(pp);
    
  int index = par[par.size()-1];
  int ind_min = index*vec->step_v;
  int ind_max = ind_min+vec->step_v;

  double sigma12 = par[1];

  vector<double> par2;
  par2.push_back(par[0]);
  par2.push_back(par[2]);
  for (unsigned int i=3; i<par.size(); i++) par2.push_back(par[i]);
  

  // convolution of the linear xi(rp,pi) with the velocity distribution function f(v)		

  double xi2D_nl = 0;
  double norm = 0.; // to avoid numerical problems when sigma12 is too small

  for (int i=ind_min; i<ind_max; i++) {
    par2[par2.size()-1] = i;
    xi2D_nl += xi2D_lin_model(vec->rp[i], vec->pi[i], pp, par2)*f_v(vec->vel[i], sigma12, vec->FV)*vec->delta_v;
    norm += f_v(vec->vel[i], sigma12, vec->FV)*vec->delta_v;
  }
  
  xi2D_nl /= norm;


  if (fabs(norm-1)>0.1) { 
    string Warn = "Attention! sigma12 = "+conv(sigma12,par::fDP2)+" ---> norm = " + conv(norm,par::fDP3) + ", the number of bins used for the convolution with f(v) should be increased!";
    WarningMsg(Warn);
    print(par);
  }
  
  return xi2D_nl;
}

/// @endcond

// ============================================================================


double cosmobl::xi2D_model (const double rp, const double pi, const double beta, const double bias, const double sigma12, const vector<double> rad_real, const vector<double> xi_real, const vector<double> xi_, const vector<double> xi__, const double var, const int FV, int index, const bool bias_nl, const double bA, const double v_min, const double v_max, const int step_v)
{
  double delta_v = (v_max-v_min)/step_v;

  double vel = v_min;

  double xi2D = 0.;


  for (int k=0; k<step_v; k++) {

    double pi_new = pi-vel*var;	

    double rr = sqrt(rp*rp+pi_new*pi_new);
    double cos = pi_new/rr;

    double xiR = (index>-1) ? xi_real[index] : -1.;
    double xiR_ = (index>-1) ? xi_[index] : -1.;
    double xiR__ = (index>-1) ? xi__[index] : -1.;
    
    if (xiR<0) {
      xiR = interpolated(rr, rad_real, xi_real, "Linear");
      xiR_ = interpolated(rr, rad_real, xi_, "Linear");
      xiR__ = interpolated(rr, rad_real, xi__, "Linear");
    }

    double Bias = bias;
    if (bias_nl) Bias *= b_nl(rr, bA);

    xi2D += xi2D_lin_model(beta, Bias, xiR, xiR_, xiR__, P_2(cos), P_4(cos))*f_v(vel, sigma12, FV)*delta_v;
  
    vel += delta_v;
    if (index>-1) index ++;
  }
 	
  return xi2D;
}


// ============================================================================


double cosmobl::f_v (const double vel, const double sigma12, const int FV) 
{
  if (FV==0) return 1./(sigma12*sqrt(2.))*exp(-sqrt(2.)*fabs(vel)/sigma12); // exponential
  
  else return 1./(sigma12*sqrt(par::pi))*exp(-(vel*vel)/(sigma12*sigma12)); // gaussian     
}


// ============================================================================


double cosmobl::f_v (const double vel, const double rp, const double pi, const double var, const double sigmav0, const double cmu, const double cs1, const double cs2)
{
  
  double sp = sqrt(pow(rp,2)+pow(pi-vel*var,2));
  double mup = 1./sp*(pi-vel*var);

  double sigmav = sigmav0*(1.+cmu*mup*mup)*(1.+cs1*exp(-cs2*rp*rp));

  return 1./(sigmav*sqrt(2.))*exp(-sqrt(2.)*fabs(vel)/sigmav); 
}


// ============================================================================


double cosmobl::f_star (const double xx, const double f_g, const double k_star) 
{
  double sigma_star = sqrt((4.*f_g+2.*f_g*f_g)/(k_star*k_star));
  
  return 1./(sigma_star*sqrt(par::pi))*exp(-xx*xx/(sigma_star*sigma_star));
}


// ============================================================================


double cosmobl::b_nl (const double rr, const double bA, const double bB, const double bC)
{
  double FF = 1./(1.+pow(rr/bB,bC));
  
  return pow(rr,bA*FF);
}
