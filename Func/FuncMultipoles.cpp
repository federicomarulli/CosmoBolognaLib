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


double cosmobl::multipole_xi0 (const int indexR, const vector<double> mu, const vector< vector<double> > xi) 
{
  double bin = mu[1]-mu[0];
  double xi0 = 0.;

  for (unsigned int j=0; j<xi[indexR].size(); j++) 
    xi0 += xi[indexR][j]*bin;

  return xi0;
}


// ============================================================================


double cosmobl::multipole_xi2 (const int indexR, const vector<double> mu, const vector< vector<double> > xi) 
{
  double bin = mu[1]-mu[0];
  double xi2 = 0.;

  for (unsigned int j=0; j<xi[indexR].size(); j++) 
    xi2 += xi[indexR][j]*P_2(mu[j])*bin;

  return 5.*xi2;
}


// ============================================================================


double cosmobl::multipole_xi4 (const int indexR, const vector<double> mu, const vector< vector<double> > xi) 
{
  double bin = mu[1]-mu[0];
  double xi4 = 0.;

  for (unsigned int j=0; j<xi[indexR].size(); j++) 
    xi4 += xi[indexR][j]*P_4(mu[j])*bin;

  return 9.*xi4;
}


// ============================================================================


double cosmobl::error_multipole_xi0 (const int indexR, const vector<double> mu, const vector<vector<double> > error) 
{
  double bin = mu[1]-mu[0];
  double err = 0.;

  for (unsigned int j=0; j<error[indexR].size(); j++) 
    err += pow(error[indexR][j],2.);
 
  return sqrt(err)*bin;
}


// ============================================================================


double cosmobl::error_multipole_xi2 (const int indexR, const vector<double> mu, const vector<vector<double> > error) 
{
  double bin = mu[1]-mu[0];
  double err = 0.;

  for (unsigned int j=0; j<error[indexR].size(); j++) 
    err += pow(error[indexR][j]*P_2(mu[j]),2.);
  
  return 5.*sqrt(err)*bin;
}


// ============================================================================


double cosmobl::error_multipole_xi4 (const int indexR, const vector<double> mu, const vector<vector<double> > error) 
{
  double bin = mu[1]-mu[0];
  double err = 0.;

  for (unsigned int j=0; j<error[indexR].size(); j++) 
    err += pow(error[indexR][j]*P_4(mu[j]),2.); 
    
  return 9.*sqrt(err)*bin;
}


// ============================================================================


double cosmobl::multipole_xi0 (const double ss, const vector<double> rp, const vector<double> pi, const vector<vector<double> > xi, const double delta_s) 
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


double cosmobl::multipole_xi2 (const double ss, const vector<double> rp, const vector<double> pi, const vector<vector<double> > xi, const double delta_s) 
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


double cosmobl::multipole_xi4 (const double ss, const vector<double> rp, const vector<double> pi, const vector<vector<double> > xi, const double delta_s) 
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


double cosmobl::error_multipole_xi0 (const double ss, const vector<double> rp, const vector<double> pi, const vector<vector<double> > error, const double delta_s) 
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


double cosmobl::error_multipole_xi2 (const double ss, const vector<double> rp, const vector<double> pi, const vector<vector<double> > error, const double delta_s) 
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


double cosmobl::error_multipole_xi4 (const double ss, const vector<double> rp, const vector<double> pi, const vector<vector<double> > error, const double delta_s) 
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

double cosmobl::multipoles (double rr, shared_ptr<void> pp, vector<double> par)
{ 
  
  int index = par[par.size()-1];
  
 
  // -----------------------------------------------------
  // ---- compute xi(rp,pi) with the dispersion model ----
  // -----------------------------------------------------

  shared_ptr<cosmobl::glob::STR_xi2D_model> vec = static_pointer_cast<cosmobl::glob::STR_xi2D_model >(pp);
  
  vector< vector<double> > Xi (vec->dim, vector<double> (vec->dim,-1.e30));
 
  int index2 = 0;
  
  for (int i=0; i<vec->dim; i++) 
    for (int j=0; j<vec->dim; j++) {
      par[par.size()-1] = index2++;
      Xi[i][j] = xi2D_model(vec->rp[i],vec->pi[j],pp,par);
    }

  /*
  // -----------------------------
  // ---- eliminate null bins ----
  // -----------------------------
 
  vector< vector<double> > XiR = Xi;
  vector<double> rpp, pii;

  for (unsigned int i=0; i<vec->rp.size(); i++) {
    rpp.push_back(vec->rp[i]); 
    pii.push_back(vec->pi[i]); 
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
    xi_cos[0].push_back(interpolated_2D(rp, pi, vec->rp, vec->pi, Xi, "Linear"));
  }


  // --------------------------------------------------
  // ---- measure the multipole of the correlation ----
  // --------------------------------------------------
  
  if (vec->type[index]==1) return multipole_xi0(0,cos_lin,xi_cos); 
  else if (vec->type[index]==2) return multipole_xi2(0,cos_lin,xi_cos); 
  else { ErrorMsg ("Error in the function multipoles of FuncMultipoles.cpp!"); return 0; } 
  
}

/// @endcond

// ============================================================================


double cosmobl::multipole_xi0_model (const double beta, const double xi_real) 
{ 
  return xi_ratio(beta)*xi_real;
}


// ============================================================================


double cosmobl::multipole_xi0_model (const double f_sigma8, const double bias_sigma8, const double sigma8z, const double xi_DM) 
{ 
  return xi_ratio(f_sigma8,bias_sigma8)*xi_DM*pow(bias_sigma8/sigma8z,2);
}


// ============================================================================

/// @cond glob

double cosmobl::multipole_xi0_model (double xx, shared_ptr<void> pp, vector<double> par) 
{ 
  // xx is not used!

  shared_ptr<cosmobl::glob::STR_xi0_model> vec = static_pointer_cast<cosmobl::glob::STR_xi0_model>(pp);

  if (par.size()==2) return multipole_xi0_model (par[0], vec->bias_sigma8, vec->sigma8z, vec->xi_DM[par[par.size()-1]]); 

  else { ErrorMsg("Error in multipole_xi0_model of FuncMultipoles.cpp!"); return 0; }
}

/// @endcond

// ============================================================================


double cosmobl::multipole_xi2_model (const double beta, const double xi_real, const double xi_) 
{ 
  return (4./3.*beta+4./7.*beta*beta)*(xi_real-xi_);
}


// ============================================================================


double cosmobl::multipole_xi4_model (const double beta, const double xi_real, const double xi_, const double xi__) 
{ 
  return (8./35.*beta*beta)*(xi_real+2.5*xi_-3.5*xi__);
}


// Multipoles from Pk


// ============================================================================


double cosmobl::Pkl_Kaiser_integrand(const double mu, void *parameters)
{
  struct cosmobl::glob::STR_Pkl_Kaiser_integrand *pp = (struct cosmobl::glob::STR_Pkl_Kaiser_integrand *) parameters;

  return pow(pp->bias+pp->f*mu*mu,2)*gsl_sf_legendre_Pl(pp->l,mu);
} 


// ============================================================================


double cosmobl::XiMultipoles_integrand(const double kk, void *parameters)
{
  struct cosmobl::glob::STR_XiMultipoles_integrand *pp = (struct cosmobl::glob::STR_XiMultipoles_integrand *) parameters;
  double pkl = pp->Pkl->operator()(kk);
  double xx = kk*pp->r;
  double k_cut = pp->k_cut;
  double cut_pow = pp->cut_pow;
  
  return kk*kk*jl(xx,pp->l)*pkl*exp(-pow(kk/k_cut,cut_pow)); 
} 


// ============================================================================


double cosmobl::XiMultipoles_from_Xi2D_integrand(const double mu, void *parameters)
{
  struct cosmobl::glob::STR_xi2D_smu_integrand *pp = (struct cosmobl::glob::STR_xi2D_smu_integrand *) parameters;
  
  return legendre_polynomial(mu,pp->order)*pp->func->operator()(mu);
} 


// ============================================================================


double cosmobl::sigma2_integrand(const double mu, void *parameters)
{

  struct cosmobl::glob::STR_sigma2_integrand *pp = (struct cosmobl::glob::STR_sigma2_integrand *) parameters;
  int l1 = pp->l1;
  int l2 = pp->l2;
  vector<int> orders = pp->orders;
  double density_inv = pp->density_inv;
  double kk = pp->kk;
  double val = 0;
   
  for(size_t i =0;i<orders.size();i++)
    val += pp->Pk_multipoles_interp[i].operator()(kk)*gsl_sf_legendre_Pl(orders[i],mu);
  return pow(val+density_inv,2)*gsl_sf_legendre_Pl(l1,mu)*gsl_sf_legendre_Pl(l2,mu);
} 


// ============================================================================


double cosmobl::covariance_XiMultipoles_integrand (const double kk, void *parameters)
{
  struct cosmobl::glob::STR_covariance_XiMultipoles_integrand *pp = (struct cosmobl::glob::STR_covariance_XiMultipoles_integrand *) parameters;
  double s2 = pp->s2->operator()(kk); 
  double jl1k = pp->jl1r1->operator()(kk);  
  double jl2k = pp->jl2r2->operator()(kk); 
  
  return kk*kk*s2*jl1k*jl2k; 
}


// ============================================================================


double cosmobl::Pkl_Kaiser_integral(const int order, const double bias, const double f)
{
  int limit_size = 1000;
  double prec = 1.e-3;
  
  cosmobl::glob::STR_Pkl_Kaiser_integrand params;
  params.l = order;
  params.bias = bias;
  params.f = f;

  gsl_function Func;
  Func.function = &Pkl_Kaiser_integrand;
  Func.params = &params;

  return 0.5*(2*order+1)*GSL_integrate_qag(Func,-1.,1.,prec,limit_size,6);
}


// ============================================================================


vector<double> cosmobl::Pk0_Kaiser(const vector<double> kk, const vector<double> Pk, const double bias, const double f)
{
  vector<double> Pk0(kk.size(),0);
  double factor = Pkl_Kaiser_integral(0,bias,f);

  for(size_t i = 0;i<kk.size();i++)
    Pk0[i]=Pk[i]*factor;

  return Pk0;
}


// ============================================================================


vector<double> cosmobl::Pk2_Kaiser(const vector<double> kk, const vector<double> Pk, const double bias, const double f)
{
  vector<double> Pk2(kk.size(),0);
  double factor = Pkl_Kaiser_integral(2,bias,f);

  for(size_t i = 0;i<kk.size();i++)
    Pk2[i]=Pk[i]*factor;

  return Pk2;
}


// ============================================================================


vector<double> cosmobl::Pk4_Kaiser(const vector<double> kk, const vector<double> Pk, const double bias, const double f)
{
  vector<double> Pk4(kk.size(),0);
  double factor = Pkl_Kaiser_integral(4,bias,f);

  for(size_t i = 0;i<kk.size();i++)
    Pk4[i]=Pk[i]*factor;

  return Pk4;
}


// ============================================================================


vector< vector<double>> cosmobl::Pkl_Kaiser(const vector<int> orders, const vector<double> kk, const vector<double> Pk, const double bias, const double f)
{
  vector<vector<double> > Pk_multipoles(orders.size(),vector<double>(kk.size(),0));
  size_t nbin_k = kk.size();

  for(size_t ll=0;ll<orders.size();ll++)
    {
      double integral= Pkl_Kaiser_integral(orders[ll],bias,f);

      for(size_t i=0;i<nbin_k;i++)
	{
	  Pk_multipoles[ll][i] = Pk[i]*integral;;
	}

    }

  return Pk_multipoles;
}


// ============================================================================



vector<double> cosmobl::Xi0(const vector<double> r, const vector<double> kk, const vector<double> Pk0, const double k_cut, const double cut_pow, const int IntegrationMethod)
{
  double f0 = 1./(2.*par::pi*par::pi); 
  int nbins = r.size();
  int nbins_k = kk.size();

  vector<double> xi0(nbins,0);

  if (IntegrationMethod==0) //Perform trapezoid integration
    { 

      for (int i=0;i<nbins;i++){
        
        vector<double> i0(kk.size(),0);

        for(int j=0;j<nbins_k;j++)
	  i0[j] =kk[j]*kk[j]*Pk0[j]*j0(kk[j]*r[i]); 

        xi0[i] = trapezoid_integration(kk,i0)*f0;
      }

    }
  else if (IntegrationMethod==1) //Perform integration with GSL
    {
      cosmobl::glob::STR_XiMultipoles_integrand params;
      int limit_size = 1000;

      gsl_function Func;
      Func.function = &XiMultipoles_integrand;
      Func.params = &params;

      double k_min = Min(kk);
      double k_max = Max(kk); 
      double prec = 1.e-3;

      classfunc::func_grid_GSL Pk0_interp(kk,Pk0,"Spline");
      params.Pkl = &Pk0_interp;
      params.l = 0;
      params.k_cut = k_cut;
      params.cut_pow = cut_pow;

      for (int i=0;i<nbins;i++){
        params.r = r[i];
	xi0[i] = GSL_integrate_qag(Func,k_min,k_max,prec,limit_size,6)*f0;
      }

    }

  return xi0;
}


// ============================================================================


vector<double> cosmobl::Xi2 (const vector<double> rr, const vector<double> kk, const vector<double> Pk2, const double k_cut, const double cut_pow, const int IntegrationMethod)
{
  double f2 = -1./(2.*par::pi*par::pi); // f4=1./(2.*par::pi*par::pi);
  int nbins = rr.size();
  int nbins_k = kk.size();

  vector<double> xi2(nbins, 0);

  if (IntegrationMethod==0) // perform trapezoid integration
    { 
      for (int i=0; i<nbins; i++) {
        
        vector<double> i2(kk.size(),0);

        for (int j=0; j<nbins_k; j++)
	  i2[j] = kk[j]*kk[j]*Pk2[j]*j2(kk[j]*rr[i]); 

        xi2[i] = trapezoid_integration(kk, i2)*f2;
      }

    }
  
  else if (IntegrationMethod==1) // perform integration with GSL
    {
      cosmobl::glob::STR_XiMultipoles_integrand params;
      int limit_size = 1000;

      gsl_function Func;
      Func.function = &XiMultipoles_integrand;
      Func.params = &params;

      double k_min = Min(kk);
      double k_max = Max(kk); 
      double prec = 1.e-3;

      classfunc::func_grid_GSL Pk2_interp(kk,Pk2,"Spline");
      params.Pkl = &Pk2_interp;
      params.l = 2;
      params.k_cut = k_cut;
      params.cut_pow = cut_pow;

      for (int i=0; i<nbins; i++) {
        params.r = rr[i];
	xi2[i] = GSL_integrate_qag(Func, k_min, k_max, prec, limit_size, 6)*f2;
      }

    }

  return xi2;
}


// ============================================================================


vector<double> cosmobl::Xi4(const vector<double> rr, const vector<double> kk, const vector<double> Pk4, const double k_cut, const double cut_pow, const int IntegrationMethod)
{
  double f4 = 1./(2.*par::pi*par::pi); 
  int nbins = rr.size();
  int nbins_k = kk.size();

  vector<double> xi4(nbins, 0);

  if (IntegrationMethod==0) //Perform trapezoid integration
    { 

      for (int i=0; i<nbins; i++) {
        
        vector<double> i4(kk.size(), 0);

        for (int j=0; j<nbins_k; j++)
	  i4[j] = kk[j]*kk[j]*Pk4[j]*j4(kk[j]*rr[i]); 
	
        xi4[i] = trapezoid_integration(kk, i4)*f4;
      }

    }
  else if (IntegrationMethod==1) // perform integration with GSL
    {
      cosmobl::glob::STR_XiMultipoles_integrand params;
      int limit_size = 1000;

      gsl_function Func;
      Func.function = &XiMultipoles_integrand;
      Func.params = &params;

      double k_min = Min(kk);
      double k_max = Max(kk); 
      double prec = 1.e-3;

      classfunc::func_grid_GSL Pk4_interp(kk, Pk4, "Spline");
      params.Pkl = &Pk4_interp;
      params.l = 4;
      params.k_cut = k_cut;
      params.cut_pow = cut_pow;

      for (int i=0; i<nbins; i++) {
        params.r = rr[i];
	xi4[i] = GSL_integrate_qag(Func, k_min, k_max, prec, limit_size, 6)*f4;
      }

    }

  return xi4;
}


// ============================================================================


vector<vector<double>> cosmobl::Xi02_AP (const double alpha_perpendicular, const double alpha_parallel, const vector<double> rr, const vector<double> rl, const vector<double> Xi0, const vector<double> Xi2)
{
  classfunc::func_grid_GSL xi0_interp(rl, Xi0, "Spline");
  classfunc::func_grid_GSL xi2_interp(rl, Xi2, "Spline");

  vector<double> xi0_new, xi2_new;

  if ((alpha_perpendicular-1)<1.e-30 && (alpha_parallel-1.)<1.e-30)
    for (size_t i=0; i<rr.size(); i++) {
      xi0_new.push_back(xi0_interp(rr[i]));
      xi2_new.push_back(xi2_interp(rr[i]));
    }

  else {

    double nbin_mu = 50.;
    vector<double> mu = linear_bin_vector(nbin_mu, 0., 1.);
    vector<double> xismu_0(nbin_mu, 0), xismu_2(nbin_mu, 0);

    for (size_t j=0; j<rr.size(); j++) {
      for (size_t i=0; i<nbin_mu; i++) {
	double alpha = sqrt(pow(alpha_parallel*mu[i], 2)+pow(alpha_perpendicular, 2)*(1.-mu[i]*mu[i]));
	double mu_new = mu[i]*alpha_parallel/alpha;
	double s_new = alpha*rr[i];
	xismu_0[i] = xi0_interp(s_new)+xi2_interp(s_new)*legendre_polynomial(mu_new, 2);
	xismu_2[i] = xismu_0[i]*legendre_polynomial(mu[i], 2);
      }
      xi0_new.push_back(trapezoid_integration(mu, xismu_0));
      xi2_new.push_back(5*trapezoid_integration(mu, xismu_2));
    }
  }

  return {xi0_new, xi2_new};
}


// ============================================================================


vector<vector<double>> cosmobl::Xi024_AP (const double alpha_perpendicular, const double alpha_parallel, const vector<double> rr, const vector<double> rl, const vector<double> Xi0, const vector<double> Xi2, const vector<double> Xi4)
{
  classfunc::func_grid_GSL xi0_interp(rl, Xi0, "Spline");
  classfunc::func_grid_GSL xi2_interp(rl, Xi2, "Spline");
  classfunc::func_grid_GSL xi4_interp(rl, Xi4, "Spline");

  vector<double> xi0_new, xi2_new, xi4_new;

  if ((alpha_perpendicular-1.)<1.e-30 && (alpha_parallel-1.)<1.e-30)
    for (size_t i=0; i<rr.size(); i++) {
      xi0_new.push_back(xi0_interp(rr[i]));
      xi2_new.push_back(xi2_interp(rr[i]));
      xi4_new.push_back(xi4_interp(rr[i]));
    }

  else {

    double nbin_mu = 50.;
    vector<double> mu = linear_bin_vector(nbin_mu, 0., 1.);
    vector<double> xismu_0(nbin_mu, 0), xismu_2(nbin_mu, 0), xismu_4(nbin_mu, 0);

    for (size_t j=0; j<rr.size(); j++) {
      for (size_t i=0; i<nbin_mu; i++) {
	double alpha = sqrt(pow(alpha_parallel*mu[i], 2)+pow(alpha_perpendicular, 2)*(1.-mu[i]*mu[i]));
	double mu_new = mu[i]*alpha_parallel/alpha;
	double s_new = alpha*rr[j];
	xismu_0[i] = xi0_interp(s_new)+xi2_interp(s_new)*legendre_polynomial(mu_new, 2)+xi4_interp(s_new)*legendre_polynomial(mu_new, 4);
	xismu_2[i] = xismu_0[i]*legendre_polynomial(mu[i], 2);
	xismu_4[i] = xismu_0[i]*legendre_polynomial(mu[i], 4);
      }
      xi0_new.push_back(trapezoid_integration(mu, xismu_0));
      xi2_new.push_back(5*trapezoid_integration(mu, xismu_2));
      xi4_new.push_back(9*trapezoid_integration(mu, xismu_4));

    }
  }

  return {xi0_new, xi2_new, xi4_new};
}


// ============================================================================


vector<vector<double>> cosmobl::XiWedges_AP (const vector<double> mu_min, const vector<double> delta_mu, const double alpha_perpendicular, const double alpha_parallel, const vector<double> rr, const vector<double> rl, const vector<double> Xi0, const vector<double> Xi2, const vector<double> Xi4)
{
  vector<vector<double> > xi_multipoles = cosmobl::Xi024_AP(alpha_perpendicular, alpha_parallel, rr, rl, Xi0, Xi2, Xi4);
  vector<vector<double> > xi_wedges(mu_min.size(), vector<double>(rr.size(), 0));

  vector<vector<double> > legendre_integral(mu_min.size(), vector<double>(3, 0));

  for(size_t j=0; j<mu_min.size(); j++) {
    double leg_int_0 = Legendre_polynomial_mu_average(0, mu_min[j], delta_mu[j]);
    double leg_int_2 = Legendre_polynomial_mu_average(2, mu_min[j], delta_mu[j]);
    double leg_int_4 = Legendre_polynomial_mu_average(4, mu_min[j], delta_mu[j]);
    for (size_t i=0; i<rr.size(); i++)
      xi_wedges[j][i] = xi_multipoles[0][i]*leg_int_0+xi_multipoles[1][i]*leg_int_2+xi_multipoles[2][i]*leg_int_4;
  }

  return xi_wedges;
}

// ============================================================================


vector< vector<double> > cosmobl::sigma2_k (const double nObjects, const double Volume, const vector<double> kk, const vector<vector<double> > Pk_multipoles, const vector<int> orders)
{
  double prec = 1.e-3;
  size_t n_orders = orders.size();
  size_t nbin_k = kk.size();

  vector< vector<double> > sigma2(n_orders*n_orders, vector<double>(nbin_k,0));
  double density_inv = Volume/nObjects;

  int limit_size = 100;

  cosmobl::glob::STR_sigma2_integrand params;

  params.orders = orders;

  vector<classfunc::func_grid_GSL> Pk_multipoles_interp;
  for(size_t i=0;i<n_orders;i++){
    Pk_multipoles_interp.push_back(classfunc::func_grid_GSL(kk,Pk_multipoles[i],"Spline"));
  }

  params.Pk_multipoles_interp = Pk_multipoles_interp;
  params.density_inv = density_inv;

  gsl_function Func;
  Func.function = &sigma2_integrand;
  Func.params = &params;

  for (size_t i=0;i<n_orders;i++){
    params.l1 = orders[i];
    for (size_t j=0;j<n_orders;j++){
      params.l2 = orders[j];
      int index = j+n_orders*i;
      for(size_t k=0;k<kk.size();k++){
	params.kk = kk[k];
	double Int = GSL_integrate_qag(Func,-1,1.,prec,limit_size,6);
	sigma2[index][k] = (2*orders[i]+1)*(2*orders[j]+1)*Int/Volume;
      }
    }
  }

  for(size_t i=0;i<n_orders;i++)
    Pk_multipoles_interp[i].free();

  return sigma2;
}


// ============================================================================


vector<vector<double> > cosmobl::Covariance_XiMultipoles (const int nbins, const double rMin, const double rMax, const double nObjects, const double Volume, const vector<double> kk, const vector<vector<double> > Pk_multipoles, const vector<int> orders)
{
  int n_orders = orders.size();
  int nbins_k = kk.size();
  vector<double> r = linear_bin_vector(nbins, rMin, rMax);
  double dr=r[1]-r[0];

  vector<vector<double>> covariance(n_orders*nbins,vector<double>(n_orders*nbins, 0));
  vector<vector<double>> sigma2 = sigma2_k(nObjects, Volume, kk, Pk_multipoles, orders);

  vector<vector<vector<double> >> jr(n_orders,vector<vector<double>>(nbins,vector<double>(nbins_k, 0)));

  for (int l=0; l<n_orders; l++) 
    for (int i=0; i<nbins; i++) 
      for (int j=0; j<nbins_k; j++) 
	jr[l][i][j] = jl_distance_average(kk[j], orders[l], r[i]-dr*0.5, r[i]+dr*0.5);

  cosmobl::glob::STR_covariance_XiMultipoles_integrand params;
  int limit_size = 1000;

  gsl_function Func;
  Func.function = &covariance_XiMultipoles_integrand;
  Func.params = &params;

  double k_min= max(1.e-4, cosmobl::Min(kk));
  double k_max = min(Max(kk),1.); //1.e0; 
  double prec = 1.e-2;
  complex<double> ii = complex<double>(0.,1.);

  // auto-correlation terms
  for (int l=0; l<n_orders; l++) {
    int index = l+n_orders*l;

    classfunc::func_grid_GSL s2(kk, sigma2[index], "Spline");
    params.s2 = &s2;
    for (int i=0; i<nbins; i++) {
      classfunc::func_grid_GSL jl1r1(kk, jr[l][i], "Spline");
      for (int j=i; j<nbins; j++) {
	classfunc::func_grid_GSL jl2r2(kk, jr[l][j], "Spline");
	params.jl1r1 = &jl1r1;
	params.jl2r2 = &jl2r2;

	double Int = GSL_integrate_qag(Func, k_min, k_max, prec, limit_size, 6);
	Int = Int/(2.*par::pi*par::pi);
	covariance[i+nbins*l][j+nbins*l] = Int;
	covariance[j+nbins*l][i+nbins*l] = Int;
	jl2r2.free();
      }
      jl1r1.free();
    }
    s2.free();
  }

  
  // cross-correlation terms

  for (int l1=0; l1<n_orders; l1++) {
    for (int l2=l1+1; l2<n_orders; l2++) {
      int index = l2+n_orders*l1;
      int sign = pow(ii,orders[l1]+orders[l2]).real();

      classfunc::func_grid_GSL s2(kk, sigma2[index], "Spline");
      params.s2 = &s2;

      for (int i=0; i<nbins; i++) {
	classfunc::func_grid_GSL jl1r1(kk, jr[l1][i], "Spline");
	for (int j=0; j<nbins; j++) {
	  classfunc::func_grid_GSL jl2r2(kk, jr[l2][j], "Spline");
	  params.jl1r1 = &jl1r1;
	  params.jl2r2 = &jl2r2;

	  double Int = GSL_integrate_qag(Func, k_min, k_max, prec, limit_size, 6);
	  Int = sign*Int/(2.*par::pi*par::pi);
	  covariance[i+nbins*l1][j+nbins*l2] = Int;
	  covariance[j+nbins*l2][i+nbins*l1] = Int;
	  jl2r2.free();
	}
	jl1r1.free();
      }
      s2.free();
    }  
  }

  return covariance;
}


// ============================================================================


vector<vector<double> > cosmobl::Covariance_XiWedges (const vector<double> mu, const vector<double> delta_mu, const int nbins, const double rMin, const double rMax, const double nObjects, const double Volume, const vector<double> kk, const vector<vector<double> > Pk_multipoles, const vector<int> orders)
{
  int n_wedges = mu.size();
  vector<int> ord = orders;
  vector<vector<double> > Pkl = Pk_multipoles;
  
  if (n_wedges>2 && ord.size()==3) {
    vector<double> Pk6(kk.size(), 0);
    ord.push_back(6);
    Pkl.push_back(Pk6);
  }

  int n_orders = ord.size();
  vector<vector<double>> covariance_wedges(n_wedges*nbins, vector<double>(n_wedges*nbins, 0));
  vector<vector<double>> covariance_multipoles = cosmobl::Covariance_XiMultipoles(nbins, rMin, rMax, nObjects, Volume, kk, Pkl, ord);

  for (int w1=0; w1<n_wedges; w1++) {
    for (int w2=0; w2<n_wedges; w2++) {

      for (int r1 =0; r1<nbins; r1++) {
	for (int r2 =0; r2<nbins; r2++) {
	  double VV = 0;

	  for (int l1=0; l1<n_orders; l1++) {
	    double leg_integral1 = Legendre_polynomial_mu_average(ord[l1], mu[w1], delta_mu[w1]);
	    for(int l2=0; l2<n_orders; l2++) {
	      double leg_integral2 = Legendre_polynomial_mu_average(ord[l2], mu[w2], delta_mu[w2]);
	      VV += covariance_multipoles[r1+nbins*l1][r2+nbins*l2]*leg_integral1*leg_integral2;
	    }
	  }
	  covariance_wedges[r1+nbins*w1][r2+nbins*w2] = VV;
	}
      }

    }
  }

  return covariance_wedges;
}


