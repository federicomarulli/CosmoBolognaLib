/********************************************************************
 *  Copyright (C) 2010 by Federico Marulli and Cosimo Fedeli        *
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
 *  @file Cosmology/Lib/NG.cpp
 *
 *  @brief Methods of the class Cosmology for cosmologies with
 *  primordial non-Gaussianity
 *
 *  This file contains the implementation of the methods of the class
 *  Cosmology used to model the large-scale structures of the Universe
 *  in cosmologies with primordial non-Gaussianity
 *
 *  @authors Federico Marulli, Cosimo Fedeli 
 *
 *  @authors federico.marulli3@unbo.it, cosimo.fedeli@oabo.inaf.it
 */

#include "Cosmology.h"
using namespace cosmobl;


// =====================================================================================


double cosmobl::Cosmology::Am (const string method_Pk, const string output_root, const int norm, const double k_min, const double k_max, const bool GSL, const double prec, const string file_par) 
{
  double kk = 1.e-4;
  bool NL = 0;
  double redshift = 0.;
  return Pk(kk, method_Pk, NL, redshift, output_root, norm, k_min, k_max, GSL, prec, file_par) / pow(kk,m_n_spec);
}


// =====================================================================================


double cosmobl::Cosmology::potential_spectral_amplitude (const string method_Pk, const string output_root, const int norm, const double k_min, const double k_max, const bool GSL, const double prec, const string file_par)
{
  return 2.78548e-14 * gsl_pow_2(m_Omega_matter) * Am(method_Pk, output_root, norm, k_min, k_max, GSL, prec, file_par);
}


// =====================================================================================


double cosmobl::Cosmology::bispectrum (const vector<double> kk, const string method_Pk, const string output_root, const int norm, const double k_min, const double k_max, const bool GSL, const double prec, const string file_par)
{
  double bs = 0.0;

  double bm = potential_spectral_amplitude (method_Pk, output_root, norm, k_min, k_max, GSL, prec, file_par);

  switch (m_type_NG)
    {
    case 1: // Local shape
      bs = (pow(kk[0]*kk[1],m_n_spec-4.0)+pow(kk[0]*kk[2],m_n_spec-4.0)+pow(kk[1]*kk[2],m_n_spec-4.0));
      bs *= 2.0*gsl_pow_2(bm);
      break;

    case 2: // Equilateral shape
      bs = pow(kk[0],(m_n_spec-4.0)/3.0)*pow(kk[1],2.0*(m_n_spec-4.0)/3.0)*pow(kk[2],m_n_spec-4.0)+pow(kk[2],(m_n_spec-4.0)/3.0)*pow(kk[0],2.0*(m_n_spec-4.0)/3.0)*pow(kk[1],m_n_spec-4.0)+
	pow(kk[1],(m_n_spec-4.0)/3.0)*pow(kk[2],2.0*(m_n_spec-4.0)/3.0)*pow(kk[0],m_n_spec-4.0)+pow(kk[1],(m_n_spec-4.0)/3.0)*pow(kk[0],2.0*(m_n_spec-4.0)/3.0)*pow(kk[2],m_n_spec-4.0)+
	pow(kk[2],(m_n_spec-4.0)/3.0)*pow(kk[1],2.0*(m_n_spec-4.0)/3.0)*pow(kk[0],m_n_spec-4.0)+pow(kk[0],(m_n_spec-4.0)/3.0)*pow(kk[2],2.0*(m_n_spec-4.0)/3.0)*pow(kk[1],m_n_spec-4.0)-
	pow(kk[0]*kk[1],m_n_spec-4.0)-pow(kk[0]*kk[2],m_n_spec-4.0)-pow(kk[1]*kk[2],m_n_spec-4.0)-2.0*pow(kk[0]*kk[1]*kk[2],2.0*(m_n_spec-4.0)/3.0); 
      bs *= 6.0*gsl_pow_2(bm);
      break;

    case 3: // Enfolded shape
      bs = pow(kk[0]*kk[1],m_n_spec-4.0)+pow(kk[0]*kk[2],m_n_spec-4.0)+pow(kk[1]*kk[2],m_n_spec-4.0)+3.0*pow(kk[0]*kk[1]*kk[2],2.0*(m_n_spec-4.0)/3.0)-
	pow(kk[0],(m_n_spec-4.0)/3.0)*pow(kk[1],2.0*(m_n_spec-4.0)/3.0)*pow(kk[2],m_n_spec-4.0)-pow(kk[2],(m_n_spec-4.0)/3.0)*pow(kk[0],2.0*(m_n_spec-4.0)/3.0)*pow(kk[1],m_n_spec-4.0)-
	pow(kk[1],(m_n_spec-4.0)/3.0)*pow(kk[2],2.0*(m_n_spec-4.0)/3.0)*pow(kk[0],m_n_spec-4.0)-pow(kk[1],(m_n_spec-4.0)/3.0)*pow(kk[0],2.0*(m_n_spec-4.0)/3.0)*pow(kk[2],m_n_spec-4.0)-
	pow(kk[2],(m_n_spec-4.0)/3.0)*pow(kk[1],2.0*(m_n_spec-4.0)/3.0)*pow(kk[0],m_n_spec-4.0)-pow(kk[0],(m_n_spec-4.0)/3.0)*pow(kk[2],2.0*(m_n_spec-4.0)/3.0)*pow(kk[1],m_n_spec-4.0);
      bs *= 6.0*gsl_pow_2(bm);
      break;

    case 4: // Orthogonal shape
      bs = 3.0*pow(kk[0],(m_n_spec-4.0)/3.0)*pow(kk[1],2.0*(m_n_spec-4.0)/3.0)*pow(kk[2],m_n_spec-4.0)+3.0*pow(kk[2],(m_n_spec-4.0)/3.0)*pow(kk[0],2.0*(m_n_spec-4.0)/3.0)*pow(kk[1],m_n_spec-4.0)+
	3.0*pow(kk[1],(m_n_spec-4.0)/3.0)*pow(kk[2],2.0*(m_n_spec-4.0)/3.0)*pow(kk[0],m_n_spec-4.0)+3.0*pow(kk[1],(m_n_spec-4.0)/3.0)*pow(kk[0],2.0*(m_n_spec-4.0)/3.0)*pow(kk[2],m_n_spec-4.0)+
	3.0*pow(kk[2],(m_n_spec-4.0)/3.0)*pow(kk[1],2.0*(m_n_spec-4.0)/3.0)*pow(kk[0],m_n_spec-4.0)+3.0*pow(kk[0],(m_n_spec-4.0)/3.0)*pow(kk[2],2.0*(m_n_spec-4.0)/3.0)*pow(kk[1],m_n_spec-4.0)-
	3.0*pow(kk[0]*kk[1],m_n_spec-4.0)-3.0*pow(kk[1]*kk[2],m_n_spec-4.0)-3.0*pow(kk[0]*kk[2],m_n_spec-4.0)-8.0*pow(kk[0]*kk[1]*kk[2],2.0*(m_n_spec-4)/3.0);
      bs *= 6.0*gsl_pow_2(bm);
      break;

    default:
      ErrorMsg("Error in the value of cosmobl::Cosmology::type_NG (in bispectrum of NG.cpp)!");
    }

  return bs;
}


// =====================================================================================


double cosmobl::Cosmology::mrk (const double kk, const double mass, const string method_Pk, const string output_root, const int norm, const double k_min, const double k_max, const bool GSL, const double prec, const string file_par) 
{
  double xx = kk * Radius(mass, m_RhoZero);

  double AA = Am(method_Pk, output_root, norm, k_min, k_max, GSL, prec, file_par);
  
  bool NL = 0;
  double redshift = 0.;
  double TT = sqrt(Pk(kk, method_Pk, NL, redshift, output_root, norm, k_min, k_max, GSL, prec, file_par) / AA / pow(kk,m_n_spec));

  return 5.99170e6 * gsl_pow_2(kk) * TopHat_WF(xx) * TT / m_Omega_matter;
}


// =====================================================================================

/// @cond TEST_NG

double cosmobl::Cosmology::bias_kernel (double xx, void *params)
{
  cosmobl::glob::GSL_f_pars *pp = (cosmobl::glob::GSL_f_pars *)params;

  int ni = 16;
  double *XX = new double[ni];
  double *Weight = new double[ni];
  gauleg (0., 1., XX, Weight, ni);

  vector<double> km(3);
  km[0] = xx;
  km[2] = pp->kt;
  
  double mass = pp->mass;
  string method_Pk = pp->method_Pk;
  string output_root = pp->output_root;
  int norm = pp->norm;
  double k_min = pp->k_min;
  double k_max = pp->k_max;
  bool GSL = pp->GSL;
  double prec = pp->prec;
  string file_par = pp->file_par;


  const double xi1 = -1.0, xi2 = 1.0;  
  double yi = 0.0;

  for (int i=0; i<ni; i++) {
    double xi = xi1+(xi2-xi1)*XX[i];
    km[1] = sqrt(gsl_pow_2(xx)+gsl_pow_2(pp->kt)+2.0*xx*pp->kt*xi); 
    yi += bispectrum(km, method_Pk, output_root, norm, k_min, k_max, GSL, prec, file_par)*mrk(km[1], mass, method_Pk, output_root, norm, k_min, k_max, GSL, prec, file_par)*Weight[i];
  }

  yi *= (xi2-xi1)*gsl_pow_2(xx)*mrk(xx, mass, method_Pk, output_root, norm, k_min, k_max, GSL, prec, file_par);

  return yi;
}


// =====================================================================================


double cosmobl::Cosmology::frk_test (const double kk, const double mass, const string method_Pk, const string output_root, const int norm, const double k_min, const double k_max, const bool GSL, const double prec, const string file_par)
{
  cosmobl::glob::GSL_f_pars pp;
  struct cosmobl::glob::GSL_f_pars *ppp = &pp;

  ppp->kt = kk;
  ppp->mass = mass;
  ppp->method_Pk = method_Pk;
  ppp->output_root = output_root;
  ppp->norm = norm;
  ppp->k_min = k_min;
  ppp->k_max = k_max;
  ppp->GSL = GSL;
  ppp->prec = prec;
  ppp->file_par = file_par;
  ppp->pt_Cosmology = this;

  gsl_function Func;

  Func.params = ppp;
 
  Func.function = &glob::GSL_bias_kernel_wrapper;

  double ibs = -1., err = -1.;

  gsl_integration_workspace *ww = gsl_integration_workspace_alloc(1000);
  gsl_integration_qagiu (&Func, 0., 1.e-8, 1e-3, 1000, ww, &ibs, &err);exit(1);
  //gsl_integration_qag (&Func, 1.e-3, 1.e3, 0., 1.e-3, 1000, 6, ww, &ibs, &err); 

  gsl_integration_workspace_free (ww);
  
  double bm = potential_spectral_amplitude (method_Pk, output_root, norm, k_min, k_max, GSL, prec, file_par);

  double redshift = 0.;
  double var = SSM_norm (mass, method_Pk, redshift, output_root, k_max, file_par); 

  ibs /= 8.0*gsl_pow_2(par::pi)*var;
  ibs /= bm*pow(kk,m_n_spec-4.0);
 
  return ibs;
}

/// @endcond

// =====================================================================================

/// @cond glob

double cosmobl::glob::bias_kernel2 (const double xx, void *params)
{
  struct cosmobl::glob::STR_NG *pp = (struct cosmobl::glob::STR_NG *) params;
  
  Cosmology cosm (pp->Omega_matter, pp->Omega_baryon, pp->Omega_neutrinos, pp->massless_neutrinos, pp->massive_neutrinos, pp->Omega_DE, pp->Omega_radiation, pp->hh, pp->scalar_amp, pp->n_spec, pp->w0, pp->wa, pp->fNL, pp->type_NG, pp->model, pp->unit);

  int ni = 16;
  double *XX = new double[ni]; 
  double *Weight = new double[ni];
  gauleg (0., 1., XX, Weight, ni);

  vector<double> km(3);
  km[0] = xx;
  km[2] = pp->kt;
  
  double mass = pp->mass;
  string method_Pk = pp->method_Pk;
  string output_root = pp->output_root;
  int norm = pp->norm;
  double k_min = pp->k_min;
  double k_max = pp->k_max;
  bool GSL = pp->GSL;
  double prec = pp->prec;
  string file_par = pp->file_par;


  const double xi1 = -1.0, xi2 = 1.0;  
  double yi = 0.0;

  for (int i=0; i<ni; i++)
    {
      double xi = xi1+(xi2-xi1)*XX[i];
      km[1] = sqrt(gsl_pow_2(xx)+gsl_pow_2(pp->kt)+2.0*xx*pp->kt*xi); 
      yi += cosm.bispectrum(km, method_Pk, output_root, norm, k_min, k_max, GSL, prec, file_par)*cosm.mrk(km[1], mass, method_Pk, output_root, norm, k_min, k_max, GSL, prec, file_par)*Weight[i];
    }

  yi *= (xi2-xi1)*gsl_pow_2(xx)*cosm.mrk(xx, mass, method_Pk, output_root, norm, k_min, k_max, GSL, prec, file_par);

  return yi;
}

/// @endcond 


// =====================================================================================


double cosmobl::Cosmology::frk (const double kk, const double mass, const string method_Pk, const string output_root, const int norm, const double k_min, const double k_max, const bool GSL, const double prec, const string file_par)
{
  string dir_grid = par::DirCosmo+"Cosmology/Tables/grid_NG/bias_kernel/unit"+conv(m_unit,par::fINT)+"/";
  string MK = "mkdir -p "+dir_grid; if (system (MK.c_str())) {};

  string Norm = (m_sigma8>0) ? "_sigma8"+conv(m_sigma8,par::fDP3) : "_scalar_amp"+conv(m_scalar_amp,par::ee3);
  string file_grid = dir_grid+"grid"+Norm+"_h"+conv(m_hh,par::fDP3)+"_OmB"+conv(m_Omega_baryon,par::fDP3)+"_OmCDM"+conv(m_Omega_CDM,par::fDP3)+"_OmL"+conv(m_Omega_DE,par::fDP3)+"_OmN"+conv(m_Omega_neutrinos,par::fDP3)+"_typeNG"+conv(m_type_NG,par::fINT)+"_k"+conv(kk,par::fDP5)+".dat";

  int bin = 100;
  double x_min = 1.e-3;
  double x_max = 1.e3;
  vector<double> xx, yy;

  cosmobl::glob::STR_NG str;
  str.Omega_matter = m_Omega_matter; 
  str.Omega_baryon = m_Omega_baryon; 
  str.Omega_neutrinos = m_Omega_neutrinos; 
  str.massless_neutrinos = m_massless_neutrinos; 
  str.massive_neutrinos = m_massive_neutrinos; 
  str.Omega_DE = m_Omega_DE; 
  str.Omega_radiation = m_Omega_radiation; 
  str.hh = m_hh; 
  str.scalar_amp = m_scalar_amp; 
  str.n_spec = m_n_spec;
  str.w0 = m_w0; 
  str.wa = m_wa; 
  str.fNL = m_fNL;
  str.type_NG = m_type_NG;
  str.model = m_model;
  str.unit = m_unit;
  str.kt = kk;
  str.mass = mass;
  str.method_Pk = method_Pk;
  str.output_root = output_root;
  str.norm = norm;
  str.k_min = k_min;
  str.k_max = k_max;
  str.GSL = GSL;
  str.prec = prec;
  str.file_par = file_par;

  bin_function (file_grid, glob::bias_kernel2, &str, bin, x_min, x_max, "loglin", xx, yy);


  cosmobl::glob::STR_grid str_grid;
  str_grid._xx = xx;
  str_grid._yy = yy;

  gsl_function Func;
  
  Func.function = func_grid_loglin;
  Func.params = &str_grid;

  double ibs = -1., err = -1.;

  gsl_integration_workspace *ww = gsl_integration_workspace_alloc(1000);
  //gsl_integration_qagiu (&Func, 0., 1.e-8, 1e-3, 1000, ww, &ibs, &err);
  gsl_integration_qag (&Func, 1.e-3, 1.e3, 0., 1.e-3, 1000, 6, ww, &ibs, &err); 

  gsl_integration_workspace_free (ww);

  double bm = potential_spectral_amplitude (method_Pk, output_root, norm, k_min, k_max, GSL, prec, file_par);

  double redshift = 0.;
  double var = SSM_norm (mass, method_Pk, redshift, output_root, k_max, file_par); 

  ibs /= 8.0*gsl_pow_2(par::pi)*var;
  ibs /= bm*pow(kk,m_n_spec-4.0);
 
  return ibs;
}


// =====================================================================================

/// @cond glob

double cosmobl::glob::GSL_bias_kernel_wrapper (const double xx, void *params) 
{
  cosmobl::glob::GSL_f_pars *pp = (cosmobl::glob::GSL_f_pars *)params;
  return pp->pt_Cosmology->bias_kernel(xx,params);
}

/// @endcond 

// =====================================================================================


double cosmobl::Cosmology::bias_correction (const double kk, const double mass, const string method_Pk, const string output_root, const int norm, const double k_min, const double k_max, const bool GSL, const double prec, const string file_par) 
{
  return m_fNL * 0.8 * frk(kk, mass, method_Pk, output_root, norm, k_min, k_max, GSL, prec, file_par) / mrk(kk, mass, method_Pk, output_root, norm, k_min, k_max, GSL, prec, file_par);
}


// =====================================================================================

/// @cond glob

double cosmobl::glob::skewness_kernel (double *kk, size_t dim, void *params)
{
  struct cosmobl::glob::STR_NG *pp = (struct cosmobl::glob::STR_NG *) params;
  
  Cosmology cosm (pp->Omega_matter, pp->Omega_baryon, pp->Omega_neutrinos, pp->massless_neutrinos, pp->massive_neutrinos, pp->Omega_DE, pp->Omega_radiation, pp->hh, pp->scalar_amp, pp->n_spec, pp->w0, pp->wa, pp->fNL, pp->type_NG, pp->output_root, pp->unit);

  int ni = 16;
  double *XX = new double[ni];
  double *Weight = new double[ni];
  gauleg (0., 1., XX, Weight, ni);

  vector<double> km(3);
  km[0] = kk[0];
  km[2] = kk[1];

  double mass = pp->mass;
  string method_Pk = pp->method_Pk;
  string output_root = pp->output_root;
  int norm = pp->norm;
  double k_min = pp->k_min;
  double k_max = pp->k_max;
  bool GSL = pp->GSL;
  double prec = pp->prec;
  string file_par = pp->file_par;


  const double xi1 = -1.0, xi2 = 1.0;  
  double xi, yi = 0.0;

  for (int i=0; i<ni; i++) {
    xi = xi1+(xi2-xi1)*XX[i];
    km[1] = sqrt(gsl_pow_2(kk[0])+gsl_pow_2(kk[1])+2.0*kk[0]*kk[1]*xi); 
    yi += cosm.bispectrum(km, method_Pk, output_root, norm, k_min, k_max, GSL, prec, file_par)*cosm.mrk(km[1], mass, method_Pk, output_root, norm, k_min, k_max, GSL, prec, file_par)*Weight[i];
  }

  yi *= (xi2-xi1);

  return cosm.mrk(kk[0], mass, method_Pk, output_root, norm, k_min, k_max, GSL, prec, file_par)*cosm.mrk(kk[1], mass, method_Pk, output_root, norm, k_min, k_max, GSL, prec, file_par)*yi*gsl_pow_2(kk[0]*kk[1]);
}

/// @endcond

// =====================================================================================



double cosmobl::Cosmology::skewness (const double mass, const string method_Pk, const string output_root, const int norm, const double k_min, const double k_max, const bool GSL, const double prec, const string file_par) 
{
  string dir_grid = par::DirCosmo+"Cosmology/Tables/grid_NG/skewness_kernel/unit"+conv(m_unit,par::fINT)+"/";
  string MK = "mkdir -p "+dir_grid; if (system (MK.c_str())) {};

  string Norm = (m_sigma8>0) ? "_sigma8"+conv(m_sigma8,par::fDP3) : "_scalar_amp"+conv(m_scalar_amp,par::ee3);
  string file_grid = dir_grid+"grid"+Norm+"_h"+conv(m_hh,par::fDP3)+"_OmB"+conv(m_Omega_baryon,par::fDP3)+"_OmCDM"+conv(m_Omega_CDM,par::fDP3)+"_OmL"+conv(m_Omega_DE,par::fDP3)+"_OmN"+conv(m_Omega_neutrinos,par::fDP3)+"_typeNG"+conv(m_type_NG,par::fINT)+"_lgMass"+conv(log10(mass),par::fDP2)+".dat";
  
  // check !!!
  int bin = 20; 
  double x1_min = 1.e-5;
  double x1_max = 1.e2;
  double x2_min = 1.e-5;
  double x2_max = 1.e2;
  vector<double> xx1, xx2;
  vector< vector<double> > yy;
 
  cosmobl::glob::STR_NG str;
  str.Omega_matter = m_Omega_matter; 
  str.Omega_baryon = m_Omega_baryon; 
  str.Omega_neutrinos = m_Omega_neutrinos; 
  str.massless_neutrinos = m_massless_neutrinos; 
  str.massive_neutrinos = m_massive_neutrinos; 
  str.Omega_DE = m_Omega_DE; 
  str.Omega_radiation = m_Omega_radiation; 
  str.hh = m_hh; 
  str.scalar_amp = m_scalar_amp; 
  str.n_spec = m_n_spec;
  str.w0 = m_w0; 
  str.wa = m_wa; 
  str.fNL = m_fNL;
  str.type_NG = m_type_NG;
  str.model = m_model;
  str.unit = m_unit;
  str.kt = -1.; // check!!!
  str.mass = mass;
  str.method_Pk = method_Pk;
  str.output_root = output_root;
  str.norm = norm;
  str.k_min = k_min;
  str.k_max = k_max;
  str.GSL = GSL;
  str.prec = prec;
  str.file_par = file_par;

  bin_function_2D (file_grid, glob::skewness_kernel, &str, bin, x1_min, x1_max, x2_min, x2_max, "loglin", xx1, xx2, yy);
  
  cosmobl::glob::STR_grid_2D str_grid_2D;
  str_grid_2D._xx1 = xx1;
  str_grid_2D._xx2 = xx2;
  str_grid_2D._yy = yy;

  
  int dim = 2;
  double *kl = new double[dim];
  double *kh = new double[dim];
  
  for (int i=0; i<dim; i++) {
    kl[i] = 0.0;
    kh[i] = 30.0/pow(Radius(mass,m_RhoZero),0.7);
  }

  
  gsl_monte_function Func;
  Func.f = func_grid_loglin_2D;
  Func.dim = dim;  
  Func.params = &str_grid_2D;
  
  
  double ibs = -1., err = -1.;

  const size_t cl = 75000;

  // random number generator
  const gsl_rng_type *TT = gsl_rng_default;
  gsl_rng *rn = gsl_rng_alloc (TT);
  
  cout <<"Please wait, I'm computing the Monte Carlo integral..."<<endl;

  gsl_monte_vegas_state *st = gsl_monte_vegas_alloc (dim); 
  gsl_monte_vegas_integrate (&Func, kl, kh, dim, cl/100, rn, st, &ibs, &err); 
  //gsl_monte_vegas_integrate (&Func, kl, kh, dim, cl, rn, st, &ibs, &err);
  gsl_monte_vegas_free (st);

  double redshift = 0.;
  double var = SSM_norm (mass, method_Pk, redshift, output_root, k_max, file_par); 

  return ibs/(8.0*gsl_pow_4(par::pi))/gsl_pow_2(var)*m_fNL;
}


// =====================================================================================


double cosmobl::Cosmology::dskewnessdM (const double mass, const string method_Pk, const string output_root, const int norm, const double k_min, const double k_max, const bool GSL, const double prec, const string file_par) 
{
  double dlogm = 0.1;
  double mInf = 6.0, mSup = 16.0;

  if (log10(mass) < mInf) ErrorMsg("Error in cosmobl::Cosmology::dskewnessdM of NG.cpp: mass should be > 10^6 Msun/h");
  
  if (log10(mass) > mSup) ErrorMsg("Error in cosmobl::Cosmology::dskewnessdM of NG.cpp: mass should be < 10^16 Msun/h");
  

  double logm1 = log10(mass)-dlogm;
  double logm2 = log10(mass)+dlogm;

  if (logm1 < mInf) logm1 = mInf;
  if (logm2 > mSup) logm2 = mSup;

  double sk = skewness(mass, method_Pk, output_root, norm, k_min, k_max, GSL, prec, file_par);
  double logsk1 = 0.0;
  double logsk2 = 0.0;

  double M1 = pow(10.0,logm1);
  double M2 = pow(10.0,logm2);

  if (m_fNL != 0.0) {
    switch (m_type_NG) {
    default:
      logsk1 = log10(skewness(M1, method_Pk, output_root, norm, k_min, k_max, GSL, prec, file_par)/m_fNL);
      logsk2 = log10(skewness(M2, method_Pk, output_root, norm, k_min, k_max, GSL, prec, file_par)/m_fNL);
      break;
    case 4:
      logsk1 = log10(-skewness(M1, method_Pk, output_root, norm, k_min, k_max, GSL, prec, file_par)/m_fNL);
      logsk2 = log10(-skewness(M2, method_Pk, output_root, norm, k_min, k_max, GSL, prec, file_par)/m_fNL);
      break;
    }
  }

  return sk*(logsk2-logsk1)/(logm2-logm1)/mass;
}


// =====================================================================================


double cosmobl::Cosmology::MF_correction (const double mass, const double redshift, const string method_Pk, const string output_root, const string interpType, const int Num, const double stepsize, const int norm, const double k_min, const double k_max, const bool GSL, const double prec, const string file_par) 
{
  double dc = deltac(redshift)*sqrt(0.8);
  double gf = DD(redshift);
  
  double z0 = 0.;
  double SSS = SSM_norm(mass, method_Pk, z0, output_root, k_max, file_par);
  double sm = sqrt(SSS);

  double DlnSigmaDlnM = dnSM(1, mass, method_Pk, z0, output_root, interpType, Num, stepsize, k_max, file_par)*(mass/(2.*SSS));
  
  double sk = skewness (mass, method_Pk, output_root, norm, k_min, k_max, GSL, prec, file_par);

  double dsk = dskewnessdM (mass, method_Pk, output_root, norm, k_min, k_max, GSL, prec, file_par);
  
  double mfc1 = DlnSigmaDlnM*(dc/gf/sm+sk*sm/6.0*(gsl_pow_4(dc/gf/sm)-2.0*gsl_pow_2(dc/gf/sm)-1.0))/mass;
  double mfc2 = dsk*sm/6.0*(gsl_pow_2(dc/gf/sm)-1.0);
  double mfc = mfc1+mfc2;
  double psg = 1.;

  if (mfc > 0.0) mfc = 0.0;
  else psg = DlnSigmaDlnM*(dc/gf/sm)/mass;

  return mfc/psg;
}
