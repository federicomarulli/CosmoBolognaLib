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
 *  @file Cosmology/Lib/PkXizSpace.cpp
 *
 *  @brief Methods of the class Cosmology used to model two-point
 *  statistics in redshift-space
 *
 *  This file contains the implementation of the methods of the class
 *  Cosmology used to model the two-point correlation function and
 *  power spectrum in redshift space
 *
 *  @author Federico Marulli 
 *
 *  @author federico.marulli3@unbo.it
 */

#include "Cosmology.h"
using namespace cosmobl;


// =====================================================================================

 
double cosmobl::Cosmology::xi0_Kaiser (double rad, double f_sigma8, double bias_sigma8, string method_Pk, double redshift, string output_root, bool xiType, double k_star, bool xiNL, int norm, double r_min, double r_max, double k_min, double k_max, double aa, bool GSL, double prec, string file_par)
{
  // ----- get the real-space DM xi(r) ----- 

  vector<double> rr, Xi;
  get_xi(rr, Xi, method_Pk, redshift, output_root, xiType, k_star, xiNL, norm, r_min, r_max, k_min, k_max, aa, GSL, prec, file_par);
  
  double err = -1, XiR = -1;  
  interpolation_extrapolation(rad,rr,Xi,"Linear",-1,&XiR,&err);

  return xi_ratio(f_sigma8, bias_sigma8)*XiR*pow(bias_sigma8/sigma8_Pk(method_Pk, redshift, output_root),2);
}


// =====================================================================================


double cosmobl::Cosmology::xi2D_DispersionModel (double rp, double pi, double f_sigma8, double bias_sigma8, double sigma12, string method_Pk, double redshift, int FV, bool NL, vector<double> rr, vector<double> Xi, vector<double> Xi_, vector<double> Xi__, string output_root, int index, bool bias_nl, double bA, bool xiType, double k_star, bool xiNL, double v_min, double v_max, int step_v, int norm, double r_min, double r_max, double k_min, double k_max, double aa, bool GSL, double prec, string file_par) 
{
  if (m_sigma8<0) { ErrorMsg("Error in cosmobl::Cosmology::xi2D_DispersionModel of PkXi_zSpace.cpp!"); return 0; }
  
  double bias = bias_sigma8/m_sigma8;
  double beta = f_sigma8/bias_sigma8;


  // ----- get the real-space xi(r) ----- 

  if (Xi.size()==0) {
    get_xi(rr, Xi, method_Pk, redshift, output_root, xiType, k_star, xiNL, norm, r_min, r_max, k_min, k_max, aa, GSL, prec, file_par);
    get_barred_xi(rr, Xi, Xi_, Xi__, method_Pk, redshift, xiType, k_star, xiNL, norm, r_min, r_max, k_min, k_max, aa, GSL, prec, file_par);
  }
 

  // ----- adding a scale dependent bias ----- 

  if (bias_nl==1) 
    for (unsigned int i=0; i<Xi.size(); i++) {
      Xi[i] *= b_nl(rr[i],bA);
      Xi_[i] *= b_nl(rr[i],bA);
      Xi__[i] *= b_nl(rr[i],bA);
    }
  

  // ----- compute xi(rp,pi) ----- 

  if (!NL) 
    return xi2D_lin_model(rp, pi, beta, bias, rr, Xi, Xi_, Xi__, index, 0, 0);
  
  else {
    double var = (1.+redshift)/HH(redshift);
    return xi2D_model(rp, pi, beta, bias, sigma12, rr, Xi, Xi_, Xi__, var, FV, index, 0, 0, v_min, v_max, step_v);
  }
}


// =====================================================================================


double cosmobl::Cosmology::xi_star (double &rr, double &redshift, string output_root, double k_star, double k_min, double k_max, bool GSL, double prec, string file_par) 
{
  string method_Pk1 = "EisensteinHu"; 
  string method_Pk2 = "CAMB";

  Pk_0 (method_Pk1, redshift, output_root, k_min, k_max, GSL, prec, file_par); 

  cosmobl::classfunc::func_xistar func (m_Omega_matter, m_Omega_baryon, m_Omega_neutrinos, m_massless_neutrinos, m_massive_neutrinos, m_Omega_DE, m_Omega_radiation, m_hh, m_scalar_amp, m_n_spec, m_w0, m_wa, m_fNL, m_type_NG, m_model, m_unit, rr, redshift, output_root, m_Pk0_EH, k_max, k_star);

  Midpnt<cosmobl::classfunc::func_xistar> q1(func,0.,1.e2); 
  Midinf<cosmobl::classfunc::func_xistar> q2(func,1.e2,1.e3);
  double Int = (rr<1) ? qromo(q1)+qromo(q2) : qromo(q1); // check!!!

  return 1./(2.*pow(par::pi,2))*Int; 
}


// =====================================================================================


double cosmobl::Cosmology::xisnl_gnw (double &rp, double &pi, double &f_sigma8, double &bias_sigma8, double &bA, double &redshift, vector<double> rr, vector<double> Xi, vector<double> Xi_, vector<double> Xi__, string output_root)
{
  string method_Pk = "EisensteinHu";
  
  return xi2D_DispersionModel (rp, pi, f_sigma8, bias_sigma8, -1, method_Pk, redshift, -1, 0, rr, Xi, Xi_, Xi__, output_root, 1, bA);
}


// =====================================================================================


double cosmobl::Cosmology::xis_gBAO (double &rp, double &pi, double &f_sigma8, double &bias_sigma8, double &redshift, vector<double> rr, vector<double> Xi, vector<double> Xi_, vector<double> Xi__, string output_root, double k_star, double x_min, double x_max, int step_x)
{
  if (m_sigma8<0) ErrorMsg("Error in cosmobl::Cosmology::xis_gBAO of PkXi_zSpace.cpp!");

  int FV = -1;
  bool NL = 0;

  string method_Pk = "CAMB";
  
  vector<double> xx = linear_bin_vector(step_x, x_min, x_max);
  double delta_x = xx[1]-xx[0];

  double xis = 0., sigma12 = -1;

  double f_g = f_sigma8/m_sigma8;

  for (unsigned int k=0; k<xx.size(); k++) {
    
    double pi_new = pi-xx[k];	

    xis += xi2D_DispersionModel(rp, pi_new, f_sigma8, bias_sigma8, sigma12, method_Pk, redshift, FV, NL, rr, Xi, Xi_, Xi__, output_root, 0, -1., 1, k_star)*f_star(xx[k], f_g, k_star)*delta_x;

  }

  return xis;
}


// =====================================================================================


double cosmobl::Cosmology::xi2D_CW (double &rp, double &pi, double &beta, double &bias_lin, double &bA, double &sigmav0, double &cmu, double &cs1, double &cs2, double &redshift, vector<double> rr1, vector<double> Xi1, vector<double> rr2, vector<double> Xi2, vector<double> Xi1_, vector<double> Xi1__, vector<double> Xi2_, vector<double> Xi2__, string output_root, bool BAO, bool xiType, double k_star, bool xiNL, double r_min, double r_max, double v_min, double v_max, int step_v, double k_min, double k_max, double x_min, double x_max, int step_x, double aa, bool GSL, double prec, string file_par)
{
  if (rr1.size()==0) {
    string method_Pk1 = "EisensteinHu"; 
    string method_Pk2 = "CAMB";
    get_xi(rr1, Xi1, method_Pk1, redshift, output_root, xiType, k_star, xiNL, 1, r_min, r_max, k_min, k_max, aa, GSL, prec, file_par);
    get_xi(rr2, Xi2, method_Pk2, redshift, output_root, xiType, k_star, xiNL, 0, r_min, r_max, k_min, k_max, aa, GSL, prec, file_par);
    get_barred_xi(rr1, Xi1, Xi1_, Xi1__, method_Pk1, redshift, xiType, k_star, xiNL, 0, r_min, r_max, k_min, k_max, aa, GSL, prec, file_par);
    get_barred_xi(rr2, Xi2, Xi2_, Xi2__, method_Pk2, redshift, xiType, k_star, xiNL, 0, r_min, r_max, k_min, k_max, aa, GSL, prec, file_par);
  }

  double var = (1.+redshift)/HH(redshift);

  double delta_v = (v_max-v_min)/step_v;

  double vel = v_min;

  double xi2D = 0.;


  for (int k=0; k<step_v; k++) {
    
    double pi_new = pi-vel*var;	

    double xi_tilde = xisnl_gnw(rp, pi_new, beta, bias_lin, bA, redshift, rr1, Xi1, Xi1_, Xi1__, output_root);
 
    if (BAO) xi_tilde += xis_gBAO(rp, pi_new, beta, bias_lin, redshift, rr2, Xi2, Xi2_, Xi2__, output_root, k_star, x_min, x_max, step_x);
    xi2D += xi_tilde*f_v(vel, rp, pi_new, var, sigmav0, cmu, cs1, cs2)*delta_v;

    vel += delta_v;
  }
 
  return xi2D;
}
