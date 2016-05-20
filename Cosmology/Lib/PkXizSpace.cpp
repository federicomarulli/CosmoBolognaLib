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

 
double cosmobl::Cosmology::xi0_Kaiser (const double rad, const double f_sigma8, const double bias_sigma8, const string method_Pk, const double redshift, const string output_root, const bool xiType, const double k_star, const bool xiNL, const int norm, const double r_min, const double r_max, const double k_min, const double k_max, const double aa, const bool GSL, const double prec, const string file_par)
{
  // ----- get the real-space DM xi(r) ----- 

  vector<double> rr, Xi;
  get_xi(rr, Xi, method_Pk, redshift, output_root, xiType, k_star, xiNL, norm, r_min, r_max, k_min, k_max, aa, GSL, prec, file_par);
  
  double XiR = interpolated(rad, rr, Xi, "Linear");

  return xi_ratio(f_sigma8, bias_sigma8)*XiR*pow(bias_sigma8/sigma8_Pk(method_Pk, redshift, output_root), 2);
}


// =====================================================================================


double cosmobl::Cosmology::xi2D_DispersionModel (const double rp, const double pi, const double f_sigma8, const double bias_sigma8, const double sigma12, const string method_Pk, const double redshift, const int FV, const bool NL, vector<double> rr, vector<double> &Xi, vector<double> &Xi_, vector<double> &Xi__, const string output_root, const int index, const bool bias_nl, const double bA, const bool xiType, const double k_star, const bool xiNL, const double v_min, const double v_max, const int step_v, const int norm, const double r_min, const double r_max, const double k_min, const double k_max, const double aa, const bool GSL, const double prec, const string file_par) 
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


double cosmobl::Cosmology::xi_star (const double rr, const double redshift, const string output_root, const double k_star, const double k_min, const double k_max, const bool GSL, const double prec, const string file_par) 
{
  string method_Pk1 = "EisensteinHu"; 
  string method_Pk2 = "CAMB";

  Pk_0(method_Pk1, redshift, output_root, k_min, k_max, GSL, prec, file_par); 

  cosmobl::classfunc::func_xistar func (m_Omega_matter, m_Omega_baryon, m_Omega_neutrinos, m_massless_neutrinos, m_massive_neutrinos, m_Omega_DE, m_Omega_radiation, m_hh, m_scalar_amp, m_n_spec, m_w0, m_wa, m_fNL, m_type_NG, m_model, m_unit, rr, redshift, output_root, m_Pk0_EH, k_max, k_star);

  Midpnt<cosmobl::classfunc::func_xistar> q1(func,0.,1.e2); 
  Midinf<cosmobl::classfunc::func_xistar> q2(func,1.e2,1.e3);
  double Int = (rr<1) ? qromo(q1)+qromo(q2) : qromo(q1); // check!!!

  return 1./(2.*pow(par::pi,2))*Int; 
}


// =====================================================================================


double cosmobl::Cosmology::xisnl_gnw (const double rp, const double pi, const double f_sigma8, const double bias_sigma8, const double bA, const double redshift, vector<double> rr, vector<double> Xi, vector<double> &Xi_, vector<double> &Xi__, const string output_root)
{
  string method_Pk = "EisensteinHu";
  
  return xi2D_DispersionModel (rp, pi, f_sigma8, bias_sigma8, -1, method_Pk, redshift, -1, 0, rr, Xi, Xi_, Xi__, output_root, 1, bA);
}


// =====================================================================================


double cosmobl::Cosmology::xis_gBAO (const double rp, const double pi, const double f_sigma8, const double bias_sigma8, const double redshift, vector<double> rr, vector<double> Xi, vector<double> &Xi_, vector<double> &Xi__, const string output_root, const double k_star, const double x_min, const double x_max, const int step_x)
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


double cosmobl::Cosmology::xi2D_CW (const double rp, const double pi, const double beta, const double bias_lin, const double bA, const double sigmav0, const double cmu, const double cs1, const double cs2, const double redshift, vector<double> rr1, vector<double> Xi1, vector<double> rr2, vector<double> Xi2, vector<double> &Xi1_, vector<double> &Xi1__, vector<double> &Xi2_, vector<double> &Xi2__, const string output_root, const bool BAO, const bool xiType, const double k_star, const bool xiNL, const double r_min, const double r_max, const double v_min, const double v_max, const int step_v, const double k_min, const double k_max, const double x_min, const double x_max, const int step_x, const double aa, const bool GSL, const double prec, const string file_par)
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
