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
 *  @author federico.marulli3@unibo.it
 */

#include "Cosmology.h"

using namespace std;

using namespace cbl;


// =====================================================================================

 
double cbl::cosmology::Cosmology::xi0_Kaiser (const double rad, const double f_sigma8, const double bias_sigma8, const std::string method_Pk, const double redshift, const bool store_output, const std::string output_root, const bool xiType, const double k_star, const bool NL, const int norm, const double r_min, const double r_max, const double k_min, const double k_max, const double aa, const bool GSL, const double prec, const std::string file_par)
{
  // ----- get the real-space DM xi(r) ----- 

  vector<double> rr, Xi;
  get_xi(rr, Xi, method_Pk, redshift, store_output, output_root, xiType, k_star, NL, norm, r_min, r_max, k_min, k_max, aa, GSL, prec, file_par);
  
  double XiR = interpolated(rad, rr, Xi, "Linear");

  return xi_ratio(f_sigma8, bias_sigma8)*XiR*pow(bias_sigma8/sigma8_Pk(method_Pk, redshift, store_output, output_root), 2);
}


// =====================================================================================

 
std::vector<double> cbl::cosmology::Cosmology::xi0_Kaiser (const std::vector<double> rad, const double bias, const std::string method_Pk, const double redshift, const std::string output_dir, const bool store_output, const std::string output_root, const bool NL, const int norm, const double k_min, const double k_max, const double prec, const std::string file_par)
{
  const vector<double> kk = logarithmic_bin_vector(100, k_min, k_max);
  const vector<double> Pk = this->Pk_DM(kk, method_Pk, NL, redshift, output_dir, store_output, output_root, norm, k_min, k_max, prec, file_par);

  vector<double> xi = wrapper::fftlog::transform_FFTlog(rad, 1, kk, Pk, 0);

  const double fact = bias*bias*xi_ratio(linear_growth_rate(redshift, 1.)/bias);
  
  for (size_t i=0; i<xi.size(); i++)
    xi[i] *= fact;

  return xi;
}


// =====================================================================================


double cbl::cosmology::Cosmology::xi2D_DispersionModel (const double rp, const double pi, const double f_sigma8, const double bias_sigma8, const double sigma12, const std::string method_Pk, const double redshift, const int FV, const bool NL, std::vector<double> rr, std::vector<double> &Xi, std::vector<double> &Xi_, std::vector<double> &Xi__, const bool store_output, const std::string output_root, const int index, const bool bias_nl, const double bA, const bool xiType, const double k_star, const bool xiNL, const double v_min, const double v_max, const int step_v, const int norm, const double r_min, const double r_max, const double k_min, const double k_max, const double aa, const bool GSL, const double prec, const std::string file_par) 
{
  if (m_sigma8<0) return ErrorCBL("sigma8<0!", "xi2D_DispersionModel", "PkXizSpace.cpp");
  
  double bias = bias_sigma8/m_sigma8; 
  double beta = f_sigma8/bias_sigma8;


  // ----- get the real-space xi(r) ----- 

  if (Xi.size()==0) {
    get_xi(rr, Xi, method_Pk, redshift, store_output, output_root, xiType, k_star, xiNL, norm, r_min, r_max, k_min, k_max, aa, GSL, prec, file_par);
    get_barred_xi(rr, Xi, Xi_, Xi__, method_Pk, redshift, xiType, k_star, xiNL, norm, r_min, r_max, k_min, k_max, aa, prec, file_par);
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


double cbl::cosmology::Cosmology::xi_star (const double rr, const double redshift, const bool store_output, const std::string output_root, const double k_star, const double k_min, const double k_max, const double prec, const std::string file_par) 
{
  string method_Pk1 = "EisensteinHu"; 
  string method_Pk2 = "CAMB";

  Pk_0(method_Pk1, redshift, store_output, output_root, k_min, k_max, prec, file_par); 

  classfunc::func_xistar func(m_Omega_matter, m_Omega_baryon, m_Omega_neutrinos, m_massless_neutrinos, m_massive_neutrinos, m_Omega_DE, m_Omega_radiation, m_hh, m_scalar_amp, m_scalar_pivot, m_n_spec, m_w0, m_wa, m_fNL, m_type_NG, m_tau, m_model, m_unit, rr, redshift, store_output, output_root, k_max, k_star);

  function<double(double)> ff = bind(&classfunc::func_xistar::operator(), func, std::placeholders::_1);

  double Int1 = wrapper::gsl::GSL_integrate_qag(ff, 0., 1.e2, 1.e-3);
  double Int2 = wrapper::gsl::GSL_integrate_qag(ff, 1.e2, 1.e3, 1.e-3);

  double Int = (rr<1) ? Int1+Int2 : Int1; // check!!!

  return 1./(2.*pow(par::pi, 2))*Int; 
}


// =====================================================================================


double cbl::cosmology::Cosmology::xisnl_gnw (const double rp, const double pi, const double f_sigma8, const double bias_sigma8, const double bA, const double redshift, std::vector<double> rr, std::vector<double> Xi, std::vector<double> &Xi_, std::vector<double> &Xi__, const bool store_output, const std::string output_root)
{
  string method_Pk = "EisensteinHu";
  
  return xi2D_DispersionModel(rp, pi, f_sigma8, bias_sigma8, -1, method_Pk, redshift, -1, 0, rr, Xi, Xi_, Xi__, store_output, output_root, 1, bA);
}


// =====================================================================================


double cbl::cosmology::Cosmology::xis_gBAO (const double rp, const double pi, const double f_sigma8, const double bias_sigma8, const double redshift, std::vector<double> rr, std::vector<double> Xi, std::vector<double> &Xi_, std::vector<double> &Xi__, const bool store_output, const std::string output_root, const double k_star, const double x_min, const double x_max, const int step_x)
{
  if (m_sigma8<0) ErrorCBL("sigma8<0!", "xis_gBAO", "PkXizSpace.cpp");

  int FV = -1;
  bool NL = 0;

  string method_Pk = "CAMB";
  
  vector<double> xx = linear_bin_vector(step_x, x_min, x_max);
  double delta_x = xx[1]-xx[0];

  double xis = 0., sigma12 = -1;

  double f_g = f_sigma8/m_sigma8;

  for (unsigned int k=0; k<xx.size(); k++) {
    
    double pi_new = pi-xx[k];	

    xis += xi2D_DispersionModel(rp, pi_new, f_sigma8, bias_sigma8, sigma12, method_Pk, redshift, FV, NL, rr, Xi, Xi_, Xi__, store_output, output_root, 0, -1., 1, k_star)*f_star(xx[k], f_g, k_star)*delta_x;

  }

  return xis;
}


// =====================================================================================


double cbl::cosmology::Cosmology::xi2D_CW (const double rp, const double pi, const double beta, const double bias_lin, const double bA, const double sigmav0, const double cmu, const double cs1, const double cs2, const double redshift, std::vector<double> rr1, std::vector<double> Xi1, std::vector<double> rr2, std::vector<double> Xi2, std::vector<double> &Xi1_, std::vector<double> &Xi1__, std::vector<double> &Xi2_, std::vector<double> &Xi2__, const bool store_output, const std::string output_root, const bool BAO, const bool xiType, const double k_star, const bool xiNL, const double r_min, const double r_max, const double v_min, const double v_max, const int step_v, const double k_min, const double k_max, const double x_min, const double x_max, const int step_x, const double aa, const bool GSL, const double prec, const std::string file_par)
{
  if (rr1.size()==0) {
    string method_Pk1 = "EisensteinHu"; 
    string method_Pk2 = "CAMB";
    get_xi(rr1, Xi1, method_Pk1, redshift, store_output, output_root, xiType, k_star, xiNL, 1, r_min, r_max, k_min, k_max, aa, GSL, prec, file_par);
    get_xi(rr2, Xi2, method_Pk2, redshift, store_output, output_root, xiType, k_star, xiNL, 0, r_min, r_max, k_min, k_max, aa, GSL, prec, file_par);
    get_barred_xi(rr1, Xi1, Xi1_, Xi1__, method_Pk1, redshift, xiType, k_star, xiNL, 0, r_min, r_max, k_min, k_max, aa, prec, file_par);
    get_barred_xi(rr2, Xi2, Xi2_, Xi2__, method_Pk2, redshift, xiType, k_star, xiNL, 0, r_min, r_max, k_min, k_max, aa, prec, file_par);
  }

  double var = (1.+redshift)/HH(redshift);

  double delta_v = (v_max-v_min)/step_v;

  double vel = v_min;

  double xi2D = 0.;


  for (int k=0; k<step_v; k++) {
    
    double pi_new = pi-vel*var;	

    double xi_tilde = xisnl_gnw(rp, pi_new, beta, bias_lin, bA, redshift, rr1, Xi1, Xi1_, Xi1__, store_output, output_root);
 
    if (BAO) xi_tilde += xis_gBAO(rp, pi_new, beta, bias_lin, redshift, rr2, Xi2, Xi2_, Xi2__, store_output, output_root, k_star, x_min, x_max, step_x);
    xi2D += xi_tilde*f_v(vel, rp, pi_new, var, sigmav0, cmu, cs1, cs2)*delta_v;

    vel += delta_v;
  }
 
  return xi2D;
}
