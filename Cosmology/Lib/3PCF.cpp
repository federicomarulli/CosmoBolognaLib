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
 *  @file Cosmology/Lib/3PCF.cpp
 *
 *  @brief Methods of the class Cosmology used to model two-point
 *  statistics
 *
 *  This file contains the implementation of the methods of the class
 *  Cosmology used to model the two-point correlation function and
 *  power spectrum
 *
 *  @author Federico Marulli
 *
 *  @author federico.marulli3@unbo.it
 */

#include "Cosmology.h"

using namespace cosmobl;
using namespace cosmology;


// =====================================================================================


double cosmobl::cosmology::Cosmology::denominator_Q (const double r1, const double r2, const double theta, const vector<double> rr, const vector<double> xi_DM) const
{
  const double r3 = sqrt(r1*r1+r2*r2-2.*r1*r2*cos(theta));
  const glob::FuncGrid interp_xi_DM(rr, xi_DM, "Spline");

  const double xi1 = interp_xi_DM(r1);
  const double xi2 = interp_xi_DM(r2);
  const double xi3 = interp_xi_DM(r3);

  return xi1*xi2+xi2*xi3+xi3*xi1;
}


// =====================================================================================


void cosmobl::cosmology::Cosmology::integrals_Q_nonLocal (vector<double> &xi_DM, vector<double> &Phi, const vector<double> rr, const vector<double> kk, const vector<double> Pk_DM, const double prec) const
{
  xi_DM = fftlog::transform_FFTlog(rr, 1, kk, Pk_DM, 0, 0, 1, 1);

  const int nk = kk.size();
  glob::FuncGrid interp_Pk = glob::FuncGrid(kk, Pk_DM, "Spline", binType::_logarithmic_);

  Phi.erase(Phi.begin(), Phi.end());
  Phi.resize(nk, 0);

  for (size_t i=0; i<rr.size(); i++) {

    auto integrand = [&] (const double _k)
    {
      const double kr = _k*rr[i];
      return pow(TopHat_WF(kr), 2)*interp_Pk(_k)*sin(kr)/kr; 
    };
    
    Phi[i] = 1./(2.*pow(par::pi, 2))*gsl::GSL_integrate_qag(integrand, 0, 100, prec);
  }
  
}


// =====================================================================================


double cosmobl::cosmology::Cosmology::Gamma_3PCF (const double r1, const double r2, const double theta, const vector<double> xi, const vector<double> dPhi) const
{
  return ((xi[0]+3*dPhi[0]/r1)*(xi[1]+3*dPhi[1]/r2))*legendre_polynomial(cos(theta), 2);
}


// =====================================================================================


double cosmobl::cosmology::Cosmology::Q_nonLocal (const double r1, const double r2, const double theta, vector<double> &rr, vector<double> &xi_DM, vector<double> &Phi, const vector<double> kk, const vector<double> Pk_DM) const
{
  if (rr.size()==0) {
    rr = linear_bin_vector(200, 1., 300.);
    integrals_Q_nonLocal(xi_DM, Phi, rr, kk, Pk_DM, 1.e-3);
  }

  glob::FuncGrid interp_xiDM(rr, xi_DM, "Spline");
  glob::FuncGrid interp_Phi(rr, Phi, "Spline");

  const double mu12 = cos(theta);
      
  const double r3 = sqrt(r1*r1+r2*r2-2.*r1*r2*mu12);
  double a12 = theta;

  double a23 = 0.;
  double a13 = 0.; 
  if (r1<=r2) {
    a23 = asin(sin(a12)*r1/r3);
    a13 = par::pi-a12-a23;
  }
  else {
    a13 = asin(sin(a12)*r2/r3);
    a23 = par::pi-a12-a13;
  }

  const double xi1 = interp_xiDM(r1);
  const double xi2 = interp_xiDM(r2);
  const double xi3 = interp_xiDM(r3);
  
  const double dPhi1 = interp_Phi.D1v(r1);
  const double dPhi2 = interp_Phi.D1v(r2);
  const double dPhi3 = interp_Phi.D1v(r3);
  
  double gamma = 0.;
  gamma += Gamma_3PCF(r1, r2, a12, {xi1, xi2}, {dPhi1, dPhi2});
  gamma += Gamma_3PCF(r2, r3, a23, {xi2, xi3}, {dPhi2, dPhi3});
  gamma += Gamma_3PCF(r3, r1, a13, {xi3, xi1}, {dPhi3, dPhi1});
  
  return 2./3*(gamma/denominator_Q(r1, r2, theta, rr, xi_DM)-1);
}


// =====================================================================================


vector<double> cosmobl::cosmology::Cosmology::Q_nonLocal (const double r1, const double r2, const vector<double> theta, const vector<double> kk, const vector<double> Pk_DM) const
{
  const int ntheta = theta.size();
  vector<double> rr, xi_DM, Phi;
  vector<double> qNL(ntheta, 0);

  for (int i=0; i<ntheta; i++)
    qNL[i] = Q_nonLocal(r1, r2, theta[i], rr, xi_DM, Phi, kk, Pk_DM);

  return qNL;
}


// =====================================================================================


void cosmobl::cosmology::Cosmology::integrals_zeta_Slepian (vector<double> &xi_DM, vector<double> &xi_DM_m1, vector<double> &xi_DM_p1, vector<double> &xi_DM_2, const vector<double> rr, const vector<double> kk, const vector<double> Pk_DM) const
{
  vector<double> Pk_DM_m1 = Pk_DM, Pk_DM_p1 = Pk_DM;
  const int nk = kk.size();

  for (int i=0; i< nk; i++) {
    Pk_DM_m1[i] *= pow(kk[i], -1);
    Pk_DM_p1[i] *= kk[i];
  }

  xi_DM = fftlog::transform_FFTlog(rr, 1, kk, Pk_DM, 0, 0, 1, 1);
  xi_DM_m1 = fftlog::transform_FFTlog(rr, 1, kk, Pk_DM_m1, 1, 0, 1, 1);
  xi_DM_p1 = fftlog::transform_FFTlog(rr, 1, kk, Pk_DM_p1, 1, 0, 1, 1);
  xi_DM_2 = fftlog::transform_FFTlog(rr, 1, kk, Pk_DM, 2, 0, 1, 1);
}


// =====================================================================================


double cosmobl::cosmology::Cosmology::zeta_precyclic_Slepian (const double r1, const double r2, const double mu, const double b1, const double b2, const glob::FuncGrid interp_xi_DM, const glob::FuncGrid interp_xi_DM_m1, const glob::FuncGrid interp_xi_DM_p1, const glob::FuncGrid interp_xi_DM_2) const
{
  const double mu12 = mu;
    
  const double r3 = sqrt(r1*r1+r2*r2-2.*r1*r2*mu12);
  const double a12 = acos(mu12);

  double a23 = 0.;
  double a13 = 0.;
  if (r1<=r2) {
    a23 = asin(sin(a12)*r1/r3);
    a13 = par::pi-a12-a23;
  }
  else {
    a13 = asin(sin(a12)*r2/r3);
    a23 = par::pi-a12-a13;
  }
  
  const double mu13 = cos(a13);
  const double mu23 = cos(a23);

  const double b1_3 = b1*b1*b1;
  const double zeta_pc_0_12 = (2.*b1*b1*b2+(34./21)*b1_3)*interp_xi_DM(r1)*interp_xi_DM(r2);
  const double zeta_pc_1_12 = -b1_3*(interp_xi_DM_m1(r1)*interp_xi_DM_p1(r2)+interp_xi_DM_m1(r2)*interp_xi_DM_p1(r1))*legendre_polynomial(mu12, 1);
  const double zeta_pc_2_12 = (8./21)*b1_3*interp_xi_DM_2(r1)*interp_xi_DM_2(r2)*legendre_polynomial(mu12, 2);

  const double zeta_pc_0_13 = (2.*b1*b1*b2+(34./21)*b1_3)*interp_xi_DM(r1)*interp_xi_DM(r3);
  const double zeta_pc_1_13 = -b1_3*(interp_xi_DM_m1(r1)*interp_xi_DM_p1(r3)+interp_xi_DM_m1(r3)*interp_xi_DM_p1(r1))*legendre_polynomial(mu13, 1);
  const double zeta_pc_2_13 = 8./21*b1_3*interp_xi_DM_2(r1)*interp_xi_DM_2(r3)*legendre_polynomial(mu13, 2);

  const double zeta_pc_0_23 = (2.*b1*b1*b2+(34./21)*b1_3)*interp_xi_DM(r2)*interp_xi_DM(r3);
  const double zeta_pc_1_23 = -b1_3*(interp_xi_DM_m1(r2)*interp_xi_DM_p1(r3)+interp_xi_DM_m1(r3)*interp_xi_DM_p1(r2))*legendre_polynomial(mu23, 1);
  const double zeta_pc_2_23 = (8./21)*b1_3*interp_xi_DM_2(r2)*interp_xi_DM_2(r3)*legendre_polynomial(mu23, 2);

  return zeta_pc_0_12+zeta_pc_1_12+zeta_pc_2_12+zeta_pc_0_13+zeta_pc_1_13+zeta_pc_2_13+zeta_pc_0_23+zeta_pc_1_23+zeta_pc_2_23;
}


// =====================================================================================


vector<double> cosmobl::cosmology::Cosmology::zeta_expansion_Slepian (const double r1, const double r2, const double b1, const double b2, vector<double> &rr, vector<double> &xi_DM, vector<double> &xi_DM_m1, vector<double> &xi_DM_p1, vector<double> &xi_DM_2, const int norders, const double prec) const
{
  glob::FuncGrid interp_xi_DM(rr, xi_DM, "Spline");
  glob::FuncGrid interp_xi_DM_m1(rr, xi_DM_m1, "Spline");
  glob::FuncGrid interp_xi_DM_p1(rr, xi_DM_p1, "Spline");
  glob::FuncGrid interp_xi_DM_2(rr, xi_DM_2, "Spline");

  vector<double> zeta_r1_r2(norders, 0);

  for (int i=0; i<norders; i++) {
    auto integrand = [&] (const double mu12) {
      return Cosmology::zeta_precyclic_Slepian(r1, r2, mu12, b1, b2, interp_xi_DM, interp_xi_DM_m1, interp_xi_DM_p1, interp_xi_DM_2)*legendre_polynomial (mu12, i);
    };
    zeta_r1_r2[i] = 0.5*(2*i+1)*gsl::GSL_integrate_qag(integrand, -1, 1, prec);
  }

  return zeta_r1_r2;
}


// =====================================================================================


double cosmobl::cosmology::Cosmology::zeta_DM_Slepian (const double r1, const double r2, const double theta, vector<double> &rr, vector<double> &xi_DM, vector<double> &xi_DM_m1, vector<double> &xi_DM_p1, vector<double> &xi_DM_2, const vector<double> kk, const vector<double> Pk_DM, const int norders, const double prec) const
{
  if (rr.size()==0) {
    rr = linear_bin_vector(200, 1., 300.);
    integrals_zeta_Slepian(xi_DM, xi_DM_m1, xi_DM_p1, xi_DM_2, rr, kk, Pk_DM);
  }

  const double mu = cos(theta);
  vector<double> z_r1_r2 = Cosmology::zeta_expansion_Slepian(r1, r2, 1, 0, rr, xi_DM, xi_DM_m1, xi_DM_p1, xi_DM_2, norders, prec);

  double zeta_r1_r2_theta = 0.;
  for (size_t i=0; i<z_r1_r2.size(); i++)
    zeta_r1_r2_theta += z_r1_r2[i]*legendre_polynomial(mu, i);

  return zeta_r1_r2_theta;
}


// =====================================================================================


double cosmobl::cosmology::Cosmology::Q_DM_Slepian (const double r1, const double r2, const double theta, vector<double> &rr, vector<double> &xi_DM, vector<double> &xi_DM_m1, vector<double> &xi_DM_p1, vector<double> &xi_DM_2, const vector<double> kk, const vector<double> Pk_DM, const int norders, const double prec) const
{
  const double zeta_DM = zeta_DM_Slepian(r1, r2, theta, rr, xi_DM, xi_DM_m1, xi_DM_p1, xi_DM_2, kk, Pk_DM, norders, prec);

  return zeta_DM/denominator_Q(r1, r2, theta, rr, xi_DM);
}


// =====================================================================================


void cosmobl::cosmology::Cosmology::integrals_zeta_BarrigaGatzanaga (vector<double> &xi_DM, vector<double> &Phi, const vector<double> rr, const vector<double> kk, const vector<double> Pk_DM) const
{
  const int nk = kk.size();

  vector<double> Pk_DM_m4 = Pk_DM;
  for (int i=0; i< nk; i++)
    Pk_DM_m4[i] *= pow(kk[i], -2);
  
  xi_DM = fftlog::transform_FFTlog(rr, 1, kk, Pk_DM, 0, 0, 1, 1);
  Phi = fftlog::transform_FFTlog(rr, 1, kk, Pk_DM_m4, 0, 0, 1, 1);
}


// =====================================================================================


double cosmobl::cosmology::Cosmology::zeta_single_BarrigaGatzanaga (const double r1, const double r2, const double theta, const vector<double> xi, const vector<double> dxi, const vector<double> dPhi) const
{
  const double mu = cos(theta);

  const double xi1 = xi[0];
  const double xi2 = xi[1];

  const double dxi1 = dxi[0];
  const double dxi2 = dxi[1];
 
  const double dPhi1 = dPhi[0];
  const double dPhi2 = dPhi[1];

  const double t1 = (10./7)*xi1*xi2;
  const double t2 = -3*dPhi1*dPhi2/(r1*r2)-xi1*dPhi2/r2-xi2*dPhi1/r1;
  const double t3 = mu*mu*( (xi1+3*dPhi1/r1)*(xi2+3*dPhi2/r2));
  const double t4 = -mu*(dxi1*dPhi2+dxi2*dPhi1);

  return t1+(4./7)*(t2+t3)+t4;
}


// =====================================================================================


double cosmobl::cosmology::Cosmology::zeta_DM_BarrigaGatzanaga (const double r1, const double r2, const double theta, vector<double> &rr, vector<double> &xi_DM, vector<double> &Phi, const vector<double> kk, const vector<double> Pk_DM) const
{
  if (rr.size()==0) {
    rr = linear_bin_vector(200, 1., 300.);
    integrals_zeta_BarrigaGatzanaga (xi_DM, Phi, rr, kk, Pk_DM);
  }

  glob::FuncGrid interp_xiDM(rr, xi_DM, "Spline");
  glob::FuncGrid interp_Phi(rr, Phi, "Spline");

  const double mu12 = cos(theta);
  
  const double r3 = sqrt(r1*r1+r2*r2-2.*r1*r2*mu12);
  const double a12 = theta;
  
  double a23 = 0.;
  double a13 = 0.; 
  if (r1<=r2) {
    a23 = asin(sin(a12)*r1/r3);
    a13 = par::pi-a12-a23;
  }
  else {
    a13 = asin(sin(a12)*r2/r3);
    a23 = par::pi-a12-a13;
  }

  const double xi1 = interp_xiDM(r1);
  const double xi2 = interp_xiDM(r2);
  const double xi3 = interp_xiDM(r3);
  
  const double dxi1 = interp_xiDM.D1v(r1);
  const double dxi2 = interp_xiDM.D1v(r2);
  const double dxi3 = interp_xiDM.D1v(r3);
  
  const double dPhi1 = interp_Phi.D1v(r1);
  const double dPhi2 = interp_Phi.D1v(r2);
  const double dPhi3 = interp_Phi.D1v(r3);

  const double t1 = zeta_single_BarrigaGatzanaga(r1, r2, a12, {xi1, xi2}, {dxi1, dxi2}, {dPhi1, dPhi2});
  const double t2 = zeta_single_BarrigaGatzanaga(r2, r3, a23, {xi2, xi3}, {dxi2, dxi3}, {dPhi2, dPhi3});
  const double t3 = zeta_single_BarrigaGatzanaga(r3, r1, a13, {xi3, xi1}, {dxi3, dxi1}, {dPhi3, dPhi1});

  const double zeta = t1+t2+t3;
  return zeta;
}


// =====================================================================================


double cosmobl::cosmology::Cosmology::Q_DM_BarrigaGatzanaga (const double r1, const double r2, const double theta, vector<double> &rr, vector<double> &xi_DM, vector<double> &Phi, const vector<double> kk, const vector<double> Pk_DM) const
{
  return zeta_DM_BarrigaGatzanaga(r1, r2, theta, rr, xi_DM, Phi, kk, Pk_DM)/denominator_Q(r1, r2, theta, rr, xi_DM);
}


// =====================================================================================


vector<double> cosmobl::cosmology::Cosmology::zeta_DM (const double r1, const double r2, const vector<double> theta, const string model, const vector<double> kk, const vector<double> Pk_DM) const
{
  const int ntheta = theta.size();
  vector<double> rr, xi_DM;
  vector<double> zDM(ntheta, 0.);

  if (model=="Slepian") {
    vector<double> xi_DM_m1, xi_DM_p1, xi_DM_2;
    for (int i=0; i<ntheta; i++)
      zDM[i] = zeta_DM_Slepian(r1, r2, theta[i], rr, xi_DM, xi_DM_m1, xi_DM_p1, xi_DM_2, kk, Pk_DM, 9, 1.e-3);
  }
  else if (model=="BarrigaGatzanaga") {
    vector<double> Phi;
    for (int i=0; i<ntheta; i++)
      zDM[i] = zeta_DM_BarrigaGatzanaga(r1, r2, theta[i], rr, xi_DM, Phi, kk, Pk_DM);
  }
  else
    ErrorCBL("Error in cosmobl::cosmology::Cosmology::zeta_DM() of 3PCF.cpp: the chosen model is not implemented!");

  return zDM;
}

// =====================================================================================


vector<double> cosmobl::cosmology::Cosmology::Q_DM (const double r1, const double r2, const vector<double> theta, const string model, const vector<double> kk, const vector<double> Pk_DM) const
{
  const int ntheta = theta.size();
  vector<double> rr, xi_DM;
  vector<double> qDM(ntheta, 0);

  if (model=="Slepian") {
    vector<double> xi_DM_m1, xi_DM_p1, xi_DM_2;
    for (int i=0; i<ntheta; i++)
      qDM[i] = Q_DM_Slepian(r1, r2, theta[i], rr, xi_DM, xi_DM_m1, xi_DM_p1, xi_DM_2, kk, Pk_DM, 9, 1.e-3);
  }
  else if (model=="BarrigaGatzanaga")
  {
    vector<double> Phi;
    for (int i=0; i<ntheta; i++)
      qDM[i] = Q_DM_BarrigaGatzanaga(r1, r2, theta[i], rr, xi_DM, Phi, kk, Pk_DM);
  }
  else
    ErrorCBL("Error in cosmobl::cosmology::Cosmology::Q_DM() of 3PCF.cpp: the chosen model is not implemented!");

  return qDM;

}


// =====================================================================================


vector<double> cosmobl::cosmology::Cosmology::zeta_halo (const double r1, const double r2, const vector<double> theta, const double b1, const double b2, const string model, const vector<double> kk, const vector<double> Pk_DM) const
{
  const int ntheta = theta.size();
  vector<double> rr, xi_DM;
  vector<double> zH(ntheta, 0);

  if (model=="Slepian") {
    vector<double> xi_DM_m1, xi_DM_p1, xi_DM_2;
    for (int i=0; i<ntheta; i++)
      zH[i] = b1*b1*b1*zeta_DM_Slepian(r1, r2, theta[i], rr, xi_DM, xi_DM_m1, xi_DM_p1, xi_DM_2, kk, Pk_DM, 9, 1.e-3)+b1*b1*b2*denominator_Q(r1, r2, theta[i], rr, xi_DM);
  }
  else if (model=="BarrigaGatzanaga")
  {
    vector<double> Phi;
    for (int i=0; i<ntheta; i++)
      zH[i] = b1*b1*b1*zeta_DM_BarrigaGatzanaga(r1, r2, theta[i], rr, xi_DM, Phi, kk, Pk_DM)+b1*b1*b2*denominator_Q(r1, r2, theta[i], rr, xi_DM);
  }
  else
    ErrorCBL("Error in cosmobl::cosmology::Cosmology::z_halo() of 3PCF.cpp: the chosen model is not implemented!");

  return zH;
}


// =====================================================================================


vector<double> cosmobl::cosmology::Cosmology::Q_halo (const double r1, const double r2, const vector<double> theta, const double b1, const double b2, const string model, const vector<double> kk, const vector<double> Pk_DM) const
{
  const int ntheta = theta.size();
  vector<double> qDM = Cosmology::Q_DM(r1, r2, theta, model, kk, Pk_DM);
  vector<double> qH(ntheta, 0);

  for (int i=0; i<ntheta; i++)
    qH[i] = qDM[i]/b1+b2/(b1*b1); 

  return qH;
}


// =====================================================================================


vector<double> cosmobl::cosmology::Cosmology::Q_halo (const double r1, const double r2, const vector<double> theta, const double b1, const double b2, const double g2, const string model, const vector<double> kk, const vector<double> Pk_DM) const
{
  const int ntheta = theta.size();
  vector<double> qH = Cosmology::Q_halo(r1, r2, theta, b1, b2, model, kk, Pk_DM);
  vector<double> qNL = Cosmology::Q_nonLocal(r1, r2, theta, kk, Pk_DM);

  for (int i=0; i<ntheta; i++)
    qH[i] += g2/b1*qNL[i]; 

  return qH;
}


// =====================================================================================


vector<double> cosmobl::cosmology::Cosmology::zeta_DM_eq (const vector<double> rr, const string model, const vector<double> kk, const vector<double> Pk_DM) const
{
  const int nr = rr.size();
  const double theta = par::pi/3;
  vector<double> _rr, xi_DM;
  vector<double> zDM(nr, 0);

  if (model=="Slepian") {
    vector<double> xi_DM_m1, xi_DM_p1, xi_DM_2;
    for (int i=0; i<nr; i++)
      zDM[i] = zeta_DM_Slepian(rr[i], rr[i], theta, _rr, xi_DM, xi_DM_m1, xi_DM_p1, xi_DM_2, kk, Pk_DM, 9, 1.e-3);
  }
  else if (model=="BarrigaGatzanaga")
  {
    vector<double> Phi;
    for (int i=0; i<nr; i++)
      zDM[i] = zeta_DM_BarrigaGatzanaga(rr[i], rr[i], theta, _rr, xi_DM, Phi, kk, Pk_DM);
  }
  else
    ErrorCBL("Error in cosmobl::cosmology::Cosmology::zeta_DM_eq() of 3PCF.cpp: the chosen model is not implemented!");

  return zDM;
}


// =====================================================================================


vector<double> cosmobl::cosmology::Cosmology::Q_DM_eq (const vector<double> rr, const string model, const vector<double> kk, const vector<double> Pk_DM) const
{
  const int nr = rr.size();
  const double theta = par::pi/3;
  vector<double> _rr, xi_DM;
  vector<double> Q_DM(nr, 0);

  if (model=="Slepian") {
    vector<double> xi_DM_m1, xi_DM_p1, xi_DM_2;
    for (int i=0; i<nr; i++)
      Q_DM[i] = Q_DM_Slepian(rr[i], rr[i], theta, _rr, xi_DM, xi_DM_m1, xi_DM_p1, xi_DM_2, kk, Pk_DM, 9, 1.e-3);
  }
  else if (model=="BarrigaGatzanaga")
  {
    vector<double> Phi;
    for (int i=0; i<nr; i++)
      Q_DM[i] = Q_DM_BarrigaGatzanaga(rr[i], rr[i], theta, _rr, xi_DM, Phi, kk, Pk_DM);
  }
  else
    ErrorCBL("Error in cosmobl::cosmology::Cosmology::Q_DM_eq() of 3PCF.cpp: the chosen model is not implemented!");

  return Q_DM;
}


// =====================================================================================


double cosmobl::cosmology::Cosmology::zeta_multipoles_covariance (const double Volume, const double nObjects, const int l, const int l_prime, const double r1, const double r2, const double r1_prime, const double r2_prime, const vector<double> kk, const vector<double> Pk, const vector<double> rr, const vector<double> Xi, const double prec)
{
  const double inverse_density = Volume/nObjects;
  const string interpType = "Spline";

  vector<double> kk_sn, pk_sn;
  for (size_t i=0; i<kk.size(); i++)
    if (kk[i]<=1) {
      kk_sn.push_back(kk[i]);
      pk_sn.push_back(Pk[i]+inverse_density);
    }

  const int nk = kk_sn.size();
  glob::FuncGrid interp_Pk(kk_sn, pk_sn, interpType);
  glob::FuncGrid interp_Xi(rr, Xi, interpType);

  // l2
  vector<int> l2;
  for (int ll=abs(l-l_prime); ll<l+l_prime+1; ll++)
    l2.push_back(ll);

  const int n_l2 = (int)l2.size();

  // integrand
  
  vector<double> pk_r1(nk, 0), pk_r2(nk, 0), pk_r1p(nk, 0), pk_r2p(nk, 0);

  for (int i=0; i<nk; i++) {
    pk_r1[i] = pk_sn[i]*jl(kk_sn[i]*r1, l);
    pk_r2[i] = pk_sn[i]*jl(kk_sn[i]*r2, l);
    pk_r1p[i] = pk_sn[i]*jl(kk_sn[i]*r1_prime, l_prime);
    pk_r2p[i] = pk_sn[i]*jl(kk_sn[i]*r2_prime, l_prime);
  }

  const double kr0 = par::pi; //check

  vector<double> f_r_r1 = fftlog::transform_FFTlog(rr, 1, kk_sn, pk_r1, l, 0, kr0, 1);
  vector<double> f_r_r2 = fftlog::transform_FFTlog(rr, 1, kk_sn, pk_r2, l, 0, kr0, 1);
  vector<double> f_r_r1p = fftlog::transform_FFTlog(rr, 1, kk_sn, pk_r1p, l_prime, 0, kr0, 1);
  vector<double> f_r_r2p = fftlog::transform_FFTlog(rr, 1, kk_sn, pk_r2p, l_prime, 0, kr0, 1);

  glob::FuncGrid interp_f_r_r1(rr, f_r_r1, interpType);
  glob::FuncGrid interp_f_r_r2(rr, f_r_r2, interpType);
  glob::FuncGrid interp_f_r_r1p(rr, f_r_r1p, interpType);
  glob::FuncGrid interp_f_r_r2p(rr, f_r_r2p, interpType);
  

  vector<double> pk_r1_r1p(nk, 0), pk_r2_r2p(nk, 0), pk_r2_r1p(nk, 0), pk_r1_r2p(nk, 0);
  for (int i=0; i<nk; i++) {
    pk_r1_r1p[i] = pk_sn[i]*jl(kk_sn[i]*r1, l)*jl(kk_sn[i]*r1_prime, l_prime);
    pk_r2_r2p[i] = pk_sn[i]*jl(kk_sn[i]*r2, l)*jl(kk_sn[i]*r2_prime, l_prime);
    pk_r2_r1p[i] = pk_sn[i]*jl(kk_sn[i]*r2, l)*jl(kk_sn[i]*r1_prime, l_prime);
    pk_r1_r2p[i] = pk_sn[i]*jl(kk_sn[i]*r1, l)*jl(kk_sn[i]*r2_prime, l_prime);
  }

  vector<shared_ptr<glob::FuncGrid>> interp_f_r_r1_r1p(n_l2), interp_f_r_r2_r2p(n_l2), interp_f_r_r1_r2p(n_l2), interp_f_r_r2_r1p(n_l2);

  for (int ll=0; ll<n_l2; ll++) {

    vector<double> f_r_r1_r1p = fftlog::transform_FFTlog(rr, 1, kk_sn, pk_r1_r1p, l2[ll], 0, kr0, 1);
    vector<double> f_r_r2_r2p = fftlog::transform_FFTlog(rr, 1, kk_sn, pk_r2_r2p, l2[ll], 0, kr0, 1);
    vector<double> f_r_r2_r1p = fftlog::transform_FFTlog(rr, 1, kk_sn, pk_r2_r1p, l2[ll], 0, kr0, 1);
    vector<double> f_r_r1_r2p = fftlog::transform_FFTlog(rr, 1, kk_sn, pk_r1_r2p, l2[ll], 0, kr0, 1);
    interp_f_r_r1_r1p[ll] = make_shared<glob::FuncGrid> (glob::FuncGrid(rr, f_r_r1_r1p, interpType));
    interp_f_r_r2_r2p[ll] = make_shared<glob::FuncGrid> (glob::FuncGrid(rr, f_r_r2_r2p, interpType));
    interp_f_r_r2_r1p[ll] = make_shared<glob::FuncGrid> (glob::FuncGrid(rr, f_r_r2_r1p, interpType));
    interp_f_r_r1_r2p[ll] = make_shared<glob::FuncGrid> (glob::FuncGrid(rr, f_r_r1_r2p, interpType));
  }

  // integrand

  auto integrand = [&] (const double rr)
  {
    double Xi_r = interp_Xi(rr);

    double f_r_r1 = interp_f_r_r1(rr);
    double f_r_r2 = interp_f_r_r2 (rr);
    double f_r_r1p = interp_f_r_r1p(rr);
    double f_r_r2p = interp_f_r_r2p(rr);
    
    double sum = 0;

    for (int ll=0; ll<n_l2; ll++) {
      double t1 = (2*l2[ll]+1)*pow(gsl_sf_coupling_3j(2*l, 2*l_prime, 2*l2[ll], 0, 0, 0),2);
      if (t1>0) {
	double f_r_r1_r1p = interp_f_r_r1_r1p[ll]->operator()(rr);
	double f_r_r2_r2p = interp_f_r_r2_r2p[ll]->operator()(rr);
	double f_r_r2_r1p = interp_f_r_r2_r1p[ll]->operator()(rr);
	double f_r_r1_r2p = interp_f_r_r1_r2p[ll]->operator()(rr);

	double t2 = pow(-1, l2[ll])*Xi_r*(f_r_r1_r1p*f_r_r2_r2p+f_r_r2_r1p*f_r_r1_r2p);
	double t3  = pow(-1, 0.5*(l+l_prime+l2[ll]))*(f_r_r1*f_r_r1p*f_r_r2_r2p+f_r_r1*f_r_r2p*f_r_r2_r1p+f_r_r2*f_r_r1p*f_r_r1_r2p+f_r_r2*f_r_r2p*f_r_r1_r1p);
	sum += t1*(t2+t3);
      }
    }

    return rr*rr*sum;
  };
  
  double fact = 4.*par::pi/Volume*(2*l+1)*(2*l_prime+1)*pow(-1, l+l_prime);
  
  return fact*gsl::GSL_integrate_qag(integrand, 0., 1000, prec);
}


// =====================================================================================


vector<vector<double>> cosmobl::cosmology::Cosmology::zeta_covariance (const double Volume, const double nObjects, const vector<double> theta, const double r1, const double r2, const vector<double> kk, const vector<double> Pk, const int norders, const double prec, const bool method, const int nExtractions, vector<double> mean, const int seed)
{
  vector<double> rr = linear_bin_vector(4100, 1.e-5, 1.e3);
  vector<double> Xi = fftlog::transform_FFTlog(rr, 1, kk, Pk, 0);

  vector<vector<double>> zeta_l1l2_covariance(norders, vector<double>(norders, 0.));

  for (int i=0; i<norders; i++)
    for (int j=i; j<norders; j++) {
      zeta_l1l2_covariance[i][j] = zeta_multipoles_covariance(Volume, nObjects, i, j, r1, r2, r1, r2, kk, Pk, rr, Xi, prec);
      zeta_l1l2_covariance[j][i] = zeta_l1l2_covariance[i][j];
    }

  const int ntheta = int(theta.size());

  vector<vector<double>> Pl_theta(ntheta, vector<double>(norders, 0));
  for (int i=0; i<ntheta; i++)
    for (int j=0; j<norders; j++)
      Pl_theta[i][j] = legendre_polynomial (cos(theta[i]), j);

  if (!method) {

    coutCBL << "method 1" << endl;
    vector<vector<double>> zeta_covariance(ntheta, vector<double>(ntheta, 0));
    for (int i=0; i<ntheta; i++)
      for (int j=0; j<ntheta; j++)
	for (int l1=0; l1<norders; l1++)
	  for (int l2=0; l2<norders; l2++)
	    zeta_covariance[i][j] += zeta_l1l2_covariance[l1][l2]*Pl_theta[i][l1]*Pl_theta[j][l2];

    return zeta_covariance;
  }
  
  else
  {
    coutCBL << "method 2" << endl;
    vector<vector<double>> zeta_covariance;
    vector<double> mean_signal(norders, 0.);
    mean_signal = (mean.size()==0) ? mean_signal : mean; 

    vector<vector<double>> zeta_l_extracted = generate_correlated_data (nExtractions, mean_signal, zeta_l1l2_covariance, seed);
    
    vector<vector<double>> zeta_extracted(nExtractions, vector<double>(ntheta, 0));
    for (int i=0; i<nExtractions; i++)
      for (int j=0; j<ntheta; j++)
	for (int l=0; l<norders; l++)
	  zeta_extracted[i][j] += zeta_l_extracted[i][l]*Pl_theta[j][l];

    covariance_matrix (zeta_extracted, zeta_covariance); 
    coutCBL <<"Done!"<<endl;

    return zeta_covariance;
  }

  vector<vector<double>> vv; return vv;
}

