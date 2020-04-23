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
 *  @author federico.marulli3@unibo.it
 */

#include "Cosmology.h"

using namespace std;

using namespace cbl;
using namespace cosmology;


// =====================================================================================


double cbl::cosmology::Cosmology::denominator_Q (const double r1, const double r2, const double theta, const vector<double> rr, const vector<double> xi_DM) const
{
  const double r3 = sqrt(r1*r1+r2*r2-2.*r1*r2*cos(theta));
  const glob::FuncGrid interp_xi_DM(rr, xi_DM, "Spline");

  const double xi1 = interp_xi_DM(r1);
  const double xi2 = interp_xi_DM(r2);
  const double xi3 = interp_xi_DM(r3);

  return xi1*xi2+xi2*xi3+xi3*xi1;
}


// =====================================================================================


void cbl::cosmology::Cosmology::integrals_Q_nonLocal (vector<double> &xi_DM, vector<double> &Phi, const vector<double> rr, const vector<double> kk, const vector<double> Pk_DM, const double prec) const
{
  xi_DM = wrapper::fftlog::transform_FFTlog(rr, 1, kk, Pk_DM, 0, 0, 1, 1);

  const int nk = kk.size();
  glob::FuncGrid interp_Pk = glob::FuncGrid(kk, Pk_DM, "Spline", BinType::_logarithmic_);

  Phi.erase(Phi.begin(), Phi.end());
  Phi.resize(nk, 0);

  for (size_t i=0; i<rr.size(); i++) {

    auto integrand = [&] (const double _k)
    {
      const double kr = _k*rr[i];
      return pow(TopHat_WF(kr), 2)*interp_Pk(_k)*sin(kr)/kr; 
    };
    
    Phi[i] = 1./(2.*pow(par::pi, 2))*wrapper::gsl::GSL_integrate_qag(integrand, 0, 100, prec);
  }
  
}


// =====================================================================================


double cbl::cosmology::Cosmology::Gamma_3PCF (const double r1, const double r2, const double theta, const vector<double> xi, const vector<double> dPhi) const
{
  return ((xi[0]+3*dPhi[0]/r1)*(xi[1]+3*dPhi[1]/r2))*legendre_polynomial(cos(theta), 2);
}


// =====================================================================================


double cbl::cosmology::Cosmology::Q_nonLocal (const double r1, const double r2, const double theta, std::vector<double> &rr, std::vector<double> &xi_DM, std::vector<double> &Phi, const std::vector<double> kk, const std::vector<double> Pk_DM) const
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


std::vector<double> cbl::cosmology::Cosmology::Q_nonLocal (const double r1, const double r2, const std::vector<double> theta, const std::vector<double> kk, const std::vector<double> Pk_DM) const
{
  const int ntheta = theta.size();
  vector<double> rr, xi_DM, Phi;
  vector<double> qNL(ntheta, 0);

  for (int i=0; i<ntheta; i++)
    qNL[i] = Q_nonLocal(r1, r2, theta[i], rr, xi_DM, Phi, kk, Pk_DM);

  return qNL;
}


// =====================================================================================


void cbl::cosmology::Cosmology::integrals_zeta_Slepian (std::vector<double> &xi_DM, std::vector<double> &xi_DM_m1, std::vector<double> &xi_DM_p1, std::vector<double> &xi_DM_2, const std::vector<double> rr, const std::vector<double> kk, const std::vector<double> Pk_DM) const
{
  vector<double> Pk_DM_m1 = Pk_DM, Pk_DM_p1 = Pk_DM;
  const int nk = kk.size();

  for (int i=0; i< nk; i++) {
    Pk_DM_m1[i] *= pow(kk[i], -1);
    Pk_DM_p1[i] *= kk[i];
  }

  xi_DM = wrapper::fftlog::transform_FFTlog(rr, 1, kk, Pk_DM, 0, 0, 1, 1);
  xi_DM_m1 = wrapper::fftlog::transform_FFTlog(rr, 1, kk, Pk_DM_m1, 1, 0, 1, 1);
  xi_DM_p1 = wrapper::fftlog::transform_FFTlog(rr, 1, kk, Pk_DM_p1, 1, 0, 1, 1);
  xi_DM_2 = wrapper::fftlog::transform_FFTlog(rr, 1, kk, Pk_DM, 2, 0, 1, 1);
}


// =====================================================================================


double cbl::cosmology::Cosmology::zeta_precyclic_Slepian (const double r1, const double r2, const double mu, const double b1, const double b2, const glob::FuncGrid interp_xi_DM, const glob::FuncGrid interp_xi_DM_m1, const glob::FuncGrid interp_xi_DM_p1, const glob::FuncGrid interp_xi_DM_2) const
{
  const double mu12 = mu;

  if (abs(mu12)>1)
    return 0;
    
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


double cbl::cosmology::Cosmology::zeta_precyclic_Slepian (const double r1, const double r2, const double r3, const double deltaR, const double b1, const double b2, const glob::FuncGrid interp_xi_DM, const glob::FuncGrid interp_xi_DM_m1, const glob::FuncGrid interp_xi_DM_p1, const glob::FuncGrid interp_xi_DM_2) const
{
  double r1Min = r1-0.5*deltaR;
  double r1Max = r1+0.5*deltaR;
  double r2Min = r2-0.5*deltaR;
  double r2Max = r2+0.5*deltaR;
  double r3Min = r3-0.5*deltaR;
  double r3Max = r3+0.5*deltaR;

  // Numerator
  auto integrandNum_r1 = [&] (const double _r1) {

    auto integrandNum_r2 = [&] (const double _r2) {

      double a = max(r3Min, _r2-_r1);
      double b = min(r3Max, _r2+_r1);
	
      if (a>b) 
	return 0.;
      else {
	auto integrandNum_r3 = [&] (const double _r3) {
	  double mu = (_r1*_r1+_r2*_r2-_r3*_r3)/(2*_r1*_r2);
	  return zeta_precyclic_Slepian (_r1, _r2, mu, b1, b2, interp_xi_DM, interp_xi_DM_m1, interp_xi_DM_p1, interp_xi_DM_2);
	};

	return wrapper::gsl::GSL_integrate_cquad(integrandNum_r3, a, b, 1.e-4);
      }
    };

    return wrapper::gsl::GSL_integrate_cquad(integrandNum_r2, r2Min, r2Max, 1.e-4);
  };

  double num = wrapper::gsl::GSL_integrate_cquad(integrandNum_r1, r1Min, r1Max, 1.e-4);

  // Denominator
  auto integrandDen_r1 = [&] (const double _r1) {

    auto integrandDen_r2 = [&] (const double _r2) {

      double a = max(r3Min, _r2-_r1);
      double b = min(r3Max, _r2+_r1);

      if (a>b) 
	return 0.;
      else {
	auto integrandDen_r3 = [&] (const double _r3) {
	  (void)_r3;
	  return 1.;
	};

	return wrapper::gsl::GSL_integrate_cquad(integrandDen_r3,  a, b, 1.e-4);
      }
    };

    return wrapper::gsl::GSL_integrate_cquad(integrandDen_r2, r2Min, r2Max, 1.e-4);
  };

  double den = wrapper::gsl::GSL_integrate_cquad(integrandDen_r1, r1Min, r1Max, 1.e-4);

  return num/den;
}


// =====================================================================================


std::vector<double> cbl::cosmology::Cosmology::zeta_expansion_Slepian (const double r1, const double r2, const double b1, const double b2, std::vector<double> &rr, std::vector<double> &xi_DM, std::vector<double> &xi_DM_m1, std::vector<double> &xi_DM_p1, std::vector<double> &xi_DM_2, const int norders, const double prec) const
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
    zeta_r1_r2[i] = 0.5*(2*i+1)*wrapper::gsl::GSL_integrate_qag(integrand, -1, 1, prec);
  }

  return zeta_r1_r2;
}


// =====================================================================================


double cbl::cosmology::Cosmology::zeta_DM_Slepian (const double r1, const double r2, const double theta, std::vector<double> &rr, std::vector<double> &xi_DM, std::vector<double> &xi_DM_m1, std::vector<double> &xi_DM_p1, std::vector<double> &xi_DM_2, const std::vector<double> kk, const std::vector<double> Pk_DM, const int norders, const double prec) const
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


double cbl::cosmology::Cosmology::Q_DM_Slepian (const double r1, const double r2, const double theta, std::vector<double> &rr, std::vector<double> &xi_DM, std::vector<double> &xi_DM_m1, std::vector<double> &xi_DM_p1, std::vector<double> &xi_DM_2, const std::vector<double> kk, const std::vector<double> Pk_DM, const int norders, const double prec) const
{
  const double zeta_DM = zeta_DM_Slepian(r1, r2, theta, rr, xi_DM, xi_DM_m1, xi_DM_p1, xi_DM_2, kk, Pk_DM, norders, prec);

  return zeta_DM/denominator_Q(r1, r2, theta, rr, xi_DM);
}


// =====================================================================================


void cbl::cosmology::Cosmology::integrals_zeta_BarrigaGatzanaga (std::vector<double> &xi_DM, std::vector<double> &Phi, const std::vector<double> rr, const std::vector<double> kk, const std::vector<double> Pk_DM) const
{
  const int nk = kk.size();

  vector<double> Pk_DM_m4 = Pk_DM;
  for (int i=0; i< nk; i++)
    Pk_DM_m4[i] *= pow(kk[i], -2);
  
  xi_DM = wrapper::fftlog::transform_FFTlog(rr, 1, kk, Pk_DM, 0, 0, 1, 1);
  Phi = wrapper::fftlog::transform_FFTlog(rr, 1, kk, Pk_DM_m4, 0, 0, 1, 1);
}


// =====================================================================================


double cbl::cosmology::Cosmology::zeta_single_BarrigaGatzanaga (const double r1, const double r2, const double theta, const std::vector<double> xi, const std::vector<double> dxi, const std::vector<double> dPhi) const
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


double cbl::cosmology::Cosmology::zeta_DM_BarrigaGatzanaga (const double r1, const double r2, const double theta, std::vector<double> &rr, std::vector<double> &xi_DM, std::vector<double> &Phi, const std::vector<double> kk, const std::vector<double> Pk_DM) const
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


double cbl::cosmology::Cosmology::Q_DM_BarrigaGatzanaga (const double r1, const double r2, const double theta, std::vector<double> &rr, std::vector<double> &xi_DM, std::vector<double> &Phi, const std::vector<double> kk, const std::vector<double> Pk_DM) const
{
  return zeta_DM_BarrigaGatzanaga(r1, r2, theta, rr, xi_DM, Phi, kk, Pk_DM)/denominator_Q(r1, r2, theta, rr, xi_DM);
}


// =====================================================================================


std::vector<double> cbl::cosmology::Cosmology::zeta_DM (const double r1, const double r2, const std::vector<double> theta, const string model, const std::vector<double> kk, const std::vector<double> Pk_DM) const
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
    ErrorCBL("the chosen model is not implemented!", "zeta_DM", "3PCF.cpp");

  return zDM;
}

// =====================================================================================


std::vector<double> cbl::cosmology::Cosmology::Q_DM (const double r1, const double r2, const std::vector<double> theta, const string model, const std::vector<double> kk, const std::vector<double> Pk_DM) const
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
    ErrorCBL("the chosen model is not implemented!", "Q_DM", "3PCF.cpp");

  return qDM;

}


// =====================================================================================


std::vector<double> cbl::cosmology::Cosmology::zeta_halo (const double r1, const double r2, const std::vector<double> theta, const double b1, const double b2, const string model, const std::vector<double> kk, const std::vector<double> Pk_DM) const
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
    ErrorCBL("the chosen model is not implemented!", "z_halo", "3PCF.cpp");

  return zH;
}


// =====================================================================================


std::vector<double> cbl::cosmology::Cosmology::Q_halo (const double r1, const double r2, const std::vector<double> theta, const double b1, const double b2, const std::string model, const std::vector<double> kk, const std::vector<double> Pk_DM) const
{
  const int ntheta = theta.size();
  vector<double> qDM = Cosmology::Q_DM(r1, r2, theta, model, kk, Pk_DM);
  vector<double> qH(ntheta, 0);

  for (int i=0; i<ntheta; i++)
    qH[i] = qDM[i]/b1+b2/(b1*b1); 

  return qH;
}


// =====================================================================================


std::vector<double> cbl::cosmology::Cosmology::Q_halo (const double r1, const double r2, const std::vector<double> theta, const double b1, const double b2, const double g2, const std::string model, const std::vector<double> kk, const std::vector<double> Pk_DM) const
{
  const int ntheta = theta.size();
  vector<double> qH = Cosmology::Q_halo(r1, r2, theta, b1, b2, model, kk, Pk_DM);
  vector<double> qNL = Cosmology::Q_nonLocal(r1, r2, theta, kk, Pk_DM);

  for (int i=0; i<ntheta; i++)
    qH[i] += g2/b1*qNL[i]; 

  return qH;
}


// =====================================================================================


std::vector<double> cbl::cosmology::Cosmology::zeta_DM_eq (const std::vector<double> rr, const std::string model, const std::vector<double> kk, const std::vector<double> Pk_DM) const
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
    ErrorCBL("the chosen model is not implemented!", "zeta_DM_eq", "3PCF.cpp");

  return zDM;
}


// =====================================================================================


std::vector<double> cbl::cosmology::Cosmology::Q_DM_eq (const std::vector<double> rr, const std::string model, const std::vector<double> kk, const std::vector<double> Pk_DM) const
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
    ErrorCBL("the chosen model is not implemented!", "Q_DM_eq", "3PCF.cpp");

  return Q_DM;
}


// =====================================================================================


double cbl::cosmology::Cosmology::zeta_multipoles_covariance (const double Volume, const double nObjects, const int l, const int l_prime, const double r1, const double r2, const double r1_prime, const double r2_prime, const double deltaR, const std::vector<double> kk, const std::vector<double> Pk, const std::vector<double> rr, const std::vector<double> Xi, const double prec)
{
  (void)prec;

  const double kr0 = 1.; //check

  size_t nk = static_cast<int>(kk.size());
  size_t nr = static_cast<int>(rr.size());
  double dlogr = log(rr[1])-log(rr[0]);

  const double density = nObjects/Volume;

  // r binning

  function<double(double, int, double)> func1; 
  function<double(double, int, int, double, double)> func2; 
  function<double(double, double)> func1_noise; 

  if (deltaR<0) {
    func1 = [&] (const double kk, const int l, const double rr) {
      return jl(kk*rr, l); };
    func2 = [&] (const double kk, const int l, const int lp, const double rr, const double rrp) {
      return jl(kk*rr, l)*jl(kk*rrp, lp); };
  } else {
    func1 = [&] (const double kk, const int l, const double rr) {
      return jl_distance_average(kk, l, rr-deltaR*0.5, rr+deltaR*0.5); };
    func2 = [&] (const double kk, const int l, const int lp, const double rr, const double rrp) {
      return jl_distance_average(kk, l, rr-deltaR*0.5, rr+deltaR*0.5)*jl_distance_average(kk, lp, rrp-deltaR*0.5, rrp+deltaR*0.5); };
  }

  func1_noise = [&] (const double rr, const double bin) {
    double rmin = bin-deltaR/2;
    double rmax = bin+deltaR/2;
    if (rr>=rmin && rr<=rmax && deltaR>0)
      return 1./(4*par::pi*pow(rr, 2));
      //return 1./(4*par::pi/3*(pow(rmax, 3)-pow(rmin, 3)));
    else 
      return 0.;
  };

  // integrand 
  
  // Terms of type  I_ll = int [(Pk+1/n) jl] jl
    
  vector<double> pk_r1 = Pk, pk_r2 = Pk;
  vector<double> pk_r1p = Pk, pk_r2p = Pk;

  for (size_t i=0; i<nk; i++) {
    pk_r1[i] *= func1(kk[i], l, r1);
    pk_r2[i] *= func1(kk[i], l, r2);
    pk_r1p[i] *= func1(kk[i], l_prime, r1_prime);
    pk_r2p[i] *= func1(kk[i], l_prime, r2_prime);
  }

  vector<double> I1_r_r1 = wrapper::fftlog::transform_FFTlog(rr, 1, kk, pk_r1, l, 0, kr0, 1);
  vector<double> I1_r_r2 = wrapper::fftlog::transform_FFTlog(rr, 1, kk, pk_r2, l, 0, kr0, 1);
  vector<double> I1_r_r1p = wrapper::fftlog::transform_FFTlog(rr, 1, kk, pk_r1p, l_prime, 0, kr0, 1);
  vector<double> I1_r_r2p = wrapper::fftlog::transform_FFTlog(rr, 1, kk, pk_r2p, l_prime, 0, kr0, 1);

  for (size_t i=0; i<nr; i++) {
    I1_r_r1[i] += func1_noise(rr[i], r1)/density;
    I1_r_r2[i] += func1_noise(rr[i], r2)/density;
    I1_r_r1p[i] += func1_noise(rr[i], r1_prime)/density;
    I1_r_r2p[i] += func1_noise(rr[i], r2_prime)/density;
  }

  // l2
  vector<int> l2;
  for (int ll=abs(l-l_prime); ll<l+l_prime+1; ll++)
    l2.push_back(ll);

  const size_t n_l2 = l2.size();

  // Terms of type  I_l1_l2_l = int [(Pk+1/n) jl1 jl2] jl
  
  vector<double> pk_r1_r1p = Pk, pk_r2_r2p = Pk, pk_r2_r1p = Pk, pk_r1_r2p = Pk;
  for (size_t i=0; i<nk; i++) {
    pk_r1_r1p[i] = (pk_r1_r1p[i]+1./density)*func2(kk[i], l, l_prime, r1, r1_prime);
    pk_r2_r2p[i] = (pk_r2_r2p[i]+1./density)*func2(kk[i], l, l_prime, r2, r2_prime);
    pk_r2_r1p[i] = (pk_r2_r1p[i]+1./density)*func2(kk[i], l, l_prime, r2, r1_prime);
    pk_r1_r2p[i] = (pk_r1_r2p[i]+1./density)*func2(kk[i], l, l_prime, r1, r2_prime);
  }

  vector<vector<double>> I2_r1_r1p(n_l2, vector<double>(nr, 0.));
  vector<vector<double>> I2_r2_r2p(n_l2, vector<double>(nr, 0.));
  vector<vector<double>> I2_r2_r1p(n_l2, vector<double>(nr, 0.));
  vector<vector<double>> I2_r1_r2p(n_l2, vector<double>(nr, 0.));

  for (size_t ll=0; ll<n_l2; ll++) {
    double wig = gsl_sf_coupling_3j(2*l, 2*l_prime, 2*l2[ll], 0, 0, 0);
    if (wig!=0) {
      I2_r1_r1p[ll] = wrapper::fftlog::transform_FFTlog(rr, 1, kk, pk_r1_r1p, l2[ll], 0, kr0, 1);
      I2_r2_r2p[ll] = wrapper::fftlog::transform_FFTlog(rr, 1, kk, pk_r2_r2p, l2[ll], 0, kr0, 1);
      I2_r2_r1p[ll] = wrapper::fftlog::transform_FFTlog(rr, 1, kk, pk_r2_r1p, l2[ll], 0, kr0, 1);
      I2_r1_r2p[ll] = wrapper::fftlog::transform_FFTlog(rr, 1, kk, pk_r1_r2p, l2[ll], 0, kr0, 1);
    }
  }

  // integrand

  double Int = 0;

  for (size_t i=0; i<nr; i++) {
    double sum = 0;

    double Xi_r = Xi[i];

    double f_r_r1 = I1_r_r1[i];
    double f_r_r2 = I1_r_r2[i];
    double f_r_r1p = I1_r_r1p[i];
    double f_r_r2p = I1_r_r2p[i];

    for (size_t ll=0; ll<n_l2; ll++) {
      int ell2 = l2[ll];
      double wig = gsl_sf_coupling_3j(2*l, 2*l_prime, 2*ell2, 0, 0, 0);
      if (wig!=0) {
	double t1 = (2*ell2+1)*pow(wig,2);
	double f_r_r1_r1p = I2_r1_r1p[ll][i];
	double f_r_r2_r2p = I2_r2_r2p[ll][i];
	double f_r_r2_r1p = I2_r2_r1p[ll][i];
	double f_r_r1_r2p = I2_r1_r2p[ll][i];

	double t2 = pow(-1, l2[ll])*Xi_r*(f_r_r1_r1p*f_r_r2_r2p+f_r_r2_r1p*f_r_r1_r2p);
	double t3  = pow(-1, 0.5*(l+l_prime+ell2))*(f_r_r1*f_r_r1p*f_r_r2_r2p+f_r_r1*f_r_r2p*f_r_r2_r1p+f_r_r2*f_r_r1p*f_r_r1_r2p+f_r_r2*f_r_r2p*f_r_r1_r1p);
	sum += t1*(t2+t3);
      }

      Int += rr[i]*rr[i]*sum*rr[i];
    }
  }
  
  double fact = 4.*par::pi/Volume*(2*l+1)*(2*l_prime+1)*pow(-1, l+l_prime);
  
  return fact*Int*dlogr;
}


// =====================================================================================


std::vector<std::vector<double>> cbl::cosmology::Cosmology::zeta_covariance (const double Volume, const double nObjects, const std::vector<double> theta, const double r1, const double r2, const double deltaR, const std::vector<double> kk, const std::vector<double> Pk, const int norders, const double prec, const bool method, const int nExtractions, std::vector<double> mean, const int seed)
{
  (void)method;
  (void)nExtractions;
  (void)mean;
  (void)seed;

  vector<double> rr, Xi;
  wrapper::fftlog::transform_FFTlog(rr, Xi, 1, kk, Pk, 0);

  vector<vector<double>> zeta_l1l2_covariance(norders, vector<double>(norders, 0.));

  for (int i=0; i<norders; i++)
    for (int j=i; j<norders; j++) {
      zeta_l1l2_covariance[i][j] = zeta_multipoles_covariance(Volume, nObjects, i, j, r1, r2, r1, r2, deltaR, kk, Pk, rr, Xi, prec);
      zeta_l1l2_covariance[j][i] = zeta_l1l2_covariance[i][j];
    }

  const int ntheta = int(theta.size());

  vector<vector<double>> Pl_theta(ntheta, vector<double>(norders, 0));
  for (int i=0; i<ntheta; i++) {
    for (int j=0; j<norders; j++)
      Pl_theta[i][j] = legendre_polynomial (cos(theta[i]), j);
  }

  vector<vector<double>> zeta_covariance(ntheta, vector<double>(ntheta, 0));
  for (int i=0; i<ntheta; i++)
    for (int j=0; j<ntheta; j++)
      for (int l1=0; l1<norders; l1++)
	for (int l2=0; l2<norders; l2++)
	  zeta_covariance[i][j] += zeta_l1l2_covariance[l1][l2]*Pl_theta[i][l1]*Pl_theta[j][l2];

  return zeta_covariance;
}


// =====================================================================================a


void cbl::cosmology::Cosmology::xi_r_n (std::vector<double> &xi_n, const std::vector<double> rr, const int nn, const std::vector<double> kk, const std::vector<double> Pk)
{
  xi_n = wrapper::fftlog::transform_FFTlog (rr, 1, kk, Pk, nn, 0, par::pi, 1);
}


// =====================================================================================


void cbl::cosmology::Cosmology::xi_r_n_pm (std::vector<double> &xi_n_p, std::vector<double> &xi_n_m, const std::vector<double> rr, const int nn, const std::vector<double> kk, const std::vector<double> Pk)
{
  vector<double> pk_p(Pk.size(), 0), pk_m(Pk.size(), 0);

  for (size_t i=0; i<Pk.size(); i++){
    pk_p[i] = kk[i]*Pk[i];
    pk_m[i] = Pk[i]/kk[i];
  }

  xi_n_p = wrapper::fftlog::transform_FFTlog (rr, 1, kk, pk_p, nn, 0, par::pi, 1);
  xi_n_m = wrapper::fftlog::transform_FFTlog (rr, 1, kk, pk_m, nn, 0, par::pi, 1);
}


// =====================================================================================


void cbl::cosmology::Cosmology::eff_l_l1 (std::vector<std::vector<double>> &eff, const std::vector<double> rr, const int l, const int l1, const std::vector<double> kk, const std::vector<double> Pk)
{
  double min_rr = Min(rr);
  double max_rr = Max(rr);
  vector<double> new_r = linear_bin_vector(rr.size(), min_rr, max_rr);
  eff.resize(rr.size());
  
  for (size_t i=0; i<rr.size(); i++)
  {
    vector<double> _pk(Pk.size(), 0);

    for (size_t j=0; j<Pk.size(); j++)
      _pk[j] = kk[j]*Pk[j]*jl(kk[j]*rr[i], l);

    eff[i] = wrapper::fftlog::transform_FFTlog (new_r, 1, kk, _pk, l1, 0, par::pi, 1);
  }
  
}


// =====================================================================================


void cbl::cosmology::Cosmology::I_ELL_ell (std::vector<std::vector<double>> &II, const std::vector<double> rr, const int ll, const int LL, const std::vector<double> kk, const std::vector<double> Pk)
{
  II.resize(rr.size(), vector<double>(rr.size(), 0));
  double min_rr = Min(rr);
  double max_rr = Max(rr);
  vector<double> new_r = linear_bin_vector(rr.size(), min_rr, max_rr);

  for (int l1 = 0; l1<=LL+ll; l1++)
  {
    if ( (LL>= fabs(l1-ll)) && (LL <= l1+ll)){
      double fact = pow(-1., l1+ll)*(2.*l1+1)*(2.*ll+1)*pow(gsl_sf_coupling_3j(2*l1, 2*ll, 2*LL, 0, 0, 0),2);

      if(fact!=0){
	vector<vector<double>> eff;
	eff_l_l1 (eff, rr, ll, l1, kk, Pk);
	for (size_t r1=0; r1<rr.size(); r1++){
	  for (size_t r2=r1; r2<rr.size(); r2++){

	    glob::FuncGrid interp_r1_eff(new_r, eff[r1], "Spline");
	    glob::FuncGrid interp_r2_eff(new_r, eff[r2], "Spline");

	    auto integrand = [&] ( const double _r) {
	      return interp_r1_eff(_r)*interp_r2_eff(_r)*_r;
	    };
	    II[r1][r2] += fact*wrapper::gsl::GSL_integrate_qag(integrand, min_rr, max_rr, 1.e-3); //Check the integral limits
	    if(r1!=r2)
	      II[r2][r1] += II[r1][r2];
	    
	  }
	}
      }
    }
  }
}


// =====================================================================================


void cbl::cosmology::Cosmology::k_ell (std::vector<std::vector<double>> &KK, const std::vector<double> rr, const int ll, const std::vector<double> kk, const std::vector<double> Pk)
{

  vector<vector<double>> I1l, I3l, I5l;

  I_ELL_ell (I1l, rr, ll, 1, kk, Pk);
  I_ELL_ell (I3l, rr, ll, 3, kk, Pk);
  I_ELL_ell (I5l, rr, ll, 5, kk, Pk);

  KK.resize(rr.size(), vector<double>(rr.size(), 0));

  const double fact = 64./77175;
  for (size_t r1=0; r1<rr.size(); r1++)
    for (size_t r2=r1; r2<rr.size(); r2++) {
      KK[r1][r2] = fact*(9.*I1l[r1][r2]-14.*I3l[r1][r2]+5.*I5l[r1][r2]);

      if(r1!=r2)
	KK[r2][r1] = KK[r1][r2];
    }

}


// =====================================================================================


double cbl::cosmology::Cosmology::zeta_ell_0_factor (const double b1, const double gamma, const double beta)
{
  return pow(b1, 3)*(34./21*(1.+4.*beta/3+1154.*beta*beta/1275+936*pow(beta, 3)/2975+21*pow(beta, 4)/425)+gamma*(1+2.*beta/3+beta*beta/9));
}


// =====================================================================================


double cbl::cosmology::Cosmology::zeta_ell_0_factor_tidal (const double gamma_t, const double beta)
{
  return 16.*beta*beta*gamma_t/675;
}


// =====================================================================================


double cbl::cosmology::Cosmology::zeta_ell_1_factor (const double b1, const double beta)
{
  return -pow(b1, 3)*(1.+4.*beta/3+82*beta*beta/75+12.*pow(beta, 3)/25+3.*pow(beta, 4)/35);
}


// =====================================================================================


double cbl::cosmology::Cosmology::zeta_ell_2_factor (const double b1, const double gamma, const double beta)
{
  return pow(b1, 3)*(8./21*(1.+4.*beta/3+52*beta*beta/21+81.*pow(beta, 3)/49+12.*pow(beta, 4)/35)+32*gamma/945*beta*beta);
}


// =====================================================================================


double cbl::cosmology::Cosmology::zeta_ell_2_factor_tidal (const double gamma_t, const double beta)
{
  return 2.5*(8./15+16*beta/45+344*beta*beta/4725)*gamma_t;
}


// =====================================================================================


double cbl::cosmology::Cosmology::zeta_ell_3_factor (const double b1, const double beta)
{
  return -pow(b1, 3)*(8*beta*beta/75+16.*pow(beta, 3)/175+8.*pow(beta, 4)/315);
}


// =====================================================================================


double cbl::cosmology::Cosmology::zeta_ell_4_factor (const double b1, const double beta)
{
  return pow(b1, 3)*(-32.*beta*beta/3675+32.*pow(beta, 3)/8575+128.*pow(beta, 4)/11025);
}


// =====================================================================================


double cbl::cosmology::Cosmology::zeta_ell_4_factor_tidal (const double gamma_t, const double beta)
{
  return 32.*beta*beta*gamma_t/525;
}


// =====================================================================================


double cbl::cosmology::Cosmology::zeta_ell_k_factor (const double b1, const double beta)
{
  return pow(b1, 3)*(7.*beta*beta+3*pow(beta, 3));
}


// =====================================================================================


double cbl::cosmology::Cosmology::zeta_ell_precyclic (const double r1, const double r2, const int ell, const double b1, const double b2, const double bt, const double beta, std::vector<std::shared_ptr<cbl::glob::FuncGrid>> interp_xi_ell, const bool use_k, std::shared_ptr<cbl::glob::FuncGrid2D> interp_k_ell)
{
  const double gamma = 2.*b2/b1;
  const double gamma_t = 2.*bt/b1;

  double fact=0;

  if(ell==0){
    fact = (zeta_ell_0_factor( b1, gamma, beta)+zeta_ell_0_factor_tidal(gamma_t, beta))*(interp_xi_ell[0]->operator()(r1)*interp_xi_ell[0]->operator()(r2));
  }
  else if (ell==1){
    fact = zeta_ell_1_factor(b1, beta)*(interp_xi_ell[0]->operator()(r1)*interp_xi_ell[1]->operator()(r2)+interp_xi_ell[0]->operator()(r2)*interp_xi_ell[1]->operator()(r1));
  }
  else if (ell==2){
    fact = (zeta_ell_2_factor( b1, gamma, beta)+zeta_ell_2_factor_tidal(gamma_t, beta))*(interp_xi_ell[0]->operator()(r1)*interp_xi_ell[0]->operator()(r2));
  }
  else if (ell==3){
    fact = zeta_ell_3_factor(b1, beta)*(interp_xi_ell[0]->operator()(r1)*interp_xi_ell[1]->operator()(r2)+interp_xi_ell[0]->operator()(r2)*interp_xi_ell[1]->operator()(r1));
  }
  else if (ell==4){
    fact = (zeta_ell_4_factor(b1, beta) + zeta_ell_4_factor_tidal(gamma_t, beta))*(interp_xi_ell[0]->operator()(r1)*interp_xi_ell[0]->operator()(r2));
  }

  return ((use_k) ? fact+zeta_ell_k_factor (b1, beta)*interp_k_ell->operator()(r1, r2) : fact);
}


// =====================================================================================


std::vector<double> cbl::cosmology::Cosmology::zeta_RSD (const double r1, const double r2, const int ntheta, const double b1, const double b2, const double bt, const double beta, const std::vector<double> rr, const std::vector<double> kk, const std::vector<double> Pk, const bool include_limits, const int max_ll, const bool use_k)
{
  (void)max_ll; (void)use_k;
  vector<vector<double>> Kl (rr.size(), vector<double>(rr.size(), 0));
  auto interp_Kl= make_shared< glob::FuncGrid2D > (glob::FuncGrid2D(rr, rr, Kl, "Linear"));

  const int nbins = (include_limits) ? ntheta+1 : ntheta;
  vector<double> zeta(ntheta, 0);
  vector<double> r3(ntheta), cos_a12(ntheta), cos_a23(ntheta), cos_a31(ntheta);

  double theta_binSize = (include_limits) ? 0. : par::pi/ntheta*0.5;
  for(int i=0; i<nbins; i++) {
    double a12 = double(i)*par::pi/ntheta+theta_binSize;
    cos_a12[i] = cos(a12);

    r3[i] = sqrt(r1*r1+r2*r2-2.*r1*r2*cos_a12[i]);

    double a23, a31;
    if (r1<=r2) {
      a23 = asin(sin(a12)*r1/r3[i]);
      a31 = par::pi-a12-a23;
    }
    else {
      a31 = asin(sin(a12)*r2/r3[i]);
      a23 = par::pi-a12-a31;
    }
    cos_a23[i] = cos(a23);
    cos_a31[i] = cos(a31);
  }
  
  //Even multipoles (0, 2, 4)

  for (int i=0; i<3; i++) {
    double ll = 2*i;
    vector<double> xil; xi_r_n(xil, rr, ll, kk, Pk);
    auto interp_xil = make_shared<glob::FuncGrid> (glob::FuncGrid (rr, xil, "Spline") );

    for(int t=0; t<nbins; t++){
      zeta[t] += zeta_ell_precyclic (r1, r2, ll, b1, b2, bt, beta, {interp_xil}, use_k, interp_Kl)*legendre_polynomial(cos_a12[t], ll);
      zeta[t] += zeta_ell_precyclic (r2, r3[t], ll, b1, b2, bt, beta, {interp_xil}, use_k, interp_Kl)*legendre_polynomial(cos_a23[t], ll);
      zeta[t] += zeta_ell_precyclic (r3[t], r1, ll, b1, b2, bt, beta, {interp_xil}, use_k, interp_Kl)*legendre_polynomial(cos_a31[t], ll);
    }
  }

  // Odd Multipoles (1, 3)
  for (int i=0; i<2; i++) {
    double ll = 2*i+1;
    vector<double> xil_p, xil_m; xi_r_n_pm (xil_p, xil_m, rr, ll, kk, Pk);
    auto interp_xil_p = make_shared<glob::FuncGrid> (glob::FuncGrid (rr, xil_p, "Spline") );
    auto interp_xil_m = make_shared<glob::FuncGrid> (glob::FuncGrid (rr, xil_m, "Spline") );

    for(int t=0; t<nbins; t++){
      zeta[t] += zeta_ell_precyclic (r1, r2, ll, b1, b2, bt, beta, {interp_xil_p, interp_xil_m}, use_k, interp_Kl)*legendre_polynomial(cos_a12[t], ll);
      zeta[t] += zeta_ell_precyclic (r2, r3[t], ll, b1, b2, bt, beta, {interp_xil_p, interp_xil_m}, use_k, interp_Kl)*legendre_polynomial(cos_a23[t], ll);
      zeta[t] += zeta_ell_precyclic (r3[t], r1, ll, b1, b2, bt, beta, {interp_xil_p, interp_xil_m}, use_k, interp_Kl)*legendre_polynomial(cos_a31[t], ll);
    }
  }

  return zeta;
}


// =====================================================================================


std::vector<double> cbl::cosmology::Cosmology::zeta_RSD (const double r1, const double r2, const int ntheta, const double b1, const double b2, const double bt, const double redshift, const std::string method_Pk, const int step_r, const int step_k, const std::string output_dir, const bool store_output, const std::string output_root, const bool force_RealSpace, const bool include_limits, const int max_ll, const bool use_k)
{
  double rmax = r1+r2;

  double beta = (force_RealSpace) ? 0 : linear_growth_rate(redshift)/b1;
  vector<double> rr = linear_bin_vector(step_r, 1., rmax);
  vector<double> kk = logarithmic_bin_vector(step_k, 1.e-4, 10.);
  vector<double> _Pk = Pk_DM(kk, method_Pk, false, redshift, output_dir, store_output, output_root);  

  return zeta_RSD (r1, r2, ntheta, b1, b2, bt, beta, rr, kk, _Pk, include_limits, max_ll, use_k);
}

