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
 *  @file CosmoBolognaLib/Func/SphericalHarmonics_Coefficients.cpp
 *
 *  @brief Methods of the class
 *  SphericalHarmonics_Coefficients used to compute the spherical
 *  harmonics coefficients at a given position on
 *
 *  This file contains the implementation of the methods of the class
 *  SphericalHarmonics_Coefficients used to measure the
 *  spherical harmonics coefficients at a given position on the unit sphere
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unibo.it, alfonso.veropalumbo@unibo.it
 */

#include "SphericalHarmonics_Coefficients.h"

using namespace std;

using namespace cbl;

// ============================================================================================

shared_ptr<cbl::glob::SphericalHarmonicsArray> cbl::glob::SphericalHarmonicsArray::factory(const SpHarMethod method, const int lMax)
{
  // get correct pointer to selected method
  unique_ptr<SphericalHarmonicsArray> ptr (nullptr);

  switch (method) {
    case (cbl::glob::SpHarMethod::_STANDARD_):
      ptr.reset(new cbl::glob::SphericalHarmonicsArray_Standard(lMax));
      break;
    case (cbl::glob::SpHarMethod::_GSL_):
      ptr.reset(new cbl::glob::SphericalHarmonicsArray_GSL(lMax));
      break;
    case (cbl::glob::SpHarMethod::_EIGEN_):
      ptr.reset(new cbl::glob::SphericalHarmonicsArray_Eigen(lMax));
      break;
    default:
      ErrorCBL("UNKNOWN or WRONG spherical harmonics method. Exiting.", "cbl::glob::SphericalHarmonicsArray::factory", "SphericalHarmonics_Coefficients.cpp");
  }

  return move(ptr);
}


// ============================================================================


void cbl::glob::SphericalHarmonicsArray_Standard::m_setAlm(const int lMax) 
{
  // write triplet bins info
  m_nOrders = lMax+1;

  // compute the number of spherical harmonics
  m_nSpH = (m_nOrders)*(m_nOrders+1)/2;

  coutCBL << "Size of Alm " << m_nSpH << endl;

  m_alm.resize(m_nSpH);

  for (size_t ell=0; ell<m_nOrders; ell++) {
    for (size_t emm=0; emm<ell+1; emm++) {
      int n=ell*(ell+1)/2+emm;

      m_alm[n].reset();
      double norm =  pow(-1, emm)*sqrt(
	  (2. * ell + 1.) / (4. * cbl::par::pi) *
	  gsl_sf_fact (ell - emm) / gsl_sf_fact ( ell + emm) );
      //double norm =  pow(-1, emm)*sqrt(
      // 	       (2. * ell + 1.) / (4. * Elements::Units::pi) *
      //	       factorial<double> (ell - emm) / factorial<double> ( ell + emm) );

      m_alm[n].Norm = norm;
      m_alm[n].ell = ell;
      m_alm[n].emm = emm;
    }
  }
}


// ============================================================================


void cbl::glob::SphericalHarmonicsArray_Standard::print ()
{
  for (auto &alm : m_alm) 
    alm.print();
}


// ============================================================================


void cbl::glob::SphericalHarmonicsArray_Standard::m_resetAlm ()
{
  for (auto &alm : m_alm) 
    alm.reset();
}


// ============================================================================


void cbl::glob::SphericalHarmonicsArray_Standard::compute (const double &xx, const double &yy, const double &zz, const double weight)
{
  double cosTheta = zz;
  double sinTheta = sqrt(1 - cosTheta*cosTheta);

  double cosPhi = xx/sinTheta;
  double sinPhi = yy/sinTheta;

  double cc_real = cosPhi;
  double cc_imag = sinPhi;

  m_alm[0].set(1., weight);
  m_alm[1].set(cosTheta, weight);

  double Plm0 = 1.;
  double Plm1 = cosTheta;
  double Plm2, Plm3;

  // P(l, 0)
  for (size_t ell = 1; ell < m_nOrders-1; ell++) {
    int q = ell+1;
    Plm2 = ((2*ell+1)*cosTheta*Plm1-ell*Plm0)/(ell+1);
    m_alm[q*(q+1)/2].set(Plm2, weight);
    Plm0 = Plm1;
    Plm1 = Plm2;
  }

  // P(m, m) -> P(m+1, m) -> P( m+1<l<lMax, m)
  Plm0 = 1.;
  int ell;
  for (int m=1; m<static_cast<int>(m_nOrders-1); m++) {
    ell = m;
    int lp1 = m+1;

    Plm1 = -(2*(m-1)+1)*Plm0*sinTheta;
    Plm0 = Plm1;
    m_alm[ell*(ell+1)/2+m].set(Plm1, weight*cc_real, weight*cc_imag);
    Plm2 = cosTheta*(2*ell+1)*Plm1;

    m_alm[lp1*(lp1+1)/2+m].set(Plm2, weight*cc_real, weight*cc_imag);

    for (int lp2=lp1+1; lp2<static_cast<int>(m_nOrders); lp2++) {
      Plm3 = ((2.*lp2-1.)*cosTheta*Plm2-(lp2+m-1)*Plm1)/(lp2-m);
      Plm1 = Plm2;
      Plm2 = Plm3;

      m_alm[lp2*(lp2+1)/2+m].set(Plm3, weight*cc_real, weight*cc_imag);
    }

    double old_cc_real = cc_real;
    double old_cc_imag = cc_imag;
    cc_real = old_cc_real*cosPhi-old_cc_imag*sinPhi;
    cc_imag = old_cc_real*sinPhi+old_cc_imag*cosPhi;
  }

  ell = m_nOrders-1;
  Plm1 = -(2*(ell-1)+1)*Plm0*sinTheta;
  m_alm[ell*(ell+1)/2+ell].set(Plm1, weight*cc_real, weight*cc_imag);
}


// ============================================================================


void cbl::glob::SphericalHarmonicsArray_GSL::compute (const double &xx, const double &yy, const double &zz, const double weight)
{
  const int n_sph = gsl_sf_legendre_array_n(m_nOrders-1);
  vector<double> Plm(n_sph);

  double phi = atan2(yy, xx);
  vector<complex<double>> pow_exp(m_nOrders+1, complex<double>(cos(phi), sin(phi)));

  for (size_t mm=0; mm<pow_exp.size(); mm++)
    pow_exp[mm] = pow(pow_exp[mm], mm);

  gsl_sf_legendre_array_e(GSL_SF_LEGENDRE_NONE, m_nOrders-1, zz, 1., Plm.data());
  for(size_t ell=0; ell<m_nOrders; ell++)
    for (size_t emm=0; emm<ell+1; emm++){
      size_t n=ell*(ell+1)/2+emm;
      m_alm[n].set(Plm[n], weight*pow_exp[emm].real(), weight*pow_exp[emm].imag());
    }
}


// ============================================================================


void cbl::glob::SphericalHarmonicsArray_Eigen::m_setAlm(const int lMax) 
{
  m_zero = {0, 0, 0, 0};
  m_one = {1, 1, 1, 1};
  m_xx = {0, 0, 1, 0};
  m_yy = {0, 1, 0, 0};
  m_zz = {1, 0, 0., 1};
  m_weight = {1, 1, 1, 1};

  // write triplet bins info
  m_nOrders = lMax+1;

  // compute the number of spherical harmonics
  m_nSpH = (m_nOrders)*(m_nOrders+1)/2;

  coutCBL << "Size of Alm " << m_nSpH << endl;

  for (size_t ell=0; ell<m_nOrders; ell++) {
    for (size_t emm=0; emm<ell+1; emm++) {
      //int n=ell*(ell+1)/2+emm;
      AlmEigen alm;
      alm.reset();

      double norm =  pow(-1, emm)*sqrt(
	  (2. * ell + 1.) / (4. * cbl::par::pi) *
	  gsl_sf_fact (ell - emm) / gsl_sf_fact ( ell + emm) );
      //double norm =  pow(-1, emm)*sqrt(
      // 	       (2. * ell + 1.) / (4. * Elements::Units::pi) *
      //	       factorial<double> (ell - emm) / factorial<double> ( ell + emm) );

      alm.Norm = {norm, norm, norm, norm};
      alm.ell = ell;
      alm.emm = emm;
      m_alm.push_back(alm);
    }
  }
}
    
// ============================================================================


void cbl::glob::SphericalHarmonicsArray_Eigen::m_buffer_reset ()
{
  m_buffer_size = 0;

  m_xx = m_zero;
  m_yy = m_zero;
  m_zz = m_zero;
  m_weight = m_zero;
}


// ============================================================================


void cbl::glob::SphericalHarmonicsArray_Eigen::m_resetAlm ()
{
  for (auto &alm : m_alm) 
    alm.reset();
}


// ============================================================================


void cbl::glob::SphericalHarmonicsArray_Eigen::print ()
{
  for (auto &alm : m_alm) 
    alm.print();
}


// ============================================================================

void cbl::glob::SphericalHarmonicsArray_Eigen::compute (const double &xx, const double &yy, const double &zz, const double weight)
{
  m_xx[m_buffer_size] = xx;
  m_yy[m_buffer_size] = yy;
  m_zz[m_buffer_size] = zz;
  m_weight[m_buffer_size] = weight;
  m_buffer_size += 1;

  if (m_buffer_size > 3) {
    compute();
    m_buffer_reset();
  }
}

// ============================================================================


void cbl::glob::SphericalHarmonicsArray_Eigen::compute ()
{
  (void)m_weight;

  dVector cosTheta = m_zz;
  dVector sinTheta = (m_one - cosTheta*cosTheta).sqrt();

  dVector cosPhi =(sinTheta>0).select(m_xx/sinTheta, m_zero);
  dVector sinPhi =(sinTheta>0).select(m_yy/sinTheta, m_zero);

  dVector cc_real = cosPhi;
  dVector cc_imag = sinPhi;

  m_alm[0].set(m_one, m_weight, m_zero);
  m_alm[1].set(cosTheta, m_weight, m_zero);
  dVector Plm0 = m_one;
  dVector Plm1 = cosTheta;
  dVector Plm2, Plm3;

  // P(l, 0)
  for (int ell = 1; ell < static_cast<int>(m_nOrders)-1; ell++) {
    int q = ell+1;
    Plm2 = ((2*ell+1)*cosTheta*Plm1-ell*Plm0)/(ell+1);
    m_alm[q*(q+1)/2].set(Plm2, m_weight, m_zero);
    Plm0 = Plm1;
    Plm1 = Plm2;
  }

  // P(m, m) -> P(m+1, m) -> P( m+1<l<lMax, m)
  Plm0 = m_one;
  int ell;
  for (int m=1; m<static_cast<int>(m_nOrders-1); m++) {
    ell = m;
    int lp1 = m+1;

    Plm1 = -(2*(m-1)+1)*Plm0*sinTheta;
    Plm0 = Plm1;
    m_alm[ell*(ell+1)/2+m].set(Plm1, m_weight*cc_real, m_weight*cc_imag);
    Plm2 = cosTheta*(2*ell+1)*Plm1;

    m_alm[lp1*(lp1+1)/2+m].set(Plm2, m_weight*cc_real, m_weight*cc_imag);

    for (int lp2=lp1+1; lp2<static_cast<int>(m_nOrders); lp2++) {
      Plm3 = ((2.*lp2-1.)*cosTheta*Plm2-(lp2+m-1)*Plm1)/(lp2-m);
      Plm1 = Plm2;
      Plm2 = Plm3;

      m_alm[lp2*(lp2+1)/2+m].set(Plm3, m_weight*cc_real, m_weight*cc_imag);
    }

    dVector old_cc_real = cc_real;
    dVector old_cc_imag = cc_imag;
    cc_real = old_cc_real*cosPhi-old_cc_imag*sinPhi;
    cc_imag = old_cc_real*sinPhi+old_cc_imag*cosPhi;
  }

  ell = m_nOrders-1;
  Plm1 = -(2*(ell-1)+1)*Plm0*sinTheta;
  m_alm[ell*(ell+1)/2+ell].set(Plm1, m_weight*cc_real, m_weight*cc_imag);
}


// ============================================================================


std::vector<double> cbl::glob:: legendre_polynomials (const SphericalHarmonicsArray& Ylm_1, const SphericalHarmonicsArray& Ylm_2)
{
  (void)Ylm_2;
  std::vector<double> p_ell(Ylm_1.nOrders());
  return p_ell;
}


// ============================================================================


void cbl::glob::SphericalHarmonics_Coefficients::initialize (const int norder, const int nbins)
{

  m_nbins = nbins;
  m_norder = norder;
  m_lmax = m_norder-1;
  m_n_sph = gsl_sf_legendre_array_n(m_lmax);

  vector<vector<complex<double>>> _alm(m_nbins, vector<complex<double>>(m_n_sph, 0));

  vector<double> normalization (m_n_sph);

  vector<double> Plm(m_n_sph, 0);
  vector<complex<double>> sph(m_n_sph);

  for (int l=0; l<m_norder; l++)
    for (int m=0; m<l+1; m++) {
      int n=l*(l+1)/2+m;
      normalization[n] = gsl_sf_fact (l - m) / gsl_sf_fact ( l + m); 
    }

  m_alm = _alm;
  m_normalization = normalization;

  m_Plm = Plm;
  m_sph = sph;
}


// ============================================================================


void cbl::glob::SphericalHarmonics_Coefficients::reset ()
{
  for (int b=0; b<m_nbins; b++)
    for (int k1=0; k1<m_n_sph; k1++)
	m_alm[b][k1] = 0;
}


// ============================================================================


std::vector<std::complex<double>> cbl::glob::SphericalHarmonics_Coefficients::alm (const double xx, const double yy, const double zz)
{
  double cosTheta = zz;

  m_sph[0].real(1);
  m_sph[0].imag(0);

  m_sph[1].real(cosTheta);
  m_sph[1].imag(0);

  double Plm0 = 1.;
  double Plm1 = cosTheta;
  double Plm2, Plm3;

  int q, l, lp1;

  // P(l, 0)
  for (int l = 1; l < m_norder-1; l++) {
    q = l+1;
    Plm2 = ((2*l+1)*cosTheta*Plm1-l*Plm0)/(l+1);
    m_sph[q*(q+1)/2].real(Plm2);
    m_sph[q*(q+1)/2].imag(0.);
    Plm0 = Plm1;
    Plm1 = Plm2;
  }

  double sinTheta = sqrt(1 - cosTheta*cosTheta);
  if (sinTheta > 0) {
    double cosPhi = xx/sinTheta;
    double sinPhi = yy/sinTheta;

    double old_cc_real = 1.;
    double old_cc_imag = 0.;
    double cc_real = cosPhi;
    double cc_imag = sinPhi;

    // P(m, m) -> P(m+1, m) -> P( m+1<l<lMax, m)
    Plm0 = 1.;
    for (int m=1; m<m_norder; m++) {
      l = m;
      lp1 = m+1;
      Plm1 = -(2*(m-1)+1)*Plm0*sinTheta;
      Plm0 = Plm1;
      m_sph[l*(l+1)/2+m].real(Plm1*cc_real);
      m_sph[l*(l+1)/2+m].imag(Plm1*cc_imag);

      Plm2 = cosTheta*(2*l+1)*Plm1;
      m_sph[lp1*(lp1+1)/2+m].real(Plm2*cc_real);
      m_sph[lp1*(lp1+1)/2+m].imag(Plm2*cc_imag);

      for (int lp2=lp1+1; lp2<m_norder; lp2++) {
	Plm3 = ((2.*lp2-1.)*cosTheta*Plm2-(lp2+m-1)*Plm1)/(lp2-m);
	Plm1 = Plm2;
	Plm2 = Plm3;

	m_sph[lp2*(lp2+1)/2+m].real(Plm3*cc_real);
	m_sph[lp2*(lp2+1)/2+m].imag(Plm3*cc_imag);
      }

      old_cc_real = cc_real;
      old_cc_imag = cc_imag;
      cc_real = old_cc_real*cosPhi-old_cc_imag*sinPhi;
      cc_imag = old_cc_real*sinPhi+old_cc_imag*cosPhi;
    }
  }
  else {
    for (int m=1; m<m_norder; m++) {
      l = m;
      lp1 = m+1;

      m_sph[l*(l+1)/2+m].real(0.);
      m_sph[l*(l+1)/2+m].imag(0.);

      m_sph[lp1*(lp1+1)/2+m].real(0.);
      m_sph[lp1*(lp1+1)/2+m].imag(0.);

      for (int lp2=lp1+1; lp2<m_norder; lp2++) {
	Plm3 = ((2.*lp2-1.)*cosTheta*Plm2-(lp2+m-1)*Plm1)/(lp2-m);

	m_sph[lp2*(lp2+1)/2+m].real(0.);
	m_sph[lp2*(lp2+1)/2+m].imag(0.);
      }
    }
  }


  return m_sph;
}


// ============================================================================


void cbl::glob::SphericalHarmonics_Coefficients::add (const std::vector<std::complex<double>> alm, const double ww, const int bin)
{
  for(int n=0; n<m_n_sph; n++)
    m_alm[bin][n] += ww*alm[n];
}


// ============================================================================


void cbl::glob::SphericalHarmonics_Coefficients::add (const double xx, const double yy, const double zz, const double ww, const int bin)
{
  vector<complex<double>> alm = this->alm (xx, yy, zz);
  for(int n=0; n<m_n_sph; n++)
    m_alm[bin][n] += ww*alm[n];
}


// ============================================================================


double cbl::glob::SphericalHarmonics_Coefficients::power (const int l, const int bin1, const int bin2)
{
  const int min_n = l*(l+1)/2;
  double power = m_normalization[min_n]*(real(min_n, bin1)*real(min_n, bin2)+imag(min_n, bin1)*imag(min_n, bin2));
  for (int m=1; m<l+1; m++) {
    const int pos = min_n+m;
    power += 2.*m_normalization[pos]*(real(pos, bin1)*real(pos, bin2)+imag(pos, bin1)*imag(pos, bin2));
  }

  return power;
}
