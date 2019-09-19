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
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#include "SphericalHarmonics_Coefficients.h"

using namespace std;

using namespace cbl;

// ============================================================================================


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
