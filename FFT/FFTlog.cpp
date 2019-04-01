/*******************************************************************
 *  Copyright (C) 2016 by Alfonso Veropalumbo                      *
 *  alfonso.veropalumbo@unibo.it                                   *
 *                                                                 *
 *  This program is free software; you can redistribute it and/or  *
 *  modify it under the terms of the GNU General Public License as *
 *  published by the Free Software Foundation; either version 2 of *
 *  the License, or (at your option) any later version.            *
 *                                                                 *
 *  This program is distributed in the hope that it will be useful,*
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the  *
 *  GNU General Public License for more details.                   *
 *                                                                 *
 *  You should have received a copy of the GNU General Public      *
 *  License along with this program; if not, write to the Free     *
 *  Software Foundation, Inc.,                                     *
 *  59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.      *
 *******************************************************************/

/**
 *  @file FFT/FFTlog.cpp
 *
 *  @brief Wrapper for fftlog wripper
 *
 *  This file contains the prototypes of the class FFTlog,
 *  wrapping functions for fast Hankel transform 
 *  implemented in the FFTlog library (Hamilton 2000)
 *
 *  @author Alfonso Veropalumbo 
 *
 *  @author alfonso.veropalumbo@unbo.it
 */

#include "FFTlog.h"

using namespace std;


// ============================================================================


vector<double> cbl::wrapper::fftlog::transform_FFTlog (const std::vector<double> yy, const int dir, const std::vector<double> xx, const std::vector<double> fx, const double mu, const double q, const double kr, const int kropt)
{
  vector<double> _yy, _fy;
  cbl::wrapper::fftlog::transform_FFTlog(_yy, _fy, dir, xx, fx, mu, q, kr, kropt);

  cbl::glob::FuncGrid interp(_yy, _fy, "Spline");

  vector<double> interpolated_values = interp.eval_func(yy);

  interp.free();

  return interpolated_values;
}


// ============================================================================


void cbl::wrapper::fftlog::transform_FFTlog (std::vector<double> &yy, std::vector<double> &fy, const int dir, const std::vector<double> xx, const std::vector<double> fx, const double mu, const double q, const double kr, const int kropt)
{
  const int NMAX = 4096;

  double fact = 1./(2*cbl::par::pi*cbl::par::pi)*sqrt(cbl::par::pi/2);

  int i_ok;
  double wsave[2*NMAX+3*(NMAX/2)+19];

  int n = fx.size();
  double _mu = mu+0.5;
  double _q = q;
  double _kr = kr;
  int _kropt = kropt;
  int _dir = dir;


  double xmin = log10(cbl::Min(xx)), xmax = log10(cbl::Max(xx));
  double dlogx = (xmax-xmin)/(n-1);
  double dlnx = dlogx*log(10.);

  double ci = double(n+1)/2;

  double logxmedian = (xmax+xmin)/2;

  double ap[NMAX];
  for (int i=0; i<n; i++)
    ap[i] = fx[i]*xx[i];

  wrapper::fftlog::fhti_(&n, &_mu, &_q, &dlnx, &_kr, &_kropt, wsave, &i_ok);

  double logymedian = log10(_kr)-logxmedian;

  double rk = pow(10, logxmedian - logymedian);

  if (i_ok==0)
    ErrorCBL("Problems in cbl::wrapper::fftlog::transform_FFTlog() of FFTlog.cpp!");

  wrapper::fftlog::fftl_(&n, ap, &rk, &_dir, wsave);

  yy.erase(yy.begin(), yy.end()); fy.erase(fy.begin(), fy.end());

  for (int i=0; i<n; i++) { 
    double _y = pow(10., logymedian + ( i+1 - ci ) * dlogx );
    yy.push_back(_y);
    fy.push_back(fact*ap[i]/_y);
  }
}
