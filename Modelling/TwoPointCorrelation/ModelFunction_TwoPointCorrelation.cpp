/*******************************************************************
 *  Copyright (C) 2016 by Federico Marulli and Alfonso Veropalumbo *
 *  federico.marulli3@unibo.it                                     *
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
 *  @file
 *  ModelFunction/TwoPointCorrelation/ModelFunction_TwoPointCorrelation.cpp
 *
 *  @brief Global functions to model two-point correlation functions
 *  of any type
 *
 *  This file contains the implementation of the used to model
 *  two-point correlation functions of any type
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#include "ModelFunction_TwoPointCorrelation.h"

using namespace cosmobl;


// ============================================================================================


double cosmobl::modelling::twopt::Pkmu (const double kk, const double mu, const double alpha_perp, const double alpha_par, const double sigmaNL_perp, const double sigmaNL_par, const double bias, const double linear_growth_rate, const double SigmaS, const shared_ptr<cosmobl::glob::FuncGrid> Pk, const shared_ptr<cosmobl::glob::FuncGrid> Pk_NW)
{
  const double beta = linear_growth_rate/bias;
  const double FF = alpha_par/alpha_perp;
  double fact = sqrt(1.+mu*mu*(pow(FF, -2)-1));

  const double kp = kk/alpha_perp*fact;
  const double mup = mu/FF/fact;

  const double KaiserBoost = pow(1.+mup*mup*beta, 2);
  const double sigmaNL2 = 0.5*((1.-mup*mup)*sigmaNL_perp*sigmaNL_perp+mup*mup*sigmaNL_par*sigmaNL_par);

  const double Fstreaming = pow(1+kp*kp*mup*mup*SigmaS*SigmaS, -2);

  const double sigmaNL = sqrt(sqrt(sigmaNL_perp*sigmaNL_perp+sigmaNL_par*sigmaNL_par));
  const double _Pk = (sigmaNL<1.e-5) ? Pk->operator()(kp) : (Pk->operator()(kp)-Pk_NW->operator()(kp))*exp(-kp*kp*sigmaNL2)+Pk_NW->operator()(kp);

  return bias*bias*KaiserBoost*Fstreaming*_Pk;
}


// ============================================================================================


double cosmobl::modelling::twopt::Pk_l (const double kk, const int order, const double alpha_perp, const double alpha_par, const double sigmaNL_perp, const double sigmaNL_par, const double bias, const double linear_growth_rate, const double SigmaS, const shared_ptr<cosmobl::glob::FuncGrid> Pk, const shared_ptr<cosmobl::glob::FuncGrid> Pk_NW, const double prec)
{
  auto integrand = [&] (const double mu)
  {
    return cosmobl::modelling::twopt::Pkmu(kk, mu, alpha_perp, alpha_par, sigmaNL_perp, sigmaNL_par, bias, linear_growth_rate, SigmaS, Pk, Pk_NW)*cosmobl::legendre_polynomial(mu, order);
  };

  return 0.5*(2*order+1)*cosmobl::gsl::GSL_integrate_qag(integrand, -1., 1., prec);
}


// ============================================================================================


vector<vector<double>> cosmobl::modelling::twopt::Xi_l(const vector<vector<double>> rr, const int nmultipoles, const double alpha_perp, const double alpha_par, const double sigmaNL_perp, const double sigmaNL_par, const double bias, const double linear_growth_rate, const double SigmaS, const vector<double> kk, const shared_ptr<cosmobl::glob::FuncGrid> Pk, const shared_ptr<cosmobl::glob::FuncGrid> Pk_NW, const double prec)
{
  checkDim(rr, nmultipoles, "rad");
  vector<vector<double>> Xil(nmultipoles);

  for(int i=0; i<nmultipoles; i++){
    int order = i*2;
    const int sign = pow(complex<double>(0.,1.), order).real();

    vector<double> Pkl;

    for(size_t j=0; j< kk.size(); j++)
      Pkl.push_back(cosmobl::modelling::twopt::Pk_l(kk[j], order, alpha_perp, alpha_par, sigmaNL_perp, sigmaNL_par, bias, linear_growth_rate, SigmaS, Pk, Pk_NW, prec));

    vector<double> xx = cosmobl::fftlog::transform_FFTlog(rr[i], 1, kk, Pkl, order);
    for(size_t j=0; j<rr[i].size(); j++)
      Xil[i].push_back(sign*xx[j]);
  }

  return Xil;
}


// ============================================================================================


vector<vector<double>> cosmobl::modelling::twopt::Xi_l(const vector<double> rr, const int nmultipoles, const double alpha_perp, const double alpha_par, const double sigmaNL_perp, const double sigmaNL_par, const double bias, const double linear_growth_rate, const double SigmaS, const vector<double> kk, const shared_ptr<cosmobl::glob::FuncGrid> Pk, const shared_ptr<cosmobl::glob::FuncGrid> Pk_NW, const double prec)
{
  vector<vector<double>> rad(nmultipoles, rr);

  return Xi_l(rad, nmultipoles, alpha_perp, alpha_par, sigmaNL_perp, sigmaNL_par, bias, linear_growth_rate, SigmaS, kk, Pk, Pk_NW, prec);
}

// ============================================================================================


vector<vector<double>> cosmobl::modelling::twopt::Xi_wedges (const vector<double> rr, const int nwedges, const double alpha_perp, const double alpha_par, const double sigmaNL_perp, const double sigmaNL_par, const double bias, const double linear_growth_rate, const double SigmaS, const vector<double> kk, const shared_ptr<cosmobl::glob::FuncGrid> Pk, const shared_ptr<cosmobl::glob::FuncGrid> Pk_NW, const double prec)
{
  vector<vector<double>> Xil = cosmobl::modelling::twopt::Xi_l(rr, 3, alpha_perp, alpha_par, sigmaNL_perp, sigmaNL_par, bias, linear_growth_rate, SigmaS, kk, Pk, Pk_NW, prec);

  vector<vector<double>> XiW(nwedges, vector<double>(rr.size(),0));

  for(int i=0; i<nwedges; i++){
    double mu_min = i*1./nwedges;
    double mu_max = (i+1)*1./nwedges;

    double f2 = 0.5*((pow(mu_max, 3)-pow(mu_min,3))/(mu_max-mu_min)-1.);
    double f4 = 0.125*( ( 7.*(pow(mu_max, 5)-pow(mu_min,5)) - 10.*(pow(mu_max, 3)-pow(mu_min,3)))/(mu_max-mu_min)+3.);

    for(size_t j=0; j<rr.size(); j++)
      XiW[i][j] = Xil[0][j]+f2*Xil[1][j]+f4*Xil[2][j];
  }

  return XiW;
}


// ============================================================================================


vector<vector<double>> cosmobl::modelling::twopt::Xi_rppi(const vector<double> rp, const vector<double> pi, const double alpha_perp, const double alpha_par, const double sigmaNL_perp, const double sigmaNL_par, const double bias, const double linear_growth_rate, const double SigmaS, const vector<double> kk, const shared_ptr<cosmobl::glob::FuncGrid> Pk, const shared_ptr<cosmobl::glob::FuncGrid> Pk_NW, const double prec)
{
  const int nr=200;
  const int nmultipoles = 3;
  vector<double> rr = linear_bin_vector(nr, min(Min(rp), Min(pi))*0.999, sqrt(Max(rp)*Max(rp)+Max(pi)*Max(pi))*1.001);
  
  vector<vector<double>> Xil = Xi_l(rr, nmultipoles, alpha_perp, alpha_par, sigmaNL_perp, sigmaNL_par, bias, linear_growth_rate, SigmaS, kk, Pk, Pk_NW, prec); 

  vector<shared_ptr<cosmobl::glob::FuncGrid>> Xil_interp(nmultipoles);

  for(size_t i=0; i<nmultipoles; i++)
    Xil_interp[i] = make_shared<cosmobl::glob::FuncGrid>(cosmobl::glob::FuncGrid(rr, Xil[i], "Spline"));

  vector<vector<double>> xi_rppi(rp.size(), vector<double>(pi.size(), 0));

  for(size_t i =0; i<rp.size(); i++)
    for(size_t j =0; j<pi.size(); j++){
      double s = sqrt(rp[i]*rp[i]+pi[j]*pi[j]);
      double mu = pi[j]/s;
      for(int l=0; l<nmultipoles; l++)
        xi_rppi[i][j] += Xil_interp[l]->operator()(s)*cosmobl::legendre_polynomial (mu, l*2);
    }

  return xi_rppi;
}


// ============================================================================================


vector<double> cosmobl::modelling::twopt::wp_from_Xi_rppi(const vector<double> rp, const double pimax, const double alpha_perp, const double alpha_par, const double sigmaNL_perp, const double sigmaNL_par, const double bias, const double linear_growth_rate, const double SigmaS, const vector<double> kk, const shared_ptr<cosmobl::glob::FuncGrid> Pk, const shared_ptr<cosmobl::glob::FuncGrid> Pk_NW, const double prec){

  vector<double> pi = linear_bin_vector(100, 1.e-4, pimax*1.001);

  vector<vector<double>> xi_rppi = cosmobl::modelling::twopt::Xi_rppi(rp, pi, alpha_perp, alpha_par, sigmaNL_perp, sigmaNL_par, bias, linear_growth_rate, SigmaS, kk, Pk, Pk_NW, prec);
  vector<double> wp(rp.size());

  for(size_t i=0; i<rp.size(); i++){
    cosmobl::glob::FuncGrid func(pi, xi_rppi[i], "Spline");
    wp[i] = func.integrate_qag(0., pimax);
  }

  return wp;
}


// ============================================================================================


vector<vector<double>> cosmobl::modelling::twopt::damped_Pk_terms (const vector<double> kk, const double linear_growth_rate, const double SigmaS, const shared_ptr<cosmobl::glob::FuncGrid> PkDM)
{
  (void)PkDM;
  double sqrt_pi = sqrt(par::pi);
  vector<vector<double>> pk(3, vector<double>(kk.size(), 0));

  for (size_t i=0; i< kk.size(); i++){
    double kSigma = kk[i]*SigmaS;
    double pk_dm = PkDM->operator()(kk[i]);
    double erf_kSigma = erf(kSigma);

    pk[0][i] = pk_dm*sqrt_pi/(2.*kSigma)*erf_kSigma;
    pk[1][i] = linear_growth_rate*pow(kSigma, -3.)*pk_dm*(0.5*sqrt_pi*erf_kSigma-kSigma*exp(-kSigma*kSigma));
    pk[2][i] = linear_growth_rate*linear_growth_rate*pow(kSigma, -5.)*pk_dm*(3.*sqrt_pi/8.*erf_kSigma-0.25*kSigma*(2*kSigma*kSigma+3)*exp(-kSigma*kSigma));
  }
  

  return pk;
}

// ============================================================================================


vector<double> cosmobl::modelling::twopt::damped_Xi (const vector<double> ss, const double bias, const double linear_growth_rate, const double SigmaS, const vector<double> kk, const shared_ptr<cosmobl::glob::FuncGrid> PkDM)
{
  vector<vector<double>> pk_terms = cosmobl::modelling::twopt::damped_Pk_terms(kk, linear_growth_rate, SigmaS, PkDM);

  vector<double> xi(ss.size(), 0);

  for (size_t i=0; i<pk_terms.size(); i++){
    vector<double> xi_term = cosmobl::fftlog::transform_FFTlog(ss, 1, kk, pk_terms[i], 0);
    for (size_t j=0; j<ss.size(); j++)
      xi[j] += pow(bias, 2-i)*xi_term[j];
    
  }

  return xi;
}
