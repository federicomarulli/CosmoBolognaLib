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

using namespace std;

using namespace cbl;


// ============================================================================================


double cbl::modelling::twopt::Pkmu_DeWiggled (const double kk, const double mu, const double alpha_perp, const double alpha_par, const double sigmaNL_perp, const double sigmaNL_par, const double linear_growth_rate, const double bias, const double SigmaS, const std::shared_ptr<cbl::glob::FuncGrid> Pk, const std::shared_ptr<cbl::glob::FuncGrid> Pk_NW)
{
  const double beta = linear_growth_rate/bias;
  const double FF = alpha_par/alpha_perp;
  const double fact = sqrt(1.+mu*mu*(pow(FF, -2)-1));

  const double kp = kk/alpha_perp*fact;
  const double mup = mu/FF/fact;

  const double KaiserBoost = pow(1.+mup*mup*beta, 2);
  const double sigmaNL2 = 0.5*((1.-mup*mup)*sigmaNL_perp*sigmaNL_perp+mu*mu*sigmaNL_par*sigmaNL_par);

  const double Fstreaming = pow(1+kp*kp*mup*mup*linear_growth_rate*linear_growth_rate*SigmaS*SigmaS, -2);

  const double sigmaNL = sqrt(sigmaNL_perp*sigmaNL_perp+sigmaNL_par*sigmaNL_par);
  const double _Pk = (sigmaNL<1.e-5) ? Pk->operator()(kp) : (Pk->operator()(kp)-Pk_NW->operator()(kp))*exp(-kk*kk*sigmaNL2)+Pk_NW->operator()(kp);

  return bias*bias*KaiserBoost*Fstreaming*_Pk;
}


// ============================================================================================


double cbl::modelling::twopt::Pkmu_ModeCoupling (const double kk, const double mu, const double alpha_perp, const double alpha_par, const double linear_growth_rate, const double bias, const double SigmaV, const double AMC, const std::shared_ptr<cbl::glob::FuncGrid> Pk, const std::shared_ptr<cbl::glob::FuncGrid> Pk_1loop)
{
  const double beta = linear_growth_rate/bias;
  const double FF = alpha_par/alpha_perp;
  const double fact = sqrt(1.+mu*mu*(pow(FF, -2)-1));

  const double kp = kk/alpha_perp*fact;
  const double mup = mu/FF/fact;

  const double KaiserBoost = pow(1.+mup*mup*beta, 2);
  const double Fstreaming = pow(1+kp*kp*mup*mup*linear_growth_rate*linear_growth_rate*SigmaV*SigmaV, -2);
  double Pk_NL = bias*bias*(Pk->operator()(kp)*exp(-kp*kp*SigmaV*SigmaV));
  if (kp<5)
    Pk_NL += bias*bias*AMC*Pk_1loop->operator()(kp)*pow(2.*par::pi, -3);

  return KaiserBoost*Fstreaming*Pk_NL;
}


// ============================================================================================


double cbl::modelling::twopt::Pkmu (const double kk, const double mu, const std::string model, const std::vector<double> parameter, const std::vector<std::shared_ptr<glob::FuncGrid>> pk_interp)
{
  if (model=="dispersion_dewiggled") {
    if (parameter.size()!=7)
      ErrorCBL("Error in cbl::modelling::twopt::Pkmu() of ModelFunction_TwoPointCorrelation.cpp: the "+model+" model has 7 parameters, while in the parameter vector in input has "+conv(parameter.size(), par::fINT)+" parameters!");
    return Pkmu_DeWiggled(kk, mu, parameter[0], parameter[1], parameter[2], parameter[3], parameter[4], parameter[5], parameter[6], pk_interp[0], pk_interp[1]);
  }

  else if (model=="dispersion_modecoupling") {
    if (parameter.size()!=6)
      ErrorCBL("Error in cbl::modelling::twopt::Pkmu() of ModelFunction_TwoPointCorrelation.cpp: the "+model+" model has 6 parameters, while in the parameter vector in input has "+conv(parameter.size(), par::fINT)+" parameters!");
    return Pkmu_ModeCoupling(kk, mu, parameter[0], parameter[1], parameter[2], parameter[3], parameter[4], parameter[5], pk_interp[0], pk_interp[1]);
  }
  
  else 
    ErrorCBL("Error in cbl::modelling::twopt::Pkmu() of ModelFunction_TwoPointCorrelation.cpp: the chosen model ("+model+") is not currently implemented!");
  
  return 0.;
}


// ============================================================================================


double cbl::modelling::twopt::Pk_l (const double kk, const int l, const std::string model, const std::vector<double> parameter, const std::vector<std::shared_ptr<glob::FuncGrid>> pk_interp, const double prec)
{
  auto integrand = [&] (const double mu)
  {
    return modelling::twopt::Pkmu(kk, mu, model, parameter, pk_interp)*legendre_polynomial(mu, l);
  };

  return 0.5*(2*l+1)*wrapper::gsl::GSL_integrate_qag(integrand, -1., 1., prec);
}

// ============================================================================================


std::vector<double> cbl::modelling::twopt::Pk_l (const std::vector<double> kk, const int l, const std::string model, const std::vector<double> parameter, const std::vector<std::shared_ptr<glob::FuncGrid>> pk_interp, const double prec)
{
  vector<double> Pkl(kk.size(), 0);

  for (size_t i=0; i<kk.size(); i++) {

    auto integrand = [&] (const double mu)
    {
      return modelling::twopt::Pkmu(kk[i], mu, model, parameter, pk_interp)*legendre_polynomial(mu, l);
    };

    Pkl[i] = 0.5*(2*l+1)*wrapper::gsl::GSL_integrate_qag(integrand, -1., 1., prec);
  }

  return Pkl;
}


// ============================================================================================


cbl::glob::FuncGrid cbl::modelling::twopt::Xil_interp (const std::vector<double> kk, const int l, const std::string model, const std::vector<double> parameter, const std::vector<std::shared_ptr<glob::FuncGrid>> pk_interp, const double prec)
{
  vector<double> Pkl = Pk_l(kk, l, model, parameter, pk_interp, prec);
  vector<double> rr, Xil;

  cbl::wrapper::fftlog::transform_FFTlog(rr, Xil, 1, kk, Pkl, l);

  cbl::glob::FuncGrid interp(rr, Xil, "Spline");

  return interp;
}


// ============================================================================================


std::vector<std::vector<double>> cbl::modelling::twopt::Xi_l (const std::vector<double> rr, const int nmultipoles, const std::string model, const std::vector<double> parameter, const std::vector<std::shared_ptr<glob::FuncGrid>> pk_interp, const double prec)
{
  vector<vector<double>> Xil(3);

  for (int i=0; i<nmultipoles; i++) {
    double sign = (i%2==0) ? 1 : -1;
    Xil[i] = Xil_interp(pk_interp[0]->x(), 2*i, model, parameter, pk_interp, prec).eval_func(rr);
    for (size_t j=0; j<rr.size(); j++)
      Xil[i][j] *= sign;
  }

  return Xil;
}


// ============================================================================================


std::vector<double> cbl::modelling::twopt::Xi_l (const std::vector<double> rr, const std::vector<int> dataset_order, const std::vector<bool> use_pole, const std::string model, const std::vector<double> parameter, const std::vector<std::shared_ptr<glob::FuncGrid>> pk_interp, const double prec)
{
  vector<cbl::glob::FuncGrid> interp_Xil(3);
  vector<double> sign={1., -1., 1.};
  
  for (size_t i=0; i<3; i++){
    if (use_pole[i])
      interp_Xil[i] = Xil_interp(pk_interp[0]->x(), 2*i, model, parameter, pk_interp, prec);
  }

  vector<double> Xil(rr.size());
  
  for (size_t i=0; i<rr.size(); i++) 
    Xil[i] = sign[dataset_order[i]]*interp_Xil[dataset_order[i]](rr[i]);

  return Xil;
}


// ============================================================================================


std::vector<std::vector<double>> cbl::modelling::twopt::Xi_rppi (const std::vector<double> rp, const std::vector<double> pi, const std::string model, const std::vector<double> parameter, const std::vector<std::shared_ptr<glob::FuncGrid>> pk_interp, const double prec)
{
  const int nr=200;
  const int nmultipoles = 3;
  vector<double> sign = {1., -1., 1.};
  vector<double> rr = linear_bin_vector(nr, min(Min(rp), Min(pi))*0.999, sqrt(Max(rp)*Max(rp)+Max(pi)*Max(pi))*1.001);

  vector<glob::FuncGrid> interp_Xil(nmultipoles);

  for (size_t i=0; i<nmultipoles; i++)
    interp_Xil[i] = Xil_interp(pk_interp[0]->x(), 2*i, model, parameter, pk_interp, prec);

  vector<vector<double>> xi_rppi(rp.size(), vector<double>(pi.size(), 0));

  for (size_t i =0; i<rp.size(); i++)
    for (size_t j =0; j<pi.size(); j++) {
      double s = sqrt(rp[i]*rp[i]+pi[j]*pi[j]);
      double mu = pi[j]/s;
      for (int l=0; l<nmultipoles; l++)
        xi_rppi[i][j] += sign[i]*interp_Xil[l](s)*legendre_polynomial (mu, l*2);
    }

  return xi_rppi;
}


// ============================================================================================


std::vector<double> cbl::modelling::twopt::wp_from_Xi_rppi (const std::vector<double> rp, const double pimax, const std::string model, const std::vector<double> parameter, const std::vector<std::shared_ptr<glob::FuncGrid>> pk_interp, const double prec) {

  vector<double> pi = linear_bin_vector(100, 1.e-4, pimax*1.001);

  vector<vector<double>> xi_rppi = modelling::twopt::Xi_rppi(rp, pi, model, parameter, pk_interp, prec);
  vector<double> wp(rp.size());

  for (size_t i=0; i<rp.size(); i++) {
    glob::FuncGrid func(pi, xi_rppi[i], "Spline");
    wp[i] = func.integrate_qag(0., pimax);
  }

  return wp;
}


// ============================================================================================


std::vector<std::vector<double>> cbl::modelling::twopt::damped_Pk_terms (const std::vector<double> kk, const double linear_growth_rate, const double SigmaS, const std::shared_ptr<cbl::glob::FuncGrid> PkDM)
{
  (void)PkDM;
  double sqrt_pi = sqrt(par::pi);
  vector<vector<double>> pk(3, vector<double>(kk.size(), 0));

  for (size_t i=0; i< kk.size(); i++) {
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


std::vector<double> cbl::modelling::twopt::damped_Xi (const std::vector<double> ss, const double bias, const double linear_growth_rate, const double SigmaS, const std::vector<double> kk, const std::shared_ptr<cbl::glob::FuncGrid> PkDM)
{
  vector<vector<double>> pk_terms = modelling::twopt::damped_Pk_terms(kk, linear_growth_rate, SigmaS, PkDM);

  vector<double> xi(ss.size(), 0);

  for (size_t i=0; i<pk_terms.size(); i++) {
    vector<double> xi_term = wrapper::fftlog::transform_FFTlog(ss, 1, kk, pk_terms[i], 0);
    for (size_t j=0; j<ss.size(); j++)
      xi[j] += pow(bias, 2-i)*xi_term[j];
    
  }

  return xi;
}
