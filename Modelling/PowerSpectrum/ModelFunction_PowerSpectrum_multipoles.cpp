/********************************************************************
 *  Copyright (C) 2016 by Federico Marulli and Alfonso Veropalumbo  *
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
 *  @file
 *  Modelling/PowerSpectrum/ModelFunction_PowerSpectrum_multipoles.cpp
 *
 *  @brief Functions to model the multipoles of the power spectrum
 *
 *  This file contains the implementation of the functions used to
 *  model the multipoles of the power spectrum
 */


#include "ModelFunction_PowerSpectrum_multipoles.h"

using namespace std;

using namespace cbl;


// ============================================================================================


std::vector<double> cbl::modelling::powspec::true_k_mu_AP (const double kk, const double mu, const double alpha_perp, const double alpha_par)
{
  const double FF = alpha_par/alpha_perp;
  const double fact = sqrt(1.+mu*mu*(1./(FF*FF)-1.));

  return {kk/alpha_perp*fact, mu/FF/fact};
}


// ============================================================================================


double cbl::modelling::powspec::Pkmu_DeWiggled (const double kk, const double mu, const double sigmaNL_perp, const double sigmaNL_par, const double linear_growth_rate, const double bias, const double SigmaS, const std::shared_ptr<cbl::glob::FuncGrid> Pk, const std::shared_ptr<cbl::glob::FuncGrid> Pk_NW)
{
  const double beta = linear_growth_rate/bias;

  const double KaiserBoost = pow(1.+mu*mu*beta, 2);
  const double sigmaNL2 = 0.5*((1.-mu*mu)*sigmaNL_perp*sigmaNL_perp+mu*mu*sigmaNL_par*sigmaNL_par);

  const double Fstreaming = pow(1+kk*kk*mu*mu*linear_growth_rate*linear_growth_rate*SigmaS*SigmaS, -2);

  const double sigmaNL = sqrt(sigmaNL_perp*sigmaNL_perp+sigmaNL_par*sigmaNL_par);
  const double _Pk = (sigmaNL<1.e-5) ? Pk->operator()(kk) : (Pk->operator()(kk)-Pk_NW->operator()(kk))*exp(-kk*kk*sigmaNL2)+Pk_NW->operator()(kk);

  return bias*bias*KaiserBoost*Fstreaming*_Pk;
}


// ============================================================================================


double cbl::modelling::powspec::Pkmu_ModeCoupling (const double kk, const double mu, const double linear_growth_rate, const double bias, const double SigmaV, const double AMC, const std::shared_ptr<cbl::glob::FuncGrid> Pk, const std::shared_ptr<cbl::glob::FuncGrid> Pk_1loop)
{
  const double beta = linear_growth_rate/bias;

  const double KaiserBoost = pow(1.+mu*mu*beta, 2);
  const double Fstreaming = pow(1+kk*kk*mu*mu*linear_growth_rate*linear_growth_rate*SigmaV*SigmaV, -2);
  double Pk_NL = bias*bias*(Pk->operator()(kk)*exp(-kk*kk*SigmaV*SigmaV));
  
  if (kk<5)
    Pk_NL += bias*bias*AMC*Pk_1loop->operator()(kk)*pow(2.*par::pi, -3);

  return KaiserBoost*Fstreaming*Pk_NL;
}


// ============================================================================================


double cbl::modelling::powspec::Pkmu_dispersion (const double kk, const double mu, const std::string DFoG, const double linear_growth_rate, const double bias, const double sigmav, const std::shared_ptr<cbl::glob::FuncGrid> Pklin)
{
  double DispFactor;

  if (DFoG=="Gaussian")  
    DispFactor = exp(-kk*kk*mu*mu*linear_growth_rate*linear_growth_rate*sigmav*sigmav);

  else if (DFoG=="Lorentzian")
    DispFactor = 1./(1.+(kk*kk*mu*mu*linear_growth_rate*linear_growth_rate*sigmav*sigmav));

  else
    ErrorCBL("the chosen DFoG ("+DFoG+") is not currently implemented!", "Pkmu_dispersion", "ModelFunction_TwoPointCorrelation.cpp");

  return DispFactor*bias*bias*pow(1+linear_growth_rate/bias*mu*mu, 2)*Pklin->operator()(kk);
}


// ============================================================================================


double cbl::modelling::powspec::Pkmu_Scoccimarro (const double kk, const double mu, const std::string DFoG, const double linear_growth_rate, const double bias, const double sigmav, const std::shared_ptr<cbl::glob::FuncGrid> Pk_DeltaDelta, const std::shared_ptr<cbl::glob::FuncGrid> Pk_DeltaTheta, const std::shared_ptr<cbl::glob::FuncGrid> Pk_ThetaTheta)
{
  const double Pk_dd = Pk_DeltaDelta->operator()(kk);
  const double Pk_dt = Pk_DeltaTheta->operator()(kk);
  const double Pk_tt = Pk_ThetaTheta->operator()(kk);

  double DispFactor;
  
  if (DFoG=="Gaussian")
    DispFactor = exp(-kk*kk*mu*mu*linear_growth_rate*linear_growth_rate*sigmav*sigmav);
  
  else if (DFoG=="Lorentzian")
    DispFactor = 1./(1.+(kk*kk*mu*mu*linear_growth_rate*linear_growth_rate*sigmav*sigmav));

  else
    ErrorCBL("the chosen DFoG ("+DFoG+") is not currently implemented!", "Pkmu_Scoccimarro", "ModelFunction_TwoPointCorrelation.cpp");

  const double beta = linear_growth_rate/bias;
  
  return DispFactor*bias*bias*(Pk_dd+2.*beta*pow(mu, 2)*Pk_dt+pow(beta, 2)*pow(mu, 4)*Pk_tt);
}


// ============================================================================================


double cbl::modelling::powspec::Pkmu_Scoccimarro_fitPezzotta (const double kk, const double mu, const std::string DFoG, const double linear_growth_rate, const double bias, const double sigmav, const double kd, const double kt, const std::shared_ptr<cbl::glob::FuncGrid> Pklin, const std::shared_ptr<cbl::glob::FuncGrid> Pknonlin)
{
  const double beta = linear_growth_rate/bias;
  double DispFactor;
  double Pk_dd = Pknonlin->operator()(kk);
  double Pk_dt = sqrt(Pknonlin->operator()(kk)*Pklin->operator()(kk)*exp(-kk/kd));
  double Pk_tt = Pklin->operator()(kk)*exp(-kk/kt);

  if (DFoG=="Gaussian")
    DispFactor = exp(-kk*kk*mu*mu*linear_growth_rate*linear_growth_rate*sigmav*sigmav);

  else if (DFoG=="Lorentzian")  
    DispFactor = 1./(1.+(kk*kk*mu*mu*linear_growth_rate*linear_growth_rate*sigmav*sigmav));

  else
    ErrorCBL("the chosen DFoG ("+DFoG+") is not currently implemented!", "Pkmu_Scoccimarro_fitPezzotta", "ModelFunction_TwoPointCorrelation.cpp");

  return DispFactor*bias*bias*(Pk_dd+2.*beta*pow(mu, 2)*Pk_dt+pow(beta, 2)*pow(mu, 4)*Pk_tt);
}


// ============================================================================================


double cbl::modelling::powspec::Pkmu_Scoccimarro_fitBel (const double kk, const double mu, const std::string DFoG, const double linear_growth_rate, const double bias, const double sigmav, const double kd, const double bb, const double a1, const double a2, const double a3, const std::shared_ptr<cbl::glob::FuncGrid> Pklin, const std::shared_ptr<cbl::glob::FuncGrid> Pknonlin)
{
  const double beta = linear_growth_rate/bias;
  double DispFactor;
  double Pk_dd = Pknonlin->operator()(kk);
  double Pk_dt = sqrt(Pknonlin->operator()(kk)*Pklin->operator()(kk))*exp(-kk/kd-bb*pow(kk, 6.0));
  double Pk_tt = Pklin->operator()(kk)*exp(-kk*(a1 + a2*kk + a3*kk*kk));

  if (DFoG=="Gaussian")
    DispFactor = exp(-kk*kk*mu*mu*linear_growth_rate*linear_growth_rate*sigmav*sigmav);
  else if (DFoG=="Lorentzian")
    DispFactor = 1./(1.+(kk*kk*mu*mu*linear_growth_rate*linear_growth_rate*sigmav*sigmav));
  else
    ErrorCBL("the chosen DFoG ("+DFoG+") is not currently implemented!", "Pkmu_Scoccimarro_fitBel", "ModelFunction_TwoPointCorrelation.cpp");

  return DispFactor*bias*bias*(Pk_dd+2.*beta*pow(mu, 2)*Pk_dt+pow(beta, 2)*pow(mu, 4)*Pk_tt);
}


// ============================================================================================


double cbl::modelling::powspec::Pkmu_TNS (const double kk, const double mu, const std::string DFoG, const double linear_growth_rate, const double bias, const double sigmav, const std::shared_ptr<cbl::glob::FuncGrid> Pk_DeltaDelta, const std::shared_ptr<cbl::glob::FuncGrid> Pk_DeltaTheta, const std::shared_ptr<cbl::glob::FuncGrid> Pk_ThetaTheta, const std::shared_ptr<cbl::glob::FuncGrid> Pk_A11, const std::shared_ptr<cbl::glob::FuncGrid> Pk_A12, const std::shared_ptr<cbl::glob::FuncGrid> Pk_A22, const std::shared_ptr<cbl::glob::FuncGrid> Pk_A23, const std::shared_ptr<cbl::glob::FuncGrid> Pk_A33, const std::shared_ptr<cbl::glob::FuncGrid> Pk_B12, const std::shared_ptr<cbl::glob::FuncGrid> Pk_B13, const std::shared_ptr<cbl::glob::FuncGrid> Pk_B14, const std::shared_ptr<cbl::glob::FuncGrid> Pk_B22, const std::shared_ptr<cbl::glob::FuncGrid> Pk_B23, const std::shared_ptr<cbl::glob::FuncGrid> Pk_B24, const std::shared_ptr<cbl::glob::FuncGrid> Pk_B33, const std::shared_ptr<cbl::glob::FuncGrid> Pk_B34, const std::shared_ptr<cbl::glob::FuncGrid> Pk_B44)
{
  const double beta = linear_growth_rate/bias;

  const double Pk_dd  = Pk_DeltaDelta->operator()(kk);
  const double Pk_dt  = Pk_DeltaTheta->operator()(kk);
  const double Pk_tt  = Pk_ThetaTheta->operator()(kk);

  const double pkA11 = Pk_A11->operator()(kk);
  const double pkA12 = Pk_A12->operator()(kk);
  const double pkA22 = Pk_A22->operator()(kk);
  const double pkA23 = Pk_A23->operator()(kk);
  const double pkA33 = Pk_A33->operator()(kk);
  const double pkB12 = Pk_B12->operator()(kk);
  const double pkB13 = Pk_B13->operator()(kk);
  const double pkB14 = Pk_B14->operator()(kk);
  const double pkB22 = Pk_B22->operator()(kk);
  const double pkB23 = Pk_B23->operator()(kk);
  const double pkB24 = Pk_B24->operator()(kk);
  const double pkB33 = Pk_B33->operator()(kk);
  const double pkB34 = Pk_B34->operator()(kk);
  const double pkB44 = Pk_B44->operator()(kk);

  double DispFactor;
  if (DFoG=="Gaussian")
    DispFactor = exp(-kk*kk*mu*mu*linear_growth_rate*linear_growth_rate*sigmav*sigmav);
  else if (DFoG=="Lorentzian")
    DispFactor = 1./(1.+(kk*kk*mu*mu*linear_growth_rate*linear_growth_rate*sigmav*sigmav));
  else
    ErrorCBL("the chosen DFoG ("+DFoG+") is not currently implemented!", "Pkmu_TNS", "ModelFunction_TwoPointCorrelation.cpp");

  const double Pk_Scoccimarro_large_scales = bias*bias*(Pk_dd+2*beta*pow(mu, 2)*Pk_dt+pow(beta, 2)*pow(mu, 4)*Pk_tt);

  const double A2 = pow(mu, 2)*(beta*pkA11+pow(beta, 2)*pkA12);
  const double A4 = pow(mu, 4)*pow(beta, 2)*(pkA22+beta*pkA23);
  const double A6 = pow(mu, 6)*pow(beta, 3)*pkA33;
  const double B2 = pow(mu, 2)*(pow(beta, 2)*pkB12+pow(beta, 3)*pkB13+pow(beta, 4)*pkB14);
  const double B4 = pow(mu, 4)*(pow(beta, 2)*pkB22+pow(beta, 3)*pkB23+pow(beta, 4)*pkB24);
  const double B6 = pow(mu, 6)*(pow(beta, 3)*pkB33+pow(beta, 4)*pkB34);
  const double B8 = pow(mu, 8)*pow(beta, 4)*pkB44;

  const double Pk_A = pow(bias, 3)*(A2+A4+A6);
  const double Pk_B = pow(bias, 4)*(B2+B4+B6+B8);

  return DispFactor*(Pk_Scoccimarro_large_scales+Pk_A+Pk_B);
}


// ============================================================================================


double cbl::modelling::powspec::Pkmu_eTNS (const double kk, const double mu, const std::string DFoG, const double linear_growth_rate, const double bias, const double bias2, const double sigmav, const double Ncorr, const std::shared_ptr<cbl::glob::FuncGrid> Pk_DeltaDelta, const std::shared_ptr<cbl::glob::FuncGrid> Pk_DeltaTheta, const std::shared_ptr<cbl::glob::FuncGrid> Pk_ThetaTheta, const std::shared_ptr<cbl::glob::FuncGrid> Pk_A11, const std::shared_ptr<cbl::glob::FuncGrid> Pk_A12, const std::shared_ptr<cbl::glob::FuncGrid> Pk_A22, const std::shared_ptr<cbl::glob::FuncGrid> Pk_A23, const std::shared_ptr<cbl::glob::FuncGrid> Pk_A33, const std::shared_ptr<cbl::glob::FuncGrid> Pk_B12, const std::shared_ptr<cbl::glob::FuncGrid> Pk_B13, const std::shared_ptr<cbl::glob::FuncGrid> Pk_B14, const std::shared_ptr<cbl::glob::FuncGrid> Pk_B22, const std::shared_ptr<cbl::glob::FuncGrid> Pk_B23, const std::shared_ptr<cbl::glob::FuncGrid> Pk_B24, const std::shared_ptr<cbl::glob::FuncGrid> Pk_B33, const std::shared_ptr<cbl::glob::FuncGrid> Pk_B34, const std::shared_ptr<cbl::glob::FuncGrid> Pk_B44, const std::shared_ptr<cbl::glob::FuncGrid> Pk_b2d, const std::shared_ptr<cbl::glob::FuncGrid> Pk_b2v, const std::shared_ptr<cbl::glob::FuncGrid> Pk_b22, const std::shared_ptr<cbl::glob::FuncGrid> Pk_bs2d, const std::shared_ptr<cbl::glob::FuncGrid> Pk_bs2v, const std::shared_ptr<cbl::glob::FuncGrid> Pk_b2s2, const std::shared_ptr<cbl::glob::FuncGrid> Pk_bs22, const std::shared_ptr<cbl::glob::FuncGrid> sigma32Pklin)
{
  const double beta = linear_growth_rate/bias;

  const double pk_dd  = Pk_DeltaDelta->operator()(kk);
  const double pk_dt  = Pk_DeltaTheta->operator()(kk);
  const double pk_tt  = Pk_ThetaTheta->operator()(kk);

  const double pkA11 = Pk_A11->operator()(kk);
  const double pkA12 = Pk_A12->operator()(kk);
  const double pkA22 = Pk_A22->operator()(kk);
  const double pkA23 = Pk_A23->operator()(kk);
  const double pkA33 = Pk_A33->operator()(kk);
  const double pkB12 = Pk_B12->operator()(kk);
  const double pkB13 = Pk_B13->operator()(kk);
  const double pkB14 = Pk_B14->operator()(kk);
  const double pkB22 = Pk_B22->operator()(kk);
  const double pkB23 = Pk_B23->operator()(kk);
  const double pkB24 = Pk_B24->operator()(kk);
  const double pkB33 = Pk_B33->operator()(kk);
  const double pkB34 = Pk_B34->operator()(kk);
  const double pkB44 = Pk_B44->operator()(kk);

  const double pk_b2d = Pk_b2d->operator()(kk);
  const double pk_b2v = Pk_b2v->operator()(kk);
  const double pk_b22 = Pk_b22->operator()(kk);
  const double pk_bs2d = Pk_bs2d->operator()(kk);
  const double pk_bs2v = Pk_bs2v->operator()(kk);
  const double pk_b2s2 = Pk_b2s2->operator()(kk);
  const double pk_bs22 = Pk_bs22->operator()(kk);
  const double sigma32pklin = sigma32Pklin->operator()(kk);
  
  double DispFactor;
  if (DFoG=="Gaussian")
    DispFactor = exp(-kk*kk*mu*mu*linear_growth_rate*linear_growth_rate*sigmav*sigmav);
  else if (DFoG=="Lorentzian")
    DispFactor = 1./(1.+(kk*kk*mu*mu*linear_growth_rate*linear_growth_rate*sigmav*sigmav));
  else
    ErrorCBL("the chosen DFoG ("+DFoG+") is not currently implemented!", "Pkmu_TNS", "ModelFunction_TwoPointCorrelation.cpp");

  const double A2 = pow(mu, 2)*(beta*pkA11+pow(beta, 2)*pkA12);
  const double A4 = pow(mu, 4)*pow(beta, 2)*(pkA22+beta*pkA23);
  const double A6 = pow(mu, 6)*pow(beta, 3)*pkA33;
  const double B2 = pow(mu, 2)*(pow(beta, 2)*pkB12+pow(beta, 3)*pkB13+pow(beta, 4)*pkB14);
  const double B4 = pow(mu, 4)*(pow(beta, 2)*pkB22+pow(beta, 3)*pkB23+pow(beta, 4)*pkB24);
  const double B6 = pow(mu, 6)*(pow(beta, 3)*pkB33+pow(beta, 4)*pkB34);
  const double B8 = pow(mu, 8)*pow(beta, 4)*pkB44;

  const double Pk_A = pow(bias, 3)*(A2+A4+A6);
  const double Pk_B = pow(bias, 4)*(B2+B4+B6+B8);

  //----- NL bias terms -----
  const double bs2  = (-4./7.)*(bias-1.);
  const double b3nl = (32./315.)*(bias-1.);

  const double Pgdd = bias*bias*pk_dd+2*bias2*bias*pk_b2d+2*bs2*bias*pk_bs2d+2*b3nl*bias*sigma32pklin+bias2*bias2*pk_b22+2*bias2*bs2*pk_b2s2+bs2*bs2*pk_bs22+Ncorr;
  const double Pgdv = bias*pk_dt+bias2*pk_b2v+bs2*pk_bs2v+b3nl*sigma32pklin;
  const double Pk_Scoccimarro_large_scales = Pgdd+2.*linear_growth_rate*pow(mu, 2)*Pgdv+pow(linear_growth_rate, 2)*pow(mu, 4)*pk_tt;
  
  return DispFactor*(Pk_Scoccimarro_large_scales+Pk_A+Pk_B);
}


// ============================================================================================


double cbl::modelling::powspec::Pkmu (const double kk, const double mu, const std::string model, const std::vector<double> parameter, const std::vector<std::shared_ptr<glob::FuncGrid>> pk_interp, const double alpha_perp, const double alpha_par)
{
  const double kp = true_k_mu_AP(kk, mu, alpha_perp, alpha_par)[0];
  const double mup = true_k_mu_AP(kk, mu, alpha_perp, alpha_par)[1];

  double pkmu = 0;

  if (model=="dispersion_dewiggled") {
    if (parameter.size()!=5)
      ErrorCBL("the "+model+" model has 5 parameters, while in the parameter vector in input has "+conv(parameter.size(), par::fINT)+" parameters!", "Pkmu", "ModelFunction_TwoPointCorrelation.cpp");
    pkmu = Pkmu_DeWiggled(kp, mup, parameter[0], parameter[1], parameter[2], parameter[3], parameter[4], pk_interp[0], pk_interp[1]);
  }

  else if (model=="dispersion_modecoupling") {
    if (parameter.size()!=4)
      ErrorCBL("the "+model+" model has 4 parameters, while in the parameter vector in input has "+conv(parameter.size(), par::fINT)+" parameters!", "Pkmu", "ModelFunction_TwoPointCorrelation.cpp");
    pkmu = Pkmu_ModeCoupling(kp, mup, parameter[0], parameter[1], parameter[2], parameter[3], pk_interp[0], pk_interp[1]);
  }
  
  else if (model=="dispersion_Gauss") {
    if (parameter.size()!=3)
      ErrorCBL("the "+model+" model has 3 parameters, while in the parameter vector in input has "+conv(parameter.size(), par::fINT)+" parameters!", "Pkmu", "ModelFunction_TwoPointCorrelation.cpp");
    pkmu = Pkmu_dispersion(kp, mup, "Gaussian", parameter[0], parameter[1], parameter[2], pk_interp[0]);
  }

  else if (model=="dispersion_Lorentz") {
    if (parameter.size()!=3)
      ErrorCBL("the "+model+" model has 3 parameters, while in the parameter vector in input has "+conv(parameter.size(), par::fINT)+" parameters!", "Pkmu", "ModelFunction_TwoPointCorrelation.cpp");
    pkmu = Pkmu_dispersion(kp, mup, "Lorentzian", parameter[0], parameter[1], parameter[2], pk_interp[0]);
  }
  
  else if (model=="Scoccimarro_Gauss") {
    if (parameter.size()!=3)
      ErrorCBL("the "+model+" model has 3 parameters, while in the parameter vector in input has "+conv(parameter.size(), par::fINT)+" parameters!", "Pkmu", "ModelFunction_TwoPointCorrelation.cpp");
    pkmu = Pkmu_Scoccimarro(kp, mup, "Gaussian", parameter[0], parameter[1], parameter[2], pk_interp[0], pk_interp[1], pk_interp[2]);
  }

  else if (model=="Scoccimarro_Lorentz") {
    if (parameter.size()!=3)
      ErrorCBL("the "+model+" model has 3 parameters, while in the parameter vector in input has "+conv(parameter.size(), par::fINT)+" parameters!", "Pkmu", "ModelFunction_TwoPointCorrelation.cpp");
    pkmu = Pkmu_Scoccimarro(kp, mup, "Lorentzian", parameter[0], parameter[1], parameter[2], pk_interp[0], pk_interp[1], pk_interp[2]);
  }
  
  else if (model=="Scoccimarro_Pezzotta_Gauss") {
    if (parameter.size()!=5)
      ErrorCBL("the "+model+" model has 5 parameters, while in the parameter vector in input has "+conv(parameter.size(), par::fINT)+" parameters!", "Pkmu", "ModelFunction_TwoPointCorrelation.cpp");
    pkmu = Pkmu_Scoccimarro_fitPezzotta(kp, mup, "Gaussian", parameter[0], parameter[1], parameter[2], parameter[3], parameter[4], pk_interp[0], pk_interp[1]);
  }

  else if (model=="Scoccimarro_Pezzotta_Lorentz") {
    if (parameter.size()!=5)
      ErrorCBL("the "+model+" model has 5 parameters, while in the parameter vector in input has "+conv(parameter.size(), par::fINT)+" parameters!", "Pkmu", "ModelFunction_TwoPointCorrelation.cpp");
    pkmu = Pkmu_Scoccimarro_fitPezzotta(kp, mup, "Lorentzian", parameter[0], parameter[1], parameter[2], parameter[3], parameter[4], pk_interp[0], pk_interp[1]);
  }

  else if (model=="Scoccimarro_Bel_Gauss") {
    if (parameter.size()!=8)
      ErrorCBL("the "+model+" model has 8 parameters, while in the parameter vector in input has "+conv(parameter.size(), par::fINT)+" parameters!", "Pkmu", "ModelFunction_TwoPointCorrelation.cpp");
    pkmu = Pkmu_Scoccimarro_fitBel(kp, mup, "Gaussian", parameter[0], parameter[1], parameter[2], parameter[3], parameter[4], parameter[5], parameter[6], parameter[7], pk_interp[0], pk_interp[1]);
  }

  else if (model=="Scoccimarro_Bel_Lorentz") {
    if (parameter.size()!=8)
      ErrorCBL("the "+model+" model has 8 parameters, while in the parameter vector in input has "+conv(parameter.size(), par::fINT)+" parameters!", "Pkmu", "ModelFunction_TwoPointCorrelation.cpp");
    pkmu = Pkmu_Scoccimarro_fitBel(kp, mup, "Lorentzian", parameter[0], parameter[1], parameter[2], parameter[3], parameter[4], parameter[5], parameter[6], parameter[7], pk_interp[0], pk_interp[1]);
  }

  else if (model=="TNS_Gauss") {
    if (parameter.size()!=3)
      ErrorCBL("the "+model+" model has 3 parameters, while in the parameter vector in input has "+conv(parameter.size(), par::fINT)+" parameters!", "Pkmu", "ModelFunction_TwoPointCorrelation.cpp");
    pkmu = Pkmu_TNS(kp, mup, "Gaussian", parameter[0], parameter[1], parameter[2], pk_interp[0], pk_interp[1], pk_interp[2], pk_interp[3], pk_interp[4], pk_interp[5], pk_interp[6], pk_interp[7], pk_interp[8], pk_interp[9], pk_interp[10], pk_interp[11], pk_interp[12], pk_interp[13], pk_interp[14], pk_interp[15], pk_interp[16]);
  }

  else if (model=="TNS_Lorentz") {
    if (parameter.size()!=3)
      ErrorCBL("the "+model+" model has 3 parameters, while in the parameter vector in input has "+conv(parameter.size(), par::fINT)+" parameters!", "Pkmu", "ModelFunction_TwoPointCorrelation.cpp");
    pkmu = Pkmu_TNS(kp, mup, "Lorentzian", parameter[0], parameter[1], parameter[2], pk_interp[0], pk_interp[1], pk_interp[2], pk_interp[3], pk_interp[4], pk_interp[5], pk_interp[6], pk_interp[7], pk_interp[8], pk_interp[9], pk_interp[10], pk_interp[11], pk_interp[12], pk_interp[13], pk_interp[14], pk_interp[15], pk_interp[16]);
  }

  else if (model=="eTNS_Gauss") {
    if (parameter.size()!=5)
      ErrorCBL("the "+model+" model has 5 parameters, while in the parameter vector in input has "+conv(parameter.size(), par::fINT)+" parameters!", "Pkmu", "ModelFunction_TwoPointCorrelation.cpp");
    pkmu = Pkmu_eTNS(kp, mup, "Gaussian", parameter[0], parameter[1], parameter[2], parameter[3], parameter[4], pk_interp[0], pk_interp[1], pk_interp[2], pk_interp[3], pk_interp[4], pk_interp[5], pk_interp[6], pk_interp[7], pk_interp[8], pk_interp[9], pk_interp[10], pk_interp[11], pk_interp[12], pk_interp[13], pk_interp[14], pk_interp[15], pk_interp[16], pk_interp[17], pk_interp[18], pk_interp[19], pk_interp[20], pk_interp[21], pk_interp[22], pk_interp[23], pk_interp[24]);
  }

  else if (model=="eTNS_Lorentz") {
    if (parameter.size()!=5)
      ErrorCBL("the "+model+" model has 5 parameters, while in the parameter vector in input has "+conv(parameter.size(), par::fINT)+" parameters!", "Pkmu", "ModelFunction_TwoPointCorrelation.cpp");
    pkmu = Pkmu_eTNS(kp, mup, "Lorentzian", parameter[0], parameter[1], parameter[2], parameter[3], parameter[4], pk_interp[0], pk_interp[1], pk_interp[2], pk_interp[3], pk_interp[4], pk_interp[5], pk_interp[6], pk_interp[7], pk_interp[8], pk_interp[9], pk_interp[10], pk_interp[11], pk_interp[12], pk_interp[13], pk_interp[14], pk_interp[15], pk_interp[16], pk_interp[17], pk_interp[18], pk_interp[19], pk_interp[20], pk_interp[21], pk_interp[22], pk_interp[23], pk_interp[24]);
  }

  else
    ErrorCBL("the chosen model ("+model+") is not currently implemented!", "Pkmu", "ModelFunction_TwoPointCorrelation.cpp");

  return pkmu/(alpha_perp*alpha_perp*alpha_par);
}


// ============================================================================================


double cbl::modelling::powspec::Pk_l (const double kk, const int l, const std::string model, const std::vector<double> parameter, const std::vector<std::shared_ptr<glob::FuncGrid>> pk_interp, const double prec, const double alpha_perp, const double alpha_par)
{
  auto integrand = [&] (const double mu)
    {
      return modelling::powspec::Pkmu(kk, mu, model, parameter, pk_interp, alpha_perp, alpha_par)*legendre_polynomial(mu, l);
    };

  return 0.5*(2*l+1)*wrapper::gsl::GSL_integrate_cquad(integrand, -1., 1., prec);
}

// ============================================================================================


std::vector<double> cbl::modelling::powspec::Pk_l (const std::vector<double> kk, const int l, const std::string model, const std::vector<double> parameter, const std::vector<std::shared_ptr<glob::FuncGrid>> pk_interp, const double prec, const double alpha_perp, const double alpha_par)
{
  vector<double> Pkl(kk.size(), 0);

  for (size_t i=0; i<kk.size(); i++) {

    auto integrand = [&] (const double mu)
      {
	return modelling::powspec::Pkmu(kk[i], mu, model, parameter, pk_interp, alpha_perp, alpha_par)*legendre_polynomial(mu, l);
      };

    Pkl[i] = 0.5*(2*l+1)*wrapper::gsl::GSL_integrate_cquad(integrand, -1., 1., prec);
  }

  return Pkl;
}


// ============================================================================================


std::vector<std::vector<double>> cbl::modelling::powspec::damped_Pk_terms (const std::vector<double> kk, const double linear_growth_rate, const double SigmaS, const std::shared_ptr<cbl::glob::FuncGrid> PkDM)
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


std::vector<double> cbl::modelling::powspec::Pl_Kaiser (const std::vector<double> kk, const std::shared_ptr<void> inputs, std::vector<double> &parameter)
{ 
 // structure contaning the required input data
  shared_ptr<STR_data_model> pp = static_pointer_cast<STR_data_model>(inputs);

  // redefine the cosmology
  auto cosmo = pp->cosmology->clone();
  auto cosmo_fid = pp->cosmology_fid->clone();
  
  // set the cosmological parameters 
  for (size_t i=0; i<pp->Cpar.size(); ++i)
    cosmo->set_parameter(pp->Cpar[i], parameter[i]);
    
  cosmology::PkXi PX(cosmo);
 
  double bias,f;
  vector<double> Pk_l;
  int kk_stop = kk.size()/(pp->redshifts.size()*pp->ll.size());
  vector<double> kk_obs(kk.begin(),kk.begin()+kk_stop);
  vector<double> Pk(kk_stop,0), Pl(kk_stop,0), qq(2,0);
  const string mod = "dispersion_Gauss"; 
  const double prec = 1.e-5;

  for (size_t i=0; i<pp->redshifts.size(); ++i)
    {    
      bias = pp->bias[i];
      f = cosmo->linear_growth_rate(pp->redshifts[i],1);
      vector<double> par={f,bias,0.0}; 
      Pk = PX.Pk_matter(kk_obs, pp->method_Pk, pp->NL, pp->redshifts[i], pp->store_output, pp->output_root, pp->norm, pp->k_min, pp->k_max, pp->prec,pp->file_par, pp->unit);
      qq = {cosmo->D_A(pp->redshifts[i],0)/cosmo_fid->D_A(pp->redshifts[i], 0), cosmo_fid->Hubble(pp->redshifts[i])/cosmo->Hubble(pp->redshifts[i])};

      for (size_t j=0; j<pp->ll.size();++j)
      {
        Pl = cbl::modelling::powspec::Pk_l(kk_obs,pp->ll[j],mod,par,{make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(kk_obs, Pk, "Spline"))},prec,qq[0],qq[1]) ;

        for (size_t h=0; h<Pl.size();++h)
          Pk_l.push_back(Pl[h]);
      }
    }
      
  return Pk_l;

}


