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
 *  @file Cosmology/Lib/MassFunction.cpp
 *
 *  @brief Methods of the class Cosmology used to model the mass
 *  function
 *
 *  This file contains the implementation of the methods of the class
 *  Cosmology used to model the mass function of dark matter haloes
 *
 *  @author Federico Marulli 
 *
 *  @author federico.marulli3@unibo.it
 */

#include "Cosmology.h"

using namespace std;

using namespace cbl;


// =====================================================================================


vector<double> cbl::cosmology::Cosmology::mass_function (const vector<double> Mass, const vector<double> Sigma, const vector<double> Dln_Sigma, const double redshift, const std::string model_MF, const bool store_output, const std::string output_root, const double Delta, const std::string interpType, const int norm, const double k_min, const double k_max, const double prec, const std::string method_SS, const std::string input_file, const bool is_parameter_file) 
{
  double fact = (m_unit) ? 1 : m_hh;
  double D_N = DN(redshift);

  vector<double> MASS(Mass.size(), 0), MF(Mass.size(), 0);
  for (size_t i=0; i<Mass.size(); i++) {
      MASS[i] = Mass[i]*fact;
      MF[i] = m_MF_generator(MASS[i], Sigma[i], Dln_Sigma[i], redshift, D_N, model_MF, Delta)*pow(fact, 4.);
  }

  if (m_fNL!=0)
      for (size_t i=0; i<Mass.size(); i++)
          MF[i] *= MF_correction(MASS[i], redshift, method_SS, store_output, output_root, interpType, norm, k_min, k_max, prec, input_file, is_parameter_file);

  return MF;
}


// =====================================================================================


double cbl::cosmology::Cosmology::m_MF_generator (const double Mass, const double Sigma, const double Dln_Sigma, const double redshift, const double D_N, const std::string model_MF, const double Delta, const bool default_delta, const double delta_t) 
{ 
  const double deltacz = (default_delta) ? deltac(redshift) : fabs(delta_t*D_N);
  const double sigmaz = Sigma*D_N;
  
  const double RHO = rho_m(0., true); 
  
  if (model_MF=="PS") // Press & Schechter
    return sqrt(2./par::pi)*RHO/(Mass*Mass)*deltacz/sigmaz*fabs(Dln_Sigma)*exp(-(deltacz*deltacz)*0.5/(sigmaz*sigmaz));
   
  else if (model_MF=="ST") { // Sheth & Tormen
    const double aa = 0.707;
    const double qq = 0.3;
    const double AA = 0.3222;
    return AA*sqrt(2.*aa/par::pi)*RHO/(Mass*Mass)*deltacz/sigmaz*(1+pow(sigmaz/(sqrt(aa)*deltacz),2.*qq))*fabs(Dln_Sigma)*exp(-(aa*deltacz*deltacz)*0.5/(sigmaz*sigmaz));
  }
  
  else if (model_MF=="Jenkins") // Jenkins et al. (2001)
    return 0.315*exp(-pow(fabs(log(1./sigmaz)+0.61),3.8))*RHO/(Mass*Mass)*fabs(Dln_Sigma);
  
  else if (model_MF=="Warren") { // Warren et al. (2006) 
    const double AA = 0.7234;
    const double aa = 1.625;
    const double bb = 0.2538;
    const double cc = 1.1982;
    return AA*(pow(sigmaz,-aa)+bb)*exp(-cc/(sigmaz*sigmaz))*RHO/(Mass*Mass)*fabs(Dln_Sigma);
  }
  
  else if (model_MF=="Reed") { // Reed et al. (2007) 
    const double cc = 1.08;
    const double aa = 0.764/cc;
    const double pp = 0.3;
    const double AA = 0.3222;
    const double Fatt = AA*sqrt(2.*aa/par::pi);
    const double n_eff = -6.*Dln_Sigma-3.;
    const double G1 = exp(-pow(log(1./sigmaz)-0.4,2)/0.72);
    const double G2 = exp(-pow(log(1./sigmaz)-0.75,2)/0.08);
    return Fatt*RHO/(Mass*Mass)*deltacz/sigmaz*(1+pow(sigmaz*sigmaz/(aa*deltacz*deltacz),pp)+0.6*G1+0.4*G2)*fabs(Dln_Sigma)*exp(-(cc*aa*deltacz*deltacz)*0.5/(sigmaz*sigmaz)-(0.03/pow(n_eff+3.,2)*pow(deltacz/sigmaz,0.6)));
  }

  else if (model_MF=="Pan") { // Pan (2007)
    const double alpha = 0.435; // check!!!
    return 4.*alpha/sqrt(2.*par::pi)*RHO/(Mass*Mass)*deltacz/pow(sigmaz,2.*alpha)*fabs(Dln_Sigma)*exp(-(deltacz*deltacz)*0.5/pow(sigmaz,4.*alpha));
  }

  else if (model_MF=="ShenH") { // halo MF, Shen et al. (2006) // check!!!
    const double alpha = -0.55;
    const double beta = -0.56;
    const double ni = pow(deltacz/sigmaz,2);
    const double ni_fni = sqrt(ni/(2.*par::pi))*exp(-ni*pow(1.-beta*pow(ni,alpha),2)*0.5)*(1.-beta/pow(ni,-alpha)*(1+alpha+(alpha*(alpha+1.))*0.5));
    return ni_fni*RHO/(Mass*Mass)*fabs(Dln_Sigma)*2;
  }
  
  else if (model_MF=="ShenF") { // filament MF, Shen et al. (2006) // check!!!
    const double alpha = -0.28;
    const double beta = -0.012;
    const double ni = pow(deltacz/sigmaz,2);
    const double ni_fni = sqrt(ni/(2.*par::pi))*exp(-ni*pow(1.-beta*pow(ni,alpha),2)*0.5)*(1.-beta/pow(ni,-alpha)*(1+alpha+(alpha*(alpha+1.))*0.5));
    return ni_fni*RHO/(Mass*Mass)*fabs(Dln_Sigma)*2;
  }

  else if (model_MF=="ShenS") { // sheet MF, Shen et al. (2006) // check!!!
    const double alpha = -0.61;
    const double beta = 0.45;
    const double ni = pow(deltacz/sigmaz,2);
    const double ni_fni = sqrt(ni/(2.*par::pi))*exp(-ni*pow(1.-beta*pow(ni,alpha),2)*0.5)*(1.-beta/pow(ni,-alpha)*(1+alpha+(alpha*(alpha+1.))*0.5));
    return fabs(ni_fni*RHO/(Mass*Mass)*fabs(Dln_Sigma)*2); // check!!!
  }

  else if (model_MF=="Tinker") { // Tinker et al. (2008)

    //if (redshift>2) WarningMsgCBL("the Tinker mass function has been tested for z<~2!", "m_MF_generator", "MassFunction.cpp");

    double A0, a0, b0, c0;

    if      (Delta==200)  {A0 = 0.186; a0 = 1.47; b0 = 2.57; c0 = 1.19;}
    else if (Delta==300)  {A0 = 0.200; a0 = 1.52; b0 = 2.25; c0 = 1.27;}
    else if (Delta==400)  {A0 = 0.212; a0 = 1.56; b0 = 2.05; c0 = 1.34;}
    else if (Delta==600)  {A0 = 0.218; a0 = 1.61; b0 = 1.87; c0 = 1.45;}
    else if (Delta==800)  {A0 = 0.248; a0 = 1.87; b0 = 1.59; c0 = 1.58;}
    else if (Delta==1200) {A0 = 0.255; a0 = 2.13; b0 = 1.51; c0 = 1.80;}
    else if (Delta==1600) {A0 = 0.260; a0 = 2.30; b0 = 1.46; c0 = 1.97;}
    else if (Delta==2400) {A0 = 0.260; a0 = 2.53; b0 = 1.44; c0 = 2.24;}
    else if (Delta==3200) {A0 = 0.260; a0 = 2.66; b0 = 1.41; c0 = 2.44;}
    else {
      A0 = (Delta<1600) ? 0.1*log10(Delta)-0.05 : 0.26;
      a0 = 1.43+pow(max(log10(Delta)-2.3,0.),1.5); // check!!!
      b0 = 1.+pow(log10(Delta)-1.6,-1.5);
      c0 = 1.2+pow(max(log10(Delta)-2.35,0.),1.6); // check!!!
    }

    const double alpha = pow(10.,-pow(max(0.75/log10(Delta/75.),0.),1.2)); // check!!!
    const double AA = A0*pow(1.+redshift,-0.14);
    const double aa = a0*pow(1.+redshift,-0.06);
    const double bb = b0*pow(1.+redshift,-alpha);
    const double cc = c0;

    return AA*(pow(sigmaz/bb,-aa)+1.)*exp(-cc/(sigmaz*sigmaz))*RHO/(Mass*Mass)*fabs(Dln_Sigma);
  }

  else if (model_MF=="Tinker08_Interp") { // Tinker et al. (2008)

    //if (redshift>2) WarningMsgCBL("the Tinker mass function has been tested for z<~2!", "m_MF_generator", "MassFunction.cpp");

    vector<double> _Delta = {200, 300, 400, 600, 800, 1200, 1600, 2400, 3200};
    vector<double> _A0 = {1.858659e-01, 1.995973e-01, 2.115659e-01, 2.184113e-01, 2.480968e-01, 2.546053e-01, 2.600000e-01, 2.600000e-01, 2.600000e-01};
    vector<double> _a0 = {1.466904, 1.521782, 1.559186, 1.614585, 1.869936, 2.128056, 2.301275, 2.529241, 2.661983};
    vector<double> _b0 = {2.571104, 2.254217, 2.048674, 1.869559, 1.588649, 1.507134, 1.464374, 1.436827, 1.405210};
    vector<double> _c0 = {1.193958, 1.270316, 1.335191, 1.446266, 1.581345, 1.795050, 1.965613, 2.237466, 2.439729};

    cbl::glob::FuncGrid interp_A0(_Delta, _A0, "Spline");
    cbl::glob::FuncGrid interp_a0(_Delta, _a0, "Spline");
    cbl::glob::FuncGrid interp_b0(_Delta, _b0, "Spline");
    cbl::glob::FuncGrid interp_c0(_Delta, _c0, "Spline");

    const double alpha = pow(10.,-pow(max(0.75/log10(Delta/75.),0.),1.2)); // check!!!
    const double AA = interp_A0(Delta)*pow(1.+redshift,-0.14);
    const double aa = interp_a0(Delta)*pow(1.+redshift,-0.06);
    const double bb = interp_b0(Delta)*pow(1.+redshift,-alpha);
    const double cc = interp_c0(Delta);

    return AA*(pow(sigmaz/bb,-aa)+1.)*exp(-cc/(sigmaz*sigmaz))*RHO/(Mass*Mass)*fabs(Dln_Sigma);
  }


  else if (model_MF=="Crocce") { // Crocce et al. (2010)
    const double AA = 0.58*pow(1.+redshift,-0.13);
    const double aa = 1.37*pow(1.+redshift,-0.15);
    const double bb = 0.3*pow(1.+redshift,-0.084);
    const double cc = 1.036*pow(1.+redshift,-0.024);
    return AA*(pow(sigmaz,-aa)+bb)*exp(-cc*pow(sigmaz,-2))*RHO/(Mass*Mass)*fabs(Dln_Sigma);
  }
  
  else if (model_MF=="Angulo_FOF") // FoF MF, Angulo et al. (2012)
    return 0.201*(pow(2.08/sigmaz,1.7)+1)*exp(-1.172/(sigmaz*sigmaz))*RHO/(Mass*Mass)*fabs(Dln_Sigma);
 
  else if (model_MF=="Angulo_Sub") // SUBFIND MF, Angulo et al. (2012)
    return 0.265*(pow(1.675/sigmaz,1.9)+1)*exp(-1.4/(sigmaz*sigmaz))*RHO/(Mass*Mass)*fabs(Dln_Sigma);

  else if (model_MF=="Bhattacharya") { // Bhattacharya et al. (2011) // check!!!
    const double aa = 0.788/pow(1+redshift,0.01);
    const double qq = 1.795;
    const double AA = 0.333/pow(1+redshift,0.11);
    const double pp = 0.807;
    const double nu = deltacz/sigmaz;
    //if (redshift>=2) aa = 0.799/pow(1+redshift,0.024); 
    return RHO/(Mass*Mass)*fabs(Dln_Sigma)*AA*sqrt(2/par::pi)*exp(-nu*nu*0.5*aa)*(1+pow(1/(nu*nu*aa),pp))*pow(nu*sqrt(aa),qq);
  }

  else if (model_MF=="Courtin") { // Courtin et al. (2010) 
    const double aa = 0.695;
    const double qq = 0.1;
    const double AA = 0.348; 
    return AA*sqrt(2.*aa/par::pi)*RHO/(Mass*Mass)*deltacz/sigmaz*(1+pow(sigmaz/(sqrt(aa)*deltacz),2.*qq))*fabs(Dln_Sigma)*exp(-(aa*deltacz*deltacz)*0.5/(sigmaz*sigmaz));
  }

  else if (model_MF=="Manera") { // Manera et al. (2010)
    
    double aa, qq;
    if (fabs(redshift-0.)<1.e-30) {
      aa = 0.709;
      qq = 0.248;
    }
    else if (fabs(redshift-0.5)<1.e-30) {
      aa = 0.724;
      qq = 0.241;
    }
    else 
      ErrorCBL("the Manera et al. mass function has been tested only for z=0 and 0.5!", "m_MF_generator", "MassFunction.cpp");
    
    return 0.3222*sqrt(2.*aa/par::pi)*RHO/(Mass*Mass)*deltacz/sigmaz*(1+pow(sigmaz/(sqrt(aa)*deltacz),2.*qq))*fabs(Dln_Sigma)*exp(-(aa*deltacz*deltacz)*0.5/(sigmaz*sigmaz));
  }

  else if (model_MF=="Watson_FOF") { // FoF MF, Watson et al. (2012) 
    const double a = 0.282;
    const double b = 1.406;
    const double c = 2.163;
    const double d = 1.210;
    return a*(pow((b/sigmaz),c)+1)*exp(-d/(sigmaz*sigmaz))*RHO/(Mass*Mass)*fabs(Dln_Sigma);
  }

  else if (model_MF=="Watson_SOH") { // Spherical Overdensity halo MF, Watson et al. (2012) 
    const double Omega = OmegaM(redshift);
    const double C = exp(0.023*(Delta/178-1));
    const double d = -0.456*Omega-0.139;
    const double p = 0.072;
    const double q = 2.13;
    const double Gamma = C*pow((Delta/178),d)*exp(p*(1-Delta/178)/pow(sigmaz,q));
    double A = 0.194;
    double alpha = 2.267;
    double beta = 1.805;
    double gamma = 1.287;
    if (redshift>0 && redshift<6) {
      A = Omega*(1.097*pow((1+redshift),-3.216)+0.074);
      alpha = Omega*(3.136*pow((1+redshift),-3.056)+2.349);
      beta = Omega*(5.907*pow((1+redshift),-3.599)+2.344);
      gamma = 1.318;
    }
    if (redshift>=6) {
      A = 0.563;
      alpha = 0.874;
      beta = 3.810;
      gamma = 1.453;
    }
    return RHO/(Mass*Mass)*fabs(Dln_Sigma)*Gamma*A*(pow((beta/sigmaz),alpha)+1.)*exp(-gamma/(sigmaz*sigmaz));
  }

  else if (model_MF=="Peacock") { // Peacock et al. (2007)
    const double aa = 1.529;
    const double bb = 0.704;
    const double cc = 0.412;
    const double nu = deltacz/sigmaz;
    return nu*exp(-cc*nu*nu)/pow(1+aa*pow(nu,bb),2)*(bb*aa*pow(nu,bb-1)+2*cc*nu*(1+aa*pow(nu,bb)))*RHO/(Mass*Mass)*fabs(Dln_Sigma);
  }

  else if (model_MF=="Despali_Z0") { // Despali et al. 2016
    const double aa = 0.794; // 0.794 ± 0.005, z=0 fit
    const double pp = 0.247; // 0.247 ± 0.009, z=0 fit
    const double AA = 0.333; // 0.333 ± 0.001, z=0 fit
    const double nu = pow(deltacz/sigmaz, 2);
    const double nup = aa*nu;
    return 2.*AA*(1.+pow(nup,-pp))*sqrt(nup/(2.*par::pi))*exp(-0.5*nup)*RHO/(Mass*Mass)*fabs(Dln_Sigma);
  }
    
  else if (model_MF=="Despali_AllZ") { // Despali et al. 2016
    const double aa = 0.7663; // 0.7663 ± 0.0013, all z fit
    const double pp = 0.2579; // 0.2579 ± 0.0026, all z fit
    const double AA = 0.3298; // 0.3298 ± 0.0003, all z fit
    const double nu = pow(deltacz/sigmaz, 2);
    const double nup = aa*nu;
    return 2.*AA*(1.+pow(nup,-pp))*sqrt(nup/(2.*par::pi))*exp(-0.5*nup)*RHO/(Mass*Mass)*fabs(Dln_Sigma);
  }
  
  else if (model_MF=="Despali_AllZAllCosmo") { // Despali et al. 2016
    const double aa = 0.7689; // 0.7689 ± 0.0011, all z fit
    const double pp = 0.2536; // 0.2536 ± 0.0026, all z fit
    const double AA = 0.3295; // 0.3295 ± 0.0003, all z fit
    const double nu = pow(deltacz/sigmaz, 2);
    const double nup = aa*nu;
    return 2.*AA*(1.+pow(nup,-pp))*sqrt(nup/(2.*par::pi))*exp(-0.5*nup)*RHO/(Mass*Mass)*fabs(Dln_Sigma);
  }
  
  else if (model_MF=="Despali_HighM") { // Despali et al. 2016
    const double aa = 0.8199; // 0.8199 ± 0.0010, M>3e13
    const double pp = 0.0; // 0, M>3e13
    const double AA = 0.3141; // 0.3141 ± 0.0006, M>3e13
    const double nu = pow(deltacz/sigmaz, 2);
    const double nup = aa*nu;
    return 2.*AA*(1.+pow(nup,-pp))*sqrt(nup/(2.*par::pi))*exp(-0.5*nup)*RHO/(Mass*Mass)*fabs(Dln_Sigma);
  }
  
  else return ErrorCBL("model_MF not allowed!", "m_MF_generator", "MassFunction.cpp");  
}
