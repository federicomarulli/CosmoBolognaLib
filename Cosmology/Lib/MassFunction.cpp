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
 *  @author federico.marulli3@unbo.it
 */

#include "Cosmology.h"
using namespace cosmobl;


// =====================================================================================


double cosmobl::Cosmology::mass_function (double &Mass, double &redshift, string &author, string &method_SS, string output_root, double Delta, string interpType, int Num, double stepsize, int norm, double k_min, double k_max, bool GSL, double prec, string file_par) 
{
  double fact = (m_unit) ? 1 : m_hh;
  double MASS = Mass*fact;

  double zero = 0.;
  double SSS = SSM_norm(MASS, method_SS, zero, output_root, k_max, file_par);
  double Sigma = sqrt(SSS);
  double Dln_Sigma = dnSM(1, MASS, method_SS, zero, output_root, interpType, Num, stepsize, k_max, file_par)*(MASS/(2.*SSS));

  double MF = MF_generator(MASS, Sigma, Dln_Sigma, redshift, author, Delta)*pow(fact,4.);

  if (m_fNL!=0) MF *= MF_correction(MASS, redshift, author, output_root, interpType, Num, stepsize, norm, k_min, k_max, GSL, prec, file_par);

  return MF;
}


// =====================================================================================


double cosmobl::Cosmology::mass_function_fast (double &Mass, double &redshift, string &author, string &method_SS, string output_root, double Delta, string interpType, int Num, double stepsize, int norm, double k_min, double k_max, bool GSL, double prec, string file_par) 
{
  double fact = (m_unit) ? 1 : m_hh;
  double MASS = Mass*fact;
  

  // ---------- read the grid file ---------- 
  
  double zero = 0.;
  string file_grid = create_grid_sigmaM (method_SS, zero, output_root, interpType, Num, stepsize, k_max);
  ifstream fin (file_grid.c_str()); checkIO (file_grid,1); 

  double MMass, Sigma, Dln_Sigma;
  vector<double> mass, sigma, dln_sigma;
  while (fin >>MMass>>Sigma>>Dln_Sigma) {
    mass.push_back(MMass);
    sigma.push_back(Sigma);
    dln_sigma.push_back(Dln_Sigma);
  } 


  // ---------- compute the MF ---------- 
   
  double err = -1., sig = -1., dlsig = -1.;
  interpolation_extrapolation(MASS, mass, sigma, "Rat", 4, &sig, &err);
  interpolation_extrapolation(MASS, mass, dln_sigma, "Rat", 4, &dlsig, &err);

  double MF = MF_generator(MASS, sig, dlsig, redshift, author, Delta)*pow(fact,4.);
  if (std::isnan(MF)) { string Err = "Error in cosmobl::Cosmology::mass_function_fast of MassFunction.cpp: MF = " + conv(MF,par::fDP3) + "!"; ErrorMsg(Err); }

  if (m_fNL!=0) MF *= MF_correction(MASS, redshift, method_SS, output_root, interpType, Num, stepsize, norm, k_min, k_max, GSL, prec, file_par);

  return MF;  
}


// =====================================================================================


double cosmobl::Cosmology::mass_function (double &Mass, double &Sigma, double &Dln_Sigma, double &redshift, string &author_MF, string output_root, double Delta, string interpType, int Num, double stepsize, int norm, double k_min, double k_max, bool GSL, double prec, string method_SS, string file_par) 
{
  double fact = (m_unit) ? 1 : m_hh;
  double MASS = Mass*fact;

  double MF = MF_generator(MASS, Sigma, Dln_Sigma, redshift, author_MF, Delta)*pow(fact,4.);

  if (m_fNL!=0) MF *= MF_correction(MASS, redshift, method_SS, output_root, interpType, Num, stepsize, norm, k_min, k_max, GSL, prec, file_par);

  return MF;
}


// =====================================================================================


double cosmobl::Cosmology::MF_generator (double &Mass, double &Sigma, double &Dln_Sigma, double &redshift, string &author, double Delta) 
{ 
  double Z0 = 0.;
  double deltacz = deltac(redshift);
  double sigmaz = Sigma*DD(redshift)/DD(Z0);

  double RHO = Rho(m_Omega_matter,m_Omega_neutrinos); 

  double dndm = -1.;
  
  if (author=="PS") // Press-Schechter
    dndm = sqrt(2./par::pi)*RHO/(Mass*Mass)*deltacz/sigmaz*fabs(Dln_Sigma)*exp(-(deltacz*deltacz)*0.5/(sigmaz*sigmaz));
   
  if (author=="ST") { // Sheth-Tormen
    double aa = 0.707;
    double qq = 0.3;
    double AA = 0.3222;
    dndm = AA*sqrt(2.*aa/par::pi)*RHO/(Mass*Mass)*deltacz/sigmaz*(1+pow(sigmaz/(sqrt(aa)*deltacz),2.*qq))*fabs(Dln_Sigma)*exp(-(aa*deltacz*deltacz)*0.5/(sigmaz*sigmaz));
  }
  
  if (author=="Jenkins") // Jenkins et al. (2001)
    dndm = 0.315*exp(-pow(fabs(log(1./sigmaz)+0.61),3.8))*RHO/(Mass*Mass)*fabs(Dln_Sigma);
  
  if (author=="Warren") { // Warren et al. (2006) 
    double AA = 0.7234;
    double aa = 1.625;
    double bb = 0.2538;
    double cc = 1.1982;
    dndm = AA*(pow(sigmaz,-aa)+bb)*exp(-cc/(sigmaz*sigmaz))*RHO/(Mass*Mass)*fabs(Dln_Sigma);
  }
  
  if (author=="Reed") { // Reed et al. (2007) 
    double aa = 0.707;
    double cc = 1.08;
    double pp = 0.3;
    double AA = 0.3222;
    double Fatt = AA*sqrt(2.*aa/par::pi);
    double n_eff = -6.*Dln_Sigma-3.;
    double G1 = exp(-pow(log(1./sigmaz)-0.4,2)/0.72);
    double G2 = exp(-pow(log(1./sigmaz)-0.75,2)/0.08);
    dndm = Fatt*RHO/(Mass*Mass)*deltacz/sigmaz*(1+pow(sigmaz*sigmaz/(aa*deltacz*deltacz),pp)+0.6*G1+0.4*G2)*fabs(Dln_Sigma)*exp(-(cc*aa*deltacz*deltacz)*0.5/(sigmaz*sigmaz)-(0.03/pow(n_eff+3.,2)*pow(deltacz/sigmaz,0.6)));
  }

  if (author=="Pan") { // Pan (2007)
    double alpha = 0.435; // check!!!
    dndm = 4.*alpha/sqrt(2.*par::pi)*RHO/(Mass*Mass)*deltacz/pow(sigmaz,2.*alpha)*fabs(Dln_Sigma)*exp(-(deltacz*deltacz)*0.5/pow(sigmaz,4.*alpha));
  }

  if (author=="ShenH") { // halo MF by Shen et al. (2006) // check!!!
    double alpha = -0.55;
    double beta = -0.56;
    double ni = pow(deltacz/sigmaz,2);
    double ni_fni = sqrt(ni/(2.*par::pi))*exp(-ni*pow(1.-beta*pow(ni,alpha),2)*0.5)*(1.-beta/pow(ni,-alpha)*(1+alpha+(alpha*(alpha+1.))*0.5));
    dndm = ni_fni*RHO/(Mass*Mass)*fabs(Dln_Sigma)*2;
  }
  
  if (author=="ShenF") { // filaments MF by Shen et al. (2006) // check!!!
    double alpha = -0.28;
    double beta = -0.012;
    double ni = pow(deltacz/sigmaz,2);
    double ni_fni = sqrt(ni/(2.*par::pi))*exp(-ni*pow(1.-beta*pow(ni,alpha),2)*0.5)*(1.-beta/pow(ni,-alpha)*(1+alpha+(alpha*(alpha+1.))*0.5));
    dndm = ni_fni*RHO/(Mass*Mass)*fabs(Dln_Sigma)*2;
  }

  if (author=="ShenS") { // sheets MF by Shen et al. (2006) // check!!!
    double alpha = -0.61;
    double beta = 0.45;
    double ni = pow(deltacz/sigmaz,2);
    double ni_fni = sqrt(ni/(2.*par::pi))*exp(-ni*pow(1.-beta*pow(ni,alpha),2)*0.5)*(1.-beta/pow(ni,-alpha)*(1+alpha+(alpha*(alpha+1.))*0.5));
    dndm = fabs(ni_fni*RHO/(Mass*Mass)*fabs(Dln_Sigma)*2); // check!!!
  }

  if (author=="Tinker") { // Tinker et al. (2008)

    if (redshift>2) WarningMsg("Attention: the Tinker mass function has been tested for z<~2!");

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
      a0 = 1.43+pow(log10(Delta)-2.3,1.5);
      b0 = 1.+pow(log10(Delta)-1.6,-1.5);
      c0 = (log10(Delta)>2.35) ? 1.2+pow(log10(Delta)-2.35,1.6) : 1.2; // check!!!
    }

    double alpha = pow(10.,-pow(0.75/log10(Delta/75.),1.2));
    double AA = A0*pow(1.+redshift,-0.14);
    double aa = a0*pow(1.+redshift,-0.06);
    double bb = b0*pow(1.+redshift,-alpha);
    double cc = c0;
   
    dndm = AA*(pow(sigmaz/bb,-aa)+1.)*exp(-cc/(sigmaz*sigmaz))*RHO/(Mass*Mass)*fabs(Dln_Sigma);
  }

  if (author=="Angulo_FOF") { // FOF MF by Angulo et al. (2012)
    double ff = 0.201*pow(2.08/sigmaz,1.7)*exp(-1.172/(sigmaz*sigmaz));
    dndm = RHO/(Mass*Mass)*(-Dln_Sigma)*ff;
  }

  if (author=="Angulo_Sub") { // SUBFIND MF by Angulo et al. (2012)
    double ff = 0.265*pow(1.675/sigmaz,1.9)*exp(-1.4/(sigmaz*sigmaz));
    dndm = RHO/(Mass*Mass)*(-Dln_Sigma)*ff;
  }

  if (dndm<0) ErrorMsg("Error in cosmobl::Cosmology::MF of MassFunction.cpp!");
  return dndm;  
}


// =====================================================================================


double cosmobl::Cosmology::n_haloes (double &Mass_min, double &Mass_max, double &z_min, double &z_max, bool &angle_rad, string &author_MF, string &method_SS, string output_root, string interpType, int Num, double stepsize, double k_max, string file_par)
{

  // ---------- read the grid file ---------- 

  double zero = 0.;
  string file_grid = create_grid_sigmaM (method_SS,zero,output_root,interpType,Num,stepsize,k_max,file_par);
  ifstream fin (file_grid.c_str()); checkIO (file_grid,1); 
  
  double Mass, Sigma, Dln_Sigma;
  vector<double> mass, sigma, dlnsigma;
  while (fin >>Mass>>Sigma>>Dln_Sigma) {
    if (Mass_min<Mass && Mass<Mass_max) {
      mass.push_back(Mass);
      sigma.push_back(Sigma);
      dlnsigma.push_back(Dln_Sigma);
    }
  }


  // ---------- compute the number of haloes ---------- 
  
  int step_z = 100;
  double delta_z = (z_max-z_min)/step_z; 
  double redshift = z_min;
  double N_haloes = 0.;

  for (int Red=0; Red<step_z; Red++) {

    double Int = 0.;
    for (unsigned int k=0; k<mass.size()-1; k++) 
      Int += mass_function(mass[k], sigma[k], dlnsigma[k], redshift, author_MF, output_root)*dV_dZdOmega(redshift, angle_rad)*(mass[k+1]-mass[k]);

    N_haloes += Int*delta_z;

    redshift += delta_z;
  }
  
  return N_haloes;
}


// =====================================================================================


double cosmobl::Cosmology::MhaloMin (int &n_halo, double &Area, bool &angle_rad, double &z_min, double &z_max, double &Mmax, double &lgM1_guess, double &lgM2_guess, string &author_MF, string &method_SS, string output_root, string interpType, int Num, double stepsize, double k_max, string file_par)
{
  cosmobl::classfunc::func_MhaloMin func (m_Omega_matter, m_Omega_baryon, m_Omega_neutrinos, m_massless_neutrinos, m_massive_neutrinos, m_Omega_DE, m_Omega_radiation, m_hh, m_scalar_amp, m_n_spec, m_w0, m_wa, m_fNL, m_type_NG, m_model, m_unit, n_halo, Area, angle_rad, z_min, z_max, Mmax, author_MF, method_SS, output_root, interpType, Num, stepsize, k_max, file_par);
  
  double prec = 0.0001;
  double lgM = zbrent (func, lgM1_guess, lgM2_guess, prec);

  return pow(10.,lgM);
}


// =====================================================================================


double cosmobl::Cosmology::unevolved_mass_function (double &mass_accr)
{
  double aa = -0.8;
  double mm = mass_accr/aa; // mass_accr = m/M
  
  double yn = log(0.3); // log(0.21) check!!!

  return aa*log(-mm)+log(exp(2.*par::pi*pow(mm,3.)))+yn;
}
