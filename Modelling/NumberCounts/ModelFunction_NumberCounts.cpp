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
 *  Modelling/NumberCounts/ModelFunction_NumberCounts.cpp
 *
 *  @brief Functions to model the number counts
 *
 *  This file contains the implementation of the functions used to
 *  model the monopole of the number counts
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unibo.it, alfonso.veropalumbo@unibo.it
 */


#include "ModelFunction_NumberCounts.h"

using namespace std;

using namespace cbl;


// ===========================================================================================


double cbl::modelling::numbercounts::Filter_sigmaR (const double kk, const double radius)
{
  return pow(TopHat_WF(kk*radius),2);
}


// ===========================================================================================


double cbl::modelling::numbercounts::Filter_dsigmaR (const double kk, const double radius)
{
  return 2.*cbl::TopHat_WF(kk*radius)*cbl::TopHat_WF_D1(kk*radius)*kk;
}


// ===========================================================================================


void cbl::modelling::numbercounts::sigmaM_dlnsigmaM (double &sigmaM, double &dlnsigmaM, const double mass, const cbl::glob::FuncGrid interp_Pk, const double kmax, const double rho)
{
  double norm =  1./(2.*pow(par::pi, 2));
  double dRdM_fact = pow(3./(4.*par::pi*rho), 1./3.);

  double RR =  Radius(mass, rho);
  double dRdM =  dRdM_fact*pow(mass, -2./3.)/3.;

  auto integrand_sigmaR = [&] (const double kk)
  {
    return kk*kk*interp_Pk(kk)*Filter_sigmaR(kk, RR);
  };

  sigmaM = norm*cbl::wrapper::gsl::GSL_integrate_cquad(integrand_sigmaR, 1.e-4, kmax, 1.e-5);

  auto integrand_dsigmaR = [&] (const double kk)
  {
    return kk*kk*interp_Pk(kk)*Filter_dsigmaR(kk, RR);
  };

  dlnsigmaM = norm*cbl::wrapper::gsl::GSL_integrate_cquad(integrand_dsigmaR, 1.e-4, kmax, 1.e-5)*dRdM*(mass/(2*sigmaM));
  sigmaM = sqrt(sigmaM);
}


// ===========================================================================================


void cbl::modelling::numbercounts::sigmaM_dlnsigmaM (std::vector<double> &sigmaM, std::vector<double> &dlnsigmaM, const std::vector<double> mass, const std::vector<double> kk, const std::vector<double> Pk, const std::string interpType, const double kmax, const double rho)
{
  double norm =  1./(2.*pow(par::pi, 2));
  double dRdM_fact = pow(3./(4.*par::pi*rho), 1./3.);
  cbl::glob::FuncGrid Pk_interp(kk, Pk, interpType);
  sigmaM.resize(mass.size(), 0);
  dlnsigmaM.resize(mass.size(), 0);

  for (size_t i=0; i<mass.size(); i++) {

    double RR =  Radius(mass[i], rho);
    double dRdM =  dRdM_fact*pow(mass[i], -2./3.)/3.;

    auto integrand_sigmaR = [&] (const double kk)
    {
      return kk*kk*Pk_interp(kk)*Filter_sigmaR(kk, RR);
    };

    sigmaM[i] = norm*cbl::wrapper::gsl::GSL_integrate_cquad(integrand_sigmaR, 1.e-4, kmax, 1.e-5);

    auto integrand_dsigmaR = [&] (const double kk)
    {
      return kk*kk*Pk_interp(kk)*Filter_dsigmaR(kk, RR);
    };

    dlnsigmaM[i] = norm*cbl::wrapper::gsl::GSL_integrate_cquad(integrand_dsigmaR, 1.e-4, kmax, 1.e-5)*dRdM*(mass[i]/(2*sigmaM[i]));
    sigmaM[i] = sqrt(sigmaM[i]);
  }
}


// ===========================================================================================


std::vector<cbl::glob::FuncGrid> cbl::modelling::numbercounts::sigmaM_dlnsigmaM (const std::vector<double> mass, cbl::cosmology::Cosmology cosmology, const std::vector<double> kk, const std::vector<double> Pk, const std::string interpType, const double kmax)
{
  const double rho = cosmology.rho_m(0., true);

  vector<double> sigmaM, dlnsigmaM; 

  cbl::modelling::numbercounts::sigmaM_dlnsigmaM (sigmaM, dlnsigmaM, mass, kk, Pk, interpType, kmax, rho);

  vector<cbl::glob::FuncGrid> interp(2);
  interp[0] = cbl::glob::FuncGrid(mass, sigmaM, interpType);
  interp[1] = cbl::glob::FuncGrid(mass, dlnsigmaM, interpType);

  return interp;
}


// ===========================================================================================


double cbl::modelling::numbercounts::mass_function (const double mass, cbl::cosmology::Cosmology cosmology, const double redshift, const std::string model_MF, const bool store_output, const double Delta, const bool isDelta_vir, const cbl::glob::FuncGrid interp_Pk, const double kmax)
{
  const double rho = cosmology.rho_m(0., true);

  double sigmaM, dlnsigmaM;

  sigmaM_dlnsigmaM (sigmaM, dlnsigmaM, ((cosmology.unit()) ? mass : mass*cosmology.hh()), interp_Pk, kmax, rho);
  double _Delta = (isDelta_vir) ? cosmology.Delta_vir(Delta, redshift) : Delta;

  return cosmology.mass_function(mass, sigmaM, dlnsigmaM, redshift, model_MF, store_output, cbl::par::defaultString, _Delta);
}


// ===========================================================================================


std::vector<double> cbl::modelling::numbercounts::mass_function (const std::vector<double> mass, cbl::cosmology::Cosmology cosmology, const double redshift, const std::string model_MF, const bool store_output, const double Delta, const bool isDelta_vir, const std::vector<double> kk, const std::vector<double> Pk, const std::string interpType, const double kmax)
{
  vector<double> mass_function(mass.size());
  
  const double rho = cosmology.rho_m(0., true);

  vector<double> _mass = mass;

  if (!cosmology.unit())
    for(size_t i=0; i<mass.size(); i++)
      _mass[i] = mass[i]*cosmology.hh();

  vector<double> sigmaM, dlnsigmaM;
  cbl::modelling::numbercounts::sigmaM_dlnsigmaM (sigmaM, dlnsigmaM, _mass, kk, Pk, interpType, kmax, rho);

  double _Delta = (isDelta_vir) ? cosmology.Delta_vir(Delta, redshift) : Delta;

  for (size_t i=0; i<mass.size(); i++) 
    mass_function[i] = cosmology.mass_function(mass[i], sigmaM[i], dlnsigmaM[i], redshift, model_MF, store_output, cbl::par::defaultString, _Delta);

  return mass_function;
}


// ===========================================================================================


std::vector<std::vector<double>> cbl::modelling::numbercounts::mass_function (const std::vector<double> redshift, const std::vector<double> mass, cbl::cosmology::Cosmology cosmology, const std::string model_MF, const bool store_output, const double Delta, const bool isDelta_vir, const std::vector<double> kk, const std::vector<double> Pk, const std::string interpType, const double kmax)
{
  vector<vector<double>> mass_function(redshift.size(), vector<double>(mass.size()));
  
  const double rho = cosmology.rho_m(0., true);

  vector<double> _mass = mass;

  if (!cosmology.unit()) 
    for(size_t i=0; i<mass.size(); i++)
      _mass[i] = mass[i]*cosmology.hh();

  vector<double> sigmaM, dlnsigmaM;
  cbl::modelling::numbercounts::sigmaM_dlnsigmaM(sigmaM, dlnsigmaM, _mass, kk, Pk, interpType, kmax, rho);

  for (size_t j=0; j<redshift.size(); j++) 
    for (size_t i=0; i<mass.size(); i++) 
      mass_function[j][i] = cosmology.mass_function(mass[i], sigmaM[i], dlnsigmaM[i], redshift[j], model_MF, store_output, par::defaultString, ((isDelta_vir) ? cosmology.Delta_vir(Delta, redshift[j]) : Delta));

  return mass_function;
}


// ===========================================================================================


double cbl::modelling::numbercounts::number_counts (const double redshift_min, const double redshift_max, const double Mass_min, const double Mass_max, cbl::cosmology::Cosmology cosmology, const double Area, const std::string model_MF, const bool store_output, const double Delta, const bool isDelta_vir, const cbl::glob::FuncGrid interp_sigmaM, const  cbl::glob::FuncGrid interp_DlnsigmaM, const int npt_redshift, const int npt_mass)
{
  
  double fact = (cosmology.unit()) ? 1 : cosmology.hh();

  double deltaz = (redshift_max-redshift_min)/(npt_redshift);
  double deltaLogM = (log10(Mass_max)-log10(Mass_min))/(npt_mass);

  vector<double> MM(npt_mass), deltaM(npt_mass), sigmaM(npt_mass), dlnsigmaM(npt_mass);
  for (int j=0; j<npt_mass; j++) {
    double M1 = pow(10., log10(Mass_min)+j*deltaLogM);
    double M2 = pow(10., log10(Mass_min)+(j+1)*deltaLogM);
    MM[j] = pow(10., log10(Mass_min)+(j+0.5)*deltaLogM);
    deltaM[j] = M2-M1;
    sigmaM[j] = interp_sigmaM(MM[j]*fact);
    dlnsigmaM[j] = interp_DlnsigmaM(MM[j]*fact);
  }

  double nc = 0.;

  for (int i=0; i<npt_redshift; i++) {
    double zz = redshift_min+(i+0.5)*deltaz;
    double _Delta = (isDelta_vir) ? cosmology.Delta_vir(Delta, zz) : Delta;
    double dV_dZ = Area*cosmology.dV_dZdOmega(zz, true);

    double Int = 0;
    for (int j=0; j<npt_mass; j++) 
      Int += cosmology.mass_function(MM[j], sigmaM[j], dlnsigmaM[j], zz, model_MF, store_output, cbl::par::defaultString, _Delta)*deltaM[j]*dV_dZ;
    nc += Int*deltaz;

  }

  return nc;
  
  
  /*
  (void)npt_redshift; (void)npt_mass;

  auto integrand_nc_z = [&] (const double redshift)
  {
    double _Delta = (isDelta_vir) ? cosmology.Delta_vir(Delta, redshift) : Delta;
    double dV_dZ = Area*cosmology.dV_dZdOmega(redshift, true);

    auto integrand_nc_m = [&] (const double mass)
    { 
      return cosmology.mass_function (mass, interp_sigmaM(mass), interp_DlnsigmaM(mass), redshift, model_MF, cbl::par::defaultString, _Delta)*dV_dZ;
    };

    return wrapper::gsl::GSL_integrate_qag(integrand_nc_m, Mass_min, Mass_max);
  };

  return wrapper::gsl::GSL_integrate_qag(integrand_nc_z, redshift_min, redshift_max);
  */
}


// ===========================================================================================

std::vector<double> cbl::modelling::numbercounts::size_function (cbl::cosmology::Cosmology cosmology, const std::vector<double> radii, const double redshift, const std::string model, const double b_eff, double slope, double offset, const double deltav_NL, const double del_c, const std::string method_Pk, const bool store_output, const std::string output_root, const std::string interpType, const double k_max, const std::string input_file, const bool is_parameter_file)
{
  
  vector<double> size_function = cosmology.size_function(radii, redshift, model, b_eff, slope, offset, deltav_NL, del_c, method_Pk, store_output, output_root, interpType, k_max, input_file, is_parameter_file);

  return size_function;
  
}
