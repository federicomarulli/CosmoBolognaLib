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
 *  @file Cosmology/Lib/Cosmology.cpp
 *
 *  @brief Generic methods of the class Cosmology  
 *
 *  This file contains the implementation of the \e generic methods of
 *  the class Cosmology
 *
 *  @author Federico Marulli
 *
 *  @author federico.marulli3@unibo.it
 */

#include "Cosmology.h"

using namespace std;

using namespace cbl;
using namespace glob;


// =====================================================================================


std::map<std::string, double> cbl::cosmology::Cosmology::parameter_value (std::string cosmological_model) const
{
  map<string, double> value;
  
  // the radiation density, given the redshift of radiation-matter equality
  auto OmegaR_zeq = [] (const double z_eq, const double Omega_matter) { return Omega_matter/(z_eq+1.); };

  // the density of massive neutrinos, given the neutrino mass
  auto Omega_neutrinos = [] (const double Mnu, const double h) { return Mnu/(93.14*pow(h, 2)); };

  // Komatsu et al. 2009: Table 1, WMAP 5 Year Mean
  if (cosmological_model=="WMAP5") {
    value["massless_neutrinos"] = 3.04;
    value["Omega_radiation"] = 2.469e-5*(1.+0.2271*value["massless_neutrinos"]); // Eq.8
    value["Omega_DE"] = 0.742;                                             // Omega_DE = 0.742 ± 0.030
    value["Omega_matter"] = 1.-value["Omega_radiation"]-value["Omega_DE"]; // flat LCDM (e.g. Eq. 11)
    value["Omega_baryon"] = 0.0441;                                        // Omega_b = 0.0441 ± 0.0030
    value["h"] = 0.719;                                                    // h = 0.719 ± 0.027 Km/s/Mpc
    value["scalar_pivot"] = 0.002;                                         // baseline
    value["scalar_amp"] = 2.41e-9;                                         // scalar amplitude = (2.41 ± 0.11)e-9 -> sigma8 = 0.796 ± 0.036
    value["n_spec"] = 0.963;                                               // n = 0.963 ± 0.015
    value["tau"] = 0.087;                                                  // tau = 0.087 ± 0.017
    value["Omega_neutrinos"] = 0.;
    value["massive_neutrinos"] = 0;
    value["w0"] = -1.;
    value["wa"] = 0.;
  }

  // Komatsu et al. 2011: Table 1, WMAP Seven-year Mean
  else if (cosmological_model=="WMAP7") {
    value["massless_neutrinos"] = 3.04;
    value["Omega_radiation"] = 2.469e-5*(1.+0.2271*value["massless_neutrinos"]);
    value["Omega_DE"] = 0.727;                                             // Omega_DE = 0.727 ± 0.030
    value["Omega_matter"] = 1.-value["Omega_radiation"]-value["Omega_DE"]; // flat LCDM
    value["Omega_baryon"] = 0.0455;                                        // Omega_b = 0.0455 ± 0.0028
    value["h"] = 0.704;                                                    // h = 0.704 ± 0.025 Km/s/Mpc
    value["scalar_amp"] = 2.43e-9;                                         // scalar amplitude = (2.43 ± 0.11)e-9 -> sigma8 = 0.811 ± 0.031
    value["scalar_pivot"] = 0.002;                                         // baseline
    value["n_spec"] = 0.968;                                               // n = 0.968 ± 0.012
    value["tau"] = 0.088;                                                  // tau = 0.088 ± 0.015
    value["Omega_neutrinos"] = 0.;
    value["massive_neutrinos"] = 0;
    value["w0"] = -1.;
    value["wa"] = 0.;
  }

  // Hinshaw et al. 2013: Table 3, WMAP-only Nine-year
  else if (cosmological_model=="WMAP9") {
    value["Omega_baryon"] = 0.0463;                                        // Omega_b = 0.0463 ± 0.0024
    value["Omega_DE"] = 0.721;                                             // Omega_DE = 0.721 ± 0.025
    value["Omega_matter"] = 0.233+value["Omega_baryon"];                   // Omega_CDM = 0.233 ± 0.023
    value["Omega_radiation"] = OmegaR_zeq(3273., value["Omega_matter"]);   // Table 2
    value["h"] = 0.70;                                                     // h = 0.700 ± 0.022 Km/s/Mpc
    value["scalar_amp"] = 2.41e-9;                                         // scalar amplitude = (2.41 ± 0.10)e-9 -> sigma8 = 0.821 ± 0.023
    value["scalar_pivot"] = 0.002;                                         // baseline
    value["n_spec"] = 0.972;                                               // n = 0.972 ± 0.013
    value["tau"] = 0.089;                                                  // tau = 0.089 ± 0.014
    value["Omega_neutrinos"] = 0.;
    value["massless_neutrinos"] = 3.04;
    value["massive_neutrinos"] = 0;
    value["w0"] = -1.;
    value["wa"] = 0.;
  }

  // Planck Collab 2013, Paper XVI: Table 2, Planck+WP
  else if (cosmological_model=="Planck13") {
    value["Omega_matter"] = 0.315;                                         // Omega_M = 0.315 ± 0.018
    value["Omega_baryon"] = 0.0487;                                        // Omega_b*h^2 = 0.02205 ± 0.00028
    value["Omega_radiation"] = OmegaR_zeq(3391., value["Omega_matter"]);   // z_eq = 3391 ± 60
    value["Omega_DE"] = 1.-value["Omega_matter"]-value["Omega_radiation"]; // assuming Omega_k = 0 (in Table 3: Omega_DE = 0.685 ± 0.018)
    value["h"] = 0.673;                                                    // h = 0.673 ± 0.012 Km/s/Mpc
    value["scalar_amp"] = 2.196e-9;                                        // scalar amplitude = (2.196 ± 0.06)e-9 -> sigma8 = 0.829 ± 0.012
    value["scalar_pivot"] = 0.05;                                          // baseline
    value["n_spec"] = 0.9603;                                              // n = 0.9603 ± 0.0073 
    value["tau"] = 0.089;                                                  // tau = 0.089 ± 0.014
    value["Omega_neutrinos"] = Omega_neutrinos(0.06, value["h"]);          // baseline (see Table 1)
    value["massless_neutrinos"] = 2.04;                                    // baseline (see Table 1)
    value["massive_neutrinos"] = 1;                                        // baseline (see Table 1)
    value["w0"] = -1.;
    value["wa"] = 0.;
  }
    
  // Planck Collab 2015, Paper XIII: Table 4, TT,TE,EE+lowP+lensing
  else if (cosmological_model=="Planck15") {
    value["Omega_matter"] = 0.3121;                                        // Omega_M = 0.3121 ± 0.0087
    value["Omega_baryon"] = 0.0488;                                        // Omega_b*h^2 = 0.02226 ± 0.00016
    value["Omega_radiation"] = OmegaR_zeq(3382., value["Omega_matter"]);   // z_eq = 3382 ± 32
    value["Omega_DE"] = 1.-value["Omega_matter"]-value["Omega_radiation"]; // assuming Omega_k = 0 (in Table 4: Omega_DE = 0.6879 ± 0.0087)
    value["h"] = 0.6751;                                                   // h = 0.6751 ± 0.0064 Km/s/Mpc
    value["scalar_amp"] = 2.13e-9;                                         // scalar amplitude = (2.130 ± 0.053)e-9 -> sigma8 = 0.8150 ± 0.0087
    value["scalar_pivot"] = 0.05;                                          // baseline
    value["n_spec"] = 0.9653;                                              // n = 0.9653 ± 0.0048
    value["tau"] = 0.063;                                                  // tau = 0.063 ± 0.014
    value["Omega_neutrinos"] = Omega_neutrinos(0.06, value["h"]);          // baseline
    value["massless_neutrinos"] = 2.04;                                    // baseline 
    value["massive_neutrinos"] = 1;                                        // baseline
    value["w0"] = -1.;
    value["wa"] = 0.;
  }
    
  // Planck Collab 2015, Paper XIII: Table 4, TT+lowP+lensing
  else if (cosmological_model=="Planck15_TT") {
    value["Omega_matter"] = 0.308;                                         // Omega_M = 0.308 ± 0.012
    value["Omega_baryon"] = 0.0484;                                        // Omega_b*h^2 = 0.02226 ± 0.0023
    value["Omega_radiation"] = OmegaR_zeq(3365., value["Omega_matter"]);   // z_eq = 3365 ± 44
    value["Omega_DE"] = 1.-value["Omega_matter"]-value["Omega_radiation"]; // assuming Omega_k = 0 (in Table 4: Omega_DE = 0.692 ± 0.012)
    value["h"] = 0.6781;                                                   // h = 0.6781 ± 0.092 Km/s/Mpc
    value["scalar_amp"] = 2.139e-9;                                        // scalar amplitude = (2.139 ± 0.063)e-9 -> sigma8 = 0.8149 ± 0.0093
    value["scalar_pivot"] = 0.05;                                          // baseline
    value["n_spec"] = 0.9677;                                              // n = 0.9677 ± 0.060
    value["tau"] = 0.066;                                                  // tau = 0.066 ± 0.016
    value["Omega_neutrinos"] = Omega_neutrinos(0.06, value["h"]);          // baseline
    value["massless_neutrinos"] = 2.04;                                    // baseline
    value["massive_neutrinos"] = 1;                                        // baseline
    value["w0"] = -1.;
    value["wa"] = 0.;
  }
    
  // Planck Collab 2018, Paper VI: Table 2, TT,TE,EE+lowE+lensing
  else if (cosmological_model=="Planck18") {
    value["Omega_matter"] = 0.3153;                                        // Omega_M = 0.3153 ± 0.0073
    value["Omega_baryon"] = 0.0486;                                        // Omega_b*h^2 = 0.02237 ± 0.00015
    value["Omega_radiation"] = OmegaR_zeq(3402., value["Omega_matter"]);   // z_eq = 3402 ± 26
    value["Omega_DE"] = 1.-value["Omega_matter"]-value["Omega_radiation"]; // assuming Omega_k = 0 (in Table 4: Omega_DE = 0.692 ± 0.012)
    value["h"] = 0.6736;                                                   // h = 0.6736 ± 0.054 Km/s/Mpc
    value["scalar_amp"] = 2.100e-9;                                        // scalar amplitude = (2.100 ± 0.030)e-9 -> sigma8 = 0.8111 ± 0.0060
    value["scalar_pivot"] = 0.05;                                          // baseline
    value["n_spec"] = 0.9649;                                              // n = 0.9649 ± 0.0042
    value["tau"] = 0.0544;                                                 // tau = 0.0544 ± 0.0073
    value["Omega_neutrinos"] = Omega_neutrinos(0.06, value["h"]);          // baseline
    value["massless_neutrinos"] = 2.04;                                    // baseline
    value["massive_neutrinos"] = 1;                                        // baseline
    value["w0"] = -1.;
    value["wa"] = 0.;
  }

  else {
    ErrorCBL("the chosen built-in cosmological model is not implemented", "Cosmology", "Cosmology.cpp");
  }
  
  return value;
}


// =====================================================================================


double cbl::cosmology::Cosmology::Omega_matter (const double redshift) const
{
  return m_Omega_matter/(EE(redshift)*EE(redshift))*pow(1.+redshift, 3);
}


// =====================================================================================


double cbl::cosmology::Cosmology::Omega_radiation (const double redshift) const 
{
  return m_Omega_radiation/(EE(redshift)*EE(redshift))*pow(1.+redshift, 4);
}


// =====================================================================================


void cbl::cosmology::Cosmology::set_parameters (const std::vector<std::string> parameter, const std::vector<double> value)
{
  for (size_t i=0; i<parameter.size(); i++)
    set_parameter(parameter[i], value[i]);
}


// =====================================================================================


double cbl::cosmology::Cosmology::Hubble (const double redshift) const 
{
  const double H0 = (m_unit) ? 100. : 100.*m_h;
  return H0*EE(redshift);
}


// =====================================================================================

double cbl::cosmology::Cosmology::DN (const double redshift, const double redshift_norm, const double prec) const 
{
  auto func = [prec, this] (const double aa) { return linear_growth_rate(1./aa-1., prec)/aa; };

  return exp(wrapper::gsl::GSL_integrate_qag(func, 1./(1.+redshift_norm), 1./(1.+redshift)));
}


// =====================================================================================


double cbl::cosmology::Cosmology::DD (const double redshift) const 
{
  const double aa = 1./(1.+redshift);
  
  auto func = [&] (const double aa) { return pow(aa*EE(1./aa-1.), -3); };
  
  return cbl::wrapper::gsl::GSL_integrate_qag(func, 0., aa)*2.5*Omega_matter()*EE(redshift);
}


// =====================================================================================


double cbl::cosmology::Cosmology::gg (const double redshift) const 
{
  return DD(redshift)*(1.+redshift);
}


// =====================================================================================


double cbl::cosmology::Cosmology::sigma8 (const double redshift) const 
{
  if (m_sigma8<0)
    ErrorCBL("sigma8 at z=0 is not set!", "sigma8", "Cosmology.cpp");

  return m_sigma8*DN(redshift);
}


// =====================================================================================


void cbl::cosmology::Cosmology::D_C_table (const std::string file_table, const double z_min, const double z_max, const int step, std::vector<double> &Redshift, std::vector<double> &dc) const
{
  cbl::Path path;
  string File_table = path.DirCosmo()+"Cosmology/Tables/dc/"+file_table;
 
  ifstream fin;
  fin.open (File_table.c_str());
  if (!fin) {

    ofstream fout(File_table.c_str()); checkIO(fout, File_table); 
    
    double delta_z = (z_max-z_min)/step;
    double z1 = z_min;
    double z2 = z_min+delta_z;

    for (int i=0; i<step; i++) {
      double zmean = (z1+z2)*0.5;
      fout << zmean << "   " << D_C(zmean) << endl;
      z1 = z2; z2 += delta_z;
    }
    
    fout.clear(); fout.close(); coutCBL <<"I wrote the file: "<<File_table<<endl;
  }
  fin.clear(); fin.close();
  
  fin.open(File_table.c_str());
  double Red, DC;
  while (fin >>Red>>DC) {
    Redshift.push_back(Red);
    dc.push_back(DC);
  }
  fin.clear(); fin.close();
}


// =====================================================================================


double cbl::cosmology::Cosmology::D_A (const double redshift) const 
{
  return D_M(redshift)/(1.+redshift);
}


// =====================================================================================


double cbl::cosmology::Cosmology::D_L (const double redshift) const 
{
  return (1.+redshift)*D_M(redshift);
}


// =====================================================================================


double cbl::cosmology::Cosmology::D_V (const double redshift) const 
{
  return pow(pow(D_M(redshift),2)*par::cc*redshift/Hubble(redshift),1./3.);
}


// =====================================================================================


double cbl::cosmology::Cosmology::F_AP (const double redshift) const 
{
  return D_M(redshift)*Hubble(redshift)/par::cc;
}


// =====================================================================================


double cbl::cosmology::Cosmology::Distance (const double redshift, const std::string distance_type) const 
{
  if (distance_type=="DC")
    return D_C(redshift);

  else if (distance_type=="DL")
    return D_L(redshift);

  else if (distance_type=="DA")
    return D_A(redshift);

  else if (distance_type=="Dv")
    return D_V(redshift);

  else if (distance_type=="Dvrs")
    return D_V(redshift)/sound_horizon_at_drag_epoch();

  else if (distance_type=="rsDv")
    return sound_horizon_at_drag_epoch()/D_V(redshift);

  else {
    ErrorCBL("no such a distance type!", "Distance", "Cosmology.cpp");
    return -1;
  }
}


// =====================================================================================


double cbl::cosmology::Cosmology::lookback_time (const double redshift) const 
{
  function<double(double)> integrand = [&] (const double redshift)
  { return 1./(1.+redshift)/EE(redshift); };
  
  double tt =  wrapper::gsl::GSL_integrate_qag(integrand, 0, redshift); 

  double Mpc = par::mega*par::pc*1.e-3; // in Km;
  double Gyr = par::giga*par::yr; // in sec

  return 1./(m_h*100.)*tt*Mpc/Gyr;
}


// =====================================================================================


double cbl::cosmology::Cosmology::cosmic_time (const double redshift) const
{
  function<double(double)> integrand = [&] (const double a)
  {
    const double redshift = 1./a-1.;
    return (1.+redshift)/EE(redshift);
  };
    
  const double aa = 1./(1.+redshift);
  const double tt = wrapper::gsl::GSL_integrate_qag(integrand, 0., aa); 

  const double Mpc = par::mega*par::pc*1.e-3; // in km
  const double Gyr = par::giga*par::yr; // in sec

  return 1./(m_h*100.)*tt*Mpc/Gyr;
}


// =====================================================================================


double cbl::cosmology::Cosmology::Hdot (const double redshift) const
{
  return -pow(Hubble(redshift), 2)*(1.+qq(redshift));
}


// =====================================================================================


double cbl::cosmology::Cosmology::rho_crit (const double redshift, const bool unit1) const
{
  // km sec^-1 Mpc^-1 -> sec^-1
  double HHc = Hubble(redshift)*pow(cbl::par::kilo*cbl::par::pc, -1);

  // force cosmological units
  if (!m_unit && unit1) HHc /= m_h;
  
  // m^3 Kg^-1 sec^-2 -> Mpc^3 Msun^-1 sec^-2
  const double GNc = cbl::par::GN*(cbl::par::Msol*pow(cbl::par::mega*cbl::par::pc, -3));
  
  return 3.*pow(HHc, 2)/(8.*par::pi*GNc);
}


// =====================================================================================


double cbl::cosmology::Cosmology::dV_dZdOmega (const double redshift, const bool angle_rad) const 
{
  // angle_rad: true -> Omega in steradians; false -> Omega in square degrees
  const double conv = (angle_rad) ? 1. : 3282.80635; 

  return D_H()*pow((1.+redshift)*D_A(redshift), 2)/EE(redshift)/conv;
}


// =====================================================================================


double cbl::cosmology::Cosmology::deltac (const double redshift) const 
{ 
  return 1.686*(1.+0.012299*log10(Omega_matter(redshift))); 
}


