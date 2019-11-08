/********************************************************************
 *  Copyright (C) 2010 by Federico Marulli and Alfonso Veropalumbo  *
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
 *  @file Cosmology/Lib/BAO.cpp
 *
 *  @brief Methods of the class Cosmology used to model the BAO
 *
 *  This file contains the implementation of the methods of the class
 *  Cosmology used to model the Baryon Acoustic Oscillations (BAO)
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @author federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#include "Cosmology.h"
#include "recombination.Recfast.h"
#include "cosmology.Recfast.h"

using namespace std;

using namespace cbl;


// =====================================================================================


//redshift at wich occurs baryon photon decoupling, see Hu & Sugiyama (1996).
double cbl::cosmology::Cosmology::z_decoupling() const
{
  double ombh2 = m_Omega_baryon*m_hh*m_hh;
  double omdmh2 = m_Omega_CDM*m_hh*m_hh;
  double g1 = 0.0783*pow(ombh2,-0.238)/(1.+39.5*pow(ombh2,0.763));
  double g2 = 0.560/(1+21.1*pow(ombh2,1.81));
  double zdec = 1048*(1.+0.00124*pow(ombh2,-0.738))*(1.+g1*pow(ombh2+omdmh2,g2));
  return zdec;

}


// =====================================================================================


double cbl::cosmology::Cosmology::z_drag() const
{
  double wb = m_Omega_baryon*m_hh*m_hh;
  double wb2 = wb*wb;
  double dNeff = m_massless_neutrinos+m_massive_neutrinos-3.046;
  double dNeff2 = dNeff*dNeff;

  double Yp = 0.2311+0.9520*wb-11.27*wb2+dNeff*(0.01356+0.008581*wb-0.1810*wb2)+dNeff2*(-0.0009795-0.001370*wb+0.01746*wb2);

  vector<double> params(14);
  int npz=10000;
  double zstart=1.0e+4, zend=0.001;

  params[0]=npz;
  params[1]=zstart;
  params[2]=zend;
  params[3] = Yp;  // Yp
  params[4] = par::TCMB;  // Temperature of CMB at z=0
  params[5] = m_Omega_matter;  // Omega matter 
  params[6] = m_Omega_baryon;  // Omega Baryons 
  params[7] = m_Omega_DE;  // Omega Lambda
  params[8] = m_Omega_k;  // Omega Curvature 
  params[9] = m_hh;  // h100
  params[10] = m_massless_neutrinos+m_massive_neutrinos; // effective number of neutrinos 
  params[11] = 1.14; // fudge-factor; normally F=1.14
 
  params[12] = 0.; // fDM [eV/s] which gives annihilation efficiency; 
              // typical value fDM=2.0e-24 eV/s (see Chluba 2010 for definitions)
  params[13] = 0; // switch on/off recombination corrections (Chluba & Thomas 2010)

  vector<double> zarr(npz), Xe_H(npz), Xe_He(npz), Xe(npz), TM(npz);

  Xe_frac(&params[0], &zarr[0], &Xe_H[0], &Xe_He[0], &Xe[0], &TM[0], 0);

  double HHc = m_hh*100*pow(cbl::par::kilo*cbl::par::pc, -1);
  double cc_m = par::cc*1000.;

  double rho_cr = 3.*pow(HHc, 2)/(8.*par::pi*par::GN);
  
  double rho_b = m_Omega_baryon*rho_cr;
  double rho_rad = 4.*par::sSB*pow(par::TCMB, 4)*pow(cc_m, -3);
  double R = 3./4*rho_b/rho_rad;

  double Np = rho_b*(1.-Yp)/par::mp;
  double kthom = 6.6524616e-29;

  vector<double> zz(npz), dtau(npz), tau(npz);
  for(int i=0; i<npz; i++){
    zz[i] = zarr[npz-1-i];
    double a = 1./(1+zz[i]);
    dtau[i] = Xe[npz-1-i]*Np*kthom/R*a/(HHc*EE(zz[i]));
  }

  glob::FuncGrid interp_dtau(zz, dtau, "Spline");

  auto integrand = [&] (double redshift) {return interp_dtau(redshift);};

  auto func = [&] (double redshift)
  {
    return wrapper::gsl::GSL_integrate_qag(integrand, 0., redshift);
  };

  return wrapper::gsl::GSL_root_brent (func, 1., 500, 2000);
  /*
  double wb = m_Omega_baryon*m_hh*m_hh;
  double wm = m_Omega_matter*m_hh*m_hh;

  double b1 = 0.313*pow(wm,-0.419)*(1+0.607*pow(wm,0.674));
  double b2 = 0.238*pow(wm,0.223);
  double zd = 1291.*pow(wm,0.251)*(1+b1*pow(wb,b2))/(1+0.659*pow(wm,0.828));
  return zd;
*/
}


// =====================================================================================

// Sound horizon at drag epoch

double cbl::cosmology::Cosmology::rs (const std::string method_Pk, const double T_CMB) const
{
  if (method_Pk=="EisensteinHu") 
    return rs_EH(T_CMB);
   
  else if (method_Pk=="CAMB")
    return rs_CAMB();

  else
    return ErrorCBL(" the input parameter method_Pk is not allowed!", "rs", "BAO.cpp");
}


// =====================================================================================

// Sound horizon at drag epoch (Eisentein & Hu 1998, Section 2.1)

double cbl::cosmology::Cosmology::rs_EH (const double T_CMB) const
{
  double Om0h2 = m_Omega_matter*pow(m_hh,2);
  double Ombh2 = m_Omega_baryon*pow(m_hh,2);
  double Tratio = T_CMB/2.7;

  double zeq = 2.5*pow(10.,4)*Om0h2*pow(Tratio,-4);    // Matter radiation Equivalence (explicit function?)
  double keq = 7.46*pow(10.,-2)*Om0h2*pow(Tratio,-2);  // Particle Horizon at zeq

  double b1 = 0.313 * pow(Om0h2,-0.419) * (1.+0.607*pow(Om0h2,0.674));
  double b2 = 0.238 * pow(Om0h2,0.223);

  double zdrag = 1291.*(pow(Om0h2,0.251)/(1.+0.659*pow(Om0h2,0.828)))*(1.+b1*pow(Ombh2,b2));

  double Rd = 31.5*Ombh2*pow(Tratio,-4)*pow(10.,3)/zdrag;
  double Req = 31.5*Ombh2*pow(Tratio,-4)*pow(10.,3)/zeq;

  double rs = 2.*pow(3.*keq,-1)*pow(6./Req,0.5)*log((pow(1.+Rd,0.5)+pow(Rd+Req,0.5))/(1.+pow(Req,0.5)));
  return ((m_unit) ? rs*m_hh : rs);
}


// =====================================================================================


double cbl::cosmology::Cosmology::rs_CAMB () const
{
  double wcb= m_Omega_matter*pow(m_hh,2);
  double wb= m_Omega_baryon*pow(m_hh,2);
  double wnu = m_Omega_neutrinos*pow(m_hh,2);

  double rd =  55.154*exp(-72.3*pow(wnu+0.0006,2))/(pow(wcb,0.25351)*pow(wb,0.12807));
  return ((m_unit) ? rd*m_hh : rd);
}


// =====================================================================================


double cbl::cosmology::Cosmology::ys (const double redshift, const std::string method_Pk, const double T_CMB) const
{
  return rs(method_Pk, T_CMB)/((m_unit) ? D_V(redshift)/m_hh : D_V(redshift));
}


// =====================================================================================


double cbl::cosmology::Cosmology::Az (const double redshift) const
{
  return ((m_unit) ? D_V(redshift)/m_hh : D_V(redshift))*1.e2*sqrt(m_Omega_matter*m_hh*m_hh)/(par::cc*redshift);
}


// =====================================================================================


double cbl::cosmology::Cosmology::sound_speed(const double redshift, const double T_CMB) const
{
  double rho_b = 3.*pow(100.*m_hh/par::cc, 2)*m_Omega_baryon; // Mpc^-2

  double cc_m = par::cc*1000.;
  double Mpc = par::pc*1.e6;
  double rho_rad = 8.*par::pi*par::GN*pow(cc_m, -2)*4.*par::sSB*pow(T_CMB, 4)*pow(cc_m, -3)*pow(Mpc, 2); // Mpc^-2
  
  double R = 3./4*rho_b/rho_rad/(1+redshift);
  double cs = 1./sqrt(3.*(1.+R));
  return par::cc*cs;
}


// =====================================================================================


double cbl::cosmology::Cosmology::rs_integrand (const double a, const double T_CMB) const
{
  double redshift=1./a-1;

  double zeq = 2.5e4*m_Omega_matter*m_hh*m_hh*pow(T_CMB/2.7,-4);
  double a_eq = 1./(1+zeq);

  double factor;
  if(m_Omega_radiation ==0)
    factor = sqrt(m_Omega_matter*(a+a_eq)+m_Omega_k*a*a+m_Omega_DE*f_DE(redshift)*pow(a,4));
  else
    factor = a*a*EE(redshift);

  return sound_speed(redshift, T_CMB)/factor;
}


// =====================================================================================


double cbl::cosmology::Cosmology::rs (const double redshift, const double T_CMB) const
{
  function<double(double)> integrand = bind(&Cosmology::rs_integrand, this, std::placeholders::_1, T_CMB);
  double a = 1./(1+redshift);
  return wrapper::gsl::GSL_integrate_qag(integrand, 0, a)/m_H0;
}


// =====================================================================================


vector<double> cbl::cosmology::Cosmology::linear_point (const double redshift, const double rmin, const double rmax, const int nbinr, const std::string interpType)
{
  vector<double> rr = linear_bin_vector(nbinr, rmin, rmax);

  vector<double> kk, Pk;
  run_CAMB(kk, Pk, false, redshift); 
  for(size_t i=0; i<kk.size(); i++){
    kk[i] = pow(10., kk[i]);
    Pk[i] = pow(10., Pk[i]);
  }

  vector<double> xi = cbl::wrapper::fftlog::transform_FFTlog(rr, 1, kk, Pk);

  cbl::glob::FuncGrid xi_interp(rr, xi, interpType);

  double rsCAMB = rs_CAMB();
  vector<double> boundaries = {rsCAMB-5, rsCAMB+5};

  
  // procedure to find the peak

  bool end = false;
  double rpeak, rdip;

  while (!end) {
    rpeak = xi_interp.root_D1v(boundaries[0], boundaries[1], 0, 1.e-10);

    if ((rpeak<boundaries[1]) && (rpeak > boundaries[0]))
      end = true;
    else if (rpeak==boundaries[1])
      boundaries[1]+=2;
    else 
      boundaries[0]-=2;
  }

  
  // procedure to find the dip

  end = false;
  boundaries[1] = boundaries[0];
  boundaries[0] = boundaries[1]-10;

  while (!end) {
    rdip = xi_interp.root_D1v(boundaries[0], boundaries[1], 0, 1.e-10);

    if ((rdip<boundaries[1]) &&(rdip > boundaries[0]))
      end = true;
    else if (rdip==boundaries[1])
      boundaries[1] += 2;
    else 
      boundaries[0] -= 2;
  }

  return {0.5*(rdip+rpeak), rdip, rpeak};
}
