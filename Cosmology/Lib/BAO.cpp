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
 *  @brief Methods of the class BAO used to model the features of the
 *  baryon acoustic oscillations
 *
 *  This file contains the implementation of the methods of the class
 *  BAO, used to model the baryon acoustic oscillations (BAO)
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @author federico.marulli3@unibo.it, alfonso.veropalumbo@unibo.it
 */

#include "BAO.h"
#include "recombination.Recfast.h"
#include "cosmology.Recfast.h"

using namespace std;

using namespace cbl;


// =====================================================================================


double cbl::cosmology::BAO::z_decoupling () const
{
  const double ombh2 = m_cosmology->Omega_baryon()*m_cosmology->little_h()*m_cosmology->little_h();
  const double omdmh2 = m_cosmology->Omega_CDM()*m_cosmology->little_h()*m_cosmology->little_h();
  const double g1 = 0.0783*pow(ombh2,-0.238)/(1.+39.5*pow(ombh2,0.763));
  const double g2 = 0.560/(1+21.1*pow(ombh2,1.81));
  const double zdec = 1048*(1.+0.00124*pow(ombh2,-0.738))*(1.+g1*pow(ombh2+omdmh2,g2));
  return zdec;
}


// =====================================================================================


double cbl::cosmology::BAO::z_drag () const
{
  const double wb = m_cosmology->Omega_baryon()*m_cosmology->little_h()*m_cosmology->little_h();
  const double wb2 = wb*wb;
  const double dNeff = m_cosmology->massless_neutrinos()+m_cosmology->massive_neutrinos()-3.046;
  const double dNeff2 = dNeff*dNeff;

  const double Yp = 0.2311+0.9520*wb-11.27*wb2+dNeff*(0.01356+0.008581*wb-0.1810*wb2)+dNeff2*(-0.0009795-0.001370*wb+0.01746*wb2);

  vector<double> params(14);
  int npz = 10000;
  const double zstart = 1.0e+4, zend = 0.001;

  params[0] = npz;
  params[1] = zstart;
  params[2] = zend;
  params[3] = Yp; 
  params[4] = par::TCMB;  
  params[5] = m_cosmology->Omega_matter(); 
  params[6] = m_cosmology->Omega_baryon(); 
  params[7] = m_cosmology->Omega_DE(); 
  params[8] = m_cosmology->Omega_k(); 
  params[9] = m_cosmology->little_h(); 
  params[10] = m_cosmology->massless_neutrinos()+m_cosmology->massive_neutrinos(); 
  params[11] = 1.14; // fudge-factor; normally F=1.14
 
  params[12] = 0.; // fDM [eV/s] which gives annihilation efficiency; 
  // typical value fDM=2.0e-24 eV/s (see Chluba 2010 for definitions)
  params[13] = 0; // switch on/off recombination corrections (Chluba & Thomas 2010)

  vector<double> zarr(npz), Xe_H(npz), Xe_He(npz), Xe(npz), TM(npz);

  Xe_frac(&params[0], &zarr[0], &Xe_H[0], &Xe_He[0], &Xe[0], &TM[0], 0);

  const double HHc = m_cosmology->little_h()*100*pow(cbl::par::kilo*cbl::par::pc, -1);
  const double cc_m = par::cc*1000.;

  const double rho_cr = 3.*pow(HHc, 2)/(8.*par::pi*par::GN);
  
  const double rho_b = m_cosmology->Omega_baryon()*rho_cr;
  const double rho_rad = 4.*par::sSB*pow(par::TCMB, 4)*pow(cc_m, -3);
  const double R = 3./4*rho_b/rho_rad;

  const double Np = rho_b*(1.-Yp)/par::mp;
  const double kthom = 6.6524616e-29;

  vector<double> zz(npz), dtau(npz), tau(npz);
  for (int i=0; i<npz; i++)
    dtau[i] = Xe[npz-1-i]*Np*kthom/R/(1+zarr[npz-1-i])/(HHc*m_cosmology->EE(zz[i]));

  glob::FuncGrid interp_dtau(zz, dtau, "Spline");

  auto integrand = [&] (double redshift) { return interp_dtau(redshift); };

  auto func = [&] (double redshift)
  {
    return wrapper::gsl::GSL_integrate_qag(integrand, 0., redshift);
  };

  return wrapper::gsl::GSL_root_brent(func, 1., 500, 2000);
  /*
    double wb = m_cosmology->Omega_baryon()*m_cosmology->little_h()*m_cosmology->little_h();
    double wm = m_cosmology->Omega_matter()*m_cosmology->little_h()*m_cosmology->little_h();

    double b1 = 0.313*pow(wm,-0.419)*(1+0.607*pow(wm,0.674));
    double b2 = 0.238*pow(wm,0.223);
    double zd = 1291.*pow(wm,0.251)*(1+b1*pow(wb,b2))/(1+0.659*pow(wm,0.828));
    return zd;
  */
}


// =====================================================================================


double cbl::cosmology::BAO::rs (const std::string method_Pk, const double T_CMB) const
{
  if (method_Pk=="EisensteinHu") 
    return m_cosmology->sound_horizon_at_drag_epoch("EisensteinHu",T_CMB);
   
  else if (method_Pk=="CAMB")
    return m_cosmology->sound_horizon_at_drag_epoch();

  else
    return ErrorCBL(" the input parameter method_Pk is not allowed!", "rs", "BAO.cpp");
}


// =====================================================================================


double cbl::cosmology::BAO::ys (const double redshift, const std::string method_Pk, const double T_CMB) const
{
  return rs(method_Pk, T_CMB)/((m_cosmology->unit()) ? m_cosmology->D_V(redshift)/m_cosmology->little_h() : m_cosmology->D_V(redshift));
}


// =====================================================================================


double cbl::cosmology::BAO::Az (const double redshift) const
{
  return ((m_cosmology->unit()) ? m_cosmology->D_V(redshift)/m_cosmology->little_h() : m_cosmology->D_V(redshift))*1.e2*sqrt(m_cosmology->Omega_matter()*m_cosmology->little_h()*m_cosmology->little_h())/(par::cc*redshift);
}


// =====================================================================================


double cbl::cosmology::BAO::sound_speed(const double redshift, const double T_CMB) const
{
  double rho_b = 3.*pow(100.*m_cosmology->little_h()/par::cc, 2)*m_cosmology->Omega_baryon(); // Mpc^-2

  double cc_m = par::cc*1000.;
  double Mpc = par::pc*1.e6;
  double rho_rad = 8.*par::pi*par::GN*pow(cc_m, -2)*4.*par::sSB*pow(T_CMB, 4)*pow(cc_m, -3)*pow(Mpc, 2); // Mpc^-2
  
  double R = 3./4*rho_b/rho_rad/(1+redshift);
  double cs = 1./sqrt(3.*(1.+R));
  return par::cc*cs;
}


// =====================================================================================


double cbl::cosmology::BAO::rs_integrand (const double a, const double T_CMB) const
{
  double redshift=1./a-1;

  double zeq = 2.5e4*m_cosmology->Omega_matter()*m_cosmology->little_h()*m_cosmology->little_h()*pow(T_CMB/2.7, -4);
  double a_eq = 1./(1+zeq);

  double factor;
  if(m_cosmology->Omega_radiation() ==0)
    factor = sqrt(m_cosmology->Omega_matter()*(a+a_eq)+m_cosmology->Omega_k()*a*a+m_cosmology->Omega_DE()*m_cosmology->f_DE(redshift)*pow(a,4));
  else
    factor = a*a*m_cosmology->EE(redshift);

  return sound_speed(redshift, T_CMB)/factor;
}


// =====================================================================================


double cbl::cosmology::BAO::rs (const double redshift, const double T_CMB) const
{
  function<double(double)> integrand = bind(&BAO::rs_integrand, this, std::placeholders::_1, T_CMB);
  double a = 1./(1+redshift);
  return wrapper::gsl::GSL_integrate_qag(integrand, 0, a)/m_cosmology->H0();
}


// =====================================================================================


vector<double> cbl::cosmology::BAO::linear_point (const double redshift, const double rmin, const double rmax, const int nbinr, const std::string interpType)
{
  vector<double> rr = linear_bin_vector(nbinr, rmin, rmax);

  vector<double> kk, Pk;

  PkXi PX(m_cosmology);
  
  PX.run_CAMB(kk, Pk, false, redshift); 
  for (size_t i=0; i<kk.size(); i++) {
    kk[i] = pow(10., kk[i]);
    Pk[i] = pow(10., Pk[i]);
  }

  vector<double> xi = cbl::wrapper::fftlog::transform_FFTlog(rr, 1, kk, Pk);

  cbl::glob::FuncGrid xi_interp(rr, xi, interpType);

  double rsCAMB = m_cosmology->sound_horizon_at_drag_epoch();
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
