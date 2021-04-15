/********************************************************************
 *  Copyright (C) 2021 by Federico Marulli and Giorgio Lesci        *
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
 *  @file Statistics/SuperSampleCovariance.cpp
 *
 *  @brief Methods of the class SuperSampleCovariance
 *
 *  This file contains the implementation of the methods of the class
 *  SuperSampleCovariance
 *
 *  @author Giorgio Lesci
 *
 *  @author giorgio.lesci2@unibo.it
 */

#include "SuperSampleCovariance.h"

// ======================================================================================


void cbl::statistics::SuperSampleCovariance::set_SSC (cbl::cosmology::Cosmology cosm, const std::vector<cbl::cosmology::CosmologicalParameter> cosmo_param, const std::vector<double> redshift_edges, const double area, const std::string method_Pk, const double delta_z, const double precision, const bool NL, const bool store_output)
{
  m_cosmo = std::make_shared<cosmology::Cosmology>(cosm);
  m_cosmo->set_unit(false); // force physical units
  
  m_cosmo_param = cosmo_param;
  m_method_Pk = method_Pk;
  m_NL = NL;
  m_store_output = store_output;
  
  m_nbins = (int)(redshift_edges.size()-1);
  m_area = area * pow(cbl::par::pi/180,2);
  m_precision = precision;
  
  m_compute_topHat_window(delta_z, redshift_edges);
}


// ======================================================================================


void cbl::statistics::SuperSampleCovariance::set_SSC (cbl::cosmology::Cosmology cosm, const std::vector<cbl::cosmology::CosmologicalParameter> cosmo_param, const double area, const std::vector<double> W_mean, const std::vector<double> W_std, const std::string method_Pk, const double delta_z, const double precision, const bool NL, const bool store_output)
{
  m_cosmo = std::make_shared<cosmology::Cosmology>(cosm);
  m_cosmo->set_unit(false); // force physical units
  
  m_cosmo_param = cosmo_param;
  m_method_Pk = method_Pk;
  m_NL = NL;
  m_store_output = store_output;
  
  m_nbins = (int)(W_mean.size());
  m_area = area * pow(cbl::par::pi/180,2);
  m_precision = precision;
  
  m_compute_gaussian_window(delta_z, W_mean, W_std);
}


// ======================================================================================


void cbl::statistics::SuperSampleCovariance::m_compute_topHat_window (const double delta_z, const std::vector<double> redshift_edges)
{
  m_nsteps = (int)((redshift_edges[redshift_edges.size()-1]-redshift_edges[0])/delta_z + 1);
  if (m_nsteps < m_nbins)
    cbl::ErrorCBL("m_nsteps can not be lower than the number of redshift bins!","m_compute_topHat_window","SuperSampleCovariance");

  m_redshifts.resize(m_nsteps);
  m_windows.resize(m_nbins);
  
  for (int i=0; i<m_nbins; i++) {
    m_windows[i].resize(m_nsteps);
    double iter = redshift_edges[0]-delta_z;
    double Delta_z = redshift_edges[i+1] - redshift_edges[i];
    for (int j=0; j<m_nsteps; j++) {
      iter += delta_z;
      m_redshifts[j] = iter;
      if ((iter>redshift_edges[i]) && (iter<=redshift_edges[i+1]))
	m_windows[i][j] = 1/Delta_z;
    }
  }
}


// ======================================================================================


void cbl::statistics::SuperSampleCovariance::m_compute_gaussian_window (const double delta_z, const std::vector<double> W_mean, const std::vector<double> W_std)
{
  double max_z = W_mean[W_mean.size()-1]+3.5*W_std[W_mean.size()-1];
  double min_z = W_mean[0]-3.5*W_std[0];
  
  m_nsteps = (int)((max_z - min_z) / delta_z + 1);
  m_redshifts.resize(m_nsteps);
  m_windows.resize(m_nbins);
  
  for (int i=0; i<m_nbins; i++) {
    m_windows[i].resize(m_nsteps);
    double iter = min_z - delta_z;
    for (int j=0; j<m_nsteps; j++) {
      iter += delta_z;
      m_redshifts[j] = iter;
      m_windows[i][j] = exp( - (m_redshifts[j] - W_mean[i]) * (m_redshifts[j] - W_mean[i]) / (2*W_std[i]*W_std[i]) ) / sqrt(2*cbl::par::pi*W_std[i]*W_std[i]);
    }
  }
}


// ======================================================================================


std::vector<std::vector<double>> cbl::statistics::SuperSampleCovariance::compute_Sij (cbl::cosmology::Cosmology cosmo) const
{
  // Compute comoving distances, volumes, growth factor
  std::vector<double> comov_dist(m_nsteps), dV(m_nsteps), growthf(m_nsteps);
  for (int i=0; i<m_nsteps; i++){
    comov_dist[i] = cosmo.D_C(m_redshifts[i]);
    dV[i] = comov_dist[i] * comov_dist[i] * cosmo.D_H() * cosmo.EE_inv(m_redshifts[i]); // D_H()*EE_inv(z) is the D_C normalized derivative
    growthf[i] = cosmo.gg(m_redshifts[i])/cosmo.gg(0)/(1+m_redshifts[i]);               // normalized scale-independent growth factor
  }
  
  // Compute normalizations
  std::vector<double> Inorm(m_nbins);
  std::vector<std::vector<double>> integrand (m_nbins,std::vector<double>(m_nsteps));
  for (int i=0; i<m_nbins; i++){
    for (int s=0; s<m_nsteps; s++){
      integrand[i][s] = dV[s] * m_windows[i][s] * m_windows[i][s] / 1.e10;
    }
    cbl::glob::FuncGrid integ(m_redshifts, integrand[i], "Spline");
    Inorm[i] = integ.integrate_cquad(m_redshifts[0], m_redshifts[m_redshifts.size()-1]) * 1.e10;
  }

  // Compute U(k)
  const double h = cosmo.hh();
  const double keq = 0.02/h; // Equality matter radiation in 1/Mpc (more or less)
  const double klogwidth = 10; // Factor of width of the integration range. 10 seems ok
  double max_comov_dist = comov_dist[comov_dist.size()-1];
  double min_comov_dist = comov_dist[0];

  double kmin = 0;
  double kmax = 0;
  if (keq<(1./max_comov_dist))
    kmin = keq/klogwidth;
  else
    kmin = (1./max_comov_dist)/klogwidth;
  if (keq>(1./min_comov_dist))
    kmax = keq*klogwidth;
  else
    kmax = min_comov_dist*klogwidth;
  
  const double nk = pow(2,m_precision); // m_precision=10 seems to be enough. Increase to test precision, reduce to speed up.
  const double logkmin = log(kmin);
  const double logkmax = log(kmax);
  std::vector<double> logk(nk);
  std::vector<double> kk(nk);
  for (size_t i=0; i<logk.size(); i++) {
    logk[i]=logkmin+i*(logkmax-logkmin)/(nk-1);
    kk[i]=exp(logk[i]);
  }
  
  std::vector<double> Pk_new = cosmo.Pk_DM(kk, m_method_Pk, m_NL, 0., m_store_output, "test", -1, kk[0], kk[1], 1.e-2, cbl::par::defaultString, false);
  std::vector<std::vector<double>> Uarr(m_nbins, std::vector<double>(logk.size()));
  std::vector<double> kr(m_nsteps);
  std::vector<std::vector<double>> integrand2(m_nbins, std::vector<double>(m_nsteps));  
  for (int i=0; i<m_nbins; i++) {
    for (size_t j=0; j<logk.size(); j++) {
      for (size_t s=0; s<m_redshifts.size(); s++) {
	kr[s] = kk[j]*comov_dist[s];
	integrand2[i][s] = dV[s] * m_windows[i][s] * m_windows[i][s] * growthf[s] * sin(kr[s]) / kr[s] / 1.e10;
      }
      cbl::glob::FuncGrid integ(m_redshifts, integrand2[i], "Spline");
      Uarr[i][j] = integ.integrate_cquad(m_redshifts[0], m_redshifts[m_redshifts.size()-1]) * 1.e10;
    }
  }
  
  // Compute S_ij
  std::vector<std::vector<double>> Cl_zero(m_nbins,std::vector<double>(m_nbins));
  std::vector<std::vector<double>> U1(m_nbins, std::vector<double>(logk.size()));
  std::vector<std::vector<double>> U2(m_nbins, std::vector<double>(logk.size()));
  std::vector<std::vector<double>> integrand3(m_nbins, std::vector<double>(logk.size()));
  for (int i=0; i<m_nbins; i++){
    for (size_t j=0; j<logk.size(); j++){
      U1[i][j] = Uarr[i][j]/Inorm[i];
    }
    for (int k=i; k<m_nbins; k++){
      for (size_t j=0; j<logk.size(); j++){
	U2[k][j] = Uarr[k][j]/Inorm[k];
	integrand3[k][j] = kk[j] * kk[j] * Pk_new[j] * U1[i][j] * U2[k][j] * 1.e10;
      }
      cbl::glob::FuncGrid integ(kk, integrand3[k], "Spline");
      Cl_zero[i][k] = 2/cbl::par::pi * integ.integrate_cquad(kk[0], kk[kk.size()-1]) / 1.e10;
    }
  }

  // Fill by symmetry
  for (int i=0; i<m_nbins; i++){
    for (int j=0; j<m_nbins; j++){
      Cl_zero[i][j] = Cl_zero[std::min(i,j)][std::max(i,j)];
    }
  }
  std::vector<std::vector<double>> Sij(m_nbins, std::vector<double>(m_nbins));
  for (int i=0; i<m_nbins; i++){
    for (int j=0; j<m_nbins; j++){
      Sij[i][j] = Cl_zero[i][j] / m_area; // S_ij is given by dividing Cl_zero/4pi by the fraction of the sky, i.e. Area(steradians)/4pi
    }
  }  

  return Sij;
}


// ======================================================================================


std::vector<std::vector<double>> cbl::statistics::SuperSampleCovariance::operator () (std::vector<double> &parameter) const
{
  // Redefine the cosmology
  cbl::cosmology::Cosmology cosmo = *m_cosmo;

  // Set the cosmological parameters
  for (size_t i=0; i<m_cosmo_param.size(); ++i)
    cosmo.set_parameter(m_cosmo_param[i], parameter[i]);

  std::vector<std::vector<double>> Sij = compute_Sij(cosmo);
  
  return Sij;
}



// ======================================================================================


std::vector<std::vector<double>> cbl::statistics::SuperSampleCovariance::get_window_function ()
{
  return m_windows;
}



// ======================================================================================


void cbl::statistics::SuperSampleCovariance::write_window_function (const std::string dir, const std::string file)
{
  // Create the directory
  std::string mkdir = "mkdir -p "+dir; if (system(mkdir.c_str())) {}
  
  // Write on file
  std::string file_out = dir+file;
  std::ofstream fout(file_out.c_str());

  int precision = 6;

  fout << "### [1] i # [2] j # [3] window[i][j] # [4] redshift " << std::endl;

  for (size_t i=0; i<m_windows.size(); ++i)
    for (size_t j=0; j<m_windows[i].size(); ++j)
      fout << setiosflags(std::ios::fixed) << std::setprecision(precision) << std::setw(15) << std::right << i
	   << "  " << setiosflags(std::ios::fixed) << std::setprecision(precision) << std::setw(15) << std::right << j
	   << "  " << setiosflags(std::ios::fixed) << std::setprecision(precision) << std::setw(15) << std::right << m_windows[i][j]
	   << "  " << setiosflags(std::ios::fixed) << std::setprecision(precision) << std::setw(15) << std::right << m_redshifts[j] << std::endl;
   
  fout.close(); std::cout << std::endl; coutCBL << "I wrote the file: " << file_out << std::endl;
}


// ======================================================================================


void cbl::statistics::SuperSampleCovariance::write_Sij (const std::string dir, const std::string file)
{
  // Compute S_ij using the cosmology given in input to set_SSC
  std::vector<std::vector<double>> Sij = compute_Sij(*m_cosmo);

  // Create the directory
  std::string mkdir = "mkdir -p "+dir; if (system(mkdir.c_str())) {}
  
  // Write on file
  std::string file_out = dir+file;
  std::ofstream fout(file_out.c_str());

  int precision = 15;

  fout << "### [1] i # [2] j # [3] S_ij " << std::endl;

  for (size_t i=0; i<Sij.size(); ++i) 
    for (size_t j=0; j<Sij[i].size(); ++j) 
      fout << setiosflags(std::ios::fixed) << std::setprecision(precision) << std::setw(15) << std::right << i
	   << "  " << setiosflags(std::ios::fixed) << std::setprecision(precision) << std::setw(15) << std::right << j
	   << "  " << setiosflags(std::ios::fixed) << std::setprecision(precision) << std::setw(15) << std::right << Sij[i][j] << std::endl;
   
  fout.close(); std::cout << std::endl; coutCBL << "I wrote the file: " << file_out << std::endl;
}
