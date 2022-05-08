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
 *******************************************************************/

/** 
 * @file LogNormal/LogNormalFull.cpp
 *
 * @brief Functions for the LogNormalFull data structure
 *
 * This file contains the implementation of the \e methods of
 * the class LogNormalFull
 *
 * @authors Federico Marulli, Alfonso Veropalumbo
 *
 * @authors federico.marulli3@unibo.it, alfonso.veropalumbo@unibo.it
 */

#include "LogNormalFull.h"

using namespace std;

using namespace cbl;
using namespace cosmology;
using namespace catalogue;


// ============================================================================


cbl::lognormal::LogNormalFull::LogNormalFull (const Cosmology cosmology, const double redshift_min, const double redshift_max, const int n_redshift_bins, const std::string author)
{
  set_cosmo_function(cosmology, redshift_min, redshift_max, n_redshift_bins, author);
}


// ============================================================================


cbl::lognormal::LogNormalFull::LogNormalFull (const double rmin, const double xMin, const double xMax, const double yMin, const double yMax, const double zMin, const double zMax, const Cosmology cosmology, const double redshift_min, const double redshift_max, const int n_redshift_bins, const std::string author)
{
  set_grid_parameters(rmin, xMin, xMax, yMin, yMax, zMin, zMax);

  set_cosmo_function(cosmology, redshift_min, redshift_max, n_redshift_bins, author);
}


// ============================================================================


cbl::lognormal::LogNormalFull::LogNormalFull (const double rmin, const std::vector<std::shared_ptr<catalogue::Catalogue>> random, const double pad, const Cosmology cosmology, const double redshift_min, const double redshift_max, const int n_redshift_bins, const std::string author)
{
  set_grid_parameters(rmin, random, pad);

  set_cosmo_function(cosmology, redshift_min, redshift_max, n_redshift_bins, author);
}


// ============================================================================


void cbl::lognormal::LogNormalFull::set_cosmo_function (const Cosmology cosmology, const double redshift_min, const double redshift_max, const int nredshift, const std::string author)
{
  m_cosmology = make_shared<Cosmology>(cosmology);

  m_author = author;

  vector<double> redshift = linear_bin_vector(nredshift, redshift_min, redshift_max);

  vector<double> dc, ff, dd, HH; 

  for (int i=0; i<nredshift;i++) {
    HH.push_back(cosmology.HH(redshift[i]));
    dc.push_back(cosmology.D_C(redshift[i]));
    ff.push_back(cosmology.linear_growth_rate(redshift[i], 0.));
    dd.push_back(cosmology.DD(redshift[i])/cosmology.DD(0.));
  }

  m_func_DC = make_shared<glob::FuncGrid>(glob::FuncGrid(redshift, dc, "Spline"));

  m_func_redshift = make_shared<glob::FuncGrid>(glob::FuncGrid(dc, redshift, "Spline"));

  m_func_HH = make_shared<glob::FuncGrid>(glob::FuncGrid(redshift, HH, "Spline"));

  m_func_growth_rate = make_shared<glob::FuncGrid>(glob::FuncGrid(redshift, ff, "Spline"));

  m_func_growth_factor = make_shared<glob::FuncGrid>(glob::FuncGrid(redshift, dd, "Spline"));
  
  int nk = 500;
  double kmin = 1.e-4;
  double kmax = 1.e2;
  vector<double> kk = logarithmic_bin_vector(nk, kmin, kmax);
  vector<double> Pk = m_cosmology->Pk_matter(kk, m_author, 0, 0.);

  m_func_pk = make_shared<glob::FuncGrid>(glob::FuncGrid(kk, Pk, "Spline"));
}


// ============================================================================


void cbl::lognormal::LogNormalFull::m_set_grid_parameters ()
{
  m_nx = (m_xMax-m_xMin)/m_rmin;
  m_nx = (m_nx%2==0) ? m_nx : m_nx+1; 

  m_ny = (m_yMax-m_yMin)/m_rmin;
  m_ny = (m_ny%2==0) ? m_ny : m_ny+1; 

  m_nz = (m_zMax-m_zMin)/m_rmin;
  m_nz = (m_nz%2==0) ? m_nz : m_nz+1; 

  m_nzF = m_nz*0.5+1;
}


// ============================================================================


void cbl::lognormal::LogNormalFull::set_grid_parameters (const double rmin, const double xMin, const double xMax, const double yMin, const double yMax, const double zMin, const double zMax)
{
  m_rmin = rmin;

  m_xMin = xMin;
  m_yMin = yMin;
  m_zMin = zMin;

  m_xMax = xMax;
  m_yMax = yMax;
  m_zMax = zMax;

  m_set_grid_parameters();
}


// ============================================================================


void cbl::lognormal::LogNormalFull::set_grid_parameters (const double rmin, const std::vector<std::shared_ptr<catalogue::Catalogue>> random, const double pad)
{
  m_rmin = rmin;

  m_random = random;

  vector<double> xMin, xMax, yMin, yMax, zMin, zMax;

  for (size_t i=0; i< m_random.size(); i++) {
    xMin.push_back(m_random[i]->Min(Var::_X_));
    yMin.push_back(m_random[i]->Min(Var::_Y_));
    zMin.push_back(m_random[i]->Min(Var::_Z_));

    xMax.push_back(m_random[i]->Max(Var::_X_));
    yMax.push_back(m_random[i]->Max(Var::_Y_));
    zMax.push_back(m_random[i]->Max(Var::_Z_));
  }

  m_xMin = Min(xMin)-pad;
  m_yMin = Min(yMin)-pad;
  m_zMin = Min(zMin)-pad;

  m_xMax = Max(xMax)+pad;
  m_yMax = Max(yMax)+pad;
  m_zMax = Max(zMax)+pad;

  m_set_grid_parameters();
}


// ============================================================================


void cbl::lognormal::LogNormalFull::m_set_fields (const bool use_random, const bool doRSD)
{
  m_clustering_signal = make_shared<data::ScalarField3D> (data::ScalarField3D(m_nx, m_ny, m_nz, m_xMin, m_xMax, m_yMin, m_yMax, m_zMin, m_zMax));

  m_densityG = make_shared<data::ScalarField3D> (data::ScalarField3D(m_nx, m_ny, m_nz, m_xMin, m_xMax, m_yMin, m_yMax, m_zMin, m_zMax));

  m_density = make_shared<data::ScalarField3D> (data::ScalarField3D(m_nx, m_ny, m_nz, m_xMin, m_xMax, m_yMin, m_yMax, m_zMin, m_zMax));

  if (doRSD) {
    m_los_velocity = make_shared<data::ScalarField3D> (data::ScalarField3D(m_nx, m_ny, m_nz, m_xMin, m_xMax, m_yMin, m_yMax, m_zMin, m_zMax));
    m_velocity = make_shared<data::VectorField3D> (data::VectorField3D(m_nx, m_ny, m_nz, m_xMin, m_xMax, m_yMin, m_yMax, m_zMin, m_zMax));
  }

  if (use_random)
    for (size_t i=0;i<m_random.size();i++) {
      m_visibility_random.push_back(make_shared<data::ScalarField3D> (data::ScalarField3D(m_nx, m_ny, m_nz, m_xMin, m_xMax, m_yMin, m_yMax, m_zMin, m_zMax)));
    }
  else
    m_visibility = make_shared<data::ScalarField3D> (data::ScalarField3D(m_nx, m_ny, m_nz, m_xMin, m_xMax, m_yMin, m_yMax, m_zMin, m_zMax));

}


// ============================================================================


void cbl::lognormal::LogNormalFull::m_reset_fields (const bool use_random, const bool doRSD)
{
  m_clustering_signal->reset();

  m_density->reset();

  if (doRSD) {
    m_los_velocity->reset();
    m_velocity->reset();
  }

  if (use_random) {
    m_visibility.reset();
    for (size_t i=0;i<m_random.size();i++) {
      m_visibility_random[i]->reset();
    }
  }
  else{
    m_visibility->reset();
    for (size_t i=0;i<m_random.size();i++) {
      m_visibility_random[i].reset();
    }
  }

}


// ============================================================================

/*
void cbl::lognormal::LogNormalFull::write_XiMultipoles(const std::vector<double> rr, const std::string file_out)
{
  vector<double> kk = logarithmic_bin_vector(200, 1.e-4, 1.e1); 
}
*/

// ============================================================================


void cbl::lognormal::LogNormalFull::m_set_clustering_signal ()
{
  double Volume = m_clustering_signal->Volume();

  for (int i=0; i<m_nx; i++) {
    double kx = m_clustering_signal->kX(i);
    for (int j=0; j<m_ny; j++) {
      double ky = m_clustering_signal->kY(j);
      for (int k=0; k<m_nzF; k++) {
	double kz = m_clustering_signal->kZ(k);

	double kk = sqrt(kx*kx+ky*ky+kz*kz);

	m_clustering_signal->set_ScalarField_FourierSpace_real(m_func_pk->operator()(kk)/Volume, i, j, k, 0);
	m_clustering_signal->set_ScalarField_FourierSpace_complex(0, i, j, k, 0);
      }
    }
  }

  m_clustering_signal->FourierAntiTransformField ();

  for (int i=0; i<m_nx; i++) { 
    for (int j=0; j<m_ny; j++) { 
      for (int k=0; k<m_nz; k++) {

	double value = log(m_clustering_signal->ScalarField(i, j, k)+1);
	m_clustering_signal->set_ScalarField(value, i, j, k);
      }
    }
  }

  m_clustering_signal->FourierTransformField ();
  
}


// ============================================================================

      
void cbl::lognormal::LogNormalFull::m_set_density_field (const double smoothing_radius)
{
  random::NormalRandomNumbers normal(0., 1., m_generator());

  for (int i=0; i<m_nx; i++) {
    double kx = m_density->kX(i);
    for (int j=0; j<m_ny; j++) {
      double ky = m_density->kY(j);
      for (int k=0; k<m_nzF; k++) {
	double kz = m_density->kZ(k);
	double kk = sqrt(kx*kx+ky*ky+kz*kz);

	if (i==0 && j==0 && k==0) {
	  m_densityG->set_ScalarField_FourierSpace_real(0, i, j, k, 0);
	  m_densityG->set_ScalarField_FourierSpace_complex(0, i, j, k, 0);
	}
	else if (i == m_nx/2 || j == m_ny/2) {

	  double smoothing = exp(-0.5*pow(kk*smoothing_radius, 2));
	  double Pkk = sqrt(max(0., 0.5*m_clustering_signal->ScalarField_FourierSpace_real(i, j, k)));

	  normal.set_mean_sigma(0., Pkk);

	  double val = normal()*smoothing;

	  while (val!=val) 
	    val = normal()*smoothing;

	  m_densityG->set_ScalarField_FourierSpace_real(val, i, j, k, 0);
	  m_densityG->set_ScalarField_FourierSpace_complex(0, i, j, k, 0);
	}
	else {
	  double smoothing = exp(-0.5*pow(kk*smoothing_radius,2));
	  double Pkk = sqrt(max(0., 0.5*m_clustering_signal->ScalarField_FourierSpace_real(i,j,k)));

	  normal.set_mean_sigma(0., Pkk);
	  double val1 = normal()*smoothing;
	  double val2 = normal()*smoothing;

	  while (val1!=val1) 
	    val1 = normal()*smoothing; 

	  while (val2!=val2) 
	    val2 = normal()*smoothing;

	  m_densityG->set_ScalarField_FourierSpace_real(val1, i, j, k, 0);
	  m_densityG->set_ScalarField_FourierSpace_complex(val2, i, j, k, 0);
	}

      }
    }
  }

  m_densityG->FourierAntiTransformField ();
  double mean = Average(m_densityG->ScalarField());
  m_sigma2G = pow(Sigma(m_densityG->ScalarField()), 2);

  coutCBL << "Average: "  << mean << " Sigma^2: " << m_sigma2G << endl;

  for (int i=0; i<m_nx; i++) {
    double xx = m_density->XX(i);
    for (int j=0; j<m_ny; j++) {
      double yy = m_density->YY(j);
      for (int k=0; k<m_nz; k++) {
	double zz = m_density->ZZ(k);

	double dc = sqrt(xx*xx+yy*yy+zz*zz);

	double redshift = m_func_redshift->operator()(dc);
	double fact = m_func_growth_factor->operator()(redshift);
	double delta = exp(fact*(m_densityG->ScalarField(i,j,k)-fact*0.5*m_sigma2G))-1;
	m_density->set_ScalarField(delta, i, j, k);
      }
    }
  }

  m_density->FourierTransformField();
}


// ============================================================================

      
void cbl::lognormal::LogNormalFull::m_set_radial_velocity ()
{
  for (int i=0; i<m_nx; i++) {
    double kx = m_clustering_signal->kX(i);
    for (int j=0; j<m_ny; j++) {
      double ky = m_clustering_signal->kY(j);
      for (int k=0; k<m_nzF; k++) {
	double kz = m_clustering_signal->kZ(k);

	double k2 = kx*kx+ky*ky+kz*kz;

	double vx_real = 0;
	double vy_real = 0;
	double vz_real = 0;

	double vx_complex = 0;
	double vy_complex = 0;
	double vz_complex = 0;


	if (i==0 && j==0 && k==0) {
	}
	else if (i == m_nx/2 || j == m_ny/2) {
	  vx_real = kx/k2*m_density->ScalarField_FourierSpace_complex(i,j,k);
	  vy_real = ky/k2*m_density->ScalarField_FourierSpace_complex(i,j,k);
	  vz_real = kz/k2*m_density->ScalarField_FourierSpace_complex(i,j,k);
	}
	else {
	  vx_real = kx/k2*m_density->ScalarField_FourierSpace_complex(i,j,k);
	  vy_real = ky/k2*m_density->ScalarField_FourierSpace_complex(i,j,k);
	  vz_real = kz/k2*m_density->ScalarField_FourierSpace_complex(i,j,k);

	  vx_complex = -kx/k2*m_density->ScalarField_FourierSpace_real(i,j,k);
	  vy_complex = -ky/k2*m_density->ScalarField_FourierSpace_real(i,j,k);
	  vz_complex = -kz/k2*m_density->ScalarField_FourierSpace_real(i,j,k);
	}

	m_velocity->set_VectorField_FourierSpace_real( {vx_real, vy_real, vz_real}, i, j, k, 0);
	m_velocity->set_VectorField_FourierSpace_complex( {vx_complex, vy_complex, vz_complex}, i, j, k, 0);
      }
    }
  }
  m_velocity->FourierAntiTransformField();

  for (int i=0; i<m_nx; i++) {
    double xx = m_velocity->XX(i);
    for (int j=0; j<m_ny; j++) {
      double yy = m_velocity->YY(j);
      for (int k=0; k<m_nz; k++) {
	double zz = m_velocity->ZZ(k);

	double rr2 = xx*xx+yy*yy+zz*zz;
	double rr = sqrt(rr2);

	double redshift = m_func_redshift->operator()(rr);
	double fact = m_func_growth_rate->operator()(redshift)*m_func_HH->operator()(redshift)*m_func_growth_factor->operator()(redshift);

	vector<double> vel = m_velocity->VectorField(i,j,k);

	double vr = fact*(xx*vel[0]+yy*vel[1]+zz*vel[2])/rr;

	m_los_velocity->set_ScalarField(vr, i, j, k, 0);
      }
    }
  }

}


// ============================================================================


void cbl::lognormal::LogNormalFull::m_set_visibility ()
{
  double mmax = max(max((m_xMax-m_xMin)/2, (m_yMax-m_yMin)/2),(m_xMax-m_xMin)/2); 
  m_nCells_eff=0;

  for (int i=0; i<m_nx; i++) {
    double xx = m_visibility->XX(i);
    for (int j=0; j<m_ny; j++) {
      double yy = m_visibility->YY(j);
      for (int k=0; k<m_nz; k++) {
	double zz = m_visibility->ZZ(k);

	double dc = sqrt(xx*xx+yy*yy+zz*zz);
	if (dc<mmax) {
	  m_nCells_eff +=1;
	  m_visibility->set_ScalarField(1., i, j, k, 0);
	}
	else {
	  m_visibility->set_ScalarField(0., i, j, k, 0);
	}
      }
    }
  }

  for (int i=0; i<m_nx; i++) 
    for (int j=0; j<m_ny; j++) 
      for (int k=0; k<m_nz; k++) {
	double vis = m_visibility->ScalarField(i, j, k)/m_nCells_eff;
	m_visibility->set_ScalarField(vis, i, j, k, 0);
      }
}


// ============================================================================


void cbl::lognormal::LogNormalFull::m_set_visibility_from_random ()
{
  coutCBL << "I'm setting the visibility from random sample..." << endl;
  int nsamples = m_random.size();
  double cell_size_inv = 1./m_rmin;

  for (int s=0; s<nsamples; s++) {
    m_visibility_random[s]->reset();

    int nObjects = m_random[s]->weightedN();  
    for (int obj = 0; obj<nObjects; obj++) {
      int i = min(int((m_random[s]->xx(obj)-m_xMin)*cell_size_inv), m_nx-1);
      int j = min(int((m_random[s]->yy(obj)-m_yMin)*cell_size_inv), m_ny-1);
      int k = min(int((m_random[s]->zz(obj)-m_zMin)*cell_size_inv), m_nz-1);

      double ww = m_random[s]->weight(obj)/nObjects;
      m_visibility_random[s]->set_ScalarField(ww, i, j, k, 1);
    }

  }

  coutCBL << "Done!" << endl;
}


// ============================================================================


void cbl::lognormal::LogNormalFull::m_set_visibility_from_random_RSD ()
{
  coutCBL << "I'm setting the visibility from random sample, using velocity field..." << endl;
  int nsamples = m_random.size();
  double cell_size_inv = 1./m_rmin;

  for (int s=0; s<nsamples; s++) {
    m_visibility_random[s]->reset();

    int nObjects = m_random[s]->nObjects();
    for (int obj = 0; obj<nObjects; obj++) {
      int i = min(int((m_random[s]->xx(obj)-m_xMin)*cell_size_inv), m_nx-1);
      int j = min(int((m_random[s]->yy(obj)-m_yMin)*cell_size_inv), m_ny-1);
      int k = min(int((m_random[s]->zz(obj)-m_zMin)*cell_size_inv), m_nz-1);

      double dc = m_func_DC->operator()(m_random[s]->redshift(obj));
      double dz = m_los_velocity->ScalarField(i, j, k)/par::cc;
      double newred = m_random[s]->redshift(obj)-dz*(1+m_random[s]->redshift(obj));
      double newdc = m_func_DC->operator()(newred);

      double newx = newdc*m_random[s]->xx(obj)/dc, newy=newdc*m_random[s]->yy(obj)/dc, newz = newdc*m_random[s]->zz(obj)/dc;

      i = min(int((newx-m_xMin)*cell_size_inv), m_nx-1);
      j = min(int((newy-m_yMin)*cell_size_inv), m_ny-1);
      k = min(int((newz-m_zMin)*cell_size_inv), m_nz-1);
      
      double ww = m_random[s]->weight(obj)/nObjects;
      m_visibility_random[s]->set_ScalarField(ww, i, j, k, 1);
    }

  }

  coutCBL << "Done!" << endl;
}


// ============================================================================


void cbl::lognormal::LogNormalFull::m_extract_points_lognormal_field (const double nObjects, const bool doRSD, const std::vector<double> redshift, const std::vector<double> bias, const std::shared_ptr<data::Field3D> visibility, const std::string file_out)
{
  ofstream fout(file_out.c_str());

  double max_redshift = m_func_redshift->operator()(sqrt(m_xMax*m_xMax+m_yMax*m_yMax+m_zMax*m_zMax));
  glob::FuncGrid func_bias(redshift, bias, "Spline");

  random::PoissonRandomNumbers poisson(1, m_generator());

  double min = -m_rmin*0.5, max = m_rmin*0.5;
  random::UniformRandomNumbers ran(min, max, m_generator());

  int nObj = 0;
  for (int i=0; i<m_nx; i++) {
    double xx = m_density->XX(i);
    for (int j=0; j<m_ny; j++) {
      double yy = m_density->YY(j);
      for (int k=0; k<m_nz; k++) {
	double zz = m_density->ZZ(k);

	double dc2 = xx*xx+yy*yy+zz*zz;
	double dc = sqrt(dc2);

	double redshift = m_func_redshift->operator()(dc);
	double bias = func_bias(redshift);

	double fact = bias*m_func_growth_factor->operator()(redshift);

	double delta = exp(fact*(m_densityG->ScalarField(i,j,k)-fact*0.5*m_sigma2G));
	double p = nObjects*visibility->ScalarField(i,j,k)*delta; 

	poisson.set_mean(p);
	int no = poisson();

	double vv = (doRSD) ? m_los_velocity->ScalarField(i, j, k)/par::cc : 0.;

	for (int n = 0; n<no; n++) {
	  comovingCoordinates comcoord;
	  observedCoordinates obscoord;
	  
	  comcoord.xx = ran()+xx; 
	  comcoord.yy = ran()+yy;
	  comcoord.zz = ran()+zz;
	  double ddcc = sqrt(comcoord.xx*comcoord.xx+comcoord.yy*comcoord.yy+comcoord.zz*comcoord.zz);

   	  polar_coord(comcoord.xx, comcoord.yy, comcoord.zz, obscoord.ra, obscoord.dec, ddcc);

	  obscoord.ra = obscoord.ra*180./par::pi;
	  obscoord.dec = obscoord.dec*180./par::pi;
	  
	  obscoord.ra = (obscoord.ra<0) ? obscoord.ra+360. : obscoord.ra;

	  obscoord.redshift =  m_func_redshift->operator()(ddcc);
	  obscoord.redshift -= vv*(1+obscoord.redshift);
	  if (obscoord.redshift >0 && obscoord.redshift<max_redshift) {
	    fout << obscoord.ra << " " << obscoord.dec << " " << obscoord.redshift << " " << endl;
	    nObj++;
	  }
	}
      }
    }
  }
  fout.clear(); fout.close();
  coutCBL << "Extracted "<< nObj << " points, written in " << file_out << endl;
}


// ============================================================================


void cbl::lognormal::LogNormalFull::generate_lognormal (const int start, const int stop, const bool doRSD, const double smoothing_radius, const std::vector<double> nObjects, const std::vector<std::vector<double>> redshift, const std::vector<std::vector<double> > bias, const std::string dir, const std::string filename, const int seed, const bool set_fields, const bool use_random)
{
  m_generator.seed(seed);

  if (set_fields)
    m_set_fields(use_random, doRSD);
  else
    m_reset_fields(use_random, doRSD);

  m_set_clustering_signal();

  int nsamples = nObjects.size();

  if (use_random)
    m_set_visibility_from_random();
  else
    m_set_visibility();

  for (int i=start; i<stop; i++) {
    m_set_density_field(smoothing_radius);

    if (doRSD) m_set_radial_velocity();

    for (int n=0;n<nsamples; n++) {
      string ff = dir+filename+"_sample_"+conv(n+1, par::fINT)+"_realization_"+conv(i+1, par::fINT);

      m_extract_points_lognormal_field(nObjects[n], doRSD, redshift[n], bias[n], (use_random) ? m_visibility_random[n] : m_visibility, ff);
    }

  }

}


// ============================================================================
/*
      
void cbl::lognormal::LogNormalFull::m_set_potential()
{
  double OmH02 = m_cosmology->OmegaM()*pow(m_cosmology->H0(),2);

  for (int i=0; i<m_nx; i++) {
    double kx = m_clustering_signal->kX(i);
    for (int j=0; j<m_ny; j++) {
      double ky = m_clustering_signal->kY(j);
      for (int k=0; k<m_nzF; k++) {
	double kz = m_clustering_signal->kZ(k);

	double k2 = kx*kx+ky*ky+kz*kz;

	double phi_real = 0;
	double phi_complex = 0;

	if (i==0 && j==0 && k==0) {
	}
	else if (i == m_nx/2 || j == m_ny/2) {
	  phi_real = -1.5*OmH02*m_density->ScalarField_FourierSpace_real(i,j,k)/k2;
	  phi_complex = 0;
	}
	else {
	  phi_real = -1.5*OmH02*m_density->ScalarField_FourierSpace_real(i,j,k)/k2;
	  phi_complex = -1.5*OmH02*m_density->ScalarField_FourierSpace_complex(i,j,k)/k2;
	}

	m_potential->set_ScalarField_FourierSpace_real( phi_real, i, j, k, 0);
	m_potential->set_ScalarField_FourierSpace_complex( phi_complex, i, j, k, 0);
      }
    }
  }
  m_potential->FourierAntiTransformField();
  
}


// ============================================================================


void cbl::lognormal::LogNormalFull::set_displacement_field ()
{

  m_density->FourierTransformField ();

  for (int i=0; i<m_nx; i++) {
    double kx = m_clustering_signal->kX(i);
    for (int j=0; j<m_ny; j++) {
      double ky = m_clustering_signal->kY(j);
      for (int k=0; k<m_nzF; k++) {
	double kz = m_clustering_signal->kZ(k);
	double kk2 = kx*kx+ky*ky+kz*kz;

	vector<double> displ_real(3,0), displ_complex(3,0.);

	if (i==0 && j==0 && k==0) {
	}
	else if (i == m_nx/2 || j == m_ny/2) {
	  displ_real[0] = kx/kk2*m_density->ScalarField_FourierSpace_complex(i,j,k);
	  displ_real[1] = ky/kk2*m_density->ScalarField_FourierSpace_complex(i,j,k);
	  displ_real[2] = kz/kk2*m_density->ScalarField_FourierSpace_complex(i,j,k);
	}
	else {
	  displ_real[0] = kx/kk2*m_density->ScalarField_FourierSpace_complex(i,j,k);
	  displ_real[1] = ky/kk2*m_density->ScalarField_FourierSpace_complex(i,j,k);
	  displ_real[2] = kz/kk2*m_density->ScalarField_FourierSpace_complex(i,j,k);

	  displ_complex[0] = -kx/kk2*m_density->ScalarField_FourierSpace_complex(i,j,k);
	  displ_complex[1] = -ky/kk2*m_density->ScalarField_FourierSpace_complex(i,j,k);
	  displ_complex[2] = -kz/kk2*m_density->ScalarField_FourierSpace_complex(i,j,k);  
	}

	m_displacement->set_VectorField_FourierSpace_real( displ_real, i, j, k, 0);
	m_displacement->set_VectorField_FourierSpace_complex( displ_complex, i, j, k, 0);


      }
    }
  }

  m_displacement->FourierAntiTransformField ();

  set_rsd_displacement_field();

  for (int i=0; i<m_nx; i++)
      for (int j=0; j<m_ny; j++)
	for (int k=0; k<m_nz; k++)
	  m_displacement->set_VectorField(m_rsd_displacement->VectorField(i,j,k), i, j, k, 0);
	
}


// ============================================================================

      
void cbl::lognormal::LogNormalFull::set_rsd_displacement_field ()
{
  for (int i=0; i<m_nx; i++) {
    double xx = m_displacement->XX(i);
      for (int j=0; j<m_ny; j++) {
	double yy = m_displacement->YY(j);
	for (int k=0; k<m_nz; k++) {
	  double zz = m_displacement->ZZ(k);

	  double dc2 = xx*xx+yy*yy+zz*zz;
	  double dc = sqrt(dc2);

	  double redshift = m_func_redshift->operator()(dc);
	  double ff = m_func_growth_rate->operator()(redshift);

	  vector<double> displ = m_displacement->VectorField(i, j, k);
	  double prod = displ[0]*xx+displ[1]*yy+displ[2]*zz;
	  double rsd_xx = ff*prod*xx/(dc2*(1.+redshift));
	  double rsd_yy = ff*prod*yy/(dc2*(1.+redshift));
	  double rsd_zz = ff*prod*zz/(dc2*(1.+redshift));
	  m_rsd_displacement->set_VectorField({rsd_xx, rsd_yy, rsd_zz}, i, j, k);
	}
      }
  }

}
*/


