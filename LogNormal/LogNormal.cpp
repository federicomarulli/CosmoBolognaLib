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
 *  @file LogNormal/LogNormal.cpp
 *
 *  @brief Functions for the LogNormal data structure
 *
 *  This file contains the implementation of the \e methods of
 *  the class LogNormal
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unibo.it, alfonso.veropalumbo@unibo.it
 */

#include "FFTlog.h"
#include "LogNormal.h"

using namespace std;

using namespace cbl;
using namespace catalogue;
using namespace lognormal;


// ============================================================================


std::shared_ptr<catalogue::Catalogue> cbl::lognormal::LogNormal::catalogue (const size_t i)
{
  if (i<m_catalogue.size())
    return m_catalogue[i];

  else {
    ErrorCBL("The log-normal catalogue "+conv(i+1, par::fINT)+" does not exist!", "LogNormal_catalogue", "LogNormal.cpp");
    return NULL;
  }
}


// ============================================================================


void cbl::lognormal::LogNormal::generate (const int n_lognormal_mocks, const std::string output_dir, const std::string filename, const int start, const int seed)
{ 
  if (n_lognormal_mocks==0)  
    ErrorCBL("set number of log-normal realizations first!", "generate", "LogNormal.cpp");
  
  random::UniformRandomNumbers ran(0., 1., seed);
  
  
  // compute the visibility mask from the input random catalogue
  
  const double minX_random = m_random.Min(Var::_X_);
  const double minY_random = m_random.Min(Var::_Y_);
  const double minZ_random = m_random.Min(Var::_Z_);

  const double DeltaX = (m_random.Max(Var::_X_)-minX_random); 
  const double DeltaY = (m_random.Max(Var::_Y_)-minY_random); 
  const double DeltaZ = (m_random.Max(Var::_Z_)-minZ_random); 

  const int nTot = m_random.nObjects();

  int nx = DeltaX/m_cell_size; // nx: number of grid cells along the x axis
  nx = (nx%2==0) ? nx : nx+1; // if nx is odd, add 1 
  const double cell_size_x = DeltaX/nx;

  int ny = DeltaY/m_cell_size;
  ny = (ny%2==0) ? ny : ny+1; 
  const double cell_size_y = DeltaY/ny;

  int nz = DeltaZ/m_cell_size;
  nz = (nz%2==0) ? nz : nz+1; 
  const double cell_size_z = DeltaZ/nz;

  const int nzp = nz*0.5+1;

  const int nRtot = nx*ny*nz;
  const int nKtot = nx*ny*nzp;
  
  double *grid, *xxii;
  fftw_complex *ppkk;
  grid = fftw_alloc_real(nRtot);
  xxii = fftw_alloc_real(nRtot);
  ppkk = fftw_alloc_complex(nKtot);

  for (int i=0; i<nRtot; i++) {
    grid[i] = 0;
    xxii[i] = 0;
  }

  for (int i=0; i<nKtot; i++) {
    ppkk[i][0] = 0;
    ppkk[i][1] = 0;
  }
  
  for (int i=0; i<nTot; i++) { // loop on the random objects to put them in the grid cells
    int i1 = min(int((m_random.xx(i)-minX_random)/cell_size_x), nx-1);
    int j1 = min(int((m_random.yy(i)-minY_random)/cell_size_y), ny-1);
    int z1 = min(int((m_random.zz(i)-minZ_random)/cell_size_z), nz-1);
    long int index = z1+nz*(j1+ny*i1);
    grid[index] += 1./nTot; // grid[index] : the fraction of random objects in the index-th grid cell 
  }

  const double fact = (m_real) ? pow(m_bias, 2) : pow(m_bias, 2)*xi_ratio(m_cosmology.linear_growth_rate(m_redshift, 0.)/m_bias);
    
  vector<double> kG = logarithmic_bin_vector(500, 1.e-4, 1.e2), PkG;
    
  for (auto &&kg : kG)
    PkG.emplace_back(fact*m_cosmology.Pk_DM(kg, m_method_Pk, m_NL, m_redshift));


  cbl::glob::FuncGrid interpPk(kG, PkG, "Spline");

  const double xfact = 2.*par::pi/(nx*cell_size_x);
  const double yfact = 2.*par::pi/(ny*cell_size_y);
  const double zfact = 2.*par::pi/(nz*cell_size_z);
  
  const double Volume = nRtot*cell_size_x*cell_size_y*cell_size_z;
  
  // interpolate the power spectrum on the 3D grid
  for (int i=0; i<nx; i++) {
    const double kx = (i<=nx/2) ? xfact*i : -(nx-i)*xfact;
    for (int j=0; j<ny; j++) {
      const double ky = (j<=ny/2) ? yfact*j : -(ny-j)*yfact;
      for (int k=0; k<nzp; k++) {
	const double kz = zfact*k;
	const double kk = pow(kx*kx+ky*ky+kz*kz, 0.5);
	const long int kindex = k+nzp*(j+ny*i);
	ppkk[kindex][0] = interpPk(kk)/Volume; //interpolated(kk, kG, PkG, "Poly")/Volume;
      }
    }
  }

  coutCBL << "I'm computing Fourier transforms..."<<endl;
  
  
  // compute the Fourier transform of the 3D power spectrum divided by the volume to get the 3D two-point correlation function
  fftw_plan plan_xi_from_pk = fftw_plan_dft_c2r_3d(nx, ny, nz, ppkk, xxii, FFTW_ESTIMATE);
  fftw_execute(plan_xi_from_pk);
  fftw_destroy_plan(plan_xi_from_pk);

  for (int i=0; i<nRtot; i++)
    xxii[i] = log(1+xxii[i]);

  // compute the Fourier transform of the 3D two-point correlation function to get the 3D power spectrum
  fftw_plan plan_pk_from_xi = fftw_plan_dft_r2c_3d(nx, ny, nz, xxii, ppkk, FFTW_ESTIMATE);
  fftw_execute(plan_pk_from_xi);
  fftw_destroy_plan(plan_pk_from_xi);
  

  for (int nn=0; nn<n_lognormal_mocks; nn++) { // loop on the log-normal mocks to be constructed

    coutCBL << "I'm constructing the log-normal mock: " << nn+1 << "..." << endl;
    
    double *densX;
    fftw_complex *densK;
    densK = fftw_alloc_complex(nKtot);
    densX = fftw_alloc_real(nRtot);

    for (int i=0; i<nRtot; i++) 
      densX[i] = 0;
    
    for (int i=0; i<nKtot; i++) {
      densK[i][0] = 0;
      densK[i][1] = 0;
    }

    
    // compute the density field in Fourier space from the 3D power spectrum
    
    random::NormalRandomNumbers rang(0., 1., seed);
    
    for (int i=0; i<nx; i++) 
      for (int j=0; j<ny; j++) 
	for (int k=0; k<nzp; k++) {

	  const int kindex = k+nzp*(j+ny*i);
	  const double dk = sqrt(max(0., ppkk[kindex][0]/nRtot)*0.5);
	  
	  const double v1 = dk*rang();
	  const double v2 = dk*rang();

	  if (i==0 && j==0 && k==0) {
	    densK[kindex][0] = 0;
	    densK[kindex][1] = 0;
	  }
	  else if (i==nx/2 || j==ny/2){
	    densK[kindex][0] = v1;
	    densK[kindex][1] = 0.;
	  }
	  else {
	    densK[kindex][0] = v1;
	    densK[kindex][1] = v2;
	  }
	}

    
    // compute the Fourier transform of the 3D density field in Fourier space to the real space one
    fftw_plan plan_dr_from_dk = fftw_plan_dft_c2r_3d(nx, ny, nz, densK, densX, FFTW_ESTIMATE);
    fftw_execute(plan_dr_from_dk);
    fftw_destroy_plan(plan_dr_from_dk);
 
    // compute the variance of the density field
    vector<double> ff;
    for (int i=0; i<nx; i++) 
      for (int j=0; j<ny; j++) 
	for (int k=0; k<nz; k++) 
	  ff.push_back(densX[k+nz*(j+ny*i)]);
    
    const double sigma = Sigma(ff);

    //coutCBL << "Log-normal mock " << nn+1 << ": average density = " << Average(ff) << ", sigma = "<< sigma << endl;


    // draw a Poisson sampling of the grid density field

    default_random_engine generator;
    vector<shared_ptr<Object>> mock_sample;
    
    for (int i=0; i<nx; i++) {
      for (int j=0; j<ny; j++) {
	for (int k=0; k<nz; k++) {
	  long int index = k+nz*(j+ny*i);
	  
	  // number of objects in the (i, j, k) grid cell
	  const double nObjects = m_nObjects*grid[index]*(exp(densX[index]-sigma*sigma*0.5)); // transform the density field to get a log-normal distribution

	  // number of objects (from Poissonian PDF) whose coordinates will be randomly extracted (flat PDF for x,y,z) in the grid cells
	  poisson_distribution<int> distribution(nObjects);
	  int no = distribution(generator); 

	  for (int nnoo=0; nnoo<no; nnoo++) {
	    comovingCoordinates coord;
	    coord.xx = cell_size_x*(i+ran())+minX_random;
	    coord.yy = cell_size_y*(j+ran())+minY_random;
	    coord.zz = cell_size_z*(k+ran())+minZ_random;
	    shared_ptr<Galaxy> SMP(new Galaxy(coord));
	    mock_sample.push_back(SMP);
	  }

	}
      }
    }
    
    shared_ptr<Catalogue> p_cat(new Catalogue{mock_sample});
    
    p_cat->write_comoving_coordinates(output_dir+filename+"_"+conv(nn+start, par::fINT)+".dat");
    
    m_catalogue.emplace_back(p_cat);
  }
}
