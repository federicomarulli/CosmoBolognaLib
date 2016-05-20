/*******************************************************************
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

/** @file LogNormal/LogNormal.cpp
 *
 *  @brief Functions for the LogNormal data structure
 *
 *  This file contains the implementation of the \e methods of
 *  the class LogNormal
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#include "LogNormal.h"

using namespace cosmobl;
using namespace catalogue;


// ============================================================================


void cosmobl::LogNormal::setCatalogues (const shared_ptr<Catalogue> data, const shared_ptr<Catalogue> random)
{
  m_data = data;
  m_random = random;
}


// ============================================================================


void cosmobl::LogNormal::setParameters_from_xi (const vector<double> rr, const vector<double> xi) 
{
  // TBD: add parameters for extrapolation
  m_rmodel = rr;
  m_ximodel = xi;
  m_withxi = 1;
}


// ============================================================================


void cosmobl::LogNormal::setParameters_from_model (const shared_ptr<Cosmology> cosmology, const double bias, const bool Real, const string author, const bool NL, const string model)
{ 
  m_cosmology = cosmology;
  m_bias = bias;
  m_Real = Real;
  m_author= author;
  m_NL = NL;
  m_model = model;
  m_withxi=0;
}


// ============================================================================


void cosmobl::LogNormal::generate_LogNormal_mock (const double rmin, const string dir, const int start, const string filename)
{ 
  if (m_nLN==0)  
    ErrorMsg("Error in cosmobl::LogNormal::generate_LogNormal_mock of LogNormal.cpp, set number of LN realization first!");

  default_random_engine gen;
  uniform_real_distribution<float> ran(0., 1.);
  
  m_rmin = rmin;

  
  // compute the visibility mask
  
  double DeltaX = (m_random->Max(Var::_X_)-m_random->Min(Var::_X_)); 
  double DeltaY = (m_random->Max(Var::_Y_)-m_random->Min(Var::_Y_)); 
  double DeltaZ = (m_random->Max(Var::_Z_)-m_random->Min(Var::_Z_)); 

  int nTot = m_random->nObjects();
  
  vector<double> stat;
  m_data->stats_var(Var::_Redshift_, stat);

  int nx = DeltaX/m_rmin;
  nx = (nx%2==0) ? nx : nx+1; 
  double xmin = DeltaX/nx;

  int ny = DeltaY/m_rmin;
  ny = (ny%2==0) ? ny : ny+1; 
  double ymin = DeltaY/ny;

  int nz = DeltaZ/m_rmin;
  nz = (nz%2==0) ? nz : nz+1; 
  double zmin = DeltaZ/nz;

  int nzp = nz*0.5+1;

  int nRtot = nx*ny*nz;
  int nKtot = nx*ny*nzp;

  double VV = nRtot*zmin*ymin*xmin; 

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

  for (int i=0;i<nTot; i++) {
    int i1 = min(int((m_random->xx(i)-m_random->Min(Var::_X_))/xmin), nx-1);
    int j1 = min(int((m_random->yy(i)-m_random->Min(Var::_Y_))/ymin), ny-1);
    int z1 = min(int((m_random->zz(i)-m_random->Min(Var::_Z_))/zmin), nz-1);
    long int index = z1+nz*(j1+ny*i1);
    grid[index] += 1./nTot;
  }

  if (!m_withxi) { // with the model xi(r)
    vector<double> PkG;
    vector<double> kG = linear_bin_vector(500, -4., 1.);

    double ff = m_cosmology->linear_growth_rate(stat[0]);
    double beta = ff/m_bias;
    double fact = pow(m_bias,2);
    fact = (m_Real) ? fact : fact*(1+2.*beta/3.+0.2*beta*beta);
    cout << fact << "   " << pow(m_bias,2) << endl;

    for (size_t i=0; i<kG.size(); i++) {
      kG[i] = pow(10,kG[i]);
      PkG.push_back(fact*m_cosmology->Pk(kG[i], m_author, m_NL, stat[0], m_model));
    }

    double xfact = 2*par::pi/(nx*xmin);
    double yfact = 2*par::pi/(ny*ymin);
    double zfact = 2*par::pi/(nz*zmin);

    for (int i=0; i<nx; i++) {
      double kx = (i<=nx/2) ? xfact*i : -(nx-i)*xfact;
      for (int j=0; j<ny; j++) {
	double ky =(j<=ny/2) ? yfact*j : -(ny-j)*yfact;
	for (int k=0; k<nzp; k++) {
	  double kz = zfact*k;
	  double kk = pow(kx*kx+ky*ky+kz*kz, 0.5);
	  long int kindex = k+nzp*(j+ny*i);
	  
	  ppkk[kindex][0] = interpolated(kk, kG, PkG, "Poly")/VV;
	}
      }
    }

     
    fftw_plan pk2xi;
    pk2xi = fftw_plan_dft_c2r_3d(nx, ny, nz, ppkk, xxii, FFTW_ESTIMATE);
    fftw_execute(pk2xi);
    fftw_destroy_plan(pk2xi);

    for (int i=0; i<nRtot; i++)
      xxii[i] = log(1+xxii[i]);
    
  }
  
  else 
    ErrorMsg("Work in progres in cosmobl::LogNormal::generate_LogNormal_mock of LogNormal.cpp");
  

  fftw_plan xi2pk;
  xi2pk = fftw_plan_dft_r2c_3d(nx, ny, nz, xxii, ppkk, FFTW_ESTIMATE);
  fftw_execute(xi2pk);
  fftw_destroy_plan(xi2pk);
   
  cout << "Ready to extract Mocks" << endl;

  for (int nn=0; nn<m_nLN; nn++) {
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

   
    Erf erf;

    for (int i=0; i<nx; i++) {
      for (int j=0; j<ny; j++) {
	for (int k=0; k<nzp; k++) {

	  int kindex = k+nzp*(j+ny*i);
	  double Pkk = max(0., ppkk[kindex][0]/nRtot);
	  Pkk = sqrt(Pkk/2);
	  double v1 = sqrt(2.)*Pkk*erf.inverf(2*ran(gen)-1);
	  double v2 = sqrt(2.)*Pkk*erf.inverf(2*ran(gen)-1);

	  if (i==0 && j==0 && k==0) {
	    densK[kindex][0] = 0;
	    densK[kindex][1] = 0;
	  }
	  else if (i == nx/2 || j == ny/2){
	    densK[kindex][0] = v1;
	    densK[kindex][1] = 0.;
	  }
	  else {
	    densK[kindex][0] = v1;
	    densK[kindex][1] = v2;
	  }

	}
      }
    }

    fftw_plan dk2dr;
    dk2dr = fftw_plan_dft_c2r_3d(nx, ny, nz, densK, densX, FFTW_ESTIMATE);
    fftw_execute(dk2dr);
    fftw_destroy_plan(dk2dr);

    vector<double> ff;
    for (int i=0; i<nx; i++) {
      for (int j=0; j<ny; j++) {
	for (int k=0; k<nz;k++) {
	  long int index = k+nz*(j+ny*i);
	  ff.push_back(densX[index]);
	}
      }
    }
    double sigma = Sigma(ff); double average = Average(ff);
    cout << "Step " << nn+1 << " " <<"Sigma = "<< sigma<< " " << "Average = "<<average << endl;
   
    for (int i=0;i<nx;i++) {
      for (int j=0;j<ny;j++) {
	for (int k=0;k<nz;k++) {
	  long int index = k+nz*(j+ny*i);
	  densX[index] = double(m_data->nObjects())*grid[index]*(exp(densX[index]-sigma*sigma/2));
	}
      }
    }
    
    int number = nn+start;
    string cat = dir+filename+conv(number, par::fINT);

    default_random_engine generator;
    vector<shared_ptr<Object>> mock_sample;
    for (int i=0; i<nx; i++) {
      for (int j=0; j<ny; j++) {
	for (int k=0; k<nz; k++) {
	  long int index = k+nz*(j+ny*i);
	  double p = densX[index]; 
	  poisson_distribution<int> distribution(p);
	  int no = distribution(generator);

	  for (int nnoo = 0; nnoo<no; nnoo++) {
	    double XX = xmin*(i+ran(gen))+m_random->Min(Var::_X_);
	    double YY = ymin*(j+ran(gen))+m_random->Min(Var::_Y_);
	    double ZZ = zmin*(k+ran(gen))+m_random->Min(Var::_Z_);
	    shared_ptr<Galaxy> SMP(new Galaxy(XX, YY, ZZ, 1.));
	    mock_sample.push_back(SMP);
	  }

	}
      }
    }
    shared_ptr<Catalogue> p_cat(new Catalogue{mock_sample});

    m_LNCat[nn] = p_cat;
     
    m_LNCat[nn]->write_comoving_coordinates(cat);
    
  }

}


// ============================================================================


void cosmobl::LogNormal::set_nLN (const int nLN)
{ 
  m_nLN = nLN;
  m_LNCat.erase(m_LNCat.begin(),m_LNCat.end());
  m_LNCat.resize(m_nLN);
}
