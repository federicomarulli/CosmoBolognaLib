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

/** 
 *  @file Field/Field3D.cpp
 *
 *  @brief Functions for the Field3D data structure
 *
 *  This file contains the implementation of the \e methods of
 *  the class Field3D
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unibo.it, alfonso.veropalumbo@unibo.it
 */

#include "Field3D.h"

using namespace std;

using namespace cbl;
using namespace data;


// ============================================================================


cbl::data::Field3D::Field3D (const double deltaR, const double minX, const double maxX, const double minY, const double maxY, const double minZ, const double maxZ)
{
  set_parameters(deltaR, minX, maxX, minY, maxY, minZ, maxZ);
}


// ============================================================================


cbl::data::Field3D::Field3D (const int nx, const int ny, const int nz, const double minX, const double maxX, const double minY, const double maxY, const double minZ, const double maxZ)
{
  set_parameters(nx, ny, nz, minX, maxX, minY, maxY, minZ, maxZ);
}


// ============================================================================


void cbl::data::Field3D::set_parameters (const double deltaR, const double minX, const double maxX, const double minY, const double maxY, const double minZ, const double maxZ)
{
  m_MinX = minX;
  m_MaxX = maxX;
  m_MinY = minY;
  m_MaxY = maxY;
  m_MinZ = minZ;
  m_MaxZ = maxZ;

  double DeltaX = m_MaxX-m_MinX;
  double DeltaY = m_MaxY-m_MinY;
  double DeltaZ = m_MaxZ-m_MinZ;

  m_nX = DeltaX/deltaR;
  m_nX = (m_nX%2==0) ? m_nX : m_nX+1; 
  m_deltaX  = DeltaX/m_nX;

  m_nY = DeltaY/deltaR;
  m_nY = (m_nY%2==0) ? m_nY : m_nY+1; 
  m_deltaY  = DeltaY/m_nY;

  m_nZ = DeltaZ/deltaR;
  m_nZ = (m_nZ%2==0) ? m_nZ : m_nZ+1; 
  m_deltaZ  = DeltaZ/m_nZ;

  m_Volume = DeltaX*DeltaY*DeltaZ;

  m_nZF = m_nZ*0.5+1;

  m_nCells = m_nX*m_nY*m_nZ;
  m_nCells_Fourier = m_nX*m_nY*m_nZF;

// Set coordinates of the grid

  double xfact = 2*par::pi/DeltaX;
  double yfact = 2*par::pi/DeltaY;
  double zfact = 2*par::pi/DeltaZ;

  for(int i=0; i<m_nX; i++){
    m_X.push_back(m_MinX+(i+0.5)*m_deltaX);
    m_kX.push_back( (i<=m_nX/2) ? xfact*i : -(m_nX-i)*xfact);
  }

  for(int i=0; i<m_nY; i++){
    m_Y.push_back(m_MinY+(i+0.5)*m_deltaY);
    m_kY.push_back( (i<=m_nY/2) ? yfact*i : -(m_nY-i)*yfact);
  }

  for(int i=0; i<m_nZ; i++)
    m_Z.push_back(m_MinZ+(i+0.5)*m_deltaZ);

  for(int i=0; i<m_nZF; i++)
    m_kZ.push_back(zfact*i);

}


// ============================================================================


void cbl::data::Field3D::set_parameters (const int nx, const int ny, const int nz, const double minX, const double maxX, const double minY, const double maxY, const double minZ, const double maxZ)
{
  m_MinX = minX;
  m_MaxX = maxX;
  m_MinY = minY;
  m_MaxY = maxY;
  m_MinZ = minZ;
  m_MaxZ = maxZ;

  m_nX = (nx%2==0) ? nx : nx+1; 
  m_nY = (ny%2==0) ? ny : ny+1; 
  m_nZ = (nz%2==0) ? nz : nz+1; 

  m_nZF = m_nZ*0.5+1;

  double DeltaX = m_MaxX-m_MinX;
  double DeltaY = m_MaxY-m_MinY;
  double DeltaZ = m_MaxZ-m_MinZ;

  m_deltaX  = DeltaX/m_nX;
  m_deltaY  = DeltaY/m_nY;
  m_deltaZ  = DeltaZ/m_nZ;

  m_Volume = DeltaX*DeltaY*DeltaZ;

  m_nZF = m_nZ*0.5+1;

  m_nCells = m_nX*m_nY*m_nZ;
  m_nCells_Fourier = m_nX*m_nY*m_nZF;

// Set coordinates of the grid

  double xfact = 2*par::pi/DeltaX;
  double yfact = 2*par::pi/DeltaY;
  double zfact = 2*par::pi/DeltaZ;

  for(int i=0; i<m_nX; i++){
    m_X.push_back(m_MinX+(i+0.5)*m_deltaX);
    m_kX.push_back( (i<=m_nX/2) ? xfact*i : -(m_nX-i)*xfact);
  }

  for(int i=0; i<m_nY; i++){
    m_Y.push_back(m_MinY+(i+0.5)*m_deltaY);
    m_kY.push_back( (i<=m_nY/2) ? yfact*i : -(m_nY-i)*yfact);
  }

  for(int i=0; i<m_nZ; i++)
    m_Z.push_back(m_MinZ+(i+0.5)*m_deltaZ);

  for(int i=0; i<m_nZF; i++)
    m_kZ.push_back(zfact*i);
}


// ============================================================================


cbl::data::ScalarField3D::ScalarField3D (const double deltaR, const double minX, const double maxX, const double minY, const double maxY, const double minZ, const double maxZ) : Field3D(deltaR, minX, maxX, minY, maxY, minZ, maxZ) 
{
  m_field = fftw_alloc_real(m_nCells);
  m_field_FourierSpace = fftw_alloc_complex(m_nCells_Fourier);

  for(int i=0;i<m_nCells;i++)
    m_field[i]=0;

  for(int i=0;i<m_nCells_Fourier;i++){
    m_field_FourierSpace[i][0] = 0;
    m_field_FourierSpace[i][1] = 0;
  }
}


// ============================================================================


cbl::data::ScalarField3D::ScalarField3D (const int nx, const int ny, const int nz, const double minX, const double maxX, const double minY, const double maxY, const double minZ, const double maxZ) : Field3D(nx, ny, nz, minX, maxX, minY, maxY, minZ, maxZ) 
{
  m_field = fftw_alloc_real(m_nCells);
  m_field_FourierSpace = fftw_alloc_complex(m_nCells_Fourier);

  for (int i=0; i<m_nCells; i++)
    m_field[i] = 0;

  for (int i=0; i<m_nCells_Fourier; i++) {
    m_field_FourierSpace[i][0] = 0;
    m_field_FourierSpace[i][1] = 0;
  }

}


// ============================================================================


void cbl::data::ScalarField3D::FourierTransformField () 
{
  
  for (int i=0; i<m_nCells_Fourier; i++) {
    m_field_FourierSpace[i][0] = 0;
    m_field_FourierSpace[i][1] = 0;
  }

  fftw_plan real2complex;
  real2complex = fftw_plan_dft_r2c_3d(m_nX, m_nY, m_nZ, m_field, m_field_FourierSpace, FFTW_ESTIMATE);
  fftw_execute(real2complex);
  fftw_destroy_plan(real2complex);

  for(int i=0;i<m_nCells_Fourier;i++){
    m_field_FourierSpace[i][0] = m_field_FourierSpace[i][0]/m_nCells;
    m_field_FourierSpace[i][1] = m_field_FourierSpace[i][1]/m_nCells;
  }

}


// ============================================================================


void cbl::data::ScalarField3D::FourierAntiTransformField () 
{
  for (int i=0; i<m_nCells; i++)
    m_field[i] = 0;

  fftw_plan complex2real;
  complex2real = fftw_plan_dft_c2r_3d(m_nX, m_nY, m_nZ, m_field_FourierSpace, m_field, FFTW_ESTIMATE);
  fftw_execute(complex2real);
  fftw_destroy_plan(complex2real);

}


// ============================================================================


void cbl::data::ScalarField3D::GaussianConvolutionField (const double kernel_size) 
{
  FourierTransformField();

  double xfact = 2*par::pi/(m_nX*m_deltaX);
  double yfact = 2*par::pi/(m_nY*m_deltaY);
  double zfact = 2*par::pi/(m_nZ*m_deltaZ);

  double ks2 = kernel_size*kernel_size;

  for (int i=0; i<m_nX; i++) {
    double kx = (i<=m_nX/2) ? xfact*i : -(m_nX-i)*xfact;
    for (int j=0; j<m_nY; j++) {
      double ky = (j<=m_nY/2) ? yfact*j : -(m_nY-j)*yfact;
      for (int k=0; k<m_nZF; k++) {
	double kz = zfact*k;
	double kk2 = kx*kx+ky*ky+kz*kz;
	long int index = inds_to_index_Fourier(i,j,k);
	m_field_FourierSpace[index][0] = m_field_FourierSpace[index][0]*exp(-0.5*kk2*ks2);
	m_field_FourierSpace[index][1] = m_field_FourierSpace[index][1]*exp(-0.5*kk2*ks2);
      }
    }
  }

  FourierAntiTransformField();

}


// ============================================================================


void cbl::data::ScalarField3D::set_ScalarField (const double value, const int i, const int j, const int k, const bool add)
{
  m_field[inds_to_index(i,j,k)] = (add) ? m_field[inds_to_index(i,j,k)] + value : value;
}


// ============================================================================


void cbl::data::ScalarField3D::set_ScalarField_FourierSpace_real (const double value, const int i, const int j, const int k, const bool add)
{
  m_field_FourierSpace[inds_to_index_Fourier(i,j,k)][0] = (add) ? m_field_FourierSpace[inds_to_index_Fourier(i,j,k)][0] + value: value;
}


// ============================================================================


void cbl::data::ScalarField3D::set_ScalarField_FourierSpace_complex (const double value, const int i, const int j, const int k, const bool add)
{
  m_field_FourierSpace[inds_to_index_Fourier(i,j,k)][1] = (add) ? m_field_FourierSpace[inds_to_index_Fourier(i,j,k)][1] + value: value;
}


// ============================================================================


double cbl::data::ScalarField3D::ScalarField (const int i, const int j, const int k) const 
{
  return m_field[inds_to_index(i,j,k)];
}


// ============================================================================


double cbl::data::ScalarField3D::ScalarField_FourierSpace_real (const int i, const int j, const int k) const 
{
  return m_field_FourierSpace[inds_to_index_Fourier(i,j,k)][0];
}


// ============================================================================


double cbl::data::ScalarField3D::ScalarField_FourierSpace_complex (const int i, const int j, const int k) const 
{
  return m_field_FourierSpace[inds_to_index_Fourier(i,j,k)][1];
}


// ============================================================================


double cbl::data::ScalarField3D::ScalarField (const std::vector<double> pos) const 
{
  int i = min(int((pos[0]-m_MinX)/m_deltaX), m_nX-1);
  int j = min(int((pos[1]-m_MinY)/m_deltaY), m_nY-1);
  int k = min(int((pos[2]-m_MinZ)/m_deltaZ), m_nZ-1);

  return m_field[inds_to_index(i,j,k)];
}


// ============================================================================


vector<double> cbl::data::ScalarField3D::ScalarField () const 
{
  vector<double> vv(m_field, m_field+m_nCells);
  return vv;
}


// ============================================================================


void cbl::data::ScalarField3D::reset () 
{
  for(int i=0;i<m_nCells;i++)
    m_field[i]=0;

  for(int i=0;i<m_nCells_Fourier;i++){
    m_field_FourierSpace[i][0] = 0;
    m_field_FourierSpace[i][1] = 0;
  }
}

// ============================================================================


cbl::data::VectorField3D::VectorField3D (const double deltaR, const double minX, const double maxX, const double minY, const double maxY, const double minZ, const double maxZ) : Field3D(deltaR, minX, maxX, minY, maxY, minZ, maxZ) 
{
  m_field.resize(3);
  
  m_field_FourierSpace.resize(3);


  m_field[0] = fftw_alloc_real(m_nCells);
  m_field[1] = fftw_alloc_real(m_nCells);
  m_field[2] = fftw_alloc_real(m_nCells);

  m_field_FourierSpace[0] = fftw_alloc_complex(m_nCells_Fourier);
  m_field_FourierSpace[1] = fftw_alloc_complex(m_nCells_Fourier);
  m_field_FourierSpace[2] = fftw_alloc_complex(m_nCells_Fourier);

  for(int i=0;i<m_nCells;i++){
    m_field[0][i]=0;
    m_field[1][i]=0;
    m_field[2][i]=0;
  }
  for(int i=0;i<m_nCells_Fourier;i++){
    m_field_FourierSpace[0][i][0] = 0;
    m_field_FourierSpace[0][i][1] = 0;
    m_field_FourierSpace[1][i][0] = 0;
    m_field_FourierSpace[1][i][1] = 0;
    m_field_FourierSpace[2][i][0] = 0;
    m_field_FourierSpace[2][i][1] = 0;
  }
}


// ============================================================================


cbl::data::VectorField3D::VectorField3D (const int nx, const int ny, const int nz, const double minX, const double maxX, const double minY, const double maxY, const double minZ, const double maxZ) : Field3D(nx, ny, nz, minX, maxX, minY, maxY, minZ, maxZ) 
{
  m_field.resize(3);

  m_field_FourierSpace.resize(3);


  m_field[0] = fftw_alloc_real(m_nCells);
  m_field[1] = fftw_alloc_real(m_nCells);
  m_field[2] = fftw_alloc_real(m_nCells);

  m_field_FourierSpace[0] = fftw_alloc_complex(m_nCells_Fourier);
  m_field_FourierSpace[1] = fftw_alloc_complex(m_nCells_Fourier);
  m_field_FourierSpace[2] = fftw_alloc_complex(m_nCells_Fourier);

  for(int i=0;i<m_nCells;i++){
    m_field[0][i]=0;
    m_field[1][i]=0;
    m_field[2][i]=0;
  }
  for(int i=0;i<m_nCells_Fourier;i++){
    m_field_FourierSpace[0][i][0] = 0;
    m_field_FourierSpace[0][i][1] = 0;
    m_field_FourierSpace[1][i][0] = 0;
    m_field_FourierSpace[1][i][1] = 0;
    m_field_FourierSpace[2][i][0] = 0;
    m_field_FourierSpace[2][i][1] = 0;
  }

}


// ============================================================================


void cbl::data::VectorField3D::FourierTransformField () 
{

  for(int i=0;i<m_nCells_Fourier;i++){
    m_field_FourierSpace[0][i][0] = 0;
    m_field_FourierSpace[0][i][1] = 0;
    m_field_FourierSpace[1][i][0] = 0;
    m_field_FourierSpace[1][i][1] = 0;
    m_field_FourierSpace[2][i][0] = 0;
    m_field_FourierSpace[2][i][1] = 0;
  }

  fftw_plan real2complex;

  real2complex = fftw_plan_dft_r2c_3d(m_nX, m_nY, m_nZ, m_field[0], m_field_FourierSpace[0], FFTW_ESTIMATE);
  fftw_execute(real2complex);
  fftw_destroy_plan(real2complex);

  real2complex = fftw_plan_dft_r2c_3d(m_nX, m_nY, m_nZ, m_field[1], m_field_FourierSpace[1], FFTW_ESTIMATE);
  fftw_execute(real2complex);
  fftw_destroy_plan(real2complex);

  real2complex = fftw_plan_dft_r2c_3d(m_nX, m_nY, m_nZ, m_field[2], m_field_FourierSpace[2], FFTW_ESTIMATE);
  fftw_execute(real2complex);
  fftw_destroy_plan(real2complex);

  for(int i=0;i<m_nCells_Fourier;i++){
    m_field_FourierSpace[0][i][0] = m_field_FourierSpace[0][i][0]/m_nCells;
    m_field_FourierSpace[0][i][1] = m_field_FourierSpace[0][i][1]/m_nCells;
    m_field_FourierSpace[1][i][0] = m_field_FourierSpace[1][i][0]/m_nCells;
    m_field_FourierSpace[1][i][1] = m_field_FourierSpace[1][i][1]/m_nCells;
    m_field_FourierSpace[2][i][0] = m_field_FourierSpace[2][i][0]/m_nCells;
    m_field_FourierSpace[2][i][1] = m_field_FourierSpace[2][i][1]/m_nCells;
  }

}


// ============================================================================


void cbl::data::VectorField3D::FourierAntiTransformField () 
{
  for (int i=0; i<m_nCells; i++) {
    m_field[0][i] = 0;
    m_field[1][i] = 0;
    m_field[2][i] = 0;
  }

  fftw_plan complex2real;

  complex2real = fftw_plan_dft_c2r_3d(m_nX, m_nY, m_nZ, m_field_FourierSpace[0], m_field[0], FFTW_ESTIMATE);
  fftw_execute(complex2real);
  fftw_destroy_plan(complex2real);

  complex2real = fftw_plan_dft_c2r_3d(m_nX, m_nY, m_nZ, m_field_FourierSpace[1], m_field[1], FFTW_ESTIMATE);
  fftw_execute(complex2real);
  fftw_destroy_plan(complex2real);
  
  complex2real = fftw_plan_dft_c2r_3d(m_nX, m_nY, m_nZ, m_field_FourierSpace[2], m_field[2], FFTW_ESTIMATE);
  fftw_execute(complex2real);
  fftw_destroy_plan(complex2real);

}


// ============================================================================


void cbl::data::VectorField3D::set_VectorField (const vector<double> value, const int i, const int j, const int k, const bool add)
{
  m_field[0][inds_to_index(i,j,k)] = (add) ? m_field[0][inds_to_index(i,j,k)]+value[0]: value[0];
  m_field[1][inds_to_index(i,j,k)] = (add) ? m_field[1][inds_to_index(i,j,k)]+value[1]: value[1];
  m_field[2][inds_to_index(i,j,k)] = (add) ? m_field[2][inds_to_index(i,j,k)]+value[2]: value[2];
}


// ============================================================================


void cbl::data::VectorField3D::set_VectorField_FourierSpace_real (const std::vector<double> value, const int i, const int j, const int k, const bool add)
{
  m_field_FourierSpace[0][inds_to_index_Fourier(i,j,k)][0] = (add) ? m_field_FourierSpace[0][inds_to_index(i,j,k)][0]+value[0] : value[0];
  m_field_FourierSpace[1][inds_to_index_Fourier(i,j,k)][0] = (add) ? m_field_FourierSpace[1][inds_to_index(i,j,k)][0]+value[1] : value[1];
  m_field_FourierSpace[2][inds_to_index_Fourier(i,j,k)][0] = (add) ? m_field_FourierSpace[2][inds_to_index(i,j,k)][0]+value[2] : value[2];
}


// ============================================================================


void cbl::data::VectorField3D::set_VectorField_FourierSpace_complex(const vector<double> value, const int i, const int j, const int k, const bool add)
{
  m_field_FourierSpace[0][inds_to_index_Fourier(i,j,k)][1] = (add) ? m_field_FourierSpace[0][inds_to_index(i,j,k)][1]+value[0] : value[0];
  m_field_FourierSpace[1][inds_to_index_Fourier(i,j,k)][1] = (add) ? m_field_FourierSpace[1][inds_to_index(i,j,k)][1]+value[1] : value[1];
  m_field_FourierSpace[2][inds_to_index_Fourier(i,j,k)][1] = (add) ? m_field_FourierSpace[2][inds_to_index(i,j,k)][1]+value[2] : value[2];
}


// ============================================================================


vector<double> cbl::data::VectorField3D::VectorField(const int i, const int j, const int k) const 
{
  return {m_field[0][inds_to_index(i,j,k)], m_field[1][inds_to_index(i,j,k)], m_field[2][inds_to_index(i,j,k)]};
}


// ============================================================================


vector<double> cbl::data::VectorField3D::VectorField_FourierSpace_real(const int i, const int j, const int k) const 
{
  return {m_field_FourierSpace[0][inds_to_index_Fourier(i,j,k)][0], m_field_FourierSpace[1][inds_to_index_Fourier(i,j,k)][0], m_field_FourierSpace[2][inds_to_index_Fourier(i,j,k)][0]};
}


// ============================================================================


vector<double> cbl::data::VectorField3D::VectorField_FourierSpace_complex(const int i, const int j, const int k) const 
{
  return {m_field_FourierSpace[0][inds_to_index_Fourier(i,j,k)][1], m_field_FourierSpace[1][inds_to_index_Fourier(i,j,k)][1], m_field_FourierSpace[2][inds_to_index_Fourier(i,j,k)][1]};
}


// ============================================================================


std::vector<double> cbl::data::VectorField3D::VectorField (const std::vector<double> pos) const 
{
  int i = min(int((pos[0]-m_MinX)/m_deltaX), m_nX-1);
  int j = min(int((pos[1]-m_MinY)/m_deltaY), m_nY-1);
  int k = min(int((pos[2]-m_MinZ)/m_deltaZ), m_nZ-1);

  return {m_field[0][inds_to_index(i,j,k)], m_field[1][inds_to_index(i,j,k)], m_field[2][inds_to_index(i,j,k)]};
}


// ============================================================================


void cbl::data::VectorField3D::reset () 
{
  for(int i=0;i<m_nCells;i++){
    m_field[0][i]=0;
    m_field[1][i]=0;
    m_field[2][i]=0;
  }
  for(int i=0;i<m_nCells_Fourier;i++){
    m_field_FourierSpace[0][i][0] = 0;
    m_field_FourierSpace[0][i][1] = 0;
    m_field_FourierSpace[1][i][0] = 0;
    m_field_FourierSpace[1][i][1] = 0;
    m_field_FourierSpace[2][i][0] = 0;
    m_field_FourierSpace[2][i][1] = 0;
  }
}
