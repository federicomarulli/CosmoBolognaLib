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
 *  @file GlobalFunc/Reconstruction.cpp
 *
 *  @brief Functions to compute displacement for the reconstructed
 *  density field of a collection of points
 *
 *  This file contains the implementation of a set of functions to
 *  perform density field reconstruction
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#include "GlobalFunc.h"

using namespace cosmobl;


// ============================================================================


void cosmobl::reconstruction_fourier_space (const catalogue::Catalogue data, const catalogue::Catalogue random, const bool random_RSD, const cosmology::Cosmology cosmology, const double redshift, const double bias, const double cell_size, const double smoothing_radius, const int interpolation_type)
{
  double ff = cosmology.linear_growth_rate(redshift);
  double beta = ff/bias;
  
  //double rsd_term = 3*ff/(7+3*ff);
  //double rsd_term = (ff-beta)/(1+beta);
  //double rsd_term= ff/(ff+1);

  double rsd_term = 3*beta/(7+3*beta);
  
data::ScalarField3D density_field = data.density_field(cell_size, random, interpolation_type, smoothing_radius);

  density_field.FourierTransformField();

  data::VectorField3D displacement_field(density_field.nx(), density_field.ny(), density_field.nz(), density_field.MinX(), density_field.MaxX(), density_field.MinY(), density_field.MaxY(), density_field.MinZ(), density_field.MaxZ()); 

  data::VectorField3D displacement_field_RSD(density_field.nx(), density_field.ny(), density_field.nz(), density_field.MinX(), density_field.MaxX(), density_field.MinY(), density_field.MaxY(), density_field.MinZ(), density_field.MaxZ()); 

  
  // convert the grid to the displacement field (in Fourier space)

  double xfact = 2.*par::pi/(density_field.deltaX()*density_field.nx());
  double yfact = 2.*par::pi/(density_field.deltaY()*density_field.ny());
  double zfact = 2.*par::pi/(density_field.deltaZ()*density_field.nz());

  for (int i=0; i<density_field.nx(); i++) {
    double kx = (i<=density_field.nx()/2) ? xfact*i : -(density_field.nx()-i)*xfact;
    for (int j=0; j<density_field.ny(); j++) {
      double ky =(j<=density_field.ny()/2) ? yfact*j : -(density_field.ny()-j)*yfact;
      for (int k=0; k<density_field.nzFourier(); k++) {
	double kz = zfact*k;
	double kk = pow(kx*kx+ky*ky+kz*kz, 0.5);

	vector<double> gridK = {density_field.ScalarField_FourierSpace_real (i, j, k)/bias, density_field.ScalarField_FourierSpace_complex (i, j, k)/bias };

	if (kk!=0){
	  displacement_field.set_VectorField_FourierSpace_real({kx*gridK[1]/(kk*kk), ky*gridK[1]/(kk*kk), kz*gridK[1]/(kk*kk)} ,i ,j ,k);

	  displacement_field.set_VectorField_FourierSpace_complex({-kx*gridK[0]/(kk*kk), -ky*gridK[0]/(kk*kk)/bias, -kz*gridK[0]/(kk*kk)} ,i , j, k);
	}

      }
    }
  }
  
  displacement_field.FourierAntiTransformField();

  
  // apply the RSD correction (in Fourier space)
  
  for (int i=0; i<density_field.nx(); i++)
    for (int j=0; j<density_field.ny(); j++)
      for (int k=0; k<density_field.nz(); k++) {

	double rx = density_field.XX(i), ry = density_field.YY(j), rz = density_field.ZZ(k); 
	double rr = (sqrt(rx*rx+ry*ry+rz*rz)==0) ? 1 : sqrt(rx*rx+ry*ry+rz*rz);

	vector<double> displ = displacement_field.VectorField(i, j, k);
	double scalar_product = (rr==0) ? 0. : (displ[0]*rx+displ[1]*ry+displ[2]*rz)/rr;

	double displ_x_rsd = displ[0]+rsd_term*scalar_product*rx/rr; 
	double displ_y_rsd = displ[1]+rsd_term*scalar_product*ry/rr; 
	double displ_z_rsd = displ[2]+rsd_term*scalar_product*rz/rr; 

	displacement_field_RSD.set_VectorField({displ_x_rsd, displ_y_rsd, displ_z_rsd}, i, j, k);
      }

  
  // set the displacement

  for (size_t i=0; i<data.nObjects(); i++) {
    int i1 = min(int((data.xx(i)-density_field.MinX())/density_field.deltaX()), density_field.nx()-1);
    int j1 = min(int((data.yy(i)-density_field.MinY())/density_field.deltaY()), density_field.ny()-1);
    int k1 = min(int((data.zz(i)-density_field.MinZ())/density_field.deltaZ()), density_field.nz()-1);

    vector<double> displacement = displacement_field_RSD.VectorField(i1, j1, k1);

    data.catalogue_object(i)->set_x_displacement(displacement[0]);
    data.catalogue_object(i)->set_y_displacement(displacement[1]);
    data.catalogue_object(i)->set_z_displacement(displacement[2]);
  }

  for (size_t i=0; i<random.nObjects(); i++) {
    int i1 = min(int((random.xx(i)-density_field.MinX())/density_field.deltaX()), density_field.nx()-1);
    int j1 = min(int((random.yy(i)-density_field.MinY())/density_field.deltaY()), density_field.ny()-1);
    int k1 = min(int((random.zz(i)-density_field.MinZ())/density_field.deltaZ()), density_field.nz()-1);

    vector<double> displacement = (random_RSD) ? displacement_field_RSD.VectorField(i1, j1, k1) : displacement_field.VectorField(i1, j1, k1);

    random.catalogue_object(i)->set_x_displacement(displacement[0]);
    random.catalogue_object(i)->set_y_displacement(displacement[1]);
    random.catalogue_object(i)->set_z_displacement(displacement[2]);
  }
}


// ============================================================================


cosmobl::catalogue::Catalogue cosmobl::displaced_catalogue (const catalogue::Catalogue input_catalogue) 
{
  catalogue::Catalogue displaced_cat(input_catalogue);

  // apply the displacements

  for (size_t i=0; i<displaced_cat.nObjects(); i++) {
    
    double xx = displaced_cat.xx(i);
    double yy = displaced_cat.yy(i);
    double zz = displaced_cat.zz(i);

    double x_displ = displaced_cat.x_displacement(i);
    double y_displ = displaced_cat.y_displacement(i);
    double z_displ = displaced_cat.z_displacement(i);

    displaced_cat.catalogue_object(i)->set_xx(xx+x_displ);
    displaced_cat.catalogue_object(i)->set_yy(yy+y_displ);
    displaced_cat.catalogue_object(i)->set_zz(zz+z_displ);
  }

  displaced_cat.computePolarCoordinates();

  return displaced_cat;
}
