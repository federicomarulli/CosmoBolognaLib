/*******************************************************************
 *  Copyright (C) 2010 by Federico Marulli                         *
 *  federico.marulli3@unibo.it                                     *
 *                                                                 *
 *  This program is free software; you can redistribute it and/or  *
 *  modify it under the terms of the GNU General Public License as *
 *  published by the Free Software Foundation; either version 2 of *
 *  the License, or (at your option) any later version.            *
 *                                                                 *
 *  This program is distributed in the hope that it will be useful,*
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the  *
 *  GNU General Public License for more details.                   *
 *                                                                 *
 *  You should have received a copy of the GNU General Public      *
 *  License along with this program; if not, write to the Free     *
 *  Software Foundation, Inc.,                                     *
 *  59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.      *
 *******************************************************************/

/**
 *  @file CosmoBolognaLib/Func/Func.cpp
 *
 *  @brief Useful generic functions
 *
 *  This file contains the implementation of a large set of useful
 *  functions of wide use
 *
 *  @author Federico Marulli
 *
 *  @author federico.marulli3@unbo.it
 */

#include "Func.h"

using namespace std;

using namespace cbl;
using namespace glob;


// ============================================================================


double cbl::closest_probability (double xx, std::shared_ptr<void> pp, std::vector<double> par)
{
  (void)par;
  shared_ptr<cbl::glob::STR_closest_probability> pars = static_pointer_cast<cbl::glob::STR_closest_probability>(pp);
  return pars->weights[(cbl::index_closest(xx, pars->values))];
}


// ============================================================================


double cbl::distribution_probability (double xx, std::shared_ptr<void> pp, std::vector<double> par)
{
  (void)par;
  shared_ptr<cbl::glob::STR_distribution_probability> pars = static_pointer_cast<cbl::glob::STR_distribution_probability>(pp);
  return pars->func->operator()(xx);
}


// ============================================================================


double cbl::Filter (const double r, const double rc)
{
  double x = pow(r/rc, 3);
  return pow(2*x*(1.-x), 2)*(0.5-x)*pow(rc, -3);
}


// ============================================================================================


double cbl::degrees (const double angle, const CoordinateUnits inputUnits)
{
  if (inputUnits==CoordinateUnits::_radians_) 
    return angle*180./par::pi;
      
  else if (inputUnits==CoordinateUnits::_arcseconds_)
    return angle/3600.;
      
  else if (inputUnits==CoordinateUnits::_arcminutes_)
    return angle/60.;

  else if (inputUnits==CoordinateUnits::_degrees_)
    return angle;
  
  else 
    { ErrorCBL("Error in cbl::degrees() in Func.cpp: inputUnits type not allowed!"); return 0; }
}


// ============================================================================================


double cbl::radians (const double angle, const CoordinateUnits inputUnits)
{
  if (inputUnits==CoordinateUnits::_degrees_)
    return angle/180.*par::pi;
  
  else if (inputUnits==CoordinateUnits::_arcseconds_)
    return angle/180.*par::pi/3600.;
  
  else if (inputUnits==CoordinateUnits::_arcminutes_)
    return angle/180.*par::pi/60.;

  else if (inputUnits==CoordinateUnits::_radians_)
    return angle;
  
  else 
    { ErrorCBL("Error in cbl::radians() in Func.cpp: inputUnits type not allowed!"); return 0; }
}


// ============================================================================================


double cbl::arcseconds (const double angle, const CoordinateUnits inputUnits)
{
  if (inputUnits==CoordinateUnits::_radians_)
    return angle*180./par::pi*3600.;
      
  else if (inputUnits==CoordinateUnits::_degrees_)
    return angle*3600.;
  
  else if (inputUnits==CoordinateUnits::_arcminutes_)
    return angle*60.;

  else if (inputUnits==CoordinateUnits::_arcseconds_)
    return angle;
      
  else 
    { ErrorCBL("Error in cbl::arcseconds() in Func.cpp: inputUnits type not allowed!"); return 0; }
}


// ============================================================================================


double cbl::arcminutes (const double angle, const CoordinateUnits inputUnits)
{
  if (inputUnits==CoordinateUnits::_radians_)
    return angle*180./par::pi*60.;
      
  else if (inputUnits==CoordinateUnits::_degrees_)
    return angle*60.;
  
  else if (inputUnits==CoordinateUnits::_arcseconds_)
    return angle/60.;

  else if (inputUnits==CoordinateUnits::_arcminutes_)
    return angle;

  else 
    { ErrorCBL("Error in cbl::arcminutes() in Func.cpp: inputUnits type not allowed!"); return 0; }
}


// ============================================================================================


double cbl::converted_angle (const double angle, const CoordinateUnits inputUnits, const CoordinateUnits outputUnits)
{
  if (outputUnits==CoordinateUnits::_radians_)
    return radians(angle, inputUnits);
      
  else if (outputUnits==CoordinateUnits::_degrees_)
    return degrees(angle, inputUnits);
  
  else if (outputUnits==CoordinateUnits::_arcseconds_)
    return arcseconds(angle, inputUnits);

  else if (outputUnits==CoordinateUnits::_arcminutes_)
    return arcminutes(angle, inputUnits);
  
  else 
    { ErrorCBL("Error in cbl::converted_angle() in Func.cpp: outputUnits type not allowed!"); return 0; }
}


// ============================================================================================


void cbl::polar_coord (const double XX, const double YY, const double ZZ, double &ra, double &dec, double &dd) 
{   
  dd = sqrt(XX*XX+YY*YY+ZZ*ZZ);
  ra = atan2(XX,YY);
  dec = asin(ZZ/dd);  
}

void cbl::cartesian_coord (const double ra, const double dec, const double dd, double &XX, double &YY, double &ZZ) 
{
  XX = dd*cos(dec)*sin(ra);
  YY = dd*cos(dec)*cos(ra);
  ZZ = dd*sin(dec);
}

void cbl::polar_coord (const std::vector<double> XX, const std::vector<double> YY, const std::vector<double> ZZ, std::vector<double> &ra, std::vector<double> &dec, std::vector<double> &dd) 
{     
  for (size_t i=0; i<XX.size(); i++) {
    dd[i] = sqrt(XX[i]*XX[i]+YY[i]*YY[i]+ZZ[i]*ZZ[i]);
    ra[i] = atan2(XX[i],YY[i]);
    dec[i] = asin(ZZ[i]/dd[i]);
  }  
}

void cbl::cartesian_coord (const std::vector<double> ra, const std::vector<double> dec, const std::vector<double> dd, std::vector<double> &XX, std::vector<double> &YY, std::vector<double> &ZZ) 
{
  for (size_t i=0; i<XX.size(); i++) {
    XX[i] = dd[i]*cos(dec[i])*sin(ra[i]);
    YY[i] = dd[i]*cos(dec[i])*cos(ra[i]);
    ZZ[i] = dd[i]*sin(dec[i]);
  }
}


// ============================================================================================


double cbl::Euclidean_distance (const double x1, const double x2, const double y1, const double y2, const double z1, const double z2)
{
  return sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2));
}


// ============================================================================================


double cbl::perpendicular_distance (const double ra1, const double ra2, const double dec1, const double dec2, const double d1, const double d2)
{
  double costheta = sin(dec1)*sin(dec2)+cos(dec1)*cos(dec2)*(cos(ra1)*cos(ra2)+sin(ra1)*sin(ra2));
  double theta = 0.;
  if (fabs(costheta)<1.-1.e-30) theta = acos(costheta);
  else if (costheta>=1.-1.e-30) theta = 0.;
  else theta = par::pi;
             
  double rp = (d1+d2)*tan(theta*0.5);
  rp *= 4.*d1*d2/((d1+d2)*(d1+d2));

  return rp;
}


// ============================================================================================


double cbl::angular_distance (const double x1, const double x2, const double y1, const double y2, const double z1, const double z2)
{
  return 2.*asin(0.5*Euclidean_distance(x1, x2, y1, y2, z1, z2));
}


// ============================================================================================


double cbl::haversine_distance (const double ra1, const double ra2, const double dec1, const double dec2)
{
  return 2*asin(sqrt(pow(sin((dec2-dec1)*0.5),2.)+cos(dec1)*cos(dec2)*pow(sin((ra2-ra1)*0.5),2.)));
}


// ============================================================================================


double cbl::MC_Int (double func(const double), const double x1, const double x2, const int seed) 
{
  int step = 100000;
  double delta_x = (x2-x1)/step;
  double xx = x1;
  vector<double> ff(step+1);
  for (int i=0; i<step+1; i++) {ff[i] = func(xx); xx += delta_x;}

  vector<double>::iterator fmin = min_element (ff.begin(),ff.end()),
    fmax = max_element (ff.begin(),ff.end());
  double f1 = *fmin;
  double f2 = *fmax;
  
  f1 = (f1>0) ? f1*0.5 : -fabs(f1)*2.;
  f2 *= 2.; 

  random::UniformRandomNumbers ran(0., 1., seed);
  
  double xt, yt, INT;
  int sub = 0, subn = 0, numTOT = 10000000;

  if (f1>0) {
    for (int i=0; i<numTOT; i++) {
      xt = ran()*(x2-x1)+x1;
      yt = ran()*(f2-f1);
      if (yt<func(xt)) sub ++;
    }
    INT = double(sub)/double(numTOT)*(x2-x1)*(f2-f1);
  }
  else {
    for (int i=0; i<numTOT; i++) {
      xt = ran()*(x2-x1)+x1;
      yt = ran()*f2;
      if (yt<func(xt)) sub ++;
    }
    for (int i=0; i<numTOT; i++) {
      xt = ran()*(x2-x1)+x1;
      yt = ran()*f1;
      if (yt>func(xt)) subn ++;
    }
    INT = (double(sub)/double(numTOT)*(x2-x1)*f2)-(double(subn)/double(numTOT)*(x2-x1)*fabs(f1));
  }
  
  return INT;
}


// ============================================================================================


double cbl::MC_Int (double func(const double, const double AA), const double AA, const double x1, const double x2, const int seed) 
{
  int step = 100000;
  double delta_x = (x2-x1)/step;
  double xx = x1;
  vector<double> ff(step+1);
  for (int i=0; i<step+1; i++) {ff[i] = func(xx,AA); xx += delta_x;}

  vector<double>::iterator fmin = min_element (ff.begin(),ff.end()),
    fmax = max_element (ff.begin(),ff.end());
  double f1 = *fmin;
  double f2 = *fmax;
  
  f1 = (f1>0) ? f1*0.5 : -fabs(f1)*2.;
  f2 *= 2.; 

  random::UniformRandomNumbers ran(0., 1., seed);

  double xt, yt, INT;
  int sub = 0, subn = 0, numTOT = 10;

  if (f1>0) {
    for (int i=0; i<numTOT; i++) {
      xt = ran()*(x2-x1)+x1;
      yt = ran()*(f2-f1);
      if (yt<func(xt, AA)) sub ++;
    }
    INT = double(sub)/double(numTOT)*(x2-x1)*(f2-f1);
  }
  else {
    for (int i=0; i<numTOT; i++) {
      xt = ran()*(x2-x1)+x1;
      yt = ran()*f2;
      if (yt<func(xt,AA)) sub ++;
    }
    for (int i=0; i<numTOT; i++) {
      xt = ran()*(x2-x1)+x1;
      yt = ran()*f1;
      if (yt>func(xt,AA)) subn ++;
    }
    INT = (double(sub)/double(numTOT)*(x2-x1)*f2)-(double(subn)/double(numTOT)*(x2-x1)*fabs(f1));
  }
  
  return INT;
}


// ============================================================================================


double cbl::MC_Int (double func(const double, const double AA, const double BB, const double CC, const double DD, const double EE), const double AA, const double BB, const double CC, const double DD, const double EE, const double x1, const double x2, const int seed) 
{
  int step = 100000;
  double delta_x = (x2-x1)/step;
  double xx = x1;
  vector<double> ff(step+1);
  for (int i=0; i<step+1; i++) {ff[i] = func(xx,AA,BB,CC,DD,EE); xx += delta_x;}
 
  vector<double>::iterator fmin = min_element (ff.begin(),ff.end()),
    fmax = max_element (ff.begin(),ff.end());
  double f1 = *fmin;
  double f2 = *fmax;
  
  f1 = (f1>0) ? f1*0.5 : -fabs(f1)*2.;
  f2 *= 2.; 

  random::UniformRandomNumbers ran(0., 1., seed);
  
  double xt, yt, INT;
  int sub = 0, subn = 0, numTOT = 100000;

  if (f1>0) {
    for (int i=0; i<numTOT; i++) {
      xt = ran()*(x2-x1)+x1;
      yt = ran()*(f2-f1);
      if (yt<func(xt,AA,BB,CC,DD,EE)) sub ++;
    }
    INT = double(sub)/double(numTOT)*(x2-x1)*(f2-f1);
  }
  else {
    for (int i=0; i<numTOT; i++) {
      xt = ran()*(x2-x1)+x1;
      yt = ran()*f2;
      if (yt<func(xt,AA,BB,CC,DD,EE)) sub ++;
    }
    for (int i=0; i<numTOT; i++) {
      xt = ran()*(x2-x1)+x1;
      yt = ran()*f1;
      if (yt>func(xt,AA,BB,CC,DD,EE)) subn ++;
    }
    INT = (double(sub)/double(numTOT)*(x2-x1)*f2)-(double(subn)/double(numTOT)*(x2-x1)*fabs(f1));
  }
  
  return INT;
}


// ============================================================================


double cbl::interpolated (const double _xx, const std::vector<double> xx, const std::vector<double> yy, const std::string type)
{
  if (xx.size()!=yy.size() || xx.size()<2) 
    ErrorCBL("Error in cbl::interpolated() in Func.cpp: "+conv(xx.size(), par::fINT)+"!");

  size_t size = xx.size();
  
  if (_xx<xx[0]) // perform a linear extrapolation
    return yy[0]+(_xx-xx[0])/(xx[1]-xx[0])*(yy[1]-yy[0]);

  else if (_xx>xx[size-1])
    return yy[size-2]+(_xx-xx[size-2])/(xx[size-1]-xx[size-2])*(yy[size-1]-yy[size-2]);
  

  gsl_interp_accel *acc = gsl_interp_accel_alloc();
  const gsl_interp_type *TT;
  
  if (type=="Linear") 
    TT = gsl_interp_linear;

  else if (type=="Poly") 
    TT = gsl_interp_polynomial;

  else if (type=="Spline") 
    TT = gsl_interp_cspline;

  else if (type=="Spline_periodic") 
    TT = gsl_interp_cspline_periodic;
  
  else if (type=="Akima") 
    TT = gsl_interp_akima;
  
  else if (type=="Akima_periodic") 
    TT = gsl_interp_akima_periodic;
  
  else if (type=="Steffen") 
    TT = gsl_interp_steffen;
  
  else 
    ErrorCBL("Error in cbl::interpolated() in Func.cpp: the value of string 'type' is not permitted!");

  gsl_interp *interp = gsl_interp_alloc(TT, size);
  gsl_interp_init(interp, xx.data(), yy.data(), size);
  
  double _yy;
  interp->type->eval(interp->state, xx.data(), yy.data(), interp->size, _xx, acc, &_yy);
  
  gsl_interp_free(interp);
  gsl_interp_accel_free(acc);

  return _yy;
}


// ============================================================================


double cbl::interpolated_2D (const double _x1, const double _x2, const std::vector<double> x1, const std::vector<double> x2, const std::vector<std::vector<double>> yy, const std::string type)
{
  bool extr = false;

  const size_t size_x1 = x1.size();
  const size_t size_x2 = x2.size();
  double *ydata = new double[size_x1*size_x2];
  
  if ((_x1>Max(x1) || Min(x1)>_x1) ||(_x2>Max(x2) || Min(x2)>_x2))
    extr = 1;

  gsl_interp_accel *x1acc = gsl_interp_accel_alloc();
  gsl_interp_accel *x2acc = gsl_interp_accel_alloc();
  const gsl_interp2d_type *TT = gsl_interp2d_bilinear;

  if (type=="Linear") 
    TT = gsl_interp2d_bilinear;

  else if (type=="Cubic") 
    TT = gsl_interp2d_bicubic;

  for (size_t i=0; i<size_x1; i++)
    for (size_t j=0; j<size_x2; j++)
      ydata[i+j*size_x1] = yy[i][j];

  gsl_interp2d *interp = gsl_interp2d_alloc(TT, size_x1, size_x2);
  gsl_interp2d_init (interp, x1.data(), x2.data(), ydata, size_x1, size_x2);

  double val;
  if (extr)
    val= gsl_interp2d_eval_extrap(interp, x1.data(), x2.data() , ydata, _x1, _x2, x1acc, x2acc);
  else
    val= gsl_interp2d_eval(interp, x1.data(), x2.data() , ydata, _x1, _x2, x1acc, x2acc);

  gsl_interp2d_free(interp);
  gsl_interp_accel_free(x1acc);
  gsl_interp_accel_free(x2acc);
  free(ydata);

  return val;
}


// ============================================================================


void cbl::read_vector (const std::string file_vector, std::vector<double> &xx, std::vector<double> &vec, const std::vector<int> col)
{
  vector<int> cols = {0, 1};
  cols = (col.size() != 2) ? cols : col; 
  const size_t max_col = Max(cols);

  xx.erase(xx.begin(), xx.end());
  vec.erase(vec.begin(), vec.end());

  ifstream fin(file_vector.c_str()); checkIO(fin, file_vector);

  string line;

  while (getline(fin, line)) {
    stringstream ss(line);
    vector<double> num; double NN = -1.e30;
    while (ss>>NN) num.push_back(NN);
    if (num.size()>=max_col) {
      xx.push_back(num[cols[0]]);
      vec.push_back(num[cols[1]]);
    }
  }

  fin.clear(); fin.close();
}


// ============================================================================


void cbl::read_matrix (const std::string file_matrix, std::vector<double> &xx, std::vector<double> &yy, std::vector<std::vector<double>> &matrix, const std::vector<int> col)
{
  vector<int> cols = {0, 1, 2};
  cols = (col.size() != 3) ? cols : col; 
  const size_t max_col = Max(cols);

  matrix.erase(matrix.begin(), matrix.end());

  ifstream fin(file_matrix.c_str()); checkIO(fin, file_matrix);

  vector<double> vv;
  matrix.push_back(vv);
  string line; size_t i = 0;

  vector<double> _xx, _yy;
  while (getline(fin, line)) {
    stringstream ss(line);
    vector<double> num; double NN = -1.e30;
    while (ss>>NN) num.push_back(NN);
    if (num.size()>=max_col) {
      _xx.push_back(num[cols[0]]);
      _yy.push_back(num[cols[1]]);
      matrix[i].push_back(num[cols[2]]);
    }
    else {i++; matrix.push_back(vv);}
  }

  xx = different_elements(_xx);
  yy = different_elements(_yy);

  fin.clear(); fin.close();

}


// ============================================================================

double cbl::determinant_matrix (const std::vector<std::vector<double>> mat)
{
  int n = mat.size();
  int s;
  gsl_matrix *mm = gsl_matrix_alloc(n, n);
  gsl_permutation *perm = gsl_permutation_alloc(n);

  for (int i=0; i<n; i++)
    for (int j=0; j<n; j++)
      gsl_matrix_set(mm, i, j, mat[i][j]);

  // make the LU decomposition of the matrix mm
  gsl_linalg_LU_decomp(mm, perm, &s);
  double det = gsl_linalg_LU_det(mm, s);

  gsl_matrix_free(mm);
  gsl_permutation_free(perm);
  return det;
}

// ============================================================================


void cbl::invert_matrix (const std::vector<std::vector<double>> mat, std::vector<std::vector<double>> &mat_inv, const double prec)
{
  int n = mat.size();
  int s;
  if (n==0)
    ErrorCBL("Error in cbl::invert_matrix() in Func.cpp. 0 size for the input matrix");

  mat_inv.erase(mat_inv.begin(), mat_inv.end());
  mat_inv = mat;

  gsl_matrix *mm = gsl_matrix_alloc(n, n);
  gsl_matrix *im = gsl_matrix_alloc(n, n);
  gsl_permutation *perm = gsl_permutation_alloc(n);

  for (int i=0; i<n; i++)
    for (int j=0; j<n; j++)
      gsl_matrix_set(mm, i, j, mat[i][j]);

  // make the LU decomposition of the matrix mm
  gsl_linalg_LU_decomp(mm, perm, &s);

  // invert the matrix mm
  gsl_linalg_LU_invert(mm, perm, im);

  for (int i=0; i<n; i++)
    for (int j=0; j<n; j++)
      mat_inv[i][j] = gsl_matrix_get(im, i, j);

  for (int i=0; i<n; i++) {
    for (int j=0; j<n; j++) {
      double fact = (i==j) ? 1 : 0;
      double prod = 0;
      for (int el=0; el<n; el++)
	prod += mat[i][el]*mat_inv[el][j];
      
      if (fabs(fact - prod) > prec)  
	WarningMsg("Exceeded precision for element "+conv(i, par::fINT)+" "+conv(j, par::fINT)+"; "+conv(fact, par::fDP4)+" "+conv(prod, par::fDP4));
    }
  }
  
  gsl_matrix_free(mm);
  gsl_matrix_free(im);
  gsl_permutation_free(perm);
}


// ============================================================================


void cbl::invert_matrix (const std::vector<std::vector<double>> mat, std::vector<std::vector<double>> &mat_inv, const int i1, const int i2, const double prec)
{
  int n = i2-i1;
  int s;
  if (n==0)
    ErrorCBL("Error in cbl::invert_matrix() in Func.cpp: 0 size for the input matrix");

  mat_inv.erase(mat_inv.begin(),mat_inv.end());
  mat_inv = mat;

  gsl_matrix *mm = gsl_matrix_alloc (n, n);
  gsl_matrix *im = gsl_matrix_alloc (n, n);
  gsl_permutation * perm = gsl_permutation_alloc (n);

  for (int i=i1; i<i2; i++)
    for (int j=i1; j<i2; j++)
      gsl_matrix_set(mm,i-i1,j-i1,mat[i][j]);

  // Make LU decomposition of matrix m
  gsl_linalg_LU_decomp (mm, perm, &s);

  // Invert the matrix m
  gsl_linalg_LU_invert (mm, perm, im);

  for (int i=0; i<n; i++) {
    for (int j=0; j<n; j++) {
      double fact = (i==j) ? 1 : 0;
      double prod = 0;
      for (int el=0; el<n; el++)
	prod += mat[i+i1][el+i1]*gsl_matrix_get(im, el, j); 
    
      if (fabs(fact - prod) > prec)  
	WarningMsg("Exceeded precision for element "+conv(i,par::fINT)+" "+conv(j,par::fINT)+"; "+conv(fact,par::fDP4)+" "+conv(prod,par::fDP4));
    }
  }

  for (size_t i=0; i<mat.size(); i++) {
    for (size_t j=0; j<mat[i].size(); j++) {
      if (int(i)<i2 && int(i)>=i1 && int(j)<i2 && int(j)>=i1)
	mat_inv[i][j] = gsl_matrix_get(im,i-i1,j-i1);
      else
	mat_inv[i][j] = 0;
    }
  }
}


// ============================================================================


void cbl::covariance_matrix (const std::vector<std::vector<double>> mat, std::vector<std::vector<double>> &cov, const bool JK) 
{  
  cov.erase(cov.begin(), cov.end());
  vector<double> vv (mat[0].size(), 0.);
  for (size_t i=0; i<mat[0].size(); i++) cov.push_back(vv);
  

  // measure the mean values at each bin

  vector<double> mean;

  for (size_t i=0; i<mat[0].size(); i++) { // loop on the i-th bin 
    
    vector<double> vect_temp;
    for (size_t j=0; j<mat.size(); j++) // loop on the j-th realisation
      vect_temp.push_back(mat[j][i]);

    if (vect_temp.size()>2) 
      mean.push_back(Average(vect_temp));
    else 
      mean.push_back(-1.e30);
    
  }
  
  
  // compute the elements of the covariance matrix
  
  for (size_t i=0; i<cov.size(); i++) 
    for (size_t j=0; j<cov.size(); j++) {
      int nm = 0;
      for (size_t k=0; k<mat.size(); k++) {
	if (mean[i]>-1.e30) {
	  cov[i][j] += (mean[i]-mat[k][i])*(mean[j]-mat[k][j]);
	  nm ++;
	}
      }
      if (nm>1) cov[i][j] = (JK) ? double(nm-1)/(nm)*cov[i][j] : cov[i][j]/(nm-1);
    }

}


// ============================================================================


void cbl::covariance_matrix (const std::vector<std::string> file, std::vector<double> &rad, std::vector<double> &mean, std::vector<std::vector<double>> &cov, const bool JK) 
{  
  vector<vector<double>> xi_mocks;

  string line; double r,xi;
  for (size_t i=0; i<file.size(); i++) {
    rad.erase(rad.begin(),rad.end());
    
    ifstream fin(file[i].c_str()); checkIO(fin, file[i]);
    getline(fin, line);

    vector<double> vv;
    while (getline(fin,line)) {
      stringstream ss(line);
      ss >> r; ss >> xi;
      rad.push_back(r);
      vv.push_back(xi);
    }
    fin.clear(); fin.close();

    xi_mocks.push_back(vv);
  }

  covariance_matrix(xi_mocks, cov, JK);

  mean.erase(mean.begin(),mean.end());
  
  for (size_t i=0; i<rad.size(); i++) {
    vector<double> vv;
    
    for (size_t j=0; j<xi_mocks.size(); j++)
      if (xi_mocks[j][i]>-1.e30)
	vv.push_back(xi_mocks[j][i]);

    mean.push_back(Average(vv));
  }
  
}

// ============================================================================


void cbl::covariance_matrix (const std::vector<std::string> file, const std::string covariance_matrix_file, const bool JK) 
{  

  vector<double> rad, mean;
  vector<vector<double>> cov;
  covariance_matrix(file, rad, mean, cov, JK);

  ofstream fout(covariance_matrix_file.c_str()); checkIO(fout, covariance_matrix_file);

  for (size_t i=0; i<rad.size(); i++) {
    for (size_t j=0; j<rad.size(); j++)
      fout << rad[i] << " " << rad[j] << " " << cov[i][j] << " " <<  cov[i][j]/sqrt(cov[i][i]*cov[j][j]) << endl;
    fout << endl;
  }

  fout.clear(); fout.close();
}


// ============================================================================

 
double cbl::Average (const std::vector<double> vect) 
{
  if (vect.size()==0) ErrorCBL("Error in cbl::Average() in Func.cpp: 0 size for the input vector");

  double aver = 0., NT = 0.;
  
#pragma omp parallel num_threads(omp_get_max_threads())
  {
    
    double averT = 0., nT = 0.;
	
#pragma omp for schedule(static, 2)
    for (size_t i=0; i<vect.size(); ++i) 
      averT += 1./(++nT)*(vect[i]-averT);
    
#pragma omp critical
    aver += ((NT+=nT)>0) ? nT/NT*(averT-aver) : 0.;
    
  }
  
  return aver;
}


// ============================================================================


double cbl::Average (const std::vector<double> vect, const std::vector<double> weight) 
{
  if (vect.size()==0 || vect.size()!=weight.size())
    ErrorCBL("Error in cbl::Average() in Func.cpp: 0 size for the input vector or vect.size()!=weight.size()");

  double aver = 0., WeightTOT = 0.;
  
#pragma omp parallel num_threads(omp_get_max_threads())
  {

    double averT = 0., WeightT = 0.;
    
#pragma omp for schedule(static, 2)
    for (size_t i=0; i<vect.size(); ++i) 
      averT += weight[i]/(WeightT+=weight[i])*(vect[i]-averT);
    
#pragma omp critical
    aver += ((WeightTOT+=WeightT)>0) ? WeightT/WeightTOT*(averT-aver) : 0.;

  }
  
  return aver;  
}


// ============================================================================


double cbl::Sigma (const std::vector<double> vect) 
{
  if (vect.size()==0) ErrorCBL("Error in cbl::Sigma() in Func.cpp: 0 size for the input vector");
  
  double aver_n1 = 0., aver_n = 0., Sn = 0., sigma = 0., NT = 0.;

#pragma omp parallel num_threads(omp_get_max_threads())
  {
    
    double aver_n1T = 0., aver_nT = 0., SnT = 0., nT = 0.;

#pragma omp for schedule(static, 2)
    for (size_t i=0; i<vect.size(); ++i) {
      aver_n1T = aver_nT;
      aver_nT += 1./(++nT)*(vect[i]-aver_nT);
      SnT += (vect[i]-aver_n1T)*(vect[i]-aver_nT);
    }
    
#pragma omp critical
    {
      NT += nT;
      if (NT>0) {
	aver_n1 = aver_n;
	aver_n += nT/NT*(aver_nT-aver_n);
	Sn += SnT+pow((aver_nT-aver_n1), 2)*nT*(NT-nT)/NT;
	sigma = sqrt(Sn/NT);
      }
    }
    
  }
  
  return sigma;
}


// ============================================================================


double cbl::Sigma (const std::vector<double> vect, const std::vector<double> weight) 
{
  if (vect.size()==0 || vect.size()!=weight.size()) ErrorCBL("Error in cbl::Sigma() in Func.cpp: 0 size for the input vector or vect.size()!=weight.size()");
  
  double aver_n1 = 0., aver_n = 0., Sn = 0., sigma = 0., WeightTOT = 0.;

#pragma omp parallel num_threads(omp_get_max_threads())
  {
    
    double aver_n1T = 0., aver_nT = 0., SnT = 0., WeightT = 0.;

#pragma omp for schedule(static, 2)
    for (size_t i=0; i<vect.size(); ++i) {
      aver_n1T = aver_nT;
      aver_nT += weight[i]/(WeightT+=weight[i])*(vect[i]-aver_nT);
      SnT += weight[i]*(vect[i]-aver_n1T)*(vect[i]-aver_nT);
    }
    
#pragma omp critical
    {
      WeightTOT += WeightT;
      if (WeightTOT>0) {
	aver_n1 = aver_n;
	aver_n += WeightT/WeightTOT*(aver_nT-aver_n);
	Sn += SnT+pow((aver_nT-aver_n1), 2)*WeightT*(WeightTOT-WeightT)/WeightTOT;
	sigma = sqrt(Sn/WeightTOT);
      }
    }
    
  }
  
  return sigma;
}


// ============================================================================


vector<double> cbl::Quartile (const std::vector<double> Vect) 
{
  vector<double> vect = Vect;
  sort(vect.begin(), vect.end()); 
  vector<double> vect1, vect2;
  
  int start, n = vect.size();
  double first = 0., second = 0., third = 0.;
  
  if (n>0) {
    if (n==1) {
      first = -1e10;
      second = vect[0];
      third = 1e10;
    }
    if (n>1) {
      if (n % 2 == 0)  // the number of elemens is even
	start = int(vect.size()*0.5);
      else 
	start = int(vect.size()*0.5)+1;
	  
      for (size_t i=0; i<vect.size()*0.5; i++)
	vect1.push_back(vect[i]);
      for (size_t i=start; i<vect.size(); i++)
	vect2.push_back(vect[i]);

      // first quartile
      n = vect1.size();
      if (n % 2 == 0) 
	first = (vect1[n*0.5-1]+vect1[(n*0.5)])*0.5;
      else 
	first = vect1[(n+1)*0.5-1];
	  
      // second quartile = median
      n = vect.size();
      if (n % 2 == 0)  
	second = (vect[n*0.5-1]+vect[(n*0.5)])*0.5;
      else
	second = vect[(n+1)*0.5-1];
	 
      // third quartile
      n = vect2.size();
      if (n % 2 == 0) 
	third = (vect2[n*0.5-1]+vect2[(n*0.5)])*0.5;
      else 
	third = vect2[(n+1)*0.5-1];
    }
  }

  return {first, second, third};
}


// ============================================================================


void cbl::Moment (const std::vector<double> data, double &ave, double &adev, double &sdev, double &var, double &skew, double &curt) 
{
  ave = gsl_stats_mean(data.data(), 1, data.size());  
  adev = gsl_stats_absdev_m(data.data(), 1, data.size(), ave);
  var = gsl_stats_variance_m(data.data(), 1, data.size(), ave);
  sdev = sqrt(var);
  skew = gsl_stats_skew_m_sd(data.data(), 1, data.size(), ave, sdev);
  curt = gsl_stats_kurtosis_m_sd(data.data(), 1, data.size(), ave, sdev);
}


// ============================================================================


double cbl::relative_error_beta (const double bias, const double Volume, const double density) // from Eq. 20 of Bianchi et al. 2012
{ 
  double n0 = 1.7e-4; // in (h/Mpc)^3
  double CC = 4.9e2;  // in (Mpc/h)^1.5

  return CC*pow(bias,0.7)/sqrt(Volume)*exp(n0/(bias*bias*density))*100.;
}


// ============================================================================


void cbl::measure_var_function (const std::vector<double> var, const int bin, const double V_min, const double V_max, const double Volume, std::vector<double> &Var, std::vector<double> &Phi, std::vector<double> &err)
{
  if (var.size()==0) ErrorCBL("Error in cbl::measure_var_functions() in Func.cpp: there are no objectes in the catalogue!");

  Var.erase(Var.begin(),Var.end());
  Phi.erase(Phi.begin(),Phi.end());
  err.erase(err.begin(),err.end());

  double delta_logV = (log10(V_max)-log10(V_min))/bin;
  double V1 = V_min;
  double V2 = V_min*pow(10.,delta_logV);

  for (int y=0; y<bin; y++) {
    double nHalo = 0.;

    for (size_t k=0; k<var.size(); k++) 
      if (V1<var[k] && var[k]<=V2) nHalo ++; 

    double PHI = nHalo/(V2-V1)/Volume;
    double ERR = sqrt(nHalo)/(V2-V1)/Volume;

    Var.push_back(pow(10.,(log10(V1)+log10(V2))*0.5));
    Phi.push_back(PHI);
    err.push_back(ERR);
  
    V1 = V2;
    V2 = V1*pow(10.,delta_logV);   
  }
}


// ============================================================================================


void cbl::bin_function (const std::string file_grid, double func(double, void *), void *par, const int bin, const double x_min, const double x_max, const std::string binning, std::vector<double> &xx, std::vector<double> &yy) 
{
  if (binning != "lin" && binning != "loglin" && binning != "log") 
    ErrorCBL("Error in cbl::bin_function() in Func.cpp: binning can only be: lin, loglin or log!");

  xx.resize(bin), yy.resize(bin);

  ifstream fin (file_grid.c_str());

  if (fin) {

    double XX, YY;    
    
    for (int i=0; i<bin; i++) {
      fin >>XX>>YY;
      xx[i] = XX;
      yy[i] = YY;
    }
    fin.clear(); fin.close();

    if (int(xx.size())!=bin) {
      string Err = "Error in cbl::bin_function() in Func.cpp: xx.size()=" + conv(xx.size(),par::fINT) + " != bin=" + conv(bin,par::fINT) + "!";
      ErrorCBL(Err);
    }

  }

  else {
    
    coutCBL <<"I'm creating the grid file: "<<file_grid<<"..."<<endl;

    fin.clear(); fin.close();

    double X_min = x_min;
    double X_max = x_max;
    
    if (binning != "lin") {
      if (x_min<0 || x_max<0) {
	string Err = "Error in cbl::bin_function() in Func.cpp: x_min=" + conv(x_min,par::fDP3) + ", x_max=" + conv(x_max,par::fDP3) + "!";
	ErrorCBL(Err);
      }
    
      X_min = log10(x_min);
      X_max = log10(x_max);
    }

    xx = linear_bin_vector(bin, X_min, X_max);

    ofstream fout(file_grid.c_str()); checkIO(fout, file_grid);

    for (size_t i=0; i<xx.size(); i++) {
      yy[i] = (binning=="lin") ? func(xx[i], par) : func(pow(10., xx[i]), par);

      if (binning=="log") {
	if (yy[i]<0) ErrorCBL("Error in cbl::bin_function() in Func.cpp: yy[i]<0!");
	else yy[i] = log10(yy[i]);
      }
      
      fout <<xx[i]<<"   "<<yy[i]<<endl; 
      coutCBL <<xx[i]<<"   "<<yy[i]<<endl; 
    }
    fout.clear(); fout.close(); coutCBL <<"I wrote the file: "<<file_grid<<endl;
  }

}


// ============================================================================================


void cbl::bin_function_2D (const std::string file_grid, double func(double *, size_t, void *), void *par, const int bin, const double x1_min, const double x1_max, const double x2_min, const double x2_max, const std::string binning, std::vector<double> &xx1, std::vector<double> &xx2, std::vector<std::vector<double>> &yy) 
{
  if (binning != "lin" && binning != "loglin" && binning != "log") 
    ErrorCBL("Error in cbl::bin_function_2D() in Func.cpp: binning can only be: lin, loglin or log !");

  xx1.resize(bin), xx2.resize(bin); yy.resize(bin);

  ifstream fin (file_grid.c_str());

  if (fin) {

    double XX1, XX2, YY;    
    
    for (int i=0; i<bin; i++) {
      for (int j=0; j<bin; j++) {
	fin >>XX1>>XX2>>YY;

	if (j==0) xx1[i] = XX1;
	if (i==0) xx2[j] = XX2;
	yy[i].push_back(YY);
      }
    }
    fin.clear(); fin.close();

    if (int(xx1.size())!=bin || int(xx2.size())!=bin) {
      string Err = "Error in cbl::bin_function_2D() in Func.cpp: xx1.size()=" + conv(xx1.size(),par::fINT) + ", xx2.size()=" + conv(xx2.size(),par::fINT) + " != bin=" + conv(bin,par::fINT) + "!";  
      ErrorCBL(Err);
    }

  }

  else {
    
    coutCBL <<"I'm creating the grid file: "<<file_grid<<"..."<<endl; 

    fin.clear(); fin.close();

    double X1_min = x1_min;
    double X1_max = x1_max;
    double X2_min = x2_min;
    double X2_max = x2_max;
    
    if (binning != "lin") {
      if (x1_min<0 || x1_max<0 || x2_min<0 || x2_max<0) {
	string Err = "Error in cbl::bin_function_2D() in Func.cpp: x1_min=" + conv(x1_min,par::fDP3) + ", x1_max=" + conv(x1_max,par::fDP3) + ", x2_min=" + conv(x2_min,par::fDP3) + ", x2_max=" + conv(x2_max,par::fDP3) + "!";
	ErrorCBL(Err);
      }
    
      X1_min = log10(x1_min);
      X1_max = log10(x1_max);
      X2_min = log10(x2_min);
      X2_max = log10(x2_max);
    }

    xx1 = linear_bin_vector(bin, X1_min, X1_max);
    xx2 = linear_bin_vector(bin, X2_min, X2_max);

    ofstream fout(file_grid.c_str()); checkIO(fout, file_grid);
    double vec[2];

    for (int i=0; i<bin; i++) {
      for (int j=0; j<bin; j++) {
	
	if (binning=="lin") {vec[0] = xx1[i]; vec[1] = xx2[j];} 
	else {vec[0] = pow(10.,xx1[i]); vec[1] = pow(10.,xx2[j]);}

	double ff = (binning=="log") ? log10(func(vec, 2, par)) : func(vec, 2, par);
	yy[i].push_back(ff);

	fout <<xx1[i]<<"   "<<xx2[j]<<"   "<<yy[i][j]<<endl; 
	coutCBL <<"--> "<<xx1[i]<<"   "<<xx2[j]<<"   "<<yy[i][j]<<endl;
      
      }
    }
    
    fout.clear(); fout.close(); coutCBL <<"I wrote the file: "<<file_grid<<endl;
  }

}


// ============================================================================================

/// @cond glob

double cbl::func_grid_lin (double xx, void *params)
{
  struct cbl::glob::STR_grid *pp = (struct cbl::glob::STR_grid *) params;

  return interpolated(xx, pp->_xx, pp->_yy, "Linear");
}


// ============================================================================================


double cbl::func_grid_loglin (double xx, void *params)
{
  struct cbl::glob::STR_grid *pp = (struct cbl::glob::STR_grid *) params;
 
  double lgx = log10(xx);

  return interpolated(lgx, pp->_xx, pp->_yy, "Linear");
}


// ============================================================================================


double cbl::func_grid_log (double xx, void *params)
{
  struct cbl::glob::STR_grid *pp = (struct cbl::glob::STR_grid *) params;
 
  double lgx = log10(xx);

  return pow(10., interpolated(lgx, pp->_xx, pp->_yy, "Linear"));
}


// ============================================================================================


double cbl::func_grid_lin_2D (double *xx, size_t dim, void *params)
{
  (void)dim;
  
  struct cbl::glob::STR_grid_2D *pp = (struct cbl::glob::STR_grid_2D *) params;
   
  return interpolated_2D(xx[0], xx[1], pp->_xx1, pp->_xx2, pp->_yy, "Linear");
}


// ============================================================================================


double cbl::func_grid_loglin_2D (double *xx, size_t dim, void *params)
{
  (void)dim;
  
  struct cbl::glob::STR_grid_2D *pp = (struct cbl::glob::STR_grid_2D *) params;
 
  double lgx1 = log10(xx[0]);
  double lgx2 = log10(xx[1]);

  return interpolated_2D(lgx1, lgx2, pp->_xx1, pp->_xx2, pp->_yy, "Linear");
}


// ============================================================================================


double cbl::func_grid_log_2D (double *xx, size_t dim, void *params)
{
  (void)dim;
  
  struct cbl::glob::STR_grid_2D *pp = (struct cbl::glob::STR_grid_2D *) params;

  double lgx1 = log10(xx[0]);
  double lgx2 = log10(xx[1]);

  return pow(10., interpolated_2D(lgx1, lgx2, pp->_xx1, pp->_xx2, pp->_yy, "Linear"));
}

/// @endcond

// ============================================================================================


void cbl::sdss_atbound(double &angle, const double minval, const double maxval)
{
  while (angle < minval)
    angle += 360.;

  while (angle > maxval)
    angle -= 360.0;
}

void cbl::sdss_atbound2(double &theta, double &phi)
{

  sdss_atbound(theta, -180.0, 180.0);

  if (fabs(theta)>90){
    theta = 180.0 - theta;
    phi += 180.0;
  }

  sdss_atbound(theta, -180.0, 180.0);
  sdss_atbound(phi, 0.0, 360.0);

  if (fabs(theta)==90.)
    phi = 0.0;
}


void cbl::eq2sdss (const std::vector<double> ra, const std::vector<double> dec, std::vector<double> &lambda, std::vector<double> &eta) 
{
  lambda.resize(ra.size());
  eta.resize(ra.size());

  double SurveyCenterRa = 185.-90, SurveyCenterDec = 32.5;
  double d2r = par::pi/180.;
  
  for (size_t i=0; i<ra.size(); i++) {
    double x = cos((ra[i]-SurveyCenterRa*d2r))*cos(dec[i]);
    double y = sin((ra[i]-SurveyCenterRa*d2r))*cos(dec[i]);
    double z = sin(dec[i]);

    lambda[i] = -asin(x)/d2r;
    eta[i] = atan2(z, y)/d2r - SurveyCenterDec;

    sdss_atbound(eta[i], -180.0, 180.0);
  }
}


// ============================================================================


void cbl::sdss2eq (const std::vector<double> lambda, const std::vector<double> eta, std::vector<double> &ra, std::vector<double> &dec) 
{
  ra.resize(lambda.size());
  dec.resize(lambda.size());

  double SurveyCenterRa = 185., SurveyCenterDec = 32.5;
  double d2r = par::pi/180.;
  
  for (size_t i=0; i<ra.size(); i++) {
    double x =  -1.0*sin(lambda[i]*d2r);
    double y = cos(lambda[i]*d2r)*cos(eta[i]*d2r+SurveyCenterDec*d2r);
    double z = cos(lambda[i]*d2r)*sin(eta[i]*d2r+SurveyCenterDec*d2r);
  
    ra[i] = atan2(y,x)/d2r + SurveyCenterRa-90;
    dec[i] = asin(z)/d2r;
    sdss_atbound2(dec[i], ra[i]);
  }
}


// ============================================================================


void cbl::sdss_stripe (const std::vector<double> eta, const std::vector<double> lambda, std::vector<int> &stripe, std::vector<int> &str_u)
{
  stripe.resize(eta.size());
  str_u.resize(eta.size());

  double stripe_sep=2.5;
  double cen = 58.75;

  for (size_t i=0; i<eta.size(); i++) {

    if (lambda[i]<90.0) stripe[i] = (eta[i]+cen)/stripe_sep;
    if (lambda[i]>90.0) stripe[i] = (eta[i]+cen+180.)/stripe_sep;

    str_u[i] = stripe[i];
  }

  sort(str_u.begin(), str_u.end());
  vector<int>::iterator it = unique(str_u.begin(), str_u.end());
  str_u.resize(distance(str_u.begin(), it));

}


// ============================================================================


double cbl::number_from_distribution (const std::vector<double> xx, const std::vector<double> fx, const double xmin, const double xmax, const int seed)
{
  return vector_from_distribution(1, xx, fx, xmin, xmax, seed)[0];
}


// ============================================================================


std::vector<double> cbl::vector_from_distribution (const int nRan, const std::vector<double> xx, const std::vector<double> fx, const double xmin, const double xmax, const int seed)
{
  random::UniformRandomNumbers ran(0., 1., seed);
  
  const int sz = 100;
  vector<double> Fx(sz, 0.);
  const vector<double> new_x = linear_bin_vector(sz, xmin, xmax); 

  glob::FuncGrid dist(xx, fx, "Spline");
  const double NN = dist.integrate_qag(xmin, xmax);

  for (int i=1; i<sz; ++i) 
    Fx[i] = dist.integrate_qag(xmin, new_x[i])/NN;
  
  glob::FuncGrid ff(Fx, new_x, "Spline");
  
  vector<double> varRandom(nRan);

  for (auto &&vR : varRandom)
    vR = ff(ran());
 
  return varRandom;
}


// ============================================================================


std::vector<size_t> cbl::minimum_maximum_indexes (const std::vector<double> xx, const double x_min, const double x_max)
{
  size_t ind1 = xx.size(), ind2 = 0;
  
  for (size_t i=0; i<xx.size(); i++) 
    if (x_min<xx[i] && xx[i]<x_max) {
      ind1 = min(i, ind1);
      ind2 = max(i, ind2);
    }
  
  ind2 ++;

  return {ind1, ind2};
}


// ============================================================================


void cbl::read_invert_covariance (const std::string filecov, std::vector< std::vector<double>> &cov, std::vector<std::vector<double>> &cov_inv, const size_t i1, const size_t i2)
{
  size_t size = i2-i1+1;

  cov.erase(cov.begin(), cov.end());
  cov_inv.erase(cov_inv.begin(), cov_inv.end());

  ifstream fin(filecov.c_str()); checkIO(fin, filecov);

  cov.erase(cov.begin(), cov.end());

  vector<double> vv;
  cov.push_back(vv);
  string line; size_t i = 0;

  while (getline(fin, line)) {
    stringstream ss(line);
    vector<double> num; double NN = -1.e30;
    while (ss>>NN) num.push_back(NN);
    if (num.size()==3 && num[2]>-1.e29) 
      cov[i].push_back(num[2]);
    else {i++; cov.push_back(vv);}
  }

  cov.erase(cov.end()-1, cov.end());
  fin.clear(); fin.close();

  cov_inv = cov;
  vector<vector<double>> cov_lim(size,vector<double>(size,0)), cov_lim_inv;

  size_t tot_size = cov.size();

  for (size_t i=0; i<tot_size; i++) {
    for (size_t j=0; j<tot_size; j++) {
      if (i>i2 || i<i1 || j>i2 || j<i1) {
      }
      else
	cov_lim[i-i1][j-i1] = cov[i][j]; 
    }
  }
  
  invert_matrix(cov_lim, cov_lim_inv);

  for (size_t i=0; i<tot_size; i++) {
    for (size_t j=0; j<tot_size; j++) {
      if (i>i2 || i<i1 || j>i2 || j<i1) 
	cov_inv[i][j] = 0.;
      else
	cov_inv[i][j] = cov_lim_inv[i-i1][j-i1]; 
    }
  }

}


// ============================================================================


void cbl::convolution (const std::vector<double> f1, const std::vector<double> f2, std::vector<double> &res, const double deltaX)
{
  size_t nn = f1.size();
  if (nn!=f2.size()) ErrorCBL("Error in cbl::convolution() in Func.cpp: the two functions have to have equal sizes!");

  double *ff1 = new double[nn];
  double *ff2 = new double[nn];
  double *rr = new double[nn];

  for (size_t i=0; i<nn; i++) {
    ff1[i] = f1[i];
    ff2[i] = f2[i];
    rr[i] = 0.;
  }
  
  gsl_fft_real_wavetable *real; 
  gsl_fft_halfcomplex_wavetable *hc;
  gsl_fft_real_workspace *work;

  work = gsl_fft_real_workspace_alloc(nn);
  real = gsl_fft_real_wavetable_alloc(nn);

  gsl_fft_real_transform(ff1, 1., nn, real, work);
  gsl_fft_real_transform(ff2, 1., nn, real, work);

  gsl_fft_real_wavetable_free(real);

  
  // Combining Fourier series coefficients (complex number product)
 
  rr[0] = ff1[0]*ff2[0];
  
  if (nn%2!=0) {
    for (unsigned int i=1; i<nn; i+=2) {
      rr[i] = ff1[i]*ff2[i]-ff1[i+1]*ff2[i+1];
      rr[i+1] = ff1[i]*ff2[i+1]+ff2[i]*ff1[i+1];
    }
  }
  else {
    for (unsigned int i=1; i<nn-1; i+=2) {
      rr[i] = ff1[i]*ff2[i]-ff1[i+1]*ff2[i+1];
      rr[i+1] = ff1[i]*ff2[i+1]+ff2[i]*ff1[i+1];
    }
    rr[nn-1] = ff1[nn-1]*ff2[nn-1];
  }


  // Back to real space

  hc = gsl_fft_halfcomplex_wavetable_alloc(nn);
  
  gsl_fft_halfcomplex_inverse(rr, 1., nn, hc, work);

  gsl_fft_halfcomplex_wavetable_free(hc);

   
  // Re-order result
  
  res.resize(nn, 0.);

  for (unsigned int i=0; i<=nn/2; i++)
    res[i] = deltaX*rr[i+nn/2];
  
  for (unsigned int i=nn/2+1; i<nn; i++)
    res[i] = deltaX*rr[i-nn/2-1];

  if (nn%2==0) res[nn/2] = res[nn/2-1];

}


// ============================================================================


void cbl::distribution (std::vector<double> &xx, std::vector<double> &fx, std::vector<double> &err, const std::vector<double> FF, const std::vector<double> WW, const int nbin, const bool linear, const std::string file_out, const double fact, const double V1, const double V2, const string bin_type, const bool conv, const double sigma)
{
  if (xx.size()>0 || fx.size()>0 || FF.size()<=0 || nbin<=0) ErrorCBL("Error in cbl::distribution() in Func.cpp:  the following conditions have to be satisfied: xx.size()<=0, fx.size()<=0, FF.size()>0 and nbin>0. The values recived are instead: xx.size() = " + cbl::conv(xx.size(), par::fINT) + ", fx.size() = " + cbl::conv(fx.size(), par::fINT) + ", FF.size() = " + cbl::conv(FF.size(), par::fINT) + "and nbin = " + cbl::conv(nbin, par::fINT) + "!");
  
  double minFF = (V1>cbl::par::defaultDouble) ? V1 : Min(FF)*0.9999;
  double maxFF = (V2>cbl::par::defaultDouble) ? V2 : Max(FF)*1.0001;

  
  // using GSL to create the histogram 

  gsl_histogram *histo = gsl_histogram_alloc(nbin);

  if (linear) gsl_histogram_set_ranges_uniform(histo, minFF, maxFF);

  else {
    vector<double> vv = logarithmic_bin_vector(nbin+1, minFF, maxFF);
    double *vvv = new double[nbin+1]; for (int i=0; i<nbin+1; i++) vvv[i] = vv[i];
    gsl_histogram_set_ranges(histo, vvv, nbin+1);
  }

  vector<double> Weight = WW;
  if (Weight.size()==0) Weight.resize(FF.size(), 1.);
  checkDim(Weight, FF.size(), "WW");
  
  for (size_t i=0; i<FF.size(); i++)
    gsl_histogram_accumulate(histo, FF[i], Weight[i]);
  
  double x1, x2;

  for (int i=0; i<nbin; i++) {

    gsl_histogram_get_range(histo, i, &x1, &x2);
    double val = gsl_histogram_get(histo, i);
    
    if (linear) xx.push_back(0.5*(x1+x2));
    else xx.push_back(pow(10., 0.5*(log10(x1)+log10(x2))));

    if (bin_type == "Linear") {
      fx.push_back(val/((x2-x1)*fact));
      err.push_back(sqrt(val)/((x2-x1)*fact));
    }

    else if (bin_type == "Log") {
      fx.push_back(val/((log(x2)-log(x1))*fact));
      err.push_back(sqrt(val)/((log(x2)-log(x1))*fact));
    }

    else if (bin_type == "Log10") {
      fx.push_back(val/((log10(x2)-log10(x1))*fact));
      err.push_back(sqrt(val)/((log10(x2)-log10(x1))*fact));
    }
    
  else 
    ErrorCBL("Error in cbl::distribution() in Func.cpp: the value of string 'bin_type' is not permitted, possible selections are 'Linear', 'Log', 'Log10'!");

    

  }

  
  // Gaussian convolution

  if (conv) {
    coutCBL << "The distribution is smoothed with a Gaussian filter" << endl;
    double *func;
    fftw_complex *func_tr;

    if (!linear) ErrorCBL("Work in progress...", ExitCode::_workInProgress_);
    int nbinN = 2*nbin;
    int i1 = nbin*0.5, i2 = 1.5*nbin;

    int nbinK = 0.5*nbinN+1;

    func = fftw_alloc_real(nbinN);
    func_tr = fftw_alloc_complex(nbinK);

    for (int i=0; i<nbinN; i++)
      func[i] = 0;
    
    for (int i=i1; i<i2; i++)
      func[i] = fx[i-i1];
    
    for (int i=0; i<nbinK; i++) {
      func_tr[i][0] = 0;
      func_tr[i][1] = 0;
    }

    fftw_plan real2complex;
    real2complex = fftw_plan_dft_r2c_1d(nbinN, func, func_tr, FFTW_ESTIMATE);
    fftw_execute(real2complex);
    fftw_destroy_plan(real2complex);

    double delta = (maxFF-minFF)/nbin;
    double SS = pow(sigma,2);

    double fact = 2*par::pi/(nbinN*delta);
    for (int i=0; i<nbinK; i++) {
      double kk = i*fact;
      func_tr[i][0] = func_tr[i][0]*exp(-0.5*kk*kk*SS);
      func_tr[i][1] = func_tr[i][1]*exp(-0.5*kk*kk*SS);
    }
    
    fftw_plan complex2real;
    complex2real = fftw_plan_dft_c2r_1d(nbinN, func_tr, func, FFTW_ESTIMATE);
    fftw_execute(complex2real);
    fftw_destroy_plan(complex2real);

    for (int i=i1; i<i2; i++)
      fx[i-i1] = func[i]/nbinN;
  }

  
  if (file_out!=par::defaultString && file_out!="") {

    ofstream fout(file_out.c_str()); checkIO(fout, file_out);
    
    for (size_t i=0; i<xx.size(); i++)
      fout << xx[i] << "   " << fx[i] << "   " << err[i] << endl;

    fout.clear(); fout.close(); coutCBL << "I wrote the file: " << file_out << endl;
  }

  gsl_histogram_free(histo);
  fftw_cleanup();

}


// ============================================================================


double cbl::legendre_polynomial (const double mu, const int l)
{
  return gsl_sf_legendre_Pl(l,mu);
}


// ============================================================================


double cbl::legendre_polynomial_integral (double mu, void *params)
{
  //  cbl::glob::STR_jl_distance_average *par = (cbl::glob::STR_jl_distance_average *)(params);
  int *l = (int *)(params);
  return legendre_polynomial(mu, *l);
}


// ============================================================================


double cbl::Legendre_polynomial_mu_average (const int ll, const double mu, const double delta_mu)
{
  gsl_function Func;
  Func.function = &cbl::legendre_polynomial_integral;
  int lll = ll;
  Func.params = &lll;

  return gsl::GSL_integrate_qag(Func, mu, mu+delta_mu, 1.e-3, 1000, 6)/delta_mu;
}


// ============================================================================


double cbl::Legendre_polynomial_mu_average (const double mu_min, const double mu_max, const int ll)
{
  auto integrand = [&] (const double mu) { return cbl::legendre_polynomial(mu, ll); };

  return gsl::GSL_integrate_qag(integrand, mu_min, mu_max, 1.e-3, 1.e-6, 1000, 6);
}


// ============================================================================


double cbl::Legendre_polynomial_theta_average (const double theta_min, const double theta_max, const int ll)
{
  auto integrand = [&] (const double theta) { return sin(theta)*cbl::legendre_polynomial(cos(theta), ll); };

  return gsl::GSL_integrate_qag(integrand, theta_min, theta_max, 1.e-3, 1.e-6, 1000, 6);
}


// ============================================================================


double cbl::Legendre_polynomial_triangles_average (const double r12_min, const double r12_max, const double r13_min, const double r13_max, const double r23_min, const double r23_max, const int ll, const double rel_err, const double abs_err, const int nevals)
{
  auto r12_integrand = [&] (const double r12) {

    auto r13_integrand = [&] (const double r13) {

      double mu_min = max(min( (r12*r12+r13*r13-r23_max*r23_max)/(2*r12*r13), 1.), -1.);
      double mu_max = min(max( (r12*r12+r13*r13-r23_min*r23_min)/(2*r12*r13), -1.), 1.);
       
      auto mu_integrand = [&] (const double mu) {
	return cbl::legendre_polynomial(mu, ll);
      };

      return ((mu_min == mu_max) ? 0 : cbl::gsl::GSL_integrate_cquad(mu_integrand, mu_min, mu_max, rel_err, abs_err, nevals));
    };

    return cbl::gsl::GSL_integrate_cquad(r13_integrand, r13_min, r13_max, rel_err, abs_err, nevals);
  };

  return cbl::gsl::GSL_integrate_cquad(r12_integrand, r12_min, r12_max, rel_err, abs_err, nevals)/((r12_max-r12_min)*(r13_max-r13_min));
}


// ============================================================================


complex<double> cbl::spherical_harmonics (const int l, const int m, const double xx, const double yy, const double zz)
{
  const double sintheta = sin(acos(zz));
  complex<double> exp_iphi(xx/(sintheta), yy/(sintheta));
  complex<double> pow_exp = pow(exp_iphi, m);

  double fact = pow(-1., m)*gsl_sf_legendre_sphPlm (l, m, zz);

  return fact*pow_exp;
}


// ============================================================================


std::vector<std::vector<complex<double>>> cbl::spherical_harmonics (const int lmax, const double xx, const double yy, const double zz)
{
  const int nl = lmax+1;

  const double sintheta = sin(acos(zz));
  complex<double> exp_iphi(xx/(sintheta), yy/(sintheta));

  vector<complex<double>> pow_exp(nl+1);
  vector<double> norm(nl+1);

  for (int mm=0; mm<nl+1; mm++){
    norm[mm] = pow(-1., mm);
    pow_exp[mm] = pow(exp_iphi, mm);
  }

  vector<vector<complex<double>>> coeff(nl);

  for(int ll=0; ll< nl; ll++){
    
    vector<complex<double>> ylm(ll+1);
    
    for (int mm=0; mm<ll+1; mm++){
      double fact = norm[mm]*gsl_sf_legendre_sphPlm (ll, mm, zz);
      ylm[mm] = fact*pow_exp[mm];
    }

    coeff[ll] = ylm;
  }

  return coeff;
}


// ============================================================================


std::vector<std::complex<double>> cbl::spherical_harmonics_array (const int lmax, const double xx, const double yy, const double zz)
{
  const int n_sph = gsl_sf_legendre_array_n(lmax);

  vector<double> Plm(n_sph, 0);
  vector<complex<double>> sph(n_sph);

  double phi = atan2(yy, xx);

  vector<complex<double>> pow_exp(lmax+2, complex<double>(cos(phi), sin(phi)));

  for (int mm=0; mm<lmax+2; mm++)
    pow_exp[mm] = pow(pow_exp[mm], mm);

  gsl_sf_legendre_array_e(GSL_SF_LEGENDRE_SPHARM, lmax, zz, 1., Plm.data());

  
  int n=0;
  for(int l=0; l<lmax+1; l++)
    for (int m=0; m<l+1; m++){
      sph[n] = pow_exp[m]*Plm[n]; 
      n++;
    }

  return sph;
}

// ============================================================================


double cbl::j0 (const double xx)
{
  return gsl_sf_bessel_jl(0, xx);
}


// ============================================================================


double cbl::j2 (const double xx)
{
  return gsl_sf_bessel_jl(2, xx);
}


// ============================================================================


double cbl::j4 (const double xx)
{
  return gsl_sf_bessel_jl(4, xx);
}


// ============================================================================


double cbl::jl (const double xx, const int order)
{
  return gsl_sf_bessel_jl(order, xx);
}


// ============================================================================


double cbl::j0_distance_average (const double kk, const double r_down, const double r_up)
{
  double volume = (pow(r_up,3)-pow(r_down,3))/3;
  double up = (sin(kk*r_up)-kk*r_up*cos(kk*r_up))*pow(kk,-3);
  double down = (sin(kk*r_down)-kk*r_down*cos(kk*r_down))*pow(kk,-3);
  return (up-down)/volume;
}


// ============================================================================


double cbl::j2_distance_average (const double kk, const double r_down, const double r_up)
{
  double volume = (pow(r_up,3)-pow(r_down,3))/3;
  double up = (3*gsl_sf_Si(kk*r_up)-4*sin(kk*r_up)+kk*r_up*cos(kk*r_up))*pow(kk,-3);
  double down = (3*gsl_sf_Si(kk*r_down)-4*sin(kk*r_down)+kk*r_down*cos(kk*r_down))*pow(kk,-3);
  return (up-down)/volume;
}


// ============================================================================


double cbl::jl_spherical_integrand (double rr, void *params)
{
  cbl::glob::STR_jl_distance_average *par = (cbl::glob::STR_jl_distance_average *)(params);
  return rr*rr*gsl_sf_bessel_jl(par->order, par->k*rr);
}


// ============================================================================


double cbl::jl_distance_average (const double kk, const int order, const double r_down, const double r_up)
{
  double volume = (pow(r_up,3)-pow(r_down,3))/3;

  cbl::glob::STR_jl_distance_average str;
  str.order = order;
  str.k = kk;

  gsl_function Func;
   
  Func.function=&cbl::jl_spherical_integrand;
  Func.params=&str;

  double prec=1.e-2;
  int limit_size = 1000;

  double Int = cbl::gsl::GSL_integrate_qag(Func, r_down, r_up, prec, limit_size, 6);

  return Int/volume;
}


// ============================================================================


std::vector<double> cbl::generate_correlated_data (const std::vector<double> mean, const std::vector<std::vector<double>> covariance, const int seed)
{
  random::UniformRandomNumbers ran(0., 1., seed);
  
  size_t sample_size = mean.size();
  vector<double> sample;
  vector<double> std;

  gsl_matrix *correlation = gsl_matrix_alloc(sample_size, sample_size);
  
  for (size_t i=0; i<sample_size; i++) {
    std.push_back(sqrt(covariance[i][i]));
    sample.push_back(ran());
    for (size_t j=0; j<sample_size; j++) {
      double corr = covariance[i][j]/sqrt(covariance[i][i]*covariance[j][j]);
      if (corr!=corr)
        ErrorCBL("Error cbl::generate_correlated_data() in Func.cpp: negative value on the covariance diagonal!");
      gsl_matrix_set(correlation,i,j, corr);
    }
  }

  
  // correlation matrix eigensystem
  
  gsl_vector *eigenvalues = gsl_vector_alloc(sample_size);
  gsl_matrix *VV = gsl_matrix_alloc(sample_size, sample_size);
  gsl_matrix_set_zero(VV);

  gsl_matrix *eigenvectors = gsl_matrix_alloc(sample_size, sample_size);

  gsl_eigen_symmv_workspace *ww = gsl_eigen_symmv_alloc(sample_size); 
  gsl_eigen_symmv (correlation, eigenvalues, eigenvectors, ww);
  gsl_eigen_symmv_free(ww);

  for (size_t j = 0; j<sample_size; j++) {
    for (size_t i = 0; i<sample_size; i++) {
      if (gsl_vector_get(eigenvalues, j)<0)
        ErrorCBL("Error cbl::generate_correlated_data() in Func.cpp: covariance matrix must be positive (semi-)definite but has at least one negative eigenvalue!");
      double v1 = gsl_matrix_get(eigenvectors, i, j);
      double v2 = sqrt(gsl_vector_get(eigenvalues, j));
      gsl_matrix_set(VV, i, j, v1*v2);
    }
  }
  
  vector<double> cov_sample;
  for (size_t i=0; i<sample_size; i++) {
    gsl_vector *row = gsl_vector_alloc(sample_size);
    gsl_matrix_get_row(row, VV, i);
    cov_sample.push_back(0);
    for (size_t j=0; j<sample_size; j++)
      cov_sample[i]+=gsl_vector_get(row, j)*sample[j];
    cov_sample[i] = std[i]*cov_sample[i]+mean[i];
  }

  return cov_sample;
}


// ============================================================================


double cbl::trapezoid_integration (const std::vector<double> xx, const std::vector<double> yy)
{
  double Int = 0.;

  for (size_t i=0; i<xx.size()-1; i++)
    Int += 0.5*(yy[i+1]+yy[i])*(xx[i+1]-xx[i]);
  
  return Int;
}


// ============================================================================


std::vector<std::vector<double>> cbl::generate_correlated_data (const int nExtractions, const std::vector<double> mean, const std::vector<std::vector<double>> covariance, const int seed)
{
  random::UniformRandomNumbers ran(0., 1., seed);
  
  size_t sample_size = mean.size();
  vector<double> std;

  gsl_matrix *correlation = gsl_matrix_alloc(sample_size,sample_size);
  
  for (size_t i=0; i<sample_size; i++) {
    std.push_back(sqrt(covariance[i][i]));
    for (size_t j=0; j<sample_size; j++) {
      double corr = covariance[i][j]/sqrt(covariance[i][i]*covariance[j][j]);
      if (corr!=corr)
        ErrorCBL("Error cbl::generate_correlated_data() in Func.cpp: negative value on the covariance diagonal!");
      gsl_matrix_set(correlation,i,j, corr);
    }
  }

  vector<vector<double>> sample;
  for (int j=0; j<nExtractions; j++) {
    vector<double> subS(sample_size, 0);
    for (size_t i=0; i<sample_size; i++)
      subS[i] = ran();
    sample.push_back(subS);
  }

  
  // correlation matrix eigensystem
  
  gsl_vector *eigenvalues = gsl_vector_alloc(sample_size);
  gsl_matrix *VV = gsl_matrix_alloc(sample_size,sample_size);
  gsl_matrix_set_zero(VV);

  gsl_matrix *eigenvectors = gsl_matrix_alloc(sample_size,sample_size);

  gsl_eigen_symmv_workspace *ww = gsl_eigen_symmv_alloc (sample_size); 
  gsl_eigen_symmv (correlation, eigenvalues, eigenvectors, ww);
  gsl_eigen_symmv_free(ww);

  for (size_t j=0; j<sample_size; j++) {
    for (size_t i=0; i<sample_size; i++) {
      if (gsl_vector_get(eigenvalues, j)<0)
        ErrorCBL("Error cbl::generate_correlated_data() in Func.cpp: covariance matrix must be positive (semi-)definite but has at least one negative eigenvalue!");
        
      double v1 = gsl_matrix_get(eigenvectors, i, j);
      double v2 = sqrt(gsl_vector_get(eigenvalues, j));
      gsl_matrix_set(VV, i, j, v1*v2);
    }
  }
  
  vector<vector<double>> cov_sample;

  for (int k=0; k<nExtractions; k++) {
    vector<double> cov_SSample(sample_size, 0);
    for (size_t i=0; i<sample_size; i++) {

      gsl_vector *row = gsl_vector_alloc(sample_size);
      gsl_matrix_get_row(row,VV,i);

      for (size_t j=0; j<sample_size; j++)
	cov_SSample[i] += gsl_vector_get(row, j)*sample[k][j];
      
      cov_SSample[i] = std[i]*cov_SSample[i]+mean[i];
    }
    cov_sample.push_back(cov_SSample);
  }
  
  return cov_sample;
}



// ============================================================================


/* ======== Cosimo Fedeli ======== */

/// @cond glob

void cbl::gauleg (const double x1, const double x2, double *x, double *w, const int n)
{
  const double eps = 3.0e-11;
  int m = (n+1)/2;
  double xm = 0.5*(x2+x1);
  double xl = 0.5*(x2-x1);
  for (int i=1; i<=m; i++) {
    double z1, pp;
    double z = cos(par::pi*(i-0.25)/(n+0.5));
    do {
      double p1 = 1.0;
      double p2 = 0.0;
      for (int j=1; j<=n; j++) {
	double p3 = p2;
	p2 = p1;
	p1 = ((2.0*j-1.0)*z*p2-(j-1.0)*p3)/j;
      }
      pp = n*(z*p1-p2)/(z*z-1.0);
      z1 = z;
      z = z1-p1/pp;
    }
    while (fabs(z-z1)>eps);
    x[i-1] = xm-xl*z;
    x[n-i] = xm+xl*z;
    w[i-1] = 2.0*xl/((1.0-z*z)*pp*pp);
    w[n-i] = w[i-1];
  }
}

/// @endcond
