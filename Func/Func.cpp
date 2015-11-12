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
 *  @file /home/frea/CosmoBolognaLib/Func/Func.cpp
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
using namespace cosmobl;


// ============================================================================


void cosmobl::unique_unsorted (vector<int> &vv) // erase all equal elements
{
  sort(vv.begin(),vv.end());
  vector<int>::iterator it;
  it = unique(vv.begin(),vv.end());
  vv.resize(it-vv.begin()); 
}

void cosmobl::unique_unsorted (vector<double> &vv) // erase all equal elements
{
  sort(vv.begin(),vv.end());
  vector<double>::iterator it;
  it = unique(vv.begin(),vv.end());
  vv.resize(it-vv.begin()); 
}

void cosmobl::unique_unsorted (vector<string> &vv) // erase all equal elements
{
  sort(vv.begin(),vv.end());
  vector<string>::iterator it;
  it = unique(vv.begin(),vv.end());
  vv.resize(it-vv.begin()); 
}


// ============================================================================

/// @cond glob
bool cosmobl::glob::operator<(const cosmobl::glob::CL &c1, const cosmobl::glob::CL &c2) {return c1.VV[0] < c2.VV[0];}
/// @endcond

void cosmobl::sort_2vectors (vector<double>::iterator p1, vector<double>::iterator p2, int dim) 
{
  int temp = 0;
  vector<cosmobl::glob::CL> ccc; 
  for (int i=0; i<dim; i++) { 
    vector<double> vect = {*p1, *p2}; 
    cosmobl::glob::CL cl(vect); 
    ccc.push_back(cl);
    if (i+1<dim) {*p1 ++; *p2 ++; temp ++;}
  }
  sort (ccc.begin(),ccc.end());
  p1 -= temp; p2 -= temp;
  for (int i=0; i<dim; i++) { 
    *p1 = ccc[i].VV[0]; *p2 = ccc[i].VV[1];
    if (i+1<dim) {*p1 ++; *p2 ++;} 
  }
}

void cosmobl::sort_3vectors (vector<double>::iterator p1, vector<double>::iterator p2, vector<double>::iterator p3, int dim) 
{
  int temp = 0;
  vector<cosmobl::glob::CL> ccc;
  for (int i=0; i<dim; i++) { 
    vector<double> vect = {*p1, *p2, *p3};
    cosmobl::glob::CL cl(vect); 
    ccc.push_back(cl);
    if (i+1<dim) {*p1 ++; *p2 ++; *p3 ++; temp ++;}
  }
  sort (ccc.begin(),ccc.end());
  p1 -= temp; p2 -= temp; p3 -= temp;
  for (int i=0; i<dim; i++) { 
    *p1 = ccc[i].VV[0]; *p2 = ccc[i].VV[1]; *p3 = ccc[i].VV[2];
    if (i+1<dim) {*p1 ++; *p2 ++; *p3 ++;} 
  }
}

void cosmobl::sort_4vectors (vector<double>::iterator p1, vector<double>::iterator p2, vector<double>::iterator p3, vector<double>::iterator p4, int dim) 
{
  int temp = 0;
  vector<cosmobl::glob::CL> ccc;
  for (int i=0; i<dim; i++) { 
    vector<double> vect = {*p1, *p2, *p3, *p4};
    cosmobl::glob::CL cl(vect); 
    ccc.push_back(cl);
    if (i+1<dim) {*p1 ++; *p2 ++; *p3 ++; *p4 ++; temp ++;}
  }
  sort (ccc.begin(),ccc.end());
  p1 -= temp; p2 -= temp; p3 -= temp; p4 -= temp;
  for (int i=0; i<dim; i++) { 
    *p1 = ccc[i].VV[0]; *p2 = ccc[i].VV[1]; *p3 = ccc[i].VV[2]; *p4 = ccc[i].VV[3];
    if (i+1<dim) {*p1 ++; *p2 ++; *p3 ++; *p4 ++;} 
  }
}


// ============================================================================================


void cosmobl::polar_coord (double XX, double YY, double ZZ, double *ra, double *dec, double *dd) 
{   
  *dd = sqrt(XX*XX+YY*YY+ZZ*ZZ);
  *ra = atan2(XX,YY);
  *dec = asin(ZZ/(*dd));  
}

void cosmobl::cartesian_coord (double ra, double dec, double dd, double *XX, double *YY, double *ZZ) 
{
  *XX = dd*cos(dec)*sin(ra);
  *YY = dd*cos(dec)*cos(ra);
  *ZZ = dd*sin(dec);

  /*
  *XX = dd*cos(ra)*sin(0.5*par::pi-dec);
  *YY = dd*sin(ra)*sin(0.5*par::pi-dec);
  *ZZ = dd*cos(0.5*par::pi-dec);
  */
}

void cosmobl::polar_coord (vector<double> XX, vector<double> YY, vector<double> ZZ, vector<double> &ra, vector<double> &dec, vector<double> &dd) 
{     
  for (unsigned int i=0; i<XX.size(); i++) {
    dd[i] = sqrt(XX[i]*XX[i]+YY[i]*YY[i]+ZZ[i]*ZZ[i]);
    ra[i] = atan2(XX[i],YY[i]);
    dec[i] = asin(ZZ[i]/dd[i]);
  }  
}

void cosmobl::cartesian_coord (vector<double> ra, vector<double> dec, vector<double> dd, vector<double> &XX, vector<double> &YY, vector<double> &ZZ) 
{
  for (unsigned int i=0; i<XX.size(); i++) {
    XX[i] = dd[i]*cos(dec[i])*sin(ra[i]);
    YY[i] = dd[i]*cos(dec[i])*cos(ra[i]);
    ZZ[i] = dd[i]*sin(dec[i]);
    /*
    XX[i] = dd[i]*cos(ra[i])*sin(0.5*par::pi-dec[i]);
    YY[i] = dd[i]*sin(ra[i])*sin(0.5*par::pi-dec[i]);
    ZZ[i] = dd[i]*cos(0.5*par::pi-dec[i]);
    */
  }
}


// ============================================================================================


double cosmobl::MC_Int (double func(double), double x1, double x2) 
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

  Ran myran(5);
  double xt, yt, INT;
  int sub = 0, subn = 0, numTOT = 10000000;
  

  if (f1>0) {
    for (int i=0; i<numTOT; i++) {
      xt = myran.doub()*(x2-x1)+x1;
      yt = myran.doub()*(f2-f1);
      if (yt<func(xt)) sub ++;
    }
    INT = double(sub)/double(numTOT)*(x2-x1)*(f2-f1);
  }
  else {
    for (int i=0; i<numTOT; i++) {
      xt = myran.doub()*(x2-x1)+x1;
      yt = myran.doub()*f2;
      if (yt<func(xt)) sub ++;
    }
    for (int i=0; i<numTOT; i++) {
      xt = myran.doub()*(x2-x1)+x1;
      yt = myran.doub()*f1;
      if (yt>func(xt)) subn ++;
    }
    INT = (double(sub)/double(numTOT)*(x2-x1)*f2)-(double(subn)/double(numTOT)*(x2-x1)*fabs(f1));
  }
  
  return INT;
}

double cosmobl::MC_Int (double func(double, double AA), double AA, double x1, double x2) 
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

  Ran myran(5);
  double xt, yt, INT;
  int sub = 0, subn = 0, numTOT = 10;

  if (f1>0) {
    for (int i=0; i<numTOT; i++) {
      xt = myran.doub()*(x2-x1)+x1;
      yt = myran.doub()*(f2-f1);
      if (yt<func(xt,AA)) sub ++;
    }
    INT = double(sub)/double(numTOT)*(x2-x1)*(f2-f1);
  }
  else {
    for (int i=0; i<numTOT; i++) {
      xt = myran.doub()*(x2-x1)+x1;
      yt = myran.doub()*f2;
      if (yt<func(xt,AA)) sub ++;
    }
    for (int i=0; i<numTOT; i++) {
      xt = myran.doub()*(x2-x1)+x1;
      yt = myran.doub()*f1;
      if (yt>func(xt,AA)) subn ++;
    }
    INT = (double(sub)/double(numTOT)*(x2-x1)*f2)-(double(subn)/double(numTOT)*(x2-x1)*fabs(f1));
  }
  
  return INT;
}

double cosmobl::MC_Int (double func(double, double AA, double BB, double CC, double DD, double EE), double AA, double BB, double CC, double DD, double EE, double x1, double x2) 
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

  Ran myran(5);
  double xt, yt, INT;
  int sub = 0, subn = 0, numTOT = 100000;
  

  if (f1>0) {
    for (int i=0; i<numTOT; i++) {
      xt = myran.doub()*(x2-x1)+x1;
      yt = myran.doub()*(f2-f1);
      if (yt<func(xt,AA,BB,CC,DD,EE)) sub ++;
    }
    INT = double(sub)/double(numTOT)*(x2-x1)*(f2-f1);
  }
  else {
    for (int i=0; i<numTOT; i++) {
      xt = myran.doub()*(x2-x1)+x1;
      yt = myran.doub()*f2;
      if (yt<func(xt,AA,BB,CC,DD,EE)) sub ++;
    }
    for (int i=0; i<numTOT; i++) {
      xt = myran.doub()*(x2-x1)+x1;
      yt = myran.doub()*f1;
      if (yt>func(xt,AA,BB,CC,DD,EE)) subn ++;
    }
    INT = (double(sub)/double(numTOT)*(x2-x1)*f2)-(double(subn)/double(numTOT)*(x2-x1)*fabs(f1));
  }
  
  return INT;
}


// ============================================================================================


short cosmobl::ShortSwap (short s)
{
  unsigned char b1, b2;
  b1 = s & 255;
  b2 = (s>>8) & 255;
  return (b1<<8) + b2;
}

int cosmobl::IntSwap (int i)
{
  unsigned char b1, b2, b3;
  b1 = i & 255;
  b2 = (i>>8) & 255;
  b3 = (i>>16) & 255;
  return ((int)b1<<16) + ((int)b2<<8) + b3;
}

long cosmobl::LongSwap (long i)
{
  unsigned char b1, b2, b3, b4;
  b1 = i & 255;
  b2 = (i>>8) & 255;
  b3 = (i>>16) & 255;
  b4 = (i>>24) & 255;
  return ((int)b1<<24) + ((int)b2<<16) + ((int)b3<<8) + b4;
}

long long cosmobl::LongLongSwap (long long i)
{
  unsigned char b1, b2, b3, b4, b5, b6, b7, b8;
  b1 = i & 255;
  b2 = (i>>8) & 255;
  b3 = (i>>16) & 255;
  b4 = (i>>24) & 255;
  b5 = (i>>32) & 255;
  b6 = (i>>40) & 255;
  b7 = (i>>48) & 255;
  b8 = (i>>56) & 255;
  return ((long long)b1<<56) + ((long long)b2<<48) + ((long long)b3<<40) + ((long long)b4<<32) + ((long long)b5<<24) + ((long long)b6<<16) + ((long long)b7<<8) + b8;
}

float cosmobl::FloatSwap (float f)
{
  union
  {
    float f;
    unsigned char b[4];
  } dat1, dat2;
  dat1.f = f;
  dat2.b[0] = dat1.b[3];
  dat2.b[1] = dat1.b[2];
  dat2.b[2] = dat1.b[1];
  dat2.b[3] = dat1.b[0];
  return dat2.f;
}

double cosmobl::DoubleSwap (double d)
{
  union
  {
    double d;
    unsigned char b[8];
  } dat1, dat2;
  dat1.d = d;
  dat2.b[0] = dat1.b[7];
  dat2.b[1] = dat1.b[6];
  dat2.b[2] = dat1.b[5];
  dat2.b[3] = dat1.b[4];
  dat2.b[4] = dat1.b[3];
  dat2.b[5] = dat1.b[2];
  dat2.b[6] = dat1.b[1];
  dat2.b[7] = dat1.b[0];
  return dat2.d;
}


// ============================================================================


double cosmobl::interpolated (double _xx, vector<double> xx, vector<double> yy, string type, int nPt, double &err)
{
  if (xx.size()!=yy.size()) 
    ErrorMsg("Error in interpolated of Func.cpp)!");

  VecDoub XX(xx.size()), YY(yy.size());
  for (int i=0; i<XX.size(); i++) {
    XX[i] = xx[i];
    YY[i] = yy[i];
  }
  
  if (type=="Linear") {
    Linear_interp func (XX, YY);
    return func.interp(_xx);
  }

  else if (type=="Poly") {
    Poly_interp func (XX, YY, nPt);
    err = func.dy;
    return func.interp(_xx);
  }

  else if (type=="Spline") {
    Spline_interp func (XX, YY);
    return func.interp(_xx); 
  }

  else if (type=="Rat") {
    Rat_interp func (XX, YY, nPt);
    err = func.dy;
    return func.interp(_xx);
  }

  else if (type=="BaryRat") {
    BaryRat_interp func (XX, YY, nPt);
    return func.interp(_xx);
  }

  else 
    ErrorMsg("Error in interpolated of Func.cpp: the value of string 'type' is not permitted!");

  return -1;
}


double cosmobl::interpolated (double _xx, vector<double> xx, vector<double> yy, string type, int nPt)
{
  double err = -1.;
  return interpolated(_xx, xx, yy, type, nPt, err);
}


// ============================================================================


double cosmobl::interpolated_2D (double _x1, double _x2, vector<double> x1, vector<double> x2, vector<vector<double> > yy, string type, int nPt)
{
  VecDoub X1(x1.size()), X2(x2.size());
  MatDoub YY(x1.size(),x2.size());  
  if (yy.size()!=x1.size() || yy[0].size()!=x2.size()) ErrorMsg("Error in interpolated_2D of Func.cpp!"); 

  for (unsigned int i=0; i<x1.size(); i++) X1[i] = x1[i];
  for (unsigned int i=0; i<x2.size(); i++) X2[i] = x2[i];
  for (unsigned int i=0; i<x1.size(); i++) 
    for (unsigned int j=0; j<x2.size(); j++) 
      YY[i][j] = yy[i][j];

  if (type=="Linear") {
    Bilin_interp func(X1, X2, YY);
    return func.interp(_x1,_x2);
  }

  else if (type=="Poly") {
    Poly2D_interp func(X1, X2, YY, nPt, nPt);
    return func.interp(_x1,_x2);
  }

  else if (type=="Spline") {
    Spline2D_interp func(X1, X2, YY);
    return func.interp(_x1,_x2); 
  }

  else 
    ErrorMsg("Error in interpolated_2D of Func.cpp: the value of string 'type' is not permitted!");

  return -1;
}


// ============================================================================


void cosmobl::checkIO (string file, bool isInput)
{
  fstream ff;

  if (isInput) {
    ff.open (file.c_str(), fstream::in); 
    if (!ff.is_open()) { string Err = "Error in opening the input file: " + file + "!"; ErrorMsg(Err); }
  }
  else {
    ff.open (file.c_str(), fstream::out);
    if (!ff.is_open()) { string Err = "Error in opening the output file: "+ file + " ! "; ErrorMsg(Err); }
  }
  
  ff.clear(); ff.close();
}


// ============================================================================


void cosmobl::invert_matrix (vector<vector<double> > mat, vector<vector<double> > &mat_inv, double prec)
{
  int n = mat.size();
  cout << n << endl;
  int s;
  if (n==0)
    ErrorMsg("Error in invert_matrix of Func.cpp. 0 size for the input matrix");

  mat_inv.erase(mat_inv.begin(), mat_inv.end());
  mat_inv = mat;

  gsl_matrix *mm = gsl_matrix_alloc(n, n);
  gsl_matrix *im = gsl_matrix_alloc(n, n);
  gsl_permutation * perm = gsl_permutation_alloc(n);

  for (int i=0; i<n; i++)
    for (int j=0; j<n; j++)
      gsl_matrix_set(mm, i, j, mat[i][j]);

  // Make LU decomposition of matrix m
  gsl_linalg_LU_decomp (mm, perm, &s);

  // Invert the matrix m
  gsl_linalg_LU_invert (mm, perm, im);

  for (int i=0; i<n; i++)
    for (int j=0; j<n; j++)
      mat_inv[i][j] = gsl_matrix_get(im,i,j);

  for(int i=0; i<n; i++) {
    for(int j=0; j<n; j++) {
      double fact = (i==j) ? 1 : 0;
      double prod = 0;
      for (int el=0; el<n; el++)
	prod += mat[i][el]*mat_inv[el][j];
      
      if (fabs(fact - prod) > prec)  
	WarningMsg("Exceeded precision for element "+conv(i,par::fINT)+" "+conv(j,par::fINT)+"; "+conv(fact,par::fDP4)+" "+conv(prod,par::fDP4));
    }
  }
}


// ============================================================================


void cosmobl::invert_matrix (vector<vector<double> > mat, vector<vector<double> > &mat_inv, int i1, int i2, double prec)
{
  int n = i2-i1;
  int s;
  if (n==0)
    ErrorMsg("Error in invert_matrix of Func.cpp. 0 size for the input matrix");

  mat_inv.erase(mat_inv.begin(),mat_inv.end());
  mat_inv = mat;

  gsl_matrix *mm = gsl_matrix_alloc (n, n);
  gsl_matrix *im = gsl_matrix_alloc (n, n);
  gsl_permutation * perm = gsl_permutation_alloc (n);

  for(int i=i1;i<i2;i++)
    for(int j=i1;j<i2;j++)
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

  for (size_t i=0; i<mat.size(); i++){
    for (size_t j=0; j<mat[i].size(); j++){
      if (int(i)<i2 && int(i)>=i1 && int(j)<i2 && int(j)>=i1)
	mat_inv[i][j] = gsl_matrix_get(im,i-i1,j-i1);
      else
	mat_inv[i][j] = 0;
    }
  }
}


// ============================================================================


void cosmobl::invert_small_matrix (vector< vector<double> > mat, vector< vector<double> > &mat_inv, double prec)
{
  mat_inv.erase(mat_inv.begin(), mat_inv.end());
  vector<double> vv (mat.size(),0.);
  for (unsigned int i=0; i<mat.size(); i++) mat_inv.push_back(vv);
  
  MatDoub Mat(mat.size(),mat.size()), Mat_inv(mat.size(),mat.size()), ID(mat.size(),mat.size());
  
  for (unsigned int i=0; i<mat.size(); i++) {
    if (mat.size()!=mat[i].size()) ErrorMsg("Error in invert_small_matrix of Func.cpp! The matrix is not quadratic!");
    for (unsigned int j=0; j<mat.size(); j++) 
      Mat[i][j] = mat[i][j];
  }

  LUdcmp alu(Mat); 

  if (fabs(alu.det())<1.e-10) { 
    string Warn = "Attention: the determinant of the matrix is: "+ conv(alu.det(),par::fDP3) + "!";
    WarningMsg(Warn);
  }

  alu.inverse(Mat_inv); 

  ID = Mat*Mat_inv;

  for (unsigned int i=0; i<mat.size(); i++) 
    for (unsigned int j=0; j<mat.size(); j++) 
      if (i!=j && ID[i][j]>prec) 
	ErrorMsg("Error in function: invert_small_matrix of Func.cpp!"); 	
      else 
	mat_inv[i][j] = Mat_inv[i][j];     
}


// ============================================================================


void cosmobl::covariance_matrix (vector< vector<double> > mat, vector< vector<double> > &cov, bool JK) 
{  
  cov.erase(cov.begin(),cov.end());
  vector<double> vv (mat[0].size(),0.);
  for (unsigned int i=0; i<mat[0].size(); i++) cov.push_back(vv);
  

  // measure the mean values at each bin

  vector<double> mean, sigma;

  for (unsigned int i=0; i<mat[0].size(); i++) { // loop on the bin i
    
    vector<double> vect_temp;
    for (unsigned int j=0; j<mat.size(); j++) // loop on the realization j
      vect_temp.push_back(mat[j][i]);

    if (vect_temp.size()>2) {
      mean.push_back(Average(vect_temp));
      sigma.push_back(Sigma(vect_temp));
    }
    else 
      mean.push_back(-1.e30);
  }
  
  
  // compute the elements of the covariance matrix
  
  for (unsigned int i=0; i<cov.size(); i++) 
    for (unsigned int j=0; j<cov.size(); j++) {
      int nm = 0;
      for (unsigned int k=0; k<mat.size(); k++) {
	cov[i][j] += (mean[i]-mat[k][i])*(mean[j]-mat[k][j]);
	nm ++;
      }
      if (nm>1) cov[i][j] = (JK) ? double(nm-1)/(nm)*cov[i][j] : cov[i][j]/(nm-1);
    }

}


// ============================================================================


void cosmobl::covariance_matrix (vector<string> file, vector<double> &rad, vector<double> &mean, vector< vector<double> > &cov) 
{  
  rad.erase(rad.begin(),rad.end());
  mean.erase(mean.begin(),mean.end());
  cov.erase(cov.begin(),cov.end());
  int n_mocks = file.size(); 
  

  // count the number of bins

  ifstream fin_check(file[0].c_str()); checkIO (file[0],1);
  int dim = 0; string line;
  while (getline(fin_check,line)) dim++;
  fin_check.clear(); fin_check.close();

  vector<double> vv(dim,0.);
  for (int i=0; i<dim; i++) cov.push_back(vv);

  vector<double> vv2 (dim,0.);
  vector< vector<double> > mat_temp;
  for (int i=0; i<n_mocks; i++) mat_temp.push_back(vv2);


  // read the input files

  for (int ff=0; ff<n_mocks; ff++) {
    ifstream fin (file[ff].c_str()); checkIO (file[ff],1);

    int index = 0;
    while (getline(fin,line)) {

      stringstream ss(line);
      vector<double> num; double NN;
      while (ss>>NN) num.push_back(NN);

      if (ff==0) rad.push_back(num[0]);
      mat_temp[ff][index++] = num[1];
    }
    
    fin.clear(); fin.close();
  }

  
  // measure the mean values at each bin

  for (int i=0; i<dim; i++) {
    
    vector<double> vect_temp;
    for (int ff=0; ff<n_mocks; ff++) 
      if (mat_temp[ff][i]>-1.e20) 
	vect_temp.push_back(mat_temp[ff][i]);

    if (vect_temp.size()>2)
      mean.push_back(Average(vect_temp));
    else 
      mean.push_back(-1.e30);
  }

  
  // compute the elements of the covariance matrix
  
  for (int i=0; i<dim; i++) 
    for (int j=0; j<dim; j++) {
      int nm = 0;
      for (int ff=0; ff<n_mocks; ff++) 
	if (mat_temp[ff][i]>-1.e29 && mat_temp[ff][j]>-1.e29 && mean[i]>-1. && mean[j]>-1.) { // check!!!
	  cov[i][j] += (mean[i]-mat_temp[ff][i])*(mean[j]-mat_temp[ff][j]);
	  nm ++;
	}
      if (nm>1) cov[i][j] *= 1./(nm-1);

      if (fabs(cov[i][j])<1.e-30) {
	cov[i][j] = 1.e30;
	string Warn = "Attention: cov[" + conv(i,par::fINT) + "][" + conv(j,par::fINT) + "] = 1.e30! (in Func.cpp:covariance_matrix)";
	WarningMsg(Warn);
      }
    }

}


// ============================================================================


void cosmobl::Moment (vector<double> data, double *ave, double *adev, double *sdev, double *var, double *skew, double *curt) 
{
  VecDoub Data(data.size());
  for (unsigned int i=0; i<data.size(); i++) Data[i] = data[i];
  moment(Data, *ave, *adev, *sdev, *var, *skew, *curt);
}


// ============================================================================


double cosmobl::relative_error_beta (double &bias, double &Volume, double &density) // from Eq. 20 of Bianchi et al. 2012
{ 
  double n0 = 1.7e-4; // in (h/Mpc)^3
  double CC = 4.9e2;  // in (Mpc/h)^1.5

  return CC*pow(bias,0.7)/sqrt(Volume)*exp(n0/(bias*bias*density))*100.;
}


// ============================================================================


void cosmobl::measure_var_function (vector<double> var, int &bin, double &V_min, double &V_max, double &Volume, vector<double> &Var, vector<double> &Phi, vector<double> &err)
{
  if (var.size()==0) ErrorMsg("Error in measure_var_functions of Func.cpp: there are no objectes in the catalogue!");

  Var.erase(Var.begin(),Var.end());
  Phi.erase(Phi.begin(),Phi.end());
  err.erase(err.begin(),err.end());

  double delta_logV = (log10(V_max)-log10(V_min))/bin;
  double V1 = V_min;
  double V2 = V_min*pow(10.,delta_logV);

  for (int y=0; y<bin; y++) {
    double nHalo = 0.;

    for (unsigned int k=0; k<var.size(); k++) 
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


// ============================================================================


void cosmobl::quad_fit (vector<double> xx, vector<double> fx, vector<double> err, double *AA, double *BB, double *CC)
{
  double (*p_quad) (double, void *, vector<double>);
  p_quad = Pol2; 

  vector<double> params;
  vector< vector<double> > par_lim;

  cosmobl::classfunc::Chi2Sigma_1D func (xx, fx, err, p_quad, &params, par_lim);

  Powell<cosmobl::classfunc::Chi2Sigma_1D> powell (func,1.e-5);
  
  VecDoub par(3); par[0] = 1.; par[1] = 1.; par[2] = 1.;
  par = powell.minimize(par);

  *AA = par[0];
  *BB = par[1];
  *CC = par[2];
}


// ============================================================================


void cosmobl::gaussian_fit (vector<double> xx, vector<double> fx, double *mean, double *sigma)
{
  double (*p_gaussian) (double, void *, vector<double>);
  p_gaussian = gaussian; 

  vector<double> params;
  vector< vector<double> > par_lim;

  cosmobl::classfunc::Chi2Model_1D func (xx, fx, p_gaussian, &params, par_lim);

  Powell<cosmobl::classfunc::Chi2Model_1D> powell (func,1.e-5);
  
  VecDoub par(2); par[0] = 1.; par[1] = 1.;
  par = powell.minimize(par);
 
  *mean = par[0];
  *sigma = par[1];
  
}


// ============================================================================


double cosmobl::gaussian_convolution (vector<double> xx, vector<double> fx, double &mean, double &sigma)
{
  cosmobl::classfunc::func_conv_gauss func(xx, fx, mean, sigma);
  Midinf<cosmobl::classfunc::func_conv_gauss> qq(func, Min(xx), Max(xx));
  return qromo(qq);
}


// ============================================================================================


void cosmobl::set_EnvVar (vector<string> Var) 
{
  for (unsigned int vv=0; vv<Var.size(); vv++) 
    putenv(&Var[0][0]);
}


// ============================================================================================


void cosmobl::check_EnvVar (string Var) 
{
  string COM = "if [ $"+Var+" ]; then touch tmp; fi";
  if (system (COM.c_str())) {};
  ifstream fin_check("tmp");
  if (!fin_check) { 
    string Warn = "Attention: the variable " + Var + " has not been defined! (see check_EnvVar of Func.cpp)";
    WarningMsg(Warn);
  }
  fin_check.clear(); fin_check.close();
  if (system ("rm -f tmp")) {};
}


// ============================================================================================


int cosmobl::used_memory (int type)
{
#ifdef LINUX
  int memory = -1;

  string mem;
  if (type==1) mem = "VmRSS:";
  else if (type==2) mem = "VmSize:";
  else ErrorMsg("Error in cosmobl::used_memory of Func.cpp: the input value of type is not allowed!");
  
  string file = "/proc/self/status";
  ifstream fin (file.c_str()); checkIO(file, 1);
  
  string line, aa;
  while (getline(fin, line)) {
    stringstream ss(line);
    vector<string> val;
    while (ss>>aa) val.push_back(aa);
    if (val.size()==3 && val[0]==mem) {
      memory = atoi(val[1].c_str());
      break;
    }
  }
  
  fin.clear(); fin.close();
  
  return memory;

#else 
  WarningMsg("Attention: used_memory of Func.cpp works only on Linux systems");
  return 1;

#endif
  
}


// ============================================================================


int cosmobl::check_memory (double frac, bool exit, string func, int type)
{
#ifdef LINUX
  struct sysinfo memInfo;
  sysinfo (&memInfo);
  
  long long freePhysMem = memInfo.freeram;
  freePhysMem *= memInfo.mem_unit;
  int used = used_memory(type);
  
  if (used > freePhysMem*0.001*frac) { // 0.001 is to convert kbytes in bytes
    string Err = "Attention: possible memory problem";
    Err += (func.empty()) ? "!\n" : " in "+func+" !\n";
    Err += "freePhysMem = "+conv((double)freePhysMem*1.e-9, par::fDP3)+" GB\n";
    Err += "memory used by the process: = "+conv((double)used*1.e-6, par::fDP3)+" GB\n";
    if (exit) ErrorMsg(Err);
    else { WarningMsg(Err); return 0; }
  }
  return 1;
  
#else
  WarningMsg("Attention: check_memory of Func.cpp works only on Linux systems");
  return 1;
  
#endif
  
}

  
// ============================================================================

  
double cosmobl::D1 (double XX, vector<double> xx, vector<double> yy, string interpType, int Num, double stepsize)
{
  cosmobl::classfunc::func_grid Func (xx, yy, interpType, Num);
  double err = -1.;
  double D1 = dfridr(Func, XX, stepsize, err);
  if (err>D1*1.e-3) ErrorMsg("Error in D1 of Func.cpp!");
  return D1;
}


// ============================================================================


double cosmobl::D2 (double XX, vector<double> xx, vector<double> yy, string interpType, int Num, double stepsize)
{
  vector<double> dy;
  for (unsigned int i=0; i<xx.size(); i++)
    dy.push_back(D1(xx[i], xx, yy, interpType, Num, stepsize));
  
  cosmobl::classfunc::func_grid Func (xx, dy, interpType, Num);
  double err = -1.;
  double D2 = dfridr(Func, XX, stepsize, err);
  if (err>D2*1.e-3) ErrorMsg("Error in D2 of Func.h!");
  return D2;
}


// ============================================================================


double cosmobl::Deriv (int nd, double XX, vector<double> xx, vector<double> yy, string interpType, int Num, double stepsize)
{
  vector<double> dy;
  for (unsigned int i=0; i<xx.size(); i++) 
    if (nd-1>0) dy.push_back(Deriv(nd-1,xx[i],xx,yy,interpType,Num,stepsize));
    else dy.push_back(yy[i]);

  cosmobl::classfunc::func_grid Func (xx, dy, interpType, Num);
  double err = -1.;
  double DD = dfridr(Func, XX, stepsize, err);

  if (err>fabs(DD)*1.e-1) { 
    double errR = err/fabs(DD)*100;
    string Warn = "Attention: the error in the derivative is = " + conv(errR,par::fDP3) + " % !";
    WarningMsg(Warn);
  }

  return DD;
}


// ============================================================================


void cosmobl::convert_map_gnuplot_sm (string &file_gnu, string &file_sm)
{
  string dir = par::DirCosmo+"Func/";
  string file1 = dir+"file1";
  string file2 = dir+"file2";

  ifstream fin (file_gnu.c_str()); checkIO(file_gnu,1); 
  
  ofstream fout (file1.c_str()); checkIO(file1,0); 
  
  string line; 
  int ind_i = 0, ind_j = 0, dim2 = 0;
  bool count = 1;

  while (getline(fin,line)) {
    double XX = -1.e30, YY, FF;
    stringstream ss(line);
    ss >>XX>>YY>>FF;
    if (XX>-1.e10) {
      fout <<ind_i<<"   "<<ind_j++<<"   "<<FF<<endl;
      if (count) dim2 = ind_j;
    }
    else {ind_i ++; count = 0;}
  }
  fin.clear(); fin.close();

  int dim1 = ind_i;
  
  string CONV = dir+"conv "+conv(dim1,par::fINT)+" "+conv(dim2,par::fINT)+" "+dir;
  if (system (CONV.c_str())) {};
  string MV = "cp "+file2+" "+file_sm+"; rm -f "+file2+" "+file1;
  if (system (MV.c_str())) {};
  cout <<"I wrote the file: "<<file_sm<<endl;
  
}


// ============================================================================


void cosmobl::random_numbers (int nRan, int idum, vector<double> xx, vector<double> fx, vector<double> &nRandom, double n_min, double n_max)
{
  vector<double> fx_temp = fx;
  sort(fx_temp.begin(), fx_temp.end());
  double fx_max = fx_temp[fx_temp.size()-1];

  double xx_min = Min(xx)*0.999;
  double xx_max = Max(xx)*1.001;
  double delta_xx = xx_max-xx_min;
  double num1, num2;

  int step = 1000;
  double delta_bin = (xx_max-xx_min)/step;
  double XXB = xx_min;
  vector<double> xx_interp, fx_interp;
  
  for (int i=0; i<step; i++) {
    
    xx_interp.push_back(XXB);
    fx_interp.push_back(interpolated(XXB, xx, fx, "Poly", 4));
    
    XXB += delta_bin;
  }
 
  double xx_interp_min = Min(xx_interp)*1.0001;
  double xx_interp_max = Max(xx_interp)*0.9999;

  
  Ran ran(idum);
  double ddd, delta_bin_inv = 1./delta_bin;
 
  for (int i=0; i<nRan; i++) {
 
    int nTry = 0;

    do {

      nTry ++;
      num1 = min(ran.doub()*delta_xx+xx_min,xx_interp_max);
      
      num1 = max(xx_interp_min,num1);
      num2 = ran.doub()*fx_max;    
      ddd = int((num1-xx_min)*delta_bin_inv);      

    } while ((num2>fx_interp[ddd] || num1<n_min || num1>n_max) && nTry<100000);
    
    if (nTry==100000) ErrorMsg("Error in random_numbers of Func.cpp!");
 
    nRandom.push_back(num1);
  }
}


// ============================================================================================


void cosmobl::bin_function (string file_grid, double func(double, void *), void *par, int bin, double x_min, double x_max, string binning, vector<double> &xx, vector<double> &yy) 
{
  if (binning != "lin" && binning != "loglin" && binning != "log") 
    ErrorMsg("Error in bin_function of Func.cpp: binning can only be: lin, loglin or log !");

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
      string Err = "Error in bin_function of Func.cpp: xx.size()=" + conv(xx.size(),par::fINT) + " != bin=" + conv(bin,par::fINT) + "!";
      ErrorMsg(Err);
    }

  }

  else {
    
    cout <<"I'm creating the grid file: "<<file_grid<<"..."<<endl;

    fin.clear(); fin.close();

    if (binning != "lin") {
      if (x_min<0 || x_max<0) {
	string Err = "Error in bin_function of Func.cpp: x_min=" + conv(x_min,par::fDP3) + ", x_max=" + conv(x_max,par::fDP3) + "!";
	ErrorMsg(Err);
      }
    
      x_min = log10(x_min);
      x_max = log10(x_max);
    }

    xx = linear_bin_vector(bin, x_min, x_max);

    ofstream fout (file_grid.c_str()); checkIO(file_grid,0);

    for (unsigned int i=0; i<xx.size(); i++) {
      yy[i] = (binning=="lin") ? func(xx[i], par) : func(pow(10.,xx[i]), par);

      if (binning=="log") {
	if (yy[i]<0) ErrorMsg("Error in bin_function of Func.cpp: yy[i]<0!");
	else yy[i] = log10(yy[i]);
      }
      
      fout <<xx[i]<<"   "<<yy[i]<<endl; 
      cout <<xx[i]<<"   "<<yy[i]<<endl; 
    }
    fout.clear(); fout.close(); cout <<"I wrote the file: "<<file_grid<<endl;
  }

}


// ============================================================================================


void cosmobl::bin_function_2D (string file_grid, double func(double *, size_t, void *), void *par, int bin, double x1_min, double x1_max, double x2_min, double x2_max, string binning, vector<double> &xx1, vector<double> &xx2, vector< vector<double> > &yy) 
{
  if (binning != "lin" && binning != "loglin" && binning != "log") 
    ErrorMsg("Error in bin_function of Func.cpp: binning can only be: lin, loglin or log !");

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
      string Err = "Error in create_grid of Func.cpp: xx1.size()=" + conv(xx1.size(),par::fINT) + ", xx2.size()=" + conv(xx2.size(),par::fINT) + " != bin=" + conv(bin,par::fINT) + "!";  
      ErrorMsg(Err);
    }

  }

  else {
    
    cout <<"I'm creating the grid file: "<<file_grid<<"..."<<endl; 

    fin.clear(); fin.close();

    if (binning != "lin") {
      if (x1_min<0 || x1_max<0 || x2_min<0 || x2_max<0) {
	string Err = "Error in create_grid of Func.cpp: x1_min=" + conv(x1_min,par::fDP3) + ", x1_max=" + conv(x1_max,par::fDP3) + ", x2_min=" + conv(x2_min,par::fDP3) + ", x2_max=" + conv(x2_max,par::fDP3) + "!";
	ErrorMsg(Err);
      }
    
      x1_min = log10(x1_min);
      x1_max = log10(x1_max);
      x2_min = log10(x2_min);
      x2_max = log10(x2_max);
    }

    xx1 = linear_bin_vector(bin, x1_min, x1_max);
    xx2 = linear_bin_vector(bin, x2_min, x2_max);

    ofstream fout (file_grid.c_str()); checkIO(file_grid,0);
    double vec[2];

    for (int i=0; i<bin; i++) {
      for (int j=0; j<bin; j++) {
	
	if (binning=="lin") {vec[0] = xx1[i]; vec[1] = xx2[j];} 
	else {vec[0] = pow(10.,xx1[i]); vec[1] = pow(10.,xx2[j]);}

	double ff = (binning=="log") ? log10(func(vec, 2, par)) : func(vec, 2, par);
	yy[i].push_back(ff);

	fout <<xx1[i]<<"   "<<xx2[j]<<"   "<<yy[i][j]<<endl; 
	cout <<"--> "<<xx1[i]<<"   "<<xx2[j]<<"   "<<yy[i][j]<<endl;
      
      }
    }
    fout.clear(); fout.close(); cout <<"I wrote the file: "<<file_grid<<endl;
  }

}


// ============================================================================================

/// @cond glob

double cosmobl::func_grid_lin (double xx, void *params)
{
  struct cosmobl::glob::STR_grid *pp = (struct cosmobl::glob::STR_grid *) params;

  return interpolated(xx, pp->_xx, pp->_yy, "Linear", -1);
}


// ============================================================================================


double cosmobl::func_grid_loglin (double xx, void *params)
{
  struct cosmobl::glob::STR_grid *pp = (struct cosmobl::glob::STR_grid *) params;
 
  double lgx = log10(xx);

  return interpolated(lgx, pp->_xx, pp->_yy, "Linear", -1);
}


// ============================================================================================


double cosmobl::func_grid_log (double xx, void *params)
{
  struct cosmobl::glob::STR_grid *pp = (struct cosmobl::glob::STR_grid *) params;
 
  double lgx = log10(xx);

  return pow(10., interpolated(lgx, pp->_xx, pp->_yy, "Linear", -1));
}


// ============================================================================================


double cosmobl::func_grid_lin_2D (double *xx, __attribute__((unused)) size_t dim, void *params)
{
  struct cosmobl::glob::STR_grid_2D *pp = (struct cosmobl::glob::STR_grid_2D *) params;
   
  return interpolated_2D(xx[0], xx[1], pp->_xx1, pp->_xx2, pp->_yy, "Linear", -1);
}


// ============================================================================================


double cosmobl::func_grid_loglin_2D (double *xx, __attribute__((unused)) size_t dim, void *params)
{
  struct cosmobl::glob::STR_grid_2D *pp = (struct cosmobl::glob::STR_grid_2D *) params;
 
  double lgx1 = log10(xx[0]);
  double lgx2 = log10(xx[1]);

  return interpolated_2D(lgx1, lgx2, pp->_xx1, pp->_xx2, pp->_yy, "Linear", -1);
}


// ============================================================================================


double cosmobl::func_grid_log_2D (double *xx, __attribute__((unused)) size_t dim, void *params)
{
  struct cosmobl::glob::STR_grid_2D *pp = (struct cosmobl::glob::STR_grid_2D *) params;

  double lgx1 = log10(xx[0]);
  double lgx2 = log10(xx[1]);

  return pow(10., interpolated_2D(lgx1, lgx2, pp->_xx1, pp->_xx2, pp->_yy, "Linear", -1));
}

/// @endcond

// ============================================================================================

/* ======== Alfonso Veropalumbo ======== */

void cosmobl::eq2sdss (vector<double> &ra, vector<double> &dec, vector<double> &lambda, vector<double> &eta) 
{
  double SurveyCenterRa = 185.-90, SurveyCenterDec = 32.5;
  double d2r = par::pi/180.;
  for (unsigned int i=0; i<ra.size(); i++){
    double x = cos((ra[i]-SurveyCenterRa*d2r))*cos(dec[i]);
    double y = sin((ra[i]-SurveyCenterRa*d2r))*cos(dec[i]);
    double z = sin(dec[i]);
  
    lambda[i] = -asin(x)/d2r;
    eta[i] = atan2(z,y)/d2r - SurveyCenterDec; 
    if (eta[i] < -180.0) eta[i] += 360.0;
    if (eta[i] > 180.0)  eta[i] -= 360.0;
  }
}


// ============================================================================


void cosmobl::sdss2eq (vector<double> &lambda, vector<double> &eta,vector<double> &ra, vector<double> &dec) 
{
  double SurveyCenterRa = 185., SurveyCenterDec = 32.5;
  double d2r = par::pi/180.;
  for (unsigned int i=0; i<ra.size(); i++){
    double x =  -1.0*sin(lambda[i]*d2r);
    double y = cos(lambda[i]*d2r)*cos(eta[i]*d2r+SurveyCenterDec*d2r);
    double z = cos(lambda[i]*d2r)*sin(eta[i]*d2r+SurveyCenterDec*d2r);
  
    ra[i] = atan2(y,x)/d2r + SurveyCenterRa-90;
    if (ra[i] < 0){ ra[i] += 360;}
    dec[i] = asin(z)/d2r;
  }
}


// ============================================================================


void cosmobl::sdss_stripe(vector<double> &eta, vector<double> &lambda, vector<int> &stripe,vector<int> &str_u)
{
  double stripe_sep=2.5;
  double cen = 58.75;

  for (unsigned int i=0;i<eta.size();i++){

    if(lambda[i]<90.0){stripe[i] = (eta[i]+cen)/stripe_sep;}
    if(lambda[i]>90.0){stripe[i] = (eta[i]+cen+180.)/stripe_sep;}

    str_u[i] = stripe[i];

  }

  sort(str_u.begin(),str_u.end());
  vector<int>::iterator it = unique(str_u.begin(),str_u.end());
  str_u.resize(distance(str_u.begin(),it));

}


// ============================================================================


void cosmobl::fill_distr (vector<double> var, vector<double> prob, double ntot, vector<double> &out_Var)
{  
  int s1 = 1242;
  int s2 = -24323;
  
  Ran ran1(s1); Ran ran2(s2);
  
  double minVar = Min(var), maxVar = Max(var);
  double minP = Min(prob), maxP = Max(prob);
  double deltaVar = maxVar-minVar;
  double deltaP = maxP-minP;
  
  for (int i=0; i<ntot; i++) {
    bool go_on = 0;
    while (go_on==0) {
      double VV = minVar+deltaVar*ran1.doub();
      double pp = 2*deltaP*ran2.doub()+minP;
      double p2 = interpolated(VV, var, prob, "Spline", 3);
      if (abs(pp-p2)/p2<=0.01) { out_Var.push_back(VV); go_on = 1; }
    }
  
  }
  
}


// ============================================================================


void cosmobl::fill_distr (vector<double> &data_var, double &delta_Var, vector<double> &var, vector<double> &nObj, double &ntot, vector<double> &random_Var)
{
  double sum = 0.;
  
  for (unsigned int i=0; i<nObj.size(); i++) sum += nObj[i];
  for (unsigned int i=0; i<nObj.size(); i++) nObj[i] = ntot*nObj[i]/sum;
  
  vector<double> bin_limit;
  vector<double> occupation;
  
  double minVar = Min(data_var);
  double maxVar = Max(data_var);
  
  bin_limit.push_back(minVar);
  occupation.push_back(0);
  
  vector<double> index;
  
  for (unsigned int i=0; i<var.size(); i++)
    if (var[i] > minVar && var[i] < maxVar) 
      index.push_back(i);


  double V1 = minVar, V2, temp = 0;

  for (int i=Min(index); i<Max(index); i++) {
  
    V2 = V1+delta_Var;
    bin_limit.push_back(V2);
    V1 = V2;
    temp += nObj[i];
    occupation.push_back(temp);
  
  }
  
  bin_limit.push_back(maxVar);
  occupation.push_back(ntot);
  
  Ran ran(21314);
  random_Var.erase(random_Var.begin(),random_Var.end());
  
  for (unsigned int i=1; i<bin_limit.size(); i++) {
    for (int j=ceil(occupation[i-1]); j<ceil(occupation[i]); j++) {
      temp = (bin_limit[i]-bin_limit[i-1])*ran.doub()+bin_limit[i-1];
      random_Var.push_back(temp);
    }
  }
}


// ============================================================================


void cosmobl::fill_distr (int nRan, vector<double> &xx, vector<double> &fx, vector<double> &varRandom, double &xmin, double &xmax, int idum)
{
  Ran ran(idum);
  varRandom.erase(varRandom.begin(), varRandom.end());

  double minVar = xmin, maxVar = xmax;

  int sz = 100;
  vector<double> Fx(sz, 0.);
  vector<double> new_x = linear_bin_vector(sz, minVar, maxVar); 

  cosmobl::classfunc::func_grid Int (xx, fx, "Poly", 4);

  Midpnt<cosmobl::classfunc::func_grid> qtot(Int,minVar, maxVar);
  double NN = qromo(qtot); // normalization of f(x)

  for (int i=1; i<sz; i++) {
    Midpnt<cosmobl::classfunc::func_grid> qq(Int, minVar, new_x[i]);
    Fx[i] = qromo(qq)/NN;
  }

  for (int i=0; i<nRan; i++) {
    double pp = ran.doub();
    varRandom.push_back(interpolated(pp, Fx, new_x, "Poly", 4));
  }

}


// ============================================================================


void cosmobl::find_index (vector<double> xx, double x_min, double x_max, int &ind1, int &ind2)
{
  ind1 = xx.size();
  ind2 = 0;
  
  for (int i=0; i<int(xx.size()); i++) {
    if (x_min < xx[i] && xx[i] < x_max) {
      ind1 = min(i,ind1);
      ind2 = max(i,ind2);
    }
  }
}


// ============================================================================


void cosmobl::read_cov (string filecov, vector< vector<double> > &cov, vector<vector<double> > &cov_inv, int i1, int i2)
{
  int size=i2-i1+1;

  cov.erase(cov.begin(),cov.end());
  cov_inv.erase(cov_inv.begin(),cov_inv.end());

  ifstream fin (filecov.c_str());
  if (!fin) {
    string Warn = "Attention: the file " + filecov + " does not exist!";
    WarningMsg(Warn);
  }

  cov.erase(cov.begin(),cov.end());

  vector<double> vv ;
  cov.push_back(vv);
  string line; int i = 0;

  while(getline(fin,line)) {
    stringstream ss(line);
    vector<double> num; double NN = -1.e30;
    while (ss>>NN) num.push_back(NN);
    if (num.size()==3 && num[2]>-1.e29){
      cov[i].push_back(num[2]);
    }
    else {i++; cov.push_back(vv);}
  }

  cov.erase(cov.end()-1,cov.end());
  fin.clear(); fin.close();

  cov_inv=cov;
  vector<vector<double> > cov_lim(size,vector<double>(size,0)), cov_lim_inv;

  int tot_size=cov.size();

  for(int i=0;i<tot_size;i++){
    for(int j=0;j<tot_size;j++){
      if(i>i2 || i<i1 || j>i2 || j<i1){
      }
      else{  cov_lim[i-i1][j-i1]=cov[i][j]; }
    }
  }

  invert_small_matrix(cov_lim,cov_lim_inv);

  for(int i=0;i<tot_size;i++){
    for(int j=0;j<tot_size;j++){
      if(i>i2 || i<i1 || j>i2 || j<i1){
	cov_inv[i][j]=0.;
      }
      else {cov_inv[i][j]=cov_lim_inv[i-i1][j-i1];}
    }
  }

}


/*
  Get convolution of two function f1, f2 in output vector res. The two
  functions need to be defined on the same X-axis range, with equal
  points number DeltaX = (Xmax-Xmin)/n_x.
*/

void cosmobl::convolution (vector<double> f1, vector<double> f2, vector<double> &res, double deltaX)
{
  size_t nn = f1.size();
  if (nn!=f2.size()) ErrorMsg("Error in convolution of Func.cpp!"); // Arrays must have equal size 

  double *ff1 = new double[nn];
  double *ff2 = new double[nn];
  double *rr = new double[nn];

  for (unsigned int i=0; i<nn; i++) {
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


void cosmobl::distribution (vector<double> &xx, vector<double> &fx, vector<double> FF, vector<double> WW, int nbin, bool linear, string file_out, double fact, double V1, double V2, bool bin_type, bool conv, double sigma)
{
  if (xx.size()>0 || fx.size()>0 || FF.size()<=0 || nbin<=0) ErrorMsg("Error in distribution of Func.cpp!");

  ofstream fout;
  if (file_out!="NULL") { fout.open (file_out.c_str()); checkIO(file_out,0); }
  
  double minFF = (V1>-1.e29) ? V1 : Min(FF)*0.9999;
  double maxFF = (V2>-1.e29) ? V2 : Max(FF)*1.0001;

  // Using GSL to create histogram 

  gsl_histogram *histo = gsl_histogram_alloc(nbin);

  if (linear) gsl_histogram_set_ranges_uniform(histo, minFF, maxFF);

  else {
    vector<double> vv = logarithmic_bin_vector(nbin+1, minFF, maxFF);
    double *vvv = new double[nbin+1]; for (int i=0; i<nbin+1; i++) vvv[i] = vv[i];
    gsl_histogram_set_ranges(histo, vvv, nbin+1);
  }
  
  for (unsigned int i=0; i<FF.size(); i++) 
    gsl_histogram_accumulate (histo, FF[i], WW[i]);
  
  double x1, x2;

  for (int i=0; i<nbin; i++) {

    gsl_histogram_get_range(histo, i, &x1, &x2);
    double val = gsl_histogram_get(histo, i);
    
    if (linear) xx.push_back(0.5*(x1+x2));
    else xx.push_back(pow(10.,0.5*(log10(x1)+log10(x2))));

    if (bin_type) fx.push_back(val/((x2-x1)*fact));
    else fx.push_back(val/((log10(x2)-log10(x1))*fact));
  }
  
  if (file_out!="NULL") {
    for (size_t i=0; i<xx.size(); i++)
      fout << xx[i] << "   " << fx[i] << endl;
    fout.clear(); fout.close(); cout <<"I wrote the file: "<<file_out<<endl;
  }
  
  
  // in case of Gaussian convolution

  if (conv) {

    if (!linear) ErrorMsg("Work in progress...");
    
    double deltaX = (maxFF-minFF)/nbin;
    sigma *= deltaX;

    vector<double> res, xx1, ff1, fx2;

    // create an extra space on the two side of the array (smoothing on the borders)
    int extra = sigma/deltaX;

    for (int i=1; i<extra; i++) {
      xx1.push_back(xx[0]-(extra-i)*deltaX);
      ff1.push_back(fx[0]*double(i)/extra);
    }

    for (int i=0; i<nbin; i++) {
      xx1.push_back(xx[i]);
      ff1.push_back(fx[i]);
    }
 
    for (int i=1; i<extra; i++) {
      xx1.push_back(xx[nbin-1]+i*deltaX);
      ff1.push_back(double(extra-i)/extra*fx[nbin-1]);
    }
       

    // Define xmean (center of the xrange)
    double mean = (Max(xx1)+Min(xx1))*0.5; // gsl_histogram_mean(histo);

    void *pp = NULL; vector<double> pars(2); pars[0] = mean; pars[1] = sigma;

    for (unsigned int i=0; i<xx1.size(); i++) 
      fx2.push_back(gaussian(xx1[i], pp, pars));
   
    // convolve
    convolution(ff1, fx2, res, deltaX);

    // delete extra bins added
    for (int i=extra-1; i<int(res.size())-extra+1; i++)
      fx[i-extra+1] = res[i];
  }
}


// ============================================================================


/* ======== Cosimo Fedeli ======== */

/// @cond glob

void cosmobl::gauleg (const double x1, const double x2, double *x, double *w, const int n)
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
