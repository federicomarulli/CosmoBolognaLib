/********************************************************************
 *  Copyright (C) 2014 by Federico Marulli and Alfonso Veropalumbo  *
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
 *  @file Chi2/Chi2.cpp
 *
 *  @brief Methods of the class Chi2 
 *
 *  This file contains the implementation of the methods of the class
 *  Chi2, used for &chi;&sup2; analyses
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#include "Chi2.h"
using namespace cosmobl;


// ======================================================================================


void cosmobl::Chi2::set_limits (double x_min, double x_max) 
{
  find_index (m_xx, x_min, x_max, m_min1, m_max1);
}


// ======================================================================================


void cosmobl::Chi2::set_limits (int m1, int m2) 
{
  m_min1 = m1;
  m_max1 = m2;
}


// ======================================================================================


void cosmobl::Chi2::set_limits (double x_min, double x_max, double y_min, double y_max) 
{
  find_index (m_xx, x_min, x_max, m_min1, m_max1);
  find_index (m_yy, y_min, y_max, m_min2, m_max2);
}


// ======================================================================================


void cosmobl::Chi2::set_limits (int m1, int m2, int m3, int m4) 
{
  m_min1 = m1;
  m_max1 = m2;
  m_min2 = m3;
  m_max2 = m4;
}


// ======================================================================================


double cosmobl::Chi2::get_chi2 (vector<double> pp) 
{
  VecDoub Par(pp.size());
  for (unsigned int i=0; i<pp.size(); i++) Par[i] = pp[i];

  if (m_dim==1) {

    if (m_fit_type==0) { // diagonal chi^2 
      cosmobl::classfunc::Chi2Sigma_1D fit (m_xx, m_fx1D, m_df1D, m_func1D, m_params, m_par_limits, m_ptype, m_min1, m_max1);
      return fit(Par);
    }

    else if (m_fit_type==1) { // diagonal log(chi^2)
      cosmobl::classfunc::Chi2log_1D fit (m_xx, m_fx1D, m_df1D, m_func1D, m_params, m_par_limits, m_ptype, m_min1, m_max1);
      return fit(Par);
    }

    else { // full covariance chi^2
      cosmobl::classfunc::Chi2Cov_1D fit (m_xx, m_fx1D, m_covinv, m_func1D, m_params, m_par_limits, m_ptype, m_min1, m_max1);
      return fit(Par);
    }
    
  }
  
  else if (m_dim==2) {

    if (m_fit_type==0) { // diagonal chi^2 
      cosmobl::classfunc::Chi2Sigma_2D fit (m_xx, m_yy, m_fx2D, m_df2D, m_func2D, m_params, m_par_limits, m_ptype, m_min1, m_min2,  m_max1, m_max2);
      return fit(Par);
    }
    
    else if (m_fit_type==1) { // diagonal log(chi^2)   
      cosmobl::classfunc::Chi2log_2D fit (m_xx, m_yy, m_fx2D, m_df2D, m_func2D, m_params, m_par_limits, m_ptype, m_min1, m_min2, m_max1, m_max2);
      return fit(Par);
    }

    else  // full covariance chi^2
      { ErrorMsg ("Work in progress... (see Chi2.cpp)"); return 0; }
  }

  else 
    { ErrorMsg ("Work in progress... (see Chi2.cpp)"); return 0; }
  
}


// ======================================================================================


void cosmobl::Chi2::get_bestfit (vector<double> &bf) 
{
  int npar = bf.size();
  VecDoub tpar(npar);
  for (int i=0; i<npar; i++) tpar[i] = bf[i];

  if (m_dim==1) {

    if (m_fit_type==0) { // diagonal chi^2 
      cosmobl::classfunc::Chi2Sigma_1D fit(m_xx, m_fx1D, m_df1D, m_func1D, m_params, m_par_limits, m_ptype, m_min1, m_max1);
      Powell<cosmobl::classfunc::Chi2Sigma_1D> pp(fit,1.e-4);
      tpar = pp.minimize(tpar);
    }

    else if (m_fit_type==1) { // diagonal log(chi^2)
      cosmobl::classfunc::Chi2log_1D fit(m_xx, m_fx1D, m_df1D, m_func1D, m_params, m_par_limits, m_ptype, m_min1, m_max1);
      Powell<cosmobl::classfunc::Chi2log_1D> pp(fit,1.e-4);
      tpar = pp.minimize(tpar);
    }

    else { // full covariance chi^2
      cosmobl::classfunc::Chi2Cov_1D fit(m_xx, m_fx1D, m_covinv, m_func1D, m_params, m_par_limits, m_ptype, m_min1, m_max1);
      Powell<cosmobl::classfunc::Chi2Cov_1D> pp(fit,1.e-4);
      tpar = pp.minimize(tpar);
    }

  }

  else if (m_dim==2) {
    if (m_fit_type==0) { // diagonal chi^2 
      cosmobl::classfunc::Chi2Sigma_2D fit(m_xx, m_yy, m_fx2D, m_df2D, m_func2D, m_params, m_par_limits, m_ptype, m_min1, m_min2, m_max1, m_max2);
      Powell<cosmobl::classfunc::Chi2Sigma_2D> pp(fit,1.e-4);
      tpar = pp.minimize(tpar);
    }
    
    else if (m_fit_type==1) { // diagonal log(chi^2)
      cosmobl::classfunc::Chi2log_2D fit(m_xx, m_yy, m_fx2D, m_df2D, m_func2D, m_params, m_par_limits, m_ptype, m_min1, m_min2, m_max1, m_max2);
      Powell<cosmobl::classfunc::Chi2log_2D> pp(fit,1.e-4);
      tpar = pp.minimize(tpar);
    }

    else // full covariance chi^2
      ErrorMsg("Work in progress... (see Chi2.cpp)");
  }
  
  m_bestfit.erase(m_bestfit.begin(),m_bestfit.end());
  for (int i=0; i<npar; i++) { m_bestfit.push_back(tpar[i]); bf[i] = tpar[i]; }
}


// ======================================================================================

// Inverting this operation ->  n=sum_{j=0}^{j=npar} par_j*pow(size,npar-1-j)

void cosmobl::Chi2::decompose_index (int nn, int npar, int size, vector<int> &inds)
{ 
  inds.erase(inds.begin(),inds.end()); 
  inds.resize(npar);

  for (int i=0; i<npar; i++) { 
    int kk = 0;
    for (int j=0; j<i; j++) 
      kk += inds[j]*pow(size,(npar-1-j));
    inds[i] = (nn-kk)*pow(size,-(npar-1-i));  
  }
}


// ======================================================================================

// Create grid of Chi^2

void cosmobl::Chi2::create_chi2grid (int size, vector<double> &bf, string fileout)
{
  int npar = bf.size();

  m_grid_par.erase(m_grid_par.begin(), m_grid_par.end());
  m_chi2_grid.erase(m_chi2_grid.begin(), m_chi2_grid.end());

  m_bestfit.erase(m_bestfit.begin(), m_bestfit.end());
  m_bestfit.resize(npar);
  
  for (int i=0; i<npar; i++) {
     vector<double> vv = linear_bin_vector(size, m_par_limits[i][0], m_par_limits[i][1]);
    m_grid_par.push_back(vv);
  }

  int tot_step = pow(size,npar);
  int nn = 0;

  double chi2 = 1.e30;
  ofstream fout(fileout.c_str()); checkIO(fileout,0);

  while(nn<tot_step){

    vector<int> inds;
    decompose_index(nn,npar,size,inds); 

    vector<double> Par(npar);

    for (int i=0; i<npar; i++) Par[i] = m_grid_par[i][inds[i]];
    double temp = get_chi2(Par);

    m_chi2_grid.push_back(temp);
      
    for (int i=0; i<npar; i++) fout << Par[i] << " " ; 
    fout << temp << endl;

    if (temp<=chi2) { // get the parameters that give the best chi square in the grid
      for (int i=0; i<npar; i++) m_bestfit[i]=Par[i];  
      chi2 = temp;
    } 

    nn ++;
  }

  cout <<"----> chi2 = "<<chi2<<endl;
  bf.erase(bf.begin(),bf.end()); bf.resize(npar);
  
  for (int i=0; i<npar; i++) { bf[i] = m_bestfit[i]; cout << bf[i] << " " ; } cout << endl;
}


// ======================================================================================


void cosmobl::Chi2::likelihood_normalization (int npar, int size)
{
  int tot_step=pow(size,npar);

  int n=npar-1;
  vector<double> ff = m_chi2_grid, fff;

  for(int i=0; i<tot_step; i++)
    ff[i] = exp(-0.5*ff[i]);

  while(n>0){

    int tt=tot_step*pow(size,-(npar-n));
    for(int i=0;i<tt;i++){
      double imin=i*size; double imax=(i+1)*size;

      vector<double> temp;
      for (int j=imin; j<imax; j++)
	temp.push_back(ff[j]);
	 
      cosmobl::classfunc::func_grid ii(m_grid_par[n],temp,"Spline",-1);
      Midpnt<cosmobl::classfunc::func_grid> q(ii,Min(m_grid_par[n]),Max(m_grid_par[n]));
      fff.push_back(qromo(q));
    }

    ff.erase(ff.begin(),ff.end());
    for (unsigned int i=0; i<fff.size(); i++) ff.push_back(fff[i]);
    fff.erase(fff.begin(), fff.end());

    n -= 1;
  }

  cosmobl::classfunc::func_grid ii(m_grid_par[0], ff, "Spline", -1);
  Midpnt<cosmobl::classfunc::func_grid> q(ii, Min(m_grid_par[0]), Max(m_grid_par[0]));

  m_lnorm = qromo(q);
}


// ======================================================================================


void cosmobl::Chi2::single_par_pdf (int par_num, vector<double> &par, vector<double> &fpar, bool Likelihood, bool norm)
{
  par.erase(par.begin(),par.end());
  fpar.erase(fpar.begin(),fpar.end());

  vector<double> grid_copy = m_chi2_grid;
  vector <vector<double> > grid_par_copy=m_grid_par;
  int npar=m_grid_par.size();
  int tot_step=m_chi2_grid.size();
  int size=ceil(pow(tot_step,1./npar));

  if (norm) likelihood_normalization(npar, size);

  if (m_bestfit.size()==0) { vector<double> bf(npar, 0); get_bestfit(bf); }
   
  //RE-ORDER, par_num as the first

  vector<double> p = grid_par_copy[0];
  grid_par_copy[0] = grid_par_copy[par_num];
  grid_par_copy[par_num] = p;

  int nn = 0;
  while (nn<tot_step) {
    vector<int> inds;
    decompose_index(nn,npar,size,inds); 

    int n2=nn+pow(size,npar-1)*(inds[par_num]-inds[0])-pow(size,npar-1-par_num)*(inds[par_num]-inds[0]);
    grid_copy[nn]= exp(-0.5*m_chi2_grid[n2])/m_lnorm;
    nn++;
  }

  int tt = tot_step;
  int nnn = npar-1;
  vector<double> ff = grid_copy, fff;

  while (nnn>0) {

    tt = tot_step*pow(size,-(npar-nnn));
    
    for (int i=0; i<tt; i++) {
      double imin = i*size; double imax = (i+1)*size;

      vector<double> temp;
      for (int j=imin; j<imax; j++)
	temp.push_back(ff[j]);
      
      cosmobl::classfunc::func_grid ii(grid_par_copy[nnn], temp, "Spline", -1);
      Midpnt<cosmobl::classfunc::func_grid> q(ii, Min(grid_par_copy[nnn]), Max(grid_par_copy[nnn]));
      fff.push_back(qromo(q));
    }

    ff.erase(ff.begin(), ff.end());
    for (unsigned int i=0; i<fff.size(); i++) ff.push_back(fff[i]);
    fff.erase(fff.begin(),fff.end());

    nnn --;
  }

  for (unsigned int i=0; i<grid_par_copy[0].size(); i++) {
    par.push_back(grid_par_copy[0][i]);
    fpar.push_back((!Likelihood) ? -2*log(ff[i]) : ff[i]);
  }

  if (!Likelihood) {
    cosmobl::classfunc::func_grid_minimum_1D func(par,fpar,"Spline",1);
    Powell<cosmobl::classfunc::func_grid_minimum_1D> ii(func,1.e-8);
    VecDoub XX(1); XX[0]=m_bestfit[par_num];
    XX=ii.minimize(XX);
    double min=func(XX);
    for (unsigned int i=0; i<grid_par_copy[0].size(); i++)
      fpar[i] -= min;
  }
}


// ======================================================================================


void cosmobl::Chi2::contour_2par (int p1, int p2, string out_file, bool norm)
{
  if (p1==p2) ErrorMsg("Error in cosmobl::Chi2::contour_2par of Chi2.cpp: parameters must be different!");
  if (p1>p2) { int bb = p2; p2 = p1; p1 = bb; }

  vector<double> grid_copy = m_chi2_grid;
  vector<vector<double>> grid_par_copy = m_grid_par;

  int npar = int(m_grid_par.size());
  int tot_step = int(m_chi2_grid.size());
  int size = ceil(pow(tot_step,1./npar));

  if (norm) likelihood_normalization(npar, size);

  vector<int> neworder(npar); 
  for (int i=0; i<npar; i++) neworder[i] = i;

  neworder[0] = p1; neworder[p1] = 0;
  int buff = neworder[1];
  neworder[1] = p2; neworder[p2] = buff;
  
  for (int i=0; i<npar; i++)
    grid_par_copy[i] = m_grid_par[neworder[i]];
  
  int nn = 0;

  while(nn<tot_step){
    vector<int> inds;

    decompose_index(nn, npar,size,inds);

    int n2 = 0;
    for (int i=0; i<npar; i++)
      n2 += inds[neworder[i]]*pow(size,npar-1-i);
      
    grid_copy[nn] = exp(-0.5*m_chi2_grid[n2])/m_lnorm;
    nn ++;
  }

  int tt = tot_step;
  int nnn = npar-1;
  vector<double> ff = grid_copy, fff;

  while (nnn>1) {

    tt = tot_step*pow(size,-(npar-nnn));
    for (int i=0; i<tt; i++) {
      double imin = i*size; double imax = (i+1)*size;

      vector<double> temp;
      for (int j=imin; j<imax; j++)
	temp.push_back(ff[j]);
      
      cosmobl::classfunc::func_grid ii(grid_par_copy[nnn], temp, "Spline", -1);
      Midpnt<cosmobl::classfunc::func_grid> qq(ii, Min(grid_par_copy[nnn]), Max(grid_par_copy[nnn]));
      fff.push_back(qromo(qq));
    }

    ff.erase(ff.begin(), ff.end());
    for (unsigned int i=0; i<fff.size(); i++) ff.push_back(fff[i]);
    fff.erase(fff.begin(), fff.end());

    nnn --;
  }
  double min;
  vector<vector<double> > grid_2D;

  for (unsigned int i=0; i<grid_par_copy[0].size(); i++) {
    vector<double> temp;
    for (unsigned int j=0; j<grid_par_copy[1].size(); j++) 
      temp.push_back(-2*log(ff[i*size+j]));
    grid_2D.push_back(temp);
  }

  cosmobl::classfunc::func_grid_minimum_2D func(grid_par_copy[0], grid_par_copy[1], grid_2D, "Linear", 1);
  Powell<cosmobl::classfunc::func_grid_minimum_2D> ii(func, 1.e-5);
  VecDoub XX(2); XX[0] = m_bestfit[p1]; XX[1] = m_bestfit[p2];

  cout << XX[0] << " " << XX[1] << endl;
  XX = ii.minimize(XX);
  min = func(XX);
  cout << XX[0] << " " << XX[1] << endl;

  ofstream fout(out_file.c_str()); checkIO(out_file,0);
  for (unsigned int i=0; i<grid_par_copy[0].size(); i++) {
    for (unsigned int j=0; j<grid_par_copy[1].size(); j++)
      fout << grid_par_copy[0][i] << " " << grid_par_copy[1][j] << " " << ( (isfinite(grid_2D[i][j])) ? grid_2D[i][j]-min : 1.e30) << endl;
    fout << endl;
  }
  fout.clear(); fout.close();

}
