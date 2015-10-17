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
 *  @file Headers/Lib/ModelTwoPointClassFunc.h
 *
 *  @brief Class functions used to model the two-point correlation
 *  function
 *
 *  This file contains some class functions used by Numerical methods
 *  to model the two-point correlation function
 *
 *  @author Federico Marulli
 *
 *  @author federico.marulli3@unbo.it
 */

#ifndef __MTWOCLASSFUNC__
#define __MTWOCLASSFUNC__ 


// =====================================================================================


namespace cosmobl {

  namespace classfunc {

    class func_chi2_xidpl
    {
    private:
      vector<double> r_log, xi_log, error_xi_log;
      int r_ind1, r_ind2;

    public:
      func_chi2_xidpl (vector<double> _r_log, vector<double> _xi_log, vector<double> _error_xi_log, int &_r_ind1, int &_r_ind2) 
	: r_log(_r_log), xi_log(_xi_log), error_xi_log(_error_xi_log), r_ind1(_r_ind1), r_ind2(_r_ind2) {}
  
      double operator() (VecDoub &XX) 
      {
	if (XX[0]<0. || XX[1]<0. || XX[2]<0.) return 1.e30;

	double chi2 = 0;
	for (int rr=r_ind1; rr<r_ind2; rr++) 
	  if (xi_log[rr]>-1.e29 && error_xi_log[rr]>0) 
	    chi2 += pow((xi_log[rr]-cosmobl::double_powerlaw(r_log[rr],XX[0],XX[1],XX[2]))/error_xi_log[rr],2);
    
	return chi2;    
      }
    };


    // =====================================================================================


    class func_chi2_xidpl_error_r0
    {
    private:
      vector<double> r_log, xi_log, error_xi_log;
      int r_ind1, r_ind2;
      double alpha, beta, chi2_min;

    public:
      func_chi2_xidpl_error_r0 (vector<double> _r_log, vector<double> _xi_log, vector<double> _error_xi_log, int &_r_ind1, int &_r_ind2, double _alpha, double _beta, double _chi2_min) 
	: r_log(_r_log), xi_log(_xi_log), error_xi_log(_error_xi_log), r_ind1(_r_ind1), r_ind2(_r_ind2), alpha(_alpha), beta(_beta), chi2_min(_chi2_min) {}
  
      double operator() (VecDoub &XX) 
      {
	if (XX[0]<0.) return 1.e30;

	double chi2 = 0;
	for (int rr=r_ind1; rr<r_ind2; rr++) 
	  if (xi_log[rr]>-1.e29 && error_xi_log[rr]>0) 
	    chi2 += pow((xi_log[rr]-cosmobl::double_powerlaw(r_log[rr],XX[0],alpha,beta))/error_xi_log[rr],2);
    
	return fabs((chi2-chi2_min)-1.);    
      }
    };


    // =====================================================================================


    class func_chi2_xidpl_error_alpha
    {
    private:
      vector<double> r_log, xi_log, error_xi_log;
      int r_ind1, r_ind2;
      double r0, beta, chi2_min;

    public:
      func_chi2_xidpl_error_alpha (vector<double> _r_log, vector<double> _xi_log, vector<double> _error_xi_log, int &_r_ind1, int &_r_ind2, double _r0, double _beta, double _chi2_min) 
	: r_log(_r_log), xi_log(_xi_log), error_xi_log(_error_xi_log), r_ind1(_r_ind1), r_ind2(_r_ind2), r0(_r0), beta(_beta), chi2_min(_chi2_min) {}
  
      double operator() (VecDoub &XX) 
      {
	if (XX[0]<0.) return 1.e30;

	double chi2 = 0;
	for (int rr=r_ind1; rr<r_ind2; rr++) 
	  if (xi_log[rr]>-1.e29 && error_xi_log[rr]>0) 
	    chi2 += pow((xi_log[rr]-cosmobl::double_powerlaw(r_log[rr],r0,XX[0],beta))/error_xi_log[rr],2);
    
	return fabs((chi2-chi2_min)-1.);    
      }
    };


    // =====================================================================================


    class func_chi2_xidpl_error_beta
    {
    private:
      vector<double> r_log, xi_log, error_xi_log;
      int r_ind1, r_ind2;
      double r0, alpha, chi2_min;

    public:
      func_chi2_xidpl_error_beta (vector<double> _r_log, vector<double> _xi_log, vector<double> _error_xi_log, int &_r_ind1, int &_r_ind2, double _r0, double _alpha, double _chi2_min) 
	: r_log(_r_log), xi_log(_xi_log), error_xi_log(_error_xi_log), r_ind1(_r_ind1), r_ind2(_r_ind2), r0(_r0), alpha(_alpha), chi2_min(_chi2_min) {}
  
      double operator() (VecDoub &XX) 
      {
	if (XX[0]<0.) return 1.e30;

	double chi2 = 0;
	for (int rr=r_ind1; rr<r_ind2; rr++) 
	  if (xi_log[rr]>-1.e29 && error_xi_log[rr]>0) 
	    chi2 += pow((xi_log[rr]-cosmobl::double_powerlaw(r_log[rr],r0,alpha,XX[0]))/error_xi_log[rr],2);
    
	return fabs((chi2-chi2_min)-1.);    
      }
    };


    // =====================================================================================


    class func_chi2cov_xidpl
    {
    private:
      vector<double> r_log, xi_log;
      vector< vector<double> > cov_inv;
      int r_ind1, r_ind2;

    public:
      func_chi2cov_xidpl (vector<double> _r_log, vector<double> _xi_log, vector< vector<double> > _cov_inv, int &_r_ind1, int &_r_ind2) 
	: r_log(_r_log), xi_log(_xi_log), cov_inv(_cov_inv), r_ind1(_r_ind1), r_ind2(_r_ind2) {}
  
      double operator() (VecDoub &XX) 
      {
	if (XX[0]<0. || XX[1]<0. || XX[2]<0.) return 1.e30;
    
	double chi2 = 0;
	for (int i=r_ind1; i<r_ind2; i++) 
	  for (int j=r_ind1; j<r_ind2; j++) 
	    if (xi_log[i]>-1.e29 && xi_log[j]>-1.e29)
	      chi2 += (cosmobl::double_powerlaw(r_log[i],XX[0],XX[1],XX[2])-xi_log[i])*cov_inv[i][j]*(cosmobl::double_powerlaw(r_log[j],XX[0],XX[1],XX[2])-xi_log[j]);
    
	if (chi2<0 || chi2!=chi2) cosmobl::ErrorMsg("Error in func_chi2cov_xidpl of ModelTwoPointClassFunc.h!");

	return chi2;    
      }
    };


    // =====================================================================================


    class func_chi2cov_xidpl_error_r0
    {
    private:
      vector<double> r_log, xi_log;
      vector< vector<double> > cov_inv;
      int r_ind1, r_ind2;
      double alpha, beta, chi2_min;

    public:
      func_chi2cov_xidpl_error_r0 (vector<double> _r_log, vector<double> _xi_log, vector< vector<double> > _cov_inv, int &_r_ind1, int &_r_ind2, double _alpha, double _beta, double _chi2_min) 
	: r_log(_r_log), xi_log(_xi_log), cov_inv(_cov_inv), r_ind1(_r_ind1), r_ind2(_r_ind2), alpha(_alpha), beta(_beta), chi2_min(_chi2_min) {}
  
      double operator() (VecDoub &XX) 
      {
	if (XX[0]<0.) return 1.e30;

	double chi2 = 0;
	for (int i=r_ind1; i<r_ind2; i++) 
	  for (int j=r_ind1; j<r_ind2; j++) 
	    if (xi_log[i]>-1.e29 && xi_log[j]>-1.e29)
	      chi2 += (cosmobl::double_powerlaw(r_log[i],XX[0],alpha,beta)-xi_log[i])*cov_inv[i][j]*(cosmobl::double_powerlaw(r_log[j],XX[0],alpha,beta)-xi_log[j]);

	if (chi2<0 || chi2!=chi2) cosmobl::ErrorMsg("Error in func_chi2cov_xidpl_error_r0 of ModelTwoPointClassFunc.h!");

	return fabs((chi2-chi2_min)-1.);    
      }
    };


    // =====================================================================================


    class func_chi2cov_xidpl_error_alpha
    {
    private:
      vector<double> r_log, xi_log;
      vector< vector<double> > cov_inv;
      int r_ind1, r_ind2;
      double r0, xi0, beta, chi2_min;

    public:
      func_chi2cov_xidpl_error_alpha (vector<double> _r_log, vector<double> _xi_log, vector< vector<double> > _cov_inv, int &_r_ind1, int &_r_ind2, double _r0, double _beta, double _chi2_min) 
	: r_log(_r_log), xi_log(_xi_log), cov_inv(_cov_inv), r_ind1(_r_ind1), r_ind2(_r_ind2), r0(_r0), beta(_beta), chi2_min(_chi2_min) {}
  
      double operator() (VecDoub &XX) 
      {
	if (XX[0]<0.) return 1.e30;

	double chi2 = 0;
	for (int i=r_ind1; i<r_ind2; i++) 
	  for (int j=r_ind1; j<r_ind2; j++) 
	    if (xi_log[i]>-1.e29 && xi_log[j]>-1.e29)
	      chi2 += (cosmobl::double_powerlaw(r_log[i],r0,XX[0],beta)-xi_log[i])*cov_inv[i][j]*(cosmobl::double_powerlaw(r_log[j],r0,XX[0],beta)-xi_log[j]);
    
	if (chi2<0 || chi2!=chi2) cosmobl::ErrorMsg("Error in func_chi2cov_xidpl_error_alpha of ModelTwoPointClassFunc.h!");

	return fabs((chi2-chi2_min)-1.);    
      }
    };


    // =====================================================================================


    class func_chi2cov_xidpl_error_beta
    {
    private:
      vector<double> r_log, xi_log;
      vector< vector<double> > cov_inv;
      int r_ind1, r_ind2;
      double r0, alpha, chi2_min;

    public:
      func_chi2cov_xidpl_error_beta (vector<double> _r_log, vector<double> _xi_log, vector< vector<double> > _cov_inv, int &_r_ind1, int &_r_ind2, double _r0, double _alpha, double _chi2_min) 
	: r_log(_r_log), xi_log(_xi_log), cov_inv(_cov_inv), r_ind1(_r_ind1), r_ind2(_r_ind2), r0(_r0), alpha(_alpha), chi2_min(_chi2_min) {}
  
      double operator() (VecDoub &XX) 
      {
	if (XX[0]<0.) return 1.e30;
    
	double chi2 = 0;
	for (int i=r_ind1; i<r_ind2; i++) 
	  for (int j=r_ind1; j<r_ind2; j++) 
	    if (xi_log[i]>-1.e29 && xi_log[j]>-1.e29)
	      chi2 += (cosmobl::double_powerlaw(r_log[i],r0,alpha,XX[0])-xi_log[i])*cov_inv[i][j]*(cosmobl::double_powerlaw(r_log[j],r0,alpha,XX[0])-xi_log[j]);
    
	if (chi2<0 || chi2!=chi2) cosmobl::ErrorMsg("Error in func_chi2cov_xidpl_error_beta of ModelTwoPointClassFunc.h!");

	return fabs((chi2-chi2_min)-1.);    
      }
    };


    // =====================================================================================


    class func_chi2_xiproj
    {
    private:
      vector<double> rp, xi_proj, error_xi_proj;
      int r_ind1, r_ind2;

    public:
      func_chi2_xiproj (vector<double> _rp, vector<double> _xi_proj, vector<double> _error_xi_proj, int &_r_ind1, int &_r_ind2) 
	: rp(_rp), xi_proj(_xi_proj), error_xi_proj(_error_xi_proj), r_ind1(_r_ind1), r_ind2(_r_ind2) {}
  
      double operator() (VecDoub &XX) 
      {
	if (XX[0]<0. || XX[1]<1.000001) return 1.e30;

	double chi2 = 0;
	for (int rr=r_ind1; rr<r_ind2; rr++) 
	  if (error_xi_proj[rr]>0)
	    chi2 += pow((xi_proj[rr]-cosmobl::xi_projected_powerlaw(rp[rr],XX[0],XX[1]))/error_xi_proj[rr],2);

	return chi2;    
      }
    };


    // =====================================================================================


    class func_chi2_xiproj_error_r0
    {
    private:
      vector<double> rp, xi_proj, error_xi_proj;
      int r_ind1, r_ind2;
      double gamma, chi2_min;

    public:
      func_chi2_xiproj_error_r0 (vector<double> _rp, vector<double> _xi_proj, vector<double> _error_xi_proj, int &_r_ind1, int &_r_ind2, double _gamma, double _chi2_min) 
	: rp(_rp), xi_proj(_xi_proj), error_xi_proj(_error_xi_proj), r_ind1(_r_ind1), r_ind2(_r_ind2), chi2_min(_chi2_min) 
      {
	if (_gamma<1) cosmobl::ErrorMsg("Error in func_chi2_xiproj_error_r0 of ModelTwoPointClassFunc.h!");
	else gamma = _gamma; 
      }
  
      double operator() (VecDoub &XX) 
      {
	if (XX[0]<0.) return 1.e30;

	double chi2 = 0;
	for (int rr=r_ind1; rr<r_ind2; rr++) 
	  if (error_xi_proj[rr]>0) 
	    chi2 += pow((xi_proj[rr]-cosmobl::xi_projected_powerlaw(rp[rr],XX[0],gamma))/error_xi_proj[rr],2);
    
	return fabs((chi2-chi2_min)-1.);    
      }
    };


    // =====================================================================================


    class func_chi2_xiproj_error_gamma
    {
    private:
      vector<double> rp, xi_proj, error_xi_proj;
      int r_ind1, r_ind2;
      double r0, chi2_min;

    public:
      func_chi2_xiproj_error_gamma (vector<double> _rp, vector<double> _xi_proj, vector<double> _error_xi_proj, int &_r_ind1, int &_r_ind2, double _r0, double _chi2_min) 
	: rp(_rp), xi_proj(_xi_proj), error_xi_proj(_error_xi_proj), r_ind1(_r_ind1), r_ind2(_r_ind2), r0(_r0), chi2_min(_chi2_min) {}
  
      double operator() (VecDoub &XX) 
      {
	if (XX[0]<1.) return 1.e30;

	double chi2 = 0;
	for (int rr=r_ind1; rr<r_ind2; rr++) 
	  if (error_xi_proj[rr]>0) 
	    chi2 += pow((xi_proj[rr]-cosmobl::xi_projected_powerlaw(rp[rr],r0,XX[0]))/error_xi_proj[rr],2);
    
	return fabs((chi2-chi2_min)-1.);    
      }
    };


    // =====================================================================================


    class func_chi2cov_xiproj
    {
    private:
      vector<double> rp, xi_proj;
      vector< vector<double> > cov_inv;
      int r_ind1, r_ind2;

    public:
      func_chi2cov_xiproj (vector<double> _rp, vector<double> _xi_proj, vector< vector<double> > _cov_inv, int &_r_ind1, int &_r_ind2) 
	: rp(_rp), xi_proj(_xi_proj), cov_inv(_cov_inv), r_ind1(_r_ind1), r_ind2(_r_ind2) {}
  
      double operator() (VecDoub &XX) 
      {
	if (XX[0]<0. || XX[1]<1.000001) return 1.e30;
    
	double chi2 = 0;
	for (int i=r_ind1; i<r_ind2; i++) 
	  for (int j=r_ind1; j<r_ind2; j++) 
	    if (xi_proj[i]>-1.e29 && xi_proj[j]>-1.e29)
	      chi2 += (cosmobl::xi_projected_powerlaw(rp[i],XX[0],XX[1])-xi_proj[i])*cov_inv[i][j]*(cosmobl::xi_projected_powerlaw(rp[j],XX[0],XX[1])-xi_proj[j]);

	if (chi2<0 || chi2!=chi2) cosmobl::ErrorMsg("Error in func_chi2cov_xiproj of ModelTwoPointClassFunc.h!");

	return chi2;    
      }
    };


    // =====================================================================================


    class func_chi2cov_xiproj_error_r0
    {
    private:
      vector<double> rp, xi_proj;
      vector< vector<double> > cov_inv;
      int r_ind1, r_ind2;
      double gamma, chi2_min;

    public:
      func_chi2cov_xiproj_error_r0 (vector<double> _rp, vector<double> _xi_proj, vector< vector<double> > _cov_inv, int &_r_ind1, int &_r_ind2, double _gamma, double _chi2_min) 
	: rp(_rp), xi_proj(_xi_proj), cov_inv(_cov_inv), r_ind1(_r_ind1), r_ind2(_r_ind2), chi2_min(_chi2_min) 
      {
	if (_gamma<1) cosmobl::ErrorMsg("Error in func_chi2_xiproj_error_r0 of ModelTwoPointClassFunc.h!");
	else gamma = _gamma; 
      }
  
      double operator() (VecDoub &XX) 
      {
	if (XX[0]<0.) return 1.e30;

	double chi2 = 0;
	for (int i=r_ind1; i<r_ind2; i++) 
	  for (int j=r_ind1; j<r_ind2; j++) 
	    if (xi_proj[i]>-1.e29 && xi_proj[j]>-1.e29)
	      chi2 += (cosmobl::xi_projected_powerlaw(rp[i],XX[0],gamma)-xi_proj[i])*cov_inv[i][j]*(cosmobl::xi_projected_powerlaw(rp[j],XX[0],gamma)-xi_proj[j]);
    
	if (chi2<0 || chi2!=chi2) cosmobl::ErrorMsg("Error in  of func_chi2_xiproj_error_r0 ModelTwoPointClassFunc.h!");
	return fabs((chi2-chi2_min)-1.);    
      }
    };


    // =====================================================================================


    class func_chi2cov_xiproj_error_gamma
    {
    private:
      vector<double> rp, xi_proj;
      vector< vector<double> > cov_inv;
      int r_ind1, r_ind2;
      double r0, chi2_min;

    public:
      func_chi2cov_xiproj_error_gamma (vector<double> _rp, vector<double> _xi_proj, vector< vector<double> > _cov_inv, int &_r_ind1, int &_r_ind2, double _r0, double _chi2_min) 
	: rp(_rp), xi_proj(_xi_proj), cov_inv(_cov_inv), r_ind1(_r_ind1), r_ind2(_r_ind2), r0(_r0), chi2_min(_chi2_min) {}
  
      double operator() (VecDoub &XX) 
      {
	if (XX[0]<1.) return 1.e30;

	double chi2 = 0;
	for (int i=r_ind1; i<r_ind2; i++) 
	  for (int j=r_ind1; j<r_ind2; j++)
	    if (xi_proj[i]>-1.e29 && xi_proj[j]>-1.e29)
	      chi2 += (cosmobl::xi_projected_powerlaw(rp[i],r0,XX[0])-xi_proj[i])*cov_inv[i][j]*(cosmobl::xi_projected_powerlaw(rp[j],r0,XX[0])-xi_proj[j]);
    
	if (chi2<0 || chi2!=chi2) cosmobl::ErrorMsg("Error in func_chi2cov_xiproj_error_gamma of ModelTwoPointClassFunc.h!");
	return fabs((chi2-chi2_min)-1.);    
      }
    };


    // =====================================================================================


    class func_bias
    {
    private:
      Cosmology *cosmology;
      double mean_bias;
      double Mmax, mean_redshift, DeltaR;
      string author_bias, author_MF, author_SS, Model;

    public:
      func_bias (Cosmology &_cosmology, double &_mean_bias, double &_Mmax, double &_mean_redshift, string &_author_bias, string &_author_MF, string &_author_SS, string &_Model, double &_DeltaR) 
	{
	  cosmology = new Cosmology; *cosmology = _cosmology;
	  mean_bias = _mean_bias; Mmax = _Mmax; mean_redshift = _mean_redshift; DeltaR = _DeltaR;
	  author_bias = _author_bias; author_MF = _author_MF; author_SS = _author_SS; Model = _Model;
	};
  
      ~func_bias () {delete cosmology;};
  
      double operator() (double lgMmin) 
      {
	double Mmin = pow(10.,lgMmin);

	double beff = cosmology->bias_eff(Mmin, Mmax, mean_redshift, author_bias, author_MF, author_SS, Model, DeltaR);
    
	return mean_bias-beff;
      }
    };

  }
}

#endif
