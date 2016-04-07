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
 *  @file Headers/Lib/Chi2ClassFunc.h
 *
 *  @brief Class functions for &chi;&sup2; analyses used by Numerical
 *  methods
 *
 *  This file contains some class functions used by Numerical methods
 *  for &chi;&sup2; analyses
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#ifndef __CHI2CLASSFUNC__
#define __CHI2CLASSFUNC__ 


// =====================================================================================

// ------------------------------------------------------------------
// ----- chi2 class function, with no errors, for a 1D function -----
// ------------------------------------------------------------------

namespace cosmobl {

  namespace classfunc {

    class Chi2Model_1D 
    {
  
    private: 

      typedef double(*model) (double, shared_ptr<void>, vector<double>); // 1D function: y (x, fixed parameters, free parameters)  
  
      // ----- data -----
      vector<double> xx, fx;
  
      // ----- model, parameters, parameter limits -----
      model func;                          // model
      shared_ptr<void> params;                        // fixed parameters of the model
      vector< vector<double> > par_limits; // prior limits on the free parameters of the model
      bool ptype;                          // ptype = 0 -> fitPar[fitPar.size()-1] = i ; ptype = 1 -> fitPar[fitPar.size()-1] = index++

      // ----- limits of the fitting region -----
      int min, max;


    public:

      Chi2Model_1D (vector<double> _xx, vector<double> _fx, model _func, shared_ptr<void> _params, vector< vector<double> > _par_limits, bool _ptype=1, int _min=0, int _max=0)
	: xx(_xx), fx(_fx), func(_func), params(_params), par_limits(_par_limits), ptype(_ptype), min(_min), max((_max>0) ? _max : xx.size()) {} 
  
      double operator() (double Par)
      {  
	if (par_limits.size()>0)
	  if (Par < par_limits[0][0] || Par > par_limits[0][1]) return 1.e30;
    
	vector<double> fitPar(2,Par);
 
	double chi2 = 0;
	int index = 0;

	for (int i=min; i<max; i++) {
	  fitPar[fitPar.size()-1] = (ptype) ? index++ : i;
	  chi2 += pow((fx[i]-func(xx[i],params,fitPar))/func(xx[i],params,fitPar),2);
	}

	return chi2;
      }


      double operator() (VecDoub Par)
      {  
	if (Par.size()!=int(par_limits.size())) 
	  ErrorMsg("Error in Chi2ClassFunc.h!");

	if (par_limits.size()>0)
	  for (int i=0; i<Par.size(); i++) 
	    if (Par[i] < par_limits[i][0] || Par[i] > par_limits[i][1]) return 1.e30;
    
	vector<double> fitPar(Par.size()+1);
	for (int i=0; i<Par.size(); i++) fitPar[i] = Par[i];
  
	double chi2 = 0;
	int index = 0;

	for (int i=min; i<max; i++) {
	  fitPar[fitPar.size()-1] = (ptype) ? index++ : i;
	  chi2 += pow((fx[i]-func(xx[i],params,fitPar))/func(xx[i],params,fitPar),2);
	}

	return chi2;
      }

    };



    // ============================================================================

    // -----------------------------------------------------------------------------
    // ----- chi2 class function, with only diagonal errors, for a 1D function -----
    // -----------------------------------------------------------------------------

    class Chi2Sigma_1D 
    {

    private: 

      typedef double(*model) (double, shared_ptr<void> , vector<double>); // 1D function: y (x, fixed parameters, free parameters)  
  
      // ----- data -----
      vector<double> xx, fx, df;

      // ----- model, parameters, parameters limits -----
      model func;                          // model
      shared_ptr<void> params;                        // fixed parameters of the model
      vector< vector<double> > par_limits; // prior limits on the free parameters of the model
      bool ptype;                          // ptype = 0 -> fitPar[fitPar.size()-1] = i ; ptype = 1 -> fitPar[fitPar.size()-1] = index++

      // ----- limits of the fitting region -----
      int min, max;


    public:

      Chi2Sigma_1D (vector<double> _xx, vector<double> _fx, vector<double> _df, model _func, shared_ptr<void> _params, vector< vector<double> > _par_limits, bool _ptype=1, int _min=0, int _max=0)
	: xx(_xx), fx(_fx), df(_df), func(_func), params(_params), par_limits(_par_limits), ptype(_ptype), min(_min), max((_max>0) ? _max : xx.size()) {} 
    
      double operator() (double Par)
      {  
	if (par_limits.size()>0)
	  if (Par < par_limits[0][0] || Par > par_limits[0][1]) return 1.e30;
    
	vector<double> fitPar(2,Par);

	double chi2 = 0.;
	int index = 0;

	for (int i=min; i<max; i++) {
	  fitPar[fitPar.size()-1] = (ptype) ? index++ : i;
	  chi2 += pow((fx[i]-func(xx[i],params,fitPar))/df[i],2);
	}

	return chi2;
      }


      double operator() (VecDoub Par)
      {  
	if (Par.size()!=int(par_limits.size())) 
	  ErrorMsg("Error in Chi2ClassFunc.h!");

	if (par_limits.size()>0)
	  for (int i=0; i<Par.size(); i++) 
	    if (Par[i] < par_limits[i][0] || Par[i] > par_limits[i][1]) return 1.e30;

	vector<double> fitPar(Par.size()+1);
	for (int i=0; i<Par.size(); i++) fitPar[i] = Par[i];

	double chi2 = 0.;
	int index = 0;

	for (int i=min; i<max; i++) {
	  fitPar[fitPar.size()-1] = (ptype) ? index++ : i;
	  chi2 += pow((fx[i]-func(xx[i],params,fitPar))/df[i],2);
	}

	return chi2;
      }

    };


    // ============================================================================

    // ----------------------------------------------------------------------------------
    // ----- log(chi2) class function, with only diagonal errors, for a 1D function -----
    // ----------------------------------------------------------------------------------

    class Chi2log_1D 
    {

    private: 

      typedef double(*model) (double, shared_ptr<void> , vector<double>); // 1D function: y (x, fixed parameters, free parameters)  
  
      // ----- data -----
      vector<double> xx, fx, df;

      // ----- model, parameters, parameters limits -----
      model func;                          // model
      shared_ptr<void> params;                        // fixed parameters of the model
      vector< vector<double> > par_limits; // prior limits on the free parameters of the model
      bool ptype;                          // ptype = 0 -> fitPar[fitPar.size()-1] = i ; ptype = 1 -> fitPar[fitPar.size()-1] = index++

      // ----- limits of the fitting region -----
      int min, max;


    public:

      Chi2log_1D (vector<double> _xx, vector<double> _fx, vector<double> _df, model _func, shared_ptr<void> _params, vector< vector<double> > _par_limits, bool _ptype=1, int _min=0, int _max=0)
	: xx(_xx), fx(_fx), df(_df), func(_func), params(_params), par_limits(_par_limits), ptype(_ptype), min(_min), max((_max>0) ? _max : xx.size()) {} 
    
      double operator() (double Par)
      {  
	if (par_limits.size()>0)
	  if (Par < par_limits[0][0] || Par > par_limits[0][1]) return 1.e30;
    
	vector<double> fitPar(2,Par);

	double chi2 = 0, num = -1., den = -1.;
	int index = 0;

	for (int i=min; i<max; i++) {
	  fitPar[fitPar.size()-1] = (ptype) ? index++ : i;
	  num = log(1.+fx[i])-log(1.+func(xx[i],params,fitPar));
	  den = log(1.+fx[i]+df[i])-log(1.+fx[i]-df[i]);
	  if (den>0) chi2 += pow(num/den,2);
	}

	return chi2;
      }


      double operator() (VecDoub Par)
      {  
	if (Par.size()!=int(par_limits.size())) 
	  ErrorMsg("Error in Chi2ClassFunc.h!");

	if (par_limits.size()>0)
	  for (int i=0; i<Par.size(); i++) 
	    if (Par[i] < par_limits[i][0] || Par[i] > par_limits[i][1]) return 1.e30;
    
	vector<double> fitPar(Par.size()+1);
	for (int i=0; i<Par.size(); i++) fitPar[i] = Par[i];
    
	double chi2 = 0, num = -1., den = -1.;
	int index = 0;

	for (int i=min; i<max; i++) {
	  fitPar[fitPar.size()-1] = (ptype) ? index++ : i;
	  num = log(1.+fx[i])-log(1.+func(xx[i],params,fitPar));
	  den = log(1.+fx[i]+df[i])-log(1.+fx[i]-df[i]);
	  if (den>0) chi2 += pow(num/den,2);
	}

	return chi2;
      }

    };


    // ============================================================================

    // -----------------------------------------------------------------------------------
    // ----- chi2 class function, with the full covariance matrix, for a 1D function -----
    // -----------------------------------------------------------------------------------

    class Chi2Cov_1D 
    {

    private: 

      typedef double(*model) (double, shared_ptr<void> , vector<double>); // 1D function: y (x, fixed parameters, free parameters)  
  
      // ----- data -----
      vector<double> xx, fx;
      vector< vector<double> > cov_inv;

      // ----- model, parameters, parameters limits -----
      model func;                          // model
      shared_ptr<void> params;                        // fixed parameters of the model
      vector< vector<double> > par_limits; // prior limits on the free parameters of the model
      bool ptype;                          // ptype = 0 -> fitPar[fitPar.size()-1] = i ; ptype = 1 -> fitPar[fitPar.size()-1] = index++

      // ----- limits in the fitted function -----
      int min, max;


    public:

      Chi2Cov_1D (vector<double> _xx, vector<double> _fx, vector< vector<double> > _cov_inv, model _func, shared_ptr<void> _params, vector< vector<double> > _par_limits, bool _ptype=1, int _min=0, int _max=0)
	: xx(_xx), fx(_fx), cov_inv(_cov_inv), func(_func), params(_params), par_limits(_par_limits), ptype(_ptype), min(_min), max((_max>0) ? _max : xx.size()) {} 
    
      double operator() (double Par)
      {
	if (Par < par_limits[0][0] || Par > par_limits[0][1]) return 1.e30;
    
	vector<double> fitPar(2,Par);
 
	double chi2 = 0;
	int index = 0;

	for (int i=min; i<max; i++)
	  for (int j=min; j<max; j++) { 
	    fitPar[fitPar.size()-1] = (ptype) ? index++ : i;
	    chi2 += (func(xx[i],params,fitPar)-fx[i])*cov_inv[i][j]*(func(xx[j],params,fitPar)-fx[j]);
	  }

	return chi2;
      }


      double operator() (VecDoub Par)
      {  
	if (Par.size()!=int(par_limits.size())) 
	  ErrorMsg("Error in Chi2ClassFunc.h!");

	for (int i=0; i<Par.size(); i++) 
	  if (Par[i] < par_limits[i][0] || Par[i] > par_limits[i][1]) return 1.e30;
    
	vector<double> fitPar(Par.size()+1);
	for (int i=0; i<Par.size(); i++) fitPar[i] = Par[i];

	double chi2 = 0;
	int index = 0;

	for (int i=min; i<max; i++)
	  for (int j=min; j<max; j++) {
	    fitPar[fitPar.size()-1] = (ptype) ? index++ : i;
	    chi2 += (func(xx[i],params,fitPar)-fx[i])*cov_inv[i][j]*(func(xx[j],params,fitPar)-fx[j]);
	  }

	return chi2;
      }

    };


    // =====================================================================================

    // ------------------------------------------------------------------
    // ----- chi2 class function, with no errors, for a 2D function -----
    // ------------------------------------------------------------------

    class Chi2Model_2D 
    {
  
    private: 

      // 2D function: y (x, y, fixed parameters, free parameters)  
      typedef double(*model) (double, double, shared_ptr<void> , vector<double>); 
  
      // ----- data -----
      vector<double> xx, yy;
      vector< vector<double> > fx;
  
      // ----- model, parameters, parameter limits -----
      model func;                          // model
      shared_ptr<void> params;                        // fixed parameters of the model
      vector< vector<double> > par_limits; // prior limits on the free parameters of the model
      bool ptype;                          // ptype = 0 -> fitPar[fitPar.size()-1] = i ; ptype = 1 -> fitPar[fitPar.size()-1] = index++

      // ----- limits of the fitting region -----
      int min1, min2, max1, max2;


    public:

      Chi2Model_2D (vector<double> _xx, vector<double> _yy, vector< vector<double> > _fx, model _func, shared_ptr<void> _params, vector< vector<double> > _par_limits, int _min1=0, int _min2=0, int _max1=0, int _max2=0)
	: xx(_xx), yy(_yy), fx(_fx), func(_func), params(_params), par_limits(_par_limits), min1(_min1), min2(_min2), max1((_max1>0) ? _max1 : xx.size()), max2((_max2>0) ? _max2 : yy.size()) {} 
  

      double operator() (double Par)
      {  
	if (par_limits.size()>0)
	  if (Par < par_limits[0][0] || Par > par_limits[0][1]) return 1.e30;
    
	vector<double> fitPar(2,Par); 
    
	double chi2 = 0;
	int index = 0;

	for (int i=min1; i<max1; i++)
	  for (int j=min2; j<max2; j++) {
	    fitPar[fitPar.size()-1] = (ptype) ? index++ : i;
	    chi2 += pow((fx[i][j]-func(xx[i],yy[i],params,fitPar))/func(xx[i],yy[i],params,fitPar),2);
	  }

	return chi2;
      }


      double operator() (VecDoub Par)
      {  
	if (Par.size()!=int(par_limits.size())) 
	  ErrorMsg("Error in Chi2ClassFunc.h!");

	if (par_limits.size()>0)
	  for (int i=0; i<Par.size(); i++) 
	    if (Par[i] < par_limits[i][0] || Par[i] > par_limits[i][1]) return 1.e30;
    
	vector<double> fitPar(Par.size()+1);
	for (int i=0; i<Par.size(); i++) fitPar[i] = Par[i];

	double chi2 = 0;
	int index = 0;

	for (int i=min1; i<max1; i++)
	  for (int j=min2; j<max2; j++) {
	    fitPar[fitPar.size()-1] = (ptype) ? index++ : i;
	    chi2 += pow((fx[i][j]-func(xx[i],yy[i],params,fitPar))/func(xx[i],yy[i],params,fitPar),2);
	  }

	return chi2;
      }

    };


    // =====================================================================================

    // -----------------------------------------------------------------------------
    // ----- chi2 class function, with only diagonal errors, for a 2D function -----
    // -----------------------------------------------------------------------------

    class Chi2Sigma_2D 
    {
  
    private: 

      // 2D function: y (x, y, fixed parameters, free parameters)  
      typedef double(*model) (double, double, shared_ptr<void> , vector<double>); 
  
      // ----- data -----
      vector<double> xx, yy;
      vector< vector<double> > fx, df;
  
      // ----- model, parameters, parameter limits -----
      model func;                          // model
      shared_ptr<void> params;                        // fixed parameters of the model
      vector< vector<double> > par_limits; // prior limits on the free parameters of the model
      bool ptype;                          // ptype = 0 -> fitPar[fitPar.size()-1] = i ; ptype = 1 -> fitPar[fitPar.size()-1] = index++

      // ----- limits of the fitting region -----
      int min1, min2, max1, max2;


    public:

      Chi2Sigma_2D (vector<double> _xx, vector<double> _yy, vector< vector<double> > _fx, vector <vector<double> > _df, model _func, shared_ptr<void> _params, vector< vector<double> > _par_limits, bool _ptype=1, int _min1=0, int _min2=0, int _max1=0, int _max2=0)
	: xx(_xx), yy(_yy), fx(_fx), df(_df), func(_func), params(_params), par_limits(_par_limits), ptype(_ptype), min1(_min1), min2(_min2), max1((_max1>0) ? _max1 : _fx.size()), max2((_max2>0) ? _max2 : _fx[0].size()) {} 
  

      double operator() (double Par)
      {  
	if (par_limits.size()>0)
	  if (Par < par_limits[0][0] || Par > par_limits[0][1]) return 1.e30;
    
	vector<double> fitPar(2,Par); 

	double chi2 = 0;
	int index = 0;

	for (int i=min1; i<max1; i++)
	  for (int j=min2; j<max2; j++) {
	    fitPar[fitPar.size()-1] = (ptype) ? index++ : i;
	    chi2 += pow((fx[i][j]-func(xx[i],yy[j],params,fitPar))/df[i][j],2);
	  }

	return chi2;
      }


      double operator() (VecDoub Par)
      {  
	if (Par.size()!=int(par_limits.size())) 
	  ErrorMsg("Error in Chi2ClassFunc.h!");

	if (par_limits.size()>0)
	  for (int i=0; i<Par.size(); i++) 
	    if (Par[i] < par_limits[i][0] || Par[i] > par_limits[i][1]) return 1.e30;

	vector<double> fitPar(Par.size()+1);
	for (int i=0; i<Par.size(); i++) fitPar[i] = Par[i];

	double chi2 = 0;
	int index = 0;

	for (int i=min1; i<max1; i++) 
	  for (int j=min2; j<max2; j++) {
	    fitPar[fitPar.size()-1] = (ptype) ? index++ : i;
	    chi2 += pow((fx[i][j]-func(xx[i],yy[j],params,fitPar))/df[i][j],2);
	  }

	return chi2;
      }

    };


    // =====================================================================================

    // ----------------------------------------------------------------------------------
    // ----- log(chi2) class function, with only diagonal errors, for a 2D function -----
    // ----------------------------------------------------------------------------------

    class Chi2log_2D 
    {
  
    private: 

      // 2D function: y (x, y, fixed parameters, free parameters)  
      typedef double(*model) (double, double, shared_ptr<void> , vector<double>); 
  
      // ----- data -----
      vector<double> xx, yy;
      vector< vector<double> > fx, df;
  
      // ----- model, parameters, parameter limits -----
      model func;                          // model
      shared_ptr<void> params;                        // fixed parameters of the model
      vector< vector<double> > par_limits; // prior limits on the free parameters of the model
      bool ptype;                          // ptype = 0 -> fitPar[fitPar.size()-1] = i ; ptype = 1 -> fitPar[fitPar.size()-1] = index++

      // ----- limits of the fitting region -----
      int min1, min2, max1, max2;


    public:

      Chi2log_2D (vector<double> _xx, vector<double> _yy, vector< vector<double> > _fx, vector <vector<double> > _df, model _func, shared_ptr<void> _params, vector< vector<double> > _par_limits, bool _ptype=1, int _min1=0, int _min2=0, int _max1=0, int _max2=0)
	: xx(_xx), yy(_yy), fx(_fx), df(_df), func(_func), params(_params), par_limits(_par_limits), ptype(_ptype), min1(_min1), min2(_min2), max1((_max1>0) ? _max1 : _fx.size()), max2((_max2>0) ? _max2 : _fx[0].size()) {} 
  

      double operator() (double Par)
      {  
	if (par_limits.size()>0)
	  if (Par < par_limits[0][0] || Par > par_limits[0][1]) return 1.e30;
    
	vector<double> fitPar(2,Par); 

	double chi2 = 0, num = -1., den = -1.;
	int index = 0;

	for (int i=min1; i<max1; i++)
	  for (int j=min2; j<max2; j++) {
	    fitPar[fitPar.size()-1] = (ptype) ? index++ : i;
	    num = log(1.+fx[i][j])-log(1.+func(xx[i],yy[i],params,fitPar));
	    den = log(1.+fx[i][j]+df[i][j])-log(1.+fx[i][j]-df[i][j]);
	    if (den>0) chi2 += pow(num/den,2);
	  }

	return chi2;
      }


      double operator() (VecDoub Par)
      {  
	if (Par.size()!=int(par_limits.size())) 
	  ErrorMsg("Error in Chi2ClassFunc.h!");

	if (par_limits.size()>0)
	  for (int i=0; i<Par.size(); i++) 
	    if (Par[i] < par_limits[i][0] || Par[i] > par_limits[i][1]) return 1.e30;
    
	vector<double> fitPar(Par.size()+1);
	for (int i=0; i<Par.size(); i++) fitPar[i] = Par[i];

	double chi2 = 0, num = -1., den = -1.;
	int index = 0;

	for (int i=min1; i<max1; i++) 
	  for (int j=min2; j<max2; j++) {
	    fitPar[fitPar.size()-1] = (ptype) ? index++ : i;
	    num = log(1.+fx[i][j])-log(1.+func(xx[i],yy[i],params,fitPar));
	    den = log(1.+fx[i][j]+df[i][j])-log(1.+fx[i][j]-df[i][j]);
	    if (den>0) chi2 += pow(num/den,2);
	  }

	return chi2;
      }

    };

  }
}

#endif
