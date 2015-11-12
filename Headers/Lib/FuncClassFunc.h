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

/** @file Headers/Lib/FuncClassFunc.h
 *
 *  @brief Generic class functions used by Numerical methods
 *
 *  This file contains some generic class functions used by Numerical
 *  methods
 *
 *  @author Federico Marulli 
 *
 *  @author federico.marulli3@unbo.it
 */

#ifndef __FUNCCLASSFUNC__
#define __FUNCCLASSFUNC__ 


// =====================================================================================


namespace cosmobl {
  
  /** @brief The namespace of the class functions of the CosmoBolognaLib  
   *  
   *  The \e classfunc namespace contains all the class functions of
   *  the CosmoBolognaLib
   */
  namespace classfunc {

    class func_grid
    {
    private: 
      vector<double> xg, yg;
      string interpType;
      int Num;

    public:
      func_grid (vector<double> _xg, vector<double> _yg, string _interpType, int _Num)
	: xg(_xg), yg(_yg), interpType(_interpType), Num(_Num) {}  
  
      double operator() (double XX) 
      {
	return interpolated(XX, xg, yg, interpType, Num);
      }
    };


    // =====================================================================================


    class func_conv_gauss
    {
    private:
      vector<double> xx, fx;
      double xX;
      double sigma;
  
    public:
      func_conv_gauss (vector<double> _xx, vector<double> _fx, double _xX, double _sigma) 
	: xx(_xx), fx(_fx), xX(_xX), sigma(_sigma) {}
  
      double operator() (const double &XX) 
      {
	vector<double> par(2); par[0] = 0.; par[1] = sigma;
	vector<double> pp;

	return interpolated(XX, xx, fx, "Poly", 4)*gaussian(xX-XX, &pp, par);
      }
    };


    // =====================================================================================


    class func_xi
    {
    private:
      vector<double> lgkk, lgPk;
      double rr, aa;

    public:
      func_xi (vector<double> _lgkk, vector<double> _lgPK, double _rr, double _aa=0) 
	: lgkk(_lgkk), lgPk(_lgPK), rr(_rr), aa(_aa) {}

      double operator() (double kk) 
      { 
	double lgk = log10(kk);

	double lgPkK = interpolated(lgk, lgkk, lgPk, "Linear", -1);

	double Int = pow(10.,lgPkK)*sin(kk*rr)*kk/rr;

	return Int * exp(-kk*kk*aa*aa); // eq. 24 of Anderson et al. 2012
      }
    };


    // =====================================================================================


    class func_Pk
    {
    private:
      vector<double> lgrr, lgxi;
      double kk;

    public:
      func_Pk (vector<double> _lgrr, vector<double> _lgxi, double _kk) 
	: lgrr(_lgrr), lgxi(_lgxi), kk(_kk) {}

      double operator() (double rr) 
      { 
	double lgr = log10(rr);

	double lgxiR = interpolated(lgr, lgrr, lgxi, "Linear", -1);

	return pow(10.,lgxiR)*sin(rr*kk)*rr/kk;
      }
    };


    // =====================================================================================


    class func_wp
    {
    private: 
      vector<double> rr, xi;
      double rp;

    public:
      func_wp (vector<double> _rr, vector<double> _xi, double _rp) 
	: rr(_rr), xi(_xi), rp(_rp) {}

      double operator() (double rrr) 
      { 
	return interpolated(rrr, rr, xi, "Linear", -1)/sqrt(rrr*rrr-rp*rp)*rrr;
      }

    };


    // =====================================================================================


    class func_sigma_xi
    {
    private:
      vector<double> rr, xi;  
      double RR;

    public:
      func_sigma_xi (vector<double> _rr, vector<double> _xi, double _RR) 
	: rr(_rr), xi(_xi), RR(_RR) {}
  
      double operator() (const double &rad) 
      {
	double xiR = interpolated(rad, rr, xi, "Poly", 4); 

	return (3.-2.25*rad/RR+0.1875*pow(rad/RR,3))*rad*rad*xiR;
      }
    };


    // =====================================================================================


    class func_sigma_wp
    {
    private:
      vector<double> rp, wp;  
      double RR;

    public:
      func_sigma_wp (vector<double> _rp, vector<double> _wp, double _RR) 
	: rp(_rp), wp(_wp), RR(_RR) {}
  
      double operator() (const double &rad) 
      {
	double wpR = interpolated(rad, rp, wp, "Poly", 4); 

	double xx = rad/RR, gg;
	if (xx<=2) 
	  gg = 1./(2.*par::pi)*(3.*par::pi-9.*xx+pow(xx,3));
	else 
	  gg = 1./(2.*par::pi)*((-pow(xx,4)+11.*pow(xx,2)-28.)/sqrt(pow(xx,2)-4.)+pow(xx,3)-9.*xx+6.*asin(2./xx));
    
	return wpR*rad*gg;
      }
    };


    // =====================================================================================


    class func_xi_
    {
    private:
      vector<double> lgrr, lgxi;
      bool type;

    public:
      func_xi_ (vector<double> _lgrr, vector<double> _lgxi, bool _type) 
	: lgrr(_lgrr), lgxi(_lgxi), type(_type) {}  
  
      double operator() (double &rr) 
      {
	double lgr = log10(rr);

	double lgxiR = interpolated(lgr, lgrr, lgxi, "Linear", -1);
    
	return (type==0) ? pow(10.,lgxiR)*pow(rr,2) : pow(10.,lgxiR)*pow(rr,4);
      }
    };


    /* Alfonso Veropalumbo */

    // Basic class_function to find minima on an interpolated function ( VecDoub instead of double, to be used with Powell)

    //1D

    class func_grid_minimum_1D
    {
    private: 
      vector<double> xg, yg;
      string interpType;
      int Num;

    public:
      func_grid_minimum_1D (vector<double> _xg, vector<double> _yg, string _interpType, int _Num)
	: xg(_xg), yg(_yg), interpType(_interpType), Num(_Num) {}  
  
      double operator() (VecDoub XX) 
      {
	if (XX[0]>Max(xg) || XX[0] <Min(xg)) return 1.e30;

	return interpolated(XX[0], xg, yg, interpType, Num);
      }
    };


    //2D

    class func_grid_minimum_2D
    {
    private: 
      vector<double> x1g, x2g;
      vector<vector<double> > yg;
      string interpType;
      int Num;

    public:
      func_grid_minimum_2D (vector<double> _x1g, vector<double> _x2g, vector< vector<double> > _yg , string _interpType, int _Num)
	: x1g(_x1g), x2g(_x2g), yg(_yg), interpType(_interpType), Num(_Num) {}  
  
      double operator() (VecDoub XX) 
      {
	if (XX[0]>Max(x1g) || XX[0] <Min(x1g)) {return 1.e30;}
	if (XX[1]>Max(x2g) || XX[1] <Min(x2g)) return 1.e30;
	
	return interpolated_2D(XX[0], XX[1], x1g, x2g, yg, interpType, Num);
      }
    };

  }
}

#endif
