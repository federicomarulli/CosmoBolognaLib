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

    public:
      func_grid (vector<double> _xg, vector<double> _yg, string _interpType)
	: xg(_xg), yg(_yg), interpType(_interpType) {}  
  
      double operator() (double XX) 
      {
	return interpolated(XX, xg, yg, interpType);
      }
    };


    // =====================================================================================


    class func_grid_GSL
    {
       private: 
          gsl_spline *m_spline;
          const gsl_interp_type *m_type;
          gsl_interp_accel *m_acc;
          size_t m_size;

       public:
          func_grid_GSL (vector<double> _xg, vector<double> _yg, string _interpType)
          {
             m_size = _xg.size();
             m_acc = gsl_interp_accel_alloc();

             if (_interpType=="Linear") 
                m_type = gsl_interp_linear;

             else if (_interpType=="Poly") 
                m_type = gsl_interp_polynomial;

             else if (_interpType=="Spline") 
                m_type = gsl_interp_cspline;

             else if (_interpType=="Spline_periodic") 
                m_type = gsl_interp_cspline_periodic;

             else if (_interpType=="Akima") 
                m_type = gsl_interp_akima;

             else if (_interpType=="Akima_periodic") 
                m_type = gsl_interp_akima_periodic;

             else if (_interpType=="Steffen") 
                m_type = gsl_interp_steffen;

             else 
                ErrorMsg("Error in interpolated of Func.cpp: the value of string 'type' is not permitted!");

             m_spline = gsl_spline_alloc (m_type, m_size); 
             gsl_spline_init (m_spline, _xg.data(), _yg.data(), m_size);
          }  

          void free()
          {
             gsl_spline_free (m_spline);
             gsl_interp_accel_free (m_acc);
          }

          double operator() (double XX) 
          {
             return gsl_spline_eval (m_spline, XX, m_acc);
          }

	  vector<double> eval_func (vector<double> XX)
	  {
	    vector<double> YY;
	    
	    for(size_t i=0;i<XX.size();i++)
	      YY.push_back(gsl_spline_eval(m_spline,XX[i],m_acc));

	    return YY;
	  }

	  double integrate_qag(const double a, const double b, const double prec=1.e-2, const int limit_size=6, const int rule = 6)
	  {
	    function<double(double)> f = bind(&func_grid_GSL::operator(), this, std::placeholders::_1);
	   
	    glob::STR_generic_integrand pp;
	    pp.f = f;

	    gsl_function Func;
	    Func.function = &generic_integrand;
	    Func.params = &pp;

	    return GSL_integrate_qag(Func,a, b, prec, limit_size, rule);
	  }

	  double integrate_qaws(const double a, const double b, const double alpha =0, const double beta =0, const int mu=0, const int nu =0, const double prec=1.e-2, const int limit_size=6, const int rule = 6)
	  {
	    function<double(double)> f = bind(&func_grid_GSL::operator(), this, std::placeholders::_1);
	   
	    glob::STR_generic_integrand pp;
	    pp.f = f;

	    gsl_function Func;
	    Func.function = &generic_integrand;
	    Func.params = &pp;

	    return GSL_integrate_qaws(Func,a, b, alpha, beta, mu, nu, prec, limit_size);
	  }

	  double root(const double x_low, const double x_up, const double fx0=0, const double prec=1.e-2)
	  {
	    function<double(double)> f = bind(&func_grid_GSL::operator(), this, std::placeholders::_1);
	   
	    glob::STR_generic_roots pp;
	    pp.f = f;
	    pp.xx0 = fx0;

	    gsl_function Func;
	    Func.function = &generic_roots;
	    Func.params = &pp;

	    return GSL_brent(Func, x_low, x_up, prec);
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

	double lgPkK = interpolated(lgk, lgkk, lgPk, "Linear");

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

	double lgxiR = interpolated(lgr, lgrr, lgxi, "Linear");

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
	return interpolated(rrr, rr, xi, "Linear")/sqrt(rrr*rrr-rp*rp)*rrr;
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
	double xiR = interpolated(rad, rr, xi, "Poly"); 

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
	double wpR = interpolated(rad, rp, wp, "Poly"); 

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

	double lgxiR = interpolated(lgr, lgrr, lgxi, "Linear");
    
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

    public:
      func_grid_minimum_1D (vector<double> _xg, vector<double> _yg, string _interpType, int _Num)
	: xg(_xg), yg(_yg), interpType(_interpType) {}  
  
      double operator() (VecDoub XX) 
      {
	if (XX[0]>Max(xg) || XX[0] <Min(xg)) return 1.e30;

	return interpolated(XX[0], xg, yg, interpType);
      }
    };


    //2D

    class func_grid_minimum_2D
    {
    private: 
      vector<double> x1g, x2g;
      vector<vector<double> > yg;
      string interpType;
      
    public:
      func_grid_minimum_2D (vector<double> _x1g, vector<double> _x2g, vector< vector<double> > _yg , string _interpType, int _Num)
	: x1g(_x1g), x2g(_x2g), yg(_yg), interpType(_interpType) {}  
  
      double operator() (VecDoub XX) 
      {
	if (XX[0]>Max(x1g) || XX[0] <Min(x1g)) {return 1.e30;}
	if (XX[1]>Max(x2g) || XX[1] <Min(x2g)) return 1.e30;
	
	return interpolated_2D(XX[0], XX[1], x1g, x2g, yg, interpType);
      }
    };

  }

  namespace glob {
     struct STR_sigma2_integrand
     {
        int l1, l2;
        double density_inv,kk;
        vector<int> orders;
        vector<classfunc::func_grid_GSL> Pk_multipoles_interp;
     };

     struct STR_XiMultipoles_integrand
     {
        double r;
        int l;
        classfunc::func_grid_GSL *Pkl;
	double k_cut;
	double cut_pow;
     };

     struct STR_xi2D_smu_integrand
     {
       classfunc::func_grid_GSL *func;
       int order;
     };

     struct STR_covariance_XiMultipoles_integrand
     {
       classfunc::func_grid_GSL *s2, *jl1r1, *jl2r2;
     };
  }
}

#endif
