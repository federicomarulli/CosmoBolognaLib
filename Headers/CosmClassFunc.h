/********************************************************************
 *  Copyright (C) 2010 by Federico Marulli                          *
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
 *  @file Headers/CosmClassFunc.h
 *
 *  @brief Class functions used by Numerical methods inside the class
 *  Cosmology
 *
 *  This file contains the cosmological class functions used by
 *  Numerical methods inside the class Cosmology
 *
 *  @author Federico Marulli 
 *
 *  @author federico.marulli3@unibo.it
 */

#ifndef __COSMCLASSFUNC__
#define __COSMCLASSFUNC__ 


// =====================================================================================


namespace cbl {

  namespace classfunc {

    class E_inv {

    private: 
      double m_Omega_matter;  
      double m_Omega_baryon;         
      double m_Omega_neutrinos;      
      double m_massless_neutrinos;   
      int m_massive_neutrinos;       
      double m_Omega_DE;        
      double m_Omega_radiation;     
      double m_hh;                    
      double m_scalar_amp;
      double m_scalar_pivot;
      double m_n_spec;               
      double m_w0; 
      double m_wa;   
      double m_fNL;
      int m_type_NG;
      double m_tau;
      std::string m_model;
      bool m_unit;      

    public:
      E_inv (double Omega_matter, double Omega_baryon, double Omega_neutrinos, double massless_neutrinos, int massive_neutrinos, double Omega_DE, double Omega_radiation, double hh, double scalar_amp, double scalar_pivot, double n_spec, double w0, double wa, double fNL, int type_NG, double tau, std::string model, bool unit)
	: m_Omega_matter(Omega_matter), m_Omega_baryon(Omega_baryon), m_Omega_neutrinos(Omega_neutrinos), m_massless_neutrinos(massless_neutrinos), m_massive_neutrinos(massive_neutrinos), m_Omega_DE(Omega_DE), m_Omega_radiation(Omega_radiation), m_hh(hh), m_scalar_amp(scalar_amp), m_scalar_pivot(scalar_pivot), m_n_spec(n_spec), m_w0(w0), m_wa(wa), m_fNL(fNL), m_type_NG(type_NG), m_tau(tau), m_model(model), m_unit(unit) {}  
  
      double operator() (double redshift) 
      {
	cosmology::Cosmology cosm(m_Omega_matter, m_Omega_baryon, m_Omega_neutrinos, m_massless_neutrinos, m_massive_neutrinos, m_Omega_DE, m_Omega_radiation, m_hh, m_scalar_amp, m_scalar_pivot, m_n_spec, m_w0, m_wa, m_fNL, m_type_NG, m_tau, m_model, m_unit);
	
	return 1./cosm.EE(redshift);   
      }
    };


    // =====================================================================================


    class E_inv2
    {
    private: 
      double m_Omega_matter;  
      double m_Omega_baryon;         
      double m_Omega_neutrinos;      
      double m_massless_neutrinos;   
      int m_massive_neutrinos;       
      double m_Omega_DE;        
      double m_Omega_radiation;     
      double m_hh;                  
      double m_scalar_amp;
      double m_scalar_pivot;           
      double m_n_spec;               
      double m_w0;
      double m_wa; 
      double m_fNL;
      int m_type_NG;
      double m_tau;
      std::string m_model;
      bool m_unit;      

    public:
      E_inv2 (double Omega_matter, double Omega_baryon, double Omega_neutrinos, double massless_neutrinos, int massive_neutrinos, double Omega_DE, double Omega_radiation, double hh, double scalar_amp, double scalar_pivot, double n_spec, double w0, double wa, double fNL, int type_NG, double tau, std::string model, bool unit)
	: m_Omega_matter(Omega_matter), m_Omega_baryon(Omega_baryon), m_Omega_neutrinos(Omega_neutrinos), m_massless_neutrinos(massless_neutrinos), m_massive_neutrinos(massive_neutrinos), m_Omega_DE(Omega_DE), m_Omega_radiation(Omega_radiation), m_hh(hh), m_scalar_amp(scalar_amp), m_scalar_pivot(scalar_pivot), m_n_spec(n_spec), m_w0(w0), m_wa(wa), m_fNL(fNL), m_type_NG(type_NG), m_tau(tau), m_model(model), m_unit(unit) {}  
  
      double operator() (double redshift) 
      {
	cosmology::Cosmology cosm (m_Omega_matter, m_Omega_baryon, m_Omega_neutrinos, m_massless_neutrinos, m_massive_neutrinos, m_Omega_DE, m_Omega_radiation, m_hh, m_scalar_amp, m_scalar_pivot, m_n_spec, m_w0, m_wa, m_fNL, m_type_NG, m_tau, m_model, m_unit);
    
	return 1./(1.+redshift)/cosm.EE(redshift);   
      }
    };


    // =====================================================================================


    class E_inv3 {

    private:   
      double m_Omega_matter;  
      double m_Omega_baryon;         
      double m_Omega_neutrinos;      
      double m_massless_neutrinos;   
      int m_massive_neutrinos;       
      double m_Omega_DE;        
      double m_Omega_radiation;     
      double m_hh;                  
      double m_scalar_amp;
      double m_scalar_pivot;           
      double m_n_spec;               
      double m_w0;
      double m_wa; 
      double m_fNL;
      int m_type_NG;
      double m_tau;
      std::string m_model;

    public:
      E_inv3 (double Omega_matter, double Omega_baryon, double Omega_neutrinos, double massless_neutrinos, int massive_neutrinos, double Omega_DE, double Omega_radiation, double hh, double scalar_amp, double scalar_pivot, double n_spec, double w0, double wa, double fNL, int type_NG, double tau, std::string model)
	: m_Omega_matter(Omega_matter), m_Omega_baryon(Omega_baryon), m_Omega_neutrinos(Omega_neutrinos), m_massless_neutrinos(massless_neutrinos), m_massive_neutrinos(massive_neutrinos), m_Omega_DE(Omega_DE), m_Omega_radiation(Omega_radiation), m_hh(hh), m_scalar_amp(scalar_amp), m_scalar_pivot(scalar_pivot), m_n_spec(n_spec), m_w0(w0), m_wa(wa), m_fNL(fNL), m_type_NG(type_NG), m_tau(tau), m_model(model) {}  

      double operator() (double aa)
      {
	cosmology::Cosmology cosm (m_Omega_matter, m_Omega_baryon, m_Omega_neutrinos, m_massless_neutrinos, m_massive_neutrinos, m_Omega_DE, m_Omega_radiation, m_hh, m_scalar_amp, m_scalar_pivot, m_n_spec, m_w0, m_wa, m_fNL, m_type_NG, m_tau, m_model, 0);
	double redshift = 1./aa-1.;
	return (1.+redshift)/cosm.EE(redshift);  
      }
    };


    // =====================================================================================


    class func_fDE { // test in CPL parameterisation

    private:
      double m_w0;
      double m_wa;

    public:
      func_fDE (double w0, double wa) 
	: m_w0(w0), m_wa(wa) {}  
  
      double operator() (double redshift) 
      {
	cosmology::Cosmology cosm;
	cosm.set_w0(m_w0); cosm.set_wa(m_wa);
    
	return (1.+cosm.w_CPL(redshift))/(1.+redshift);    
      }
    };


    // =====================================================================================


    class func_z {

    private: 
      double m_Omega_matter;  
      double m_Omega_baryon;         
      double m_Omega_neutrinos;      
      double m_massless_neutrinos;   
      int m_massive_neutrinos;       
      double m_Omega_DE;        
      double m_Omega_radiation;     
      double m_hh;                  
      double m_scalar_amp;
      double m_scalar_pivot;           
      double m_n_spec;               
      double m_w0, m_wa;   
      double m_fNL;
      int m_type_NG;
      double m_tau;
      std::string m_model;
      bool m_unit;  
      double m_dd;

    public:
      func_z (double Omega_matter, double Omega_baryon, double Omega_neutrinos, double massless_neutrinos, int massive_neutrinos, double Omega_DE, double Omega_radiation, double hh, double scalar_amp, double scalar_pivot, double n_spec, double w0, double wa, double fNL, int type_NG, double tau, std::string model, bool unit, double dd)
	: m_Omega_matter(Omega_matter), m_Omega_baryon(Omega_baryon), m_Omega_neutrinos(Omega_neutrinos), m_massless_neutrinos(massless_neutrinos), m_massive_neutrinos(massive_neutrinos), m_Omega_DE(Omega_DE), m_Omega_radiation(Omega_radiation), m_hh(hh), m_scalar_amp(scalar_amp), m_scalar_pivot(scalar_pivot), m_n_spec(n_spec), m_w0(w0), m_wa(wa), m_fNL(fNL), m_type_NG(type_NG), m_tau(tau), m_model(model), m_unit(unit), m_dd(dd) {}  

      double operator() (double redshift) 
      {
	cosmology::Cosmology cosm(m_Omega_matter, m_Omega_baryon, m_Omega_neutrinos, m_massless_neutrinos, m_massive_neutrinos, m_Omega_DE, m_Omega_radiation, m_hh, m_scalar_amp, m_scalar_pivot, m_n_spec, m_w0, m_wa, m_fNL, m_type_NG, m_tau, m_model, m_unit);

	return cosm.D_C(redshift)-m_dd;
      }
    };


    // =====================================================================================


    class func_V {

    private: 
      double m_Omega_matter;  
      double m_Omega_baryon;         
      double m_Omega_neutrinos;      
      double m_massless_neutrinos;   
      int m_massive_neutrinos;       
      double m_Omega_DE;        
      double m_Omega_radiation;     
      double m_hh;                  
      double m_scalar_amp;
      double m_scalar_pivot;           
      double m_n_spec;               
      double m_w0, m_wa;
      double m_fNL;
      int m_type_NG;
      double m_tau;
      std::string m_model;
      bool m_unit;  
      double m_z_min;
      double m_Area;
      double m_VV;

    public:
      func_V (double Omega_matter, double Omega_baryon, double Omega_neutrinos, double massless_neutrinos, int massive_neutrinos, double Omega_DE, double Omega_radiation, double hh, double scalar_amp, double scalar_pivot, double n_spec, double w0, double wa, double fNL, int type_NG, double tau, std::string model, bool unit, double z_min, double Area, double VV)
	: m_Omega_matter(Omega_matter), m_Omega_baryon(Omega_baryon), m_Omega_neutrinos(Omega_neutrinos), m_massless_neutrinos(massless_neutrinos), m_massive_neutrinos(massive_neutrinos), m_Omega_DE(Omega_DE), m_Omega_radiation(Omega_radiation), m_hh(hh), m_scalar_amp(scalar_amp), m_scalar_pivot(scalar_pivot), m_n_spec(n_spec), m_w0(w0), m_wa(wa), m_fNL(fNL), m_type_NG(type_NG), m_tau(tau), m_model(model), m_unit(unit), m_z_min(z_min), m_Area(Area), m_VV(VV) {}  

      double operator() (double z_max) 
      {
	cosmology::Cosmology cosm (m_Omega_matter, m_Omega_baryon, m_Omega_neutrinos, m_massless_neutrinos, m_massive_neutrinos, m_Omega_DE, m_Omega_radiation, m_hh, m_scalar_amp, m_scalar_pivot, m_n_spec, m_w0, m_wa, m_fNL, m_type_NG, m_tau, m_model, m_unit);
    
	return cosm.Volume(m_z_min, z_max, m_Area)-m_VV;
      }
    };


    // =====================================================================================


    class func_zt {

    private: 
      double m_Omega_matter;  
      double m_Omega_baryon;         
      double m_Omega_neutrinos;      
      double m_massless_neutrinos;   
      int m_massive_neutrinos;       
      double m_Omega_DE;        
      double m_Omega_radiation;     
      double m_hh;                  
      double m_scalar_amp;
      double m_scalar_pivot;           
      double m_n_spec;               
      double m_w0;
      double m_wa;   
      double m_fNL;
      int m_type_NG;
      double m_tau;
      std::string m_model;
      bool m_unit;  
      double m_tt;

    public:
      func_zt (double Omega_matter, double Omega_baryon, double Omega_neutrinos, double massless_neutrinos, int massive_neutrinos, double Omega_DE, double Omega_radiation, double hh, double scalar_amp, double scalar_pivot, double n_spec, double w0, double wa, double fNL, int type_NG, double tau, std::string model, bool unit, double tt)
	: m_Omega_matter(Omega_matter), m_Omega_baryon(Omega_baryon), m_Omega_neutrinos(Omega_neutrinos), m_massless_neutrinos(massless_neutrinos), m_massive_neutrinos(massive_neutrinos), m_Omega_DE(Omega_DE), m_Omega_radiation(Omega_radiation), m_hh(hh), m_scalar_amp(scalar_amp), m_scalar_pivot(scalar_pivot), m_n_spec(n_spec), m_w0(w0), m_wa(wa), m_fNL(fNL), m_type_NG(type_NG), m_tau(tau), m_model(model), m_unit(unit), m_tt(tt) {}  

      double operator() (double redshift) 
      {
	cosmology::Cosmology cosm (m_Omega_matter, m_Omega_baryon, m_Omega_neutrinos, m_massless_neutrinos, m_massive_neutrinos, m_Omega_DE, m_Omega_radiation, m_hh, m_scalar_amp, m_scalar_pivot, m_n_spec, m_w0, m_wa, m_fNL, m_type_NG, m_tau, m_model, m_unit);

	return cosm.cosmic_time(redshift)-m_tt;
      }
    };
    

    // =====================================================================================


    class func_xistar {

    private:
      double m_Omega_matter;  
      double m_Omega_baryon;         
      double m_Omega_neutrinos;      
      double m_massless_neutrinos;   
      int m_massive_neutrinos;       
      double m_Omega_DE;        
      double m_Omega_radiation;     
      double m_hh;                  
      double m_scalar_amp;
      double m_scalar_pivot;           
      double m_n_spec;               
      double m_w0;
      double m_wa; 
      double m_fNL;
      int m_type_NG;
      double m_tau;
      std::string m_model;
      bool m_unit;   
      double m_rr;
      double m_redshift;
      bool m_store_output;
      std::string m_output_root;
      double m_kmax;
      double m_k_star;

    public:
      func_xistar (double Omega_matter, double Omega_baryon, double Omega_neutrinos, double massless_neutrinos, int massive_neutrinos, double Omega_DE, double Omega_radiation, double hh, double scalar_amp, double scalar_pivot, double n_spec, double w0, double wa, double fNL, int type_NG, double tau, std::string model, bool unit, double rr, double redshift, bool store_output, std::string output_root, double kmax, double k_star)
	: m_Omega_matter(Omega_matter), m_Omega_baryon(Omega_baryon), m_Omega_neutrinos(Omega_neutrinos), m_massless_neutrinos(massless_neutrinos), m_massive_neutrinos(massive_neutrinos), m_Omega_DE(Omega_DE), m_Omega_radiation(Omega_radiation), m_hh(hh), m_scalar_amp(scalar_amp), m_scalar_pivot(scalar_pivot), m_n_spec(n_spec), m_w0(w0), m_wa(wa), m_fNL(fNL), m_type_NG(type_NG), m_tau(tau), m_model(model), m_unit(unit), m_rr(rr), m_redshift(redshift), m_store_output(store_output), m_output_root(output_root), m_kmax(kmax), m_k_star(k_star) {}  
  
      double operator() (double kk) 
      { 
	cosmology::Cosmology cosm(m_Omega_matter, m_Omega_baryon, m_Omega_neutrinos, m_massless_neutrinos, m_massive_neutrinos, m_Omega_DE, m_Omega_radiation, m_hh, m_scalar_amp, m_scalar_pivot, m_n_spec, m_w0, m_wa, m_fNL, m_type_NG, m_tau, m_model, m_unit);
    
	std::string method_PkC = "CAMB";
	std::string method_PkEH = "EisensteinHu";
	bool NL = 0;
	int norm = 0;

	double Plin_BAO = cosm.Pk_matter({kk}, method_PkC, NL, m_redshift, m_store_output, m_output_root, norm, m_kmax)[0]-cosm.Pk_matter({kk}, method_PkEH, NL, m_redshift, m_store_output)[0];
	Plin_BAO *= exp(-kk*kk*0.5/(m_k_star*m_k_star));
    
	return Plin_BAO*sin(kk*m_rr)*kk/m_rr;  
      }
    };


    // =====================================================================================


    class func_V2 {

    private:    
      double m_Omega_matter;  
      double m_Omega_baryon;         
      double m_Omega_neutrinos;      
      double m_massless_neutrinos;   
      int m_massive_neutrinos;       
      double m_Omega_DE;        
      double m_Omega_radiation;     
      double m_hh;                   
      double m_scalar_amp;
      double m_scalar_pivot;           
      double m_n_spec;               
      double m_w0, m_wa;  
      double m_fNL;
      int m_type_NG;
      double m_tau;
      std::string m_model;
      bool m_unit;  
      std::string m_method_Pk;
      double m_rr;
      double m_redshift;
      bool m_store_output;

    public:
      func_V2 (double Omega_matter, double Omega_baryon, double Omega_neutrinos, double massless_neutrinos, int massive_neutrinos, double Omega_DE, double Omega_radiation, double hh, double scalar_amp, double scalar_pivot, double n_spec, double w0, double wa, double fNL, int type_NG, double tau, std::string model, bool unit, std::string method_Pk, double rr, double redshift, bool store_output)
	: m_Omega_matter(Omega_matter), m_Omega_baryon(Omega_baryon), m_Omega_neutrinos(Omega_neutrinos), m_massless_neutrinos(massless_neutrinos), m_massive_neutrinos(massive_neutrinos), m_Omega_DE(Omega_DE), m_Omega_radiation(Omega_radiation), m_hh(hh), m_scalar_amp(scalar_amp), m_scalar_pivot(scalar_pivot), m_n_spec(n_spec), m_w0(w0), m_wa(wa), m_fNL(fNL), m_type_NG(type_NG), m_tau(tau), m_model(model), m_unit(unit), m_method_Pk(method_Pk), m_rr(rr), m_redshift(redshift), m_store_output(store_output) {}  

      double operator() (double kk) 
      { 
	cosmology::Cosmology cosm (m_Omega_matter, m_Omega_baryon, m_Omega_neutrinos, m_massless_neutrinos, m_massive_neutrinos, m_Omega_DE, m_Omega_radiation, m_hh, m_scalar_amp, m_scalar_pivot, m_n_spec, m_w0, m_wa, m_fNL, m_type_NG, m_tau, m_model, m_unit);
	
	if (m_method_Pk=="EisensteinHu") 
	  return pow(cosm.linear_growth_rate(m_redshift, kk),2)*cosm.Pk_matter({kk}, m_method_Pk, false, m_redshift, m_store_output)[0]*pow(cbl::TopHat_WF(kk*m_rr),2);   

	else return cbl::ErrorCBL("", "func_V2", "CosmClassFunc.h"); 
      }
    };


    // =====================================================================================


    class func_V2_Table {

    private:   
      double m_Omega_matter;  
      double m_Omega_baryon;         
      double m_Omega_neutrinos;      
      double m_massless_neutrinos;   
      int m_massive_neutrinos;       
      double m_Omega_DE;        
      double m_Omega_radiation;     
      double m_hh;                  
      double m_scalar_amp;
      double m_scalar_pivot;           
      double m_n_spec;               
      double m_w0;
      double m_wa;  
      double m_fNL;
      int m_type_NG;
      double m_tau;
      std::string m_model;
      bool m_unit;  
      std::vector<double> m_lgkk, m_lgPk;
      double m_rr;
      double m_redshift;

    public:
      func_V2_Table (double Omega_matter, double Omega_baryon, double Omega_neutrinos, double massless_neutrinos, int massive_neutrinos, double Omega_DE, double Omega_radiation, double hh, double scalar_amp, double scalar_pivot, double n_spec, double w0, double wa, double fNL, int type_NG, double tau, std::string model, bool unit, std::vector<double> lgkk, std::vector<double> lgPk, double rr, double redshift)
	: m_Omega_matter(Omega_matter), m_Omega_baryon(Omega_baryon), m_Omega_neutrinos(Omega_neutrinos), m_massless_neutrinos(massless_neutrinos), m_massive_neutrinos(massive_neutrinos), m_Omega_DE(Omega_DE), m_Omega_radiation(Omega_radiation), m_hh(hh), m_scalar_amp(scalar_amp), m_scalar_pivot(scalar_pivot), m_n_spec(n_spec), m_w0(w0), m_wa(wa), m_fNL(fNL), m_type_NG(type_NG), m_tau(tau), m_model(model), m_unit(unit), m_lgkk(lgkk), m_lgPk(lgPk), m_rr(rr), m_redshift(redshift) {}  
  
      double operator() (double kk) 
      { 
	double fact = (m_unit) ? 1. : m_hh;
	double lgk = log10(kk/fact);

	cosmology::Cosmology cosm (m_Omega_matter, m_Omega_baryon, m_Omega_neutrinos, m_massless_neutrinos, m_massive_neutrinos, m_Omega_DE, m_Omega_radiation, m_hh, m_scalar_amp, m_scalar_pivot, m_n_spec, m_w0, m_wa, m_fNL, m_type_NG, m_tau, m_model, m_unit);

	double lgPkK = cbl::interpolated(lgk, m_lgkk, m_lgPk, "Linear");

	return pow(cosm.linear_growth_rate(m_redshift, kk),2)*pow(10.,lgPkK)/pow(fact, m_n_spec)*pow(cbl::TopHat_WF(kk*m_rr),2);  
      }
    };


    // =====================================================================================


    class func_sigma2 {

    private:     
      double m_Omega_matter;  
      double m_Omega_baryon;         
      double m_Omega_neutrinos;      
      double m_massless_neutrinos;   
      int m_massive_neutrinos;       
      double m_Omega_DE;        
      double m_Omega_radiation;     
      double m_hh;                  
      double m_scalar_amp;
      double m_scalar_pivot;           
      double m_n_spec;               
      double m_w0;
      double m_wa;   
      double m_fNL;
      int m_type_NG;
      double m_tau;
      std::string m_model;
      bool m_unit;  
      std::string m_method_Pk;
      double m_rr;
      double m_redshift;
      bool m_store_output;

    public:
      func_sigma2 (double Omega_matter, double Omega_baryon, double Omega_neutrinos, double massless_neutrinos, int massive_neutrinos, double Omega_DE, double Omega_radiation, double hh, double scalar_amp, double scalar_pivot, double n_spec, double w0, double wa, double fNL, int type_NG, double tau, std::string model, bool unit, std::string method_Pk, double rr, double redshift, bool store_output)
	: m_Omega_matter(Omega_matter), m_Omega_baryon(Omega_baryon), m_Omega_neutrinos(Omega_neutrinos), m_massless_neutrinos(massless_neutrinos), m_massive_neutrinos(massive_neutrinos), m_Omega_DE(Omega_DE), m_Omega_radiation(Omega_radiation), m_hh(hh), m_scalar_amp(scalar_amp), m_scalar_pivot(scalar_pivot), m_n_spec(n_spec), m_w0(w0), m_wa(wa), m_fNL(fNL), m_type_NG(type_NG), m_tau(tau), m_model(model), m_unit(unit), m_method_Pk(method_Pk), m_rr(rr), m_redshift(redshift), m_store_output(store_output) {}  

      double operator() (double kk) 
      { 
	cosmology::Cosmology cosm (m_Omega_matter, m_Omega_baryon, m_Omega_neutrinos, m_massless_neutrinos, m_massive_neutrinos, m_Omega_DE, m_Omega_radiation, m_hh, m_scalar_amp, m_scalar_pivot, m_n_spec, m_w0, m_wa, m_fNL, m_type_NG, m_tau, m_model, m_unit);
    
	if (m_method_Pk=="EisensteinHu") 
	  return pow(cosm.linear_growth_rate(m_redshift, kk),2)*cosm.Pk_matter({kk}, m_method_Pk, false, m_redshift, m_store_output)[0]*(1.-pow(cbl::TopHat_WF(kk*m_rr),2));
	
	else return cbl::ErrorCBL("", "func_sigma2", "CosmClassFunc.h"); 
      }
    };


    // =====================================================================================


    class func_sigma2_Table {

    private:  
      double m_Omega_matter;  
      double m_Omega_baryon;         
      double m_Omega_neutrinos;      
      double m_massless_neutrinos;   
      int m_massive_neutrinos;       
      double m_Omega_DE;        
      double m_Omega_radiation;     
      double m_hh;                  
      double m_scalar_amp;
      double m_scalar_pivot;           
      double m_n_spec;               
      double m_w0;
      double m_wa; 
      double m_fNL;
      int m_type_NG;
      double m_tau;
      std::string m_model;
      bool m_unit;  
      std::vector<double> m_lgkk, m_lgPk;
      double m_rr;
      double m_redshift;

    public:
      func_sigma2_Table (double Omega_matter, double Omega_baryon, double Omega_neutrinos, double massless_neutrinos, int massive_neutrinos, double Omega_DE, double Omega_radiation, double hh, double scalar_amp, double scalar_pivot, double n_spec, double w0, double wa, double fNL, int type_NG, double tau, std::string model, bool unit, std::vector<double> lgkk, std::vector<double> lgPk, double rr, double redshift)
	: m_Omega_matter(Omega_matter), m_Omega_baryon(Omega_baryon), m_Omega_neutrinos(Omega_neutrinos), m_massless_neutrinos(massless_neutrinos), m_massive_neutrinos(massive_neutrinos), m_Omega_DE(Omega_DE), m_Omega_radiation(Omega_radiation), m_hh(hh), m_scalar_amp(scalar_amp), m_scalar_pivot(scalar_pivot), m_n_spec(n_spec), m_w0(w0), m_wa(wa), m_fNL(fNL), m_type_NG(type_NG), m_tau(tau), m_model(model), m_unit(unit), m_lgkk(lgkk), m_lgPk(lgPk), m_rr(rr), m_redshift(redshift) {}  
  
      double operator() (double kk) 
      { 
	double fact = (m_unit) ? 1. : m_hh;
	double lgk = log10(kk/fact);

	cosmology::Cosmology cosm (m_Omega_matter, m_Omega_baryon, m_Omega_neutrinos, m_massless_neutrinos, m_massive_neutrinos, m_Omega_DE, m_Omega_radiation, m_hh, m_scalar_amp, m_scalar_pivot, m_n_spec, m_w0, m_wa, m_fNL, m_type_NG, m_tau, m_model, m_unit);
 
	double lgPkK = cbl::interpolated(lgk, m_lgkk, m_lgPk, "Linear");

	return pow(cosm.linear_growth_rate(m_redshift, kk),2)*pow(10.,lgPkK)/pow(fact, m_n_spec)*(1.-pow(cbl::TopHat_WF(kk*m_rr),2));  
      }
    };


    // =====================================================================================


    class func_MhaloMin {

    private:
      double m_Omega_matter;  
      double m_Omega_baryon;         
      double m_Omega_neutrinos;      
      double m_massless_neutrinos;   
      int m_massive_neutrinos;       
      double m_Omega_DE;        
      double m_Omega_radiation;     
      double m_hh;                  
      double m_scalar_amp;
      double m_scalar_pivot;           
      double m_n_spec;               
      double m_w0;
      double m_wa;   
      double m_fNL;
      int m_type_NG;
      double m_tau;
      std::string m_model;
      bool m_unit;   
      bool m_angle_rad;
      double m_n_halo, m_Mmax, m_z_min, m_z_max, m_Area;
      std::string m_model_MF, m_method_SS;
      bool m_store_output;
      std::string m_output_root;
      double m_Delta;
      std::string m_interpType; 
      double m_kmax; 
      std::string m_input_file;
      bool m_is_parameter_file;
      
    public:
      func_MhaloMin (double Omega_matter, double Omega_baryon, double Omega_neutrinos, double massless_neutrinos, int massive_neutrinos, double Omega_DE, double Omega_radiation, double hh, double scalar_amp, double scalar_pivot, double n_spec, double w0, double wa, double fNL, int type_NG, double tau, std::string model, bool unit, double n_halo, double Area, bool angle_rad, double z_min, double z_max, double Mmax, std::string model_MF, std::string method_SS, bool store_output, std::string output_root, const double Delta, std::string interpType, double kmax, std::string input_file, bool is_parameter_file)  
	: m_Omega_matter(Omega_matter), m_Omega_baryon(Omega_baryon), m_Omega_neutrinos(Omega_neutrinos), m_massless_neutrinos(massless_neutrinos), m_massive_neutrinos(massive_neutrinos), m_Omega_DE(Omega_DE), m_Omega_radiation(Omega_radiation), m_hh(hh), m_scalar_amp(scalar_amp), m_scalar_pivot(scalar_pivot), m_n_spec(n_spec), m_w0(w0), m_wa(wa), m_fNL(fNL), m_type_NG(type_NG), m_tau(tau), m_model(model), m_unit(unit), m_angle_rad(angle_rad), m_n_halo(n_halo), m_Mmax(Mmax), m_z_min(z_min), m_z_max(z_max), m_Area(Area), m_model_MF(model_MF), m_method_SS(method_SS), m_store_output(store_output), m_output_root(output_root), m_Delta(Delta), m_interpType(interpType), m_kmax(kmax), m_input_file(input_file), m_is_parameter_file(is_parameter_file) {}
  
      double operator() (double lgMmin) 
      {
	cosmology::Cosmology cosm (m_Omega_matter, m_Omega_baryon, m_Omega_neutrinos, m_massless_neutrinos, m_massive_neutrinos, m_Omega_DE, m_Omega_radiation, m_hh, m_scalar_amp, m_scalar_pivot, m_n_spec, m_w0, m_wa, m_fNL, m_type_NG, m_tau, m_model, m_unit);

	double Mmin = pow(10., lgMmin);

	double n_halo_expected = cosm.n_haloes(Mmin, m_Mmax, m_z_min, m_z_max, m_angle_rad, m_model_MF, m_method_SS, m_store_output, m_output_root, m_Delta, m_interpType, m_kmax, m_input_file, m_is_parameter_file)*m_Area;
   
	return m_n_halo-n_halo_expected;
      }
    };


    // =====================================================================================


    class func_kstar {

    private:
      double m_hh;
      bool m_unit;
      std::vector<double> m_lgkk, m_lgPk;

    public:
      func_kstar (double hh, bool unit, std::vector<double> lgkk, std::vector<double> lgPk)
	: m_hh(hh), m_unit(unit), m_lgkk(lgkk), m_lgPk(lgPk) {}
  
      double operator() (double kk) 
      { 
	double fact = (m_unit) ? 1. : m_hh;
	double lgk = log10(kk/fact);

	double lgPkK = cbl::interpolated(lgk, m_lgkk, m_lgPk, "Linear");

	return pow(10.,lgPkK); 
      }
    };

  }
}

#endif
