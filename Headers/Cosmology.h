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
 *  @file Headers/Cosmology.h
 *
 *  @brief The class Cosmology 
 *
 *  This file defines the interface of the class Cosmology, used for
 *  all kind of cosmological calculations
 *
 *  @author Federico Marulli 
 *
 *  @author federico.marulli3@unbo.it
 */

#ifndef __COSM__
#define __COSM__


#include "Likelihood.h"
#include "EisensteinHu.h"


// ===================================================================================================


namespace cbl {

  /**
   *  @brief The namespace of the functions and classes used for
   *  <B> cosmological calculations </B> 
   *  
   *  The \e cosmology namespace contains all the functions and classes
   *  used for cosmological calculations
   */
  namespace cosmology {

    /**
     *  @enum CosmologicalModel
     *
     *  @brief built-in cosmological models
     */
    enum class CosmologicalModel {
      
      /// Komatsu et al. 2009: Table 1, WMAP 5 Year Mean
      _WMAP5_,
      
      /// Komatsu et al. 2011: Table 1, WMAP Seven-year Mean
      _WMAP7_,

      /// Hinshaw et al. 2013: Table 3, WMAP-only Nine-year
      _WMAP9_,

      /// Planck collaboration 2013, paper XVI: Table 3, Planck+WP
      _Planck13_,
      
      /// Planck collaboration 2015, paper XIII: Table 4, TT,TE,EE+lowP+lensing
      _Planck15_,
      
      /// Planck collaboration 2015, paper XIII: Table 4, TT+lowP+lensing
      _Planck15_TT_,

      /// Planck collaboration 2018, Paper VI: Table 2, TT,TE,EE+lowE+lensing
      _Planck18_	
      
    };

    /**
     * @brief return a vector containing the
     * CosmologicalModel names
     * @return a vector containing the
     * CosmologicalModel names
     */
    inline std::vector<std::string> CosmologicalModelNames () { return {"WMAP5", "WMAP7", "WMAP9", "Planck13", "Planck15", "Planck15_TT", "Planck18"}; }

    /**
     *
     * @brief cast an enum of type CosmologicalModel
     * from its index
     * @param cosmologicalModelIndex the cosmologicalModel index
     * @return object of class CosmologicalModel
     */
    inline CosmologicalModel CosmologicalModelCast (const int cosmologicalModelIndex) { return castFromValue<CosmologicalModel>(cosmologicalModelIndex); }

    /**
     * @brief cast an enum of type CosmologicalModel
     * from its name
     * @param cosmologicalModelName the cosmologicalModel name
     * @return object of class CosmologicalModel
     */
    inline CosmologicalModel CosmologicalModelCast (const std::string cosmologicalModelName) { return castFromName<CosmologicalModel>(cosmologicalModelName, CosmologicalModelNames()); }

    /**
     * @brief cast an enum of type CosmologicalModel
     * from indeces
     * @param cosmologicalModelIndeces the cosmologicalModel indeces
     * @return object of class CosmologicalModel
     */
    inline std::vector<CosmologicalModel> CosmologicalModelCast (const std::vector<int> cosmologicalModelIndeces) { return castFromValues<CosmologicalModel>(cosmologicalModelIndeces); } 

    /**
     * @brief cast an enum of type CosmologicalModel
     * from thier names
     * @param cosmologicalModelNames the cosmologicalModel names
     * @return vector of CosmologicalModel enums
     */
    inline std::vector<CosmologicalModel> CosmologicalModelCast (const std::vector<std::string> cosmologicalModelNames) { return castFromNames<CosmologicalModel>(cosmologicalModelNames, CosmologicalModelNames()); }

    /**
     *  @enum CosmologicalParameter
     *  @brief the cosmological parameters
     */
    enum class CosmologicalParameter {

      /// \f$\Omega_M\f$: the density of baryons, cold dark matter and massive neutrinos (in units of the critical density) at z=0 in the LCDM case
      _Omega_matter_LCDM_,

      /// \f$\Omega_M\f$: the density of baryons, cold dark matter and massive neutrinos (in units of the critical density) at z=0
      _Omega_matter_,
  
      /// \f$\Omega_b\f$: the baryon density at z=0
      _Omega_baryon_,         
    
      /// \f$\Omega_b h^2\f$: the baryon density times \f$h^2\f$ at z=0
      _Omega_baryon_h2_,        

      /// \f$\Omega_\nu\f$: the density of massive neutrinos at z=0
      _Omega_neutrinos_,      

      /// \f$N_{eff}\f$: the effective number (for QED + non-instantaneous decoupling)
      _massless_neutrinos_,   

      /// the number of degenerate massive neutrino species 
      _massive_neutrinos_,

      /// the total neutrino mass
      _neutrino_mass_,

      /// \f$\Omega_{DE}\f$: the dark energy density at z=0
      _Omega_DE_,             

      /// \f$\Omega_{rad}\f$: the radiation density at z=0
      _Omega_radiation_,     
             
      /// \f$H_0\f$: the Hubble constant at z=0 [km/sec/Mpc] 
      _H0_,
    
      /// \f$h\f$: the Hubble constant at z=0 divided by 100 
      _hh_,           

      /// \f$A_s\f$: the initial scalar amplitude of the power spectrum
      _scalar_amp_,           

      /// \f$\ln(10^{10}A_s)\f$: the logarithm of 1e10 times the initial scalar amplitude of the power spectrum
      _ln_scalar_amp_,

      /// the scalar pivot k in \f$Mpc^{-1}\f$
      _scalar_pivot_,
      
      /// \f$n_{spec}\f$: the primordial spectral index
      _n_spec_,               

      /// \f$w_0\f$: the parameter of the dark energy equation of state (CPL parameterisation)
      _w0_,

      /// \f$w_a\f$: the parameter of the dark energy equation of state (CPL parameterisation)
      _wa_,             

      /// \f$f_{NL}\f$: the non-Gaussian amplitude
      _fNL_,                  

      /// \f$\sigma_8\f$: the power spectrum normalisation
      _sigma8_,
      
      /// \f$\tau\f$: Thomson scattering optical depth due to reionization
      _tau_,

      /// sound horizon
      _rs_
    };

    /**
     * @brief return a vector containing the
     * CosmologicalParameter names
     * @return a vector containing the
     * CosmologicalParameter names
     */
    inline std::vector<std::string> CosmologicalParameterNames () { return {"Omega_matter_LCDM", "Omega_matter", "Omega_baryon", "Omega_baryon_h2", "Omega_matter", "massless_neutrinos", "massive_neutrinos", "neutrino_mass", "Omega_DE", "Omega_radiation", "H0", "hh", "ln_scalar_amp", "scalar_amp", "scalar_pivot", "n_spec", "w0", "wa", "fNL", "sigma8", "tau", "rs"}; }
    
    /**
     *
     * @brief cast an enum of type CosmologicalParameter
     * from its index
     * @param cosmologicalParameterIndex the cosmologicalParameter index
     * @return object of class CosmologicalParameter
     */
    inline CosmologicalParameter CosmologicalParameterCast (const int cosmologicalParameterIndex) { return castFromValue<CosmologicalParameter>(cosmologicalParameterIndex); }

    /**
     * @brief cast an enum of type CosmologicalParameter
     * from its name
     * @param cosmologicalParameterName the cosmologicalParameter name
     * @return object of class CosmologicalParameter
     */
    inline CosmologicalParameter CosmologicalParameterCast (const std::string cosmologicalParameterName) { return castFromName<CosmologicalParameter>(cosmologicalParameterName, CosmologicalParameterNames()); }

    /**
     * @brief cast an enum of type CosmologicalParameter
     * from indeces
     * @param cosmologicalParameterIndeces the cosmologicalParameter indeces
     * @return object of class CosmologicalParameter
     */
    inline std::vector<CosmologicalParameter> CosmologicalParameterCast (const std::vector<int> cosmologicalParameterIndeces) { return castFromValues<CosmologicalParameter>(cosmologicalParameterIndeces); } 

    /**
     * @brief cast an enum of type CosmologicalParameter
     * from thier names
     * @param cosmologicalParameterNames the cosmologicalParameter names
     * @return vector of CosmologicalParameter enums
     */
    inline std::vector<CosmologicalParameter> CosmologicalParameterCast (const std::vector<std::string> cosmologicalParameterNames) { return castFromNames<CosmologicalParameter>(cosmologicalParameterNames, CosmologicalParameterNames()); }

    /**
     *  @brief name of the cosmological parameter
     *
     *  @param parameter the cosmological parameter
     *
     *  @return a std::string containing the name of the cosmological
     *  parameter provided in input
     */
    std::string CosmologicalParameter_name (const CosmologicalParameter parameter);

    /**
     *  @class Cosmology Cosmology.h "Headers/Cosmology.h"
     *
     *  @brief The class Cosmology
     *
     *  This class is used to handle objects of type <EM> cosmology
     *  </EM>. It can be used to estimate i) any kind of cosmic
     *  distances and volumes, ii) the halo mass function and bias, iii)
     *  the real-space and redshift-space dark matter power spectrum and
     *  correlation function (with both CAMB, CLASS, MPTbreeze and the
     *  analytic fitting forms by Eisenstein & Hu), iv) the accreted
     *  mass function, v) several BAO parameters, and vi) the halo mass
     *  function and bias in different non-Gaussian cosmological
     *  frameworks.
     */
    class Cosmology {

    private:
      
      /// \f$\Omega_M\f$: the density of baryons, cold dark matter and massive neutrinos (in units of the critical density) at z=0 in the LCDM case
      double m_Omega_matter;         
   
      /// \f$\Omega_b\f$: the baryon density at z=0
      double m_Omega_baryon;         

      /// \f$\Omega_\nu\f$: the density of massive neutrinos at z=0
      double m_Omega_neutrinos;      

      /// \f$N_{eff}\f$: the effective number (for QED + non-instantaneous decoupling)
      double m_massless_neutrinos;   

      /// the number of degenerate massive neutrino species 
      int m_massive_neutrinos;    

      /// \f$\Omega_{DE}\f$: the dark energy density at z=0
      double m_Omega_DE;             

      /// \f$\Omega_{rad}\f$: the radiation density at z=0 
      double m_Omega_radiation;     

      /// \f$\Omega_k\f$: the density of curvature energy
      double m_Omega_k;              

      /// \f$Omega_{CDM}\f$: the cold dark matter density at z=0
      double m_Omega_CDM;   
             
      /// \f$H_0\f$: the Hubble constant at z=0 [km/sec/Mpc] 
      double m_H0;

      /// \f$h\f$: the Hubble parameter, \f$H_0/100\f$
      double m_hh;

      /// \f$t_H\f$: the Hubble time
      double m_t_H;                  

      /// \f$D_H\f$: the Hubble distance
      double m_D_H;                  

      /// \f$sigma_8\f$: the power spectrum normalisation
      double m_sigma8;               

      /// \f$A_s\f$: the initial scalar amplitude of the power spectrum
      double m_scalar_amp;           

      /// the scalar pivot k in \f$Mpc^{-1}\f$
      double m_scalar_pivot;
      
      /// \f$n_{spec}\f$: the primordial spectral index
      double m_n_spec;               

      /// \f$w_0\f$: the parameter of the dark energy equation of state (CPL parameterisation)
      double m_w0;

      /// \f$w_a\f$: the parameter of the dark energy equation of state (CPL parameterisation)
      double m_wa;             
    
      /// \f$\rho_0\f$: the mean density of the Universe at z=0 [Msun*Mpc^-3]
      double m_RhoZero;             

      /// \f$f_{NL}\f$: the non-Gaussian amplitude
      double m_fNL;                  

      /// the non-Gaussian shape (type=1 local, type=2 equilateral, type=3 enfolded, type=4 orthogonal)
      int m_type_NG;    

      /// \f$\tau\f$: Thomson scattering optical depth due to reionization
      double m_tau;

      /// \f$r_s\f$ the sound horizon
      double m_rs = -1;
      
      /// the normalisation of the power spectrum for Eisenstein & Hu [http://background.uchicago.edu/~whu/transfer/transferpage.html]
      double m_Pk0_EH;

      /// the normalisation of the power spectrum for CAMB [http://camb.info/]
      double m_Pk0_CAMB;

      /// the normalisation of the power spectrum for MPTbreeze [http://arxiv.org/abs/1207.1465]
      double m_Pk0_MPTbreeze; 

      /// the normalisation of the power spectrum for CLASS [http://class-code.net/]
      double m_Pk0_CLASS; 
        
      /// the cosmologial model used to compute distances
      std::string m_model;                

      /// false \f$\rightarrow\f$ phyical units; true \f$\rightarrow\f$ cosmological units (i.e. without \e h)
      bool m_unit;

      
      /**
       *  @name Auxiliary functions of internal usage
       */
      ///@{

      /**
       *  @brief internal function to set default values
       *
       *  @return none
       */
      void set_default ();
      
      /**
       *  @brief function to compute the not-yet-normalised mass
       *  variances and their derivatives
       *
       *  this function computes the not-yet-normalised mass variances
       *  and their derivatives:
       *
       *  \f[ \sigma^2(R) = \frac{1}{2\pi^2}\int_0^\infty {\rm d}k\,
       *  k^2 P_{lin}(k, z) F^2(k, R)\f]
       *
       *  where \f$F(x)\f$ is a generic filter
       *
       *  @param method_Pk method used to compute the power spectrum;
       *  valid choices for method_Pk are: CAMB [http://camb.info/],
       *  classgal_v1 [http://class-code.net/], MPTbreeze-v1
       *  [http://arxiv.org/abs/1207.1465], EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       *
       *  @param redshift the redshift
       *
       *  @param output_root the output_root parameter of the
       *  parameter file used to compute the power spectrum; it can be
       *  any name
       * 
       *  @param interpType method to interpolate the power spectrum
       *
       *  @param k_max maximum wave vector module up to which the
       *  power spectrum is computed
       *
       *  @param input_file either the parameter file or the power
       *  spectrum file; if a parameter file is provided,
       *  i.e. input_file!=NULL and is_parameter_file=true, it will be
       *  used to compute the power spectrum; if a power spectrum file
       *  is provided, i.e. input_file!=NULL and
       *  is_parameter_file=false, then the provided power spectrum
       *  will be used directly; in both cases &sigma;<SUP>2</SUP>(M)
       *  is computed by integrating the computed/provided power
       *  spectrum ignoring the cosmological parameters of the object
       *
       *  @param is_parameter_file true \f$\rightarrow\f$ the
       *  input_file is a parameter file, used to compute the power
       *  spectrum with the method specified by method_Pk; false
       *  \f$\rightarrow\f$ the input_file is a file containing the
       *  power spectrum
       *
       *  @param filter the filter
       * 
       *  @param unit1 true \f$\rightarrow\f$ force cosmological units
       *
       *  @return the funciton to compute the not-yet-normalised 
       *  mass variances
       */
      double m_func_sigma (const std::string method_Pk, const double redshift, const std::string output_root="test", const std::string interpType="Linear", const double k_max=100., const std::string input_file=par::defaultString, const bool is_parameter_file=true, std::function<double(double)> filter={}, const bool unit1=false) const;

      /**
       *  @brief the not-yet-normalised mass variance,
       *  \f$\sigma^2(R)\f$
       *
       *  this function computes the not-yet-normalised variance of
       *  the linear density field:
       *
       *  \f[ \sigma^2(R) = \frac{1}{2\pi^2}\int_0^\infty {\rm d}k\,
       *  k^2 P_{lin}(k, z) W^2(k, R)\f]
       *
       *  where \f$W(x)=(3/x)^3(\sin x-x\cos x)\f$ and
       *  \f$R=(3M/4\pi\rho_m)^{1/3}\f$
       *
       *  @param radius the radius, \f$R\f$
       *
       *  @param method_Pk method used to compute the power spectrum;
       *  valid choices for method_Pk are: CAMB [http://camb.info/],
       *  classgal_v1 [http://class-code.net/], MPTbreeze-v1
       *  [http://arxiv.org/abs/1207.1465], EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       *
       *  @param redshift the redshift
       *
       *  @param output_root the output_root parameter of the
       *  parameter file used to compute the power spectrum; it can be
       *  any name
       * 
       *  @param interpType method to interpolate the power spectrum
       *
       *  @param k_max maximum wave vector module up to which the power
       *  spectrum is computed
       *
       *  @param input_file either the parameter file or the power
       *  spectrum file; if a parameter file is provided,
       *  i.e. input_file!=NULL and is_parameter_file=true, it will be
       *  used to compute the power spectrum; if a power spectrum file
       *  is provided, i.e. input_file!=NULL and
       *  is_parameter_file=false, then the provided power spectrum
       *  will be used directly; in both cases &sigma;<SUP>2</SUP>(M)
       *  is computed by integrating the computed/provided power
       *  spectrum ignoring the cosmological parameters of the object
       *
       *  @param is_parameter_file true \f$\rightarrow\f$ the
       *  input_file is a parameter file, used to compute the power
       *  spectrum with the method specified by method_Pk; false
       *  \f$\rightarrow\f$ the input_file is a file containing the
       *  power spectrum
       *
       *  @param unit1 true \f$\rightarrow\f$ force cosmological units
       *
       *  @return the not-yet-normalised \f$\sigma^2(R)\f$
       */
      double m_sigma2R_notNormalised (const double radius, const std::string method_Pk, const double redshift, const std::string output_root="test", const std::string interpType="Linear", const double k_max=100., const std::string input_file=par::defaultString, const bool is_parameter_file=true, const bool unit1=false) const;
      
      /**
       *  @brief the not-yet-normalised mass variance,
       *  \f$\sigma^2(M)\f$
       *
       *  this function computes the not-yet-normalised variance of
       *  the linear density field:
       *
       *  \f[ \sigma^2(M) = \frac{1}{2\pi^2}\int_0^\infty {\rm d}k\,
       *  k^2 P_{lin}(k, z) W^2(k, R)\f]
       *
       *  where \f$W(x)=(3/x)^3(\sin x-x\cos x)\f$ and
       *  \f$R=(3M/4\pi\rho_m)^{1/3}\f$
       *
       *  @param mass the mass
       *
       *  @param method_Pk method used to compute the power spectrum;
       *  valid choices for method_Pk are: CAMB [http://camb.info/],
       *  classgal_v1 [http://class-code.net/], MPTbreeze-v1
       *  [http://arxiv.org/abs/1207.1465], EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       *
       *  @param redshift the redshift
       *
       *  @param output_root the output_root parameter of the
       *  parameter file used to compute the power spectrum; it can be
       *  any name
       *
       *  @param interpType method to interpolate the power spectrum
       *
       *  @param k_max maximum wave vector module up to which the power
       *  spectrum is computed
       *
       *  @param input_file either the parameter file or the power
       *  spectrum file; if a parameter file is provided,
       *  i.e. input_file!=NULL and is_parameter_file=true, it will be
       *  used to compute the power spectrum; if a power spectrum file
       *  is provided, i.e. input_file!=NULL and
       *  is_parameter_file=false, then the provided power spectrum
       *  will be used directly; in both cases &sigma;<SUP>2</SUP>(M)
       *  is computed by integrating the computed/provided power
       *  spectrum ignoring the cosmological parameters of the object
       *
       *  @param is_parameter_file true \f$\rightarrow\f$ the input_file
       *  is a parameter file, used to compute the power spectrum with
       *  the method specified by method_Pk; false \f$\rightarrow\f$
       *  the input_file is a file containing the power spectrum
       *
       *  @param unit1 true \f$\rightarrow\f$ force cosmological units
       *
       *  @return the not-yet-normalised \f$\sigma^2(M)\f$
       */
      double m_sigma2M_notNormalised (const double mass, const std::string method_Pk, const double redshift, const std::string output_root="test", const std::string interpType="Linear", const double k_max=100., const std::string input_file=par::defaultString, const bool is_parameter_file=true, const bool unit1=false) const; 
      
      /**
       *  @brief auxiliary function to compute the mass function of
       *  dark matter haloes (filaments and sheets)
       *
       *  @author Alfonso Veropalumbo, Jacopo Neri (and Federico
       *  Marulli)
       *
       *  @author alfonso.veropalumbo@unibo.it, jacopo.neri6@gmail.com
       *  (and federico.marulli3@unibo.it)
       *
       *  @param Mass mass
       *
       *  @param mass_function_params function to a container of the
       *  mass function parameters
       *
       *  @return the mass function, d&Phi;/dM=dn(M)/dM
       */
      double m_mass_function (const double Mass, std::shared_ptr<void> mass_function_params);
      
      /**
       *  @brief auxiliary function to compute the mass function
       *
       *  @author Alfonso Veropalumbo, Jacopo Neri (and Federico
       *  Marulli)
       *
       *  @author alfonso.veropalumbo@unibo.it, jacopo.neri6@gmail.com
       *  (and federico.marulli3@unibo.it)
       *
       *  @param Mass mass
       *
       *  @param Sigma &sigma;(mass): the mass variance
       *
       *  @param Dln_Sigma dln&sigma;/dM: the derivative of the mass
       *  variance
       *
       *  @param redshift the redshift
       *
       *  @param model_MF author(s) who proposed the mass function;
       *  valid authors are: PS (Press & Schechter), ST (Sheth &
       *  Tormen), Jenkins (Jenkins et al. 2001), Warren (Warren et
       *  al. 2006), Reed (Reed et al. 2007), Pan (Pan 2007), ShenH
       *  (halo MF, Shen et al. 2006), ShenF (filament MF, Shen et
       *  al. 2006), ShenS (sheet MF, Shen et al. 2006), Tinker
       *  (Tinker et al. 2008), Crocce (Crocce et al. 2010),
       *  Angulo_FOF (FoF MF, Angulo et al. 2012), Angulo_Sub (SUBFIND
       *  MF, Angulo et al. 2012), Watson_FOF (FoF MF, Watson et
       *  al. 2012), Watson_SOH (Spherical Overdensity halo MF, Watson
       *  et al. 2012), Manera (Manera et al. 2010), Bhattacharya
       *  (Bhattacharya et al. 2011), Courtin (Courtin et al. 2010),
       *  Peacock (Peacock at al. 2007), Despali (Despali et al. 2016)
       *
       *  @param Delta \f$\Delta\f$: the overdensity, defined as the
       *  mean interior density relative to the background
       *
       *  @param default_delta true = using function
       *  cbl::cosmology::deltac; false = using delta_t*growth
       *  factor
       *  
       *  @param delta_t user defined density contrast at \f$z = 0\f$
       *
       *  @return the mass function, d&Phi;/dM=dn(M)/dM
       *
       *  @warning the input parameter \f$\Delta\f$ is the background
       *  overdensity, not the critical overdensity
       *  \f$\Delta_{crit}\f$; the function
       *  cbl::Cosmology::Delta_vir can be used to convert
       *  \f$\Delta_{crit}\f$ into \f$\Delta\f$
       *
       *  @warning the mass function by Manera et al. (2010) has been
       *  tested only for z=0 and z=0.5; the mass function by Despali
       *  et al. (2016) is currently implemented only for virial
       *  masses and at \f$z<1.25\f$
       */
      double m_MF_generator (const double Mass, const double Sigma, const double Dln_Sigma, const double redshift, const std::string model_MF, const double Delta=200., const bool default_delta=true, const double delta_t=1.686); 

      /**
       *  @brief auxiliary function to compute the halo bias
       *
       *  @param Sigma &sigma;(mass, z=0): the mass variance at z=0
       *
       *  @param redshift the redshift
       *
       *  @param author author(s) who proposed the bias; valid authors
       *  are: ST99 (Sheth & Tormen 1999), SMT01 (Sheth, Mo & Tormen
       *  2001), SMT01_WL04 (Sheth, Mo & Tormen 2001 with the correction
       *  of Warren 2004), Tinker (Tinker et al. 2010)
       *
       *  @param Delta \f$\Delta\f$: the overdensity, defined as the
       *  mean interior density relative to the background
       *
       *  @return the halo bias
       *
       *  @warning the input parameter \f$\Delta\f$ is the background
       *  overdensity, not the critical overdensity
       *  \f$\Delta_{crit}\f$; the function
       *  cbl::Cosmology::Delta_vir can be used to convert
       *  \f$\Delta_{crit}\f$ into \f$\Delta\f$
       */
      double m_bias_halo_generator (const double Sigma, const double redshift, const std::string author, const double Delta=200.) const;    

      /**
       *  @brief the incomplete elliptic integral 
       *
       *  this method computes the incomplete elliptic integral
       *  of the first kind F(&phi;|m), with shape parameter fixed to
       *  m=(2+3<SUP>0.5</SUP>)/4; it is used to compute the relation
       *  between comoving distance and redshift (see
       *  Numer. Math. (2010) 116:687–719, T. Fukushima "Fast
       *  computation of incomplete elliptic integral of first kind by
       *  half argument transformation")
       *
       *  @author Mauro Roncarelli
       *  @author mauro.roncarelli@unibo.it
       *
       *  @param phi &phi; scalar, argument
       *
       *  @return F(&phi;|m)
       *
       *  @warning this method works only for a fixed value of the shape
       *  parameter m=(2+3<SUP>0.5</SUP>)/4
       */
      double m_elf_dz (const double phi) const;

      /**
       *  @brief the inverse cosine amplitude of the Jacobian elliptic
       *  function
       *
       *  this method computes the inverse cosine amplitude of the
       *  Jacobian elliptic function cn<SUP>-1</SUP>(c|m), with shape
       *  parameter fixed to m=(2+3<SUP>0.5</SUP>)/4; namely
       *  cn<SUP>-1</SUP>(s|m) = F(arccos(c)|m); it is used to compute
       *  the relation between comoving distance and redshift (see
       *  Numer. Math. (2010) 116:687–719, T. Fukushima "Fast
       *  computation of incomplete elliptic integral of first kind by
       *  half argument transformation")
       *
       *  @author Mauro Roncarelli
       *  @author mauro.roncarelli@unibo.it
       *
       *  @param cc scalar, argument
       *
       *  @return cn<SUP>-1</SUP>(c|m)
       *
       *  @warning this method works only for a fixed value of the shape
       *  parameter m=(2+3<SUP>0.5</SUP>)/4; the argument c must be in
       *  the range 0<c<1
       */
      double m_acn_dz (const double cc) const;

      /**
       *  @brief the inverse sine amplitude of the Jacobian elliptic
       *  function
       *
       *  this method computes the inverse sine amplitude of the
       *  Jacobian elliptic function sn<SUP>-1</SUP>(s|m), with shape
       *  parameter fixed to m=(2+3<SUP>0.5</SUP>)/4, namely
       *  sn<SUP>-1</SUP>(s|m) = F(arcsin(s)|m); it is used to compute
       *  the relation between comoving distance and redshift (see
       *  Numer. Math. (2010) 116:687–719, T. Fukushima "Fast
       *  computation of incomplete elliptic integral of first kind by
       *  half argument transformation")
       *
       *  @author Mauro Roncarelli
       *  @author mauro.roncarelli@unibo.it
       *
       *  @param ss scalar, argument
       *
       *  @return sn<SUP>-1</SUP>(s|m)
       *
       *  @warning this method works only for a fixed value of the shape
       *  parameter m=(2+3<SUP>0.5</SUP>)/4; the argument s must be in
       *  the range 0<s<1
       */
      double m_asn_dz (const double ss) const;

      /**
       *  @brief the inverse truncated series necessary to compute
       *  sn<SUP>-1</SUP>(s|m) in ASN_DZ
       *
       *  this method computes the inverse truncated series necessary to
       *  compute sn<SUP>-1</SUP>(s|m) in ASN_DZ; the shape parameter is
       *  fixed to m=(2+3<SUP>0.5</SUP>)/4; it is used to compute the
       *  relation between comoving distance and redshift (see
       *  Numer. Math. (2010) 116:687–719, T. Fukushima "Fast
       *  computation of incomplete elliptic integral of first kind by
       *  half argument transformation"; the output of the function
       *  corresponds to the second factor of eq. (23))
       *
       *  @author Mauro Roncarelli
       *  @author mauro.roncarelli@unibo.it
       *
       *  @param yy scalar, argument
       *
       *  @return \f$\sum_{l=0}^L u_l(m)y^l\f$, with L fixed to 9
       *
       *  @warning this method works only for a fixed value of the shape
       *  parameter m=(2+3<SUP>0.5</SUP>)/4
       */
      double m_serf_dz (const double yy) const;

      ///@}

      
      // -----------------------------------------------------------------------


    public:
    
      /**
       *  @name Constructors/destructors
       */
      ///@{
    
      /**
       *  @brief constructor
       *
       *  @param Omega_matter \f$\Omega_M\f$: the density of baryons,
       *  cold dark matter and massive neutrinos (in units of the
       *  critical density) at z=0
       *
       *  @param Omega_baryon \f$\Omega_b\f$: the density of baryons
       *  at z=0
       *
       *  @param Omega_neutrinos \f$\Omega_\nu\f$: the density of
       *  massive neutrinos at z=0
       *
       *  @param massless_neutrinos the effective number (for QED +
       *  non-instantaneous decoupling)
       *
       *  @param massive_neutrinos the number of degenerate massive
       *  neutrino species
       *
       *  @param Omega_DE \f$\Omega_{DE}\f$: the density of dark
       *  energy at z=0
       *
       *  @param Omega_radiation \f$\Omega_{rad}\f$: the density of
       *  radiation at z=0
       *
       *  @param hh \e h: the Hubble parameter, \f$H_0/100\f$
       *
       *  @param scalar_amp \f$A_s\f$: the initial scalar amplitude of
       *  the power spectrum
       *
       *  @param scalar_pivot the scalar pivot k in \f$Mpc^{-1}\f$
       *
       *  @param n_spec \f$n_{spec}\f$: the primordial spectral index
       *
       *  @param w0 \f$w_0\f$: one of the two parameters of the dark
       *  energy equation of state (CPL parameterisation)
       *
       *  @param wa \f$w_a\f$: one of the two parameters of the
       *  dark energy equation of state (CPL parameterisation)
       *
       *  @param fNL \f$f_{NL}\f$: the non-Gaussian amplitude
       *
       *  @param type_NG the non-Gaussian shape (type=1 local, type=2
       *  equilateral, type=3 enfolded, type=4 orthogonal)
       *
       *  @param tau \f$\tau\f$: Thomson scattering optical depth due
       *  to reionization
       *
       *  @param model the cosmologial model used to compute distances
       *
       *  @param unit false \f$\rightarrow\f$ phyical units; true
       *  \f$\rightarrow\f$ cosmological units (i.e. without \e h)
       *
       *  @return none
       *
       *  @warning by default: \f$\Omega_k =
       *  1-\Omega_M-\Omega_{rad}-\Omega_{DM}\f$, \f$\Omega_{CDM} =
       *  \Omega_M-\Omega_b-\Omega_\nu\f$, \f$t_H = 1-H_0\f$, \f$D_H =
       *  c/t_H\f$, \f$\rho_0 = \rho_0(\Omega_M, \Omega_\nu)\f$
       */
      Cosmology (const double Omega_matter=0.27, const double Omega_baryon=0.046, const double Omega_neutrinos=0., const double massless_neutrinos=3.04, const int massive_neutrinos=0, const double Omega_DE=0.73, const double Omega_radiation=0., const double hh=0.7, const double scalar_amp=2.46e-9, const double scalar_pivot=0.05, const double n_spec=0.96, const double w0=-1., const double wa=0., const double fNL=0., const int type_NG=1, const double tau=0.09, const std::string model="LCDM", const bool unit=true);

      /**
       *  @brief constructor using built-in cosmological models
       *
       *  @param cosmoModel the built-in cosmological model
       *
       *  @param model the cosmologial model used to compute distances
       *
       *  @param unit 0 \f$\rightarrow\f$ phyical units; 1
       *  \f$\rightarrow\f$ cosmological units (i.e. without \e h)
       *
       *  @return object of class Cosmology; by default:
       *  \f$\Omega_k = 1-\Omega_M-\Omega_{rad}-\Omega_{DM}\f$,
       *  \f$\Omega_{CDM} = \Omega_M-\Omega_b-\Omega_\nu\f$,
       *  \f$t_H = 1-H_0\f$
       *  \f$D_H = c/t_H\f$
       *  \f$\rho_0 = \rho_0(\Omega_M, \Omega_nu)\f$
       */
      Cosmology (const CosmologicalModel cosmoModel, const std::string model="LCDM", const bool unit=true);    

      /**
       *  @brief default destructor
       *  @return none
       */
      ~Cosmology () = default;
    
      ///@}


      /**
       *  @name Functions to get the private members of the class
       */
      ///@{

      /**
       *  @brief get the private member specified by the enum CosmologicalParameter
       *  
       *  @param parameter the cosmological parameter
       *
       *  @return value of the cosmological parameter
       */
      double value (const CosmologicalParameter parameter) const;

      /**
       *  @brief get the private member Cosmology::m_Omega_matter
       *
       *  @return &Omega;<SUB>M</SUB>: the matter density, i.e. the
       *  density of baryons, cold dark matter, and massive neutrinos
       *  (in units of the critical density)
       */
      double Omega_matter () const { return m_Omega_matter; }; 

      /**
       *  @brief get the private member Cosmology::m_Omega_baryon
       *
       *  @return &Omega;<SUB>b</SUB>: the baryons density
       */
      double Omega_baryon () const { return m_Omega_baryon; };

      /**
       *  @brief get the private member Cosmology::m_Omega_neutrinos
       *
       *  @return &Omega;<SUB>&nu;</SUB>: the density of massive
       *  neutrinos
       */
      double Omega_neutrinos () const { return m_Omega_neutrinos; };

      /**
       *  @brief get the private member Cosmology::m_massless_neutrinos
       *
       *  @return N<SUB>eff</SUB>: the effective number (for QED +
       *  non-instantaneous decoupling)
       */
      double massless_neutrinos () const { return m_massless_neutrinos; };

      /**
       *  @brief get the private member Cosmology::m_massive_neutrinos
       *
       *  @return the number of degenerate massive neutrino species
       */
      int massive_neutrinos () const { return m_massive_neutrinos; };

      /**
       *  @brief get the private member Cosmology::m_Omega_DE
       *
       *  @return &Omega;<SUB>DE</SUB>: the dark energy density
       */
      double Omega_DE () const { return m_Omega_DE; };

      /**
       *  @brief get the private member Cosmology::m_Omega_radiation
       *
       *  @return &Omega;<SUB>rad</SUB>: the radiation density
       */
      double Omega_radiation () const { return m_Omega_radiation; };

      /**
       *  @brief get the private member Cosmology::m_Omega_k
       *
       *  @return &Omega;<SUB>k</SUB>: the density of curvature energy
       */
      double Omega_k () const { return m_Omega_k; };

      /**
       *  @brief get the private member Cosmology::m_Omega_CDM
       *
       *  @return &Omega;<SUB>CDM</SUB>: the cold dark matter density
       */
      double Omega_CDM () const { return m_Omega_CDM; };

      /**
       *  @brief get the private member Cosmology::m_H0
       *
       *  @return H<SUB>0</SUB>: the Hubble constant [km/sec/Mpc]
       */
      double H0 () const { return m_H0; };

      /**
       *  @brief get the private member Cosmology::m_hh
       *
       *  @return \e h: the Hubble parameter, H<SUB>0</SUB>/100
       */
      double hh () const { return m_hh; };

      /**
       *  @brief get the private member Cosmology::m_t_H
       *
       *  @return t<SUB>H</SUB>: the Hubble time
       */
      double t_H () const { return m_t_H; };

      /**
       *  @brief get the private member Cosmology::m_D_H
       *
       *  @return D<SUB>H</SUB>: the Hubble distance
       */
      double D_H () const { return m_D_H; };

      /**
       *  @brief get the private member Cosmology::m_sigma8
       *
       *  @return &sigma;<SUB>8</SUB>: the power spectrum normalisation
       */
      double sigma8 () const { return m_sigma8; };

      /**
       *  @brief get the private member Cosmology::m_scalar_amp
       *
       *  @return \f$A_s\f$: the initial scalar amplitude of the
       *  power spectrum
       */
      double scalar_amp () const { return m_scalar_amp; };

      /**
       *  @brief get the private member Cosmology::m_scalar_pivot
       *
       *  @return the scalar pivot k in \f$Mpc^{-1}\f$
       */
      double scalar_pivot () const { return m_scalar_pivot; };

      /**
       *  @brief get the private member Cosmology::m_n_spec
       *
       *  @return n<SUB>spec</SUB>: the primordial spectral index
       */
      double n_spec () const { return m_n_spec; };

      /**
       *  @brief get the private member Cosmology::m_w0
       *
       *  @return w<SUB>0</SUB>: one of the parameters of the dark
       *  energy equation of state (CPL parameterisation)
       */
      double w0 () const { return m_w0; };

      /**
       *  @brief get the private member Cosmology::m_wa
       *
       *  @return w<SUB>a</SUB>: one of the parameters of the dark
       *  energy equation of state (CPL parameterisation)
       */
      double wa () const { return m_wa; };

      /**
       *  @brief get the private member Cosmology::m_RhoZero
       *
       *  @return &rho;<SUB>0</SUB>: the mean density of the Universe at
       *  z=0 [Msun*Mpc^-3]
       */
      double RhoZero () const { return m_RhoZero; };

      /**
       *  @brief get the private member Cosmology::m_fNL
       *
       *  @return f<SUB>NL</SUB>: the non-Gaussian amplitude
       */
      double fNL () const { return m_fNL; };

      /**
       *  @brief get the private member Cosmology::m_type_NG
       *
       *  @return the non-Gaussian shape (type=1 local, type=2
       *  equilateral, type=3 enfolded, type=4 orthogonal)
       */
      int type_NG () const { return m_type_NG; };

      /**
       *  @brief get the private member Cosmology::m_tau
       *
       *  @return &tau; the Thomson scattering optical depth due to
       *  reionization
       */
      double tau () const { return m_tau; }; 

      /**
       *  @brief get the sound horizon at recombination
       *
       *  @return \f$r_s\f$; the sound horizon at recombination
       *  epoch
       */
      double rs () const { return ((m_rs==-1) ? rs_CAMB() : m_rs); };

      /**
       *  @brief get the private member Cosmology::m_Pk0_EH
       *
       *  @return the normalisation of the power spectrum for
       *  Eisenstein & Hu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       */
      double Pk0_EH () const { return m_Pk0_EH; };

      /**
       *  @brief get the private member Cosmology::m_Pk0_CAMB
       *
       *  @return the normalisation of the power spectrum for CAMB
       *  [http://camb.info/]
       */
      double Pk0_CAMB () const { return m_Pk0_CAMB; };

      /**
       *  @brief get the private member Cosmology::m_Pk0_MPTbreeze 
       *
       *  @return the normalisation of the power spectrum for MPTbreeze
       *  [http://arxiv.org/abs/1207.1465]
       */
      double Pk0_MPTbreeze () const { return m_Pk0_MPTbreeze; };

      /**
       *  @brief get the private member Cosmology::m_Pk0_CLASS
       *
       *  @return the normalisation of the power spectrum for CLASS [http://class-code.net/]
       */
      double Pk0_CLASS () const { return m_Pk0_CLASS; };

      /**
       *  @brief get the private member Cosmology::m_model
       *
       *  @return the cosmologial model used to compute distances
       */
      std::string model () const { return m_model; }

      /**
       *  @brief get the private member Cosmology::m_unit
       *
       *  @return unit: false \f$\rightarrow\f$ phyical units; true
       *  \f$\rightarrow\f$ cosmological units (i.e. without \e h)
       */
      bool unit () const { return m_unit; };

      /**
       *  @brief print the values of the private members on the screen
       *
       *  @return none
       */
      void print_parameters () const;
      
      ///@}


      /**
       *  @name Functions to set the private members of the class
       */
      ///@{

      /**
       *  @brief set the value of one cosmological paramter
       *  
       *  @param parameter cosmological parameter to set
       *  @param value the new value for the parameter 
       *
       *  @return none
       */
      void set_parameter (const CosmologicalParameter parameter, const double value);

      /**
       *  @brief set the value of some cosmological paramters 
       *  
       *  @param parameter vector containing the cosmological parameters
       *  to set
       *  @param value vector containing the new values for the
       *  parameters
       *
       *  @return none
       */
      void set_parameters (const std::vector<CosmologicalParameter> parameter, const std::vector<double> value);

      /**
       *  @brief set the value of &Omega;<SUB>M</SUB>, keeping
       *  &Omega;<SUB>DE</SUB>=1-&Omega;<SUB>M</SUB>-&Omega;<SUB>rad</SUB>-&Omega;<SUB>k</SUB>
       *
       *  @param Omega_matter &Omega;<SUB>M</SUB>: density of baryons,
       *  cold dark matter and massive neutrinos (in units of the
       *  critical density)
       *
       *  @return none
       */
      void set_Omega (const double Omega_matter=0.27) {
	m_Omega_matter = Omega_matter; 
	m_Omega_DE = 1.-m_Omega_matter-m_Omega_radiation-m_Omega_k;
	m_Omega_CDM = m_Omega_matter-m_Omega_baryon-m_Omega_neutrinos;
      };

      /**
       * @brief set the value of &Omega;<SUB>b</SUB>, keeping
       * &Omega;<SUB>CDM</SUB>=&Omega;<SUB>M</SUB>-&Omega;<SUB>b</SUB>
       *
       *  @param Omega_baryon &Omega;<SUB>b</SUB>: density of baryons,
       *  (in units of the critical density)
       *
       *  @return none
       */
      void set_OmegaB (const double Omega_baryon=0.046)
      {
	m_Omega_baryon = Omega_baryon; 
	m_Omega_CDM = m_Omega_matter-m_Omega_baryon-m_Omega_neutrinos;
      };

      /**
       *  @brief set the value of &Omega;<SUB>b</SUB>, keeping
       *  &Omega;<SUB>CDM</SUB>=&Omega;<SUB>M</SUB>-&Omega;<SUB>b</SUB>
       *
       *  @param Omega_baryonh2 &Omega;<SUB>b</SUB>h<SUP>2</SUP>:
       *  density of baryons, (in units of the critical density) times
       *  h<SUP>2</SUP>
       *
       *  @return none
       */
      void set_OmegaB_h2 (const double Omega_baryonh2=0.0222)
      {
	m_Omega_baryon = Omega_baryonh2*pow(m_hh,-2); 
	m_Omega_CDM = m_Omega_matter-m_Omega_baryon-m_Omega_neutrinos;
      };

      /**
       *  @brief set the value of &Omega;<SUB>M</SUB>
       *
       *  @param Omega_matter &Omega;<SUB>M</SUB>: density of baryons,
       *  cold dark matter and massive neutrinos (in units of the
       *  critical density)
       *
       *  @return none
       */
      void set_OmegaM (const double Omega_matter=0.27)
      {
	m_Omega_matter = Omega_matter; 
	m_Omega_k = 1.-m_Omega_matter-m_Omega_radiation-m_Omega_DE;
	m_Omega_CDM = m_Omega_matter-m_Omega_baryon-m_Omega_neutrinos;
      };

      /**
       *  @brief set the value of &Omega;<SUB>DE</SUB>
       *
       *  @param Omega_DE &Omega;<SUB>DE</SUB>: density of dark energy
       *
       *  @return none
       */
      void set_OmegaDE (const double Omega_DE=0.73)
      {
	m_Omega_DE = Omega_DE; 
	m_Omega_k = 1.-m_Omega_matter-m_Omega_radiation-m_Omega_DE;
	m_Omega_CDM = m_Omega_matter-m_Omega_baryon-m_Omega_neutrinos;
      };

      /**
       *  @brief set the value of &Omega;<SUB>&nu;</SUB>
       *  @param Omega_neutrinos &Omega;<SUB>&nu;</SUB>: density of
       *  massive neutrinos
       *  @param massless_neutrinos N<SUB>eff</SUB>: the effective
       *  number (for QED + non-instantaneous decoupling)
       *  @param massive_neutrinos the number of degenerate massive
       *  neutrino species
       *  @return none
       */
      void set_OmegaNu (const double Omega_neutrinos=0., const double massless_neutrinos=3.04, const int massive_neutrinos=0)
      {
	m_Omega_neutrinos = Omega_neutrinos; 
	m_massless_neutrinos = massless_neutrinos;
	m_massive_neutrinos = massive_neutrinos;
	if (m_Omega_neutrinos>0 && m_massive_neutrinos==0)
	  { m_massive_neutrinos = 1; m_massless_neutrinos = 2.04; }
	m_Omega_CDM = m_Omega_matter-m_Omega_baryon-m_Omega_neutrinos;
      };

      /**
       *  @brief set the private member Cosmology::m_Omega_radiation
       *
       *  @param Omega_radiation \f$\Omega_{rad}\f$: the radiation
       *  density
       *
       *  @return none
       */
      void set_Omega_radiation (const double Omega_radiation)
      { 
	m_Omega_radiation = Omega_radiation;
	m_Omega_k = 1.-m_Omega_matter-m_Omega_radiation-m_Omega_DE;
      };

      /**
       *  @brief set the value of h 
       *
       *  @param hh the Hubble constant H0/100
       *
       *  @param warn true \f$\rightarrow\f$ print a warning message
       *  if m_unit is true
       *
       *  @return none
       */
      void set_hh (const double hh=0.7, const bool warn=true)
      {
	if (m_unit && warn) WarningMsg("if unit=true then H0=100 (by internal definition)");
	m_hh = hh; 
	m_H0 = (m_unit) ? 100. : m_hh*100.;   
	m_t_H = 1./m_H0; 
	m_D_H = par::cc*m_t_H;
      };

      /**
       *  @brief set the value of H<SUB>0</SUB>
       *
       *  @param H0 H<SUB>0</SUB>: Hubble constant [km/sec/Mpc]
       *
       *  @param warn true \f$\rightarrow\f$ print a warning message if m_unit is
       *  true
       *
       *  @return none
       */
      void set_H0 (const double H0=70., const bool warn=true)
      {
	if (m_unit && warn) WarningMsg("if unit=true then H0=100 (by internal definition)");
	m_hh = H0/100.; 
	m_H0 = (m_unit) ? 100. : m_hh*100.;   
	m_t_H = 1./m_H0; 
	m_D_H = par::cc*m_t_H;
      };

      /**
       *  @brief set the value of &sigma;<SUB>8</SUB>
       *
       *  @param sigma8 &sigma;<SUB>8</SUB>: power spectrum normalisation
       *
       *  @return none
       */
      void set_sigma8 (const double sigma8=-1.) { m_sigma8 = sigma8; };

      /**
       *  @brief set the value of A<SUB>s</SUB>
       *
       *  @param scalar_amp \f$A_s\f$: initial scalar amplitude of
       *  the power spectrum
       *
       *  @return none
       */
      void set_scalar_amp (const double scalar_amp=2.46e-9) { m_scalar_amp = scalar_amp; }; 

      /**
       *  @brief set the value of the scalar pivot
       *
       *  @param scalar_pivot the scalar pivot k in \f$Mpc^{-1}\f$
       *
       *  @return none
       */
      void set_scalar_pivot (const double scalar_pivot=0.05) { m_scalar_pivot = scalar_pivot; }; 
      
      /**
       *  @brief set the value of n<SUB>spec</SUB>
       *
       *  @param n_spec n<SUB>spec</SUB>: the primordial spectral index
       *
       *  @return none
       */
      void set_n_spec (const double n_spec) { m_n_spec = n_spec; }; 

      /**
       *  @brief set the value of w<SUB>0</SUB>
       *
       *  @param w0 w<SUB>0</SUB>: parameter of the dark energy equation
       *  of state (CPL parameterisation)
       *
       *  @return none 
       */
      void set_w0 (const double w0=-1.) { m_w0 = w0; };  
    
      /**
       *  @brief set the value of w<SUB>a</SUB>
       *
       *  @param wa w<SUB>a</SUB>: parameter of the dark energy
       *  equation of state (CPL parameterisation)
       *
       *  @return none
       */
      void set_wa (const double wa=0.) { m_wa = wa; };  

      /**
       *  @brief set the value of &rho;<SUB>0</SUB>
       *
       *  @param RhoZero: the mean density of the Universe at z=0
       *  [Msun*Mpc^-3]
       *
       *  @return none
       */
      void set_RhoZero (const double RhoZero=7.5e10) { m_RhoZero = RhoZero; };  
    
      /**
       *  @brief set the value of f<SUB>NL</SUB>
       *
       *  @param fNL f<SUB>NL</SUB>: the non-Gaussian amplitude
       *
       *  @return none
       */
      void set_fNL (const double fNL=0.) { m_fNL = fNL; };  
    
      /**
       *  @brief set the value of the non-Gaussian shape
       *
       *  @param type_NG the non-Gaussian shape (type=1 local, type=2
       *  equilateral, type=3 enfolded, type=4 orthogonal)
       *
       *  @return none
       */
      void set_type_NG (const int type_NG=1) { m_type_NG = type_NG; };  

      /**
       *  @brief set the value of the &tau;
       *
       *  @param &tau; the Thomson scattering optical depth due to
       *  reionization
       *
       *  @return none
       */
      void set_tau (const double tau=0.09) { m_tau = tau; };  
    
      /**
       *  @brief set the value of the \f$r_s\f$;
       *
       *  @param rs the sound horizon
       *
       *  @return none
       */
      void set_rs (const double rs=-1) { m_rs = (rs==-1) ? rs_CAMB() : rs; };  

      /**
       *  @brief set the cosmologial model used to compute distances
       *
       *  @param model the cosmologial model used to compute distances
       *
       *  @return none
       */
      void set_model (const std::string model="LCDM") { m_model = model; };  

      /**
       *  @brief set the value of unit
       *
       *  @param unit false \f$\rightarrow\f$ phyical units; true \f$\rightarrow\f$
       *  cosmological units (i.e. without \e h)
       *
       *  @return none
       */
      void set_unit (const bool unit=true) { m_unit = unit; set_H0(100*m_hh, false); }

      ///@}

    
      /**
       *  @name Functions to estimate general cosmological parameters
       */
      ///@{

      /**
       *  @brief the matter density at a given redshift
       *
       *  @param redshift the redshift
       *
       *  @return &Omega;<SUB>M</SUB>
       */
      double OmegaM (const double redshift=0.) const; 

      /**
       *  @brief the dark energy density at a given redshift
       *
       *  @param redshift the redshift
       *
       *  @return &Omega;<SUB>DE</SUB>
       */
      double OmegaDE (const double redshift=0.) const;

      /**
       *  @brief the radiation density at a given redshift
       *
       *  @param redshift the redshift
       *
       *  @return &Omega;<SUB>rad</SUB>
       */
      double OmegaR (const double redshift=0.) const;

      /**
       *  @brief the radiation density, as a function of the redshift
       *  of radiation-matter equality
       *
       *  @param z_eq the redshift of radiation-matter equality
       *
       *  @return &Omega;<SUB>rad</SUB>
       */
      double OmegaR_zeq (const double z_eq=3395.) const;

      /**
       *  @brief the neutrino density at a given redshift
       *
       *  @param redshift the redshift
       *
       *  @return &Omega;<SUB>&nu;</SUB>
       */
      double OmegaNu (const double redshift=0.) const; 

      /**
       *  @brief the density of curvature energy at a given redshift
       *
       *  @param redshift the redshift
       *
       *  @return &Omega;<SUB>k</SUB>
       */
      double OmegaK (const double redshift=0.) const;

      /**
       *  @brief the cosmic density at a given redshift
       *
       *  @param redshift the redshift
       *
       *  @return \f$Omega\f$
       */
      double Omega (const double redshift=0.) const;

      /**
       *  @brief the density of massive neutrinos, given the neutrino
       *  mass
       *
       *  this function computes the density of massive neutrinos as
       *  follows:
       *
       *  \f[\Omega_nu = \frac{\sum m_\nu}{93.8h^2 eV}\f]
       *
       *  @param Mnu \f$sum m_\nu\f$ in eV
       *  
       *  @return \f$\Omega_nu\f$
       */
      double Omega_neutrinos (const double Mnu) const;

      /**
       *  @brief the total neutrino mass
       *
       *  this function computes the neutrino mass as follows:
       *
       *  \f[\sum m_\nu = \Omega_\nu\cdot94\, h^2 eV\f]
       *  
       *  @return \f$\sum m_\nu\f$
       */
      double neutrino_mass () const;

      /**
       *  @brief the critical cosmic density
       *
       *  this function computes the critical cosmic density at a
       *  given redshift:
       *
       *  \f[\rho_{crit}(z)=\frac{3H^2(z)}{8\pi G}\f]
       *
       *  @param redshift the redshift
       *
       *  @param unit1 true \f$\rightarrow\f$ force cosmological units
       *
       *  @return the critical cosmic density [Msun*Mpc^-3(*h^2)]
       */
      double rho_crit (const double redshift, const bool unit1=false) const;
      
      /**
       *  @brief the mean cosmic background density
       *
       *  \f[\rho_m(z) = \rho_{crit}(z)\Omega_M(z) =
       *  \frac{3H^2(z)}{8\pi G}\Omega_M(z)\f]
       *
       *  @param redshift the redshift
       *
       *  @param unit1 true \f$\rightarrow\f$ force cosmological units
       *
       *  @param nu true \f$\rightarrow\f$ compute \f$\rho_m(z) =
       *  \rho_{crit}(z)[\Omega_M(z)-\Omega_\nu(z)]\f$
       *
       *  @return &rho;<SUB>mean</SUB>: the mean cosmic background
       *  density [Msun*Mpc^-3(*h^2)]
       */
      double rho_m (const double redshift=0., const bool unit1=false, const bool nu=false) const;  
 
      /**
       *  @brief the critical overdensity
       *
       *  this function computes the critical overdensity,
       *  \f$\Delta_c(z)\equiv\Delta^{vir}_c(z)\f$, using approximated
       *  equations valid only for a flat Universe (see e.g. Coe
       *  2010), where \f$\rho_{vir} = \Delta^{vir}_c\rho_c =
       *  \Delta^{vir}_b\rho_m = \Delta^{vir}_b\Omega_M\rho_c\f$:
       *  
       *  - Bryan & Norman (1998)
       *
       *  \f[\Delta_c\simeq18\pi^2+60x-32x^2\, \mbox{for}\,
       *  \Omega_\Lambda=0\f]
       *
       *  \f[\Delta_c\simeq18\pi^2+82x-39x^2\, \mbox{for}\,
       *  \Omega_k=0\f]
       *  
       *  \f[x=\Omega_M(z)-1\f]
       *
       *  - Eke et al. (1998)
       * 
       *  \f[\Delta_c\simeq178\Omega_M(z)^{0.3}\, \mbox{for}\,
       *  \Omega_\Lambda=0\f]
       *
       *  \f[\Delta_c\simeq178\Omega_M(z)^{0.45}\, \mbox{for}\,
       *  \Omega_k=0\f]
       *
       *  - Nakamura & Suto (1998)
       *  
       *  \f[\Delta_c\simeq18\pi^2(1+0.4093x^{2.7152})\Omega_M(z)\f]
       *
       *  \f[x=(1-\Omega_{M,0})^{1/3}(1+z)^{-1}\f]
       *
       *  @param redshift the redshift 
       *
       *  @param author the author of the equation implemented;
       *  available options are: "BryanNorman", "Eke", "NakamuraSuto"
       *
       *  @return \f$\Delta_c\f$
       *
       *  @warning The implemented functions are approximated, and
       *  valid only for a restricted range of cosmological parameters
       *  (see the reported references)
       */
      double Delta_c (const double redshift, const std::string author="BryanNorman") const;

      /**
       *  @brief the virial overdensity given a critical overdensity
       *
       *  this function converts a given critical overdensity
       *  \f$\Delta_c(z)\equiv\Delta^{vir}_c(z)\f$ into the
       *  correspondent virial overdensity, at a given redshift
       *
       *  \f[\Delta_{vir}(z) \equiv \Delta^{vir}_b(z) =
       *  \frac{\Delta^{vir}_c(z)}{\Omega_M(z)}\f]
       *
       *  where \f$\rho_{vir} = \Delta^{vir}_c\rho_c =
       *  \Delta^{vir}_b\rho_m = \Delta^{vir}_b\Omega_M\rho_c\f$ (see
       *  e.g. Coe 2010)
       *
       *  @param Delta_c \f$\Delta_{crit}\f$: critical overdensity
       *
       *  @param redshift the redshift
       *
       *  @return \f$\Delta_{vir}\f$
       */
      double Delta_vir (const double Delta_c, const double redshift) const;
      
      /**
       *  @brief the virial overdensity
       *
       *  this function computes the virial overdensity:
       *
       *  \f[\Delta_{vir}(z) \equiv \Delta^{vir}_b(z) =
       *  \frac{\Delta^{vir}_c(z)}{\Omega_M(z)}\f]
       *
       *  where \f$\rho_{vir} = \Delta^{vir}_c\rho_c =
       *  \Delta^{vir}_b\rho_m = \Delta^{vir}_b\Omega_M\rho_c\f$ (see
       *  e.g. Coe 2010), and \f$\Delta^{vir}_c(z)\f$ is computed by
       *  cbl::cosmology::Cosmology::Delta_c
       *  
       *  @param redshift the redshift
       *
       *  @param author the author of the equation implemented;
       *  available options are: "BryanNorma", "Eke", "NakamuraSuto"
       *
       *  @return \f$\Delta_{vir}\f$
       */
      double Delta_vir (const double redshift, const std::string author="BryanNorman") const; 

      /**
       *  @brief the virial mass, given the virial radius and the
       *  redshift
       *
       *  this function computes the virial halo mass as follows:
       *
       *  \f[M_{vir}(z) = \frac{4}{3}\pi
       *  r_{vir}^3\Delta_c(z)\rho_{crit}(z) =
       *  \frac{r_{vir}^3\Delta_c(z)H^2(z)}{2G}\f]
       *
       *  where \f$\Delta_c(z)\f$ is computed by
       *  cbl::cosmology::Cosmology::Delta_c and
       *  \f$\rho_{crit}(z)\f$ is computed by
       *  cbl::cosmology::Cosmology::rho_crit
       *
       *  @param r_vir the virial radius
       * 
       *  @param redshift the redshift
       *
       *  @param author the author of the \f$\Delta_c(z)\f$
       *  equation (see cbl::cosmology::Cosmology::Delta_c)
       *
       *  @param unit1 true \f$\rightarrow\f$ force cosmological units
       *
       *  @return \f$M_{vir}\f$ 
       */
      double M_vir (const double r_vir, const double redshift, const std::string author="BryanNorman", const bool unit1=false) const;

      /**
       *  @brief the virial radius, given the virial mass and the
       *  redshift
       *
       *  this function computes the virial halo radius as follows:
       *
       *  \f[r_{vir}(z) = \left(\frac{3
       *  M_{vir}}{4\pi\Delta_c(z)\rho_{crit}(z)}\right)^{1/3}\f]
       *
       *  where \f$\Delta_c(z)\f$ is computed by
       *  cbl::cosmology::Cosmology::Delta_c and
       *  \f$\rho_{crit}(z)\f$ is computed by
       *  cbl::cosmology::Cosmology::rho_crit
       *
       *  @param M_vir the virial mass
       * 
       *  @param redshift the redshift
       *
       *  @param author the author of the \f$\Delta_c(z)\f$
       *  equation (see cbl::cosmology::Cosmology::Delta_c)
       *
       *  @param unit1 true \f$\rightarrow\f$ force cosmological units
       *
       *  @return \f$r_{vir}\f$ 
       */
      double r_vir (const double M_vir, const double redshift, const std::string author="BryanNorman", const bool unit1=false) const;
      
      /**
       *  @brief the DE equation of state in the CPL parameterisation,
       *  as a function of redshift
       *
       *  @param redshift the redshift
       *
       *  @return w: the DE equation of state in the CPL
       *  parameterisation
       */
      double w_CPL (const double redshift=0.) const;

      /**
       *  @brief auxiliary function used to compute the Hubble function
       *
       *  this function returns f<SUB>DE</SUB>(z), a parameter that
       *  multiplies &Omega;<SUB>DE</SUB> in the Hubble function (see
       *  e.g. Bassett & Hlozek 2010)
       *
       *  @param redshift the redshift
       *
       *  @return f<SUB>DE</SUB>
       */
      double f_DE (const double redshift=0.) const;
   
      /**
       *  @brief auxiliary function used to compute the Hubble function
       *  @param redshift the redshift
       *  @return E=H/H<SUB>0</SUB> 
       */
      double EE (const double redshift=0.) const;

      /**
       *  @brief inverse of the auxiliary function used to compute the
       *  Hubble function integrand of the comoving distance
       *  @param redshift the redshift
       *  @return 1./E=H<SUB>0</SUB>/H
       */
      double EE_inv (const double redshift=0.) const  
      { return 1./EE(redshift); }

      /**
       *  @brief inverse of the auxiliary function used to compute the
       *  Hubble function, integrand of the lookback time
       *  @param redshift the redshift
       *  @return 1/(1+z)/E=H<SUB>0</SUB>/H
       */
      double EE_inv2 (const double redshift=0.) const  
      { return 1./(1.+redshift)/EE(redshift); }

      /**
       *  @brief inverse of the auxiliary function used to compute the
       *  Hubble function integrand of the cosmic time
       *  @param aa 1/(1+z)
       *  @return (1+z)/E=H<SUB>0</SUB>/H
       */
      double EE_inv3 (const double aa) const
      {
	double redshift = 1./aa-1.;
	return (1.+redshift)/EE(redshift);
      }

      /**
       *  @brief the Hubble function
       *  @param redshift the redshift
       *  @return H(z)
       */
      double HH (const double redshift=0.) const;

      /**
       *  @brief the linear growth factor at a given redshift,
       *  \f$g(z)\f$
       *
       *  this function computes the following quantity (e.g. Eq.1 by
       *  Hamilton 2001):
       *
       *  \f[ g(z) \equiv D(z)(1+z) = \frac{5 \Omega_{M,0} H(a)}{2a}
       *  \int_0^a \frac{{\rm d}\,a'}{a^{'3} H^3(a')} \f]
       *
       *  where \f$a=1/(1+z)\f$
       *
       *  @param redshift the redshift
       *
       *  @return the linear growth factor
       *
       *  @warning the current implementation is valid in a
       *  Friedmann-Robertson-Walker Universe containing only matter
       *  and vacuum energy (e.g. Hamilton 2001 and references
       *  therein)
       */
      double gg (const double redshift=0.) const;   

      /**
       *  @brief the amplitude of the growing mode at a given
       *  redshift, \f$D(z)\f$
       *
       *  this function computes the following quantity (e.g. Eq.1
       *  by Hamilton 2001):
       *
       *  \f[ D(z) = \frac{g(z)}{1+z} = \frac{5 \Omega_{M,0} H(a)}{2}
       *  \int_0^a \frac{{\rm d}\,a'}{a^{'3} H^3(a')} \f]
       *
       *  where \f$a=1/(1+z)\f$ and \f$g(z)\f$ is computed by
       *  cbl::cosmology::Cosmology::gg
       *
       *  @param redshift the redshift 
       *
       *  @return the amplitude of the growing mode
       *
       *  @warning the current implementation is valid in a
       *  Friedmann-Robertson-Walker Universe containing only matter
       *  and vacuum energy (e.g. Hamilton 2001 and references
       *  therein)
       */
      double DD (const double redshift=0.) const;   

      /**
       *  @brief &sigma;<SUB>8</SUB> at a given redshift
       *
       *  @param redshift the redshift
       *
       *  @return &sigma;<SUB>8</SUB>
       */
      double sigma8 (const double redshift) const; 

      /**
       *  @brief lookback time at a given redshift
       *  @param redshift the redshift
       *  @return t<SUB>lookback</SUB> [Gyr]
       */
      double lookback_time (const double redshift=0.) const; 

      /**
       *  @brief cosmic time at a given redshift
       *  @param redshift the redshift
       *  @return t<SUB>cosmic</SUB> [Gyr]
       */
      double cosmic_time (const double redshift=0.) const; 

      /**
       *  @brief auxiliary function used to compute the deceleration
       *  parameter
       * 
       *  see e.g. de Araujo 2005
       *  
       *  @param redshift the redshift
       *  @return E<SUB>2</SUB>
       */
      double EE2 (const double redshift=0.) const;
    
      /**
       *  @brief the deceleration parameter at a given redshift
       *  @param redshift the redshift
       *  @return q
       */
      double qq (const double redshift=0.) const;

      /**
       *  @brief derivative of the Hubble function at a given redshift
       *  @param redshift the redshift
       *  @return dH/dz
       */
      double Hdot (const double redshift=0.) const;

      /**
       *  @brief redshift at wich occurs baryon photon decoupling
       *
       *  see Hu & Sugiyama (1996)
       *
       *  @return z<SUB>dec</SUB>
       */
      double z_decoupling () const;

      /**
       *  @brief redshift of drag epoch 
       *
       *  see Hu & Sugiyama (1996).
       *
       *  @return z<SUB>dec</SUB>
       */
      double z_drag () const;   

      /**
       *  @brief redshift at which the Universe begins to accelerate 
       *
       *  see e.g. de Araujo 2005
       *
       *  @return z<SUB>acc</SUB>
       */
      double z_acc () const; 

      /**
       *  @brief redshift of matter-dark energy equality
       *
       *  see e.g. de Araujo 2005
       *
       *  @return z<SUB>acc</SUB>
       */
      double z_eq () const;

      /**
       *  @brief the sound speed
       *
       *  @param redshift the redshift
       *
       *  @param T_CMB the temperature of the Cosmic Microwave
       *  Background
       *
       *  @return the sound speed
       */
      double sound_speed (const double redshift, const double T_CMB=2.7255) const;
    
      /**
       *  @brief the sound horizon integrand
       *
       *  @param redshift the redshift
       *
       *  @param T_CMB the temperature of the Cosmic Microwave
       *  Background
       *
       *  @return the sound horizon integrand
       */
      double rs_integrand (const double redshift, const double T_CMB=2.7255) const;

      /**
       *  @brief the sound horizon
       *
       *  @param redshift the redshift
       *
       *  @param T_CMB the temperature of the Cosmic Microwave
       *  Background
       *
       *  @return the sound horizon
       */
      double rs (const double redshift, const double T_CMB=2.7255) const;

      /**
       *  @brief maximum absolute magnitude to have a volume-limited
       *  catalogue
       *
       *  @param z_max maximum redshift
       *  @param mag_lim magnitude limit
       *  @return z<SUB>acc</SUB>
       */
      double Mag_Volume_limited (const double z_max=1., const double mag_lim=-20.) const;

      /**
       *  @brief bolometric luminosity
       *  @param redshift the redshift
       *  @param flux flux
       *  @return L<SUB>bol</SUB>
       */
      double Lum_bol (const double redshift=0., const double flux=1.) const;  

      /**
       *  @brief redshift at a given comoving distance
       *  
       *  this method provides the redshift for a given comoving
       *  distance
       *
       *  @param d_c line-of-sight comoving distance
       *
       *  @param z1_guess minimum prior on the redshift
       *
       *  @param z2_guess maximum prior on the redshift
       *
       *  @param prec precision of the computation ( prec =
       *  min(prec,1.e-5) )
       *
       *  @return redshift
       */
      double Redshift (const double d_c=1., const double z1_guess=0., const double z2_guess=10., const double prec=0.0001) const; 

      /**
       *  @brief redshift at a given comoving distance
       *
       *  this method provides the redshift for a given comoving
       *  distance; the iteration process exploits the properties of the
       *  f=d(z) function (i.e. f'=1/E(z), f"<0)
       *     
       *  @author Mauro Roncarelli
       *  @author mauro.roncarelli@unibo.it
       *
       *  @param d_c line-of-sight comoving distance
       *
       *  @param z1_guess minimum redshift of the region explored to
       *  search the redshift
       *
       *  @param z2_guess maximum redshift used to search the redshift
       *
       *  @param go_fast 0 \f$\rightarrow\f$ the method uses the
       *  function Cosmology::D_C; 1 \f$\rightarrow\f$ the method uses
       *  the function Cosmology::D_C_LCDM (much faster than D_C)
       *
       *  @param prec precision of the computation; ( prec =
       *  min(prec,1.e-5) )
       *
       *  @return redshift
       *
       *  @warning if go_fast=1, the method works correctly only for a
       *  flat &Lambda;CDM universe at low enough redshifts where
       *  &Omega;<SUB>r</SUB> is negligible;
       *
       *  @warning wrong redshift interval limits do not lead to an
       *  error, but just slow down the computation
       */
      double Redshift_LCDM (const double d_c=1., const double z1_guess=0., const double z2_guess=10., const bool go_fast=1, const double prec=0.0001) const; 

      /**
       *  @brief redshift at a given wf
       *
       *  this routine estimates the redshift from wf, given the
       *  parent halo mass at z=z', z', and its assembled fraction f
       *
       *  @author Carlo Giocoli
       *  @author cgiocoli@gmail.com
       *  @param mm mass
       *  @param redshift the redshift
       *  @param ff assembled fraction
       *
       *  @param method_SS method used to compute the power spectrum
       *  and &sigma;(mass) const; valid method_SS are: CAMB
       *  [http://camb.info/], classgal_v1 [http://class-code.net/],
       *  EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       *
       *  @param wwf rescaled variable w as in Lacey and Coles 1993
       *
       *  @param output_root output_root of the parameter file used to
       *  compute the power spectrum and &sigma;(mass); it can be any
       *  name
       *
       *  @return redshift
       */
      double Redshift (const double mm, const double redshift, const double ff, const std::string method_SS, const double wwf, const std::string output_root="test") const; 

      /**
       *  @brief redshift at a given cosmic time
       *  @param time cosmic time
       *  @param z1_guess minimum redshift used to search the redshift
       *  @param z2_guess maximum redshift used to search the redshift
       *  @return redshift
       */
      double Redshift_time (const double time, const double z1_guess, const double z2_guess) const;

      /**
       *  @brief spherical collapse density threshold at a given
       *  redshift
       *
       *  this function computes the spherical collapse density
       *  threshold, \f$\delta_c\f$, by using the approximated
       *  equation (C.28) provided by Nakamura & Suto (1997)
       *
       *  \f[\delta_c(z) \simeq \frac{3}{20}(12\pi)^{2/3}
       *  \{1.+0.012299\log[\Omega_M(z)]\}\f]
       *
       *  @param redshift the redshift
       *  @return \f$\delta_c\f$
       */
      double deltac (const double redshift) const;  
    
      /**
       *  @brief Linear (under)density contrast
       *
       *  @author Tommaso Ronconi
       *  @author tommaso.ronconi@studio.unibo.it
       *
       *  @param bias the bias of the sample
       *
       *  @param rho_vm the non linear density contrast:
       *  \f$\rho_v/\rho_m\f$ (default value set to \f$0.205\f$)
       *  
       *  @return The linear density contrast used as second barrier
       *  in the excursion set formalism for voids, based on the fit
       *  by Bernardeu (1994): \f$\delta_v^L \equiv \frac{\rho_v -
       *  \rho_m}{\rho_m} \approx C [1 - (\rho_v/\rho_m)^{- 1/C}]\f$
       *  where \f$\rho_v =\ \f$ average void density, \f$\rho_m =\
       *  \f$ average density of the surrounding Universe and \f$C =
       *  1.594\f$, a costant.
       */
      double deltav_L (const double bias=1., const double rho_vm=0.205) const;
    
      /**
       *  @brief Non-Linear (under)density contrast
       *
       *  @author Tommaso Ronconi
       *  @author tommaso.ronconi@studio.unibo.it
       *
       *  @param deltav the linear density contrast: \f$\delta_v\f$
       *  (default value set to \f$-2.71\f$)
       *  
       *  @return The non linear density contrast to be used in
       *  samples of tracers, independently of the bias value, based
       *  on the fit by Bernardeu (1994): \f$ \delta_v^{NL} \equiv
       *  \frac{\rho_v}{\rho_m} - 1 \approx \bigl(1 -
       *  C^{-1}\delta_v^L\bigr)^{-C} - 1\f$ where \f$\rho_v =\ \f$
       *  average void density, \f$\rho_m =\ \f$ average density of
       *  the surrounding Universe and \f$C = 1.594\f$, a costant.
       */
      double deltav_NL (const double deltav=-2.71) const;

      /**
       *  @brief expansion factor
       *
       *  @author Tommaso Ronconi
       *  @author tommaso.ronconi@studio.unibo.it
       *
       *  @param deltav the linear density contrast: \f$\delta_v^L\f$
       *  (default value set to \f$-2.71\f$)
       *  
       *  @return the expansion factor: \f$\frac{r}{r_L} = \bigl((1 -
       *  C^{-1}\delta_v^L\bigr)^{C/3}\f$ where \f$C = 1.594\f$
       */
      double r_rL (const double deltav = -2.71) const;

      ///@}
    

      /**
       *  @name Functions to estimate cosmological distances and volumes
       */
      ///@{

      /**
       *  @brief the comoving line-of-sight distance at a given redshift
       *
       *  for demonstration, see \ref distances.cpp
       *
       *  @param redshift the redshift
       *  @return D<SUB>C</SUB>
       */
      double D_C (const double redshift) const;  

      /**
       *  @brief the comoving line-of-sight distance at a given redshift
       *
       *  this method provides the comoving distance for a given
       *  redshift using elliptic integrals (see Numer. Math. (2010)
       *  116:687–719, T. Fukushima "Fast computation of incomplete
       *  elliptic integral of first kind by half argument
       *  transformation")
       *
       *  d(z) is computed using the following formula, taken from
       *  A. Meszaros & J. Ripa 2013, arXiv:1306.4736,
       *
       *  d(z) = C*[F(&phi;<SUB>0</SUB>,m)-F(&phi;<SUB>1</SUB>,m)]
       *
       *  where
       *
       *  C = 3<SUP>-1/4</SUP> /
       *  (&Omega;<Sub>M</SUB><SUP>1/3</SUP>*&Omega;<SUB>&Lambda;</SUB><SUP>1/6</SUP>)
       *
       *  &phi;<SUB>0</SUB> = arccos(
       *  [1+(1-3<SUP>0.5</SUP>)*(&Omega;<SUB>&Lambda;</SUB>/&Omega;<Sub>M</SUB>)<SUP>1/3</SUP>]
       *  /
       *  [1+(1+3<SUP>0.5</SUP>)*(&Omega;<SUB>&Lambda;</SUB>/&Omega;<Sub>M</SUB>)<SUP>1/3</SUP>]
       *  )
       *
       *  &phi;<SUB>1</SUB> = arccos(
       *  [1+(1-3<SUP>0.5</SUP>)*(&Omega;<SUB>&Lambda;</SUB>/&Omega;<Sub>M</SUB>)<SUP>1/3</SUP>*(1+z)]
       *  /
       *  [1+(1+3<SUP>0.5</SUP>)*(&Omega;<SUB>&Lambda;</SUB>/&Omega;<Sub>M</SUB>)<SUP>1/3</SUP>*(1+z)]
       *  )
       *
       *  F(&phi;, m) is the incomplete elliptic integral of the first
       *  kind, with argument &phi;
       *
       *  m = (2+3<SUP>0.5</SUP>)/4 is the elliptic integral shape
       *  parameter
       *  
       *  for demonstration, see \ref distances.cpp
       *
       *  @author Mauro Roncarelli
       *  @author mauro.roncarelli@unibo.it
       *
       *  @param redshift the redshift
       *  @return D<SUB>C</SUB>
       *
       *  @warning this method works only for a flat &Lambda;CDM
       *  universe at low enough redshifts where &Omega;<SUB>r</SUB> is
       *  negligible; it does not work for non-standard dark energy or
       *  non-flat models
       */
      double D_C_LCDM (const double redshift) const;  


      // table of redshift -- comoving line-of-sight distance

      /**
       *  @brief create a table of [redshift, comoving line-of-sight
       *  distance]
       *
       *  this function is used to create a table of [redshift, comoving
       *  line-of-sight distance], useful to speed up the analysis
       *
       *  @param [in] file_table name of the file where the table is
       *  stored
       *  @param [in] z_min minimum redshift of the table
       *  @param [in] z_max maximum redshift of the table
       *  @param [in] step redshift step
       *  @param [out] Redshift vector of redshifts
       *  @param [out] dc vector of comoving line-of-sight distances
       *  @return none
       */
      void D_C_table (const std::string file_table, const double z_min, const double z_max, const int step, std::vector<double> &Redshift, std::vector<double> &dc) const;

      /**
       *  @brief the comoving transverse distance at a given redshift
       *  @param redshift the redshift
       *  @return D<SUB>M</SUB>
       */
      double D_M (const double redshift) const;

      /**
       *  @brief the angular diameter distance at a given redshift
       *  @param redshift the redshift
       *  @return D<SUB>A</SUB>
       */
      double D_A (const double redshift) const; 

      /**
       *  @brief the angular diameter distance between objects at two
       *  redshifts
       *
       *  this function provides the angular diameter distance
       *  \f$D_{A,12}\f$ between two objects at redshifts \f$z_1\f$
       *  and \f$z_2\f$, as follows:
       *
       *  \f[ D_{A,12} = \frac{1}{1+z_2} \left[
       *  D_M(z_2)\sqrt{1+\Omega_k\frac{D_M^2(z_1)}{D_H^2}} -
       *  D_M(z_1)\sqrt{1+\Omega_k\frac{D_M^2(z_2)}{D_H^2}} \right]
       *  \f]
       *
       *  (e.g. Eq.19 of Hogg 2000)
       *
       *  @param z1 the redshift of the first object
       *  @param z2 the redshift of the second object
       *
       *  @return D<SUB>A</SUB> between z1 and z2
       *
       *  @warning the implemented formula is not correct for
       *  \f$\Omega_k<0\f$
       */
      double D_A (const double z1, const double z2) const;
      
      /**
       *  @brief the luminosity distance at a given redshift
       *  @param redshift the redshift
       *  @return D<SUB>L</SUB>
       */
      double D_L (const double redshift) const; 
  
      /**
       *  @brief the average distance at a given redshift, used to
       rescale the correlation function
       *  @param redshift the redshift
       *  @return D<SUB>V</SUB>
       */
      double D_V (const double redshift) const;

      /**
       *  @brief F_AP, the ALCOCK-PACZYNSKI distortion parameter
       *  @param redshift the redshift
       *  @return return F_AP
       */
      double F_AP (const double redshift) const;

      /**
       *  @brief the distance at a given redshift. Distance available
       *  are:
       *  D<SUB>C</SUB>,D<SUB>L</SUB>,D<SUB>A</SUB>,D<SUB>V</SUB>,
       *  D<SUB>V</SUB>/r<SUB>s</SUB>, r<SUB>s</SUB>/D<SUB>V</SUB>
       *  
       *  @author Alfonso Veropalumbo
       *  @author alfonso.veropalumbo@unibo.it
       *  @param redshift the redshift
       *  @param distance_type the type of distance to return 
       *  @return Distance
       */
      double Distance (const double redshift, const std::string distance_type) const;

      /**
       *  @brief comoving volume for a given redshift range and sky area
       *
       *  this function provides an approximated expression valid only
       *  for a LCDM model
       *
       *  @param z1 minimum redshift
       *  @param z2 maximum redshift
       *  @param Area sky area
       *  @return comoving volume
       */
      double Volume (const double z1, const double z2, const double Area) const; 

      /**
       *  @brief total comoving volume from z=0 to z
       *
       *  from Hogg 2000, Eq. 29
       *
       *  @author Alfonso Veropalumbo
       *  @author alfonso.veropalumbo@unibo.it
       *  @param zz redshift
       *  @return comoving volume
       */
      double Volume (const double zz) const;
  
      /**
       *  @brief maximum redshift for a given volume, sky area and
       *  minimum redshift
       *  @param Volume volume
       *  @param Area sky area
       *  @param z_min minimum redshift
       *  @return redshift
       */
      double max_redshift (const double Volume, const double Area, const double z_min) const; 

      /**
       *  @brief the derivative of the comoving volume,
       *  d<SUP>2</SUP>V/(dz*d&Omega;) at a given redshift
       *
       *  @param redshift the redshift
       *
       *  @param angle_rad false \f$\rightarrow\f$ &Omega; in square
       *  degrees; true \f$\rightarrow\f$ &Omega; in steradians
       *
       *  @return d<SUP>2</SUP>V/(dz*d&Omega;)
       */
      double dV_dZdOmega (const double redshift, const bool angle_rad) const; 

      ///@}


      /**
       *  @name Functions to estimate the halo mass function and related quantities
       */
      
      ///@{
      /**
       *  @brief auxiliary function to create a grid file with
       *  &sigma;(M) 
       *
       *  @param method_SS method used to compute the power spectrum
       *  and &sigma;(mass); valid method_SS are: CAMB
       *  [http://camb.info/], classgal_v1 [http://class-code.net/],
       *  EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       *
       *  @param redshift the redshift
       *
       *  @param output_root output_root of the parameter file used to
       *  compute the power spectrum and &sigma;(mass); it can be any
       *  name
       *
       *  @param interpType method to interpolate the power spectrum
       *     
       *  @param k_max maximum wave vector module up to which the power
       *  spectrum is computed
       *
       *  @param input_file either the parameter file or the power
       *  spectrum file; if a parameter file is provided,
       *  i.e. input_file!=NULL and is_parameter_file=true, it will be
       *  used to compute the power spectrum; if a power spectrum file
       *  is provided, i.e. input_file!=NULL and
       *  is_parameter_file=false, then the provided power spectrum
       *  will be used directly; in both cases &sigma;<SUP>2</SUP>(M)
       *  is computed by integrating the computed/provided power
       *  spectrum ignoring the cosmological parameters of the object
       *
       *  @param is_parameter_file true \f$\rightarrow\f$ the input_file
       *  is a parameter file, used to compute the power spectrum with
       *  the method specified by method_Pk; false \f$\rightarrow\f$
       *  the input_file is a file containing the power spectrum
       *
       *  @return file_grid name of the file where the grid is stored
       */
      std::string create_grid_sigmaM (const std::string method_SS, const double redshift, const std::string output_root="test", const std::string interpType="Linear", const double k_max=100., const std::string input_file=par::defaultString, const bool is_parameter_file=true) const;         

      /**
       *  @brief the mass function of dark matter haloes (filaments
       *  and sheets)
       *
       *  @author Alfonso Veropalumbo, Jacopo Neri (and Federico
       *  Marulli)
       *
       *  @author alfonso.veropalumbo@unibo.it, jacopo.neri6@gmail.com
       *  (and federico.marulli3@unibo.it)
       *
       *  @param Mass mass
       *
       *  @param redshift the redshift
       *
       *  @param model_MF author(s) who proposed the mass function;
       *  valid authors are: PS (Press & Schechter), ST (Sheth &
       *  Tormen), Jenkins (Jenkins et al. 2001), Warren (Warren et
       *  al. 2006), Reed (Reed et al. 2007), Pan (Pan 2007), ShenH
       *  (halo MF, Shen et al. 2006), ShenF (filament MF, Shen et
       *  al. 2006), ShenS (sheet MF, Shen et al. 2006), Tinker
       *  (Tinker et al. 2008), Crocce (Crocce et al. 2010),
       *  Angulo_FOF (FoF MF, Angulo et al. 2012), Angulo_Sub (SUBFIND
       *  MF, Angulo et al. 2012), Watson_FOF (FoF MF, Watson et
       *  al. 2012), Watson_SOH (Spherical Overdensity halo MF, Watson
       *  et al. 2012), Manera (Manera et al. 2010), Bhattacharya
       *  (Bhattacharya et al. 2011), Courtin (Courtin et al. 2010),
       *  Peacock (Peacock at al. 2007), Despali (Despali et al. 2016)
       *
       *  @param method_SS method used to compute the power spectrum
       *  and &sigma;(mass); valid method_SS are: CAMB
       *  [http://camb.info/], classgal_v1 [http://class-code.net/],
       *  EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       *
       *  @param output_root output_root of the parameter file used to
       *  compute the power spectrum and &sigma;(mass); it can be any
       *  name
       *
       *  @param Delta \f$\Delta\f$: the overdensity, defined as the mean
       *  interior density relative to the background
       *
       *  @param interpType method to interpolate the power spectrum
       *
       *  @param norm 0 \f$\rightarrow\f$ don't normalise the power
       *  spectrum; 1 \f$\rightarrow\f$ normalise the power spectrum;
       *  -1 \f$\rightarrow\f$ normalise only if sigma8 is set
       *
       *  @param k_min minimum wave vector module up to which the power
       *  spectrum is computed
       *
       *  @param k_max maximum wave vector module up to which the power
       *  spectrum is computed
       *
       *  @param prec accuracy of the integration 
       *
       *  @param input_file either the parameter file or the power
       *  spectrum file; if a parameter file is provided,
       *  i.e. input_file!=NULL and is_parameter_file=true, it will be
       *  used to compute the power spectrum; if a power spectrum file
       *  is provided, i.e. input_file!=NULL and
       *  is_parameter_file=false, then the provided power spectrum
       *  will be used directly; in both cases &sigma;<SUP>2</SUP>(M)
       *  is computed by integrating the computed/provided power
       *  spectrum ignoring the cosmological parameters of the object
       *
       *  @param is_parameter_file true \f$\rightarrow\f$ the input_file
       *  is a parameter file, used to compute the power spectrum with
       *  the method specified by method_Pk; false \f$\rightarrow\f$
       *  the input_file is a file containing the power spectrum
       *
       *  @param default_delta true = using function
       *  cbl::cosmology::deltac; false = using delta_t*growth
       *  factor
       *  
       *  @param delta_t user defined density contrast at \f$z = 0\f$
       *
       *  @return the mass function, d&Phi;/dM=dn(M)/dM
       *
       *  @warning the input parameter \f$\Delta\f$ is the background
       *  overdensity, not the critical overdensity
       *  \f$\Delta_{crit}\f$; the function
       *  cbl::Cosmology::Delta_vir can be used to convert
       *  \f$\Delta_{crit}\f$ into \f$\Delta\f$
       *
       *  @warning the mass function by Manera et al. (2010) has been
       *  tested only for z=0 and z=0.5; the mass function by Despali
       *  et al. (2016) is currently implemented only for virial
       *  masses and at \f$z<1.25\f$
       */
      double mass_function (const double Mass, const double redshift, const std::string model_MF, const std::string method_SS, const std::string output_root="test", const double Delta=200., const std::string interpType="Linear", const int norm=-1, const double k_min=0., const double k_max=100., const double prec=1.e-2, const std::string input_file=par::defaultString, const bool is_parameter_file=true, const bool default_delta=true, const double delta_t=1.686);

      /**
       *  @brief the mass function of dark matter haloes (filaments and
       *  sheets) computed quickly using a grid
       *
       *  @author Alfonso Veropalumbo, Jacopo Neri (and Federico
       *  Marulli)
       *
       *  @author alfonso.veropalumbo@unibo.it, jacopo.neri6@gmail.com
       *  (and federico.marulli3@unibo.it)
       *
       *  @param Mass mass
       *
       *  @param redshift the redshift
       *
       *  @param model_MF author(s) who proposed the mass function;
       *  valid authors are: PS (Press & Schechter), ST (Sheth &
       *  Tormen), Jenkins (Jenkins et al. 2001), Warren (Warren et
       *  al. 2006), Reed (Reed et al. 2007), Pan (Pan 2007), ShenH
       *  (halo MF, Shen et al. 2006), ShenF (filament MF, Shen et
       *  al. 2006), ShenS (sheet MF, Shen et al. 2006), Tinker
       *  (Tinker et al. 2008), Crocce (Crocce et al. 2010),
       *  Angulo_FOF (FoF MF, Angulo et al. 2012), Angulo_Sub (SUBFIND
       *  MF, Angulo et al. 2012), Watson_FOF (FoF MF, Watson et
       *  al. 2012), Watson_SOH (Spherical Overdensity halo MF, Watson
       *  et al. 2012), Manera (Manera et al. 2010), Bhattacharya
       *  (Bhattacharya et al. 2011), Courtin (Courtin et al. 2010),
       *  Peacock (Peacock at al. 2007), Despali (Despali et al. 2016)
       *
       *  @param method_SS method used to compute the power spectrum
       *  and &sigma;(mass); valid method_SS are: CAMB
       *  [http://camb.info/], classgal_v1 [http://class-code.net/],
       *  EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       *
       *  @param output_root output_root of the parameter file used to
       *  compute the power spectrum and &sigma;(mass); it can be any
       *  name
       *
       *  @param Delta \f$\Delta\f$: the overdensity, defined as the mean
       *  interior density relative to the background
       *
       *  @param interpType method to interpolate the power spectrum
       *
       *  @param norm 0 \f$\rightarrow\f$ don't normalise the power
       *  spectrum; 1 \f$\rightarrow\f$ normalise the power spectrum;
       *  -1 \f$\rightarrow\f$ normalise only if sigma8 is set
       *
       *  @param k_min minimum wave vector module up to which the power
       *  spectrum is computed
       *
       *  @param k_max maximum wave vector module up to which the power
       *  spectrum is computed
       *
       *  @param prec accuracy of the integration 
       *
       *  @param input_file either the parameter file or the power
       *  spectrum file; if a parameter file is provided,
       *  i.e. input_file!=NULL and is_parameter_file=true, it will be
       *  used to compute the power spectrum; if a power spectrum file
       *  is provided, i.e. input_file!=NULL and
       *  is_parameter_file=false, then the provided power spectrum
       *  will be used directly; in both cases &sigma;<SUP>2</SUP>(M)
       *  is computed by integrating the computed/provided power
       *  spectrum ignoring the cosmological parameters of the object
       *
       *  @param is_parameter_file true \f$\rightarrow\f$ the input_file
       *  is a parameter file, used to compute the power spectrum with
       *  the method specified by method_Pk; false \f$\rightarrow\f$
       *  the input_file is a file containing the power spectrum
       *
       *  @return the mass function, d&Phi;/dM=dn(M)/dM
       *
       *  @warning the input parameter \f$\Delta\f$ is the background
       *  overdensity, not the critical overdensity
       *  \f$\Delta_{crit}\f$; the function
       *  cbl::Cosmology::Delta_vir can be used to convert
       *  \f$\Delta_{crit}\f$ into \f$\Delta\f$
       *
       *  @warning the mass function by Manera et al. (2010) has been
       *  tested only for z=0 and z=0.5; the mass function by Despali
       *  et al. (2016) is currently implemented only for virial
       *  masses and at \f$z<1.25\f$
       */
      double mass_function_fast (const double Mass, const double redshift, const std::string model_MF, const std::string method_SS, const std::string output_root="test", const double Delta=200., const std::string interpType="Linear", const int norm=-1, const double k_min=0., const double k_max=100., const double prec=1.e-2, const std::string input_file=par::defaultString, const bool is_parameter_file=true); 

      /**
       *  @brief the mass function of dark matter haloes (filaments and
       *  sheets) computed quickly passing directly the mass variance
       *  and its derivative as inputs
       *
       *  @author Alfonso Veropalumbo, Jacopo Neri (and Federico
       *  Marulli)
       *
       *  @author alfonso.veropalumbo@unibo.it, jacopo.neri6@gmail.com
       *  (and federico.marulli3@unibo.it)
       *
       *  @param Mass mass
       *
       *  @param Sigma &sigma;(mass): the mass variance
       *
       *  @param Dln_Sigma dln&sigma;/dM: the derivative of the mass
       *  variance
       *
       *  @param redshift the redshift
       *
       *  @param model_MF author(s) who proposed the mass function;
       *  valid authors are: PS (Press&Schechter), ST (Sheth&Tormen),
       *  Jenkins (Jenkins et al. 2001), Warren (Warren et al. 2006),
       *  Reed, (Reed et al. 2007), Pan (Pan 2007), ShenH (halo MF by
       *  Shen et al. 2006), ShenF (filaments MF by Shen et al. 2006),
       *  ShenS (sheets MF by Shen et al. 2006), Tinker (Tinker et
       *  al. 2008), Crocce (Crocce et al. 2010), Angulo_FOF (FOF MF by
       *  Angulo et al. 2012), Angulo_Sub (SUBFIND MF by Angulo et
       *  al. 2012), Watson_FOF(FOF MF by Watson et al. 2012),
       *  Watson_SOH (MF for Spherical Overdensity Haloes by Watson 
       *  et al. 2012), Manera (Manera et al. 2010), Bhattacharya 
       *  (Bhattacharya et al. 2011), Courtin (Courtin et al. 2010),
       *  Peacock (by Peacock at al. 2007)
       * 
       *  @param output_root output_root of the parameter file used to
       *  compute the power spectrum and &sigma;(mass); it can be any
       *  name
       *
       *  @param Delta \f$\Delta\f$: the overdensity, defined as the mean
       *  interior density relative to the background
       *
       *  @param interpType method to interpolate the power spectrum
       *
       *  @param norm 0 \f$\rightarrow\f$ don't normalise the power
       *  spectrum; 1 \f$\rightarrow\f$ normalise the power spectrum;
       *  -1 \f$\rightarrow\f$ normalise only if sigma8 is set
       *
       *  @param k_min minimum wave vector module up to which the power
       *  spectrum is computed
       *
       *  @param k_max maximum wave vector module up to which the power
       *  spectrum is computed
       *
       *  @param prec accuracy of the integration 
       *
       *  @param method_SS method used to compute the power spectrum
       *  and &sigma;(mass); valid method_SS are: CAMB
       *  [http://camb.info/], classgal_v1 [http://class-code.net/],
       *  EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       *
       *  @param input_file either the parameter file or the power
       *  spectrum file; if a parameter file is provided,
       *  i.e. input_file!=NULL and is_parameter_file=true, it will be
       *  used to compute the power spectrum; if a power spectrum file
       *  is provided, i.e. input_file!=NULL and
       *  is_parameter_file=false, then the provided power spectrum
       *  will be used directly; in both cases &sigma;<SUP>2</SUP>(M)
       *  is computed by integrating the computed/provided power
       *  spectrum ignoring the cosmological parameters of the object
       *
       *  @param is_parameter_file true \f$\rightarrow\f$ the input_file
       *  is a parameter file, used to compute the power spectrum with
       *  the method specified by method_Pk; false \f$\rightarrow\f$
       *  the input_file is a file containing the power spectrum
       *
       *  @return the mass function, d&Phi;/dM=dn(M)/dM
       *
       *  @warning the input parameter \f$\Delta\f$ is the background
       *  overdensity, not the critical overdensity
       *  \f$\Delta_{crit}\f$; the function
       *  cbl::Cosmology::Delta_vir can be used to convert
       *  \f$\Delta_{crit}\f$ into \f$\Delta\f$
       *
       *  @warning the mass function by Manera et al. (2010) has been
       *  tested only for z=0 and z=0.5; the mass function by Despali
       *  et al. (2016) is currently implemented only for virial
       *  masses and at \f$z<1.25\f$
       */
      double mass_function (const double Mass, const double Sigma, const double Dln_Sigma, const double redshift, const std::string model_MF, const std::string output_root="test", const double Delta=200., const std::string interpType="Linear", const int norm=-1, const double k_min=0., const double k_max=100., const double prec=1.e-2, const std::string method_SS="CAMB", const std::string input_file=par::defaultString, const bool is_parameter_file=true); 
    
      /**
       *  @brief number of dark matter haloes per steradian or square
       *  degree, for a given redshift range
       *
       *  @author Alfonso Veropalumbo, Jacopo Neri (and Federico
       *  Marulli)
       *
       *  @author alfonso.veropalumbo@unibo.it, jacopo.neri6@gmail.com
       *  (and federico.marulli3@unibo.it)
       *
       *  @param Mass_min minimum mass
       *
       *  @param Mass_max maximum mass
       *
       *  @param z_min minimum redshift
       *
       *  @param z_max maximum redshift
       *
       *  @param angle_rad 0 \f$\rightarrow\f$ &Omega; in square
       *  degrees; 1 \f$\rightarrow\f$ &Omega; in steradians
       *
       *  @param model_MF author(s) who proposed the mass function;
       *  valid authors are: PS (Press&Schechter), ST (Sheth&Tormen),
       *  Jenkins (Jenkins et al. 2001), Warren (Warren et al. 2006),
       *  Reed, (Reed et al. 2007), Pan (Pan 2007), ShenH (halo MF by
       *  Shen et al. 2006), ShenF (filaments MF by Shen et al. 2006),
       *  ShenS (sheets MF by Shen et al. 2006), Tinker (Tinker et
       *  al. 2008), Crocce (Crocce et al. 2010), Angulo_FOF (FOF MF by
       *  Angulo et al. 2012), Angulo_Sub (SUBFIND MF by Angulo et
       *  al. 2012), Watson_FOF(FOF MF by Watson et al. 2012),
       *  Watson_SOH (MF for Spherical Overdensity Haloes by Watson 
       *  et al. 2012), Manera (Manera et al. 2010), Bhattacharya 
       *  (Bhattacharya et al. 2011), Courtin (Courtin et al. 2010),
       *  Peacock (by Peacock at al. 2007)
       *
       *  @param method_SS method used to compute the power spectrum and
       *  &sigma;(mass); valid method_SS are: CAMB [http://camb.info/],
       *  classgal_v1 [http://class-code.net/], EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       *
       *  @param output_root output_root of the parameter file used to
       *  compute the power spectrum and &sigma;(mass); it can be any
       *  name
       *       
       *  @param Delta \f$\Delta\f$: the overdensity, defined as the
       *  mean interior density relative to the background
       *
       *  @param interpType method to interpolate the power spectrum
       *
       *  @param k_max maximum wave vector module up to which the power
       *  spectrum is computed
       *
       *  @param input_file either the parameter file or the power
       *  spectrum file; if a parameter file is provided,
       *  i.e. input_file!=NULL and is_parameter_file=true, it will be
       *  used to compute the power spectrum; if a power spectrum file
       *  is provided, i.e. input_file!=NULL and
       *  is_parameter_file=false, then the provided power spectrum
       *  will be used directly; in both cases &sigma;<SUP>2</SUP>(M)
       *  is computed by integrating the computed/provided power
       *  spectrum ignoring the cosmological parameters of the object
       *
       *  @param is_parameter_file true \f$\rightarrow\f$ the input_file
       *  is a parameter file, used to compute the power spectrum with
       *  the method specified by method_Pk; false \f$\rightarrow\f$
       *  the input_file is a file containing the power spectrum
       *
       *  @return n<SUB>haloes</SUB>: the number density of dark matter
       *  haloes (per steradian or square degree)
       *
       *  @warning the input parameter \f$\Delta\f$ is the background
       *  overdensity, not the critical overdensity
       *  \f$\Delta_{crit}\f$; the function
       *  cbl::Cosmology::Delta_vir can be used to convert
       *  \f$\Delta_{crit}\f$ into \f$\Delta\f$
       *
       *  @warning the mass function by Manera et al. (2010) has been
       *  tested only for z=0 and z=0.5; the mass function by Despali
       *  et al. (2016) is currently implemented only for virial
       *  masses and at \f$z<1.25\f$
       */
      double n_haloes (const double Mass_min, const double Mass_max, const double z_min, const double z_max, const bool angle_rad, const std::string model_MF, const std::string method_SS, const std::string output_root="test", const double Delta=200, const std::string interpType="Linear", const double k_max=100., const std::string input_file=par::defaultString, const bool is_parameter_file=true);
      
      /**
       *  @brief number of dark matter haloes per volume at fixed
       *  redshift 
       *
       *  this function computes the number of dark matter haloes per
       *  volume at fixed redshift as follows:
       *
       *  \f[ N_h = \int_{M_{min}}^{M_{max}} d M \Phi(M)\f]
       *
       *  @author Alfonso Veropalumbo, Jacopo Neri (and Federico
       *  Marulli)
       *
       *  @author alfonso.veropalumbo@unibo.it, jacopo.neri6@gmail.com
       *  (and federico.marulli3@unibo.it)
       *
       *  @param Mass_min minimum mass
       *
       *  @param Mass_max maximum mass
       *
       *  @param Volume the volume
       *
       *  @param redshift the redshift
       *
       *  @param model_MF author(s) who proposed the mass function;
       *  valid authors are: PS (Press & Schechter), ST (Sheth &
       *  Tormen), Jenkins (Jenkins et al. 2001), Warren (Warren et
       *  al. 2006), Reed, (Reed et al. 2007), Pan (Pan 2007), ShenH
       *  (halo MF by Shen et al. 2006), ShenF (filaments MF by Shen
       *  et al. 2006), ShenS (sheets MF by Shen et al. 2006), Tinker
       *  (Tinker et al. 2008), Crocce (Crocce et al. 2010),
       *  Angulo_FOF (FOF MF by Angulo et al. 2012), Angulo_Sub
       *  (SUBFIND MF by Angulo et al. 2012), Watson_FOF(FOF MF by
       *  Watson et al. 2012), Watson_SOH (MF for Spherical Overdensity
       *  Haloes by Watson et al. 2012), Manera (Manera et al. 2010),
       *  Bhattacharya (Bhattacharya et al. 2011), Courtin (Courtin 
       *  et al. 2010), Peacock (by Peacock at al. 2007)
       *
       *  @param method_SS method used to compute the power spectrum
       *  and &sigma;(mass); valid method_SS are: CAMB
       *  [http://camb.info/], classgal_v1 [http://class-code.net/],
       *  EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       *
       *  @param nbin_mass number of bin for the mass function
       *  computation
       *
       *  @param output_root output_root of the parameter file used to
       *  compute the power spectrum and &sigma;(mass); it can be any
       *  name
       *
       *  @param Delta \f$\Delta\f$: the overdensity, defined as the
       *  mean interior density relative to the background
       *
       *  @param interpType method to interpolate the power spectrum
       *
       *  @param norm 0 \f$\rightarrow\f$ don't normalise the power
       *  spectrum; 1 \f$\rightarrow\f$ normalise the power spectrum;
       *  -1 \f$\rightarrow\f$ normalise only if sigma8 is set
       *
       *  @param k_min minimum wave vector module up to which the power
       *  spectrum is computed
       *
       *  @param k_max maximum wave vector module up to which the power
       *  spectrum is computed
       *
       *  @param prec accuracy of the integration 
       *
       *  @param input_file either the parameter file or the power
       *  spectrum file; if a parameter file is provided,
       *  i.e. input_file!=NULL and is_parameter_file=true, it will be
       *  used to compute the power spectrum; if a power spectrum file
       *  is provided, i.e. input_file!=NULL and
       *  is_parameter_file=false, then the provided power spectrum
       *  will be used directly; in both cases &sigma;<SUP>2</SUP>(M)
       *  is computed by integrating the computed/provided power
       *  spectrum ignoring the cosmological parameters of the object
       *
       *  @param is_parameter_file true \f$\rightarrow\f$ the input_file
       *  is a parameter file, used to compute the power spectrum with
       *  the method specified by method_Pk; false \f$\rightarrow\f$
       *  the input_file is a file containing the power spectrum
       *
       *  @param default_delta true = using function
       *  cbl::cosmology::deltac; false = using delta_t*growth
       *  factor
       *  
       *  @param delta_t user defined density contrast at \f$z = 0\f$
       *
       *  @return the mass function, d&Phi;/dM=dn(M)/dM
       *
       *  @warning the input parameter \f$\Delta\f$ is the background
       *  overdensity, not the critical overdensity
       *  \f$\Delta_{crit}\f$; the function
       *  cbl::Cosmology::Delta_vir can be used to convert
       *  \f$\Delta_{crit}\f$ into \f$\Delta\f$
       *
       *  @warning the mass function by Manera et al. (2010) has been
       *  tested only for z=0 and z=0.5; the mass function by Despali
       *  et al. (2016) is currently implemented only for virial
       *  masses and at \f$z<1.25\f$
       */
      double n_haloes (const double Mass_min, const double Mass_max, const double Volume, const double redshift, const std::string model_MF, const std::string method_SS, const int nbin_mass=0, const std::string output_root="test", const double Delta=200., const std::string interpType="Linear", const int norm=-1, const double k_min=0., const double k_max=100., const double prec=1.e-2, const std::string input_file=par::defaultString, const bool is_parameter_file=true, const bool default_delta=true, const double delta_t=1.686);
      
      /**
       *  @brief number of dark matter haloes per steradian or square
       *  degree, for a given redshift range and with selection function
       *  defined on a grid
       *
       *  @author Alfonso Veropalumbo, Jacopo Neri (and Federico
       *  Marulli)
       *
       *  @author alfonso.veropalumbo@unibo.it, jacopo.neri6@gmail.com
       *  (and federico.marulli3@unibo.it)
       *
       *  @param Mass_min minimum mass
       *
       *  @param Mass_max maximum mass
       *
       *  @param z_min minimum redshift
       *
       *  @param z_max maximum redshift
       *
       *  @param angle_rad 0 \f$\rightarrow\f$ &Omega; in square
       *  degrees; 1 \f$\rightarrow\f$ &Omega; in steradians
       *
       *  @param model_MF author(s) who proposed the mass function;
       *  valid authors are: PS (Press&Schechter), ST (Sheth&Tormen),
       *  Jenkins (Jenkins et al. 2001), Warren (Warren et al. 2006),
       *  Reed, (Reed et al. 2007), Pan (Pan 2007), ShenH (halo MF by
       *  Shen et al. 2006), ShenF (filaments MF by Shen et al. 2006),
       *  ShenS (sheets MF by Shen et al. 2006), Tinker (Tinker et
       *  al. 2008), Crocce (Crocce et al. 2010), Angulo_FOF (FOF MF by
       *  Angulo et al. 2012), Angulo_Sub (SUBFIND MF by Angulo et
       *  al. 2012), Watson_FOF(FOF MF by Watson et al. 2012),
       *  Watson_SOH (MF for Spherical Overdensity Haloes by Watson 
       *  et al. 2012), Manera (Manera et al. 2010), Bhattacharya 
       *  (Bhattacharya et al. 2011), Courtin (Courtin et al. 2010),
       *  Peacock (by Peacock at al. 2007)
       *
       *  @param method_SS method used to compute the power spectrum and
       *  &sigma;(mass); valid method_SS are: CAMB [http://camb.info/],
       *  classgal_v1 [http://class-code.net/], EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       *
       *  @param selection_function_file input file where the selection
       *  function is stored
       *
       *  @param [in] column the columns to be read
       *
       *  @param output_root output_root of the parameter file used to
       *  compute the power spectrum and &sigma;(mass); it can be any
       *  name
       *
       *  @param Delta \f$\Delta\f$: the overdensity, defined as the
       *  mean interior density relative to the background
       *
       *  @param isDelta_vir \f$\rightarrow\f$ \f$\Delta\f$ is the
       *  virial overdensity
       *
       *  @param interpType method to interpolate the power spectrum
       *
       *  @param k_max maximum wave vector module up to which the power
       *  spectrum is computed
       *
       *  @param input_file either the parameter file or the power
       *  spectrum file; if a parameter file is provided,
       *  i.e. input_file!=NULL and is_parameter_file=true, it will be
       *  used to compute the power spectrum; if a power spectrum file
       *  is provided, i.e. input_file!=NULL and
       *  is_parameter_file=false, then the provided power spectrum
       *  will be used directly; in both cases &sigma;<SUP>2</SUP>(M)
       *  is computed by integrating the computed/provided power
       *  spectrum ignoring the cosmological parameters of the object
       *
       *  @param is_parameter_file true \f$\rightarrow\f$ the input_file
       *  is a parameter file, used to compute the power spectrum with
       *  the method specified by method_Pk; false \f$\rightarrow\f$
       *  the input_file is a file containing the power spectrum
       *
       *  @return n<SUB>haloes</SUB>: the number density of dark matter
       *  haloes (per steradian or square degree)
       *
       *  @warning the input parameter \f$\Delta\f$ is the background
       *  overdensity, not the critical overdensity
       *  \f$\Delta_{crit}\f$; the function
       *  cbl::Cosmology::Delta_vir can be used to convert
       *  \f$\Delta_{crit}\f$ into \f$\Delta\f$
       *
       *  @warning the mass function by Manera et al. (2010) has been
       *  tested only for z=0 and z=0.5; the mass function by Despali
       *  et al. (2016) is currently implemented only for virial
       *  masses and at \f$z<1.25\f$
       */
      double n_haloes_selection_function (const double Mass_min, const double Mass_max, const double z_min, const double z_max, const bool angle_rad, const std::string model_MF, const std::string method_SS, const std::string selection_function_file, const std::vector<int> column={}, const std::string output_root="test", const double Delta=200, const bool isDelta_vir=false, const std::string interpType="Linear", const double k_max=100., const std::string input_file=par::defaultString, const bool is_parameter_file=true);

      /**
       *  @brief mass function for a range of masses
       *
       *  @author Alfonso Veropalumbo
       *
       *  @author alfonso.veropalumbo@unibo.it
       *
       *  @param mass vector of mass
       *
       *  @param z_min minimum redshift
       *
       *  @param z_max maximum redshift
       *
       *  @param model_MF author(s) who proposed the mass function;
       *  valid authors are: PS (Press&Schechter), ST (Sheth&Tormen),
       *  Jenkins (Jenkins et al. 2001), Warren (Warren et al. 2006),
       *  Reed, (Reed et al. 2007), Pan (Pan 2007), ShenH (halo MF by
       *  Shen et al. 2006), ShenF (filaments MF by Shen et al. 2006),
       *  ShenS (sheets MF by Shen et al. 2006), Tinker (Tinker et
       *  al. 2008), Crocce (Crocce et al. 2010), Angulo_FOF (FOF MF by
       *  Angulo et al. 2012), Angulo_Sub (SUBFIND MF by Angulo et
       *  al. 2012), Watson_FOF(FOF MF by Watson et al. 2012),
       *  Watson_SOH (MF for Spherical Overdensity Haloes by Watson 
       *  et al. 2012), Manera (Manera et al. 2010), Bhattacharya 
       *  (Bhattacharya et al. 2011), Courtin (Courtin et al. 2010),
       *  Peacock (by Peacock at al. 2007)
       *
       *  @param method_SS method used to compute the power spectrum and
       *  &sigma;(mass); valid method_SS are: CAMB [http://camb.info/],
       *  classgal_v1 [http://class-code.net/], EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       *
       *  @param output_root output_root of the parameter file used to
       *  compute the power spectrum and &sigma;(mass); it can be any
       *  name
       *
       *  @param Delta \f$\Delta\f$: the overdensity, defined as the
       *  mean interior density relative to the background
       *
       *  @param isDelta_vir \f$\rightarrow\f$ \f$\Delta\f$ is the
       *  virial overdensity
       *
       *  @param interpType method to interpolate the power spectrum
       *
       *  @param k_max maximum wave vector module up to which the power
       *  spectrum is computed
       *
       *  @param input_file either the parameter file or the power
       *  spectrum file; if a parameter file is provided,
       *  i.e. input_file!=NULL and is_parameter_file=true, it will be
       *  used to compute the power spectrum; if a power spectrum file
       *  is provided, i.e. input_file!=NULL and
       *  is_parameter_file=false, then the provided power spectrum
       *  will be used directly; in both cases &sigma;<SUP>2</SUP>(M)
       *  is computed by integrating the computed/provided power
       *  spectrum ignoring the cosmological parameters of the object
       *
       *  @param is_parameter_file true \f$\rightarrow\f$ the input_file
       *  is a parameter file, used to compute the power spectrum with
       *  the method specified by method_Pk; false \f$\rightarrow\f$
       *  the input_file is a file containing the power spectrum
       *
       *  @return n<SUB>haloes</SUB>: the number density of dark matter
       *  haloes (per steradian or square degree)
       *
       *  @warning the input parameter \f$\Delta\f$ is the background
       *  overdensity, not the critical overdensity
       *  \f$\Delta_{crit}\f$; the function
       *  cbl::Cosmology::Delta_vir can be used to convert
       *  \f$\Delta_{crit}\f$ into \f$\Delta\f$
       *
       *  @warning the mass function by Manera et al. (2010) has been
       *  tested only for z=0 and z=0.5; the mass function by Despali
       *  et al. (2016) is currently implemented only for virial
       *  masses and at \f$z<1.25\f$
       */
      std::vector<double> mass_function (const std::vector<double> mass, const double z_min, const double z_max, const std::string model_MF, const std::string method_SS, const std::string output_root="test", const double Delta=200, const bool isDelta_vir=false, const std::string interpType="Linear", const double k_max=100., const std::string input_file=par::defaultString, const bool is_parameter_file=true);

      /**
       *  @brief mass function given a selection function
       *
       *  @author Alfonso Veropalumbo
       *
       *  @author alfonso.veropalumbo@unibo.it
       *
       *  @param mass vector of mass
       *
       *  @param z_min minimum redshift
       *
       *  @param z_max maximum redshift
       *
       *  @param model_MF author(s) who proposed the mass function;
       *  valid authors are: PS (Press&Schechter), ST (Sheth&Tormen),
       *  Jenkins (Jenkins et al. 2001), Warren (Warren et al. 2006),
       *  Reed, (Reed et al. 2007), Pan (Pan 2007), ShenH (halo MF by
       *  Shen et al. 2006), ShenF (filaments MF by Shen et al. 2006),
       *  ShenS (sheets MF by Shen et al. 2006), Tinker (Tinker et
       *  al. 2008), Crocce (Crocce et al. 2010), Angulo_FOF (FOF MF by
       *  Angulo et al. 2012), Angulo_Sub (SUBFIND MF by Angulo et
       *  al. 2012), Watson_FOF(FOF MF by Watson et al. 2012),
       *  Watson_SOH (MF for Spherical Overdensity Haloes by Watson 
       *  et al. 2012), Manera (Manera et al. 2010), Bhattacharya 
       *  (Bhattacharya et al. 2011), Courtin (Courtin et al. 2010),
       *  Peacock (by Peacock at al. 2007)
       *
       *  @param method_SS method used to compute the power spectrum and
       *  &sigma;(mass); valid method_SS are: CAMB [http://camb.info/],
       *  classgal_v1 [http://class-code.net/], EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       *
       *  @param selection_function_file input file where the selection
       *  function is stored
       *
       *  @param column the columns to be read
       *
       *  @param output_root output_root of the parameter file used to
       *  compute the power spectrum and &sigma;(mass); it can be any
       *  name
       *
       *  @param Delta \f$\Delta\f$: the overdensity, defined as the
       *  mean interior density relative to the background
       *
       *  @param isDelta_vir \f$\rightarrow\f$ \f$\Delta\f$ is the
       *  virial overdensity
       *
       *  @param interpType method to interpolate the power spectrum
       *
       *  @param k_max maximum wave vector module up to which the power
       *  spectrum is computed
       *
       *  @param input_file either the parameter file or the power
       *  spectrum file; if a parameter file is provided,
       *  i.e. input_file!=NULL and is_parameter_file=true, it will be
       *  used to compute the power spectrum; if a power spectrum file
       *  is provided, i.e. input_file!=NULL and
       *  is_parameter_file=false, then the provided power spectrum
       *  will be used directly; in both cases &sigma;<SUP>2</SUP>(M)
       *  is computed by integrating the computed/provided power
       *  spectrum ignoring the cosmological parameters of the object
       *
       *  @param is_parameter_file true \f$\rightarrow\f$ the input_file
       *  is a parameter file, used to compute the power spectrum with
       *  the method specified by method_Pk; false \f$\rightarrow\f$
       *  the input_file is a file containing the power spectrum
       *
       *  @return n<SUB>haloes</SUB>: the number density of dark matter
       *  haloes (per steradian or square degree)
       *
       *  @warning the input parameter \f$\Delta\f$ is the background
       *  overdensity, not the critical overdensity
       *  \f$\Delta_{crit}\f$; the function
       *  cbl::Cosmology::Delta_vir can be used to convert
       *  \f$\Delta_{crit}\f$ into \f$\Delta\f$
       *
       *  @warning the mass function by Manera et al. (2010) has been
       *  tested only for z=0 and z=0.5; the mass function by Despali
       *  et al. (2016) is currently implemented only for virial
       *  masses and at \f$z<1.25\f$
       */
      std::vector<double> mass_function_selection_function_vector (const std::vector<double> mass, const double z_min, const double z_max, const std::string model_MF, const std::string method_SS, const std::string selection_function_file, const std::vector<int> column={}, const std::string output_root="test", const double Delta=200, const bool isDelta_vir=false, const std::string interpType="Linear", const double k_max=100., const std::string input_file=par::defaultString, const bool is_parameter_file=true);

      /**
       *  @brief redshift distribution of dark matter haloes
       *
       *  @author Alfonso Veropalumbo
       *
       *  @author alfonso.veropalumbo@unibo.it
       *
       *  @param z_min minimum redshift
       *
       *  @param z_max maximum redshift
       *
       *  @param step_z redshift step
       *
       *  @param Area_degrees the survey area, in degrees
       *
       *  @param Mass_min minimum halo mass
       *
       *  @param Mass_max maximum halo mass
       *
       *  @param model_MF author(s) who proposed the mass function;
       *  valid authors are: PS (Press&Schechter), ST (Sheth&Tormen),
       *  Jenkins (Jenkins et al. 2001), Warren (Warren et al. 2006),
       *  Reed, (Reed et al. 2007), Pan (Pan 2007), ShenH (halo MF by
       *  Shen et al. 2006), ShenF (filaments MF by Shen et al. 2006),
       *  ShenS (sheets MF by Shen et al. 2006), Tinker (Tinker et
       *  al. 2008), Crocce (Crocce et al. 2010), Angulo_FOF (FOF MF
       *  by Angulo et al. 2012), Angulo_Sub (SUBFIND MF by Angulo et
       *  al. 2012), Watson_FOF(FOF MF by Watson et al. 2012),
       *  Watson_SOH (MF for Spherical Overdensity Haloes by Watson et
       *  al. 2012), Manera (Manera et al. 2010), Bhattacharya
       *  (Bhattacharya et al. 2011), Courtin (Courtin et al. 2010),
       *  Peacock (by Peacock at al. 2007)
       *
       *  @param method_SS method used to compute the power spectrum
       *  and &sigma;(mass); valid method_SS are: CAMB
       *  [http://camb.info/], classgal_v1 [http://class-code.net/],
       *  EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       *
       *  @param output_root output_root of the parameter file used to
       *  compute the power spectrum and &sigma;(mass); it can be any
       *  name
       *
       *  @param Delta \f$\Delta\f$: the overdensity, defined as the
       *  mean interior density relative to the background
       *
       *  @param isDelta_vir \f$\rightarrow\f$ \f$\Delta\f$ is the
       *  virial overdensity
       *
       *  @param interpType method to interpolate the power spectrum
       *
       *  @param k_max maximum wave vector module up to which the power
       *  spectrum is computed
       *
       *  @param input_file either the parameter file or the power
       *  spectrum file; if a parameter file is provided,
       *  i.e. input_file!=NULL and is_parameter_file=true, it will be
       *  used to compute the power spectrum; if a power spectrum file
       *  is provided, i.e. input_file!=NULL and
       *  is_parameter_file=false, then the provided power spectrum
       *  will be used directly; in both cases &sigma;<SUP>2</SUP>(M)
       *  is computed by integrating the computed/provided power
       *  spectrum ignoring the cosmological parameters of the object
       *
       *  @param is_parameter_file true \f$\rightarrow\f$ the input_file
       *  is a parameter file, used to compute the power spectrum with
       *  the method specified by method_Pk; false \f$\rightarrow\f$
       *  the input_file is a file containing the power spectrum
       *
       *  @return the redshift distribution of dark matter haloes
       *
       *  @warning the input parameter \f$\Delta\f$ is the background
       *  overdensity, not the critical overdensity
       *  \f$\Delta_{crit}\f$; the function
       *  cbl::Cosmology::Delta_vir can be used to convert
       *  \f$\Delta_{crit}\f$ into \f$\Delta\f$
       *
       *  @warning the mass function by Manera et al. (2010) has been
       *  tested only for z=0 and z=0.5; the mass function by Despali
       *  et al. (2016) is currently implemented only for virial
       *  masses and at \f$z<1.25\f$
       */
      std::vector<double> redshift_distribution_haloes (const double z_min, const double z_max, const int step_z, const double Area_degrees, const double Mass_min, const double Mass_max, const std::string model_MF, const std::string method_SS, const std::string output_root="test", const double Delta=200, const bool isDelta_vir=false, const std::string interpType="Linear", const double k_max=100., const std::string input_file=par::defaultString, const bool is_parameter_file=true);

      /**
       *  @brief redshift distribution of dark matter haloes, given a
       *  selection function
       *
       *  @author Alfonso Veropalumbo
       *
       *  @author alfonso.veropalumbo@unibo.it
       *
       *  @param redshift vector containing the redshift at which the
       *  halo distribution will be computed
       *
       *  @param Area_degrees the survey area, in degrees
       *
       *  @param Mass_min minimum halo mass
       *
       *  @param Mass_max maximum halo mass
       *
       *  @param model_MF author(s) who proposed the mass function;
       *  valid authors are: PS (Press&Schechter), ST (Sheth&Tormen),
       *  Jenkins (Jenkins et al. 2001), Warren (Warren et al. 2006),
       *  Reed, (Reed et al. 2007), Pan (Pan 2007), ShenH (halo MF by
       *  Shen et al. 2006), ShenF (filaments MF by Shen et al. 2006),
       *  ShenS (sheets MF by Shen et al. 2006), Tinker (Tinker et
       *  al. 2008), Crocce (Crocce et al. 2010), Angulo_FOF (FOF MF
       *  by Angulo et al. 2012), Angulo_Sub (SUBFIND MF by Angulo et
       *  al. 2012), Watson_FOF(FOF MF by Watson et al. 2012),
       *  Watson_SOH (MF for Spherical Overdensity Haloes by Watson et
       *  al. 2012), Manera (Manera et al. 2010), Bhattacharya
       *  (Bhattacharya et al. 2011), Courtin (Courtin et al. 2010),
       *  Peacock (by Peacock at al. 2007)
       *
       *  @param method_SS method used to compute the power spectrum
       *  and &sigma;(mass); valid method_SS are: CAMB
       *  [http://camb.info/], classgal_v1 [http://class-code.net/],
       *  EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       *
       *  @param selection_function_file input file where the selection
       *  function is stored
       *
       *  @param column the columns to be read
       *
       *  @param output_root output_root of the parameter file used to
       *  compute the power spectrum and &sigma;(mass); it can be any
       *  name
       *
       *  @param Delta \f$\Delta\f$: the overdensity, defined as the
       *  mean interior density relative to the background
       *
       *  @param isDelta_vir \f$\rightarrow\f$ \f$\Delta\f$ is the
       *  virial overdensity
       *
       *  @param interpType method to interpolate the power spectrum
       *
       *  @param k_max maximum wave vector module up to which the power
       *  spectrum is computed
       *
       *  @param input_file either the parameter file or the power
       *  spectrum file; if a parameter file is provided,
       *  i.e. input_file!=NULL and is_parameter_file=true, it will be
       *  used to compute the power spectrum; if a power spectrum file
       *  is provided, i.e. input_file!=NULL and
       *  is_parameter_file=false, then the provided power spectrum
       *  will be used directly; in both cases &sigma;<SUP>2</SUP>(M)
       *  is computed by integrating the computed/provided power
       *  spectrum ignoring the cosmological parameters of the object
       *
       *  @param is_parameter_file true \f$\rightarrow\f$ the input_file
       *  is a parameter file, used to compute the power spectrum with
       *  the method specified by method_Pk; false \f$\rightarrow\f$
       *  the input_file is a file containing the power spectrum
       *
       *  @return the redshift distribution of dark matter haloes
       *
       *  @warning the input parameter \f$\Delta\f$ is the background
       *  overdensity, not the critical overdensity
       *  \f$\Delta_{crit}\f$; the function
       *  cbl::Cosmology::Delta_vir can be used to convert
       *  \f$\Delta_{crit}\f$ into \f$\Delta\f$
       *
       *  @warning the mass function by Manera et al. (2010) has been
       *  tested only for z=0 and z=0.5; the mass function by Despali
       *  et al. (2016) is currently implemented only for virial
       *  masses and at \f$z<1.25\f$
       */
      std::vector<double> redshift_distribution_haloes_selection_function (const std::vector<double> redshift, const double Area_degrees, const double Mass_min, const double Mass_max, const std::string model_MF, const std::string method_SS, const std::string selection_function_file, const std::vector<int> column={}, const std::string output_root="test", const double Delta=200, const bool isDelta_vir=false, const std::string interpType="Linear", const double k_max=100., const std::string input_file=par::defaultString, const bool is_parameter_file=true);

      /**
       *  @brief the mean redshift of a dark matter haloe sample,
       *  given a selection function
       *
       *  @author Alfonso Veropalumbo
       *
       *  @author alfonso.veropalumbo@unibo.it
       *
       *  @param z_min minimum redshift
       *
       *  @param z_max maximum redshift
       *
       *  @param Mass_min minimum halo mass
       *
       *  @param Mass_max maximum halo mass
       *
       *  @param model_MF author(s) who proposed the mass function;
       *  valid authors are: PS (Press&Schechter), ST (Sheth&Tormen),
       *  Jenkins (Jenkins et al. 2001), Warren (Warren et al. 2006),
       *  Reed, (Reed et al. 2007), Pan (Pan 2007), ShenH (halo MF by
       *  Shen et al. 2006), ShenF (filaments MF by Shen et al. 2006),
       *  ShenS (sheets MF by Shen et al. 2006), Tinker (Tinker et
       *  al. 2008), Crocce (Crocce et al. 2010), Angulo_FOF (FOF MF
       *  by Angulo et al. 2012), Angulo_Sub (SUBFIND MF by Angulo et
       *  al. 2012), Watson_FOF(FOF MF by Watson et al. 2012),
       *  Watson_SOH (MF for Spherical Overdensity Haloes by Watson et
       *  al. 2012), Manera (Manera et al. 2010), Bhattacharya
       *  (Bhattacharya et al. 2011), Courtin (Courtin et al. 2010),
       *  Peacock (by Peacock at al. 2007)
       *
       *  @param method_SS method used to compute the power spectrum
       *  and &sigma;(mass); valid method_SS are: CAMB
       *  [http://camb.info/], classgal_v1 [http://class-code.net/],
       *  EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       *
       *  @param selection_function_file input file where the selection
       *  function is stored
       *
       *  @param column the columns to be read
       *
       *  @param output_root output_root of the parameter file used to
       *  compute the power spectrum and &sigma;(mass); it can be any
       *  name
       *
       *  @param Delta \f$\Delta\f$: the overdensity, defined as the
       *  mean interior density relative to the background
       *
       *  @param isDelta_vir \f$\rightarrow\f$ \f$\Delta\f$ is the
       *  virial overdensity
       *
       *  @param interpType method to interpolate the power spectrum
       *
       *  @param k_max maximum wave vector module up to which the power
       *  spectrum is computed
       *
       *  @param input_file either the parameter file or the power
       *  spectrum file; if a parameter file is provided,
       *  i.e. input_file!=NULL and is_parameter_file=true, it will be
       *  used to compute the power spectrum; if a power spectrum file
       *  is provided, i.e. input_file!=NULL and
       *  is_parameter_file=false, then the provided power spectrum
       *  will be used directly; in both cases &sigma;<SUP>2</SUP>(M)
       *  is computed by integrating the computed/provided power
       *  spectrum ignoring the cosmological parameters of the object
       *
       *  @param is_parameter_file true \f$\rightarrow\f$ the input_file
       *  is a parameter file, used to compute the power spectrum with
       *  the method specified by method_Pk; false \f$\rightarrow\f$
       *  the input_file is a file containing the power spectrum
       *
       *  @return the mean redshift of dark matter haloes
       *
       *  @warning the input parameter \f$\Delta\f$ is the background
       *  overdensity, not the critical overdensity
       *  \f$\Delta_{crit}\f$; the function
       *  cbl::Cosmology::Delta_vir can be used to convert
       *  \f$\Delta_{crit}\f$ into \f$\Delta\f$
       *
       *  @warning the mass function by Manera et al. (2010) has been
       *  tested only for z=0 and z=0.5; the mass function by Despali
       *  et al. (2016) is currently implemented only for virial
       *  masses and at \f$z<1.25\f$
       */
      double mean_redshift_haloes_selection_function (const double z_min, const double z_max, const double Mass_min, const double Mass_max, const std::string model_MF, const std::string method_SS, const std::string selection_function_file, const std::vector<int> column={}, const std::string output_root="test", const double Delta=200, const bool isDelta_vir=false, const std::string interpType="Linear", const double k_max=100., const std::string input_file=par::defaultString, const bool is_parameter_file=true);
  
      /**
       *  @brief minimum halo mass, given the number of haloes in a
       *  given region of sky
       *
       *  @author Alfonso Veropalumbo, Jacopo Neri (and Federico
       *  Marulli)
       *
       *  @author alfonso.veropalumbo@unibo.it, jacopo.neri6@gmail.com
       *  (and federico.marulli3@unibo.it)
       *
       *  @param n_halo number density of dark matter haloes
       *
       *  @param Area sky area
       *
       *  @param angle_rad 0 \f$\rightarrow\f$ &Omega; in square
       *  degrees; 1 \f$\rightarrow\f$ &Omega; in steradians
       *
       *  @param z_min minimum redshift
       * 
       *  @param z_max maximum redshift
       *
       *  @param Mmax maximum mass
       *  
       *  @param lgM1_guess logarithm of the minimum mass used by the
       *  root finder
       *
       *  @param lgM2_guess logarithm of the maximum mass used by the
       *  root finder
       *  
       *  @param model_MF author(s) who proposed the mass function;
       *  valid authors are: PS (Press&Schechter), ST (Sheth&Tormen),
       *  Jenkins (Jenkins et al. 2001), Warren (Warren et al. 2006),
       *  Reed, (Reed et al. 2007), Pan (Pan 2007), ShenH (halo MF by
       *  Shen et al. 2006), ShenF (filaments MF by Shen et al. 2006),
       *  ShenS (sheets MF by Shen et al. 2006), Tinker (Tinker et
       *  al. 2008), Crocce (Crocce et al. 2010), Angulo_FOF (FOF MF by
       *  Angulo et al. 2012), Angulo_Sub (SUBFIND MF by Angulo et
       *  al. 2012), Watson_FOF(FOF MF by Watson et al. 2012),
       *  Watson_SOH (MF for Spherical Overdensity Haloes by Watson 
       *  et al. 2012), Manera (Manera et al. 2010), Bhattacharya 
       *  (Bhattacharya et al. 2011), Courtin (Courtin et al. 2010),
       *  Peacock (by Peacock at al. 2007)
       *
       *  @param method_SS method used to compute the power spectrum
       *  and &sigma;(mass); valid method_SS are: CAMB
       *  [http://camb.info/], classgal_v1 [http://class-code.net/],
       *  EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       *
       *  @param output_root output_root of the parameter file used to
       *  compute the power spectrum and &sigma;(mass); it can be any
       *  name
       *
       *  @param Delta \f$\Delta\f$: the overdensity, defined as the
       *  mean interior density relative to the background
       *
       *  @param interpType method to interpolate the power spectrum
       *
       *  @param k_max maximum wave vector module up to which the power
       *  spectrum is computed
       *
       *  @param input_file either the parameter file or the power
       *  spectrum file; if a parameter file is provided,
       *  i.e. input_file!=NULL and is_parameter_file=true, it will be
       *  used to compute the power spectrum; if a power spectrum file
       *  is provided, i.e. input_file!=NULL and
       *  is_parameter_file=false, then the provided power spectrum
       *  will be used directly; in both cases &sigma;<SUP>2</SUP>(M)
       *  is computed by integrating the computed/provided power
       *  spectrum ignoring the cosmological parameters of the object
       *
       *  @param is_parameter_file true \f$\rightarrow\f$ the input_file
       *  is a parameter file, used to compute the power spectrum with
       *  the method specified by method_Pk; false \f$\rightarrow\f$
       *  the input_file is a file containing the power spectrum
       *
       *  @return minimum halo mass
       *
       *  @warning the input parameter \f$\Delta\f$ is the background
       *  overdensity, not the critical overdensity
       *  \f$\Delta_{crit}\f$; the function
       *  cbl::Cosmology::Delta_vir can be used to convert
       *  \f$\Delta_{crit}\f$ into \f$\Delta\f$
       *
       *  @warning the mass function by Manera et al. (2010) has been
       *  tested only for z=0 and z=0.5; the mass function by Despali
       *  et al. (2016) is currently implemented only for virial
       *  masses and at \f$z<1.25\f$
       */
      double MhaloMin (const int n_halo, const double Area, const bool angle_rad, const double z_min, const double z_max, const double Mmax, const double lgM1_guess, const double lgM2_guess, const std::string model_MF, const std::string method_SS, const std::string output_root="test", const double Delta=200, const std::string interpType="Linear", const double k_max=100., const std::string input_file=par::defaultString, const bool is_parameter_file=true) const;

      /**
       *  @brief convert a cluster mass estimated in a different
       *  cosmology
       *
       *  this function converts a cluster mass estimated assuming a
       *  different cosmology, following Eq.C4 of Sereno & Ettori 2015
       *  (https://arxiv.org/abs/1407.7868)
       *
       *  \f[ M^{(2)} = M^{(1)} \left( \frac{D_{ds}^{(1)}}{D_s^{(1)}}
       *  \right)^{3/2} H^{(1)} \left(
       *  \frac{D_{ds}^{(2)}}{D_s^{(2)}}\right)^{-3/2} H^{(2) -1} \f]
       *
       *  where and \f$D_s\f$ and \f$D_{ds}\f$ are the source and the
       *  lens-source angular diameter distances, respectively
       *  
       *  @param mass the cluster mass to be converted (since
       *  estimated in a different cosmology), \f$M^{(1)}\f$
       *
       *  @param cosmology the cosmology assumed to measure the
       *  cluster mass
       * 
       *  @param redshift redshift of the cluster
       *
       *  @param redshift_source redshift of the source, if the
       *  cluster mass is estimated from weak lensing, -1 otherwise
       *
       *  @return the cluster mass converted in this cosmology,
       *  \f$M^{(2)}\f$
       */
      double converted_mass (const double mass, const cosmology::Cosmology cosmology, const double redshift, const double redshift_source=-1.) const;
      
      ///@}


      /**
       *  @name Functions to estimate the cosmic mass accretion history
       */
      ///@{

      /**
       *  @brief differential distribution 
       *
       *  this function provides the differential rescaled and
       *  generalized formation redshift distribution
       *
       *  @author Carlo Giocoli
       *  @author cgiocoli@gmail.com
       *
       *  @param ww rescaled variable w as in Lacey and Coles 1993
       *  @param ff assembled fraction
       *  @param author valid authors are: NS (Nusser and Sheth), GTS
       *  (Giocoli et al. 2012)
       *
       *  @return p(w): differential distribution 
       */
      double pw (const double ww, const double ff, const std::string author) const; 

      /**
       *  @brief formation probability
       *
       *  this function provides the probability that a halo of a given
       *  mass m0 at redshift z0 make a mass fraction f at redshift z
       *
       *  @author Carlo Giocoli
       *  @author cgiocoli@gmail.com
       *  @param m0 halo mass
       *  @param z0 redshift when the halo has a mass m0
       *  @param frac mass fraction
       *  @param redshift the redshift
       *
       *  @param model_model valid authors are: NS (Nusser and Sheth),
       *  GTS (Giocoli et al. 2012)
       *
       *  @param method_SS method used to compute the power spectrum and
       *  &sigma;(mass); valid method_SS are: CAMB [http://camb.info/],
       *  classgal_v1 [http://class-code.net/], EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       *     
       *  @param output_root output_root of the parameter file used to
       *  compute the power spectrum and &sigma;(mass); it can be any
       *  name
       *
       *  @return p(z): formation probability
       */
      double pz (const double m0, const double z0, const double frac, const double redshift, const std::string model_model, const std::string method_SS, const std::string output_root="test") const; 

      /**
       *  @brief cumulative distribution 
       *
       *  this function provides the cumulative rescaled and generalized
       *  formation redshift distribution
       *
       *  @author Carlo Giocoli
       *  @author cgiocoli@gmail.com
       *  @param ww rescaled variable w as in Lacey and Coles 1993
       *  @param ff assembled fraction
       *  @param author valid authors are: NS (Nusser and Sheth), GTS
       *  (Giocoli et al. 2012)
       *  @return P(w): cumulative distribution
       */
      double cumPw (const double ww, const double ff, const std::string author) const; 
    
      /**
       *  @brief median formation w     
       *  @author Carlo Giocoli
       *  @author cgiocoli@gmail.com
       *  @param [in] ff assembled fraction
       *  @param [in] model_model valid authors are: NS (Nusser and Sheth), GTS
       *  (Giocoli et al. 2012)
       *  @param [out] wf vector of w(f)
       *  @return none
       */
      void medianwf (const double ff, const std::string model_model, std::vector<double> &wf) const; 

      /**
       *  @brief median formation z 
       *  @author Carlo Giocoli
       *  @author cgiocoli@gmail.com
       *  @param [in] ff assembled fraction
       *  @param [in] mass halo mass
       *  @param [in] z0 redshift when the halo has a mass mass
       *
       *  @param [in] model_model valid authors are: NS (Nusser and
       *  Sheth), GTS (Giocoli et al. 2012)
       *
       *  @param method_SS method used to compute the power spectrum and
       *  &sigma;(mass); valid method_SS are: CAMB [http://camb.info/],
       *  classgal_v1 [http://class-code.net/], EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       *
       *  @param [out] zf vector of z(f)
       *      
       *  @param output_root output_root of the parameter file used to
       *  compute the power spectrum and &sigma;(mass); it can be any
       *  name
       *
       *  @return none
       */
      void medianzf (const double ff, const double mass, const double z0, const std::string model_model, const std::string method_SS, std::vector<double> &zf, const std::string output_root="test") const; 
  
      /**
       *  @brief rescaled variable w as in Lacey and Coles 1993
       *
       *  this functions provides the conditional variable
       *  w=[&delta;<SUB>c</SUB>(zf) - &delta;<SUB>c</SUB>(z)] / &radic;
       *  [s(fm)-s(m)] where&delta;<SUB>c</SUB>(z) =
       *  &delta;<SUB>c0</SUB>(z)/D+(z)
       *
       *  @author Carlo Giocoli
       *  @author cgiocoli@gmail.com
       *  @param mm halo mass
       *  @param redshift the redshift 
       *  @param ff assembled fraction
       *  @param zf redshift at which the mass is accreted
       *
       *  @param method_SS method used to compute the power spectrum and
       *  &sigma;(mass); valid method_SS are: CAMB [http://camb.info/],
       *  classgal_v1 [http://class-code.net/], EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       *     
       *  @param output_root output_root of the parameter file used to
       *  compute the power spectrum and &sigma;(mass); it can be any
       *  name
       *
       *  @return the conditional variable w
       */
      double wf (const double mm, const double redshift, const double ff, const double zf, const std::string method_SS, const std::string output_root="test") const; 

      /**
       *  @brief the unevolved mass function
       *  @author Carlo Giocoli
       *  @author cgiocoli@gmail.com
       *  @param mass_accr mass accreted
       *  @return the unevolved mass function
       */
      double unevolved_mass_function (const double mass_accr) const; 

      ///@}


      /**
       *  @name Functions to estimate the power spectrum and related quantities
       */
      ///@{

      /**
       *  @brief amplitude of the curvature perturbations
       *  
       *  this function provides an approximate value of A<SUB>s</SUB>,
       *  for a fiven value of &sigma;<SUB>8</SUB>, from Vikhlinin et
       *  al. 2009, Eq.3; it is valid only without massive neutrinos! 
       *  (see also Hu & Jain 2004)
       *
       *  @param sigma8 &sigma;<SUB>8</SUB> the power spectrum normalisation
       *  @return A<SUB>s</SUB>
       */
      double As (const double sigma8) const; 

      /**
       *  @brief &sigma;<SUB>8</SUB>
       *  
       *  this function provides an approximate value of
       *  &sigma;<SUB>8</SUB>, at a given redshift, from Aubourg et
       *  al. 2015, eq.(32)
       *  
       *  @warning it is valid only without massive neutrinos! To
       *  account for non-zero neutrino mass authors multiplies by an
       *  extra factor of 0.995
       *
       *  @param redshift the redshift
       *  @return &sigma;<SUB>8</SUB>
       */
      double sigma8_interpolated (const double redshift) const; 

      /**
       *  @brief return the path to the power spectrum output
       *
       *  @param code method used to compute the power spectrum
       * 
       *  @param NL 0 \f$\rightarrow\f$ linear power spectrum; 1
       *  \f$\rightarrow\f$ non-linear power spectrum
       * 
       *  @param redshift the redshift
       *
       *  @param run true \f$\rightarrow\f$ write or read the table
       *  where the dark matter power spectrum is stored
       *
       *  @param output_root the output_root parameter of the
       *  parameter file used to compute the power spectrum; it can be
       *  any name
       *
       *  @param k_max the maximum wave vector module up to which the
       *  power spectrum is computed
       *
       *  @param file_par name of the parameter file; if a parameter
       *  file is provided (i.e. file_par!=NULL), it will use be used,
       *  ignoring the cosmological parameters of the object
       *
       *  @return the path to the power spectrum output
       */
      std::string Pk_output_file (const std::string code, const bool NL, const double redshift, const bool run=0, const std::string output_root="test", const double k_max=100., const std::string file_par=par::defaultString);

      /**
       *  @brief run CAMB [http://camb.info/]
       *  
       *  this function runs CAMB [http://camb.info/], after editing the parameter file
       *  appropriately (if file_par=NULL)
       *
       *  @param NL 0 \f$\rightarrow\f$ linear power spectrum; 1 \f$\rightarrow\f$ non-linear power
       *  spectrum
       *
       *  @param redshift the redshift
       *
       *  @param output_root the output_root parameter of the
       *  parameter file used to compute the power spectrum; it can be
       *  any name. If NULL, the 
       *  output will be deleted after running CAMB
       *
       *  @param output_dir std::string containing the output directory
       *
       *  @param k_max maximum wave vector module up to which the power
       *  spectrum is computed
       *
       *  @param file_par name of the parameter file; if a parameter
       *  file is provided (i.e. file_par!=NULL), it will use be used,
       *  ignoring the cosmological parameters of the object
       *
       *  @return none
       */
      void run_CAMB (const bool NL, const double redshift, const std::string output_root=par::defaultString, const std::string output_dir=par::defaultString, const double k_max=100., const std::string file_par=par::defaultString) const; 

      /**
       *  @brief run CAMB [http://camb.info/] and read the matter
       *  power spectrum
       *  
       *  this function runs CAMB [http://camb.info/], after editing
       *  the parameter file appropriately (if file_par=NULL) and
       *  store the matter power spectrum in two vectors
       *
       *  @param [out] lgkk vector of log(k)
       *
       *  @param [out] lgPk vector of log(P(k))
       *
       *  @param NL 0 \f$\rightarrow\f$ linear power spectrum; 1 \f$\rightarrow\f$ non-linear power
       *  spectrum
       *
       *  @param redshift the redshift
       *
       *  @param output_root the output_root parameter of the
       *  parameter file used to compute the power spectrum; it can be
       *  any name. If NULL, the 
       *  output will be deleted after running CAMB
       *
       *  @param output_dir std::string containing the output directory
       *
       *  @param k_max maximum wave vector module up to which the power
       *  spectrum is computed
       *
       *  @param file_par name of the parameter file; if a parameter
       *  file is provided (i.e. file_par!=NULL), it will use be used,
       *  ignoring the cosmological parameters of the object
       *
       *  @return none
       */
      void run_CAMB (std::vector<double> &lgkk, std::vector<double> &lgPk, const bool NL, const double redshift, const std::string output_root=par::defaultString, const std::string output_dir=par::defaultString, const double k_max=100., const std::string file_par=par::defaultString) const;

      /**
       *  @brief write or read the table where the dark matter power
       *  spectrum is stored
       *
       *  @param [in] code method used to compute the power spectrum;
       *  valid codes are: CAMB [http://camb.info/], classgal_v1
       *  [http://class-code.net/], MPTbreeze-v1
       *  [http://arxiv.org/abs/1207.1465], EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       *
       *  @param [in] NL 0 \f$\rightarrow\f$ linear power spectrum; 1 \f$\rightarrow\f$ non-linear power spectrum
       *
       *  @param [out] lgkk vector of log(k)
       *
       *  @param [out] lgPk vector of log(P(k))
       *
       *  @param [in] redshift redshift
       *
       *  @param [in] output_root output_root of the parameter file used
       *  to compute the power spectrum; it can be any name
       *
       *  @param [in] k_max maximum wave vector module up to which the
       *  power spectrum is computed
       *
       *  @param [in] file_par name of the parameter file; if a
       *  parameter file is provided (i.e. file_par!=NULL), it will be
       *  used, ignoring the cosmological parameters of the object
       *
       *  @return none
       */
      void Table_PkCodes (const std::string code, const bool NL, std::vector<double> &lgkk, std::vector<double> &lgPk, const double redshift, const std::string output_root="test", const double k_max=100., std::string file_par=par::defaultString) const;

      /**
       *  @brief write or read the table where the dark matter two-point
       *  correlation function is stored
       *
       *  @param [in] code method used to compute the power spectrum;
       *  valid codes are: CAMB [http://camb.info/], classgal_v1
       *  [http://class-code.net/], MPTbreeze-v1
       *  [http://arxiv.org/abs/1207.1465], EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       *
       *  @param [in] NL 0 \f$\rightarrow\f$ linear power spectrum; 1 \f$\rightarrow\f$ non-linear power spectrum
       *
       *  @param [out] rr vector of comoving separations
       *
       *  @param [out] xi vector of the binned values of &xi;(r)
       *
       *  @param [in] redshift redshift
       *
       *  @param [in] output_root output_root of the parameter file used
       *  to compute the power spectrum; it can be any name
       *
       *  @param [in] k_max maximum wave vector module up to which the
       *  power spectrum is computed
       *
       *  @param [in] file_par name of the parameter file; if a
       *  parameter file is provided (i.e. file_par!=NULL), it will be
       *  used, ignoring the cosmological parameters of the object
       *
       *  @return none
       */
      void Table_XiCodes (const std::string code, const bool NL, std::vector<double> &rr, std::vector<double> &xi, const double redshift, const std::string output_root, const double k_max, std::string file_par) const;

      /**
       *  @brief normalisation of the power spectrum
       *
       *  this function sets the value of the private member m_Pk0_*,
       *  i.e. the normalisation of the power spectrum
       *
       *  @param method_Pk method used to compute the power spectrum;
       *  valid choices for method_Pk are: CAMB [http://camb.info/],
       *  classgal_v1 [http://class-code.net/], MPTbreeze-v1
       *  [http://arxiv.org/abs/1207.1465], EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       *
       *  @param redshift the redshift
       *
       *  @param output_root the output_root parameter of the
       *  parameter file used to compute the power spectrum; it can be
       *  any name
       *
       *  @param k_min minimum wave vector module up to which the power
       *  spectrum is computed
       *
       *  @param k_max maximum wave vector module up to which the power
       *  spectrum is computed
       *
       *  @param prec accuracy of the integration 
       *
       *  @param file_par name of the parameter file; if a parameter
       *  file is provided (i.e. file_par!=NULL), it will be used,
       *  ignoring the cosmological parameters of the object
       *
       *  @return none
       */
      void Pk_0 (const std::string method_Pk, const double redshift, const std::string output_root="test", const double k_min=0., const double k_max=100., const double prec=1.e-2, const std::string file_par=par::defaultString); 

      /**
       *  @brief normalised power spectrum
       *
       *  this function provides the power spectrum P(k); it can use
       *  either CAMB, CLASS, MPTbreeze or the analytic
       *  approximation by Eisenstein & Hu
       *
       *  @param kk the wave vector module
       *
       *  @param method_Pk method used to compute the power spectrum;
       *  valid choices for method_Pk are: CAMB [http://camb.info/],
       *  classgal_v1 [http://class-code.net/], MPTbreeze-v1
       *  [http://arxiv.org/abs/1207.1465], EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       *
       *  @param NL 0 \f$\rightarrow\f$ linear power spectrum; 1 \f$\rightarrow\f$
       *  non-linear power spectrum
       *
       *  @param redshift the redshift
       *
       *  @param output_root the output_root parameter of the
       *  parameter file used to compute the power spectrum; it can be
       *  any name
       *
       *  @param norm 0 \f$\rightarrow\f$ don't normalise the power
       *  spectrum; 1 \f$\rightarrow\f$ normalise the power spectrum;
       *  -1 \f$\rightarrow\f$ normalise only if sigma8 is set
       *
       *  @param k_min minimum wave vector module up to which the power
       *  spectrum is computed
       *
       *  @param k_max maximum wave vector module up to which the power
       *  spectrum is computed
       *
       *  @param prec accuracy of the integration 
       *
       *  @param file_par name of the parameter file; if a parameter
       *  file is provided (i.e. file_par!=NULL), it will be used,
       *  ignoring the cosmological parameters of the object
       *
       *  @param unit1 true \f$\rightarrow\f$ force cosmological units
       *
       *  @return P(k)
       */
      double Pk (const double kk, const std::string method_Pk, const bool NL, const double redshift, const std::string output_root="test", const int norm=-1, const double k_min=0., const double k_max=100., const double prec=1.e-2, const std::string file_par=par::defaultString, const bool unit1=false); 

      /**
       *  @brief normalised power spectrum
       *
       *  this function provides the power spectrum P(k); it can use
       *  either CAMB, CLASS, MPTbreeze or the analytic
       *  approximation by Eisenstein & Hu
       *
       *  @param kk the wave vector module
       *
       *  @param method_Pk method used to compute the power spectrum;
       *  valid choices for method_Pk are: CAMB [http://camb.info/],
       *  classgal_v1 [http://class-code.net/], MPTbreeze-v1
       *  [http://arxiv.org/abs/1207.1465], EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       *
       *  @param NL 0 \f$\rightarrow\f$ linear power spectrum; 1
       *  \f$\rightarrow\f$ non-linear power spectrum
       *
       *  @param redshift the redshift
       *
       *  @param output_dir the output_dir directory
       *  where the output of external codes are written
       *  
       *  @param output_root the output_root parameter of the
       *  parameter file used to compute the power spectrum; it can be
       *  any name
       *
       *  @param norm 0 \f$\rightarrow\f$ don't normalise the power
       *  spectrum; 1 \f$\rightarrow\f$ normalise the power spectrum;
       *  -1 \f$\rightarrow\f$ normalise only if sigma8 is set
       *
       *  @param k_min minimum wave vector module up to which the power
       *  spectrum is computed
       *
       *  @param k_max maximum wave vector module up to which the power
       *  spectrum is computed
       *
       *  @param prec accuracy of the integration 
       *
       *  @param file_par name of the parameter file; if a parameter
       *  file is provided (i.e. file_par!=NULL), it will be used,
       *  ignoring the cosmological parameters of the object
       *
       *  @param unit1 true \f$\rightarrow\f$ force cosmological units
       *
       *  @return P(k)
       */
      std::vector<double> Pk (const std::vector<double> kk, const std::string method_Pk, const bool NL, const double redshift, const std::string output_dir, const std::string output_root="test", const int norm=-1, const double k_min=0., const double k_max=100., const double prec=1.e-2, const std::string file_par=par::defaultString, const bool unit1=false); 

      /**
       *  @brief return power spectrum first three multipoles using linear kaiser model
       *
       *  this function provides the power spectrum P(k); it can use
       *  either CAMB, CLASS, MPTbreeze or the analytic
       *  approximation by Eisenstein & Hu
       *
       *  @param Pk0 the monopole of the power spectrum
       *
       *  @param Pk2 the quadrupole of the power spectrum
       *
       *  @param Pk4 the hexadecapole of the power spectrum
       *
       *  @param kk the scale at which compute the power spectrum
       *
       *  @param method_Pk method used to compute the power spectrum;
       *  valid choices for method_Pk are: CAMB [http://camb.info/],
       *  classgal_v1 [http://class-code.net/], MPTbreeze-v1
       *  [http://arxiv.org/abs/1207.1465], EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       *
       *  @param NL 0 \f$\rightarrow\f$ linear power spectrum; 1
       *  \f$\rightarrow\f$ non-linear power spectrum
       *
       *  @param redshift the redshift
       * 
       *  @param bias the linear bias
       *
       *  @param sigma_NL the BAO damping parameter
       *
       *  @param output_root the output_root parameter of the
       *  parameter file used to compute the power spectrum; it can be
       *  any name
       *
       *  @param norm 0 \f$\rightarrow\f$ don't normalise the power
       *  spectrum; 1 \f$\rightarrow\f$ normalise the power spectrum;
       *  -1 \f$\rightarrow\f$ normalise only if sigma8 is set
       *
       *  @param k_min minimum wave vector module up to which the power
       *  spectrum is computed
       *
       *  @param k_max maximum wave vector module up to which the power
       *  spectrum is computed
       *
       *  @param prec accuracy of the integration 
       *
       *  @param file_par name of the parameter file; if a parameter
       *  file is provided (i.e. file_par!=NULL), it will be used,
       *  ignoring the cosmological parameters of the object
       *
       *  @return none
       */
      void Pk_Kaiser_multipoles (std::vector<double> &Pk0, std::vector<double> &Pk2, std::vector<double> &Pk4, const std::vector<double> kk, const std::string method_Pk, const bool NL, const double redshift, const double bias, const double sigma_NL = 0., const std::string output_root="test", const int norm=-1, const double k_min=0., const double k_max=100., const double prec=1.e-2, const std::string file_par=par::defaultString); 

      /**
       *  @brief the dark matter power spectrum, de-wiggled (see
       *  e.g. Anderson et al 2014)
       *
       *  this function provides the De-Wiggled dark matter power
       *  spectrum
       *
       *  @author Alfonso Veropalumbo
       *  @author alfonso.veropalumbo@unibo.it
       *  
       *  @param kk the wave vector module
       *
       *  @param redshift the redshift
       *
       *  @param sigma_NL the non linear BAO damping
       *
       *  @param output_root output_root of the parameter file used to
       *  compute the power spectrum and &sigma;(mass); it can be any
       *  name
       *
       *  @param norm 0 \f$\rightarrow\f$ don't normalise the power
       *  spectrum; 1 \f$\rightarrow\f$ normalise the power spectrum;
       *  -1 \f$\rightarrow\f$ normalise only if sigma8 is set
       *
       *  @param k_min minimum wave vector module up to which the power
       *  spectrum is computed
       *
       *  @param k_max maximum wave vector module up to which the power
       *  spectrum is computed
       *
       *  @param aa parameter \e a of Eq. 24 of Anderson et al. 2012
       *
       *  @param prec accuracy of the integration 
       *
       *  @return P;<SUB>DW</SUB>(k): the De-Wiggled power
       *  spectrum of dark matter
       */
      double Pk_DeWiggle (const double kk, const double redshift, const double sigma_NL, const std::string output_root = "test", const bool norm=1, const double k_min=0., const double k_max=100., const double aa=1., const double prec=1.e-2);
 
      /**
       *  @brief the mass variance, \f$\sigma^2(R)\f$
       *
       *  this function computes the variance of the linear density
       *  field:
       *
       *  \f[ \sigma^2(R)=\frac{1}{2\pi^2}\int_0^\infty {\rm d}k\, k^2
       *  P_{lin}(k, z) W^2(k, R)\f]
       *
       *  where \f$W(x)=(3/x)^3(\sin x-x\cos x)\f$
       *
       *  @param radius the radius, \f$R\f$
       *
       *  @param method_Pk method used to compute the power spectrum;
       *  valid choices for method_Pk are: CAMB [http://camb.info/],
       *  classgal_v1 [http://class-code.net/], MPTbreeze-v1
       *  [http://arxiv.org/abs/1207.1465], EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       *
       *  @param redshift the redshift
       *
       *  @param output_root the output_root parameter of the
       *  parameter file used to compute the power spectrum; it can be
       *  any name
       *
       *  @param interpType method to interpolate the power spectrum
       *
       *  @param k_max maximum wave vector module up to which the power
       *  spectrum is computed
       *
       *  @param input_file either the parameter file or the power
       *  spectrum file; if a parameter file is provided,
       *  i.e. input_file!=NULL and is_parameter_file=true, it will be
       *  used to compute the power spectrum; if a power spectrum file
       *  is provided, i.e. input_file!=NULL and
       *  is_parameter_file=false, then the provided power spectrum
       *  will be used directly; in both cases &sigma;<SUP>2</SUP>(M)
       *  is computed by integrating the computed/provided power
       *  spectrum ignoring the cosmological parameters of the object
       *
       *  @param is_parameter_file true \f$\rightarrow\f$ the input_file
       *  is a parameter file, used to compute the power spectrum with
       *  the method specified by method_Pk; false \f$\rightarrow\f$
       *  the input_file is a file containing the power spectrum
       *
       *  @param unit1 true \f$\rightarrow\f$ force cosmological units
       *
       *  @return \f$\sigma^2(R)\f$
       */
      double sigma2R (const double radius, const std::string method_Pk, const double redshift, const std::string output_root="test", const std::string interpType="Linear", const double k_max=100., const std::string input_file=par::defaultString, const bool is_parameter_file=true, const bool unit1=false) const; 

      /**
       *  @brief the mass variance, \f$\sigma^2(M)\f$
       *
       *  this function computes the variance of the
       *  linear density field:
       *
       *  \f[ \sigma^2(M) = \frac{1}{2\pi^2}\int_0^\infty {\rm d}k\,
       *  k^2 P_{lin}(k, z) W^2(k, R)\f]
       *
       *  where \f$W(x)=(3/x)^3(\sin x-x\cos x)\f$ and
       *  \f$R=(3M/4\pi\rho_m)^{1/3}\f$
       *
       *  @param mass the mass 
       *
       *  @param method_Pk the method used to compute the power
       *  spectrum; valid choices for method_Pk are: CAMB
       *  [http://camb.info/], classgal_v1 [http://class-code.net/],
       *  MPTbreeze-v1 [http://arxiv.org/abs/1207.1465], EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       *
       *  @param redshift the redshift
       *
       *  @param output_root the output_root parameter of the
       *  parameter file used to compute the power spectrum; it can be
       *  any name
       *
       *  @param interpType method to interpolate the power spectrum
       *
       *  @param k_max maximum wave vector module up to which the power
       *  spectrum is computed
       *
       *  @param input_file either the parameter file or the power
       *  spectrum file; if a parameter file is provided,
       *  i.e. input_file!=NULL and is_parameter_file=true, it will be
       *  used to compute the power spectrum; if a power spectrum file
       *  is provided, i.e. input_file!=NULL and
       *  is_parameter_file=false, then the provided power spectrum
       *  will be used directly; in both cases &sigma;<SUP>2</SUP>(M)
       *  is computed by integrating the computed/provided power
       *  spectrum ignoring the cosmological parameters of the object
       *
       *  @param is_parameter_file true \f$\rightarrow\f$ the input_file
       *  is a parameter file, used to compute the power spectrum with
       *  the method specified by method_Pk; false \f$\rightarrow\f$
       *  the input_file is a file containing the power spectrum
       *
       *  @param unit1 true \f$\rightarrow\f$ force cosmological units
       *
       *  @return \f$\sigma^2(M)\f$
       */
      double sigma2M (const double mass, const std::string method_Pk, const double redshift, const std::string output_root="test", const std::string interpType="Linear", const double k_max=100., const std::string input_file=par::defaultString, const bool is_parameter_file=true, const bool unit1=false) const; 
 
      /**
       *  @brief the nth-order derivative of the mass variance,
       *  \f${\rm d}^n\sigma^2(R)/{\rm d}R^n\f$
       *
       *  @param nd the derivative order, \f$n\f$
       *
       *  @param radius the radius, \f$R\f$
       *
       *  @param method_Pk method used to compute the power spectrum;
       *  valid choices for method_Pk are: CAMB [http://camb.info/],
       *  classgal_v1 [http://class-code.net/], MPTbreeze-v1
       *  [http://arxiv.org/abs/1207.1465], EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       *
       *  @param redshift the redshift
       *
       *  @param output_root the output_root parameter of the
       *  parameter file used to compute the power spectrum; it can be
       *  any name
       *
       *  @param interpType method to interpolate the power spectrum
       *
       *  @param k_max maximum wave vector module up to which the power
       *  spectrum is computed
       *
       *  @param input_file either the parameter file or the power
       *  spectrum file; if a parameter file is provided,
       *  i.e. input_file!=NULL and is_parameter_file=true, it will be
       *  used to compute the power spectrum; if a power spectrum file
       *  is provided, i.e. input_file!=NULL and
       *  is_parameter_file=false, then the provided power spectrum
       *  will be used directly; in both cases &sigma;<SUP>2</SUP>(M)
       *  is computed by integrating the computed/provided power
       *  spectrum ignoring the cosmological parameters of the object
       *
       *  @param is_parameter_file true \f$\rightarrow\f$ the input_file
       *  is a parameter file, used to compute the power spectrum with
       *  the method specified by method_Pk; false \f$\rightarrow\f$
       *  the input_file is a file containing the power spectrum
       *
       *  @param unit1 true \f$\rightarrow\f$ force cosmological units
       *
       *  @return the nth-order derivative of the mass variance
       *
       *  @warning the current implementation computes the derivative
       *  using the simplest numerical approximation, with fixed
       *  incremental step; it is computationally efficient, but the
       *  accuracy might be lowt
       */
      double dnsigma2R (const int nd, const double radius, const std::string method_Pk, const double redshift, const std::string output_root="test", const std::string interpType="Linear", const double k_max=100., const std::string input_file=par::defaultString, const bool is_parameter_file=true, const bool unit1=false) const; 

      /**
       *  @brief the first derivative of the mass variance, \f${\rm
       *  d}^n\sigma^2(M)/{\rm d}M^n\f$
       * 
       *  @param nd the derivative order, \f$n\f$
       *
       *  @param mass the mass, \f$M\f$
       *
       *  @param method_Pk method used to compute the power spectrum;
       *  valid choices for method_Pk are: CAMB [http://camb.info/],
       *  classgal_v1 [http://class-code.net/], MPTbreeze-v1
       *  [http://arxiv.org/abs/1207.1465], EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       *
       *  @param redshift the redshift
       *
       *  @param output_root the output_root parameter of the
       *  parameter file used to compute the power spectrum; it can be
       *  any name
       *
       *  @param interpType method to interpolate the power spectrum
       *
       *  @param k_max maximum wave vector module up to which the power
       *  spectrum is computed
       *
       *  @param input_file either the parameter file or the power
       *  spectrum file; if a parameter file is provided,
       *  i.e. input_file!=NULL and is_parameter_file=true, it will be
       *  used to compute the power spectrum; if a power spectrum file
       *  is provided, i.e. input_file!=NULL and
       *  is_parameter_file=false, then the provided power spectrum
       *  will be used directly; in both cases &sigma;<SUP>2</SUP>(M)
       *  is computed by integrating the computed/provided power
       *  spectrum ignoring the cosmological parameters of the object
       *
       *  @param is_parameter_file true \f$\rightarrow\f$ the input_file
       *  is a parameter file, used to compute the power spectrum with
       *  the method specified by method_Pk; false \f$\rightarrow\f$
       *  the input_file is a file containing the power spectrum
       *
       *  @param unit1 true \f$\rightarrow\f$ force cosmological units
       *
       *  @return the first derivative of the mass variance
       *
       *  @warning the current implementation computes the derivative
       *  using the simplest numerical approximation, with fixed
       *  incremental step; it is computationally efficient, but the
       *  accuracy might be low
       */
      double dnsigma2M (const int nd, const double mass, const std::string method_Pk, const double redshift, const std::string output_root="test", const std::string interpType="Linear", const double k_max=100., const std::string input_file=par::defaultString, const bool is_parameter_file=true, const bool unit1=false) const; 

      ///@}


      /**
       *  @name Functions to estimate the halo density profile
       */
      ///@{
      
      /**
       *  @brief the concentration-mass relation
       * 
       *  this function computes the concentration of a dark matter
       *  halo of a given a mass, at a given redshift; the models
       *  implemented are the following:
       *
       *  - Duffy et al. 2008:
       *  \f[c(M_h, z) = A(M_h/M_{pivot})^B\,(1+z)^C\f]
       *
       *  @param Mass the halo mass
       *
       *  @param redshift the redshift
       *
       *  @param author the author(s) of the relation; available
       *  options are: "Duffy" (Duffy et al. 2008)
       *
       *  @param profile the density profile; available options are:
       *  "NFW" \f$\rightarrow\f$ Navarro-Frenk-White profile;
       *  "Einasto" \f$\rightarrow\f$ Einasto profile
       *
       *  @param halo_def the halo definition; available options are:
       *  "vir" \f$\rightarrow\f$ all matter withing the radius
       *  \f$r_{vir}\f$ for which the mean internal density is
       *  \f$\Delta\f$ times the critical density
       *  \f$\rho_{crit}=3H^2/8\pi G\f$; "200" \f$\rightarrow\f$ all
       *  matter withing the radius \f$r_{200}\f$ for which the mean
       *  internal density is 200 times the critical density; "mean"
       *  \f$\rightarrow\f$ all matter withing the radius
       *  \f$r_{200}\f$ for which the mean internal density is 200
       *  times the critical mean background density
       *
       *  @return the halo concentration
       *
       *  @warning the Duffy et al. concentrantion-mass relation
       *  refers to the 0<z<2 redshift range, obtained from their full
       *  samples (see Table 1 of Duffy et al. 2008); actually, the
       *  current implementation does not depend on cosmology, it is
       *  implemented here for possible future implementations (Zhao
       *  et al. 2009, Giocoli, Tormen & Sheth 2012)
       */
      double concentration (const double Mass, const double redshift, const std::string author="Duffy", const std::string profile="NFW", const std::string halo_def="vir") const;
          
      /**
       *  @brief compute the halo concentration
       *  @author Carlo Giocoli
       *  @author cgiocoli@gmail.com
       *  @param Vmax V<SUB>max</SUB>
       *  @param Rmax R<SUB>max</SUB>
       *  @return the halo concentration
       */
      double concentration2 (const double Vmax, const double Rmax) const;

      /**
       *  @brief virial halo concentration given \f$c_{200}\f$
       *  
       *  this function provides an approximate conversion to compute
       *  \f$c_{vir}\f$ from \f$c_{200}\f$ (from Coe 2010):
       *
       *  \f[c_{vir}\simeq a\,c_{200}+b\f]
       * 
       *  \f[a\simeq-1.119\log\Delta_c(z)+3.537\f]
       *
       *  \f[b\simeq-0.967\log\Delta_c(z)+2.181\f]
       *
       *  where \f$\Delta_c(z)\f$ is computed by
       *  cbl::cosmology::Cosmology::Delta_c
       *
       *  @param c200 \f$c_{200}\f$
       *
       *  @param redshift the redshift
       *
       *  @param author the author of the \f$\Delta_c(z)\f$
       *  equation (see cbl::cosmology::Cosmology::Delta_c)
       *
       *  @return \f$c_{vir}\f$
       */
      double c_vir (const double c200, const double redshift, const std::string author="BryanNorman") const;
      
      /**
       *  @brief the normalised halo density profile
       *
       *  this function computes the normalised density distribution
       *  of dark matter haloes; the Navarro-Frenk-White profile is
       *  the only one currently implemented (see e.g. eq. 74 of
       *  Cooray & Sheth 2002, eq. 68 of van den Bosch et
       *  al. 2012eq. 38 of Coe 2010):
       *
       *  \f[u_h(r, M_h, z) = \frac{\rho_h(r, M_h, z)}{M_h} =
       *  \frac{\rho_s}{(r/r_s)(1+r/r_s)^2}\f]
       *
       *  where
       *
       *  \f[\rho_s =
       *  \frac{\rho_{crit}\Delta_c}{3}\frac{c^3}{\ln(1+c)-c/(1+c)}\f]
       *
       *  the relation between the halo concentration,
       *  \f$c=c_{vir}=r_{vir}/r_s\f$, and halo mass, \f$M_h\f$, is
       *  computed by cbl::modelling::twopt::concentration;
       *  \f$\Delta_c(z)\f$ is computed by
       *  cbl::cosmology::Cosmology::Delta_c and
       *  \f$\rho_{crit}(z)\f$ is computed by
       *  cbl::cosmology::Cosmology::rho_crit; \f$r_{vir}(M_h,
       *  z)\f$ is computed by cbl::cosmology::Cosmology::r_vir
       *
       *  @param rad the scale
       *
       *  @param Mass the dark matter halo mass
       *
       *  @param redshift the redshift
       *
       *  @param model_cM the author(s) of the concentration-mass
       *  relation (see cbl::modelling::twopt::concentration)
       *
       *  @param profile the density profile (see
       *  cbl::modelling::twopt::concentration)
       *
       *  @param halo_def the halo definition (see
       *  cbl::modelling::twopt::concentration)
       *
       *  @return the halo density profile
       */
      double density_profile (const double rad, const double Mass, const double redshift, const std::string model_cM="Duffy", const std::string profile="NFW", const std::string halo_def="vir") const;

      /**
       *  @brief the Fourier transform of the normalised halo density
       *  profile
       *
       *  this function computes the Fourier transform of the
       *  normalised density distribution of dark matter haloes; the
       *  Navarro-Frenk-White profile is the only one currently
       *  implemented (see e.g. eq. 81 of Cooray & Sheth 2002 and 70
       *  of van den Bosch et al. 2012)
       *
       *  \f[\tilde{u}_h(k, M_h, z) = \frac{4\pi\rho_sr_s^3}{M_h}
       *  \left[\cos\mu\left[{\rm Ci}(\mu+\mu c) - {\rm
       *  Ci}(\mu)\right] + \sin\mu\left[{\rm Si}(\mu+\mu c) - {\rm
       *  Si}(\mu)\right] - \frac{\sin\mu c}{\mu+\mu c}\right]\f]
       *
       *  where 
       *
       *  \f[\mu\equiv kr_s\,,\f]
       *
       *  \f[\rho_s = \frac{\rho_{crit}\Delta_c}{3}
       *  \frac{c^3}{\ln(1+c)-c/(1+c)}\,,\f]
       *
       *  \f[{\rm Ci}(x)=-\int_x^\infty\frac{\cos t}{t}\,{\rm d}t\,,\f]
       *
       *  \f[{\rm Si}(x)=-\int_x^\infty\frac{\sin t}{t}\,{\rm d}t\f]
       *
       *  the relation between the halo concentration,
       *  \f$c=c_{vir}=r_{vir}/r_s\f$, and halo mass, \f$M_h\f$, is
       *  computed by cbl::modelling::twopt::concentration;
       *  \f$\Delta_c(z)\f$ is computed by
       *  cbl::cosmology::Cosmology::Delta_c and
       *  \f$\rho_{crit}(z)\f$ is computed by
       *  cbl::cosmology::Cosmology::rho_crit; \f$r_{vir}(M_h,
       *  z)\f$ is computed by cbl::cosmology::Cosmology::r_vir
       *
       *  @param kk the wave vector module at which the model is
       *  computed
       *
       *  @param Mass the halo mass
       *
       *  @param redshift the redshift
       *
       *  @param model_cM the author(s) of the concentration-mass
       *  relation (see cbl::modelling::twopt::concentration)
       *
       *  @param profile the density profile, see
       *  cbl::modelling::twopt::concentration
       *
       *  @param halo_def the halo definition, see
       *  cbl::modelling::twopt::concentration
       *
       *  @return the halo density profile
       */
      double density_profile_FourierSpace (const double kk, const double Mass, const double redshift, const std::string model_cM="Duffy", const std::string profile="NFW", const std::string halo_def="vir") const;

      /**
       *  @brief The halo mass converted to a different value of
       *  \f$\Delta\f$, assuming the Navarro-Frenk-White density
       *  profile
       *
       *  This function converts a given input mass \f$M_\Delta\f$ to
       *  \f$M_{\Delta^{new}}\f$ (e.g. \f$M_{500} \rightarrow
       *  M_{200}\f$).
       *
       *  Specifically, the algorithm currently implemented can be
       *  derived as follows. Given the Navarro-Frenk-White profile:
       *
       *  \f[ \rho(r) = \frac{\rho_0}{\frac{r}{R_s}
       *  \left(1+\frac{r}{R_s} \right)^2} \f]
       *
       *  the total halo mass contained within a radius \f$R_\Delta\f$
       *  is:
       *  
       *  \f[ M_\Delta = \int_0^{R_\Delta} 4\pi r^2\rho(r)dr = 4\pi
       *  \rho_0 R_s^3 \left[
       *  \ln(1+c_\Delta)-\frac{c_\Delta}{1+c_\Delta} \right] \f]
       *
       *  where the concentration is defined as \f$c_\Delta\equiv
       *  R_\Delta/R_s\f$. Thus, we can write:
       *
       *  \f[ \frac{\ln(1+c_\Delta) - \frac{c_\Delta}{1+c_\Delta}}
       *  {M_\Delta} = \frac{\ln(1+c_{\Delta^{new}}) -
       *  \frac{c_{\Delta^{new}}}{1+c_{\Delta^{new}}}}
       *  {M_{\Delta^{new}}} \f]
       *
       *  \f$ M_\Delta \f$ can be written as follows:
       *
       *  \f[ M_\Delta = \frac{4}{3}\pi\Delta\rho_{crit}R_\Delta^3 =
       *  \frac{\Delta}{\Delta^{new}}
       *  \left(\frac{R_\Delta}{R_{\Delta^{new}}} \right)^3
       *  M_{\Delta^{new}} = \frac{\Delta}{\Delta^{new}} x^3
       *  M_{\Delta^{new}} \f]
       *
       *  where \f$ x\equiv R_\Delta/R_{\Delta^{new}} \f$. Thus we
       *  have:
       *
       *  \f[ \frac{M_\Delta}{M_{\Delta^{new}}} \left[
       *  \ln(1+c_{\Delta^{new}}) -
       *  \frac{c_{\Delta^{new}}}{1+c_{\Delta^{new}}} \right] - \left[
       *  \ln(1+c_\Delta) - \frac{c_\Delta}{1+c_\Delta} \right] = 0
       *  \f]
       *
       *  where \f$ c_{\Delta^{new}} \equiv R_{\Delta^{new}}/R_s =
       *  c_\Delta R_{\Delta^{new}}/R_\Delta = c_\Delta/x \f$. The
       *  algorithm solves the above equation as a function of \f$ x
       *  \f$, providing in output:
       *
       *  \f[ M_{\Delta^{new}} =
       *  \frac{1}{x^3}\frac{\Delta^{new}}{\Delta} M_\Delta \f]
       *
       *  @param Mass the input mass \f$M_\Delta\f$
       *  (e.g. \f$M_{500}\f$)
       *
       *  @param Delta_in the input \f$\Delta\f$
       *  (e.g. \f$\Delta=500\f$)
       *
       *  @param Delta_out the output \f$\Delta^{new}\f$
       *  (e.g. \f$\Delta^{new}=200\f$)
       *
       *  @param conc the concentration related to either the input
       *  \f$\Delta\f$, or the output \f$\Delta^{new}\f$ (this is
       *  specified by is_input_conc)
       *
       *  @param is_input_conc true \f$\rightarrow\f$ the given
       *  concentration is related to the input \f$\Delta\f$; false
       *  \f$\rightarrow\f$ the given concentration is related to the
       *  output \f$\Delta^{new}\f$
       * 
       *  @param rRmin_guess the minimum guess value of
       *  \f$r_\Delta/r_{\Delta^{new}}\f$ used by the gsl minimisation
       *  function
       *
       *  @param rRmax_guess the maximum guess value of
       *  \f$r_\Delta/r_{\Delta^{new}}\f$ used by the gsl minimisation
       *  function
       *
       *  @return \f$M_{\Delta^{new}}\f$
       *
       *  @warning the current implementation assumes the
       *  Navarro-Frenk-White profile
       */
      double Mass_Delta (const double Mass, const double Delta_in, const double Delta_out, const double conc, const bool is_input_conc, const double rRmin_guess=1.e-3, const double rRmax_guess=10.) const;
      
      ///@}

 
      /**
       *  @name Functions to estimate the two-point correlation function, bias and related quantities
       */
      ///@{

      /**
       *  @brief the dark matter two-point correlation function
       *
       *  this function provides the dark matter correlation function,
       *  obtained by Fourier transforming the matter power spectrum
       *
       *  @param rr the module of the comoving separation
       *
       *  @param method_Pk method used to compute the power spectrum;
       *  valid choices for method_Pk are: CAMB [http://camb.info/],
       *  classgal_v1 [http://class-code.net/], MPTbreeze-v1
       *  [http://arxiv.org/abs/1207.1465], EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       *
       *  @param redshift the redshift
       *
       *  @param output_root output_root of the parameter file used to
       *  compute the power spectrum and &sigma;(mass); it can be any
       *  name
       *
       *  @param NL 0 \f$\rightarrow\f$ linear power spectrum; 1 \f$\rightarrow\f$
       *  non-linear power spectrum
       *
       *  @param norm 0 \f$\rightarrow\f$ don't normalise the power
       *  spectrum; 1 \f$\rightarrow\f$ normalise the power spectrum;
       *  -1 \f$\rightarrow\f$ normalise only if sigma8 is set
       *
       *  @param k_min minimum wave vector module up to which the power
       *  spectrum is computed
       *
       *  @param k_max maximum wave vector module up to which the power
       *  spectrum is computed
       *
       *  @param aa parameter \e a of Eq. 24 of Anderson et al. 2012
       *
       *  @param GSL false \f$\rightarrow\f$ FFTlog is used; true
       *  \f$\rightarrow\f$ the GSL libraries are used
       *
       *  @param prec accuracy of the integration 
       *
       *  @param file_par name of the parameter file; if a
       *  parameter file is provided (i.e. file_par!=NULL), it will be
       *  used, const ignoring the cosmological parameters of the object
       *
       *  @return &xi;<SUB>DM</SUB>(r): the spherically
       *  averaged (monopole) of the two-point correlation function of
       *  dark matter
       */
      double xi_DM (const double rr, const std::string method_Pk, const double redshift, const std::string output_root="test", const bool NL=true, const int norm=-1, const double k_min=0., const double k_max=100., const double aa=0., const bool GSL=false, const double prec=1.e-2, const std::string file_par=par::defaultString);

      /**
       *  @brief the dark matter angular two-point correlation function
       *
       *  this function provides the dark matter angular correlation function,
       *  obtained by integrating the 2PCF using Limber approximation 
       *  (see Salazar et al. 2014, eqs. 5-8):
       *
       *  \f[ w(\theta) = \int_{z_{min}}^{z_{max}} \mathrm{d}z_1
       *  \int_{z_{min}}^{z_{max}} \mathrm{d}z_2 \phi(z_1) \phi(z_2)
       *  \xi(r,z) \f]
       *
       *  where \f$ r = \sqrt{\chi(z_1)^2+\chi(z_2)^2
       *  -2\chi(z_1)\chi(z_2)\cos(\theta)} \f$, \f$\chi\f$ is the
       *  comoving distance and \f$z = (z_1+z_2)/2 \f$
       *
       *  @param theta the angular separation
       *
       *  @param zz the redshift range
       *
       *  @param phiz the number density
       *
       *  @param interpolationMethod the method in interpolation
       *
       *  @param coordUnits the angular separation units
       *
       *  @param GSL false \f$\rightarrow\f$ FFTlog is used; true
       *  \f$\rightarrow\f$ the GSL libraries are used
       *
       *  @param method_Pk method used to compute the power spectrum;
       *  valid choices for method_Pk are: CAMB [http://camb.info/],
       *  classgal_v1 [http://class-code.net/], MPTbreeze-v1
       *  [http://arxiv.org/abs/1207.1465], EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       *
       *  @param NL 0 \f$\rightarrow\f$ linear power spectrum; 1
       *  \f$\rightarrow\f$ non-linear power spectrum
       *
       *  @param output_root output_root of the parameter file used to
       *  compute the power spectrum and &sigma;(mass); it can be any
       *  name
       *
       *  @param norm 0 \f$\rightarrow\f$ don't normalize the power
       *  spectrum; 1 \f$\rightarrow\f$ normalize the power spectrum;
       *  -1 \f$\rightarrow\f$ normalise only if sigma8 is set
       *
       *  @param k_min minimum wave vector module up to which the power
       *  spectrum is computed
       *
       *  @param k_max maximum wave vector module up to which the power
       *  spectrum is computed
       *
       *  @param prec accuracy of the integration 
       *
       *  @param file_par name of the parameter file; if a parameter
       *  file is provided (i.e. file_par!=NULL), it will be used,
       *  ignoring the cosmological parameters of the object
       *
       *  @return \f$w_{DM}(\theta)\f$: the angular two point
       *  correlation function of dark matter
       */
      double wtheta_DM (const double theta, const std::vector<double> zz, const std::vector<double> phiz, const std::string interpolationMethod, const CoordinateUnits coordUnits = CoordinateUnits::_degrees_, const bool GSL=false, const std::string method_Pk="CAMB", const bool NL=false, const std::string output_root="test", const int norm=-1, const double k_min=1.e-4, const double k_max=100, const double prec=1.e-2, const std::string file_par=par::defaultString);

      /**
       *  @brief the dark matter angular two-point correlation function
       *
       *  this function provides the dark matter angular correlation function,
       *  obtained by integrating the 2PCF using Limber approximation 
       *  (see Sawangwit et al. 2011, eqs. 14,15):
       *
       *  \f[ w(\theta) = \frac{\int_{z_{min}}^{z_{max}} \mathrm{d}z_1
       *  \int_{z_{min}}^{z_{max}} \mathrm{d}z_2 f(z_1) f(z_2)
       *  \xi(r,z)} { \left\{ \int_{z_{min}}^{z_{max}} f(z) dz
       *  \right\}^2} \f]
       *
       *  where \f$ r = \sqrt{\chi(z_1)^2+\chi(z_2)^2
       *  -2\chi(z_1)\chi(z_2)\cos(\theta)} \f$, \f$\chi\f$ is the
       *  comoving distance and \f$z = (z_1+z_2)/2 \f$.
       *
       *  The function \f$f(z)\f$ is the number of objects per unitar
       *  volume:
       *
       *  \f[ f(z) = \frac{\mathrm{d}V}{\mathrm{d}z \mathrm{d\Omega}}
       *     n(z) \phi(z) \f]
       *
       *  where \f$\frac{\mathrm{d}V}{\mathrm{d}z \mathrm{d\Omega}}\f$
       *  is the comoving volume element, \f$n(z)\f$ is the comoving
       *  number density and \f$\phi(z)\f$ is the selection function
       *
       *  @param theta the angular separation
       *
       *  @param kk the wave vector module
       *
       *  @param Pk linear power spectrum
       *
       *  @param zz the redshift range
       *
       *  @param nz the comoving number density
       *
       *  @param phiz the selection function
       *
       *  @param interpolationType the method in interpolation
       *
       *  @param coordUnits the angular separation units
       *
       *  @param GSL false \f$\rightarrow\f$ FFTlog is used; true
       *  \f$\rightarrow\f$ the GSL libraries are used
       *
       *  @param redshift_Pk the redshift of the input power spectrum
       *
       *  @return \f$w_{DM}(\theta)\f$: the angular two point
       *  correlation function of dark matter
       */
      double wtheta_DM (const double theta, const std::vector<double> kk, const std::vector<double> Pk, const std::vector<double> zz, const std::vector<double> nz, const std::vector<double> phiz, const std::string interpolationType="Spline", const CoordinateUnits coordUnits = CoordinateUnits::_degrees_, const bool GSL=false, const double redshift_Pk=0);

      /**
       * @brief the dark matter angular linear power spectrum
       * \f$C_l\f$.
       *
       * this function provides the angular linear power spectrum
       * using the Limber approximation up to a given \f$l_{max}\f$:
       * \f[ C_l = \int_{z_{min}}^{z_{max}}\phi^2(z) D^2(z) P(\frac{l}{D_C(z)}) \frac{H(z)}{c \cdot D_c(z)} dz \f]
       *
       * where l is the multipole order, \f$\phi(z)\f$ is the tracers redshift distribution, \f$P\f$ is the \f$z=0\f$ power
       * spectrum and \f$D(z)\f$ is the growth factor. \f$D_C(z)\f$ and \f$H(z)\f$ are
       * respectively the comoving function and the Hubble parameter.
       *
       * @param lmax the maximum multipole order
       *
       * @param zz the redshift range
       *
       * @param phiz the number density
       *
       * @param interpolationMethod the method in interpolation
       *
       * @param method_Pk method used to compute the power spectrum;
       * valid choices for method_Pk are: CAMB [http://camb.info/],
       * classgal_v1 [http://class-code.net/], MPTbreeze-v1
       * [http://arxiv.org/abs/1207.1465], EisensteinHu
       * [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       *
       * @param output_root output_root of the parameter file used to
       * compute the power spectrum and &sigma;(mass); it can be any
       * name
       *
       * @param norm 0 \f$\rightarrow\f$ don't normalize the power
       * spectrum; 1 \f$\rightarrow\f$ normalize the power spectrum;
       * -1 \f$\rightarrow\f$ normalise only if sigma8 is set
       *
       * @param k_min minimum wave vector module up to which the power
       * spectrum is computed
       *
       * @param k_max maximum wave vector module up to which the power
       * spectrum is computed
       *
       * @param prec accuracy of the integration
       *
       * @param file_par name of the parameter file; if a parameter
       * file is provided (i.e. file_par!=NULL), it will be used,
       * ignoring the cosmological parameters of the object
       *
       * @return vector containing the angular linear power spectrum up to \f$l_{max}\f$
       *
       */
      std::vector<double> C_l_DM (const int lmax, const std::vector<double> zz, const std::vector<double> phiz, const std::string interpolationMethod, const std::string method_Pk="CAMB", const std::string output_root="test", const int norm=-1, const double k_min=1.e-4, const double k_max=100, const double prec=1.e-2, const std::string file_par=par::defaultString);

      /**
       *  @brief the dark matter two-point correlation function,
       *  de-wiggled (see e.g. Anderson et al 2014)
       *
       *  this function provides the dark matter correlation function,
       *  obtained by Fourier transforming the De-Wiggled matter power
       *  spectrum
       *
       *  @author Alfonso Veropalumbo
       *  @author alfonso.veropalumbo@unibo.it
       *  
       *  @param rr the module of the comoving separation
       *
       *  @param redshift the redshift
       *
       *  @param sigma_NL the non linear BAO damping
       *
       *  @param output_root output_root of the parameter file used to
       *  compute the power spectrum and &sigma;(mass); it can be any
       *  name
       *
       *  @param norm 0 \f$\rightarrow\f$ don't normalise the power
       *  spectrum; 1 \f$\rightarrow\f$ normalise the power spectrum;
       *  -1 \f$\rightarrow\f$ normalise only if sigma8 is set
       *
       *  @param k_min minimum wave vector module up to which the power
       *  spectrum is computed
       *
       *  @param k_max maximum wave vector module up to which the power
       *  spectrum is computed
       *
       *  @param aa parameter \e a of Eq. 24 of Anderson et al. 2012
       *
       *  @param prec accuracy of the integration 
       *
       *  @return &xi;<SUB>DW</SUB>(r): the De-Wiggled spherically
       *  averaged (monopole) of the two-point correlation function of
       *  dark matter
       */
      double xi_DM_DeWiggle (const double rr, const double redshift, const double sigma_NL, const std::string output_root = "test", const bool norm=1, const double k_min=0., const double k_max=100., const double aa=1., const double prec=1.e-2);

      /**
       *  @brief get the dark matter two-point correlation function
       *
       *  this function provides the dark matter correlation function,
       *  obtained by Fourier transforming the matter power spectrum
       *
       *  @param [out] rr vector of r, the module of the comoving
       *  separation
       *
       *  @param [out] Xi vector of &xi;(r), the two-point correlation
       *  function of dark matter
       *
       *  @param [in] method_Pk method used to compute the power
       *  spectrum; valid choices for method_Pk are: CAMB
       *  [http://camb.info/], classgal_v1 [http://class-code.net/],
       *  MPTbreeze-v1 [http://arxiv.org/abs/1207.1465], EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       *
       *  @param [in] redshift redshift
       *
       *  @param [in] output_root output_root of the parameter file used
       *  to compute the power spectrum and &sigma;(mass); it can be any
       *  name
       *
       *  @param [in] xiType 0 \f$\rightarrow\f$ standard; 1 \f$\rightarrow\f$ Chuang & Wang
       *  model
       *
       *  @param [in] k_star k<SUB>*</SUB> of the Chuang & Wang model
       *
       *  @param [in] xiNL 0 \f$\rightarrow\f$ linear power spectrum; 1 \f$\rightarrow\f$
       *  non-linear power spectrum
       *
       *  @param [in] norm 0 \f$\rightarrow\f$ don't normalise the
       *  power spectrum; 1 \f$\rightarrow\f$ normalise the power
       *  spectrum; -1 \f$\rightarrow\f$ normalise only if sigma8 is
       *  set
       *
       *  @param [in] r_min minimum separation up to which the
       *  correlation function is computed
       *
       *  @param [in] r_max maximum separation up to which the
       *  correlation function is computed
       *
       *  @param [in] k_min minimum wave vector module up to which the
       *  power spectrum is computed
       *
       *  @param [in] k_max maximum wave vector module up to which the
       *  power spectrum is computed
       *
       *  @param [in] aa parameter \e a of Eq. 24 of Anderson et
       *  al. 2012
       *
       *  @param [in] GSL false \f$\rightarrow\f$ FFTlog is used; true
       *  \f$\rightarrow\f$ the GSL libraries are used
       *
       *  @param [in] prec accuracy of the integration
       *
       *  @param [in] file_par name of the parameter file; if a
       *  parameter file is provided (i.e. file_par!=NULL), it will be
       *  used, const ignoring the cosmological parameters of the object
       *
       *  @return none
       */
      void get_xi (std::vector<double> &rr, std::vector<double> &Xi, const std::string method_Pk, const double redshift, const std::string output_root="test", const bool xiType=0, const double k_star=-1., const bool xiNL=0, const int norm=-1, const double r_min=0.1, const double r_max=150., const double k_min=0., const double k_max=100., const double aa=0., const bool GSL=false, const double prec=1.e-2, const std::string file_par=par::defaultString);
  
      /**
       *  @brief get the barred dark matter correlation functions
       *
       *  this function provides the dark matter \e barred correlation
       *  functions, used to model the two-point correlation function
       *  in redshift-space
       *
       *  @param [in] rr vector of r, the module of the comoving
       *  separation
       *
       *  @param [in] Xi vector of &xi;(r), the two-point correlation
       *  function of dark matter
       *
       *  @param [out] Xi_ vector of barred &xi;(r),
       *
       *  @param [out] Xi__ vector of double-barred &xi;(r)
       *
       *  @param [in] method_Pk method used to compute the power
       *  spectrum; valid choices for method_Pk are: CAMB
       *  [http://camb.info/], classgal_v1 [http://class-code.net/],
       *  MPTbreeze-v1 [http://arxiv.org/abs/1207.1465], EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       *
       *  @param [in] redshift redshift
       *
       *  @param [in] xiType 0 \f$\rightarrow\f$ standard; 1 \f$\rightarrow\f$ Chuang & Wang
       *  model
       *
       *  @param [in] k_star k<SUB>*</SUB> of the Chuang & Wang model
       *
       *  @param [in] xiNL 0 \f$\rightarrow\f$ linear power spectrum; 1 \f$\rightarrow\f$
       *  non-linear power spectrum
       *
       *  @param [in] norm 0 \f$\rightarrow\f$ don't normalise the
       *  power spectrum; 1 \f$\rightarrow\f$ normalise the power
       *  spectrum; -1 \f$\rightarrow\f$ normalise only if sigma8 is
       *  set
       *
       *  @param [in] r_min minimum separation up to which the
       *  correlation function is computed
       *
       *  @param [in] r_max maximum separation up to which the
       *  correlation function is computed
       *
       *  @param [in] k_min minimum wave vector module up to which the
       *  power spectrum is computed
       *
       *  @param [in] k_max maximum wave vector module up to which the
       *  power spectrum is computed
       *
       *  @param [in] aa parameter \e a of Eq. 24 of Anderson et
       *  al. 2012
       *
       *  @param [in] prec accuracy of the integration
       *
       *  @param [in] file_par name of the parameter file; if a
       *  parameter file is provided (i.e. file_par!=NULL), it will be
       *  used, ignoring the cosmological parameters of the object
       *
       *  @return none
       */
      void get_barred_xi (std::vector<double> rr, std::vector<double> Xi, std::vector<double> &Xi_, std::vector<double> &Xi__, const std::string method_Pk, const double redshift, const bool xiType=0, const double k_star=-1., const bool xiNL=0, const int norm=-1, const double r_min=0.1, const double r_max=150., const double k_min=0., const double k_max=100., const double aa=0., const double prec=1.e-2, const std::string file_par=par::defaultString) const;

      /**
       *  @brief the dark matter projected correlation function
       *
       *  this function provides the dark matter projected correlation
       *  functions, obtained by Fourier transforming the matter power
       *  spectrum
       *
       *  @param rp r<SUB>p</SUB>: projected separation
       *
       *  @param method_Pk method used to compute the power spectrum;
       *  valid choices for method_Pk are: CAMB [http://camb.info/],
       *  classgal_v1 [http://class-code.net/], MPTbreeze-v1
       *  [http://arxiv.org/abs/1207.1465], EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       *
       *  @param redshift the redshift
       *
       *  @param pimax the upper limit of the line-of-sight integration
       *
       *  @param output_root output_root of the parameter file used to
       *  compute the power spectrum and &sigma;(mass); it can be any
       *  name
       *
       *  @param NL 0 \f$\rightarrow\f$ linear power spectrum; 1 \f$\rightarrow\f$
       *  non-linear power spectrum
       *
       *  @param norm 0 \f$\rightarrow\f$ don't normalise the power
       *  spectrum; 1 \f$\rightarrow\f$ normalise the power spectrum;
       *  -1 \f$\rightarrow\f$ normalise only if sigma8 is set
       *
       *  @param r_min minimum separation up to which the
       *  correlation function is computed
       *
       *  @param r_max maximum separation up to which the
       *  correlation function is computed
       *
       *  @param k_min minimum wave vector module up to which the
       *  power spectrum is computed
       *
       *  @param k_max maximum wave vector module up to which the
       *  power spectrum is computed
       *
       *  @param aa parameter \e a of Eq. 24 of Anderson et
       *  al. 2012
       *
       *  @param GSL false \f$\rightarrow\f$ FFTlog is used; true
       *  \f$\rightarrow\f$ the GSL libraries are used
       *
       *  @param prec accuracy of the integration
       *
       *  @param file_par name of the parameter file; if a
       *  parameter file is provided (i.e. file_par!=NULL), it will be
       *  used, ignoring the cosmological parameters of the object
       *
       *  @return w<SUB>p,DM</SUB>(&theta;): the projected correlation
       *  function of dark matter
       */
      double wp_DM (const double rp, const std::string method_Pk, const double redshift, const double pimax, const std::string output_root="test", const bool NL=1, const int norm=-1, const double r_min=1.e-3, const double r_max=350., const double k_min=0., const double k_max=100., const double aa=0., const bool GSL=false, const double prec=1.e-2, const std::string file_par=cbl::par::defaultString);

      /**
       *  @brief the k<SUB>*</SUB> parameter 
       *
       *  this function provides the k<SUB>*</SUB> parameter used to
       *  model the BAO (see e.g. Chuang & Wang 2012, Crocce &
       *  Scoccimarro2006, Matsubara 2008)
       *
       *  @param method_Pk method used to compute the power spectrum;
       *  valid choices for method_Pk are: CAMB [http://camb.info/],
       *  classgal_v1 [http://class-code.net/], MPTbreeze-v1
       *  [http://arxiv.org/abs/1207.1465], EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       *
       *  @param redshift the redshift
       *
       *  @param output_root output_root of the parameter file used to
       *  compute the power spectrum and &sigma;(mass); it can be any
       *  name
       *
       *  @param k_max maximum wave vector module up to which the
       *  power spectrum is computed
       *
       *  @param file_par name of the parameter file; if a
       *  parameter file is provided (i.e. file_par!=NULL), it will be
       *  used, ignoring the cosmological parameters of the object
       *
       *  @return k<SUB>*</SUB>
       */
      double k_star (const std::string method_Pk, const double redshift, const std::string output_root="test", const double k_max=100., const std::string file_par=par::defaultString) const; 


      /**
       *  @brief the dark matter rms mass fluctuation
       *
       *  @param RR radius inside which the dark matter rms mass
       *  fluctuation is computed
       *
       *  @param corrType 0 \f$\rightarrow\f$ the projected correlation function,
       *  w(&theta;), is used; 1 \f$\rightarrow\f$ the spherically averaged
       *  correlation function, &xi;(r), is used
       *
       *  @param method_Pk method used to compute the power spectrum;
       *  valid choices for method_Pk are: CAMB [http://camb.info/],
       *  classgal_v1 [http://class-code.net/], MPTbreeze-v1
       *  [http://arxiv.org/abs/1207.1465], EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       *
       *  @param redshift the redshift
       *
       *  @param pimax the upper limit of the line-of-sight integration
       *
       *  @param output_root output_root of the parameter file used to
       *  compute the power spectrum and &sigma;(mass); it can be any
       *  name
       *  
       *  @param NL 0 \f$\rightarrow\f$ linear power spectrum; 1 \f$\rightarrow\f$
       *  non-linear power spectrum
       *
       *  @param norm 0 \f$\rightarrow\f$ don't normalise the power
       *  spectrum; 1 \f$\rightarrow\f$ normalise the power spectrum;
       *  -1 \f$\rightarrow\f$ normalise only if sigma8 is set
       *
       *  @param r_min minimum separation up to which the
       *  correlation function is computed
       *
       *  @param r_max maximum separation up to which the
       *  correlation function is computed
       *
       *  @param k_min minimum wave vector module up to which the
       *  power spectrum is computed
       *
       *  @param k_max maximum wave vector module up to which the
       *  power spectrum is computed
       *
       *  @param aa parameter \e a of Eq. 24 of Anderson et
       *  al. 2012
       *
       *  @param GSL false \f$\rightarrow\f$ FFTlog is used; true
       *  \f$\rightarrow\f$ the GSL libraries are used
       *
       *  @param prec accuracy of the integration
       *
       *  @param file_par name of the parameter file; if a
       *  parameter file is provided (i.e. file_par!=NULL), it will be
       *  used, ignoring the cosmological parameters of the object
       *
       *  @return &sigma;<SUB>R</SUB>: the dark matter rms mass
       *  fluctuation
       */
      double sigmaR_DM (const double RR, const int corrType, const std::string method_Pk, const double redshift, const double pimax=40, const std::string output_root="test", const bool NL=1, const int norm=-1, const double r_min=1.e-3, const double r_max=350., const double k_min=0., const double k_max=100., const double aa=0., const bool GSL=false, const double prec=1.e-2, const std::string file_par=par::defaultString); 

      /**
       *  @brief the dark matter rms mass fluctuation within 8 Mpc/h
       *
       *  this function provides the rms mass fluctuation within 8
       *  Mpc/h, estimated directly from the power spectrum
       *
       *  @param method_Pk method used to compute the power spectrum;
       *  valid choices for method_Pk are: CAMB [http://camb.info/],
       *  classgal_v1 [http://class-code.net/], MPTbreeze-v1
       *  [http://arxiv.org/abs/1207.1465], EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       *
       *  @param redshift the redshift
       *
       *  @param output_root output_root of the parameter file used to
       *  compute the power spectrum and &sigma;(mass); it can be any
       *  name
       *  
       *  @param NL 0 \f$\rightarrow\f$ linear power spectrum; 1 \f$\rightarrow\f$
       *  non-linear power spectrum
       *
       *  @param k_min minimum wave vector module up to which the
       *  power spectrum is computed
       *
       *  @param k_max maximum wave vector module up to which the
       *  power spectrum is computed
       *
       *  @param prec accuracy of the integration
       *
       *  @param file_par name of the parameter file; if a
       *  parameter file is provided (i.e. file_par!=NULL), it will be
       *  used, ignoring the cosmological parameters of the object
       *
       *  @return &sigma;<SUB>8</SUB>: the dark matter rms mass
       *  fluctuation within 8 Mpc/h
       */
      double sigma8_Pk (const std::string method_Pk, const double redshift, const std::string output_root="test", const bool NL=0, const double k_min=0., const double k_max=100., const double prec=1.e-2, const std::string file_par=par::defaultString) const; 

      /**
       *  @brief bias of dark matter haloes
       *
       *  @param Mass halo mass
       *
       *  @param redshift the redshift
       *
       *  @param author author(s) who proposed the bias; valid authors
       *  are: ST99 (Sheth & Tormen 1999), SMT01 (Sheth, Mo & Tormen
       *  2001), SMT01_WL04 (Sheth, Mo & Tormen 2001 with the correction
       *  of Warren 2004), Tinker (Tinker et al. 2010)
       *
       *  @param method_SS method used to compute the power spectrum
       *  and &sigma;(mass); valid method_SS are: CAMB
       *  [http://camb.info/], classgal_v1 [http://class-code.net/],
       *  EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       *
       *  @param output_root output_root of the parameter file used to
       *  compute the power spectrum and &sigma;(mass); it can be any
       *  name
       *
       *  @param interpType method to interpolate the power spectrum
       *
       *  @param Delta \f$\Delta\f$: the overdensity, defined as the
       *  mean interior density relative to the background
       *  
       *  @param kk wave vector module
       *
       *  @param norm 0 \f$\rightarrow\f$ don't normalise the power
       *  spectrum; 1 \f$\rightarrow\f$ normalise the power spectrum;
       *  -1 \f$\rightarrow\f$ normalise only if sigma8 is set
       *
       *  @param k_min minimum wave vector module up to which the
       *  power spectrum is computed
       *
       *  @param k_max maximum wave vector module up to which the
       *  power spectrum is computed
       *
       *  @param prec accuracy of the integration
       *
       *  @param input_file either the parameter file or the power
       *  spectrum file; if a parameter file is provided,
       *  i.e. input_file!=NULL and is_parameter_file=true, it will be
       *  used to compute the power spectrum; if a power spectrum file
       *  is provided, i.e. input_file!=NULL and
       *  is_parameter_file=false, then the provided power spectrum
       *  will be used directly; in both cases &sigma;<SUP>2</SUP>(M)
       *  is computed by integrating the computed/provided power
       *  spectrum ignoring the cosmological parameters of the object
       *
       *  @param is_parameter_file true \f$\rightarrow\f$ the input_file
       *  is a parameter file, used to compute the power spectrum with
       *  the method specified by method_Pk; false \f$\rightarrow\f$
       *  the input_file is a file containing the power spectrum
       *
       *  @return b<SUB>halo</SUB>: the dark matter bias
       *
       *  @warning the input parameter \f$\Delta\f$ is the background
       *  overdensity, not the critical overdensity
       *  \f$\Delta_{crit}\f$; the function
       *  cbl::Cosmology::Delta_vir can be used to convert
       *  \f$\Delta_{crit}\f$ into \f$\Delta\f$
       */
      double bias_halo (const double Mass, const double redshift, const std::string author, const std::string method_SS, const std::string output_root="test", const std::string interpType="Linear", const double Delta=200., const double kk=-1., const int norm=-1, const double k_min=0., const double k_max=100., const double prec=1.e-2, const std::string input_file=par::defaultString, const bool is_parameter_file=true); 

      /**
       *  @brief bias of dark matter haloes
       *
       *  @param Mass halo mass
       *
       *  @param Sigma &sigma;(mass, z=0): the mass variance at z=0
       *
       *  @param redshift the redshift
       *
       *  @param model_bias author(s) who proposed the bias; valid
       *  authors are: ST99 (Sheth & Tormen 1999), SMT01 (Sheth, Mo &
       *  Tormen 2001), SMT01_WL04 (Sheth, Mo & Tormen 2001 with the
       *  correction of Warren 2004), Tinker (Tinker et al. 2010)
       *
       *  @param output_root output_root of the parameter file used to
       *  compute the power spectrum and &sigma;(mass); it can be any
       *  name
       *
       *  @param interpType method to interpolate the power spectrum
       *
       *  @param Delta \f$\Delta\f$: the overdensity, defined as the
       *  mean interior density relative to the background
       *  
       *  @param kk wave vector module
       *
       *  @param norm 0 \f$\rightarrow\f$ don't normalise the power
       *  spectrum; 1 \f$\rightarrow\f$ normalise the power spectrum;
       *  -1 \f$\rightarrow\f$ normalise only if sigma8 is set
       *
       *  @param k_min minimum wave vector module up to which the
       *  power spectrum is computed
       *
       *  @param k_max maximum wave vector module up to which the
       *  power spectrum is computed
       *
       *  @param prec accuracy of the integration
       *
       *  @param method_SS method used to compute the power spectrum
       *  and &sigma;(mass); valid method_SS are: CAMB
       *  [http://camb.info/], classgal_v1 [http://class-code.net/],
       *  EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       *    
       *  @param input_file either the parameter file or the power
       *  spectrum file; if a parameter file is provided,
       *  i.e. input_file!=NULL and is_parameter_file=true, it will be
       *  used to compute the power spectrum; if a power spectrum file
       *  is provided, i.e. input_file!=NULL and
       *  is_parameter_file=false, then the provided power spectrum
       *  will be used directly; in both cases &sigma;<SUP>2</SUP>(M)
       *  is computed by integrating the computed/provided power
       *  spectrum ignoring the cosmological parameters of the object
       *
       *  @param is_parameter_file true \f$\rightarrow\f$ the input_file
       *  is a parameter file, used to compute the power spectrum with
       *  the method specified by method_Pk; false \f$\rightarrow\f$
       *  the input_file is a file containing the power spectrum
       *
       *  @return b<SUB>halo</SUB>: the dark matter bias
       *
       *  @warning the input parameter \f$\Delta\f$ is the background
       *  overdensity, not the critical overdensity
       *  \f$\Delta_{crit}\f$; the function
       *  cbl::Cosmology::Delta_vir can be used to convert
       *  \f$\Delta_{crit}\f$ into \f$\Delta\f$
       */
      double bias_halo (const double Mass, const double Sigma, const double redshift, const std::string model_bias, const std::string output_root="test", const std::string interpType="Linear", const double Delta=200., const double kk=-1., const int norm=-1, const double k_min=0., const double k_max=100., const double prec=1.e-2, const std::string method_SS="CAMB", const std::string input_file=par::defaultString, const bool is_parameter_file=true);
  
      /**
       *  @brief the effective bias of dark matter haloes, with masses
       *  in a given range and at a given mean redshift
       *
       *  this function computes the effective bias of dark matter
       *  haloes:
       *
       *  \f[ b_{eff}(z) = \frac{\int_{M_{min}}^{M_{max}} {\rm d}M\,
       *  b(M, z) \Phi(M, z)}{\int_{M_{min}}^{M_{max}} {\rm
       *  d}M\,\Phi(M, z)} \f]
       *
       *  in the current implementation, the integral is actually
       *  replaced by the Riemann sum, as follows:
       *
       *  \f[ b_{eff}(z) \simeq \frac{\sum_{i_{min}}^{i_{max}} b(M, z)
       *  \Phi(M, z) (M_{i+1}-M_i)}{\sum_{i_{min}}^{i_{max}}
       *  \Phi(M, z) (M_{i+1}-M_i)} \f]
       *
       *  where the halo mass function, \f$\Phi(M, z)\f$, is computed
       *  by cbl::cosmology::Cosmology::mass_function and the
       *  linear bias, \f$b(M, z)\f$, is computed by
       *  cbl::cosmology::Cosmology::bias_halo
       *
       *  @param Mass_min minimum halo mass
       *
       *  @param Mass_max maximum halo mass
       *
       *  @param redshift the redshift
       *
       *  @param model_bias author(s) who proposed the bias; valid
       *  authors are: ST99 (Sheth & Tormen 1999), SMT01 (Sheth, Mo &
       *  Tormen 2001), SMT01_WL04 (Sheth, Mo & Tormen 2001 with the
       *  correction of Warren 2004), Tinker (Tinker et al. 2010)
       *
       *  @param model_MF author(s) who proposed the mass function;
       *  valid authors are: PS (Press&Schechter), ST (Sheth&Tormen),
       *  Jenkins (Jenkins et al. 2001), Warren (Warren et al. 2006),
       *  Reed, (Reed et al. 2007), Pan (Pan 2007), ShenH (halo MF by
       *  Shen et al. 2006), ShenF (filaments MF by Shen et al. 2006),
       *  ShenS (sheets MF by Shen et al. 2006), Tinker (Tinker et
       *  al. 2008), Crocce (Crocce et al. 2010), Angulo_FOF (FOF MF
       *  by Angulo et al. 2012), Angulo_Sub (SUBFIND MF by Angulo et
       *  al. 2012), Watson_FOF(FOF MF by Watson et al. 2012),
       *  Watson_SOH (MF for Spherical Overdensity Haloes by Watson et
       *  al. 2012), Manera (Manera et al. 2010), Bhattacharya
       *  (Bhattacharya et al. 2011), Courtin (Courtin et al. 2010),
       *  Peacock (by Peacock at al. 2007)
       *
       *  @param method_SS method used to compute the power spectrum and
       *  &sigma;(mass); valid method_SS are: CAMB [http://camb.info/],
       *  classgal_v1 [http://class-code.net/], EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       *
       *  @param output_root output_root of the parameter file used to
       *  compute the power spectrum and &sigma;(mass); it can be any
       *  name
       *
       *  @param Delta \f$\Delta\f$: the overdensity, defined as the
       *  mean interior density relative to the background
       *  
       *  @param kk wave vector module
       *
       *  @param interpType method to interpolate the power spectrum
       *
       *  @param norm 0 \f$\rightarrow\f$ don't normalise the power
       *  spectrum; 1 \f$\rightarrow\f$ normalise the power spectrum;
       *  -1 \f$\rightarrow\f$ normalise only if sigma8 is set
       *
       *  @param k_min minimum wave vector module up to which the
       *  power spectrum is computed
       *
       *  @param k_max maximum wave vector module up to which the
       *  power spectrum is computed
       *
       *  @param prec accuracy of the integration 
       *
       *  @param input_file either the parameter file or the power
       *  spectrum file; if a parameter file is provided,
       *  i.e. input_file!=NULL and is_parameter_file=true, it will be
       *  used to compute the power spectrum; if a power spectrum file
       *  is provided, i.e. input_file!=NULL and
       *  is_parameter_file=false, then the provided power spectrum
       *  will be used directly; in both cases &sigma;<SUP>2</SUP>(M)
       *  is computed by integrating the computed/provided power
       *  spectrum ignoring the cosmological parameters of the object
       *
       *  @param is_parameter_file true \f$\rightarrow\f$ the input_file
       *  is a parameter file, used to compute the power spectrum with
       *  the method specified by method_Pk; false \f$\rightarrow\f$
       *  the input_file is a file containing the power spectrum
       *
       *  @return b<SUB>eff</SUB>: the effective dark matter bias
       *
       *  @warning the input parameter \f$\Delta\f$ is the background
       *  overdensity, not the critical overdensity
       *  \f$\Delta_{crit}\f$; the function
       *  cbl::Cosmology::Delta_vir can be used to convert
       *  \f$\Delta_{crit}\f$ into \f$\Delta\f$
       *
       *  @warning the mass function by Manera et al. (2010) has been
       *  tested only for z=0 and z=0.5; the mass function by Despali
       *  et al. (2016) is currently implemented only for virial
       *  masses and at \f$z<1.25\f$
       */
      double bias_eff (const double Mass_min, const double Mass_max, const double redshift, const std::string model_bias, const std::string model_MF, const std::string method_SS, const std::string output_root="test", const double Delta=200., const double kk=-1., const std::string interpType="Linear", const int norm=-1, const double k_min=0., const double k_max=100., const double prec=1.e-2, const std::string input_file=par::defaultString, const bool is_parameter_file=true);
 
      /**
       *  @brief effective bias of dark matter haloes, computed by
       *  averaging over the bias of a given set of haloes
       *
       *  this function computes the effective bias of dark matter
       *  haloes as follows:
       *
       *  \f[ b_{eff}(z) = \frac{\int_{M_{min}}^{M_{max}} {\rm d}M\,
       *  b(M, z) \Phi(M, z)}{\int_{M_{min}}^{M_{max}} {\rm
       *  d}M\,\Phi(M, z)} \f]
       *
       *  in the current implementation, the integral is actually
       *  replaced by the Riemann sum, as follows:
       *
       *  \f[ b_{eff}(z) \simeq \frac{\sum_{i_{min}}^{i_{max}} b(M, z)
       *  \Phi(M, z) (M_{i+1}-M_i)}{\sum_{i_{min}}^{i_{max}}
       *  \Phi(M, z) (M_{i+1}-M_i)} \f]
       *
       *  where the halo mass function is provided in input and the
       *  linear bias, \f$b(M, z)\f$, is computed by
       *  cbl::cosmology::Cosmology::bias_halo
       *
       *  @param MM vector of halo masses
       *
       *  @param MF vector of mass function values, d&Phi;/dM=dn(M)/dM
       *
       *  @param redshift the redshift
       *
       *  @param model_bias author(s) who proposed the bias; valid
       *  authors are: ST99 (Sheth & Tormen 1999), SMT01 (Sheth, Mo &
       *  Tormen 2001), SMT01_WL04 (Sheth, Mo & Tormen 2001 with the
       *  correction of Warren 2004), Tinker (Tinker et al. 2010)
       *
       *  @param method_SS method used to compute the power spectrum
       *  and &sigma;(mass); valid method_SS are: CAMB
       *  [http://camb.info/], classgal_v1 [http://class-code.net/],
       *  EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       *
       *  @param output_root output_root of the parameter file used to
       *  compute the power spectrum and &sigma;(mass); it can be any
       *  name
       *
       *  @param Delta \f$\Delta\f$: the overdensity, defined as the
       *  mean interior density relative to the background
       *  
       *  @param kk wave vector module
       *
       *  @param interpType method to interpolate the power spectrum
       *
       *  @param norm 0 \f$\rightarrow\f$ don't normalise the power
       *  spectrum; 1 \f$\rightarrow\f$ normalise the power spectrum;
       *  -1 \f$\rightarrow\f$ normalise only if sigma8 is set
       *
       *  @param k_min minimum wave vector module up to which the
       *  power spectrum is computed
       *
       *  @param k_max maximum wave vector module up to which the
       *  power spectrum is computed
       *
       *  @param prec accuracy of the integration 
       *
       *  @param input_file either the parameter file or the power
       *  spectrum file; if a parameter file is provided,
       *  i.e. input_file!=NULL and is_parameter_file=true, it will be
       *  used to compute the power spectrum; if a power spectrum file
       *  is provided, i.e. input_file!=NULL and
       *  is_parameter_file=false, then the provided power spectrum
       *  will be used directly; in both cases &sigma;<SUP>2</SUP>(M)
       *  is computed by integrating the computed/provided power
       *  spectrum ignoring the cosmological parameters of the object
       *
       *  @param is_parameter_file true \f$\rightarrow\f$ the input_file
       *  is a parameter file, used to compute the power spectrum with
       *  the method specified by method_Pk; false \f$\rightarrow\f$
       *  the input_file is a file containing the power spectrum
       *
       *  @return b<SUB>eff</SUB>: the effective dark matter bias
       *
       *  @warning the input parameter \f$\Delta\f$ is the background
       *  overdensity, not the critical overdensity
       *  \f$\Delta_{crit}\f$; the function
       *  cbl::Cosmology::Delta_vir can be used to convert
       *  \f$\Delta_{crit}\f$ into \f$\Delta\f$
       */
      double bias_eff (const std::vector<double> MM, const std::vector<double> MF, const double redshift, const std::string model_bias, const std::string method_SS, const std::string output_root="test", const double Delta=200., const double kk=-1., const std::string interpType="Linear", const int norm=-1, const double k_min=0., const double k_max=100., const double prec=1.e-2, const std::string input_file=par::defaultString, const bool is_parameter_file=true);


      /**
       *  @brief effective bias of dark matter haloes, computed by
       *  averaging the bias of a set of haloes, with the mass variance
       *  estimated from a grid
       *
       *  this function computes the effective bias of dark matter
       *  haloes by either averaging the bias of a set of haloes with
       *  a given mass:
       *
       *  \f[b_{eff}(z) = \frac{1}{N_{halo}}\sum_{i=1}^{N_{halo}}
       *  b(M_i, z_i) \; , \; (1)\f]
       *
       *  or by averaging over halo pairs:
       *
       *  \f[b_{eff}(z) = \sqrt{ \frac{2}{N_{halo}(N_{halo}-1)}
       *  \sum_{i=1}^{N_{halo}}\sum_{j=i+1}^{N_{halo}} b(M_i,
       *  z_i)b(M_j, z_j)} \; , \; (2)\f]
       *
       *  where the linear bias of the \f$i\f$-th halo, \f$b^{i}(M,
       *  z)\f$, is computed by
       *  cbl::cosmology::Cosmology::bias_halo; the mass variance
       *  is estimated from a grid by
       *  cbl::cosmology::Cosmology::create_grid_sigmaM
       *
       *  @param MM vector containing the halo masses
       *
       *  @param redshift vector containing the redshifts; if it has
       *  size=1, it will be considered as the main redshift
       *
       *  @param model_bias author(s) who proposed the bias; valid
       *  authors are: ST99 (Sheth & Tormen 1999), SMT01 (Sheth, Mo &
       *  Tormen 2001), SMT01_WL04 (Sheth, Mo & Tormen 2001 with the
       *  correction of Warren 2004), Tinker (Tinker et al. 2010)
       *
       *  @param method_SS method used to compute the power spectrum
       *  and &sigma;(mass); valid method_SS are: CAMB
       *  [http://camb.info/], classgal_v1 [http://class-code.net/],
       *  EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       *
       *  @param meanType meanType="mean_bias" \f$\rightarrow\f$ the
       *  effective bias is computed with Eq.(1);
       *  meanType="mean_pair_bias" \f$\rightarrow\f$ the effective
       *  bias is computed with Eq.(2)
       * 
       *  @param output_root output_root of the parameter file used to
       *  compute the power spectrum and &sigma;(mass); it can be any
       *  name
       *
       *  @param Delta_crit \f$\Delta_{crit}\f$: the critical
       *  overdensity
       *  
       *  @param kk wave vector module
       *
       *  @param interpType method to interpolate the power spectrum
       *
       *  @param norm 0 \f$\rightarrow\f$ don't normalise the power
       *  spectrum; 1 \f$\rightarrow\f$ normalise the power spectrum;
       *  -1 \f$\rightarrow\f$ normalise only if sigma8 is set
       *
       *  @param k_min minimum wave vector module up to which the
       *  power spectrum is computed
       *
       *  @param k_max maximum wave vector module up to which the
       *  power spectrum is computed
       *
       *  @param prec accuracy of the integration 
       *
       *  @param input_file either the parameter file or the power
       *  spectrum file; if a parameter file is provided,
       *  i.e. input_file!=NULL and is_parameter_file=true, it will be
       *  used to compute the power spectrum; if a power spectrum file
       *  is provided, i.e. input_file!=NULL and
       *  is_parameter_file=false, then the provided power spectrum
       *  will be used directly; in both cases &sigma;<SUP>2</SUP>(M)
       *  is computed by integrating the computed/provided power
       *  spectrum ignoring the cosmological parameters of the object
       *
       *  @param is_parameter_file true \f$\rightarrow\f$ the input_file
       *  is a parameter file, used to compute the power spectrum with
       *  the method specified by method_Pk; false \f$\rightarrow\f$
       *  the input_file is a file containing the power spectrum
       *
       *  @return a vector containing the mean and standard deviation
       *  of the effective dark matter bias
       */
      std::vector<double> bias_eff_mass_grid (const std::vector<double> MM, const std::vector<double> redshift, const std::string model_bias, const std::string method_SS, const std::string meanType="mean_bias", const std::string output_root="test", const double Delta_crit=200., const double kk=-1., const std::string interpType="Linear", const int norm=-1, const double k_min=0., const double k_max=100., const double prec=1.e-2, const std::string input_file=par::defaultString, const bool is_parameter_file=true);

      /**
       *  @brief effective bias of dark matter haloes, computed by
       *  averaging the bias of a set of haloes
       *
       *  this function computes the effective bias of dark matter
       *  haloes by either averaging the bias of a set of haloes with
       *  a given mass:
       *
       *  \f[b_{eff}(z) = \frac{1}{N_{halo}}\sum_{i=1}^{N_{halo}}
       *  b(M_i, z_i) \; , \; (1)\f]
       *
       *  or by averaging over halo pairs:
       *
       *  \f[b_{eff}(z) = \sqrt{ \frac{2}{N_{halo}(N_{halo}-1)}
       *  \sum_{i=1}^{N_{halo}}\sum_{j=i+1}^{N_{halo}} b(M_i,
       *  z_i)b(M_j, z_j)} \; , \; (2)\f]
       *
       *  where the linear bias of the \f$i\f$-th halo, \f$b^{i}(M,
       *  z)\f$, is computed by
       *  cbl::cosmology::Cosmology::bias_halo; the mass variance
       *  is computed by cbl::cosmology::Cosmology::sigma2M
       *
       *  @param MM vector containing the halo masses
       *
       *  @param redshift vector containing the redshifts; if it has
       *  size=1, it will be considered as the main redshift
       *
       *  @param model_bias author(s) who proposed the bias; valid
       *  authors are: ST99 (Sheth & Tormen 1999), SMT01 (Sheth, Mo &
       *  Tormen 2001), SMT01_WL04 (Sheth, Mo & Tormen 2001 with the
       *  correction of Warren 2004), Tinker (Tinker et al. 2010)
       *
       *  @param method_SS method used to compute the power spectrum
       *  and &sigma;(mass); valid method_SS are: CAMB
       *  [http://camb.info/], classgal_v1 [http://class-code.net/],
       *  EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       *
       *  @param meanType meanType="mean_bias" \f$\rightarrow\f$ the
       *  effective bias is computed with Eq.(1);
       *  meanType="mean_pair_bias" \f$\rightarrow\f$ the effective
       *  bias is computed with Eq.(2)
       *
       *  @param output_root output_root of the parameter file used to
       *  compute the power spectrum and &sigma;(mass); it can be any
       *  name
       *
       *  @param Delta_crit \f$\Delta_{crit}\f$: the critical
       *  overdensity
       *  
       *  @param kk wave vector module
       *
       *  @param interpType method to interpolate the power spectrum
       *
       *  @param norm 0 \f$\rightarrow\f$ don't normalise the power
       *  spectrum; 1 \f$\rightarrow\f$ normalise the power spectrum;
       *  -1 \f$\rightarrow\f$ normalise only if sigma8 is set
       *
       *  @param k_min minimum wave vector module up to which the
       *  power spectrum is computed
       *
       *  @param k_max maximum wave vector module up to which the
       *  power spectrum is computed
       *
       *  @param prec accuracy of the integration 
       *
       *  @param input_file either the parameter file or the power
       *  spectrum file; if a parameter file is provided,
       *  i.e. input_file!=NULL and is_parameter_file=true, it will be
       *  used to compute the power spectrum; if a power spectrum file
       *  is provided, i.e. input_file!=NULL and
       *  is_parameter_file=false, then the provided power spectrum
       *  will be used directly; in both cases &sigma;<SUP>2</SUP>(M)
       *  is computed by integrating the computed/provided power
       *  spectrum ignoring the cosmological parameters of the object
       *
       *  @param is_parameter_file true \f$\rightarrow\f$ the input_file
       *  is a parameter file, used to compute the power spectrum with
       *  the method specified by method_Pk; false \f$\rightarrow\f$
       *  the input_file is a file containing the power spectrum
       *
       *  @return a vector containing the mean and standard deviation
       *  of the effective dark matter bias
       */
      std::vector<double> bias_eff_mass (const std::vector<double> MM, const std::vector<double> redshift, const std::string model_bias, const std::string method_SS, const std::string meanType="mean_bias", const std::string output_root="test", const double Delta_crit=200., const double kk=-1., const std::string interpType="Linear", const int norm=-1, const double k_min=0., const double k_max=100., const double prec=1.e-2, const std::string input_file=par::defaultString, const bool is_parameter_file=true);
      
      /**
       *  @brief effective bias of dark matter haloes, computed by
       *  averaging the bias of a set of haloes, interpolating the
       *  mass variance on a grid
       *
       *  this function computes the effective bias of dark matter
       *  haloes by either averaging the bias of a set of haloes with
       *  a given mass:
       *
       *  \f[b_{eff}(z) = \frac{1}{N_{halo}}\sum_{i=1}^{N_{halo}}
       *  b(M_i, z_i) \; , \; (1)\f]
       *
       *  or by averaging over halo pairs:
       *
       *  \f[b_{eff}(z) = \sqrt{ \frac{2}{N_{halo}(N_{halo}-1)}
       *  \sum_{i=1}^{N_{halo}}\sum_{j=i+1}^{N_{halo}} b(M_i,
       *  z_i)b(M_j, z_j)} \; , \; (2)\f]
       *
       *  where the linear bias of the \f$i\f$-th halo, \f$b^{i}(M,
       *  z)\f$, is computed by
       *  cbl::cosmology::Cosmology::bias_halo; the mass variance
       *  is computed by cbl::cosmology::Cosmology::sigma2M
       *
       *  @param mass vector containing the halo masses
       *
       *  @param mass_grid vector containing the halo masses on the
       *  grid used to interpolate the mass variance
       *
       *  @param redshift vector containing the redshifts; if it has
       *  size=1, it will be considered as the main redshift
       *
       *  @param model_bias author(s) who proposed the bias; valid
       *  authors are: ST99 (Sheth & Tormen 1999), SMT01 (Sheth, Mo &
       *  Tormen 2001), SMT01_WL04 (Sheth, Mo & Tormen 2001 with the
       *  correction of Warren 2004), Tinker (Tinker et al. 2010)
       *
       *  @param method_SS method used to compute the power spectrum
       *  and &sigma;(mass); valid method_SS are: CAMB
       *  [http://camb.info/], classgal_v1 [http://class-code.net/],
       *  EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       *
       *  @param meanType meanType="mean_bias" \f$\rightarrow\f$ the
       *  effective bias is computed with Eq.(1);
       *  meanType="mean_pair_bias" \f$\rightarrow\f$ the effective
       *  bias is computed with Eq.(2)
       *
       *  @param output_root output_root of the parameter file used to
       *  compute the power spectrum and &sigma;(mass); it can be any
       *  name
       *
       *  @param Delta_crit \f$\Delta_{crit}\f$: the critical
       *  overdensity
       *  
       *  @param kk wave vector module
       *
       *  @param interpType method to interpolate the power spectrum
       *
       *  @param norm 0 \f$\rightarrow\f$ don't normalise the power
       *  spectrum; 1 \f$\rightarrow\f$ normalise the power spectrum;
       *  -1 \f$\rightarrow\f$ normalise only if sigma8 is set
       *
       *  @param k_min minimum wave vector module up to which the
       *  power spectrum is computed
       *
       *  @param k_max maximum wave vector module up to which the
       *  power spectrum is computed
       *
       *  @param prec accuracy of the integration 
       *
       *  @param input_file either the parameter file or the power
       *  spectrum file; if a parameter file is provided,
       *  i.e. input_file!=NULL and is_parameter_file=true, it will be
       *  used to compute the power spectrum; if a power spectrum file
       *  is provided, i.e. input_file!=NULL and
       *  is_parameter_file=false, then the provided power spectrum
       *  will be used directly; in both cases &sigma;<SUP>2</SUP>(M)
       *  is computed by integrating the computed/provided power
       *  spectrum ignoring the cosmological parameters of the object
       *
       *  @param is_parameter_file true \f$\rightarrow\f$ the input_file
       *  is a parameter file, used to compute the power spectrum with
       *  the method specified by method_Pk; false \f$\rightarrow\f$
       *  the input_file is a file containing the power spectrum
       *
       *  @return a vector containing the mean and standard deviation
       *  of the effective dark matter bias
       */
      std::vector<double> bias_eff_mass (const std::vector<double> mass, const std::vector<double> mass_grid, const std::vector<double> redshift, const std::string model_bias, const std::string method_SS, const std::string meanType="mean_bias", const std::string output_root="test", const double Delta_crit=200., const double kk=-1., const std::string interpType="Linear", const int norm=-1, const double k_min=0., const double k_max=100., const double prec=1.e-2, const std::string input_file=par::defaultString, const bool is_parameter_file=true);
      
      /**
       *  @brief compute the effective bias of dark matter haloes, by
       *  averaging the bias of a set of haloes, interpolating the
       *  mass variance on a grid of masses and of one input
       *  cosmological parameter; this function is used when modelling
       *  the two-point correlation function
       *
       *  this function computes the effective bias of dark matter
       *  haloes by either averaging the bias of a set of haloes with
       *  a given mass:
       *
       *  \f[b_{eff}(z) = \frac{1}{N_{halo}}\sum_{i=1}^{N_{halo}}
       *  b(M_i, z_i) \; , \; (1)\f]
       *
       *  or by averaging over halo pairs:
       *
       *  \f[b_{eff}(z) = \sqrt{ \frac{2}{N_{halo}(N_{halo}-1)}
       *  \sum_{i=1}^{N_{halo}}\sum_{j=i+1}^{N_{halo}} b(M_i,
       *  z_i)b(M_j, z_j)} \; , \; (2)\f]
       *
       *  where the linear bias of the \f$i\f$-th halo, \f$b^{i}(M,
       *  z)\f$, is computed by
       *  cbl::cosmology::Cosmology::bias_halo; the mass variance
       *  is computed by cbl::cosmology::Cosmology::sigma2M
       *
       *  @param parameter vector containing the grid of the
       *  cosmological parameters on which the effective bias grid is
       *  computed
       *
       *  @param bias_eff vector containing the effective bias grid
       *
       *  @param dir_output the directory where the effective bias
       *  grid is stored
       *
       *  @param file_bias_eff_grid the file there the effective bias
       *  grid is stored
       *
       *  @param cosmoPar the cosmological parameter for which the
       *  effective bias grid is computed
       *
       *  @param min_par the minimum value for the
       *  parameter where the effective bias is computed
       *  
       *  @param max_par the maximum value for the
       *  parameter where the effective bias is computed
       *
       *  @param nbin_par the number of points for the
       *  parameter where the effective bias is computed
       *
       *  @param mass vector containing the halo masses
       *
       *  @param mass_grid vector containing the halo masses on the
       *  grid used to interpolate the mass variance
       *
       *  @param redshift vector containing the redshifts; if it has
       *  size=1, it will be considered as the main redshift
       *
       *  @param model_bias author(s) who proposed the bias; valid
       *  authors are: ST99 (Sheth & Tormen 1999), SMT01 (Sheth, Mo &
       *  Tormen 2001), SMT01_WL04 (Sheth, Mo & Tormen 2001 with the
       *  correction of Warren 2004), Tinker (Tinker et al. 2010)
       *
       *  @param method_SS method used to compute the power spectrum
       *  and &sigma;(mass); valid method_SS are: CAMB
       *  [http://camb.info/], classgal_v1 [http://class-code.net/],
       *  EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       *
       *  @param meanType meanType="mean_bias" \f$\rightarrow\f$ the
       *  effective bias is computed with Eq.(1);
       *  meanType="mean_pair_bias" \f$\rightarrow\f$ the effective
       *  bias is computed with Eq.(2)
       *
       *  @param output_root output_root of the parameter file used to
       *  compute the power spectrum and &sigma;(mass); it can be any
       *  name
       *
       *  @param Delta_crit \f$\Delta_{crit}\f$: the critical
       *  overdensity
       *  
       *  @param kk wave vector module
       *
       *  @param interpType method to interpolate the power spectrum
       *
       *  @param norm 0 \f$\rightarrow\f$ don't normalise the power
       *  spectrum; 1 \f$\rightarrow\f$ normalise the power spectrum;
       *  -1 \f$\rightarrow\f$ normalise only if sigma8 is set
       *
       *  @param k_min minimum wave vector module up to which the
       *  power spectrum is computed
       *
       *  @param k_max maximum wave vector module up to which the
       *  power spectrum is computed
       *
       *  @param prec accuracy of the integration 
       *
       *  @param input_file either the parameter file or the power
       *  spectrum file; if a parameter file is provided,
       *  i.e. input_file!=NULL and is_parameter_file=true, it will be
       *  used to compute the power spectrum; if a power spectrum file
       *  is provided, i.e. input_file!=NULL and
       *  is_parameter_file=false, then the provided power spectrum
       *  will be used directly; in both cases &sigma;<SUP>2</SUP>(M)
       *  is computed by integrating the computed/provided power
       *  spectrum ignoring the cosmological parameters of the object
       *
       *  @param is_parameter_file true \f$\rightarrow\f$ the input_file
       *  is a parameter file, used to compute the power spectrum with
       *  the method specified by method_Pk; false \f$\rightarrow\f$
       *  the input_file is a file containing the power spectrum
       *
       *  @param cosmology_mass cosmology used to measure the cluster
       *  masses
       *
       *  @param redshift_source vector containing the redshifts of
       *  the source galaxies, in case the cluster masses are
       *  estimated from weak lensing
       *
       *  @return a vector containing the mean and standard deviation
       *  of the effective dark matter bias
       */
      void generate_bias_eff_grid_one_cosmopar (std::vector<double> &parameter, std::vector<double> &bias_eff, const std::string dir_output, const std::string file_bias_eff_grid, const cbl::cosmology::CosmologicalParameter cosmoPar, const double min_par, const double max_par, const int nbin_par, const std::vector<double> mass, const std::vector<double> mass_grid, const std::vector<double> redshift, const std::string model_bias, const std::string method_SS, const std::string meanType="mean_bias", const std::string output_root="test", const double Delta_crit=200., const double kk=-1., const std::string interpType="Linear", const int norm=-1, const double k_min=0., const double k_max=100., const double prec=1.e-2, const std::string input_file=par::defaultString, const bool is_parameter_file=true, const cbl::cosmology::Cosmology cosmology_mass={}, const std::vector<double> redshift_source={});
      
      /**
       *  @brief effective bias of dark matter haloes, computed by
       *  weighting on the selection function on a grid of one input
       *  cosmological parameter; this function is used when modelling
       *  the two-point correlation function
       *
       *  this function computes the effective bias of dark matter
       *  haloes:
       *
       *  \f[ b_{eff}(z) = \frac{\int_{M_{min}}^{M_{max}} {\rm d}M\,
       *  b(M, z) \Phi(M, z) f(M, z)}{\int_{M_{min}}^{M_{max}} {\rm
       *  d}M\,\Phi(M, z) f(M, z)} \f]
       *
       *  where the linear bias of the \f$i\f$-th halo, \f$b^{i}(M,
       *  z)\f$, is computed by
       *  cbl::cosmology::Cosmology::bias_halo; the mass variance
       *  is computed by cbl::cosmology::Cosmology::sigma2M
       *
       *  @param parameter vector containing the grid of the
       *  cosmological parameters on which the effective bias grid is
       *  computed
       *
       *  @param bias_eff vector containing the effective bias grid
       *
       *  @param dir_output the directory where the effective bias
       *  grid is stored
       *
       *  @param file_bias_eff_grid the file there the effective bias
       *  grid is stored
       *
       *  @param cosmoPar the cosmological parameter for which the
       *  effective bias grid is computed
       *
       *  @param min_par the minimum value for the
       *  parameter where the effective bias is computed
       *  
       *  @param max_par the maximum value for the
       *  parameter where the effective bias is computed
       *
       *  @param nbin_par the number of points for the
       *  parameter where the effective bias is computed
       *
       *  @param redshift vector containing the redshifts; if it has
       *  size=1, it will be considered as the main redshift
       *
       *  @param Mass_min minimum cluster mass
       *
       *  @param Mass_max maximum cluster mass
       *
       *  @param model_bias author(s) who proposed the bias; valid
       *  authors are: ST99 (Sheth & Tormen 1999), SMT01 (Sheth, Mo &
       *  Tormen 2001), SMT01_WL04 (Sheth, Mo & Tormen 2001 with the
       *  correction of Warren 2004), Tinker (Tinker et al. 2010)
       * 
       *  @param model_MF author(s) who proposed the mass function;
       *  valid authors are: PS (Press&Schechter), ST (Sheth&Tormen),
       *  Jenkins (Jenkins et al. 2001), Warren (Warren et al. 2006),
       *  Reed, (Reed et al. 2007), Pan (Pan 2007), ShenH (halo MF by
       *  Shen et al. 2006), ShenF (filaments MF by Shen et al. 2006),
       *  ShenS (sheets MF by Shen et al. 2006), Tinker (Tinker et
       *  al. 2008), Crocce (Crocce et al. 2010), Angulo_FOF (FOF MF
       *  by Angulo et al. 2012), Angulo_Sub (SUBFIND MF by Angulo et
       *  al. 2012), Watson_FOF(FOF MF by Watson et al. 2012),
       *  Watson_SOH (MF for Spherical Overdensity Haloes by Watson et
       *  al. 2012), Manera (Manera et al. 2010), Bhattacharya
       *  (Bhattacharya et al. 2011), Courtin (Courtin et al. 2010),
       *  Peacock (by Peacock at al. 2007)
       *
       *  @param method_SS method used to compute the power spectrum
       *  and &sigma;(mass); valid method_SS are: CAMB
       *  [http://camb.info/], classgal_v1 [http://class-code.net/],
       *  EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       *
       *  @param alpha the \f$\alpha\f$ parameter of the cluster mass
       *  scaling relation
       *
       *  @param selection_function_file the input selection function
       *  file
       *
       *  @param column vector containing the columns with {mass,
       *  redshift, selection function}
       *
       *  @param output_root output_root of the parameter file used to
       *  compute the power spectrum and &sigma;(mass); it can be any
       *  name
       *
       *  @param Delta_crit \f$\Delta_{crit}\f$: the critical
       *  overdensity
       *  
       *  @param kk wave vector module
       *
       *  @param interpType method to interpolate the power spectrum
       *
       *  @param norm 0 \f$\rightarrow\f$ don't normalise the power
       *  spectrum; 1 \f$\rightarrow\f$ normalise the power spectrum;
       *  -1 \f$\rightarrow\f$ normalise only if sigma8 is set
       *
       *  @param k_min minimum wave vector module up to which the
       *  power spectrum is computed
       *
       *  @param k_max maximum wave vector module up to which the
       *  power spectrum is computed
       *
       *  @param prec accuracy of the integration 
       *
       *  @param input_file either the parameter file or the power
       *  spectrum file; if a parameter file is provided,
       *  i.e. input_file!=NULL and is_parameter_file=true, it will be
       *  used to compute the power spectrum; if a power spectrum file
       *  is provided, i.e. input_file!=NULL and
       *  is_parameter_file=false, then the provided power spectrum
       *  will be used directly; in both cases &sigma;<SUP>2</SUP>(M)
       *  is computed by integrating the computed/provided power
       *  spectrum ignoring the cosmological parameters of the object
       *
       *  @param is_parameter_file true \f$\rightarrow\f$ the input_file
       *  is a parameter file, used to compute the power spectrum with
       *  the method specified by method_Pk; false \f$\rightarrow\f$
       *  the input_file is a file containing the power spectrum
       *
       *  @return a vector containing the mean and standard deviation
       *  of the effective dark matter bias
       */
      void generate_bias_eff_grid_one_cosmopar (std::vector<double> &parameter, std::vector<double> &bias_eff, const std::string dir_output, const std::string file_bias_eff_grid, const cbl::cosmology::CosmologicalParameter cosmoPar, const double min_par, const double max_par, const int nbin_par, const double redshift, const double Mass_min, const double Mass_max, const std::string model_bias, const std::string model_MF, const std::string method_SS, const std::string selection_function_file, const std::vector<int> column={}, const double alpha=1., const std::string output_root="test", const double Delta_crit=200., const double kk=-1., const std::string interpType="Linear", const int norm=-1, const double k_min=0., const double k_max=100., const double prec=1.e-2, const std::string input_file=par::defaultString, const bool is_parameter_file=true);
      
      /**
       *  @brief effective bias of dark matter haloes, computed by
       *  averaging the bias of a set of haloes, interpolating the
       *  mass variance on a grid of masses and two input cosmological
       *  parameters; this function is used when modelling the
       *  two-point correlation function
       *
       *  this function computes the effective bias of dark matter
       *  haloes by either averaging the bias of a set of haloes with
       *  a given mass:
       *
       *  \f[b_{eff}(z) = \frac{1}{N_{halo}}\sum_{i=1}^{N_{halo}}
       *  b(M_i, z_i) \; , \; (1)\f]
       *
       *  or by averaging over halo pairs:
       *
       *  \f[b_{eff}(z) = \sqrt{ \frac{2}{N_{halo}(N_{halo}-1)}
       *  \sum_{i=1}^{N_{halo}}\sum_{j=i+1}^{N_{halo}} b(M_i,
       *  z_i)b(M_j, z_j)} \; , \; (2)\f]
       *
       *  where the linear bias of the \f$i\f$-th halo, \f$b^{i}(M,
       *  z)\f$, is computed by
       *  cbl::cosmology::Cosmology::bias_halo; the mass variance
       *  is computed by cbl::cosmology::Cosmology::sigma2M
       *
       *  @param parameter1 vector containing the grid of the first
       *  cosmological parameters on which the effective bias grid is
       *  computed
       *
       *  @param parameter2 vector containing the grid of the second
       *  cosmological parameters on which the effective bias grid is
       *  computed
       *
       *  @param bias_eff vector containing the effective bias grid
       *
       *  @param dir_output the directory where the effective bias
       *  grid is stored
       *
       *  @param file_bias_eff_grid the file there the effective bias
       *  grid is stored
       *
       *  @param cosmoPar1 the first cosmological parameter for which
       *  the effective bias grid is computed
       *
       *  @param min_par1 the minimum value for the first
       *  parameter where the effective bias is computed
       *  
       *  @param max_par1 the maximum value for the first
       *  parameter where the effective bias is computed
       *
       *  @param nbin_par1 the number of points for the first
       *  parameter where the effective bias is computed
       *
       *  @param cosmoPar2 the second cosmological parameter for which
       *  the effective bias grid is computed
       *
       *  @param min_par2 the minimum value for the second
       *  parameter where the effective bias is computed
       *  
       *  @param max_par2 the maximum value for the second
       *  parameter where the effective bias is computed
       *
       *  @param nbin_par2 the number of points for the second
       *  parameter where the effective bias is computed
       *
       *  @param mass vector containing the halo masses
       *
       *  @param mass_grid vector containing the halo masses on the
       *  grid used to interpolate the mass variance
       *
       *  @param redshift vector containing the redshifts; if it has
       *  size=1, it will be considered as the main redshift
       *
       *  @param model_bias author(s) who proposed the bias; valid
       *  authors are: ST99 (Sheth & Tormen 1999), SMT01 (Sheth, Mo &
       *  Tormen 2001), SMT01_WL04 (Sheth, Mo & Tormen 2001 with the
       *  correction of Warren 2004), Tinker (Tinker et al. 2010)
       *
       *  @param method_SS method used to compute the power spectrum
       *  and &sigma;(mass); valid method_SS are: CAMB
       *  [http://camb.info/], classgal_v1 [http://class-code.net/],
       *  EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       *
       *  @param meanType meanType="mean_bias" \f$\rightarrow\f$ the
       *  effective bias is computed with Eq.(1);
       *  meanType="mean_pair_bias" \f$\rightarrow\f$ the effective
       *  bias is computed with Eq.(2)
       *
       *  @param output_root output_root of the parameter file used to
       *  compute the power spectrum and &sigma;(mass); it can be any
       *  name
       *
       *  @param Delta \f$\Delta\f$: the overdensity, defined as the
       *  mean interior density relative to the background
       *  
       *  @param kk wave vector module
       *
       *  @param interpType method to interpolate the power spectrum
       *
       *  @param norm 0 \f$\rightarrow\f$ don't normalise the power
       *  spectrum; 1 \f$\rightarrow\f$ normalise the power spectrum;
       *  -1 \f$\rightarrow\f$ normalise only if sigma8 is set
       *
       *  @param k_min minimum wave vector module up to which the
       *  power spectrum is computed
       *
       *  @param k_max maximum wave vector module up to which the
       *  power spectrum is computed
       *
       *  @param prec accuracy of the integration 
       *
       *  @param input_file either the parameter file or the power
       *  spectrum file; if a parameter file is provided,
       *  i.e. input_file!=NULL and is_parameter_file=true, it will be
       *  used to compute the power spectrum; if a power spectrum file
       *  is provided, i.e. input_file!=NULL and
       *  is_parameter_file=false, then the provided power spectrum
       *  will be used directly; in both cases &sigma;<SUP>2</SUP>(M)
       *  is computed by integrating the computed/provided power
       *  spectrum ignoring the cosmological parameters of the object
       *
       *  @param is_parameter_file true \f$\rightarrow\f$ the input_file
       *  is a parameter file, used to compute the power spectrum with
       *  the method specified by method_Pk; false \f$\rightarrow\f$
       *  the input_file is a file containing the power spectrum
       *
       *  @param cosmology_mass cosmology used to measure the cluster
       *  masses
       *
       *  @param redshift_source vector containing the redshifts of
       *  the source galaxies, in case the cluster masses are estimated
       *  from weak lensing
       *
       *  @return a vector containing the mean and standard deviation
       *  of the effective dark matter bias
       */
      void generate_bias_eff_grid_two_cosmopars (std::vector<double> &parameter1, std::vector<double> &parameter2, std::vector<std::vector<double>> &bias_eff, const std::string dir_output, const std::string file_bias_eff_grid, const cbl::cosmology::CosmologicalParameter cosmoPar1, const double min_par1, const double max_par1, const int nbin_par1, const cbl::cosmology::CosmologicalParameter cosmoPar2, const double min_par2, const double max_par2, const int nbin_par2, const std::vector<double> mass, const std::vector<double> mass_grid, const std::vector<double> redshift, const std::string model_bias, const std::string method_SS, const std::string meanType="mean_bias", const std::string output_root="test", const double Delta=200., const double kk=-1., const std::string interpType="Linear", const int norm=-1, const double k_min=0., const double k_max=100., const double prec=1.e-2, const std::string input_file=par::defaultString, const bool is_parameter_file=true, const cbl::cosmology::Cosmology cosmology_mass={}, const std::vector<double> redshift_source={});

      /**
       *  @brief effective bias of dark matter haloes, computed using
       *  a given selection function; &sigma;(mass) and dln&sigma;/dM
       *  are provided in input
       *
       *  this function computes the effective bias of dark matter
       *  haloes:
       *
       *  \f[ b_{eff}(z) = \frac{\int_{M_{min}}^{M_{max}} {\rm d}M\,
       *  b(M, z) \Phi(M, z) f(M, z)}{\int_{M_{min}}^{M_{max}} {\rm
       *  d}M\,\Phi(M, z) f(M, z)} \f]
       *
       *  in the current implementation, the integral is actually
       *  replaced by the Riemann sum, as follows:
       *
       *  \f[ b_{eff}(z) \simeq \frac{\sum_{i_{min}}^{i_{max}} b(M, z)
       *  \Phi(M, z) f(M, z) (M_{i+1}-M_i)}{\sum_{i_{min}}^{i_{max}}
       *  \Phi(M, z) f(M, z) (M_{i+1}-M_i)} \f]
       *
       *  where the halo mass function, \f$\Phi(M, z)\f$, is computed
       *  by cbl::cosmology::Cosmology::mass_function, the linear
       *  bias, \f$b(M, z)\f$, is computed by
       *  cbl::cosmology::Cosmology::bias_halo, and f(M, z) is the
       *  selection function
       *
       *  @param interp_sigma FuncGrid object containing the values
       *  of &sigma;(mass) computed on a grid, used by interpolation
       *
       *  @param interp_DnSigma FuncGrid object containing the values
       *  of dln&sigma;/dM computed on a grid, used by interpolation
       *
       *  @param interp_SF FuncGrid object containing the values of
       *  the selection function computed on a grid in mass, at the
       *  mean redshift, used by interpolation
       *
       *  @param Mass_min minimum halo mass
       *
       *  @param Mass_max maximum halo mass
       *
       *  @param redshift vector containing the input redshifts
       *
       *  @param model_bias author(s) who proposed the bias; valid
       *  authors are: ST99 (Sheth & Tormen 1999), SMT01 (Sheth, Mo &
       *  Tormen 2001), SMT01_WL04 (Sheth, Mo & Tormen 2001 with the
       *  correction of Warren 2004), Tinker (Tinker et al. 2010)
       *
       *  @param model_MF author(s) who proposed the mass function;
       *  valid authors are: PS (Press&Schechter), ST (Sheth&Tormen),
       *  Jenkins (Jenkins et al. 2001), Warren (Warren et al. 2006),
       *  Reed, (Reed et al. 2007), Pan (Pan 2007), ShenH (halo MF by
       *  Shen et al. 2006), ShenF (filaments MF by Shen et al. 2006),
       *  ShenS (sheets MF by Shen et al. 2006), Tinker (Tinker et
       *  al. 2008), Crocce (Crocce et al. 2010), Angulo_FOF (FOF MF
       *  by Angulo et al. 2012), Angulo_Sub (SUBFIND MF by Angulo et
       *  al. 2012), Watson_FOF(FOF MF by Watson et al. 2012),
       *  Watson_SOH (MF for Spherical Overdensity Haloes by Watson et
       *  al. 2012), Manera (Manera et al. 2010), Bhattacharya
       *  (Bhattacharya et al. 2011), Courtin (Courtin et al. 2010),
       *  Peacock (by Peacock at al. 2007)
       *
       *  @param method_SS method used to compute the power spectrum
       *  and &sigma;(mass); valid method_SS are: CAMB
       *  [http://camb.info/], classgal_v1 [http://class-code.net/],
       *  EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       *
       *  @param alpha the \f$\alpha\f$ parameter of the cluster mass
       *  scaling relation
       *
       *  @param output_root output_root of the parameter file used to
       *  compute the power spectrum and &sigma;(mass); it can be any
       *  name
       *
       *  @param Delta_crit \f$\Delta_{crit}\f$: the critical
       *  overdensity
       *  
       *  @param kk wave vector module
       *
       *  @param interpType method to interpolate the power spectrum
       *
       *  @param norm 0 \f$\rightarrow\f$ don't normalise the power
       *  spectrum; 1 \f$\rightarrow\f$ normalise the power spectrum;
       *  -1 \f$\rightarrow\f$ normalise only if sigma8 is set
       *
       *  @param k_min minimum wave vector module up to which the
       *  power spectrum is computed
       *
       *  @param k_max maximum wave vector module up to which the
       *  power spectrum is computed
       *
       *  @param prec accuracy of the integration 
       *
       *  @param input_file either the parameter file or the power
       *  spectrum file; if a parameter file is provided,
       *  i.e. input_file!=NULL and is_parameter_file=true, it will be
       *  used to compute the power spectrum; if a power spectrum file
       *  is provided, i.e. input_file!=NULL and
       *  is_parameter_file=false, then the provided power spectrum
       *  will be used directly; in both cases &sigma;<SUP>2</SUP>(M)
       *  is computed by integrating the computed/provided power
       *  spectrum ignoring the cosmological parameters of the object
       *
       *  @param is_parameter_file true \f$\rightarrow\f$ the
       *  input_file is a parameter file, used to compute the power
       *  spectrum with the method specified by method_Pk; false
       *  \f$\rightarrow\f$ the input_file is a file containing the
       *  power spectrum
       *
       *  @return b<SUB>eff</SUB>: the effective dark matter bias
       *
       *  @warning interp_sigma and interp_DnSigma have to be
       *  computed at the same cosmology of the object. They are
       *  provided as an input just to improve the performances in
       *  some applications (e.g. MCMC) where these quantities can be
       *  computed once
       */
      std::vector<double> bias_eff_selection_function (const glob::FuncGrid interp_sigma, const glob::FuncGrid interp_DnSigma, const glob::FuncGrid interp_SF, const double Mass_min, const double Mass_max, const std::vector<double> redshift, const std::string model_bias, const std::string model_MF, const std::string method_SS, const double alpha=1., const std::string output_root="test", const double Delta_crit=200., const double kk=-1., const std::string interpType="Linear", const int norm=-1, const double k_min=0., const double k_max=100., const double prec=1.e-2, const std::string input_file=par::defaultString, const bool is_parameter_file=true);
       
      /**
       *  @brief effective bias of dark matter haloes, computed using
       *  a given selection function; &sigma;(mass) and dln&sigma;/dM
       *  are provided in input
       *
       *  this function computes the effective bias of dark matter
       *  haloes:
       *
       *  \f[ b_{eff}(z) = \frac{\int_{M_{min}}^{M_{max}} {\rm d}M\,
       *  b(M, z) \Phi(M, z) f(M, z)}{\int_{M_{min}}^{M_{max}} {\rm
       *  d}M\,\Phi(M, z) f(M, z)} \f]
       *
       *  in the current implementation, the integral is actually
       *  replaced by the Riemann sum, as follows:
       *
       *  \f[ b_{eff}(z) \simeq \frac{\sum_{i_{min}}^{i_{max}} b(M, z)
       *  \Phi(M, z) f(M, z) (M_{i+1}-M_i)}{\sum_{i_{min}}^{i_{max}}
       *  \Phi(M, z) f(M, z) (M_{i+1}-M_i)} \f]
       *
       *  where the halo mass function, \f$\Phi(M, z)\f$, is computed
       *  by cbl::cosmology::Cosmology::mass_function, the linear
       *  bias, \f$b(M, z)\f$, is computed by
       *  cbl::cosmology::Cosmology::bias_halo, and f(M, z) is the
       *  selection function
       *
       *  @param interp_sigma FuncGrid object containing the values
       *  of &sigma;(mass) computed on a grid, used by interpolation
       *
       *  @param interp_DnSigma FuncGrid object containing the values
       *  of dln&sigma;/dM computed on a grid, used by interpolation
       *
       *  @param interp_SF FuncGrid2D object containing the values of
       *  the selection function computed on a grid in mass and
       *  redshift, used by interpolation
       *
       *  @param Mass_min minimum halo mass
       *
       *  @param Mass_max maximum halo mass
       *
       *  @param redshift vector containing the input redshifts
       *
       *  @param model_bias author(s) who proposed the bias; valid
       *  authors are: ST99 (Sheth & Tormen 1999), SMT01 (Sheth, Mo &
       *  Tormen 2001), SMT01_WL04 (Sheth, Mo & Tormen 2001 with the
       *  correction of Warren 2004), Tinker (Tinker et al. 2010)
       *
       *  @param model_MF author(s) who proposed the mass function;
       *  valid authors are: PS (Press&Schechter), ST (Sheth&Tormen),
       *  Jenkins (Jenkins et al. 2001), Warren (Warren et al. 2006),
       *  Reed, (Reed et al. 2007), Pan (Pan 2007), ShenH (halo MF by
       *  Shen et al. 2006), ShenF (filaments MF by Shen et al. 2006),
       *  ShenS (sheets MF by Shen et al. 2006), Tinker (Tinker et
       *  al. 2008), Crocce (Crocce et al. 2010), Angulo_FOF (FOF MF
       *  by Angulo et al. 2012), Angulo_Sub (SUBFIND MF by Angulo et
       *  al. 2012), Watson_FOF(FOF MF by Watson et al. 2012),
       *  Watson_SOH (MF for Spherical Overdensity Haloes by Watson et
       *  al. 2012), Manera (Manera et al. 2010), Bhattacharya
       *  (Bhattacharya et al. 2011), Courtin (Courtin et al. 2010),
       *  Peacock (by Peacock at al. 2007)
       *
       *  @param method_SS method used to compute the power spectrum
       *  and &sigma;(mass); valid method_SS are: CAMB
       *  [http://camb.info/], classgal_v1 [http://class-code.net/],
       *  EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       *
       *  @param alpha the \f$\alpha\f$ parameter of the cluster mass
       *  scaling relation
       *
       *  @param output_root output_root of the parameter file used to
       *  compute the power spectrum and &sigma;(mass); it can be any
       *  name
       *
       *  @param Delta_crit \f$\Delta_{crit}\f$: the critical
       *  overdensity
       *  
       *  @param kk wave vector module
       *
       *  @param interpType method to interpolate the power spectrum
       *
       *  @param norm 0 \f$\rightarrow\f$ don't normalise the power
       *  spectrum; 1 \f$\rightarrow\f$ normalise the power spectrum;
       *  -1 \f$\rightarrow\f$ normalise only if sigma8 is set
       *
       *  @param k_min minimum wave vector module up to which the
       *  power spectrum is computed
       *
       *  @param k_max maximum wave vector module up to which the
       *  power spectrum is computed
       *
       *  @param prec accuracy of the integration 
       *
       *  @param input_file either the parameter file or the power
       *  spectrum file; if a parameter file is provided,
       *  i.e. input_file!=NULL and is_parameter_file=true, it will be
       *  used to compute the power spectrum; if a power spectrum file
       *  is provided, i.e. input_file!=NULL and
       *  is_parameter_file=false, then the provided power spectrum
       *  will be used directly; in both cases &sigma;<SUP>2</SUP>(M)
       *  is computed by integrating the computed/provided power
       *  spectrum ignoring the cosmological parameters of the object
       *
       *  @param is_parameter_file true \f$\rightarrow\f$ the
       *  input_file is a parameter file, used to compute the power
       *  spectrum with the method specified by method_Pk; false
       *  \f$\rightarrow\f$ the input_file is a file containing the
       *  power spectrum
       *
       *  @return b<SUB>eff</SUB>: the effective dark matter bias
       *
       *  @warning interp_sigma and interp_DnSigma have to be
       *  computed at the same cosmology of the object. They are
       *  provided as an input just to improve the performances in
       *  some applications (e.g. MCMC) where these quantities can be
       *  computed once
       */
      std::vector<double> bias_eff_selection_function (const glob::FuncGrid interp_sigma, const glob::FuncGrid interp_DnSigma, const glob::FuncGrid2D interp_SF, const double Mass_min, const double Mass_max, const std::vector<double> redshift, const std::string model_bias, const std::string model_MF, const std::string method_SS, const double alpha=1., const std::string output_root="test", const double Delta_crit=200., const double kk=-1., const std::string interpType="Linear", const int norm=-1, const double k_min=0., const double k_max=100., const double prec=1.e-2, const std::string input_file=par::defaultString, const bool is_parameter_file=true);

      /**
       *  @brief effective bias of dark matter haloes, computed using
       *  a given selection function
       *
       *  this function computes the effective bias of dark matter
       *  haloes:
       *
       *  \f[ b_{eff}(z) = \frac{\int_{M_{min}}^{M_{max}} {\rm d}M\,
       *  b(M, z) \Phi(M, z) f(M, z)}{\int_{M_{min}}^{M_{max}} {\rm
       *  d}M\,\Phi(M, z) f(M, z)} \f]
       *
       *  in the current implementation, the integral is actually
       *  replaced by the Riemann sum, as follows:
       *
       *  \f[ b_{eff}(z) \simeq \frac{\sum_{i_{min}}^{i_{max}} b(M, z)
       *  \Phi(M, z) f(M, z) (M_{i+1}-M_i)}{\sum_{i_{min}}^{i_{max}}
       *  \Phi(M, z) f(M, z) (M_{i+1}-M_i)} \f]
       *
       *  where the halo mass function, \f$\Phi(M, z)\f$, is computed
       *  by cbl::cosmology::Cosmology::mass_function, the linear
       *  bias, \f$b(M, z)\f$, is computed by
       *  cbl::cosmology::Cosmology::bias_halo, and f(M, z) is the
       *  selection function
       *
       *  @param Mass_min minimum halo mass
       *
       *  @param Mass_max maximum halo mass
       *
       *  @param redshift vector containing the input redshifts
       *
       *  @param model_bias author(s) who proposed the bias; valid
       *  authors are: ST99 (Sheth & Tormen 1999), SMT01 (Sheth, Mo &
       *  Tormen 2001), SMT01_WL04 (Sheth, Mo & Tormen 2001 with the
       *  correction of Warren 2004), Tinker (Tinker et al. 2010)
       *
       *  @param model_MF author(s) who proposed the mass function;
       *  valid authors are: PS (Press&Schechter), ST (Sheth&Tormen),
       *  Jenkins (Jenkins et al. 2001), Warren (Warren et al. 2006),
       *  Reed, (Reed et al. 2007), Pan (Pan 2007), ShenH (halo MF by
       *  Shen et al. 2006), ShenF (filaments MF by Shen et al. 2006),
       *  ShenS (sheets MF by Shen et al. 2006), Tinker (Tinker et
       *  al. 2008), Crocce (Crocce et al. 2010), Angulo_FOF (FOF MF
       *  by Angulo et al. 2012), Angulo_Sub (SUBFIND MF by Angulo et
       *  al. 2012), Watson_FOF(FOF MF by Watson et al. 2012),
       *  Watson_SOH (MF for Spherical Overdensity Haloes by Watson et
       *  al. 2012), Manera (Manera et al. 2010), Bhattacharya
       *  (Bhattacharya et al. 2011), Courtin (Courtin et al. 2010),
       *  Peacock (by Peacock at al. 2007)
       *
       *  @param method_SS method used to compute the power spectrum
       *  and &sigma;(mass); valid method_SS are: CAMB
       *  [http://camb.info/], classgal_v1 [http://class-code.net/],
       *  EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       *
       *  @param selection_function_file input file with the selection
       *  functon
       *
       *  @param column vector containing the three columns of the
       *  selection function file to be read
       *
       *  @param alpha the \f$\alpha\f$ parameter of the cluster mass
       *  scaling relation
       *
       *  @param output_root output_root of the parameter file used to
       *  compute the power spectrum and &sigma;(mass); it can be any
       *  name
       *
       *  @param Delta_crit \f$\Delta_{crit}\f$: the critical
       *  overdensity
       *  
       *  @param kk wave vector module
       *
       *  @param interpType method to interpolate the power spectrum
       *
       *  @param norm 0 \f$\rightarrow\f$ don't normalise the power
       *  spectrum; 1 \f$\rightarrow\f$ normalise the power spectrum;
       *  -1 \f$\rightarrow\f$ normalise only if sigma8 is set
       *
       *  @param k_min minimum wave vector module up to which the
       *  power spectrum is computed
       *
       *  @param k_max maximum wave vector module up to which the
       *  power spectrum is computed
       *
       *  @param prec accuracy of the integration 
       *
       *  @param input_file either the parameter file or the power
       *  spectrum file; if a parameter file is provided,
       *  i.e. input_file!=NULL and is_parameter_file=true, it will be
       *  used to compute the power spectrum; if a power spectrum file
       *  is provided, i.e. input_file!=NULL and
       *  is_parameter_file=false, then the provided power spectrum
       *  will be used directly; in both cases &sigma;<SUP>2</SUP>(M)
       *  is computed by integrating the computed/provided power
       *  spectrum ignoring the cosmological parameters of the object
       *
       *  @param is_parameter_file true \f$\rightarrow\f$ the
       *  input_file is a parameter file, used to compute the power
       *  spectrum with the method specified by method_Pk; false
       *  \f$\rightarrow\f$ the input_file is a file containing the
       *  power spectrum
       *
       *  @return b<SUB>eff</SUB>: the effective dark matter bias
       */
      std::vector<double> bias_eff_selection_function (const double Mass_min, const double Mass_max, const std::vector<double> redshift, const std::string model_bias, const std::string model_MF, const std::string method_SS, const std::string selection_function_file, const std::vector<int> column={}, const double alpha=1., const std::string output_root="test", const double Delta_crit=200., const double kk=-1., const std::string interpType="Linear", const int norm=-1, const double k_min=0., const double k_max=100., const double prec=1.e-2, const std::string input_file=par::defaultString, const bool is_parameter_file=true);
       
      ///@}


      /**
       *  @name Functions to model redshift-space distortions
       */
      ///@{

      /**
       *  @brief the linear growth rate at a given redshift,
       *  \f$f(z)\f$
       *
       *  this function computes the following function:
       *
       *  \f[ f(z) = \frac{{\rm d}\,\ln D}{{\rm d}\,\ln a} \f]
       *
       *  using approximated functions provided by Wang & Steinhardt
       *  1998, Kiakotou, Elgarøy & Lahav 2008, Gong et al. 2009
       *
       *  @param redshift the redshift
       *  @param kk wave vector module
       *  @return the linear growth rate
       *
       *  @warning the current implementation is not correct if w_a is
       *  different than 0
       */
      double linear_growth_rate (const double redshift, const double kk=-1.) const;

      /**
       *  @brief f*&sigma;<SUB>8</SUB>: the linear growth rate times
       *  the dark matter rms mass fluctuation within 8 Mpc/h
       *
       *  @param redshift the redshift
       *
       *  @param method_Pk method used to compute the power spectrum and
       *  &sigma;(mass); valid choices for method_Pk are: CAMB [http://camb.info/],
       *  classgal_v1 [http://class-code.net/], MPTbreeze-v1
       *  [http://arxiv.org/abs/1207.1465], EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       *  
       *  @param output_root output_root of the parameter file used to compute
       *  the power spectrum and &sigma;(mass); it can be any name
       *
       *  @param kk wave vector module 
       *
       *  @param NL 0 \f$\rightarrow\f$ linear power spectrum; 1 \f$\rightarrow\f$ non-linear power
       *  spectrum
       *
       *  @param k_min minimum wave vector module up to which the power
       *  spectrum is computed
       *
       *  @param k_max maximum wave vector module up to which the power
       *  spectrum is computed
       *   
       *  @param prec accuracy of the integration
       *
       *  @param file_par name of the parameter file; if a
       *  parameter file is provided (i.e. file_par!=NULL), it will be
       *  used, ignoring the cosmological parameters of the object
       *
       *  @return f*&sigma;<SUB>8</SUB>
       */
      double fsigma8 (const double redshift, const std::string method_Pk, const std::string output_root="test", const double kk=1., const bool NL=0, const double k_min=0., const double k_max=100., const double prec=1.e-2, const std::string file_par=par::defaultString) const;

      /**
       *  @brief the specific growth rate &beta;
       *  @param redshift the redshift
       *  @param bias bias
       *  @param kk wave vector module
       *  @return &beta;=f/b, where f is the linear growth rate and b is
       *  the bias
       */
      double beta (const double redshift, const double bias, const double kk=-1.) const;

      /**
       *  @brief the error on the specific growth rate &beta;
       *  @param redshift the redshift
       *  @param bias bias
       *  @param err_bias error on the bias
       *  @param kk wave vector module
       *  @return error on &beta;=f/b, where f is the linear growth rate
       *  and b is the bias
       */
      double error_beta (const double redshift, const double bias, const double err_bias, const double kk=-1.) const;

      /**
       *  @brief the error on the specific growth rate &beta;
       *
       *  @param Mass_min minimum halo mass
       * 
       *  @param Mass_max maximum halo mass
       *
       *  @param redshift the redshift
       *
       *  @param model_bias author(s) who proposed the bias; valid
       *  authors are: ST99 (Sheth & Tormen 1999), SMT01 (Sheth, Mo &
       *  Tormen 2001), SMT01_WL04 (Sheth, Mo & Tormen 2001 with the
       *  correction of Warren 2004), Tinker (Tinker et al. 2010)
       *
       *  @param model_MF author(s) who proposed the mass function;
       *  valid authors are: PS (Press&Schechter), ST (Sheth&Tormen),
       *  Jenkins (Jenkins et al. 2001), Warren (Warren et al. 2006),
       *  Reed, (Reed et al. 2007), Pan (Pan 2007), ShenH (halo MF by
       *  Shen et al. 2006), ShenF (filaments MF by Shen et al. 2006),
       *  ShenS (sheets MF by Shen et al. 2006), Tinker (Tinker et
       *  al. 2008), Crocce (Crocce et al. 2010), Angulo_FOF (FOF MF by
       *  Angulo et al. 2012), Angulo_Sub (SUBFIND MF by Angulo et
       *  al. 2012), Watson_FOF(FOF MF by Watson et al. 2012),
       *  Watson_SOH (MF for Spherical Overdensity Haloes by Watson 
       *  et al. 2012), Manera (Manera et al. 2010), Bhattacharya 
       *  (Bhattacharya et al. 2011), Courtin (Courtin et al. 2010),
       *  Peacock (by Peacock at al. 2007)
       *
       *  @param method_SS method used to compute the power spectrum
       *  and &sigma;(mass); valid method_SS are: CAMB
       *  [http://camb.info/], classgal_v1 [http://class-code.net/],
       *  EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       *
       *  @param output_root output_root of the parameter file used to
       *  compute the power spectrum and &sigma;(mass); it can be any
       *  name
       *
       *  @param Delta \f$\Delta\f$: the overdensity, defined as the
       *  mean interior density relative to the background
       *  
       *  @param kk wave vector module
       *
       *  @param interpType method to interpolate the power spectrum
       *
       *  @param norm 0 \f$\rightarrow\f$ don't normalise the power
       *  spectrum; 1 \f$\rightarrow\f$ normalise the power spectrum;
       *  -1 \f$\rightarrow\f$ normalise only if sigma8 is set
       *
       *  @param k_min minimum wave vector module up to which the
       *  power spectrum is computed
       *
       *  @param k_max maximum wave vector module up to which the
       *  power spectrum is computed
       *
       *  @param prec accuracy of the integration 
       *
       *  @param input_file either the parameter file or the power
       *  spectrum file; if a parameter file is provided,
       *  i.e. input_file!=NULL and is_parameter_file=true, it will be
       *  used to compute the power spectrum; if a power spectrum file
       *  is provided, i.e. input_file!=NULL and
       *  is_parameter_file=false, then the provided power spectrum
       *  will be used directly; in both cases &sigma;<SUP>2</SUP>(M)
       *  is computed by integrating the computed/provided power
       *  spectrum ignoring the cosmological parameters of the object
       *
       *  @param is_parameter_file true \f$\rightarrow\f$ the input_file
       *  is a parameter file, used to compute the power spectrum with
       *  the method specified by method_Pk; false \f$\rightarrow\f$
       *  the input_file is a file containing the power spectrum
       *
       *  @return &beta;=f/b, where f is the linear growth rate and b is
       *  the bias
       *
       *  @warning the input parameter \f$\Delta\f$ is the background
       *  overdensity, not the critical overdensity
       *  \f$\Delta_{crit}\f$; the function
       *  cbl::Cosmology::Delta_vir can be used to convert
       *  \f$\Delta_{crit}\f$ into \f$\Delta\f$
       *
       *  @warning the mass function by Manera et al. (2010) has been
       *  tested only for z=0 and z=0.5; the mass function by Despali
       *  et al. (2016) is currently implemented only for virial
       *  masses and at \f$z<1.25\f$
       */
      double beta (const double Mass_min, const double Mass_max, const double redshift, const std::string model_bias, const std::string model_MF, const std::string method_SS, const std::string output_root="test", const double Delta=200., const double kk=-1., const std::string interpType="Linear", const int norm=-1, const double k_min=0., const double k_max=100., const double prec=1.e-2, const std::string input_file=par::defaultString, const bool is_parameter_file=true);

      /**
       *  @brief the specific growth rate &beta;
       *
       *  @param Mass_min minimum halo mass
       * 
       *  @param Mass_max maximum halo mass
       *
       *  @param redshift the redshift
       *
       *  @param model_bias author(s) who proposed the bias; valid
       *  authors are: ST99 (Sheth & Tormen 1999), SMT01 (Sheth, Mo &
       *  Tormen 2001), SMT01_WL04 (Sheth, Mo & Tormen 2001 with the
       *  correction of Warren 2004), Tinker (Tinker et al. 2010)
       *
       *  @param model_MF author(s) who proposed the mass function;
       *  valid authors are: PS (Press&Schechter), ST (Sheth&Tormen),
       *  Jenkins (Jenkins et al. 2001), Warren (Warren et al. 2006),
       *  Reed, (Reed et al. 2007), Pan (Pan 2007), ShenH (halo MF by
       *  Shen et al. 2006), ShenF (filaments MF by Shen et al. 2006),
       *  ShenS (sheets MF by Shen et al. 2006), Tinker (Tinker et
       *  al. 2008), Crocce (Crocce et al. 2010), Angulo_FOF (FOF MF by
       *  Angulo et al. 2012), Angulo_Sub (SUBFIND MF by Angulo et
       *  al. 2012), Watson_FOF(FOF MF by Watson et al. 2012),
       *  Watson_SOH (MF for Spherical Overdensity Haloes by Watson 
       *  et al. 2012), Manera (Manera et al. 2010), Bhattacharya 
       *  (Bhattacharya et al. 2011), Courtin (Courtin et al. 2010),
       *  Peacock (by Peacock at al. 2007)
       *
       *  @param method_SS method used to compute the power spectrum and
       *  &sigma;(mass); valid method_SS are: CAMB [http://camb.info/],
       *  classgal_v1 [http://class-code.net/], EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       *
       *  @param err_bias error on the bias
       *
       *  @param output_root output_root of the parameter file used to
       *  compute the power spectrum and &sigma;(mass); it can be any
       *  name
       *
       *  @param Delta \f$\Delta\f$: the overdensity, defined as the
       *  mean interior density relative to the background
       *  
       *  @param kk wave vector module
       *
       *  @param interpType method to interpolate the power spectrum
       *
       *  @param norm 0 \f$\rightarrow\f$ don't normalise the power
       *  spectrum; 1 \f$\rightarrow\f$ normalise the power spectrum;
       *  -1 \f$\rightarrow\f$ normalise only if sigma8 is set
       *
       *  @param k_min minimum wave vector module up to which the
       *  power spectrum is computed
       *
       *  @param k_max maximum wave vector module up to which the
       *  power spectrum is computed
       *
       *  @param prec accuracy of the integration 
       *
       *  @param input_file either the parameter file or the power
       *  spectrum file; if a parameter file is provided,
       *  i.e. input_file!=NULL and is_parameter_file=true, it will be
       *  used to compute the power spectrum; if a power spectrum file
       *  is provided, i.e. input_file!=NULL and
       *  is_parameter_file=false, then the provided power spectrum
       *  will be used directly; in both cases &sigma;<SUP>2</SUP>(M)
       *  is computed by integrating the computed/provided power
       *  spectrum ignoring the cosmological parameters of the object
       *
       *  @param is_parameter_file true \f$\rightarrow\f$ the input_file
       *  is a parameter file, used to compute the power spectrum with
       *  the method specified by method_Pk; false \f$\rightarrow\f$
       *  the input_file is a file containing the power spectrum
       *
       *  @return error on &beta;=f/b, where f is the linear growth
       *  rate and b is the bias
       *
       *  @warning the input parameter \f$\Delta\f$ is the background
       *  overdensity, not the critical overdensity
       *  \f$\Delta_{crit}\f$; the function
       *  cbl::Cosmology::Delta_vir can be used to convert
       *  \f$\Delta_{crit}\f$ into \f$\Delta\f$
       *
       *  @warning the mass function by Manera et al. (2010) has been
       *  tested only for z=0 and z=0.5; the mass function by Despali
       *  et al. (2016) is currently implemented only for virial
       *  masses and at \f$z<1.25\f$
       */
      double error_beta (const double Mass_min, const double Mass_max, const double redshift, const std::string model_bias, const std::string model_MF, const std::string method_SS, const double err_bias, const std::string output_root="test", const double Delta=200., const double kk=-1., const std::string interpType="Linear", const int norm=-1, const double k_min=0., const double k_max=100., const double prec=1.e-2, const std::string input_file=par::defaultString, const bool is_parameter_file=true); 
  
      /**
       *  @brief the specific growth rate &beta;
       *
       *  @param MM vector of halo masses
       *
       *  @param MF vector of mass function values, d&Phi;/dM=dn(M)/dM
       *
       *  @param redshift the redshift
       *
       *  @param model_bias author(s) who proposed the bias; valid
       *  authors are: ST99 (Sheth & Tormen 1999), SMT01 (Sheth, Mo &
       *  Tormen 2001), SMT01_WL04 (Sheth, Mo & Tormen 2001 with the
       *  correction of Warren 2004), Tinker (Tinker et al. 2010)
       *
       *  @param method_SS method used to compute the power spectrum
       *  and &sigma;(mass); valid method_SS are: CAMB
       *  [http://camb.info/], classgal_v1 [http://class-code.net/],
       *  EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       *
       *  @param output_root output_root of the parameter file used to
       *  compute the power spectrum and &sigma;(mass); it can be any
       *  name
       *
       *  @param Delta \f$\Delta\f$: the overdensity, defined as the
       *  mean interior density relative to the background
       *  
       *  @param kk wave vector module
       *
       *  @param interpType method to interpolate the power spectrum
       *
       *  @param norm 0 \f$\rightarrow\f$ don't normalise the power
       *  spectrum; 1 \f$\rightarrow\f$ normalise the power spectrum;
       *  -1 \f$\rightarrow\f$ normalise only if sigma8 is set
       *
       *  @param k_min minimum wave vector module up to which the
       *  power spectrum is computed
       *
       *  @param k_max maximum wave vector module up to which the
       *  power spectrum is computed
       *
       *  @param prec accuracy of the integration 
       *      
       *  @param input_file either the parameter file or the power
       *  spectrum file; if a parameter file is provided,
       *  i.e. input_file!=NULL and is_parameter_file=true, it will be
       *  used to compute the power spectrum; if a power spectrum file
       *  is provided, i.e. input_file!=NULL and
       *  is_parameter_file=false, then the provided power spectrum
       *  will be used directly; in both cases &sigma;<SUP>2</SUP>(M)
       *  is computed by integrating the computed/provided power
       *  spectrum ignoring the cosmological parameters of the object
       *
       *  @param is_parameter_file true \f$\rightarrow\f$ the input_file
       *  is a parameter file, used to compute the power spectrum with
       *  the method specified by method_Pk; false \f$\rightarrow\f$
       *  the input_file is a file containing the power spectrum
       *
       *  @return &beta;=f/b, where f is the linear growth rate and b is
       *  the bias
       *
       *  @warning the input parameter \f$\Delta\f$ is the background
       *  overdensity, not the critical overdensity
       *  \f$\Delta_{crit}\f$; the function
       *  cbl::Cosmology::Delta_vir can be used to convert
       *  \f$\Delta_{crit}\f$ into \f$\Delta\f$
       *
       *  @warning the mass function by Manera et al. (2010) has been
       *  tested only for z=0 and z=0.5; the mass function by Despali
       *  et al. (2016) is currently implemented only for virial
       *  masses and at \f$z<1.25\f$
       */
      double beta (const std::vector<double> MM, const std::vector<double> MF, const double redshift, const std::string model_bias, const std::string method_SS, const std::string output_root="test", const double Delta=200., const double kk=-1., const std::string interpType="Linear", const int norm=-1, const double k_min=0., const double k_max=100., const double prec=1.e-2, const std::string input_file=par::defaultString, const bool is_parameter_file=true);

      /**
       *  @brief the error on the specific growth rate &beta;
       *
       *  @param MM vector of halo masses
       *
       *  @param MF vector of mass function values, d&Phi;/dM=dn(M)/dM
       *
       *  @param redshift the redshift
       *
       *  @param model_bias author(s) who proposed the bias; valid authors
       *  are: ST99 (Sheth & Tormen 1999), SMT01 (Sheth, Mo & Tormen
       *  2001), SMT01_WL04 (Sheth, Mo & Tormen 2001 with the correction
       *  of Warren 2004), Tinker (Tinker et al. 2010)
       *
       *  @param method_SS method used to compute the power spectrum and
       *  &sigma;(mass); valid method_SS are: CAMB [http://camb.info/],
       *  classgal_v1 [http://class-code.net/], EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       *
       *  @param err_bias error on the bias
       *
       *  @param output_root output_root of the parameter file used to
       *  compute the power spectrum and &sigma;(mass); it can be any
       *  name
       *
       *  @param Delta \f$\Delta\f$: the overdensity, defined as the
       *  mean interior density relative to the background
       *  
       *  @param kk wave vector module
       *
       *  @param interpType method to interpolate the power spectrum
       *
       *  @param norm 0 \f$\rightarrow\f$ don't normalise the power
       *  spectrum; 1 \f$\rightarrow\f$ normalise the power spectrum;
       *  -1 \f$\rightarrow\f$ normalise only if sigma8 is set
       *
       *  @param k_min minimum wave vector module up to which the
       *  power spectrum is computed
       *
       *  @param k_max maximum wave vector module up to which the
       *  power spectrum is computed
       *
       *  @param prec accuracy of the integration 
       *
       *  @param input_file either the parameter file or the power
       *  spectrum file; if a parameter file is provided,
       *  i.e. input_file!=NULL and is_parameter_file=true, it will be
       *  used to compute the power spectrum; if a power spectrum file
       *  is provided, i.e. input_file!=NULL and
       *  is_parameter_file=false, then the provided power spectrum
       *  will be used directly; in both cases &sigma;<SUP>2</SUP>(M)
       *  is computed by integrating the computed/provided power
       *  spectrum ignoring the cosmological parameters of the object
       *
       *  @param is_parameter_file true \f$\rightarrow\f$ the input_file
       *  is a parameter file, used to compute the power spectrum with
       *  the method specified by method_Pk; false \f$\rightarrow\f$
       *  the input_file is a file containing the power spectrum
       *
       *  @return error on &beta;=f/b, where f is the linear growth
       *  rate and b is the bias
       *
       *  @warning the input parameter \f$\Delta\f$ is the background
       *  overdensity, not the critical overdensity
       *  \f$\Delta_{crit}\f$; the function
       *  cbl::Cosmology::Delta_vir can be used to convert
       *  \f$\Delta_{crit}\f$ into \f$\Delta\f$
       *
       *  @warning the mass function by Manera et al. (2010) has been
       *  tested only for z=0 and z=0.5; the mass function by Despali
       *  et al. (2016) is currently implemented only for virial
       *  masses and at \f$z<1.25\f$
       */
      double error_beta (const std::vector<double> MM, const std::vector<double> MF, const double redshift, const std::string model_bias, const std::string method_SS, const double err_bias, const std::string output_root="test", const double Delta=200., const double kk=-1., const std::string interpType="Linear", const int norm=-1, const double k_min=0., const double k_max=100., const double prec=1.e-2, const std::string input_file=par::defaultString, const bool is_parameter_file=true);
 
      /**
       *  @brief the error on the specific growth rate &beta; from
       * Bianchi et al. 2012
       *
       *  @param Volume comoving volume 
       *
       *  @param density comoving density
       *
       *  @param Mass_min minimum halo mass
       * 
       *  @param Mass_max maximum halo mass
       *
       *  @param redshift the redshift
       *
       *  @param model_bias author(s) who proposed the bias; valid
       *  authors are: ST99 (Sheth & Tormen 1999), SMT01 (Sheth, Mo &
       *  Tormen 2001), SMT01_WL04 (Sheth, Mo & Tormen 2001 with the
       *  correction of Warren 2004), Tinker (Tinker et al. 2010)
       *
       *  @param model_MF author(s) who proposed the mass function;
       *  valid authors are: PS (Press&Schechter), ST (Sheth&Tormen),
       *  Jenkins (Jenkins et al. 2001), Warren (Warren et al. 2006),
       *  Reed, (Reed et al. 2007), Pan (Pan 2007), ShenH (halo MF by
       *  Shen et al. 2006), ShenF (filaments MF by Shen et al. 2006),
       *  ShenS (sheets MF by Shen et al. 2006), Tinker (Tinker et
       *  al. 2008), Crocce (Crocce et al. 2010), Angulo_FOF (FOF MF by
       *  Angulo et al. 2012), Angulo_Sub (SUBFIND MF by Angulo et
       *  al. 2012), Watson_FOF(FOF MF by Watson et al. 2012),
       *  Watson_SOH (MF for Spherical Overdensity Haloes by Watson 
       *  et al. 2012), Manera (Manera et al. 2010), Bhattacharya 
       *  (Bhattacharya et al. 2011), Courtin (Courtin et al. 2010),
       *  Peacock (by Peacock at al. 2007)
       *
       *  @param method_SS method used to compute the power spectrum
       *  and &sigma;(mass); valid method_SS are: CAMB
       *  [http://camb.info/], classgal_v1 [http://class-code.net/],
       *  EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       *
       *  @param output_root output_root of the parameter file used to
       *  compute the power spectrum and &sigma;(mass); it can be any
       *  name
       *
       *  @param Delta \f$\Delta\f$: the overdensity, defined as the
       *  mean interior density relative to the background
       *  
       *  @param kk wave vector module
       *
       *  @param interpType method to interpolate the power spectrum
       *
       *  @param norm 0 \f$\rightarrow\f$ don't normalise the power
       *  spectrum; 1 \f$\rightarrow\f$ normalise the power spectrum;
       *  -1 \f$\rightarrow\f$ normalise only if sigma8 is set
       *
       *  @param k_min minimum wave vector module up to which the
       *  power spectrum is computed
       *
       *  @param k_max maximum wave vector module up to which the
       *  power spectrum is computed
       *
       *  @param prec accuracy of the integration 
       *
       *  @param input_file either the parameter file or the power
       *  spectrum file; if a parameter file is provided,
       *  i.e. input_file!=NULL and is_parameter_file=true, it will be
       *  used to compute the power spectrum; if a power spectrum file
       *  is provided, i.e. input_file!=NULL and
       *  is_parameter_file=false, then the provided power spectrum
       *  will be used directly; in both cases &sigma;<SUP>2</SUP>(M)
       *  is computed by integrating the computed/provided power
       *  spectrum ignoring the cosmological parameters of the object
       *
       *  @param is_parameter_file true \f$\rightarrow\f$ the input_file
       *  is a parameter file, used to compute the power spectrum with
       *  the method specified by method_Pk; false \f$\rightarrow\f$
       *  the input_file is a file containing the power spectrum
       *
       *  @return error on &beta;=f/b, where f is the linear growth rate
       *  and b is the bias
       *
       *  @warning the input parameter \f$\Delta\f$ is the background
       *  overdensity, not the critical overdensity
       *  \f$\Delta_{crit}\f$; the function
       *  cbl::Cosmology::Delta_vir can be used to convert
       *  \f$\Delta_{crit}\f$ into \f$\Delta\f$
       *
       *  @warning the mass function by Manera et al. (2010) has been
       *  tested only for z=0 and z=0.5; the mass function by Despali
       *  et al. (2016) is currently implemented only for virial
       *  masses and at \f$z<1.25\f$
       */
      double error_beta_measured (const double Volume, const double density, const double Mass_min, const double Mass_max, const double redshift, const std::string model_bias, const std::string model_MF, const std::string method_SS, const std::string output_root="test", const double Delta=200., const double kk=-1., const std::string interpType="Linear", const int norm=-1, const double k_min=0., const double k_max=100., const double prec=1.e-2, const std::string input_file=par::defaultString, const bool is_parameter_file=true); 

      /**
       *  @brief the normalised quadrupole Q
       *
       *  @param Mass_min minimum halo mass
       * 
       *  @param Mass_max maximum halo mass
       *
       *  @param redshift the redshift
       *
       *  @param model_bias author(s) who proposed the bias; valid
       *  authors are: ST99 (Sheth & Tormen 1999), SMT01 (Sheth, Mo &
       *  Tormen 2001), SMT01_WL04 (Sheth, Mo & Tormen 2001 with the
       *  correction of Warren 2004), Tinker (Tinker et al. 2010)
       *
       *  @param model_MF author(s) who proposed the mass function;
       *  valid authors are: PS (Press&Schechter), ST (Sheth&Tormen),
       *  Jenkins (Jenkins et al. 2001), Warren (Warren et al. 2006),
       *  Reed, (Reed et al. 2007), Pan (Pan 2007), ShenH (halo MF by
       *  Shen et al. 2006), ShenF (filaments MF by Shen et al. 2006),
       *  ShenS (sheets MF by Shen et al. 2006), Tinker (Tinker et
       *  al. 2008), Crocce (Crocce et al. 2010), Angulo_FOF (FOF MF by
       *  Angulo et al. 2012), Angulo_Sub (SUBFIND MF by Angulo et
       *  al. 2012), Watson_FOF(FOF MF by Watson et al. 2012),
       *  Watson_SOH (MF for Spherical Overdensity Haloes by Watson 
       *  et al. 2012), Manera (Manera et al. 2010), Bhattacharya 
       *  (Bhattacharya et al. 2011), Courtin (Courtin et al. 2010),
       *  Peacock (by Peacock at al. 2007)
       *
       *  @param method_SS method used to compute the power spectrum
       *  and &sigma;(mass); valid method_SS are: CAMB
       *  [http://camb.info/], classgal_v1 [http://class-code.net/],
       *  EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       *
       *  @param output_root output_root of the parameter file used to
       *  compute the power spectrum and &sigma;(mass); it can be any
       *  name
       *
       *  @param Delta \f$\Delta\f$: the overdensity, defined as the
       *  mean interior density relative to the background
       *  
       *  @param kk wave vector module
       *
       *  @param interpType method to interpolate the power spectrum
       *
       *  @param norm 0 \f$\rightarrow\f$ don't normalise the power
       *  spectrum; 1 \f$\rightarrow\f$ normalise the power spectrum;
       *  -1 \f$\rightarrow\f$ normalise only if sigma8 is set
       *
       *  @param k_min minimum wave vector module up to which the
       *  power spectrum is computed
       *
       *  @param k_max maximum wave vector module up to which the
       *  power spectrum is computed
       *
       *  @param prec accuracy of the integration 
       *
       *  @param input_file either the parameter file or the power
       *  spectrum file; if a parameter file is provided,
       *  i.e. input_file!=NULL and is_parameter_file=true, it will be
       *  used to compute the power spectrum; if a power spectrum file
       *  is provided, i.e. input_file!=NULL and
       *  is_parameter_file=false, then the provided power spectrum
       *  will be used directly; in both cases &sigma;<SUP>2</SUP>(M)
       *  is computed by integrating the computed/provided power
       *  spectrum ignoring the cosmological parameters of the object
       *
       *  @param is_parameter_file true \f$\rightarrow\f$ the input_file
       *  is a parameter file, used to compute the power spectrum with
       *  the method specified by method_Pk; false \f$\rightarrow\f$
       *  the input_file is a file containing the power spectrum
       *
       *  @return Q: the normalised quadrupole
       *
       *  @warning the input parameter \f$\Delta\f$ is the background
       *  overdensity, not the critical overdensity
       *  \f$\Delta_{crit}\f$; the function
       *  cbl::Cosmology::Delta_vir can be used to convert
       *  \f$\Delta_{crit}\f$ into \f$\Delta\f$
       *
       *  @warning the mass function by Manera et al. (2010) has been
       *  tested only for z=0 and z=0.5; the mass function by Despali
       *  et al. (2016) is currently implemented only for virial
       *  masses and at \f$z<1.25\f$
       */
      double quadrupole (const double Mass_min, const double Mass_max, const double redshift, const std::string model_bias, const std::string model_MF, const std::string method_SS, const std::string output_root="test", const double Delta=200., const double kk=-1., const std::string interpType="Linear", const int norm=-1, const double k_min=0., const double k_max=100., const double prec=1.e-2, const std::string input_file=par::defaultString, const bool is_parameter_file=true); 

      /**
       *  @brief the normalised quadrupole Q
       *
       *  @param MM vector of halo masses
       *
       *  @param MF vector of mass function values, d&Phi;/dM=dn(M)/dM
       *
       *  @param redshift the redshift
       *
       *  @param model_bias author(s) who proposed the bias; valid
       *  authors are: ST99 (Sheth & Tormen 1999), SMT01 (Sheth, Mo &
       *  Tormen 2001), SMT01_WL04 (Sheth, Mo & Tormen 2001 with the
       *  correction of Warren 2004), Tinker (Tinker et al. 2010)
       *
       *  @param method_SS method used to compute the power spectrum
       *  and &sigma;(mass); valid method_SS are: CAMB
       *  [http://camb.info/], classgal_v1 [http://class-code.net/],
       *  EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       *
       *  @param output_root output_root of the parameter file used to
       *  compute the power spectrum and &sigma;(mass); it can be any
       *  name
       *
       *  @param Delta \f$\Delta\f$: the overdensity, defined as the
       *  mean interior density relative to the background
       *  
       *  @param kk wave vector module
       *
       *  @param interpType method to interpolate the power spectrum
       *
       *  @param norm 0 \f$\rightarrow\f$ don't normalise the power
       *  spectrum; 1 \f$\rightarrow\f$ normalise the power spectrum;
       *  -1 \f$\rightarrow\f$ normalise only if sigma8 is set
       *
       *  @param k_min minimum wave vector module up to which the
       *  power spectrum is computed
       *
       *  @param k_max maximum wave vector module up to which the
       *  power spectrum is computed
       *
       *  @param prec accuracy of the integration 
       *
       *  @param input_file either the parameter file or the power
       *  spectrum file; if a parameter file is provided,
       *  i.e. input_file!=NULL and is_parameter_file=true, it will be
       *  used to compute the power spectrum; if a power spectrum file
       *  is provided, i.e. input_file!=NULL and
       *  is_parameter_file=false, then the provided power spectrum
       *  will be used directly; in both cases &sigma;<SUP>2</SUP>(M)
       *  is computed by integrating the computed/provided power
       *  spectrum ignoring the cosmological parameters of the object
       *
       *  @param is_parameter_file true \f$\rightarrow\f$ the input_file
       *  is a parameter file, used to compute the power spectrum with
       *  the method specified by method_Pk; false \f$\rightarrow\f$
       *  the input_file is a file containing the power spectrum
       *
       *  @return Q: the normalised quadrupole
       *
       *  @warning the input parameter \f$\Delta\f$ is the background
       *  overdensity, not the critical overdensity
       *  \f$\Delta_{crit}\f$; the function
       *  cbl::Cosmology::Delta_vir can be used to convert
       *  \f$\Delta_{crit}\f$ into \f$\Delta\f$
       */
      double quadrupole (const std::vector<double> MM, const std::vector<double> MF, const double redshift, const std::string model_bias, const std::string method_SS, const std::string output_root="test", const double Delta=200., const double kk=-1., const std::string interpType="Linear", const int norm=-1, const double k_min=0., const double k_max=100., const double prec=1.e-2, const std::string input_file=par::defaultString, const bool is_parameter_file=true);

      /**
       *  @brief the mean square bulk flow
       *
       *  @param rr comoving radius 
       *
       *  @param k_int_min minimum wave vector module up to which the
       *  integral is computed
       *
       *  @param method_Pk method used to compute the power spectrum
       *  and &sigma;(mass); valid choices for method_Pk are: CAMB
       *  [http://camb.info/], classgal_v1 [http://class-code.net/],
       *  EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       *
       *  @param redshift the redshift
       *
       *  @param output_root output_root of the parameter file used to
       *  compute the power spectrum and &sigma;(mass); it can be any
       *  name
       *
       *  @param k_min minimum wave vector module up to which the
       *  power spectrum is computed
       *
       *  @param k_max maximum wave vector module up to which the
       *  power spectrum is computed
       *
       *  @param prec accuracy of the integration
       *
       *  @param file_par name of the parameter file; if a
       *  parameter file is provided (i.e. file_par!=NULL), it will be
       *  used, ignoring the cosmological parameters of the object
       *
       *  @return the mean square bulk flow
       */
      double square_bulk_flow (const double rr, const double k_int_min, const std::string method_Pk, const double redshift, const std::string output_root="test", const double k_min=0., const double k_max=100., const double prec=1.e-2, const std::string file_par=par::defaultString);

      /**
       *  @brief the mean square bulk flow
       *
       *  @param rr comoving radius 
       *
       *  @param k_int_min minimum wave vector module up to which the
       *  integral is computed
       *
       *  @param lgkk vector of log(k)
       *
       *  @param lgPk vector of log(P(k))
       *
       *  @param redshift the redshift
       *
       *  @return the mean square bulk flow
       */
      double square_bulk_flow_Table (const double rr, const double k_int_min, const std::vector<double> lgkk, const std::vector<double> lgPk, const double redshift) const; 

      /**
       *  @brief the mean square velocity dispersion
       *
       *  @param rr comoving radius 
       *
       *  @param k_int_min minimum wave vector module up to which the
       *  integral is computed
       *
       *  @param method_Pk method used to compute the power spectrum
       *  and &sigma;(mass); valid choices for method_Pk are: CAMB
       *  [http://camb.info/], classgal_v1 [http://class-code.net/],
       *  EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       *
       *  @param redshift the redshift
       *
       *  @param output_root output_root of the parameter file used to
       *  compute the power spectrum and &sigma;(mass); it can be any
       *  name
       *
       *  @param k_min minimum wave vector module up to which the
       *  power spectrum is computed
       *
       *  @param k_max maximum wave vector module up to which the
       *  power spectrum is computed
       *
       *  @param prec accuracy of the integration
       *
       *  @param file_par name of the parameter file; if a
       *  parameter file is provided (i.e. file_par!=NULL), it will be
       *  used, ignoring the cosmological parameters of the object
       *
       *  @return the mean square velocity dispersion
       */
      double square_velocity_dispersion (const double rr, const double k_int_min, const std::string method_Pk, const double redshift, const std::string output_root="test", const double k_min=0., const double k_max=100., const double prec=1.e-2, const std::string file_par=par::defaultString);
    
      /**
       *  @brief the Cosmic Mach Number
       *
       *  @param rr comoving radius 
       *
       *  @param k_int_min minimum wave vector module up to which the
       *  integral is computed
       *
       *  @param method_Pk method used to compute the power spectrum
       *  and &sigma;(mass); valid choices for method_Pk are: CAMB
       *  [http://camb.info/], classgal_v1 [http://class-code.net/],
       *  EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       *
       *  @param redshift the redshift
       *
       *  @param output_root output_root of the parameter file used to
       *  compute the power spectrum and &sigma;(mass); it can be any
       *  name
       *
       *  @param k_max maximum wave vector module up to which the
       *  power spectrum is computed
       *
       *  @param file_par name of the parameter file; if a
       *  parameter file is provided (i.e. file_par!=NULL), it will be
       *  used, ignoring the cosmological parameters of the object
       *
       *  @return the Cosmic Mach Number
       */
      double CMN (const double rr, const double k_int_min, const std::string method_Pk, const double redshift, const std::string output_root="test", const double k_max=100., const std::string file_par=par::defaultString) const;

      /**
       *  @brief the hierarchical moments S<SUB>n</SUB>
       *
       *  this function provides the hierarchical moments S<SUB>n</SUB>
       *  given by the perturbation theory (see e.g. Juszkiewicz et
       *  al. 1993, Bernardeau 1994, Wolk 2013)
       *
       *  @param nn order of the moment
       *
       *  @param RR comoving separation
       *
       *  @param method_SS method used to compute the power spectrum
       *  and &sigma;(mass); valid method_SS are: CAMB
       *  [http://camb.info/], classgal_v1 [http://class-code.net/],
       *  EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       *
       *  @param output_root output_root of the parameter file used to
       *  compute the power spectrum and &sigma;(mass); it can be any
       *  name
       *
       *  @param interpType method to interpolate the power spectrum
       *
       *  @param k_max maximum wave vector module up to which the power
       *  spectrum is computed
       *
       *  @param input_file either the parameter file or the power
       *  spectrum file; if a parameter file is provided,
       *  i.e. input_file!=NULL and is_parameter_file=true, it will be
       *  used to compute the power spectrum; if a power spectrum file
       *  is provided, i.e. input_file!=NULL and
       *  is_parameter_file=false, then the provided power spectrum
       *  will be used directly; in both cases &sigma;<SUP>2</SUP>(M)
       *  is computed by integrating the computed/provided power
       *  spectrum ignoring the cosmological parameters of the object
       *
       *  @param is_parameter_file true \f$\rightarrow\f$ the input_file
       *  is a parameter file, used to compute the power spectrum with
       *  the method specified by method_Pk; false \f$\rightarrow\f$
       *  the input_file is a file containing the power spectrum
       *
       *  @return the hierarchical moments, S<SUB>n</SUB>, given by the
       *  perturbation theory
       */
      double Sn_PT (const int nn, const double RR, const std::string method_SS, const std::string output_root="test", const std::string interpType="Linear", const double k_max=100., const std::string input_file=par::defaultString, const bool is_parameter_file=true) const;
  
      /**
       *  @brief the deprojected hierarchical moments
       *  &Sigma;<SUB>n</SUB>
       *
       *  this function provides the deprojected hierarchical moments
       *  &Sigma;<SUB>n</SUB> given by the perturbation theory (see
       *  e.g. Juszkiewicz et al. 1993, Bernardeau 1994, Wolk 2013)
       *
       *  @param nn order of the moment
       *
       *  @param RR comoving separation
       *
       *  @param method_SS method used to compute the power spectrum
       *  and &sigma;(mass); valid method_SS are: CAMB
       *  [http://camb.info/], classgal_v1 [http://class-code.net/],
       *  EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       *
       *  @param output_root output_root of the parameter file used to
       *  compute the power spectrum and &sigma;(mass); it can be any
       *  name
       *
       *  @param interpType method to interpolate the power spectrum
       *
       *  @param k_max maximum wave vector module up to which the
       *  power spectrum is computed
       *
       *  @param input_file either the parameter file or the power
       *  spectrum file; if a parameter file is provided,
       *  i.e. input_file!=NULL and is_parameter_file=true, it will be
       *  used to compute the power spectrum; if a power spectrum file
       *  is provided, i.e. input_file!=NULL and
       *  is_parameter_file=false, then the provided power spectrum
       *  will be used directly; in both cases &sigma;<SUP>2</SUP>(M)
       *  is computed by integrating the computed/provided power
       *  spectrum ignoring the cosmological parameters of the object
       *
       *  @param is_parameter_file true \f$\rightarrow\f$ the input_file
       *  is a parameter file, used to compute the power spectrum with
       *  the method specified by method_Pk; false \f$\rightarrow\f$
       *  the input_file is a file containing the power spectrum
       *
       *  @return the deprojected hierarchical moments,
       *  &Sigma;<SUB>n</SUB>, given by the perturbation theory
       */
      double Sigman_PT (const int nn, const double RR, const std::string method_SS, const std::string output_root="test", const std::string interpType="Linear", const double k_max=100., const std::string input_file=par::defaultString, const bool is_parameter_file=true) const;
    
      /**
       *  @brief 1D monopole in the Kaiser limit
       *
       *  this function provides the monopole &xi;<SUB>0</SUB>(r)
       *  predicted at large scales, in the Kaiser limit
       *
       *  @param rad comoving separation
       *
       *  @param f_sigma8 f*&sigma;<SUB>8</SUB>
       *
       *  @param bias_sigma8 b*&sigma;<SUB>8</SUB>
       *
       *  @param method_Pk method used to compute the power spectrum
       *  and &sigma;(mass); valid choices for method_Pk are: CAMB
       *  [http://camb.info/], classgal_v1 [http://class-code.net/],
       *  MPTbreeze-v1 [http://arxiv.org/abs/1207.1465], EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       *
       *  @param redshift the redshift
       *
       *  @param output_root output_root of the parameter file used to
       *  compute the power spectrum and &sigma;(mass); it can be any
       *  name
       *
       *  @param xiType 0 \f$\rightarrow\f$ standard; 1
       *  \f$\rightarrow\f$ Chuang & Wang model
       *
       *  @param k_star k<SUB>*</SUB> of the Chuang & Wang model
       *
       *  @param xiNL 0 \f$\rightarrow\f$ linear power spectrum; 1
       *  \f$\rightarrow\f$ non-linear power spectrum
       *
       *  @param norm 0 \f$\rightarrow\f$ don't normalise the power
       *  spectrum; 1 \f$\rightarrow\f$ normalise the power spectrum;
       *  -1 \f$\rightarrow\f$ normalise only if sigma8 is set
       *
       *  @param r_min minimum separation up to which the
       *  correlation function is computed
       *
       *  @param r_max maximum separation up to which the
       *  correlation function is computed
       *
       *  @param k_min minimum wave vector module up to which the
       *  power spectrum is computed
       *
       *  @param k_max maximum wave vector module up to which the
       *  power spectrum is computed
       *
       *  @param aa parameter \e a of Eq. 24 of Anderson et
       *  al. 2012
       *
       *  @param GSL true \f$\rightarrow\f$ the GSL libraries are
       *  used
       *
       *  @param prec accuracy of the integration
       *
       *  @param file_par name of the parameter file; if a
       *  parameter file is provided (i.e. file_par!=NULL), it will be
       *  used, ignoring the cosmological parameters of the object
       *
       *  @return &xi;<SUB>0</SUB>
       *
       */
      double xi0_Kaiser (const double rad, const double f_sigma8, const double bias_sigma8, const std::string method_Pk, const double redshift, const std::string output_root="test", const bool xiType=0, const double k_star=-1., const bool xiNL=0, const int norm=-1, const double r_min=0.1, const double r_max=150., const double k_min=0., const double k_max=100., const double aa=0., const bool GSL=false, const double prec=1.e-2, const std::string file_par=par::defaultString);
    
      /**
       *  @brief 1D monopole in the Kaiser limit
       *
       *  this function provides the monopole &xi;<SUB>0</SUB>(r)
       *  predicted at large scales, in the Kaiser limit
       *
       *  @param rad comoving separations
       *
       *  @param bias b
       *
       *  @param method_Pk method used to compute the power spectrum
       *  and &sigma;(mass); valid choices for method_Pk are: CAMB
       *  [http://camb.info/], classgal_v1 [http://class-code.net/],
       *  MPTbreeze-v1 [http://arxiv.org/abs/1207.1465], EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       *
       *  @param NL 0 \f$\rightarrow\f$ linear power spectrum; 1 \f$\rightarrow\f$ non-linear power
       *  spectrum
       *
       *  @param redshift the redshift
       *
       *  @param output_dir the output_dir directory
       *  where the output of external codes are written
       *
       *  @param output_root output_root of the parameter file used to
       *  compute the power spectrum and &sigma;(mass); it can be any
       *  name
       *
       *  @param norm 0 \f$\rightarrow\f$ don't normalise the power
       *  spectrum; 1 \f$\rightarrow\f$ normalise the power spectrum;
       *  -1 \f$\rightarrow\f$ normalise only if sigma8 is set
       *
       *  @param k_min minimum wave vector module up to which the
       *  power spectrum is computed
       *
       *  @param k_max maximum wave vector module up to which the
       *  power spectrum is computed
       *
       *  @param step number of steps
       *
       *  @param prec accuracy of the integration
       *
       *  @param file_par name of the parameter file; if a
       *  parameter file is provided (i.e. file_par!=NULL), it will be
       *  used, ignoring the cosmological parameters of the object
       *
       *  @return &xi;<SUB>0</SUB>
       *
       */
      std::vector<double> xi0_Kaiser (const std::vector<double> rad, const double bias, const std::string method_Pk, const bool NL, const double redshift, const std::string output_dir, const std::string output_root, const int norm, const double k_min, const double k_max, const int step, const double prec, const std::string file_par);

      /**
       *  @brief 2D correlation function, &xi;(r<SUB>p</SUB>,&pi;),
       *  predicted by the dispersion model
       *
       *  @param rp r<SUB>p</SUB>: the comoving separation perpendicular
       *  to the line-of-sight
       *
       *  @param pi &pi;: the comoving separation parallel to the
       *  line-of-sight
       *
       *  @param f_sigma8 f*&sigma;<SUB>8</SUB>
       *
       *  @param bias_sigma8 b*&sigma;<SUB>8</SUB>
       *
       *  @param sigma12 &sigma;<SUB>12</SUB>: pairwise peculiar
       *  velocity dispersion
       *
       *  @param method_Pk method used to compute the power spectrum
       *  and &sigma;(mass); valid choices for method_Pk are: CAMB
       *  [http://camb.info/], classgal_v1 [http://class-code.net/],
       *  MPTbreeze-v1 [http://arxiv.org/abs/1207.1465], EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       *
       *  @param redshift the redshift
       *
       *  @param FV 0 \f$\rightarrow\f$ exponential form for f(v); 1
       *  \f$\rightarrow\f$ Gaussian form for f(v); where f(v) is the
       *  velocity distribution function
       *
       *  @param NL 0 \f$\rightarrow\f$ linear power spectrum; 1
       *  \f$\rightarrow\f$ non-linear power spectrum
       *
       *  @param rr vector of r, the module of the comoving
       *  separation
       *
       *  @param Xi vector of &xi;(r), the two-point correlation
       *  function of dark matter
       *
       *  @param Xi_ vector of barred &xi;(r),
       *
       *  @param Xi__ vector of double-barred &xi;(r)
       *
       *  @param output_root output_root of the parameter file used to
       *  compute the power spectrum and &sigma;(mass); it can be any
       *  name
       *
       *  @param index internal parameter used when minimizing the
       *  &chi;<SUB>2</SUB>
       *
       *  @param bias_nl 0 \f$\rightarrow\f$ linear bias; 1
       *  \f$\rightarrow\f$ non-linear bias
       *
       *  @param bA b<SUB>a</SUB> non-linear bias parameter
       *
       *  @param xiType 0 \f$\rightarrow\f$ standard; 1
       *  \f$\rightarrow\f$ Chuang & Wang model
       *
       *  @param k_star k<SUB>*</SUB> of the Chuang & Wang model
       *
       *  @param xiNL 0 \f$\rightarrow\f$ linear power spectrum; 1
       *  \f$\rightarrow\f$ non-linear power spectrum
       *
       *  @param v_min minimum velocity used in the convolution of the
       *  correlation function
       *
       *  @param v_max maximum velocity used in the convolution of the
       *  correlation function
       *
       *  @param step_v number of steps used in the convolution of the
       *  correlation function
       *
       *  @param norm 0 \f$\rightarrow\f$ don't normalise the power
       *  spectrum; 1 \f$\rightarrow\f$ normalise the power spectrum;
       *  -1 \f$\rightarrow\f$ normalise only if sigma8 is set
       *
       *  @param r_min minimum separation up to which the
       *  correlation function is computed
       *
       *  @param r_max maximum separation up to which the
       *  correlation function is computed
       *
       *  @param k_min minimum wave vector module up to which the
       *  power spectrum is computed
       *
       *  @param k_max maximum wave vector module up to which the
       *  power spectrum is computed
       *
       *  @param aa parameter \e a of Eq. 24 of Anderson et
       *  al. 2012
       *
       *  @param GSL true \f$\rightarrow\f$ the GSL libraries are
       *  used
       *
       *  @param prec accuracy of the integration
       *
       *  @param file_par name of the parameter file; if a
       *  parameter file is provided (i.e. file_par!=NULL), it will be
       *  used, ignoring the cosmological parameters of the object
       *
       *  @return &xi;(r<SUB>p</SUB>,&pi;)
       */
      double xi2D_DispersionModel (const double rp, const double pi, const double f_sigma8, const double bias_sigma8, const double sigma12, const std::string method_Pk, const double redshift, const int FV, const bool NL, std::vector<double> rr, std::vector<double> &Xi, std::vector<double> &Xi_, std::vector<double> &Xi__, const std::string output_root="test", const int index=-1, const bool bias_nl=0, const double bA=-1., const bool xiType=0, const double k_star=-1., const bool xiNL=0, const double v_min=-3000., const double v_max=3000., const int step_v=500, const int norm=-1, const double r_min=0.1, const double r_max=150., const double k_min=0., const double k_max=100., const double aa=0., const bool GSL=false, const double prec=1.e-2, const std::string file_par=par::defaultString);

      /**
       *  @brief the function &xi;<SUB>*</SUB> of the Chuang & Wang 2012
       *  model
       *
       *  see Chuang & Wang 2012, 1209.0210
       *
       *  @param rr comoving separation
       *
       *  @param redshift the redshift
       *
       *  @param output_root output_root of the parameter file used to
       *  compute the power spectrum and &sigma;(mass); it can be any
       *  name
       *
       *  @param k_star k<SUB>*</SUB> of the Chuang & Wang model
       *
       *  @param k_min minimum wave vector module up to which the
       *  power spectrum is computed
       *
       *  @param k_max maximum wave vector module up to which the
       *  power spectrum is computed
       *
       *  @param prec accuracy of the integration
       *
       *  @param file_par name of the parameter file; if a
       *  parameter file is provided (i.e. file_par!=NULL), it will be
       *  used, ignoring the cosmological parameters of the object
       *
       *  @return &xi;<SUB>*</SUB>
       */
      double xi_star (const double rr, const double redshift, const std::string output_root="test", const double k_star=-1., const double k_min=0., const double k_max=100., const double prec=1.e-2, const std::string file_par=par::defaultString);
  
      /**
       *  @brief the function &xi;<SUB>g,nw</SUB>(s) of the Chuang &
       *  Wang 2012 model
       *
       *  see Chuang & Wang 2012, 1209.0210
       *
       *  @param rp r<SUB>p</SUB>: the comoving separation perpendicular
       *  to the line-of-sight
       *
       *  @param pi &pi;: the comoving separation parallel to the
       *  line-of-sight
       *      
       *  @param f_sigma8 f*&sigma;<SUB>8</SUB>
       *
       *  @param bias_sigma8 b*&sigma;<SUB>8</SUB>
       *
       *  @param bA b<SUB>a</SUB> non-linear bias parameter
       *
       *  @param redshift the redshift
       *    
       *  @param rr vector of r, the module of the comoving
       *  separation
       *
       *  @param Xi vector of &xi;(r), the two-point correlation
       *  function of dark matter
       *
       *  @param Xi_ vector of barred &xi;(r),
       *
       *  @param Xi__ vector of double-barred &xi;(r)
       *
       *  @param output_root output_root of the parameter file used to
       *  compute the power spectrum and &sigma;(mass); it can be any
       *  name
       *
       *  @return &xi;<SUB>g,nw</SUB>(s)
       */
      double xisnl_gnw (const double rp, const double pi, const double f_sigma8, const double bias_sigma8, const double bA, const double redshift, std::vector<double> rr, std::vector<double> Xi, std::vector<double> &Xi_, std::vector<double> &Xi__, const std::string output_root="test");
 
      /**
       *  @brief the function &xi;<SUB>g,BAO</SUB>(s) of the Chuang &
       *  Wang 2012 model
       *
       *  see Chuang & Wang 2012, 1209.0210
       *
       *  @param rp r<SUB>p</SUB>: the comoving separation perpendicular
       *  to the line-of-sight
       *
       *  @param pi &pi;: the comoving separation parallel to the
       *  line-of-sight
       *      
       *  @param f_sigma8 f*&sigma;<SUB>8</SUB>
       *
       *  @param bias_sigma8 b*&sigma;<SUB>8</SUB>
       *
       *  @param redshift the redshift
       *
       *  @param rr vector of r, the module of the comoving
       *  separation
       *
       *  @param Xi vector of &xi;(r), the two-point correlation
       *  function of dark matter
       *
       *  @param Xi_ vector of barred &xi;(r),
       *
       *  @param Xi__ vector of double-barred &xi;(r)
       *
       *  @param output_root output_root of the parameter file used to
       *  compute the power spectrum and &sigma;(mass); it can be any
       *  name
       *
       *  @param k_star k<SUB>*</SUB> of the Chuang & Wang model
       *
       *  @param x_min minimum velocity used in the integral of the
       *  Chuang & Wang model
       *
       *  @param x_max maximum velocity used in the integral of the
       *  Chuang & Wang model
       *
       *  @param step_x number of steps in the integral of the Chuang &
       *  Wang model
       *
       *  @return &xi;<SUB>g,BAO</SUB>(s)
       */
      double xis_gBAO (const double rp, const double pi, const double f_sigma8, const double bias_sigma8, const double redshift, std::vector<double> rr, std::vector<double> Xi, std::vector<double> &Xi_, std::vector<double> &Xi__, const std::string output_root="test", const double k_star=-1., const double x_min=-3000., const double x_max=3000., const int step_x=500);
 
      /**
       *  @brief 2D correlation function, &xi;(r<SUB>p</SUB>,&pi;),
       *  predicted by the Chuang & Wang model
       *
       *  @param rp r<SUB>p</SUB>: the comoving separation perpendicular
       *  to the line-of-sight
       *
       *  @param pi &pi;: the comoving separation parallel to the
       *  line-of-sight
       *
       *  @param beta &beta;=f/b, where f is the linear growth rate
       *  and b is the bias
       *
       *  @param bias_lin linear bias
       *
       *  @param bA b<SUB>a</SUB> non-linear bias parameter
       *
       *  @param sigmav0 &sigma;<SUB>0</SUB>(v): parameter of the
       *  velocity distribution function, f(v)
       *
       *  @param cmu parameter of the velocity distribution function,
       *  f(v)
       *
       *  @param cs1 parameter of the velocity distribution function,
       *  f(v)
       *
       *  @param cs2 parameter of the velocity distribution function,
       *  f(v)
       *
       *  @param redshift the redshift
       *
       *  @param rr1 vector of r, the module of the comoving separation
       *
       *  @param Xi1 vector of &xi;(r), the two-point correlation
       *  function of dark matter
       *
       *  @param rr2 vector of r, the module of the comoving separation
       *
       *  @param Xi2 vector of &xi;(r), the two-point correlation
       *  function of dark matter
       *
       *  @param Xi1_ vector of barred &xi;(r),
       *
       *  @param Xi1__ vector of double-barred &xi;(r)
       *   
       *  @param Xi2_ vector of barred &xi;(r),
       *
       *  @param Xi2__ vector of double-barred &xi;(r)
       *
       *  @param output_root output_root of the parameter file used to
       *  compute the power spectrum and &sigma;(mass); it can be any
       *  name
       *
       *  @param BAO 0 \f$\rightarrow\f$ no BAO convolution; 1 \f$\rightarrow\f$ BAO
       *  convolution
       *
       *  @param xiType 0 \f$\rightarrow\f$ standard; 1 \f$\rightarrow\f$ Chuang & Wang model
       *
       *  @param k_star k<SUB>*</SUB> of the Chuang & Wang model
       *
       *  @param xiNL 0 \f$\rightarrow\f$ linear power spectrum; 1 \f$\rightarrow\f$
       *  non-linear power spectrum
       *
       *  @param r_min minimum separation up to which the
       *  correlation function is computed
       *
       *  @param r_max maximum separation up to which the
       *  correlation function is computed
       *
       *  @param v_min minimum velocity used in the convolution of the
       *  correlation function
       *
       *  @param v_max maximum velocity used in the convolution of the
       *  correlation function
       *
       *  @param step_v number of steps used in the convolution of the
       *  correlation function
       *
       *  @param k_min minimum wave vector module up to which the
       *  power spectrum is computed
       *
       *  @param k_max maximum wave vector module up to which the
       *  power spectrum is computed
       *
       *  @param x_min minimum velocity used in the integral of the
       *  Chuang & Wang model
       *
       *  @param x_max maximum velocity used in the integral of the
       *  Chuang & Wang model
       *
       *  @param step_x number of steps in the integral of the Chuang &
       *  Wang model
       *
       *  @param aa parameter \e a of Eq. 24 of Anderson et
       *  al. 2012
       *
       *  @param GSL true \f$\rightarrow\f$ the GSL libraries are
       *  used
       *
       *  @param prec accuracy of the integration
       *
       *  @param file_par name of the parameter file; if a
       *  parameter file is provided (i.e. file_par!=NULL), it will be
       *  used, ignoring the cosmological parameters of the object
       *
       *  @return &xi;(r<SUB>p</SUB>,&pi;)
       */
      double xi2D_CW (const double rp, const double pi, const double beta, const double bias_lin, const double bA, const double sigmav0, const double cmu, const double cs1, const double cs2, const double redshift, std::vector<double> rr1, std::vector<double> Xi1, std::vector<double> rr2, std::vector<double> Xi2, std::vector<double> &Xi1_, std::vector<double> &Xi1__, std::vector<double> &Xi2_, std::vector<double> &Xi2__, const std::string output_root="test", const bool BAO=1, const bool xiType=0, const double k_star=-1, const bool xiNL=0, const double r_min=0.1, const double r_max=150., const double v_min=-3000., const double v_max=3000., const int step_v=500, const double k_min=0., const double k_max=100., const double x_min=-3000., const double x_max=3000., const int step_x=500, const double aa=0., const bool GSL=false, const double prec=1.e-2, const std::string file_par=par::defaultString);

      ///@}


      /**
       *  @name Functions to model baryon acoustic oscillations
       */
      ///@{
    
      /**
       *  @brief the sound horizon at the drag epoch
       *  r<SUB>s</SUB>(z<SUB>d</SUB>), valid choices for method_Pk
       *  are: EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html],
       *  CAMB [http://camb.info/]
       *
       *  @author Alfonso Veropalumbo
       *  @author alfonso.veropalumbo@unibo.it
       *
       *  @param method_Pk the method to compute the sound horizon
       *  
       *  @param T_CMB T<SUB>CMB</SUB>: the present day CMB temperature [K]
       *     
       *  @return r<SUB>s</SUB>
       */
      double rs (const std::string method_Pk, const double T_CMB=par::TCMB) const;

      /**
       *  @brief the sound horizon at the drag epoch predicted by
       *  Eisentein & Hu 1998 
       *
       *  see Eisentein & Hu 1998, Section 2.1
       *
       *  @author Alfonso Veropalumbo 
       *  @author alfonso.veropalumbo@unibo.it 
       *  @param T_CMB CMB temperature
       *  @return r<SUB>s</SUB>
       */
      double rs_EH (const double T_CMB=par::TCMB) const;

      /**
       *  @brief the sound horizon at the drag epoch estimated with CAMB [http://camb.info/],
       *  analytical formula by Aubourg et al. 2014
       *
       *  see Anderson et al 2014, Eq. 16
       *  
       *  @author Alfonso Veropalumbo
       *  @author alfonso.veropalumbo@unibo.it
       *  @return r<SUB>s</SUB>
       */
      double rs_CAMB () const;
  
      /**
       *  @brief the fiducial cosmology independent ratio
       *  r<SUB>s</SUB>/D<SUB>V</SUB>,  valid choices for method_Pk are: EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html],
       *  CAMB [http://camb.info/]
       *
       *  both r<SUB>s</SUB> and D<SUB>V</SUB> are in Mpc
       *
       *  @author Alfonso Veropalumbo
       *  @author alfonso.veropalumbo@unibo.it
       *
       *  @param redshift the redshift
       *  @param method_Pk method used to compute the sound horizon;
       *  @param T_CMB CMB temperature
       *
       *  @return y<SUB>s</SUB>
       */
      double ys (const double redshift, const std::string method_Pk, const double T_CMB=par::TCMB) const;
  
      /**
       *  @brief the acoustic parameter 
       *
       *  see Eisenstein 2005 
       *
       *  @author Alfonso Veropalumbo
       *  @author alfonso.veropalumbo@unibo.it
       *  @param redshift the redshift
       *  @return the acoustic parameter 
       */
      double Az (const double redshift) const;
  
      /**
       *  @brief the linear point  
       *
       *  see Anselmi et al. 2016 
       *
       *  @author Alfonso Veropalumbo
       *  @author alfonso.veropalumbo@unibo.it
       *
       *  @param redshift the redshift
       *
       *  @param rmin the minimum scale
       *
       *  @param rmax the maximum scale
       *
       *  @param nbinr the number of scale bins
       *
       *  @param interpType the interpolation type
       *
       *  @return vector containing the linear point, the dip 
       *  and the BAO peak for the correlation function at the redshift
       *  provided
       *
       */
      std::vector<double> linear_point (const double redshift, const double rmin=60., const double rmax=150., const int nbinr=100, const std::string interpType="Spline");

      ///@}


      /**
       *  @name Functions to model cosmological quantities in non-Gaussian cosmologies
       */
      ///@{

      /**
       *  @brief the amplitude of the matter power spectrum
       *
       *  @author Cosimo Fedeli
       *  @author cosimo.fedeli@oabo.inaf.it
       *
       *  @param method_Pk method used to compute the power spectrum;
       *  valid choices for method_Pk are: CAMB [http://camb.info/],
       *  classgal_v1 [http://class-code.net/], MPTbreeze-v1
       *  [http://arxiv.org/abs/1207.1465], EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       * 
       *  @param output_root the output_root parameter of the
       *  parameter file used to compute the power spectrum; it can be
       *  any name
       *
       *  @param norm 0 \f$\rightarrow\f$ don't normalise the power
       *  spectrum; 1 \f$\rightarrow\f$ normalise the power spectrum;
       *  -1 \f$\rightarrow\f$ normalise only if sigma8 is set
       *
       *  @param k_min minimum wave vector module up to which the power
       *  spectrum is computed
       *
       *  @param k_max maximum wave vector module up to which the power
       *  spectrum is computed
       *
       *  @param prec accuracy of the integration 
       *       
       *  @param file_par name of the parameter file; if a parameter
       *  file is provided (i.e. file_par!=NULL), it will use be used,
       *  ignoring the cosmological parameters of the object
       *
       *  @return A<SUB>m</SUB>
       */
      double Am (const std::string method_Pk, const std::string output_root="test", const int norm=-1, const double k_min=0., const double k_max=100., const double prec=1.e-2, const std::string file_par=par::defaultString); 

      /**
       *  @brief the potential spectral amplitude 
       *
       *  @author Cosimo Fedeli
       *  @author cosimo.fedeli@oabo.inaf.it
       *
       *  @param method_Pk method used to compute the power spectrum;
       *  valid choices for method_Pk are: CAMB [http://camb.info/],
       *  classgal_v1 [http://class-code.net/], MPTbreeze-v1
       *  [http://arxiv.org/abs/1207.1465], EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       * 
       *  @param output_root the output_root parameter of the
       *  parameter file used to compute the power spectrum; it can be
       *  any name
       *
       *  @param norm 0 \f$\rightarrow\f$ don't normalise the power
       *  spectrum; 1 \f$\rightarrow\f$ normalise the power spectrum;
       *  -1 \f$\rightarrow\f$ normalise only if sigma8 is set
       *
       *  @param k_min minimum wave vector module up to which the power
       *  spectrum is computed
       *
       *  @param k_max maximum wave vector module up to which the power
       *  spectrum is computed
       *
       *  @param prec accuracy of the integration 
       *
       *  @param file_par name of the parameter file; if a parameter
       *  file is provided (i.e. file_par!=NULL), it will use be used,
       *  ignoring the cosmological parameters of the object
       *
       *  @return the potential spectral amplitude
       */
      double potential_spectral_amplitude (const std::string method_Pk, const std::string output_root="test", const int norm=-1, const double k_min=0., const double k_max=100., const double prec=1.e-2, const std::string file_par=par::defaultString);

      /**
       *  @brief the bispectrum
       *
       *  @author Cosimo Fedeli
       *  @author cosimo.fedeli@oabo.inaf.it
       *
       *  @param kk wave vector module
       *  
       *  @param method_Pk method used to compute the power spectrum;
       *  valid choices for method_Pk are: CAMB [http://camb.info/],
       *  classgal_v1 [http://class-code.net/], MPTbreeze-v1
       *  [http://arxiv.org/abs/1207.1465], EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       * 
       *  @param output_root the output_root parameter of the
       *  parameter file used to compute the power spectrum; it can be
       *  any name
       *
       *  @param norm 0 \f$\rightarrow\f$ don't normalise the power
       *  spectrum; 1 \f$\rightarrow\f$ normalise the power spectrum;
       *  -1 \f$\rightarrow\f$ normalise only if sigma8 is set
       *
       *  @param k_min minimum wave vector module up to which the power
       *  spectrum is computed
       *
       *  @param k_max maximum wave vector module up to which the power
       *  spectrum is computed
       *
       *  @param prec accuracy of the integration 
       *
       *  @param file_par name of the parameter file; if a parameter
       *  file is provided (i.e. file_par!=NULL), it will use be used,
       *  ignoring the cosmological parameters of the object
       *
       *  @return the potential spectral amplitude
       */
      double bispectrum (const std::vector<double> kk, const std::string method_Pk, const std::string output_root="test", const int norm=-1, const double k_min=0., const double k_max=100., const double prec=1.e-2, const std::string file_par=par::defaultString);
    
      /**
       *  @brief auxiliary function to estimate cosmological quantities
       *  in non-Gaussian cosmologies
       *
       *  @author Cosimo Fedeli
       *  @author cosimo.fedeli@oabo.inaf.it
       *
       *  @param kk wave vector module
       *
       *  @param mass halo mass
       *  
       *  @param method_Pk method used to compute the power spectrum; 
       *  valid choices for method_Pk are: CAMB [http://camb.info/], classgal_v1
       *  [http://class-code.net/], MPTbreeze-v1
       *  [http://arxiv.org/abs/1207.1465], EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       * 
       *  @param output_root the output_root parameter of the
       *  parameter file used to compute the power spectrum; it can be
       *  any name
       *
       *  @param norm 0 \f$\rightarrow\f$ don't normalise the power
       *  spectrum; 1 \f$\rightarrow\f$ normalise the power spectrum;
       *  -1 \f$\rightarrow\f$ normalise only if sigma8 is set
       *
       *  @param k_min minimum wave vector module up to which the power
       *  spectrum is computed
       *
       *  @param k_max maximum wave vector module up to which the power
       *  spectrum is computed
       *
       *  @param prec accuracy of the integration 
       *
       *  @param file_par name of the parameter file; if a parameter
       *  file is provided (i.e. file_par!=NULL), it will use be used,
       *  ignoring the cosmological parameters of the object
       *
       *  @return mrk
       */
      double mrk (const double kk, const double mass, const std::string method_Pk, const std::string output_root="test", const int norm=-1, const double k_min=0., const double k_max=100., const double prec=1.e-2, const std::string file_par=par::defaultString);

      /**
       *  @brief auxiliary function to estimate cosmological quantities
       *  in non-Gaussian cosmologies
       *
       *  @author Cosimo Fedeli
       *  @author cosimo.fedeli@oabo.inaf.it
       *
       *  @param kk wave vector module
       *
       *  @param mass halo mass
       *  
       *  @param method_Pk method used to compute the power spectrum;
       *  valid choices for method_Pk are: CAMB [http://camb.info/],
       *  classgal_v1 [http://class-code.net/], MPTbreeze-v1
       *  [http://arxiv.org/abs/1207.1465], EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       * 
       *  @param output_root the output_root parameter of the
       *  parameter file used to compute the power spectrum; it can be
       *  any name
       *
       *  @param interpType method to interpolate the power spectrum
       *
       *  @param norm 0 \f$\rightarrow\f$ don't normalise the power
       *  spectrum; 1 \f$\rightarrow\f$ normalise the power spectrum;
       *  -1 \f$\rightarrow\f$ normalise only if sigma8 is set
       *
       *  @param k_min minimum wave vector module up to which the power
       *  spectrum is computed
       *
       *  @param k_max maximum wave vector module up to which the power
       *  spectrum is computed
       *
       *  @param prec accuracy of the integration 
       *
       *  @param input_file either the parameter file or the power
       *  spectrum file; if a parameter file is provided,
       *  i.e. input_file!=NULL and is_parameter_file=true, it will be
       *  used to compute the power spectrum; if a power spectrum file
       *  is provided, i.e. input_file!=NULL and
       *  is_parameter_file=false, then the provided power spectrum
       *  will be used directly; in both cases &sigma;<SUP>2</SUP>(M)
       *  is computed by integrating the computed/provided power
       *  spectrum ignoring the cosmological parameters of the object
       *
       *  @param is_parameter_file true \f$\rightarrow\f$ the input_file
       *  is a parameter file, used to compute the power spectrum with
       *  the method specified by method_Pk; false \f$\rightarrow\f$
       *  the input_file is a file containing the power spectrum
       *
       *  @return frk
       */
      double frk (const double kk, const double mass, const std::string method_Pk, const std::string output_root="test", const std::string interpType="Linear", const int norm=-1, const double k_min=0., const double k_max=100., const double prec=1.e-2, const std::string input_file=par::defaultString, const bool is_parameter_file=true);

      /// @cond TEST_NG
      double bias_kernel (const double, void *); 

      double frk_test (const double, const double, const std::string, const std::string output_root="test", const std::string interpType="Linear", const int norm=-1, const double k_min=0., const double k_max=100., const double prec=1.e-2, const std::string input_file=par::defaultString, const bool is_parameter_file=true);
      /// @endcond


      /**
       *  @brief correction to the halo bias in non-Gaussian cosmologies
       *
       *  @author Cosimo Fedeli
       *  @author cosimo.fedeli@oabo.inaf.it
       *
       *  @param kk wave vector module
       *
       *  @param mass halo mass
       *  
       *  @param method_Pk method used to compute the power spectrum;
       *  valid choices for method_Pk are: CAMB [http://camb.info/],
       *  classgal_v1 [http://class-code.net/], MPTbreeze-v1
       *  [http://arxiv.org/abs/1207.1465], EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       * 
       *  @param output_root the output_root parameter of the
       *  parameter file used to compute the power spectrum; it can be
       *  any name
       *
       *  @param interpType method to interpolate the power spectrum
       *
       *  @param norm 0 \f$\rightarrow\f$ don't normalise the power
       *  spectrum; 1 \f$\rightarrow\f$ normalise the power spectrum;
       *  -1 \f$\rightarrow\f$ normalise only if sigma8 is set
       *
       *  @param k_min minimum wave vector module up to which the power
       *  spectrum is computed
       *
       *  @param k_max maximum wave vector module up to which the power
       *  spectrum is computed
       *
       *  @param prec accuracy of the integration 
       *
       *  @param input_file either the parameter file or the power
       *  spectrum file; if a parameter file is provided,
       *  i.e. input_file!=NULL and is_parameter_file=true, it will be
       *  used to compute the power spectrum; if a power spectrum file
       *  is provided, i.e. input_file!=NULL and
       *  is_parameter_file=false, then the provided power spectrum
       *  will be used directly; in both cases &sigma;<SUP>2</SUP>(M)
       *  is computed by integrating the computed/provided power
       *  spectrum ignoring the cosmological parameters of the object
       *
       *  @param is_parameter_file true \f$\rightarrow\f$ the input_file
       *  is a parameter file, used to compute the power spectrum with
       *  the method specified by method_Pk; false \f$\rightarrow\f$
       *  the input_file is a file containing the power spectrum
       *
       *  @return bias correction
       */
      double bias_correction (const double kk, const double mass, const std::string method_Pk, const std::string  output_root="test", const std::string interpType="Linear", const int norm=-1, const double k_min=0., const double k_max=100., const double prec=1.e-2, const std::string input_file=par::defaultString, const bool is_parameter_file=true);

      /**
       *  @brief the skewness
       *
       *  @author Cosimo Fedeli
       *  @author cosimo.fedeli@oabo.inaf.it
       *
       *  @param mass halo mass
       *  
       *  @param method_Pk method used to compute the power spectrum;
       *  valid choices for method_Pk are: CAMB [http://camb.info/],
       *  classgal_v1 [http://class-code.net/], MPTbreeze-v1
       *  [http://arxiv.org/abs/1207.1465], EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       * 
       *  @param output_root the output_root parameter of the
       *  parameter file used to compute the power spectrum; it can be
       *  any name
       *
       *  @param interpType method to interpolate the power spectrum
       *
       *  @param norm 0 \f$\rightarrow\f$ don't normalise the power
       *  spectrum; 1 \f$\rightarrow\f$ normalise the power spectrum;
       *  -1 \f$\rightarrow\f$ normalise only if sigma8 is set
       *
       *  @param k_min minimum wave vector module up to which the power
       *  spectrum is computed
       *
       *  @param k_max maximum wave vector module up to which the power
       *  spectrum is computed
       *
       *  @param prec accuracy of the integration 
       *
       *  @param input_file either the parameter file or the power
       *  spectrum file; if a parameter file is provided,
       *  i.e. input_file!=NULL and is_parameter_file=true, it will be
       *  used to compute the power spectrum; if a power spectrum file
       *  is provided, i.e. input_file!=NULL and
       *  is_parameter_file=false, then the provided power spectrum
       *  will be used directly; in both cases &sigma;<SUP>2</SUP>(M)
       *  is computed by integrating the computed/provided power
       *  spectrum ignoring the cosmological parameters of the object
       *
       *  @param is_parameter_file true \f$\rightarrow\f$ the input_file
       *  is a parameter file, used to compute the power spectrum with
       *  the method specified by method_Pk; false \f$\rightarrow\f$
       *  the input_file is a file containing the power spectrum
       *
       *  @return skewness
       */
      double skewness (const double mass, const std::string method_Pk, const std::string output_root="test", const std::string interpType="Linear", const int norm=-1, const double k_min=0., const double k_max=100., const double prec=1.e-2, const std::string input_file=par::defaultString, const bool is_parameter_file=true);

      /**
       *  @brief the derivative of the skewness, ds/dM
       *
       *  @author Cosimo Fedeli
       *  @author cosimo.fedeli@oabo.inaf.it
       *
       *  @param mass halo mass
       *  
       *  @param method_Pk method used to compute the power spectrum; valid
       *  choices for method_Pk are: CAMB [http://camb.info/], classgal_v1
       *  [http://class-code.net/], MPTbreeze-v1
       *  [http://arxiv.org/abs/1207.1465], EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       * 
       *  @param output_root the output_root parameter of the
       *  parameter file used to compute the power spectrum; it can be
       *  any name
       *
       *  @param interpType method to interpolate the power spectrum
       *
       *  @param norm 0 \f$\rightarrow\f$ don't normalise the power
       *  spectrum; 1 \f$\rightarrow\f$ normalise the power spectrum;
       *  -1 \f$\rightarrow\f$ normalise only if sigma8 is set
       *
       *  @param k_min minimum wave vector module up to which the power
       *  spectrum is computed
       *
       *  @param k_max maximum wave vector module up to which the power
       *  spectrum is computed
       *
       *  @param prec accuracy of the integration 
       *
       *  @param input_file either the parameter file or the power
       *  spectrum file; if a parameter file is provided,
       *  i.e. input_file!=NULL and is_parameter_file=true, it will be
       *  used to compute the power spectrum; if a power spectrum file
       *  is provided, i.e. input_file!=NULL and
       *  is_parameter_file=false, then the provided power spectrum
       *  will be used directly; in both cases &sigma;<SUP>2</SUP>(M)
       *  is computed by integrating the computed/provided power
       *  spectrum ignoring the cosmological parameters of the object
       *
       *  @param is_parameter_file true \f$\rightarrow\f$ the input_file
       *  is a parameter file, used to compute the power spectrum with
       *  the method specified by method_Pk; false \f$\rightarrow\f$
       *  the input_file is a file containing the power spectrum
       *
       *  @return derivative of the skewness
       */
      double dskewnessdM (const double mass, const std::string method_Pk, const std::string output_root="test", const std::string interpType="Linear", const int norm=-1, const double k_min=0., const double k_max=100., const double prec=1.e-2, const std::string input_file=par::defaultString, const bool is_parameter_file=true);

      /**
       *  @brief correction to the halo mass in non-Gaussian cosmologies
       *
       *  @author Cosimo Fedeli
       *  @author cosimo.fedeli@oabo.inaf.it
       *
       *  @param mass the halo mass
       *  
       *  @param redshift the redshift
       *
       *  @param method_Pk method used to compute the power spectrum;
       *  valid choices for method_Pk are: CAMB [http://camb.info/],
       *  classgal_v1 [http://class-code.net/], MPTbreeze-v1
       *  [http://arxiv.org/abs/1207.1465], EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       * 
       *  @param output_root the output_root parameter of the
       *  parameter file used to compute the power spectrum; it can be
       *  any name
       *
       *  @param interpType method to interpolate the power spectrum
       *   
       *  @param norm 0 \f$\rightarrow\f$ don't normalise the power
       *  spectrum; 1 \f$\rightarrow\f$ normalise the power spectrum
       *
       *  @param k_min minimum wave vector module up to which the power
       *  spectrum is computed
       *
       *  @param k_max maximum wave vector module up to which the power
       *  spectrum is computed
       *
       *  @param prec accuracy of the integration 
       *
       *  @param input_file either the parameter file or the power
       *  spectrum file; if a parameter file is provided,
       *  i.e. input_file!=NULL and is_parameter_file=true, it will be
       *  used to compute the power spectrum; if a power spectrum file
       *  is provided, i.e. input_file!=NULL and
       *  is_parameter_file=false, then the provided power spectrum
       *  will be used directly; in both cases &sigma;<SUP>2</SUP>(M)
       *  is computed by integrating the computed/provided power
       *  spectrum ignoring the cosmological parameters of the object
       *
       *  @param is_parameter_file true \f$\rightarrow\f$ the input_file
       *  is a parameter file, used to compute the power spectrum with
       *  the method specified by method_Pk; false \f$\rightarrow\f$
       *  the input_file is a file containing the power spectrum
       *
       *  @return bias correction
       */
      double MF_correction (const double mass, const double redshift, const std::string method_Pk, const std::string output_root="test", const std::string interpType="Linear", const int norm=-1, const double k_min=0., const double k_max=100., const double prec=1.e-2, const std::string input_file=par::defaultString, const bool is_parameter_file=true);

      ///@}


      /**
       *  @name Functions to estimate the void size function
       */
      ///@{

      /**
       *  @brief \f$f_{\ln \sigma}(\sigma)\f$ (approximation)
       *
       *  @author Tommaso Ronconi
       *  @author tommaso.ronconi@studio.unibo.it
       *
       *  @param SS variance of the linear density field
       *  (\f$\sigma^2(R)\f$)
       *
       *  @param del_v linear density contrast defining a void
       *
       *  @param del_c critical value of the linear density field
       *  
       *  @return the fraction of trajectories that evolve into voids,
       *  as given in equation (8) of Jennings et al. (2013)
       */
      double f_nu (const double SS, const double del_v, const double del_c) const;

      /**
       *  @brief the void size function
       *
       *  @author Tommaso Ronconi
       *  @author tommaso.ronconi@studio.unibo.it
       *
       *  @param RV radius
       *
       *  @param redshift the redshift
       *
       *  @param del_v linear density contrast defining a void
       *
       *  @param del_c critical value of the linear density field
       *
       *  @param model size function model name; valid choices for
       *  model name are SvdW (Sheth and van de Weygaert, 2004),
       *  linear and Vdn (Jennings et al., 2013)
       *
       *  @param method_Pk method used to compute the power spectrum;
       *  valid choices for method_Pk are: CAMB [http://camb.info/],
       *  classgal_v1 [http://class-code.net/], MPTbreeze-v1
       *  [http://arxiv.org/abs/1207.1465], EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       *
       *  @param output_root output_root of the parameter file used to
       *  compute the power spectrum and &sigma;(mass); it can be any
       *  name
       *
       *  @param interpType method to interpolate the power spectrum
       *
       *  @param k_max maximum wave vector module up to which the power
       *  spectrum is computed
       *           
       *  @param input_file either the parameter file or the power
       *  spectrum file; if a parameter file is provided,
       *  i.e. input_file!=NULL and is_parameter_file=true, it will be
       *  used to compute the power spectrum; if a power spectrum file
       *  is provided, i.e. input_file!=NULL and
       *  is_parameter_file=false, then the provided power spectrum
       *  will be used directly; in both cases &sigma;<SUP>2</SUP>(M)
       *  is computed by integrating the computed/provided power
       *  spectrum ignoring the cosmological parameters of the object
       *
       *  @param is_parameter_file true \f$\rightarrow\f$ the input_file
       *  is a parameter file, used to compute the power spectrum with
       *  the method specified by method_Pk; false \f$\rightarrow\f$
       *  the input_file is a file containing the power spectrum
       *
       *  @return the number density of voids as a function of radius.
       *  Volume Conserving Model, equation (17) from Jennings et
       *  al.(2013)
       */
      double size_function (const double RV, const double redshift, const double del_v, const double del_c, const std::string model, const std::string method_Pk="CAMB", const std::string output_root="test", const std::string interpType="Linear", const double k_max=100., const std::string input_file=par::defaultString, const bool is_parameter_file=true) const;

      /**
       *  @brief the void size function
       *
       *  @author Tommaso Ronconi
       *  @author tommaso.ronconi@studio.unibo.it
       *
       *  @param RV radius
       *
       *  @param redshift the redshift
       *
       *  @param model_mf author(s) who proposed the mass function;
       *  valid authors are: PS (Press & Schechter), ST (Sheth &
       *  Tormen), Jenkins (Jenkins et al. 2001), Warren (Warren et
       *  al. 2006), Reed, (Reed et al. 2007), Pan (Pan 2007), ShenH
       *  (halo MF by Shen et al. 2006), ShenF (filaments MF by Shen
       *  et al. 2006), ShenS (sheets MF by Shen et al. 2006), Tinker
       *  (Tinker et al. 2008), Crocce (Crocce et al. 2010),
       *  Angulo_FOF (FOF MF by Angulo et al. 2012), Angulo_Sub
       *  (SUBFIND MF by Angulo et al. 2012), Watson_FOF(FOF MF by
       *  Watson et al. 2012), Watson_SOH (MF for Spherical Overdensity
       *  Haloes by Watson et al. 2012), Manera (Manera et al. 2010),
       *  Bhattacharya (Bhattacharya et al. 2011), Courtin (Courtin 
       *  et al. 2010), Peacock (by Peacock at al. 2007)
       *
       *  @param del_v linear density contrast defining a void
       *
       *  @param model_sf size function model name; valid choices for
       *  model name are SvdW (Sheth and van de Weygaert, 2004),
       *  linear and Vdn (Jennings et al., 2013)
       *
       *  @param method_Pk method used to compute the power spectrum;
       *  valid choices for method_Pk are: CAMB [http://camb.info/],
       *  classgal_v1 [http://class-code.net/], MPTbreeze-v1
       *  [http://arxiv.org/abs/1207.1465], EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       *
       *  @param output_root output_root of the parameter file used to
       *  compute the power spectrum and &sigma;(mass); it can be any
       *  name
       *
       *  @param Delta \f$\Delta\f$: the overdensity, defined as the
       *  mean interior density relative to the background
       *
       *  @param interpType method to interpolate the power spectrum
       *
       *  @param norm 0 \f$\rightarrow\f$ don't normalise the power
       *  spectrum; 1 \f$\rightarrow\f$ normalise the power spectrum;
       *  -1 \f$\rightarrow\f$ normalise only if sigma8 is set
       *
       *  @param k_min minimum wave vector module up to which the
       *  power spectrum is computed
       *
       *  @param k_max maximum wave vector module up to which the power
       *  spectrum is computed
       *
       *  @param prec accuracy of the integration
       *
       *  @param input_file either the parameter file or the power
       *  spectrum file; if a parameter file is provided,
       *  i.e. input_file!=NULL and is_parameter_file=true, it will be
       *  used to compute the power spectrum; if a power spectrum file
       *  is provided, i.e. input_file!=NULL and
       *  is_parameter_file=false, then the provided power spectrum
       *  will be used directly; in both cases &sigma;<SUP>2</SUP>(M)
       *  is computed by integrating the computed/provided power
       *  spectrum ignoring the cosmological parameters of the object
       *
       *  @param is_parameter_file true \f$\rightarrow\f$ the
       *  input_file is a parameter file, used to compute the power
       *  spectrum with the method specified by method_Pk; false
       *  \f$\rightarrow\f$ the input_file is a file containing the
       *  power spectrum
       *  
       *  @return the number density of voids as a function of radius.
       *  Volume Conserving Model, equation (17) from Jennings et
       *  al.(2013)
       */
      double size_function (const double RV, const double redshift, const std::string model_mf, const double del_v, const std::string model_sf, const std::string method_Pk="CAMB", const std::string output_root="test", const double Delta=200., const std::string interpType="Linear", const int norm=-1, const double k_min=0., const double k_max=100., const double prec=1.e-2, const std::string input_file=par::defaultString, const bool is_parameter_file=true);
      
      ///@}

      
      /**
       *  @name Functions to estimate the multipoles/wedges covariance matrix
       */
      ///@{

      /**
       *  @brief the first three non-null multipoles of the two-point
       *  correlation function
       *
       *  @param nbins the number of bins of the two-point
       *  correlation function multipoles
       *
       *  @param rMin the minimum scale
       *
       *  @param rMax the maximum scale 
       *
       *  @param kk vector containing the wave vector modules
       *
       *  @param Pk0 vector containing the monopole of the power
       *  spectrum
       *
       *  @param Pk2 vector containing the quadrupole of the power
       *  spectrum
       * 
       *  @param Pk4 vector containing the hexadecapole of the power
       *  spectrum
       *
       *  @param IntegrationMethod the integration method
       *
       *  @return the matrix containing the first three non-null
       *  multipoles of the two-point correlation function
       *
       */
      std::vector<std::vector<double>> XiMultipoles (const int nbins, const double rMin, const double rMax, const std::vector<double> kk, const std::vector<double> Pk0, const std::vector<double> Pk2, const std::vector<double> Pk4, const int IntegrationMethod=1);

      /**
       *  @brief the covariance matrix of the first three non-null
       *  multipoles of the two-point correlation function
       *
       *  @param nbins the number of bins of the two-point
       *  correlation function multipoles
       *
       *  @param rMin the minimum scale
       *
       *  @param rMax the maximum scale 
       *
       *  @param nn order of the moment
       *
       *  @param Volume the volume
       *
       *  @param kk vector containing the wave vector modules
       *
       *  @param Pk0 vector containing the monopole of the power
       *  spectrum
       *
       *  @param IntegrationMethod the integration method
       *
       *  @return the covariance matrix of the first three non-null
       *  multipoles of the two-point correlation function
       */
      std::vector<std::vector<double>> XiMonopole_covariance (const int nbins, const double rMin, const double rMax, const double nn, const double Volume, const std::vector<double> kk, const std::vector<double> Pk0, const int IntegrationMethod=1);
       
      /**
       *  @brief the covariance matrix of the first three non-null
       *  multipole moments of the two-point correlation function
       *
       *  @param nbins the number of bins of the two-point
       *  correlation function multipoles
       *
       *  @param rMin the minimum scale
       *
       *  @param rMax the maximum scale 
       *
       *  @param nn order of the moment
       *
       *  @param Volume the volume
       *
       *  @param kk vector containing the wave vector modules
       *
       *  @param Pk0 vector containing the monopole of the power
       *  spectrum
       *
       *  @param Pk2 vector containing the quadrupole of the power
       *  spectrum
       * 
       *  @param Pk4 vector containing the hexadecapole of the power
       *  spectrum
       *
       *  @param IntegrationMethod the integration method
       *
       *  @return the covariance matrix of the first three non-null
       *  multipole moments of the two-point correlation function
       */
      std::vector<std::vector<double>> XiMultipoles_covariance (const int nbins, const double rMin, const double rMax, const double nn, const double Volume, const std::vector<double> kk, const std::vector<double> Pk0, const std::vector<double> Pk2, const std::vector<double> Pk4, const int IntegrationMethod=1);

      ///@}
      /**
       *  @name Functions to estimate the non linear power spectrum
       */
      ///@{

      double F2 (const double k, const double q, const double kq);

      double G2 (const double k, const double q, const double kq);

      double f_k (const double k, const std::shared_ptr<cbl::glob::FuncGrid> PkLin, const double qmin, const double qmax, const double prec=1.e-3);

      double g_k (const double k, const std::shared_ptr<cbl::glob::FuncGrid> PkLin, const double qmin, const double qmax, const double prec=1.e-3);

      double Pk_1loop (const double kk, const std::shared_ptr<cbl::glob::FuncGrid> PkLin, const int corrtype, const double qmin, const double qmax, const double prec=1.e-3);

      /**
       * @brief compute the Delta-Delta non linear power spectrum 
       * at 1-loop following MPTbreeze scheme (Crocce & Scoccimarro 2012)
       *
       * \f[
       *    P_{\delta \delta} (k) = \exp(f(k))^2 (P_L(k)+P_{1loop}(k))
       * \f]
       *
       * where \f$ f(k) \f$ is the second order correction of the non-linear
       * propagator, \f$\P_L(k)\f$ is the linear power spectrum
       * and \f$ P_{1loop} \f$ is the one loop power spectrum correction,
       * computed by cbl::cosmology::Cosmology::Pk_1loop.
       *
       * @param kk the wavevector module
       *
       * @param Pk pointer to a FuncGrid object to interpolate
       * the linear power spectrum
       *
       * @param qmin the lower integration limit
       *
       * @param qmax the upper integration limit
       *
       * @param prec the integral precision
       *
       * @return the Delta-Delta non-linear power spectrum
       */
      double Pk_DeltaDelta (const double kk, const std::shared_ptr<cbl::glob::FuncGrid> Pk, const double qmin, const double qmax, const double prec=1.e-3);

      /**
       * @brief compute the Delta-Delta non linear power spectrum 
       * at 1-loop following MPTbreeze scheme (Crocce & Scoccimarro 2012)
       *
       * \f[
       *    P_{\delta \theta} (k) = \exp(f(k))^2 (P_L(k)+P_{1loop}(k))
       * \f]
       *
       * where \f$ f(k) \f$ is the second order correction of the non-linear
       * propagator, \f$\P_L(k)\f$ is the linear power spectrum
       * and \f$ P_{1loop} \f$ is the one loop power spectrum correction,
       * computed by cbl::cosmology::Cosmology::Pk_1loop.
       *
       * @param kk the wavevector module
       *
       * @param Pk pointer to a FuncGrid object to interpolate
       * the linear power spectrum
       *
       * @param qmin the lower integration limit
       *
       * @param qmax the upper integration limit
       *
       * @param prec the integral precision
       *
       * @return the Delta-Delta non-linear power spectrum
       */
      double Pk_DeltaTheta (const double kk, const std::shared_ptr<cbl::glob::FuncGrid> Pk, const double qmin, const double qmax, const double prec=1.e-3);

      /**
       * @brief compute the Delta-Delta non linear power spectrum 
       * at 1-loop following MPTbreeze scheme (Crocce & Scoccimarro 2012)
       *
       * \f[
       *    P_{\theta \theta} (k) = \exp(g(k))^2 (P_L(k)+P_{1loop}(k))
       * \f]
       *
       * where \f$ g(k) \f$ is the second order correction of the non-linear
       * propagator, \f$\P_L(k)\f$ is the linear power spectrum
       * and \f$ P_{1loop} \f$ is the one loop power spectrum correction,
       * computed by cbl::cosmology::Cosmology::Pk_1loop.
       *
       * @param kk the wavevector module
       *
       * @param Pk pointer to a FuncGrid object to interpolate
       * the linear power spectrum
       *
       * @param qmin the lower integration limit
       *
       * @param qmax the upper integration limit
       *
       * @param prec the integral precision
       *
       * @return the Delta-Delta non-linear power spectrum
       */
      double Pk_ThetaTheta (const double kk, const std::shared_ptr<cbl::glob::FuncGrid> Pk, const double qmin, const double qmax, const double prec=1.e-3);

      std::vector<double> Pk_DeltaDelta (const std::vector<double> kk, const double redshift, const std::string method_Pk, const std::string output_dir, const std::string output_root="test", const int norm=-1, const double k_min=0., const double k_max=100., const double prec=1.e-2, const std::string file_par=par::defaultString, const bool unit1=false);
      
      std::vector<double> Pk_DeltaTheta (const std::vector<double> kk, const double redshift, const std::string method_Pk, const std::string output_dir, const std::string output_root="test", const int norm=-1, const double k_min=0., const double k_max=100., const double prec=1.e-2, const std::string file_par=par::defaultString, const bool unit1=false);

      /**
       *  @brief compute the Delta-Delta non linear power spectrum 
       *  at 1-loop following MPTbreeze scheme (Crocce & Scoccimarro 2012)
       *
       *  \f[ P_{\theta \theta} (k) = \exp(g(k))^2
       *    (P_L(k)+P_{1loop}(k)) \f]
       *
       *  where \f$ g(k) \f$ is the second order correction of the
       *  non-linear propagator, \f$\P_L(k)\f$ is the linear power
       *  spectrum and \f$ P_{1loop} \f$ is the one loop power
       *  spectrum correction, computed by
       *  cbl::cosmology::Cosmology::Pk_1loop.
       *
       *  @param kk vector of wavevector modules
       *
       *  @param redshift the redshift
       * 
       *  @param method_Pk method used to compute the power spectrum;
       *  valid choices for method_Pk are: CAMB [http://camb.info/],
       *  classgal_v1 [http://class-code.net/], MPTbreeze-v1
       *  [http://arxiv.org/abs/1207.1465], EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       *
       *  @param output_dir the output directory
       *
        *  @param norm 0 \f$\rightarrow\f$ don't normalise the power
       *  spectrum; 1 \f$\rightarrow\f$ normalise the power spectrum;
       *  -1 \f$\rightarrow\f$ normalise only if sigma8 is set
       *
       *  @param k_min minimum wave vector module up to which the power
       *  spectrum is computed
       *
       *  @param k_max maximum wave vector module up to which the power
       *  spectrum is computed
       *
       *  @param prec accuracy of the integration 
       *
       *  @param output_root the output_root parameter of the
       *  parameter file used to compute the power spectrum; it can be
       *  any name
       *
       *  @param file_par name of the parameter file; if a parameter
       *  file is provided (i.e. file_par!=NULL), it will use be used,
       *  ignoring the cosmological parameters of the object
       *
       *  @param unit1 true \f$\rightarrow\f$ force cosmological units
       *
       *  @return the \f$\Theta-\Theta\f$ non-linear power spectrum
       */
      std::vector<double> Pk_ThetaTheta (const std::vector<double> kk, const double redshift, const std::string method_Pk, const std::string output_dir, const std::string output_root="test", const int norm=-1, const double k_min=0., const double k_max=100., const double prec=1.e-2, const std::string file_par=par::defaultString, const bool unit1=false);

      ///@}
            
      /**
       *  @name Functions to estimate the three-point correlation function
       */
      ///@{

      /**
       *  @brief the normalization factor for reduced three-point
       *  correlation function
       *
       *  this function computes the normalization factor for reduced
       *  three-point correlation function:
       *
       *  \f[ \xi(r_1)\cdot\xi(r_2) + \xi(r_2)\cdot\xi(r_3) +
       *  \xi(r_3)\cdot\xi(r_1) \f]
       *
       *  with \f$ r_3 = \sqrt{r_1^2+r_2^2-2 r_1 r_2 \cos(\theta)}
       *  \f$
       *
       *  @param r1 the first side of the triangle
       *
       *  @param r2 the second side of the triangle
       *
       *  @param theta the angle between r1 and r2
       *
       *  @param rr vector containing the scales at which the
       *  two-point correlation function is computed
       *  
       *  @param xi_DM vector containing the dark matter two-point
       *  correlation function values, estimated at the scales given
       *  in rr
       *
       *  @return the normalization factor for reduced three-point
       *  correlation function
       */
      double denominator_Q (const double r1, const double r2, const double theta, const std::vector<double> rr, const std::vector<double> xi_DM) const;

      /**
       *  @brief integral functions for the three-point correlation
       *  model
       *
       *  this function computes and store functons used to model the
       *  three-point correlation model; specifically, it implements
       *  Eq. 21, in polar coordinates, of Bel et al. 2015, MNRAS,
       *  453, 259):
       *
       *  \f[ \xi_{DM}(r) = \frac{1}{2\pi^2}\int_0^\infty \mathrm{d}
       *  k\, k^2 P_{DM}(k) j_0(k r), \\\
       *
       *  \Phi(r) = \frac{1}{2\pi^2}\int_0^\infty \mathrm{d} k\,
       *  P_{DM}(k) W^2(kr) j_0(k r), \\ \f]
       *
       *  where \f$j_0(k r)=\sin(k r)/(kr)\f$ is the l=0 spherical
       *  Bessel function, and \f$W(kr)\f$ is the top-hat window
       *  function computed by cbl::TopHat_WF
       *
       *  @param [out] xi_DM vector containing the dark matter
       *  two-point correlation function values
       *
       *  @param [out] Phi vector containing the \f$ \Phi(r)\f$
       *  values, estimated at the scales given in rr
       *
       *  @param [in] rr vector or scales at which the dark matter
       *  two-point correlation function (xi_DM) will be computed
       *
       *  @param [in] kk vector of the wave vector modules at which the
       *  power spectrum is computed
       *
       *  @param [in] Pk_DM vector of containing the dark matter power
       *  spectrum values, estimated at the wave vector modules given
       *  in kk
       *
       *  @param prec the integral precision
       *
       *  @return none
       */
      void integrals_Q_nonLocal (std::vector<double> &xi_DM, std::vector<double> &Phi, const std::vector<double> rr, const std::vector<double> kk, const std::vector<double> Pk_DM, const double prec) const;

      /**
       *  @brief function to compute non-local contribution to
       *  three-point correlation function; specifically, it
       *  implements Eq. 20 of Bel el et al. 2015, MNRAS, 453, 259:
       *
       *  \f[ \Gamma_{123} = \left[
       *  \xi(r_1)+3\frac{\Phi^\prime(r_1)}{r1}\right] \left[
       *  \xi(r_2)+3\frac{\Phi^\prime(r_2)}{r_2}\right]P_2(\cos\theta)
       *  \f]
       *
       *  where the prime indicates the derivative with respect to
       *  \f$r\f$, \f$P_2\f$ is the second Legandre polynomial
       *  computed by cbl::legendre_polynomial, and \f$\xi(r),
       *  \Phi(r)\f$ are the integrals of the power spectrum computed
       *  by cbl::cosmology::Cosmology::integrals_Q_nonLocal
       *
       *  @param r1 the first side of the triangle
       *
       *  @param r2 the second side of the triangle
       *
       *  @param theta the angle betwee r1 and r2
       *
       *  @param xi vector containing the value of xi at r1, r2
       *
       *  @param dPhi vector containing the value of the derivative
       *  of Phi at r1, r2
       *
       *  @return the value of the \f$\Gamma_{123}\f$
       */
      double Gamma_3PCF (const double r1, const double r2, const double theta, const std::vector<double> xi, const std::vector<double> dPhi) const;

      /**
       *  @brief the non-local contribution to the reduced dark
       *  matter three-point correlation function
       *
       *  this function computes the non-local contribution to
       *  three-point correlation function; specifically, it
       *  implements Eq. 22 of Bel el et al. 2015, MNRAS, 453, 259:
       *
       *  \f[ Q_{non-local}(r_1, r_2, \theta) = \frac{2}{3} \left(
       *  \frac{\Gamma_{123} + \Gamma_{312} + \Gamma_{231}}
       *  {\xi(r_1)\cdot\xi(r_2) + \xi(r_2)\cdot\xi(r_3) +
       *  \xi(r_3)\cdot\xi(r_1)}-1 \right) \f]
       *
       *  where the prime indicates the derivative with respect to
       *  \f$r\f$, and \f$\xi(r), \Phi(r)\f$ are the integrals of the
       *  power spectrum computed by
       *  cbl::cosmology::Cosmology::integrals_Q_nonLocal
       *
       *  @param [in] r1 the first side of the triangle
       *
       *  @param [in] r2 the second side of the triangle
       *
       *  @param [in] theta the angle betwee r1 and r2
       *
       *  @param [out] rr vector or scales at which the dark matter
       *  two-point correlation function is computed
       *
       *  @param [out] xi_DM vector containing the dark matter
       *  two-point correlation function values
       *
       *  @param [out] Phi vector containing the \f$ \Phi(r)\f$
       *  values, estimated at the scales given in rr
       *
       *  @param [in] kk vector of the wave vector modules at which
       *  the power spectrum is computed
       *
       *  @param [in] Pk_DM vector of containing the dark matter
       *  power spectrum values, estimated at the wave vector modules
       *  given in kk
       *
       *  @return the value of non-local Q
       */
      double Q_nonLocal (const double r1, const double r2, const double theta, std::vector<double> &rr, std::vector<double> &xi_DM, std::vector<double> &Phi, const std::vector<double> kk, const std::vector<double> Pk_DM) const;

      /**
       *  @brief all the non-local contribution terms of the reduced
       *  dark matter three-point correlation function
       *  
       *  this function computes all the the non-local contribution
       *  terms of the reduced three-point correlation function,
       *  computed by cbl::cosmology::Cosmology::Q_nonLocal
       *
       *  @param r1 the first side of the triangle
       *
       *  @param r2 the second side of the triangle
       *
       *  @param theta the angle betwee r1 and r2
       *
       *  @param kk vector of the wave vector modules at which the
       *  power spectrum is computed
       *
       *  @param Pk_DM vector of containing the dark matter power
       *  spectrum values, estimated at the wave vector modules given
       *  in kk
       *
       *  @return vector containing the DM reduced three-point
       *  correlation function
       */
      std::vector<double> Q_nonLocal (const double r1, const double r2, const std::vector<double> theta, const std::vector<double> kk, const std::vector<double> Pk_DM) const;

      /**
       *  @brief integrals used to compute the Slepian et al. 2015
       *  three-point correlation function model
       *
       *  this function computes the integrals used to
       *  model the three-point correlation as described in Slepian
       *  et. al 2015:
       *
       *  \f[ \xi_{DM}(r) = \frac{1}{2\pi^2}\int_0^\infty
       *  \mathrm{d} k k^2 P_{DM}(k) j_0(k r), \\\
       *  \xi^{[1\pm]}_{DM}(r) = \frac{1}{2\pi^2}\int_0^\infty
       *  \mathrm{d} k k^2 P_{DM}(k) k^{\pm 1} j_1(k r), \\
       *  \xi^{[2]}_{DM}(r) = \frac{1}{2\pi^2}\int_0^\infty
       *  \mathrm{d} k k^2 P_{DM}(k) j_2(k r) . \f]
       *
       *  @param [out] xi_DM vector containing the dark matter
       *  two-point correlation function
       *
       *  @param [out] xi_DM_m1 vector containing
       *  \f$\xi^{[1-]}_{DM}(r)\f$
       *
       *  @param [out] xi_DM_p1 vector containing
       *  \f$\xi^{[1+]}_{DM}(r)\f$
       *
       *  @param [out] xi_DM_2 vector containing
       *  \f$\xi^{[2]}_{DM}(r)\f$
       *
       *  @param [in] rr vector or scales
       *
       *  @param [in] kk vector of the wave vector modules
       *
       *  @param [in] Pk_DM the dark matter power spectrum
       *
       *  @return none
       */
      void integrals_zeta_Slepian (std::vector<double> &xi_DM, std::vector<double> &xi_DM_m1, std::vector<double> &xi_DM_p1, std::vector<double> &xi_DM_2, const std::vector<double> rr, const std::vector<double> kk, const std::vector<double> Pk_DM) const;

      /**
       *  @brief the pre-cyclic three-point correlation function as
       *  described in Slepian et al. 2015
       *
       *  this function computes the pre-cyclic three-point
       *  correlation function as described in Slepian et al. 2015, as
       *  follows:
       *
       *  \f[ \zeta_{pc} = \sum_{l=0}^2 \zeta_{pc l}(r_1, r_2)
       *		   P_l(\hat{r_1}\cdot \hat{r_2})+ \sum_{l=0}^2
       *		   \zeta_{pc l}(r_2, r_3) P_l(\hat{r_2}\cdot
       *		   \hat{r_3})+ \sum_{l=0}^2 \zeta_{pc l}(r_3,
       *		   r_1) P_l(\hat{r_3}\cdot \hat{r_1}) \f]
       *
       *  with 
       *
       *  \f[r_3 = \sqrt{r_1^2+r_2^2-2 r_1 r_2 \cos(\theta)}\f]
       *
       *  and
       *
       *  \f[ \zeta_{pc0}(r_i, r_j) = \left[ 2 b_1^2 b_2 +
       *  \frac{34}{21} b_1^3 \right] \xi(r_i) \xi(r_j), \f] \f[
       *  \zeta_{pc1}(r_i, r_j) = -b_1^3\left[ \xi^{[1-]}(r_i)
       *  \xi^{[1+]}(r_j) + \xi^{[1-]}(r_j) \xi^{[1+]}(r_i) \right] ,
       *  \f] \f[ \zeta_{pc2}(r_i, r_j) = \frac{8}{21}
       *  b_1^3\xi^{[2]}(r_i)\xi^{[2]}(r_j) . \f]
       *
       *  where \f$ b_1, b_2 \f$ are the linear and non-linear bias,
       *  respectively, and \f$\xi_{DM}(r), \xi^{[1\pm]}_{DM}(r),
       *  \xi^{[2]}_{DM}(r)\f$ are the integrals of the dark matter
       *  power spectrum computed by
       *  cbl::cosmology::Cosmology::integrals_zeta_Slepian
       *
       *  @param r1 the first side
       *
       *  @param r2 the second side
       *
       *  @param mu the cosine of the angle between r1 and r2
       *
       *  @param b1 the linear bias
       *
       *  @param b2 the non-linear bias
       *
       *  @param interp_xi_DM interpolating function
       *  for \f$\xi_DM\f$
       *
       *  @param interp_xi_DM_m1 interpolating function
       *  for \f$\xi^{[1-]}_{DM}(r)\f$      
       *
       *  @param interp_xi_DM_p1 interpolating function
       *  for \f$\xi^{[1+]}_{DM}(r)\f$
       *  
       *  @param interp_xi_DM_2 interpolating function
       *  for \f$\xi^{[2]}_{DM}(r)\f$
       *
       *  @return the pre-cyclic three-point
       *  correlation function
       */
      double zeta_precyclic_Slepian (const double r1, const double r2, const double mu, const double b1, const double b2, const glob::FuncGrid interp_xi_DM, const glob::FuncGrid interp_xi_DM_m1, const glob::FuncGrid interp_xi_DM_p1, const glob::FuncGrid interp_xi_DM_2) const;

      /**
       *  @brief the terms of the \f$\zeta(r_1, r_2)\f$ expansion 
       *
       *  this function computes the terms of the \f$\zeta(r_1,
       *  r_2)\f$ expansion up to an arbitrary order \f$l\f$ (the
       *  default value is \f$l_{max}=9\f$), as described in Slepian
       *  et al. 2015:
       * 
       *  \f[\zeta_l(r_1, r_2) = \frac{2l+1}{2} \int_{-1}^{1}
       *  \mathrm{d}\mu_{12} \left[\zeta_{pc}(r_1, r_2,
       *  \mu_{12})+\zeta_{pc}(r_2, r_3, \mu_{23})+ \zeta_{pc}(r_3,
       *  r_1, \mu_{31})\right] P_l(\mu_{12}) .\f] 
       * 
       *  the terms in square brackets is computed by
       *  cbl::cosmology::Cosmology::zeta_precyclic_Slepian
       *
       *  @param [in] r1 the first side of the triangle
       *
       *  @param [in] r2 the second side of the triangle
       *
       *  @param [in] b1 the linear bias of the triangle
       *
       *  @param [in] b2 the non-linear bias
       *
       *  @param [out] rr vector or scales
       *
       *  @param [out] xi_DM vector containing the dark matter
       *  two-point correlation function
       *
       *  @param [out] xi_DM_m1 vector containing
       *  \f$\xi^{[1-]}_{DM}(r)\f$
       *
       *  @param [out] xi_DM_p1 vector containing
       *  \f$\xi^{[1+]}_{DM}(r)\f$
       *
       *  @param [out] xi_DM_2 vector containing
       *  \f$\xi^{[2]}_{DM}(r)\f$
       *
       *  @param [in] norders the maximum numbers of orders
       *
       *  @param [in] prec the integral precision
       *
       *  @return vector containing the terms of legendre expansion
       */
      std::vector<double> zeta_expansion_Slepian (const double r1, const double r2, const double b1, const double b2, std::vector<double> &rr, std::vector<double> &xi_DM, std::vector<double> &xi_DM_m1, std::vector<double> &xi_DM_p1, std::vector<double> &xi_DM_2, const int norders=9, const double prec=1.e-3) const;

      /**
       *  @brief the dark matter three-point correlation function
       *  model by Slepian et al. 2015
       *
       *  this function computes \f$\zeta_{DM} (r_1, r_2, \hat{r_1}
       *  \cdot \hat{r_2})\f$, as described in Slepian et al. 2015:
       *
       *  \f[ \zeta_{DM} (r_1, r_2, \hat{r_1} \cdot \hat{r_2}) = 
       *  \sum_l \zeta_l(r_1, r_2) P_l(\hat{r_1} \cdot \hat{r_2}) .\f]
       *
       *  The coefficients of the expansion are computed by
       *  cbl::cosmology::Cosmology::zeta_expansion_Slepian
       *
       *  @param [in] r1 the first side of the triangle
       *
       *  @param [in] r2 the second side of the triangle
       *
       *  @param [in] theta the angle between r1 and r2
       *
       *  @param [out] rr vector or scales
       *
       *  @param [out] xi_DM vector containing the dark matter
       *  two-point correlation function
       *
       *  @param [out] xi_DM_m1 vector containing
       *  \f$\xi^{[1-]}_{DM}(r)\f$
       *
       *  @param [out] xi_DM_p1 vector containing
       *  \f$\xi^{[1+]}_{DM}(r)\f$
       *
       *  @param [out] xi_DM_2 vector containing
       *  \f$\xi^{[2]}_{DM}(r)\f$
       *
       *  @param [in] kk vector of the wave vector modules
       *
       *  @param [in] Pk_DM the dark matter power spectrum
       *
       *  @param [in] norders the maximum number of orders
       *
       *  @param [in] prec the integral precision
       *
       *  @return the connected dark matter three-point correlation
       *  function
       */
      double zeta_DM_Slepian (const double r1, const double r2, const double theta, std::vector<double> &rr, std::vector<double> &xi_DM, std::vector<double> &xi_DM_m1, std::vector<double> &xi_DM_p1, std::vector<double> &xi_DM_2, const std::vector<double> kk, const std::vector<double> Pk_DM, const int norders=9, const double prec=1.e-3) const;

      /**
       *  @brief the dark matter reduced three-point correlation
       *  function model by Slepian et al. 2015
       *
       *  this function computes \f$Q_{DM} (r_1, r_2, \hat{r_1} \cdot
       *  \hat{r_2})\f$ as described in Slepian et al. 2015:
       *
       *  \f[ Q_{DM} (r_1, r_2, \hat{r_1} \cdot \hat{r_2}) =
       *  \frac{\zeta_{DM} (r_1, r_2, \hat{r_1} \cdot \hat{r_2})}
       *  {\left(
       *  \xi(r_1)\xi(r_2)+\xi(r_2)\xi(r_3)+\xi(r_3)\xi(r_1)\right)}
       *  \f]
       *
       *  see cbl::cosmology::Cosmology::zeta_DM_Slepian for
       *  more details
       *
       *  @param [in] r1 the first side of the triangle
       *
       *  @param [in] r2 the second side of the triangle
       *
       *  @param [in] theta the angle between r1 and r2
       *
       *  @param [out] rr vector or scales
       *
       *  @param [out] xi_DM vector containing the dark matter
       *  two-point orrelation function
       *
       *  @param [out] xi_DM_m1 vector containing
       *  \f$\xi^{[1-]}_{DM}(r)\f$
       *
       *  @param [out] xi_DM_p1 vector containing
       *  \f$\xi^{[1+]}_{DM}(r)\f$
       *
       *  @param [out] xi_DM_2 vector containing
       *  \f$\xi^{[2]}_{DM}(r)\f$
       *
       *  @param [out] kk vector of the wave vector modules
       *
       *  @param [out] Pk_DM the dark matter power spectrum
       *
       *  @param [in] norders the maximum numbers of orders
       *
       *  @param [in] prec the integral precision
       *
       *  @return the dark matter reduced three-point correlation
       *  function
       */
      double Q_DM_Slepian (const double r1, const double r2, const double theta, std::vector<double> &rr, std::vector<double> &xi_DM, std::vector<double> &xi_DM_m1, std::vector<double> &xi_DM_p1, std::vector<double> &xi_DM_2, const std::vector<double> kk, const std::vector<double> Pk_DM, const int norders=9, const double prec=1.e-3) const;

      /**
       *  @brief integrals used to compute the Barriga & Gatzanaga
       *  al. 2002 three-point correlation function model
       * 
       *  this function computes the integrals used to model the
       *  three-point correlation as described in Barriga & Gatzanaga
       *  2002:
       *
       *  \f[ \xi_{DM}(r) = \frac{1}{2\pi^2}\int_0^\infty \mathrm{d}
       *  k\, k^2 P_{DM}(k) j_0(k r), \f] \f[ \Phi(r) =
       *  \frac{1}{2\pi^2}\int_0^\infty \mathrm{d} k\, P_{DM}(k)
       *  j_0(k r). \f]
       *
       *  @param [out] xi_DM vector containing the dark matter
       *  two-point correlation function
       *
       *  @param [out] Phi vector containing \f$ \Phi(r)\f$
       *
       *  @param [in] rr vector or scales
       *
       *  @param [in] kk vector of the wave vector modules
       *
       *  @param [in] Pk_DM the dark matter power spectrum
       *
       *  @return none
       */
      void integrals_zeta_BarrigaGatzanaga (std::vector<double> &xi_DM, std::vector<double> &Phi, const std::vector<double> rr, const std::vector<double> kk, const std::vector<double> Pk_DM) const;

      /**
       *  @brief the single term of the dark matter three-point
       *  correlation function model by Barriga & Gatzanaga et
       *  al. 2002
       *
       *  this function computes the single term of the dark matter
       *  three-point correlation function, following Barriga &
       *  Gatzanaga et al.  2002:
       *
       *  \f[ f(r_1, r_2) = \frac{10}{7}\xi(r_1) \xi(r_2)+\frac{4}{7}
       *  \left\{ -3 \frac{\Phi^\prime(r_1) \Phi^\prime(r_2)}{r_1
       *  r_2} -\frac{\xi(r_1) \Phi^\prime(r_2)}{r_2}-\frac{\xi(r_2)
       *  \Phi^\prime(r_1)}{r_1} +\mu^2\left[
       *  \xi(r_1)+3\frac{\Phi^\prime(r_1)}{r1}\right]\left[
       *  \xi(r_2)+3\frac{\Phi^\prime(r_2)}{r_3}\right] \right\}
       *  -\mu\left[ \xi^\prime(r_1)\Phi^\prime(r_2) +
       *  \xi^\prime(r_2)\Phi^\prime(r_1)\right] \f]
       *
       *  where the prime indicates the derivative with respect to
       *  \f$r\f$, and \f$\xi(r), \Phi(r)\f$ are the integrals of the
       *  power spectrum computed by
       *  cbl::cosmology::Cosmology::integrals_zeta_BarrigaGatzanaga
       *
       *  @param r1 the first side of the triangle
       *
       *  @param r2 the second side of the triangle
       *
       *  @param theta the angle betwee r1 and r2
       *
       *  @param xi vector containing the value of xi at r1, r2
       *
       *  @param dxi vector containing the value of the derivative of
       *  xi at r1, r2
       *
       *  @param dPhi vector containing the value of the derivative
       *  of Phi at r1, r2
       *
       *  @return the dark matter reduced three-point correlation
       *  function
       */
      double zeta_single_BarrigaGatzanaga (const double r1, const double r2, const double theta, const std::vector<double> xi, const std::vector<double> dxi, const std::vector<double> dPhi) const;

      /**
       *  @brief the dark matter three-point correlation function
       *  model by Barriga & Gatzanaga et al. 2002
       *
       *  this functions computes the dark matter three-point
       *  correlation function model by Barriga & Gatzanaga et al
       *  2002:
       *
       *  \f[ f(r_1, r_2) = \frac{10}{7}\xi(r_1) \xi(r_2)+\frac{4}{7}
       *  \left\{ -3 \frac{\Phi^\prime(r_1) \Phi^\prime(r_2)}{r_1
       *  r_2} -\frac{\xi(r1) \Phi^\prime(r_2)}{r_2}-\frac{\xi(r2)
       *  \Phi^\prime(r_1)}{r_1} +\mu^2\left[
       *  \xi(r_1)+3\frac{\Phi^\prime(r_1)}{r1}\right]\left[
       *  \xi(r_2)+3\frac{\Phi^\prime(r_2)}{r_3}\right] \right\}
       *  -\mu\left[ \xi^\prime(r_1)\Phi^\prime(r_2) +
       *  \xi^\prime(r_2)\Phi^\prime(r_1)\right] +
       *  \mathrm{permutations} \f]
       *
       *  where the prime indicates the derivative with respect to
       *  \f$r\f$, and \f$\xi(r), \Phi(r)\f$ are the integrals of the
       *  power spectrum computed by
       *  cbl::cosmology::Cosmology::integrals_zeta_BarrigaGatzanaga.
       *
       *  @param [in] r1 the first side of the triangle
       *
       *  @param [in] r2 the second side of the triangle
       *
       *  @param [in] theta the angle between r1 and r2
       *
       *  @param [out] rr vector or scales
       *
       *  @param [out] xi_DM vector containing the dark matter
       *  two-point correlation function
       *
       *  @param [out] Phi vector containing \f$ \Phi(r)\f$
       *
       *  @param [in] kk vector of the wave vector modules
       *
       *  @param [in] Pk_DM the dark matter power spectrum
       *
       *  @return the dark matter three-point correlation function
       */
      double zeta_DM_BarrigaGatzanaga (const double r1, const double r2, const double theta, std::vector<double> &rr, std::vector<double> &xi_DM, std::vector<double> &Phi, const std::vector<double> kk, const std::vector<double> Pk_DM) const;

      /**
       *  @brief the dark matter reduced three-point correlation
       *  function model by Barriga & Gatzanaga et al. 2002
       *
       *  this functions computes \f$Q_{DM} (r_1, r_2, \hat{r_1}
       *  \cdot \hat{r_2})\f$, as described in Barriga & Gatzanaga et
       *  al. 2002:
       *
       *  \f[ Q_{DM} (r_1, r_2, \hat{r_1} \cdot \hat{r_2}) =
       *  \frac{\zeta_{DM} (r_1, r_2, \hat{r_1} \cdot \hat{r_2})}
       *  {\left(
       *  \xi(r_1)\xi(r_2)+\xi(r_2)\xi(r_3)+\xi(r_3)\xi(r_1)\right)}
       *  \f]
       *
       *  see cbl::cosmology::Cosmology::zeta_DM_BarrigaGatzanaga
       *  for more details
       *
       *  @param [in] r1 the first side of the triangle
       *
       *  @param [in] r2 the second side of the triangle
       *
       *  @param [in] theta the angle between r1 and r2
       * 
       *  @param [out] rr vector or scales
       * 
       *  @param [out] xi_DM vector containing the dark matter
       *  two-point correlation function
       * 
       *  @param [out] Phi vector containing \f$ \Phi(r)\f$
       * 
       *  @param [in] kk vector of the wave vector modules
       * 
       *  @param [in] Pk_DM the dark matter power spectrum
       *
       *  @return the dark matter reduced three-point correlation
       *  function
       */
      double Q_DM_BarrigaGatzanaga (const double r1, const double r2, const double theta, std::vector<double> &rr, std::vector<double> &xi_DM, std::vector<double> &Phi, const std::vector<double> kk, const std::vector<double> Pk_DM) const;

      /**
       *  @brief the dark matter three-point correlation function
       *
       *  this function computes the dark matter three-point
       *  correlation function with either the Slepian et al 2015 or
       *  the Barriga & Gatzagnaga 2002 model
       *
       *  @param r1 the first side of the triangle
       *
       *  @param r2 the second side of the triangle
       *
       *  @param theta vector containing angles between r1 and r2, in
       *  radians
       *
       *  @param model the model to compute the three-point
       *  correlation function, can be "Slepian" or
       *  "BarrigaGatzanaga"
       *
       *  @param kk vector of the wave vector modules
       *
       *  @param Pk_DM the dark matter power spectrum
       *
       *  @return vector containing the dark matter three-point
       *  correlation function
       */
      std::vector<double> zeta_DM (const double r1, const double r2, const std::vector<double> theta, const std::string model, const std::vector<double> kk, const std::vector<double> Pk_DM) const;

      /**
       *  @brief the dark matter reduced three-point correlation
       *  function
       *
       *  this function computes the dark matter reduced three-point
       *  reduced correlation function with either the Slepian et al.
       *  2015 or the Barriga & Gatzagnaga 2002 model
       *
       *  @param r1 the first side of the triangle 
       *
       *  @param r2 the second side of the triangle
       *
       *  @param theta vector containing angles between r1 and r2, in
       *  radians
       *
       *  @param model the model to compute the three-point
       *  correlation function, can be "Slepian" or
       *  "BarrigaGatzanaga"
       *
       *  @param kk vector of the wave vector modules
       *
       *  @param Pk_DM the dark matter power spectrum
       *
       *  @return vector containing the dark matter reduced
       *  three-point correlation function
       */
      std::vector<double> Q_DM (const double r1, const double r2, const std::vector<double> theta, const std::string model, const std::vector<double> kk, const std::vector<double> Pk_DM) const;

      /**
       *  @brief the local-bias model of the three-point correlation
       *  function of dark matter haloes
       *
       *  this function computes the three-point correlation function
       *  of dark matter haloes with either the Slepian et al.  2015
       *  or the Barriga & Gatzagnaga 2002 model, as follows:
       *
       *  \f[ \zeta_h (r_1, r_2, \hat{r_1} \cdot \hat{r_2}) = b_1^3
       *  \zeta_{DM}(r_1, r_2, \hat{r_1} \cdot \hat{r_2}) + b_1^2 b_2
       *  \left[ \xi(r_1)\cdot\xi(r_2) + \xi(r_2)\cdot\xi(r_3) +
       *  \xi(r_3)\cdot\xi(r_1) \right] \f]
       *
       *  with \f$r_3 = \sqrt{r_1^2+r_2^2-2 r_1 r_2 \cos(\theta)}\f$
       *  and \f$b_1, b_2\f$ the linear and non-linear halo bias,
       *  respectively; \f$\zeta_{DM}\f$ is compute by
       *  cbl::cosmology::Cosmology::zeta_DM
       *
       *  @param r1 the first side of the triangle
       *
       *  @param r2 the second side of the triangle
       *  
       *  @param theta vector containing angles between r1 and r2, in
       *  radians 
       *
       *  @param b1 the linear bias
       *
       *  @param b2 the non-linear bias
       *
       *  @param model the model to compute the three-point
       *  correlation function, can be "Slepian" or
       *  "BarrigaGatzanaga"
       *
       *  @param kk vector of the wave vector modules
       *
       *  @param Pk_DM the dark matter power spectrum
       *
       *  @return vector containing the three-point correlation
       *  function of dark matter haloes
       */
      std::vector<double> zeta_halo (const double r1, const double r2, const std::vector<double> theta, const double b1, const double b2, const std::string model, const std::vector<double> kk, const std::vector<double> Pk_DM) const;

      /**
       *  @brief the local-bias model of the reduced three-point
       *  correlation function of dark matter haloes
       *
       *  this function computes the reduced three-point correlation
       *  function of dark matter haloes with either the Slepian et
       *  al.  2015 or the Barriga & Gatzagnaga 2002 model, as
       *  follows:
       *
       *  \f[ Q_h (r_1, r_2, \hat{r_1} \cdot \hat{r_2}) =
       *  \frac{Q_{DM}(r_1, r_2, \hat{r_1} \cdot
       *  \hat{r_2})}{b_1}+\frac{b_2}{b_1^2}\f]
       *
       *  \f$Q_{DM}\f$ is compute by
       *  cbl::cosmology::Cosmology::Q_DM
       *
       *  @param r1 the first side of the triangle
       *
       *  @param r2 the second side of the triangle
       *
       *  @param theta vector containing angles between r1 and r2, in
       *  radians
       *
       *  @param b1 the linear bias
       *
       *  @param b2 the non-linear bias
       *
       *  @param model the model to compute the three-point
       *  correlation function, can be "Slepian" or
       *  "BarrigaGatzanaga"
       *
       *  @param kk vector of the wave vector modules
       *
       *  @param Pk_DM the dark matter power spectrum
       *
       *  @return vector containing the reduced three-point
       *  correlation function of dark matter haloes
       */
      std::vector<double> Q_halo (const double r1, const double r2, const std::vector<double> theta, const double b1, const double b2, const std::string model, const std::vector<double> kk, const std::vector<double> Pk_DM) const;

      /**
       *  @brief the non-local-bias model of the three-point
       *  correlation function of dark matter haloes
       *
       *  this function computes the reduced three-point correlation
       *  function of dark matter haloes, with non-local bias
       *  corrections, with either the Slepian et al. 2015 or the
       *  Barriga & Gatzagnaga 2002 model, as follows:
       *
       *  \f[ Q_h (r_1, r_2, \hat{r_1} \cdot \hat{r_2}) =
       *  \frac{Q_{DM}(r_1, r_2, \hat{r_1} \cdot \hat{r_2})}{b_1}+
       *  \frac{b_2}{b_1^2}+\frac{g_2}{b_1}Q_{non-local} \f]
       *
       *  \f$Q_{DM}\f$ is compute by
       *  cbl::cosmology::Cosmology::Q_DM and \f$Q_{non-local}\f$
       *  is the non-local contirbuion term, computed by
       *  cbl::cosmology::Cosmology::Q_nonLocal
       *
       *  @param r1 the first side of the triangle
       *
       *  @param r2 the second side of the triangle
       *
       *  @param theta vector containing angles between r1 and r2, in
       *  radians
       *
       *  @param b1 the linear bias
       *
       *  @param b2 the non-linear bias
       *
       *  @param g2 the non-local bias
       *
       *  @param model the model to compute the three-point
       *  correlation function, can be "Slepian" or
       *  "BarrigaGatzanaga"
       *
       *  @param kk vector of the wave vector modules 
       *
       *  @param Pk_DM the dark matter power spectrum
       *
       *  @return vector containing the reduced three-point
       *  correlation function of dark matter haloes
       */
      std::vector<double> Q_halo (const double r1, const double r2, const std::vector<double> theta, const double b1, const double b2, const double g2, const std::string model, const std::vector<double> kk, const std::vector<double> Pk_DM) const;

      /**
       *  @brief the dark matter equilateral three-point correlation
       *  function
       *
       *  this function computes the dark matter equilateral
       *  three-point correlation function with either the Slepian et
       *  al 2015 or the Barriga & Gatzagnaga 2002 model
       *
       *  @param rr vector of sides
       *
       *  @param model the model to compute the three-point
       *  correlation function, can be "Slepian" or "BarrigaGatzanaga"
       *
       *  @param kk vector of the wave vector modules
       *
       *  @param Pk_DM the dark matter power spectrum
       *
       *  @return vector containing the dark matter three-point
       *  correlation function
       */
      std::vector<double> zeta_DM_eq (const std::vector<double> rr, const std::string model, const std::vector<double> kk, const std::vector<double> Pk_DM) const;

      /**
       *  @brief the dark matter equilateral reduced three-point
       *  correlation function
       *
       *  this function computes the dark matter equilateral reduced
       *  three-point correlation function with either the Slepian et
       *  al 2015 or the Barriga & Gatzagnaga 2002 model
       *
       *  @param rr vector of sides
       *
       *  @param model the model to compute the three-point
       *  correlation function, can be "Slepian" or
       *  "BarrigaGatzanaga"
       *
       *  @param kk vector of the wave vector modules
       *
       *  @param Pk_DM the dark matter power spectrum
       *
       *  @return vector containing the dark matter three-point
       *  correlation function
       */
      std::vector<double> Q_DM_eq (const std::vector<double> rr, const std::string model, const std::vector<double> kk, const std::vector<double> Pk_DM) const;

      /**
       *  @brief the dark matter three-point correlation function 
       *  multipoles covariance model, by Slepian et al. 2015
       *
       *  \f[
       *  C_{{\rm GRF}, ll'}(r_{1},r_{2};r_{1}',r_{2}')=\frac{4\pi}{V}(2l+1)(2l'+1)(-1)^{l+l'} \times\int r^{2}dr\sum_{l_{2}}(2l_{2}+1)
       *  \left(\begin{array}{ccc} 
       *  l & l' & l_{2}\\
       *  0 & 0 & 0
       *   \end{array}\right)^2\nonumber\\
       *   \times\bigg\{(-1)^{l_2}\xi(r)\bigg[f_{l_{2}ll'}(r;r_{1},r_{1}')f_{l_{2}ll'}(r;r_{2},r_{2}')
       *   +f_{l_{2}ll'}(r;r_{2},r_{1}')f_{l_{2}ll'}(r;r_{1},r_{2}')\bigg]+(-1)^{(l+l'+l_{2})/2}
       *   \times\bigg[f_{ll}(r;r_{1})f_{l'l'}(r;r_{1}')f_{l_{2}ll'}(r;r_{2},r_{2}')
       *   +f_{ll}(r;r_{1})f_{l'l'}(r;r_{2}')f_{l_{2}ll'}(r;r_{2},r_{1}')
       *   +f_{ll}(r;r_{2})f_{l'l'}(r;r_{1}')f_{l_{2}ll'}(r;r_{1},r_{2}')
       *   +f_{ll}(r;r_{2})f_{l'l'}(r;r_{2}')f_{l_{2}ll'}(r;r_{1},r_{1}')\bigg]\bigg\}
       *  \f]
       *
       *  with:
       *  \f[
       *   f_{ll}(r;r_{1})=\int\frac{k^{2}dk}{2\pi^{2}}\left[P(k)+\frac{1}{n}\right]j_{l}(kr_{1})j_{l}(kr)
       *  \f]
       *   and:
       *   \f[
       *   f_{l_{2}ll'}(r;r_{1},r_{1}')=\int\frac{k^{2}dk}{2\pi^{2}} \left[P(k)+\frac{1}{n}\right] j_{l}(kr_{1})j_{l'}(kr_{1}')j_{l_{2}}(kr),
       *   \f]
       *   where \f$V\f$ is the effective survey volume, and  \f$n\f$ is the effective number density of the survey.
       *
       *  @param Volume the volume
       *
       *  @param nObjects the number of objects
       *
       *  @param l the order l of the multipoles expansion
       *
       *  @param l_prime the order \f$l\textprime\f$ of the multipoles expansion
       *
       *  @param r1 the scale \f$r_1\f$
       *  
       *  @param r2 the scale \f$r_2\f$
       *  
       *  @param r1_prime the scale \f$r_1\textprime\f$
       *
       *  @param r2_prime the scale \f$r_2\textprime\f$
       *
       *  @param kk vector of the wave vector modules
       *
       *  @param Pk the pdark matter ower spectrum
       *
       *  @param rr vector of scales
       *
       *  @param Xi vector containing the two-point correlation
       *  function
       *
       *  @param prec the integral precision
       *
       *  @return the covariance of the multipole expansion of 
       *  the dark matter three-point correlation function
       */
      double zeta_multipoles_covariance (const double Volume, const double nObjects, const int l, const int l_prime, const double r1, const double r2, const double r1_prime, const double r2_prime, const std::vector<double> kk, const std::vector<double> Pk, const std::vector<double> rr, const std::vector<double> Xi, const double prec=1.e-3);

      /**
       *  @brief the dark matter three-point correlation function 
       *  covariance model
       *
       *  this function computes the dark matter three-point
       *  correlation function covariance model, by Slepian et
       *  al. 2015, as a function of \f$\theta = r_1 \cdot r_2\f$:
       *
       *  \f[ C(r_1, r_2, \vec{r_1}\cdot\vec{r_2} \equiv \cos(\theta)
       *  = \sum_{l=0}^{l=max_l} \sum_{l^'=0}^{l^'=max_l} C_{l,
       *  l^'}(r_1, r_2) P_l(\cos(\theta) P_{l^'}(\cos(theta) \f]
       *
       *  where \f$C_{l, l^'}(r_1, r_2)\f$ is computed by
       *  cbl::cosmology::Cosmology::zeta_multipoles_covariance
       *
       *  @param Volume the volume
       *
       *  @param nObjects the number of objects
       *
       *  @param theta vector of angles at which the covariance is
       *  computed
       *
       *  @param r1 the scale \f$r_1\f$
       *  
       *  @param r2 the scale \f$r_2\f$
       *
       *  @param kk vector of the wave vector modules
       *
       *  @param Pk the pdark matter ower spectrum
       *
       *  @param norders the maximum number of orders of multipoles
       *  of the three point correlation function expansion
       *
       *  @param prec the integral precision
       *
       *  @param method false \f$\rightarrow\f$ apply method 1; true
       *  \f$\rightarrow\f$ apply method 2
       *
       *  @param nExtractions the number of mock extraction
       *  from zeta multipoles coefficient covariance matrix
       *
       *  @param mean vector containing the mean values
       *
       *  @param seed random number generator seed
       *
       *  @return the covariance of the dark matter 
       *   three-point correlation function
       */
      std::vector<std::vector<double>> zeta_covariance (const double Volume, const double nObjects, const std::vector<double> theta, const double r1, const double r2, const std::vector<double> kk, const std::vector<double> Pk, const int norders=10, const double prec=1.e-3, const bool method=false, const int nExtractions=10000, const std::vector<double> mean={}, const int seed=543);
       
      /**
       * @brief compute the  power spectrum integral transform
       *
       * this function computes the power spectrum integral transform:
       *
       * \f[ 
       *   \xi^{[n]} (r) = \int \frac{k^2\mathrm{d}k}{2\pi^2} P(k) j_n(kr).
       * \f]
       *
       * where n is the order of the transform.
       *
       * @param xi_n the power spectrum transform \f$\xi^{[n]} (r)\f$
       *
       * @param rr vector of scales
       *
       * @param nn the order of the transform
       *
       * @param kk vector of wavevector modules
       *
       * @param Pk dark matter power spectrum
       *
       * @return none
       */
      void xi_r_n (std::vector<double> &xi_n, const std::vector<double> rr, const int nn, const std::vector<double> kk, const std::vector<double> Pk);

      /**
       * @brief compute the  power spectrum integral transform
       *
       * this function computes the power spectrum integral transform:
       *
       * \f[ 
       *   \xi^{[n\pm]} (r) = \int \frac{k^2\mathrm{d}k}{2\pi^2} k^{\pm1} P(k) j_n(kr)
       * \f]
       *
       * where n is the order of the transform.
       *
       * @param xi_n_p the power spectrum transform \f$\xi^{[n+]} (r)\f$
       * 
       * @param xi_n_m the power spectrum transform \f$\xi^{[n-]} (r)\f$
       *
       * @param rr vector of scales
       *
       * @param nn the order of the transform
       *
       * @param kk vector of wavevector modules
       *
       * @param Pk dark matter power spectrum
       *
       * @return none
       */
      void xi_r_n_pm (std::vector<double> &xi_n_p, std::vector<double> &xi_n_m, const std::vector<double> rr, const int nn, const std::vector<double> kk, const std::vector<double> Pk);

      /**
       * @brief compute the  power spectrum integral transform
       *
       * this function computes the power spectrum integral transform:
       *
       * \f[ 
       *    f_{l, l_1} (r_i;r)= \int \frac{k^2\mathrm{d}k}{2\pi^2} j_l(kr_i) j_{l_1} (kr) k P(k)
       * \f]
       *
       * where \f$l, l_1\f$ are the orders of the transform.
       *
       * @param eff the power spectrum transform \f$ f_{l, l_1} (r_i;r) \f$
       * 
       * @param rr vector of scales
       *
       * @param l the order \f$l\f$
       *
       * @param l1 the order \f$l_1\f$
       *
       * @param kk vector of wavevector modules
       *
       * @param Pk dark matter power spectrum
       *
       * @return none
       */
      void eff_l_l1 (std::vector<std::vector<double>> &eff, const std::vector<double> rr, const int l, const int l1, const std::vector<double> kk, const std::vector<double> Pk);

      /**
       * @brief compute the quantity \f$ I_{\mathcal{L} l} (r_1, r_2)\f$
       *
       * This function computes the quantity \f$ I_{\mathcal{L} l} (r_1, r_2)\f$:
       *
       * \f[
       *   I_{\mathcal{L} l} (r_1, r_2) = \sum_l_1 (-1)^{l_1+l}(2l_1+1)(2l+1) \begin{pmatrix} l_1 & l & \mathcal{L} \\ 0 & 0 & 0 \end{pmatrix}^2 
       *   \\ \times \int r \mathrm{d} r f_{l, l_1}(r_1;r) f_{l, l_1} (r_2; r)
       * \f]
       *
       * where \f$ f_{l, l_1} (r_i;r) \f$ is computed by cbl::cosmology::Cosmology::eff_l_l1
       * This quantity is  used the compute the tree-level theoretical 
       * prediction  for the biased and redshift space the three-point correlation function, 
       * following Slepian&Eisenstein, 2017
       *
       * @param II the quantity \f$ I_{\mathcal{L} l} (r_1, r_2) \f$
       * 
       * @param rr vector of scales
       *
       * @param ll the order \f$l\f$
       *
       * @param LL the order \f$ \mathcal{L} \f$
       *
       * @param kk vector of wavevector modules
       *
       * @param Pk dark matter power spectrum
       *
       * @return none
       */
      void I_ELL_ell (std::vector<std::vector<double>> &II, const std::vector<double> rr, const int ll, const int LL, const std::vector<double> kk, const std::vector<double> Pk);

      /**
       * @brief compute the quantity \f$ k_l (r_1, r_2) \f$
       *
       * This function computes the quantity \f$ k_l (r_1, r_2) \f$:
       *
       * \f[
       *  k_l (r_1, r_2) = \frac{64}{77175} \left[9I_{1,l}(r_1, r_2)-
       *  14I_{3,l}(r_1, r_2) +5I_{5,l}(r_1, r_2) \right] 
       * \f]
       *
       * where \f$ I_{\mathcal{L}l} \f$ is computed by cbl::cosmology::Cosmology::I_ELL_ell
       * This quantity is  used the compute the tree-level theoretical 
       * prediction for the biased and redshift space three-point correlation function,
       * following Slepian&Eisenstein, 2017
       *
       * @param KK the quantity \f$ k_l (r_1, r_2) \f$
       * 
       * @param rr vector of scales
       *
       * @param ll the order \f$l\f$
       *
       * @param kk vector of wavevector modules
       *
       * @param Pk dark matter power spectrum
       *
       * @return none
       */
      void k_ell (std::vector<std::vector<double>> &KK, const std::vector<double> rr, const int ll, const std::vector<double> kk, const std::vector<double> Pk);

      /**
       * @brief the multiplicative factor for \f$ \zeta_0 \f$, with local bias
       *
       * This function computes the multiplicative factor for \f$ \zeta_0 \f$, with  local bias:
       *
       * \f[
       *  l = 0 : b_1^3 \left( \frac{34}{21} \left[ 1+\frac{4}{3}\beta+\frac{1154}{1275}\beta^2+
       *  \frac{936}{2975}\beta^3+\frac{21}{425}\beta^4\right]+
       *  \gamma\left[ 1+\frac{2}{3}\beta+\frac{1}{9}\beta^2 \right ]  \right)
       * \f]
       *
       * with \f$b_1\f$ the linear bias, \f$\gamma\f$ the ratio of quadratic and linear bias, 
       * \f$\gamma= 2 b_2 / b_1 \f$ and \f$ \beta = f/b_1 \f$ with \f$ f \f$ the linear growth rate
       *
       * @param b1 the linear bias
       *
       * @param gamma the ratio of quadratic and linear bias, \f$\gamma= 2 b_2 / b_1 \f$, 
       *
       * @param beta the kaiser factor \f$ \beta = f/b_1 \f$ with \f$ f \f$ the linear growth rate
       *
       * @return  the multiplicative factor for \f$ \zeta_0 \f$, with local bias
       */
      double zeta_ell_0_factor (const double b1, const double gamma, const double beta);

      
      /**
       * @brief the multiplicative factor for \f$ \zeta_1 \f$, with local bias
       *
       * This function computes the multiplicative factor for \f$ \zeta_1 \f$, with local bias:
       *
       * \f[
       *   l = 1 : -b_1^3 \left[ 1+\frac{4}{3}\beta+\frac{82}{75}\beta^2+\frac{12}{25}\beta^3+\frac{3}{35}\beta^5 \right]
       * \f]
       *
       * with \f$b_1\f$ the linear bias and \f$ \beta = f/b_1 \f$ with \f$ f \f$ the linear growth rate
       *
       * @param b1 the linear bias
       *
       * @param beta the kaiser factor \f$ \beta = f/b_1 \f$ with \f$ f \f$ the linear growth rate
       *
       * @return  the multiplicative factor for \f$ \zeta_1 \f$, with local bias
       */
      double zeta_ell_1_factor (const double b1, const double beta);
      
      /**
       * @brief the multiplicative factor for \f$ \zeta_2 \f$, with local bias
       *
       * This function computes the multiplicative factor for \f$ \zeta_2 \f$, with local bias:
       *
       * \f[
       *  l = 2 : b_1^3 \left( \frac{8}{21} \left[ 1+\frac{4}{3}\beta+\frac{52}{21}\beta^2+
       *  \frac{81}{49}\beta^3+\frac{12}{35}\beta^4\right]+\frac{32}{945}\gamma \beta^2  \right)
       * \f]
       *
       * with \f$b_1\f$ the linear bias, \f$\gamma\f$ the ratio of quadratic and linear bias, 
       * \f$\gamma= 2 b_2 / b_1 \f$ and \f$ \beta = f/b_1 \f$ with \f$ f \f$ the linear growth rate
       *
       * @param b1 the linear bias
       *
       * @param gamma the ratio of quadratic and linear bias, \f$\gamma= 2 b_2 / b_1 \f$, 
       *
       * @param beta the kaiser factor \f$ \beta = f/b_1 \f$ with \f$ f \f$ the linear growth rate
       *
       * @return  the multiplicative factor for \f$ \zeta_2 \f$, with local bias
       */
      double zeta_ell_2_factor (const double b1, const double gamma, const double beta);

      /**
       * @brief the multiplicative factor for \f$ \zeta_3 \f$, with local bias
       *
       * This function computes the multiplicative factor for \f$ \zeta_3 \f$, with local bias:
       *
       * \f[
       *    l=3: -b_1^3 \left[ \frac{8}{75}\beta^2+\frac{16}{175}\beta^3+\frac{8}{315}\beta^4 \right]
       * \f]
       *
       * with \f$b_1\f$ the linear bias and \f$ \beta = f/b_1 \f$ with \f$ f \f$ the linear growth rate
       *
       * @param b1 the linear bias
       *
       * @param beta the kaiser factor \f$ \beta = f/b_1 \f$ with \f$ f \f$ the linear growth rate
       *
       * @return  the multiplicative factor for \f$ \zeta_3 \f$, with local bias
       */
      double zeta_ell_3_factor (const double b1, const double beta);

      /**
       * @brief the multiplicative factor for \f$ \zeta_4 \f$, with local bias
       *
       * This function computes the multiplicative factor for \f$ \zeta_4 \f$, with local bias:
       *
       * \f[
       *   l=4: -b_1^3 \left[ -\frac{32}{3675}\beta^2+\frac{32}{8575}\beta^3+\frac{128}{11025}\beta^4 \right]
       * \f]
       *
       * with \f$ b_1\f$ the linear bias and \f$ \beta = f/b_1 \f$ with \f$ f \f$ the linear growth rate
       *
       * @param b1 the linear bias
       *
       * @param beta the kaiser factor \f$ \beta = f/b_1 \f$ with \f$ f \f$ the linear growth rate
       *
       * @return  the multiplicative factor for \f$ \zeta_4 \f$, with local bias
       */
      double zeta_ell_4_factor (const double b1, const double beta);

      /**
       * @brief the multiplicative factor for \f$ \zeta_l, l>4 \f$, with local bias
       *
       * This function computes the multiplicative factor for \f$ \zeta_k, l>4 \f$, with  local bias:
       *
       * \f[
       *   l>4: b_1^3(7\beta^2+3\beta^3)
       * \f]
       *
       * with \f$b_1\f$ the linear bias and \f$ \beta = f/b_1 \f$ with \f$ f \f$ the linear growth rate
       *
       * @param b1 the linear bias
       *
       * @param beta the kaiser factor \f$ \beta = f/b_1 \f$ with \f$ f \f$ the linear growth rate
       *
       * @return  the multiplicative factor for \f$ \zeta_l, l>4 \f$, with local bias
       */
      double zeta_ell_k_factor (const double b1, const double beta);

      /**
       * @brief the multiplicative factor for \f$ \zeta_l, l=0 \f$, with non-local bias
       *
       * This function computes the multiplicative factor for \f$ \zeta_0 \f$, with non-local bias:
       *
       * \f[
       *   l=0: \frac{16 \beta^2 \gamma_t} {675} 
       * \f]
       *
       * with \f$\gamma_t\ = b_t/b_1\f$ the ratio between the tidal and the linear bias 
       * and \f$ \beta = f/b_1 \f$ with \f$ f \f$ the linear growth rate
       *
       * @param gamma_t the ratio between the tidal and the linear bias
       *
       * @param beta the kaiser factor \f$ \beta = f/b_1 \f$ with \f$ f \f$ the linear growth rate
       *
       * @return  the multiplicative factor for \f$ \zeta_l, l=0 \f$, with non-local bias
       */
      double zeta_ell_0_factor_tidal (const double gamma_t, const double beta);

      /**
       * @brief the multiplicative factor for \f$ \zeta_l, l=2 \f$, with non-local bias
       *
       * This function computes the multiplicative factor for \f$ \zeta_2 \f$, with non-local bias:
       *
       * \f[
       *   l=2: \frac{5}{2} \left(\frac{8}{15}+\frac{16\beta}{45}+\frac{344\beta^2}{4725} \right)\gamma_t 
       * \f]
       *
       * with \f$\gamma_t\ = b_t/b_1\f$ the ratio between the tidal and the linear bias 
       * and \f$ \beta = f/b_1\f$ with \f$ f \f$ the linear growth rate
       *
       * @param gamma_t the ratio between the tidal and the linear bias
       *
       * @param beta the kaiser factor \f$ \beta = f/b_1 \f$ with \f$ f \f$ the linear growth rate
       *
       * @return  the multiplicative factor for \f$ \zeta_l, l=2 \f$, with non-local bias
       */
      double zeta_ell_2_factor_tidal (const double gamma_t, const double beta);

      /**
       * @brief the multiplicative factor for \f$ \zeta_l, l=4 \f$, with non-local bias
       *
       * This function computes the multiplicative factor for \f$ \zeta_4 \f$, with non-local bias:
       *
       * \f[
       *   l=4: \frac{32 \beta^2 \gamma_t} {525} 
       * \f]
       *
       * with \f$\gamma_t\ = b_t/b_1\f$ the ratio between the tidal and the linear bias 
       * and \f$ \beta = f/b_1 \f$ with \f$ f \f$ the linear growth rate
       *
       * @param gamma_t the ratio between the tidal and the linear bias
       *
       * @param beta the kaiser factor \f$ \beta = f/b_1\f$ with \f$ f \f$ the linear growth rate
       *
       * @return the multiplicative factor for \f$ \zeta_l, l=4 \f$, with non-local bias
       */
      double zeta_ell_4_factor_tidal (const double gamma_t, const double beta);

      /**
       * @brief the pre-cyclic \f$ \zeta_l \f$
       *
       * This function computes the legendre multipoles coefficients
       * of the three-point correlation function using the tree-level expansion order of the biased 
       * and redshift space density field, as derived by Slepian&Eisenstein (2017) in configuration
       * space, based on the work presented by Scoccimarro et al (1999) for the bispectrum.
       *
       * The model is computed for a specific triangle configurations and order of the expansion
       * given the biasing parameters and \f$\beta\f$.
       *
       * Integrals and coefficients are computed by:
       * cbl::cosmology::Cosmology::xi_r_n,
       * cbl::cosmology::Cosmology::xi_r_n_pm,
       * cbl::cosmology::Cosmology::k_ell,
       * cbl::cosmology::Cosmology::zeta_ell_0_factor,
       * cbl::cosmology::Cosmology::zeta_ell_1_factor,
       * cbl::cosmology::Cosmology::zeta_ell_2_factor,
       * cbl::cosmology::Cosmology::zeta_ell_3_factor,
       * cbl::cosmology::Cosmology::zeta_ell_4_factor,
       * cbl::cosmology::Cosmology::zeta_ell_k_factor,
       * cbl::cosmology::Cosmology::zeta_ell_0_factor_tidal, 
       * cbl::cosmology::Cosmology::zeta_ell_2_factor_tidal, 
       * cbl::cosmology::Cosmology::zeta_ell_4_factor_tidal
       *
       * @param r1 the first triangle side
       *
       * @param r2 the second triangle side
       *
       * @param ell the order of the expansion
       *
       * @param b1 the linear bias
       *
       * @param b2 the quadratic bias
       *
       * @param bt the tidal bias
       *
       * @param beta the kaiser factor
       *
       * @param interp_xi_ell vector of interpolating function for terms \f$\xi^{[n]}, \xi^{[n\pm]}\f$
       * 
       * @param use_k if true, use the \f$k_l\f$ part of the model \f$O(\beta^2)\f$
       *
       * @param interp_k_ell interpolationg function for \f$ k_l\f$.
       *
       * @return the pre-cyclic \f$ \zeta_l \f$
       */
      double zeta_ell_precyclic (const double r1, const double r2, const int ell, const double b1, const double b2, const double bt, const double beta, std::vector<glob::FuncGrid> interp_xi_ell, const bool use_k, glob::FuncGrid2D interp_k_ell);

      /**
       * @brief the \f$ \zeta (r_1, r_2, \theta) \f$
       *
       * This function computes tree-level prediction for three-point correlation function
       * of haloes in redshift space, as derived by Slepian&Eisenstein (2017) in configuration
       * space, based on the work presented by Scoccimarro et al (1999) for the bispectrum.
       *
       * @param r1 the first triangle side
       *
       * @param r2 the second triangle side
       *
       * @param ntheta the number of \f$\theta\f$ bins
       *
       * @param b1 the linear bias
       *
       * @param b2 the quadratic bias
       *
       * @param bt the tidal bias
       *
       * @param beta the kaiser factor
       *
       * @param rr vector of scales
       *
       * @param kk vector of wavevector modules
       *
       * @param Pk dark matter power spectrum
       *
       * @param include_limits include the \f$\theta\f$ limits
       *
       * @param max_ll maximum order in the model
       *
       * @param use_k if true, use the \f$k_l\f$ part of the model \f$O(\beta^2)\f$
       *
       * @return the halo redshift space three-point correlation function
       */
      std::vector<double> zeta_RSD (const double r1, const double r2, const int ntheta, const double b1, const double b2, const double bt, const double beta, const std::vector<double> rr, const std::vector<double> kk, const std::vector<double> Pk, const bool include_limits=false, const int max_ll=4, const bool use_k=false);

      /**
       * @brief the \f$ \zeta (r_1, r_2, \theta) \f$
       *
       * This function computes tree-level prediction for three-point correlation function
       * of haloes in redshift space, as derived by Slepian&Eisenstein (2017) in configuration
       * space, based on the work presented by Scoccimarro et al (1999) for the bispectrum.
       *
       * @param r1 the first triangle side
       *
       * @param r2 the second triangle side
       *
       * @param ntheta the number of \f$\theta\f$ bins
       *
       * @param b1 the linear bias
       *
       * @param b2 the quadratic bias
       *
       * @param bt the tidal bias
       *
       * @param redshift the redshif
       *
       * @param method_Pk method used to compute the power spectrum;
       * valid choices for method_Pk are: CAMB [http://camb.info/],
       * classgal_v1 [http://class-code.net/], MPTbreeze-v1
       * [http://arxiv.org/abs/1207.1465], EisensteinHu
       * [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       *
       * @param step_r the number of bins for r
       *
       * @param step_k the number of bins for k
       *
       * @param force_RealSpace if true, force the computation
       * to be in real space
       *
       * @param include_limits include the \f$\theta\f$ limits
       *
       * @param max_ll maximum order in the model
       *
       * @param use_k if true, use the \f$k_l\f$ part of the model \f$O(\beta^2)\f$
       *
       * @return the halo redshift space three-point correlation function
       */
      std::vector<double> zeta_RSD (const double r1, const double r2, const int ntheta, const double b1, const double b2, const double bt, const double redshift, const std::string method_Pk, const int step_r, const int step_k, const bool force_RealSpace = false, const bool include_limits=false, const int max_ll=4, const bool use_k=false);

      ///@}
       
    };
    
  }


  // =====================================================================================
  
  
  namespace glob {
    
    /// @cond glob
    double GSL_bias_kernel_wrapper (double, void *);
    
    double func_xi_EH_GSL (double, void *);
    double func_sigma2M_EH_GSL (double, void *);
    
    double bias_kernel2 (double, void *); // non-Gaussian cosmologies
    double skewness_kernel (double *, size_t, void *); // non-Gaussian cosmologies
    /// @endcond 
    
    struct GSL_f_pars {
      double kt;
      double mass;
      std::string method_Pk;
      std::string output_root;
      int norm;
      double k_min;
      double k_max;
      bool GSL;
      double prec;
      std::string file_par;
      cosmology::Cosmology *pt_Cosmology;
    };

    struct STR_xi_EH
    {
      double Omega_matter;
      double Omega_baryon;
      double Omega_neutrinos;
      double massless_neutrinos;
      double massive_neutrinos;
      double Omega_DE;
      double Omega_radiation;
      double hh;
      double scalar_amp;
      double scalar_pivot;
      double n_spec;
      double w0;
      double wa;
      double fNL;
      int type_NG;
      double tau;
      std::string model;
      bool unit;
      std::string method_Pk;
      double rr;
      double redshift;
      double aa;
    };

    struct STR_sigma2M_EH
    {
      double Omega_matter;
      double Omega_baryon;
      double Omega_neutrinos;
      double massless_neutrinos;
      double massive_neutrinos;
      double Omega_DE;
      double Omega_radiation;
      double hh;
      double scalar_amp;
      double scalar_pivot;
      double n_spec;
      double w0;
      double wa;
      double fNL;
      int type_NG;
      double tau;
      std::string model;
      bool unit;
      std::string method_Pk;
      double redshift;
      double mass;
      double rho;
    };

    struct STR_MF
    {
      double redshift;
      std::string model_MF;
      std::string method_SS;
      std::string output_root;
      double Delta;
      std::string interpType;
      int norm; 
      double k_min; 
      double k_max;
      double prec;
      std::string input_file;
      bool is_parameter_file;
      bool default_delta;
      double delta_t;
    };

    struct STR_NG /* Cosimo Fedeli */
    {
      double Omega_matter;
      double Omega_baryon;
      double Omega_neutrinos;
      double massless_neutrinos;
      double massive_neutrinos;
      double Omega_DE;
      double Omega_radiation;
      double hh;
      double scalar_amp;
      double scalar_pivot;
      double n_spec;
      double w0;
      double wa;
      double fNL;
      int type_NG;
      double tau;
      std::string model;
      bool unit;
      double kt;
      double mass;
      std::string method_Pk;
      std::string output_root;
      int norm;
      double k_min;
      double k_max;
      bool GSL;
      double prec;
      std::string file_par;
    };
    
  }

}

#include "CosmClassFunc.h"

#endif
