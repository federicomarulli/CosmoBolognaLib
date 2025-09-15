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
 *  This file defines the interface of the abstract class Cosmology,
 *  which defines the members and methods of the derived classes used
 *  to for numerical calculations in different cosmological models
 *
 *  @author Federico Marulli 
 *
 *  @author federico.marulli3@unibo.it
 */

#ifndef __COSM__
#define __COSM__

#include <boost/math/special_functions/lambert_w.hpp>
#include "Likelihood.h"
#include "EisensteinHu.h"
#include "CAMB.h"


// ===================================================================================================


namespace cbl {

  /**
   *  @brief The namespace of the functions and classes used for <B>
   *  cosmological calculations </B>
   *  
   *  The \e cosmology namespace contains all the functions and
   *  classes used for cosmological calculations
   */
  namespace cosmology {

    /**
     *  @class Cosmology Cosmology.h "Headers/Cosmology.h"
     *
     *  @brief The class Cosmology
     *
     *  This abstract class is used to define the interface of the
     *  derived classes used for numerical calculations in different
     *  cosmological models
     */
    class Cosmology {

    protected:
      
      /// name used to label the cosmological model
      std::string m_model;
      
      /// \f$\Omega_m\f$: the total matter density (in units of the critical density) at z=0
      double m_Omega_matter;
      
      /// \f$\Omega_{\rm DE}\f$: the dark energy density (in units of the critical density) at z=0
      double m_Omega_DE;
      
      /// \f$\Omega_k\f$: the density of curvature energy (in units of the critical density) at z=0
      double m_Omega_k = 0;
      
      /// \f$\Omega_\nu\f$: the density of massive neutrinos (in units of the critical density) at z=0
      double m_Omega_neutrinos = 0;     

      /// \f$N_{eff}\f$: the effective number (for QED + non-instantaneous decoupling)
      double m_massless_neutrinos = 3.04;   

      /// the number of degenerate massive neutrino species
      int m_massive_neutrinos = 0.;
      
      /// \f$\Omega_{\rm CDM}\f$: the cold dark matter density (in units of the critical density) at z=0
      double m_Omega_CDM;
      
      /// \f$\Omega_b\f$: the baryon density (in units of the critical density) at z=0
      double m_Omega_baryon;           

      /// \f$\Omega_{\rm rad}\f$: the radiation density (in units of the critical density) at z=0 
      double m_Omega_radiation;      

      /// \f$h\f$: the dimensionless Hubble parameter, \f$H_0/100\f$
      double m_h;                    

      /// \f$A_s\f$: the initial scalar amplitude of the matter power spectrum
      double m_scalar_amp;           

      /// the scalar pivot k in \f$Mpc^{-1}\f$
      double m_scalar_pivot;
      
      /// \f$\sigma_8\f$: the matter power spectrum normalisation
      double m_sigma8 = -1;   
      
      /// \f$n_{\rm spec}\f$: the primordial spectral index
      double m_n_spec;      

      /// \f$\tau\f$: the Thomson scattering optical depth due to reionization
      double m_tau;

      /// false \f$\rightarrow\f$ quantities are provided in phyical units; true \f$\rightarrow\f$ quantities are provided in cosmological units (i.e. in units of \e h, e.g. Mpc/h)
      bool m_unit;
      
      /// \f$w_0\f$: the parameter of the dark energy equation of state (CPL parameterisation)
      double m_w0 = -1;
      
      /// \f$w_a\f$: the parameter of the dark energy equation of state (CPL parameterisation)
      double m_wa = 0;
      
      /// \f$f_{NL}\f$: the non-Gaussian amplitude
      double m_fNL = 0.;
      
      /// the non-Gaussian shape (type=1 local, type=2 equilateral, type=3 enfolded, type=4 orthogonal)
      int m_type_NG = 1;

      
      // -----------------------------------------------------------------------
   
      
    public:
      
      /**
       *  @name Constructors/destructors
       */
      ///@{

      /**
       *  @brief Default constuctor
       */
      Cosmology () = default;

      /**
       *  @brief Default destructor
       */
      virtual ~Cosmology () = default;

      ///@}
      

      /**
       *  @name Functions to get the private members of the class
       */
      ///@{

      /**
       * @brief Pointer to the current class instance
       *
       * @return a pointer to the current class instance
       */
      virtual std::shared_ptr<Cosmology> clone () const = 0;
      
      /**
       * @brief Vector containing the names of the cosmological model
       * parameters
       *
       * @return a vector containing the cosmological parameter names
       */
      virtual std::vector<std::string> cosmological_parameter_names () = 0;
      
      /**
       *  @brief Print the values of the cosmological parameters on
       *  the screen
       */
      virtual void print_parameters () const = 0;
      
      /**
       *  @brief The value of a chosen cosmological model parameter
       *  
       *  @param parameter the cosmological parameter
       *
       *  @return value of the cosmological parameter
       */
      virtual double value (const std::string parameter) const = 0;

      /**
       *  @brief The name of the cosmological model
       *
       *  @return the name used to label the cosmological model
       */
      std::string model () const { return m_model; }
      
      /**
       *  @brief The total matter density at z=0
       *
       *  @return \f$\Omega_m\f$: the total matter density (in units
       *  of the critical density) at z=0
       */
      inline double Omega_matter () const { return m_Omega_matter; };
      
      /**
       *  @brief The cold dark matter density at z=0
       *
       *  @return \f$Omega_{CDM}\f$: the cold dark matter density at
       *  z=0
       */
      inline double Omega_CDM () const { return m_Omega_CDM; };
      
      /**
       *  @brief The baryon density at z=0
       *
       *  @return \f$\Omega_b\f$: the baryon density (in units of the
       *  critical density) at z=0
       */
      inline double Omega_baryon () const { return m_Omega_baryon; };

      /**
       *  @brief The radiation density at z=0
       *
       *  @return \f$\Omega_{\rm rad}\f$: the radiation density (in
       *  units of the critical density) at z=0
       */
      inline double Omega_radiation () const { return m_Omega_radiation; };   

      /**
       *  @brief The dimensionless Hubble parameter
       *
       *  @return \f$h\f$: the dimensionless Hubble parameter,
       *  \f$H_0/100\f$
       */
      inline double little_h () const { return m_h; };

      /**
       *  @brief The matter power spectrum normalisation
       *
       *  @return \f$\sigma_8\f$: the matter power spectrum
       *  normalisation
       */
      inline double sigma8 () const { return m_sigma8; };

      /**
       *  @brief The initial scalar amplitude of the matter power
       *  spectrum
       *
       *  @return \f$A_s\f$: the initial scalar amplitude of the
       *  matter power spectrum
       */
      inline double scalar_amp () const { return m_scalar_amp; };

      /**
       *  @brief The scalar pivot 
       *
       *  @return the scalar pivot k in \f$Mpc^{-1}\f$
       */
      inline double scalar_pivot () const { return m_scalar_pivot; };

      /**
       *  @brief The primordial spectral index
       *
       *  @return \f$n_{\rm spec}\f$: the primordial spectral index
       */
      inline double n_spec () const { return m_n_spec; };

      /**
       *  @brief The Thomson scattering optical depth
       *
       *  @return \f$\tau\f$: the Thomson scattering optical depth due
       *  to reionization
       */
      inline double tau () const { return m_tau; };

      /**
       *  @brief The \e unit variable, used to chose the
       *  cosmological units
       *
       *  @return unit: false \f$\rightarrow\f$ quantities are
       *  provided in phyical units; true \f$\rightarrow\f$ quantities
       *  are provided in cosmological units (i.e. in units of \e h,
       *  e.g. Mpc/h)
       */
      inline bool unit () const { return m_unit; };

      /**
       *  @brief The density of curvature energy at z=0
       *
       *  @return \f$\Omega_k\f$: the density of curvature energy at
       *  z=0
       */ 
      virtual inline double Omega_k () const
      { return m_Omega_k; };

      /**
       *  @brief The density of massive neutrinos at z=0
       *
       *  @return \f$\Omega_\nu\f$: the density of massive neutrinos
       *  at z=0
       */
      virtual inline double Omega_neutrinos () const
      { return m_Omega_neutrinos; };
      
      /**
       *  @brief The effective number of neutrino species
       *
       *  @return \f$N_{\rm eff}\f$: the effective number of neutrino
       *  species (for QED + non-instantaneous decoupling)
       */
      virtual inline double massless_neutrinos () const
      { return m_massless_neutrinos; };

      /**
       *  @brief The number of degenerate massive neutrino species
       *
       *  @return the number of degenerate massive neutrino species
       */
      virtual inline int massive_neutrinos () const
      { return m_massive_neutrinos; };
     
      /**
       *  @brief Get \f$w_0\f$
       *
       *  @return \f$w_0\f$: one of the parameters of the dark
       *  energy equation of state (CPL parameterisation)
       */
      virtual inline double w0 () const
      { return m_w0; };

      /**
       *  @brief Get \f$w_a\f$
       *
       *  @return \f$w_a\f$: one of the parameters of the dark
       *  energy equation of state (CPL parameterisation)
       */
      virtual inline double wa () const 
      { return m_wa; };

      /**
       *  @brief The non-Gaussian amplitude
       *
       *  @return f<SUB>NL</SUB>: the non-Gaussian amplitude
       */
      virtual inline double fNL () const 
      { return m_fNL; };
      
      /**
       *  @brief The non-Gaussian shape
       *
       *  @return the non-Gaussian shape (type=1 local, type=2
       *  equilateral, type=3 enfolded, type=4 orthogonal)
       */
      virtual inline int type_NG () const 
      { return m_type_NG; };
	
      ///@}


      /**
       *  @name Functions to set the private members of the class
       */
      ///@{

      /**
       *  @brief Set the value of one cosmological paramter
       *  
       *  @param parameter cosmological parameter to set
       *  @param value the new value for the parameter 
       */
      virtual void set_parameter (const std::string parameter, const double value) = 0;

      /**
       *  @brief Set the value of some cosmological paramters
       *  
       *  @param parameter vector containing the cosmological
       *  parameters to set
       *
       *  @param value vector containing the new values for the
       *  parameters
       */
      void set_parameters (const std::vector<std::string> parameter, const std::vector<double> value);

      /**
       *  @brief Set the value of unit
       *
       *  @param unit false \f$\rightarrow\f$ phyical units; true
       *  \f$\rightarrow\f$ cosmological units (i.e. in units of \e h)
       */
      inline void set_unit (const bool unit=true) { m_unit = unit; }
    
      /**
       *  @brief Set the value of f<SUB>NL</SUB>
       *
       *  @param fNL f<SUB>NL</SUB>: the non-Gaussian amplitude
       */
      virtual void set_fNL (const double fNL)
      { (void)fNL; cbl::ErrorCBL("This cosmological model does not have this parameter", "set_fNL", "Cosmology.h"); };
      
      /**
       *  @brief Set the value of the non-Gaussian shape
       *
       *  @param type_NG the non-Gaussian shape (type=1 local, type=2
       *  equilateral, type=3 enfolded, type=4 orthogonal)
       */
      virtual void set_type_NG (const int type_NG)
      { (void)type_NG; cbl::ErrorCBL("This cosmological model does not have this parameter", "set_type_NG", "Cosmology.h"); };

      ///@}

    
      /**
       *  @name Functions to estimate general cosmological parameters
       */
      ///@{
      
      /**
       *  @brief The best-fit cosmological parameter values estimated
       *  by key past cosmological missions (e.g WMAP, Planck)
       *
       *  This member function provides the values of the main
       *  cosmological parameters from different cosmological probes
       *
       *  @param cosmological_model the name used to label the
       *  cosmological model parameters from a given probe
       *  (e.g. Planck18)
       *
       *  @return a dictionary containing the cosmological model
       *  parameters
       */
      std::map<std::string, double> parameter_value (std::string cosmological_model) const;

      /**
       *  @brief The matter density parameter
       *
       *  This member function computes the matter density parameter
       *  at a given redshift (see e.g. Percival 2005, Eq. 6):
       *
       *  \f[ \Omega_m(z) = \Omega_m\frac{(1+z)^3}{E^2(z)} \f]
       *
       *  where \f$\Omega_m\equiv\Omega_m(z=0)\f$ and \f$E(z)\f$ is
       *  computed by cosmology::Cosmology::EE
       *
       *  @param redshift the redshift
       *
       *  @return \f$\Omega_m(z)\f$: the matter density parameter
       */
      double Omega_matter (const double redshift) const; 
      
      /**
       *  @brief The radiation density parameter
       *
       *  This member function computes the radiation density
       *  parameter at a given redshift:
       *
       *  \f[ \Omega_{\rm rad}(z) = \Omega_{\rm rad}\frac{(1+z)^4}{E^2(z)} \f]
       *
       *  where \f$\Omega_{\rm rad}\equiv\Omega_{\rm rad}(z=0)\f$ and
       *  \f$E(z)\f$ is computed by cosmology::Cosmology::EE
       *
       *  @param redshift the redshift
       *
       *  @return \f$\Omega_{\rm rad}(z)\f$: the radiation density
       *  parameter
       */
      double Omega_radiation (const double redshift) const;

      /**
       *  @brief The dimensionless dark energy density
       *
       *  This member function computes the dimensionless dark energy
       *  density, which is used e.g. to estimate the Hubble function
       *  (see e.g. Eq.(1.2) of Bassett & Hlozek 2010):
       *
       *  \f[ f_{\rm DE}(z) = \frac{\rho_{\rm DE}(z)}{\rho_{\rm
       *  DE}(z=0)} \f]
       *
       *  @param redshift the redshift
       *
       *  @return the dimensionless dark energy density
       */
      virtual double f_DE (const double redshift) const
      { (void)redshift; cbl::ErrorCBL("This cosmological model does not have this parameter", "f_DE", "Cosmology.h"); return 0.; }; 

      /**
       *  @brief The dark energy density parameter
       *
       *  This member function computes the dark energy density
       *  parameter at a given redshift:
       *
       *  \f[ \Omega_{\rm DE}(z) = \Omega_{\rm DE}\frac{f_{\rm
       *  DE}(z)}{E^2(z)} \f]
       *
       *  where \f$\Omega_{\rm DE}\equiv\Omega_{\rm DE}(z=0)\f$,
       *  \f$f_{\rm DE}\f$ is computed by cosmology::Cosmology::f_DE
       *  and \f$E(z)\f$ is computed by cosmology::Cosmology::EE
       *
       *  @param redshift the redshift
       *
       *  @return \f$\Omega_{\rm DE}(z)\f$: the dark energy density
       *  parameter
       */
      virtual double Omega_DE (const double redshift=0) const = 0;
      
      /**
       *  @brief The neutrino density at a given redshift
       *
       *  @param redshift the redshift
       *
       *  @return &Omega;<SUB>&nu;</SUB>
       */
      virtual double Omega_Nu (const double redshift) const 
      { (void)redshift; cbl::ErrorCBL("This cosmological model does not have this parameter", "Omega_Nu", "Cosmology.h"); return 0.; }; 
    
      
      /**
       *  @brief The dimensionless Hubble function
       *
       *  This member function computes the dimensionless Hubble
       *  function:
       *
       *  \f[ E(z) = \frac{H(z)}{H_0} \f]
       *
       *  @param redshift the redshift
       *
       *  @return the dimensionless Hubble function
       */
      virtual double EE (const double redshift) const = 0;
      
      /**
       *  @brief The linear growth rate
       *
       *  This member function computes the linear growth rate at a
       *  given redshift:
       *
       *  \f[ f(z) = \frac{d\,\ln D(z)}{d\,\ln a(z)} \f]
       *
       *  where \f$ D(z) \f$ is the linear growth factor.
       *
       *  @param redshift the redshift at which the linear growth rate
       *  is evaluated
       *
       *  @param prec precision used for the resolution of the
       *  differential equation, used in the cosmological models in
       *  which an analytic approximation is not accurate enough
       *
       *  @return the linear growth rate
       */
      virtual double linear_growth_rate (const double redshift, const double prec=1.e-4) const = 0;
      
      /**
       *  @brief The Hubble function
       *
       *  This member function computes the Hubble function,
       *  \f$H(z)\f$:
       *
       *  @param redshift the redshift
       *
       *  @return the Hubble function
       */
      virtual double Hubble (const double redshift) const;
      
      /**
       *  @brief The normalised amplitude of the growing mode at a
       *  given redshift
       *
       *  This member function computes the following quantity:
       *
       *  \f[ \frac{D(z)}{D(z_{norm})} = \exp \left[
       *  \int_{a_{norm}(z_{norm})}^{a(z)} \frac{{\rm d} a'}{a'} f(a')
       *  \right] \f]
       *
       *  where \f$a=1/(1+z)\f$, \f$a_{norm}=1/(1+z_{norm})\f$, and
       *  \f$f\equiv\frac{{\rm d}\ln D(a)}{{\rm d}\ln a}\f$ is the
       *  linear growth rate computed by
       *  cbl::cosmology::Cosmology::linear_growth_rate
       *
       *  @param redshift the redshift 
       *
       *  @param redshift_norm the redshift at with the amplitude of
       *  the growing mode is normalised
       *
       *  @param prec precision used for the resolution of the
       *  differential equation when the the normalised growth factor
       *  is computed
       *
       *  @return the amplitude of the growing mode
       */
      virtual double DN (const double redshift, const double redshift_norm=0., const double prec=1.e-4) const;

      /**
       *  @brief The amplitude of the growing mode at a given
       *  redshift
       *
       *  This member function computes the following quantity (see e.g. Eq.(15)
       *  by Percival (2005), A&A 443, 819):
       *
       *  \f[ D(z) = \frac{5 \Omega_m }{2} E(a) \int_0^a
       *  \frac{{\rm d}\,a'}{[a' E(a')]^3} \f]
       *
       *  where \f$a=1/(1+z)\f$
       *
       *  @param redshift the redshift 
       *
       *  @return the amplitude of the growing mode
       *
       *  @warning the current implementation is valid only for
       *  &Lambda;CDM cosmologies (e.g. Percival (2005) A&A 443, 819,
       *  and references therein)
       */
      virtual double DD (const double redshift) const;
      
      /**
       *  @brief The linear growth factor at a given redshift
       *
       *  This member function computes the following quantity (see e.g. Eq.(1) by
       *  Hamilton 2001):
       *
       *  \f[ g(z) \equiv \tilde{D}(z)(1+z) \f]
       *
       *  where \f$a=1/(1+z)\f$ and \f$D(z)\f$ is computed by
       *  cbl::cosmology::Cosmology::DD
       *
       *  @param redshift the redshift
       *
       *  @return the linear growth factor
       */
      virtual double gg (const double redshift) const;
      
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
       *
       *  @param redshift the redshift
       *
       *  @return t<SUB>lookback</SUB> [Gyr]
       */
      double lookback_time (const double redshift) const; 

      /**
       *  @brief cosmic time at a given redshift
       *  @param redshift the redshift
       *  @return t<SUB>cosmic</SUB> [Gyr]
       */
      double cosmic_time (const double redshift) const;
    
      /**
       *  @brief the deceleration parameter at a given redshift
       *  @param redshift the redshift
       *  @return q
       */
      virtual double qq (const double redshift) const = 0;

      /**
       *  @brief derivative of the Hubble function at a given redshift
       *  @param redshift the redshift
       *  @return dH/dz
       */
      double Hdot (const double redshift) const;

      /**
       *  @brief redshift at which the Universe begins to accelerate 
       *
       *  see e.g. de Araujo 2005
       *
       *  @return z<SUB>acc</SUB>
       */
      virtual double z_acc () const = 0; 

      /**
       *  @brief redshift of matter-dark energy equality
       *
       *  see e.g. de Araujo 2005
       *
       *  @return z<SUB>eq</SUB>
       */
      virtual double z_eq () const = 0;

      /**
       *  @brief redshift of matter-radiation equality
       *
       *  @param T_CMB the temperature of the CMB
       *
       *  @return z<SUB>eq</SUB>
       */
      virtual double z_eq_rad (const double T_CMB=2.7255) const = 0;
    
      /**
       *  @brief maximum absolute magnitude to have a volume-limited
       *  catalogue
       *
       *  @param z_max maximum redshift
       *  @param mag_lim magnitude limit
       *  @return r<SUB>s</SUB>
       */
      virtual double Mag_Volume_limited (const double z_max=1., const double mag_lim=-20.) const = 0;

      /**
       *  @brief bolometric luminosity
       *  @param redshift the redshift
       *  @param flux flux
       *  @return L<SUB>bol</SUB>
       */
      virtual double Lum_bol (const double redshift=0., const double flux=1.) const = 0;  

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
      virtual double Redshift (const double d_c=1., const double z1_guess=0., const double z2_guess=10., const double prec=0.0001) const = 0;
      
      /**
       *  @brief redshift at a given cosmic time
       *  @param time cosmic time
       *  @param z1_guess minimum redshift used to search the redshift
       *  @param z2_guess maximum redshift used to search the redshift
       *  @return redshift
       */
      virtual double Redshift_time (const double time, const double z1_guess, const double z2_guess) const = 0;

      /**
       *  @brief spherical collapse density threshold at a given
       *  redshift
       *
       *  this function computes the spherical collapse density
       *  threshold, \f$\delta_c\f$, by using the approximated
       *  equation (C.28) provided by Nakamura & Suto (1997)
       *
       *  \f[\delta_c(z) \simeq \frac{3}{20}(12\pi)^{2/3}
       *  \{1.+0.012299\log[\Omega_m(z)]\}\f]
       *
       *  @param redshift the redshift
       *  @return \f$\delta_c\f$
       */
      double deltac (const double redshift) const;  

      /**
       *  @brief The density of curvature energy 
       *
       *  This member function computes the density of curvature
       *  energy at a given redshift:
       *
       *  \f[ \Omega_k(z) = \Omega_k\frac{(1+z)^2}{E^2(z)} \f]
       *
       *  where \f$\Omega_k\equiv\Omega_k(z=0)\f$ and
       *  \f$E(z)\f$ is computed by cosmology::Cosmology::EE
       *
       *  @param redshift the redshift
       *
       *  @return \f$\Omega_k(z)\f$: the radiation density
       *  parameter
       */
      virtual double Omega_k (const double redshift) const
      { (void)redshift; cbl::ErrorCBL("This cosmological model does not have this parameter", "Omega_k", "Cosmology.h"); return 0.; };

      /**
       *  @brief The total neutrino mass
       *
       *  this function computes the neutrino mass as follows:
       *
       *  \f[\sum m_\nu = \Omega_\nu\cdot94\, h^2 eV\f]
       *  
       *  @return \f$\sum m_\nu\f$
       */
      virtual double neutrino_mass () const
      { cbl::ErrorCBL("This cosmological model does not have this parameter", "neutrino_mass", "Cosmology.h"); return 0.; };

      /**
       *  @brief The Hubble constant at z=0
       *
       *  This function computes the Hubble constant at z=0 in units
       *  of km/sec/Mpc as:
       *
       *  \f[ H_0 = 100\,h \f]
       *
       *  @return \f$H_0\f$: the Hubble constant at z=0
       */
      inline double H0 () const { return (m_unit) ? 100. : 100.*m_h; };

      /**
       *  @brief The Hubble time
       *
       *  This function computes:
       *
       *   \f[ t_H = \frac{1}{H_0} \f]
       *
       *  @return \f$t_H\f$: the Hubble time
       */
      inline double t_H () const { return 1./H0(); };

      /**
       *  @brief The Hubble distance
       *
       *  This function computes: 
       *
       *  \f[ D_H = \frac{c}{H_0} \f]
       *
       *  @return \f$D_H\f$: the Hubble distance
       */
      inline double D_H () const { return par::cc/H0(); };
      
      /**
       *  @brief The critical cosmic density
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
       *  @brief The mean cosmic background density
       *
       *  \f[\rho_m(z) = \rho_{crit}(z)\Omega_m(z) =
       *  \frac{3H^2(z)}{8\pi G}\Omega_m(z)\f]
       *
       *  @param redshift the redshift
       *
       *  @param unit1 true \f$\rightarrow\f$ force cosmological units
       *
       *  @return &rho;<SUB>mean</SUB>: the mean cosmic background
       *  density [Msun*Mpc^-3(*h^2)]
       */
      virtual double rho_m (const double redshift=0., const bool unit1=false) const = 0;  
 
      /**
       *  @brief The critical overdensity
       *
       *  this function computes the critical overdensity,
       *  \f$\Delta_c(z)\equiv\Delta^{vir}_c(z)\f$, using approximated
       *  equations valid only for a flat Universe (see e.g. Coe
       *  2010), where \f$\rho_{vir} = \Delta^{vir}_c\rho_c =
       *  \Delta^{vir}_b\rho_m = \Delta^{vir}_b\Omega_m\rho_c\f$:
       *  
       *  - Bryan & Norman (1998)
       *
       *  \f[\Delta_c\simeq18\pi^2+60x-32x^2\, \mbox{for}\,
       *  \Omega_\Lambda=0\f]
       *
       *  \f[\Delta_c\simeq18\pi^2+82x-39x^2\, \mbox{for}\,
       *  \Omega_k=0\f]
       *  
       *  \f[x=\Omega_m(z)-1\f]
       *
       *  - Eke et al. (1998)
       * 
       *  \f[\Delta_c\simeq178\Omega_m(z)^{0.3}\, \mbox{for}\,
       *  \Omega_\Lambda=0\f]
       *
       *  \f[\Delta_c\simeq178\Omega_m(z)^{0.45}\, \mbox{for}\,
       *  \Omega_k=0\f]
       *
       *  - Nakamura & Suto (1998)
       *  
       *  \f[\Delta_c\simeq18\pi^2(1+0.4093x^{2.7152})\Omega_m(z)\f]
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
      virtual double Delta_c (const double redshift, const std::string author="BryanNorman") const = 0;
      
      /**
       *  @brief The dark energy equation of state in the CPL
       *  parameterisation
       *
       *  This member function computes the dark energy equation of
       *  state in the Chevallier-Polarski-Linder (CPL)
       *  parameterisation, as a function of redshift:
       *
       *  \f[ w(z) = \frac{p_{\rm DE}(z)}{\rho_{\rm DE}(z)} = w_0+w_a
       *  \frac{z}{1+z} \f]
       * 
       *  @param redshift the redshift
       *
       *  @return \f$w(z)\f$: the dark energy equation of state
       */
      virtual double w_CPL (const double redshift) const
      { (void)redshift; cbl::ErrorCBL("This cosmological model does not have these parameters", "w_CPL", "Cosmology.h"); return 0.; };
    
      ///@}
    

      /**
       *  @name Functions to estimate cosmological distances and
       *  volumes
       */
      ///@{

      /**
       *  @brief The comoving line-of-sight distance at a given redshift
       *
       *  for demonstration, see \ref distances.cpp
       *
       *  @param redshift the redshift
       *  @return D<SUB>C</SUB>
       */
      virtual double D_C (const double redshift) const = 0;  

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
       */
      void D_C_table (const std::string file_table, const double z_min, const double z_max, const int step, std::vector<double> &Redshift, std::vector<double> &dc) const;

      /**
       *  @brief The comoving transverse distance at a given redshift
       *  @param redshift the redshift
       *  @return D<SUB>M</SUB>
       */
      virtual double D_M (const double redshift) const = 0;

      /**
       *  @brief The angular diameter distance at a given redshift
       *
       *  @param redshift the redshift
       *
       *  @return \f$D_A\f$
       */
      double D_A (const double redshift) const; 

      /**
       *  @brief The angular diameter distance between objects at two
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
      virtual double D_A (const double z1, const double z2) const = 0;
      
      /**
       *  @brief The luminosity distance at a given redshift
       *  @param redshift the redshift
       *  @return D<SUB>L</SUB>
       */
      double D_L (const double redshift) const; 
  
      /**
       *  @brief The average distance at a given redshift, used to
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
       *  @brief The distance at a given redshift 
       *
       *  The distances available are:
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
       *  @brief comoving volume for a given redshift range and sky
       *  area
       *
       *  this function provides an approximated expression of the
       *  comoving volume valid only for a &Lambda;CDM model
       *
       *  @param z1 minimum redshift
       *
       *  @param z2 maximum redshift
       *
       *  @param Area sky area [square degrees]
       *
       *  @return comoving volume [m_unit=true \f$\rightarrow
       *  (\mbox{Mpc}/h)^3\f$; m_unit=false \f$\rightarrow
       *  \mbox{Mpc}^3\f$]
       */
      virtual double Volume (const double z1, const double z2, const double Area) const = 0; 

      /**
       *  @brief comoving volume for a given redshift range and
       *  R.A.-Dec limits
       *
       *  this function provides an approximated expression of the
       *  comoving volume valid only for a &Lambda;CDM model
       *
       *  @param z1 minimum redshift
       *
       *  @param z2 maximum redshift
       *
       *  @param RA_min minimum Right ascension [radians]
       *
       *  @param RA_max maximum Right ascension [radians]
       *
       *  @param Dec_min minimum Declination [radians]
       *
       *  @param Dec_max minimum Declination [radians]
       *
       *  @return comoving volume [m_unit=true \f$\rightarrow
       *  (\mbox{Mpc}/h)^3\f$; m_unit=false \f$\rightarrow
       *  \mbox{Mpc}^3\f$]
       */
      virtual double Volume (const double z1, const double z2, const double RA_min, const double RA_max, const double Dec_min, const double Dec_max) const = 0; 
      
      /**
       *  @brief total comoving volume from z=0 to z
       *
       *  from Hogg 2000, Eq. 29
       *
       *  @author Alfonso Veropalumbo
       *  @author alfonso.veropalumbo@unibo.it
       *
       *  @param zz redshift
       *
       *  @return comoving volume [m_unit=true \f$\rightarrow
       *  (\mbox{Mpc}/h)^3\f$; m_unit=false \f$\rightarrow
       *  \mbox{Mpc}^3\f$]
       */
      virtual double Volume (const double zz) const = 0;
  
      /**
       *  @brief maximum redshift for a given volume, sky area and
       *  minimum redshift
       *  @param Volume volume
       *  @param Area sky area
       *  @param z_min minimum redshift
       *  @return redshift
       */
      virtual double max_redshift (const double Volume, const double Area, const double z_min) const = 0; 

      /**
       *  @brief The derivative of the comoving volume,
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
      
      /**
       *  @brief The sound horizon at the drag epoch predicted by
       *  Eisenstein & Hu 1998 or by CAMB [http://camb.info/],
       *  analytical formula by Aubourg et al. 2014
       *
       *  see Eisenstein & Hu 1998, Section 2.1
       *  see Anderson et al 2014, Eq. 16
       * 
       *  @author Alfonso Veropalumbo 
       *  @author alfonso.veropalumbo@unibo.it 
       *  @param T_CMB CMB temperature
       *  @param method CAMB or EisensteinHu method 
       *  @return r<SUB>s</SUB>
       */
      virtual double sound_horizon_at_drag_epoch (std::string method="CAMB", const double T_CMB=par::TCMB) const = 0;

      /**
       *  @brief Amplitude of the curvature perturbations
       *  
       *  This function provides an approximate value of
       *  A<SUB>s</SUB>, for a fiven value of &sigma;<SUB>8</SUB>,
       *  from Eq.(3) by Vikhlinin et al. 2009, ApJ, 692, 1060; it is
       *  valid only without massive neutrinos!  (see also Hu & Jain
       *  2004)
       *
       *  @param sigma8 &sigma;<SUB>8</SUB> the power spectrum normalisation
       *  @return A<SUB>s</SUB>
       */
      virtual double As (const double sigma8) const = 0;

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
      virtual double sigma8_interpolated (const double redshift) const = 0;

      ///@}

    };   
  }
}

#endif
