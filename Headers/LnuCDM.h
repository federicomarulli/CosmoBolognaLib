/********************************************************************
 *  Copyright (C) 2023 by Federico Marulli                          *
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
 *  @file Headers/LnuCDM.h
 *
 *  @brief The class LnuCDM 
 *
 *  This file defines the interface of the class LnuCDM, used for
 *  cosmological calculations in the flat Lambda-cold dark matter
 *  model with massive neutrinos
 *
 *  @author Federico Marulli 
 *
 *  @author federico.marulli3@unibo.it
 */

#ifndef __LnuCDM__
#define __LnuCDM__


#include "LCDM.h"


// ===================================================================================================


namespace cbl {

  namespace cosmology {
    
    /**
     *  @class LnuCDM LnuCDM.h "Headers/LnuCDM.h"
     *
     *  @brief The class LnuCDM
     *
     *  This class is used to handle objects of type <EM> LnuCDM
     *  </EM>. They are used for cosmological calculations in the flat
     *  Lambda-cold dark matter model with massive neutrinos
     */
    class LnuCDM : public virtual LCDM {
  
    public:
      
      /**
       *  @name Constructors/destructors
       */
      ///@{

      /**
       *  @brief Default constructor of the class LnuCDM
       */
      LnuCDM () { m_model = "LnuCDM"; }

      /**
       *  @brief Main Constructor of the class LnuCDM
       *
       *  @param Omega_matter \f$\Omega_m\f$: the density of baryons,
       *  cold dark matter and massive neutrinos (in units of the
       *  critical density) at z=0
       *
       *  @param Omega_baryon \f$\Omega_b\f$: the density of baryons
       *  at z=0
       *
       *  @param Omega_radiation \f$\Omega_{\rm rad}\f$: the density of
       *  radiation at z=0
       *
       *  @param little_h \e h: the Hubble parameter, \f$H_0/100\f$
       *
       *  @param scalar_amp \f$A_s\f$: the initial scalar amplitude of
       *  the power spectrum
       *
       *  @param scalar_pivot the scalar pivot k in \f$Mpc^{-1}\f$
       *
       *  @param n_spec \f$n_{spec}\f$: the primordial spectral index
       *
       *  @param tau \f$\tau\f$: Thomson scattering optical depth due
       *  to reionization
       *
       *  @param Omega_neutrinos \f$\Omega_\nu\f$: the density of
       *  massive neutrinos at z=0
       *
       *  @param massless_neutrinos \f$N_{eff}\f$: the effective
       *  number (for QED + non-instantaneous decoupling)
       *
       *  @param massive_neutrinos the number of degenerate massive
       *  neutrino species
       *
       *  @param unit false \f$\rightarrow\f$ quantities are provided
       *  in phyical units; true \f$\rightarrow\f$ quantities are
       *  provided in cosmological units (i.e. in units of \e h,
       *  e.g. Mpc/h)
       */
      LnuCDM (const double Omega_matter, const double Omega_baryon, const double Omega_radiation, const double little_h, const double scalar_amp, const double scalar_pivot, const double n_spec, const double tau, const double Omega_neutrinos, const double massless_neutrinos, const int massive_neutrinos, const bool unit=true);

      /**
       *  @brief Constructor of the class LnuCDM that uses built-in
       *  cosmological model parameter sets
       *
       *  @param cosmological_model the built-in cosmological model
       *
       *  @param unit false \f$\rightarrow\f$ quantities are provided
       *  in phyical units; true \f$\rightarrow\f$ quantities are
       *  provided in cosmological units (i.e. in units of \e h,
       *  e.g. Mpc/h)
       */
      LnuCDM (const std::string cosmological_model, const bool unit=true);

      /**
       *  @brief Default destructor
       */
      virtual ~LnuCDM () = default;

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
      std::shared_ptr<Cosmology> clone() const override 
      { return std::make_shared<LnuCDM>(*this); }

      /**
       *  @brief The density of massive neutrinos at z=0
       *
       *  @return \f$\Omega_\nu\f$: the density of massive neutrinos
       *  at z=0
       */
      double Omega_neutrinos () const { return m_Omega_neutrinos; };
      
      /**
       *  @brief The density of massive neutrinos at a given redshift
       *
       *  @param redshift the redshift
       *
       *  @return \f$\Omega_\nu\f$: the density of massive neutrinos
       *  at a given redshift
       */
      double Omega_neutrinos (const double redshift) const { return m_Omega_neutrinos/pow(EE(redshift),2)*(1.+redshift)*(1.+redshift)*(1.+redshift); };

      /**
       * @brief Vector containing the names of the cosmological model
       * parameters
       *
       * @return a vector containing the cosmological parameter names
       */
      inline virtual std::vector<std::string> cosmological_parameter_names  () override
      { return {"Omega_matter", "Omega_CDM", "Omega_baryon", "Omega_radiation", "little_h", "scalar_amp", "scalar_pivot", "sigma8", "n_spec", "tau", "Omega_neutrinos", "massless_neutrinos", "massive_neutrinos"}; }
      
       /**
       *  @brief Print the values of the cosmological parameters on
       *  the screen
       */
      virtual void print_parameters () const override;
      
      /**
       *  @brief The value of a chosen cosmological model parameter
       *  
       *  @param parameter the cosmological parameter
       *
       *  @return value of the cosmological parameter
       */
      virtual double value (const std::string parameter) const override;

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
      void set_parameter (const std::string parameter, const double value) override;

      ///@}

    
      /**
       *  @name Functions to estimate general cosmological parameters
       */
      ///@{
             
      /**
       *  @brief The total neutrino mass 
       *
       *  This member function computes the total neutrino mass:
       *
       *  \f[\sum m_\nu = \Omega_\nu\cdot93.12, h^2 eV\f]
       *
       *  @return the total neutrinos mass 
       */
      virtual double neutrino_mass () const override;

      /**
       *  @brief the neutrino density at a given redshift
       *
       *  This member function computes the neutrino density parameter
       *  at a given redshift :
       *
       *  \f[ \Omega_nu(z) = \Omega_nu\frac{(1+z)^3}{E^2(z)} \f]
       *
       *  where \f$\Omega_nu\equiv\Omega_nu(z=0)\f$ and \f$E(z)\f$ is
       *  computed by cosmology::Cosmology::EE
       *
       *  @param redshift the redshift
       *
       *  @return &Omega;<SUB>&nu;</SUB>
       */
      virtual double Omega_Nu (const double redshift) const override; 

      /**
       *  @brief The dark energy density parameter
       *
       *  This member function computes the dark energy density
       *  parameter at a given redshift (see e.g. Percival 2005,
       *  Eq. 6):
       *
       *  \f[ \Omega_{\rm DE}(z) = \left( 1-\Omega_m-\Omega_r \right)
       *  \frac{1}{E^2(z)} \f]
       *
       *  where \f$\Omega_{\rm DE} \equiv\Omega_{\rm DE}(z=0)\f$,
       *  \f$E(z)\f$ is computed by cosmology::Cosmology::EE and
       *  \f$f_{\rm DE}(z)\f$ is computed by
       *  cosmology::Cosmology::f_DE
       *
       *  @param redshift the redshift
       *
       *  @return \f$\Omega_{\rm DE}(z)\f$: the dark energy density
       *  parameter
       */
      virtual double Omega_DE (const double redshift=0) const override;
      
      /**
       *  @brief The dimensionless Hubble function
       *
       *  This member function computes the dimensionless Hubble
       *  function:
       *
       *  \f[ E(z) = \frac{H(z)}{H_0} = \sqrt{\Omega_m(1+z)^3 +
       *  (1-\Omega_m-\Omega_r) + \Omega_r(1+z)^4} \f]
       *
       *  @param redshift the redshift 
       *
       *  @return the dimensionless Hubble function
       */
      virtual double EE (const double redshift) const override;
    
      /**
       *  @brief the deceleration parameter at a given redshift
       *  @param redshift the redshift
       *  @return q
       */
      virtual double qq (const double redshift) const override;

      /**
       *  @brief redshift at which the Universe begins to accelerate 
       *
       *  see e.g. de Araujo 2005
       *
       *  @return z<SUB>acc</SUB>
       */
      virtual double z_acc () const override; 

      /**
       *  @brief redshift of matter-dark energy equality
       *
       *  see e.g. de Araujo 2005
       *
       *  @return z<SUB>eq</SUB>
       */
      virtual double z_eq () const override;

      /**
       *  @brief redshift of matter-radiation equality
       *
       *  @param T_CMB the temperature of the CMB
       *
       *  @return z<SUB>eq</SUB>
       */
      virtual double z_eq_rad (const double T_CMB=2.7255) const override;
    
      /**
       *  @brief maximum absolute magnitude to have a volume-limited
       *  catalogue
       *
       *  @param z_max maximum redshift
       *  @param mag_lim magnitude limit
       *  @return r<SUB>s</SUB>
       */
      virtual double Mag_Volume_limited (const double z_max=1., const double mag_lim=-20.) const override;

      /**
       *  @brief bolometric luminosity
       *  @param redshift the redshift
       *  @param flux flux
       *  @return L<SUB>bol</SUB>
       */
      virtual double Lum_bol (const double redshift=0., const double flux=1.) const override;

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
      virtual double Redshift (const double d_c=1., const double z1_guess=0., const double z2_guess=10., const double prec=0.0001) const override;
      
      /**
       *  @brief redshift at a given cosmic time
       *  @param time cosmic time
       *  @param z1_guess minimum redshift used to search the redshift
       *  @param z2_guess maximum redshift used to search the redshift
       *  @return redshift
       */
      virtual double Redshift_time (const double time, const double z1_guess, const double z2_guess) const override;
      
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
      virtual double D_C (const double redshift) const override;  

      /**
       *  @brief the comoving transverse distance at a given redshift
       *  @param redshift the redshift
       *  @return D<SUB>M</SUB>
       */
      virtual double D_M (const double redshift) const override;

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
      virtual double D_A (const double z1, const double z2) const override;

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
      virtual double Volume (const double z1, const double z2, const double Area) const override; 

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
      virtual double Volume (const double z1, const double z2, const double RA_min, const double RA_max, const double Dec_min, const double Dec_max) const override; 
      
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
      virtual double Volume (const double zz) const override;
  
      /**
       *  @brief maximum redshift for a given volume, sky area and
       *  minimum redshift
       *
       *  @param Volume volume
       *
       *  @param Area sky area
       *
       *  @param z_min minimum redshift
       *
       *  @return redshift
       */
      virtual double max_redshift (const double Volume, const double Area, const double z_min) const override;
      
      /**
       *  @brief the mean cosmic background density
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
      virtual double rho_m (const double redshift=0., const bool unit1=false) const override;

      /**
       *  @brief the critical overdensity
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
      virtual double Delta_c (const double redshift, const std::string author="BryanNorman") const override;
      
      /**
       *  @brief The linear growth rate at a given redshift,
       *  \f$f(z)\f$
       *
       *  this function computes the following function:
       *
       *  \f[ f(z) = \frac{{\rm d}\,\ln D}{{\rm d}\,\ln a} \f]
       *
       *  @param redshift the redshift at which the linear growth rate
       *  is evaluated
       *
       *  @param prec precision used for the resolution of the
       *  differential equation in the case w<SUB>a</SUB> different
       *  than 0
       *
       *  @return the linear growth rate
       *
       *  @warning for w<SUB>a</SUB> different than 0 the current
       *  implementation does not take into account the precence of
       *  massive neutrinos
       */
      virtual double linear_growth_rate (const double redshift, const double prec=1.e-4) const override;

      /**
       *  @brief the sound horizon at the drag epoch predicted by
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
      virtual double sound_horizon_at_drag_epoch (std::string method="CAMB", const double T_CMB=par::TCMB) const override;
      
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
      virtual double As (const double sigma8) const override;

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
      virtual double sigma8_interpolated (const double redshift) const override;

      ///@}
      
    };
  }
}

#endif
