/********************************************************************
 *  Copyright (C) 2022 by Federico Marulli and Giorgio Lesci        *
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
 *  @file Headers/HaloProfile.h
 *
 *  @brief The class HaloProfile 
 *
 *  This file defines the interface of the class HaloProfile, used to
 *  handle objects of type <EM> cluster of galaxies </EM>
 *
 *  @authors Giorgio Lesci, Federico Marulli
 *
 *  @authors giorgio.lesci2@unibo.it, federico.marulli3@unibo.it
 */

#ifndef __HALOPROFILE__
#define __HALOPROFILE__ 

#include "Cosmology.h"


// ============================================================================================


namespace cbl {
  
  namespace cosmology {

    /**
     *  @class HaloProfile HaloProfile.h "Headers/HaloProfile.h"
     *
     *  @brief The class HaloProfile
     */
    class HaloProfile { 

    private :
    
      /// pointer to the input cosmology
      std::shared_ptr<cosmology::Cosmology> m_cosmology = NULL;
      
      /// cluster density profile author
      std::string m_profile_author = par::defaultString;
      
      /// halo (overdensity) definition
      std::string m_halo_def;
      
      /// overdensity factor
      double m_Delta;
      
      /// function returning the overdensity factor
      std::function<double(const double, cbl::cosmology::Cosmology, const double)> m_Delta_func;
      
      /// truncation factor  $F_t$ defining the truncation radius, that is \f$r_t = F_t*r_{\Delta}\f$ 
      double m_trunc_fact;
      
      /// if true, account for the miscentering
      bool m_miscentering;
      
      /// rms of the off-centered cluster distribution
      double m_sigma_off;
      
      /// fraction of off-centered clusters
      double m_f_off;
      
      /// lower threshold for sigma_off
      double m_sigma_off_threshold;
      
      /// lower threshold for f_off
      double m_f_off_threshold;
      
      /// minimum value where the off-centered Sigma(R) is computed
      double m_min_Smis;
      
      /// maxmimum value where the off-centered Sigma(R) is computed
      double m_max_Smis;
      
      /// if true, the miscetering is related to a single cluster. Otherwise, a population of clusters is considered
      bool m_single_profile;
      
      /// cluster 3D density profile function
      std::vector<double> (HaloProfile::*m_rho_ptr) (const std::vector<double>);
      
      /// function returning the characteristic density of the cluster profile 
      double (HaloProfile::*m_rho_s_ptr) ();
      
      /// function returning the miscentering contribution to the total cluster surface density for a given radius
      std::vector<double> (HaloProfile::*m_Sigma_mis_ptr) (const std::vector<double>);
      
      /// function returning the centered cluster surface density for a given radius
      std::vector<double> (HaloProfile::*m_Sigma_cen_ptr) (const std::vector<double>);
      
      /// function returning the centered mean cluster surface density for a given radius
      std::vector<double> (HaloProfile::*m_Sigma_mean_cen_ptr) (const std::vector<double>);
      
      /// function returning the cluster surface density for a given radius
      std::vector<double> (HaloProfile::*m_Sigma_ptr) (const std::vector<double>);
      
      /// function returning the miscentering contribution to the total cluster excess surface density for a given radius
      std::vector<double> (HaloProfile::*m_DeltaSigma_mis_ptr) (const std::vector<double>);
      
      /// function returning the centered cluster excess surface density for a given radius
      std::vector<double> (HaloProfile::*m_DeltaSigma_cen_ptr) (const std::vector<double>);
      
      /// function returning the excess cluster surface density for a given radius
      std::vector<double> (HaloProfile::*m_DeltaSigma_ptr) (const std::vector<double>);
    
      /// if true, the concentration-mass relation is set
      bool m_isSet_cM_relation = false;
      
      /// pointer returning the concentration, the one set by the user or from a concentration-mass relation function
      double (HaloProfile::*m_return_concentration) ();
      
      /// parameter A in the c-M relation
      double m_A = par::defaultDouble;
      
      /// parameter B in the c-M relation
      double m_B = par::defaultDouble;
      
      /// parameter C in the c-M relation
      double m_C = par::defaultDouble;

      /// halo mass
      double m_mass = par::defaultDouble;
      
      /// halo concentration
      double m_concentration = par::defaultDouble;
      
      /// halo redshift
      double m_redshift = par::defaultDouble;
      
      /**
       *  @name Protected member functions 
       */
      ///@{
      
      /**
       *  @brief private function that sets the halo profile parameters. 
       *  If the miscentering is
       *  considered, \f$\Sigma(R)\f$ and \f$\Delta\Sigma(R)\f$ 
       *  are expressed as (see e.g. Bellagamba et al. 2019):
       *
       *  \f$\Sigma_{\rm 1h}(R)=(1-f_{\rm off})\Sigma_{\rm cen}(R)+
       *  f_{\rm off}\Sigma_{\rm off}(R),\,\,\,\,(1)\f$
       *
       *  \f$\Delta\Sigma_{\rm 1h}(R)=(1-f_{\rm off})\Delta\Sigma_{\rm cen}(R)+
       *  f_{\rm off}\Delta\Sigma_{\rm off}(R).\,\,\,\,(2)\f$
       *
       *  @param cosmology the cosmology
       *
       *  @param redshift the redshift
       *
       *  @param conc halo concentration
       *
       *  @param Mass halo mass
       *
       *  @param Delta overdensity factor which needs to be multiplied to the critical
       *  density in order to define an overdensity
       *
       *  @param profile_author the density profile author(s); available options are:
       *  "NFW" \f$\rightarrow\f$ Navarro-Frenk-White profile;
       *  "NFW_trunc" \f$\rightarrow\f$ truncated Navarro-Frenk-White profile 
       *
       *  @param halo_def the halo definition; available options are:
       *  "vir" \f$\rightarrow\f$ all matter within the radius
       *  \f$r_{vir}\f$ for which the mean internal density is
       *  cbl::cosmology::Cosmology::Delta_c times the critical density
       *  \f$\rho_{crit}=3H^2/8\pi G\f$; "critical" \f$\rightarrow\f$ all
       *  matter within the radius \f$r_{200}\f$ for which the mean
       *  internal density is 200 times the critical density; "mean"
       *  \f$\rightarrow\f$ all matter within the radius
       *  \f$r_{\Delta}\f$ for which the mean internal density is 
       *  \f$\Delta\Omega_m\f$
       *  times the critical mean background density
       *
       *  @param trunc_fact factor \f$F_t\f$ defining the truncation
       *  radius, that is \f$r_t = F_tr_{\Delta}\f$
       *
       *  @param miscentering if true, account for the miscentering
       *
       *  @param single_profile if true, the miscetering is related
       *  to a single cluster. Otherwise, a population of clusters is considered
       *  (the latter is the case of stacked weak lensing analyses)
       *
       *  @param sigma_off if single_profile=false, this is the standard deviation of the 
       *  miscentered cluster population, \f$\sigma_{\rm off}\f$, in Mpc\f$/h\f$.
       *  Otherwise, it is the miscentering of a single cluster profile
       *
       *  @param f_off fraction of miscentered clusters, \f$f_{\rm off}\f$.
       *  Used only if single_profile=false
       *
       *  @warning If the miscentering is considered, the methods HaloProfile::Sigma()
       *  and HaloProfile::DeltaSigma() include the miscentering contribution 
       *  (i.e. Eq. 1 and Eq. 2).
       *
       *  @warning If in the cbl::cosmology::Cosmology object the \f$h\f$ units are
       *  set, both the input quantities (e.g. mass, radius) and the output quantities
       *  (e.g. \f$\Sigma\f$, \f$\Delta\Sigma\f$) are in units of \f$h\f$
       *
       */
      void m_set_profile (const cbl::cosmology::Cosmology cosmology, const double redshift, const double conc, const double Mass, const double Delta, const std::string profile_author, const std::string halo_def, const double trunc_fact=3., const bool miscentering=false, const bool single_profile=false, const double sigma_off=0.1, const double f_off=1.);
	
      /**
       *  @brief Return the concentration set by the user
       *
       *  @return the concentration
       *
       */
      double m_return_set_concentration ()
      { return m_concentration; }
      
      /**
       *  @brief The concentration-mass relation by Duffy et
       *  al. (2008): \f[c(M_h, z) = A(M_h/M_{pivot})^B\,(1+z)^C\f]
       *
       *  @return the concentration
       *
       */
      double m_concentration_Duffy ();
      
      /**
       *  @brief Characteristic density of the NFW profile
       *
       *  @return the NFW characteristic density
       *
       */
      double m_rho_s_NFW ();
      
      /**
       *  @brief Characteristic density of the truncated NFW profile
       *
       *  @return the truncated NFW characteristic density
       *
       */
      double m_rho_s_NFW_trunc ();
      
      /**
       *  @brief The NFW 3D density profile
       *
       *  @param rad radius
       *
       *  @return the NFW density profile
       *
       */
      std::vector<double> m_rho_NFW (const std::vector<double> rad);
      
      /**
       *  @brief The truncated NFW 3D density profile
       *
       *  @param rad radius
       *
       *  @return the truncated NFW density profile
       *
       */
      std::vector<double> m_rho_NFW_trunc (const std::vector<double> rad);
      
      /**
       *  @brief Surface density of the NFW profile
       *
       *  @param rad radius
       *
       *  @return the NFW surface density
       *
       */
      std::vector<double> m_Sigma_NFW (const std::vector<double> rad);
      
      /**
       *  @brief Mean surface density of the NFW profile
       *
       *  @param rad radius
       *
       *  @return the NFW surface density
       *
       */
      std::vector<double> m_Sigma_mean_NFW (const std::vector<double> rad);
      
      /**
       *  @brief Surface density of the truncated NFW profile
       *
       *  @param rad radius
       *
       *  @return the truncated NFW surface density
       *
       */
      std::vector<double> m_Sigma_NFW_trunc (const std::vector<double> rad);
      
      /**
       *  @brief Mean surface density of the truncated NFW profile
       *
       *  @param rad radius
       *
       *  @return the truncated NFW surface density
       *
       */
      std::vector<double> m_Sigma_mean_NFW_trunc (const std::vector<double> rad);
      
      /**
       *  @brief the mis-centered surface density profile of the halo
       *  (see e.g. Bellagamba et al. 2019).
       *
       *  @param rad radius
       *
       *  @return the halo surface density profile in \f$M_{\odot} / \f$pc^2
       */
      std::vector<double> m_Sigma_mis (const std::vector<double> rad);
      
      /**
       *  @brief Surface density
       *  including the miscentering contribution
       *
       *  @param rad radius
       *
       *  @return the surface density
       *
       */
      std::vector<double> m_Sigma_including_miscentering (const std::vector<double> rad);
      
      /**
       *  @brief Excess surface density
       *  without the miscentering contribution
       *
       *  @param rad radius
       *
       *  @return the excess surface density
       *
       */
      std::vector<double> m_DeltaSigma_cen (const std::vector<double> rad);
      
      /**
       *  @brief the mis-centered excess surface density profile of the halo
       *  (see e.g. Bellagamba et al. 2019).
       *
       *  @param rad radius
       *
       *  @return the halo excess surface density profile in \f$M_{\odot} / \f$pc^2
       */
      std::vector<double> m_DeltaSigma_mis (const std::vector<double> rad);
      
      /**
       *  @brief Excess surface density
       *  including the miscentering contribution
       *
       *  @param rad radius
       *
       *  @return the excess surface density
       *
       */
      std::vector<double> m_DeltaSigma_including_miscentering (const std::vector<double> rad);
      
      /**
       *  @brief set the concentration-mass relation. 
       *  A function is set, computing the concentration of a dark matter
       *  halo of a given a mass, at a given redshift; the models
       *  implemented are the following:
       *
       *  - Duffy et al. 2008:
       *  \f[c(M_h, z) = A(M_h/M_{pivot})^B\,(1+z)^C\f]
       *
       *
       *  @param cM_author author(s) who proposed the 
       *  concentration-mass relation. Possibilities are:
       *  "Duffy" (Duffy et al. 2008)
       *
       *  @warning the Duffy et al. concentrantion-mass relation
       *  refers to the 0<z<2 redshift range, obtained from their full
       *  samples (see Table 1 of Duffy et al. 2008).
       *
       */
      void m_set_cM_relation (const std::string cM_author="Duffy");

      ///@}
      
      /// 
      
    public:
    
      /**
       *  @name Constructors/destructors
       */
      ///@{

      /**
       *  @brief default constructor
       *  
       */
      HaloProfile () = default;

      /**
       *  @brief constructor that sets the halo profile parameters. 
       *  If the miscentering is
       *  considered, \f$\Sigma(R)\f$ and \f$\Delta\Sigma(R)\f$ 
       *  are expressed as (see e.g. Bellagamba et al. 2019):
       *
       *  \f$\Sigma_{\rm 1h}(R)=(1-f_{\rm off})\Sigma_{\rm cen}(R)+
       *  f_{\rm off}\Sigma_{\rm off}(R),\,\,\,\,(1)\f$
       *
       *  \f$\Delta\Sigma_{\rm 1h}(R)=(1-f_{\rm off})\Delta\Sigma_{\rm cen}(R)+
       *  f_{\rm off}\Delta\Sigma_{\rm off}(R).\,\,\,\,(2)\f$
       *
       *  @param cosmology the cosmology
       *
       *  @param redshift the redshift
       *
       *  @param conc halo concentration
       *
       *  @param Mass halo mass
       *
       *  @param Delta overdensity factor which needs to be multiplied to the critical
       *  density in order to define an overdensity
       *
       *  @param profile_author the density profile author(s); available options are:
       *  "NFW" \f$\rightarrow\f$ Navarro-Frenk-White profile;
       *  "NFW_trunc" \f$\rightarrow\f$ truncated Navarro-Frenk-White profile 
       *
       *  @param halo_def the halo definition; available options are:
       *  "vir" \f$\rightarrow\f$ all matter within the radius
       *  \f$r_{vir}\f$ for which the mean internal density is
       *  cbl::cosmology::Cosmology::Delta_c times the critical density
       *  \f$\rho_{crit}=3H^2/8\pi G\f$; "critical" \f$\rightarrow\f$ all
       *  matter within the radius \f$r_{200}\f$ for which the mean
       *  internal density is 200 times the critical density; "mean"
       *  \f$\rightarrow\f$ all matter within the radius
       *  \f$r_{\Delta}\f$ for which the mean internal density is 
       *  \f$\Delta\Omega_m\f$
       *  times the critical mean background density
       *
       *  @param trunc_fact factor \f$F_t\f$ defining the truncation
       *  radius, that is \f$r_t = F_tr_{\Delta}\f$
       *
       *  @param miscentering if true, account for the miscentering
       *
       *  @param single_profile if true, the miscetering is related
       *  to a single cluster. Otherwise, a population of clusters is considered
       *  (the latter is the case of stacked weak lensing analyses)
       *
       *  @param sigma_off if single_profile=false, this is the standard deviation of the 
       *  miscentered cluster population, \f$\sigma_{\rm off}\f$, in Mpc\f$/h\f$.
       *  Otherwise, it is the miscentering of a single cluster profile
       *
       *  @param f_off fraction of miscentered clusters, \f$f_{\rm off}\f$.
       *  Used only if single_profile=false
       *
       *  @warning If the miscentering is considered, the methods HaloProfile::Sigma()
       *  and HaloProfile::DeltaSigma() include the miscentering contribution 
       *  (i.e. Eq. 1 and Eq. 2).
       *
       *  @warning If in the cbl::cosmology::Cosmology object the \f$h\f$ units are
       *  set, both the input quantities (e.g. mass, radius) and the output quantities
       *  (e.g. \f$\Sigma\f$, \f$\Delta\Sigma\f$) are in units of \f$h\f$
       *
       */
      HaloProfile (const cbl::cosmology::Cosmology cosmology, const double redshift, const double conc, const double Mass, const double Delta, const std::string profile_author, const std::string halo_def, const double trunc_fact=3., const bool miscentering=false, const bool single_profile=false, const double sigma_off=0.1, const double f_off=1.);
      
      /**
       *  @brief constructor that sets the halo profile parameters,
       *  assuming a
       *  concentration-mass relation. If the miscentering is
       *  considered, \f$\Sigma(R)\f$ and \f$\Delta\Sigma(R)\f$ 
       *  are expressed as (see e.g. Bellagamba et al. 2019):
       *
       *  \f$\Sigma_{\rm 1h}(R)=(1-f_{\rm off})\Sigma_{\rm cen}(R)+
       *  f_{\rm off}\Sigma_{\rm off}(R),\,\,\,\,(1)\f$
       *
       *  \f$\Delta\Sigma_{\rm 1h}(R)=(1-f_{\rm off})\Delta\Sigma_{\rm cen}(R)+
       *  f_{\rm off}\Delta\Sigma_{\rm off}(R).\,\,\,\,(2)\f$
       *
       *  @param cosmology the cosmology
       *
       *  @param redshift the redshift
       *
       *  @param cM_author author(s) who proposed the 
       *  concentration-mass relation. Possibilities are:
       *  "Duffy" (Duffy et al. 2008)
       *
       *  @param Mass halo mass
       *
       *  @param Delta overdensity factor which needs to be multiplied to the critical
       *  density in order to define an overdensity
       *
       *  @param profile_author the density profile author(s); available options are:
       *  "NFW" \f$\rightarrow\f$ Navarro-Frenk-White profile;
       *  "NFW_trunc" \f$\rightarrow\f$ truncated Navarro-Frenk-White profile 
       *
       *  @param halo_def the halo definition; available options are:
       *  "vir" \f$\rightarrow\f$ all matter within the radius
       *  \f$r_{vir}\f$ for which the mean internal density is
       *  cbl::cosmology::Cosmology::Delta_c times the critical density
       *  \f$\rho_{crit}=3H^2/8\pi G\f$; "critical" \f$\rightarrow\f$ all
       *  matter within the radius \f$r_{200}\f$ for which the mean
       *  internal density is 200 times the critical density; "mean"
       *  \f$\rightarrow\f$ all matter within the radius
       *  \f$r_{\Delta}\f$ for which the mean internal density is 
       *  \f$\Delta\Omega_m\f$
       *  times the critical mean background density
       *
       *  @param trunc_fact factor \f$F_t\f$ defining the truncation
       *  radius, that is \f$r_t = F_tr_{\Delta}\f$
       *
       *  @param miscentering if true, account for the miscentering
       *
       *  @param single_profile if true, the miscetering is related
       *  to a single cluster. Otherwise, a population of clusters is considered
       *  (the latter is the case of stacked weak lensing analyses)
       *
       *  @param sigma_off if single_profile=false, this is the standard deviation of the 
       *  miscentered cluster population, \f$\sigma_{\rm off}\f$, in Mpc\f$/h\f$.
       *  Otherwise, it is the miscentering of a single cluster profile
       *
       *  @param f_off fraction of miscentered clusters, \f$f_{\rm off}\f$.
       *  Used only if single_profile=false
       *
       *  @warning the Duffy et al. concentrantion-mass relation
       *  refers to the 0<z<2 redshift range, obtained from their full
       *  samples (see Table 1 of Duffy et al. 2008).
       *
       *  @warning If the miscentering is considered, the methods HaloProfile::Sigma()
       *  and HaloProfile::DeltaSigma() include the miscentering contribution 
       *  (i.e. Eq. 1 and Eq. 2).
       *
       *  @warning If in the cbl::cosmology::Cosmology object the \f$h\f$ units are
       *  set, both the input quantities (e.g. mass, radius) and the output quantities
       *  (e.g. \f$\Sigma\f$, \f$\Delta\Sigma\f$) are in units of \f$h\f$
       *
       */
      HaloProfile (const cbl::cosmology::Cosmology cosmology, const double redshift, const std::string cM_author, const double Mass, const double Delta, const std::string profile_author, const std::string halo_def, const double trunc_fact=3., const bool miscentering=false, const bool single_profile=false, const double sigma_off=0.1, const double f_off=1.);
      
      /**
       *  @brief default destructor
       */
      ~HaloProfile () = default;

      ///@}
  
  
      /**
       *  @name Member functions used to get the private members or
       *  to compute the cluster halo profiles
       */
      ///@{
      
      /**
       *  @brief get the concentration
       *
       *  @return the concentration of the halo
       *
       *  @warning if a concentration-mass relation
       *  is set, the output concentration is derived
       *  from such relation
       */
      double concentration ()
      { return (this->*m_return_concentration)(); }
	
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
       *  @brief get the overdensity
       *
       *  @return the overdensity which needs to be multiplied 
       *  to the critical density in order to define an overdensity. 
       *  Given an overdensity factor \f$\Delta\f$, this function
       *  returns \f$\Delta\f$ for a critical overdensity, 
       *  \f$\Delta\Omega_{\rm m}(z)\f$ for a mean overdensity,
       *  and cbl::cosmology::Cosmology::Delta\_c for a
       *  overdensity defined with respect to the critical density.
       *
       */
      double Delta ()
      { return m_Delta_func(m_Delta,*m_cosmology,m_redshift); }
      
      /**
       *  @brief the 3D halo density profile. 
       *
       *  @param rad radius
       *
       *  @return the halo density profile
       *
       *  @warning If in the cbl::cosmology::Cosmology object the \f$h\f$ units are
       *  set, both the input quantities (e.g. mass, radius) and the output quantities
       *  are in units of \f$h\f$
       */
      std::vector<double> density_profile_3D (const std::vector<double> rad);
      
      /**
       *  @brief the characteristic density of the
       *  3D halo density profile. 
       *
       *  @return the halo density profile
       */
      double rho_s ();
      
      /**
       *  @brief the total surface density profile of the halo. 
       *
       *  @param rad radius
       *
       *  @return the halo surface density profile in \f$M_{\odot} / \f$pc^2
       *
       *  @warning If in the cbl::cosmology::Cosmology object the \f$h\f$ units are
       *  set, both the input quantities (e.g. mass, radius) and the output quantities
       *  are in units of \f$h\f$
       */
      std::vector<double> Sigma (const std::vector<double> rad);
      
      /**
       *  @brief the miscentering contribution to
       *  the total surface density profile of the halo.
       *
       *  @param rad radius
       *
       *  @return the halo surface density profile in \f$M_{\odot} / \f$pc^2
       *
       *  @warning If in the cbl::cosmology::Cosmology object the \f$h\f$ units are
       *  set, both the input quantities (e.g. mass, radius) and the output quantities
       *  are in units of \f$h\f$
       */
      std::vector<double> Sigma_mis (const std::vector<double> rad);
      
      /**
       *  @brief the centered contribution to
       *  the total surface density profile of the halo.
       *
       *  @param rad radius
       *
       *  @return the halo surface density profile in \f$M_{\odot} / \f$pc^2
       *
       *  @warning If in the cbl::cosmology::Cosmology object the \f$h\f$ units are
       *  set, both the input quantities (e.g. mass, radius) and the output quantities
       *  are in units of \f$h\f$
       */
      std::vector<double> Sigma_cen (const std::vector<double> rad);     
      
      /**
       *  @brief the 2-halo contribution to
       *  the total surface density profile of the halo, computed
       *  by assuming a halo bias model.
       *
       *  @param rad radius
       *
       *  @param bias_author author(s) who proposed the bias; valid authors are: 
       *  ST99 (Sheth & Tormen 1999), SMT01 (Sheth, Mo & Tormen 2001), 
       *  SMT01_WL04 (Sheth, Mo & Tormen 2001 with the correction of Warren 2004), 
       *  Tinker (Tinker et al. 2010)
       *
       *  @param method_Pk method used for the computation of the power spectrum.
       *  Valid choices for method_Pk
       *  are: CAMB [http://camb.info/], CLASS
       *  [http://class-code.net/], MPTbreeze-v1
       *  [http://arxiv.org/abs/1207.1465], EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       *
       *  @param interp_type method to interpolate the power spectrum.
       *  Possibilities are: "Linear", "Poly", "Spline", "Spline_periodic", 
       *  "Akima", "Akima_periodic", "Steffen"
       *
       *  @param NL false \f$\rightarrow\f$ linear power spectrum;
       *  true \f$\rightarrow\f$ non-linear power spectrum
       *
       *  @return the halo surface density profile in \f$M_{\odot} / \f$pc^2
       *
       *  @warning If in the cbl::cosmology::Cosmology object the \f$h\f$ units are
       *  set, both the input quantities (e.g. mass, radius) and the output quantities
       *  are in units of \f$h\f$
       */
      std::vector<double> Sigma_2h (const std::vector<double> rad, const std::string bias_author, const std::string method_Pk="EisensteinHu", const std::string interp_type="Linear", const bool NL=false);     
      
      /**
       *  @brief the 2-halo contribution to
       *  the total surface density profile of the halo.
       *
       *  @param rad radius
       *
       *  @param bias the halo bias value
       *
       *  @param method_Pk method used for the computation of the power spectrum.
       *  Valid choices for method_Pk
       *  are: CAMB [http://camb.info/], CLASS
       *  [http://class-code.net/], MPTbreeze-v1
       *  [http://arxiv.org/abs/1207.1465], EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       *
       *  @param interp_type method to interpolate the power spectrum.
       *  Possibilities are: "Linear", "Poly", "Spline", "Spline_periodic", 
       *  "Akima", "Akima_periodic", "Steffen"
       *
       *  @param NL false \f$\rightarrow\f$ linear power spectrum;
       *  true \f$\rightarrow\f$ non-linear power spectrum
       *
       *  @return the halo surface density profile in \f$M_{\odot} / \f$pc^2
       *
       *  @warning If in the cbl::cosmology::Cosmology object the \f$h\f$ units are
       *  set, both the input quantities (e.g. mass, radius) and the output quantities
       *  are in units of \f$h\f$
       */
      std::vector<double> Sigma_2h (const std::vector<double> rad, const double bias, const std::string method_Pk="EisensteinHu", const std::string interp_type="Linear", const bool NL=false);
      
      /**
       *  @brief the total excess surface density profile of the halo. 
       *
       *  @param rad radius
       *
       *  @return the halo excess surface density profile in \f$M_{\odot} / \f$pc^2
       *
       *  @warning If in the cbl::cosmology::Cosmology object the \f$h\f$ units are
       *  set, both the input quantities (e.g. mass, radius) and the output quantities
       *  are in units of \f$h\f$
       *
       */
      std::vector<double> DeltaSigma (const std::vector<double> rad);
      
      /**
       *  @brief the miscentering contribution to
       *  the excess surface density profile of the halo. 
       *
       *  @param rad radius
       *
       *  @return the halo excess surface density profile in \f$M_{\odot} / \f$pc^2
       *
       *  @warning If in the cbl::cosmology::Cosmology object the \f$h\f$ units are
       *  set, both the input quantities (e.g. mass, radius) and the output quantities
       *  are in units of \f$h\f$
       */
      std::vector<double> DeltaSigma_mis (const std::vector<double> rad);
      
      /**
       *  @brief the centered contribution to
       *  the excess surface density profile of the halo. 
       *
       *  @param rad radius
       *
       *  @return the halo excess surface density profile in \f$M_{\odot} / \f$pc^2
       *
       *  @warning If in the cbl::cosmology::Cosmology object the \f$h\f$ units are
       *  set, both the input quantities (e.g. mass, radius) and the output quantities
       *  are in units of \f$h\f$
       */
      std::vector<double> DeltaSigma_cen (const std::vector<double> rad);
      
      /**
       *  @brief the 2-halo contribution to
       *  the excess surface density profile of the halo, computed
       *  by assuming a halo bias model.
       *
       *  @param rad radius
       *
       *  @param bias_author author(s) who proposed the bias; valid authors are: 
       *  ST99 (Sheth & Tormen 1999), SMT01 (Sheth, Mo & Tormen 2001), 
       *  SMT01_WL04 (Sheth, Mo & Tormen 2001 with the correction of Warren 2004), 
       *  Tinker (Tinker et al. 2010)
       *
       *  @param method_Pk method used for the computation of the power spectrum.
       *  Valid choices for method_Pk
       *  are: CAMB [http://camb.info/], CLASS
       *  [http://class-code.net/], MPTbreeze-v1
       *  [http://arxiv.org/abs/1207.1465], EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       *
       *  @param interp_type method to interpolate the power spectrum.
       *  Possibilities are: "Linear", "Poly", "Spline", "Spline_periodic", 
       *  "Akima", "Akima_periodic", "Steffen"
       *
       *  @param NL false \f$\rightarrow\f$ linear power spectrum;
       *  true \f$\rightarrow\f$ non-linear power spectrum
       *
       *  @return the halo excess surface density profile in \f$M_{\odot} / \f$pc^2
       *
       *  @warning If in the cbl::cosmology::Cosmology object the \f$h\f$ units are
       *  set, both the input quantities (e.g. mass, radius) and the output quantities
       *  are in units of \f$h\f$
       */
      std::vector<double> DeltaSigma_2h (const std::vector<double> rad, const std::string bias_author, const std::string method_Pk="EisensteinHu", const std::string interp_type="Linear", const bool NL=false);     
      
      /**
       *  @brief the 2-halo contribution to
       *  the excess surface density profile of the halo.
       *
       *  @param rad radius
       *
       *  @param bias the halo bias value
       *
       *  @param method_Pk method used for the computation of the power spectrum.
       *  Valid choices for method_Pk
       *  are: CAMB [http://camb.info/], CLASS
       *  [http://class-code.net/], MPTbreeze-v1
       *  [http://arxiv.org/abs/1207.1465], EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       *
       *  @param interp_type method to interpolate the power spectrum.
       *  Possibilities are: "Linear", "Poly", "Spline", "Spline_periodic", 
       *  "Akima", "Akima_periodic", "Steffen"
       *
       *  @param NL false \f$\rightarrow\f$ linear power spectrum;
       *  true \f$\rightarrow\f$ non-linear power spectrum
       *
       *  @return the halo excess surface density profile in \f$M_{\odot} / \f$pc^2
       *
       *  @warning If in the cbl::cosmology::Cosmology object the \f$h\f$ units are
       *  set, both the input quantities (e.g. mass, radius) and the output quantities
       *  are in units of \f$h\f$
       */
      std::vector<double> DeltaSigma_2h (const std::vector<double> rad, const double bias, const std::string method_Pk="EisensteinHu", const std::string interp_type="Linear", const bool NL=false);
      
      
      /**
       *  @brief the Fourier transform of the normalised halo density
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
       *  @return the halo density profile
       */
      double density_profile_FourierSpace (const double kk);
      
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
       *  @name Member functions used to set the private members 
       */
      ///@{
      
      /**
       *  @brief set the cosmological model
       *
       *  @param cosmology the cosmology
       *
       */
      void set_cosmology (const cbl::cosmology::Cosmology cosmology);
      
      /**
       *  @brief set the private member \e m_mass
       *  @param mass the mass of the cluster
       */
      void set_mass (const double mass=par::defaultDouble)
      { m_mass = mass; }
      
      /**
       *  @brief set the private member \e m_concentration
       *  @param conc concentration
       */
      void set_concentration (const double conc=par::defaultDouble)
      { m_concentration = conc; }
      
      /**
       *  @brief set the private member \e m_f_off
       *  @param f_off f_off
       */
      void set_f_off (const double f_off=par::defaultDouble)
      { m_f_off = f_off; }
      
      /**
       *  @brief set the private member \e m_sigma_off
       *  @param sigma_off sigma_off
       */
      void set_sigma_off (const double sigma_off=par::defaultDouble)
      { m_sigma_off = sigma_off; }
      
      /**
       *  @brief set the private member \e m_trunc_fact
       *  @param trunc_fact trunc_fact
       */
      void set_trunc_fact (const double trunc_fact=par::defaultDouble)
      { m_trunc_fact = trunc_fact; }
      
      ///@}      
    
    };
  }
}

#endif

