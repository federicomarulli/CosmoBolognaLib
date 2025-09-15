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

#include "Bias.h"


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
    
      /// Newton's constant in [Mpc^3 Msun^-1 s^-2]
      double m_G_newt = cbl::par::GN / pow(3.086e+16 * 1.e6, 3.) * 1.98855e30;
      
      /// light speed in [Mpc s^-1]
      double m_c_light = cbl::par::cc * 1.e3 / (3.086e+16 * 1.e6);
      
      /// pointer to the input cosmology
      std::shared_ptr<cosmology::Cosmology> m_cosmology = NULL;
      
      /// if true, include the 2-halo 
      bool m_two_halo;
      
      /// cluster density profile author
      std::string m_profile_author = par::defaultString;
      
      /// halo (overdensity) definition
      std::string m_halo_def;
      
      /// overdensity factor
      double m_Delta;
      
      /// function returning the overdensity factor
      std::function<double(const double, std::shared_ptr<cosmology::Cosmology>, const double)> m_Delta_func;
      
      /// truncation factor  $F_t$ defining the truncation radius, that is \f$r_t = F_t*r_{\Delta}\f$ 
      double m_trunc_fact;
      
      /// if true, account for the miscentering
      bool m_miscentering;
      
      /// rms of the off-centered cluster distribution
      double m_sigma_off;
      
      /// fraction of off-centered clusters
      double m_f_off;
      
      /// anisotropic boost factor multiplying the 2-halo
      double m_AB_fact;
      
      /// alpha_0 in Einasto's alpha profile
      double m_alpha_0;
      
      /// alpha_nu in Einasto's alpha profile
      double m_alpha_nu;
      
      /// beta in the DK14 profile
      double m_beta;
      
      /// gamma_0 in the DK14 profile
      double m_gamma_0;
      
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
      std::vector<double> (HaloProfile::*m_rho_ptr) (const std::vector<double>, const double);
      
      /// function returning the centered cluster surface density for a given radius
      std::vector<double> (HaloProfile::*m_Sigma_cen_ptr) (const std::vector<double>);
      
      /// function returning the cluster surface density for a given radius
      std::vector<double> (HaloProfile::*m_Sigma_ptr) (const std::vector<double>);
      
      /// function returning the excess cluster surface density for a given radius
      std::vector<double> (HaloProfile::*m_DeltaSigma_ptr) (const std::vector<double>);
      
      /// pointer to the interpolated unbiased 2-halo Sigma(R)
      std::shared_ptr<glob::FuncGrid> m_Sigma_2h_unbiased_interp = NULL;
    
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
      
      /// normalisation of the outer term in DK14
      double m_b_e;
      
      /// slope of the outer term in DK14
      double m_s_e;
      
      /// halo bias author
      std::string m_bias_author;
      
      /// halo bias
      double m_bias;
      
      /// method for the power spectrum computation
      std::string m_method_Pk;
      
      /**
       *  @name Protected member functions 
       */
      ///@{
      
      /**
       *  @brief private function that sets the halo profile
       *  parameters.  If the miscentering is considered,
       *  \f$\Sigma(R)\f$ and \f$\Delta\Sigma(R)\f$ are expressed as
       *  (see e.g. Bellagamba et al. 2019):
       *
       *  \f$\Sigma_{\rm 1h}(R)=(1-f_{\rm off})\Sigma_{\rm cen}(R)+
       *  f_{\rm off}\Sigma_{\rm off}(R),\,\,\,\,(1)\f$
       *
       *  \f$\Delta\Sigma_{\rm 1h}(R)=(1-f_{\rm off})\Delta\Sigma_{\rm cen}(R)+
       *  f_{\rm off}\Delta\Sigma_{\rm off}(R).\,\,\,\,(2)\f$
       *
       *  @param two_halo if true, include the 2-halo contribution
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
       *  @param profile_author the density profile author(s);
       *  available options are: "NFW" \f$\rightarrow\f$
       *  Navarro-Frenk-White profile (Navarro et al. 1997);
       *  "NFW_trunc" \f$\rightarrow\f$ truncated Navarro-Frenk-White
       *  profile, aka BMO (Baltz et al. 2009); "Einasto"
       *  \f$\rightarrow\f$ Einasto 1965; "DK14" \f$\rightarrow\f$
       *  Diemer & Kravtsov 2014
       *
       *  @param halo_def the halo definition; available options are:
       *  "vir" \f$\rightarrow\f$ all matter within the radius
       *  \f$r_{vir}\f$ for which the mean internal density is
       *  cbl::cosmology::Cosmology::Delta_c times the critical
       *  density \f$\rho_{crit}=3H^2/8\pi G\f$; "critical"
       *  \f$\rightarrow\f$ all matter within the radius \f$r_{\Delta}\f$
       *  for which the mean internal density is \f$\Delta\f$ times the
       *  critical density; "mean" \f$\rightarrow\f$ all matter within
       *  the radius \f$r_{\Delta}\f$ for which the mean internal
       *  density is \f$\Delta\Omega_m\f$ times the critical mean
       *  background density
       *
       *  @param trunc_fact factor \f$F_t\f$ defining the truncation
       *  radius, that is \f$r_t = F_tr_{\Delta}\f$
       *
       *  @param miscentering if true, account for the miscentering
       *
       *  @param single_profile if true, the miscetering is related to
       *  a single cluster. Otherwise, a population of clusters is
       *  considered (the latter is the case of stacked weak lensing
       *  analyses)
       *
       *  @param sigma_off if single_profile=false, this is the
       *  standard deviation of the miscentered cluster population,
       *  \f$\sigma_{\rm off}\f$, in Mpc\f$/h\f$.  Otherwise, it is
       *  the miscentering of a single cluster profile
       *
       *  @param f_off fraction of miscentered clusters, \f$f_{\rm
       *  off}\f$.  Used only if single_profile=false
       *
       *  @param b_e normalisation of the DK14 outer profile
       *
       *  @param s_e slope of the DK14 outer profile
       *
       *  @param AB_fact anisotropic boost factor, \f$\sigma_{\rm
       *  AB}\f$, affecting the 2-halo term and/or the infall term in
       *  DK14. Specifically, such terms are multiplied by
       *  \f$(1+\sigma_{\rm AB})\f$
       *
       *  @param bias author(s) who proposed the bias; valid authors
       *  are: ST99 (Sheth & Tormen 1999), SMT01 (Sheth, Mo & Tormen
       *  2001), SMT01_WL04 (Sheth, Mo & Tormen 2001 with the
       *  correction of Warren 2004), Tinker (Tinker et al. 2010)
       *
       *  @param method_Pk method used for the computation of the
       *  power spectrum.  Valid choices for method_Pk are: CAMB
       *  [http://camb.info/], CLASS [http://class-code.net/],
       *  MPTbreeze-v1 [http://arxiv.org/abs/1207.1465], EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       *
       *  @param alpha_0 \f$\alpha_0\f$ parameter in the equation
       *  \f$ \alpha(\nu) = \alpha_0 + \alpha_\nu\nu^2 \f$, where \f$\nu\f$ 
       *  is the virial peak height. \f$\alpha(\nu)\f$ appears in both
       *  Einasto and DK14 profiles. Following Gao et al. (2008), 
       *  \f$\alpha_0=0.155\f$ by default.
       *
       *  @param alpha_nu \f$\alpha_\nu\f$ parameter in the equation
       *  \f$ \alpha(\nu) = \alpha_0 + \alpha_\nu\nu^2 \f$, where \f$\nu\f$ 
       *  is the virial peak height. \f$\alpha(\nu)\f$ appears in both
       *  Einasto and DK14 profiles. Following Gao et al. (2008), 
       *  \f$\alpha_\nu=0.0095\f$ by default.
       *
       *  @param beta \f$\beta\f$ parameter in the DK14 profile. Following
       *  DK14, \f$\beta=4\f$ by default.
       *
       *  @param gamma_0 \f$\gamma_0\f$ parameter in the expression
       *  \f$\gamma=\gamma_0\nu\f$ appearing in the DK14 profile. Following
       *  DK14, \f$\gamma_0=4\f$ by default.
       *
       *  @warning If the miscentering is considered, the methods
       *  HaloProfile::Sigma() and HaloProfile::DeltaSigma() include
       *  the miscentering contribution (i.e. Eq. 1 and Eq. 2).
       *
       *  @warning If in the cbl::cosmology::Cosmology object the
       *  \f$h\f$ units are set, both the input quantities (e.g. mass,
       *  radius) and the output quantities (e.g. \f$\Sigma\f$,
       *  \f$\Delta\Sigma\f$) are in units of \f$h\f$
       *
       */
      void m_set_profile (const bool two_halo, const double redshift, const double conc, const double Mass, const double Delta, const std::string profile_author, const std::string halo_def, const double trunc_fact=3., const bool miscentering=false, const bool single_profile=false, const double sigma_off=0.1, const double f_off=1., const double b_e=1., const double s_e=1.5, const double AB_fact=0., const std::string bias="Tinker", const std::string method_Pk="EisensteinHu", const double alpha_0=0.155, const double alpha_nu=0.0095, const double beta=4., const double gamma_0=4.);
	
      /**
       *  @brief Return the concentration set by the user
       *
       *  @return the concentration
       *
       */
      double m_return_set_concentration ()
      { return m_concentration; }
      
      /**
       *  @brief The virial peak height
       *
       *  @return the virial peak height
       *
       */
      double m_nu_vir ();
      
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
       *  @brief Characteristic density of the Einasto profile
       *
       *  @param nu peak height
       *
       *  @return the Einasto characteristic density
       *
       */
      double m_rho_s_Einasto (const double nu=0.);
      
      /**
       *  @brief Characteristic density of the DK14 profile
       *
       *  @param nu peak height
       *
       *  @return the DK14 characteristic density
       *
       */
      double m_rho_s_DK14 (const double nu=0.);
      
      /**
       *  @brief The NFW 3D density profile
       *
       *  @param rad radius
       *
       *  @param nuisance unused parameter
       *
       *  @return the NFW density profile
       *
       */
      std::vector<double> m_rho_NFW (const std::vector<double> rad, const double nuisance=-1);
      
      /**
       *  @brief The truncated NFW 3D density profile
       *
       *  @param rad radius
       *
       *  @param nuisance unused parameter
       *
       *  @return the truncated NFW density profile
       *
       */
      std::vector<double> m_rho_NFW_trunc (const std::vector<double> rad, const double nuisance=-1);
      
      /**
       *  @brief The Einasto 3D density profile
       *
       *  @param rad radius
       *
       *  @param nu_in peak height
       *
       *  @return the Einasto density profile
       *
       */
      std::vector<double> m_rho_Einasto (const std::vector<double> rad, const double nu_in=-1);
      
      /**
       *  @brief The DK14 3D density profile
       *
       *  @param rad radius
       *
       *  @param nu_in peak height
       *
       *  @return the DK14 density profile
       *
       */
      std::vector<double> m_rho_DK14 (const std::vector<double> rad, const double nu_in=-1);
      
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
       *  @brief Surface density of the truncated NFW profile
       *
       *  @param rad radius
       *
       *  @return the truncated NFW surface density
       *
       */
      std::vector<double> m_Sigma_NFW_trunc (const std::vector<double> rad);
      
      /**
       *  @brief Generic surface density profile, obtained by
       *  integrating the 3D density
       *
       *  @param rad radius
       *
       *  @return the surface density
       *
       */
      std::vector<double> m_Sigma_generic (const std::vector<double> rad);
      
      /**
       *  @brief Centered mean surface density
       *
       *  @param rad radius
       *
       *  @return the centered mean surface density
       *
       */
      std::vector<double> m_Sigma_mean_cen (const std::vector<double> rad);
      
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
       *  @brief unbiased 2-halo surface density
       *
       *  @param rad radius
       *
       *  @return the 2-halo surface density
       *
       */
      std::vector<double> m_Sigma_2h_unbiased (const std::vector<double> rad);
      
      /**
       *  @brief set the concentration-mass relation.  A function is
       *  set, computing the concentration of a dark matter halo of a
       *  given a mass, at a given redshift; the models implemented
       *  are the following:
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
       *  @brief constructor
       *
       *  @param cosmology pointer to an object of class Cosmology
       */
      HaloProfile (std::shared_ptr<cbl::cosmology::Cosmology> cosmology)
      : m_cosmology(std::move(cosmology)) {}

      /**
       *  @brief constructor that sets the halo profile parameters.
       *  If the miscentering is considered, \f$\Sigma(R)\f$ and
       *  \f$\Delta\Sigma(R)\f$ are expressed as (see e.g. Bellagamba
       *  et al. 2019):
       *
       *  \f$\Sigma(R)=(1-f_{\rm off})\Sigma_{\rm cen}(R)+
       *  f_{\rm off}\Sigma_{\rm off}(R),\,\,\,\,(1)\f$
       *
       *  \f$\Delta\Sigma(R)=(1-f_{\rm off})\Delta\Sigma_{\rm
       *  cen}(R)+ f_{\rm off}\Delta\Sigma_{\rm off}(R).\,\,\,\,(2)\f$
       *
       *  where \f$\Sigma_{\rm cen}\f$ and \f$\Sigma_{\rm cen}\f$ are the
       *  centred and miscentred surface densities, respectively. The
       *  miscentring probability follows a Rayleigh distribution.
       * 
       *  The miscentring also affects the 2-halo term. Indeed, the
       *  general expression of \f$\Sigma_{\rm cen}\f$ is:
       *
       *  \f$\Sigma_{\rm cen} = \Sigma_{\rm 1h} + \Sigma_{\rm 2h},\f$
       *
       *  where \f$\Sigma_{\rm 1h}\f$ and \f$\Sigma_{\rm 2h}\f$ are the
       *  1-halo and 2-halo surface densities, respectively.
       *  \f$\Sigma_{\rm 2h}\f$ is expressed as follows (Oguri & Takada 2011):
       *
       *  \f$\Sigma_{\rm 2h}(R) = \int\,\frac{l{\rm d}l}{2\pi}
       *  J_0(l\theta)\frac{\bar{\rho}_{\rm m}(z)b(M,z)}{(1+z)^3D_{\rm l}^2(z)}
       *  P(k_l,z) \f$
       *
       *  @param cosmology the cosmology
       *
       *  @param two_halo if true, include the 2-halo contribution
       *
       *  @param redshift the redshift
       *
       *  @param conc halo concentration
       *
       *  @param Mass halo mass
       *
       *  @param Delta overdensity factor which needs to be multiplied
       *  to the critical density in order to define an overdensity
       *
       *  @param profile_author the density profile author(s);
       *  available options are: "NFW" \f$\rightarrow\f$
       *  Navarro-Frenk-White profile (Navarro et al. 1997);
       *  "NFW_trunc" \f$\rightarrow\f$ truncated Navarro-Frenk-White
       *  profile, aka BMO (Baltz et al. 2009); "Einasto"
       *  \f$\rightarrow\f$ Einasto 1965; "DK14" \f$\rightarrow\f$
       *  Diemer & Kravtsov 2014
       *
       *  @param halo_def the halo definition; available options are:
       *  "vir" \f$\rightarrow\f$ all matter within the radius
       *  \f$r_{vir}\f$ for which the mean internal density is
       *  cbl::cosmology::Cosmology::Delta_c times the critical
       *  density \f$\rho_{crit}=3H^2/8\pi G\f$; "critical"
       *  \f$\rightarrow\f$ all matter within the radius \f$r_{\Delta}\f$
       *  for which the mean internal density is \f$\Delta\f$ times the
       *  critical density; "mean" \f$\rightarrow\f$ all matter within
       *  the radius \f$r_{\Delta}\f$ for which the mean internal
       *  density is \f$\Delta\Omega_m\f$ times the critical mean
       *  background density
       *
       *  @param trunc_fact factor \f$F_t\f$ defining the truncation
       *  radius, that is \f$r_t = F_tr_{\Delta}\f$. The default value, 3,
       *  corresponds to a reasonable choice in the case of a truncated
       *  NFW profile
       *
       *  @param miscentering if true, account for the miscentering
       *
       *  @param single_profile if true, the miscetering is related to
       *  a single cluster. Otherwise, a population of clusters is
       *  considered (the latter is the case of stacked weak lensing
       *  analyses)
       *
       *  @param sigma_off if single_profile=false, this is the
       *  standard deviation of the miscentered cluster population,
       *  \f$\sigma_{\rm off}\f$, in Mpc\f$/h\f$.  Otherwise, it is
       *  the miscentering of a single cluster profile
       *
       *  @param f_off fraction of miscentered clusters, \f$f_{\rm off}\f$.
       *  Used only if single_profile=false
       *
       *  @param b_e normalisation of the DK14 outer profile
       *
       *  @param s_e slope of the DK14 outer profile
       *
       *  @param AB_fact anisotropic boost factor, \f$\sigma_{\rm
       *  AB}\f$, affecting the 2-halo term and/or the infall term in
       *  DK14. Specifically, such terms are multiplied by
       *  \f$(1+\sigma_{\rm AB})\f$
       *
       *  @param bias author(s) who proposed the bias; valid authors are: 
       *  ST99 (Sheth & Tormen 1999), SMT01 (Sheth, Mo & Tormen 2001), 
       *  SMT01_WL04 (Sheth, Mo & Tormen 2001 with the correction of Warren 2004), 
       *  Tinker (Tinker et al. 2010)
       *
       *  @param method_Pk method used for the computation of the
       *  power spectrum.  Valid choices for method_Pk are: CAMB
       *  [http://camb.info/], CLASS [http://class-code.net/],
       *  MPTbreeze-v1 [http://arxiv.org/abs/1207.1465], EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       *
       *  @param alpha_0 \f$\alpha_0\f$ parameter in the equation
       *  \f$ \alpha(\nu) = \alpha_0 + \alpha_\nu\nu^2 \f$, where \f$\nu\f$ 
       *  is the virial peak height. \f$\alpha(\nu)\f$ appears in both
       *  Einasto and DK14 profiles. Following Gao et al. (2008), 
       *  \f$\alpha_0=0.155\f$ by default.
       *
       *  @param alpha_nu \f$\alpha_\nu\f$ parameter in the equation
       *  \f$ \alpha(\nu) = \alpha_0 + \alpha_\nu\nu^2 \f$, where \f$\nu\f$ 
       *  is the virial peak height. \f$\alpha(\nu)\f$ appears in both
       *  Einasto and DK14 profiles. Following Gao et al. (2008), 
       *  \f$\alpha_\nu=0.0095\f$ by default.
       *
       *  @param beta \f$\beta\f$ parameter in the DK14 profile. Following
       *  DK14, \f$\beta=4\f$ by default.
       *
       *  @param gamma_0 \f$\gamma_0\f$ parameter in the expression
       *  \f$\gamma=\gamma_0\nu\f$ appearing in the DK14 profile. Following
       *  DK14, \f$\gamma_0=4\f$ by default.
       *
       *  @warning If the miscentering is considered, the methods
       *  HaloProfile::Sigma() and HaloProfile::DeltaSigma() include
       *  the miscentering contribution (i.e. Eq. 1 and Eq. 2). The
       *  same is valid for the 2-halo term.
       *
       *  @warning If in the cbl::cosmology::Cosmology object the
       *  \f$h\f$ units are set, both the input quantities (e.g. mass,
       *  radius) and the output quantities (e.g. \f$\Sigma\f$,
       *  \f$\Delta\Sigma\f$) are in units of \f$h\f$.
       *
       *  @warning Both in the Einasto and DK14 profiles, the Einasto
       *  \f$\alpha\f$ parameter is expressed as \f$ \alpha(\nu) =
       *  \alpha_0 + \alpha_\nu\nu^2 \f$, where \f$\nu\f$ is the virial peak
       *  height.
       *
       *  @warning By default, the DK14 implementation has \f$\beta=4\f$,
       *  while \f$\gamma\f$ varies with the virial peak height \f$\nu\f$,
       *  namely \f$\gamma=4\nu\f$, and \f$F_t\f$ is a constant. In
       *  this case, DK14 found \f$F_t=1.495\f$.
       *
       *  @warning If the DK14 profile is used, the outer term is
       *  always included and it replaces the 2-halo
       *  term. Specifically, the outer term is expressed as
       *
       *  \f$ \rho_{\rm out} = \rho_{\rm m} \left[ b_{\rm e} \left(
       *  \frac{r}{5R_\Delta} \right)^{-s_{\rm e}} + 1 \right] \f$
       *
       */
      HaloProfile (std::shared_ptr<cbl::cosmology::Cosmology> cosmology, const bool two_halo, const double redshift, const double conc, const double Mass, const double Delta, const std::string profile_author, const std::string halo_def, const double trunc_fact=3., const bool miscentering=false, const bool single_profile=false, const double sigma_off=0.1, const double f_off=1., const double b_e=1., const double s_e=1.5, const double AB_fact=0., const std::string bias="Tinker", const std::string method_Pk="EisensteinHu", const double alpha_0=0.155, const double alpha_nu=0.0095, const double beta=4., const double gamma_0=4.);
      
      /**
       *  @brief constructor that sets the halo profile parameters,
       *  assuming a concentration-mass relation. 
       *  If the miscentering is considered, \f$\Sigma(R)\f$ and
       *  \f$\Delta\Sigma(R)\f$ are expressed as (see e.g. Bellagamba
       *  et al. 2019):
       *
       *  \f$\Sigma(R)=(1-f_{\rm off})\Sigma_{\rm cen}(R)+
       *  f_{\rm off}\Sigma_{\rm off}(R),\,\,\,\,(1)\f$
       *
       *  \f$\Delta\Sigma(R)=(1-f_{\rm off})\Delta\Sigma_{\rm
       *  cen}(R)+ f_{\rm off}\Delta\Sigma_{\rm off}(R).\,\,\,\,(2)\f$
       *
       *  where \f$\Sigma_{\rm cen}\f$ and \f$\Sigma_{\rm cen}\f$ are the
       *  centred and miscentred surface densities, respectively. The
       *  miscentring probability follows a Rayleigh distribution.
       * 
       *  The miscentring also affects the 2-halo term. Indeed, the
       *  general expression of \f$\Sigma_{\rm cen}\f$ is:
       *
       *  \f$\Sigma_{\rm cen} = \Sigma_{\rm 1h} + \Sigma_{\rm 2h},\f$
       *
       *  where \f$\Sigma_{\rm 1h}\f$ and \f$\Sigma_{\rm 2h}\f$ are the
       *  1-halo and 2-halo surface densities, respectively.
       *  \f$\Sigma_{\rm 2h}\f$ is expressed as follows (Oguri & Takada 2011):
       *
       *  \f$\Sigma_{\rm 2h}(R) = \int\,\frac{l{\rm d}l}{2\pi}
       *  J_0(l\theta)\frac{\bar{\rho}_{\rm m}(z)b(M,z)}{(1+z)^3D_{\rm l}^2(z)}
       *  P(k_l,z) \f$
       *
       *  @param cosmology the cosmology
       *
       *  @param two_halo if true, include the 2-halo contribution
       *
       *  @param redshift the redshift
       *
       *  @param cM_author author(s) who proposed the 
       *  concentration-mass relation. Possibilities are:
       *  "Duffy" (Duffy et al. 2008)
       *
       *  @param Mass halo mass
       *
       *  @param Delta overdensity factor which needs to be multiplied
       *  to the critical density in order to define an overdensity
       *
       *  @param profile_author the density profile author(s);
       *  available options are: "NFW" \f$\rightarrow\f$
       *  Navarro-Frenk-White profile (Navarro et al. 1997);
       *  "NFW_trunc" \f$\rightarrow\f$ truncated Navarro-Frenk-White
       *  profile, aka BMO (Baltz et al. 2009); "Einasto"
       *  \f$\rightarrow\f$ Einasto 1965; "DK14" \f$\rightarrow\f$
       *  Diemer & Kravtsov 2014
       *
       *  @param halo_def the halo definition; available options are:
       *  "vir" \f$\rightarrow\f$ all matter within the radius
       *  \f$r_{vir}\f$ for which the mean internal density is
       *  cbl::cosmology::Cosmology::Delta_c times the critical density
       *  \f$\rho_{crit}=3H^2/8\pi G\f$; "critical" \f$\rightarrow\f$ all
       *  matter within the radius \f$r_{\Delta}\f$ for which the mean
       *  internal density is \f$\Delta\f$ times the critical density; "mean"
       *  \f$\rightarrow\f$ all matter within the radius
       *  \f$r_{\Delta}\f$ for which the mean internal density is 
       *  \f$\Delta\Omega_m\f$
       *  times the critical mean background density
       *
       *  @param trunc_fact factor \f$F_t\f$ defining the truncation
       *  radius, that is \f$r_t = F_tr_{\Delta}\f$. The default value, 3,
       *  corresponds to a reasonable choice in the case of a truncated
       *  NFW profile
       *
       *  @param miscentering if true, account for the miscentering
       *
       *  @param single_profile if true, the miscetering is related to
       *  a single cluster. Otherwise, a population of clusters is
       *  considered (the latter is the case of stacked weak lensing
       *  analyses)
       *
       *  @param sigma_off if single_profile=false, this is the
       *  standard deviation of the miscentered cluster population,
       *  \f$\sigma_{\rm off}\f$, in Mpc\f$/h\f$.  Otherwise, it is
       *  the miscentering of a single cluster profile
       *
       *  @param f_off fraction of miscentered clusters, \f$f_{\rm off}\f$.
       *  Used only if single_profile=false
       *
       *  @param b_e normalisation of the DK14 outer profile
       *
       *  @param s_e slope of the DK14 outer profile
       *
       *  @param AB_fact anisotropic boost factor, \f$\sigma_{\rm
       *  AB}\f$, affecting the 2-halo term and/or the infall term in
       *  DK14. Specifically, such terms are multiplied by
       *  \f$(1+\sigma_{\rm AB})\f$
       *
       *  @param bias author(s) who proposed the bias; valid authors
       *  are: ST99 (Sheth & Tormen 1999), SMT01 (Sheth, Mo & Tormen
       *  2001), SMT01_WL04 (Sheth, Mo & Tormen 2001 with the
       *  correction of Warren 2004), Tinker (Tinker et al. 2010)
       *
       *  @param method_Pk method used for the computation of the
       *  power spectrum.  Valid choices for method_Pk are: CAMB
       *  [http://camb.info/], CLASS [http://class-code.net/],
       *  MPTbreeze-v1 [http://arxiv.org/abs/1207.1465], EisensteinHu
       *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
       *
       *  @param alpha_0 \f$\alpha_0\f$ parameter in the equation
       *  \f$ \alpha(\nu) = \alpha_0 + \alpha_\nu\nu^2 \f$, where \f$\nu\f$ 
       *  is the virial peak height. \f$\alpha(\nu)\f$ appears in both
       *  Einasto and DK14 profiles. Following Gao et al. (2008), 
       *  \f$\alpha_0=0.155\f$ by default.
       *
       *  @param alpha_nu \f$\alpha_\nu\f$ parameter in the equation
       *  \f$ \alpha(\nu) = \alpha_0 + \alpha_\nu\nu^2 \f$, where \f$\nu\f$ 
       *  is the virial peak height. \f$\alpha(\nu)\f$ appears in both
       *  Einasto and DK14 profiles. Following Gao et al. (2008), 
       *  \f$\alpha_\nu=0.0095\f$ by default.
       *
       *  @param beta \f$\beta\f$ parameter in the DK14 profile. Following
       *  DK14, \f$\beta=4\f$ by default.
       *
       *  @param gamma_0 \f$\gamma_0\f$ parameter in the expression
       *  \f$\gamma=\gamma_0\nu\f$ appearing in the DK14 profile. Following
       *  DK14, \f$\gamma_0=4\f$ by default.
       *
       *  @warning the Duffy et al. concentrantion-mass relation
       *  refers to the 0<z<2 redshift range, obtained from their full
       *  samples (see Table 1 of Duffy et al. 2008).
       *
       *  @warning If the miscentering is considered, the methods
       *  HaloProfile::Sigma() and HaloProfile::DeltaSigma() include
       *  the miscentering contribution (i.e. Eq. 1 and Eq. 2). The
       *  same is valid for the 2-halo term.
       *
       *  @warning If in the cbl::cosmology::Cosmology object the
       *  \f$h\f$ units are set, both the input quantities (e.g. mass,
       *  radius) and the output quantities (e.g. \f$\Sigma\f$,
       *  \f$\Delta\Sigma\f$) are in units of \f$h\f$.
       *
       *  @warning Both in the Einasto and DK14 profiles, the Einasto
       *  \f$\alpha\f$ parameter is expressed as \f$ \alpha(\nu) =
       *  \alpha_0 + \alpha_\nu\nu^2 \f$, where \f$\nu\f$ is the virial peak
       *  height.
       *
       *  @warning By default, the DK14 implementation has \f$\beta=4\f$,
       *  while \f$\gamma\f$ varies with the virial peak height \f$\nu\f$,
       *  namely \f$\gamma=4\nu\f$, and \f$F_t\f$ is a constant. In
       *  this case, DK14 found \f$F_t=1.495\f$.
       *
       *  @warning If the DK14 profile is used, the outer term is
       *  always included and it replaces the 2-halo
       *  term. Specifically, the outer term is expressed as
       *
       *  \f$ \rho_{\rm out} = \rho_{\rm m} \left[ b_{\rm e} \left(
       *  \frac{r}{5R_\Delta} \right)^{-s_{\rm e}} + 1 \right] \f$
       *
       */
      HaloProfile (std::shared_ptr<cbl::cosmology::Cosmology> cosmology, const bool two_halo, const double redshift, const std::string cM_author, const double Mass, const double Delta, const std::string profile_author, const std::string halo_def, const double trunc_fact=3., const bool miscentering=false, const bool single_profile=false, const double sigma_off=0.1, const double f_off=1., const double b_e=1., const double s_e=1.5, const double AB_fact=0., const std::string bias="Tinker", const std::string method_Pk="EisensteinHu", const double alpha_0=0.155, const double alpha_nu=0.0095, const double beta=4., const double gamma_0=4.);
      
      /**
       *  @brief default destructor
       */
      ~HaloProfile () = default;

      ///@}
  
      /**
       *  @name Member functions used to set the private members 
       */
      ///@{
      
      /**
       *  @brief update the 2-halo term. It must be
       *  used if e.g. cosmological parameters or halo mass changes.
       *
       *  @param update_bias if true, update also the halo bias value
       */
      void update_2halo (const bool update_bias=true);
      
      /**
       *  @brief set the cosmological model
       *
       *  @param cosmology the cosmology
       *
       *  @param update_2halo if true, update the 2-halo term and the
       *  halo bias (only if the 2-halo is requested in the
       *  constructor)
       *
       */
      void set_cosmology (const std::shared_ptr<cbl::cosmology::Cosmology> cosmology, const bool update_2halo=true);
      
      /**
       *  @brief set the private member \e m_mass
       *  @param mass the mass of the cluster
       */
      void set_mass (const double mass)
      { m_mass = mass; }
      
      /**
       *  @brief set the private member \e m_concentration
       *  @param conc concentration
       */
      void set_concentration (const double conc)
      { m_concentration = conc; }
      
      /**
       *  @brief set the private member \e m_f_off
       *  @param f_off f_off
       */
      void set_f_off (const double f_off)
      { m_f_off = f_off; }
      
      /**
       *  @brief set the private member \e m_sigma_off
       *  @param sigma_off sigma_off
       */
      void set_sigma_off (const double sigma_off)
      { m_sigma_off = sigma_off; }
      
      /**
       *  @brief set the private member \e m_trunc_fact
       *  @param trunc_fact trunc_fact
       */
      void set_trunc_fact (const double trunc_fact)
      { m_trunc_fact = trunc_fact; }
      
      /**
       *  @brief set the private member \e m_AB_fact
       *  @param AB_fact AB_fact
       */
      void set_AB_fact (const double AB_fact)
      { m_AB_fact = AB_fact; }
      
      /**
       *  @brief set the private member \e m_b_e
       *  @param b_e b_e
       */
      void set_b_e (const double b_e)
      { m_b_e = b_e; }
      
      /**
       *  @brief set the private member \e m_s_e
       *  @param s_e s_e
       */
      void set_s_e (const double s_e)
      { m_s_e = s_e; }
      
      /**
       *  @brief set the private member \e m_bias
       *  @param bias the halo bias
       */
      void set_bias (const double bias)
      { m_bias = bias; }
      
      /**
       *  @brief set the \f$ \alpha(\nu) = \alpha_0 + \alpha_\nu\nu^2 \f$
       *  parameter of the Einasto and DK14 profiles
       *  @param alpha_0 \f$\alpha_0\f$
       *  @param alpha_nu \f$\alpha_\nu\f$
       */
      void set_alpha (const double alpha_0, const double alpha_nu)
      { m_alpha_0 = alpha_0; m_alpha_nu = alpha_nu; }
      
      /**
       *  @brief set the \f$\beta\f$
       *  parameter of the DK14 profile
       *  @param beta \f$\beta\f$
       */
      void set_beta (const double beta)
      { m_beta = beta; }
      
      /**
       *  @brief set \f$\gamma_0\f$ appearing in the
       *  \f$\gamma=\gamma_0\nu\f$ parameter of the DK14 profile
       *  @param gamma_0 \f$\gamma_0\f$
       */
      void set_gamma_0 (const double gamma_0)
      { m_gamma_0 = gamma_0; }
      
      ///@}

      
      /**
       *  @name Member functions used to get the private members or
       *  to compute the cluster halo profiles
       */
      ///@{
      
      /**
       *  @brief Get the private member m_cosmology
       *
       *  @return the cosmological model
       */
      std::shared_ptr<cosmology::Cosmology> cosmology ()
      { return move(m_cosmology); }
      
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
      { return m_Delta_func(m_Delta, m_cosmology, m_redshift); }
      
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
      std::vector<double> rho (const std::vector<double> rad);
      
      /**
       *  @brief The total surface density profile of the halo,
       *  including miscentering and 2-halo
       *
       *  @param rad radius
       *
       *  @return the halo surface density profile in \f$M_{\odot}\f$
       *  pc^2
       *
       *  @warning If in the cbl::cosmology::Cosmology object the
       *  \f$h\f$ units are set, both the input quantities (e.g. mass,
       *  radius) and the output quantities are in units of \f$h\f$
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
       *  @return the halo surface density profile in \f$M_{\odot} / \f$pc^2
       *
       *  @warning If in the cbl::cosmology::Cosmology object the \f$h\f$ units are
       *  set, both the input quantities (e.g. mass, radius) and the output quantities
       *  are in units of \f$h\f$
       */
      std::vector<double> Sigma_2h (const std::vector<double> rad);
      
      /**
       *  @brief the total excess surface density profile of the halo, including
       *  miscentering and 2-halo. 
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
       *  @return the halo excess surface density profile in \f$M_{\odot} / \f$pc^2
       *
       *  @warning If in the cbl::cosmology::Cosmology object the \f$h\f$ units are
       *  set, both the input quantities (e.g. mass, radius) and the output quantities
       *  are in units of \f$h\f$
       */
      std::vector<double> DeltaSigma_2h (const std::vector<double> rad);
      
      /**
       *  @brief the total reduced shear of the halo, including
       *  miscentering and 2-halo. The integration over \f$N(z_{\rm s})\f$
       *  is performed for \f$z_{\rm s}>z_{\rm l}\f$, where \f$z_{\rm s}\f$
       *  is the source redshift and \f$z_{\rm l}\f$ is the lens redshift. 
       *  In particular, the reduced shear is expressed as follows:
       *
       *  \f$ g = \frac{\gamma}{1-\kappa}\,, \f$
       *
       *  where (Seitz and Schneider 1997)
       *
       *  \f$\gamma = \Delta\Sigma\langle\Sigma_{\rm crit}^{-1}\rangle\,, \f$
       *
       *  \f$\kappa = \Sigma\langle\Sigma_{\rm crit}^{-1}\rangle^{-1}\langle\Sigma_{\rm crit}^{-2}\rangle\,, \f$
       *
       *  \f$ \langle\Sigma_{\rm crit}^{-n}\rangle = 
       *  \frac{\sum_i \Sigma_{\rm crit}^{-n}(z_{\rm l},z_{\rm s}^i) N(z_{\rm s}^i)}
       *  {\sum_i N(z_{\rm s}^i)}\,. \f$
       *
       *  @param rad radius
       *
       *  @param z_s central bin values where the binned \f$N(z_{\rm s})\f$ distribution
       *  is computed, where 
       *  \f$z_{\rm s}\f$ is the background source redshift.
       *
       *  @param N_zs the binned \f$N(z_{\rm s})\f$ distribution
       *
       *  @return the reduced shear
       *
       *  @warning If in the cbl::cosmology::Cosmology object the \f$h\f$ units are
       *  set, both the input quantities (e.g. mass, radius) and the output quantities
       *  are in units of \f$h\f$
       *
       */
      std::vector<double> g (const std::vector<double> rad, const std::vector<double> z_s, const std::vector<double> N_zs);
      
      
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
       *  computed by cbl::cosmology::HaloProfile::concentration;
       *  \f$\Delta_c(z)\f$ is computed by
       *  cbl::cosmology::Cosmology::Delta_c and \f$\rho_{crit}(z)\f$
       *  is computed by cbl::cosmology::Cosmology::rho_crit;
       *  \f$r_{vir}(M_h, z)\f$ is computed by
       *  cbl::cosmology::HaloProfile::r_vir
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

       /**
       *  @brief the halo concentration-mass relation for NFW profile
       *  and Duffy model
       * 
       *  This member function computes the concentration of a dark
       *  matter halo of a given a mass, at a given redshift following
       *  Duffy et al. 2008 model: \f[c(M_h, z) =
       *  A(M_h/M_{pivot})^B\,(1+z)^C\f]
       *
       *  @param Mass the halo mass
       *
       *  @param redshift the redshift (must be < 2)
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
       *  @return the halo concentration for NFW profile and Duffy
       *  model
       *
       *  @warning the Duffy et al. concentrantion-mass relation
       *  refers to the 0<z<2 redshift range, obtained from their full
       *  samples (see Table 1 of Duffy et al. 2008); actually, the
       *  current implementation does not depend on cosmology
       */
      double concentration_NFW_Duffy (const double Mass, const double redshift, const std::string halo_def="vir") const;

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
       *  @brief The virial radius, given the virial mass and the
       *  redshift
       *
       *  This member function computes the virial halo radius as
       *  follows:
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
      
      ///@}
    
    };
  }
}

#endif

