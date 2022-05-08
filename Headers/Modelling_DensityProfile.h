/********************************************************************
 *  Copyright (C) 2021 by Federico Marulli and Giorgio Lesci        *
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
 *  @file Headers/Modelling_DensityProfile.h
 *
 *  @brief The class Modelling_DensityProfile
 *
 *  This file defines the interface of the class Modelling_DensityProfile, used to
 *  model the density profile of galaxy clusters
 *
 *  @authors Giorgio Lesci (and Federico Marulli)
 *
 *  @authors giorgio.lesci2@unibo.it (and federico.marulli3@unibo.it)
 */

#ifndef __MODELLINGDPROFILE__
#define __MODELLINGDPROFILE__

#include "HaloProfile.h"
#include "StackedDensityProfile.h"
#include "Modelling_MassObservableRelation.h"


// ===================================================================================================


namespace cbl {

  namespace modelling {

    /**
     *  @brief The namespace of the <B> cluster projected density profile
     *  modelling </B>
     *  
     *  The \e modelling::densityprofile namespace contains all the functions
     *  and classes to model the surface density excess profile of galaxy clusters
     */
    namespace densityprofile {
      
      /**
       *  @struct STR_Profile_data_model
       *  @brief the structure STR_Profile_data_model
       *
       *  This structure contains the data used for statistical
       *  analyses of cluster surface density excess profiles
       */
      struct STR_Profile_data_model {
      
        /// function returning the concentration 
	std::function<double(const double, cbl::cosmology::HaloProfile)> conc_func;
	
	/// function returning the 2-halo term 
	std::function<std::vector<double>(const std::vector<double>, cbl::cosmology::HaloProfile, const std::string, const std::string, const std::string)> two_halo_func;
	
	/// function returning the 2-halo term, taking the halo bias in input
	std::function<std::vector<double>(const std::vector<double>, cbl::cosmology::HaloProfile, const double, const std::string, const std::string)> two_halo_func_fast;
      
	/// Fiducial cosmology pointer
	std::shared_ptr<cosmology::Cosmology> cosmology;
	
	/// Cosmological parameters
	std::vector<cosmology::CosmologicalParameter> Cpar;
	
	/// Modelling_MassObservableRelation object pointer
	std::shared_ptr<modelling::massobsrel::Modelling_MassObservableRelation> scaling_relation;
	
	/// HaloProfile object pointer
	std::shared_ptr<cosmology::HaloProfile> halo_profile;
	
	/// Redshift
	double redshift;
	
	/// Mass proxy
	double mass_proxy;
	
	/// Redshift pivot in the scaling relation
	double redshift_pivot;
	
	/// Mass proxy pivot in the scaling relation
	double proxy_pivot;
	
	/// The density contrast with respect to the critical density (i.e. if equal to 200, \f$M_{200}\f$ is considered)
	double contrast;	
	
	/// Base of the mass logarithm
	double logM_base;
	
	/// The mass pivot
	double mass_pivot;
	
	/// author(s) who proposed the bias function
	std::string bias_author;
	
	/// method used for the computation of the power spectrum
	std::string method_Pk;
	
	/// interpolation type for the power spectrum
	std::string interp_type;

	/**
	 *  @brief default constructor
	 */
	STR_Profile_data_model () = default;
      };
    
      /**
       *  @class Modelling_DensityProfile
       *  Modelling_DensityProfile.h
       *  "Headers/Modelling_DensityProfile.h"
       *
       *  @brief The class Modelling_DensityProfile
       *
       *  This file defines the interface of the base class
       *  Modelling_DensityProfile, used for modelling
       *  cluster surface density profile excess measurements, 
       *  i.e. \f$\Delta\Sigma(r)\f$ [\f$h\f$ M\f$_\odot\f$/pc\f$^2\f$].
       *  Cosmological units are forced.
       *
       *  In particular, the excess surface mass density 
       *  is expressed as  (Sheldon et al. 2004)
       *  \f[ \Delta\Sigma (R) \equiv \overline{\Sigma}(< R) - \Sigma (R), \f]
       *
       *  where \f$\Sigma(R)\f$ is the surface mass density, while
       *  \f$\overline{\Sigma}(<R)\f$ is its mean value within
       *  the projected radius R. 
       *
       */
      class Modelling_DensityProfile : public Modelling {
      
      protected:

	/// the container of parameters for the density model computation
	STR_Profile_data_model m_data_model;
	
	/// if true, consider the 2-halo contribution
	bool m_2halo;
	
	/// if true, the mass is a parameter derived from the scaling relation
	bool m_mass_is_derived;
	
	/// author(s) of the cluster density profile
        std::string m_profile_author;
        
        /// halo definition
        std::string m_halo_def;
        
        /// overdensity factor
        double m_Delta;


      public:

	/**
	 *  @name Constructors/destructors
	 */
	///@{

	/**
	 *  @brief default constuctor
	 *  _DensityProfile
	 */
	Modelling_DensityProfile () = default;
	
	/**
	 *  @brief constuctor for the modelling of
	 *  stacked density excess profile of galaxy clusters. Cosmological units are forced
	 *
	 *  @param profile the object of stacked density profile to model
	 *
	 *  @param profile_author author(s) of the cluster density profile.
	 *  See available options in cbl::cosmology::HaloProfile
	 *
	 *  @param _2halo if true, include the 2-halo contribution
	 *
	 *  @param halo_def the halo definition; available options are:
	 *  "critical", "vir", "mean"
	 *
	 *  @param Delta overdensity factor which needs to be multiplied 
	 *  to the critical density in order to define an overdensity
	 *
	 */
	Modelling_DensityProfile (const std::shared_ptr<cbl::measure::stackprofile::StackedDensityProfile> profile, const std::string profile_author="NFW_trunc", const bool _2halo=false, const std::string halo_def="critical", const double Delta=200.);
	
	/**
	 *  @brief constuctor for the modelling of
	 *  stacked density excess profile of galaxy cluster. Cosmological units are forced.
	 *
	 *  @param dataset cluster profile dataset
	 *
	 *  @param profile_author author(s) of the cluster density profile.
	 *  See available options in cbl::cosmology::HaloProfile
	 *
	 *  @param _2halo if true, include the 2-halo contribution
	 *
	 *  @param halo_def the halo definition; available options are:
	 *  "critical", "vir", "mean"
	 *
	 *  @param Delta overdensity factor which needs to be multiplied 
	 *  to the critical density in order to define an overdensity
	 *
	 */
	Modelling_DensityProfile (const std::shared_ptr<cbl::data::Data> dataset, const std::string profile_author="NFW_trunc", const bool _2halo=false, const std::string halo_def="critical", const double Delta=200.);
	
	/**
	 *  @brief default destructor
	 *  
	 */
	virtual ~Modelling_DensityProfile () = default;

      ///@}

	/**
	 *  @name Member functions used to set the protected members of the class
	 */
	///@{

	/**
	 * @brief get the member \e m_data_model
	 * @return the container of parameters for cluster density profile
	 * model computation
	 */
	STR_Profile_data_model data_model () { return m_data_model; }
	
	/**
	 *  @name Member functions used to set the model parameters
	 */
	///@{

	/**
	 *  @brief Set the data used to construct generic models of
	 *  halo profiles.
	 *  
	 *  @param cosmology the cosmological model
	 *
	 *  @param redshift cluster redshift
	 *
	 *  @param contrast density contrast with respect 
	 *  to the critical density (i.e. if equal to 200, \f$M_{200}\f$ is considered)
	 *
	 *  @param logM_base base of the mass logarithm
	 *
	 *  @param mass_pivot the mass pivot
	 *
	 *  @param bias_author author(s) who proposed the bias; valid authors are: 
	 *  ST99 (Sheth & Tormen 1999), SMT01 (Sheth, Mo & Tormen 2001), 
	 *  SMT01_WL04 (Sheth, Mo & Tormen 2001 with the correction of Warren 2004), 
	 *  Tinker (Tinker et al. 2010)
	 *
	 *  @param method_Pk method used for the computation of the power spectrum.
	 *
	 *  @param interp_type method to interpolate the power spectrum.
	 *  
	 */
	void set_data_model (const cbl::cosmology::Cosmology cosmology, const double redshift, const double contrast, const double logM_base, const double mass_pivot, const std::string bias_author="Tinker", const std::string method_Pk="EisensteinHu", std::string interp_type="Linear");
	
	/**
	 *  @brief Set the data used to construct generic models of
	 *  halo profiles, where the mass is derived from a scaling
	 *  relation.
	 *  
	 *  @param cosmology the cosmological model
	 *
	 *  @param redshift cluster redshift
	 *
	 *  @param mass_proxy cluster mass proxy
	 *
	 *  @param redshift_pivot redshift pivot in the scaling relation
	 *
	 *  @param proxy_pivot proxy pivot in the scaling relation
	 *
	 *  @param contrast density contrast with respect 
	 *  to the critical density (i.e. if equal to 200, \f$M_{200}\f$ is considered)
	 *
	 *  @param logM_base base of the mass logarithm
	 *
	 *  @param mass_pivot the mass pivot
	 *
	 *  @param Nclusters number of clusters in the bin used for the stacking
	 *
	 *  @param bias_author author(s) who proposed the bias; valid authors are: 
	 *  ST99 (Sheth & Tormen 1999), SMT01 (Sheth, Mo & Tormen 2001), 
	 *  SMT01_WL04 (Sheth, Mo & Tormen 2001 with the correction of Warren 2004), 
	 *  Tinker (Tinker et al. 2010)
	 *
	 *  @param method_Pk method used for the computation of the power spectrum.
	 *
	 *  @param interp_type method to interpolate the power spectrum.
	 *  
	 */
	void set_data_model (const cosmology::Cosmology cosmology, const double redshift, const double mass_proxy, const double redshift_pivot, const double proxy_pivot, const double contrast, const double logM_base, const double mass_pivot, const double Nclusters, const std::string bias_author="Tinker", const std::string method_Pk="EisensteinHu", std::string interp_type="Linear");

	/**
	 *  @brief Set the profile and cosmological parameters used to model the 
	 *  cluster density profile. 
	 *
	 *  The 1-halo term has the following functional form, including the
	 *  contribution of centered and off-centered populations
	 *  of galaxy clusters (e.g. Bellagamba et al. 2019, Johnston et al. 2007, Yang et al. 2006),
	 *
	 *  \f$\Delta\Sigma_{\rm 1h}(R)=
	 *  (1-f_{\rm off})\Delta\Sigma_{\rm cen}(R)+f_{\rm off}\Delta\Sigma_{\rm off}(R),\f$
	 *
	 *  while the 2-halo term is expressed as (e.g. Bellagamba et al. 2019)
	 *
	 *  \f$\Delta\Sigma_{\rm 2h}(R) = \int\,\frac{l{\rm d}l}{2\pi}
	 *  J_2(l\theta)\frac{\bar{\rho}_{\rm m}(z)b(M,z)}{(1+z)^3D_{\rm l}^2(z)}
	 *  P(k_l,z). \f$
	 *
	 *  The mass is
	 *  expressed e.g. in \f$10^{14}\f$ M\f$_\odot\f$ \f$h^{-1}\f$, if
	 *  mass_pivot set through set_data_model is \f$10^{14}\f$.
	 *
	 *  @param cosmo_param vector of enums containing cosmological
	 *  parameters
	 *
	 *  @param cosmo_prior vector containing the priors for
	 *  the cosmological parameters
	 *
	 *  @param Rt_prior prior on the NFW truncation factor \f$F_t\f$ defining the truncation
         *  radius, that is \f$r_t = F_tr_{\Delta}\f$ 
	 *
	 *  @param concentration_prior prior on the cluster concentration
	 *
	 *  @param logM_prior prior on the mass logarithm (where the mass is
	 *  expressed e.g. in \f$10^{14}\f$ M\f$_\odot\f$ \f$h^{-1}\f$, if
	 *  mass_pivot set through set_data_model is \f$10^{14}\f$). The base of
	 *  the logarithm is set through set_data_model
	 *
	 *  @param f_off_prior prior on the fraction of miscentered clusters. 
	 *  This parameter makes sense only if the user models a stacked profile, 
	 *  not a single cluster profile. If a single profile is modelled, set a constant prior equal
	 *  to 0 or 1 for f_off.
	 *
	 *  @param sigma_off_prior prior on the rms of the miscentered cluster population.
	 *
	 *  @param anisotropic_boost_prior prior on the anisotropic boost factor,
	 *  internally called "AB_fact", **entering the 
	 *  model only if the 2-halo term is considered**. In particular, the 2-halo excess surface
	 *  density is expressed as \f$\Delta\Sigma_{\rm 2h,\,correct}=
	 *  \Delta\Sigma_{\rm 2h}(1+\sigma_{\rm AB})\f$, where \f$\sigma_{\rm AB}\f$ is
	 *  the parameter set through this prior.
	 *
	 *  @param orientation_boost_prior prior on the orientation boost factor,
	 *  internally called "OB_fact". In particular, the profile mass is expressed
	 *  as \f$M_{\rm correct}=M(1+\sigma_{\rm OB})\f$, where \f$\sigma_{\rm OB}\f$ is
	 *  the parameter set through this prior.
	 *
	 *  @warning \f$F_t\f$ is used only if a truncated NFW is assumed!
	 *
	 *  @warning The off-centering is related to stacks of clusters. For
	 *  details, see cbl::cosmology::HaloProfile
	 *
	 */
	void set_model_DensityProfile_cosmology (const std::vector<cbl::cosmology::CosmologicalParameter> cosmo_param, const std::vector<statistics::PriorDistribution> cosmo_prior, const statistics::PriorDistribution Rt_prior, const statistics::PriorDistribution concentration_prior, const statistics::PriorDistribution logM_prior, const statistics::PriorDistribution f_off_prior, const statistics::PriorDistribution sigma_off_prior, const statistics::PriorDistribution anisotropic_boost_prior, const statistics::PriorDistribution orientation_boost_prior);
	
	/**
	 *  @brief Set the profile and cosmological parameters used to model the 
	 *  cluster density profile. The concentration is a derived parameter,
	 *  computed through a concentration-mass relation.
	 *
	 *  The 1-halo term has the following functional form, including the
	 *  contribution of centered and off-centered populations
	 *  of galaxy clusters (e.g. Bellagamba et al. 2019, Johnston et al. 2007, Yang et al. 2006),
	 *
	 *  \f$\Delta\Sigma_{\rm 1h}(R)=
	 *  (1-f_{\rm off})\Delta\Sigma_{\rm cen}(R)+f_{\rm off}\Delta\Sigma_{\rm off}(R),\f$
	 *
	 *  while the 2-halo term is expressed as (e.g. Bellagamba et al. 2019)
	 *
	 *  \f$\Delta\Sigma_{\rm 2h}(R) = \int\,\frac{l{\rm d}l}{2\pi}
	 *  J_2(l\theta)\frac{\bar{\rho}_{\rm m}(z)b(M,z)}{(1+z)^3D_{\rm l}^2(z)}
	 *  P(k_l,z). \f$
	 *
	 *  The mass is
	 *  expressed e.g. in \f$10^{14}\f$ M\f$_\odot\f$ \f$h^{-1}\f$, if
	 *  mass_pivot set through set_data_model is \f$10^{14}\f$.
	 *
	 *  @param cosmo_param vector of enums containing cosmological
	 *  parameters
	 *
	 *  @param cosmo_prior vector containing the priors for
	 *  the cosmological parameters
	 *
	 *  @param Rt_prior prior on the NFW truncation factor \f$F_t\f$ defining the truncation
         *  radius, that is \f$r_t = F_tr_{\Delta}\f$ 
	 *
	 *  @param cM_author author(s) who proposed the 
	 *  concentration-mass relation. For details, see
	 *  cbl::cosmology::HaloProfile
	 *
	 *  @param logM_prior prior on the mass logarithm (where the mass is
	 *  expressed e.g. in \f$10^{14}\f$ M\f$_\odot\f$ \f$h^{-1}\f$, if
	 *  mass_pivot set through set_data_model is \f$10^{14}\f$). The base of
	 *  the logarithm is set through set_data_model
	 *
	 *  @param f_off_prior prior on the fraction of miscentered clusters. 
	 *  This parameter makes sense only if the user models a stacked profile, 
	 *  not a single cluster profile. If a single profile is modelled, set a constant prior equal
	 *  to 0 or 1 for f_off.
	 *
	 *  @param sigma_off_prior prior on the rms of the miscentered cluster population.
	 *
	 *  @param anisotropic_boost_prior prior on the anisotropic boost factor,
	 *  internally called "AB_fact", **entering the 
	 *  model only if the 2-halo term is considered**. In particular, the 2-halo excess surface
	 *  density is expressed as \f$\Delta\Sigma_{\rm 2h,\,correct}=
	 *  \Delta\Sigma_{\rm 2h}(1+\sigma_{\rm AB})\f$, where \f$\sigma_{\rm AB}\f$ is
	 *  the parameter set through this prior.
	 *
	 *  @param orientation_boost_prior prior on the orientation boost factor,
	 *  internally called "OB_fact". In particular, the profile mass is expressed
	 *  as \f$M_{\rm correct}=M(1+\sigma_{\rm OB})\f$, where \f$\sigma_{\rm OB}\f$ is
	 *  the parameter set through this prior.
	 *
	 *  @warning \f$F_t\f$ is used only if a truncated NFW is assumed!
	 *
	 *  @warning The off-centering is related to stacks of clusters. For
	 *  details, see cbl::cosmology::HaloProfile
	 *
	 */
	void set_model_DensityProfile_cosmology (const std::vector<cbl::cosmology::CosmologicalParameter> cosmo_param, const std::vector<statistics::PriorDistribution> cosmo_prior, const statistics::PriorDistribution Rt_prior, const std::string cM_author, const statistics::PriorDistribution logM_prior, const statistics::PriorDistribution f_off_prior, const statistics::PriorDistribution sigma_off_prior, const statistics::PriorDistribution anisotropic_boost_prior, const statistics::PriorDistribution orientation_boost_prior);
	
	/**
	 *  @brief Set the profile and cosmological parameters used to model the 
	 *  cluster density profile. 
	 *
	 *  The 1-halo term has the following functional form, including the
	 *  contribution of centered and off-centered populations
	 *  of galaxy clusters (e.g. Bellagamba et al. 2019, Johnston et al. 2007, Yang et al. 2006),
	 *
	 *  \f$\Delta\Sigma_{\rm 1h}(R)=
	 *  (1-f_{\rm off})\Delta\Sigma_{\rm cen}(R)+f_{\rm off}\Delta\Sigma_{\rm off}(R),\f$
	 *
	 *  while the 2-halo term is expressed as (e.g. Bellagamba et al. 2019)
	 *
	 *  \f$\Delta\Sigma_{\rm 2h}(R) = \int\,\frac{l{\rm d}l}{2\pi}
	 *  J_2(l\theta)\frac{\bar{\rho}_{\rm m}(z)b(M,z)}{(1+z)^3D_{\rm l}^2(z)}
	 *  P(k_l,z). \f$
	 *
	 *  In particular, with this function the mass
	 *  is not a base parameter. Instead, the cluster mass-mass proxy scaling
	 *  relation is used to derive the mass. The mass is
	 *  expressed e.g. in \f$10^{14}\f$ M\f$_\odot\f$ \f$h^{-1}\f$, if
	 *  mass_pivot set through set_data_model is \f$10^{14}\f$.
	 *
	 *  Specifically, the model is the expectation value of the
	 *  total excess density profile, \f$\Delta\Sigma_{\rm tot} \f$,
	 *  expressed as
	 *
	 *  \f$\Delta\Sigma_{\rm tot} = \Delta\Sigma_{\rm 1h}+\Delta\Sigma_{\rm 2h}, \f$
	 *
	 *  derived through the following formula:
	 *
	 *  \f$\bar{\Delta\Sigma}_{\rm tot}(R_{\rm eff},\lambda_{\rm eff},z_{\rm eff}) =
	 *  \int_0^\infty\,{\rm d}M 
	 *  \Delta\Sigma_{\rm tot}(R_{\rm eff},\lambda_{\rm eff},z_{\rm eff}|M)
	 *  P(M|\lambda_{\rm eff},z_{\rm eff}), \f$
	 *
	 *  where \f$P(M|\lambda_{\rm eff},z_{\rm eff})\f$ is a log-normal whose
	 *  mean is given by the proxy-mass relation and whose rms is given by
	 *  the intrinsic scatter of such relation
	 *  (for details, see cbl::modelling::massobsrel::Modelling\_MassObservableRelation).
	 *
	 *  In particular, the intrinsic scatter is divided by the square root
	 *  of the number of clusters in the bin where the stacked profile is measured.
	 *
	 *  @param cosmo_param vector of enums containing cosmological
	 *  parameters
	 *
	 *  @param cosmo_prior vector containing the priors for
	 *  the cosmological parameters
	 *
	 *  @param z_evo redshift evolution function in the scaling relation.
	 *  Possibilities are: "Ez" (\f$ f(z_{\rm eff})=E(z_{\rm eff})/E(z_{piv})\f$),
	 *  "direct" (\f$ f(z_{\rm eff})=(1+z_{\rm eff})/(1+z_{piv}) \f$).
	 *
	 *  @param Rt_prior prior on the NFW \f$F_t\f$ defining the truncation
         *  radius, that is \f$r_t = F_tr_{\Delta}\f$ 
	 *
	 *  @param concentration_prior prior on the cluster concentration
	 *
	 *  @param f_off_prior prior on the fraction of miscentered clusters. 
	 *  This parameter makes sense only if the user models a stacked profile, 
	 *  not a single cluster profile. If a single profile is modelled, set a constant prior equal
	 *  to 0 or 1 for f_off.
	 *
	 *  @param sigma_off_prior prior on the rms of the miscentered cluster population.
	 *
	 *  @param anisotropic_boost_prior prior on the anisotropic boost factor,
	 *  internally called "AB_fact", **entering the 
	 *  model only if the 2-halo term is considered**. In particular, the 2-halo excess surface
	 *  density is expressed as \f$\Delta\Sigma_{\rm 2h,\,correct}=
	 *  \Delta\Sigma_{\rm 2h}(1+\sigma_{\rm AB})\f$, where \f$\sigma_{\rm AB}\f$ is
	 *  the parameter set through this prior.
	 *
	 *  @param orientation_boost_prior prior on the orientation boost factor,
	 *  internally called "OB_fact". In particular, the profile mass is expressed
	 *  as \f$M_{\rm correct}=M(1+\sigma_{\rm OB})\f$, where \f$\sigma_{\rm OB}\f$ is
	 *  the parameter set through this prior.
	 *
	 *  @param alpha_prior prior on the scaling relation normalization
	 *
	 *  @param beta_prior prior on the scaling relation slope
	 *
	 *  @param gamma_prior prior on the redshift evolution factor of the scaling relation
	 *
	 *  @param scatter0_prior prior on the 
	 *  constant term of the intrinsic scatter, \f$ \sigma_0 \f$
	 *
	 *  @param scatterM_prior prior on the factor in the
	 *  proxy-dependent term of the intrinsic scatter, \f$ \sigma_{\lambda} \f$
	 *
	 *  @param scatterM_exponent_prior prior on the exponent in the
	 *  proxy-dependent term of the intrinsic scatter, \f$ e_{\lambda} \f$
	 *
	 *  @param scatterz_prior prior on the factor in the
	 *  redshift-dependent term of the intrinsic scatter, \f$ \sigma_z \f$
	 *
	 *  @param scatterz_exponent_prior prior on the exponent in the
	 *  redshift-dependent term of the intrinsic scatter, \f$ e_z \f$
	 *
	 *  @warning \f$F_t\f$ is used only if a truncated NFW is assumed!
	 *
	 *  @warning The off-centering is related to stacks of clusters. For
	 *  details, see cbl::cosmology::HaloProfile
	 *
	 */
	void set_model_DensityProfile_cosmology (const std::vector<cbl::cosmology::CosmologicalParameter> cosmo_param, const std::vector<statistics::PriorDistribution> cosmo_prior, const std::string z_evo, const statistics::PriorDistribution Rt_prior, const statistics::PriorDistribution concentration_prior, const statistics::PriorDistribution f_off_prior, const statistics::PriorDistribution sigma_off_prior, const statistics::PriorDistribution anisotropic_boost_prior, const statistics::PriorDistribution orientation_boost_prior, const statistics::PriorDistribution alpha_prior, const statistics::PriorDistribution beta_prior, const statistics::PriorDistribution gamma_prior, const statistics::PriorDistribution scatter0_prior, const statistics::PriorDistribution scatterM_prior, const statistics::PriorDistribution scatterM_exponent_prior, const statistics::PriorDistribution scatterz_prior, const statistics::PriorDistribution scatterz_exponent_prior);

	///@}
     };

      
      /**
       * @brief Compute the excess density profile model
       * in all the radial bins
       *
       * @param radius the radius array
       *
       * @param inputs model inputs
       *
       * @param parameter model parameters
       *
       * @return halo model in each radial bin
       *
       */
      std::vector<double> model_density (const std::vector<double> radius, const std::shared_ptr<void> inputs, std::vector<double> &parameter);
      
      /**
       * @brief Compute the excess density profile model
       * in all the radial bins. The mass is derived from the 
       * scaling relation.
       *
       * @param radius the radius array
       *
       * @param inputs model inputs
       *
       * @param parameter model parameters
       *
       * @return halo model in each radial bin
       *
       */
      std::vector<double> model_density_scaling_relation (const std::vector<double> radius, const std::shared_ptr<void> inputs, std::vector<double> &parameter);
    
    }
  }
}

#endif
