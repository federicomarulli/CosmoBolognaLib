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

#include "Cosmology.h"
#include "StackedDensityProfile.h"
#include "Modelling_MassObservableRelation.h"
#include "Modelling.h"
#include "Cluster.h"


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
	std::function<double(const double, cbl::catalogue::Cluster)> conc_func;
      
	/// Fiducial cosmology pointer
	std::shared_ptr<cosmology::Cosmology> cosmology;
	
	/// Cosmological parameters
	std::vector<cosmology::CosmologicalParameter> Cpar;
	
	/// Modelling_MassObservableRelation object pointer
	std::shared_ptr<modelling::massobsrel::Modelling_MassObservableRelation> scaling_relation;
	
	/// Cluster object pointer
	std::shared_ptr<catalogue::Cluster> cluster;
	
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
	
	/// Truncation factor
	double trunc_fact;	
	
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
	 *  Available options are: "NFW" (trucated NFW)
	 *
	 *  @param _2halo if true, compute the 2-halo contribution
	 *
	 *  @param halo_def the halo definition; available options are:
	 *  "200" \f$\rightarrow\f$ all
	 *  matter within the radius \f$r_{200}\f$ for which the mean
	 *  internal density is 200 times the critical density.
	 *
	 */
	Modelling_DensityProfile (const std::shared_ptr<cbl::measure::stackprofile::StackedDensityProfile> profile, const std::string profile_author="NFW", const bool _2halo=false, const std::string halo_def="200");
	
	/**
	 *  @brief constuctor for the modelling of
	 *  stacked density excess profile of galaxy cluster. Cosmological units are forced.
	 *
	 *  @param dataset cluster profile dataset
	 *
	 *  @param profile_author author(s) of the cluster density profile.
	 *  Available options are: "NFW" (trucated NFW)
	 *
	 *  @param _2halo if true, compute the 2-halo contribution
	 *
	 *  @param halo_def the halo definition; available options are:
	 *  "200" \f$\rightarrow\f$ all
	 *  matter within the radius \f$r_{200}\f$ for which the mean
	 *  internal density is 200 times the critical density.
	 *
	 */
	Modelling_DensityProfile (const std::shared_ptr<cbl::data::Data> dataset, const std::string profile_author="NFW", const bool _2halo=false, const std::string halo_def="200");
	
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
	 *  @brief set the data used to construct generic models of
	 *  cluster profiles
	 *  
	 *  @param cosmology the cosmological model
	 *
	 *  @param cluster the cluster object
	 *
	 *  @param redshift cluster redshift
	 *
	 *  @param contrast density contrast with respect 
	 *  to the critical density (i.e. if equal to 200, \f$M_{200}\f$ is considered)
	 *
	 *  @param trunc_fact truncation factor
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
	void set_data_model (const cbl::cosmology::Cosmology cosmology, const cbl::catalogue::Cluster, const double redshift, const double contrast, const double trunc_fact, const double logM_base, const double mass_pivot, const std::string bias_author="Tinker", const std::string method_Pk="EisensteinHu", std::string interp_type="Linear");
	
	/**
	 *  @brief set the data used to construct generic models of
	 *  cluster profiles, where the mass is derived from a scaling
	 *  relation expressed as follows:
	 *
	 *  \f$\log M = \alpha + \beta 
	 *  \log (\lambda/\lambda_{\rm piv}) + \gamma \log (f(z)).\f$
	 *  
	 *  @param cosmology the cosmological model
	 *
	 *  @param cluster the cluster object
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
	 *  @param trunc_fact truncation factor
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
	void set_data_model (const cosmology::Cosmology cosmology, const catalogue::Cluster cluster, const double redshift, const double mass_proxy, const double redshift_pivot, const double proxy_pivot, const double contrast, const double trunc_fact, const double logM_base, const double mass_pivot, const std::string bias_author="Tinker", const std::string method_Pk="EisensteinHu", std::string interp_type="Linear");

	/**
	 *  @brief Set the profile and cosmological parameters used to model the 
	 *  cluster density profile. 
	 *
	 *  The 1-halo term has the following functional form, including the
	 *  contribution of centered and off-centered populations
	 *  of galaxy clusters (e.g. Bellagamba et al. 2019, Johnston et al. 2007, Yang et al. 2006):
	 *
	 *  \f$\Sigma_{\rm 1h}(R)=(1-f_{\rm off})\Sigma_{\rm cen}(R)+f_{\rm off}\Sigma_{\rm off}(R)\f$
	 *
	 *  @param cosmo_param vector of enums containing cosmological
	 *  parameters
	 *
	 *  @param cosmo_prior vector containing the priors for
	 *  the cosmological parameters
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
	 */
	void set_model_DensityProfile_cosmology (const std::vector<cbl::cosmology::CosmologicalParameter> cosmo_param, const std::vector<statistics::PriorDistribution> cosmo_prior, const statistics::PriorDistribution concentration_prior, const statistics::PriorDistribution logM_prior, const statistics::PriorDistribution f_off_prior, const statistics::PriorDistribution sigma_off_prior);
	
	/**
	 *  @brief Set the profile and cosmological parameters used to model the 
	 *  cluster density profile. The concentration is a derived parameter,
	 *  computed through a concentration-mass relation.
	 *
	 *  The 1-halo term has the following functional form, including the
	 *  contribution of centered and off-centered populations
	 *  of galaxy clusters (e.g. Bellagamba et al. 2019, Johnston et al. 2007, Yang et al. 2006):
	 *
	 *  \f$\Sigma_{\rm 1h}(R)=(1-f_{\rm off})\Sigma_{\rm cen}(R)+f_{\rm off}\Sigma_{\rm off}(R)\f$
	 *
	 *  @param cosmo_param vector of enums containing cosmological
	 *  parameters
	 *
	 *  @param cosmo_prior vector containing the priors for
	 *  the cosmological parameters
	 *
	 *  @param cM_author author(s) who proposed the 
	 *  concentration-mass relation. Possibilities are:
	 *  "Duffy" (Duffy et al. 2008)
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
	 */
	void set_model_DensityProfile_cosmology (const std::vector<cbl::cosmology::CosmologicalParameter> cosmo_param, const std::vector<statistics::PriorDistribution> cosmo_prior, const std::string cM_author, const statistics::PriorDistribution logM_prior, const statistics::PriorDistribution f_off_prior, const statistics::PriorDistribution sigma_off_prior);
	
	/**
	 *  @brief Set the profile and cosmological parameters used to model the 
	 *  cluster density profile. 
	 *
	 *  In particular, with this function the mass
	 *  is a parameter derived from the cluster mass-mass proxy scaling
	 *  relation, with the following functional form:
	 *
	 *  \f$\log M = \alpha + \beta 
	 *  \log (\lambda_{\rm eff}/\lambda_{\rm piv}) + \gamma \log (f(z_{\rm eff})).\f$
	 *
	 *  WARNING: the only way to have a dependency on the intrinsic scatter
	 *  parameters is to define a user-defined likelihood, whose covariance
	 *  depends on the intrinsic scatter. In particular the intrinsic
	 *  scatter, \f$\sigma_{\rm intr}\f$, is expressed as
	 *
	 *  \f$ \sigma_{\rm intr} = \sigma_0 + \sigma_{\lambda} 
	 *  \log (\lambda/\lambda_{\rm piv})^{e_{\lambda}} + \sigma_z \log (f(z))^{e_z}.\f$
	 *
	 *  The 1-halo term has the following functional form, including the
	 *  contribution of centered and off-centered populations
	 *  of galaxy clusters (e.g. Bellagamba et al. 2019, Johnston et al. 2007, Yang et al. 2006):
	 *
	 *  \f$\Sigma_{\rm 1h}(R)=(1-f_{\rm off})\Sigma_{\rm cen}(R)+f_{\rm off}\Sigma_{\rm off}(R)\f$
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
	 *  @param z_evo redshift evolution function in the scaling relation.
	 *  Possibilities are: "Ez" (\f$ f(z_{\rm eff})=E(z_{\rm eff})/E(z_{piv})\f$),
	 *  "direct" (\f$ f(z_{\rm eff})=(1+z_{\rm eff})/(1+z_{piv}) \f$).
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
	 */
	void set_model_DensityProfile_cosmology (const std::vector<cbl::cosmology::CosmologicalParameter> cosmo_param, const std::vector<statistics::PriorDistribution> cosmo_prior, const std::string z_evo, const statistics::PriorDistribution concentration_prior, const statistics::PriorDistribution f_off_prior, const statistics::PriorDistribution sigma_off_prior, const statistics::PriorDistribution alpha_prior, const statistics::PriorDistribution beta_prior, const statistics::PriorDistribution gamma_prior, const statistics::PriorDistribution scatter0_prior, const statistics::PriorDistribution scatterM_prior, const statistics::PriorDistribution scatterM_exponent_prior, const statistics::PriorDistribution scatterz_prior, const statistics::PriorDistribution scatterz_exponent_prior);
	
	/**
	 *  @brief Set the profile and cosmological parameters used to model the 
	 *  cluster density profile. The concentration is a derived parameter,
	 *  computed through a concentration-mass relation.
	 *
	 *  In particular, with this function the mass
	 *  is a parameter derived from the cluster mass-mass proxy scaling
	 *  relation, with the following functional form:
	 *
	 *  \f$\log M = \alpha + \beta 
	 *  \log (\lambda_{\rm eff}/\lambda_{\rm piv}) + \gamma \log (f(z_{\rm eff})).\f$
	 *
	 *  WARNING: the only way to have a dependency on the intrinsic scatter
	 *  parameters is to define a user-defined likelihood, whose covariance
	 *  depends on the intrinsic scatter. In particular the intrinsic
	 *  scatter, \f$\sigma_{\rm intr}\f$, is expressed as
	 *
	 *  \f$ \sigma_{\rm intr} = \sigma_0 + \sigma_{\lambda} 
	 *  \log (\lambda/\lambda_{\rm piv})^{e_{\lambda}} + \sigma_z \log (f(z))^{e_z}.\f$
	 *
	 *  The 1-halo term has the following functional form, including the
	 *  contribution of centered and off-centered populations
	 *  of galaxy clusters (e.g. Bellagamba et al. 2019, Johnston et al. 2007, Yang et al. 2006):
	 *
	 *  \f$\Sigma_{\rm 1h}(R)=(1-f_{\rm off})\Sigma_{\rm cen}(R)+f_{\rm off}\Sigma_{\rm off}(R)\f$
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
	 *  @param z_evo redshift evolution function in the scaling relation.
	 *  Possibilities are: "Ez" (\f$ f(z_{\rm eff})=E(z_{\rm eff})/E(z_{piv})\f$),
	 *  "direct" (\f$ f(z_{\rm eff})=(1+z_{\rm eff})/(1+z_{piv}) \f$).
	 *
	 *  @param cM_author author(s) who proposed the 
	 *  concentration-mass relation. Possibilities are:
	 *  "Duffy" (Duffy et al. 2008)
	 *
	 *  @param f_off_prior prior on the fraction of miscentered clusters. 
	 *  This parameter makes sense only if the user models a stacked profile, 
	 *  not a single cluster profile. If a single profile is modelled, set a constant prior equal
	 *  to 0 or 1 for f_off.
	 *
	 *  @param sigma_off_prior prior on the rms of the miscentered cluster population.
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
	 */
	void set_model_DensityProfile_cosmology (const std::vector<cbl::cosmology::CosmologicalParameter> cosmo_param, const std::vector<statistics::PriorDistribution> cosmo_prior, const std::string z_evo, const std::string cM_author, const statistics::PriorDistribution f_off_prior, const statistics::PriorDistribution sigma_off_prior, const statistics::PriorDistribution alpha_prior, const statistics::PriorDistribution beta_prior, const statistics::PriorDistribution gamma_prior, const statistics::PriorDistribution scatter0_prior, const statistics::PriorDistribution scatterM_prior, const statistics::PriorDistribution scatterM_exponent_prior, const statistics::PriorDistribution scatterz_prior, const statistics::PriorDistribution scatterz_exponent_prior);

	///@}
     };

      
      /**
       * @brief Compute the truncated NFW density profile as a function
       * of radius and redshift
       *
       * @param cosmology the cosmology
       *
       * @param cluster the cluster object 
       *
       * @param r radius where the profile is computed (Mpc/h)
       *
       * @param redshift cluster redshift
       *
       *  @param contrast density contrast with respect 
       *  to the critical density (i.e. if equal to 200, M_200 is considered)
       *
       * @param trunc_fact truncation factor
       *
       * @return the model for \f$\Delta\Sigma(r)\f$, expressed in units
       * of \f$h\f$ M\f$_\odot\f$/pc\f$^2\f$
       *
       */
      double NFW_truncated (cosmology::Cosmology cosmology, catalogue::Cluster cluster, const double r, const double redshift, const double contrast, const double trunc_fact);
      
      /**
       * @brief Compute the 2-halo term
       *
       * @param cosm the cosmology
       *
       * @param radius radius where the profile is computed (Mpc/h)
       *
       * @param mass cluster mass
       *
       * @param redshift cluster redshift
       *
       *  @param contrast density contrast with respect 
       *  to the critical density (i.e. if equal to 200, M_200 is considered)
       *
       * @param bias_author author(s) who proposed the bias; valid authors are: 
       *  ST99 (Sheth & Tormen 1999), SMT01 (Sheth, Mo & Tormen 2001), 
       *  SMT01_WL04 (Sheth, Mo & Tormen 2001 with the correction of Warren 2004), 
       *  Tinker (Tinker et al. 2010)
       *
       *  @param method_Pk method used for the computation of the power spectrum.
       *
       *  @param interp_type method to interpolate the power spectrum.
       *
       * @return the 2-halo term for \f$\Delta\Sigma(r)\f$, expressed in units
       * of \f$h\f$ M\f$_\odot\f$/pc\f$^2\f$
       *
       */
      double two_halo (cosmology::Cosmology cosm, const double radius, const double mass, const double redshift, const double contrast, const std::string bias_author, const std::string method_Pk, std::string interp_type);
      
      /**
       * @brief Compute the truncated NFW density profile model
       * in all the radial bins
       *
       * @param radius the radius array
       *
       * @param inputs model inputs
       *
       * @param parameter model parameters
       *
       * @return the truncated NFW 1-halo model in each radial bin
       *
       */
      std::vector<double> model_NFW_truncated (const std::vector<double> radius, const std::shared_ptr<void> inputs, std::vector<double> &parameter);
      
      /**
       * @brief Compute the truncated NFW density profile model
       * + the 2halo term (see Oguri & Takada 2011; Sereno et al. 2017)
       * in all the radial bins
       *
       * @param radius the radius array
       *
       * @param inputs model inputs
       *
       * @param parameter model parameters
       *
       * @return the truncated NFW 1-halo + 2-halo model in each radial bin
       *
       */
      std::vector<double> model_NFW_truncated_2halo (const std::vector<double> radius, const std::shared_ptr<void> inputs, std::vector<double> &parameter);
      
      /**
       * @brief Compute the truncated NFW density profile model
       * in all the radial bins. The mass is derived from the 
       * scaling relation.
       *
       * @param radius the radius array
       *
       * @param inputs model inputs
       *
       * @param parameter model parameters
       *
       * @return the truncated NFW 1-halo model in each radial bin
       *
       */
      std::vector<double> model_NFW_truncated_scaling_relation (const std::vector<double> radius, const std::shared_ptr<void> inputs, std::vector<double> &parameter);
      
      /**
       * @brief Compute the truncated NFW density profile model
       * + the 2halo term (see Oguri & Takada 2011; Sereno et al. 2017) 
       * in all the radial bins. The mass is derived from the 
       * scaling relation.
       *
       * @param radius the radius array
       *
       * @param inputs model inputs
       *
       * @param parameter model parameters
       *
       * @return the truncated NFW 1-halo model in each radial bin
       *
       */
      std::vector<double> model_NFW_truncated_2halo_scaling_relation (const std::vector<double> radius, const std::shared_ptr<void> inputs, std::vector<double> &parameter);
    
    }
  }
}

#endif
