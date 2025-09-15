/********************************************************************
 *  Copyright (C) 2024 by Federico Marulli and Giorgio Lesci        *
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
 *  @file Headers/Modelling_DensityProfile_Combined.h
 *
 *  @brief The class Modelling_DensityProfile_Combined
 *
 *  This file defines the interface of the class Modelling_DensityProfile_Combined,
 *  used to model the density profile of galaxy clusters, defined by
 *  the class Modelling_DensityProfile, in different bins and in a combined way
 *
 *  @authors Giorgio Lesci (and Federico Marulli)
 *
 *  @authors giorgio.lesci2@unibo.it (and federico.marulli3@unibo.it)
 */

#ifndef __MODELLINGDPROFILECOMB__
#define __MODELLINGDPROFILECOMB__

#include "Modelling_DensityProfile.h"


// ===================================================================================================


namespace cbl {

  namespace modelling {

    /**
     *  @brief The namespace of the <B> cluster projected density profile
     *  modelling </B>
     *  
     *  The \e modelling::densityprofilecomb namespace contains all the functions
     *  and classes to model the surface density excess profile of galaxy clusters
     */
    namespace densityprofilecomb {
      
      /**
       *  @struct STR_ProfileComb_data_model
       *  @brief the structure STR_ProfileComb_data_model
       *
       *  This structure contains the data used for statistical
       *  analyses of cluster surface density excess profiles
       */
      struct STR_ProfileComb_data_model {
        
        /// vector of Modelling_DensityProfile objects
        std::vector<cbl::modelling::densityprofile::Modelling_DensityProfile> single_object;
      
	/// Fiducial cosmology pointer
	std::shared_ptr<cosmology::Cosmology> cosmology;

	/// Cosmological parameters
	std::vector<std::string> Cpar;
	
	/// parameter types
	std::vector<statistics::ParameterType> Par_type;

	/// parameter names
	std::vector<std::string> Par_string;

	/// parameter priors
	std::vector<statistics::PriorDistribution> param_prior;
	
	/// Minimum of the interval of cluster redshifts
	double z_min;
	
	/// Maximum of the interval of cluster redshifts
	double z_max;
	
	/// Minimum of the interval of cluster mass proxies
	double mass_proxy_min;
	
	/// Maximum of the interval of cluster mass proxies
	double mass_proxy_max;

	/**
	 *  @brief default constructor
	 */
	STR_ProfileComb_data_model () = default;
      };
    
      /**
       *  @class Modelling_DensityProfile_Combined
       *  Modelling_DensityProfile_Combined.h
       *  "Headers/Modelling_DensityProfile_Combined.h"
       *
       *  @brief The class Modelling_DensityProfile_Combined
       *
       *  This file defines the interface of the base class
       *  Modelling_DensityProfile_Combined, used for modelling
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
      class Modelling_DensityProfile_Combined : public Modelling {
      
      protected:

	/// the container of parameters for the density model computation
	STR_ProfileComb_data_model m_data_model;


      public:

	/**
	 *  @name Constructors/destructors
	 */
	///@{

	/**
	 *  @brief default constuctor
	 *  _DensityProfile_Combined
	 */
	Modelling_DensityProfile_Combined () = default;
	
	/**
	 *  @brief constuctor for the modelling of
	 *  stacked density excess profile of galaxy cluster. Cosmological units are forced.
	 *
	 *  Aside from the probe parameter, all the other inputs are vectors whose
	 *  most external size equals the number of Modelling_DensityProfile objects
	 *  that will be combined
	 *
	 *  @param dataset cluster profile dataset
	 *
	 *  @param probe "DeltaSigma" \f$\rightarrow\f$ excess surface density; 
	 *  "g" \f$\rightarrow\f$ reduced shear
	 *
	 *  @param zl_bin_edges_for_N_zs lens (i.e. cluster) redshift **bin edges**. 
	 *  For each lens redshift bin defined this way,
	 *  an \f$N(z_{\rm s})\f$ distribution must be defined.
	 *  This parameter is used only in the case of reduced shear modelling
	 *
	 *  @param zs **central bin values** where the binned \f$N(z_{\rm s})\f$
	 *  distributions are computed, where 
	 *  \f$z_{\rm s}\f$ is the background source redshift. 
	 *  This parameter is used only in the case of reduced shear modelling. For
	 *  details, see cbl::cosmology::HaloProfile::g
	 *
	 *  @param N_zs the binned \f$N(z_{\rm s})\f$ distributions for each
	 *  lens redshift bin.
	 *  This parameter is used only in the case of reduced shear modelling.
	 *
	 *  @param purity_background purity of the background sample, \f$P_{\rm bkg}\f$,
	 *  defined in the lens redshift bins defined above. 
	 *  It must range between 0 and 1. The reduced shear models are multiplied
	 *  by \f$P_{\rm bkg}\f$.
	 *  This parameter is used only in the case of reduced shear modelling.
	 *
	 */
	Modelling_DensityProfile_Combined (const std::vector<std::shared_ptr<cbl::data::Data>> dataset, const std::string probe, const std::vector<std::vector<double>> zl_bin_edges_for_N_zs={}, const std::vector<std::vector<double>> zs={}, const std::vector<std::vector<std::vector<double>>> N_zs={}, const std::vector<std::vector<double>> purity_background={});
	
	/**
	 *  @brief default destructor
	 *  
	 */
	virtual ~Modelling_DensityProfile_Combined () = default;

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
	STR_ProfileComb_data_model data_model () { return m_data_model; }
	
	/**
	 *  @brief Set the data used to construct models of
	 *  average halo profiles in bins of redshift and mass proxy, 
	 *  where the mass is derived from a scaling relation and the expected
	 *  mean profile is weighted over the observed distribution of clusters.
	 *
	 *  To properly define redshift and mass proxy arrays, use
	 *  the files generated via 
	 *  cbl::measure::stackprofile::StackedDensityProfile::measure.
	 *  In these files, indeed, only the clusters that contribute
	 *  to the stacked lensing signal are listed.
	 *
	 *  **Compared to cbl::modelling::densityprofile::Modelling_DensityProfile::set_data_model**,
	 *  redshift_obs, mass_proxy_obs, Pz_params, Pproxy_params have an additional dimension
	 *  
	 *  @param cosmology the cosmological model
	 *
	 *  @param profile_author author(s) of the cluster density profile.
	 *  See available options in cbl::cosmology::HaloProfile
	 *
	 *  @param two_halo if true, include the 2-halo contribution
	 *
	 *  @param halo_def the halo definition; available options are:
	 *  "critical", "vir", "mean"
	 *
	 *  @param Delta overdensity factor which needs to be multiplied 
	 *  to the critical density in order to define an overdensity
	 *
	 *  @param redshift_obs observed redshift values of the clusters contributing to
	 *  the stacked signal
	 *
	 *  @param mass_proxy_obs observed mass proxy values of the clusters contributing to
	 *  the stacked signal
	 *
	 *  @param Pz_params three parameters must be provided, namely mean, standard deviation,
	 *  and \f$\alpha\f$ (in this order). These are the parameters of the skew-normal
	 *  distribution \f$P(x|z_{\rm ob},\lambda_{\rm ob})\f$, where
	 *  \f$x = (z_{\rm tr}-z_{\rm ob})\,/\,(1+z_{\rm ob})\f$. Estimates of the true redshift values
	 *  are extracted from this distribution. To check the extracted numbers, use
	 *  cbl::modelling::densityprofile::Modelling_DensityProfile::get_extracted_true_observables().
	 *  Note that if \f$\alpha=0\f$ the distribution is Gaussian. In addition,
	 *  **if the mean is zero and the standard deviation is zero, the observed values are used**.
	 *
	 *  @param Pproxy_params three parameters must be provided, namely mean, standard deviation,
	 *  and \f$\alpha\f$ (in this order). These are the parameters of the skew-normal
	 *  distribution \f$P(x|z_{\rm ob},\lambda_{\rm ob})\f$, where
	 *  \f$x = (\lambda_{\rm tr}-\lambda_{\rm ob})\,/\,\lambda_{\rm ob}\f$. Estimates of the true mass proxy values
	 *  are extracted from this distribution. To check the extracted numbers, use
	 *  cbl::modelling::densityprofile::Modelling_DensityProfile::get_extracted_true_observables().
	 *  Note that if \f$\alpha=0\f$ the distribution is Gaussian. In addition,
	 *  **if the mean is zero and the standard deviation is zero, the observed values are used**.
	 *
	 *  @param redshift_pivot redshift pivot in the scaling relation
	 *
	 *  @param proxy_pivot proxy pivot in the scaling relation
	 *
	 *  @param logM_base base of the mass logarithm
	 *
	 *  @param mass_pivot the mass pivot
	 *
	 *  @param bias_author author(s) who proposed the bias model.
	 *  See cbl::cosmology::Bias::Bias.
	 *
	 *  @param method_Pk method used for the computation of the
	 *  power spectrum.
	 *
	 *  @param interp_type method to interpolate the power
	 *  spectrum.
	 *
	 *  In case the 2-halo term is considered, or if the DK14 profile is
	 *  assumed, the code also computes the splashback radius (Rsp) and
	 *  the splashback radius divided by \f$r_\Delta\f$ (Rsp_rDelta). 
	 *  Both are provided as derived parameters of the MCMC.
	 *  Specifically, the splashback radius is computed as the minimum
	 *  of \f$ {\rm d}\log\rho \, / \, {\rm d}\log R \f$, where
	 *  \f$\log\f$ is the base 10 logarithm.
	 *  
	 */
	void set_data_model (const std::shared_ptr<cosmology::Cosmology> cosmology, const std::string profile_author, const bool two_halo, const std::string halo_def, const double Delta, const std::vector<std::vector<double>> redshift_obs, const std::vector<std::vector<double>> mass_proxy_obs, const std::vector<std::vector<double>> Pz_params, const std::vector<std::vector<double>> Pproxy_params, const double redshift_pivot, const double proxy_pivot, const double logM_base, const double mass_pivot, const std::string bias_author="Tinker", const std::string method_Pk="EisensteinHu", std::string interp_type="Linear");
	
	/**
	 *  @brief Analogous to
	 *  cbl::modelling::densityprofile::Modelling_DensityProfile::set_model_DensityProfile_cosmology.
	 *  The only difference is that here a vector of purity priors must
	 *  be provided.
	 *
	 *  @param cosmo_param vector of enums containing cosmological
	 *  parameters
	 *
	 *  @param cosmo_prior vector containing the priors for the
	 *  cosmological parameters
	 *
	 *  @param z_evo redshift evolution function in the scaling
	 *  relation.  Possibilities are: "E_z" (\f$ f(z)=E(z)/E(z_{piv})\f$), 
	 *  "direct" (\f$f(z)=(1+z)/(1+z_{piv}) \f$).
	 *
	 *  @param Rt_prior prior on the NFW \f$F_t\f$ defining the
         *  truncation radius, that is \f$r_t = F_tr_{\Delta}\f$
	 *
	 *  @param c0_prior prior on \f$ c_0 \f$
	 *
	 *  @param cM_prior prior on \f$ c_M \f$
	 *
	 *  @param cz_prior prior on \f$ c_z \f$
	 *
	 *  @param f_off0_prior prior on \f$ f_{\rm off,0} \f$
	 *
	 *  @param f_offM_prior prior on \f$ f_{\rm off,M} \f$
	 *
	 *  @param f_offz_prior prior on \f$ f_{\rm off,z} \f$
	 *
	 *  @param sigma_off0_prior prior on \f$ \sigma_{\rm off,0}
	 *  \f$
	 *
	 *  @param sigma_offM_prior prior on \f$ \sigma_{\rm
	 *  off,M} \f$
	 *
	 *  @param sigma_offz_prior prior on \f$ \sigma_{\rm off,z}
	 *  \f$
	 *
	 *  @param anisotropic_boost_prior prior on the anisotropic
	 *  boost factor, internally called "AB_fact", **entering the
	 *  model only if the 2-halo term is considered**. In
	 *  particular, the 2-halo excess surface density is expressed
	 *  as \f$\Sigma_{\rm 2h,\,correct}= \Sigma_{\rm
	 *  2h}(1+\sigma_{\rm AB})\f$, where \f$\sigma_{\rm AB}\f$ is
	 *  the parameter set through this prior.
	 *
	 *  @param orientation_boost_prior prior on the orientation boost factor, \f$\mathcal{B}_O\f$,
	 *  internally called "OB_fact". In particular, the average profile is multiplied
	 *  by \f$(1+\mathcal{B}_O)\f$.
	 *
	 *  @param alpha_prior prior on the scaling relation normalization
	 *
	 *  @param beta_prior prior on the scaling relation slope
	 *
	 *  @param gamma_prior prior on the redshift evolution factor
	 *  of the scaling relation
	 *
	 *  @param scatter0_prior prior on the constant term of the
	 *  intrinsic scatter, \f$ \sigma_0 \f$
	 *
	 *  @param scatterM_prior prior on the factor in the
	 *  proxy-dependent term of the intrinsic scatter, \f$
	 *  \sigma_{\lambda} \f$
	 *
	 *  @param scatterM_exponent_prior prior on the exponent in the
	 *  proxy-dependent term of the intrinsic scatter, \f$ e_{\lambda} \f$
	 *
	 *  @param scatterz_prior prior on the factor in the
	 *  redshift-dependent term of the intrinsic scatter, \f$
	 *  \sigma_z \f$
	 *
	 *  @param scatterz_exponent_prior prior on the exponent in
	 *  the redshift-dependent term of the intrinsic scatter, \f$
	 *  e_z \f$
	 *
	 *  @param purity_prior prior on the purity of the cluster samples, defined
	 *  between 0 and 1
	 *
	 *  @warning \f$F_t\f$ is used only if a truncated NFW is
	 *  assumed!
	 *
	 *  @warning The off-centering is related to stacks of
	 *  clusters. For details, see cbl::cosmology::HaloProfile.
	 */
	void set_model_DensityProfile_cosmology (const std::vector<std::string> cosmo_param, const std::vector<statistics::PriorDistribution> cosmo_prior, const std::string z_evo, const statistics::PriorDistribution Rt_prior, const statistics::PriorDistribution c0_prior, const statistics::PriorDistribution cM_prior, const statistics::PriorDistribution cz_prior, const statistics::PriorDistribution f_off0_prior, const statistics::PriorDistribution f_offM_prior, const statistics::PriorDistribution f_offz_prior, const statistics::PriorDistribution sigma_off0_prior, const statistics::PriorDistribution sigma_offM_prior, const statistics::PriorDistribution sigma_offz_prior, const statistics::PriorDistribution anisotropic_boost_prior, const statistics::PriorDistribution orientation_boost_prior, const statistics::PriorDistribution alpha_prior, const statistics::PriorDistribution beta_prior, const statistics::PriorDistribution gamma_prior, const statistics::PriorDistribution scatter0_prior, const statistics::PriorDistribution scatterM_prior, const statistics::PriorDistribution scatterM_exponent_prior, const statistics::PriorDistribution scatterz_prior, const statistics::PriorDistribution scatterz_exponent_prior, const std::vector<statistics::PriorDistribution> purity_prior);
	
	/**
       * @brief write the model computing 16th, 50th and 84th 
       * percentiles from the MCMC
       *
       * @param output_dir output directory
       *
       * @param output_file tail of the output files
       *
       * @param start minimum chain position
       *
       * @param thin step used for diluting the chain
       *
       * @param xx x points where the model is computed. If not 
       * provided, the x points set for the MCMC are considered
       *
       */
       void write_model_from_combined_chain (const std::string output_dir, const std::string output_file, const int start, const int thin, const std::vector<double> xx={});
	
	///@}

      };
      
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
      std::vector<double> combined_model_density_scaling_relation_evolving_concentration_offcentering (const std::vector<double> radius, const std::shared_ptr<void> inputs, std::vector<double> &parameter);
      
    }
  }
}

#endif
