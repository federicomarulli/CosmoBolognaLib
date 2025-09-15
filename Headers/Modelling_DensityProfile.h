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

#include "StackedDensityProfile.h"
#include "Modelling_DensityProfile_Auxiliary.h"


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
        
        /// the probe
        std::string probe;
        
        /// central bin values where the binned \f$N(z_{\rm s})\f$ distribution is computed
        std::vector<double> zs;
        
        /// binned \f$N(z_{\rm s})\f$ distributions
        std::vector<std::vector<double>> N_zs;
        
        /// lens redshift bin edges where the \f$N(z_{\rm s})\f$ distributions are defined
        std::vector<double> zl_bin_edges_for_N_zs;
        
        /// in which bin the cluster redshifts lie
        std::vector<int> N_zs_index;
        
        /// the purity of the background source catalogue as a function of the lens redshift
        std::vector<double> purity_background;
        
        /// the purity of the background source catalogue interpolated as a function of the lens redshift
        cbl::glob::FuncGrid purity_background_interp;
        
        /// cluster completeness
        cbl::glob::FuncGrid2D completeness_interp;

	/// parameter types
	std::vector<statistics::ParameterType> Par_type;

	/// parameter names
	std::vector<std::string> Par_string;

	/// parameter priors
	std::vector<statistics::PriorDistribution> param_prior;
	
	/// boolean for the two-halo contribution
	bool two_halo;
      
        /// function returning the concentration 
	std::function<double(const double, cbl::cosmology::HaloProfile)> conc_func;
	
	/// function returning the concentration from a log(concentration)-log(mass) relation
	std::function<double(const double, const double, const double, const double, const double)> conc_scaling_relation_func;
      
	/// Fiducial cosmology pointer
	std::shared_ptr<cosmology::Cosmology> cosmology;
	
	/// Cosmological parameters
	std::vector<std::string> Cpar;
	
	/// Modelling_MassObservableRelation object pointer
	std::shared_ptr<modelling::massobsrel::Modelling_MassObservableRelation> scaling_relation;
	
	/// HaloProfile object pointer
	std::shared_ptr<cosmology::HaloProfile> halo_profile;
	
	/// overdensity
	double Delta;
	
	/// author(s) of the cluster density profile
        std::string profile_author;
        
        /// halo definition
        std::string halo_def;
        
        /// minimum value an extracted quantity (redshift or mass proxy) can assume
        double min_extracted;
	
	/// Redshift
	double redshift;
	
	/// Redshift vector
	std::vector<double> redshifts;
	
	/// Mass proxy
	double mass_proxy;
	
	/// Mass proxy vector
	std::vector<double> mass_proxies;
	
	/// Minimum of the interval of observed cluster redshifts
	double z_min;
	
	/// Maximum of the interval of observed cluster redshifts
	double z_max;
	
	/// Minimum of the interval of observed cluster mass proxies
	double mass_proxy_min;
	
	/// Maximum of the interval of observed cluster mass proxies
	double mass_proxy_max;
	
	/// Redshift error
	double z_error;
	
	/// Function computing the mean of \f$P(\lambda_{\rm ob}|\lambda,z)\f$
	std::function<double(const double, const double, const double, const double, const double)> Plambda_mean_fc;
	
	/// Function computing the standard deviation of \f$P(\lambda_{\rm ob}|\lambda,z)\f$
	std::function<double(const double, const double, const double, const double, const double)> Plambda_std_fc;
	
	/// \f$A_\mu\f$, see cbl::modelling::densityprofile::Modelling_DensityProfile::set_P_proxy
	double Plambda_A_mu;
	
	/// \f$B_\mu\f$, see cbl::modelling::densityprofile::Modelling_DensityProfile::set_P_proxy
	double Plambda_B_mu;
	
	/// \f$C_\mu\f$, see cbl::modelling::densityprofile::Modelling_DensityProfile::set_P_proxy
	double Plambda_C_mu;
	
	/// \f$A_\sigma\f$, see cbl::modelling::densityprofile::Modelling_DensityProfile::set_P_proxy
	double Plambda_A_sigma;
	
	/// \f$B_\sigma\f$, see cbl::modelling::densityprofile::Modelling_DensityProfile::set_P_proxy
	double Plambda_B_sigma;
	
	/// \f$C_\sigma\f$, see cbl::modelling::densityprofile::Modelling_DensityProfile::set_P_proxy
	double Plambda_C_sigma;
	
	/// Redshift pivot in the scaling relation
	double redshift_pivot;
	
	/// Mass proxy pivot in the scaling relation
	double proxy_pivot;	
	
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
	
	/// author(s) of the theoretical halo mass function
	std::string MF_author;
	
	/// if true, the profiles are averaged over the theoretical mass function
	bool isTheoretical_MF = false;
	
	/// minimum k for P(k)
	double k_min;
	
	/// maximum k for P(k)
	double k_max;
	
	/// k vector for the P(k) interpolation
	std::vector<double> kk;
	
	/// if true the output files created by the Boltzmann solver are stored; if false the output files are removed
	bool store_output;
	
	/// output_root of the parameter file used to compute the power spectrum and σ(mass); it can be any name
	std::string output_root;
	
	/// 0: don't normalize the power spectrum; 1: normalize the power spectrum; -1: normalize only if sigma8 is set
	int norm;
	
	/// the precision for P(k) interpolation
	double prec;
	
	/// mass vector used to interpolate sigma(M)
	std::vector<double> Mass_vector;

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
	
	/// if true, the mass is a parameter derived from the scaling relation
	bool m_mass_is_derived;
        
        /// the probe
        std::string m_probe;
        
        /**
	 *  @brief Set the probe's model. 
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
	 *  This parameter is used only in the case of reduced shear modelling.
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
	void m_set_probe_model (const std::string probe, const std::vector<double> zl_bin_edges_for_N_zs, const std::vector<double> zs, const std::vector<std::vector<double>> N_zs, const std::vector<double> purity_background);


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
	 *  @param dataset the object of cbl::measure::stackprofile::StackedDensityProfile
	 *  to be modelled
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
	Modelling_DensityProfile (const std::shared_ptr<cbl::measure::stackprofile::StackedDensityProfile> dataset, const std::string probe, const std::vector<double> zl_bin_edges_for_N_zs={}, const std::vector<double> zs={}, const std::vector<std::vector<double>> N_zs={}, const std::vector<double> purity_background={});
	
	/**
	 *  @brief constuctor for the modelling of
	 *  stacked density excess profile of galaxy cluster. Cosmological units are forced.
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
	Modelling_DensityProfile (const std::shared_ptr<cbl::data::Data> dataset, const std::string probe, const std::vector<double> zl_bin_edges_for_N_zs={}, const std::vector<double> zs={}, const std::vector<std::vector<double>> N_zs={}, const std::vector<double> purity_background={});
	
	/**
	 *  @brief default destructor
	 *  
	 */
	virtual ~Modelling_DensityProfile () = default;

	///@}

	/**
	 *  @name Member functions
	 */
	///@{

	/**
	 * @brief get the member \e m_data_model
	 * @return the container of parameters for cluster density profile
	 * model computation
	 */
	STR_Profile_data_model data_model () { return m_data_model; }

	/**
	 *  @brief Set the data used to construct generic models of
	 *  halo profiles, where the mass is a base parameter (i.e. it
	 *  is not derived from a mass-observable relation). 
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
	 *  @param redshift cluster redshift
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
	void set_data_model (const std::shared_ptr<cosmology::Cosmology> cosmology, const std::string profile_author, const bool two_halo, const std::string halo_def, const double Delta, const double redshift, const double logM_base, const double mass_pivot, const std::string bias_author="Tinker", const std::string method_Pk="EisensteinHu", std::string interp_type="Linear");
	
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
	void set_data_model (const std::shared_ptr<cosmology::Cosmology> cosmology, const std::string profile_author, const bool two_halo, const std::string halo_def, const double Delta, const std::vector<double> redshift_obs, const std::vector<double> mass_proxy_obs, const std::vector<double> Pz_params, const std::vector<double> Pproxy_params, const double redshift_pivot, const double proxy_pivot, const double logM_base, const double mass_pivot, const std::string bias_author="Tinker", const std::string method_Pk="EisensteinHu", std::string interp_type="Linear");
	
	/**
	 *  @brief Set the data used to construct models of
	 *  average halo profiles in bins of redshift and mass proxy, 
	 *  where the mass is derived from a scaling relation and the mean profile
	 *  is weighted over the theoretical halo mass function.
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
	 *  @param redshift_min minimum observed redshift
	 *
	 *  @param redshift_max maximum observed redshift
	 *
	 *  @param mass_proxy_min minimum observed mass proxy
	 *
	 *  @param mass_proxy_max maximum observed mass proxy
	 *
	 *  @param z_error \f$\sigma_{z,0}\f$ in the equation 
	 *  \f$\sigma_z=\sigma_{z,0}(1+z)\f$, where \f$z\f$ is the true redshift and 
	 *  \f$\sigma_z\f$ represents the total statistical uncertainty on redshift.
	 *  In particular, \f$\sigma_z\f$ is the rms of the Gaussian PDF 
	 *  \f$P(z_{\rm ob}|z)\f$
	 *
	 *  @param proxy_rel_error relative error on the mass proxy, namely
	 *  \f$\tilde{\sigma}_{\lambda,0}\f$ in the expression
	 *  \f$\sigma_\lambda=\tilde{\sigma}_{\lambda,0}\lambda\f$, where
	 *  \f$\lambda\f$ is the true mass proxy and \f$\sigma_\lambda\f$
	 *  represents the total statistical uncertainty on the mass proxy.
	 *  In particular, \f$\sigma_\lambda\f$ is the rms of the Gaussian PDF 
	 *  \f$P(\lambda_{\rm ob}|\lambda)\f$.
	 *  For more complex expressions for both the rms and the mean of
	 *  \f$P(\lambda_{\rm ob}|\lambda)\f$, use
	 *  cbl::modelling::densityprofile::Modelling_DensityProfile::set_P_proxy
	 *
	 *  @param redshift_pivot redshift pivot in the scaling relation
	 *
	 *  @param proxy_pivot proxy pivot in the scaling relation
	 *
	 *  @param logM_base base of the mass logarithm
	 *
	 *  @param mass_pivot the mass pivot
	 *
	 *  @param redshift_points_completeness **true** cluster redshift points
	 *  where the completeness is measured. Note that this range should
	 *  be large enough given the uncertainties on the observed redshifts.
	 *  If not provided, the completeness is 1 in each point
	 *
	 *  @param proxy_points_completeness **true** cluster mass proxy points
	 *  where the completeness is measured. Note that this range should
	 *  be large enough given the uncertainties on the observed proxy.
	 *  If not provided, the completeness is 1 in each point
	 *
	 *  @param completeness completeness as a function of redshift (first
	 *  dimension) and mass proxy (second dimension). This quantity is interpolated.
	 *  To check the quality of the interpolation, use
	 *  cbl::modelling::densityprofile::Modelling_DensityProfile::completeness_interp().
	 *  If not provided, the completeness is 1 in each point
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
	 *  @param MF_author author(s) who proposed the mass function model.
	 *  See cbl::cosmology::MassFunction::m_MF_generator
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
	void set_data_model (const std::shared_ptr<cosmology::Cosmology> cosmology, const std::string profile_author, const bool two_halo, const std::string halo_def, const double Delta, const double redshift_min, const double redshift_max, const double mass_proxy_min, const double mass_proxy_max, const double z_error, const double proxy_rel_error, const double redshift_pivot, const double proxy_pivot, const double logM_base, const double mass_pivot, const std::vector<double> redshift_points_completeness={}, const std::vector<double> proxy_points_completeness={}, const std::vector<std::vector<double>> completeness={}, const std::string bias_author="Tinker", const std::string method_Pk="EisensteinHu", std::string interp_type="Linear", const std::string MF_author="Tinker");
	
	/**
	 *  @brief Set the functional forms for the mean and rms of
	 *  the Gaussian PDF \f$P(\lambda_{\rm ob}|\lambda,z)\f$.
	 *
	 *  The mean is expressed as
	 *
	 *  \f$ \mu = \lambda + A_\mu \lambda \exp[- \lambda 
	 *  (B_\mu + C_\mu \, z)]\,, \f$
	 *
	 *  while the rms has the following expression
	 *
	 *  \f$ \sigma = A_\sigma \lambda \exp[- \lambda 
	 *  (B_\sigma + C_\sigma \, z)]\,. \f$
	 *  
	 *  @param A_mu the parameter \f$A_\mu\f$
	 *  
	 *  @param B_mu the parameter \f$B_\mu\f$
	 *  
	 *  @param C_mu the parameter \f$C_\mu\f$
	 *  
	 *  @param A_sigma the parameter \f$A_\sigma\f$
	 *  
	 *  @param B_sigma the parameter \f$B_\sigma\f$
	 *  
	 *  @param C_sigma the parameter \f$C_\sigma\f$
	 *  
	 */
	void set_P_proxy (const double A_mu, const double B_mu, const double C_mu, const double A_sigma, const double B_sigma, const double C_sigma);

	/**
	 *  @brief Set the profile and cosmological parameters used to model the 
	 *  cluster density profile, by assuming a DK14 profile. For details
	 *  on its functional form see cbl::cosmology::HaloProfile.
	 *
	 *  For details on the base profiles, see cbl::cosmology::HaloProfile.
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
	 *  @param Rt_prior prior on the truncation factor \f$F_t\f$ defining the truncation
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
	 *  internally called "AB_fact", affecting the 2-halo term and/or
	 *  the infall term in the DK14 profile. 
	 *  The 2-halo excess surface density
	 *  is expressed as \f$\Sigma_{\rm 2h,\,correct}=
	 *  \Sigma_{\rm 2h}(1+\sigma_{\rm AB})\f$, where \f$\sigma_{\rm AB}\f$ is
	 *  the parameter set through this prior. This is analogous for the DK14
	 *  infall term
	 *
	 *  @param orientation_boost_prior prior on the orientation boost factor,
	 *  internally called "OB_fact". In particular, the profile mass is expressed
	 *  as \f$M_{\rm correct}=M(1+\sigma_{\rm OB})\f$, where \f$\sigma_{\rm OB}\f$ is
	 *  the parameter set through this prior.
	 *
	 *  @param be_prior prior on the normalisation of the DK14 outer profile
	 *
	 *  @param se_prior prior on the slope of the DK14 outer profile
	 *
	 *  @param alpha_0_prior prior on the \f$\alpha_0\f$ parameter in the equation
	 *  \f$ \alpha(\nu) = \alpha_0 + \alpha_\nu\nu^2 \f$, where \f$\nu\f$ 
	 *  is the virial peak height. Following Gao et al. (2008), \f$\alpha_0=0.155\f$.
	 *
	 *  @param alpha_nu_prior prior on the \f$\alpha_\nu\f$ parameter in the equation
	 *  \f$ \alpha(\nu) = \alpha_0 + \alpha_\nu\nu^2 \f$, where \f$\nu\f$ 
	 *  is the virial peak height. Following Gao et al. (2008), \f$\alpha_\nu=0.0095\f$.
	 *
	 *  @param beta_prior prior on the \f$\beta\f$ parameter in the DK14 profile. Following
	 *  DK14, \f$\beta=4\f$.
	 *
	 *  @param gamma_0_prior prior on the \f$\gamma_0\f$ parameter in the 
	 *  expression \f$\gamma=\gamma_0\nu\f$ appearing in the DK14 profile.
	 *  Following DK14, \f$\gamma_0=4\f$.
	 *
	 *  @warning The off-centering is related to stacks of clusters. For
	 *  details, see cbl::cosmology::HaloProfile.
	 *
	 */
	void set_model_DensityProfile_cosmology_DK14 (const std::vector<std::string> cosmo_param, const std::vector<statistics::PriorDistribution> cosmo_prior, const statistics::PriorDistribution Rt_prior, const statistics::PriorDistribution concentration_prior, const statistics::PriorDistribution logM_prior, const statistics::PriorDistribution f_off_prior, const statistics::PriorDistribution sigma_off_prior, const statistics::PriorDistribution anisotropic_boost_prior, const statistics::PriorDistribution orientation_boost_prior, const statistics::PriorDistribution be_prior, const statistics::PriorDistribution se_prior, const statistics::PriorDistribution alpha_0_prior, const statistics::PriorDistribution alpha_nu_prior, const statistics::PriorDistribution beta_prior, const statistics::PriorDistribution gamma_0_prior);
	
	/**
	 *  @brief Set the profile and cosmological parameters used to model the 
	 *  cluster density profile, by assuming a DK14 profile. For details
	 *  on its functional form see cbl::cosmology::HaloProfile.
	 *
	 *  For details on the base profiles, see cbl::cosmology::HaloProfile.
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
	 *  @param Rt_prior prior on the truncation factor \f$F_t\f$ defining the truncation
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
	 *  internally called "AB_fact", affecting the 2-halo term and/or
	 *  the infall term in the DK14 profile. 
	 *  The 2-halo excess surface density
	 *  is expressed as \f$\Sigma_{\rm 2h,\,correct}=
	 *  \Sigma_{\rm 2h}(1+\sigma_{\rm AB})\f$, where \f$\sigma_{\rm AB}\f$ is
	 *  the parameter set through this prior. This is analogous for the DK14
	 *  infall term
	 *
	 *  @param orientation_boost_prior prior on the orientation boost factor,
	 *  internally called "OB_fact". In particular, the profile mass is expressed
	 *  as \f$M_{\rm correct}=M(1+\sigma_{\rm OB})\f$, where \f$\sigma_{\rm OB}\f$ is
	 *  the parameter set through this prior.
	 *
	 *  @param be_prior prior on the normalisation of the DK14 outer profile
	 *
	 *  @param se_prior prior on the slope of the DK14 outer profile
	 *
	 *  @param alpha_0_prior prior on the \f$\alpha_0\f$ parameter in the equation
	 *  \f$ \alpha(\nu) = \alpha_0 + \alpha_\nu\nu^2 \f$, where \f$\nu\f$ 
	 *  is the virial peak height. Following Gao et al. (2008), \f$\alpha_0=0.155\f$.
	 *
	 *  @param alpha_nu_prior prior on the \f$\alpha_\nu\f$ parameter in the equation
	 *  \f$ \alpha(\nu) = \alpha_0 + \alpha_\nu\nu^2 \f$, where \f$\nu\f$ 
	 *  is the virial peak height. Following Gao et al. (2008), \f$\alpha_\nu=0.0095\f$.
	 *
	 *  @param beta_prior prior on the \f$\beta\f$ parameter in the DK14 profile. Following
	 *  DK14, \f$\beta=4\f$.
	 *
	 *  @param gamma_0_prior prior on the \f$\gamma_0\f$ parameter in the 
	 *  expression \f$\gamma=\gamma_0\nu\f$ appearing in the DK14 profile.
	 *  Following DK14, \f$\gamma_0=4\f$.
	 *
	 *  @warning The off-centering is related to stacks of clusters. For
	 *  details, see cbl::cosmology::HaloProfile.
	 *
	 */
	void set_model_DensityProfile_cosmology_DK14 (const std::vector<std::string> cosmo_param, const std::vector<statistics::PriorDistribution> cosmo_prior, const statistics::PriorDistribution Rt_prior, const std::string cM_author, const statistics::PriorDistribution logM_prior, const statistics::PriorDistribution f_off_prior, const statistics::PriorDistribution sigma_off_prior, const statistics::PriorDistribution anisotropic_boost_prior, const statistics::PriorDistribution orientation_boost_prior, const statistics::PriorDistribution be_prior, const statistics::PriorDistribution se_prior, const statistics::PriorDistribution alpha_0_prior, const statistics::PriorDistribution alpha_nu_prior, const statistics::PriorDistribution beta_prior, const statistics::PriorDistribution gamma_0_prior);
	
	/**
	 *  @brief Set the profile and cosmological parameters used to model the 
	 *  cluster density profile. 
	 *
	 *  For details on the base profiles, see cbl::cosmology::HaloProfile.
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
	 *  density is expressed as \f$\Sigma_{\rm 2h,\,correct}=
	 *  \Sigma_{\rm 2h}(1+\sigma_{\rm AB})\f$, where \f$\sigma_{\rm AB}\f$ is
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
	void set_model_DensityProfile_cosmology (const std::vector<std::string> cosmo_param, const std::vector<statistics::PriorDistribution> cosmo_prior, const statistics::PriorDistribution Rt_prior, const statistics::PriorDistribution concentration_prior, const statistics::PriorDistribution logM_prior, const statistics::PriorDistribution f_off_prior, const statistics::PriorDistribution sigma_off_prior, const statistics::PriorDistribution anisotropic_boost_prior, const statistics::PriorDistribution orientation_boost_prior);
	
	/**
	 *  @brief Set the profile and cosmological parameters used to model the 
	 *  cluster density profile. The concentration is a derived parameter,
	 *  computed through a concentration-mass relation.
	 *
	 *  For details on the base profiles, see cbl::cosmology::HaloProfile.
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
	 *  density is expressed as \f$\Sigma_{\rm 2h,\,correct}=
	 *  \Sigma_{\rm 2h}(1+\sigma_{\rm AB})\f$, where \f$\sigma_{\rm AB}\f$ is
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
	void set_model_DensityProfile_cosmology (const std::vector<std::string> cosmo_param, const std::vector<statistics::PriorDistribution> cosmo_prior, const statistics::PriorDistribution Rt_prior, const std::string cM_author, const statistics::PriorDistribution logM_prior, const statistics::PriorDistribution f_off_prior, const statistics::PriorDistribution sigma_off_prior, const statistics::PriorDistribution anisotropic_boost_prior, const statistics::PriorDistribution orientation_boost_prior);
	
	/**
	 *  @brief Set the profile and cosmological parameters used to model the 
	 *  cluster density profile, assuming a redshift and mass proxy evolution
	 *  of the concentration-mass relation and off-centering parameters.
	 *
	 *  For details on the base profiles, see cbl::cosmology::HaloProfile.
	 *
	 *  In particular, with this function the mass
	 *  is not a base parameter. Instead, the cluster mass-mass proxy scaling
	 *  relation is used to derive the mass. The mass is
	 *  expressed e.g. in \f$10^{14}\f$ M\f$_\odot\f$ \f$h^{-1}\f$, if
	 *  mass_pivot set through set_data_model is \f$10^{14}\f$.
	 *
	 *  Specifically, the model is the expectation value of the
	 *  total excess density profile, expressed as:
	 *
	 *  \f$\langle\Delta\Sigma_{\rm tot}(
	 *  R|\Delta\lambda_{\rm ob}^i,\Delta z_{\rm ob}^j)\rangle =
	 *  \mathcal{P}^{ij} \,
	 *  \left[ (1 - f_{\rm off})
	 *  \langle\Delta\Sigma_{\rm cen}(
	 *  R|\Delta\lambda_{\rm ob}^i,\Delta z_{\rm ob}^j)\rangle + 
	 *  f_{\rm off} \langle\Delta\Sigma_{\rm off}(
	 *  R|\Delta\lambda_{\rm ob}^i,\Delta z_{\rm ob}^j)\rangle
	 *  \right] \f$
	 *
	 *  Depending on the cbl::modelling::densityprofile::Modelling_DensityProfile::set_data_model()
	 *  used, the model can be expressed in alternative ways.
	 *
	 *  CASE 1: 
	 *
	 *  \f$\langle\Delta\Sigma_{\rm cen}(
	 *  R|\Delta\lambda_{\rm ob}^i,\Delta z_{\rm ob}^j)\rangle = 
	 *  \frac{1}{N^{ij}} \, \sum_{k=1}^{N^{ij}} \,
	 *  \int_0^\infty\,{\rm d}M\,
	 *  \Delta\Sigma_{\rm cen}(R|M,\hat{z})
	 *  P(M|\hat{\lambda},\hat{z})\,, \f$
	 *
	 *  where \f$\hat{\lambda}\f$ is extracted from \f$P(\lambda|\lambda_{\rm ob},z_{\rm ob})\f$
	 *  and \f$\hat{z}\f$ is extracted from \f$P(z|\lambda_{\rm ob},z_{\rm ob})\f$,
	 *  \f$N^{ij}\f$ is the number of clusters in the stack,
	 *  \f$\mathcal{P}^{ij}\f$ is the purity of the cluster sample.
	 *  If the probe is the reduced shear, the model is multiplied also by
	 *  the purity of the background sample.
	 *
	 *  CASE 2:
	 *
	 *  \f$\langle\Delta\Sigma_{\rm cen}(
	 *  R|\Delta\lambda_{\rm ob}^i,\Delta z_{\rm ob}^j)\rangle = 
	 *  \frac{\Omega}{\langle N(\Delta\lambda_{\rm ob}^i,\Delta z_{\rm ob}^j) \rangle}
	 *  \int_0^\infty{\rm d}z_{\rm tr}\,\frac{{\rm d}^2V}{{\rm d}\Omega{\rm d}z_{\rm tr}}
	 *  \int_0^\infty{\rm d}\lambda_{\rm tr}\,\mathcal{C}_{\rm clu}(\lambda_{\rm tr},z_{\rm tr})
	 *  \frac{{\rm d}n(M,z_{\rm tr})}{{\rm d}M}
	 *  \int_0^\infty{\rm d}M\, \Delta\Sigma_{\rm cen}(R|M,z) P(\lambda|M,z)
	 *  \int_{\Delta\lambda_{\rm ob}}{\rm d}\lambda_{\rm ob}P(\lambda_{\rm ob}|\lambda_{\rm tr},z_{\rm tr})
	 *  \int_{\Delta z_{\rm ob}}{\rm d}z_{\rm ob}P(z_{\rm ob}|z_{\rm tr},\lambda_{\rm tr}) \,. \f$
	 *
	 *  We note that in any case \f$\Sigma_{\rm off}\f$ is expressed as
	 *
	 *  \f$\Sigma_{\rm off}(R) = \int P(R_{\rm s})\, 
	 *  \Sigma (R | R_{\rm s})\, {\rm d} R_{\rm s}\,,\f$
	 *
	 *  where \f$R_{\rm s}\f$ is the projected displacement distance, and
	 *
	 *  \f$ P(R_{\text{s}}) = \frac {R_{\text{s}}} {\sigma_{\text{off}}^2}
	 *  \exp \bigg[-\frac 1 2 \bigg( \frac {R_{\text{s}}} 
	 *  {\sigma_{\text{off}}} \bigg)^2 \bigg]\,,\f$
	 *
	 *  \f$ \Sigma (R | R_{\text{s}}) = \frac 1 {2 \pi} \int_0^{2 \pi} 
	 *  \Sigma_{\text{cen}} \, \bigg( \sqrt {R^2 + R_{\text{s}}^2 + 
	 *  2 R R_{\text{s}} \cos \theta} \bigg) \, {\rm d} \theta\,.\f$
	 *
	 *  \f$P(M|\lambda,z)\f$ is a log-normal whose
	 *  mean is given by the proxy-mass relation and whose rms is given by
	 *  the intrinsic scatter of such relation
	 *  (for details, see cbl::modelling::massobsrel::Modelling\_MassObservableRelation).
	 *
	 *  The standard deviations of the PDFs
	 *  \f$P(z_{{\rm ob},k}|z)\f$ and \f$P(\lambda_{{\rm ob},k}|\lambda)\f$
	 *  are given by the uncertainties set through 
	 *  cbl::modelling::densityprofile::Modelling_DensityProfile::set_data_model.
	 *
	 *  The concentration-mass relation is logarithmic, with logarithmic base
	 *  10, and it is expressed as follows:
	 *
	 *  \f$\log c = c_0 + c_M\, \log_{10}\frac{M}{M_{\rm piv}} +
	 *  c_z\, \log_{10}\frac{1+z}{1+z_{\rm piv}}. \f$
	 *
	 *  The off-centering parameters, namely \f$f_{\rm off}\f$ and
	 *  \f$\sigma_{\rm off}\f$, are expressed as follows:
	 *
	 *  \f$f_{\rm off} = f_{\rm off,0}
	 *  \left(\frac{\bar{\lambda}_{\rm ob}}{\lambda_{\rm piv}}\right)^{f_{\rm off,M}}
	 *  \left(\frac{1+\bar{z}_{\rm ob}}{1+z_{\rm piv}}\right)^{f_{\rm off,z}} \f$,
	 *
	 *  \f$\sigma_{\rm off} = \sigma_{\rm off,0}
	 *  \left(\frac{M}{M_{\rm piv}}\right)^{\sigma_{\rm
	 *  off,M}} \left(\frac{1+z}{1+z_{\rm piv}}\right)^{\sigma_{\rm
	 *  off,z}} \f$.
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
	 *  @param purity_prior prior on the purity of the cluster sample, defined
	 *  between 0 and 1
	 *
	 *  @warning \f$F_t\f$ is used only if a truncated NFW is
	 *  assumed!
	 *
	 *  @warning The off-centering is related to stacks of
	 *  clusters. For details, see cbl::cosmology::HaloProfile.
	 */
	void set_model_DensityProfile_cosmology (const std::vector<std::string> cosmo_param, const std::vector<statistics::PriorDistribution> cosmo_prior, const std::string z_evo, const statistics::PriorDistribution Rt_prior, const statistics::PriorDistribution c0_prior, const statistics::PriorDistribution cM_prior, const statistics::PriorDistribution cz_prior, const statistics::PriorDistribution f_off0_prior, const statistics::PriorDistribution f_offM_prior, const statistics::PriorDistribution f_offz_prior, const statistics::PriorDistribution sigma_off0_prior, const statistics::PriorDistribution sigma_offM_prior, const statistics::PriorDistribution sigma_offz_prior, const statistics::PriorDistribution anisotropic_boost_prior, const statistics::PriorDistribution orientation_boost_prior, const statistics::PriorDistribution alpha_prior, const statistics::PriorDistribution beta_prior, const statistics::PriorDistribution gamma_prior, const statistics::PriorDistribution scatter0_prior, const statistics::PriorDistribution scatterM_prior, const statistics::PriorDistribution scatterM_exponent_prior, const statistics::PriorDistribution scatterz_prior, const statistics::PriorDistribution scatterz_exponent_prior, const statistics::PriorDistribution purity_prior);
	
	/**
	*  @brief Mass expected values divided by the mass pivot,
	*  weighted over the theoretical halo mass function. For details,
	*  see cbl::modelling::densityprofileaux::compute_mass_expected_value
	*
	*  @param z_ob vector of cluster observed redshifts
	*
	*  @param mass_proxy_ob vector of cluster observed mass proxies
	*
	*  @param cosmology cosmological model
	*
	*  @param alpha normalisation of the scaling relation
	*
	*  @param beta slope of the scaling relation
	*
	*  @param gamma scaling relation redshift evolution
	*  parameter
	*
	*  @param scatter0 scatter normalisation
	*
	*  @param scatterM scatter mass proxy dependence parameter
	*
	*  @param scatterM_exp scatter mass proxy exponent parameter
	*
	*  @param scatterz scatter redshift dependence parameter
	*
	*  @param scatterz_exp scatter redshift exponent parameter
	*
	*  @return mass expected values divided by the mass pivot
	*
	**/
	std::vector<double> mass_expected_value (const std::vector<double> z_ob, const std::vector<double> mass_proxy_ob, const std::shared_ptr<cosmology::Cosmology> cosmology, const double alpha, const double beta, const double gamma, const double scatter0, const double scatterM, const double scatterM_exp, const double scatterz, const double scatterz_exp);
	
	/**
	*  @brief interpolated cluster sample completeness
	*
	*  @param z_tr true redshift
	*
	*  @param proxy_tr true mass proxy
	*
	*  @return completeness
	*
	**/
	double completeness_interp(const double z_tr, const double proxy_tr) { return m_data_model.completeness_interp(z_tr, proxy_tr); };

	///@}
	
	/**
	 *  @name Member functions used to get the protected members of the class
	 */
	///@{
	
	/**
	 *  @brief Get the extracted \f$z_{\rm tr}\f$ and \f$\lambda_{\rm tr}\f$
	 *
	 *  @return the extracted \f$z_{\rm tr}\f$ (first elements of the vector)
	 *  and \f$\lambda_{\rm tr}\f$ (second element)
	 */
	std::vector<std::vector<double>> get_extracted_true_observables() { return {m_data_model.redshifts, m_data_model.mass_proxies}; };
	
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
       * @brief Compute the DK14 excess density profile model
       * in all the radial bins
       *
       * @param radius the radius array
       *
       * @param inputs model inputs
       *
       * @param parameter model parameters
       *
       * @return DK14 halo model in each radial bin
       *
       */
      std::vector<double> model_density_DK14 (const std::vector<double> radius, const std::shared_ptr<void> inputs, std::vector<double> &parameter);
      
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
      std::vector<double> model_density_scaling_relation_evolving_concentration_offcentering (const std::vector<double> radius, const std::shared_ptr<void> inputs, std::vector<double> &parameter);
    
    }
  }
}

#endif
