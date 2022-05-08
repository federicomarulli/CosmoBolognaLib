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
 *  @file Headers/Modelling_MassObservableRelation.h
 *
 *  @brief The class Modelling_MassObservableRelation
 *
 *  This file defines the interface of the class Modelling_MassObservableRelation, used to
 *  model the mass - mass proxy scaling relation for galaxy clusters.
 *
 *  @authors Giorgio Lesci (and Federico Marulli)
 *
 *  @authors giorgio.lesci2@unibo.it (and federico.marulli3@unibo.it)
 */

#ifndef __MODELLINGMORELATION__
#define __MODELLINGMORELATION__

#include "Cosmology.h"
#include "Measure.h"
#include "Modelling.h"
#include "Cluster.h"


// ===================================================================================================


namespace cbl {

  namespace modelling {

    /**
     *  @brief The namespace of the <B> cluster mass - mass proxy scaling relation
     *  modelling </B>
     *  
     *  The \e modelling::massobsrel namespace contains all the functions
     *  and classes to model the mass - mass proxy scaling relation of galaxy clusters
     */
    namespace massobsrel {
      
      /**
       *  @struct STR_MOrelation_data_model
       *  @brief the structure STR_MOrelation_data_model
       *
       *  This structure contains the data used for statistical
       *  analyses of cluster mass - mass proxy scaling relation
       */
      struct STR_MOrelation_data_model {
      
	/// Fiducial cosmology
	std::shared_ptr<cosmology::Cosmology> cosmology;
	
	/// Cosmological parameters
	std::vector<cosmology::CosmologicalParameter> Cpar;
	
	/// Cluster object
	std::shared_ptr<catalogue::Cluster> cluster;
	
	/// pointer to the redshift evolution function in the scaling relation
	double (*fz)(std::vector<double>, std::shared_ptr<void>);
	
	/// pointer to the redshift error
	double (*z_error)(std::vector<double>);
	
	/// pointer to the mass proxy error
	double (*proxy_error)(std::vector<double>);
	
	/// Array of redshift values where the scaling relation is evaluated
	std::vector<double> redshift;
	
	/// Redshift pivot
	double redshift_pivot;
	
	/// Proxy (or mass) pivot
	double proxy_or_mass_pivot;
	
	/// Proxy pivot
	double proxy_pivot;
	
	/// Mass) pivot
	double mass_pivot;
	
	/// logarithmic base
	double log_base;
	
	/// the mass proxy bin edges
	std::vector<std::vector<double>> edges_proxy;
	
	/// the redshift bin edges
	std::vector<double> edges_z;
	
	/// method to compute the dark matter power spectrum
	std::string method_Pk;
	
	/// minimum wave vector module up to which the power spectrum is computed
	double k_min;

	/// maximum wave vector module up to which the power spectrum is computed
	double k_max;

	/// number of steps used to compute the binned dark matter correlation function
	int step;

	/// vector of wave vector modules
	std::vector<double> kk;
	
	/// 0 &rarr; don't normalize the power spectrum; 1 &rarr; normalize the power spectrum
	int norm;
	
	/// true \f$\rightarrow\f$ the output files created by the Boltzmann solver are stored; false \f$\rightarrow\f$ the output files are removed
	bool store_output;
      
	/// output root of the parameter file used to compute the dark matter power spectrum
	std::string output_root;
	
	/// name of the parameter file
	std::string file_par;
	
	/// accuracy of the GSL integration
	double prec;

	/// &Delta;: the overdensity, defined as the mean interior density relative to the background
	double Delta;

	/// isDelta_Vir
	bool isDelta_Vir;

	/// author(s) who proposed the mass function
	std::string model_MF;
	
	/// vector of masses
	std::vector<double> Mass_vector;
	
	/// the survey aperture Area
	double area_rad;
	
	/// true &rarr; sigma8 is a free parameter; false &rarr; sigma8 can be considered a derived parameter
	bool is_sigma8_free;

	/**
	 *  @brief default constructor
	 */
	STR_MOrelation_data_model () = default;
      };
    
      /**
       *  @class Modelling_MassObservableRelation
       *  Modelling_MassObservableRelation.h
       *  "Headers/Modelling_MassObservableRelation.h"
       *
       *  @brief The class Modelling_MassObservableRelation
       *
       *  This file defines the interface of the base class
       *  Modelling_MassObservableRelation, used for modelling the
       *  cluster mass - mass proxy scaling relation.
       *  Cosmological units are forced.
       *
       */
      class Modelling_MassObservableRelation : public Modelling {
      
      protected:

	/// the container of parameters for the model computation
	STR_MOrelation_data_model m_data_model;

      public:

	/**
	 *  @name Constructors/destructors
	 */
	///@{

	/**
	 *  @brief default constuctor
	 *  _MassObservableRelation
	 */
	Modelling_MassObservableRelation () = default;
	
	/**
	 *  @brief constuctor for the modelling of the
	 *  cluster mass - mass proxy relation. Cosmological units are forced
	 *  @param dataset the dataset containing x, data and errors
	 */
	Modelling_MassObservableRelation (const std::shared_ptr<cbl::data::Data> dataset)
	{ m_data = dataset; }
	
	/**
	 *  @brief default destructor
	 *  
	 */
	virtual ~Modelling_MassObservableRelation () = default;

      ///@}

	/**
	 *  @name Member functions used to set the protected members of the class
	 */
	///@{

	/**
	 * @brief get the member \e m_data_model
	 * @return the container of parameters for cluster mass - mass proxy
	 * scaling relation model computation
	 */
	STR_MOrelation_data_model data_model () { return m_data_model; }
	
	
	/**
	 *  @name Member functions used to set the model parameters
	 */
	///@{

	/**
	 *  @brief Set the data used to construct the scaling relation,
	 *  written as:
	 * 
	 *  \f$\log M = \alpha + \beta 
	 *  \log (\lambda/\lambda_{\rm piv}) + \gamma \log (f(z)),\f$
	 *
	 *  or
	 *
	 *  \f$\log \lambda = \alpha + \beta 
	 *  \log (M/M_{\rm piv}) + \gamma \log (f(z)),\f$
	 *
	 *  where \f$\lambda\f$ is the mass proxy.
	 *  
	 *  @param cosmology the cosmological model
	 *
	 *  @param cluster the cluster object
	 *
	 *  @param redshift redshift array
	 *
	 *  @param redshift_pivot redshift pivot value
	 *
	 *  @param proxy_or_mass_pivot proxy or mass pivot value
	 *
	 *  @param log_base base of the mass and proxy logarithms
	 *  
	 */
	void set_data_model (const cbl::cosmology::Cosmology cosmology, const cbl::catalogue::Cluster cluster, const std::vector<double> redshift, const double redshift_pivot, const double proxy_or_mass_pivot, const double log_base);
	
	/**
	 *  @brief set the data used to construct a model of masses derived from
	 *  number counts as a function of a mass proxy, here 
	 *  expressed as \f$\lambda\f$, with the following
	 *  functional form:
	 *  
	 *  \f$ \langle \bar{M}(\Delta{\lambda_{\text{ob},i}},\Delta z_{\text{ob},j})\rangle =
	 *  \frac { \langle M^{\rm tot}(\Delta{\lambda_{\text{ob},i}},\Delta z_{\text{ob},j})\rangle }
	 *  { \langle N(\Delta{\lambda_{\text{ob},i}},\Delta z_{\text{ob},j})\rangle } \f$
	 *
	 *  where
	 *
	 *  \f$ \langle N(\Delta{\lambda_{\text{ob},i}},\Delta z_{\text{ob},j})\rangle
	 *  = \Omega 
	 *  \int_{0}^{\infty} {\rm d} z_{\rm tr}\,\,
	 *  \frac{{\rm d} V}{{\rm d} z_{\rm tr}{\rm d}\Omega}\int_{0}^{\infty} 
	 *  {\rm d} M_{\rm tr} \,\,\frac{{\rm d} n(M_{\rm tr},z_{\rm tr})}{{\rm d} M_{\rm tr}}\,\, 
	 *  \int_{0}^{\infty}{\rm d}\lambda_{\rm tr}\,\,
	 *  P(\lambda_{\rm tr}| M_{\rm tr},z_{\rm tr})\,
	 *  \int_{\Delta z_{\text{ob},j}}{\rm d} z_{\rm ob} 
	 *  \,\,P(z_{\rm ob}|z_{\rm tr})\,
	 *  \int_{\Delta\lambda_{\text{ob},i}}{\rm d} \lambda_{\rm ob} 
	 *  \,\,P(\lambda_{\rm ob}|\lambda_{\rm tr}), \f$
	 *
	 *  \f$ \langle M^{\rm tot}(\Delta{\lambda_{\text{ob},i}},\Delta z_{\text{ob},j})\rangle
	 *  = \Omega 
	 *  \int_{0}^{\infty} {\rm d} z_{\rm tr}\,\,
	 *  \frac{{\rm d} V}{{\rm d} z_{\rm tr}{\rm d}\Omega}\int_{0}^{\infty} 
	 *  {\rm d} M_{\rm tr} \,\,M\frac{{\rm d} n(M_{\rm tr},z_{\rm tr})}{{\rm d} M_{\rm tr}}\,\, 
	 *  \int_{0}^{\infty}{\rm d}\lambda_{\rm tr}\,\,
	 *  P(\lambda_{\rm tr}| M_{\rm tr},z_{\rm tr})\,
	 *  \int_{\Delta z_{\text{ob},j}}{\rm d} z_{\rm ob} 
	 *  \,\,P(z_{\rm ob}|z_{\rm tr})\,
	 *  \int_{\Delta\lambda_{\text{ob},i}}{\rm d} \lambda_{\rm ob} 
	 *  \,\,P(\lambda_{\rm ob}|\lambda_{\rm tr}). \f$
	 *
	 *  The distirbution \f$P(\lambda_{\rm tr}| M_{\rm tr},z_{\rm tr})\f$ is 
	 *  a log-normal whose mean is given by the mass-mass proxy
	 *  relation, i.e. 
	 *
	 *  \f$\log (\lambda/\lambda_{\rm piv}) = \alpha + \beta 
	 *  \log (M/M_{\rm piv}) + \gamma \log (f(z)),\f$
	 *
	 *  and whose standard deviation is given by the intrinsic scatter, \f$ \sigma_{\rm intr} \f$, expressed as
	 *
	 *  \f$ \sigma_{\rm intr} = \sigma_0 + \sigma_{M} 
	 *  \log (M/M_{\rm piv})^{e_{M}} + \sigma_z \log (f(z))^{e_z}. \f$
	 *
	 *  Finally, \f$ P(z_{\rm ob}|z_{\rm tr})\f$ and \f$ P(\lambda_{\rm ob}|\lambda_{\rm tr})\f$
	 *  are Gaussian distributions.
	 *  
	 *  @param cosmology the cosmological model
	 *
	 *  @param cluster object of the class Cluster
	 *
	 *  @param z_edges redshift bin edges
	 *
	 *  @param proxy_edges proxy bin edges
	 *
	 *  @param z_pivot the redshift pivot in the scaling relation
	 *
	 *  @param proxy_pivot the mass proxy pivot in the scaling relation
	 *
	 *  @param mass_pivot the mass pivot; for example, if the scaling relation
	 *  is written as \f$ \log [M/(10^{14}M_\odot /h)] = 
	 *  \alpha+\beta \log(\lambda/\lambda_{piv})+\gamma\log f(z) \f$,
	 *  then mass_pivot = 1.e14
	 *
	 *  @param log_base the base of the logarithm used in the 
	 *  mass-mass proxy scaling relation
	 *
	 *  @param method_Pk method used to compute the power spectrum
	 *  (i.e. the Boltzmann solver); valid choices for method_Pk
	 *  are: CAMB [http://camb.info/], CLASS
	 *  [http://class-code.net/], MPTbreeze-v1
	 *  [http://arxiv.org/abs/1207.1465], EisensteinHu
	 *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
	 *
	 *  @param store_output if true the output files created by
	 *  the Boltzmann solver are stored; if false the output files
	 *  are removed
	 *
	 *  @param norm 0 &rarr; don't normalize the power spectrum; 1
	 *  &rarr; normalize the power spectrum
	 *
	 *  @param Delta \f$\Delta\f$: the overdensity, defined as the
	 *  mean interior density relative to the background
	 *
	 *  @param isDelta_vir \f$\rightarrow\f$ \f$\Delta\f$ is the
	 *  virial overdensity
	 *
	 *  @param model_MF author(s) who proposed the mass function
	 *
	 *  @param area_degrees the area in degrees
	 *
	 *  @param prec the precision
	 *
	 *  
	 */
	void set_data_model (const cosmology::Cosmology cosmology, const catalogue::Cluster cluster, const std::vector<double> z_edges, const std::vector<std::vector<double>> proxy_edges, const double z_pivot, const double proxy_pivot, const double mass_pivot, const double log_base, const std::string method_Pk, const bool store_output=true, const int norm=-1, const double Delta=200., const bool isDelta_vir=true, const std::string model_MF="Tinker", const double area_degrees=par::defaultDouble, const double prec=1.e-4);

	/**
	 *  @brief Set the scaling relation and cosmological parameters,
	 *  where the scaling relation is written, e.g., as:
	 *
	 *  \f$\log M = \alpha + \beta 
	 *  \log (\lambda/\lambda_{\rm piv}) + \gamma \log (f(z)),\f$
	 *
	 *  with the redshift evolution function \f$ f(z)=E(z)/E(z_{piv}) \f$.
	 *
	 *  WARNING: the only way to have a dependency on the intrinsic scatter
	 *  parameters is to define a user-defined likelihood, whose covariance
	 *  depends on the intrinsic scatter. In particular the intrinsic
	 *  scatter, \f$\sigma_{\rm intr}\f$, is expressed as
	 *
	 *  \f$ \sigma_{\rm intr} = \sigma_0 + \sigma_{\lambda} 
	 *  \log (\lambda/\lambda_{\rm piv})^{e_{\lambda}} + \sigma_z \log (f(z))^{e_z}.\f$
	 *
	 *  @param cosmo_param vector of enums containing cosmological
	 *  parameters
	 *
	 *  @param cosmo_prior vector containing the priors for
	 *  the cosmological parameters
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
	void set_model_MassObservableRelation_cosmology (const std::vector<cbl::cosmology::CosmologicalParameter> cosmo_param, const std::vector<statistics::PriorDistribution> cosmo_prior, const statistics::PriorDistribution alpha_prior, const statistics::PriorDistribution beta_prior, const statistics::PriorDistribution gamma_prior, const statistics::PriorDistribution scatter0_prior, const statistics::PriorDistribution scatterM_prior, const statistics::PriorDistribution scatterM_exponent_prior, const statistics::PriorDistribution scatterz_prior, const statistics::PriorDistribution scatterz_exponent_prior);
	
	/**
	 *  @brief Set the scaling relation parameters,
	 *  where the scaling relation is written, e.g., as:
	 *
	 *  \f$\log M = \alpha + \beta 
	 *  \log (\lambda/\lambda_{\rm piv}) + \gamma \log (f(z)),\f$
	 *
	 *  with the redshift evolution function \f$ f(z)=(1+z)/(1+z_{piv}) \f$.
	 *
	 *  WARNING: the only way to have a dependency on the intrinsic scatter
	 *  parameters is to define a user-defined likelihood, whose covariance
	 *  depends on the intrinsic scatter. In particular the intrinsic
	 *  scatter, \f$\sigma_{\rm intr}\f$, is expressed as
	 *
	 *  \f$ \sigma_{\rm intr} = \sigma_0 + \sigma_{\lambda} 
	 *  \log (\lambda/\lambda_{\rm piv})^{e_{\lambda}} + \sigma_z \log (f(z))^{e_z}.\f$
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
	void set_model_MassObservableRelation_cosmology (const statistics::PriorDistribution alpha_prior, const statistics::PriorDistribution beta_prior, const statistics::PriorDistribution gamma_prior, const statistics::PriorDistribution scatter0_prior, const statistics::PriorDistribution scatterM_prior, const statistics::PriorDistribution scatterM_exponent_prior, const statistics::PriorDistribution scatterz_prior, const statistics::PriorDistribution scatterz_exponent_prior);
	
	/**
	 *  @brief Set the cosmological parameters used to model the masses,
	 *  where the logarithmic base, \f$b\f$, is defined by the log_base parameter, from
	 *  number counts as a function of a mass proxy, here written
	 *  as \f$\lambda\f$, with the model expressed as follows:
	 *
	 *  \f$ \langle \log_b\bar{M}(\Delta{\lambda_{\text{ob},i}},\Delta z_{\text{ob},j})\rangle =
	 *  \frac { \langle M^{\rm tot}(\Delta{\lambda_{\text{ob},i}},\Delta z_{\text{ob},j})\rangle }
	 *  { \langle N(\Delta{\lambda_{\text{ob},i}},\Delta z_{\text{ob},j})\rangle } \f$
	 *
	 *  where
	 *
	 *  \f$ \langle N(\Delta{\lambda_{\text{ob},i}},\Delta z_{\text{ob},j})\rangle
	 *  = \Omega 
	 *  \int_{0}^{\infty} {\rm d} z_{\rm tr}\,\,
	 *  \frac{{\rm d} V}{{\rm d} z_{\rm tr}{\rm d}\Omega}\int_{0}^{\infty} 
	 *  {\rm d} M_{\rm tr} \,\,\frac{{\rm d} n(M_{\rm tr},z_{\rm tr})}{{\rm d} M_{\rm tr}}\,\, 
	 *  \int_{0}^{\infty}{\rm d}\lambda_{\rm tr}\,\,
	 *  P(\lambda_{\rm tr}| M_{\rm tr},z_{\rm tr})\,
	 *  \int_{\Delta z_{\text{ob},j}}{\rm d} z_{\rm ob} 
	 *  \,\,P(z_{\rm ob}|z_{\rm tr})\,
	 *  \int_{\Delta\lambda_{\text{ob},i}}{\rm d} \lambda_{\rm ob} 
	 *  \,\,P(\lambda_{\rm ob}|\lambda_{\rm tr}), \f$
	 *
	 *  \f$ \langle M^{\rm tot}(\Delta{\lambda_{\text{ob},i}},\Delta z_{\text{ob},j})\rangle
	 *  = \Omega 
	 *  \int_{0}^{\infty} {\rm d} z_{\rm tr}\,\,
	 *  \frac{{\rm d} V}{{\rm d} z_{\rm tr}{\rm d}\Omega}\int_{0}^{\infty} 
	 *  {\rm d} M_{\rm tr} \,\,M\frac{{\rm d} n(M_{\rm tr},z_{\rm tr})}{{\rm d} M_{\rm tr}}\,\, 
	 *  \int_{0}^{\infty}{\rm d}\lambda_{\rm tr}\,\,
	 *  P(\lambda_{\rm tr}| M_{\rm tr},z_{\rm tr})\,
	 *  \int_{\Delta z_{\text{ob},j}}{\rm d} z_{\rm ob} 
	 *  \,\,P(z_{\rm ob}|z_{\rm tr})\,
	 *  \int_{\Delta\lambda_{\text{ob},i}}{\rm d} \lambda_{\rm ob} 
	 *  \,\,P(\lambda_{\rm ob}|\lambda_{\rm tr}). \f$
	 *
	 *  The distirbution \f$P(\lambda_{\rm tr}| M_{\rm tr},z_{\rm tr})\f$ is 
	 *  a log-normal whose mean is given by the mass-mass proxy
	 *  relation, i.e. 
	 *
	 *  \f$\log (\lambda/\lambda_{\rm piv}) = \alpha + \beta 
	 *  \log (M/M_{\rm piv}) + \gamma \log (f(z)),\f$
	 *
	 *  and whose standard deviation is given by the intrinsic scatter, \f$ \sigma_{\rm intr} \f$, expressed as
	 *
	 *  \f$ \sigma_{\rm intr} = \sigma_0 + \sigma_{M} 
	 *  \log (M/M_{\rm piv})^{e_{M}} + \sigma_z \log (f(z))^{e_z}. \f$
	 *
	 *  Finally, \f$ P(z_{\rm ob}|z_{\rm tr})\f$ and \f$ P(\lambda_{\rm ob}|\lambda_{\rm tr})\f$
	 *  are Gaussian distributions.
	 *
	 *
	 *  @param scalrel_z_evo functional form of the redshift evolution
	 *  function in the scaling relation: "E_z" \f$\rightarrow\f$ 
	 *  \f$ f(z)=E(z)/E(z_{piv}) \f$, "direct" \f$\rightarrow\f$ 
	 *  \f$ f(z)=(1+z)/(1+z_{piv}) \f$
	 *
	 *  @param z_error_type if "absolute", set the absolute error
	 *  on redshift in the model; if "relative", set the relative error
	 *
	 *  @param proxy_error_type if "absolute", set the absolute error
	 *  on the mass proxy in the model; if "relative", set the relative error
	 *
	 *  @param cosmo_param vector of enums containing cosmological
	 *  parameters
	 *
	 *  @param cosmo_param_prior vector containing the priors for
	 *  the cosmological parameters
	 *
	 *  @param alpha_prior prior on the normalization 
	 *  of the mass-mass proxy scaling relation, \f$\alpha\f$
	 *
	 *  @param beta_prior prior on the slope 
	 *  of the mass-mass proxy scaling relation, \f$\beta\f$
	 *
	 *  @param gamma_prior prior on the redshift evolution factor 
	 *  of the mass-mass proxy scaling relation, \f$\gamma\f$
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
	 *  @param z_bias_prior prior on the redshif bias,
	 *  governing the position of the mean of the Gaussian \f$ P(z_{\rm ob}|z_{\rm tr}) \f$.
	 *  For example, let's say that the redshift bias in the data is \f$\Delta z=0.02\f$:
	 *  the data redshifts would be corrected as \f$z_{\rm corr} = z - 0.02 (1+z)\f$, 
	 *  where \f$z\f$ is the not corrected redshift. If you correct the data in this way,
	 *  set this prior as constant and equal to zero. Otherwise you can choose to not
	 *  correct the data and set a non-zero prior on the redshift bias: in this case,
	 *  the mean of \f$ P(z_{\rm ob}|z_{\rm tr}) \f$ will be given by 
	 *  \f$z_{\rm corr} = z + \Delta z (1+z)\f$.
	 *
	 *  @param proxy_bias_prior prior on the mass proxy bias, analogous to z_bias_prior
	 *
	 *  @param z_error_prior prior on the redshift error, governing the distribution
	 *  \f$ P(z_{\rm ob}|z_{\rm tr}) \f$. Depending on the input parameter
	 *  z_error_type, this prior is related to the absolute error or to the relative error
	 *
	 *  @param proxy_error_prior prior on the mass proxy error, governing the distribution
	 *  \f$ P(\lambda_{\rm ob}|\lambda_{\rm tr}) \f$. Depending on the input parameter
	 *  proxy_error_type, this prior is related to the absolute error or to the relative error
	 *  
	 */
	void set_model_MassObservableRelation_cosmology (const std::string scalrel_z_evo, const std::string z_error_type, const std::string proxy_error_type, const std::vector<cbl::cosmology::CosmologicalParameter> cosmo_param, const std::vector<statistics::PriorDistribution> cosmo_param_prior, const statistics::PriorDistribution alpha_prior, const statistics::PriorDistribution beta_prior, const statistics::PriorDistribution gamma_prior, const statistics::PriorDistribution scatter0_prior, const statistics::PriorDistribution scatterM_prior, const statistics::PriorDistribution scatterM_exponent_prior, const statistics::PriorDistribution scatterz_prior, const statistics::PriorDistribution scatterz_exponent_prior, const statistics::PriorDistribution z_bias_prior, const statistics::PriorDistribution proxy_bias_prior, const statistics::PriorDistribution z_error_prior, const statistics::PriorDistribution proxy_error_prior);

	///@}
     };

    
    /**
       * @brief compute the mass - mass proxy scaling relation
       *
       * @param cluster the cluster object
       *
       * @param div_proxy_or_mass the logarithm of the mass proxy (or mass) divided by the pivot
       *
       * @param f_z the logarithm of the redshift evolution function
       *
       * @return \f$\log M\f$ or \f$\log \lambda\f$. E.g. for the mass the relation is 
       *
       * \f$\log M = \alpha+\beta \log\frac{\lambda}{\lambda_{piv}}+\gamma\log f(z) \f$,
       *
       * where \f$\lambda_{piv}\f$ and \f$z_{piv}\f$ are, respectively, the mass proxy and the redshift pivots. 
       * The functional form of the redshift evolution function, \f$ f(z) \f$, depends on the 
       * set_model_MassObservableRelation_cosmology used. It can have the following functional forms:
       *
       * \f$ f(z) = \frac{E(z)}{E(z_{piv})} \equiv \frac{H(z)/H0}{H(z_{piv})/H0}\f$
       *
       * \f$ f(z) = \frac{1+z}{1+z_{piv}}\f$
       *
       */
      double scaling_relation(cbl::catalogue::Cluster cluster, const double div_proxy_or_mass, const double f_z);
      
      /**
       *  @brief function used to model the logarithm of the masses,
       *    where the logarithmic base, \f$b\f$, is defined by the log_base parameter, from
	 *  number counts as a function of a mass proxy, here written
	 *  as \f$\lambda\f$, with the model expressed as follows:
	 *
	 *  \f$ \langle \log_b\bar{M}(\Delta{\lambda_{\text{ob},i}},\Delta z_{\text{ob},j})\rangle =
	 *  \frac { \langle M^{\rm tot}(\Delta{\lambda_{\text{ob},i}},\Delta z_{\text{ob},j})\rangle }
	 *  { \langle N(\Delta{\lambda_{\text{ob},i}},\Delta z_{\text{ob},j})\rangle } \f$
	 *
	 *  where
	 *
	 *  \f$ \langle N(\Delta{\lambda_{\text{ob},i}},\Delta z_{\text{ob},j})\rangle
	 *  = \Omega 
	 *  \int_{0}^{\infty} {\rm d} z_{\rm tr}\,\,
	 *  \frac{{\rm d} V}{{\rm d} z_{\rm tr}{\rm d}\Omega}\int_{0}^{\infty} 
	 *  {\rm d} M_{\rm tr} \,\,\frac{{\rm d} n(M_{\rm tr},z_{\rm tr})}{{\rm d} M_{\rm tr}}\,\, 
	 *  \int_{0}^{\infty}{\rm d}\lambda_{\rm tr}\,\,
	 *  P(\lambda_{\rm tr}| M_{\rm tr},z_{\rm tr})\,
	 *  \int_{\Delta z_{\text{ob},j}}{\rm d} z_{\rm ob} 
	 *  \,\,P(z_{\rm ob}|z_{\rm tr})\,
	 *  \int_{\Delta\lambda_{\text{ob},i}}{\rm d} \lambda_{\rm ob} 
	 *  \,\,P(\lambda_{\rm ob}|\lambda_{\rm tr}), \f$
	 *
	 *  \f$ \langle M^{\rm tot}(\Delta{\lambda_{\text{ob},i}},\Delta z_{\text{ob},j})\rangle
	 *  = \Omega 
	 *  \int_{0}^{\infty} {\rm d} z_{\rm tr}\,\,
	 *  \frac{{\rm d} V}{{\rm d} z_{\rm tr}{\rm d}\Omega}\int_{0}^{\infty} 
	 *  {\rm d} M_{\rm tr} \,\,M\frac{{\rm d} n(M_{\rm tr},z_{\rm tr})}{{\rm d} M_{\rm tr}}\,\, 
	 *  \int_{0}^{\infty}{\rm d}\lambda_{\rm tr}\,\,
	 *  P(\lambda_{\rm tr}| M_{\rm tr},z_{\rm tr})\,
	 *  \int_{\Delta z_{\text{ob},j}}{\rm d} z_{\rm ob} 
	 *  \,\,P(z_{\rm ob}|z_{\rm tr})\,
	 *  \int_{\Delta\lambda_{\text{ob},i}}{\rm d} \lambda_{\rm ob} 
	 *  \,\,P(\lambda_{\rm ob}|\lambda_{\rm tr}). \f$
	 *
	 *  The distirbution \f$P(\lambda_{\rm tr}| M_{\rm tr},z_{\rm tr})\f$ is 
	 *  a log-normal whose mean is given by the mass-mass proxy
	 *  relation, i.e. 
	 *
	 *  \f$\log (\lambda/\lambda_{\rm piv}) = \alpha + \beta 
	 *  \log (M/M_{\rm piv}) + \gamma \log (f(z)),\f$
	 *
	 *  and whose standard deviation is given by the intrinsic scatter, \f$ \sigma_{\rm intr} \f$, expressed as
	 *
	 *  \f$ \sigma_{\rm intr} = \sigma_0 + \sigma_{M} 
	 *  \log (M/M_{\rm piv})^{e_{M}} + \sigma_z \log (f(z))^{e_z}. \f$
	 *
	 *  Finally, \f$ P(z_{\rm ob}|z_{\rm tr})\f$ and \f$ P(\lambda_{\rm ob}|\lambda_{\rm tr})\f$
	 *  are Gaussian distributions.
       *
       * @param fz redshift evolution function in the scaling relation
       *
       * @param z_error relative or absolute error on the redshift
       *
       * @param proxy_error relative or absolute error on the mass proxy
       *
       * @param redshift_min minimum redshift
       *
       * @param redshift_max maximum redshift
       *
       * @param proxy_min minimum mass proxy
       *
       * @param proxy_max maximum mass proxy
       *
       * @param cosmology the cosmology 
       *
       * @param cluster the cluster object
       *
       * @param Area the area in degrees
       *
       * @param model_MF author(s) who proposed the mass function;
       * valid authors are: PS (Press & Schechter), ST (Sheth &
       * Tormen), Jenkins (Jenkins et al. 2001), Warren (Warren et
       * al. 2006), Reed, (Reed et al. 2007), Pan (Pan 2007), ShenH
       * (halo MF by Shen et al. 2006), ShenF (filament MF by Shen et
       * al. 2006), ShenS (sheet MF by Shen et al. 2006), Tinker
       * (Tinker et al. 2008), Crocce (Crocce et al. 2010),
       * Angulo_FOF (FoF MF by Angulo et al. 2012), Angulo_Sub
       * (SUBFIND MF by Angulo et al. 2012), Watson_FOF (FoF MF by
       * Watson et al. 2012), Watson_SOH (Spherical Overdensity halo
       * MF by Watson et al. 2012), Manera (Manera et al. 2010),
       * Bhattacharya (Bhattacharya et al. 2011), Courtin (Courtin et
       * al. 2010), Peacock (by Peacock at al. 2007)
       *
       *  @param store_output if true the output files created by the
       *  Boltzmann solver are stored; if false the output files are
       *  removed
       *
       * @param Delta \f$\Delta\f$: the overdensity, defined as the
       * mean interior density relative to the background
       *
       * @param isDelta_vir \f$\rightarrow\f$ \f$\Delta\f$ is the
       * virial overdensity
       *
       * @param interp_sigmaM interpolating function of \f$
       * \sigma(M)\f$
       *
       * @param interp_DlnsigmaM interpolating function of \f$
       * \mathrm{d} \ln(\sigma(M)) / \mathrm{d} M \f$
       *
       * @param proxy_pivot mass proxy pivot in the scaling relation
       *
       * @param z_pivot redshift pivot in the scaling relation
       *
       * @param mass_pivot mass pivot in the scaling relation
       *
       * @param log_base logarithmic base used in the scaling relation
       *
       * @return values of the mass function as a function of redshift
       * and mass proxy
       */
      double mass_from_counts(double (*fz)(std::vector<double>, std::shared_ptr<void>), double (*z_error)(std::vector<double>), double (*proxy_error)(std::vector<double>), const double redshift_min, const double redshift_max, const double proxy_min, const double proxy_max, cbl::cosmology::Cosmology cosmology, cbl::catalogue::Cluster cluster, const double Area, const std::string model_MF, const bool store_output, const double Delta, const bool isDelta_vir, const cbl::glob::FuncGrid interp_sigmaM, const  cbl::glob::FuncGrid interp_DlnsigmaM, const double proxy_pivot, const double z_pivot, const double mass_pivot, const double log_base);
      
      /**
       * @brief set the cluster scaling relation model
       *
       * @param proxy_or_mass the mass proxy (or mass) array
       *
       * @param inputs model inputs
       *
       * @param parameter model parameters
       *
       * @return the scaling relation as a function of
       * the mass proxy in each bin
       *
       */
      std::vector<double> model_scaling_relation (const std::vector<double> proxy_or_mass, const std::shared_ptr<void> inputs, std::vector<double> &parameter);
      
      /**
       * @brief set the cluster scaling relation model
       *
       * @param proxy the mass proxy array
       *
       * @param inputs model inputs
       *
       * @param parameter model parameters
       *
       * @return the scaling relation as a function of
       * the mass proxy in each bin
       *
       */
      std::vector<double> model_mass_from_counts (const std::vector<double> proxy, const std::shared_ptr<void> inputs, std::vector<double> &parameter);
      
      /**
       * @brief redshift evolution function in the scaling relation, 
       * with the functional form \f$ f(z)=E(z)/E(z_{\rm piv}) = H(z)/H(z_{\rm piv}) \f$
       *
       * @param x vector containing the redhsift and the redshift pivot, in this order
       *
       * @param cosmo the cosmological model
       *
       * @return the value of the redshift evolution function
       *
       */
      double Fz_Ez (const std::vector<double> x, const std::shared_ptr<void> cosmo);
      
      /**
       * @brief redshift evolution function in the scaling relation, 
       * with the functional form \f$ f(z)=(1+z)/(1+z_{\rm piv})\f$
       *
       * @param x vector containing the redhsift and the redshift pivot, in this order
       *
       * @param cosmo the cosmological model
       *
       * @return the value of the redshift evolution function
       *
       */
      double Fz_direct (const std::vector<double> x, const std::shared_ptr<void> cosmo);
      
      /**
       * @brief return the absolute error
       *
       * @param x vector containing the absolute error
       *
       * @return the absolute error
       *
       */
      double ReturnAbsoluteError (const std::vector<double> x);
      
      /**
       * @brief return the absolute error starting from the relative error
       *
       * @param x vector containing the relative error and the measure
       *
       * @return the absolute error
       *
       */
      double AbsoluteFromRelativeError (const std::vector<double> x);
      
    }
  }
}

#endif
