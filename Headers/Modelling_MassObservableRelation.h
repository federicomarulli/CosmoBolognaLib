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
	
	/// the redshift evolution function in the scaling relation
	std::function<double(const double, const double, const std::shared_ptr<void>)> fz;
	
	/// function returning the redshift error (absolute or relative)
	std::function<double(const double, const double)> z_error;
	
	/// function returning the mass proxy error (absolute or relative)
	std::function<double(const double, const double)> proxy_error;
	
	/// Array of redshift values where the scaling relation is evaluated
	std::vector<double> redshift;
	
	/// Redshift pivot
	double redshift_pivot;
	
	/// Proxy pivot
	double proxy_pivot;
	
	/// Mass pivot
	double mass_pivot;
	
	/// logarithmic base
	double log_base;
	
	/// Vector of redshift error values
	std::vector<double> z_eff_err;
	
	/// Vector of proxy error values
	std::vector<double> proxy_eff_err;
	
	/// Vector of number of clusters in the stacked bin
	std::vector<double> Nclusters;
	
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

	/// &Delta;, the overdensity
	double Delta;

	/// isDelta_critical
	bool isDelta_critical;

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
	
	/// check if the number of clusters in the bin is set
	bool m_isSet_Nclusters = false;

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
	
	///@}
	
	
	/**
	 *  @name Member functions used to set the model parameters
	 */
	///@{
	
	/**
	 *  @brief Set the mass pivot
	 *
	 *  @param mass_pivot the mass pivot
	 *
	 */	 
	void set_mass_pivot (const double mass_pivot) { m_data_model.mass_pivot = mass_pivot; }

	/**
	 *  @brief Set the data used to construct the scaling relation,
	 *  written as:
	 * 
	 *  \f$\log M = \alpha + \beta 
	 *  \log (\lambda/\lambda_{\rm piv}) + \gamma \log (f(z)),\f$
	 *
	 *  where \f$\lambda\f$ is the mass proxy.
	 *  
	 *  @param cosmology the cosmological model
	 *
	 *  @param redshift redshift array
	 *
	 *  @param redshift_pivot redshift pivot value
	 *
	 *  @param proxy_pivot proxy or mass pivot value
	 *
	 *  @param log_base base of the mass and proxy logarithms
	 *  
	 */
	void set_data_model (const cbl::cosmology::Cosmology cosmology, const std::vector<double> redshift, const double redshift_pivot, const double proxy_pivot, const double log_base);
	
	/**
	 *  @brief Set the data used to construct the scaling relation,
	 *  written as:
	 * 
	 *  \f$\log M = \int_0^\infty {\rm d}\log M_{\rm tr}\,\,
	 *  \log M_{\rm tr} \, P(\log M_{\rm tr}|\lambda_{\rm eff},z_{\rm eff}) \,,\f$ 
	 *
	 *  where \f$\lambda\f$ is the mass proxy.
	 *  The distirbution \f$P(M| \lambda,z)\f$ is 
	 *  a log-normal whose mean is given by the mass-mass proxy
	 *  relation, i.e. 
	 *
	 *  \f$\log (M/M_{\rm piv}) = \alpha + \beta 
	 *  \log (\lambda/\lambda_{\rm piv}) + \gamma \log (f(z)),\f$
	 *
	 *  and whose standard deviation is given by the intrinsic scatter
	 *  computed in the \f$j\f$-th bin of proxy and redshift, 
	 *  \f$ \sigma_{{\rm intr},j} \f$, expressed as
	 *
	 *  \f$ \sigma_{{\rm intr},j} = \frac{1}{N_{{\rm cl},j}}\left[\sigma_0 
	 *  + \sigma_{\lambda} \log (\lambda/\lambda_{\rm piv})^{e_{\lambda}} 
	 *  + \sigma_z \log (f(z))^{e_z}\right]\,, \f$
	 *
	 *  where \f$N_{{\rm cl},j}\f$ is the number of clusters 
	 *  used for the stacking in the bin.
	 *  
	 *  @param cosmology the cosmological model
	 *
	 *  @param redshift redshift array
	 *
	 *  @param redshift_pivot redshift pivot value
	 *
	 *  @param proxy_pivot proxy or mass pivot value
	 *
	 *  @param log_base base of the mass and proxy logarithms
	 *
	 *  @param Nclusters number of clusters in the bin
	 *  
	 */
	void set_data_model (const cbl::cosmology::Cosmology cosmology, const std::vector<double> redshift, const double redshift_pivot, const double proxy_pivot, const double log_base, const std::vector<double> Nclusters);

	/**
	 *  @brief Set the scaling relation and cosmological parameters,
	 *  where the scaling relation is written as
	 *
	 *  \f$\log M = \alpha + \beta 
	 *  \log (\lambda_{\rm eff}/\lambda_{\rm piv}) 
	 *  + \gamma \log (f(z_{\rm eff})),\,\,\,\,(1)\f$
	 *
	 *  or as 
	 * 
	 *  \f$\log M = \log \left[ \int_0^\infty {\rm d}M_{\rm tr}\,\,
	 *  M_{\rm tr} \, P(M_{\rm tr}|\lambda_{\rm eff},z_{\rm eff})
	 *  \right] \,,\,\,\,\,(2)\f$ 
	 *
	 *  depending on the 
	 *  cbl::modelling::massobsrel::Modelling_MassObservableRelation::set_data_model() used.
	 *  In Eq. (2), the distirbution \f$P(M| \lambda,z)\f$ is 
	 *  a log-normal whose mean is given by the mass-mass proxy
	 *  relation, i.e. 
	 *
	 *  \f$\log (M/M_{\rm piv}) = \alpha + \beta 
	 *  \log (\lambda/\lambda_{\rm piv}) + \gamma \log (f(z)),\f$
	 *
	 *  and whose standard deviation is given by the intrinsic scatter, 
	 *  computed in the \f$j\f$-th bin of proxy and redshift, 
	 *  \f$ \sigma_{{\rm intr},j} \f$, expressed as
	 *
	 *  \f$ \sigma_{{\rm intr},j} = \frac{1}{N_{{\rm cl},j}}\left[\sigma_0 
	 *  + \sigma_{\lambda} \log (\lambda/\lambda_{\rm piv})^{e_{\lambda}} 
	 *  + \sigma_z \log (f(z))^{e_z}\right]\,, \f$
	 *
	 *  where \f$N_{{\rm cl},j}\f$ is the number of clusters 
	 *  used for the stacking in the bin.
	 *
	 *  WARNING: in case of Eq. (1),
	 *  the only way to have a dependency on the intrinsic scatter
	 *  parameters is to define a user-defined likelihood, whose covariance
	 *  depends on the intrinsic scatter. In particular the intrinsic
	 *  scatter, \f$\sigma_{\rm intr}\f$, is expressed as
	 *
	 *  \f$ \sigma_{\rm intr} = \sigma_0 + \sigma_{\lambda} 
	 *  \log (\lambda/\lambda_{\rm piv})^{e_{\lambda}} + \sigma_z \log (f(z))^{e_z}.\f$
	 *
	 *  @param z_evo functional form of the redshift evolution
	 *  function in the scaling relation: "E_z" \f$\rightarrow\f$ 
	 *  \f$ f(z)=E(z)/E(z_{piv}) \f$, "direct" \f$\rightarrow\f$ 
	 *  \f$ f(z)=(1+z)/(1+z_{piv}) \f$
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
	void set_model_MassObservableRelation_cosmology (const std::string z_evo, const std::vector<cbl::cosmology::CosmologicalParameter> cosmo_param, const std::vector<statistics::PriorDistribution> cosmo_prior, const statistics::PriorDistribution alpha_prior, const statistics::PriorDistribution beta_prior, const statistics::PriorDistribution gamma_prior, const statistics::PriorDistribution scatter0_prior, const statistics::PriorDistribution scatterM_prior, const statistics::PriorDistribution scatterM_exponent_prior, const statistics::PriorDistribution scatterz_prior, const statistics::PriorDistribution scatterz_exponent_prior);

	///@}
     };

    
    /**
       * @brief compute the mass - mass proxy scaling relation
       *
       * @param alpha normalisation
       *
       * @param beta slope
       *
       * @param gamma redshift evolution factor
       *
       * @param div_proxy the logarithm of the mass proxy divided by the pivot
       *
       * @param f_z the logarithm of the redshift evolution function
       *
       * @return \f$\log M\f$.
       *
       */
      double scaling_relation(const double alpha, const double beta, const double gamma, const double div_proxy, const double f_z);
      
      /**
       * @brief set the cluster scaling relation model
       *
       * @param proxy the mass proxy (or mass) array
       *
       * @param inputs model inputs
       *
       * @param parameter model parameters
       *
       * @return the scaling relation as a function of
       * the mass proxy in each bin
       *
       */
      std::vector<double> model_scaling_relation (const std::vector<double> proxy, const std::shared_ptr<void> inputs, std::vector<double> &parameter);
      
      /**
       * @brief set the cluster scaling relation model
       *
       * @param proxy the mass proxy (or mass) array
       *
       * @param inputs model inputs
       *
       * @param parameter model parameters
       *
       * @return the scaling relation as a function of
       * the mass proxy in each bin
       *
       */
      std::vector<double> model_scaling_relation_int (const std::vector<double> proxy, const std::shared_ptr<void> inputs, std::vector<double> &parameter);
      
    }
  }
}

#endif
