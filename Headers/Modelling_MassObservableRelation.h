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
	
	/// Array of redshift values where the scaling relation is evaluated
	std::vector<double> redshift;
	
	/// Redshift pivot
	double redshift_pivot;
	
	/// Proxy (or mass) pivot
	double proxy_or_mass_pivot;
	
	/// logarithmic base
	double log_base;

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
       * @brief set the cluster scaling relation model, for which
       * the redshift evolution function is \f$E(z)/E(z_{piv})\f$
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
      std::vector<double> model_E_z (const std::vector<double> proxy_or_mass, const std::shared_ptr<void> inputs, std::vector<double> &parameter); 
      
      /**
       * @brief set the cluster scaling relation model, for which
       * the redshift evolution function is \f$(1+z)/(1+z_{piv})\f$
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
      std::vector<double> model_direct_z (const std::vector<double> proxy_or_mass, const std::shared_ptr<void> inputs, std::vector<double> &parameter);
      
    }
  }
}

#endif
