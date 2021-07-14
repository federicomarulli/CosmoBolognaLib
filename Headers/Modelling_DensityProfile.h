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
     *  and classes to model the surface density profile of galaxy clusters
     */
    namespace densityprofile {
      
      /**
       *  @struct STR_Profile_data_model
       *  @brief the structure STR_Profile_data_model
       *
       *  This structure contains the data used for statistical
       *  analyses of cluster surface density profiles
       */
      struct STR_Profile_data_model {
      
	/// Fiducial cosmology
	std::shared_ptr<cosmology::Cosmology> cosmology;
	
	/// Cosmological parameters
	std::vector<cosmology::CosmologicalParameter> Cpar;
	
	/// Cluster object
	std::shared_ptr<catalogue::Cluster> cluster;
	
	/// Redshift
	double redshift;
	
	/// The density contrast with respect to the critical density (i.e. if equal to 200, \f$M_{200}\f$ is considered)
	double contrast;
	
	/// Truncation factor
	double trunc_fact;	
	
	/// Base of the mass logarithm
	double logM_base;

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
       *  cluster surface density profile measurements, 
       *  i.e. \f$\Delta\Sigma(r)\f$ [\f$h\f$ M\f$_\odot\f$/pc\f$^2\f$].
       *  Cosmological units are forced.
       *
       */
      class Modelling_DensityProfile : public Modelling {
      
      protected:

	/// the container of parameters for the density model computation
	STR_Profile_data_model m_data_model;
	
	/// if true, consider the 2-halo contribution
	bool m_2halo;


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
	 *  stacked density profiles of galaxy clusters. Cosmological units are forced
	 *  @param profile the object of stacked density profile to model
	 *  @param _2halo if true, compute the 2-halo contribution
	 */
	Modelling_DensityProfile (const std::shared_ptr<cbl::measure::stackprofile::StackedDensityProfile> profile, const bool _2halo=false)
	{ m_data = profile->dataset(); m_2halo = _2halo; }
	
	/**
	 *  @brief constuctor for the modelling of
	 *  a generic cluster density profile. Cosmological units are forced
	 *  @param dataset cluster profile dataset
	 *  @param _2halo if true, compute the 2-halo contribution
	 */
	Modelling_DensityProfile (const std::shared_ptr<cbl::data::Data> dataset, const bool _2halo=false)
	{ m_data = dataset; m_2halo = _2halo; }
	
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
	 */
	void set_data_model (const cbl::cosmology::Cosmology cosmology, const cbl::catalogue::Cluster, const double redshift, const double contrast, const double trunc_fact, const double logM_base);

	/**
	 *  @brief Set the profile and cosmological parameters used to model the 
	 *  cluster density profile. 
	 *
	 *  The 1-halo term has the following functional form, including the
	 *  contribution of centered and off-centered populations
	 *  of galaxy clusters (see Bellagamba et al. 2019):
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
	 *  expressed in \f$10^{14}\f$ M\f$_\odot\f$ \f$h^{-1}\f$). The base of
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

	///@}
     };

    
    /**
       * @brief compute the density profile as a function
       * of mass proxy and redshift
       *
       * @param cosmology the cosmology
       *
       * @param cluster the cluster object 
       *
       * @param radius radius where the profile is computed (Mpc/h)
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
      double nfw_1halo(cosmology::Cosmology cosmology, catalogue::Cluster cluster, const double radius, const double redshift, const double contrast, const double trunc_fact);
      
      /**
       * @brief set the truncated NFW density profile model
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
      std::vector<double> model_nfw_1halo (const std::vector<double> radius, const std::shared_ptr<void> inputs, std::vector<double> &parameter);
    
    }
  }
}

#endif
