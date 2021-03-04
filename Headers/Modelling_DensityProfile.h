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
     *  @enum ProfileParameter
     *  @brief the density profile parameters
     */
    enum class ProfileParameter {
  
  	/// Concentration parameter
        _concentration_,
        
        /// Base 10 logarithm of the cluster mass (expressed in \f$10^{14}\f$ \f$M_\odot\f$)
        _LogM_,
        
        /// Fraction of halos that belong to the miscentred population
        _f_off_,
        
        /// rms of the distribution of the misplacement of the halos on the plane of the sky
        _sigma_off_
  
      };    
      
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
	
	/// Pointer to the modelling object
	std::shared_ptr<modelling::Modelling> profile;
	
	/// Redshift
	double redshift;
	
	/// The density contrast with respect to the critical density (i.e. if equal to 200, \f$M_{200}\f$ is considered)
	double contrast;
	
	/// Truncation factor
	double trunc_fact;	

	/// Cosmological parameters
	std::vector<cosmology::CosmologicalParameter> Cpar;
	
	/// Profile parameters
	std::vector<ProfileParameter> Ppar;
	
	/// Profile paramater name strings
	std::vector<std::string> Ppar_string;

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
       *
       */
      class Modelling_DensityProfile : public Modelling {
      
      protected:

	/// the container of parameters for the density model computation
	STR_Profile_data_model m_data_model;
	
	/// if true, consider the 2-halo contribution
	bool m_2halo;
	
	/// Concentration parameter
        double m_concentration;
        
        /// Base 10 logarithm of the cluster mass (expressed in \f$10^{14}\f$ \f$M_\odot\f$)
        double m_LogM;
        
        /// Fraction of halos that belong to the miscentred population
        double m_f_off;
        
        /// rms of the distribution of the misplacement of the halos on the plane of the sky
        double m_sigma_off;	


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
	 *  @brief constuctor for the modelling of the
	 *  stacked density profile of clusters. Cosmological units are forced
	 *  @param profile the object of stacked density profile to model
	 *  @param _2halo if true, compute the 2-halo contribution
	 */
	Modelling_DensityProfile (const std::shared_ptr<cbl::measure::stackprofile::StackedDensityProfile> profile, const bool _2halo=false)
	{ m_data = profile->dataset(); m_2halo = _2halo; m_concentration = 5; m_LogM = 1.e14; m_f_off = 0.2; m_sigma_off = 0.2; }
	
	/**
	 *  @brief default destructor
	 *  
	 */
	virtual ~Modelling_DensityProfile () = default;

      /**
       *  @name Functions to get the private members of the class
       */
      ///@{
      
      /**
       *  @brief get the private member specified by the enum ProfileParameter
       *  
       *  @param parameter the profile parameter
       *
       *  @return value of the profile parameter
       */
      double value (const ProfileParameter parameter) const;
      
      /**
       *  @brief get the private member Modelling_DensityProfile::m_concentration
       *
       *  @return galaxy cluster concentration
       */
      double concentration () const { return m_concentration; };
      
      /**
       *  @brief get the private member Modelling_DensityProfile::m_LogM
       *
       *  @return base 10 log of the galaxy cluster mass
       */
      double LogM () const { return m_LogM; };
      
      /**
       *  @brief get the private member Modelling_DensityProfile::m_f_off
       *
       *  @return Fraction of halos that belong to the miscentred population
       */
      double f_off () const { return m_f_off; };
      
      /**
       *  @brief get the private member Modelling_DensityProfile::m_sigma_off
       *
       *  @return rms of the distribution of the misplacement of the halos on the plane of the sky
       */
      double sigma_off () const { return m_sigma_off; };

      ///@}

	///@}

	/**
	 *  @name Member functions used to set the protected members of the class
	 */
	///@{
	
       /**
       *  @brief set the value of the concentration
       *
       *  @param concentration
       */
      void set_concentration (const double concentration) {
	m_concentration = concentration;
      };
      
      /**
       *  @brief set the value of the base 10 log of the cluster mass
       *
       *  @param LogM
       */
      void set_LogM (const double LogM) {
	m_LogM = LogM;
      };
      
      /**
       *  @brief set the value of m_f_off
       *
       *  @param f_off
       */
      void set_f_off (const double f_off) {
	m_f_off = f_off;
      };
      
      /**
       *  @brief set the value of m_sigma_off
       *
       *  @param sigma_off
       */
      void set_sigma_off (const double sigma_off) {
	m_sigma_off = sigma_off;
      };
      
      /**
       *  @brief set the value of one profile paramter
       *  
       *  @param parameter profile parameter to set
       *  @param value the new value for the parameter 
       */
      void set_parameter (const ProfileParameter parameter, const double value);
      
      /**
       *  @brief set the value of one profile parameter 
       *  providing its name string
       *  
       *  @param parameter profile parameter name to set
       *  @param value the new value for the parameter 
       */
      void set_parameter_from_string (const std::string parameter, const double value) override;
      
      /**
       *  @brief get the value of a profile parameter 
       *  providing its name string
       *  
       *  @param parameter parameter name to get
       *  @return the parameter value
       */
      double get_parameter_from_string (const std::string parameter) const override;

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
	 *  number counts
         *  
	 *  @param cosmology the cosmological model
	 *
	 *  @param redshift the cluster redshift
	 *
	 *  @param contrast The density contrast with respect 
	 *  to the critical density (i.e. if equal to 200, \f$M_200\f$ is considered)
	 *
	 *  @param trunc_fact truncation factor
	 *  
	 */
	void set_data_model (const cbl::cosmology::Cosmology cosmology, const double redshift, const double contrast, const double trunc_fact);

	/**
	 *  @brief set the profile and cosmological parameters used to model the 
	 *  cluster density profile function
	 *
	 *  @param cosmo_param vector of enums containing cosmological
	 *  parameters
	 *
	 *  @param profile_param vector of enums containing the profile
	 *  parameters
	 *
	 *  @param cosmo_prior vector containing the priors for
	 *  the cosmological parameters
	 *
	 *  @param profile_prior vector containing the priors for
	 *  the profile parameters
	 *
	 */
	void set_model_DensityProfile_cosmology (const std::vector<cbl::cosmology::CosmologicalParameter> cosmo_param={}, std::vector<ProfileParameter> profile_param={}, const std::vector<statistics::PriorDistribution> cosmo_prior={}, const std::vector<statistics::PriorDistribution> profile_prior={});

	///@}
     };
      
     /**
     * @brief return a vector containing the ProfileParameter
     * names
     *
     * @return a vector containing the ProfileParameter names
     *
     */
      inline std::vector<std::string> ProfileParameterNames () { return {"concentration", "LogM", "f_off", "sigma_off"}; }      
     
     /**
     * @brief cast an enum of type ProfileParameter from its
     * index
     *
     * @param profileParameterIndex the profileParameter
     * index
     *
     * @return object of class ProfileParameter
     *
     */
    inline ProfileParameter ProfileParameterCast (const int profileParameterIndex) { return castFromValue<ProfileParameter>(profileParameterIndex); }

    /**
     * @brief cast an enum of type ProfileParameter from its name
     *
     * @param profileParameterName the profileParameter name
     *
     * @return object of class ProfileParameter
     *
     */
    inline ProfileParameter ProfileParameterCast (const std::string profileParameterName) { return castFromName<ProfileParameter>(profileParameterName, ProfileParameterNames()); }
    
    /**
     * @brief cast an enum of type ProfileParameter from indeces
     *
     * @param profileParameterIndeces the profileParameter
     * indeces
     *
     * @return vector of objects of class ProfileParameter
     *
     */
    inline std::vector<ProfileParameter> ProfileParameterCast (const std::vector<int> profileParameterIndeces) { return castFromValues<ProfileParameter>(profileParameterIndeces); } 

    /**
     * @brief cast an enum of type ProfileParameter from thier
     * names
     *
     * @param profileParameterNames the profileParameter
     * names
     *
     * @return vector of objects of class ProfileParameter
     *
     */
    inline std::vector<ProfileParameter> ProfileParameterCast (const std::vector<std::string> profileParameterNames) { return castFromNames<ProfileParameter>(profileParameterNames, ProfileParameterNames()); }

    /**
     *  @brief name of the profile parameter
     *
     *  @param parameter the profile parameter
     *
     *  @return a std::string containing the name of the profile
     *  parameter provided in input
     *
     */
    std::string ProfileParameter_name (const ProfileParameter parameter);
    
    /**
       * @brief compute the density profile as a function
       * of mass proxy and redshift
       *
       * @param cosmology the cosmology 
       *
       * @param radius the radius where the profile is computed (Mpc/h)
       *
       * @param redshift the cluster redshift
       *
       *  @param contrast density contrast with respect 
       *  to the critical density (i.e. if equal to 200, M_200 is considered)
       *
       * @param trunc_fact truncation factor
       *
       * @param concentration cluster concentration
       * 
       * @param mass cluster mass (in solar masses)
       *
       * @param f_off Fraction of halos that belong to the miscentred population
       *
       * @param sigma_off rms of the distribution of the misplacement 
       * of the halos on the plane of the sky
       *
       * @return the model for \f$\Delta\Sigma(r)\f$, expressed in units
       * of \f$h\f$ M\f$_\odot\f$/pc\f$^2\f$
       *
       */
      double nfw_1halo(cosmology::Cosmology cosmology, const double radius, const double redshift, const double contrast, const double trunc_fact, const double concentration, const double mass, const double f_off, const double sigma_off);
      
      /**
       * @brief set the density profile model
       *
       * @param radius the radius array
       *
       * @param inputs model inputs
       *
       * @param parameter model parameters
       *
       * @return the 1-halo model in each radius bin
       *
       */
      std::vector<double> nfw_1halo_allBins (const std::vector<double> radius, const std::shared_ptr<void> inputs, std::vector<double> &parameter);
    
    }
  }
}

#endif
