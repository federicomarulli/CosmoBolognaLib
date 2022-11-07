/********************************************************************
 *  Copyright (C) 2016 by Federico Marulli and Alfonso Veropalumbo  *
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
 *  @file Headers/Modelling_Cosmology.h
 *
 *  @brief The class Modelling_Cosmology
 *
 *  This file defines the interface of the class Modelling, used for
 *  modelling cosmological measurements
 *
 *  @author Federico Marulli, Alfonso Veropalumbo
 *
 *  @author federico.marulli3@unibo.it, alfonso.veropalumbo@unibo.it
 */

#ifndef __MODELLINGCOSM__
#define __MODELLINGCOSM__


#include "Cosmology.h"
#include "Modelling.h"
#include "ModelFunction_Cosmology.h"


// ===================================================================================================


namespace cbl {

  namespace modelling {

    /**
     *  @brief The namespace of the <B> modelling of cosmological data
     *  </B>
     *  
     *  The \e modelling::cosmology namespace contains all the
     *  functions and classes to model cosmological data
     */
    namespace cosmo {

      /**
       *  @class Modelling_Cosmology Modelling_Cosmology.h
       *  "Headers/Modelling_Cosmology.h"
       *
       *  @brief The class Modelling_Cosmology
       *
       *  This file defines the interface of the base class Modelling_Cosmology,
       *  used for modelling cosmological measurements
       *
       */
      class Modelling_Cosmology : public Modelling {

      protected:

	/// cosmology
	std::shared_ptr<cbl::cosmology::Cosmology> m_cosmology;

	/// map with cosmological parameters
	std::map<cbl::cosmology::CosmologicalParameter, int> m_map_cosmoPar;

	/// the cosmological measurements
	std::vector<std::string> m_data_type;

	/// the container of parameters for model computation
	STR_data_model_cosmology m_data_model;

      public:

	/**
	 *  @name Constructors/destructors
	 */
	///@{
	
	/**
	 *  @brief default constuctor
	 */
	Modelling_Cosmology () = default;

	/**
	 *  @brief constuctor
	 *
	 *  @param dataset the dataset containing the input
	 *  measurements to be modelled
	 *
	 *  @param data_type the type of cosmological measurements
	 *  to be modelled
	 */
	Modelling_Cosmology (const std::shared_ptr<cbl::data::Data> dataset, const std::vector<std::string> data_type);

	/**
	 *  @brief default destructor
	 */
	~Modelling_Cosmology () = default;

	///@}
	
	/**
	 *  @brief set the fiducial cosmology
	 *
	 *  @param cosmology the fiducial cosmology
	 */
	void set_fiducial_cosmology (const cbl::cosmology::Cosmology cosmology);

	/**
	 *  @brief set the fiducial cosmology
	 *
	 *  @return pointer to object of class Cosmology
	 */
	std::shared_ptr<cbl::cosmology::Cosmology> fiducial_cosmology ()
	{ return m_cosmology; }

	/**
	 *  @brief set the cosmological parameters
	 *
	 *  @param cosmoPar_name vector of cosmological
	 *  parameter names
	 *
	 *  @param cosmoPar_prior vector of cosmological
	 *  parameter priors
	 *
	 *  @param distance_prior the CMB distance prior name
	 *
	 *  @param external_dataset external dataset file name
	 */
	void set_cosmological_parameters (const std::vector<cbl::cosmology::CosmologicalParameter> cosmoPar_name, const std::vector<cbl::statistics::PriorDistribution> cosmoPar_prior, const std::string distance_prior=par::defaultString, const std::vector<std::string> external_dataset={});

      };
    }
  }
}

#endif
