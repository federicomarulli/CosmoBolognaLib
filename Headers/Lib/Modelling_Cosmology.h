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
 *  @file Headers/Lib/Modelling_Cosmology.h
 *
 *  @brief The class Modelling_Cosmology
 *
 *  This file defines the interface of the class Modelling, used for
 *  modelling cosmological measurements
 *
 *  @author Federico Marulli, Alfonso Veropalumbo
 *
 *  @author federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#ifndef __MODELLINGCOSM__
#define __MODELLINGCOSM__


#include "Cosmology.h"
#include "Modelling.h"
#include "ModelFunction_Cosmology.h"


// ===================================================================================================


namespace cosmobl {

  namespace modelling {

    /**
     *  @brief The namespace of the <B> modelling of cosmological data
     *  </B>
     *  
     *  The \e modelling::cosmology namespace contains all the
     *  functions and classes to model cosmological data
     */
    namespace cosmology {

      /**
       *  @class Modelling_Cosmology Modelling_Cosmology.h
       *  "Headers/Lib/Modelling_Cosmology.h"
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
	  shared_ptr<cosmobl::cosmology::Cosmology> m_cosmology;

	  /// map with cosmological parameters
	  std::map<cosmobl::cosmology::CosmoPar, shared_ptr<cosmobl::statistics::Parameter>> m_map_cosmoPar;

	  /// the cosmological measurements
	  vector<string> m_data_type;

	  /// the container of parameters for model computation
	  STR_data_model_cosmology m_data_model;

	public:

	  /**
	   *  @name Constructors/destructors
	   */
	  ///@{

	  /**
	   *  @brief default constuctor
	   *  @return object of class Modelling
	   */
	  Modelling_Cosmology () {}

	  /**
	   *  @brief constuctor
	   *
	   *  @param dataset the dataset containing the input
	   *  measurements to be modelled
	   *
	   *  @param data_type the type of cosmological measurements
	   *  to be modelled
	   *
	   *  @return object of class Modelling
	   */
	  Modelling_Cosmology (const shared_ptr<cosmobl::data::Data> dataset, const vector<string> data_type);

	  /**
	   *  @brief default destructor
	   *  @return none
	   */
	  ~Modelling_Cosmology () {}

	  ///@}

	  /**
	   *  @brief set the fiducial cosmology
	   *
	   *  @param cosmology the fiducial cosmology
	   * 
	   *  @return none
	   */
	  void set_fiducial_cosmology (const cosmobl::cosmology::Cosmology cosmology);

	  /**
	   *  @brief set the fiducial cosmology
	   *
	   * 
	   *  @return none
	   */
	  shared_ptr<cosmobl::cosmology::Cosmology> fiducial_cosmology () {return m_cosmology;}

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
	   * 
	   *  @return none
	   */
	  void set_cosmological_parameters (const vector<cosmobl::cosmology::CosmoPar> cosmoPar_name, const vector<cosmobl::statistics::Prior> cosmoPar_prior, const string distance_prior=par::defaultString, const vector<string> external_dataset={});

	  /**
	   *  @brief return the cosmological parameter
	   *
	   *  @param cosmoPar_name name of the cosmological
	   *  parameter
	   * 
	   *  @return the cosmological parameter
	   */
	  shared_ptr<statistics::Parameter> cosmological_parameter (const cosmobl::cosmology::CosmoPar cosmoPar_name) {return m_map_cosmoPar[cosmoPar_name];}
      };
    }
  }
}

#endif
