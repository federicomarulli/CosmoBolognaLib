/********************************************************************
 *  Copyright (C) 2017 by Federico Marulli                          *
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
 *  @file
 *  Headers/Modelling_ThreePointCorrelation_comoving_reduced.h
 *
 *  @brief The class Modelling_ThreePointCorrelation_comoving_reduced
 *
 *  This file defines the interface of the class
 *  Modelling_ThreePointCorrelation_comoving_reduced, that contains
 *  all the methods to model the connected three-point correlation
 *  function in angular coordinates
 *
 *  @author Federico Marulli, Michele Moresco
 *
 *  @author federico.marulli3@unbo.it, michele.moresco@unibo.it
 */

#ifndef __MODELLINGTHREEPOINTCOMRED__
#define __MODELLINGTHREEPOINTCOMRED__


#include "Modelling_ThreePointCorrelation_comoving_connected.h"
#include "ModelFunction_ThreePointCorrelation_comoving_reduced.h"


// ===================================================================================================


namespace cbl {

  namespace modelling {
    
    namespace threept {
      
      /**
       *  @class Modelling_ThreePointCorrelation_comoving_reduced
       *  Modelling_ThreePointCorrelation_comoving_reduced.h
       *  "Headers/Modelling_ThreePointCorrelation_comoving_reduced.h"
       *
       *  @brief The class
       *  Modelling_ThreePointCorrelation_comoving_reduced
       *
       *  This file defines the interface of the class
       *  Modelling_ThreePointCorrelation_comoving_reduced, that
       *  contains all the methods to model the connected three-point
       *  correlation function in angular coordinates
       *
       */
      class Modelling_ThreePointCorrelation_comoving_reduced : public Modelling_ThreePointCorrelation_comoving_connected {

      public:

	/**
	 *  @name Constructors/destructors
	 */
	///@{

	/**
	 *  @brief default constuctor
	 *  @return object of class
	 *  Modelling_ThreePointCorrelation_comoving_reduced
	 */
	Modelling_ThreePointCorrelation_comoving_reduced () = default;

	/**
	 *  @brief constructor
	 *  
	 *  @param threep the three-point correlation function to
	 *  model
	 *
	 *  @return object of type
	 *  Modelling_ThreePointCorrelation_comoving_reduced
	 */
	Modelling_ThreePointCorrelation_comoving_reduced (const std::shared_ptr<cbl::measure::threept::ThreePointCorrelation> threep)
	  : Modelling_ThreePointCorrelation_comoving_connected(threep) {}

	/**
	 *  @brief constructor
	 *  
	 *  @param threep_dataset the dataset containing the
	 *  three-point correlation function to model
	 *
	 *  @return object of type
	 *  Modelling_ThreePointCorrelation_comoving_reduced
	 */
	Modelling_ThreePointCorrelation_comoving_reduced (const std::shared_ptr<data::Data> threep_dataset)
	  : Modelling_ThreePointCorrelation_comoving_connected() { set_data(threep_dataset); }
	
	/**
	 *  @brief default destructor
	 *  @return none
	 */
	virtual ~Modelling_ThreePointCorrelation_comoving_reduced () = default;
	
	///@}

	/**
	 *  @brief set the parameters used to model the reduced
	 *  three-point correlation function in comoving coordinates
	 *
	 *  @param bias1_prior prior for the parameter \f$\b_1\f$
	 *
	 *  @param bias2_prior prior for the parameter \f$\b_2\f$
	 *
	 *  @return none
	 */
	void set_model_nonlinear_localbias (const statistics::PriorDistribution bias1_prior={}, const statistics::PriorDistribution bias2_prior={});

	/**
	 *  @brief set the parameters used to model the reduced
	 *  three-point correlation function in comoving coordinates
	 *  with non-local contributions
	 *
	 *  @param bias1_prior prior for the parameter \f$\b_1\f$
	 *
	 *  @param bias2_prior prior for the parameter \f$\b_2\f$
	 *
	 *  @param g2_prior prior for the parameter \f$\g_2\f$
	 *
	 *  @return none
	 */
	void set_model_nonlinear_nonlocalbias (const statistics::PriorDistribution bias1_prior={}, const statistics::PriorDistribution bias2_prior={}, const statistics::PriorDistribution g2_prior={});

	/**
	 *  @brief set the parameters used to model the reduced
	 *  three-point correlation function in comoving coordinates
	 *  with non-local contributions
	 *
	 *  @param bias1_prior prior for the parameter \f$\b_1\f$
	 *
	 *  @param bias2_prior prior for the parameter \f$\b_2\f$
	 *
	 *  @param g2_prior prior for the parameter \f$\g_2\f$
	 *
	 *  @param alpha_prior prior for the parameter \f$\alpha\f$
	 *
	 *  @return none
	 */
	void set_model_nonlinear_nonlocalbias_alpha (const statistics::PriorDistribution bias1_prior={}, const statistics::PriorDistribution bias2_prior={}, const statistics::PriorDistribution g2_prior={}, const statistics::PriorDistribution alpha_prior={});

      };
    }
  }
}

#endif
