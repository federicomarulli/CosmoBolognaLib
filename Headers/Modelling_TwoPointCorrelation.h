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
 *  @file Headers/Modelling_TwoPointCorrelation.h
 *
 *  @brief The class Modelling_TwoPointCorrelation
 *
 *  This file defines the interface of the class Modelling, used to
 *  model two-point correlation functions of any kind
 *
 *  @author Federico Marulli, Alfonso Veropalumbo
 *
 *  @author federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#ifndef __MODELLINGTWOPOINT__
#define __MODELLINGTWOPOINT__


#include "TwoPointCorrelation.h"
#include "Modelling.h"
#include "ModelFunction_TwoPointCorrelation.h"


// ===================================================================================================


namespace cbl {

  namespace modelling {

    /**
     *  @brief The namespace of the <B> two-point correlation function
     *  modelling </B>
     *  
     *  The \e modelling::twopt namespace contains all the functions
     *  and classes to model the two-point correlation function
     */
    namespace twopt {
    
      /**
       *  @class Modelling_TwoPointCorrelation
       *  Modelling_TwoPointCorrelation.h
       *  "Headers/Modelling_TwoPointCorrelation.h"
       *
       *  @brief The class Modelling_TwoPointCorrelation
       *
       *  This file defines the interface of the base class
       *  Modelling_TwoPointCorrelation, used for modelling any kind of
       *  two-point correlation function measurements
       *
       */
      class Modelling_TwoPointCorrelation
      {
      
      protected:
	
	/// the two-point correlation function type
	measure::twopt::TwoPType m_twoPType;

	/// the container of parameters for two-point correlation function model computation
	std::shared_ptr<modelling::twopt::STR_data_model> m_data_model;

	
      public:

	/**
	 *  @name Constructors/destructors
	 */
	///@{

	/**
	 *  @brief default constuctor
	 *  @return object of class ModellingTwoPointCorrelation
	 */
	Modelling_TwoPointCorrelation () = default;
	
	/**
	 *  @brief default destructor
	 *  @return none
	 */
	virtual ~Modelling_TwoPointCorrelation () = default;

	/**
	 *  @brief static factory used to construct modelling of 
	 *  two-point correlation functions of any type
	 *
	 *  @param twop the two-point correlation function to model
	 *
	 *  @return a pointer to an object of class 
	 *  Modelling_TwoPointCorrelation of a given type
	 */
	static std::shared_ptr<Modelling_TwoPointCorrelation> Create (const std::shared_ptr<measure::twopt::TwoPointCorrelation> twop);

	/**
	 *  @brief static factory used to construct modelling of 
	 *  two-point correlation functions of any type
	 *
	 *  @param twoPType type of the two-point correlation function
	 *
	 *  @param twop_dataset the dataset containing the two-point
	 *  correlation function to model
	 *
	 *  @return a pointer to an object of class 
	 *  Modelling_TwoPointCorrelation of a given type
	 */
	static std::shared_ptr<Modelling_TwoPointCorrelation> Create (const measure::twopt::TwoPType twoPType, const std::shared_ptr<data::Data> twop_dataset);

	///@}


	/**
	 *  @name Member functions used to get the protected members of the class
	 */
	///@{
	
	/**
	 * @brief get the member \e m_twoPType
	 * @return the two-point correlation function type
	 */
	measure::twopt::TwoPType twoPType () { return m_twoPType; }

	/**
	 * @brief get the member \e m_data_model
	 * @return the container of parameters for two-point
	 * correlation function model computation
	 */
	std::shared_ptr<modelling::twopt::STR_data_model> data_model () { return m_data_model; }
	
	///@}

      };
    }
  }
}

#endif
