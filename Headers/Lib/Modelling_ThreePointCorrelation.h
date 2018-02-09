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
 *  @file Headers/Lib/Modelling_ThreePointCorrelation.h
 *
 *  @brief The class Modelling_ThreePointCorrelation
 *
 *  This file defines the interface of the class
 *  Modelling_ThreePointCorrelation, used to model three-point
 *  correlation functions of any kind
 *
 *  @author Federico Marulli, Michele Moresco
 *
 *  @author federico.marulli3@unbo.it, michele.moresco@unibo.it
 */

#ifndef __MODELLINGTHREEPOINT__
#define __MODELLINGTHREEPOINT__


#include "ThreePointCorrelation.h"
#include "Modelling.h"
#include "ModelFunction_ThreePointCorrelation.h"


// ===================================================================================================


namespace cosmobl {

  namespace modelling {

    /**
     *  @brief The namespace of the <B> three-point correlation
     *  function modelling </B>
     *  
     *  The \e modelling::threept namespace contains all the functions
     *  and classes to model the three-point correlation function
     */
    namespace threept {
    
      /**
       *  @class Modelling_ThreePointCorrelation
       *  Modelling_ThreePointCorrelation.h
       *  "Headers/Lib/Modelling_ThreePointCorrelation.h"
       *
       *  @brief The class Modelling_ThreePointCorrelation
       *
       *  This file defines the interface of the base class
       *  Modelling_ThreePointCorrelation, used for modelling any kind
       *  of three-point correlation function measurements
       *
       */
      class Modelling_ThreePointCorrelation : public Modelling 
      {
      
      protected:
	
	/// the three-point correlation function type
	measure::threept::ThreePType m_threePType;

	/// the container of parameters for three-point correlation function model computation
	modelling::threept::STR_data_model_threept m_data_model;


      public:

	/**
	 *  @name Constructors/destructors
	 */
	///@{

	/**
	 *  @brief default constuctor
	 *  @return object of class ModellingThreePointCorrelation
	 */
	Modelling_ThreePointCorrelation () = default;
	
	/**
	 *  @brief constuctor
	 *  @param threep the three-point correlation function to model
	 *  @return object of class Modelling_ThreePointCorrelation
	 */
	Modelling_ThreePointCorrelation (const shared_ptr<cosmobl::measure::threept::ThreePointCorrelation> threep)
	  { m_data = threep->dataset(); }
	
	/**
	 *  @brief default destructor
	 *  @return none
	 */
	virtual ~Modelling_ThreePointCorrelation () = default;

	/**
	 *  @brief static factory used to construct modelling of 
	 *  three-point correlation functions of any type
	 *
	 *  @param threep the three-point correlation function to model
	 *
	 *  @return a pointer to an object of class 
	 *  Modelling_ThreePointCorrelation of a given type
	 */
	static shared_ptr<Modelling_ThreePointCorrelation> Create (const shared_ptr<measure::threept::ThreePointCorrelation> threep);

	/**
	 *  @brief static factory used to construct modelling of 
	 *  three-point correlation functions of any type
	 *
	 *  @param threePType type of the three-point correlation
	 *  function
	 *
	 *  @param threept_dataset the dataset containing the
	 *  three-point correlation function to model
	 *
	 *  @return a pointer to an object of class 
	 *  Modelling_ThreePointCorrelation of a given type
	 */
	static shared_ptr<Modelling_ThreePointCorrelation> Create (const measure::threept::ThreePType threePType, const shared_ptr<data::Data> threept_dataset);

	///@}

	
	/**
	 *  @brief return the type of correlation function
	 *  @return the type of correlation function
	 */
	measure::threept::ThreePType threePType () { return m_threePType; }


	// ============================================================================================

	
	/**
	 *  @brief set the data model for the three-point correlation
	 *  function 
	 *
	 *  @param Q_DM vector contaning the DM reduced three-point
	 *  correlation function
	 *
	 *  @return none
	 */
	void set_data_model (const vector<double> Q_DM);

	/**
	 *  @brief set the data model for the three-point correlation
	 *  function with non-local contributions
	 *
	 *  @param cosmology
	 *
	 *  @param r1
	 *
	 *  @param r2
	 *
	 *  @param theta
	 *
	 *  @param model
	 *
	 *  @param kk
	 *
	 *  @param Pk_DM
	 *
	 *  @return none
	 */
	void set_data_Q_nonlocal (const cosmology::Cosmology cosmology, const double r1, const double r2, const vector<double> theta, const string model, const vector<double> kk, const vector<double> Pk_DM);

	/**
	 *  @name Member functions used to write the outputs
	 */
	///@{
      
	/**
	 *  @brief compute and write the model
	 *
	 *  @param dir the output directory
	 *
	 *  @param file the output file
	 *
	 *  @param xx vector of points at which the model is computed
	 *
	 *  @param parameter vector containing the input parameters
	 *  used to compute the model; if it is not provided, best-fit
	 *  parameters will be used
	 *
	 *  @return none
	 */
	void write_model (const string dir, const string file, const vector<double> xx, const vector<double> parameter) const
	{ vector<double> par = parameter; m_model->write_model(dir, file, xx, par); }
      
	///@}
	
      };
    }
  }
}

#endif
