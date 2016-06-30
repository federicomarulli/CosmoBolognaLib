/********************************************************************
 *  Copyright (C) 2010 by Federico Marulli                          *
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
 *  @file Headers/Lib/Modelling_TwoPointCorrelation_projected.h
 *
 *  @brief The class Modelling_TwoPointCorrelation_projected
 *
 *  This file defines the interface of the class
 *  Modelling_TwoPointCorrelation_projected, used for modelling projected 2pcf
 *
 *  @author Federico Marulli
 *
 *  @author federico.marulli3@unbo.it
 */

#ifndef __MODELLINGPROJ__
#define __MODELLINGPROJ__


#include "Modelling_TwoPointCorrelation.h"


// ===================================================================================================


namespace cosmobl {

  namespace modelling {
    
    /**
     *  @class Modelling_TwoPointCorrelation_projected
     *  Modelling_TwoPointCorrelation_projected.h
     *  "Headers/Lib/Modelling_TwoPointCorrelation_projected.h"
     *
     *  @brief The class Modelling_TwoPointCorrelation_projected
     *
     *  This file defines the interface of the base class
     *  Modelling_TwoPointCorrelation_projected, used to model the
     *  projected of the two-point correlation function
     *
     */
    class Modelling_TwoPointCorrelation_projected : public Modelling_TwoPointCorrelation {

      public:

	/**
	 *  @name Constructors/destructors
	 */
	///@{

	/**
	 *  @brief default constuctor
	 *  @return object of class ModellingTwoPointCorrelation_projected
	 */
	Modelling_TwoPointCorrelation_projected () {}

	/**
	 *  @brief default destructor
	 *  @return none
	 */
	virtual ~Modelling_TwoPointCorrelation_projected () {}

	/**
	 *  @brief constructor of the Modelling_TwoPointCorrelation_projected
	 *  
	 *  @param twop the two-point correlation function to model
	 *
	 *  @return object of type Modelling_TwoPointCorrelation_projected
	 */
	Modelling_TwoPointCorrelation_projected(const shared_ptr<cosmobl::twopt::TwoPointCorrelation> twop);
	
	///@}

	/**
	 * @brief set the fiducial model for dark matter 
	 * two point correlation function
	 *
	 *  @return none
	 */
	void set_fiducial_twop() override;

	/**
	 * @brief fit the bias of the measured
	 * two point correlation function
	 *
	 * @param bias_prior prior for the bias
	 * @param pT_bias the parameter type of the bias: it can be
	 * either _free_ or _fixed_
	 *
	 * @return none
	 */
	void set_model_bias(const statistics::Prior bias_prior, const statistics::ParameterType pT_bias=statistics::_free_) override;

	/**
	 * @brief compute and write the model using the stored 
	 * parameter values
	 *
	 * @param xx vector of point at which the model
	 * is computed
	 * @param dir_model the output directory of the model
	 * @param file_model the name of the file
	 *
	 * @return none
	 */
	void write_model(const vector<double> xx, const string dir_model, const string file_model)
	{ m_model->write_model(xx, dir_model, file_model); }

	/**
	 * @brief compute and write the model using the stored 
	 * parameter values
	 *
	 * @param xx vector of point at which the model
	 * is computed
	 * @param parameters vector of parameters values
	 * at which the model is computed
	 * @param dir_model the output directory of the model
	 * @param file_model the name of the file
	 *
	 * @return none
	 */
	virtual void write_model_parameters(const vector<double> xx, const vector<double> parameters, const string dir_model, const string file_model)
	{ m_model->write_model(xx, parameters, dir_model, file_model); }
    };
  }
}

#endif
