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
 *  @file Headers/Lib/Modelling_TwoPointCorrelation_deprojected.h
 *
 *  @brief The class Modelling_TwoPointCorrelatoin_deprojected
 *
 *  This file defines the interface of the class
 *  Modelling_TwoPointCorrelation_deprojected, used for modelling deprojected 2pcf
 *
 *  @author Federico Marulli
 *
 *  @author federico.marulli3@unbo.it
 */

#ifndef __MODELLINGDEP__
#define __MODELLINGDEP__


#include "Modelling_TwoPointCorrelation.h"


// ===================================================================================================


namespace cosmobl {
  
  namespace modelling {
    
    /**
     *  @class Modelling_TwoPointCorrelation_deprojected Modelling_TwoPointCorrelation_deprojected.h
     *  "Headers/Lib/Modelling_TwoPointCorrelation_deprojected.h"
     *
     *  @brief The class Modelling_TwoPointCorrelation_deprojected
     *
     *  This file defines the interface of the base class Modelling_TwoPointCorrelation_deprojected,
     *  used for modelling deprojected 2pcf
     *
     */
    class Modelling_TwoPointCorrelation_deprojected : public Modelling_TwoPointCorrelation {

      public:

	/**
	 *  @name Constructors/destructors
	 */
	///@{

	/**
	 *  @brief default constuctor
	 *  @return object of class Modelling_TwoPointCorrelation_deprojected
	 */
	Modelling_TwoPointCorrelation_deprojected () {}

	/**
	 *  @brief default destructor
	 *  @return none
	 */
	virtual ~Modelling_TwoPointCorrelation_deprojected () {}

	/**
	 *  @brief constructor of the Modelling_TwoPointCorrelation_deprojected
	 *  
	 *  @param twop the two-point correlation function to model
	 *
	 *  @return object of type Modelling_TwoPointCorrelation_deprojected
	 */
	Modelling_TwoPointCorrelation_deprojected(const shared_ptr<cosmobl::twopt::TwoPointCorrelation> twop);
	
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
	void set_model_bias (const statistics::Prior bias_prior, const statistics::ParameterType pT_bias=statistics::_free_) override;

	/**
	 * @brief fit the Alcock-Paczynski effect for the two point
	 * correlation function. 
	 * \f$\xi(s)= b^2  \xi_{DM}(\alpha s)\f$.
	 * &xi<SUB>DM</SUB> is computed at the fiducial cosmology.
	 *
	 * @param bias_prior prior for the bias
	 *
	 * @param alpha_prior prior for &alpha;
	 *
	 * @param pT_bias the parameter type of the bias: it can be
	 * either _free_ or _fixed_
	 *
	 * @param pT_alpha the parameter type of &alpha;: it can be
	 * either _free_ or _fixed_
	 *
	 * @return none
	 */
	void set_model_bias_AP_isotropic(const statistics::Prior bias_prior, const statistics::Prior alpha_prior, const statistics::ParameterType pT_bias=statistics::_free_, const statistics::ParameterType pT_alpha=statistics::_free_) override;

	/**
	 * @brief fit the Alcock-Paczynski effect for the two point
	 * correlation function. 
	 * \f$\xi(s)= B^2  \xi_{DM}(\alpha s)\ + A_0 + A_1/s +A_2/s^2\f$.
	 * &xi<SUB>DM</SUB> is computed at the fiducial cosmology.
	 * B, A<SUB>0</SUB>, A<SUB>1</SUB>, A<SUB>2</SUB> are nuisance
	 *
	 * @param alpha_prior the prior for &alpha;
	 *
	 * @param B_prior the prior for B
	 *
	 * @param A0_prior the prior for A0
	 *
	 * @param A1_prior the prior for A1
	 *
	 * @param A2_prior the prior for A2
	 *
	 * @param pT_alpha the parameter type of &alpha;: it can be
	 * either _free_ or _fixed_
	 *
	 * @param pT_B the parameter type of B: it can be either _free_
	 * or _fixed_
	 *
	 * @param pT_A0 the parameter type of A0: it can be either
	 * _free_ or _fixed_
	 *
	 * @param pT_A1 the parameter type of A1: it can be either
	 * _free_ or _fixed_
	 *
	 * @param pT_A2 the parameter type of A2: it can be either
	 * _free_ or _fixed_
	 *
	 * @return none
	 */
	void set_model_AP_isotropic (const statistics::Prior alpha_prior, const statistics::Prior B_prior, const statistics::Prior A0_prior, const statistics::Prior A1_prior, const statistics::Prior A2_prior, const statistics::ParameterType pT_alpha=statistics::_free_, const statistics::ParameterType pT_B=statistics::_free_, const statistics::ParameterType pT_A0=statistics::_free_, const statistics::ParameterType pT_A1=statistics::_free_, const statistics::ParameterType pT_A2=statistics::_free_) override;

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
