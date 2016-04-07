/********************************************************************
 *  Copyright (C) 2014 by Federico Marulli and Alfonso Veropalumbo  *
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
 *  @file Headers/Lib/Chi2.h
 *
 *  @brief The class Chi2 
 *
 *  This file defines the interface of the class Chi2, used for
 *  &chi;&sup2; analyses
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#ifndef __CHI2__
#define __CHI2__

#include "Model.h"


// ============================================================================================


namespace cosmobl {

  /**
   *  @brief The namespace of functions and classes used for statistical
   *  analysis
   *  
   * The \e statistic namespace contains all the functions and classes
   * used for statistical analyis
   */
  namespace statistics {

    struct STR_params{
      shared_ptr<Data> data;
      shared_ptr<Model> model;

      STR_params(shared_ptr<Data> _data, shared_ptr<Model> _model) :
	data(_data), model(_model) {}
    };

    typedef function<double(double , shared_ptr<void>)> chi2_1par;
    typedef function<double(vector<double> , shared_ptr<void> )> chi2_npar;


    double chi2_1D_model_1par (double model_parameters, const shared_ptr<void> fixed_parameters);
    double chi2_1D_error_1par (double model_parameters, const shared_ptr<void> fixed_parameters);
    double chi2_1D_covariance_1par (double model_parameters, const shared_ptr<void> fixed_parameters);
    double chi2_2D_error_1par (double model_parameters, const shared_ptr<void> fixed_parameters);

    double chi2_1D_model_npar (vector<double> model_parameters, const shared_ptr<void> fixed_parameters);
    double chi2_1D_error_npar (vector<double> model_parameters, const shared_ptr<void> fixed_parameters);
    double chi2_1D_covariance_npar (vector<double> model_parameters, const shared_ptr<void> fixed_parameters);
    double chi2_2D_error_npar (vector<double> model_parameters, const shared_ptr<void> fixed_parameters);

    /**
     *  @class Chi2 Chi2.h "Headers/Lib/Chi2.h"
     *
     *  @brief The class Chi2
     *
     *  This class is used to handle objects of type &chi;&sup2;. It is
     *  used for all kind of &chi;&sup2; analyses, such as &chi;&sup2; minimisation
     *  and the estimation of confidence contours
     */

    class Chi2
    {
      protected:

	/**
	 *  @name Data and Model
	 */
	///@{

	/// data containers
	shared_ptr<Data> m_data;

	/// model to test
	shared_ptr<Model> m_model; 

	///@}

      public:

	/**
	 *  @name Constructors/destructors
	 */
	///@{

	/**
	 *  @brief default constructor
	 *  @return object of class Chi2
	 */
	Chi2 () {}

	/**
	 *  @brief constructor
	 *  @param data pointers to the data container
	 *  @param model pointers to the model 
	 *  @return object of class Chi2
	 */
	Chi2 (const shared_ptr<Data> data, const shared_ptr<Model> model)
	  : m_data(data), m_model(model) {}

	/**
	 *  @brief default destructor
	 *  @return none
	 */
	~Chi2 () {}

	///@}

	/**
	 *  @brief funciton that minimize chi square whit one free parameter,
	 * find best fit parameters and store them in model
	 *  @param parameter starting value of the parameter 
	 *  @param type chi2 function to be used
	 *  @param dim dimension of data; can be 1 or 2
	 *  @param max_iter maximum number of iteration 
	 *  @param min minumum value for minima finding
	 *  @param max maximum value for minima finding
	 *  @return none
	 */
	void minimize (double parameter, const string type="model", const int dim=1, const unsigned int max_iter=100, const double min=-1.e30, const double max=1.e30);

	/**
	 *  @brief funciton that minimize chi square, find best fit
	 *  parameters and store them in model
	 *  @param parameters vector containing parameters starting values
	 *  @param type chi2 function to be used
	 *  @param dim dimension of data; can be 1 or 2
	 *  @param max_iter maximum number of iteration 
	 *  @param tol the tolerance for minima finding convergence
	 *  @return none
	 */   
	void minimize (const vector<double> parameters, const string type="model", const int dim=1, const unsigned int max_iter=100, const double tol=1.e-6); 

	/**
	 *  @brief funciton that minimize chi square, find best fit
	 *  parameters and store them in model
	 *  @param parameter starting value of the parameter 
	 *  @param f function of type chi2_1par
	 *  @param max_iter maximum number of iteration 
	 *  @param min minumum value for minima finding
	 *  @param max maximum value for minima finding
	 *  @return none
	 */
	void minimize (const double parameter, const chi2_1par f, const unsigned int max_iter=100, const double min=-1.e30, const double max=1.e30);

	/**
	 *  @brief funciton that minimize chi square, find best fit
	 *  parameters and store them in model
	 *  @param parameters vector containing parameters starting values
	 *  @param f function of type chi2_npar
	 *  @param max_iter maximum number of iteration 
	 *  @param tol the tolerance for minima finding convergence  
	 *  @return none
	 */
	void minimize (const vector<double> parameters, const chi2_npar f, const unsigned int max_iter=100, const double tol=1.e-6); 

    };

  }
}

#endif
