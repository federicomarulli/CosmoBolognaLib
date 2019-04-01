/*******************************************************************
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

/** @file Statistics/Likelihood.cpp
 *
 *  @brief Methods of the class Likelihood 
 *
 *  This file contains the implementation of the methods of the class
 *  Likelihood
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#include "Chi2.h"

using namespace std;

using namespace cbl;


// ============================================================================================


double cbl::statistics::Chi2::operator() (vector<double> &pp) const
{
  if (m_likelihood_type==statistics::LikelihoodType::_NotSet_)
    ErrorCBL("Error in cbl::statistics::Chi2::operator() of Chi2.cpp: you should provide a dataset!");

  return -2*m_log_likelihood_function(pp, m_likelihood_inputs);
}


// ============================================================================================


cbl::statistics::Chi2::Chi2 (const std::shared_ptr<data::Data> data, const std::shared_ptr<Model> model, const bool use_covariance, const std::vector<size_t> x_index, const int w_index) 
{
  set_data(data);
  set_model(model);
  m_x_index = x_index;
  m_w_index = w_index;

  if (use_covariance) 
    set_function(LikelihoodType::_Gaussian_Covariance_, m_x_index, m_w_index);
  else
    set_function(LikelihoodType::_Gaussian_Error_, m_x_index, m_w_index);

  m_likelihood_inputs = make_shared<STR_likelihood_inputs>(STR_likelihood_inputs(m_data, m_model, m_x_index, m_w_index));
}


// ============================================================================================


void cbl::statistics::Chi2::minimize (const vector<double> start, vector<vector<double>> parameter_limits, const unsigned int max_iter, const double tol, const double epsilon)
{
  if (m_likelihood_type==statistics::LikelihoodType::_NotSet_)
    ErrorCBL("Error in cbl::statistics::Chi2::minimize() of Chi2.cpp: a dataset should be provided!");

  unsigned int npar = m_model->parameters()->nparameters_free();

  if (npar==0)
    ErrorCBL("Error in cbl::statistics::Chi2::minimize() of Likelihood.cpp: there is no parameter free to vary!");
  if (start.size() != npar && parameter_limits.size()!=0)
    ErrorCBL("Error in cbl::statistics::Chi2::minimize() of Likelihood.cpp: wrong size for the vector of starting parameters!");
  if (parameter_limits.size() != npar && parameter_limits.size()!=0)
    ErrorCBL("Error in cbl::statistics::Chi2::minimize() of Likelihood.cpp: wrong size for the vector of parameter limits!");


  function<double(vector<double> &)> func = [this](vector<double> & pp) { 
    return this->operator()(pp); 
  };

  coutCBL << "Minimizing..." << endl;
  vector<double> result = cbl::wrapper::gsl::GSL_minimize_nD(func, start, parameter_limits, max_iter, tol, epsilon);
  coutCBL << "Done!" << endl << endl;

  m_model->parameters()->set_bestfit_values(result);

  m_model->parameters()->write_bestfit_info();
  coutCBL << "Chi2 = " << this->operator()(result) << endl << endl;
}
