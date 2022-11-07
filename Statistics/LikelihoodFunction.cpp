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

/** @file Statistics/LikelihoodFunction.cpp
 *
 *  @brief Implementation of likelihood functions 
 *
 *  This file contains the implementation of some likelihood functions
 *  Likelihood, used for bayesian analyses
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unibo.it, alfonso.veropalumbo@unibo.it
 */

#include "LikelihoodFunction.h"

using namespace std;

using namespace cbl;

// ============================================================================================


cbl::statistics::STR_likelihood_inputs::STR_likelihood_inputs (const std::shared_ptr<data::Data> input_data, const std::shared_ptr<Model> input_model, const vector<size_t> input_x_index, const int input_w_index) : data(input_data), model(input_model) 
{																							
  switch (data->dataType()) {
  case(data::DataType::_1D_):
    xx = data->xx();
    weights1D.resize(data->ndata(), 1.);
    break;

  case (data::DataType::_1D_extra_):

    if (input_x_index.size()==0)
      xx = data->xx();
    else 
      for (int i=0; i<data->ndata(); i++) // using extra info
	xx.push_back(data->extra_info(input_x_index[0], i));
      
    if (input_w_index<0)
      weights1D.resize(data->ndata(), 1.);
    else 
      for (int i=0; i<data->ndata(); i++) // using extra info
	weights1D.push_back(data->extra_info(input_w_index, i));
    break;

  case (data::DataType::_2D_):
    xx = data->xx();
    yy = data->yy();
    if (input_w_index<0) 
      weights2D.resize(data->xsize(), vector<double>(data->ysize(), 1.));
    break;

  case (data::DataType::_2D_extra_):
    if (input_x_index.size()==0) {
      xx = data->xx();
      yy = data->yy();
    }
    else {
      for (int i=0; i<data->ndata(); i++) { // using extra info
	xx.push_back(data->extra_info(input_x_index[0], i));
	yy.push_back(data->extra_info(input_x_index[1], i));
      }
      xx = different_elements(xx);
      yy = different_elements(yy);
    }

    if (input_w_index<0) 
      weights2D.resize(data->xsize(), vector<double>(data->ysize(), 1.));
    else {
      weights2D.resize(data->xsize(), vector<double>(data->ysize(), 0));
      for (int i=0; i<data->xsize(); i++)
	for (int j=0; j<data->ysize(); j++)
	  weights2D[i][j] = data->extra_info(input_w_index, i*data->ysize()+j);
    }
    break;
  default:
    ErrorCBL("wrong dataType!", "STR_likelihood_inputs", "LikelihoodFunction.cpp");
  }
}


// ============================================================================================


double statistics::LogLikelihood_1D_interpolated (std::vector<double> &likelihood_parameter, const std::shared_ptr<void> fixed_parameter)
{
  // ----- extract the parameters ----- 
  shared_ptr<statistics::STR_likelihood_inputs> pp = static_pointer_cast<statistics::STR_likelihood_inputs>(fixed_parameter);

  double ll = pp->interp_function1D->operator()(likelihood_parameter[0]);
  pp->model->parameters()->full_parameter(likelihood_parameter);

  return ll;
}


// ============================================================================================


double statistics::LogLikelihood_2D_interpolated (std::vector<double> &likelihood_parameter, const std::shared_ptr<void> fixed_parameter)
{
  // ----- extract the parameters ----- 
  shared_ptr<statistics::STR_likelihood_inputs> pp = static_pointer_cast<statistics::STR_likelihood_inputs>(fixed_parameter);

  double ll = pp->interp_function2D->operator()(likelihood_parameter[0], likelihood_parameter[1]);
  pp->model->parameters()->full_parameter(likelihood_parameter);

  return ll;
}


// ============================================================================================


double statistics::LogLikelihood_Gaussian_1D_error (std::vector<double> &likelihood_parameter, const std::shared_ptr<void> fixed_parameter)
{
  // ----- extract the parameters ----- 
  shared_ptr<statistics::STR_likelihood_inputs> pp = static_pointer_cast<statistics::STR_likelihood_inputs>(fixed_parameter);

  // ----- compute the model values ----- 

  vector<double> computed_model = pp->model->operator()(pp->xx, likelihood_parameter);

  
  // ----- estimate the Gaussian log-likelihood -----
  
  double LogLikelihood = 0.;
  for (int i=0; i<pp->data->ndata(); i++)  
    LogLikelihood += pow((pp->data->data(i)-computed_model[i])/pp->data->error(i), 2);

  return -0.5*LogLikelihood;
}


// ============================================================================================


double statistics::LogLikelihood_Gaussian_1D_covariance (std::vector<double> &likelihood_parameter, const std::shared_ptr<void> fixed_parameter)
{
  // ----- extract the parameters ----- 
  
  shared_ptr<statistics::STR_likelihood_inputs> pp = static_pointer_cast<statistics::STR_likelihood_inputs>(fixed_parameter);


  // ----- compute the model values ----- 

  vector<double> computed_model = pp->model->operator()(pp->xx, likelihood_parameter);

  // ----- compute the difference between model and data at each bin ----- 
  
  vector<double> diff(pp->data->ndata(), 0);
  for (int i=0; i<pp->data->ndata(); i++)
    diff[i] = pp->data->data(i)-computed_model[i];
  

  // ----- estimate the Gaussian log-likelihood ----- 
  
  double LogLikelihood = 0.;
  for (int i=0; i<pp->data->ndata(); i++)
    for (int j=0; j<pp->data->ndata(); j++)
      LogLikelihood += diff[i]*pp->data->inverse_covariance(i, j)*diff[j];    

  return -0.5*LogLikelihood;
}


// ============================================================================================


double statistics::LogLikelihood_Gaussian_2D_error (std::vector<double> &likelihood_parameter, const std::shared_ptr<void> fixed_parameter)
{
  // ----- extract the parameters ----- 
  shared_ptr<statistics::STR_likelihood_inputs> pp = static_pointer_cast<statistics::STR_likelihood_inputs>(fixed_parameter);
  
  // ----- compute the model values -----
  vector<vector<double>> computed_model = pp->model->operator()(pp->xx, pp->yy, likelihood_parameter);
  
  // ----- estimate the Gaussian log-likelihood -----

  double LogLikelihood = 0.;
  for (int i=0; i<pp->data->xsize(); i++)
    for (int j=0; j<pp->data->ysize(); j++)
      LogLikelihood += pow((pp->data->data(i, j)-computed_model[i][j])/pp->data->error(i, j), 2);

  return -0.5*LogLikelihood;
}


// ============================================================================================


double statistics::LogLikelihood_Poissonian_1D_ (std::vector<double> &likelihood_parameter, const std::shared_ptr<void> fixed_parameter)
{
  // ----- extract the parameters ----- 
  shared_ptr<statistics::STR_likelihood_inputs> pp = static_pointer_cast<statistics::STR_likelihood_inputs>(fixed_parameter);

  // ----- compute the model values ----- 
  vector<double> computed_model = pp->model->operator()(pp->xx, likelihood_parameter);
 
  // ----- estimate the Poissonian log-likelihood -----
  
  double LogLikelihood = 0.;

  for (int i=0; i<pp->data->ndata(); i++)
    LogLikelihood += pp->data->data(i)*cbl::Ln(computed_model[i],1.e-50)-computed_model[i]-gsl_sf_lnfact(int(pp->data->data(i)));

  return LogLikelihood;
}


// ============================================================================================


double statistics::LogLikelihood_Poissonian_2D_ (std::vector<double> &likelihood_parameter, const std::shared_ptr<void> fixed_parameter)
{
  // ----- extract the parameters ----- 
  shared_ptr<statistics::STR_likelihood_inputs> pp = static_pointer_cast<statistics::STR_likelihood_inputs>(fixed_parameter);
  
  // ----- compute the model values -----
  vector<vector<double>> computed_model = pp->model->operator()(pp->xx, pp->yy, likelihood_parameter);
  
  // ----- estimate the Poissonian log-likelihood -----

  double LogLikelihood = 0.;

  for (int i=0; i<pp->data->xsize(); i++)
    for (int j=0; j<pp->data->ysize(); j++)
      LogLikelihood += pp->data->data(i,j)*cbl::Ln(computed_model[i][j],1.e-50)-computed_model[i][j]-gsl_sf_lnfact(int(pp->data->data(i, j)));

  return LogLikelihood;
}

