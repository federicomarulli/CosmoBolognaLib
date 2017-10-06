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
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#include "LikelihoodFunction.h"

using namespace cosmobl;


// ============================================================================================


double statistics::LogLikelihood_Gaussian_1D_error (vector<double> &likelihood_parameters, const shared_ptr<void> fixed_parameters)
{
  // ----- extract the parameters ----- 
  shared_ptr<statistics::STR_likelihood_parameters> pp = static_pointer_cast<statistics::STR_likelihood_parameters>(fixed_parameters);

  // ----- compute the model values ----- 

  vector<double> xx; 

  vector<double> computed_model;
  
  if (pp->data->dataType()==data::_1D_data_){
    pp->data->xx(xx);
  }

  else if (pp->data->dataType()==data::_1D_data_extra_) {
    int x_index = pp->x_index[0];
    for (int i=0; i<pp->data->ndata(); i++) // using extra info
      xx.push_back(pp->data->extra_info(x_index, i));
  }

  else ErrorCBL("Error in LogLikelihood_Gaussian_1D_error of LikelihoodFunction.cpp: wrong dataType!");  

  computed_model = pp->model->operator()(xx, likelihood_parameters);
 
  // ----- estimate the Gaussian log-likelihood -----
  
  double LogLikelihood = 0.;
  for (int i=0; i<pp->data->ndata(); i++) 
    LogLikelihood += pow((pp->data->data(i)-computed_model[i])/pp->data->error(i), 2);

  return -0.5*LogLikelihood;
}


// ============================================================================================


double statistics::LogLikelihood_Gaussian_1D_covariance (vector<double> &likelihood_parameters, const shared_ptr<void> fixed_parameters)
{
  // ----- extract the parameters ----- 
  shared_ptr<statistics::STR_likelihood_parameters> pp = static_pointer_cast<statistics::STR_likelihood_parameters>(fixed_parameters);

  // ----- compute the model values ----- 

  vector<double> xx; 

  vector<double> computed_model;
  
  if (pp->data->dataType()==data::_1D_data_){
    pp->data->xx(xx);
  }
  else if (pp->data->dataType()==data::_1D_data_extra_) {
    int x_index = pp->x_index[0];
    for (int i=0; i<pp->data->ndata(); i++) // using extra info
      xx.push_back(pp->data->extra_info(x_index, i));
  }
  else ErrorCBL("Error in LogLikelihood_Gaussian_1D_covariance of LikelihoodFunction.cpp: wrong dataType!");  

  computed_model = pp->model->operator()(xx, likelihood_parameters);
  

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


double statistics::LogLikelihood_Gaussian_2D_error (vector<double> &likelihood_parameters, const shared_ptr<void> fixed_parameters)
{
  // ----- extract the parameters ----- 
  shared_ptr<statistics::STR_likelihood_parameters> pp = static_pointer_cast<statistics::STR_likelihood_parameters>(fixed_parameters);
  
  // ----- compute the model values -----

  vector<double> xx, yy;

  if (pp->data->dataType()==data::_2D_data_){
    pp->data->xx(xx);
    pp->data->yy(yy);
  }

  else if (pp->data->dataType()==data::_2D_data_extra_) {

    int scaleD1_index = pp->x_index[0];
    for (int i=0; i<pp->data->xsize(); ++i) 
      xx.push_back(pp->data->extra_info(scaleD1_index, i));

    int scaleD2_index = pp->x_index[1];
    for (int i=0; i<pp->data->ysize(); ++i) 
      yy.push_back(pp->data->extra_info(scaleD2_index, i));
    
  }

  else ErrorCBL("Error in LogLikelihood_Gaussian_2D_Error of LikelihoodFunction.cpp: wrong dataType!");  

  vector<vector<double>> computed_model = pp->model->operator()(xx, yy, likelihood_parameters);
  
  // ----- estimate the Gaussian log-likelihood -----

  double LogLikelihood = 0.;
  for (int i=0; i<pp->data->xsize(); i++)
    for (int j=0; j<pp->data->ysize(); j++)
      LogLikelihood += pow((pp->data->data(i, j)-computed_model[i][j])/pp->data->error(i, j), 2);

  return -0.5*LogLikelihood;
}

