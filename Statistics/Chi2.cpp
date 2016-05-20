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
 *  @file Chi2.cpp
 *
 *  @brief Methods of the class Chi2 
 *
 *  This file contains the implementation of the methods of the class
 *  Chi2, used for &chi;&sup2; analyses
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#include "Chi2.h"
using namespace cosmobl;


// ============================================================================================


double cosmobl::statistics::chi2_1D_model_1par (double model_parameter, const shared_ptr<void> fixed_parameters)
{  
  shared_ptr<Data> data = static_pointer_cast<STR_params>(fixed_parameters)->data;
  shared_ptr<Model> model = static_pointer_cast<STR_params>(fixed_parameters)->model;

  model_parameter = model->update_parameter(model_parameter);
  vector<double> pars = {model_parameter};

  vector<double> computed_model(data->ndata(),0);
  for (int i=data->x_down(); i< data->x_up(); i++)
    computed_model[i]=model->operator()(data->xx(i),pars);

  double c2=0;
  for (int i=data->x_down(); i< data->x_up(); i++)
    c2+=pow((data->fx(i)-computed_model[i])/computed_model[i],2);

  return c2;

}


// ============================================================================================


double cosmobl::statistics::chi2_1D_model_npar (vector<double> model_parameters, const shared_ptr<void> fixed_parameters)
{
  shared_ptr<Data> data = static_pointer_cast<STR_params>(fixed_parameters)->data;
  shared_ptr<Model> model = static_pointer_cast<STR_params>(fixed_parameters)->model;

  model_parameters = model->update_parameters(model_parameters);

  vector<double> computed_model(data->ndata(),0);
  for (int i=data->x_down(); i< data->x_up(); i++)
    computed_model[i]=model->operator()(data->xx(i),model_parameters);

  double c2=0;
  for (int i=data->x_down(); i< data->x_up(); i++)
    c2+=pow((data->fx(i)-computed_model[i])/computed_model[i],2);

  return c2;
}

// ============================================================================================


double cosmobl::statistics::chi2_1D_error_1par (double model_parameter, const shared_ptr<void> fixed_parameters)
{  
  shared_ptr<Data> data = static_pointer_cast<STR_params>(fixed_parameters)->data;
  shared_ptr<Model> model = static_pointer_cast<STR_params>(fixed_parameters)->model;

  model_parameter = model->update_parameter(model_parameter);
  vector<double> pars = {model_parameter};

  double c2 = 0;
  
  for (int i=data->x_down(); i<data->x_up(); i++)
    c2 += pow((data->fx(i)-model->operator()(data->xx(i),pars))/data->error_fx(i),2);

  return c2;
}


// ============================================================================================


double cosmobl::statistics::chi2_1D_error_npar (vector<double> model_parameters, const shared_ptr<void> fixed_parameters)
{
  shared_ptr<Data> data = static_pointer_cast<STR_params>(fixed_parameters)->data;
  shared_ptr<Model> model = static_pointer_cast<STR_params>(fixed_parameters)->model;

  model_parameters = model->update_parameters(model_parameters);

  double c2 = 0;
  
  for (int i=data->x_down(); i<data->x_up(); i++)
    c2 += pow((data->fx(i)-model->operator()(data->xx(i),model_parameters))/data->error_fx(i),2);

  return c2;
}


// ============================================================================================


double cosmobl::statistics::chi2_1D_covariance_1par (double model_parameter, const shared_ptr<void> fixed_parameters)
{
  shared_ptr<Data> data = static_pointer_cast<STR_params>(fixed_parameters)->data;
  shared_ptr<Model> model = static_pointer_cast<STR_params>(fixed_parameters)->model;

  model_parameter = model->update_parameter(model_parameter);
  vector<double> pars = {model_parameter};

  vector<double> computed_model(data->ndata(),0);
  for (int i=data->x_down(); i< data->x_up(); i++)
    computed_model[i]=model->operator()(data->xx(i),pars);

  double c2 = 0;
  
  for (int i=data->x_down(); i<data->x_up(); i++)
    for (int j=data->x_down(); j<data->x_up(); j++)
      c2 += (data->fx(i)-computed_model[i])*data->inverse_covariance_fx(i,j)*(data->fx(j)-computed_model[j]);    

  return c2;
}


// ============================================================================================


double cosmobl::statistics::chi2_1D_covariance_npar (vector<double> model_parameters, const shared_ptr<void> fixed_parameters)
{
  shared_ptr<Data> data = static_pointer_cast<STR_params>(fixed_parameters)->data;
  shared_ptr<Model> model = static_pointer_cast<STR_params>(fixed_parameters)->model;

  model_parameters = model->update_parameters(model_parameters);

  vector<double> computed_model(data->ndata(),0);
  for (int i=data->x_down(); i< data->x_up(); i++)
    computed_model[i]=model->operator()(data->xx(i),model_parameters);

  double c2 = 0;
  
  for (int i=data->x_down(); i<data->x_up(); i++)
    for (int j=data->x_down(); j<data->x_up(); j++)
      c2 += (data->fx(i)-computed_model[i])*data->inverse_covariance_fx(i,j)*(data->fx(j)-computed_model[j]);
    
  return c2;
}


// ============================================================================================


double cosmobl::statistics::chi2_2D_error_1par (double model_parameter, const shared_ptr<void> fixed_parameters)
{  
  shared_ptr<Data> data = static_pointer_cast<STR_params>(fixed_parameters)->data;
  shared_ptr<Model> model = static_pointer_cast<STR_params>(fixed_parameters)->model;

  vector<double> pars = {model_parameter};
  pars = model->update_parameters(pars);

  double c2 = 0;
  
  for (int i=data->x_down(); i<data->x_up(); i++)
    for (int j=data->y_down(); j<data->y_up(); j++)
      c2 += pow((data->fxy(i,j)-model->operator()(data->xx(i),data->yy(j),pars))/data->error_fxy(i,j),2);
    
  return c2;
}


// ============================================================================================


double cosmobl::statistics::chi2_2D_error_npar (vector<double> model_parameters, const shared_ptr<void> fixed_parameters)
{
  shared_ptr<Data> data = static_pointer_cast<STR_params>(fixed_parameters)->data;
  shared_ptr<Model> model = static_pointer_cast<STR_params>(fixed_parameters)->model;

  model_parameters = model->update_parameters(model_parameters);

  double c2 = 0;
  for (int i=data->x_down(); i<data->x_up(); i++)
    for (int j=data->y_down(); j<data->y_up(); j++)
      c2 += pow((data->fxy(i,j)-model->operator()(data->xx(i),data->yy(j),model_parameters))/data->error_fxy(i,j),2);
 
  return c2;
}


// ============================================================================================


void cosmobl::statistics::Chi2::minimize (const double parameter, const string type, const int dim, const unsigned int max_iter, const double min, const double max)
{
  if (type == "error" && dim ==1)
    minimize(parameter, chi2_1D_error_1par, max_iter, min, max);
  
  else if (type =="model" && dim ==1)
    minimize(parameter, chi2_1D_model_1par, max_iter, min, max);
  
  else if (type =="covariance" && dim ==1)
    minimize(parameter, chi2_1D_covariance_1par, max_iter, min, max);
  
  else if (type =="error" && dim ==2)
    minimize(parameter, chi2_2D_error_1par, max_iter, min, max);
  
  else
    ErrorMsg("Error in minimize of Chi2, no such type or dim");
}

// ============================================================================================


void cosmobl::statistics::Chi2::minimize (const vector<double> parameters, const string type, const int dim, const unsigned int max_iter, const double tol)
{
  if (type == "error" && dim ==1)
    minimize(parameters, chi2_1D_error_npar, max_iter, tol);
  
  else if (type =="model" && dim ==1)
    minimize(parameters, chi2_1D_model_npar, max_iter, tol);
  
  else if (type =="covariance" && dim ==1)
    minimize(parameters, chi2_1D_covariance_npar, max_iter, tol);
  
  else if (type =="error" && dim ==2)
    minimize(parameters, chi2_2D_error_npar, max_iter, tol);
  
  else
    ErrorMsg("Error in minimize of Chi2, no such type or dim");
}


// ============================================================================================


void cosmobl::statistics::Chi2::minimize (double parameter, const chi2_1par f, const unsigned int max_iter, const double min, const double max)
{
  unsigned int npar = m_model->npar_eff();

  if (npar!=1)
    ErrorMsg("Error in minimize of Chi2, wrong number of parameters");

  double Parameter = parameter;
  
  if (!m_model->parameter(0)->prior()->isIncluded(parameter)) { 
    string msg = "Warning, starting value for parameter '"+m_model->parameter(0)->name()+"' is out of range, it will be changed";
    WarningMsg(msg);
    Parameter = 0.5*(m_model->parameter(0)->prior()->xmax()+m_model->parameter(0)->prior()->xmin());
  }

  double Min = (min>-1.e29) ? min : m_model->parameter(0)->prior()->xmin();
  double Max = (max<1.e29) ? max : m_model->parameter(0)->prior()->xmax();

  auto fixed_parameters = make_shared<STR_params>(STR_params(m_data,m_model));

  GSLfunction_1D_1 func(f, fixed_parameters);
  func.minimize(Parameter, max_iter, Min, Max);
  cout << "Done" << endl;
}


// ============================================================================================


void cosmobl::statistics::Chi2::minimize (const vector<double> parameters, const chi2_npar f, const unsigned int max_iter, const double tol)
{
  unsigned int npar = m_model->npar_eff();
  if(npar==0)
    ErrorMsg("Error in minimize of Chi2, there is no parameter to vary");

  auto fixed_parameters = make_shared<STR_params>(STR_params(m_data,m_model));

  vector<double> pars ;
  if(parameters.size() == npar){
    pars = parameters;
  }
  else if(parameters.size() == m_model->npar()){
    for(unsigned int i=0;i<m_model->npar();i++){
      if(!m_model->parameter(i)->isFreezed())
        pars.push_back(parameters[i]);
    }
  }
  else{ErrorMsg("Error in minimize of Chi2, unrecognized number of parameters");}

  vector<double> step_size(npar,1);
  int nn = 0;

  for (unsigned int i=0; i<m_model->npar(); i++) {
    if (!m_model->parameter(i)->isFreezed()) {
      
      if(!m_model->parameter(i)->prior()->isIncluded(pars[nn])) { 
        string msg = "Warning, starting value for parameter '"+m_model->parameter(i)->name()+"' is out of range, it will be changed";
        WarningMsg(msg);
        pars[nn] = 0.5*(m_model->parameter(i)->prior()->xmax()+m_model->parameter(i)->prior()->xmin());
      }
    
      if ( m_model->parameter(i)->interval_size()<1.e29)
        step_size[nn] = 1.e-1*m_model->parameter(i)->interval_size();
     
      nn ++;
    }
  }

  GSLfunction_nD_1 func(npar, f, fixed_parameters);
  func.minimize(pars, step_size, max_iter, tol);
  cout << "Done" << endl;
}
