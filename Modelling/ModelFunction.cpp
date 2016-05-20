/*******************************************************************
 *  Copyright (C) 2010 by Federico Marulli and Alfonso Veropalumbo *
 *  federico.marulli3@unibo.it                                     *
 *                                                                 *
 *  This program is free software; you can redistribute it and/or  *
 *  modify it under the terms of the GNU General Public License as *
 *  published by the Free Software Foundation; either version 2 of *
 *  the License, or (at your option) any later version.            *
 *                                                                 *
 *  This program is distributed in the hope that it will be useful,*
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the  *
 *  GNU General Public License for more details.                   *
 *                                                                 *
 *  You should have received a copy of the GNU General Public      *
 *  License along with this program; if not, write to the Free     *
 *  Software Foundation, Inc.,                                     *
 *  59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.      *
 *******************************************************************/

/** 
 *  @file Modelling/ModelFunction.cpp
 *
 *  @brief implementation of the function for data modelling
 *
 *  This file contains the implementation of a set of functions 
 *  to model the data
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#include "ModelFunction.h"


// ============================================================================================


double cosmobl::glob::wp_bias_cosmology(double r, shared_ptr<void> parameters, vector<double> model_parameters){

  shared_ptr<STR_twop_model> pp = static_pointer_cast<STR_twop_model>(parameters);

  double bias = model_parameters[0];
  return bias*bias*pp->cosmology->wp_DM(r, pp->method, pp->redshift, pp->output_root, pp->NL, pp->norm, pp->r_min, pp->r_max, pp-> k_min, pp->k_max, pp->aa, pp->GSL, pp->prec, pp->file_par);
}


// ============================================================================================


vector<double> cosmobl::glob::wp_bias_cosmology_vector(vector<double> r, shared_ptr<void> parameters, vector<double> model_parameters)
{

  shared_ptr<STR_twop_model> pp = static_pointer_cast<STR_twop_model>(parameters);

  double bias = model_parameters[0];
  vector<double> wp(r.size(),0);
  for(size_t i=0;i<wp.size();i++)
    wp[i] =bias*bias*pp->cosmology->wp_DM(r[i], pp->method, pp->redshift, pp->output_root, pp->NL, pp->norm, pp->r_min, pp->r_max, pp-> k_min, pp->k_max, pp->aa, pp->GSL, pp->prec, pp->file_par);

  return wp;
}


// ============================================================================================


double cosmobl::glob::wp_bias(double r, shared_ptr<void> parameters, vector<double> model_parameters){
  
  shared_ptr<STR_twop_model> pp = static_pointer_cast<STR_twop_model>(parameters);
  double bias = model_parameters[0];
  return bias*bias*pp->func_xi->operator()(r);
}


// ============================================================================================


vector<double> cosmobl::glob::wp_bias_vector(vector<double> r, shared_ptr<void> parameters, vector<double> model_parameters)
{
  
  shared_ptr<STR_twop_model> pp = static_pointer_cast<STR_twop_model>(parameters);
  double bias2 = pow(model_parameters[0],2);

  vector<double> wp = pp->func_xi->eval_func(r);

  for(size_t i=0;i<wp.size();i++)
    wp[i] *= bias2;

  return wp;
}


// ============================================================================================


double cosmobl::glob::xi_bias_cosmology(double r, shared_ptr<void> parameters, vector<double> model_parameters)
{

  shared_ptr<STR_twop_model> pp = static_pointer_cast<STR_twop_model>(parameters);
  double bias = model_parameters[0];
  for (size_t i =1;i<model_parameters.size();i++)
    pp->cosmology->set_parameter(pp->Cpar[i-1],model_parameters[i]);

  return bias*bias*pp->cosmology->xi_DM(r,pp->method,pp->redshift, pp->output_root, pp->NL, pp->norm, pp->k_min, pp->k_max, pp->aa, pp->GSL, pp->prec, pp->file_par);
}


// ============================================================================================


vector<double> cosmobl::glob::xi_bias_cosmology_vector(vector<double> r, shared_ptr<void> parameters, vector<double> model_parameters)
{
  shared_ptr<STR_twop_model> pp = static_pointer_cast<STR_twop_model>(parameters);
  double bias = model_parameters[0];
  for (size_t i =1;i<model_parameters.size();i++)
    pp->cosmology->set_parameter(pp->Cpar[i-1],model_parameters[i]);

  vector<double> xi(r.size(),0);
  for (size_t i=0; i<xi.size(); i++)
    xi[i] = bias*bias*pp->cosmology->xi_DM(r[i],pp->method,pp->redshift, pp->output_root, pp->NL, pp->norm, pp->k_min, pp->k_max, pp->aa, pp->GSL, pp->prec, pp->file_par);

  return xi;

}


// ============================================================================================


double cosmobl::glob::xi_bias(double r, shared_ptr<void> parameters, vector<double> model_parameters){

  shared_ptr<STR_twop_model> pp = static_pointer_cast<STR_twop_model>(parameters);
  double bias = model_parameters[0];
  return bias*bias*pp->func_xi->operator()(r);
}


// ============================================================================================


vector<double> cosmobl::glob::xi_bias_vector(vector<double> r, shared_ptr<void> parameters, vector<double> model_parameters){

  shared_ptr<STR_twop_model> pp = static_pointer_cast<STR_twop_model>(parameters);
  double bias2 = (model_parameters[0],2);

  vector<double> xi = pp->func_xi->eval_func(r);

  for (size_t i=0; i<xi.size(); i++)
    xi[i] *= bias2;

  return xi;
}


// ============================================================================================


double cosmobl::glob::xi0_bias_cosmology(double r, shared_ptr<void> parameters, vector<double> model_parameters)
{
  shared_ptr<STR_twop_model> pp = static_pointer_cast<STR_twop_model>(parameters);
  double bias = model_parameters[0];
  for (size_t i =1;i<model_parameters.size();i++)
    pp->cosmology->set_parameter(pp->Cpar[i-1],model_parameters[i]);

  double beta = pp->cosmology->beta(pp->redshift,bias);
  double fact = bias*bias*(1+2./3*beta+beta*beta/5);
  return fact*pp->cosmology->xi_DM(r,pp->method,pp->redshift, pp->output_root, pp->NL, pp->norm, pp->k_min, pp->k_max, pp->aa, pp->GSL, pp->prec, pp->file_par);

}


// ============================================================================================


vector<double> cosmobl::glob::xi0_bias_cosmology_vector(vector<double> r, shared_ptr<void> parameters, vector<double> model_parameters)
{
  shared_ptr<STR_twop_model> pp = static_pointer_cast<STR_twop_model>(parameters);
  double bias = model_parameters[0];
  for (size_t i =1;i<model_parameters.size();i++)
    pp->cosmology->set_parameter(pp->Cpar[i-1],model_parameters[i]);

  double beta = pp->cosmology->beta(pp->redshift,bias);
  double fact = bias*bias*(1+2./3*beta+beta*beta/5);

  vector<double> xi(r.size(),0);

  for(size_t i=0;i<xi.size();i++)
    xi[i] = fact*pp->cosmology->xi_DM(r[i],pp->method,pp->redshift, pp->output_root, pp->NL, pp->norm, pp->k_min, pp->k_max, pp->aa, pp->GSL, pp->prec, pp->file_par);

  return xi;

}


// ============================================================================================


double cosmobl::glob::xi0_bias(double r, shared_ptr<void> parameters, vector<double> model_parameters)
{
  shared_ptr<STR_twop_model> pp = static_pointer_cast<STR_twop_model>(parameters);

  double bias = model_parameters[0];
  double beta = pp->f/bias;
  double fact = bias*bias*(1+2./3*beta+beta*beta/5);
  return fact*pp->func_xi->operator()(r);
}


// ============================================================================================


vector<double> cosmobl::glob::xi0_bias_vector(vector<double> r, shared_ptr<void> parameters, vector<double> model_parameters)
{
  shared_ptr<STR_twop_model> pp = static_pointer_cast<STR_twop_model>(parameters);

  double bias = model_parameters[0];
  double beta = pp->f/bias;
  double fact = bias*bias*(1+2./3*beta+beta*beta/5);

  vector<double> xi = pp->func_xi->eval_func(r);

  for(size_t i=0;i<xi.size();i++)
    xi[i] *= fact;

  return xi;
}


// ============================================================================================


double cosmobl::glob::xi_alpha_B(double r, shared_ptr<void> parameters, vector<double> model_parameters)
{
  shared_ptr<STR_twop_model> pp = static_pointer_cast<STR_twop_model>(parameters);

  double new_r = model_parameters[1]*r;
//  cout << "I'm in model " << new_r << endl;
  
  return model_parameters[0]*pp->func_xi->operator()(new_r); 
}


// ============================================================================================


vector<double> cosmobl::glob::xi_alpha_B_vector(vector<double> r, shared_ptr<void> parameters, vector<double> model_parameters)
{
  shared_ptr<STR_twop_model> pp = static_pointer_cast<STR_twop_model>(parameters);
  
  vector<double> xi(r.size(),0);

  for(size_t i=0;i<xi.size();i++)
    xi[i] = model_parameters[0]*pp->func_xi->operator()(model_parameters[1]*r[i]); 

  return xi;
}


// ============================================================================================


double cosmobl::glob::xi_alpha_B_poly(double r, shared_ptr<void> parameters, vector<double> model_parameters)
{
  shared_ptr<STR_twop_model> pp = static_pointer_cast<STR_twop_model>(parameters);

  double new_r = model_parameters[1]*r;
  double poly = model_parameters[2]/(r*r)+model_parameters[3]/r+model_parameters[4];
  return model_parameters[0]*pp->func_xi->operator()(new_r)+poly; 
}


// ============================================================================================


vector<double> cosmobl::glob::xi_alpha_B_poly_vector(vector<double> r, shared_ptr<void> parameters, vector<double> model_parameters)
{
  shared_ptr<STR_twop_model> pp = static_pointer_cast<STR_twop_model>(parameters);

  vector<double> xi(r.size(),0);

  for(size_t i=0;i<xi.size();i++){
    double poly = model_parameters[2]/(r[i]*r[i])+model_parameters[3]/r[i]+model_parameters[4];
    xi[i] =model_parameters[0]*pp->func_xi->operator()(model_parameters[1]*r[i])+poly; 
  }

  return xi;
}


// ============================================================================================


double cosmobl::glob::xi0_alpha_bias(double r, shared_ptr<void> parameters, vector<double> model_parameters)
{
  shared_ptr<STR_twop_model> pp = static_pointer_cast<STR_twop_model>(parameters);

  double bias = model_parameters[0];
  double beta = pp->f/bias;
  double fact = bias*bias*(1+2./3*beta+beta*beta/5);
  double new_r = model_parameters[1]*r;

  return fact*pp->func_xi->operator()(new_r); 
}


// ============================================================================================


vector<double> cosmobl::glob::xi0_alpha_bias_vector(vector<double> r, shared_ptr<void> parameters, vector<double> model_parameters)
{
  shared_ptr<STR_twop_model> pp = static_pointer_cast<STR_twop_model>(parameters);

  double bias = model_parameters[0];
  double beta = pp->f/bias;
  double fact = bias*bias*(1+2./3*beta+beta*beta/5);

  vector<double> xi(r.size(),0);

  for(size_t i=0;i<xi.size();i++)
    xi[i] = fact*pp->func_xi->operator()(model_parameters[1]*r[i]); 

  return xi;
}


// ============================================================================================


double cosmobl::glob::xi0_alpha_bias_cosmology(double r, shared_ptr<void> parameters, vector<double> model_parameters)
{
  shared_ptr<STR_twop_model> pp = static_pointer_cast<STR_twop_model>(parameters);
  double bias = model_parameters[0];
  double new_r = model_parameters[1]*r;

  for (size_t i =2;i<model_parameters.size();i++)
    pp->cosmology->set_parameter(pp->Cpar[i-2],model_parameters[i]);

  double beta = pp->cosmology->beta(pp->redshift,bias);
  double fact = bias*bias*(1+2./3*beta+beta*beta/5);
  return fact*pp->cosmology->xi_DM(new_r,pp->method,pp->redshift, pp->output_root, pp->NL, pp->norm, pp->k_min, pp->k_max, pp->aa, pp->GSL, pp->prec, pp->file_par);
}


// ============================================================================================


vector<double> cosmobl::glob::xi0_alpha_bias_cosmology_vector(vector<double> r, shared_ptr<void> parameters, vector<double> model_parameters)
{
  shared_ptr<STR_twop_model> pp = static_pointer_cast<STR_twop_model>(parameters);
  double bias = model_parameters[0];

  for (size_t i =2;i<model_parameters.size();i++)
    pp->cosmology->set_parameter(pp->Cpar[i-2],model_parameters[i]);

  double beta = pp->cosmology->beta(pp->redshift,bias);
  double fact = bias*bias*(1+2./3*beta+beta*beta/5);

  vector<double> xi(r.size(),0);
  for(size_t i=0;i<xi.size();i++)
    xi[i] = fact*pp->cosmology->xi_DM(model_parameters[1]*r[i],pp->method,pp->redshift, pp->output_root, pp->NL, pp->norm, pp->k_min, pp->k_max, pp->aa, pp->GSL, pp->prec, pp->file_par);

  return xi;
}
