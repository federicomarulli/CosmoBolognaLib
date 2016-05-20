/********************************************************************
 *  Copyright (C) 2015 by Federico Marulli and Alfonso Veropalumbo  *
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
 *  @file Modelling/ModelBias.cpp
 *
 *  @brief Methods of the class ModelBias, used for modelling bias in two
 *  point correlation function
 *
 *  This file contains the implementation of the methods of the class
 *  ModelBias
 *
 *  @authors Federico Marulli
 *
 *  @authors federico.marulli3@unbo.it
 */


#include "ModelBias.h"

using namespace cosmobl;


// ============================================================================================


void cosmobl::ModelBias::set_xi_parameters (const vector<double> r, const shared_ptr<Cosmology> cosmology, const double redshift, const string type, const int nPt, const string method, const string output_root, const bool NL, const int norm, const double k_min, const double k_max, const double aa, const bool GSL, const double prec, const string file_par)
{
  auto model_parameters = make_shared<cosmobl::glob::STR_twop_model>(cosmobl::glob::STR_twop_model());

  vector<double> xi(r.size(),0);

  for (size_t i=0; i<r.size(); i++)
    xi[i] = cosmology->xi_DM(r[i], method, redshift, output_root, 1, norm, k_min, k_max, aa, GSL, prec, file_par);

  model_parameters->func_xi = make_shared<classfunc::func_grid_GSL>(classfunc::func_grid_GSL(r,xi,type));
  m_model_parameters = move(model_parameters);

  m_model = &cosmobl::glob::xi_bias;
}


// ============================================================================================


void cosmobl::ModelBias::set_xi_parameters_cosmology (const shared_ptr<Cosmology> cosmology, const vector<CosmoPar> cosmo_parameters, const double redshift, const string method, const string output_root, const bool NL, const int norm, const double k_min, const double k_max, const double aa, const bool GSL, const double prec, const string file_par)
{

  auto model_parameters = make_shared<cosmobl::glob::STR_twop_model>(cosmobl::glob::STR_twop_model());

  model_parameters->cosmology=cosmology;
  model_parameters->Cpar=cosmo_parameters;
  model_parameters->redshift=redshift;
  model_parameters->method=method;
  model_parameters->output_root=output_root;
  model_parameters->NL=NL; 
  model_parameters->norm=norm; 
  model_parameters->k_min=k_min;  
  model_parameters->k_max=k_max;
  model_parameters->aa=aa;  
  model_parameters->GSL=GSL;
  model_parameters->prec=prec; 
  model_parameters->file_par=file_par;

  m_model_parameters = move(model_parameters);

  m_model = &cosmobl::glob::xi_bias_cosmology;
}


// ============================================================================================


void cosmobl::ModelBias::set_wp_parameters (const vector<double> r, const shared_ptr<Cosmology> cosmology, const double redshift, const double pi_max, const string type, const int nPt, const string method, const string output_root, const bool NL, const int norm, const double r_min, const double r_max, const double k_min, const double k_max, const double aa, const bool GSL, const double prec, const string file_par)
{
  auto model_parameters = make_shared<cosmobl::glob::STR_twop_model>(cosmobl::glob::STR_twop_model());

  vector<double> wp(r.size(),0);

  for (size_t i=0; i<r.size(); i++)
    wp[i] = cosmology->wp_DM(r[i], method, redshift, output_root, NL, norm, r_min, r_max, k_min, k_max, aa, GSL, prec, file_par);
  

  model_parameters->func_xi = make_shared<classfunc::func_grid_GSL>(classfunc::func_grid_GSL(r,wp,type));
  m_model_parameters = move(model_parameters);

  m_model = &cosmobl::glob::wp_bias;
}


// ============================================================================================


void cosmobl::ModelBias::set_wp_parameters_cosmology (const shared_ptr<Cosmology> cosmology, const vector<CosmoPar> cosmo_parameters, const double redshift, const double pi_max, const string method, const string output_root, const int norm, const double r_min, const double r_max, const double k_min, const double k_max, const double aa, const bool GSL, const double prec, const string file_par)
{

  auto model_parameters = make_shared<cosmobl::glob::STR_twop_model>(cosmobl::glob::STR_twop_model());

  model_parameters->cosmology=cosmology;
  model_parameters->Cpar=cosmo_parameters;
  model_parameters->redshift=redshift;
  model_parameters->pi_max=pi_max;
  model_parameters->method=method;
  model_parameters->output_root=output_root;
  model_parameters->norm=norm; 
  model_parameters->r_min=r_min;
  model_parameters->r_max=r_max;
  model_parameters->k_min=k_min; 
  model_parameters->k_max=k_max;
  model_parameters->aa=aa;
  model_parameters->GSL=GSL;
  model_parameters->prec=prec; 
  model_parameters->file_par=file_par;

  m_model_parameters = move(model_parameters);

  m_model = &cosmobl::glob::wp_bias_cosmology;
}


// ============================================================================================


void cosmobl::ModelBias::set_xi0_parameters (const vector<double> r, const shared_ptr<Cosmology> cosmology, const double redshift, const string type, const int nPt, const string method, const string output_root, const bool NL, const int norm, const double k_min, const double k_max, const double aa, const bool GSL, const double prec, const string file_par)
{
  auto model_parameters = make_shared<cosmobl::glob::STR_twop_model>(cosmobl::glob::STR_twop_model());

  vector<double> xi(r.size(),0);

  for (size_t i=0; i<r.size(); i++)
    xi[i] = cosmology->xi_DM(r[i], method, redshift, output_root, 1, norm, k_min, k_max, aa, GSL, prec, file_par);

  model_parameters->func_xi = make_shared<classfunc::func_grid_GSL>(classfunc::func_grid_GSL(r,xi,type));
  m_model_parameters = move(model_parameters);

  m_model = &cosmobl::glob::xi0_bias;
}


// ============================================================================================


void cosmobl::ModelBias::set_xi0_parameters_cosmology (const shared_ptr<Cosmology> cosmology, const vector<CosmoPar> cosmo_parameters, const double redshift, const string method, const string output_root, const bool NL, const int norm, const double k_min, const double k_max, const double aa, const bool GSL, const double prec, const string file_par)
{

  auto model_parameters = make_shared<cosmobl::glob::STR_twop_model>(cosmobl::glob::STR_twop_model());

  model_parameters->cosmology=cosmology;
  model_parameters->Cpar=cosmo_parameters;
  model_parameters->redshift=redshift;
  model_parameters->method=method;
  model_parameters->output_root=output_root;
  model_parameters->NL=NL; 
  model_parameters->norm=norm; 
  model_parameters->k_min=k_min;  
  model_parameters->k_max=k_max;
  model_parameters->aa=aa;  
  model_parameters->GSL=GSL;
  model_parameters->prec=prec; 
  model_parameters->file_par=file_par;

  m_model_parameters = move(model_parameters);

  m_model = &cosmobl::glob::xi0_bias_cosmology;
}


// ============================================================================================


cosmobl::ModelBias::ModelBias (const double bias_value, const statistics::Prior bias_prior) : Model1D()
{
  m_parameters.resize(1); 
  m_parameters[0] = make_shared<statistics::Parameter>(bias_value,bias_prior,0,"bias");
  m_npar=m_parameters.size();
}


// ============================================================================================


cosmobl::ModelBias::ModelBias (const double bias_value, const statistics::Prior bias_prior, const vector<CosmoPar> cosmo_parameters, const vector<double> cosmo_parameters_values, const vector<statistics::Prior> cosmo_parameters_priors, const vector<string> cpar_name) : Model1D()
{

  m_npar = cosmo_parameters.size()+1;
  m_parameters.resize(m_npar); 
  m_parameters[0] = make_shared<statistics::Parameter>(bias_value,bias_prior,0,"bias");

  for (unsigned int i=0; i<m_npar-1; i++) {
    string par_name = (cpar_name.size()!=m_npar-1) ? "par"+conv(i+1, par::fINT) : cpar_name[i];
    m_parameters[i+1] = make_shared<statistics::Parameter>(cosmo_parameters_values[i], cosmo_parameters_priors[i],0, par_name);
  }
  m_npar=m_parameters.size();
}
