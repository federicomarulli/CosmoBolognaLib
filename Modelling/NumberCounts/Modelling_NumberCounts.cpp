/********************************************************************
 *  Copyright (C) 2016 by Federico Marulli and Alfonso Veropalumbo  *
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
 *  @file
 *  Modelling/NumberCounts/Modelling_NumberCounts.cpp
 *
 *  @brief Methods of the class Modelling_NumberCounts
 *
 *  This file contains the implementation of the methods of the class
 *  Modelling_NumberCounts, i.e. the common functions to model
 *  the number counts of any kind
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unibo.it, alfonso.veropalumbo@unibo.it
 */


#include "Modelling_NumberCounts.h"
#include "Data1D.h"

using namespace std;

using namespace cbl;


// ===========================================================================================


void cbl::modelling::numbercounts::Modelling_NumberCounts::set_P_proxy (const double A_mu, const double B_mu, const double C_mu, const double A_sigma, const double B_sigma, const double C_sigma)
{
  if (m_model != NULL)
    ErrorCBL("This function must be called before the definition of the model!", "set_P_proxy", "Modelling_DensityProfile.cpp");
  
  m_data_model.Plambda_A_mu = A_mu;
  m_data_model.Plambda_B_mu = B_mu;
  m_data_model.Plambda_C_mu = C_mu;
  m_data_model.Plambda_A_sigma = A_sigma;
  m_data_model.Plambda_B_sigma = B_sigma;
  m_data_model.Plambda_C_sigma = C_sigma;
  
  m_data_model.Plambda_mean_fc = [] (const double proxy_tr, const double z_tr, const double A_mu, const double B_mu, const double C_mu)
				 {
				   return proxy_tr + A_mu * proxy_tr * exp(- proxy_tr * (B_mu + C_mu * z_tr));
				 };

  m_data_model.Plambda_std_fc = [] (const double proxy_tr, const double z_tr, const double A_sigma, const double B_sigma, const double C_sigma)
				 {
				   return A_sigma * proxy_tr * exp(- proxy_tr * (B_sigma + C_sigma * z_tr));
				 };

  m_data_model.isSet_P_proxy = true;
}


// ===========================================================================================


void cbl::modelling::numbercounts::Modelling_NumberCounts::set_data_model (const std::shared_ptr<cosmology::Cosmology> cosmology, const std::vector<double> redshift_points_completeness, const std::vector<double> proxy_points_completeness, const std::vector<std::vector<double>> completeness, const std::vector<double> purity, const double z_error_normalised, const std::string scalrel_z_evo, const double z_pivot, const double proxy_pivot, const double mass_pivot, const double log_base, const std::string method_Pk, const bool store_output, const int norm, const double Delta, const bool isDelta_critical, const std::string model_MF, const std::string model_bias, const double z_min, const double z_max, const double area_degrees, const double prec)
{
  if (m_data_model.isSet_P_proxy == false) ErrorCBL("You must set the P(proxy_ob|proxy_tr,z_tr) PDF first! Use set_P_proxy().","set_data_model","Modelling_NumberCounts.cpp");
  
  m_data_model.isSnapshot = false;

  if (m_fit_range==false) ErrorCBL("You must set the fit range (through set_fit_range) first!","set_data_model","Modelling_NumberCounts.cpp");
  m_data_model.edges_x = m_data_fit->edges_xx();

  for (size_t i=0; i<m_data_model.edges_x.size(); i++)
    if (m_data_model.edges_x[i] <= 0)
      ErrorCBL("The values of the proxy edges cannot be <= 0.","set_data_model","Modelling_NumberCounts.cpp");

  m_data_model.cosmology = move(cosmology);

  if ( (redshift_points_completeness.size() > 0) && (proxy_points_completeness.size() > 0) )
    m_data_model.completeness_interp = cbl::glob::FuncGrid2D(redshift_points_completeness, proxy_points_completeness, completeness, "Cubic");  
  else {    
    std::vector<double> dummy_values = cbl::linear_bin_vector(10, 0., 10000.);
    std::vector<std::vector<double>> dummy_completeness(dummy_values.size(), std::vector<double>(dummy_values.size(), 1.));
    m_data_model.completeness_interp = cbl::glob::FuncGrid2D(dummy_values, dummy_values, dummy_completeness, "Cubic");
  }

  if (purity.size() != m_data_model.edges_x.size() - 1)
    ErrorCBL("The purity vector must have the same size of the observed mass proxy bin vector.","set_data_model","Modelling_NumberCounts.cpp");
  m_data_model.purity = purity;

  m_data_model.z_error_normalised = z_error_normalised;
  
  m_data_model.z_pivot = z_pivot;
  m_data_model.proxy_pivot = proxy_pivot;
  m_data_model.mass_pivot = mass_pivot;
  m_data_model.log_base = log_base;
  
  m_data_model.method_Pk = method_Pk;
  m_data_model.kk = logarithmic_bin_vector(500, 1.e-4, 100.);
  m_data_model.norm = norm;
  m_data_model.store_output = store_output;
  m_data_model.output_root = "test";
  m_data_model.file_par = par::defaultString;

  m_data_model.isDelta_critical = isDelta_critical;
  m_data_model.Delta = Delta;
  m_data_model.model_MF = model_MF;
  m_data_model.model_bias = model_bias;

  m_data_model.Mass_vector = logarithmic_bin_vector(200, 1.e10, 1.e16);

  m_data_model.z_min = z_min;
  m_data_model.z_max = z_max;

  m_data_model.prec = prec;

  m_data_model.area_rad = area_degrees*pow(par::pi/180., 2);

  m_data_model.is_sigma8_free = false;

  // build a dummy dataset for the scaling relation Modelling object, useful only to avoid internal errors
  std::vector<double> dummy_vec = {1.};
  std::shared_ptr<cbl::data::Data> dataset = std::make_shared<cbl::data::Data1D>(cbl::data::Data1D(dummy_vec, dummy_vec, dummy_vec));

  // set the scaling relation model
  modelling::massobsrel::Modelling_MassObservableRelation scaling_relation (dataset);
  m_data_model.scaling_relation = make_shared<modelling::massobsrel::Modelling_MassObservableRelation>(scaling_relation);
  (m_data_model.scaling_relation)->set_data_model(cosmology, {0.}, z_pivot, proxy_pivot, log_base);

  m_data_model.z_evo = scalrel_z_evo;
}


// ===========================================================================================


void cbl::modelling::numbercounts::Modelling_NumberCounts::set_data_model (const std::shared_ptr<cosmology::Cosmology> cosmology, const double redshift, const std::string method_Pk, const double k_min, const double k_max, const int step, const bool store_output, const int norm, const double Delta, const bool isDelta_critical, const std::string model_MF, const std::string selection_function_file, const std::vector<int> selection_function_column, const double z_min, const double z_max, const int z_step, const double Mass_min, const double Mass_max, const int Mass_step, const double area_degrees, const double prec)
{
  m_data_model.isSnapshot = false;

  if (m_fit_range==false) ErrorCBL("You must set the fit range (through set_fit_range) first!","set_data_model","Modelling_NumberCounts.cpp");
  m_data_model.edges_x = m_data_fit->edges_xx();
  m_data_model.edges_y = m_data_fit->edges_yy();

  m_data_model.cosmology = move(cosmology);
  m_data_model.redshift = redshift;
  m_data_model.method_Pk = method_Pk;
  m_data_model.k_min = k_min;
  m_data_model.k_max = k_max;
  m_data_model.step = step;
  m_data_model.kk = logarithmic_bin_vector(step, k_min, k_max);
  m_data_model.norm = norm;
  
  m_data_model.store_output = store_output;
  m_data_model.output_root = "test";
  m_data_model.file_par = par::defaultString;

  m_data_model.isDelta_critical = isDelta_critical;
  m_data_model.Delta = Delta;
  m_data_model.model_MF = model_MF;

  m_data_model.Mass_min = Mass_min;
  m_data_model.Mass_max = Mass_max;
  m_data_model.Mass_step = Mass_step;
  m_data_model.Mass_vector = logarithmic_bin_vector(200, 1.e10, 1.e16);

  m_data_model.z_min = z_min;
  m_data_model.z_max = z_max;
  m_data_model.z_step = z_step;
  m_data_model.z_vector = linear_bin_vector(z_step, z_min, z_max);

  m_data_model.prec = prec;

  m_data_model.area_rad = area_degrees*pow(par::pi/180., 2);
  if (m_data_model.z_min>0)
    m_data_model.Volume = cosmology->Volume(z_min, z_max, area_degrees);

  if (selection_function_file!=par::defaultString) {
    m_data_model.use_SF = true;
    std::vector<double> mass, redshift;
    std::vector<std::vector<double>> SF;
    read_matrix(selection_function_file, redshift, mass, SF, selection_function_column);
    m_data_model.interp_SelectionFunction = make_shared<glob::FuncGrid2D> (glob::FuncGrid2D(redshift, mass, SF, "Cubic"));
  }
  else
    m_data_model.use_SF = false;

  m_data_model.is_sigma8_free = false;
}



// ===========================================================================================


void cbl::modelling::numbercounts::Modelling_NumberCounts::set_data_model (const std::shared_ptr<cosmology::Cosmology> cosmology, const std::vector<double> SF_weights, const double z_pivot, const double proxy_pivot, const double mass_pivot, const double log_base, const std::string method_Pk, const bool store_output, const int norm, const double Delta, const bool isDelta_critical, const std::string model_MF, const std::string model_bias, const double z_min, const double z_max, const double area_degrees, const double prec)
{
  m_data_model.isSnapshot = false;

  if (m_fit_range==false) ErrorCBL("You must set the fit range (through set_fit_range) first!","set_data_model","Modelling_NumberCounts.cpp");
  m_data_model.edges_x = m_data_fit->edges_xx();

  for (size_t i=0; i<m_data_model.edges_x.size(); i++)
    if (m_data_model.edges_x[i] <= 0)
      ErrorCBL("The values of the proxy edges cannot be <= 0.","set_data_model","Modelling_NumberCounts.cpp");

  m_data_model.cosmology = move(cosmology);

  if (SF_weights.size() == (m_data_fit->edges_xx()).size()-1)
    m_data_model.SF_weights = SF_weights;
  else
    ErrorCBL("The weights vector must have the same size of the x vector!","set_data_model","Modelling_NumberCounts.cpp");
  m_data_model.z_pivot = z_pivot;
  m_data_model.proxy_pivot = proxy_pivot;
  m_data_model.mass_pivot = mass_pivot;
  m_data_model.log_base = log_base;
  
  m_data_model.method_Pk = method_Pk;
  m_data_model.k_min = 1.e-4;
  m_data_model.k_max = 100.;
  m_data_model.step = 500;
  m_data_model.kk = logarithmic_bin_vector(500, 1.e-4, 100.);
  m_data_model.norm = norm;
  
  m_data_model.store_output = store_output;
  m_data_model.output_root = "test";
  m_data_model.file_par = par::defaultString;

  m_data_model.isDelta_critical = isDelta_critical;
  m_data_model.Delta = Delta;
  m_data_model.model_MF = model_MF;
  m_data_model.model_bias = model_bias;

  m_data_model.Mass_vector = logarithmic_bin_vector(200, 1.e10, 1.e16);

  m_data_model.z_min = z_min;
  m_data_model.z_max = z_max;

  m_data_model.prec = prec;

  m_data_model.area_rad = area_degrees*pow(par::pi/180., 2);
  if (m_data_model.z_min>0)
    m_data_model.Volume = cosmology->Volume(z_min, z_max, area_degrees);

  m_data_model.is_sigma8_free = false;
}



// ===========================================================================================


void cbl::modelling::numbercounts::Modelling_NumberCounts::set_data_model_SF (const std::shared_ptr<cosmology::Cosmology> cosmology, const std::vector<double> radii, const double redshift, const std::string model_SF, const double b_eff, double slope, double offset, const double deltav_NL, const double del_c, const std::string method_Pk, const bool store_output, const std::string output_root, const std::string interpType, const double k_max, const std::string input_file, const bool is_parameter_file)
{
  m_data_model_SF.cosmology = move(cosmology);
  m_data_model_SF.radii = radii;
  m_data_model_SF.redshift = redshift;
  m_data_model_SF.model_SF = model_SF;
  m_data_model_SF.b_eff = b_eff;
  m_data_model_SF.b_slope = slope;
  m_data_model_SF.b_offset = offset;
  m_data_model_SF.deltav_NL = deltav_NL;
  m_data_model_SF.delta_c = del_c;
  m_data_model_SF.method_Pk = method_Pk;
  m_data_model_SF.store_output = store_output;
  m_data_model_SF.output_root = output_root;
  m_data_model_SF.interpType = interpType;
  m_data_model_SF.k_max = k_max;
  m_data_model_SF.input_file = input_file;
  m_data_model_SF.is_parameter_file = is_parameter_file;  

}
