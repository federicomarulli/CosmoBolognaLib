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

using namespace std;

using namespace cbl;


// ===========================================================================================


void cbl::modelling::numbercounts::Modelling_NumberCounts::set_data_model (const cosmology::Cosmology cosmology, const double redshift, const std::string method_Pk, const double k_min, const double k_max, const int step, const bool store_output, const int norm, const double Delta, const bool isDelta_vir, const std::string model_MF, const std::string selection_function_file, const std::vector<int> selection_function_column, const double z_min, const double z_max, const int z_step, const double Mass_min, const double Mass_max, const int Mass_step, const double area_degrees, const double prec)
{
  m_data_model.isSnapshot = false;

  if (m_fit_range==false) ErrorCBL("You must set the fit range (through set_fit_range) first!","set_data_model","Modelling_NumberCounts.cpp");
  m_data_model.edges_x = m_data->edges_xx();
  m_data_model.edges_y = m_data->edges_yy();

  m_data_model.cosmology = make_shared<cosmology::Cosmology>(cosmology);
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

  m_data_model.isDelta_Vir = isDelta_vir;
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
    m_data_model.Volume = cosmology.Volume(z_min, z_max, area_degrees);

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


void cbl::modelling::numbercounts::Modelling_NumberCounts::set_data_model_SF (const cosmology::Cosmology cosmology, const std::vector<double> radii, const double redshift, const std::string model_SF, const double b_eff, double slope, double offset, const double deltav_NL, const double del_c, const std::string method_Pk, const double k_Pk_ratio, const bool store_output, const std::string output_root, const std::string interpType, const double k_max, const std::string input_file, const bool is_parameter_file)
{
  m_data_model_SF.cosmology = make_shared<cosmology::Cosmology>(cosmology);
  m_data_model_SF.radii = radii;
  m_data_model_SF.redshift = redshift;
  m_data_model_SF.model_SF = model_SF;
  m_data_model_SF.b_eff = b_eff;
  m_data_model_SF.b_slope = slope;
  m_data_model_SF.b_offset = offset;
  m_data_model_SF.deltav_NL = deltav_NL;
  m_data_model_SF.delta_c = del_c;
  m_data_model_SF.method_Pk = method_Pk;
  m_data_model_SF.k_Pk_ratio = k_Pk_ratio;
  m_data_model_SF.store_output = store_output;
  m_data_model_SF.output_root = output_root;
  m_data_model_SF.interpType = interpType;
  m_data_model_SF.k_max = k_max;
  m_data_model_SF.input_file = input_file;
  m_data_model_SF.is_parameter_file = is_parameter_file;  

}
