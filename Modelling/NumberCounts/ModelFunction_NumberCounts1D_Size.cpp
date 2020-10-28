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
 *  Modelling/NumberCounts/ModelFunction_NumberCounts1D_Size.cpp
 *
 *  @brief Functions to model the size number counts
 *
 *  This file contains the implementation of the functions used to
 *  model mass number counts
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unibo.it, alfonso.veropalumbo@unibo.it
 */


#include "ModelFunction_NumberCounts1D_Size.h"

using namespace std;

using namespace cbl;

// ===========================================================================================

std::vector<double> cbl::modelling::numbercounts::size_function_model (const std::vector<double> radii, const std::shared_ptr<void> inputs, std::vector<double> &parameter)
{
  // structure contaning the required input data
  shared_ptr<STR_NCSF_data_model> pp = static_pointer_cast<STR_NCSF_data_model>(inputs);

  // redefine the cosmology
  cbl::cosmology::Cosmology cosmo = *pp->cosmology;
  
  for (size_t i=0; i<pp->Cpar.size(); ++i)
    cosmo.set_parameter(pp->Cpar[i], parameter[i]);

  double bias_eff;
  double bias_slope;
  double bias_offset;
  if (parameter.size()==pp->Cpar.size()+3){
    bias_eff = parameter[parameter.size()-3];
    bias_slope = parameter[parameter.size()-2];
    bias_offset = parameter[parameter.size()-1];}
  else {
    bias_eff = pp->b_eff;
    bias_slope = pp->b_slope;
    bias_offset = pp->b_offset;}

  return cbl::modelling::numbercounts::size_function(cosmo, radii, pp->redshift, pp->model_SF, bias_eff, bias_slope, bias_offset, pp->deltav_NL, pp->delta_c, pp->method_Pk, pp->k_Pk_ratio, pp->store_output, pp->output_root, pp->interpType, pp->k_max, pp->input_file, pp->is_parameter_file);
  
}



