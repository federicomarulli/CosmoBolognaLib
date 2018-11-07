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
 *  @file Statistics/Model.cpp
 *
 *  @brief Methods of the class Model
 *
 *  This file contains the implementation of the methods of the class
 *  Model
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#include "Model.h"

using namespace std;

using namespace cbl;

// ======================================================================================


vector<vector<double>> cbl::statistics::Model::operator() (const std::vector<std::vector<double>> xx, std::vector<double> &parameters) const
{
  parameters = m_parameters->full_parameters(parameters);
  return m_function(xx, m_inputs, parameters);
}


// ======================================================================================


void cbl::statistics::Model::set_parameters (const size_t nparameters, std::vector<ParameterType> parameterTypes, std::vector<std::string> parameterNames)
{
  m_parameters = make_shared<cbl::statistics::ModelParameters>(cbl::statistics::ModelParameters(nparameters, parameterTypes, parameterNames));
}


// ======================================================================================


void cbl::statistics::Model::stats_from_chains (const std::vector<std::vector<double>> xx, std::vector<std::vector<double>> &median_model, std::vector<std::vector<double>> &low_model, std::vector<std::vector<double>> &up_model, const int start, const int thin) 
{
  int sz1, sz2;
  if (m_dimension==Dim::_1D_) {
    sz1 = 1;
    sz2 = int(xx[0].size());
  }
  else if (m_dimension==Dim::_2D_) {
    sz1 = int(xx[0].size());
    sz2 = int(xx[1].size());
  }
  else
    ErrorCBL("Error in cbl::statistics::Model::stats_from_chains() of Model.cpp: wrong size for xx vector!");
  
  vector<double> _median(sz1*sz2, 0);
  vector<double> _low(sz1*sz2, 0);
  vector<double> _up(sz1*sz2, 0);

  vector<vector<double>> models;

  for (size_t j=start; j<m_parameters->chain_size(); j+=thin) {
    for (size_t i=0; i<m_parameters->chain_nwalkers(); i++) {
      vector<double> parameters;

      for (size_t k=0; k<m_parameters->nparameters(); k++) 
	parameters.push_back(m_parameters->chain_value(k, j, i));

      models.push_back(flatten(this->operator()(xx, parameters)));
    }
  }

  vector<vector<double>> tr_models = transpose(models);
  
  for (size_t i=0; i<tr_models.size(); i++) {
    vector<double> vv = tr_models[i];
    sort(vv.begin(), vv.end());
    int low = vv.size()*0.16;
    int up = vv.size()*0.84;
    int median = vv.size()*0.5;
    _median[i] = vv[median];
    _low[i] = vv[low];
    _up[i] = vv[up];
  }

  median_model = reshape(_median, sz1, sz2);
  low_model = reshape(_low, sz1, sz2);
  up_model = reshape(_up, sz1, sz2);
}
