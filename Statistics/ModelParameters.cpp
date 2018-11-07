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
 *  @file Statistics/ModelParameters.cpp
 *
 *  @brief Methods of the class ModelParameters 
 *
 *  This file contains the implementation of the methods of the class
 *  ModelParameters
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo 
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#include "ModelParameters.h"

using namespace std;

using namespace cbl;


// ============================================================================================


vector<double> cbl::statistics::ModelParameters::full_parameters (const vector<double> parameter_values) const
{
  if(parameter_values.size() == m_nparameters_base){
    vector<double> all_parameters(m_nparameters, 0);

    for(size_t i=0; i<m_nparameters_base; i++)
      all_parameters[m_base_parameters[i]] = parameter_values[i];

    for(size_t i=0; i<m_nparameters_derived; i++)
      all_parameters[m_derived_parameters[i]] = 0.;

    return all_parameters;
  }
  else if (parameter_values.size() == m_nparameters)
    return parameter_values;
  else
    ErrorCBL("Error in cbl::statistics::ModelParameters::full_parameters() of ModelParameters.cpp, provided vector has the wrong size!");

  vector<double> vv;
  return vv;
}


// ============================================================================================


void cbl::statistics::ModelParameters::m_set_parameter_type ()
{
  m_nparameters_base = 0;
  m_nparameters_derived = 0;

  m_base_parameters.erase(m_base_parameters.begin(), m_base_parameters.end());
  m_derived_parameters.erase(m_derived_parameters.begin(), m_derived_parameters.end());

  for (size_t i=0; i<m_nparameters; i++) {
    switch (m_parameter_type[i]) {
      case statistics::ParameterType::_Base_:
	m_nparameters_base +=1;
	m_base_parameters.push_back(i);
	break;

      case statistics::ParameterType::_Derived_:
	m_nparameters_derived += 1;
	m_derived_parameters.push_back(i);
	break;

      default:
	ErrorCBL("Error in cbl::statistics::ModelParameters::m_set_parameter_type() of ModelParameters.cpp: no such kind of parameter!");
    }
  }
}


// ============================================================================================


cbl::statistics::ModelParameters::ModelParameters (const size_t nparameters, std::vector<ParameterType> parameterTypes, std::vector<std::string> parameterNames)
{
  set_parameters(nparameters, parameterTypes, parameterNames);
}


// ============================================================================================

  
size_t cbl::statistics::ModelParameters::nparameters () const
{
  return m_nparameters;
}


// ============================================================================================


size_t cbl::statistics::ModelParameters::nparameters_base () const
{
  return m_nparameters_base;
}


// ============================================================================================


size_t cbl::statistics::ModelParameters::nparameters_derived () const
{
  return m_nparameters_derived;
}


// ============================================================================================


void cbl::statistics::ModelParameters::set_parameters (const size_t nparameters, std::vector<ParameterType> parameterTypes, std::vector<std::string> parameterNames)
{
  //Check parameterTypes size
  
  if (nparameters==0)
    ErrorCBL("Error in cbl::statistics::ModelParameters::set_parameters() of ModelParameters.cpp: nparameters should be > 0.");

  if ((parameterTypes.size()!=nparameters) && (parameterTypes.size()!=0))
    ErrorCBL("Error in cbl::statistics::ModelParameters::set_parameters() of ModelParameters.cpp: wrong size for the vector parameterTypes.");

  if ((parameterNames.size()!=nparameters) && (parameterNames.size()!=0))
    ErrorCBL("Error in cbl::statistics::ModelParameters::set_parameters() of ModelParameters.cpp: wrong size for the vector parameterNames.");


  if ((parameterTypes.size()==nparameters) && (parameterNames.size()==nparameters)) {
    m_nparameters=nparameters;
    m_parameter_type = parameterTypes;
    m_parameter_name = parameterNames;
  }
  else if ((parameterTypes.size()==0) && (parameterNames.size()==0)){
    m_nparameters=nparameters;
    vector<ParameterType> pTypes(m_nparameters);
    vector<string> pNames(m_nparameters);
    for(size_t i=0; i<m_nparameters; i++){
      pTypes[i] = ParameterType::_Base_;
      pNames[i] = "par_"+conv(i+1, par::fINT);
    }
    m_parameter_type = pTypes;
    m_parameter_name = pNames;
  }
  else if ((parameterTypes.size()==nparameters) && (parameterNames.size()==0)){
    m_nparameters=nparameters;
    vector<string> pNames(m_nparameters);
    for(size_t i=0; i<m_nparameters; i++)
      pNames[i] = "par_"+conv(i+1, par::fINT);
    
    m_parameter_type = parameterTypes;
    m_parameter_name = pNames;
  }
  else if ((parameterTypes.size()==0) && (parameterNames.size()==0)){
    m_nparameters=nparameters;
    vector<ParameterType> pTypes(m_nparameters);
    for(size_t i=0; i<m_nparameters; i++)
      pTypes[i] = ParameterType::_Base_;
    
    m_parameter_type = pTypes;
    m_parameter_name = parameterNames;
  }

  m_set_parameter_type();
}
