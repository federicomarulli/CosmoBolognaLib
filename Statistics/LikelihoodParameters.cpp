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
 *  @file Statistics/LikelihoodParameters.cpp
 *
 *  @brief Methods of the class LikelihoodParameters 
 *
 *  This file contains the implementation of the methods of the class
 *  LikelihoodParameters
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo 
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#include "LikelihoodParameters.h"

using namespace std;

using namespace cbl;

// ============================================================================================


void cbl::statistics::LikelihoodParameters::reset()
{
  set_parameters(m_nparameters, m_parameter_type, m_parameter_name);
}


// ============================================================================================


vector<double> cbl::statistics::LikelihoodParameters::full_parameters (const vector<double> parameter_values) const
{
  if(parameter_values.size() == m_nparameters_free){
    vector<double> all_parameters(m_nparameters, 0);

    for(size_t i=0; i<m_nparameters_free; i++)
      all_parameters[m_free_parameters[i]] = parameter_values[i];

    for(size_t i=0; i<m_nparameters_fixed; i++)
      all_parameters[m_fixed_parameters[i]] = m_parameter_fixed_value[m_fixed_parameters[i]];

    for(size_t i=0; i<m_nparameters_derived; i++)
      all_parameters[m_derived_parameters[i]] = 0.;

    return all_parameters;
  }
  else if (parameter_values.size() == m_nparameters)
    return parameter_values;
  else
    ErrorCBL("Error in parameters of LikelihoodParameters.cpp, the vector of values of free parameters has the wrong size!");

  vector<double> vv;
  return vv;
}


// ============================================================================================


void cbl::statistics::LikelihoodParameters::m_set_parameter_type ()
{
  m_nparameters_free = 0;
  m_nparameters_fixed = 0;
  m_nparameters_base = 0;
  m_nparameters_derived = 0;

  m_base_parameters.erase(m_base_parameters.begin(), m_base_parameters.end());
  m_fixed_parameters.erase(m_fixed_parameters.begin(), m_fixed_parameters.end());
  m_free_parameters.erase(m_free_parameters.begin(), m_free_parameters.end());
  m_derived_parameters.erase(m_derived_parameters.begin(), m_derived_parameters.end());

  for (size_t i=0; i<m_nparameters; i++) {
    switch (m_parameter_type[i]) {
      case statistics::ParameterType::_Base_:
	if (m_parameter_isFixed[i]) {
	  m_nparameters_fixed +=1;
	  m_fixed_parameters.push_back(i);
	}
	else {
	  m_nparameters_free +=1;
	  m_free_parameters.push_back(i);
	}
	m_nparameters_base +=1;
	m_base_parameters.push_back(i);
	break;

      case statistics::ParameterType::_Derived_:
	m_nparameters_derived += 1;
	m_derived_parameters.push_back(i);
	break;

      default:
	ErrorCBL("Error in cbl::statistics::LikelihoodParameters of LikelihoodParameters.cpp: no such kind of parameter!");
    }
  }
}


// ============================================================================================


cbl::statistics::LikelihoodParameters::LikelihoodParameters (const size_t nparameters, std::vector<ParameterType> parameterTypes, std::vector<std::string> parameterNames) 
{
  set_parameters(nparameters, parameterTypes, parameterNames);
}


// ============================================================================================


size_t cbl::statistics::LikelihoodParameters::nparameters_free () const
{
  return m_nparameters_free;
}


// ============================================================================================


size_t cbl::statistics::LikelihoodParameters::nparameters_fixed () const
{
  return m_nparameters_fixed;
}


// ============================================================================================


void cbl::statistics::LikelihoodParameters::set_parameters (const size_t nparameters, std::vector<ParameterType> parameterTypes, std::vector<std::string> parameterNames)
{
  //Check parameterTypes size
  
  if (nparameters==0)
    ErrorCBL("Error in cbl::statistics::LikelihoodParameters::set_parameters of ModelParameters.cpp! nparameters should be > 0.");

  if ((parameterTypes.size()!=nparameters) && (parameterTypes.size()!=0))
    ErrorCBL("Error in cbl::statistics::LikelihoodParameters::set_parameters of ModelParameters.cpp! Wrong size for the vector parameterTypes.");

  if ((parameterNames.size()!=nparameters) && (parameterNames.size()!=0))
    ErrorCBL("Error in cbl::statistics::LikelihoodParameters::set_parameters of ModelParameters.cpp! Wrong size for the vector parameterNames.");


  if ((parameterTypes.size()==nparameters) && (parameterNames.size()==nparameters)){
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

  m_parameter_bestfit_value.erase(m_parameter_bestfit_value.begin(), m_parameter_bestfit_value.end());
  m_parameter_isFixed.resize(m_nparameters, false);
  m_parameter_fixed_value.resize(m_nparameters, 0);
  m_set_parameter_type();
}


// ============================================================================================


void cbl::statistics::LikelihoodParameters::free (const int p)
{
  switch (m_parameter_type[p]) {

    case ParameterType::_Base_:
      m_parameter_isFixed[p] = false;
      m_set_parameter_type();
      break;

    case ParameterType::_Derived_:
      WarningMsg("Warning in fix of LikelihoodParameters, "+m_parameter_name[p]+" is a derived parameter");
      break;

    default:
      ErrorCBL("Error in cbl::statistics::fix of LikelihoodParameters.cpp: no such kind of parameter!");
  }
}


// ============================================================================================


void cbl::statistics::LikelihoodParameters::fix (const int p, const double value)
{
  switch (m_parameter_type[p]) {

    case statistics::ParameterType::_Base_:
      m_parameter_isFixed[p]=true;
      m_parameter_fixed_value[p] = value;
      m_set_parameter_type();
      break;

    case statistics::ParameterType::_Derived_:
      WarningMsg("Warning in fix of LikelihoodParameters, "+m_parameter_name[p]+" is a derived parameter");
      break;

    default:
      ErrorCBL("Error in cbl::statistics::fix of LikelihoodParameters.cpp: no such kind of parameter!");
  }
}


// ============================================================================================


void cbl::statistics::LikelihoodParameters::fix_at_bestfit (const int p)
{
  fix(p, m_parameter_bestfit_value[p]);
}


// ============================================================================================


double cbl::statistics::LikelihoodParameters::bestfit_value (const int p) const
{
  if (m_parameter_bestfit_value.size() == 0) 
    ErrorCBL("Error in bestfit_values of LikelihoodParameters! Can't found best fit values!"); 

  return m_parameter_bestfit_value[p];
}

// ============================================================================================


vector<double> cbl::statistics::LikelihoodParameters::bestfit_values () const
{
  if (m_parameter_bestfit_value.size() == 0) 
    ErrorCBL("Error in bestfit_values of LikelihoodParameters! Can't found best fit values!"); 

  return m_parameter_bestfit_value;
}


// ============================================================================================


void cbl::statistics::LikelihoodParameters::set_bestfit_value (const vector<double> bestfit_value)
{
  if (bestfit_value.size() != m_nparameters)
    ErrorCBL("Error in set_bestfit_value of LikelihoodParameters! Wrong size for the input vector!"); 

  m_parameter_bestfit_value.erase(m_parameter_bestfit_value.begin(), m_parameter_bestfit_value.end());
  for (size_t i=0; i<m_nparameters; i++) 
    m_parameter_bestfit_value.push_back(bestfit_value[i]);
}


// ============================================================================================


void cbl::statistics::LikelihoodParameters::write_bestfit_info ()
{
  if (m_parameter_bestfit_value.size() == m_nparameters) {
    for (size_t i=0; i<m_nparameters; i++) {

      switch (m_parameter_type[i]) {
	case statistics::ParameterType::_Base_:
	if (m_parameter_isFixed[i]) 
	  coutCBL << "Parameter: " << par::col_yellow << m_parameter_name[i] << par::col_default << " --> status: " << par::col_purple << "FIXED" << endl;
	else 
	  coutCBL << "Parameter: " << par::col_yellow << m_parameter_name[i] << par::col_default << " --> status: " <<  par::col_green << "FREE" << endl;
	break;

	case statistics::ParameterType::_Derived_:
	coutCBL << "Parameter: " << par::col_yellow << m_parameter_name[i] << par::col_default << " --> status: "  << par::col_bred << "OUTPUT" << endl;
	break;

	default:
	ErrorCBL("Error in write_bestfit_info() of LikelihoodParameters.cpp: no such kind of parameter!");
      }

      coutCBL << "value = " << m_parameter_bestfit_value[i] << endl;
      cout << endl;
    }
  }
  else
    ErrorCBL("Error in write_bestfit_info of LikelihoodParameters! Can't found best fit values!"); 

}
