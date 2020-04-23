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
 *  @authors federico.marulli3@unibo.it, alfonso.veropalumbo@unibo.it
 */

#include "LikelihoodParameters.h"

using namespace std;

using namespace cbl;


// ============================================================================================


void cbl::statistics::LikelihoodParameters::reset ()
{
  set_parameters(m_nparameters, m_parameter_type, m_parameter_name);
}

// ============================================================================================


string cbl::statistics::LikelihoodParameters::status (const int p) const
{
  string stat;

  switch (m_parameter_type[p]) {
    case statistics::ParameterType::_Base_:
      if (m_parameter_isFixed[p]) 
	stat = "FIXED";
      else 
	stat = "FREE";
      break;

    case statistics::ParameterType::_Derived_:
      stat = "OUTPUT";
      break;

    default:
      ErrorCBL("no such kind of parameter!", "status", "LikelihoodParameters.cpp");
  }

  return stat;
}


// ============================================================================================


vector<string> cbl::statistics::LikelihoodParameters::status () const
{
  vector<string> stat;
  for (size_t i=0; i<m_nparameters; i++) {
    switch (m_parameter_type[i]) {
      case statistics::ParameterType::_Base_:
	if (m_parameter_isFixed[i]) 
	  stat.push_back("FIXED");
	else 
	  stat.push_back("FREE");
	break;

      case statistics::ParameterType::_Derived_:
	stat.push_back("OUTPUT");
	break;

      default:
	ErrorCBL("no such kind of parameter!", "status", "LikelihoodParameters.cpp");
    }
  }
  return stat;
}


// ============================================================================================


vector<double> cbl::statistics::LikelihoodParameters::full_parameter (const vector<double> parameter_value) const
{
  if (parameter_value.size() == m_nparameters_free) {
    vector<double> all_parameter(m_nparameters, 0);

    for(size_t i=0; i<m_nparameters_free; i++)
      all_parameter[m_free_parameter[i]] = parameter_value[i];

    for(size_t i=0; i<m_nparameters_fixed; i++)
      all_parameter[m_fixed_parameter[i]] = m_parameter_fixed_value[m_fixed_parameter[i]];

    for(size_t i=0; i<m_nparameters_derived; i++)
      all_parameter[m_derived_parameter[i]] = 0.;

    return all_parameter;
  }
  else if (parameter_value.size()==m_nparameters)
    return parameter_value;
  else
    ErrorCBL("the size of the vector of free parameters is incorrect!", "full_parameters", "LikelihoodParameters.cpp");

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

  m_base_parameter.erase(m_base_parameter.begin(), m_base_parameter.end());
  m_fixed_parameter.erase(m_fixed_parameter.begin(), m_fixed_parameter.end());
  m_free_parameter.erase(m_free_parameter.begin(), m_free_parameter.end());
  m_derived_parameter.erase(m_derived_parameter.begin(), m_derived_parameter.end());

  for (size_t i=0; i<m_nparameters; i++) {
    switch (m_parameter_type[i]) {
      case statistics::ParameterType::_Base_:
	if (m_parameter_isFixed[i]) {
	  m_nparameters_fixed ++;
	  m_fixed_parameter.push_back(i);
	}
	else {
	  m_nparameters_free ++;
	  m_free_parameter.push_back(i);
	}
	m_nparameters_base ++;
	m_base_parameter.push_back(i);
	break;

      case statistics::ParameterType::_Derived_:
	m_nparameters_derived ++;
	m_derived_parameter.push_back(i);
	break;

      default:
	ErrorCBL("no such kind of parameter!", "m_set_parameter_type", "LikelihoodParameters.cpp");
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
    ErrorCBL("nparameters must be >0!", "set_parameters", "LikelihoodParameters.cpp");

  if ((parameterTypes.size()!=nparameters) && (parameterTypes.size()!=0))
    ErrorCBL("the size of parameterTypes is incorrect!", "set_parameters", "LikelihoodParameters.cpp");

  if ((parameterNames.size()!=nparameters) && (parameterNames.size()!=0))
    ErrorCBL("the size of parameterNames is incorrect!", "set_parameters", "LikelihoodParameters.cpp");


  if ((parameterTypes.size()==nparameters) && (parameterNames.size()==nparameters)) {
    m_nparameters=nparameters;
    m_parameter_type = parameterTypes;
    m_parameter_name = parameterNames;
  }
  else if ((parameterTypes.size()==0) && (parameterNames.size()==0)) {
    m_nparameters=nparameters;
    vector<ParameterType> pTypes(m_nparameters);
    vector<string> pNames(m_nparameters);
    for(size_t i=0; i<m_nparameters; i++) {
      pTypes[i] = ParameterType::_Base_;
      pNames[i] = "par_"+conv(i+1, par::fINT);
    }
    m_parameter_type = pTypes;
    m_parameter_name = pNames;
  }
  else if ((parameterTypes.size()==nparameters) && (parameterNames.size()==0)) {
    m_nparameters=nparameters;
    vector<string> pNames(m_nparameters);
    for(size_t i=0; i<m_nparameters; i++)
      pNames[i] = "par_"+conv(i+1, par::fINT);
    
    m_parameter_type = parameterTypes;
    m_parameter_name = pNames;
  }
  else if ((parameterTypes.size()==0) && (parameterNames.size()==0)) {
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
      WarningMsgCBL(m_parameter_name[p]+" is a derived parameter!", "free", "LikelihoodParameters.cpp");
      break;

    default:
      ErrorCBL("no such kind of parameter!", "free", "LikelihoodParameters.cpp");
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
      WarningMsgCBL(m_parameter_name[p]+" is a derived parameter!", "fix", "LikelihoodParameters.cpp");
      break;

    default:
      ErrorCBL("no such kind of parameter!", "fix", "LikelihoodParameters.cpp");
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
    ErrorCBL("the best-fit values have not been computed!", "bestfit_value", "LikelihoodParameters.cpp"); 

  return m_parameter_bestfit_value[p];
}

// ============================================================================================


vector<double> cbl::statistics::LikelihoodParameters::bestfit_value () const
{
  if (m_parameter_bestfit_value.size()==0) 
    ErrorCBL("the best-fit values have not been computed!", "bestfit_value", "LikelihoodParameters.cpp"); 

  return m_parameter_bestfit_value;
}


// ============================================================================================


void cbl::statistics::LikelihoodParameters::set_bestfit_values (const vector<double> bestfit_value)
{
  if (bestfit_value.size() != m_nparameters)
    ErrorCBL("the size of the input vector is incorrect!", "set_bestfit_values", "LikelihoodParameters.cpp"); 

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
	ErrorCBL("no such kind of parameter!", "write_bestfit_info", "LikelihoodParameters.cpp");
      }

      Print(m_parameter_bestfit_value[i], 5, 10, "value = ", "\n", true, std::cout); 
      cout << endl;
    }
  }
  else
    ErrorCBL("the best-fit values have not been computed!", "write_bestfit_info", "LikelihoodParameters.cpp"); 

}
