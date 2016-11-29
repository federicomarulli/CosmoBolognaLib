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
using namespace cosmobl;


// ======================================================================================


cosmobl::statistics::Model::Model (const vector<Parameter> parameters, const shared_ptr<void> fixed_parameters)
  : m_fixed_parameters(fixed_parameters)
{
  for (size_t i=0; i<parameters.size(); i++)
    m_parameters.push_back(move(make_shared<Parameter>(parameters[i])));

  m_npar = m_parameters.size();
  m_npar_eff = npar_eff();
}


// ======================================================================================


vector<double> cosmobl::statistics::Model::parameter_values ()
{  
  vector<double> parameters;
  for (unsigned int i=0; i<m_npar; i++)
    parameters.push_back(m_parameters[i]->value());
  return parameters;
}


// ======================================================================================


vector<double> cosmobl::statistics::Model::parameter_values_from_chain (const int chain, const int position)
{  
  vector<double> parameters;
  
  for (unsigned int i=0; i<m_npar; i++)
    parameters.push_back(parameter(i)->chain(chain)->chain_value(position));

  return parameters;
}


// ======================================================================================


void cosmobl::statistics::Model::set_parameter_values (const vector<double> parameter_values)
{  
  if (parameter_values.size() == m_npar)
    for (unsigned int i=0; i<m_npar; i++) 
      m_parameters[i]->set_value(parameter_values[i]);
  else
    ErrorCBL("Error in set_parameter_values of Model.cpp, size of parameter vector doesn't match the number of model parameters");
}


// ======================================================================================


void cosmobl::statistics::Model::set_chains (const int nchains, const int chain_size)
{
  for (unsigned int i=0; i<m_npar; i++) 
    m_parameters[i]->set_chains(nchains, chain_size);
}


// ======================================================================================


void cosmobl::statistics::Model1D::write_model (const vector<double> xx, const string output_dir, const string output_file, const vector<double> parameter)
{
  if (parameter.size()!=0 && parameter.size()!=m_npar)
    ErrorCBL("Error in Model1D::write_model of Model.cpp: the number of parameters does not match with the internal variable m_npar");
  
  string file_out = output_dir+output_file;
  ofstream fout(file_out.c_str()); checkIO(fout, file_out);
  
  if (parameter.size()==0)
    for (size_t i=0; i<xx.size(); ++i)
      fout << xx[i] << " " << this->operator()(xx[i]) << endl;
  
  else
    for (size_t i=0; i<xx.size(); i++)
      fout << xx[i] << " " << this->operator()(xx[i], parameter) << endl;
  
  fout.clear(); fout.close();
  coutCBL << "I wrote the file: " << file_out << endl;
}


// ======================================================================================


void cosmobl::statistics::Model2D::write_model (const vector<double> xx, vector<double> yy, const string output_dir, const string output_file, const vector<double> parameter)
{
  if (parameter.size()!=0 && parameter.size()!=m_npar)
    ErrorCBL("Error in Model2D::write_model of Model.cpp: the number of parameters does not match with the internal variable m_npar");
  
  string file_out = output_dir+output_file;
  ofstream fout(file_out.c_str()); checkIO(fout, file_out);

  if (parameter.size()==0)
    for (size_t i=0; i<xx.size(); i++)
      for (size_t j=0; j<yy.size(); j++)
	fout << xx[i] << " " << yy[j] << " " << this->operator()(xx[i], yy[j]) << endl;

  else 
    for (size_t i=0; i<xx.size(); i++)
      for (size_t j=0; j<yy.size(); j++)
	fout << xx[i] << " " << yy[j] << " " << this->operator()(xx[i], yy[j], parameter) << endl;

  fout.clear(); fout.close();
  coutCBL << "I wrote the file: " << file_out << endl;
}
