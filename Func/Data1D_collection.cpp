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
 *  @file Func/Data1D_collection.cpp
 *
 *  @brief Methods of the class Data1D_collection
 *
 *  This file contains the implementation of the methods of the class
 *  Data1D_collection
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#include "Data1D_collection.h"

using namespace cosmobl;
using namespace data;


// ======================================================================================


cosmobl::data::Data1D_collection::Data1D_collection (const int ndata)
{
  m_n_data = ndata;
  m_data.resize(ndata);
}


// ======================================================================================


cosmobl::data::Data1D_collection::Data1D_collection (const vector<Data1D> data, const vector<double> x_min, const vector<double> x_max)
{
  m_n_data = data.size();
  m_data = data;

  if (x_min.size() == x_max.size() && int(x_min.size())==m_n_data)
    for (int i=0; i<m_n_data; i++)
      m_data[i].set_limits(x_min[i], x_max[i]);
}


// ======================================================================================


cosmobl::data::Data1D_collection::Data1D_collection (const vector<Data1D> data, const vector<vector<double > > covariance_matrix, const vector<double> x_min, const vector<double> x_max)
{
  m_n_data = data.size();
  m_data = data;
  m_covariance_matrix = covariance_matrix;

  for (int i=0; i<m_n_data; i++)
    m_data[i].set_limits(x_min[i], x_max[i]);
}


// ======================================================================================

      
vector<double> cosmobl::data::Data1D_collection::xx () const
{
  vector<double> xx;

  for (int i=0; i< m_n_data; i++)
    for(int j=0;j< m_data[i].ndata(); j++)
      xx.push_back(m_data[i].xx(j));

  return xx;

}


// ======================================================================================
      
      
vector<double> cosmobl::data::Data1D_collection::fx () const
{
  vector<double> fx;

  for (int i=0; i< m_n_data; i++)
    for(int j=0;j< m_data[i].ndata(); j++)
      fx.push_back(m_data[i].fx(j));

  return fx;
}

// ======================================================================================
      
      
vector<double> cosmobl::data::Data1D_collection::error_fx () const
{
  vector<double> error_fx;

  for (int i=0; i< m_n_data; i++)
    for(int j=0;j< m_data[i].ndata(); j++)
      error_fx.push_back(m_data[i].error_fx(j));

  return error_fx;
}


// ======================================================================================
      
      
void cosmobl::data::Data1D_collection::set_covariance (const string filename, const int skipped_lines)
{
  m_covariance_matrix.erase(m_covariance_matrix.begin(), m_covariance_matrix.end());

  ifstream fin(filename.c_str()); checkIO(fin, filename);

  vector<double> vv;
  m_covariance_matrix.push_back(vv);
  string line; int i = 0;

  for (int i=0; i<skipped_lines; ++i) getline(fin, line);
 
  while (getline(fin, line)) {
    
    stringstream ss(line);
    vector<double> num; double NN = -1.e30;
    while (ss>>NN) num.push_back(NN);

    if (num.size()>=3 && num[2]>-1.e29)
      m_covariance_matrix[i].push_back(num[2]);
    else { i++; m_covariance_matrix.push_back(vv); }

  }

  fin.clear(); fin.close();
  
  m_covariance_matrix.erase(m_covariance_matrix.end()-1, m_covariance_matrix.end());
  
}


// ======================================================================================


void cosmobl::data::Data1D_collection::set_covariance ()
{
  m_covariance_matrix.erase( m_covariance_matrix.begin(), m_covariance_matrix.end());
  m_covariance_matrix.resize(ndata(),vector<double>(ndata(),0));


  for(int d1=0; d1<m_n_data; d1++)
    for(int i=x_down(d1); i<x_up(d1); i++)
      for(int j=x_down(d1); j<x_up(d1); j++)
	m_covariance_matrix[i+m_n_data*d1][j+m_n_data*d1] = m_data[d1].covariance(i,j);
}


// ======================================================================================


void cosmobl::data::Data1D_collection::invert_covariance()
{
  vector<vector<double> > temporary_covariance(ndata_eff(),vector<double>(ndata_eff(),0)), temporary_inverted_covariance(ndata_eff(),vector<double>(ndata_eff(),0));

  for(int d1=0; d1<m_n_data; d1++)
    for(int d2=0; d2<m_n_data; d2++)
      for(int i=x_down(d1); i<x_up(d1); i++)
	for(int j=x_down(d2); j<x_up(d2); j++)
	  temporary_covariance[i-x_down(d1)][j-x_down(d2)] = m_covariance_matrix[i][j];

  invert_matrix (temporary_covariance, temporary_inverted_covariance);  

  m_inverse_covariance_matrix.erase( m_inverse_covariance_matrix.begin(), m_inverse_covariance_matrix.end());
  m_inverse_covariance_matrix.resize(ndata(),vector<double>(ndata(),0));

  for(int d1=0; d1<m_n_data; d1++)
    for(int d2=0; d2<m_n_data; d2++)
      for(int i=x_down(d1); i<x_up(d1); i++)
	for(int j=x_down(d2); j<x_up(d2); j++)
	  m_inverse_covariance_matrix[i][j] = temporary_inverted_covariance[i-x_down(d1)][j-x_down(d2)];
}


// ======================================================================================


int cosmobl::data::Data1D_collection::ndata_eff () const 
{ 
  int ndata_eff = 0; 

  for (int i=0; i<m_n_data; i++)
    ndata_eff += m_data[i].ndata_eff();
    
  return ndata_eff; 
} 


// ======================================================================================


int cosmobl::data::Data1D_collection::ndata () const 
{ 
  int ndata = 0; 

  for (int i=0; i<m_n_data; i++)
    ndata += m_data[i].ndata();
    
  return ndata; 
}  
