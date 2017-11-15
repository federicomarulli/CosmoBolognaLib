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
 *  @file Func/Data.cpp
 *
 *  @brief Methods of the class Data
 *
 *  This file contains the implementation of the methods of the class
 *  Data
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#include "Data.h"
#include "Data1D.h"
#include "Data2D.h"
#include "Data1D_collection.h"
#include "Data1D_extra.h"
#include "Data2D_extra.h"

using namespace cosmobl;
using namespace data;


// ======================================================================================


shared_ptr<Data> cosmobl::data::Data::Create (const DataType dataType)
{
  if (dataType==DataType::_1D_data_) return move(unique_ptr<Data1D>(new Data1D()));
  else if (dataType==DataType::_2D_data_) return move(unique_ptr<Data2D>(new Data2D()));
  else if (dataType==DataType::_1D_data_collection_) return move(unique_ptr<Data1D_collection>(new Data1D_collection()));
  else if (dataType==DataType::_1D_data_extra_) return move(unique_ptr<Data1D_extra>(new Data1D_extra()));
  else if (dataType==DataType::_2D_data_extra_) return move(unique_ptr<Data2D_extra>(new Data2D_extra()));
  else ErrorCBL("Error in cosmobl::data::Data::Create of Data.cpp: no such type of object, or error in the input parameters!!");

  return NULL;
}


// ======================================================================================


void cosmobl::data::Data::reset (const int ndata) 
{
  m_ndata = ndata;
  m_data.resize(m_ndata, 0);
  m_error.resize(m_ndata, 0);
  m_covariance.resize(m_ndata, vector<double>(m_ndata, 0));
  m_inverse_covariance.resize(m_ndata, vector<double>(m_ndata, 0));
}


// ======================================================================================


cosmobl::data::Data::Data (const DataType dataType, const vector<double> data) 
{
  m_dataType = dataType;
  reset(data.size());

  set_data(data);
}


// ======================================================================================


cosmobl::data::Data::Data (const DataType dataType, const vector<double> data, const vector<double> error) 
{
  m_dataType = dataType;
  reset(data.size());

  set_data(data);
  set_error(error);
  set_covariance(error);
}


// ======================================================================================


cosmobl::data::Data::Data (const DataType dataType, const vector<double> data, const vector<vector<double>> covariance) 
{
  m_dataType = dataType;
  reset(data.size());

  set_data(data);
  set_covariance(covariance);
  set_error(covariance);
}


// ======================================================================================


vector<vector<double>> cosmobl::data::Data::correlation() const
{
  vector<vector<double>> corr(m_ndata, vector<double>(m_ndata,0));

  for (int i=0; i<m_ndata; i++)
    for (int j=0; j<m_ndata; j++)
      corr[i][j] = m_covariance[i][j]/sqrt(m_covariance[i][i]*m_covariance[j][j]);

  return corr;
}


// ======================================================================================


void cosmobl::data::Data::set_data(const vector<double> data)
{
  checkDim(data, m_ndata, "data");
  m_data = data;
}


// ======================================================================================


void cosmobl::data::Data::set_error(const vector<double> error)
{
  checkDim(error, m_ndata, "error");
  m_error = error;
}


// ======================================================================================


void cosmobl::data::Data::set_error(const vector<vector<double>> covariance)
{
  checkDim(covariance, m_ndata, m_ndata, "covariance");

  for (int i=0; i<m_ndata; i++)
    m_error[i]=sqrt(covariance[i][i]);
}



// ======================================================================================


void cosmobl::data::Data::set_covariance (const vector<double> error)
{
  checkDim(error, m_ndata, "error");

  for (int i=0; i<m_ndata; i++)
    m_covariance[i][i]=pow(m_error[i],2);

  invert_covariance();
}


// ======================================================================================


void cosmobl::data::Data::set_covariance (const vector<vector<double>> covariance)
{
  checkDim(covariance, m_ndata, m_ndata, "covariance");

  m_covariance = covariance;
  set_error(covariance);

  invert_covariance();
}


// ======================================================================================


void cosmobl::data::Data::set_covariance (const string filename, const int cov_col, const int skipped_lines)
{  
  ifstream fin(filename.c_str()); checkIO(fin, filename);
  string line;

  for (int i=0; i<skipped_lines; ++i) getline(fin, line);

  vector<double> covariance;
  
  while (getline(fin, line)) {
    
    stringstream ss(line);
    vector<double> num; double NN = par::defaultDouble;
    while (ss>>NN) num.push_back(NN);
    
    if (int(num.size())>=cov_col && num[cov_col]>par::defaultDouble) 
      covariance.push_back(num[cov_col]);
  }

  fin.clear(); fin.close();

  m_covariance = reshape(covariance, m_ndata, m_ndata);
  
  set_covariance(m_covariance);
  set_error(m_covariance);

  invert_covariance();
}


// ======================================================================================


void cosmobl::data::Data::cut(const vector<bool> mask, vector<double> &data, vector<double> &error, vector<vector<double>> &covariance_matrix) const
{
  checkDim (mask, m_ndata, "mask");

  int ndata_eff = 0;

  for (int i=0; i<m_ndata; i++)
    if (mask[i])
      ndata_eff +=1;

  if (ndata_eff <1)
    ErrorCBL("Error in cut of Data, no elements left");

  data.resize(ndata_eff, 0);
  error.resize(ndata_eff, 0);
  covariance_matrix.resize(ndata_eff, vector<double>(ndata_eff, 0));

  int index1 = 0;
  for (int i=0; i<m_ndata; i++) {
    if (mask[i]) {
      data[index1] = m_data[i];
      error[index1] = m_data[i];

      int index2 = 0;
      for (int j=0; j<m_ndata; j++) {
	if (mask[j]) {
	  covariance_matrix[index1][index2] = m_covariance[i][j];
	  index2++;
	}
      }
      index1++;
    }
  }
}


// ============================================================================


shared_ptr<data::Data> cosmobl::data::join_dataset (vector<shared_ptr<data::Data>> dataset)
{
  if (dataset.size()<2)
    cosmobl::ErrorCBL("Error in join_dataset(). You must provide at least 2 dataset");

  cosmobl::data::DataType dt = dataset[0]->dataType();

  for (size_t i=0; i<dataset.size(); i++)
    if (dt!=dataset[i]->dataType())
      cosmobl::ErrorCBL("Error in join_dataset(). Dataset must be equal");

  switch(dt) {

    case cosmobl::data::DataType::_1D_data_:
      return join_dataset_1D(dataset);
      break;

    case cosmobl::data::DataType::_1D_data_extra_:
      return join_dataset_1D_extra(dataset);
      break;

    default:
      ErrorCBL("Error in join_dataset(). Work in progress!");

  }

  return NULL;
}


// ============================================================================


shared_ptr<data::Data> cosmobl::data::join_dataset_1D(vector<shared_ptr<data::Data>> dataset)
{
  int ndataset = (int)dataset.size();
  int ndata = 0;
  vector<vector<int>> data_index;

  int nn = 0;
  for (int i=0; i<ndataset; i++) {
    ndata+=dataset[i]->ndata();
    vector<int> _dd;
    for (int j=0; j<dataset[i]->ndata(); j++) {
      _dd.push_back(nn);
      nn++;
    }
    data_index.push_back(_dd);
  }

  vector<double> xx(ndata, 0), data(ndata, 0);
  vector<vector<double>> covariance(ndata, vector<double>(ndata, 0));

  for (int i=0; i<ndataset;i++)
    for (int j=0; j<dataset[i]->ndata(); j++) {
      xx[data_index[i][j]] = dataset[i]->xx(j);
      data[data_index[i][j]] = dataset[i]->data(j);

      for (int k=0; k<dataset[i]->ndata(); k++)
	covariance[data_index[i][j]][data_index[i][k]] = dataset[i]->covariance(j, k);
    }

  return move(unique_ptr<cosmobl::data::Data1D>(new cosmobl::data::Data1D(xx, data, covariance)));
}


// ============================================================================


shared_ptr<data::Data> data::join_dataset_1D_extra (vector<shared_ptr<data::Data>> dataset)
{
  ErrorCBL("Error in join_dataset_1D_extra, work in progress!");
  int ndataset = (int)dataset.size();
  int ndata = 0;
  vector<int> n_data(ndataset);
  int nextra = dataset[0]->extra_info().size();

  for (int i=0; i<ndataset; i++)
    ndata += dataset[i]->ndata();

  vector<double> xx(ndata, 0), data(ndata, 0);
  vector<vector<double>> extra_info (nextra, vector<double>(ndata, 0));
  vector<vector<double>> covariance(ndata, vector<double>(ndata, 0));

  auto merged_dataset = make_shared<cosmobl::data::Data1D_extra>(cosmobl::data::Data1D_extra(xx, data, covariance, extra_info));
  
  return merged_dataset;
}
