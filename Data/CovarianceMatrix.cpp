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
 *  @file Data/CovarianceMatrix.cpp
 *
 *  @brief Methods of the class CovarianceMatrix
 *
 *  This file contains the implementation of the methods of the class
 *  CovarianceMatrix
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unibo.it, alfonso.veropalumbo@unibo.it
 */

#include "Data.h"
#include "EigenWrapper.h"
#include "CovarianceMatrix.h"

using namespace std;

using namespace cbl;
using namespace data;

// ======================================================================================


void cbl::data::CovarianceMatrix::m_set_default ()
{
  m_order = 0;
  m_hartlap_factor = 1;

  m_matrix.resize(0, 0);
  m_precision.resize(0, 0);
  m_correlation.resize(0, 0);
  m_variance.resize(0);
  m_std.resize(0);
  
}


// ======================================================================================


void cbl::data::CovarianceMatrix::m_set (const vector<double> covariance, const double nmeasures, const double prec)
{
  (void)prec;

  m_order = sqrt(covariance.size());

  m_hartlap_factor = hartlap_factor(m_order, nmeasures);

  m_matrix = cbl::wrapper::eigen::SquareMatrixToEigen(covariance);
  m_precision = m_matrix.inverse();
  
  m_variance = m_matrix.diagonal();
  m_std = m_variance.cwiseSqrt();

  m_correlation = m_matrix.cwiseQuotient((m_std * m_std.transpose()));
}


// ======================================================================================


vector<vector<double>> cbl::data::CovarianceMatrix::operator () () const
{
  return cbl::wrapper::eigen::EigenToMatrix(m_matrix);
}


// ======================================================================================


vector<vector<double>> cbl::data::CovarianceMatrix::correlation () const
{
  return cbl::wrapper::eigen::EigenToMatrix(m_correlation);
}


// ======================================================================================


vector<vector<double>> cbl::data::CovarianceMatrix::precision () const
{
  return cbl::wrapper::eigen::EigenToMatrix(m_precision);
}


// ======================================================================================


vector<vector<double>> cbl::data::CovarianceMatrix::precision_hartlap () const 
{
  return cbl::wrapper::eigen::EigenToMatrix(m_hartlap_factor*m_precision);
}


// ======================================================================================


vector<double> cbl::data::CovarianceMatrix::standard_deviation () const 
{
  return cbl::wrapper::eigen::EigenToVector(m_std);
}


// ======================================================================================


vector<double> cbl::data::CovarianceMatrix::variance () const 
{
  return cbl::wrapper::eigen::EigenToVector(m_variance);
}


// ======================================================================================


void cbl::data::CovarianceMatrix::set_from_standard_deviation (const vector<double> standard_deviation, const double nmeasures) 
{
  size_t order = standard_deviation.size();
  vector<double> matrix(order*order);

  for (size_t i=0; i<order; i++)
    matrix[i*order+1] = pow(standard_deviation[i], 2);

  m_set(matrix, nmeasures);
}


// ======================================================================================


void cbl::data::CovarianceMatrix::read (const std::string filename, const int cov_col, const int skipped_lines, const double nmeasures, const double prec) 
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
  
  m_set(covariance, nmeasures, prec);
}


// ======================================================================================


void cbl::data::CovarianceMatrix::write (const string dir, const string file, const int precision, const int rank) const 
{
  (void)rank;
  string file_out = dir+file;
  ofstream fout(file_out.c_str()); checkIO(fout, file_out);

  fout << "### [1] i # [2] j # [3] covariance # [4] correlation " << endl;

  for (size_t i=0; i<m_order; ++i) 
    for (size_t j=0; j<m_order; ++j) 
      fout << setiosflags(ios::fixed) << setprecision(precision) << setw(15) << right << i
	   << "  " << setiosflags(ios::fixed) << setprecision(precision) << setw(15) << right << j
	   << "  " << setiosflags(ios::fixed) << setprecision(precision) << setw(15) << right << m_matrix(i, j)
	   << "  " << setiosflags(ios::fixed) << setprecision(precision) << setw(15) << right << m_correlation(i, j) <<  endl;
   
  fout.close(); cout << endl; coutCBL << "I wrote the file: " << file_out << endl;
}


// ======================================================================================


void cbl::data::CovarianceMatrix::measure (const std::vector<std::shared_ptr<Data>> dataset, const double normalization, const double prec)
{
  const size_t nbins = dataset[0]->ndata();
  const size_t nmocks = dataset.size();

  vector<double> mean(nbins);

  for (size_t i=0; i<nbins; i++)
    for (size_t j=0; j<nmocks; j++)
      mean[i] = dataset[j]->data(i)/nmocks;

  double norm = normalization/(nmocks-1);

  vector<double> covariance;
  for (size_t i=0; i<nbins; i++)
    for (size_t j=0; j<nbins; j++) {
      double val = 0;
      for (size_t n=0; n<nmocks; n++)
	val += (dataset[n]->data(i)-mean[i])*(dataset[n]->data(j)-mean[j])*norm;
      covariance.push_back(val);
    }

  m_set(covariance, nmocks, prec);
}

// ======================================================================================


void cbl::data::CovarianceMatrix::measure(const std::vector<std::vector<std::shared_ptr<Data>>> dataset, const double normalization, const double prec)
{
  size_t nmocks = dataset[0].size();
  size_t nbins = 0;

  for (size_t i=0; i<dataset.size(); i++)
    nbins += dataset[0][i]->ndata();

  m_hartlap_factor = 1-(nbins+1)/(nmocks-1);

  vector<vector<double>> table(nbins, vector<double>(nmocks, 0));
  vector<double> mean(nbins);

  int start=0;
  for (size_t i=0; i<dataset.size(); i++){
    for (size_t j=0; j<dataset[i].size(); j++) {
      for (int k=0; k<dataset[i][j]->ndata(); k++){
	table[start+k][j] = dataset[i][j]->data(k);
	mean[start+k] += dataset[i][j]->data(k)/nmocks;
      }
    }
    start += dataset[i][0]->ndata();
  }

  double norm = normalization/(nmocks-1);

  vector<double> covariance;
  for (size_t i=0; i<nbins; i++)
    for (size_t j=0; j<nbins; j++) {
      double val = 0;
      for (size_t n=0; n<nmocks; n++)
	val += (table[i][n]-mean[i])*(table[j][n]-mean[j])*norm;
      covariance.push_back(val);
    }

  m_set(covariance, nmocks, prec);
}


// ======================================================================================


CovarianceMatrix cbl::data::CovarianceMatrix::cut (const std::vector<bool> mask) const
{
  if (mask.size() != m_order)
    ErrorCBL("size of the mask doesn't match covariance matrix order!", "cut", "CovarianceMatrix.cpp");

  vector<double> matrix;

  for (size_t i=0; i<m_order; i++)
    for (size_t j=0; j<m_order; j++) 
      if (mask[i] and mask[j])
	matrix.push_back(this->operator()(i, j));

  size_t size = sqrt(matrix.size());
  return CovarianceMatrix(reshape(matrix, size, size));
}



// ======================================================================================


cbl::data::CovarianceMatrix cbl::data::CovarianceMatrix::operator += (const cbl::data::CovarianceMatrix covariance) const
{
  const size_t new_size = m_order+covariance.order();
  vector<vector<double>> matrix(new_size, vector<double>(new_size, 0));

  for (size_t i=0; i<m_order; i++)
    for (size_t j=0; j<m_order; j++)
      matrix[i][j] = m_matrix (i, j);

  for (size_t i=m_order; i<new_size; i++)
    for (size_t j=m_order; j<new_size; j++)
      matrix[i][j] = covariance(i-m_order, j-m_order);

  return CovarianceMatrix(matrix);
}


// ======================================================================================


cbl::data::CovarianceMatrix cbl::data::CovarianceMatrix::operator += (const std::shared_ptr<CovarianceMatrix> covariance) const
{
  const size_t new_size = m_order+covariance->order();
  vector<vector<double>> matrix(new_size, vector<double>(new_size, 0));

  for (size_t i=0; i<m_order; i++)
    for (size_t j=0; j<m_order; j++)
      matrix[i][j] = m_matrix (i, j);

  for (size_t i=m_order; i<new_size; i++)
    for (size_t j=m_order; j<new_size; j++)
      matrix[i][j] = covariance->operator()(i-m_order, j-m_order);

  return CovarianceMatrix(matrix);
}


// ======================================================================================


cbl::data::CovarianceMatrix CovarianceMatrix::operator += (const std::vector<CovarianceMatrix> covariance) const 
{
  size_t new_size = m_order;

  for (size_t i=0; i<covariance.size(); i++)
    new_size += covariance[i].order();

  vector<vector<double>> matrix(new_size, vector<double>(new_size, 0));

  for (size_t i=0; i<m_order; i++)
    for (size_t j=0; j<m_order; j++)
      matrix[i][j] = this->operator() (i, j);

  size_t start = m_order;
  size_t stop = m_order;

  for (size_t n=0; n<covariance.size(); n++) {
    stop += covariance[n].order();
    for (size_t i=start; i<stop; i++)
      for (size_t j=start; j<stop; j++)
	matrix[i][j] = covariance[n](i-start, j-start);
    start = stop;
  }

  return CovarianceMatrix(matrix);
}



// ======================================================================================


cbl::data::CovarianceMatrix cbl::data::CovarianceMatrix::operator += (const std::vector<std::shared_ptr<CovarianceMatrix>> covariance) const
{
  size_t new_size = m_order;

  for (size_t i=0; i<covariance.size(); i++)
    new_size += covariance[i]->order();

  vector<vector<double>> matrix(new_size, vector<double>(new_size, 0));

  for (size_t i=0; i<m_order; i++)
    for (size_t j=0; j<m_order; j++)
      matrix[i][j] = this->operator() (i, j);

  size_t start = m_order;
  size_t stop = m_order;

  for (size_t n=0; n<covariance.size(); n++) {
    stop += covariance[n]->order();
    for (size_t i=start; i<stop; i++)
      for (size_t j=start; j<stop; j++)
	matrix[i][j] = covariance[n]->operator()(i-start, j-start);
    start = stop;
  }

  return CovarianceMatrix(matrix);
}
