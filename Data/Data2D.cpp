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
 *  @file Data/Data2D.cpp
 *
 *  @brief Methods of the class Data2D
 *
 *  This file contains the implementation of the methods of the class
 *  Data2D
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unibo.it, alfonso.veropalumbo@unibo.it
 */

#include "Data2D.h"

using namespace std;

using namespace cbl;
using namespace data;


// ======================================================================================


cbl::data::Data2D::Data2D (const std::vector<double> x, const std::vector<double> y, const std::vector<std::vector<double>> data, const std::vector<double> bin_edges_x, const std::vector<double> bin_edges_y) : Data(DataType::_2D_)
{
  m_x = x;
  m_y = y;
  if (bin_edges_x.size()>0) 
    m_edges_xx = bin_edges_x;
  if (bin_edges_y.size()>0) 
    m_edges_yy = bin_edges_y;

  m_xsize = m_x.size();
  m_ysize = m_y.size();

  checkDim(data, m_xsize, "data");
  
  for (int i=0; i<m_xsize; i++) {
    checkDim(data[i], m_ysize, "data["+conv(i,par::fINT)+"]");
    for (int j=0; j<m_ysize; j++)
      m_data.push_back(data[i][j]);
  }

  m_ndata = m_data.size();
  m_error.resize(m_ndata,0);

  m_covariance.resize(m_ndata, vector<double>(m_ndata,0));
}


// ======================================================================================


cbl::data::Data2D::Data2D (const std::vector<double> x, const std::vector<double> y, const std::vector<std::vector<double>> data, const std::vector<std::vector<double>> error, const std::vector<double> bin_edges_x, const std::vector<double> bin_edges_y) : Data(DataType::_2D_)
{
  m_x = x;
  m_y = y;
  if (bin_edges_x.size()>0) 
    m_edges_xx = bin_edges_x;
  if (bin_edges_y.size()>0) 
    m_edges_yy = bin_edges_y;

  m_xsize = m_x.size();
  m_ysize = m_y.size();

  checkDim(data, m_xsize, "data");

  checkDim(error, m_xsize, "error");
  
  for (int i=0; i<m_xsize; i++) {
    checkDim(data[i], m_ysize, "data["+conv(i, par::fINT)+"]");
    checkDim(error[i], m_ysize, "error["+conv(i, par::fINT)+"]");
    for (int j=0; j<m_ysize; j++) {
      m_data.push_back(data[i][j]);
      m_error.push_back(error[i][j]);
    }
  }

  m_ndata = m_data.size();

  m_covariance.resize(m_ndata, vector<double>(m_ndata, 0));
  for (int i=0; i<m_ndata; i++)
    m_covariance[i][i] = pow(m_error[i], 2);
}


// ======================================================================================


cbl::data::Data2D::Data2D (const std::vector<double> x, const std::vector<double> y, const std::vector<double> data, const std::vector<double> error, const std::vector<double> bin_edges_x, const std::vector<double> bin_edges_y) : Data(DataType::_2D_)
{
  m_x = x;
  m_y = y;
  if (bin_edges_x.size()>0) 
    m_edges_xx = bin_edges_x;
  if (bin_edges_y.size()>0) 
    m_edges_yy = bin_edges_y;

  m_xsize = m_x.size();
  m_ysize = m_y.size();

  checkDim(data, m_xsize*m_ysize, "data");
  checkDim(error, m_xsize*m_ysize, "error");

  m_data = data;
  m_error = error;

  m_ndata = m_data.size();

  m_covariance.resize(m_ndata, vector<double>(m_ndata,0));
  for (int i=0; i<m_ndata; i++)
    m_covariance[i][i] = pow(m_error[i], 2);
}


// ======================================================================================


cbl::data::Data2D::Data2D (const std::vector<double> x, const std::vector<double> y, const std::vector<double> data, const std::vector<std::vector<double>> covariance_matrix, const std::vector<double> bin_edges_x, const std::vector<double> bin_edges_y) : Data(DataType::_2D_)
{
  m_x = x;
  m_y = y;
  if (bin_edges_x.size()>0) 
    m_edges_xx = bin_edges_x;
  if (bin_edges_y.size()>0) 
    m_edges_yy = bin_edges_y;

  m_xsize = m_x.size();
  m_ysize = m_y.size();
  
  checkDim(data, m_xsize*m_ysize, "data");
  
  m_data = data;

  m_ndata = m_data.size();
  m_error.resize(m_ndata,0);

  checkDim(covariance_matrix, m_ndata, "covariance_matrix");
  for (int i=0;i<m_ndata;i++)
    checkDim(covariance_matrix[i], m_ndata, "covariance_matrix["+conv(i,par::fINT)+"]");

  m_covariance = covariance_matrix;

  for (int i=0;i<m_ndata;i++)
    m_error[i] = sqrt(m_covariance[i][i]);
}


// ======================================================================================


void cbl::data::Data2D::get_data (std::vector<std::vector<double>> &data) const
{
  data.erase(data.begin(), data.end());
  data.resize(m_xsize, vector<double>(m_ysize,0));

  for (int i=0; i<m_xsize; i++)
    for (int j=0; j< m_ysize; j++)
      data[i][j] = this->data(i,j);
}


// ======================================================================================


void cbl::data::Data2D::get_error (std::vector<std::vector<double>> &error) const
{
  error.erase(error.begin(), error.end());
  error.resize(m_xsize, vector<double>(m_ysize,0));

  for (int i=0; i<m_xsize; i++)
    for (int j=0; j< m_ysize; j++) 
      error[i][j] = this->error(i,j);
}


// ======================================================================================


void cbl::data::Data2D::read (const string input_file, const int skip_nlines, const vector<int> column, const vector<int> column_data, const vector<int> column_errors, const vector<int> column_edges)
{
  if (column_data.size()!=column_errors.size()) ErrorCBL("the sizes of column_data and columt_errors must be equal!", "read", "Data2D.cpp");

  // default column of input data values in the y dimension
  int column_data_default = column[1]+1;
  
  // default column of input error values
  int column_error_default = column[1]+2;

  // default first column of input edge values
  int column_edges_default = column[1]+3;
  
  ifstream fin(input_file.c_str()); checkIO(fin, input_file);
  string line;
  
  // loop on the data (and, possibly, error) vectors to be read; if
  // more than one data vectors are read, the data (and, possibly,
  // error) vectors will be added one after the other, in order to
  // construct a single data object
  const int cl_max = (column_data.size()<2) ? 1 : column_data.size();
  std::vector<double> dummy_edges_xx, dummy_edges_yy;
  for (int cl=0; cl<cl_max; ++cl) {
    
    // skip the first nlines (in case of header lines)
    if (skip_nlines>0)
      for (int i=0; i<skip_nlines; ++i)
	getline(fin, line);

    
    // read the file lines

    double last_bin_edge_xx = par::defaultDouble;
    double last_bin_edge_yy = par::defaultDouble;
    
    while (getline(fin, line)) {
      
      stringstream ss(line);
      vector<double> num; double NUM;
      while (ss>>NUM) num.emplace_back(NUM);

      if (num.size()<3) ErrorCBL("there are less than 3 columns in the input file!", "read", "Data2D.cpp");
      
      // if column[0] is larger than 0, the column of x coordinates is
      // the one specified in input
      size_t ind_x = max(0, column[0]-1);
      checkDim(num, ind_x+1, "num", false);
      m_x.emplace_back(num[ind_x]);

      // if column[1] is larger than 0, the column of y coordinates is
      // the one specified in input
      size_t ind_y = max(0, column[1]-1);
      checkDim(num, ind_y+1, "num", false);
      m_y.emplace_back(num[ind_y]);

      // if the size of the column_data vector is 0, the columns of
      // data values are the ones specified in input
      const size_t ind_data = (column_data.size()==0) ? column_data_default-1 : column_data[cl]-1;
      checkDim(num, ind_data+1, "num", false);
      m_data.emplace_back(num[ind_data]);

      // if the size of the column_errors vector is 0, the columns of
      // error values are the ones specified in input; if the error
      // column is not present, the errors will be set to 1, by
      // default
      const size_t ind_error = (column_errors.size()==0) ? column_error_default-1 : column_errors[cl]-1;
      if (num.size()<ind_error+1 && column_errors.size()>1)
	WarningMsgCBL("the errors cannot be retrieved from the provided input file, and will be set to 1", "read", "Data1D.cpp");
      m_error.emplace_back((num.size()<ind_error+1) ? 1 : num[ind_error]);

      // if the size of the column_edges vector is 0, the columns of
      // bin edges values are the ones specified in input
      const size_t ind_edges = (column_edges.size()==0) ? column_edges_default-1 : column_edges[cl]-1;
      if (num.size()>ind_edges+1) {
	dummy_edges_xx.emplace_back(num[ind_edges]);
	last_bin_edge_xx = num[ind_edges+1];
      }
      
      // if the size of the column_edges vector is 0, the columns of
      // bin edges values are the ones specified in input
      const size_t ind_edges_yy = (column_edges.size()==0) ? column_edges_default-1+2 : column_edges[cl]-1+2;
      if (num.size()>ind_edges_yy+1) {
	dummy_edges_yy.emplace_back(num[ind_edges_yy]);
	last_bin_edge_yy = num[ind_edges_yy+1];
      }
      
    }

    std::sort( dummy_edges_xx.begin(), dummy_edges_xx.end() );
    dummy_edges_xx.erase( std::unique( dummy_edges_xx.begin(), dummy_edges_xx.end() ), dummy_edges_xx.end() );
    std::sort( dummy_edges_yy.begin(), dummy_edges_yy.end() );
    dummy_edges_yy.erase( std::unique( dummy_edges_yy.begin(), dummy_edges_yy.end() ), dummy_edges_yy.end() );

    for (size_t i=0; i<dummy_edges_xx.size(); i++)
      m_edges_xx.emplace_back(dummy_edges_xx[i]);
    m_edges_xx.emplace_back(last_bin_edge_xx);
    for (size_t i=0; i<dummy_edges_yy.size(); i++)
      m_edges_yy.emplace_back(dummy_edges_yy[i]);
    m_edges_yy.emplace_back(last_bin_edge_yy);

    // go back at the beginning of the file
    fin.clear(); fin.seekg(ios::beg);

    column_data_default += 2;
    column_error_default += 2;
  }
  
  fin.clear(); fin.close();

  
  // set the data size and diagonal covariance

  unique_unsorted(m_x);
  unique_unsorted(m_y);
  
  m_xsize = m_x.size();
  m_ysize = m_y.size();
  m_ndata = m_data.size();

  m_covariance.resize(m_ndata, vector<double>(m_ndata, 0));
  for (int i=0; i<m_ndata; i++)
    m_covariance[i][i] = pow(m_error[i], 2);
}


// ======================================================================================


void cbl::data::Data2D::Print (const int precision) const 
{
  for (int i=0; i<m_xsize; ++i)
    for (int j=0; j<m_ysize; ++j) {
      int index = j+m_ysize*i;
      coutCBL << setprecision(precision) << setw(8) << right << m_x[i]
	      << "   " << setprecision(precision) << setw(8) << right << m_y[j]
	      << "   " << setprecision(precision) << setw(8) << right << m_data[index]
	      << "   " << setprecision(precision) << setw(8) << right << m_error[index] << endl;
    }
}


// ======================================================================================


void cbl::data::Data2D::write (const string dir, const string file, const string header, const bool full, const int prec, const int ww, const int rank) const 
{
  (void)rank;
  
  string file_out = dir+file;
  ofstream fout(file_out.c_str()); checkIO(fout, file_out);

  if (header!=par::defaultString)
    fout << "### " << header << " ###" << endl;

  const int bp = cout.precision();
  
  for (int i=0; i<m_xsize; ++i)
    for (int j=0; j<m_ysize; ++j) {
      int index = j+m_ysize*i;
      cbl::Print(m_x[i], prec, ww, "", "  ", false, fout);
      cbl::Print(m_y[j], prec, ww, "", "  ", false, fout);
      cbl::Print(m_data[index], prec, ww, "", "  ", false, fout);
      cbl::Print(m_error[index], prec, ww, "", "\n", false, fout);
    }


  if (full) { // duplicate the information in the other 3 quadrants

    for (int i=0; i<m_xsize; ++i)
      for (int j=0; j<m_ysize; ++j) {
	int index = j+m_ysize*i;
	cbl::Print(m_x[i], prec, ww, "", "  ", false, fout);
	cbl::Print(-m_y[j], prec, ww, "", "  ", false, fout);
	cbl::Print(m_data[index], prec, ww, "", "  ", false, fout);
	cbl::Print(m_error[index], prec, ww, "", "\n", false, fout);
      }

    for (int i=0; i<m_xsize; ++i)
      for (int j=0; j<m_ysize; ++j) {
	int index = j+m_ysize*i;
	cbl::Print(-m_x[i], prec, ww, "", "  ", false, fout);
	cbl::Print(-m_y[j], prec, ww, "", "  ", false, fout);
	cbl::Print(m_data[index], prec, ww, "", "  ", false, fout);
	cbl::Print(m_error[index], prec, ww, "", "\n", false, fout);
      }

    for (int i=0; i<m_xsize; ++i)
      for (int j=0; j<m_ysize; ++j) {
	int index = j+m_ysize*i;
	cbl::Print(-m_x[i], prec, ww, "", "  ", false, fout);
	cbl::Print(m_y[j], prec, ww, "", "  ", false, fout);
	cbl::Print(m_data[index], prec, ww, "", "  ", false, fout);
	cbl::Print(m_error[index], prec, ww, "", "\n", false, fout);
      }
    
  }
  
  cout.precision(bp);
  fout.close(); cout << endl; coutCBL << "I wrote the file: " << file_out << endl << endl;
}


// ======================================================================================


void cbl::data::Data2D::write_covariance (const string dir, const string file, const int precision) const
{
  checkDim(m_covariance, m_ndata, m_ndata, "covariance", false);
  
  string file_out = dir+file;
  ofstream fout(file_out.c_str()); checkIO(fout, file_out);

  fout << "### [1] r1 # [2] r2 # [3] r1 # [4] r2 # [5] covariance # [6] correlation # [7] index1 # [8] index2 ### " << endl;

  for (int i=0; i<m_xsize; ++i) {
    for (int j=0; j<m_ysize; ++j) { 
      for (int k=0; k<m_xsize; ++k) {
	for (int l=0; l<m_ysize; ++l) {    
	  int index1 = j+m_ysize*i;
	  int index2 = l+m_ysize*k;
	  fout << setiosflags(ios::fixed) << setprecision(precision) << setw(15) << right << m_x[i]
	       << "  " << setiosflags(ios::fixed) << setprecision(precision) << setw(15) << right << m_y[j]
	       << "  " << setiosflags(ios::fixed) << setprecision(precision) << setw(15) << right << m_x[k]
	       << "  " << setiosflags(ios::fixed) << setprecision(precision) << setw(15) << right << m_y[l]
	       << "  " << setiosflags(ios::fixed) << setprecision(precision) << setw(15) << right << m_covariance[index1][index2]
	       << "  " << setiosflags(ios::fixed) << setprecision(precision) << setw(15) << right << m_covariance[index1][index2]/sqrt(m_covariance[index1][index1]*m_covariance[index2][index2])
	       << "  " << setiosflags(ios::fixed) << setprecision(precision) << setw(5) << right << index1 
	       << "  " << setiosflags(ios::fixed) << setprecision(precision) << setw(5) << right << index2 <<  endl;
	}
      }
    }
  }

  fout.close(); cout << endl; coutCBL << "I wrote the file: " << file_out << endl;

}


// ======================================================================================


shared_ptr<Data> cbl::data::Data2D::cut (const double xmin, const double xmax, const double ymin, const double ymax) const
{
  vector<bool> mask(m_ndata, true);
  vector<double> xx, yy, edges_xx, edges_yy;

  for (int i=0; i<m_xsize; i++) {
    for (int j=0; j<m_ysize; j++) {
      if ( ((m_x[i] < xmin) || (m_x[i]>xmax)) || ((m_y[j] < ymin) || (m_y[j]>ymax)))
	mask[j+m_ysize*i] = false;
    }
  }

  int last_index = 0;
  for (int i=0; i<m_xsize; i++){
    if ( !((m_x[i] < xmin) || (m_x[i]>xmax))){
      xx.push_back(m_x[i]);
      if (m_edges_xx.size()>0)
	edges_xx.push_back(m_edges_xx[i]);
      last_index = i;
    }
  }
  if (m_edges_xx.size()>0)
    edges_xx.push_back(m_edges_xx[last_index+1]);

  last_index = 0;
  for (int i=0; i<m_ysize; i++){
    if (!((m_y[i] < ymin) || (m_y[i]>ymax))){
      yy.push_back(m_y[i]);
      if (m_edges_yy.size()>0)
	edges_yy.push_back(m_edges_yy[i]);
      last_index = i;
    }
  }
  if (m_edges_yy.size()>0)
    edges_yy.push_back(m_edges_yy[last_index+1]);

  vector<double> data, error;
  vector<vector<double>> covariance;
  Data::cut(mask, data, error, covariance);

  shared_ptr<Data> dd = make_shared<Data2D>(Data2D(xx, yy, data, covariance, edges_xx, edges_yy));

  return dd;
}


// ======================================================================================


shared_ptr<Data> cbl::data::Data2D::as_factory ()
{
  shared_ptr<Data> dd = make_shared<Data2D>(Data2D(m_x, m_data, m_covariance));

  return dd;
}
