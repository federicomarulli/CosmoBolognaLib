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
 *  @file Data/Data1D.cpp
 *
 *  @brief Methods of the class Data1D
 *
 *  This file contains the implementation of the methods of the class
 *  Data1D
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unibo.it, alfonso.veropalumbo@unibo.it
 */

#include "Data1D.h"

using namespace std;

using namespace cbl;
using namespace data;


// ======================================================================================


void cbl::data::Data1D::set_xx (const vector<double> x)
{
  checkDim(x, ndata(), "x");
  m_x = x;
  m_xsize = ndata();
}


// ======================================================================================


void cbl::data::Data1D::read (const string input_file, const int skip_nlines, const vector<int> column, const vector<int> column_data, const vector<int> column_errors, const vector<int> column_edges)
{
  if (column_data.size()!=column_errors.size()) ErrorCBL("the sizes of column_data and columt_errors must be equal!", "read", "Data1D.cpp");
  
  // default column of input data values
  int column_data_default = column[0]+1;
  
  // default column of input error values
  int column_error_default = column[0]+2;

  // default first column of input edge values
  int column_edges_default = column[0]+3;
  
  ifstream fin(input_file.c_str()); checkIO(fin, input_file);
  string line;
  
  // loop on the data (and, possibly, error) vectors to be read; if
  // more than one data vectors are read, the data (and, possibly,
  // error) vectors will be added one after the other, in order to
  // construct a single data object
  const int cl_max = (column_data.size()<2) ? 1 : column_data.size();
  for (int cl=0; cl<cl_max; ++cl) {
    
    // skip the first nlines (in case of header lines)
    if (skip_nlines>0)
      for (int i=0; i<skip_nlines; ++i)
	getline(fin, line);
   
    // read the file lines
    double last_bin_edge = par::defaultDouble;
    while (getline(fin, line)) {
      
      stringstream ss(line);
      vector<double> num; double NUM;
      while (ss>>NUM) num.emplace_back(NUM);

      if (num.size()<2) ErrorCBL("there are less than 2 columns in the input file!", "read", "Data1D.cpp");
      
      // if column[0] is larger than 0, the column of x coordinates is
      // the one specified in input
      size_t ind_x = max(0, column[0]-1);
      checkDim(num, ind_x+1, "num", false);
      m_x.emplace_back(num[ind_x]);
      
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
      // bin edge values are the ones specified in input
      const size_t ind_edges = (column_edges.size()==0) ? column_edges_default-1 : column_edges[cl]-1;
      if (num.size()>ind_edges+1) {
	m_edges_xx.emplace_back(num[ind_edges]);
	last_bin_edge = num[ind_edges+1];
      }
      
    }

    // add the last bin edge to m_edges_xx
    if (last_bin_edge>par::defaultDouble)
      m_edges_xx.emplace_back(last_bin_edge);

    // go back at the beginning of the file
    fin.clear(); fin.seekg(ios::beg);

    column_data_default += 2;
    column_error_default += 2;
  }
  
  fin.clear(); fin.close();


  // set the data size and diagonal covariance

  m_xsize = m_data.size();
  m_ndata = m_data.size();

  m_covariance.resize(m_ndata, vector<double>(m_ndata, 0));
  for (int i=0; i<m_ndata; i++)
    m_covariance[i][i] = pow(m_error[i], 2);

}


// ======================================================================================


void cbl::data::Data1D::Print (const int precision) const
{
  for (size_t i=0; i<m_x.size(); i++)
    coutCBL << std::fixed << setprecision(precision) << setw(8) << right << m_x[i] 
	 << "   " << std::fixed << setprecision(precision) << setw(8) << right << m_data[i] 
	 << "   " << std::fixed << setprecision(precision) << setw(8) << right << m_error[i] << endl;
}

// ======================================================================================


void cbl::data::Data1D::write (const string dir, const string file, const string header, const int prec, const int ww, int rank) const 
{
  (void)rank;
  
  string file_out = dir+file;
  ofstream fout(file_out.c_str()); checkIO(fout, file_out);

  if (header!=par::defaultString)
    fout << "### " << header << " ###" << endl;

  const int bp = std::cout.precision();

  for (size_t i=0; i<m_x.size(); i++) {
    cbl::Print(m_x[i], prec, ww, "", "  ", false, fout);
    cbl::Print(m_data[i], prec, ww, "", "  ", false, fout);
    cbl::Print(m_error[i], prec, ww, "", "\n", false, fout);
  }
  cout.precision(bp);
  fout.close(); cout << endl; coutCBL << "I wrote the file: " << file_out << endl;
}


// ======================================================================================


void cbl::data::Data1D::write_covariance (const string dir, const string file, const int precision) const 
{
  checkDim(m_covariance, m_ndata, m_ndata, "covariance", false);
  
  string file_out = dir+file;
  ofstream fout(file_out.c_str()); checkIO(fout, file_out);

  fout << "### [1] r1 # [2] r2 # [3] covariance # [4] correlation # [5] index1 # [6] index2 ### " << endl;

  for (int i=0; i<m_ndata; ++i) 
    for (int j=0; j<m_ndata; ++j) 
      fout << setiosflags(ios::fixed) << setprecision(precision) << setw(15) << right << m_x[i]
	   << "  " << setiosflags(ios::fixed) << setprecision(precision) << setw(15) << right << m_x[j]
	   << "  " << setiosflags(ios::fixed) << setprecision(precision) << setw(15) << right << m_covariance[i][j]
	   << "  " << setiosflags(ios::fixed) << setprecision(precision) << setw(15) << right << m_covariance[i][j]/sqrt(m_covariance[i][i]*m_covariance[j][j])
	   << "  " << setiosflags(ios::fixed) << setprecision(precision) << setw(5) << right << i
	   << "  " << setiosflags(ios::fixed) << setprecision(precision) << setw(5) << right << j <<  endl;
   
  fout.close(); cout << endl; coutCBL << "I wrote the file: " << file_out << endl;
}


// ======================================================================================


std::shared_ptr<Data> cbl::data::Data1D::cut (const std::vector<bool> mask) const
{
  checkDim(mask, m_ndata, "mask");

  vector<double> xx, edges_xx;
  int last_index = 0;
  for (size_t i=0; i<mask.size(); i++){
    if (mask[i]){
      xx.push_back(m_x[i]);
      if (m_edges_xx.size() > 0) 
	edges_xx.push_back(m_edges_xx[i]);
      last_index = i;
    }
  }
  if (m_edges_xx.size() > 0) 
    edges_xx.push_back(m_edges_xx[last_index+1]);

  vector<double> data, error;
  vector<vector<double>> covariance;
  Data::cut(mask, data, error, covariance);

  shared_ptr<Data> dd = make_shared<Data1D>(Data1D(xx, data, covariance, -1, edges_xx));

  return dd;
}



// ======================================================================================


shared_ptr<Data> cbl::data::Data1D::cut (const double xmin, const double xmax) const
{
  vector<bool> mask(m_ndata, true);
  for (int i=0; i<m_ndata; i++)
    if ((m_x[i] < xmin) || (m_x[i]>xmax))
      mask[i] = false;

  return cut(mask);
}


// ======================================================================================


shared_ptr<Data> cbl::data::Data1D::as_factory ()
{
  shared_ptr<Data> dd = make_shared<Data1D>(Data1D(m_x, m_data, m_covariance));

  return dd;
}

