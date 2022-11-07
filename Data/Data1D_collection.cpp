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
 *  @file Data/Data1D_collection.cpp
 *
 *  @brief Methods of the class Data1D_collection
 *
 *  This file contains the implementation of the methods of the class
 *  Data1D_collection
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unibo.it, alfonso.veropalumbo@unibo.it
 */

#include "Data1D_collection.h"

using namespace std;

using namespace cbl;
using namespace data;


// ======================================================================================


cbl::data::Data1D_collection::Data1D_collection (const std::string input_file, const int skip_header)
{
  set_dataType(DataType::_1D_collection_);
  read(input_file, skip_header);
}


// ======================================================================================


cbl::data::Data1D_collection::Data1D_collection (const std::vector<std::string> input_files, const int skip_header)
{
  set_dataType(DataType::_1D_collection_);
  read(input_files, skip_header);
}


// ======================================================================================


cbl::data::Data1D_collection::Data1D_collection (const std::vector<std::vector<double>> xx, const std::vector<std::vector<double>> data)
{
  set_dataType(DataType::_1D_collection_);

  m_x = xx;
  m_ndataset = m_x.size();
  m_xsize.resize(m_ndataset,0);

  checkDim(data, m_ndataset, "data");

  for (int i=0; i<m_ndataset; i++) {
    m_xsize[i] = m_x[i].size();
    checkDim(data[i], m_xsize[i], "data["+conv(i,par::fINT)+"]");
    for (int j=0; j<m_xsize[i]; j++)
      m_data.push_back(data[i][j]);
  }

  m_ndata = m_data.size();
  m_error.resize(m_ndata,0);

  int index = 0;
  for (int i=0; i<m_ndataset; i++) {
    vector<int> vv;
    for (int j=0; j<m_xsize[i]; j++) {
      vv.push_back(index);
      index++;
    }
    m_index.push_back(vv);
  }
}


// ======================================================================================


cbl::data::Data1D_collection::Data1D_collection (const std::vector<std::vector<double>> xx, const std::vector<std::vector<double>> data, const std::vector<std::vector<double>> covariance_matrix)
{
  set_dataType(DataType::_1D_collection_);

  m_x = xx;
  m_ndataset = m_x.size();
  m_xsize.resize(m_ndataset,0);

  checkDim(data, m_ndataset, "data");

  for (int i=0; i<m_ndataset; i++) {
    m_xsize[i] = m_x[i].size();
    checkDim(data[i], m_xsize[i], "data["+conv(i, par::fINT)+"]");
    for (int j=0; j<m_xsize[i]; j++)
      m_data.push_back(data[i][j]);
  }

  m_ndata = m_data.size();
  m_error.resize(m_ndata,0);

  checkDim (covariance_matrix, m_ndata, "covariance");
  for (int i=0; i<m_ndata; i++)
    checkDim(covariance_matrix[i], m_ndata, "covariance["+conv(i,par::fINT)+"]");

  m_covariance = covariance_matrix;

  for (int i=0; i<m_ndata; i++)
    m_error[i] = sqrt(m_covariance[i][i]);

  int index = 0;
  for (int i=0; i<m_ndataset; i++) {
    vector<int> vv;
    for (int j=0; j<m_xsize[i]; j++) {
      vv.push_back(index);
      index++;
    }
    m_index.push_back(vv);
  }
}


// ======================================================================================


cbl::data::Data1D_collection::Data1D_collection (const std::vector<std::vector<double>> xx, const std::vector<double> data, const std::vector<double> error)
{
  set_dataType(DataType::_1D_collection_);

  m_x = xx;
  m_ndataset = m_x.size();
  m_xsize.resize(m_ndataset,0);

  int nn = 0;
  
  for (int i=0; i<m_ndataset; i++) {
    m_xsize[i] = m_x[i].size();
    nn += m_x[i].size();
  }

  checkDim(data, nn, "data");

  m_data = data;
  m_error = error;
  m_ndata = m_data.size();

  m_covariance.resize(m_ndata, vector<double>(m_ndata, 0));
  for (int i=0; i<m_ndata; i++)
    m_covariance[i][i] = pow(m_error[i],2);
}


// ======================================================================================


cbl::data::Data1D_collection::Data1D_collection (const std::vector<std::vector<double>> xx, const std::vector<double> data, const std::vector<std::vector<double>> covariance_matrix)
{
  set_dataType(DataType::_1D_collection_);

  m_x = xx;
  m_ndataset = m_x.size();
  m_xsize.resize(m_ndataset,0);

  int nn = 0;
  
  for (int i=0; i<m_ndataset; i++) {
    m_xsize[i] = m_x[i].size();
    nn += m_x[i].size();
  }

  checkDim(data, nn, "data");

  m_data = data;
  m_covariance = covariance_matrix;

  m_ndata = m_data.size();


  m_error.resize(m_ndata, 0);
  for (int i=0; i<m_ndata; i++)
    m_error[i] = sqrt(m_covariance[i][i]);

}


// ======================================================================================


void cbl::data::Data1D_collection::data (std::vector<std::vector<double>> &data) const
{
  data.erase(data.begin(), data.end());

  for (int i=0; i<m_ndataset; i++) {
    vector<double> vv(m_xsize[i], 0);
    for (int j=0; j< m_xsize[i]; j++) {
      vv[j] = this->data(i,j);
    }
    data.push_back(vv);
  }
}


// ======================================================================================


void cbl::data::Data1D_collection::error (std::vector<std::vector<double>> &error) const
{
  error.erase(error.begin(), error.end());

  for (int i=0; i<m_ndataset; i++) {
    vector<double> vv(m_xsize[i], 0);
    for (int j=0; j< m_xsize[i]; j++) {
      vv[j] = this->error(i, j);
    }
    error.push_back(vv);
  }
}


// ======================================================================================


void cbl::data::Data1D_collection::read (const std::string input_file, const int skip_nlines, const std::vector<int> column, const std::vector<int> column_data, const std::vector<int> column_errors, const std::vector<int> column_edges)
{
  (void)column; (void)column_data; (void)column_errors; (void)column_edges; 
  WarningMsgCBL("column, column_data, column_errors, column_edges are not used in the current implementation", "read", "Data1D_collection.cpp");
  
  ifstream fin(input_file.c_str()); checkIO(fin, input_file);
  string line;

  if (skip_nlines>0)
    for (int i=0; i<skip_nlines; ++i)
      getline(fin, line);

  vector<vector<double>> table, tr_table;
  while (getline(fin, line)) {
    stringstream ss(line); double NUM=par::defaultDouble;
    vector<double> vv;
    while(NUM>par::defaultDouble)
      vv.push_back(NUM);
    table.push_back(vv);
  }

  fin.clear(); fin.close(); 

  tr_table = table;

  for (size_t i=0; i<table.size(); i++)
    for (size_t j=0; j<table[i].size(); j++)
      tr_table[i][j] = table[j][i];

  m_ndataset = (tr_table.size()-1)/2;
  m_xsize.resize(m_ndataset,0);

  for (int i=0; i<m_ndataset; i++) {
    m_x.push_back(tr_table[0]);
    m_xsize[i]=m_x[i].size();

    for (int j=0; j<m_xsize[i]; j++) {
      m_data.push_back(tr_table[2*i+1][j]);
      m_error.push_back(tr_table[2*i+2][j]);
    }
  }

  m_ndata = m_data.size();

  m_covariance.resize(m_ndata, vector<double>(m_ndata,0));
  for (int i=0; i<m_ndata; i++)
    m_covariance[i][i] = pow(m_error[i],2);

  int index = 0;
  for (int i=0; i<m_ndataset; i++) {
    vector<int> vv;
    for (int j=0; j<m_xsize[i]; j++) {
      vv.push_back(index);
      index++;
    }
    m_index.push_back(vv);
  }
}


// ======================================================================================


void cbl::data::Data1D_collection::read (const std::vector<std::string> input_file, const int skip_nlines, const std::vector<int> column, const std::vector<int> column_data, const std::vector<int> column_errors, const std::vector<int> column_edges)
{
  (void)column; (void)column_data; (void)column_errors; (void)column_edges;
  WarningMsgCBL("column, column_data, column_errors, column_edges are not used in the current implementation", "read", "Data1D_collection.cpp");
  
  m_ndataset = input_file.size();
  m_xsize.resize(m_ndataset,0);

  for (int i=0; i<m_ndataset; i++) {
    ifstream fin(input_file[i].c_str()); checkIO(fin, input_file[i]);
    string line;

    if (skip_nlines>0)
      for (int j=0; j<skip_nlines; ++j)
	getline(fin, line);

    vector<double> xx, data, error;

    while (getline(fin, line)) {
      stringstream ss(line); double NUM;
      ss>>NUM; xx.push_back(NUM);
      ss>>NUM; data.push_back(NUM);
      ss>>NUM; error.push_back(NUM);
    }

    m_x.push_back(xx);
    for (size_t j=0; j<data.size(); j++) {
      m_data.push_back(data[j]);
      m_error.push_back(error[j]);
    }
    fin.clear(); fin.close(); 

    m_xsize[i]=m_x.size();
  }

  m_ndata = m_data.size();

  m_covariance.resize(m_ndata, vector<double>(m_ndata,0));
  for (int i=0; i<m_ndata; i++)
    m_covariance[i][i] = pow(m_error[i],2);

  int index = 0;
  for (int i=0; i<m_ndataset; i++) {
    vector<int> vv;
    for (int j=0; j<m_xsize[i]; j++) {
      vv.push_back(index);
      index++;
    }
    m_index.push_back(vv);
  }
}


// ======================================================================================


void cbl::data::Data1D_collection::Print (const int precision) const
{
  const int ndata = m_xsize[0];
  
  for (int i=1; i<m_ndataset; i++)
    checkEqual(m_x[0], m_x[i]);

  for (int i=0; i<ndata; i++) {
    coutCBL << setprecision(precision) << setw(8) << right << m_x[0][i] << "   ";
    for (int j=0; j<m_ndataset; j++)
      coutCBL << setprecision(precision) << setw(8) << right << m_data[m_index[j][i]]
	      << "   " << setprecision(precision) << setw(8) << right << m_error[m_index[j][i]] << "   ";
    cout << endl;
  }
}


// ======================================================================================


void cbl::data::Data1D_collection::write (const std::string dir, const std::string file, const std::string header, const int prec, const int ww, const int rank) const
{
  (void)rank;

  const int ndata = m_xsize[0];

  for (int i=1; i<m_ndataset; i++)
    checkEqual(m_x[0], m_x[i]);

  string file_out = dir+file;

  ofstream fout(file_out.c_str()); checkIO(fout, file_out);
  if (header!=par::defaultString)
    fout << "### "<< header <<" ###" << endl;

  const int bp = cout.precision();
  
  for (int i=0; i<ndata; i++) {
    cbl::Print(m_x[0][i], prec, ww, "", "  ", false, fout);
    for (int j=0; j<m_ndataset; j++) {
      cbl::Print(m_data[m_index[j][i]], prec, ww, "", "  ", false, fout);
      cbl::Print(m_error[m_index[j][i]], prec, ww, "", "\n", false, fout);
    }
  }

  cout.precision(bp);
  fout.close(); cout << endl; coutCBL << "I wrote the file: " << file_out << endl;
}


// ======================================================================================


void cbl::data::Data1D_collection::write (const std::string dir, const std::vector<std::string> files, const std::string header, const int prec, const int ww, const int rank) const
{
  (void)rank;

  checkDim(files, m_ndataset, "files");

  for (int i=0; i<m_ndataset; i++) {
    string file_out =dir+files[i];
    ofstream fout(file_out.c_str()); checkIO(fout, file_out);

    if (header!=par::defaultString)
      fout << "### "<< header <<" ###" << endl;

    const int bp = cout.precision();

    for (int j=0; j<m_xsize[i]; j++){
      cbl::Print(m_x[i][j], prec, ww, "", "  ", false, fout);
      cbl::Print(m_data[m_index[i][j]], prec, ww, "", "  ", false, fout);
      cbl::Print(m_error[m_index[i][j]], prec, ww, "", "\n", false, fout);
    }

    cout.precision(bp);
    fout.close(); cout << endl; coutCBL << "I wrote the file: " << file_out << endl;
  }

}


// ======================================================================================


shared_ptr<Data> cbl::data::Data1D_collection::cut (const int dataset, const double xmin, const double xmax) const
{
  vector<double> vxmin(m_ndataset, -1), vxmax(m_ndataset, -1);

  vxmin[dataset] = xmin; vxmax[dataset] = xmax;

  vector<bool> mask(m_ndata, true);
  vector<vector<double>> xx;

  for (int i=0; i<m_ndataset; i++) {
    vector<double> vv;
    for (int j=0; j<m_xsize[i]; j++) {
      if ( (m_x[i][j] < vxmin[i]) || (m_x[i][j]>vxmax[i]))
	mask[m_index[i][j]] = false;
      else
	vv.push_back(m_x[i][j]);
    }
    if (vv.size()>0)
      xx.push_back(vv);
  }

  checkDim(xx, 1, "x");

  vector<double> data, error;
  vector<vector<double>> covariance;
  Data::cut(mask, data, error, covariance);

  checkDim(xx[0], data.size(), "x");

  shared_ptr<Data> dd = make_shared<Data1D>(Data1D(xx[0], data, covariance));

  return dd;
}


// ======================================================================================


shared_ptr<Data> cbl::data::Data1D_collection::cut (const std::vector<double> xmin, const std::vector<double> xmax) const
{
  checkDim(xmin, m_ndataset, "xmin");
  checkDim(xmax, m_ndataset, "xmax");

  vector<bool> mask(m_ndata, true);
  vector<vector<double>> xx;

  for (int i=0; i<m_ndataset; i++) {
    vector<double> vv;
    for (int j=0; j<m_xsize[i]; j++) {
      if ( (m_x[i][j] < xmin[i]) || (m_x[i][j]>xmax[i]))
	mask[m_index[i][j]] = false;
      else
	vv.push_back(m_x[i][j]);
    }
    if (vv.size()>0)
      xx.push_back(vv);
  }

  vector<double> data, error;
  vector<vector<double>> covariance;
  Data::cut(mask, data, error, covariance);

  vector<vector<double>> data2;

  int index = 0;
  for (size_t i=0; i<xx.size(); i++) {
    vector<double> vv(xx[i].size(),0);
    for (size_t j=0; j<xx[i].size(); j++) {
      vv[j] = data[index]; 
      index++;
    }
    data2.push_back(vv);
  }

  shared_ptr<Data> dd = make_shared<Data1D_collection>(Data1D_collection(xx, data2, covariance));

  return dd;
}


// ======================================================================================


void cbl::data::Data1D_collection::write_covariance (const string dir, const string file, const int precision) const 
{
  checkDim(m_covariance, m_ndata, m_ndata, "covariance", false);
  
  string file_out = dir+file;
  ofstream fout(file_out.c_str()); checkIO(fout, file_out);

  fout << "### # [1] dataset1 # [2] dataset2 # [3] r1 # [4] r2 # [5] covariance # [6] correlation # [7] index1 # [8] index2 ### " << endl;

  int index1 = 0, index2= 0;

  for (int i=0; i<m_ndataset; i++)
    for (int j=0; j<m_ndataset; j++)
      for (int m=0; m<m_xsize[i]; m++)
	for (int n=0; n<m_xsize[j]; n++) {
	  index1 = m_index[i][m];
	  index2 = m_index[j][n];
	  fout << setiosflags(ios::fixed) << setprecision(precision) << setw(15) << right << i
	       << "  " << setiosflags(ios::fixed) << setprecision(precision) << setw(15) << right << j
	       << "  " << setiosflags(ios::fixed) << setprecision(precision) << setw(15) << right << m_x[i][m]
	       << "  " << setiosflags(ios::fixed) << setprecision(precision) << setw(15) << right << m_x[j][n]
	       << "  " << setiosflags(ios::fixed) << setprecision(precision) << setw(15) << right << m_covariance[index1][index2]
	       << "  " << setiosflags(ios::fixed) << setprecision(precision) << setw(15) << right << m_covariance[index1][index2]/sqrt(m_covariance[index1][index1]*m_covariance[index2][index2])
	       << "  " << setiosflags(ios::fixed) << setprecision(precision) << setw(5) << right << index1
	       << "  " << setiosflags(ios::fixed) << setprecision(precision) << setw(5) << right << index2 <<  endl;
	}

  fout.close(); cout << endl; coutCBL << "I wrote the file: " << file_out << endl;
}
