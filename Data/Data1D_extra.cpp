/********************************************************************
 *  Copyright (C) 2016 by Federico Marulli                          *
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
 *  @file Data/Data1D_extra.cpp
 *
 *  @brief Methods of the class Data1D_extra
 *
 *  This file contains the implementation of the methods of the class
 *  Data1D
 *
 *  @authors Federico Marulli
 *
 *  @authors federico.marulli3@unibo.it
 */

#include "Data1D_extra.h"

using namespace std;

using namespace cbl;
using namespace data;


// ======================================================================================


void cbl::data::Data1D_extra::read (const string input_file, const int skip_nlines, const int column_x, const vector<int> column_data, const vector<int> column_errors)
{
  if (column_data.size()!=column_errors.size()) ErrorCBL("the sizes of column_data and columt_errors must be equal!", "read", "Data1D.cpp");
  
  // default column of input data values
  int column_data_default = column_x+1;
  
  // default column of input error values
  int column_error_default = column_x+2;

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
    while (getline(fin, line)) {
      
      stringstream ss(line);
      vector<double> num; double NUM;
      while (ss>>NUM) num.emplace_back(NUM);
      
      // if column_x is larger than 0, the column of x coordinates is
      // the one specified in input
      const int ind_x = max(0, column_x-1);
      checkDim(num, ind_x+1, "num", false);
      m_x.emplace_back(num[ind_x]);

      // if the size of the column_data vector is larger than 1, the
      // columns of data values are the ones specified in input
      const int ind_data = (column_data.size()<1) ? column_data_default : column_data[cl]-1;
      checkDim(num, ind_data+1, "num", false);
      m_data.emplace_back(num[ind_data]);

      // if the size of the column_errors vector is larger than 1, the
      // columns of error values are the ones specified in input; if
      // the error column is not present, the errors will be set to 1,
      // by default
      const size_t ind_error = (column_errors.size()<1) ? column_error_default : column_errors[cl]-1;
      if (num.size()<ind_error+1 && column_errors.size()>1)
	WarningMsgCBL("the errors cannot be retrieved from the provided input file, and will be set to 1", "read", "Data1D.cpp");
      m_error.emplace_back((num.size()<ind_error+1) ? 1 : num[ind_error]);

      // get the extra info
      int ind = 0;
      for (size_t ex=num.size()-m_extra_info.size(); ex<num.size(); ++ex) 
	m_extra_info[ind++].emplace_back(num[ex]);
    }

    // go back at the beginning of the file
    fin.clear(); fin.seekg(ios::beg);

    column_data_default += 2;
    column_error_default += 2;
  }
  
  fin.clear(); fin.close();


  // set the data size and diagonal covariance

  m_ndata = m_data.size();

  m_covariance.resize(m_ndata, vector<double>(m_ndata, 0));
  for (int i=0; i<m_ndata; i++)
    m_covariance[i][i] = pow(m_error[i], 2);

}


// ======================================================================================


void cbl::data::Data1D_extra::Print (const int precision) const 
{
  for (size_t i=0; i<m_x.size(); ++i) {
    coutCBL << setprecision(precision) << setw(8) << right << m_x[i] 
	 << "   " << setprecision(precision) << setw(8) << right << m_data[i] 
	 << "   " << setprecision(precision) << setw(8) << right << m_error[i];
    
    for (size_t ex=0; ex<m_extra_info.size(); ++ex)
      coutCBL << "   " << setprecision(precision) << setw(8) << m_extra_info[ex][i];
    cout << endl;
  }
}


// ======================================================================================


void cbl::data::Data1D_extra::write (const string dir, const string file, const string header, const int prec, const int ww, const int rank) const 
{
  (void)rank;
  
  string file_out = dir+file;
  ofstream fout(file_out.c_str()); checkIO(fout, file_out);
  
  fout << "### "<< header <<" ###" << endl;

  for (size_t i=0; i<m_x.size(); ++i) {
    cbl::Print(m_x[i], prec, ww, "", "  ", false, fout);
    cbl::Print(m_data[i], prec, ww, "", "  ", false, fout);
    cbl::Print(m_error[i], prec, ww, "", "  ", false, fout);
    for (size_t ex=0; ex<m_extra_info.size(); ++ex)
      cbl::Print(m_extra_info[ex][i], prec, ww, "", "  ", false, fout);
    fout << endl;
  }
  
  fout.close(); cout << endl ; coutCBL << "I wrote the file: " << file_out << endl;
}


// ======================================================================================


std::shared_ptr<Data> cbl::data::Data1D_extra::cut (const std::vector<bool> mask) const
{
  checkDim(mask, m_ndata, "mask");

  vector<double> xx;
  vector<double> data, error;
  vector<vector<double>> covariance;
  vector<vector<double>> extra_info;

  int ndata_eff = 0;

  for (int i=0; i<m_ndata; i++)
    if (mask[i])
      ndata_eff +=1;

  if (ndata_eff <1)
    ErrorCBL("no elements left!", "cut", "Data1D_extra.cpp");

  xx.resize(ndata_eff, 0);
  data.resize(ndata_eff, 0);
  error.resize(ndata_eff, 0);
  covariance.resize(ndata_eff, vector<double>(ndata_eff, 0));
  extra_info.resize(m_extra_info.size(), vector<double>(ndata_eff, 0));

  int index1 = 0;
  for (int i=0; i<m_ndata; i++) {
    if (mask[i]) {
      xx[index1] = m_x[i];
      data[index1] = m_data[i];
      error[index1] = m_error[i];

      for (size_t j=0; j<m_extra_info.size(); j++) {
	extra_info[j][index1] = m_extra_info[j][i];
      }

      int index2 = 0;
      for (int j=0; j<m_ndata; j++) {
	if (mask[j]) {
	  covariance[index1][index2] = m_covariance[i][j];
	  index2++;
	}
      }
      index1++;
    }
  }

  shared_ptr<Data> dd = make_shared<Data1D_extra>(Data1D_extra(xx, data, covariance, extra_info));

  return dd;
}


// ======================================================================================


shared_ptr<Data> cbl::data::Data1D_extra::cut (const double xmin, const double xmax) const
{
  vector<bool> mask(m_ndata, true);
  vector<double> xx;
  for (int i=0; i<m_ndata; i++) {
    if ((m_x[i] < xmin) || (m_x[i]>xmax))
      mask[i] = false;
    else
      xx.push_back(m_x[i]);
  }

  return cut(mask);
}
