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
 *  @file Data/Data2D_extra.cpp
 *
 *  @brief Methods of the class Data2D_extra
 *
 *  This file contains the implementation of the methods of the class
 *  Data2D
 *
 *  @authors Federico Marulli
 *
 *  @authors federico.marulli3@unibo.it
 */

#include "Data2D_extra.h"

using namespace std;

using namespace cbl;
using namespace data;
using namespace glob;


// ======================================================================================


void cbl::data::Data2D_extra::read (const string input_file, const int skip_nlines, const vector<int> column, const vector<int> column_data, const vector<int> column_errors, const std::vector<int> column_edges)
{ 
  if (column_data.size()!=column_errors.size()) ErrorCBL("the sizes of column_data and columt_errors must be equal!", "read", "Data2D_extra.cpp");
  
  // default column of input data values
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
      
      // if column[0] is larger than 0, the column of x coordinates is
      // the one specified in input
      const int ind_x = max(0, column[0]-1);
      checkDim(num, ind_x+1, "num", false);
      m_x.emplace_back(num[ind_x]);

      // if column[1] is larger than 0, the column of y coordinates is
      // the one specified in input
      const int ind_y = max(0, column[1]-1);
      checkDim(num, ind_y+1, "num", false);
      m_y.emplace_back(num[ind_y]);

      // if the size of the column_data vector is 0, the columns of
      // data values are the ones specified in input
      const int ind_data = (column_data.size()==0) ? column_data_default-1 : column_data[cl]-1;
      checkDim(num, ind_data+1, "num", false);
      m_data.emplace_back(num[ind_data]);

      // if the size of the column_errors vector is 0, the columns of
      // error values are the ones specified in input; if the error
      // column is not present, the errors will be set to 1, by
      // default
      const size_t ind_error = (column_errors.size()==0) ? column_error_default-1 : column_errors[cl]-1;
      if (num.size()<ind_error+1 && column_errors.size()>1)
	WarningMsgCBL("the errors cannot be retrieved from the provided input file, and will be set to 1", "read", "Data2D_extra.cpp");
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
      
      // get the extra info
      int ind = 0;
      for (size_t ex=num.size()-m_extra_info.size(); ex<num.size(); ++ex) 
	m_extra_info[ind++].emplace_back(num[ex]);
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


void cbl::data::Data2D_extra::Print (const int precision) const 
{
  for (int i=0; i<m_xsize; ++i)
    for (int j=0; j<m_ysize; ++j) {
      int index = j+m_ysize*i;
      coutCBL << setprecision(precision) << setw(8) << right << m_x[i]
	      << "   " << setprecision(precision) << setw(8) << right << m_y[j]
	      << "   " << setprecision(precision) << setw(8) << right << m_data[index]
	      << "   " << setprecision(precision) << setw(8) << right << m_error[index];
      
      for (size_t ex=0; ex<m_extra_info.size(); ++ex)
        coutCBL << "   " << setprecision(precision) << setw(8) << right << m_extra_info[ex][index];
      cout << endl;
    }
}


// ======================================================================================


void cbl::data::Data2D_extra::write (const string dir, const string file, const string header, const bool full, const int prec, const int ww, const int rank) const 
{
  (void)rank;

  string file_out = dir+file;
  ofstream fout(file_out.c_str()); checkIO(fout, file_out);

  if (header!=par::defaultString)
    fout << "### " << header << " ###" << endl;

  const int bp = std::cout.precision();

  for (int i=0; i<m_xsize; ++i)
    for (int j=0; j<m_ysize; ++j) {
      int index = j+m_ysize*i;
      cbl::Print(m_x[i], prec, ww, "", "  ", false, fout);
      cbl::Print(m_y[j], prec, ww, "", "  ", false, fout);
      cbl::Print(m_data[index], prec, ww, "", "  ", false, fout);
      cbl::Print(m_error[index], prec, ww, "", "  ", false, fout);
      for (size_t ex=0; ex<m_extra_info.size(); ++ex) 
	cbl::Print(m_extra_info[ex][index], prec, ww, "", "  ", false, fout);
      fout << endl;
    }


  if (full) { // duplicate the information in the other 3 quadrants

    for (int i=0; i<m_xsize; ++i)
      for (int j=0; j<m_ysize; ++j) {
        int index = j+m_ysize*i;
	cbl::Print(m_x[i], prec, ww, "", "  ", false, fout);
	cbl::Print(-m_y[j], prec, ww, "", "  ", false, fout);
	cbl::Print(m_data[index], prec, ww, "", "  ", false, fout);
	cbl::Print(m_error[index], prec, ww, "", "  ", false, fout);
	for (size_t ex=0; ex<m_extra_info.size(); ++ex) 
	  cbl::Print(m_extra_info[ex][index], prec, ww, "", "  ", false, fout);
	fout << endl;
      }
    
    for (int i=0; i<m_xsize; ++i)
      for (int j=0; j<m_ysize; ++j) {
        int index = j+m_ysize*i;
	cbl::Print(-m_x[i], prec, ww, "", "  ", false, fout);
	cbl::Print(-m_y[j], prec, ww, "", "  ", false, fout);
	cbl::Print(m_data[index], prec, ww, "", "  ", false, fout);
	cbl::Print(m_error[index], prec, ww, "", "  ", false, fout);
	for (size_t ex=0; ex<m_extra_info.size(); ++ex) 
	  cbl::Print(m_extra_info[ex][index], prec, ww, "", "  ", false, fout);
	fout << endl;
      }

    for (int i=0; i<m_xsize; ++i)
      for (int j=0; j<m_ysize; ++j) {
        int index = j+m_ysize*i;
	cbl::Print(-m_x[i], prec, ww, "", "  ", false, fout);
	cbl::Print(m_y[j], prec, ww, "", "  ", false, fout);
	cbl::Print(m_data[index], prec, ww, "", "  ", false, fout);
	cbl::Print(m_error[index], prec, ww, "", "  ", false, fout);
	for (size_t ex=0; ex<m_extra_info.size(); ++ex) 
	  cbl::Print(m_extra_info[ex][index], prec, ww, "", "  ", false, fout);
	fout << endl;
      } 
  }

  std::cout.precision(bp);
  fout.close(); cout << endl; coutCBL << "I wrote the file: " << file_out << endl << endl;
}


// ======================================================================================


shared_ptr<Data> cbl::data::Data2D_extra::cut (const double xmin, const double xmax, const double ymin, const double ymax) const
{
  vector<bool> mask(m_ndata, true);
  vector<double> xx, yy, edges_xx, edges_yy;
  
  for (int i=0; i<m_xsize; i++) 
    for (int j=0; j<m_ysize; j++) 
      if (((m_x[i]<xmin) || (m_x[i]>xmax)) || ((m_y[j]<ymin) || (m_y[j]>ymax)))
	mask[j+m_ysize*i] = false;
  
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
  vector<vector<double>> extra_info;

  int ndata_eff = 0;

  for (int i=0; i<m_ndata; i++)
    if (mask[i])
      ndata_eff ++;

  if (ndata_eff<1)
    ErrorCBL("no elements left!", "cut", "Data2D_extra.cpp");
  
  data.resize(ndata_eff, 0);
  error.resize(ndata_eff, 0);
  covariance.resize(ndata_eff, vector<double>(ndata_eff, 0));
  extra_info.resize(m_extra_info.size(), vector<double>(ndata_eff, 0));
  
  int index1 = 0;
  for (int i=0; i<m_ndata; i++) {
    if (mask[i]) {
      data[index1] = m_data[i];
      error[index1] = m_error[i];

      for (size_t j=0; j<m_extra_info.size(); j++) 
	extra_info[j][index1] = m_extra_info[j][i];
 
      int index2 = 0;
      for (int j=0; j<m_ndata; j++) {
	if (mask[j]) {
	  covariance[index1][index2] = m_covariance[i][j];
	  index2 ++;
	}
      }
      index1 ++;
    }
  }
  
  shared_ptr<Data> dd = make_shared<Data2D_extra>(Data2D_extra(xx, yy, data, covariance, extra_info, edges_xx, edges_yy));
  
  return dd;
}

