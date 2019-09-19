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
 *  @authors federico.marulli3@unbo.it
 */

#include "Data2D_extra.h"

using namespace std;

using namespace cbl;
using namespace data;
using namespace glob;


// ======================================================================================


void cbl::data::Data2D_extra::read (const string input_file, const int skip_nlines, const int column_x, const vector<int> column_data, const vector<int> column_errors)
{
  (void)input_file; (void)skip_nlines; (void)column_x; (void)column_data; (void)column_errors;

  ErrorCBL("", "read", "Data2D_extra.cpp", ExitCode::_workInProgress_);
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


void cbl::data::Data2D_extra::write (const string dir, const string file, const string header, const bool full, const int precision, const int rank) const 
{
  (void)rank;

  string file_out = dir+file;
  ofstream fout(file_out.c_str()); checkIO(fout, file_out);

  if (header!=par::defaultString)
    fout << "### " << header << " ###" << endl;

  for (int i=0; i<m_xsize; ++i)
    for (int j=0; j<m_ysize; ++j) {
      int index = j+m_ysize*i;
      fout << setprecision(precision) << setw(15) << right << m_x[i]
	   << "  " << setprecision(precision) << setw(15) << right << m_y[j]
	   << "  " << setprecision(precision) << setw(15) << right << m_data[index]
	   << "  " << setprecision(precision) << setw(15) << right << m_error[index];

      for (size_t ex=0; ex<m_extra_info.size(); ++ex)
        fout << "  " << setprecision(precision) << setw(15) << right << m_extra_info[ex][index];
      fout << endl;
    }


  if (full) { // duplicate the information in the other 3 quadrants

    for (int i=0; i<m_xsize; ++i)
      for (int j=0; j<m_ysize; ++j) {
        int index = j+m_ysize*i;
        fout << setprecision(precision) << setw(15) << right << m_x[i]
	     << "  " << setprecision(precision) << setw(15) << right << -m_y[j]
	     << "  " << setprecision(precision) << setw(15) << right << m_data[index]
	     << "  " << setprecision(precision) << setw(15) << right << m_error[index];

        for (size_t ex=0; ex<m_extra_info.size(); ++ex)
          fout  << "  " << setprecision(precision) << setw(15) << m_extra_info[ex][index];
        fout << endl;
      }

    for (int i=0; i<m_xsize; ++i)
      for (int j=0; j<m_ysize; ++j) {
        int index = j+m_ysize*i;
        fout << setprecision(precision) << setw(15) << right << -m_x[i]
	     << "  " << setprecision(precision) << setw(15) << right << -m_y[j]
	     << "  " << setprecision(precision) << setw(15) << right << m_data[index]
	     << "  " << setprecision(precision) << setw(15) << right << m_error[index];

        for (size_t ex=0; ex<m_extra_info.size(); ++ex)
          fout  << "  " << setprecision(precision) << setw(15) << right << m_extra_info[ex][index];
        fout << endl;
      }

    for (int i=0; i<m_xsize; ++i)
      for (int j=0; j<m_ysize; ++j) {
        int index = j+m_ysize*i;
        fout << setprecision(precision) << setw(15) << right << -m_x[i]
	     << "  " << setprecision(precision) << setw(15) << right << m_y[j]
	     << "  " << setprecision(precision) << setw(15) << right << m_data[index]
	     << "  " << setprecision(precision) << setw(15) << right << m_error[index];

        for (size_t ex=0; ex<m_extra_info.size(); ++ex)
          fout << "  " << setprecision(precision) << setw(15) << right << m_extra_info[ex][index];
        fout << endl;
      }

  }

  fout.close(); cout << endl; coutCBL << "I wrote the file: " << file_out << endl << endl;
}


// ======================================================================================


shared_ptr<Data> cbl::data::Data2D_extra::cut (const double xmin, const double xmax, const double ymin, const double ymax) const
{
  vector<bool> mask(m_ndata, true);
  vector<double> xx, yy;

  for (int i=0; i<m_xsize; i++) 
    for (int j=0; j<m_ysize; j++) 
      if (((m_x[i] < xmin) || (m_x[i]>xmax)) || ((m_y[j] < ymin) || (m_y[j]>ymax)))
	mask[j+m_ysize*i] = false;

  for (int i=0; i<m_xsize; i++)
    if (!((m_x[i] < xmin) || (m_x[i]>xmax)))
      xx.push_back(m_x[i]);

  for (int i=0; i<m_ysize; i++)
    if (!((m_y[i] < ymin) || (m_y[i]>ymax)))
      yy.push_back(m_y[i]);

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

  shared_ptr<Data> dd = make_shared<Data2D_extra>(Data2D_extra(xx, yy, data, covariance, extra_info));

  return dd;
}

