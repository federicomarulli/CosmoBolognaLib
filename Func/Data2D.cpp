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
 *  @file Func/Data2D.cpp
 *
 *  @brief Methods of the class Data2D
 *
 *  This file contains the implementation of the methods of the class
 *  Data2D
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#include "Data2D.h"

using namespace cosmobl;

using namespace data;


// ======================================================================================


cosmobl::data::Data2D::Data2D (const vector<double> x, const vector<double> y, const vector<vector<double> > data) : Data(DataType::_2D_data_)
{
  m_x = x;
  m_y = y;

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


cosmobl::data::Data2D::Data2D (const vector<double> x, const vector<double> y, const vector<vector<double> > data, const vector<vector<double> > error) : Data(DataType::_2D_data_)
{
  m_x = x;
  m_y = y;

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

  invert_covariance();
}


// ======================================================================================


cosmobl::data::Data2D::Data2D (const vector<double> x, const vector<double> y, const vector<double> data, const vector<double> error) : Data(DataType::_2D_data_)
{
  m_x = x;
  m_y = y;

  m_xsize = m_x.size();
  m_ysize = m_y.size();

  checkDim(data, m_xsize*m_ysize, "data");
  checkDim(error, m_xsize*m_ysize, "error");

  m_data = data;
  m_error = error;

  m_ndata = m_data.size();

  m_covariance.resize(m_ndata, vector<double>(m_ndata,0));
  for (int i=0; i<m_ndata; i++)
    m_covariance[i][i] = pow(m_error[i],2);
  invert_covariance();
}


// ======================================================================================


cosmobl::data::Data2D::Data2D (const vector<double> x, const vector<double> y, const vector<double> data, const vector<vector<double> > covariance_matrix) : Data(DataType::_2D_data_)
{
  m_x = x;
  m_y = y;

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
  invert_covariance();

  for (int i=0;i<m_ndata;i++)
    m_error[i] = sqrt(m_covariance[i][i]);
}


// ======================================================================================


void cosmobl::data::Data2D::data(vector<vector<double>> &data) const
{
  data.erase(data.begin(), data.end());
  data.resize(m_xsize, vector<double>(m_ysize,0));

  for (int i=0; i<m_xsize; i++)
    for (int j=0; j< m_ysize; j++)
      data[i][j] = this->data(i,j);
}


// ======================================================================================


void cosmobl::data::Data2D::error(vector<vector<double>> &error) const
{
  error.erase(error.begin(), error.end());
  error.resize(m_xsize, vector<double>(m_ysize,0));

  for (int i=0; i<m_xsize; i++)
    for (int j=0; j< m_ysize; j++) {
      error[i][j] = this->error(i,j);
    }
}


// ======================================================================================


void cosmobl::data::Data2D::read (const string input_file, const int skip_nlines)
{
  (void)skip_nlines;
  
  ErrorCBL("Error in cosmobl::data::Data2D::read : work in progress!", glob::ExitCode::_workInProgress_);

  ifstream fin(input_file.c_str()); checkIO(fin, input_file);
  
  fin.clear(); fin.close();
}


// ======================================================================================


void cosmobl::data::Data2D::write (const string dir, const string file, const string header, const bool full, const int precision, const int rank) const 
{
  (void)rank;
  
  string file_out = dir+file;
  ofstream fout(file_out.c_str()); checkIO(fout, file_out);

  if (header!=par::defaultString)
    fout << "### " << header << " ###" << endl;

  for (int i=0; i<m_xsize; ++i)
    for (int j=0; j<m_ysize; ++j) {
      int index = j+m_ysize*i;
      fout << setiosflags(ios::fixed) << setprecision(precision) << setw(8) << m_x[i] << "  " << setw(8) << m_y[j] << "  " << setw(8) << m_data[index] << "  " << setw(8) << m_error[index] << endl;
    }

 
  if (full) { // duplicate the information in the other 3 quadrants

    for (int i=0; i<m_xsize; ++i)
      for (int j=0; j<m_ysize; ++j) {
	int index = j+m_ysize*i;
	fout << setiosflags(ios::fixed) << setprecision(precision) << setw(8) << m_x[i] << "  " << setw(8) << -m_y[j] << "  " << setw(8) << m_data[index] << "  " << setw(8) << m_error[index]<< endl;
      }

    for (int i=0; i<m_xsize; ++i)
      for (int j=0; j<m_ysize; ++j) {
	int index = j+m_ysize*i;
	fout << setiosflags(ios::fixed) << setprecision(precision) << setw(8) << -m_x[i] << "  " << setw(8) << -m_y[j] << "  " << setw(8) << m_data[index] << "  " << setw(8) << m_error[index] << endl;
      }

    for (int i=0; i<m_xsize; ++i)
      for (int j=0; j<m_ysize; ++j) {
	int index = j+m_ysize*i;
	fout << setiosflags(ios::fixed) << setprecision(precision) << setw(8) << -m_x[i] << "  " << setw(8) << m_y[j] << "  " << setw(8) << m_data[index] << "  " << setw(8) << m_error[index]<< endl;
      }
    
  }

  fout.close(); cout << endl; coutCBL << "I wrote the file: " << file_out << endl << endl;
}


// ======================================================================================


void cosmobl::data::Data2D::write_covariance (const string dir, const string file, const int precision) const
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
	  fout << setiosflags(ios::fixed) << setprecision(precision) << setw(8) << m_x[i] << "  " << setw(8) << m_y[j] << "  " << setw(8) << m_x[k] << "  " <<  setw(8) << m_y[l] << "  " <<  setw(8) <<  m_covariance[index1][index2] << " " << m_covariance[index1][index2]/sqrt(m_covariance[index1][index1]*m_covariance[index2][index2]) << " " << index1 << " " << index2 <<  endl;
	}
      }
    }
  }

  fout.close(); cout << endl; coutCBL << "I wrote the file: " << file_out << endl;

}


// ======================================================================================


shared_ptr<Data> cosmobl::data::Data2D::cut(const double xmin, const double xmax, const double ymin, const double ymax) const
{
  vector<bool> mask(m_ndata, true);
  vector<double> xx, yy;

  for (int i=0; i<m_xsize; i++) {
    for (int j=0; j<m_ysize; j++) {
      if ( ((m_x[i] < xmin) || (m_x[i]>xmax)) || ((m_y[j] < ymin) || (m_y[j]>ymax)))
	mask[j+m_ysize*i] = false;
    }
  }

  for (int i=0; i<m_xsize; i++)
    if ( !((m_x[i] < xmin) || (m_x[i]>xmax)))
      xx.push_back(m_x[i]);

  for (int i=0; i<m_ysize; i++)
    if (!((m_y[i] < ymin) || (m_y[i]>ymax)))
      yy.push_back(m_y[i]);

  vector<double> data, error;
  vector<vector<double>> covariance;
  Data::cut(mask, data, error, covariance);

  shared_ptr<Data> dd = make_shared<Data2D>(Data2D(xx, yy, data, covariance));

  return dd;
}


// ======================================================================================


shared_ptr<Data> cosmobl::data::Data2D::as_factory()
{
  shared_ptr<Data> dd = make_shared<Data2D>(Data2D(m_x, m_data, m_covariance));

  return dd;
}
