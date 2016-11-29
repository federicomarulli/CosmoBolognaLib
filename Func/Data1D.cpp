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
 *  @file Func/Data1D.cpp
 *
 *  @brief Methods of the class Data1D
 *
 *  This file contains the implementation of the methods of the class
 *  Data1D
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#include "Data1D.h"

using namespace cosmobl;
using namespace data;


// ======================================================================================


cosmobl::data::Data1D::Data1D (const string input_file, const int skip_nlines, const double xmin, const double xmax, const DataType dataType) : Data(dataType)
{
  read(input_file, skip_nlines);
  
  set_limits((xmin>par::defaultDouble) ? xmin : Min(m_x)-0.001, (xmax<-par::defaultDouble) ? xmax : Max(m_x)+0.001);
}


// ======================================================================================


cosmobl::data::Data1D::Data1D (const vector<double> x, const vector<double> fx, const double xmin, const double xmax, const DataType dataType) : Data(dataType)
{
  m_x = x;
  m_fx = fx;
  
  set_limits((xmin>par::defaultDouble) ? xmin : Min(m_x)-0.001, (xmax<-par::defaultDouble) ? xmax : Max(m_x)+0.001);
}


// ======================================================================================


cosmobl::data::Data1D::Data1D (const vector<double> x, const vector<double> fx, const vector<double> error_fx, const double xmin, const double xmax, const DataType dataType) : Data(dataType)
{
  m_x = x;
  m_fx = fx;
  m_error_fx = error_fx;

  set_limits((xmin>par::defaultDouble) ? xmin : Min(m_x)-0.001, (xmax<-par::defaultDouble) ? xmax : Max(m_x)+0.001);
}


// ======================================================================================


cosmobl::data::Data1D::Data1D (const vector<double> x, const vector<double> fx, const vector<vector<double> > covariance, const double xmin, const double xmax, const DataType dataType) : Data(dataType)
{
  m_x = x;
  m_fx = fx;
  m_covariance = covariance;

  set_limits((xmin>par::defaultDouble) ? xmin : Min(m_x)-0.001, (xmax<-par::defaultDouble) ? xmax : Max(m_x)+0.001);
  
  for (size_t i=0; i<m_covariance.size(); i++)
    m_error_fx.push_back(sqrt(m_covariance[i][i]));
}


// ======================================================================================


vector<double> cosmobl::data::Data1D::xx () const
{
  if (isSet(m_x_down) && isSet(m_x_up)) {
    vector<double> xx;
    for (int i=m_x_down; i<m_x_up; i++)
      xx.push_back(m_x[i]);
    return xx;
  }
  else
    return m_x;
}  


// ======================================================================================


vector<double> cosmobl::data::Data1D::fx () const 
{
  if (isSet(m_x_down) && isSet(m_x_up)) {
    vector<double> fx;
    for (int i=m_x_down; i<m_x_up; i++)
      fx.push_back(m_fx[i]);
    return fx;
  }
  else
    return m_fx;
}  


// ======================================================================================


vector<double> cosmobl::data::Data1D::error_fx () const
{
  if (isSet(m_x_down) && isSet(m_x_up)) {
    vector<double> efx;
    for (int i=m_x_down; i<m_x_up; i++)
      efx.push_back(m_error_fx[i]);
    return efx;
  }
  else
    return m_error_fx;
}  


// ======================================================================================


vector<vector<double> > cosmobl::data::Data1D::covariance () const
{
  if (m_covariance.size() == 0)
    ErrorCBL("Error in covariance of Data1D, covariance matrix is not set");

  if (isSet(m_x_down) && isSet(m_x_up)) {
    vector<vector<double>> cm;
    for (int i=m_x_down; i<m_x_up; i++) {
      vector<double> vv;
      for (int j=m_x_down; j<m_x_up; j++)
	vv.push_back(m_covariance[i][j]);
      cm.push_back(vv);
    }
    return cm;
  }
  
  else
    return m_covariance;
}  


// ======================================================================================


vector<vector<double> > cosmobl::data::Data1D::inverse_covariance () const
{
  if (m_inverse_covariance.size() == 0)
    ErrorCBL("Error in inverse_covariance of Data1D, inverted covariance matrix is not set. Run invert_covariance() first");

  if (isSet(m_x_down) && isSet(m_x_up)) {
    vector<vector<double>> icm;
    for (int i=m_x_down; i<m_x_up; i++) {
      vector<double> vv;
      for (int j=m_x_down; j<m_x_up; j++)
	vv.push_back(m_inverse_covariance[i][j]);
      icm.push_back(vv);
    }
    return icm;
  }
  else
    return m_inverse_covariance;
}  


// ======================================================================================


void cosmobl::data::Data1D::invert_covariance () 
{
  vector<vector<double> > cov = covariance(), icov;
  invert_matrix(cov, icov);

  vector<vector<double> > inverse_covariance(ndata(),vector<double>(ndata(),0));

  for (int i=m_x_down; i<m_x_up; i++)
    for (int j=m_x_down; j<m_x_up; j++)
      inverse_covariance[i][j] = icov[i-m_x_down][j-m_x_down];

  m_inverse_covariance = inverse_covariance;
}


// ======================================================================================


void cosmobl::data::Data1D::set_covariance (const string filename)
{  
  m_covariance.erase(m_covariance.begin(), m_covariance.end());
  m_error_fx.erase(m_error_fx.begin(), m_error_fx.end());

  ifstream fin(filename.c_str()); checkIO(fin, filename);

  vector<double> vv;
  m_covariance.push_back(vv);
  string line; int i = 0;

  while (getline(fin, line)) {
    stringstream ss(line);
    vector<double> num; double NN = -1.e30;
    while (ss>>NN) num.push_back(NN);
    
    if (num.size()>=3 && num[2]>-1.e29) 
      m_covariance[i].push_back(num[2]);
    else { i++; m_covariance.push_back(vv); }
  }

  m_covariance.erase(m_covariance.end()-1, m_covariance.end());
  fin.clear(); fin.close();

  for (size_t i=0; i<m_covariance.size(); i++)
    m_error_fx.push_back(sqrt(m_covariance[i][i]));
}


// ======================================================================================


void cosmobl::data::Data1D::set_covariance (const vector<vector<double> > covariance)
{
  m_error_fx.erase(m_error_fx.begin(), m_error_fx.end());
  m_covariance = covariance;

  for (size_t i=0; i<m_covariance.size(); i++)
    m_error_fx.push_back(sqrt(m_covariance[i][i]));
}


// ======================================================================================


void cosmobl::data::Data1D::set_limits (const double xmin, const double xmax)
{
  find_index(m_x, xmin, xmax, m_x_down, m_x_up);
}


// ======================================================================================


void cosmobl::data::Data1D::write_covariance (const string dir, const string file, const string xname, const string fxname) const 
{
  string file_out = dir+file;
  ofstream fout(file_out.c_str()); checkIO(fout, file_out);

  fout << "### [1] "<< xname << " # [2] " << xname << " # [3] covariance # [4] " << fxname << " # [5] index1 # [6] index2 ### " << endl;

  int cntr1 = 0, cntr2 = 0;
  for (size_t i=0; i<m_x.size(); ++i) {
    cntr2 = 0;
    for (size_t j=0; j<m_x.size(); ++j) { 
      fout << setiosflags(ios::fixed) << setprecision(4) << setw(8) << m_x[i] << "  " << setw(8) << m_x[j] << "  " << setw(8) << m_covariance[i][j] << " " << m_covariance[i][j]/sqrt(m_covariance[i][i]*m_covariance[j][j]) << " " << cntr1 << " " << cntr2 <<  endl;
      cntr2 ++;
    }
    cntr1 ++;
  }
   
  fout.close(); cout << endl; coutCBL << "I wrote the file: " << file_out << endl;
}


// ======================================================================================


void cosmobl::data::Data1D::read (const string input_file, const int skip_nlines)
{
  ifstream fin(input_file.c_str()); checkIO(fin, input_file);
  string line;

  if (skip_nlines>0)
    for (int i=0; i<skip_nlines; ++i)
      getline(fin, line);

  while (getline(fin, line)) {
    stringstream ss(line); double NUM;
    ss>>NUM; m_x.push_back(NUM);
    ss>>NUM; m_fx.push_back(NUM);
    ss>>NUM; m_error_fx.push_back(NUM);
  }

  fin.clear(); fin.close(); 
}


// ======================================================================================


void cosmobl::data::Data1D::write (const string dir, const string file, const string header, const int rank) const 
{
  (void)rank;
  
  string file_out = dir+file;
  ofstream fout(file_out.c_str()); checkIO(fout, file_out);

  fout << "### "<< header <<" ###" << endl;

  for (size_t i=0; i<m_x.size(); i++) 
    fout << setiosflags(ios::fixed) << setprecision(4) << setw(8) << m_x[i] << "  " << setw(8) << m_fx[i] << "  " << setw(8) << m_error_fx[i] << endl;
   
  fout.close(); cout << endl; coutCBL << "I wrote the file: " << file_out << endl;
}
