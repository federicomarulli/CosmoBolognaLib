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


// ======================================================================================


cosmobl::Data1D::Data1D (const string input_file, const double xmin, const double xmax)
{
  read(input_file);
  set_limits(xmin, xmax);
}


// ======================================================================================


cosmobl::Data1D::Data1D (const vector<double> x, const vector<double> fx, const double xmin, const double xmax)
{
  m_x = x;
  m_fx = fx;
  double xMin = (xmin>-par::defaultDouble) ? xmin : Min(m_x)-0.001; 
  double xMax = (xmin<par::defaultDouble) ? xmin : Max(m_x)+0.001; 

  find_index(m_x, xMin, xMax, m_x_down, m_x_up);
}


// ======================================================================================


cosmobl::Data1D::Data1D (const vector<double> x, const vector<double> fx, const vector<double> error_fx, const double xmin, const double xmax)
{
  m_x = x;
  m_fx = fx;
  m_error_fx = error_fx;
  double xMin = (xmin>-par::defaultDouble) ? xmin : Min(m_x)-0.001; 
  double xMax = (xmin<par::defaultDouble) ? xmin : Max(m_x)+0.001; 

  find_index(m_x, xMin, xMax, m_x_down, m_x_up);
}


// ======================================================================================


cosmobl::Data1D::Data1D (const vector<double> x, const vector<double> fx, const vector<vector<double> > covariance_fx, const double xmin, const double xmax)
{
  m_x = x;
  m_fx = fx;
  m_covariance_fx = covariance_fx;
  double xMin = (xmin>-par::defaultDouble) ? xmin : Min(m_x)-0.001; 
  double xMax = (xmin<par::defaultDouble) ? xmin : Max(m_x)+0.001; 

  for(size_t i=0;i<m_covariance_fx.size();i++)
    m_error_fx.push_back(sqrt(m_covariance_fx[i][i]));

  find_index(m_x, xMin, xMax, m_x_down, m_x_up);
}


// ======================================================================================


vector<double> cosmobl::Data1D::xx () const
{
  if(isSet(m_x_down) && isSet(m_x_up)){
    vector<double> xx;
    for(int i=m_x_down; i<m_x_up;i++)
      xx.push_back(m_x[i]);
    return xx;
  }
  else
    return m_x;
}  


// ======================================================================================


vector<double> cosmobl::Data1D::fx () const 
{
  if(isSet(m_x_down) && isSet(m_x_up)){
    vector<double> fx;
    for(int i=m_x_down; i<m_x_up;i++)
      fx.push_back(m_fx[i]);
    return fx;
  }
  else
    return m_fx;
}  


// ======================================================================================


vector<double> cosmobl::Data1D::error_fx () const
{
  if(isSet(m_x_down) && isSet(m_x_up)){
    vector<double> efx;
    for(int i=m_x_down; i<m_x_up;i++)
      efx.push_back(m_error_fx[i]);
    return efx;
  }
  else
    return m_error_fx;
}  


// ======================================================================================


vector<vector<double> > cosmobl::Data1D::covariance_fx () const
{
  if(isSet(m_x_down) && isSet(m_x_up)){
    vector<vector<double>> cm;
    for(int i=m_x_down; i<m_x_up;i++){
      vector<double> vv;
      for(int j=m_x_down; j<m_x_up;j++)
	vv.push_back(m_covariance_fx[i][j]);
      cm.push_back(vv);
    }
    return cm;
  }
  else
    return m_covariance_fx;
}  


// ======================================================================================


vector<vector<double> > cosmobl::Data1D::inverse_covariance_fx () const
{
  if(isSet(m_x_down) && isSet(m_x_up)){
    vector<vector<double>> icm;
    for(int i=m_x_down; i<m_x_up;i++){
      vector<double> vv;
      for(int j=m_x_down; j<m_x_up;j++)
	vv.push_back(m_inverse_covariance_fx[i][j]);
      icm.push_back(vv);
    }
    return icm;
  }
  else
    return m_inverse_covariance_fx;
}  



// ======================================================================================


void cosmobl::Data1D::set_covariance_fx (const string filename)
{  
  m_covariance_fx.erase(m_covariance_fx.begin(), m_covariance_fx.end());
  m_error_fx.erase(m_error_fx.begin(), m_error_fx.end());

  ifstream fin (filename.c_str());
  if (!fin) {
    string Warn = "Attention: the file " + filename + " does not exist!";
    WarningMsg (Warn);
  }

  vector<double> vv ;
  m_covariance_fx.push_back(vv);
  string line; int i = 0;

  while (getline(fin, line)) {
    stringstream ss(line);
    vector<double> num; double NN = -1.e30;
    while (ss>>NN) num.push_back(NN);
    if (num.size()>=3 && num[2]>-1.e29){
      m_covariance_fx[i].push_back(num[2]);
    }
    else {i++; m_covariance_fx.push_back(vv);}
  }

  m_covariance_fx.erase(m_covariance_fx.end()-1, m_covariance_fx.end());
  fin.clear(); fin.close();

  for(size_t i=0;i<m_covariance_fx.size();i++)
    m_error_fx.push_back(sqrt(m_covariance_fx[i][i]));

  vector<vector<double> > cov = covariance_fx(),icov;
  invert_matrix(cov,icov);

  m_inverse_covariance_fx.resize(m_covariance_fx.size(),vector<double>(m_covariance_fx.size(),0));

  for(int i=x_down();i<x_up();i++)
    for(int j=x_down();j<x_up();j++)
      m_inverse_covariance_fx[i][j] = icov[i-x_down()][j-x_down()];

}


// ======================================================================================


void cosmobl::Data1D::set_covariance_fx (const vector<vector<double> > covariance_fx)
{
  m_error_fx.erase(m_error_fx.begin(), m_error_fx.end());
  m_covariance_fx=covariance_fx;

  for (size_t i=0; i<m_covariance_fx.size(); i++)
    m_error_fx.push_back(sqrt(m_covariance_fx[i][i]));

}


// ======================================================================================


void cosmobl::Data1D::read (const string input_file)
{
  ifstream fin(input_file.c_str());
  string line;
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


void cosmobl::Data1D::write (const string dir, const string file, const string xname, const string fxname, const int rank) const 
{
  string file_out = dir+file;
  ofstream fout (file_out.c_str()); checkIO(file_out, 0);

  fout << "### "<< xname << " " << fxname << " error ###" << endl;

  for (size_t i=0; i<m_x.size(); i++) 
      fout << setiosflags(ios::fixed) << setprecision(4) << setw(8) << m_x[i] << "  " << setw(8) << m_fx[i] << "  " << setw(8) << m_error_fx[i] << endl;
   
  fout.close(); cout << endl << "I wrote the file: " << file_out << endl;
}


// ======================================================================================


void cosmobl::Data1D::set_limits (const double xmin, const double xmax)
{
  find_index(m_x, xmin, xmax, m_x_down, m_x_up);
}
