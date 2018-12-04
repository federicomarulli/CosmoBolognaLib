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
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
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


void cbl::data::Data1D::read (const string input_file, const int skip_nlines)
{
  ifstream fin(input_file.c_str()); checkIO(fin, input_file);
  string line;

  bool moreColumns = true;
  unsigned int col = 3, ind1 = 1, ind2 = 2;

  // enter here the first time, and continue to read in case of more the 3 columns (e.g. for clustering multipoles)
  while (moreColumns) {

    // skip the first nlines (in case of header lines)
    if (skip_nlines>0)
      for (int i=0; i<skip_nlines; ++i)
	getline(fin, line);

    // read the file lines
    while (getline(fin, line)) {
      
      stringstream ss(line);
      vector<double> num; double NUM;
      while (ss>>NUM) num.emplace_back(NUM);

      if ((int)num.size()<2) ErrorCBL("Error in cbl::data::Data1D::read(): the input file has less than 2 columns!");
      
      // if there are more than 3 columns then the next columns will
      // be read and added to the first one
      if (num.size()<=col) moreColumns = false;
      
      m_x.emplace_back(num[0]);
      m_data.emplace_back(num[ind1]);
      m_error.emplace_back((num.size()>ind2) ? num[ind2] : 1);
       
    }

    // go back at the beginning of the file
    fin.clear(); fin.seekg(ios::beg);
 
    col += 2;
    ind1 += 2;
    ind2 += 2;
  }
  
  fin.clear(); fin.close();


  // set the data type and diagonal covariance

  m_ndata = m_data.size();

  m_covariance.resize(m_ndata, vector<double>(m_ndata, 0));
  for (int i=0; i<m_ndata; i++)
    m_covariance[i][i] = pow(m_error[i], 2);
}


// ======================================================================================


void cbl::data::Data1D::write (const string dir, const string file, const string header, const int precision, const int rank) const 
{
  (void)rank;
  
  string file_out = dir+file;
  ofstream fout(file_out.c_str()); checkIO(fout, file_out);

  if (header!=par::defaultString)
    fout << "### " << header << " ###" << endl;

  for (size_t i=0; i<m_x.size(); i++)
    fout << setprecision(precision) << setw(15) << right << m_x[i] 
	 << "  " << setprecision(precision) << setw(15) << right << m_data[i] 
	 << "  " << setprecision(precision) << setw(15) << right << m_error[i] << endl;
   
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
      fout << setprecision(precision) << setw(15) << right << m_x[i]
	   << "  " << setprecision(precision) << setw(15) << right << m_x[j]
	   << "  " << setprecision(precision) << setw(15) << right << m_covariance[i][j]
	   << "  " << setprecision(precision) << setw(15) << right << m_covariance[i][j]/sqrt(m_covariance[i][i]*m_covariance[j][j])
	   << "  " << setprecision(precision) << setw(5) << right << i
	   << "  " << setprecision(precision) << setw(5) << right << j <<  endl;
   
  fout.close(); cout << endl; coutCBL << "I wrote the file: " << file_out << endl;
}


// ======================================================================================


std::shared_ptr<Data> cbl::data::Data1D::cut (const std::vector<bool> mask) const
{
  checkDim(mask, m_ndata, "mask");

  vector<double> xx;
  for (size_t i=0; i<mask.size(); i++) 
    if (mask[i])
      xx.push_back(m_x[i]);

  vector<double> data, error;
  vector<vector<double>> covariance;
  Data::cut(mask, data, error, covariance);

  shared_ptr<Data> dd = make_shared<Data1D>(Data1D(xx, data, covariance));

  return dd;
}



// ======================================================================================


shared_ptr<Data> cbl::data::Data1D::cut (const double xmin, const double xmax) const
{
  vector<bool> mask(m_ndata, true);
  vector<double> xx;
  for (int i=0; i<m_ndata; i++) {
    if ((m_x[i] < xmin) || (m_x[i]>xmax))
      mask[i] = false;
    else
      xx.push_back(m_x[i]);
  }

  vector<double> data, error;
  vector<vector<double>> covariance;
  Data::cut(mask, data, error, covariance);


  shared_ptr<Data> dd = make_shared<Data1D>(Data1D(xx, data, covariance));

  return dd;
}


// ======================================================================================


shared_ptr<Data> cbl::data::Data1D::as_factory ()
{
  shared_ptr<Data> dd = make_shared<Data1D>(Data1D(m_x, m_data, m_covariance));

  return dd;
}

