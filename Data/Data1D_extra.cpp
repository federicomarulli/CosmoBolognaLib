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
 *  @authors federico.marulli3@unbo.it
 */

#include "Data1D_extra.h"

using namespace std;

using namespace cbl;
using namespace data;


// ======================================================================================


void cbl::data::Data1D_extra::read (const string input_file, const int skip_nlines)
{
  ifstream fin(input_file.c_str()); checkIO(fin, input_file);
  string line;
  
  bool moreColumns = true;
  int col = 3;
  int ind1 = 1, ind2 = 2;
  
  // continue to read in case of more the 3 columns (e.g. for clustering multipoles)
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

      // if there are more than 3 columns (e.g. for clustering
      // multipoles) then the next columns will be read and added to
      // the first one
      if (num.size()==(size_t)col+m_extra_info.size()) moreColumns = false;
    
      m_x.emplace_back(num[0]);
      m_data.emplace_back(num[ind1]);
      m_error.emplace_back(num[ind2]);
      
      int ind = 0;
      for (size_t ex=num.size()-m_extra_info.size(); ex<num.size(); ++ex) 
	m_extra_info[ind++].emplace_back(num[ex]);

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


void cbl::data::Data1D_extra::write (const string dir, const string file, const string header, const int precision, const int rank) const 
{
  (void)rank;
  
  string file_out = dir+file;
  ofstream fout(file_out.c_str()); checkIO(fout, file_out);
  
  fout << "### "<< header <<" ###" << endl;

  for (size_t i=0; i<m_x.size(); ++i) {
    fout << setiosflags(ios::fixed) << setprecision(precision) << setw(15) << right << m_x[i] 
	 << "  " << setiosflags(ios::fixed) << setprecision(precision) << setw(15) << right << m_data[i] 
	 << "  " << setiosflags(ios::fixed) << setprecision(precision) << setw(15) << right << m_error[i];
    
    for (size_t ex=0; ex<m_extra_info.size(); ++ex)
      fout << "  " << setiosflags(ios::fixed) << setprecision(precision) << setw(15) << m_extra_info[ex][i];
    fout << endl;
  }
  
  fout.close(); cout << endl ; coutCBL << endl << "I wrote the file: " << file_out << endl;
}
