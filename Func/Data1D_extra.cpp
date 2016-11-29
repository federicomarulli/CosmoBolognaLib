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
 *  @file Func/Data1D_extra.cpp
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

using namespace cosmobl;
using namespace data;


// ======================================================================================


void cosmobl::data::Data1D_extra::read (const string input_file, const int skip_nlines)
{
  ifstream fin(input_file.c_str()); checkIO(fin, input_file);
  string line;

  if (skip_nlines>0)
    for (int i=0; i<skip_nlines; ++i)
      getline(fin, line);

  bool first = true;
  
  while (getline(fin, line)) {
    
    stringstream ss(line); 
    vector<double> num; double NUM;
    while (ss>>NUM) num.push_back(NUM);
    
    m_x.push_back(num[0]);
    m_fx.push_back(num[1]);
    m_error_fx.push_back(num[2]);

    if (first) m_extra_info.resize(num.size()-3);
    
    for (size_t ex=3; ex<num.size(); ++ex) 
      m_extra_info[ex].emplace_back(num[ex]);

    first = false;
  }

  fin.clear(); fin.close(); 
}


// ======================================================================================


void cosmobl::data::Data1D_extra::write (const string dir, const string file, const string header, const int rank) const 
{
  (void)rank;
  
  string file_out = dir+file;
  ofstream fout(file_out.c_str()); checkIO(fout, file_out);
  
  fout << "### "<< header <<" ###" << endl;

  for (size_t i=0; i<m_x.size(); ++i) { 
    fout << setiosflags(ios::fixed) << setprecision(4) << setw(8) << m_x[i] << "  " << setw(8) << m_fx[i] << "  " << setw(8) << m_error_fx[i];
    for (size_t ex=0; ex<m_extra_info.size(); ++ex)
      fout << "  " << setw(8) << m_extra_info[ex][i];
    fout << endl;
  }
  
  fout.close(); cout << endl ; coutCBL << endl << "I wrote the file: " << file_out << endl;
}
