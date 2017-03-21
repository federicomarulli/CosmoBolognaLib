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
 *  @file Func/Data2D_extra.cpp
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

using namespace cosmobl;
using namespace data;
using namespace glob;


// ======================================================================================


void cosmobl::data::Data2D_extra::read (const string input_file, const int skipped_lines)
{
  (void)skipped_lines;
  
  ErrorCBL("Error in cosmobl::data::Data2D_extra::read : work in progress!", ExitCode::_workInProgress_);

  ifstream fin(input_file.c_str()); checkIO(fin, input_file);
  
  fin.clear(); fin.close();
}


// ======================================================================================


void cosmobl::data::Data2D_extra::write (const string dir, const string file, const string header, const bool full, const int rank) const 
{
  (void)rank;
  
  string file_out = dir+file;
  ofstream fout(file_out.c_str()); checkIO(fout, file_out);

  fout << "### " << header << " ###" << endl;

  int index = 0;
  for (size_t i=0; i<m_x.size(); ++i)
    for (size_t j=0; j<m_y.size(); ++j) {
      fout << setiosflags(ios::fixed) << setprecision(4) << setw(8) << m_x[i] << "  " << setw(8) << m_y[j] << "  " << setw(8) << m_fxy[i][j] << "  " << setw(8) << m_error_fxy[i][j];
      for (size_t ex=0; ex<m_extra_info.size(); ++ex)
	fout << "  " << setw(8) << m_extra_info[ex][index];
      fout << endl;
      index ++;
    }
  
  if (full) { // duplicate the information in the other 3 quadrants

    index = 0;
    for (size_t i=0; i<m_x.size(); ++i)
      for (size_t j=0; j<m_y.size(); ++j) {
	fout << setiosflags(ios::fixed) << setprecision(4) << setw(8) << m_x[i] << "  " << setw(8) << -m_y[j] << "  " << setw(8) << m_fxy[i][j] << "  " << setw(8) << m_error_fxy[i][j] << "  ";
	for (size_t ex=0; ex<m_extra_info.size(); ++ex)
	  fout << setw(8) << m_extra_info[ex][index] << "  ";
	fout << endl;
	index ++;
      }
    
    index = 0;
    for (size_t i=0; i<m_x.size(); ++i)
      for (size_t j=0; j<m_y.size(); ++j) { 
	fout << setiosflags(ios::fixed) << setprecision(4) << setw(8) << -m_x[i] << "  " << setw(8) << -m_y[j] << "  " << setw(8) << m_fxy[i][j] << "  " << setw(8) << m_error_fxy[i][j] << "  ";
	for (size_t ex=0; ex<m_extra_info.size(); ++ex)
	  fout << setw(8) << m_extra_info[ex][index] << "  ";
	fout << endl;
	index ++;
      }

    index = 0;
    for (size_t i=0; i<m_x.size(); ++i)
      for (size_t j=0; j<m_y.size(); ++j) {
	fout << setiosflags(ios::fixed) << setprecision(4) << setw(8) << -m_x[i] << "  " << setw(8) << m_y[j] << "  " << setw(8) << m_fxy[i][j] << "  " << setw(8) << m_error_fxy[i][j] << "  ";
	for (size_t ex=0; ex<m_extra_info.size(); ++ex)
	  fout << setw(8) << m_extra_info[ex][index] << "  ";
	fout << endl;
	index ++;
      }
  }

  fout.close(); cout << endl; coutCBL << "I wrote the file: " << file_out << endl << endl;

}
