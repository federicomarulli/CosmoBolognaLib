/*******************************************************************
 *  Copyright (C) 2010 by Federico Marulli                         *
 *  federico.marulli3@unibo.it                                     *
 *                                                                 *
 *  This program is free software; you can redistribute it and/or  *
 *  modify it under the terms of the GNU General Public License as *
 *  published by the Free Software Foundation; either version 2 of *
 *  the License, or (at your option) any later version.            *
 *                                                                 *
 *  This program is distributed in the hope that it will be useful,*
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the  *
 *  GNU General Public License for more details.                   *
 *                                                                 *
 *  You should have received a copy of the GNU General Public      *
 *  License along with this program; if not, write to the Free     *
 *  Software Foundation, Inc.,                                     *
 *  59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.      *
 *******************************************************************/

/**
 *  @file Wrappers/FITSwrapper.cpp
 *
 *  @brief functions that wrap FITS routines to
 *  manage .fits files
 *
 *  This file contains the implementation of 
 *  wrappers of FITS routines to
 *  manage .fits files
 *
 *  @author Alfonso Veropalumbo
 *
 *  @author alfonso.veropalumbo@unibo.it
 */

#include "FITSwrapper.h"
#include "CCfits/CCfits"

using namespace std;

using namespace cbl;


// ============================================================================


vector<vector<double>> cbl::wrapper::ccfits::read_table_fits (const std::string input_fits, const std::vector<std::string> column_names, const int next, const double fill_value)
{
  ifstream check_input(input_fits); checkIO(check_input, input_fits); check_input.close();
  
  unique_ptr<CCfits::FITS> pInfile(new CCfits::FITS(input_fits, CCfits::Read, next));

  CCfits::ExtHDU& table = pInfile->currentExtension();

  const long nrows = table.rows();
  size_t no_col = 0;

  if (nrows==0)
    ErrorCBL("no rows in the selected table extension!", "read_table_fits", "FITSwrapper.cpp");

  vector<vector<double>> cc;

  for (size_t i=0; i<column_names.size(); i++) {
    vector<double> vv;
    try {
      table.column(column_names[i]).read(vv, 0, nrows);
    }
    catch (CCfits::Table::NoSuchColumn) {
      if (fill_value==par::defaultDouble)
	ErrorCBL("no column "+column_names[i]+"!", "read_table_fits", "FITSwrapper.cpp");
      else {
	vv.erase(vv.begin(), vv.end());
	vv.resize(nrows, fill_value);
	no_col ++;
      }
    }
    cc.push_back(vv);
  }

  if (no_col==column_names.size())
    ErrorCBL("no column found!", "read_table_fits", "FITSwrapper.cpp");

  return cc;
}


// ============================================================================


void cbl::wrapper::ccfits::write_table_fits (const std::string output_dir, const std::string file_fits, const std::vector<std::string> column_names, const std::vector<std::vector<double>> table, const std::vector<std::string> column_units)
{
  const string file_name = output_dir+file_fits;

  if (system(("rm -rf "+file_name).c_str())) {}
  
  shared_ptr<CCfits::FITS> pFits(0);

  try
  {                
    pFits.reset( new CCfits::FITS(file_name, CCfits::Write) );
  }
  catch(CCfits::FITS::CantOpen)
  {
    ErrorCBL("the file "+file_name+" cannot be opened!", "write_table_fits", "FITSwrapper.cpp");
  }

  const size_t ncolumns = column_names.size();
  size_t rows = table[0].size();     
  string hduName = "TABLE_BINARY";   

  vector<string> colName = column_names;
  vector<string> colUnit = (column_units.size()==0) ? vector<string>(ncolumns, "") :  column_units;
  vector<string> colForm(ncolumns, "1D");

  CCfits::Table* newTable = pFits->addTable(hduName, rows, colName, colForm, colUnit);
    
  // numbers is a scalar column
  for (size_t i=0; i<ncolumns; i++)
    newTable->column(colName[i]).write(table[i], 1);  

  coutCBL << "I wrote the file: " << file_name << endl;
}
