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
 *  @file Data/Table.cpp
 *
 *  @brief Methods of the class Table
 *
 *  This file contains the implementation of the methods of the class
 *  Table
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#include "Table.h"

using namespace std;

using namespace cbl;
using namespace data;


// ======================================================================================


void cbl::data::Table::m_set (const std::vector<std::string> names, const std::vector<std::vector<double>> values)
{
  this->insert(names, values);
}


// ======================================================================================


vector<string> cbl::data::Table::m_keys ()
{
  vector<string> keys;

  for(table_map::iterator it = m_table.begin(); it != m_table.end(); ++it) 
    keys.push_back(it->first);

  return keys;
}


// ======================================================================================


cbl::data::Table::Table (const std::vector<std::string> names, const std::vector<std::vector<double>> values)
{
  m_set(names, values);
}


// ======================================================================================


cbl::data::Table::Table (const std::string input_dir, const std::string input_file, const std::vector<std::string> names, const vector<size_t> use_cols, const size_t header_lines_to_skip)
{
  this->read(input_dir, input_file, names, use_cols, header_lines_to_skip);
}


// ======================================================================================


vector<double> cbl::data::Table::operator[] (const std::string name)
{
  table_map::iterator it = m_table.find(name);

  if (it == m_table.end())
    ErrorCBL("Column "+name+" has not been found, exiting.", "cbl::data::Table::operator[]", "Data/Table.cpp");

  return m_table.at(name);
}


// ======================================================================================


vector<vector<double>> cbl::data::Table::operator[] (const std::vector<std::string> names)
{
  vector<vector<double>> columns;

  for (size_t i=0; i<names.size(); i++)
    columns.push_back(this->operator[](names[i]));

  return columns;
}


// ======================================================================================


void cbl::data::Table::insert (const std::string name, const std::vector<double> values, const bool replace)
{
  table_map::iterator it = m_table.find(name);

  if (it == m_table.end())
    m_table.emplace(name, values);
  else {
    if (replace) {
      m_table.erase(name);
      m_table.emplace(name, values);
    }
    else
      cbl::ErrorCBL("Column "+name+" exists and replace is false. Check your inputs. Exiting.", "cbl::data::Table::insert", "Data/Table.cpp");
  }
}


// ======================================================================================


void cbl::data::Table::insert (const std::vector<std::string> names, const std::vector<std::vector<double>> values, const bool replace)
{
  for (size_t i=0; i<names.size(); i++)
    this->insert(names[i], values[i], replace);
}


// ======================================================================================


void cbl::data::Table::read (const std::string input_dir, const std::string input_file, const std::vector<std::string> names, const vector<size_t> use_cols, const size_t header_lines_to_skip)
{
  string file = input_dir+input_file;
  coutCBL << "Reading the table " << file << endl;

  const size_t ncols = names.size();

  vector<vector<double>> values;

  if (names.size() != use_cols.size() && use_cols.size()>0)
    ErrorCBL("Wrong number of columns, check your inputs!", "cbl::data::Table::read", "Data/Table.cpp");
  	
  std::function<vector<double>(const vector<double> params)> assign;

  if (use_cols.size()==0) {
    assign = [&] (const vector<double> params) {
      return params;
    };
  }
  else {
    checkDim(use_cols, ncols, "columns");
    assign = [&] (const vector<double> params) {
      vector<double> pp(ncols);
      for (size_t i=0; i<ncols; i++)
	pp[i] = params[use_cols[i]];
      return pp;
    };
  }

  ifstream fin(file.c_str()); checkIO(fin, file);

  string line;
  for (size_t i=0; i<header_lines_to_skip; i++)
    getline(fin, line);

  while (getline(fin, line)) {
    stringstream ss(line);
    double NUM;
    vector<double> ll;

    while (ss>>NUM) ll.push_back(NUM);

    values.push_back(assign(ll));
  }

  fin.clear(); fin.close();

  m_set(names, cbl::transpose(values));

  coutCBL << "Done!" << endl << endl;
}


// ======================================================================================


void cbl::data::Table::write (const std::string output_dir, const std::string output_file, const std::vector<std::string> names)
{
  vector<string> keys = (names.size()==0) ? this->m_keys() : names;
  vector<vector<double>> values = transpose(this->operator[](keys)); 

  string MK = "mkdir -p "+output_dir;
  if (system(MK.c_str())) {}
  
  string file = output_dir+output_file; 

  ofstream fout(file.c_str()); checkIO(fout, file);
  fout.precision(10);

  // Write the header
  fout << "#";
  for (size_t i=0; i<keys.size(); i++)
    fout << " " << keys[i];
  fout << endl;

  for (size_t i=0; i<values.size(); i++) {
    for (size_t j=0; j<values[i].size(); j++)
      cbl::Print(values[i][j], 5, 14, "", "  ", false, fout);
    fout << endl;
  }

  fout.clear(); fout.close();

  coutCBL << "I wrote the file: " << file << endl;
}
