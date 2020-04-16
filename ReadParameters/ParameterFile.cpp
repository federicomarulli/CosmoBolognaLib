/********************************************************************
 *  Copyright (C) 2010 by Federico Marulli and Tommaso Ronconi      *
 *  federico.marulli3@unibo.it, tommaso.ronconi@outlook.it          *
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
 *  @file ParameterFile/ParameterFile.cpp
 *
 *  @brief Methods of the class ParameterFile used to read parameter
 *  files
 *
 *  This file contains the implementation of the methods of the class
 *  ParameterFile used to read a generic parameter file
 *
 *  @authors Tommaso Ronconi, Federico Marulli
 *
 *  @authors tommaso.ronconi@outlook.it, federico.marulli3@unibo.it
 */

#include "ParameterFile.h"

using namespace std;

using namespace cbl;
using namespace glob;


// ============================================================================


std::string cbl::glob::ParameterFile::m_trim (const string inStr)
{
  if (inStr.find_first_not_of(' ') >= inStr.find_last_not_of(' ')+1)
    return ""; //CHECK
  return inStr.substr(inStr.find_first_not_of(' '), inStr.find_last_not_of(' ')+1);
}


// ============================================================================


vector<string> cbl::glob::ParameterFile::m_trim_vect (const string inStr)
{
  string str = inStr;
  str = m_trim(str);
  vector<string> vect;
  str = str.substr(str.find_first_not_of('{'), str.find_last_not_of('}'));
  while (str.find(str.front()) != string::npos) {
    string val;
    if (str.find(',') != string::npos) {
      val = m_trim(str.substr(str.find(str.front()), str.find_first_of(',')));
      str.erase(str.find(str.front()), str.find_first_of(',')+1);
    }
    else {
      val = m_trim(str.substr(str.find(str.front()), string::npos));
      str.erase(str.find(str.front()), string::npos);
    }
    vect.emplace_back(val);
  }
  return vect;
}


// ============================================================================


cbl::glob::ParameterFile::ParameterFile (const std::string parameter_file)
{      
  read(parameter_file);
}


// ============================================================================


void cbl::glob::ParameterFile::read (const std::string parameter_file)
{
  // open the input parameter file
  ifstream fin(parameter_file.c_str()); checkIO(fin, parameter_file);
  
  // read all lines and skip comments (#)
  string line;
  while (getline(fin, line)) {

    // remove potential CR or NL or comments at end of line #
    line = line.substr(0, line.find("\r"));
    line = line.substr(0, line.find("\n"));
    line = line.substr(0, line.find("#"));

    if (line.size()>0) {
      string::size_type eqpos = line.find('=');

      if (eqpos!=string::npos) {
	string key = m_trim(line.substr(0, eqpos));
	string val = line.substr(eqpos+1, string::npos);

	m_keys.emplace_back(key);
	    
	// trim and store
	if (val.find('{') == string::npos) 
	{
	  std::string value = m_trim(val);
	  if (m_trim(val) == "")
	    m_parameters[key] = {};
	  else
	    m_parameters[key] = {m_trim(val)};
	}
	else {
	  vector<string> vval = m_trim_vect(val);
	  m_parameters[key] = vval;
	}
	  
      }
    }
  }
  fin.clear(); fin.close();
}


// ============================================================================


void cbl::glob::ParameterFile::write (const std::string parameter_file)
{
  // open the input parameter file
  ofstream fout(parameter_file.c_str()); checkIO(fout, parameter_file);

  for (auto &&key : m_keys) 
  {
    vector<string> values = m_parameters.at(key);

    fout << key << " = ";

    if (values.size() == 0) 
      fout << endl;
    else if (values.size() == 1)
      fout << values[0] << endl;
    else {
      fout << "{ " << values[0];
      for(size_t i=1; i<values.size(); i++)
	fout << ", " << values[i];
      fout << " }" << endl;
    }
  }

  fout.clear(); fout.close();
}


// ============================================================================


std::vector<std::string> & cbl::glob::ParameterFile::operator[] (const std::string key)
{
  // Add a new entry
  if (m_parameters.find(key) == m_parameters.end()) 
    m_keys.emplace_back(key);

  return m_parameters[key];
}


// ============================================================================


void cbl::glob::ParameterFile::set_key (const std::string key, std::string value, const size_t pos)
{
  if (m_parameters.find(key) == m_parameters.end()) 
    ErrorCBL("Key "+key+" has not been found, exiting.", "cbl::cbl::glob::ParameterFile::operator[]", "ReadParameters/ParameterFile.cpp");
  
  vector<string> &vv = m_parameters.at(key);

  if (vv.size()>pos)
    m_parameters[key][pos] = value;
  else if (vv.size()==pos) {
    WarningMsgCBL("Found "+conv(pos, par::fINT)+" values for key "+key+" == "+conv(pos, par::fINT)+". I will add one parameter.", "cbl::cbl::glob::ParameterFile::set_key()", "ReadParameters/ParameterFile.cpp");
    m_parameters[key].push_back(value);
  }
  else
    ErrorCBL("Found "+conv(vv.size(), par::fINT)+" values for key "+key+" < "+conv(pos, par::fINT)+". Exiting.", "cbl::cbl::glob::ParameterFile::set_key()", "ReadParameters/ParameterFile.cpp");
}


// ============================================================================


void cbl::glob::ParameterFile::set_key (const std::string key, const std::vector<std::string> values)
{
  if (m_parameters.find(key) == m_parameters.end()) 
    ErrorCBL("Key "+key+" has not been found, exiting.", "cbl::cbl::glob::ParameterFile::operator[]", "ReadParameters/ParameterFile.cpp");
  
  m_parameters[key].erase(m_parameters[key].begin(), m_parameters[key].end());
  m_parameters[key] = values;
}


// ============================================================================


std::string cbl::glob::ParameterFile::get_key (const std::string key, const std::string default_value, const size_t pos) const
{
  if (m_parameters.find(key) == m_parameters.end()) {
    WarningMsgCBL("Key "+key+" has not been found, I will return default_value.", "cbl::cbl::glob::ParameterFile::get_key", "ReadParameters/ParameterFile.cpp");
    return default_value;
  }

  const vector<string> &vv = m_parameters.at(key);

  if (vv.size() <= pos) {
    WarningMsgCBL("Found "+conv(pos, par::fINT)+" values for key "+key+" <= "+conv(pos, par::fINT)+", I will return default_values.", "cbl::cbl::glob::ParameterFile::get_key()", "ReadParameters/ParameterFile.cpp");
    return default_value;
  }

  return vv[pos];
}


// ============================================================================


std::vector<std::string> cbl::glob::ParameterFile::get_key (const std::string key, const std::vector<std::string> default_values) const
{
  if (m_parameters.find(key) == m_parameters.end()) {
    WarningMsgCBL("Key "+key+" has not been found, I will return default_values.", "cbl::cbl::glob::ParameterFile::operator[]", "ReadParameters/ParameterFile.cpp");
    return default_values;
  }

  return m_parameters.at(key);
}
