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
 *  @file ReadParameters/ReadParameters.cpp
 *
 *  @brief Methods of the class ReadParameters used to read parameter
 *  files
 *
 *  This file contains the implementation of the methods of the class
 *  ReadParameters used to read a generic parameter file
 *
 *  @authors Tommaso Ronconi, Federico Marulli
 *
 *  @authors tommaso.ronconi@outlook.it, federico.marulli3@unibo.it
 */

#include "ReadParameters.h"

using namespace cosmobl;
using namespace glob;


// ============================================================================


cosmobl::glob::ReadParameters::ReadParameters (const string parameter_file)
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
	string key = line.substr(0, eqpos);
	string val = line.substr(eqpos+1, string::npos);
	    
	// trim and store
	if (val.find('{') == string::npos) m_parameters[m_trim(key)] = m_trim(val);
	else {
	  vector<string> vval = m_trim_vect(val);
	  m_vectors[m_trim(key)] = vval;
	}
	  
      }
    }
  }

  fin.clear(); fin.close();
}


// ============================================================================


string cosmobl::glob::ReadParameters::m_trim (const string inStr)
{
  return inStr.substr(inStr.find_first_not_of(' '), inStr.find_last_not_of(' ')+1);
}


// ============================================================================

vector<string> cosmobl::glob::ReadParameters::m_trim_vect (const string inStr)
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

