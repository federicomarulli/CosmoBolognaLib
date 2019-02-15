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
 *  @file Statistics/Model2D.cpp
 *
 *  @brief Methods of the class Model2D
 *
 *  This file contains the implementation of the methods of the class
 *  Model2D
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#include "Model2D.h"

using namespace std;

using namespace cbl;


// ======================================================================================


void cbl::statistics::Model2D::set_function (const model_function_2D function)
{
  m_function = [this, function](vector<vector<double>> xx, shared_ptr<void> inputs, vector<double> &parameters) 
    {
      return function(xx[0], xx[1], inputs, parameters);
    };
}


// ======================================================================================


void cbl::statistics::Model2D::stats_from_chains (const vector<double> xx, const vector<double> yy, vector<vector<double>> &median_model, vector<vector<double>> &low_model, vector<vector<double>> &up_model, const int start, const int thin) 
{
  Model::stats_from_chains({xx, yy}, median_model, low_model, up_model, start, thin);
}


// ======================================================================================


void cbl::statistics::Model2D::write (const string output_dir, const string output_file, const vector<double> xx, const vector<double> yy, const vector<double> parameters)
{
  vector<double> pp = parameters;
  vector<vector<double>> model = this->operator()(xx, yy, pp);

  const string mkdir = "mkdir -p "+output_dir;
  if (system(mkdir.c_str())) {}

  const string file = output_dir+output_file;
  
  ofstream fout(file.c_str()); checkIO(fout, file);

  fout << "### [1] x # [2] y # [3] model(x,y) ###";
  
  for (size_t i=0; i<xx.size(); i++)
    for (size_t j=0; j<yy.size(); j++)
      fout << setprecision(5) << setw(10) << right << xx[i] << "  "
	   << setprecision(5) << setw(10) << right << yy[j] << "  "
	   << setprecision(5) << setw(10) << right << model[i][j] << endl;
  
  fout.clear(); fout.close(); coutCBL << "I wrote the file: " << file << endl << endl;
}


// ======================================================================================


void cbl::statistics::Model2D::write_at_bestfit (const string output_dir, const string output_file, const vector<double> xx, const vector<double> yy)
{
  vector<double> bf = m_parameters->bestfit_values();
  write(output_dir, output_file, xx, yy, bf);  
}


// ======================================================================================


void cbl::statistics::Model2D::write_from_chains (const string output_dir, const string output_file, const vector<double> xx, const vector<double> yy, const int start, const int thin)
{
  vector<vector<double>> median_model, low_model, up_model;
  stats_from_chains(xx, yy, median_model, low_model, up_model, start, thin);

  const string mkdir = "mkdir -p "+output_dir;
  if (system(mkdir.c_str())) {}
  
  const string file = output_dir+output_file;
  
  ofstream fout(file.c_str()); checkIO(fout, file);

  fout << "### [1] x # [2] y # [3] median model(x,y) # [4] 16% percentile model(x,y) # [5] 84% percentile ###";
  
  for (size_t i=0; i<xx.size(); i++)
    for (size_t j=0; j<yy.size(); j++)
      fout << setprecision(5) << setw(10) << right << xx[i] << "  "
	   << setprecision(5) << setw(10) << right << yy[j] << "  "
	   << setprecision(5) << setw(10) << right << median_model[i][j] << "  "
	   << setprecision(5) << setw(10) << right << low_model[i][j] << "  "
	   << setprecision(5) << setw(10) << right << up_model[i][j] << endl;
 
  fout.clear(); fout.close(); coutCBL << "I wrote the file: " << file << endl;
}

