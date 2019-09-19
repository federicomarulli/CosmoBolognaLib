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
 *  @file Statistics/Model1D.cpp
 *
 *  @brief Methods of the class Model1D
 *
 *  This file contains the implementation of the methods of the class
 *  Model1D
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#include "Model1D.h"

using namespace std;

using namespace cbl;


// ======================================================================================


void cbl::statistics::Model1D::set_function (const model_function_1D function)
{
  m_function = [this, function](vector<vector<double>> xx, shared_ptr<void> inputs, vector<double> &parameters) 
    {
      vector<vector<double>> res(1, function(xx[0], inputs, parameters));
      return res;
    };
}


// ======================================================================================


void cbl::statistics::Model1D::set_function (const std::vector<double> (*function)(const std::vector<double> xx, std::vector<double> &val))
{
  m_function = [this, function](vector<vector<double>> xx, shared_ptr<void> inputs, vector<double> &parameters) 
    {
      (void)inputs;
      vector<vector<double>> res(1, function(xx[0], parameters));
      return res;
    };
}


// ======================================================================================


void cbl::statistics::Model1D::stats_from_chains (const vector<double> xx, vector<double> &median_model, vector<double> &low_model, vector<double> &up_model, const int start, const int thin) 
{
  vector<vector<double>> _median, _low, _up;
  
  Model::stats_from_chains({1, xx}, _median, _low, _up, start, thin);

  median_model = _median[0];
  low_model = _low[0];
  up_model = _up[0];
}


// ======================================================================================


void cbl::statistics::Model1D::write (const string output_dir, const string output_file, const vector<double> xx, const vector<double> parameters)
{
  vector<double> pp = parameters;
  vector<double> xx_unique = different_elements(xx);

  if (xx.size() % xx_unique.size()!=0)
    ErrorCBL("model.size() is not a multiple of xx.size()!", "write", "Model1D.cpp");

  const int nmodels = xx.size()/xx_unique.size();
  vector<double> model = this->operator()(xx, pp);

  const string mkdir = "mkdir -p "+output_dir;
  if (system(mkdir.c_str())) {}
  
  const string file = output_dir+output_file;

  ofstream fout(file.c_str()); checkIO(fout, file);

  fout << "### [1] x";
  for (int nn=0; nn<nmodels; nn++) fout << " # [" << nn+2 << "] y_" << nn << "(x)";
  fout << " ###" << endl;
  
  for (size_t i=0; i<xx_unique.size(); ++i) {
    fout << setprecision(5) << setw(10) << right << xx_unique[i]; 
    for (int nn=0; nn<nmodels; nn++)
      fout << "  " << setprecision(5) << setw(10) << right << model[i+nn*xx_unique.size()];
    fout << endl;
  }
  
  fout.clear(); fout.close(); coutCBL << "I wrote the file: " << file << endl << endl;
}


// ======================================================================================


void cbl::statistics::Model1D::write_at_bestfit (const string output_dir, const string output_file, const vector<double> xx)
{
  vector<double> bf = m_parameters->bestfit_values();
  write(output_dir, output_file, xx, bf);  
}


// ======================================================================================


void cbl::statistics::Model1D::write_from_chains (const string output_dir, const string output_file, const vector<double> xx, const int start, const int thin)
{
  vector<double> xx_unique = different_elements(xx);

  if (xx.size() % xx_unique.size()!=0)
    ErrorCBL("model.size() is not a multiple of xx.size()!", "write_from_chains", "Model1D.cpp");

  const int nmodels = xx.size()/xx_unique.size();

  vector<double> median_model, low_model, up_model;
  stats_from_chains(xx, median_model, low_model, up_model, start, thin);

  const string mkdir = "mkdir -p "+output_dir;
  if (system(mkdir.c_str())) {}
  
  const string file = output_dir+output_file;

  ofstream fout(file.c_str()); checkIO(fout, file);
  
  fout << "### [1] x";
  for (int nn=0; nn<nmodels; nn++) fout << " # [" << nn+2 << "] median y(x) # [" << nn+3 << "] 16% percentile y(x)# [" << nn+4 << "] 84% percentile y(x)";
  fout << " ###" << endl;
  
  for (size_t i=0; i<xx_unique.size(); i++) {
    fout << setprecision(5) << setw(10) << right << xx_unique[i] << "  "; 
    for (int nn=0; nn<nmodels; nn++)
      fout << setprecision(5) << setw(10) << right << median_model[i+nn*xx_unique.size()] << "  "
	   << setprecision(5) << setw(10) << right << low_model[i+nn*xx_unique.size()] << "  "
	   << setprecision(5) << setw(10) << right << up_model[i+nn*xx_unique.size()] << "  ";
    fout << endl;
  }
  
  fout.clear(); fout.close(); coutCBL << "I wrote the file: " << file << endl;
}
