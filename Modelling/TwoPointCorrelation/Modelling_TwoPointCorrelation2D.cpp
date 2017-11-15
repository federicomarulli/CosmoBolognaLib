/********************************************************************
 *  Copyright (C) 2016 by Federico Marulli and Alfonso Veropalumbo  *
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
 *  @file
 *  Modelling/TwoPointCorrelation/Modelling_TwoPointCorrelation2D.cpp
 *
 *  @brief Methods of the class Modelling_TwoPointCorrelation2D
 *
 *  This file contains the implementation of the methods of the class
 *  Modelling_TwoPointCorrelation2D, i.e. the common functions to
 *  model the 2D two-point correlation functions of any type
 *
 *  @authors Federico Marulli
 *
 *  @authors federico.marulli3@unbo.it
 */


#include "Modelling_TwoPointCorrelation2D.h"

using namespace cosmobl;


// ============================================================================================


void cosmobl::modelling::twopt::Modelling_TwoPointCorrelation2D::set_data_model (const cosmology::Cosmology cosmology, const double redshift, const string method_Pk, const double sigmaNL, const bool NL, const int FV, const string output_root, const bool bias_nl, const double bA, const bool xiType, const double k_star, const bool xiNL, const double v_min, const double v_max, const int step_v, const int norm, const double r_min, const double r_max, const double k_min, const double k_max, const int step, const double aa, const bool GSL, const double prec, const string file_par)
{
  m_data_model.cosmology = make_shared<cosmology::Cosmology>(cosmology);
  m_data_model.redshift = redshift;
  m_data_model.method_Pk = method_Pk;
  m_data_model.sigmaNL = sigmaNL;
  m_data_model.NL = NL;
  m_data_model.FV = FV;
  m_data_model.output_root = output_root;
  m_data_model.bias_nl = bias_nl;
  m_data_model.bA = bA;
  m_data_model.xiType = xiType;
  m_data_model.k_star = k_star;
  m_data_model.xiNL = xiNL;
  m_data_model.v_min = v_min;
  m_data_model.v_max = v_max;
  m_data_model.step_v = step_v;
  m_data_model.norm = norm;
  m_data_model.r_min = r_min;
  m_data_model.r_max = r_max;
  m_data_model.k_min = k_min;
  m_data_model.k_max = k_max;
  m_data_model.step = step;
  m_data_model.aa = aa;
  m_data_model.GSL = GSL;
  m_data_model.prec = prec;
  m_data_model.file_par = file_par;

  m_data_model.sigma8_z = m_data_model.cosmology->sigma8_Pk(m_data_model.method_Pk, m_data_model.redshift, m_data_model.output_root);
  m_data_model.linear_growth_rate_z = m_data_model.cosmology->linear_growth_rate(m_data_model.redshift);
  m_data_model.var = (1.+ m_data_model.redshift)/m_data_model.cosmology->HH(m_data_model.redshift);
}


// ============================================================================================


void cosmobl::modelling::twopt::Modelling_TwoPointCorrelation2D::write_model (const string dir, const string file, const vector<double> xx, const vector<double> yy, const vector<double> parameter, const int start, const int thin)
{
  vector<double> r1 = (xx.size()==0) ? logarithmic_bin_vector(100, 0.1, 100.) : xx;
  vector<double> r2 = (yy.size()==0) ? logarithmic_bin_vector(100, 0.1, 100.) : yy;


  // compute the model with best-fit parameters
  
  if (parameter.size()==0) {
    
    vector<vector<double>> median_model, low_model, up_model;
    compute_model_from_chains(r1, r2, start, thin, median_model, low_model, up_model); // check!!!
    
  
    string mkdir = "mkdir -p "+dir; if (system(mkdir.c_str())) {}
    string file_out = dir+file;
    ofstream fout(file_out.c_str()); checkIO(fout, file_out);

    fout << "### scale in the first dimension # scale in the second dimension # median model # model at 16th percentile # model at 84th percentile ###" << endl;
    
    for (size_t i=0; i<r1.size(); i++)
      for (size_t j=0; j<r2.size(); j++) 
	fout << setw(10) << setiosflags(ios::fixed) << setprecision(5) << r1[i] << "  " << r2[j] << "  " << median_model[i][j] << "  " << low_model[i][j] << "  " << up_model[i][j] << endl;
    
    fout.clear(); fout.close();
    coutCBL << "I wrote the file: " << dir+file << endl;
  }

  
  // compute the model with input parameters
  
  else{
    vector<double> par = parameter;
    m_model->write_model(dir, file, r1, r2, par);
  }
}


// ============================================================================================


void cosmobl::modelling::twopt::Modelling_TwoPointCorrelation2D::compute_model_from_chains (const vector<double> xx, const vector<double> yy, const int start, const int thin, vector<vector<double>> &median_model, vector<vector<double>> &low_model, vector<vector<double>> &up_model)
{
  (void)xx; (void)yy; (void)start; (void)thin; (void)median_model; (void)low_model; (void)up_model;
  ErrorCBL("Work in progress in compute_model_from_chains() of ModellingTwoPointCorrelation2D.cpp...", glob::_workInProgress_);
}

