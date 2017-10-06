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
 *  Modelling/TwoPointCorrelation/ModelFunction_TwoPointCorrelation_wedges.cpp
 *
 *  @brief Functions to model the wedges of the two-point correlation
 *  function
 *
 *  This file contains the implementation of the functions used to
 *  model the wedges of the two-point correlation function
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#include "ModelFunction_TwoPointCorrelation.h"
#include "ModelFunction_TwoPointCorrelation_wedges.h"

using namespace cosmobl;


// ============================================================================================


vector<double> cosmobl::modelling::twopt::xiWedges (const vector<double> rad, const shared_ptr<void> inputs, vector<double> &parameter)
{
  // structure contaning the required input data
  shared_ptr<STR_data_model> pp = static_pointer_cast<STR_data_model>(inputs);

  const int nrad = rad.size()/pp->nwedges;

  vector<double> new_rad(rad.begin(), rad.begin()+nrad);

  // input parameters

  // AP parameter that contains the distance information
  double alpha_perpendicular = parameter[0];

  // AP parameter that contains the distance information
  double alpha_parallel = parameter[1];

  // f(z)*sigma8(z)
  double fsigma8 = parameter[2];
  
  // bias(z)*sigma8(z)
  double bsigma8 = parameter[3];

  // streaming scale
  double SigmaS = parameter[4];

  // bias
  double bias = bsigma8/pp->sigma8_z;

  // f(z)
  double linear_growth_rate = fsigma8/pp->sigma8_z;

  vector<double> Xi;

  vector<vector<double>> Xiw = Xi_wedges( new_rad, pp->nwedges, alpha_perpendicular, alpha_parallel, pp->sigmaNL_perp, pp->sigmaNL_par, bias, linear_growth_rate, SigmaS, pp->kk, pp->func_Pk, pp->func_Pk_NW, pp->prec);

  for(int i=0; i<pp->nwedges; i++)
    for(int j=0; j<nrad; j++)
      Xi.push_back(Xiw[i][j]);

  return Xi;
}


// ============================================================================================


vector<double> cosmobl::modelling::twopt::xiWedges_BAO (const vector<double> rad, const shared_ptr<void> inputs, vector<double> &parameter)
{
  // structure contaning the required input data
  shared_ptr<STR_data_model> pp = static_pointer_cast<STR_data_model>(inputs);

  vector<vector<double>> new_rad(pp->nwedges);

  for(size_t i=0; i<rad.size(); i++)
    new_rad[pp->dataset_order[i]].push_back(rad[i]);

  // input parameters

  // AP parameter that contains the distance information
  double alpha_perpendicular = parameter[0];

  // AP parameter that contains the distance information
  double alpha_parallel = parameter[1];

  vector<vector<double>> xiww = cosmobl::XiWedges_AP ({0., 0.5}, {0.5, 0.5}, alpha_perpendicular, alpha_parallel,  pp->rr, pp->func_multipoles[0], pp->func_multipoles[1], pp->func_multipoles[2]);

  cosmobl::glob::FuncGrid xiperp(pp->rr, xiww[0], "Spline");
  cosmobl::glob::FuncGrid xipar(pp->rr, xiww[1], "Spline");

  vector<vector<double>> Xiw(2);

  for(size_t i=0; i<new_rad[0].size(); i++)
    Xiw[0].push_back(parameter[2]*parameter[2]*xiperp(new_rad[0][i])+parameter[4]+parameter[6]/new_rad[0][i]+parameter[8]/(new_rad[0][i]*new_rad[0][i]));

  for(size_t i=0; i<new_rad[1].size(); i++)
    Xiw[1].push_back(parameter[3]*parameter[3]*xipar(new_rad[1][i])+parameter[5]+parameter[7]/new_rad[1][i]+parameter[9]/(new_rad[1][i]*new_rad[1][i]));

  vector<double> Xi;

  for (size_t i=0; i<Xiw.size(); i++)
    for (size_t j=0; j<Xiw[i].size(); j++)
      Xi.push_back(Xiw[i][j]);

  return Xi;
}

/*
vector<double> cosmobl::modelling::twopt::xiWedges_BAO (const vector<double> rad, const shared_ptr<void> inputs, vector<double> &parameter)
{
  // structure contaning the required input data
  shared_ptr<STR_data_model> pp = static_pointer_cast<STR_data_model>(inputs);

  vector<vector<double>> new_rad(pp->nwedges);

  for(size_t i=0; i<rad.size(); i++)
      new_rad[pp->dataset_order[i]].push_back(rad[i]);
  
  // input parameters

  // AP parameter that contains the distance information
  double alpha_perpendicular = parameter[0];

  // AP parameter that contains the distance information
  double alpha_parallel = parameter[1];

  vector<vector<double>> Xiw(2);

  for(size_t i=0; i<new_rad[0].size(); i++){

    auto xi_mu = [&] (const double mu)
    {
      double fact =sqrt(mu*mu*alpha_parallel*alpha_parallel+(1-mu*mu)*alpha_perpendicular*alpha_perpendicular);
      double sp = new_rad[0][i]*fact;
      double mup = mu*alpha_parallel/fact;
      double val=0;
      for(int j=0; j<pp->nwedges; j++)
	val += pp->func_multipoles[j]->operator()(sp)*cosmobl::legendre_polynomial(mup, j*2);
      return val;
    };

    Xiw[0].push_back(0.5*parameter[2]*cosmobl::gsl::GSL_integrate_qag(xi_mu, 0, 0.5)+parameter[4]+parameter[6]/new_rad[0][i]+parameter[8]/(new_rad[0][i]*new_rad[0][i]));
  }

  for(size_t i=0; i<new_rad[1].size(); i++){

    auto xi_mu = [&] (const double mu)
    {
      double fact =sqrt(mu*mu*alpha_parallel*alpha_parallel+(1-mu*mu)*alpha_perpendicular*alpha_perpendicular);
      double sp = new_rad[1][i]*fact;
      double mup = mu*alpha_parallel/fact;
      double val=0;
      for(int j=0; j<pp->nwedges; j++)
	val += pp->func_multipoles[j]->operator()(sp)*cosmobl::legendre_polynomial(mup, j*2);
      return val;
    };

    Xiw[1].push_back(0.5*parameter[3]*cosmobl::gsl::GSL_integrate_qag(xi_mu, 0.5, 1)+parameter[5]+parameter[7]/new_rad[1][i]+parameter[9]/(new_rad[1][i]*new_rad[1][i]));
  }

  vector<double> Xi;

  for (size_t i=0; i<Xiw.size(); i++)
    for (size_t j=0; j<Xiw[i].size(); j++)
      Xi.push_back(Xiw[i][j]);

  return Xi;
}
*/
