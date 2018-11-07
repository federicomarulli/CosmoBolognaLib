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

using namespace std;

using namespace cbl;


// ============================================================================================


std::vector<double> cbl::modelling::twopt::xiWedges (const std::vector<double> rad, const std::shared_ptr<void> inputs, std::vector<double> &parameter)
{
  // structure contaning the required input data
  shared_ptr<STR_data_model> pp = static_pointer_cast<STR_data_model>(inputs);

  // input parameters
  vector<double> pars;
  vector<shared_ptr<glob::FuncGrid>> pk_interp(2);

  if (pp->Pk_mu_model == 0) { // de-wiggled model 
    pk_interp[0] = pp->func_Pk;
    pk_interp[1] = pp->func_Pk_NW;
  }
  else if (pp->Pk_mu_model == 1) { // mode-coupling model 
    pk_interp[0] = pp->func_Pk;
    pk_interp[1] = pp->func_Pk1loop;
  }

  vector<double> Xiw = Xi_wedges(rad, pp->dataset_order, pp->nwedges, pp->Pk_mu_model, parameter, pk_interp, pp->prec);

  const double norm = pow(pp->sigma8_z, -2);
  
  for (size_t i=0; i<Xiw.size(); i++)
      Xiw[i] *= norm;

  return Xiw;

}


// ============================================================================================


std::vector<double> cbl::modelling::twopt::xiWedges_BAO (const std::vector<double> rad, const std::shared_ptr<void> inputs, std::vector<double> &parameter)
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

  vector<vector<double>> xiww = cbl::XiWedges_AP ({0., 0.5}, {0.5, 0.5}, alpha_perpendicular, alpha_parallel,  pp->rr, pp->func_multipoles[0], pp->func_multipoles[1], pp->func_multipoles[2]);

  cbl::glob::FuncGrid xiperp(pp->rr, xiww[0], "Spline");
  cbl::glob::FuncGrid xipar(pp->rr, xiww[1], "Spline");

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
