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


std::vector<double> cbl::modelling::twopt::xi_Wedges (const std::vector<double> rr, const std::vector<int> dataset_wedge, const int nwedges, const std::string model, const std::vector<double> parameter, const std::vector<std::shared_ptr<glob::FuncGrid>> pk_interp, const double prec)
{
  vector<cbl::glob::FuncGrid> interp_Xil(3);
  for (size_t i=0; i<3; i++)
    interp_Xil[i] = Xil_interp(pk_interp[0]->x(), 2*i, model, parameter, pk_interp, prec);

  vector<double> XiW(rr.size());

  for (size_t i=0; i<rr.size(); i++) {
    double mu_min = double(dataset_wedge[i])/nwedges;
    double mu_max = double(dataset_wedge[i]+1)/nwedges;

    double f2 = -0.5*((pow(mu_max, 3)-pow(mu_min, 3))/(mu_max-mu_min)-1.);
    double f4 = 0.125*((7.*(pow(mu_max, 5)-pow(mu_min, 5))-10.*(pow(mu_max, 3)-pow(mu_min, 3)))/(mu_max-mu_min)+3.);

    XiW[i] = (interp_Xil[0](rr[i])+f2*interp_Xil[1](rr[i])+f4*interp_Xil[2](rr[i]));
  }
  
  return XiW;
}


// ============================================================================================


std::vector<std::vector<double>> cbl::modelling::twopt::xi_Wedges (const std::vector<double> rr, const int nwedges, const std::string model, const std::vector<double> parameter, const std::vector<std::shared_ptr<glob::FuncGrid>> pk_interp, const double prec)
{
  vector<vector<double>> Xil = Xi_l(rr, 3, model, parameter, pk_interp, prec);

  vector<vector<double>> XiW(nwedges, vector<double>(rr.size(), 0));

  for (int i=0; i<nwedges; i++) {
    double mu_min = double(i)/nwedges;
    double mu_max = double(i+1)/nwedges;

    double f2 = 0.5*((pow(mu_max, 3)-pow(mu_min, 3))/(mu_max-mu_min)-1.);
    double f4 = 0.125*((7.*(pow(mu_max, 5)-pow(mu_min, 5))-10.*(pow(mu_max, 3)-pow(mu_min, 3)))/(mu_max-mu_min)+3.);

    for (size_t j=0; j<rr.size(); j++) 
      XiW[i][j] = Xil[0][j]+f2*Xil[1][j]+f4*Xil[2][j];
  }

  return XiW;
}


// ============================================================================================


std::vector<double> cbl::modelling::twopt::xiWedges (const std::vector<double> rad, const std::shared_ptr<void> inputs, std::vector<double> &parameter)
{
  // structure contaning the required input data
  shared_ptr<STR_data_model> pp = static_pointer_cast<STR_data_model>(inputs);

  // input parameters
  vector<double> pars;
  vector<shared_ptr<glob::FuncGrid>> pk_interp(2);

  if (pp->Pk_mu_model=="dispersion_dewiggled") {
    pk_interp[0] = pp->func_Pk;
    pk_interp[1] = pp->func_Pk_NW;
  }
  else if (pp->Pk_mu_model=="dispersion_modecoupling") { 
    pk_interp[0] = pp->func_Pk;
    pk_interp[1] = pp->func_Pk1loop;
  }
  else ErrorCBL("Error in cbl::modelling::twopt::xiWedges() of ModelFunction_TwoPointCorrelation_wedges.cpp: the chosen model ("+pp->Pk_mu_model+") is not currently implemented!");
  
  return xi_Wedges(rad, pp->dataset_order, pp->nwedges, pp->Pk_mu_model, { parameter[0], parameter[1], parameter[2], parameter[3], parameter[4]/pp->sigma8_z, parameter[5]/pp->sigma8_z, parameter[6] }, pk_interp, pp->prec);
  
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

  vector<vector<double>> xiww = cbl::XiWedges_AP({0., 0.5}, {0.5, 0.5}, alpha_perpendicular, alpha_parallel,  pp->rr, pp->func_multipoles[0], pp->func_multipoles[1], pp->func_multipoles[2]);

  cbl::glob::FuncGrid xiperp(pp->rr, xiww[0], "Spline");
  cbl::glob::FuncGrid xipar(pp->rr, xiww[1], "Spline");

  vector<vector<double>> Xiw(2);

  for (size_t i=0; i<new_rad[0].size(); i++)
    Xiw[0].push_back(parameter[2]*parameter[2]*xiperp(new_rad[0][i])+parameter[4]+parameter[6]/new_rad[0][i]+parameter[8]/(new_rad[0][i]*new_rad[0][i]));

  for (size_t i=0; i<new_rad[1].size(); i++)
    Xiw[1].push_back(parameter[3]*parameter[3]*xipar(new_rad[1][i])+parameter[5]+parameter[7]/new_rad[1][i]+parameter[9]/(new_rad[1][i]*new_rad[1][i]));

  vector<double> Xi;

  for (size_t i=0; i<Xiw.size(); i++)
    for (size_t j=0; j<Xiw[i].size(); j++)
      Xi.push_back(Xiw[i][j]);

  return Xi;
}
