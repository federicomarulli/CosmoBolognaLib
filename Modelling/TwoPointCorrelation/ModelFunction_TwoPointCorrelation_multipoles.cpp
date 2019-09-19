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
 *  Modelling/TwoPointCorrelation/ModelFunction_TwoPointCorrelation_multipoles.cpp
 *
 *  @brief Functions to model the multipoles of the two-point
 *  correlation function
 *
 *  This file contains the implementation of the functions used to
 *  model the multipoles of the two-point correlation function
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */


#include "ModelFunction_TwoPointCorrelation.h"
#include "ModelFunction_TwoPointCorrelation_multipoles.h"

using namespace std;

using namespace cbl;


// ============================================================================================


std::vector<double> cbl::modelling::twopt::xiMultipoles (const std::vector<double> rad, const std::shared_ptr<void> inputs, std::vector<double> &parameter)
{
  // structure contaning the required input data
  shared_ptr<STR_data_model> pp = static_pointer_cast<STR_data_model>(inputs);
  
  // input parameters
  vector<shared_ptr<glob::FuncGrid>> pk_interp(2);
  vector<shared_ptr<glob::FuncGrid>> pk_interp_Scoccimarro_CPT(3);
  vector<shared_ptr<glob::FuncGrid>> pk_interp_TNS_CPT(17);
  vector<shared_ptr<glob::FuncGrid>> pk_interp_Dispersion(1);
  std::vector<double> Xi_ll;

  if (pp->Pk_mu_model=="dispersion_dewiggled") { 
    pk_interp[0] = pp->func_Pk;
    pk_interp[1] = pp->func_Pk_NW;
    Xi_ll = Xi_l(rad, pp->dataset_order, pp->use_pole, pp->Pk_mu_model, { parameter[0], parameter[1], parameter[2], parameter[3], parameter[4]/pp->sigma8_z, parameter[5]/pp->sigma8_z, parameter[6] }, pk_interp, pp->prec);
  }
  else if (pp->Pk_mu_model=="dispersion_modecoupling") { 
    pk_interp[0] = pp->func_Pk;
    pk_interp[1] = pp->func_Pk1loop;
    Xi_ll = Xi_l(rad, pp->dataset_order, pp->use_pole, pp->Pk_mu_model, { parameter[0], parameter[1], parameter[2], parameter[3], parameter[4]/pp->sigma8_z, parameter[5]/pp->sigma8_z, parameter[6] }, pk_interp, pp->prec);
  }
  else if (pp->Pk_mu_model=="DispersionGauss" || pp->Pk_mu_model=="DispersionLorentz") {
    pk_interp_Dispersion[0] = pp->func_Pk;
    Xi_ll = Xi_l(rad, pp->dataset_order, pp->use_pole, pp->Pk_mu_model, { parameter[0]/pp->sigma8_z, parameter[1]/pp->sigma8_z, parameter[2] }, pk_interp_Dispersion, pp->prec);
  }  
  else if (pp->Pk_mu_model=="ScoccimarroPezzottaGauss" || pp->Pk_mu_model=="ScoccimarroPezzottaLorentz") {
    pk_interp[0] = pp->func_Pk;
    pk_interp[1] = pp->func_Pk_nonlin;
    Xi_ll = Xi_l(rad, pp->dataset_order, pp->use_pole, pp->Pk_mu_model, { parameter[0]/pp->sigma8_z, parameter[1]/pp->sigma8_z, parameter[2], parameter[3], parameter[4] }, pk_interp, pp->prec);
  }
  else if (pp->Pk_mu_model=="ScoccimarroBelGauss" || pp->Pk_mu_model=="ScoccimarroBelLorentz") {
    pk_interp[0] = pp->func_Pk;
    pk_interp[1] = pp->func_Pk_nonlin;
    Xi_ll = Xi_l(rad, pp->dataset_order, pp->use_pole, pp->Pk_mu_model, { parameter[0]/pp->sigma8_z, parameter[1]/pp->sigma8_z, parameter[2], parameter[3], parameter[4], parameter[5], parameter[6], parameter[7] }, pk_interp, pp->prec);
  }  
  else if (pp->Pk_mu_model=="ScoccimarroGauss" || pp->Pk_mu_model=="ScoccimarroLorentz") {
    pk_interp_Scoccimarro_CPT[0] = pp->func_Pk_DeltaDelta;
    pk_interp_Scoccimarro_CPT[1] = pp->func_Pk_DeltaTheta;
    pk_interp_Scoccimarro_CPT[2] = pp->func_Pk_ThetaTheta;
    Xi_ll = Xi_l(rad, pp->dataset_order, pp->use_pole, pp->Pk_mu_model, { parameter[0]/pp->sigma8_z, parameter[1]/pp->sigma8_z, parameter[2] }, pk_interp_Scoccimarro_CPT, pp->prec);
  }
  else if (pp->Pk_mu_model=="TaruyaGauss" || pp->Pk_mu_model=="TaruyaLorentz") {
    pk_interp_TNS_CPT[0]  = pp->func_Pk_DeltaDelta;
    pk_interp_TNS_CPT[1]  = pp->func_Pk_DeltaTheta;
    pk_interp_TNS_CPT[2]  = pp->func_Pk_ThetaTheta;
    pk_interp_TNS_CPT[3]  = pp->func_Pk_A11;
    pk_interp_TNS_CPT[4]  = pp->func_Pk_A12;
    pk_interp_TNS_CPT[5]  = pp->func_Pk_A22;
    pk_interp_TNS_CPT[6]  = pp->func_Pk_A23;
    pk_interp_TNS_CPT[7]  = pp->func_Pk_A33;
    pk_interp_TNS_CPT[8]  = pp->func_Pk_B12;
    pk_interp_TNS_CPT[9]  = pp->func_Pk_B13;
    pk_interp_TNS_CPT[10] = pp->func_Pk_B14;
    pk_interp_TNS_CPT[11] = pp->func_Pk_B22;
    pk_interp_TNS_CPT[12] = pp->func_Pk_B23;
    pk_interp_TNS_CPT[13] = pp->func_Pk_B24;
    pk_interp_TNS_CPT[14] = pp->func_Pk_B33;
    pk_interp_TNS_CPT[15] = pp->func_Pk_B34;
    pk_interp_TNS_CPT[16] = pp->func_Pk_B44;

    Xi_ll = Xi_l(rad, pp->dataset_order, pp->use_pole, pp->Pk_mu_model, { parameter[0]/pp->sigma8_z, parameter[1]/pp->sigma8_z, parameter[2] }, pk_interp_TNS_CPT, pp->prec);    
  }

  else ErrorCBL("the chosen model ("+pp->Pk_mu_model+") is not currently implemented!", "xiMultipoles", "ModelFunction_TwoPointCorrelation_multipoles.cpp");
  
  return Xi_ll;
}


// ============================================================================================


std::vector<double> cbl::modelling::twopt::xiMultipoles_sigma8_bias (const std::vector<double> rad, const std::shared_ptr<void> inputs, std::vector<double> &parameter)
{
  // structure contaning the required input data
  shared_ptr<STR_data_model> pp = static_pointer_cast<STR_data_model>(inputs);

  // input parameters

  // sigma8(z)
  double sigma8 = parameter[0];

  // bias
  double bias = parameter[1];

  // AP parameter that contains the distance information
  double alpha_perpendicular = 1.;

  // AP parameter that contains the distance information
  double alpha_parallel = 1.;

  // f(z)*sigma8(z)
  double fsigma8 = pp->linear_growth_rate_z*sigma8;
  
  // bias(z)*sigma8(z)
  double bsigma8 = bias*sigma8;

  // streaming scale
  double SigmaS = 0.;

  return Xi_l(rad, pp->dataset_order, pp->use_pole, pp->Pk_mu_model, {alpha_perpendicular, alpha_parallel, pp->sigmaNL_perp, pp->sigmaNL_par, fsigma8/pp->sigma8_z, bsigma8/pp->sigma8_z, SigmaS}, {pp->func_Pk, pp->func_Pk_NW}, pp->prec);
}


// ============================================================================================


std::vector<double> cbl::modelling::twopt::xiMultipoles_BAO (const std::vector<double> rad, const std::shared_ptr<void> inputs, std::vector<double> &parameter)
{
  // structure contaning the required input data
  shared_ptr<STR_data_model> pp = static_pointer_cast<STR_data_model>(inputs);

  vector<vector<double>> new_rad(pp->nmultipoles);

  for(size_t i=0; i<rad.size(); i++)
    new_rad[pp->dataset_order[i]].push_back(rad[i]);
  
  // input parameters

  // AP parameter that contains the distance information
  double alpha_perpendicular = parameter[0];

  // AP parameter that contains the distance information
  double alpha_parallel = parameter[1];

  vector<vector<double>> Xil(2);

  for(size_t i=0; i<new_rad[0].size(); i++) {

    auto xi_mu_0 = [&] (const double mu)
      {
	double alpha = sqrt(mu*mu*alpha_parallel*alpha_parallel+(1-mu*mu)*alpha_perpendicular*alpha_perpendicular);
	double sp = new_rad[0][i]*alpha;
	double mup = mu*alpha_parallel/alpha;
	double val=0;
	for (int j=0; j<pp->nmultipoles; j++)
	  val += pp->func_multipoles[j]->operator()(sp)*cbl::legendre_polynomial(mup, j*2);
	return val;
      };

    double xi0 = cbl::wrapper::gsl::GSL_integrate_qag(xi_mu_0, 0, 1);

    Xil[0].push_back(parameter[2]*xi0+parameter[4]+parameter[6]/new_rad[0][i]+parameter[8]/(new_rad[0][i]*new_rad[0][i]));

  }

  for(size_t i=0; i<new_rad[1].size(); i++) {

    auto xi_mu_0 = [&] (const double mu)
      {
	double alpha = sqrt(mu*mu*alpha_parallel*alpha_parallel+(1-mu*mu)*alpha_perpendicular*alpha_perpendicular);
	double sp = new_rad[1][i]*alpha;
	double mup = mu*alpha_parallel/alpha;
	double val = 0;
	for (int j=0; j<pp->nmultipoles; j++)
	  val += pp->func_multipoles[j]->operator()(sp)*cbl::legendre_polynomial(mup, j*2);
	return val;
      };

    auto xi_mu_2 = [&] (const double mu)
      {
	double alpha =sqrt(mu*mu*alpha_parallel*alpha_parallel+(1-mu*mu)*alpha_perpendicular*alpha_perpendicular);
	double sp = new_rad[1][i]*alpha;
	double mup = mu*alpha_parallel/alpha;
	double val=0;
	for (int j=0; j<pp->nmultipoles; j++)
	  val += pp->func_multipoles[j]->operator()(sp)*cbl::legendre_polynomial(mup, j*2);
      
	return 3*val*mu*mu;
      };

    double xi0 = cbl::wrapper::gsl::GSL_integrate_qag(xi_mu_0, 0, 1);
    double ximu2 = cbl::wrapper::gsl::GSL_integrate_qag(xi_mu_2, 0, 1);

    Xil[1].push_back(2.5*(parameter[3]*ximu2-parameter[2]*xi0)+parameter[5]+parameter[7]/new_rad[1][i]+parameter[9]/(new_rad[1][i]*new_rad[1][i]));
  }

  vector<double> Xi;

  for (size_t i=0; i<Xil.size(); i++)
    for (size_t j=0; j<Xil[i].size(); j++)
      Xi.push_back(Xil[i][j]);

  return Xi;
}

