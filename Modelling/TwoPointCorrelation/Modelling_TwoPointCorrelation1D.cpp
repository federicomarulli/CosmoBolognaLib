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
 *  Modelling/TwoPointCorrelation/Modelling_TwoPointCorrelation1D.cpp
 *
 *  @brief Methods of the class Modelling_TwoPointCorrelation1D
 *
 *  This file contains the implementation of the methods of the class
 *  Modelling_TwoPointCorrelation1D, i.e. the common functions to
 *  model the 1D two-point correlation functions of any type
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unibo.it, alfonso.veropalumbo@unibo.it
 */

#include "MassFunction.h"
#include "Bias.h"
#include "Data1D_extra.h"
#include "Modelling_TwoPointCorrelation1D.h"
#include "BAO.h"

using namespace std;

using namespace cbl;


// ============================================================================================


cbl::modelling::twopt::Modelling_TwoPointCorrelation1D::Modelling_TwoPointCorrelation1D (const std::shared_ptr<cbl::measure::twopt::TwoPointCorrelation> twop)
{
  m_data = twop->dataset();
  m_twoPType = twop->twoPType();
}


// ============================================================================================


cbl::modelling::twopt::Modelling_TwoPointCorrelation1D::Modelling_TwoPointCorrelation1D (const std::shared_ptr<cbl::data::Data> dataset, const measure::twopt::TwoPType twoPType)
{
  m_data = dataset;
  m_twoPType = twoPType;
}

// ============================================================================================

void cbl::modelling::twopt::Modelling_TwoPointCorrelation1D::set_data_model (const std::shared_ptr<cosmology::Cosmology> cosmology, const double redshift, const std::string method_Pk, const bool NL, const bool store_output, const double k_min, const double k_max, const std::string output_dir, const std::string output_root, const int norm, const double prec, const std::string file_par)
{
  if (file_par!=par::defaultString)
    WarningMsgCBL("check the consistency between the parameters of the object cosmology, provided in input, and the ones in the parameter file", "set_data_model", "Modelling_TwoPointCorrelation1D.cpp");
  
  m_data_model = make_shared<STR_data_model>(STR_data_model());

  m_data_model->cosmology = move(cosmology);
  m_data_model->redshift = redshift;
  m_data_model->method_Pk = method_Pk;
  m_data_model->NL = NL;
  m_data_model->store_output = store_output;
  m_data_model->k_min = k_min;
  m_data_model->k_max = k_max;
  m_data_model->output_dir = output_dir;
  m_data_model->output_root = output_root;
  m_data_model->norm = norm;
  m_data_model->prec = prec;
  m_data_model->file_par = file_par;
}

// ============================================================================================


void cbl::modelling::twopt::Modelling_TwoPointCorrelation1D::set_data_model (const std::shared_ptr<cosmology::Cosmology> cosmology, const double redshift, const std::string method_Pk, const double sigmaNL_perp, const double sigmaNL_par, const bool NL, const double bias, const double pimax, const double r_min, const double r_max, const double k_min, const double k_max, const int step,  const std::string output_dir, const std::string output_root, const int norm, const double aa, const bool GSL, const double prec, const std::string file_par, const double Delta, const bool isDelta_critical, const std::vector<double> cluster_redshift, const std::vector<double> cluster_mass_proxy, const std::vector<double> cluster_mass_proxy_error, const std::string model_bias, const std::string meanType, const int seed, std::shared_ptr<cosmology::Cosmology> cosmology_mass, const std::vector<double> redshift_source)
{
  if (file_par!=par::defaultString)
    WarningMsgCBL("check the consistency between the parameters of the object cosmology, provided in input, and the ones in the parameter file", "set_data_model", "Modelling_TwoPointCorrelation1D.cpp");
  
  m_data_model = make_shared<STR_data_model>(STR_data_model());

  m_data_model->cosmology = move(cosmology);
  m_data_model->redshift = redshift;
  m_data_model->method_Pk = method_Pk;
  m_data_model->sigmaNL_perp = sigmaNL_perp;
  m_data_model->sigmaNL_par = sigmaNL_par;
  m_data_model->sigmaNL = sqrt(0.5*(sigmaNL_par*sigmaNL_par+sigmaNL_perp*sigmaNL_perp));
  m_data_model->NL = NL;
  m_data_model->pi_max = pimax;
  m_data_model->r_min = r_min;
  m_data_model->r_max = r_max;
  m_data_model->k_min = k_min;
  m_data_model->k_max = k_max;
  m_data_model->step = step;
  m_data_model->output_dir = output_dir;
  m_data_model->output_root = output_root;
  m_data_model->norm = norm;
  m_data_model->aa = aa;
  m_data_model->GSL = GSL;
  m_data_model->prec = prec;
  m_data_model->file_par = file_par;
  m_data_model->bias = bias;
  m_data_model->Delta_input = Delta;
  m_data_model->isDelta_critical = isDelta_critical;
  m_data_model->Delta = (isDelta_critical) ? Delta/cosmology->Omega_matter(redshift) : Delta;
  m_data_model->model_bias = model_bias;
  m_data_model->meanType = meanType;
  m_data_model->gau_ran = make_shared<random::NormalRandomNumbers>(random::NormalRandomNumbers(0., 1., seed));

  cosmology::PkXi PX(cosmology);
  cosmology::Sigma SG(cosmology);
  
  if (cosmology->sigma8()>0) 
    m_data_model->sigma8_z = cosmology->sigma8(redshift);
  
  else { 
    coutCBL << "sigma8 is not set, it will be computed from the power spectrum with " << method_Pk << endl;
    
    m_data_model->sigma8_z = PX.sigma8_Pk(method_Pk, redshift, true, output_root, NL, k_min, k_max, prec, file_par);
    
    coutCBL << "--> sigma8(z=" << redshift << ") = " << m_data_model->sigma8_z << endl << endl;
  }
  m_data_model->linear_growth_rate_z = cosmology->linear_growth_rate(redshift, 1.);

  m_data_model->DVfid = cosmology->D_V(redshift);
  m_data_model->DAfid = cosmology->D_A(redshift);
  m_data_model->HHfid = cosmology->Hubble(redshift);

  vector<double> sigma2M = SG.sigma2M(cluster_mass_proxy, method_Pk, 0., true, output_root, "Spline", k_max);

  if (cluster_mass_proxy.size() > 0)
  {
    vector<double> cl_zz(cluster_mass_proxy.size(), redshift), cl_sigma(cluster_mass_proxy.size());
    if (cluster_redshift.size()==cluster_mass_proxy.size())
      cl_zz = cluster_redshift;
    for (size_t i=0; i<cl_sigma.size(); ++i)
      cl_sigma[i] = sqrt(sigma2M[i]);

    m_data_model->cluster_mass_proxy = make_shared<data::Data1D_extra>(data::Data1D_extra(cl_zz, cluster_mass_proxy, cluster_mass_proxy_error, {cl_sigma}));

    m_data_model->cosmology_mass = cosmology_mass;
    m_data_model->redshift_source = redshift_source;
  }
}


// ============================================================================================


void cbl::modelling::twopt::Modelling_TwoPointCorrelation1D::set_data_model (const std::shared_ptr<cosmology::Cosmology> cosmology, const double redshift, const std::vector<double> cluster_redshift, const std::vector<double> cluster_mass_proxy, const double redshift_pivot, const double proxy_pivot, double mass_pivot, const double log_base, const std::string method_Pk, const bool NL, const double r_min, const double r_max, const double k_min, const double k_max, const int step,  const std::string output_dir, const std::string output_root, const int norm, const double prec, const std::string file_par, const double Delta, const bool isDelta_critical, const std::string model_bias)
{
  if (file_par!=par::defaultString)
    WarningMsgCBL("check the consistency between the parameters of the object cosmology, provided in input, and the ones in the parameter file", "set_data_model", "Modelling_TwoPointCorrelation1D.cpp");
  
  m_data_model = make_shared<STR_data_model>(STR_data_model());

  m_data_model->cosmology = move(cosmology);
  m_data_model->redshift = redshift;
  m_data_model->method_Pk = method_Pk;
  m_data_model->NL = NL;
  m_data_model->r_min = r_min;
  m_data_model->r_max = r_max;
  m_data_model->k_min = k_min;
  m_data_model->k_max = k_max;
  m_data_model->step = step;
  m_data_model->output_dir = output_dir;
  m_data_model->output_root = output_root;
  m_data_model->norm = norm;
  m_data_model->prec = prec;
  m_data_model->file_par = file_par;
  m_data_model->Delta_input = Delta;
  m_data_model->isDelta_critical = isDelta_critical;
  m_data_model->Delta = (isDelta_critical) ? Delta/cosmology->Omega_matter(redshift) : Delta;
  m_data_model->model_bias = model_bias;

  cosmology::PkXi PX(cosmology);
  
  if (cosmology->sigma8()>0) 
    m_data_model->sigma8_z = cosmology->sigma8(redshift);  
  else { 
    coutCBL << "sigma8 is not set, it will be computed from the power spectrum with " << method_Pk << endl;    
    m_data_model->sigma8_z = PX.sigma8_Pk(method_Pk, redshift, false, output_root, NL, k_min, k_max, prec, file_par);    
    coutCBL << "--> sigma8(z=" << redshift << ") = " << m_data_model->sigma8_z << endl << endl;
  }
  
  if (cluster_mass_proxy.size() == 0 || cluster_redshift.size() == 0 || cluster_mass_proxy.size() != cluster_redshift.size())
    ErrorCBL("The mass proxy and redshift vectors must contain at least one object each, and must have the same size.", "set_data_model", "Modelling_TwoPointCorrelation1D.cpp");

  // Build a dummy dataset for the scaling relation Modelling object.
  // Only the x values of the dataset (the mass proxy values) are used in the computation
  std::shared_ptr<cbl::data::Data> dataset = std::make_shared<cbl::data::Data1D>(cbl::data::Data1D(cluster_mass_proxy, cluster_mass_proxy, cluster_mass_proxy));
  
  // Build the scaling relation Modelling object
  modelling::massobsrel::Modelling_MassObservableRelation scaling_relation (dataset);
  m_data_model->scaling_relation = make_shared<modelling::massobsrel::Modelling_MassObservableRelation>(scaling_relation);
  
  (m_data_model->scaling_relation)->set_data_model(cosmology, cluster_redshift, redshift_pivot, proxy_pivot, log_base);
  (m_data_model->scaling_relation)->set_mass_pivot(mass_pivot);
}


//===================================================================================


void cbl::modelling::twopt::Modelling_TwoPointCorrelation1D::set_data_model_bias(const double z_abs_err, const std::shared_ptr<cosmology::Cosmology> cosmology, const std::vector<double> cluster_mass, const std::vector<double> cluster_redshift, const double z_min, const double z_max, const std::string method_Pk, const bool NL, const int norm, const double k_min, const double k_max, const std::string output_dir, const std::string output_root, const bool redshift_space, const string interpType, const double Delta, const bool isDelta_critical, const std::string model_MF, const std::string model_bias, const int mass_step, const bool bias_from_sample)
{
  m_data_model = make_shared<STR_data_model>(STR_data_model());
  m_data_model->z_abs_err = z_abs_err;
  m_data_model->cosmology = move(cosmology);
  m_data_model->z_min=z_min;
  m_data_model->z_max=z_max;
  double redshift=cbl::Average(cluster_redshift); 
  m_data_model->r_min = 10;
  m_data_model->r_max = 200;
  m_data_model->step = 500;
  m_data_model->method_Pk = method_Pk;
  m_data_model->NL = NL;
  m_data_model->norm = norm;
  m_data_model->prec = 1.e-3;
  m_data_model->k_min = k_min;
  m_data_model->k_max = k_max;
  m_data_model->output_dir = output_dir;
  m_data_model->output_root = output_root;
  m_data_model->interpType=interpType;
  m_data_model->redshift_space=redshift_space;
  m_data_model->Delta_input = Delta;
  m_data_model->isDelta_critical=isDelta_critical;
  m_data_model->Delta = (isDelta_critical) ? Delta/cosmology->Omega_matter(redshift) : Delta;
  m_data_model->model_MF=model_MF;
  m_data_model->model_bias=model_bias;
  m_data_model->mass_step=mass_step;
  m_data_model->DVfid = cosmology->D_V(redshift);
  m_data_model->DAfid = cosmology->D_A(redshift);
  m_data_model->HHfid = cosmology->Hubble(redshift);
  cbl::cosmology::BAO bao(cosmology);
  m_data_model->rsfid = bao.rs(method_Pk);
  m_data_model->bias_from_sample = bias_from_sample;

  if(norm!=0) {
    cosmology::PkXi PX(cosmology);
    double s8 = PX.sigma8_Pk(method_Pk, 0., false, "test", NL, k_min, k_max);
    m_data_model->cosmology->set_parameter("sigma8", s8);
  }

  if (cluster_mass.size() == 0 || cluster_redshift.size() == 0 || cluster_mass.size() != cluster_redshift.size())
    ErrorCBL("The mass and redshift vectors must contain at least one object each, and must have the same size.", "set_data_model_bias", "Modelling_TwoPointCorrelation1D.cpp");
  m_data_model->mass=cluster_mass;
  m_data_model->cluster_redshift=cluster_redshift;

  m_data_model->Mass_min=*std::min_element(cluster_mass.begin(), cluster_mass.end());
  m_data_model->Mass_max=*std::max_element(cluster_mass.begin(), cluster_mass.end());
}


// ============================================================================================


void cbl::modelling::twopt::Modelling_TwoPointCorrelation1D::set_data_model (const double z_abs_err, const double proxy_rel_err, const std::shared_ptr<cosmology::Cosmology> cosmology, const double redshift, const std::vector<double> cluster_redshift, const std::vector<double> cluster_mass_proxy, const double redshift_pivot, const double proxy_pivot, double mass_pivot, const double log_base, const std::string method_Pk, const bool NL, const double r_min, const double r_max, const double k_min, const double k_max, const int step,  const std::string output_dir, const std::string output_root, const int norm, const double prec, const std::string file_par, const double Delta, const bool isDelta_critical, const std::string model_bias, const std::string Pk_mu_model)
{
  set_data_model(cosmology, redshift, cluster_redshift, cluster_mass_proxy, redshift_pivot, proxy_pivot, mass_pivot, log_base, method_Pk, NL, r_min, r_max, k_min, k_max, step, output_dir, output_root, norm, prec, file_par, Delta, isDelta_critical, model_bias);

  if (z_abs_err < 0 || proxy_rel_err < 0)
    ErrorCBL("The errors on z and proxy must be > 0.", "set_data_model", "Modelling_TwoPointCorrelation1D.cpp");

  m_data_model->z_abs_err = z_abs_err;
  m_data_model->proxy_rel_err = proxy_rel_err;
  m_data_model->Pk_mu_model = Pk_mu_model;

}


// ============================================================================================


void cbl::modelling::twopt::Modelling_TwoPointCorrelation1D::set_data_HOD (const std::shared_ptr<cosmology::Cosmology> cosmology, const double redshift, const string model_MF, const string model_bias, const double Mh_min, const double Mh_max, const double pi_max, const double r_max_int, const double r_min, const double r_max, const double k_min, const double k_max, const int step, const double m_min, const double m_max, const int m_step, const string method_Pk, const bool NL, const string output_root, const double Delta, const double kk, const string interpType, const int norm, const double prec, const string input_file, const bool is_parameter_file, const string model_cM, const string profile, const string halo_def)
{
  m_data_HOD = make_shared<STR_data_HOD>(STR_data_HOD());
  
  m_data_HOD->cosmology = move(cosmology);
  m_data_HOD->redshift = redshift;
  m_data_HOD->model_MF = model_MF;
  m_data_HOD->model_bias = model_bias;
  m_data_HOD->Mh_min = Mh_min;
  m_data_HOD->Mh_max = Mh_max;
  m_data_HOD->pi_max = pi_max;
  m_data_HOD->r_max_int = r_max_int;
  m_data_HOD->r_min = r_min;
  m_data_HOD->r_max = r_max;
  
  //  creation of the mass vector for interpolating pk 
  m_data_HOD->k_min = k_min;
  m_data_HOD->k_max = k_max;
  m_data_HOD->step = step;
  m_data_HOD->kkvec = logarithmic_bin_vector(step, k_min, k_max);
  
  //  creation of the mass vector for interpolating mass functiom and bias
  m_data_HOD->m_min = m_min;
  m_data_HOD->m_max = m_max;
  m_data_HOD->m_step = m_step;
  m_data_HOD->massvec = logarithmic_bin_vector(m_step, m_min, m_max);

  //  creation of vectors containing mass function and bias values for interpolation
  vector<double> mass_function_vec;
  vector<double> bias_vec;

  cosmology::MassFunction mass_function(m_data_HOD->cosmology);
  cosmology::Bias bias(m_data_HOD->cosmology);
  
  for (int i=0; i<m_step; i++) {
    mass_function_vec.emplace_back(mass_function.mass_function(m_data_HOD->massvec[i], redshift, model_MF, "CAMB"));
    bias_vec.emplace_back(bias.bias_halo(m_data_HOD->massvec[i], redshift, model_MF, "CAMB"));
  }  
  m_data_HOD->mass_function_vec = mass_function_vec;
  m_data_HOD->bias_vec = bias_vec;

  //  creation of vector containing pk values for interpolation
  cosmology::PkXi PX(m_data_HOD->cosmology);
  vector<double> pk_vec = PX.Pk_matter(m_data_HOD->kkvec, method_Pk, NL, redshift, false, output_root, norm, k_min, k_max, prec, input_file); 
  
  m_data_HOD->pk_vec = pk_vec;

  //  creation of FuncGrid objects for interpolation
  glob::FuncGrid MF_grid(m_data_HOD->massvec, mass_function_vec, "Spline");
  m_data_HOD->interpMF = MF_grid;
  glob::FuncGrid BI_grid(m_data_HOD->massvec, bias_vec, "Spline");
  m_data_HOD->interpBias = BI_grid;
  glob::FuncGrid PK_grid(m_data_HOD->kkvec, pk_vec, "Spline");
  m_data_HOD->interpPk = PK_grid;

  m_data_HOD->method_Pk = method_Pk;
  m_data_HOD->NL = NL;
  m_data_HOD->output_root = output_root;
  m_data_HOD->Delta = Delta;
  m_data_HOD->kk = kk;
  m_data_HOD->interpType = interpType;
  m_data_HOD->norm = norm;
  m_data_HOD->prec = prec;
  m_data_HOD->input_file = input_file;
  m_data_HOD->is_parameter_file = is_parameter_file;
  m_data_HOD->model_cM = model_cM;
  m_data_HOD->profile = profile;
  m_data_HOD->halo_def = halo_def;
}


// ============================================================================================


void cbl::modelling::twopt::Modelling_TwoPointCorrelation1D::set_data_model_cluster_selection_function (const std::shared_ptr<cosmology::Cosmology> cosmology, const std::shared_ptr<cosmology::Cosmology> test_cosmology, const double mean_redshift, const string model_MF, const string model_bias, const string selection_function_file, const std::vector<int> selection_function_column, const double z_min, const double z_max, const double Mass_min, const double Mass_max, const string file_par, const double Delta, const bool isDelta_critical, const string method_Pk, const string output_dir, const double k_min, const double k_max, const double prec, const int step, const int mass_step)
{
  if (file_par!=par::defaultString)
    WarningMsgCBL("check the consistency between the parameters of the object cosmology, provided in input, and the ones in the parameter file", "set_data_model_cluster_selection_function", "Modelling_TwoPointCorrelation1D.cpp");
  
  m_data_model = make_shared<STR_data_model>(STR_data_model());
  
  m_data_model->cosmology = move(cosmology);
  m_data_model->test_cosmology = move(test_cosmology);
  m_data_model->redshift = mean_redshift;
  m_data_model->method_Pk = method_Pk;
  m_data_model->k_min = k_min;
  m_data_model->k_max = k_max;
  m_data_model->step = step;
  m_data_model->prec = prec;
  m_data_model->output_root = "test";
  m_data_model->Delta = Delta;
  m_data_model->isDelta_critical = isDelta_critical;
  m_data_model->model_MF = model_MF;
  m_data_model->model_bias = model_bias;
  m_data_model->kk = logarithmic_bin_vector(step, k_min, k_max);
  m_data_model->mass = logarithmic_bin_vector(mass_step, Mass_min, Mass_max);
  m_data_model->file_par = par::defaultString;
  m_data_model->output_dir = output_dir;

  vector<double> mass_grid, redshift_grid;
  vector<vector<double>> selection_function;
  read_matrix(selection_function_file, mass_grid, redshift_grid, selection_function, selection_function_column);
  m_data_model->interp_SelectionFunction = make_shared<glob::FuncGrid2D>(glob::FuncGrid2D(mass_grid, redshift_grid, selection_function, "Linear"));

  vector<double> sf_monoz(mass_grid.size());
  for (size_t i=0; i<mass_grid.size(); i++)
    sf_monoz[i] = m_data_model->interp_SelectionFunction->operator()(mass_grid[i], mean_redshift);
  m_data_model->interp_SelectionFunction_cut = make_shared<glob::FuncGrid>(glob::FuncGrid(mass_grid, sf_monoz, "Spline"));

  m_data_model->z_min = (z_min>par::defaultDouble) ? z_min : Min(redshift_grid);
  m_data_model->z_max = (z_max>par::defaultDouble) ? z_max : Max(redshift_grid);
  m_data_model->Mass_min = (Mass_min>par::defaultDouble) ? Mass_min : Min(mass_grid);
  m_data_model->Mass_max = (Mass_max>par::defaultDouble) ? Mass_max : Max(mass_grid);

  m_data_model->DVfid = cosmology->D_V(mean_redshift);
}
