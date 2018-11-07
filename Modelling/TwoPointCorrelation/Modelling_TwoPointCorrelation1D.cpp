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
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */


#include "Data1D_extra.h"
#include "Modelling_TwoPointCorrelation1D.h"

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


void cbl::modelling::twopt::Modelling_TwoPointCorrelation1D::set_data_model (const cosmology::Cosmology cosmology, const double redshift, const string method_Pk, const double sigmaNL_perp, const double sigmaNL_par, const bool NL, const double bias, const double pimax, const double r_min, const double r_max, const double k_min, const double k_max, const int step, const string output_dir, const string output_root, const int norm, const double aa, const bool GSL, const double prec, const string file_par, const double Delta, const bool isDelta_vir, const std::vector<double> cluster_redshift, const std::vector<double> cluster_mass_proxy, const std::vector<double> cluster_mass_proxy_error, const string model_bias, const string meanType, const int seed, const cosmology::Cosmology cosmology_mass, const std::vector<double> redshift_source)
{
  m_data_model = make_shared<STR_data_model>(STR_data_model());

  m_data_model->cosmology = make_shared<cosmology::Cosmology>(cosmology);
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
  m_data_model->Delta = (isDelta_vir) ? m_data_model->cosmology->Delta_vir(Delta, redshift) : Delta;
  m_data_model->model_bias = model_bias;
  m_data_model->meanType = meanType;
  m_data_model->gau_ran = make_shared<random::NormalRandomNumbers>(random::NormalRandomNumbers(0., 1., seed));

  if (m_data_model->cosmology->sigma8()>0) {
    m_data_model->sigma8_z = m_data_model->cosmology->sigma8(m_data_model->redshift);
  }
  else { 
    coutCBL << "sigma8 is not set, it will be computed from the power spectrum with " << m_data_model->method_Pk << endl;
    m_data_model->sigma8_z = m_data_model->cosmology->sigma8_Pk(m_data_model->method_Pk, m_data_model->redshift, m_data_model->output_root);
    coutCBL << "--> sigma8(z=" << m_data_model->redshift << ") = " << m_data_model->sigma8_z << endl << endl;
  }
  m_data_model->linear_growth_rate_z = m_data_model->cosmology->linear_growth_rate(m_data_model->redshift, 1.);


  m_data_model->DVfid = m_data_model->cosmology->D_V(m_data_model->redshift);
  m_data_model->DAfid = m_data_model->cosmology->D_A(m_data_model->redshift);
  m_data_model->HHfid = m_data_model->cosmology->HH(m_data_model->redshift);
  
  if (cluster_mass_proxy.size()>0) {
    vector<double> cl_zz(cluster_mass_proxy.size(), m_data_model->redshift), cl_sigma(cluster_mass_proxy.size());
    if (cluster_redshift.size()==cluster_mass_proxy.size())
      cl_zz = cluster_redshift;
    for (size_t i=0; i<cl_sigma.size(); ++i)
      cl_sigma[i] = sqrt(m_data_model->cosmology->sigma2M(cluster_mass_proxy[i], m_data_model->method_Pk, 0., m_data_model->output_root, "Spline", m_data_model->k_max));
    
    m_data_model->cluster_mass_proxy = make_shared<data::Data1D_extra>(data::Data1D_extra(cl_zz, cluster_mass_proxy, cluster_mass_proxy_error, {cl_sigma}));

    m_data_model->cosmology_mass = cosmology_mass;
    m_data_model->redshift_source = redshift_source;
  }
}


// ============================================================================================


void cbl::modelling::twopt::Modelling_TwoPointCorrelation1D::set_data_HOD (const cosmology::Cosmology cosmology, const double redshift, const string model_MF, const string model_bias, const double Mh_min, const double Mh_max, const double pi_max, const double r_max_int, const double r_min, const double r_max, const double k_min, const double k_max, const int step, const string method_Pk, const bool NL, const string output_root, const double Delta, const double kk, const string interpType, const int norm, const double prec, const string input_file, const bool is_parameter_file, const string model_cM, const string profile, const string halo_def)
{
  m_data_HOD.cosmology = make_shared<cosmology::Cosmology>(cosmology);
  m_data_HOD.redshift = redshift;
  m_data_HOD.model_MF = model_MF;
  m_data_HOD.model_bias = model_bias;
  m_data_HOD.Mh_min = Mh_min;
  m_data_HOD.Mh_max = Mh_max;
  m_data_HOD.pi_max = pi_max;
  m_data_HOD.r_max_int = r_max_int;
  m_data_HOD.r_min = r_min;
  m_data_HOD.r_max = r_max;
  m_data_HOD.k_min = k_min;
  m_data_HOD.k_max = k_max;
  m_data_HOD.step = step;
  m_data_HOD.method_Pk = method_Pk;
  m_data_HOD.NL = NL;
  m_data_HOD.output_root = output_root;
  m_data_HOD.Delta = Delta;
  m_data_HOD.kk = kk;
  m_data_HOD.interpType = interpType;
  m_data_HOD.norm = norm;
  m_data_HOD.prec = prec;
  m_data_HOD.input_file = input_file;
  m_data_HOD.is_parameter_file = is_parameter_file;
  m_data_HOD.model_cM = model_cM;
  m_data_HOD.profile = profile;
  m_data_HOD.halo_def = halo_def;
}


// ============================================================================================



void cbl::modelling::twopt::Modelling_TwoPointCorrelation1D::set_data_model_cluster_selection_function (const cosmology::Cosmology cosmology, const cosmology::Cosmology test_cosmology, const double mean_redshift, const string model_MF, const string model_bias, const string selection_function_file, const std::vector<int> selection_function_column, const double z_min, const double z_max, const double Mass_min, const double Mass_max, const double Delta, const bool isDelta_vir, const string method_Pk, const string output_dir, const double k_min, const double k_max, const double prec, const int step, const int mass_step)

{
  m_data_model->cosmology = make_shared<cosmology::Cosmology>(cosmology);
  m_data_model->test_cosmology = make_shared<cosmology::Cosmology>(test_cosmology);
  m_data_model->redshift = mean_redshift;
  m_data_model->method_Pk = method_Pk;
  m_data_model->k_min = k_min;
  m_data_model->k_max = k_max;
  m_data_model->step = step;
  m_data_model->prec = prec;
  m_data_model->output_root = "test";
  m_data_model->Delta = Delta;
  m_data_model->isDelta_Vir = isDelta_vir;
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
    sf_monoz[i] = m_data_model->interp_SelectionFunction->operator()(mass_grid[i], m_data_model->redshift);
  m_data_model->interp_SelectionFunction_cut = make_shared<glob::FuncGrid>(glob::FuncGrid(mass_grid, sf_monoz, "Spline"));

  m_data_model->z_min = (z_min>par::defaultDouble) ? z_min : Min(redshift_grid);
  m_data_model->z_max = (z_max>par::defaultDouble) ? z_max : Max(redshift_grid);
  m_data_model->Mass_min = (Mass_min>par::defaultDouble) ? Mass_min : Min(mass_grid);
  m_data_model->Mass_max = (Mass_max>par::defaultDouble) ? Mass_max : Max(mass_grid);

  m_data_model->DVfid = m_data_model->cosmology->D_V(m_data_model->redshift);
}
