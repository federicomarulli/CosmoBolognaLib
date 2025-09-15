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
 *  Modelling/TwoPointCorrelation/Modelling_TwoPointCorrelation_multipoles.cpp
 *
 *  @brief Methods of the class
 *  Modelling_TwoPointCorrelation_multipoles
 *
 *  This file contains the implementation of the methods of the class
 *  Modelling_TwoPointCorrelation_multipoles, used to model the
 *  multipoles of the two-point correlation function
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo, Massimiliano Romanello
 *
 *  @authors federico.marulli3@unibo.it, alfonso.veropalumbo@unibo.it, massimilia.romanell2@unibo.it
 */


#include "Data1D.h"
#include "PkXiNonLinear.h"
#include "Modelling_TwoPointCorrelation_multipoles.h"

using namespace std;

using namespace cbl;


// ============================================================================================


cbl::modelling::twopt::Modelling_TwoPointCorrelation_multipoles::Modelling_TwoPointCorrelation_multipoles (const std::shared_ptr<cbl::measure::twopt::TwoPointCorrelation> twop)
  : Modelling_TwoPointCorrelation1D_monopole(twop), m_nmultipoles(3)
{
  m_ModelIsSet = false;

  m_multipoles_order.erase(m_multipoles_order.begin(), m_multipoles_order.end());

  int size = m_data->ndata()/m_nmultipoles;

  for (int j=0; j<m_nmultipoles; j++)
    for (int i=0; i<size; i++)
      m_multipoles_order.push_back(j);

  m_use_pole.resize(3, true);
}


// ============================================================================================


cbl::modelling::twopt::Modelling_TwoPointCorrelation_multipoles::Modelling_TwoPointCorrelation_multipoles (const std::shared_ptr<data::Data> twop_dataset, const int nmultipoles)
  : Modelling_TwoPointCorrelation1D_monopole(twop_dataset), m_nmultipoles(nmultipoles)
{
  m_ModelIsSet = false;

  m_multipoles_order.erase(m_multipoles_order.begin(), m_multipoles_order.end());
  m_use_pole.resize(3, false);

  for (int j=0; j<m_nmultipoles; j++) {
    m_use_pole[j] = true;
    for (int i=0; i<m_data->ndata()/m_nmultipoles; i++)
      m_multipoles_order.push_back(j);
  }

}

// ============================================================================================

cbl::modelling::twopt::Modelling_TwoPointCorrelation_multipoles::Modelling_TwoPointCorrelation_multipoles (const std::shared_ptr<data::Data> twop_dataset, const std::vector<bool> use_pole)
  : Modelling_TwoPointCorrelation1D_monopole(twop_dataset) 
{
  m_ModelIsSet = false;

  m_multipoles_order.erase(m_multipoles_order.begin(), m_multipoles_order.end());

  m_use_pole=use_pole;
  m_nmultipoles=use_pole.size(); 
  
  for (int j=0; j<m_nmultipoles; j++) 
    for (int i=0; i<m_data->ndata()/m_nmultipoles; i++)
      m_multipoles_order.push_back(j);
    
}


// ============================================================================================


cbl::modelling::twopt::Modelling_TwoPointCorrelation_multipoles::Modelling_TwoPointCorrelation_multipoles (const std::shared_ptr<data::Data> twop_dataset, const std::vector<int> dataset_order)
  : Modelling_TwoPointCorrelation1D_monopole(twop_dataset) 
{
  m_ModelIsSet = false;

  m_multipoles_order = dataset_order;
  for(size_t i = 0; i < dataset_order.size(); i++)
    m_use_pole.emplace_back(true);
  
  m_nmultipoles=dataset_order.size(); 
      
}


// ============================================================================================


void cbl::modelling::twopt::Modelling_TwoPointCorrelation_multipoles::write_model (const std::string output_dir, const std::string output_file, const int nmultipoles, const std::vector<double> xx, const std::vector<double> parameters)
{
  if (m_likelihood==NULL) ErrorCBL("this function requires the likelihood to be defined (with the function set_likelihood)!", "write_model", "Modelling_TwoPointCorrelation_multipoles.cpp");

  int nmultipoles_original = m_data_model->nmultipoles;
  vector<int> dataset_order_original = m_data_model->dataset_order;

  m_data_model->nmultipoles = nmultipoles;

  vector<bool> new_use_pole(3, false);
  vector<int> new_dataset_order;
  vector<double> new_xx;

  if (xx.size()==0)
    new_xx = m_data_fit->xx();
  else
    for (int n=0; n<m_data_model->nmultipoles; n++) {
      new_use_pole[n] = true;
      for (size_t i=0; i<xx.size(); i++) {
	new_xx.push_back(xx[i]);
	new_dataset_order.push_back(n);
      }
    }

  m_data_model->dataset_order = new_dataset_order;
  m_data_model->use_pole = new_use_pole;

  m_likelihood->write_model(output_dir, output_file, parameters, new_xx);

  m_data_model->dataset_order = dataset_order_original;
  m_data_model->nmultipoles = nmultipoles_original;
  m_data_model->use_pole = m_use_pole;

}


// ============================================================================================


void cbl::modelling::twopt::Modelling_TwoPointCorrelation_multipoles::write_model (const std::string output_dir, const std::string output_file, const std::vector<bool> use_pole, const std::vector<double> xx, const std::vector<double> parameters)
{
  if (m_likelihood==NULL) ErrorCBL("this function requires the likelihood to be defined (with the function set_likelihood)!", "write_model", "Modelling_TwoPointCorrelation_multipoles.cpp");

  int nmultipoles_original = m_data_model->nmultipoles;
  vector<int> dataset_order_original = m_data_model->dataset_order;

  m_data_model->nmultipoles = use_pole.size();

  vector<int> new_dataset_order;
  vector<double> new_xx;

  if (xx.size()==0)
    new_xx = m_data_fit->xx();
  else
    for (int n=0; n<m_data_model->nmultipoles; n++) {
      if(use_pole[n])
	for (size_t i=0; i<xx.size(); i++) {
	  new_xx.push_back(xx[i]);
	  new_dataset_order.push_back(n);
	}
    }
  
  m_data_model->dataset_order = new_dataset_order;
  
  m_likelihood->write_model(output_dir, output_file, parameters, new_xx);
  
  m_data_model->dataset_order = dataset_order_original;
  m_data_model->nmultipoles = nmultipoles_original;
  m_data_model->use_pole = m_use_pole;
  
}


// ============================================================================================


void cbl::modelling::twopt::Modelling_TwoPointCorrelation_multipoles::write_model_at_bestfit (const std::string output_dir, const std::string output_file, const int nmultipoles, const std::vector<double> xx)
{
  if (m_posterior==NULL)
    ErrorCBL("no posterior found: run maximize_posterior() first!", "write_model_at_best_fit", "Modelling_TwoPointCorrelation_multipoles.cpp");

  int nmultipoles_original = m_data_model->nmultipoles;
  vector<int> dataset_order_original = m_data_model->dataset_order;

  m_data_model->nmultipoles=nmultipoles;

  vector<bool> new_use_pole(3, false);
  vector<int> new_dataset_order;
  vector<double> new_xx;

  if (xx.size()==0)
    new_xx = m_data_fit->xx();
  else
    for (int n=0; n<m_data_model->nmultipoles; n++) {
      new_use_pole[n] = true;
      for (size_t i=0; i<xx.size(); i++) {
	new_xx.push_back(xx[i]);
	new_dataset_order.push_back(n);
      }
    }

  m_data_model->dataset_order = new_dataset_order;
  m_data_model->use_pole = new_use_pole;
  m_posterior->write_model_at_bestfit(output_dir, output_file, new_xx);

  m_data_model->dataset_order = dataset_order_original;
  m_data_model->nmultipoles = nmultipoles_original;
  m_data_model->use_pole = m_use_pole;
}


// ============================================================================================


void cbl::modelling::twopt::Modelling_TwoPointCorrelation_multipoles::write_model_at_bestfit (const std::string output_dir, const std::string output_file, const std::vector<bool> use_pole, const std::vector<double> xx)
{
  if (m_posterior==NULL)
    ErrorCBL("no posterior found: run maximize_posterior() first!", "write_model_at_best_fit", "Modelling_TwoPointCorrelation_multipoles.cpp");

  int nmultipoles_original = m_data_model->nmultipoles;
  vector<int> dataset_order_original = m_data_model->dataset_order;

  m_data_model->nmultipoles = use_pole.size();

  vector<int> new_dataset_order;
  vector<double> new_xx;

  if (xx.size()==0)
    new_xx = m_data_fit->xx();
  else
    for (int n=0; n<m_data_model->nmultipoles; n++) {
      if(use_pole[n])
	for (size_t i=0; i<xx.size(); i++) {
	  new_xx.push_back(xx[i]);
	  new_dataset_order.push_back(n);
	}
    }

  m_data_model->dataset_order = new_dataset_order;
  m_posterior->write_model_at_bestfit(output_dir, output_file, new_xx);

  m_data_model->dataset_order = dataset_order_original;
  m_data_model->nmultipoles = nmultipoles_original;
  m_data_model->use_pole = m_use_pole;
}


// ============================================================================================


void cbl::modelling::twopt::Modelling_TwoPointCorrelation_multipoles::write_model_from_chains (const std::string output_dir, const std::string output_file, const int nmultipoles, const std::vector<double> xx, const int start, const int thin)
{
  if (m_posterior==NULL)
    ErrorCBL("no posterior found: run sample_posterior() first!", "write_model_from_chains", "Modelling_TwoPointCorrelation_multipoles.cpp");

  int nmultipoles_original = m_data_model->nmultipoles;
  vector<int> dataset_order_original = m_data_model->dataset_order;

  m_data_model->nmultipoles=nmultipoles;

  vector<bool> new_use_pole(3, false);
  vector<int> new_dataset_order;
  vector<double> new_xx;

  if (xx.size()==0)
    new_xx = m_data_fit->xx();
  else
    for (int n=0; n<m_data_model->nmultipoles; n++) {
      new_use_pole[n] = true;
      for (size_t i=0; i<xx.size(); i++) {
	new_xx.push_back(xx[i]);
	new_dataset_order.push_back(n);
      }
    }

  m_data_model->dataset_order = new_dataset_order;
  m_data_model->use_pole = new_use_pole;
  m_posterior->write_model_from_chain(output_dir, output_file, new_xx, {}, start, thin);

  m_data_model->dataset_order = dataset_order_original;
  m_data_model->nmultipoles = nmultipoles_original;
  m_data_model->use_pole = m_use_pole;
}


// ============================================================================================


void cbl::modelling::twopt::Modelling_TwoPointCorrelation_multipoles::write_model_from_chains (const std::string output_dir, const std::string output_file, const std::vector<bool> use_pole, const std::vector<double> xx, const int start, const int thin)
{
  if (m_posterior==NULL)
    ErrorCBL("no posterior found: run sample_posterior() first!", "write_model_from_chains", "Modelling_TwoPointCorrelation_multipoles.cpp");

  int nmultipoles_original = m_data_model->nmultipoles;
  vector<int> dataset_order_original = m_data_model->dataset_order;

  m_data_model->nmultipoles=use_pole.size();

  vector<int> new_dataset_order;
  vector<double> new_xx;

  if (xx.size()==0)
    new_xx = m_data_fit->xx();
  else
    for (int n=0; n<m_data_model->nmultipoles; n++) {
      if(use_pole[n])
	for (size_t i=0; i<xx.size(); i++) {
	  new_xx.push_back(xx[i]);
	  new_dataset_order.push_back(n);
	}
    }

  m_data_model->dataset_order = new_dataset_order;
  m_posterior->write_model_from_chain(output_dir, output_file, new_xx, {}, start, thin);

  m_data_model->dataset_order = dataset_order_original;
  m_data_model->nmultipoles = nmultipoles_original;
  m_data_model->use_pole = m_use_pole;
}


// ============================================================================================


void cbl::modelling::twopt::Modelling_TwoPointCorrelation_multipoles::set_fit_range (const double xmin, const double xmax, const int nmultipoles)
{
  if (m_likelihood != NULL)
    ErrorCBL("You must set the fit range before the declaration of the likelihood function!", "set_fit_range", "Modelling_TwoPointCorrelation_multipoles.cpp");
  
  vector<vector<double>> fr(m_nmultipoles, vector<double>(2, -1.));

  int mp = (0<nmultipoles && nmultipoles<m_nmultipoles) ? nmultipoles : m_nmultipoles;

  for (int i=0; i<mp; i++) {
    fr[i][0] = xmin;
    fr[i][1] = xmax;
  }

  set_fit_range(fr);
}


// ============================================================================================


void cbl::modelling::twopt::Modelling_TwoPointCorrelation_multipoles::set_fit_range (const std::vector<std::vector<double>> fit_range)
{
  if (m_likelihood != NULL)
    ErrorCBL("You must set the fit range before the declaration of the likelihood function!", "set_fit_range", "Modelling_TwoPointCorrelation_multipoles.cpp");
  
  if ((int)fit_range.size()!=m_nmultipoles)
    ErrorCBL("the dimension input matrix must be equal to the number of multipoles to be fitted, i.e."+conv(m_nmultipoles, par::fINT)+"!", "set_fit_range", "Modelling_TwoPointCorrelation_multipoles.cpp");

  m_use_pole.resize(3);
  m_use_pole[0] = false;
  m_use_pole[1] = false;
  m_use_pole[2] = false;

  m_multipoles_order.erase(m_multipoles_order.begin(), m_multipoles_order.end());

  const int size = m_data->ndata()/m_nmultipoles;
  vector<bool> mask(m_data->ndata(), false);
  vector<double> xx;

  for (int j=0; j<m_nmultipoles; j++) {
    for (int i=0; i<size; i++) {
      if (fit_range[j][0]<m_data->xx(i+j*size) && m_data->xx(i+j*size)<fit_range[j][1]) {
        m_multipoles_order.push_back(j);
        xx.push_back(m_data->xx(i+j*size));
        m_use_pole[j] = true;
        mask[i+j*size] = true;
      }
    }
  }

  m_data_fit = m_data->cut(mask);

  m_fit_range = true;

  if (m_ModelIsSet)
    m_data_model->dataset_order = m_multipoles_order;
}


// ============================================================================================


void cbl::modelling::twopt::Modelling_TwoPointCorrelation_multipoles::set_fit_range (const std::vector<std::vector<double>> fit_range, const std::vector<int> dataset_order)
{
  if (m_likelihood != NULL)
    ErrorCBL("You must set the fit range before the declaration of the likelihood function!", "set_fit_range", "Modelling_TwoPointCorrelation_multipoles.cpp");

  // Find the number of different multipoles, with occurrence
  std::vector<int> unique_elements;
  std::vector<int> counts;
  
  // Iterate over each element in dataset_order
  for (int elem : dataset_order) {
    // Find the element in unique_elements
    auto it = std::find(unique_elements.begin(), unique_elements.end(), elem);
    
    if (it != unique_elements.end()) {
      // If element is found, increment the corresponding count
      int index = std::distance(unique_elements.begin(), it);
      counts[index]++;
    } else {
      // If element is not found, add it to unique_elements and start its count at 1
      unique_elements.push_back(elem);
      counts.push_back(1);
    }
  }

  m_nmultipoles = unique_elements.size();
  
  if ((int)fit_range.size()!=m_nmultipoles)
    ErrorCBL("the dimension input matrix must be equal to the number of multipoles to be fitted, i.e."+conv(m_nmultipoles, par::fINT)+"!", "set_fit_range", "Modelling_TwoPointCorrelation_multipoles.cpp");

  vector<bool> mask(m_data->ndata(), false);
  vector<double> xx;
  m_use_pole.resize(3);
  m_use_pole[0] = false;
  m_use_pole[1] = false;
  m_use_pole[2] = false;

  for(int j=0; j<m_nmultipoles; j++){
    m_use_pole[j] = true;
    for (int i=0; i<counts[j]; i++) {
      if (fit_range[j][0]<m_data->xx(i+j*counts[j]) && m_data->xx(i+j*counts[j])<fit_range[j][1]) {
	xx.push_back(m_data->xx(i+j*counts[j]));
        mask[i+j*counts[j]] = true;
	
      }
    }
  }
      
  m_data_fit = m_data->cut(mask);

  m_fit_range = true;

  if (m_ModelIsSet)
    m_data_model->dataset_order = m_multipoles_order;
}


// ============================================================================================


void cbl::modelling::twopt::Modelling_TwoPointCorrelation_multipoles::set_fiducial_PkDM ()
{
  m_data_model->nmultipoles = m_nmultipoles;

  m_data_model->kk = logarithmic_bin_vector(m_data_model->step, max(m_data_model->k_min, 1.e-4), min(m_data_model->k_max, 500.));

  cosmology::PkXi PX(m_data_model->cosmology);
  
  vector<double> Pk = PX.Pk_matter(m_data_model->kk, m_data_model->method_Pk, false, m_data_model->redshift, m_data_model->store_output, m_data_model->output_root, m_data_model->norm, m_data_model->k_min, m_data_model->k_max, m_data_model->prec, m_data_model->file_par);

  if (m_data_model->Pk_mu_model=="dispersion_dewiggled") {
    vector<double> PkNW = PX.Pk_matter(m_data_model->kk, "EisensteinHu", false, m_data_model->redshift, m_data_model->store_output, m_data_model->output_root, m_data_model->norm, m_data_model->k_min, m_data_model->k_max, m_data_model->prec, m_data_model->file_par);
    m_data_model->func_Pk = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(m_data_model->kk, Pk, "Spline"));
    m_data_model->func_Pk_NW = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(m_data_model->kk, PkNW, "Spline"));

    m_data_model->funcs_pk.push_back(m_data_model->func_Pk);
    m_data_model->funcs_pk.push_back(m_data_model->func_Pk_NW);
  }

  else if (m_data_model->Pk_mu_model=="dispersion_modecoupling") {
    vector<double> kk_1loop, Pk_1loop;
    cbl::cosmology::PkXiNonLinear PXN(m_data_model->cosmology);
    for (size_t i=0; i<(size_t)m_data_model->step; i++) {
      if (m_data_model->kk[i] < par::pi) {
	kk_1loop.push_back(m_data_model->kk[i]);
	Pk_1loop.push_back(PXN.Pk_1loop(m_data_model->kk[i], m_data_model->func_Pk, 0, m_data_model->k_min, 5., m_data_model->prec));
      }
    }
    m_data_model->func_Pk = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(m_data_model->kk, Pk, "Spline"));
    m_data_model->func_Pk1loop = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(kk_1loop, Pk_1loop, "Spline"));

    m_data_model->funcs_pk.push_back(m_data_model->func_Pk);
    m_data_model->funcs_pk.push_back(m_data_model->func_Pk1loop);
  }

  else if (m_data_model->Pk_mu_model=="dispersion_Gauss" || m_data_model->Pk_mu_model=="dispersion_Lorentz") {
    m_data_model->func_Pk = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(m_data_model->kk, Pk, "Spline"));
    m_data_model->funcs_pk.push_back(m_data_model->func_Pk);
  }

  else if (m_data_model->Pk_mu_model=="Scoccimarro_Pezzotta_Gauss" || m_data_model->Pk_mu_model=="Scoccimarro_Pezzotta_Lorentz" || m_data_model->Pk_mu_model=="Scoccimarro_Bel_Gauss" || m_data_model->Pk_mu_model=="Scoccimarro_Bel_Lorentz") {
    vector<double> Pknonlin = PX.Pk_matter(m_data_model->kk, m_data_model->method_Pk, true, m_data_model->redshift, m_data_model->store_output, m_data_model->output_root, m_data_model->norm, m_data_model->k_min, m_data_model->k_max, m_data_model->prec, m_data_model->file_par);
    m_data_model->func_Pk = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(m_data_model->kk, Pk, "Spline"));
    m_data_model->func_Pk_nonlin = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(m_data_model->kk, Pknonlin, "Spline"));

    m_data_model->funcs_pk.push_back(m_data_model->func_Pk);
    m_data_model->funcs_pk.push_back(m_data_model->func_Pk_nonlin);
  }

  else if (m_data_model->Pk_mu_model=="Scoccimarro_Gauss" || m_data_model->Pk_mu_model=="Scoccimarro_Lorentz") {
    cbl::cosmology::PkXiNonLinear PXN(m_data_model->cosmology);
    vector<vector<double>> Pk_terms = PXN.Pk_TNS_dd_dt_tt(m_data_model->kk, m_data_model->method_Pk, m_data_model->redshift, false, m_data_model->output_root, m_data_model->norm, m_data_model->k_min, m_data_model->k_max, m_data_model->prec);

    m_data_model->func_Pk_DeltaDelta = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(m_data_model->kk, Pk_terms[0], "Spline"));
    m_data_model->func_Pk_DeltaTheta = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(m_data_model->kk, Pk_terms[1], "Spline"));
    m_data_model->func_Pk_ThetaTheta = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(m_data_model->kk, Pk_terms[2], "Spline"));

    m_data_model->funcs_pk.push_back(m_data_model->func_Pk_DeltaDelta);
    m_data_model->funcs_pk.push_back(m_data_model->func_Pk_DeltaTheta);
    m_data_model->funcs_pk.push_back(m_data_model->func_Pk_ThetaTheta);
  }

  else if (m_data_model->Pk_mu_model=="TNS_Gauss" || m_data_model->Pk_mu_model=="TNS_Lorentz") {
    cbl::cosmology::PkXiNonLinear PXN(m_data_model->cosmology);
    vector<vector<double>> Pk_terms = PXN.Pk_TNS_dd_dt_tt(m_data_model->kk, m_data_model->method_Pk, m_data_model->redshift, false, m_data_model->output_root, m_data_model->norm, m_data_model->k_min, m_data_model->k_max, m_data_model->prec);
    vector<vector<double>> Pk_AB = PXN.Pk_TNS_AB_terms_1loop(m_data_model->kk, m_data_model->method_Pk, m_data_model->redshift, false, m_data_model->output_root, m_data_model->norm, m_data_model->k_min, m_data_model->k_max, m_data_model->prec);

    m_data_model->func_Pk_DeltaDelta = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(m_data_model->kk, Pk_terms[0], "Spline"));
    m_data_model->func_Pk_DeltaTheta = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(m_data_model->kk, Pk_terms[1], "Spline"));
    m_data_model->func_Pk_ThetaTheta = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(m_data_model->kk, Pk_terms[2], "Spline"));

    m_data_model->func_Pk_A11 = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(m_data_model->kk, Pk_AB[0], "Spline"));
    m_data_model->func_Pk_A12 = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(m_data_model->kk, Pk_AB[1], "Spline"));
    m_data_model->func_Pk_A22 = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(m_data_model->kk, Pk_AB[2], "Spline"));
    m_data_model->func_Pk_A23 = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(m_data_model->kk, Pk_AB[3], "Spline"));
    m_data_model->func_Pk_A33 = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(m_data_model->kk, Pk_AB[4], "Spline"));
    m_data_model->func_Pk_B12 = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(m_data_model->kk, Pk_AB[5], "Spline"));
    m_data_model->func_Pk_B13 = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(m_data_model->kk, Pk_AB[6], "Spline"));
    m_data_model->func_Pk_B14 = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(m_data_model->kk, Pk_AB[7], "Spline"));
    m_data_model->func_Pk_B22 = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(m_data_model->kk, Pk_AB[8], "Spline"));
    m_data_model->func_Pk_B23 = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(m_data_model->kk, Pk_AB[9], "Spline"));
    m_data_model->func_Pk_B24 = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(m_data_model->kk, Pk_AB[10], "Spline"));
    m_data_model->func_Pk_B33 = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(m_data_model->kk, Pk_AB[11], "Spline"));
    m_data_model->func_Pk_B34 = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(m_data_model->kk, Pk_AB[12], "Spline"));
    m_data_model->func_Pk_B44 = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(m_data_model->kk, Pk_AB[13], "Spline"));

    m_data_model->funcs_pk.push_back(m_data_model->func_Pk_DeltaDelta);
    m_data_model->funcs_pk.push_back(m_data_model->func_Pk_DeltaTheta);
    m_data_model->funcs_pk.push_back(m_data_model->func_Pk_ThetaTheta);
    m_data_model->funcs_pk.push_back(m_data_model->func_Pk_A11);
    m_data_model->funcs_pk.push_back(m_data_model->func_Pk_A12);
    m_data_model->funcs_pk.push_back(m_data_model->func_Pk_A22);
    m_data_model->funcs_pk.push_back(m_data_model->func_Pk_A23);
    m_data_model->funcs_pk.push_back(m_data_model->func_Pk_A33);
    m_data_model->funcs_pk.push_back(m_data_model->func_Pk_B12);
    m_data_model->funcs_pk.push_back(m_data_model->func_Pk_B13);
    m_data_model->funcs_pk.push_back(m_data_model->func_Pk_B14);
    m_data_model->funcs_pk.push_back(m_data_model->func_Pk_B22);
    m_data_model->funcs_pk.push_back(m_data_model->func_Pk_B23);
    m_data_model->funcs_pk.push_back(m_data_model->func_Pk_B24);
    m_data_model->funcs_pk.push_back(m_data_model->func_Pk_B33);
    m_data_model->funcs_pk.push_back(m_data_model->func_Pk_B34);
    m_data_model->funcs_pk.push_back(m_data_model->func_Pk_B44);
  }

  else if (m_data_model->Pk_mu_model=="eTNS_Gauss" || m_data_model->Pk_mu_model=="eTNS_Lorentz") {
    cbl::cosmology::PkXiNonLinear PXN(m_data_model->cosmology);
    vector<vector<double>> Pk_AB = PXN.Pk_TNS_AB_terms_1loop(m_data_model->kk, m_data_model->method_Pk, m_data_model->redshift, false, m_data_model->output_root, m_data_model->norm, m_data_model->k_min, m_data_model->k_max, m_data_model->prec);
    vector<vector<double>> Pk_eTNS_terms = PXN.Pk_eTNS_terms_1loop(m_data_model->kk, m_data_model->method_Pk, m_data_model->redshift, false, m_data_model->output_root, m_data_model->norm, m_data_model->k_min, m_data_model->k_max, m_data_model->prec);

    m_data_model->func_Pk_DeltaDelta = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(m_data_model->kk, Pk_eTNS_terms[0], "Spline"));
    m_data_model->func_Pk_DeltaTheta = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(m_data_model->kk, Pk_eTNS_terms[1], "Spline"));
    m_data_model->func_Pk_ThetaTheta = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(m_data_model->kk, Pk_eTNS_terms[2], "Spline"));

    m_data_model->func_Pk_A11 = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(m_data_model->kk, Pk_AB[0], "Spline"));
    m_data_model->func_Pk_A12 = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(m_data_model->kk, Pk_AB[1], "Spline"));
    m_data_model->func_Pk_A22 = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(m_data_model->kk, Pk_AB[2], "Spline"));
    m_data_model->func_Pk_A23 = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(m_data_model->kk, Pk_AB[3], "Spline"));
    m_data_model->func_Pk_A33 = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(m_data_model->kk, Pk_AB[4], "Spline"));
    m_data_model->func_Pk_B12 = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(m_data_model->kk, Pk_AB[5], "Spline"));
    m_data_model->func_Pk_B13 = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(m_data_model->kk, Pk_AB[6], "Spline"));
    m_data_model->func_Pk_B14 = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(m_data_model->kk, Pk_AB[7], "Spline"));
    m_data_model->func_Pk_B22 = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(m_data_model->kk, Pk_AB[8], "Spline"));
    m_data_model->func_Pk_B23 = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(m_data_model->kk, Pk_AB[9], "Spline"));
    m_data_model->func_Pk_B24 = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(m_data_model->kk, Pk_AB[10], "Spline"));
    m_data_model->func_Pk_B33 = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(m_data_model->kk, Pk_AB[11], "Spline"));
    m_data_model->func_Pk_B34 = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(m_data_model->kk, Pk_AB[12], "Spline"));
    m_data_model->func_Pk_B44 = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(m_data_model->kk, Pk_AB[13], "Spline"));

    m_data_model->func_Pk_b2d = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(m_data_model->kk, Pk_eTNS_terms[3], "Spline"));
    m_data_model->func_Pk_b2v = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(m_data_model->kk, Pk_eTNS_terms[4], "Spline"));
    m_data_model->func_Pk_b22 = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(m_data_model->kk, Pk_eTNS_terms[5], "Spline"));
    m_data_model->func_Pk_bs2d = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(m_data_model->kk, Pk_eTNS_terms[6], "Spline"));
    m_data_model->func_Pk_bs2v = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(m_data_model->kk, Pk_eTNS_terms[7], "Spline"));
    m_data_model->func_Pk_b2s2 = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(m_data_model->kk, Pk_eTNS_terms[8], "Spline"));
    m_data_model->func_Pk_bs22 = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(m_data_model->kk, Pk_eTNS_terms[9], "Spline"));
    m_data_model->func_sigma32Pklin = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(m_data_model->kk, Pk_eTNS_terms[10], "Spline"));

    //-------------------------

    m_data_model->funcs_pk.push_back(m_data_model->func_Pk_DeltaDelta);
    m_data_model->funcs_pk.push_back(m_data_model->func_Pk_DeltaTheta);
    m_data_model->funcs_pk.push_back(m_data_model->func_Pk_ThetaTheta);
    m_data_model->funcs_pk.push_back(m_data_model->func_Pk_A11);
    m_data_model->funcs_pk.push_back(m_data_model->func_Pk_A12);
    m_data_model->funcs_pk.push_back(m_data_model->func_Pk_A22);
    m_data_model->funcs_pk.push_back(m_data_model->func_Pk_A23);
    m_data_model->funcs_pk.push_back(m_data_model->func_Pk_A33);
    m_data_model->funcs_pk.push_back(m_data_model->func_Pk_B12);
    m_data_model->funcs_pk.push_back(m_data_model->func_Pk_B13);
    m_data_model->funcs_pk.push_back(m_data_model->func_Pk_B14);
    m_data_model->funcs_pk.push_back(m_data_model->func_Pk_B22);
    m_data_model->funcs_pk.push_back(m_data_model->func_Pk_B23);
    m_data_model->funcs_pk.push_back(m_data_model->func_Pk_B24);
    m_data_model->funcs_pk.push_back(m_data_model->func_Pk_B33);
    m_data_model->funcs_pk.push_back(m_data_model->func_Pk_B34);
    m_data_model->funcs_pk.push_back(m_data_model->func_Pk_B44);

    m_data_model->funcs_pk.push_back(m_data_model->func_Pk_b2d);
    m_data_model->funcs_pk.push_back(m_data_model->func_Pk_b2v);
    m_data_model->funcs_pk.push_back(m_data_model->func_Pk_b22);
    m_data_model->funcs_pk.push_back(m_data_model->func_Pk_bs2d);
    m_data_model->funcs_pk.push_back(m_data_model->func_Pk_bs2v);
    m_data_model->funcs_pk.push_back(m_data_model->func_Pk_b2s2);
    m_data_model->funcs_pk.push_back(m_data_model->func_Pk_bs22);
    m_data_model->funcs_pk.push_back(m_data_model->func_sigma32Pklin);
  }

  else ErrorCBL("the chosen model ("+m_data_model->Pk_mu_model+") is not currently implemented!", "set_fiducial_PkDM", "Modelling_TwoPointCorrelation_multipoles.cpp");
}


// ============================================================================================


void cbl::modelling::twopt::Modelling_TwoPointCorrelation_multipoles::set_fiducial_xiDM ()
{
  cout << endl; coutCBL << "Setting up the fiducial two-point correlation function model" << endl;

  m_data_model->nmultipoles = 3;

  const vector<double> rad = linear_bin_vector(m_data_model->step, m_data_model->r_min, m_data_model->r_max);

  m_data_model->rr = rad;
  m_data_model->kk = logarithmic_bin_vector(m_data_model->step, max(m_data_model->k_min, 1.e-4), min(m_data_model->k_max, 500.));

  cosmology::PkXi PX(m_data_model->cosmology);
  
  vector<double> Pk = PX.Pk_matter(m_data_model->kk, m_data_model->method_Pk, false, m_data_model->redshift, m_data_model->store_output, m_data_model->output_root, m_data_model->norm, m_data_model->k_min, m_data_model->k_max, m_data_model->prec, m_data_model->file_par);
  vector<double> PkNW = PX.Pk_matter(m_data_model->kk, "EisensteinHu", false, m_data_model->redshift, m_data_model->store_output, m_data_model->output_root, m_data_model->norm, m_data_model->k_min, m_data_model->k_max, m_data_model->prec, m_data_model->file_par);

  m_data_model->func_Pk = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(m_data_model->kk, Pk, "Spline"));
  m_data_model->func_Pk_NW = make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(m_data_model->kk, PkNW, "Spline"));

  std::vector<double> template_parameters = {m_data_model->sigmaNL_perp, m_data_model->sigmaNL_par, m_data_model->linear_growth_rate_z, m_data_model->bias, 0.};
  cbl::Print(template_parameters);

  vector<vector<double>> xil = Xi_l(rad, m_data_model->nmultipoles, "dispersion_dewiggled", {m_data_model->sigmaNL_perp, m_data_model->sigmaNL_par, m_data_model->linear_growth_rate_z, m_data_model->bias, 0.}, {m_data_model->func_Pk, m_data_model->func_Pk_NW}, m_data_model->prec, 1, 1);

  m_data_model->func_multipoles.erase(m_data_model->func_multipoles.begin(), m_data_model->func_multipoles.end());
  for (int i=0; i< m_data_model->nmultipoles; i++)
    m_data_model->func_multipoles.push_back(make_shared<cbl::glob::FuncGrid>(cbl::glob::FuncGrid(rad, xil[i], "Spline")));
}


// ============================================================================================


void cbl::modelling::twopt::Modelling_TwoPointCorrelation_multipoles::set_model_fullShape_sigma8_bias (const statistics::PriorDistribution sigma8_prior, const statistics::PriorDistribution bias_prior)
{
  // compute the fiducial dark matter two-point correlation function
  set_fiducial_PkDM();

  m_data_model->nmultipoles = m_nmultipoles;
  m_data_model->dataset_order = m_multipoles_order;
  m_data_model->use_pole = m_use_pole;

  // set the model parameters
  const int nparameters = 2;

  vector<statistics::ParameterType> parameterType(nparameters, statistics::ParameterType::_Base_);

  vector<string> parameterName(nparameters);
  parameterName[0] = "sigma8";
  parameterName[1] = "bias";

  vector<statistics::PriorDistribution> priors = {sigma8_prior, bias_prior};

  //set the priors
  m_set_prior(priors);

  // construct the model
  m_model = make_shared<statistics::Model1D>(statistics::Model1D(&xiMultipoles_sigma8_bias, nparameters, parameterType, parameterName, m_data_model));
  m_ModelIsSet = true;
}


// ============================================================================================


void cbl::modelling::twopt::Modelling_TwoPointCorrelation_multipoles::set_model_BAO (const statistics::PriorDistribution alpha_perpendicular_prior, const statistics::PriorDistribution alpha_parallel_prior, const statistics::PriorDistribution B0_prior, const statistics::PriorDistribution B2_prior, const statistics::PriorDistribution A00_prior, const statistics::PriorDistribution A20_prior, const statistics::PriorDistribution A01_prior, const statistics::PriorDistribution A21_prior, const statistics::PriorDistribution A02_prior, const statistics::PriorDistribution A22_prior, const bool compute_XiTemplate, const bool isRealSpace)
{

  // compute the fiducial dark matter two-point correlation function

  if (m_nmultipoles>2)
    ErrorCBL("BAO modelling can be done only with two multipoles!", "set_model_BAO", "Modelling_TwoPointCorrelation_multipoles");

  if (isRealSpace) {
    double lgf = m_data_model->linear_growth_rate_z;
    m_data_model->linear_growth_rate_z = 0.;

    set_fiducial_xiDM();

    m_data_model->linear_growth_rate_z = lgf;
  }
  else if (compute_XiTemplate)
    set_fiducial_xiDM();

  m_data_model->nmultipoles = m_nmultipoles;
  m_data_model->dataset_order = m_multipoles_order;
  m_data_model->use_pole = m_use_pole;

  // set the model parameters
  const int nparameters = 10;

  vector<statistics::ParameterType> parameterType(nparameters, statistics::ParameterType::_Base_);

  vector<string> parameterName(nparameters);
  parameterName[0] = "alpha_perpendicular";
  parameterName[1] = "alpha_parallel";
  parameterName[2] = "B0";
  parameterName[3] = "B2";
  parameterName[4] = "A00";
  parameterName[5] = "A20";
  parameterName[6] = "A01";
  parameterName[7] = "A21";
  parameterName[8] = "A02";
  parameterName[9] = "A22";

  vector<statistics::PriorDistribution> priors = {alpha_perpendicular_prior, alpha_parallel_prior, B0_prior, B2_prior, A00_prior, A20_prior, A01_prior, A21_prior, A02_prior, A22_prior};

  //set the priors
  m_set_prior(priors);

  // construct the model
  m_model = make_shared<statistics::Model1D>(statistics::Model1D(&xiMultipoles_BAO, nparameters, parameterType, parameterName, m_data_model));
  m_ModelIsSet = true;
}


// ============================================================================================


void cbl::modelling::twopt::Modelling_TwoPointCorrelation_multipoles::set_model_fullShape_DeWiggled (const statistics::PriorDistribution alpha_perpendicular_prior, const statistics::PriorDistribution alpha_parallel_prior, const statistics::PriorDistribution SigmaNL_perpendicular_prior, const statistics::PriorDistribution SigmaNL_parallel_prior, const statistics::PriorDistribution fsigma8_prior, const statistics::PriorDistribution bsigma8_prior, const statistics::PriorDistribution SigmaS_prior, const bool compute_PkDM)
{
  m_data_model->Pk_mu_model = "dispersion_dewiggled";

  // compute the fiducial dark matter two-point correlation function
  if (compute_PkDM) set_fiducial_PkDM();

  m_data_model->nmultipoles = m_nmultipoles;
  m_data_model->dataset_order = m_multipoles_order;
  m_data_model->use_pole = m_use_pole;

  // set the model parameters
  const int nparameters = 7;

  vector<statistics::ParameterType> parameterType(nparameters, statistics::ParameterType::_Base_);

  vector<string> parameterName(nparameters);
  parameterName[0] = "alpha_perpendicular";
  parameterName[1] = "alpha_parallel";
  parameterName[2] = "SigmaNL_perpendicular";
  parameterName[3] = "SigmaNL_parallel";
  parameterName[4] = "f*sigma8";
  parameterName[5] = "b*sigma8";
  parameterName[6] = "Sigma_S";

  vector<statistics::PriorDistribution> priors = {alpha_perpendicular_prior, alpha_parallel_prior, SigmaNL_perpendicular_prior, SigmaNL_parallel_prior, fsigma8_prior, bsigma8_prior, SigmaS_prior};

  //set the priors
  m_set_prior(priors);

  // construct the model
  m_model = make_shared<statistics::Model1D>(statistics::Model1D(&xiMultipoles, nparameters, parameterType, parameterName, m_data_model));
  m_ModelIsSet = true;
}


// ============================================================================================


void cbl::modelling::twopt::Modelling_TwoPointCorrelation_multipoles::set_model_fullShape_ModeCoupling (const statistics::PriorDistribution alpha_perpendicular_prior, const statistics::PriorDistribution alpha_parallel_prior, const statistics::PriorDistribution fsigma8_prior, const statistics::PriorDistribution bsigma8_prior, const statistics::PriorDistribution SigmaV_prior, const statistics::PriorDistribution AMC_prior, const bool compute_PkDM)
{
  m_data_model->Pk_mu_model = "dispersion_modecoupling";

  // compute the fiducial dark matter two-point correlation function
  if (compute_PkDM) set_fiducial_PkDM();

  m_data_model->nmultipoles = m_nmultipoles;
  m_data_model->dataset_order = m_multipoles_order;
  m_data_model->use_pole = m_use_pole;

  // set the model parameters
  const int nparameters = 6;

  vector<statistics::ParameterType> parameterType(nparameters, statistics::ParameterType::_Base_);

  vector<string> parameterName(nparameters);
  parameterName[0] = "alpha_perpendicular";
  parameterName[1] = "alpha_parallel";
  parameterName[2] = "f*sigma8";
  parameterName[3] = "b*sigma8";
  parameterName[4] = "sigma_v";
  parameterName[5] = "AMC";

  vector<statistics::PriorDistribution> priors = {alpha_perpendicular_prior, alpha_parallel_prior, fsigma8_prior, bsigma8_prior, SigmaV_prior, AMC_prior};

  //set the priors
  m_set_prior(priors);

  // construct the model
  m_model = make_shared<statistics::Model1D>(statistics::Model1D(&xiMultipoles, nparameters, parameterType, parameterName, m_data_model));
  m_ModelIsSet = true;
}


// ============================================================================================


void cbl::modelling::twopt::Modelling_TwoPointCorrelation_multipoles::set_model_dispersion (const statistics::PriorDistribution fsigma8_prior, const statistics::PriorDistribution bsigma8_prior, const statistics::PriorDistribution sigmav_prior, const statistics::PriorDistribution alpha_perpendicular_prior, const statistics::PriorDistribution alpha_parallel_prior, const bool DFoG, const bool compute_PkDM)
{
  if (DFoG) m_data_model->Pk_mu_model = "dispersion_Gauss";
  else m_data_model->Pk_mu_model = "dispersion_Lorentz";

  // compute the fiducial dark matter two-point correlation function
  if (compute_PkDM) set_fiducial_PkDM();

  m_data_model->nmultipoles = m_nmultipoles;
  m_data_model->dataset_order = m_multipoles_order;
  m_data_model->use_pole = m_use_pole;

  // set the model parameters
  const int nparameters = 5;

  vector<statistics::ParameterType> parameterType(nparameters, statistics::ParameterType::_Base_);

  vector<string> parameterName(nparameters);
  parameterName[0] = "f*sigma8";
  parameterName[1] = "b*sigma8";
  parameterName[2] = "sigmav";
  parameterName[3] = "alpha_perpendicular";
  parameterName[4] = "alpha_parallel";
  
  vector<statistics::PriorDistribution> priors = {fsigma8_prior, bsigma8_prior, sigmav_prior, alpha_perpendicular_prior, alpha_parallel_prior};

  //set the priors
  m_set_prior(priors);

  // construct the model
  m_model = make_shared<statistics::Model1D>(statistics::Model1D(&xiMultipoles, nparameters, parameterType, parameterName, m_data_model));
  m_ModelIsSet = true;
}


// ============================================================================================


void cbl::modelling::twopt::Modelling_TwoPointCorrelation_multipoles::set_model_Scoccimarro (const statistics::PriorDistribution fsigma8_prior, const statistics::PriorDistribution bsigma8_prior, const statistics::PriorDistribution sigmav_prior, const statistics::PriorDistribution alpha_perpendicular_prior, const statistics::PriorDistribution alpha_parallel_prior, const bool DFoG, const bool compute_PkDM)
{
  if (DFoG) m_data_model->Pk_mu_model = "Scoccimarro_Gauss";
  else m_data_model->Pk_mu_model = "Scoccimarro_Lorentz";

  // compute the fiducial dark matter two-point correlation function
  if (compute_PkDM) set_fiducial_PkDM();

  m_data_model->nmultipoles = m_nmultipoles;
  m_data_model->dataset_order = m_multipoles_order;
  m_data_model->use_pole = m_use_pole;

  // set the model parameters
  const int nparameters = 5;

  vector<statistics::ParameterType> parameterType(nparameters, statistics::ParameterType::_Base_);

  vector<string> parameterName(nparameters);
  parameterName[0] = "f*sigma8";
  parameterName[1] = "b*sigma8";
  parameterName[2] = "sigmav";  
  parameterName[3] = "alpha_perpendicular";
  parameterName[4] = "alpha_parallel";

  vector<statistics::PriorDistribution> priors = {fsigma8_prior, bsigma8_prior, sigmav_prior, alpha_perpendicular_prior, alpha_parallel_prior};

  //set the priors
  m_set_prior(priors);

  // construct the model
  m_model = make_shared<statistics::Model1D>(statistics::Model1D(&xiMultipoles, nparameters, parameterType, parameterName, m_data_model));
  m_ModelIsSet = true;
}


// ============================================================================================


void cbl::modelling::twopt::Modelling_TwoPointCorrelation_multipoles::set_model_Scoccimarro_fitPezzotta (const statistics::PriorDistribution fsigma8_prior, const statistics::PriorDistribution bsigma8_prior, const statistics::PriorDistribution sigmav_prior, const statistics::PriorDistribution kd_prior, const statistics::PriorDistribution kt_prior, const statistics::PriorDistribution alpha_perpendicular_prior, const statistics::PriorDistribution alpha_parallel_prior, const bool DFoG, const bool compute_PkDM)
{
  if (DFoG) m_data_model->Pk_mu_model = "Scoccimarro_Pezzotta_Gauss";
  else m_data_model->Pk_mu_model = "Scoccimarro_Pezzotta_Lorentz";

  // compute the fiducial dark matter two-point correlation function
  if (compute_PkDM) set_fiducial_PkDM();

  m_data_model->nmultipoles = m_nmultipoles;
  m_data_model->dataset_order = m_multipoles_order;
  m_data_model->use_pole = m_use_pole;

  // set the model parameters
  const int nparameters = 7;

  vector<statistics::ParameterType> parameterType(nparameters, statistics::ParameterType::_Base_);

  vector<string> parameterName(nparameters);
  parameterName[0] = "f*sigma8";
  parameterName[1] = "b*sigma8";
  parameterName[2] = "sigmav";
  parameterName[3] = "kd";
  parameterName[4] = "kt";
  parameterName[5] = "alpha_perpendicular";
  parameterName[6] = "alpha_parallel";

  vector<statistics::PriorDistribution> priors = {fsigma8_prior, bsigma8_prior, sigmav_prior, kd_prior, kt_prior, alpha_perpendicular_prior, alpha_parallel_prior};

  //set the priors
  m_set_prior(priors);

  // construct the model
  m_model = make_shared<statistics::Model1D>(statistics::Model1D(&xiMultipoles, nparameters, parameterType, parameterName, m_data_model));
  m_ModelIsSet = true;
}


// ============================================================================================


void cbl::modelling::twopt::Modelling_TwoPointCorrelation_multipoles::set_model_Scoccimarro_fitBel (const statistics::PriorDistribution fsigma8_prior, const statistics::PriorDistribution bsigma8_prior, const statistics::PriorDistribution sigmav_prior, const statistics::PriorDistribution kd_prior, const statistics::PriorDistribution bb_prior, const statistics::PriorDistribution a1_prior, const statistics::PriorDistribution a2_prior, const statistics::PriorDistribution a3_prior, const statistics::PriorDistribution alpha_perpendicular_prior, const statistics::PriorDistribution alpha_parallel_prior, const bool DFoG, const bool compute_PkDM)
{
  if (DFoG) m_data_model->Pk_mu_model = "Scoccimarro_Bel_Gauss";
  else m_data_model->Pk_mu_model = "Scoccimarro_Bel_Lorentz";

  // compute the fiducial dark matter two-point correlation function
  if (compute_PkDM) set_fiducial_PkDM();

  m_data_model->nmultipoles = m_nmultipoles;
  m_data_model->dataset_order = m_multipoles_order;
  m_data_model->use_pole = m_use_pole;

  // set the model parameters
  const int nparameters = 10;

  vector<statistics::ParameterType> parameterType(nparameters, statistics::ParameterType::_Base_);

  vector<string> parameterName(nparameters);
  parameterName[0] = "f*sigma8";
  parameterName[1] = "b*sigma8";
  parameterName[2] = "sigmav";
  parameterName[3] = "kd";
  parameterName[4] = "a0";
  parameterName[5] = "a1";
  parameterName[6] = "a2";
  parameterName[7] = "a3";
  parameterName[8] = "alpha_perpendicular";
  parameterName[9] = "alpha_parallel";

  vector<statistics::PriorDistribution> priors = {fsigma8_prior, bsigma8_prior, sigmav_prior, kd_prior, bb_prior, a1_prior, a2_prior, a3_prior, alpha_perpendicular_prior, alpha_parallel_prior};

  //set the priors
  m_set_prior(priors);

  // construct the model
  m_model = make_shared<statistics::Model1D>(statistics::Model1D(&xiMultipoles, nparameters, parameterType, parameterName, m_data_model));
  m_ModelIsSet = true;
}


// ============================================================================================


void cbl::modelling::twopt::Modelling_TwoPointCorrelation_multipoles::set_model_TNS (const statistics::PriorDistribution fsigma8_prior, const statistics::PriorDistribution bsigma8_prior, const statistics::PriorDistribution sigmav_prior, const statistics::PriorDistribution alpha_perpendicular_prior, const statistics::PriorDistribution alpha_parallel_prior, const bool DFoG, const bool compute_PkDM)
{
  if (DFoG) m_data_model->Pk_mu_model = "TNS_Gauss";
  else m_data_model->Pk_mu_model = "TNS_Lorentz";

  // compute the fiducial dark matter two-point correlation function
  if (compute_PkDM) set_fiducial_PkDM();

  m_data_model->nmultipoles = m_nmultipoles;
  m_data_model->dataset_order = m_multipoles_order;
  m_data_model->use_pole = m_use_pole;

  // set the model parameters
  const int nparameters = 5;

  vector<statistics::ParameterType> parameterType(nparameters, statistics::ParameterType::_Base_);

  vector<string> parameterName(nparameters);
  parameterName[0] = "f*sigma8";
  parameterName[1] = "b*sigma8";
  parameterName[2] = "sigmav";
  parameterName[3] = "alpha_perpendicular";
  parameterName[4] = "alpha_parallel";

  vector<statistics::PriorDistribution> priors = {fsigma8_prior, bsigma8_prior, sigmav_prior, alpha_perpendicular_prior, alpha_parallel_prior};

  //set the priors
  m_set_prior(priors);

  // construct the model
  m_model = make_shared<statistics::Model1D>(statistics::Model1D(&xiMultipoles, nparameters, parameterType, parameterName, m_data_model));
  m_ModelIsSet = true;
}


// ============================================================================================


void cbl::modelling::twopt::Modelling_TwoPointCorrelation_multipoles::set_model_eTNS (const statistics::PriorDistribution fsigma8_prior, const statistics::PriorDistribution b1sigma8_prior, const statistics::PriorDistribution b2sigma8_prior, const statistics::PriorDistribution sigmav_prior, const statistics::PriorDistribution alpha_perpendicular_prior, const statistics::PriorDistribution alpha_parallel_prior, const bool DFoG, const bool compute_PkDM)
{
  if (DFoG) m_data_model->Pk_mu_model = "eTNS_Gauss";
  else m_data_model->Pk_mu_model = "eTNS_Lorentz";

  // compute the fiducial dark matter two-point correlation function
  if (compute_PkDM) set_fiducial_PkDM();

  m_data_model->nmultipoles = m_nmultipoles;
  m_data_model->dataset_order = m_multipoles_order;
  m_data_model->use_pole = m_use_pole;

  // set the model parameters
  const int nparameters = 6;

  vector<statistics::ParameterType> parameterType(nparameters, statistics::ParameterType::_Base_);

  vector<string> parameterName(nparameters);
  parameterName[0] = "f*sigma8";
  parameterName[1] = "b1*sigma8";
  parameterName[2] = "b2*sigma8";
  parameterName[3] = "sigmav";
  parameterName[4] = "alpha_perpendicular";
  parameterName[5] = "alpha_parallel";

  vector<statistics::PriorDistribution> priors = {fsigma8_prior, b1sigma8_prior, b2sigma8_prior, sigmav_prior, alpha_perpendicular_prior, alpha_parallel_prior};

  //set the priors
  m_set_prior(priors);

  // construct the model
  m_model = make_shared<statistics::Model1D>(statistics::Model1D(&xiMultipoles, nparameters, parameterType, parameterName, m_data_model));
  m_ModelIsSet = true;
}


// ============================================================================================


void cbl::modelling::twopt::Modelling_TwoPointCorrelation_multipoles::set_model_linear_theoretical_bias_BAO (const std::vector<std::string> cosmo_param, const std::vector<statistics::PriorDistribution> cosmo_param_prior, const statistics::PriorDistribution sigmaNL_prior)
{
  //this is a linear model, without FoG effect. The dispersion model is used without damping in real and redshift space. We account for gaussian damping only in presence of photometric errors
  m_data_model->Pk_mu_model = "dispersion_Gauss";
  
  m_data_model->nmultipoles = m_nmultipoles;
  m_data_model->dataset_order = m_multipoles_order;
  m_data_model->use_pole = m_use_pole;

  m_data_model->Cpar = cosmo_param;

  // set the model parameters
  const int nparameters = (cosmo_param.size())+2;

  vector<statistics::ParameterType> parameterType(nparameters, statistics::ParameterType::_Base_);
  parameterType[nparameters-1] = statistics::ParameterType::_Derived_;

  vector<string> parameterName(nparameters);

  for (size_t i=0; i<cosmo_param.size(); i++) 
    parameterName[i] = cosmo_param[i];
  parameterName[cosmo_param.size()] = "sigmaNL";
  parameterName[cosmo_param.size()+1]="bias";

  vector<statistics::PriorDistribution> priors;
  for (size_t i=0; i<cosmo_param.size(); i++)
    priors.emplace_back(cosmo_param_prior[i]); 
  priors.emplace_back(sigmaNL_prior);

  //set the priors
  m_set_prior(priors);

  // construct the model
  m_model = make_shared<statistics::Model1D>(statistics::Model1D(&xiMultipoles_linear_theoretical_bias_BAO, nparameters, parameterType, parameterName, m_data_model));

}


// ============================================================================================


void cbl::modelling::twopt::Modelling_TwoPointCorrelation_multipoles::set_model_scaling_relation_sigmaz_cosmology (const std::vector<std::string> cosmo_param, const std::vector<statistics::PriorDistribution> cosmo_prior, const statistics::PriorDistribution alpha_prior, const statistics::PriorDistribution beta_prior, const statistics::PriorDistribution gamma_prior, const statistics::PriorDistribution scatter0_prior, const statistics::PriorDistribution scatterM_prior, const statistics::PriorDistribution scatterM_exponent_prior, const statistics::PriorDistribution scatterz_prior, const statistics::PriorDistribution scatterz_exponent_prior, const statistics::PriorDistribution sigmaz_prior, const std::string z_evo)
{
  // compute the fiducial dark matter two-point correlation function
  m_data_model->nmultipoles = m_nmultipoles;
  m_data_model->dataset_order = m_multipoles_order;
  m_data_model->use_pole = m_use_pole;


  // Set the scaling relation Modelling object
  if ((m_data_model->scaling_relation)->data()->xx().size() == 0)
    ErrorCBL("The mass-observable relation is not set! Use the correct set_data_model().", "set_model_scaling_relation_sigmaz_cosmology", "Modelling_TwoPointCorrelation1D_monopole.cpp");

  (m_data_model->scaling_relation)->set_model_MassObservableRelation_cosmology(z_evo, cosmo_param, cosmo_prior, alpha_prior, beta_prior, gamma_prior, scatter0_prior, scatterM_prior, scatterM_exponent_prior, scatterz_prior, scatterz_exponent_prior);

  (m_data_model->scaling_relation)->set_likelihood(cbl::statistics::LikelihoodType::_Gaussian_Error_, {}); // Set the likelihood for the scaling relation (only to avoid internal errors, of course it is not used)

  // Set the parameter names and priors
  m_data_model->Cpar = cosmo_param;

  const size_t nParams_scaling_relation = (m_data_model->scaling_relation)->likelihood()->parameters()->nparameters();
  const size_t nParams_base = 1 + nParams_scaling_relation; // sigmaz + scaling relation parameters (which include the cosmological parameters)
  const size_t nParams_derived = 1; // the effective bias is derived from the scaling relation
  
  const size_t nParams = nParams_base + nParams_derived;

  vector<statistics::ParameterType> Par_type (nParams, statistics::ParameterType::_Base_);
  Par_type[nParams-1] = statistics::ParameterType::_Derived_;

  vector<string> Par_string (nParams);
  std::vector<statistics::PriorDistribution> param_prior (nParams_base);

  // Cosmological and scaling relation parameters
  for (size_t i=0; i<nParams_base-1; i++) {
    Par_string[i] = (m_data_model->scaling_relation)->likelihood()->parameters()->name(i);
    param_prior[i] = *(m_data_model->scaling_relation)->get_prior(i);
  }

  // sigmaz
  Par_string[nParams_base-1] = "sigmaz";
  param_prior[nParams_base-1] = sigmaz_prior;

  // bias
  Par_string[nParams-1] = "bias";

  // set prior
  m_set_prior(param_prior);

  // construct the model
  m_model = make_shared<statistics::Model1D>(statistics::Model1D(&xil_damped_scaling_relation_sigmaz_cosmology, nParams, Par_type, Par_string, m_data_model));
}
