/********************************************************************
 *  Copyright (C) 2017 by Federico Marulli                          *
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
 *  Modelling/ThreePointCorrelation/Modelling_ThreePointCorrelation.cpp
 *
 *  @brief Methods of the class Modelling_ThreePointCorrelation, used to
 *  model three-point correlation functions of any kind
 *
 *  This file contains the implementation of the methods of the class
 *  Modelling_ThreePointCorrelation
 *
 *  @authors Federico Marulli, Michele Moresco
 *
 *  @authors federico.marulli3@unbo.it, michele.moresco@unibo.it
 */


#include "Modelling_ThreePointCorrelation.h"
#include "Modelling_ThreePointCorrelation_angular_connected.h"
#include "Modelling_ThreePointCorrelation_angular_reduced.h"
#include "Modelling_ThreePointCorrelation_comoving_connected.h"
#include "Modelling_ThreePointCorrelation_comoving_reduced.h"

using namespace std;

using namespace cbl;
using namespace measure;
using namespace measure::threept;
using namespace modelling;
using namespace measure::threept;


// ============================================================================================


std::shared_ptr<cbl::modelling::threept::Modelling_ThreePointCorrelation> cbl::modelling::threept::Modelling_ThreePointCorrelation::Create (const std::shared_ptr<measure::threept::ThreePointCorrelation> threept)
{
  if (threept->threePType()==measure::threept::ThreePType::_angular_connected_)
    return move(unique_ptr<Modelling_ThreePointCorrelation_angular_connected> (new Modelling_ThreePointCorrelation_angular_connected(threept)));
  
  else if (threept->threePType()==measure::threept::ThreePType::_angular_reduced_)
    return move(unique_ptr<modelling::threept::Modelling_ThreePointCorrelation_angular_reduced> (new Modelling_ThreePointCorrelation_angular_reduced(threept)));

  else if (threept->threePType()==measure::threept::ThreePType::_comoving_connected_)
    return move(unique_ptr<Modelling_ThreePointCorrelation_comoving_connected> (new Modelling_ThreePointCorrelation_comoving_connected(threept)));

  else if (threept->threePType()==measure::threept::ThreePType::_comoving_reduced_)
    return move(unique_ptr<Modelling_ThreePointCorrelation_comoving_reduced> (new Modelling_ThreePointCorrelation_comoving_reduced(threept)));

  else ErrorCBL("no such type of object, or error in the input parameters!", "Create", "Modelling_ThreePointCorrelation.cpp");
 
  return NULL;
}


// ============================================================================================


std::shared_ptr<cbl::modelling::threept::Modelling_ThreePointCorrelation> cbl::modelling::threept::Modelling_ThreePointCorrelation::Create (const measure::threept::ThreePType threePType, const std::shared_ptr<data::Data> threept_dataset)
{
  if (threePType==measure::threept::ThreePType::_angular_connected_)
    return move(unique_ptr<Modelling_ThreePointCorrelation_angular_connected> (new Modelling_ThreePointCorrelation_angular_connected(threept_dataset)));

  else if (threePType==measure::threept::ThreePType::_angular_reduced_)
    return move(unique_ptr<Modelling_ThreePointCorrelation_angular_reduced> (new Modelling_ThreePointCorrelation_angular_reduced(threept_dataset)));

  else if (threePType==measure::threept::ThreePType::_comoving_connected_)
    return move(unique_ptr<Modelling_ThreePointCorrelation_comoving_connected> (new Modelling_ThreePointCorrelation_comoving_connected(threept_dataset)));

  else if (threePType==measure::threept::ThreePType::_comoving_reduced_)
    return move(unique_ptr<Modelling_ThreePointCorrelation_comoving_reduced> (new Modelling_ThreePointCorrelation_comoving_reduced(threept_dataset)));
  
  else ErrorCBL("no such type of object, or error in the input parameters!", "Create", "Modelling_ThreePointCorrelation.cpp");

  return NULL;
}

// ============================================================================================


void cbl::modelling::threept::Modelling_ThreePointCorrelation::set_data_model (const std::vector<double> Q_DM)
{
  m_data_model.Q_DM = Q_DM;
}


// ============================================================================================


void cbl::modelling::threept::Modelling_ThreePointCorrelation::set_data_Q_nonlocal (const cosmology::Cosmology cosmology, const double r1, const double r2, const std::vector<double> theta, const string model, const std::vector<double> kk, const std::vector<double> Pk_DM)
{
  m_data_model.cosmology = make_shared<cosmology::Cosmology>(cosmology);
  m_data_model.r1 = r1;
  m_data_model.r2 = r2;
  m_data_model.theta = theta;
  m_data_model.model = model;
  m_data_model.kk = kk;
  m_data_model.Pk_DM = Pk_DM;
}


// ============================================================================================


void cbl::modelling::threept::Modelling_ThreePointCorrelation::set_data_model_zeta_RSD (const double r1, const double r2, const cosmology::Cosmology cosmology, const double redshift, const string method_Pk, const bool NL, const int max_ll, const double k_min, const double k_max, const int step_k, const double r_min, const double r_max, const int step_r, const bool force_realSpace, const bool use_k, const string output_dir, const bool store_output, const string output_root, const int norm, const double prec)
{
  m_data_model.cosmology = make_shared<cosmology::Cosmology>(cosmology);

  m_data_model.r1 = r1;
  m_data_model.r2 = r2;
  m_data_model.redshift = redshift;
  m_data_model.method_Pk = method_Pk;
  m_data_model.NL = NL;
  m_data_model.k_min = k_min;
  m_data_model.k_max = k_max;
  m_data_model.step_k = step_k;
  m_data_model.kk = logarithmic_bin_vector(m_data_model.step_k, m_data_model.k_min, m_data_model.k_max);
  m_data_model.r_min = r_min;
  m_data_model.r_max = r_max;
  m_data_model.step_r = step_r;
  m_data_model.rr = linear_bin_vector(m_data_model.step_r, m_data_model.r_min, m_data_model.r_max);
  m_data_model.output_dir = output_dir;
  m_data_model.store_output = store_output;
  m_data_model.output_root = output_root;
  m_data_model.norm = norm;
  m_data_model.prec = prec;
  m_data_model.force_realSpace = force_realSpace;
  m_data_model.max_ll = max_ll; 
  m_data_model.use_k = use_k;
  

  try {
    m_data_model.sigma8_z = m_data_model.cosmology->sigma8(m_data_model.redshift);
  }
  catch(cbl::glob::Exception &exc) { 
    coutCBL << "sigma8 is not set, computing from the power spectrum, method_Pk = "+m_data_model.method_Pk << endl; 
    m_data_model.sigma8_z = m_data_model.cosmology->sigma8_Pk(m_data_model.method_Pk, m_data_model.redshift,m_data_model.store_output,  m_data_model.output_root);
  }
  m_data_model.linear_growth_rate_z = m_data_model.cosmology->linear_growth_rate(m_data_model.redshift, 1.);

  m_data_model.Pk_DM = m_data_model.cosmology->Pk(m_data_model.kk, m_data_model.method_Pk, m_data_model.NL, m_data_model.redshift, m_data_model.output_dir, m_data_model.store_output, m_data_model.output_root, m_data_model.norm, m_data_model.k_min, m_data_model.k_max, m_data_model.prec);
}
