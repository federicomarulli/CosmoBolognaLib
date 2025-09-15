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
 *  Modelling/ShearShearLensing/Modelling_PowerSpectrum_Shear.cpp
 *
 *  @brief Methods of the class
 *  Modelling_PowerSpectrum_Shear
 *
 *  This file contains the implementation of the methods of the class
 *  Modelling_PowerSpectrum_Shear, used to model the
 *  shear-shear angular power spectrum
 *
 */


#include "Modelling_PowerSpectrum_Shear.h"
#include "Data1D.h"

using namespace std;

using namespace cbl;

  
// ============================================================================================

   void cbl::modelling::sslensing::Modelling_PowerSpectrum_shear::set_data_model(const std::shared_ptr<cbl::cosmology::Cosmology> cosmology, const double z_min, const double z_max, const std::vector<std::shared_ptr<glob::FuncGrid>> Nz_all, const vector<double> pars_IA, const cbl::glob::FuncGrid L_ratio, const std::string method_Pk, const bool NL, const int norm, const double k_min, const double k_max, const std::string output_root, const std::string file_par, const double prec, const bool unit)
  {
    m_data_model = make_shared<STR_data_model>(STR_data_model());
  
    m_data_model->cosmology = move(cosmology);
    m_data_model->z_min = z_min;
    m_data_model->z_max = z_max;
    m_data_model->Nz_all = Nz_all;
    m_data_model->pars_IA = pars_IA;
    m_data_model->L_ratio = L_ratio;
    m_data_model->method_Pk = method_Pk;
    m_data_model->NL = NL;
    m_data_model->norm = norm;
    m_data_model->k_min = k_min;
    m_data_model->k_max = k_max;
    m_data_model->output_root = output_root;
    m_data_model->file_par = file_par;
    m_data_model->prec = prec;
    m_data_model->unit = unit;
  }

//=============================================================================================

   void cbl::modelling::sslensing::Modelling_PowerSpectrum_shear::set_model_Cij_LL_eNLA (const std::vector<std::string> cosmo_param, const std::vector<statistics::PriorDistribution> cosmo_prior)
   {
     m_data_model->Cpar = cosmo_param;

     const int nparameters = cosmo_param.size(); 
  
     vector<statistics::ParameterType> parameterType (nparameters, statistics::ParameterType::_Base_);

     vector<string> parameterName(nparameters);
     vector<statistics::PriorDistribution> prior(nparameters); // check!!!
  
     // Set the names and priors of the cosmological parameters
     for (size_t i=0; i<cosmo_param.size(); i++)
     {
       parameterName[i] = cosmo_param[i];
       prior[i] = cosmo_prior[i];
     }

     //set the prior
     m_set_prior(prior);
  
     // Construct the model
     m_model = make_shared<statistics::Model1D>(statistics::Model1D(&C_ij_LL_eNLA, nparameters, parameterType, parameterName, m_data_model));
   }


