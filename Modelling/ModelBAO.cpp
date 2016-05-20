/********************************************************************
 *  Copyright (C) 2015 by Federico Marulli and Alfonso Veropalumbo  *
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
 *  @file Modelling/ModelBAO.cpp
 *
 *  @brief Methods of the class ModelBAO, used for modelling the baryon
 *  acoustic oscillation peak in two point correlation function
 *
 *  This file contains the implementation of the methods of the class
 *  ModelBAO
 *
 *  @authors Federico Marulli
 *
 *  @authors federico.marulli3@unbo.it
 */


#include "ModelBAO.h"

using namespace cosmobl;


// ============================================================================================


cosmobl::ModelBAO::ModelBAO (const double bias_value, const statistics::Prior bias_prior, const double alpha_value, const statistics::Prior alpha_prior, const bool AddPoly, const vector<double> r, const vector<double> xi, const string method) : Model1D()
{

  if(AddPoly){
    m_npar = 5;
    m_parameters.resize(m_npar); 
    m_parameters[0] = make_shared<statistics::Parameter>(bias_value, bias_prior, 0, "bias");
    m_parameters[1] = make_shared<statistics::Parameter>(alpha_value, alpha_prior, 0, "alpha");

    vector<double> vv={};
    m_parameters[2] = make_shared<statistics::Parameter>(0., -1.e3, 1.e3, 0, vv, "A0");
    m_parameters[3] = make_shared<statistics::Parameter>(0., -1.e2, 1.e2, 0, vv, "A1");
    m_parameters[4] = make_shared<statistics::Parameter>(0., -1.e2, 1.e2, 0, vv, "A2");

    m_model = &cosmobl::glob::xi_alpha_B_poly;
    m_model_vector = &cosmobl::glob::xi_alpha_B_poly_vector;

  }
  else{
    m_npar = 2;
    m_parameters.resize(m_npar); 
    m_parameters[0] = make_shared<statistics::Parameter>(bias_value, bias_prior, 0, "bias");
    m_parameters[1] = make_shared<statistics::Parameter>(alpha_value, alpha_prior, 0, "alpha");

    m_model = &cosmobl::glob::xi_alpha_B;
    m_model_vector = &cosmobl::glob::xi_alpha_B_vector;

  }

  auto model_parameters = make_shared<cosmobl::glob::STR_twop_model>(cosmobl::glob::STR_twop_model());
  model_parameters->func_xi = make_shared<classfunc::func_grid_GSL>(classfunc::func_grid_GSL(r,xi,method));
  m_model_parameters = move(model_parameters);

  classfunc::func_grid_GSL f_xi (r,xi,method);

}
