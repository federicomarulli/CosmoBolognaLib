/*******************************************************************
 *  Copyright (C) 2010 by Federico Marulli                         *
 *  federico.marulli3@unibo.it                                     *
 *                                                                 *
 *  This program is free software; you can redistribute it and/or  *
 *  modify it under the terms of the GNU General Public License as *
 *  published by the Free Software Foundation; either version 2 of *
 *  the License, or (at your option) any later version.            *
 *                                                                 *
 *  This program is distributed in the hope that it will be useful,*
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the  *
 *  GNU General Public License for more details.                   *
 *                                                                 *
 *  You should have received a copy of the GNU General Public      *
 *  License along with this program; if not, write to the Free     *
 *  Software Foundation, Inc.,                                     *
 *  59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.      *
 *******************************************************************/

/**
 *  @file Headers/Lib/ModelFunction.h
 *
 *  @brief Functions to model data 
 *
 *  This file contains the prototypes of a set of functions to model
 *  data
 *  
 *  @author Federico Marulli 
 *
 *  @author federico.marulli3@unbo.it
 */

#ifndef __MODFUNC__
#define __MODFUNC__


// ============================================================================

#include "Cosmology.h"


namespace cosmobl {

  namespace glob{

    struct STR_twop_model{

      shared_ptr<Cosmology> cosmology;
      vector<CosmoPar> Cpar;

      double redshift;
      double f;
      string method;
      string output_root;
      bool NL;
      int norm;
      double k_min;
      double k_max;
      double aa;
      bool GSL;
      double prec;
      string file_par;
      double r_min;
      double r_max;

      shared_ptr<classfunc::func_grid_GSL> func_xi;

      double pi_max;

      STR_twop_model() {}

    };

    double wp_bias(double r, shared_ptr<void> parameters, vector<double> model_parameters); 

    double wp_bias_cosmology(double r, shared_ptr<void> parameters, vector<double> model_parameters); 

    double xi_bias(double r, shared_ptr<void> parameters, vector<double> model_parameters); 

    double xi_bias_cosmology(double r, shared_ptr<void> parameters, vector<double> model_parameters); 

    double xi0_bias(double r, shared_ptr<void> parameters, vector<double> model_parameters); 

    double xi0_bias_cosmology(double r, shared_ptr<void> parameters, vector<double> model_parameters); 

    double xi_alpha_B(double r, shared_ptr<void> parameters, vector<double> model_parameters);

    double xi_alpha_B_poly(double r, shared_ptr<void> parameters, vector<double> model_parameters);

    double xi0_alpha_bias(double r, shared_ptr<void> parameters, vector<double> model_parameters);

    double xi0_alpha_bias_cosmology(double r, shared_ptr<void> parameters, vector<double> model_parameters);

    vector<double> wp_bias_vector(vector<double> r, shared_ptr<void> parameters, vector<double> model_parameters); 

    vector<double> wp_bias_cosmology_vector(vector<double> r, shared_ptr<void> parameters, vector<double> model_parameters); 

    vector<double> xi_bias_vector(vector<double> r, shared_ptr<void> parameters, vector<double> model_parameters); 

    vector<double> xi_bias_cosmology_vector(vector<double> r, shared_ptr<void> parameters, vector<double> model_parameters); 

    vector<double> xi0_bias_vector(vector<double> r, shared_ptr<void> parameters, vector<double> model_parameters); 

    vector<double> xi0_bias_cosmology_vector(vector<double> r, shared_ptr<void> parameters, vector<double> model_parameters); 

    vector<double> xi_alpha_B_vector(vector<double> r, shared_ptr<void> parameters, vector<double> model_parameters);

    vector<double> xi_alpha_B_poly_vector(vector<double> r, shared_ptr<void> parameters, vector<double> model_parameters);

    vector<double> xi0_alpha_bias_vector(vector<double> r, shared_ptr<void> parameters, vector<double> model_parameters);

    vector<double> xi0_alpha_bias_cosmology_vector(vector<double> r, shared_ptr<void> parameters, vector<double> model_parameters);

  }
}

#endif
