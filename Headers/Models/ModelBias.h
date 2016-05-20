/********************************************************************
 *  Copyright (C) 2014 by Federico Marulli and Alfonso Veropalumbo  *
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
 *  @file Headers/Models/ModelBias.h
 *
 *  @brief The class ModelBias
 *
 *  This file defines the interface of the class ModelBias
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#ifndef __MODELBIAS__
#define __MODELBIAS__

#include "ModelFunction.h"

namespace cosmobl{

  class ModelBias : public statistics::Model1D
  {
    public:
      ModelBias() {}

      ModelBias (const double bias_value, const statistics::Prior bias_prior);

      ModelBias (const double bias_value, const statistics::Prior bias_prior, const vector<CosmoPar> cosmo_parameters, const vector<double> cosmo_parameters_values, const vector<statistics::Prior> cosmo_parameters_priors, const vector<string> cpar_name={});

      void set_xi_parameters (const vector<double> r,  const shared_ptr<Cosmology> cosmology, const double redshift, const string type="Poly", const int nPt=4, const string method="CAMB", const string output_root="test", const bool NL=1, const int norm=-1, const double k_min=0., const double k_max=100., const double aa=0., const bool GSL=1, const double prec=1.e-2, const string file_par=par::defaultString);

      void set_xi_parameters_cosmology (const shared_ptr<Cosmology> cosmology, const vector<CosmoPar> cosmo_parameters, const double redshift, const string method="CAMB", const string output_root="test", const bool NL=1, const int norm=-1, const double k_min=0., const double k_max=100., const double aa=0., const bool GSL=1, const double prec=1.e-2, const string file_par=par::defaultString);

      void set_wp_parameters (const vector<double> r, const shared_ptr<Cosmology> cosmology, const double redshift, const double pi_max, const string type="Poly", const int nPt=4, const string method="CAMB", const string output_root="test", const bool NL=1, const int norm=-1, const double r_min=1.e-3, const double r_max=350., const double k_min=0., const double k_max=100., const double aa=0., const bool GSL=1, const double prec=1.e-2, const string file_par=par::defaultString);

      void set_wp_parameters_cosmology (const shared_ptr<Cosmology> cosmology, const vector<CosmoPar> cosmo_parameters, const double redshift, const double pi_max, const string method="CAMB", const string output_root="test", const int norm=-1, const double r_min=1.e-3, const double r_max=350., const double k_min=0., const double k_max=100., const double aa=0., const bool GSL=1, const double prec=1.e-2, const string file_par=par::defaultString);
      
      void set_xi0_parameters (const vector<double> r,  const shared_ptr<Cosmology> cosmology, const double redshift, const string type="Poly", const int nPt=4, const string method="CAMB", const string output_root="test", const bool NL=1, const int norm=-1, const double k_min=0., const double k_max=100., const double aa=0., const bool GSL=1, const double prec=1.e-2, const string file_par=par::defaultString);

      void set_xi0_parameters_cosmology (const shared_ptr<Cosmology> cosmology, const vector<CosmoPar> cosmo_parameters, const double redshift, const string method="CAMB", const string output_root="test", const bool NL=1, const int norm=-1, const double k_min=0., const double k_max=100., const double aa=0., const bool GSL=1, const double prec=1.e-2, const string file_par=par::defaultString);

      shared_ptr<cosmobl::glob::STR_twop_model> model_parameters() 
      { shared_ptr<cosmobl::glob::STR_twop_model> pp = static_pointer_cast<cosmobl::glob::STR_twop_model>(m_model_parameters); return pp;}


  };
}
#endif
