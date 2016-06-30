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
 *  @file Headers/Models/ModelBAO.h
 *
 *  @brief The class ModelBAO
 *
 *  This file defines the interface of the class ModelBAO
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#ifndef __MODELBAO__
#define __MODELBAO__

#include "ModelFunction.h"

namespace cosmobl{

  class ModelBAO : public statistics::Model1D
  {
    public:
      ModelBAO() {}

      ModelBAO (const double bias_value, const statistics::Prior bias_prior, const double alpha_value, const statistics::Prior alpha_prior, const bool AddPoly);

      ~ModelBAO() {}

      void set_xi_parameters (const vector<double> r,  const shared_ptr<cosmology::Cosmology> cosmology, const double redshift, const string type="Poly", const string method="CAMB", const string output_root="test", const bool NL=1, const double sigmaNL=0, const int norm=-1, const double k_min=0., const double k_max=100., const double aa=0., const double prec=1.e-2, const string file_par=par::defaultString);

  };
}
#endif
