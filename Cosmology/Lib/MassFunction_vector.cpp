/********************************************************************
 *  Copyright (C) 2010 by Federico Marulli                          *
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
 *  @file Cosmology/Lib/MassFunction.cpp
 *
 *  @brief Methods of the class Cosmology used to model the mass
 *  function
 *
 *  This file contains the implementation of the methods of the class
 *  Cosmology used to model the mass function of dark matter haloes
 *
 *  @author Federico Marulli 
 *
 *  @author federico.marulli3@unibo.it
 */

#include "Cosmology.h"

using namespace std;

using namespace cbl;


// =====================================================================================


vector<double> cbl::cosmology::Cosmology::mass_function (const std::vector<double> Mass, const std::vector<double> Sigma, const std::vector<double> Dln_Sigma, const double redshift, const std::string model_MF, const bool store_output, const std::string output_root, const double Delta, const std::string interpType, const int norm, const double k_min, const double k_max, const double prec, const std::string method_SS, const std::string input_file, const bool is_parameter_file) 
{
  double fact = (m_unit) ? 1 : m_hh;
  double D_N = DN(redshift);

  vector<double> MASS(Mass.size(), 0), MF(Mass.size(), 0);
  for (size_t i=0; i<Mass.size(); i++) {
      MASS[i] = Mass[i]*fact;
      MF[i] = m_MF_generator(MASS[i], Sigma[i], Dln_Sigma[i], redshift, D_N, model_MF, Delta)*pow(fact, 4.);
  }

  if (m_fNL!=0)
      for (size_t i=0; i<Mass.size(); i++)
          MF[i] *= MF_correction(MASS[i], redshift, method_SS, store_output, output_root, interpType, norm, k_min, k_max, prec, input_file, is_parameter_file);

  return MF;
}
