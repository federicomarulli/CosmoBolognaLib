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
 *  @file Cosmology/Lib/Bias.cpp
 *
 *  @brief Methods of the class Cosmology used to model the bias
 *
 *  This file contains the implementation of the methods of the class
 *  Cosmology used to model the bias
 *
 *  @author Federico Marulli 
 *
 *  @author federico.marulli3@unibo.it
 */

#include "Cosmology.h"

using namespace std;

using namespace cbl;


// =====================================================================================


vector<double> cbl::cosmology::Cosmology::bias_halo (const std::vector<double> Mass, const std::vector<double> Sigma, const double redshift, const std::string model_bias, const bool store_output, const std::string output_root, const std::string interpType, const double Delta, const double kk, const int norm, const double k_min, const double k_max, const double prec, const std::string method_SS, const std::string input_file, const bool is_parameter_file) 
{
  double D_N = DN(redshift);
  vector<double> bias(Mass.size());
  for (size_t i=0; i<Mass.size(); i++)
      bias[i] = m_bias_halo_generator(Sigma[i], redshift, D_N, model_bias, Delta); 

  if (m_fNL!=0)
      for (size_t i=0; i<Mass.size(); i++)
          bias[i] += bias_correction(kk, Mass[i], method_SS, store_output, output_root, interpType, norm, k_min, k_max, prec, input_file, is_parameter_file)*Sigma[i]*pow(bias[i]-1, 2); // check!!!

  return bias; 
}
