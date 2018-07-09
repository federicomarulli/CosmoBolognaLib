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
 *  @file Modelling/Global/Modelling1D.cpp
 *
 *  @brief Methods of the class Modelling1D, used for modelling any kind
 *  of measurements
 *
 *  This file contains the implementation of the methods of the class
 *  Modelling
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#include "Modelling1D.h"

using namespace std;

using namespace cbl;



// ============================================================================================


void cbl::modelling::Modelling1D::set_fit_range (const double xmin, const double xmax)
{
  m_data_fit = m_data->cut(xmin, xmax); 
  m_fit_range = true;
}

// ============================================================================================


void cbl::modelling::Modelling1D::write_model (const std::string output_dir, const std::string output_file, const std::vector<double> xx, const std::vector<double> parameters)
{
  m_posterior->write_model(output_dir, output_file, parameters, xx);
}

// ============================================================================================


void cbl::modelling::Modelling1D::write_model_at_bestfit (const std::string output_dir, const std::string output_file, const std::vector<double> xx)
{
  m_posterior->write_model_at_bestfit(output_dir, output_file, xx);
}

// ============================================================================================


void cbl::modelling::Modelling1D::write_model_from_chains (const std::string output_dir, const std::string output_file, const std::vector<double> xx, const int start, const int thin)
{
  m_posterior->write_model_from_chain(output_dir, output_file, xx, {}, start, thin);
}
