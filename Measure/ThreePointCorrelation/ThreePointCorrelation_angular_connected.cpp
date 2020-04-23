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
 *  @file
 *  CatalogueAnalysis/ThreePointCorrelation/ThreePointCorrelation_angular_connected.cpp
 *
 *  @brief Methods of the class
 *  ThreePointCorrelation_angular_connected used to measure the
 *  connected three-point correlation function in angular coordinates
 *
 *  This file contains the implementation of the methods of the class
 *  ThreePointCorrelation_angular_connected used to measure the
 *  connected three-point correlation function in angular coordinates
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unibo.it, alfonso.veropalumbo@unibo.it
 */


#include "ThreePointCorrelation_angular_connected.h"

using namespace std;

using namespace cbl;
using namespace catalogue;
using namespace triplets;
using namespace measure;
using namespace threept;
using namespace glob;


// ============================================================================================


void cbl::measure::threept::ThreePointCorrelation_angular_connected::set_parameters (const double side_s, const double side_u, const double perc_increase, const int nbins) 
{
  (void)side_s; (void)side_u; (void)perc_increase; (void)nbins;
  ErrorCBL("", "set_parameters", "ThreePointCorrelation_angular_connected.cpp", ExitCode::_workInProgress_);
}


// ============================================================================================


void cbl::measure::threept::ThreePointCorrelation_angular_connected::set_parameters (const double r12, const double r12_binSize, const double r13, const double r13_binSize, const int nbins) 
{
  (void)r12; (void)r12_binSize; (void)r13; (void)r13_binSize; (void)nbins;
  ErrorCBL("", "set_parameters", "ThreePointCorrelation_angular_connected.cpp", ExitCode::_workInProgress_);
}


// ============================================================================


void cbl::measure::threept::ThreePointCorrelation_angular_connected::measure (const string dir_output_triplets, const vector<string> dir_input_triplets, const bool count_ddd, const bool count_rrr, const bool count_ddr, const bool count_drr, const bool tcount, const int seed) 
{
  (void)dir_output_triplets; (void)dir_input_triplets; (void)count_ddd; (void)count_rrr; (void)count_ddr; (void)count_drr; (void)tcount; (void)seed;
  ErrorCBL("", "measure", "ThreePointCorrelation_angular_connected.cpp", ExitCode::_workInProgress_);
}


// ============================================================================


void cbl::measure::threept::ThreePointCorrelation_angular_connected::write (const string dir, const string file) const
{
  (void)dir; (void)file;
  ErrorCBL("", "write", "ThreePointCorrelation_angular_connected.cpp", ExitCode::_workInProgress_);
}  


// ============================================================================


void cbl::measure::threept::ThreePointCorrelation_angular_connected::write_covariance (const string dir, const string file) const
{
  m_dataset->write_covariance(dir, file);
}
