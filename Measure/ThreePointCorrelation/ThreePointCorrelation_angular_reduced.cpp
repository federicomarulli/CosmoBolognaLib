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
 *  Measure/ThreePointCorrelation/ThreePointCorrelation_angular_reduced.cpp
 *
 *  @brief Methods of the class ThreePointCorrelation_angular_reduced
 *  used to measure the reduced three-point correlation function in
 *  angular coordinates
 *
 *  This file contains the implementation of the methods of the class
 *  ThreePointCorrelation_angular_reduced used to measure the recuded
 *  three-point correlation function in angular coordinates
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */


#include "ThreePointCorrelation_angular_reduced.h"

using namespace cosmobl;
using namespace catalogue;
using namespace triplets;
using namespace measure;
using namespace threept;
using namespace glob;


// ============================================================================


void cosmobl::measure::threept::ThreePointCorrelation_angular_reduced::measure (const string dir_output_triplets, const string dir_output_2pt, const vector<string> dir_input_triplets, const bool count_ddd, const bool count_rrr, const bool count_ddr, const bool count_drr, const bool tcount) 
{
  (void)dir_output_triplets; (void)dir_output_2pt; (void)dir_input_triplets; (void)count_ddd; (void)count_rrr; (void)count_ddr; (void)count_drr; (void)tcount; 
  ErrorCBL("Error in threept::ThreePointCorrelation_angular_reduced::measure() of ThreePointCorrelation_angular_reduced.cpp", ExitCode::_workInProgress_);
}


// ============================================================================


void cosmobl::measure::threept::ThreePointCorrelation_angular_reduced::write (const string dir, const string file, const bool connected) const
{      
  (void)dir; (void)file; (void)connected;
  ErrorCBL("Error in threept::ThreePointCorrelation_angular_reduced::write() of ThreePointCorrelation_angular_reduced.cpp", ExitCode::_workInProgress_);
}  
