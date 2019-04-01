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
 *  Modelling/TwoPointCorrelation/Modelling_TwoPointCorrelation.cpp
 *
 *  @brief Methods of the class Modelling_TwoPointCorrelation
 *
 *  This file contains the implementation of the methods of the class
 *  Modelling_TwoPointCorrelation, i.e. the common functions to model
 *  the two-point correlation functions of any kind
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */


#include "Modelling_TwoPointCorrelation.h"
#include "Modelling_TwoPointCorrelation1D_monopole.h"
#include "Modelling_TwoPointCorrelation2D_cartesian.h"
#include "Modelling_TwoPointCorrelation_projected.h"
#include "Modelling_TwoPointCorrelation_deprojected.h"

using namespace std;

using namespace cbl;
using namespace modelling;
using namespace modelling::twopt;

// ============================================================================================


std::shared_ptr<Modelling_TwoPointCorrelation> modelling::twopt::Modelling_TwoPointCorrelation::Create (const std::shared_ptr<measure::twopt::TwoPointCorrelation> twop)
{
  if (twop->twoPType()==measure::twopt::TwoPType::_monopole_)
    return move(unique_ptr<Modelling_TwoPointCorrelation1D_monopole> (new Modelling_TwoPointCorrelation1D_monopole(twop)));
  
  else if (twop->twoPType()==measure::twopt::TwoPType::_2D_Cartesian_)
    return move(unique_ptr<Modelling_TwoPointCorrelation2D_cartesian> (new Modelling_TwoPointCorrelation2D_cartesian(twop)));

  else if (twop->twoPType()==measure::twopt::TwoPType::_projected_)
    return move(unique_ptr<Modelling_TwoPointCorrelation_projected> (new Modelling_TwoPointCorrelation_projected(twop)));

  else if (twop->twoPType()==measure::twopt::TwoPType::_deprojected_)
    return move(unique_ptr<Modelling_TwoPointCorrelation_deprojected> (new Modelling_TwoPointCorrelation_deprojected(twop)));

  else ErrorCBL("Error in cbl::modelling::twopt::Modelling_TwoPointCorrelation::Create of Modelling_TwoPointCorrelation.cpp: no such type of object, or error in the input parameters!");
 
  return NULL;
}


// ============================================================================================


std::shared_ptr<Modelling_TwoPointCorrelation> modelling::twopt::Modelling_TwoPointCorrelation::Create (const measure::twopt::TwoPType twoPType, const std::shared_ptr<data::Data> twop_dataset)
{
  if (twoPType==measure::twopt::TwoPType::_monopole_)
    return move(unique_ptr<Modelling_TwoPointCorrelation1D_monopole> (new Modelling_TwoPointCorrelation1D_monopole(twop_dataset)));

  else if (twoPType==measure::twopt::TwoPType::_2D_Cartesian_)
    return move(unique_ptr<Modelling_TwoPointCorrelation2D_cartesian> (new Modelling_TwoPointCorrelation2D_cartesian(twop_dataset)));

  else if (twoPType==measure::twopt::TwoPType::_projected_)
    return move(unique_ptr<Modelling_TwoPointCorrelation_projected> (new Modelling_TwoPointCorrelation_projected(twop_dataset)));

  else if (twoPType==measure::twopt::TwoPType::_deprojected_)
    return move(unique_ptr<Modelling_TwoPointCorrelation_deprojected> (new Modelling_TwoPointCorrelation_deprojected(twop_dataset)));
  
  else ErrorCBL("Error in cbl::modelling::twopt::Modelling_TwoPointCorrelation::Create of Modelling_TwoPointCorrelation.cpp: no such type of object, or error in the input parameters!");

  return NULL;
}



