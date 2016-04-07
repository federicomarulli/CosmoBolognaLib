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
 *  @file Modelling/Modelling_TwoPointCorrelation.cpp
 *
 *  @brief Methods of the class Modelling_TwoPointCorrelation,
 *  used for modelling any kind of measured two point correlation function
 *
 *  This file contains the implementation of the methods of the class
 *  Modelling_TwoPointCorrelation
 *
 *  @authors Federico Marulli
 *
 *  @authors federico.marulli3@unbo.it
 */


#include "Modelling_TwoPointCorrelation.h"
#include "Modelling_TwoPointCorrelation_monopole.h"
#include "Modelling_TwoPointCorrelation_projected.h"
#include "Modelling_TwoPointCorrelation_deprojected.h"

using namespace cosmobl;


// ============================================================================================


shared_ptr<modelling::Modelling_TwoPointCorrelation> modelling::Modelling_TwoPointCorrelation::Create (const shared_ptr<twopt::TwoPointCorrelation> twop, const double redshift, const Cosmology cosmology)
{
  if (twop->twoPType()==twopt::TwoPType::_1D_monopole_) return move(unique_ptr<Modelling_TwoPointCorrelation_monopole> (new Modelling_TwoPointCorrelation_monopole(twop, redshift, cosmology)));
  else if (twop->twoPType()==twopt::TwoPType::_1D_projected_) return move(unique_ptr<Modelling_TwoPointCorrelation_projected> (new Modelling_TwoPointCorrelation_projected(twop, redshift, cosmology)));
  else if (twop->twoPType()==twopt::TwoPType::_1D_deprojected_) return move(unique_ptr<Modelling_TwoPointCorrelation_deprojected> (new Modelling_TwoPointCorrelation_deprojected(twop, redshift, cosmology)));
  else ErrorMsg("Error in cosmobl::modelling::Modelling_TwoPointCorrelation::Create of Modelling_TwoPointCorrelation.cpp: no such type of object, or error in the input parameters!");
  return NULL;
}

