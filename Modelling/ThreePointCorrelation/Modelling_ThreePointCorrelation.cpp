/********************************************************************
 *  Copyright (C) 2017 by Federico Marulli                          *
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
 *  Modelling/ThreePointCorrelation/Modelling_ThreePointCorrelation.cpp
 *
 *  @brief Methods of the class Modelling_ThreePointCorrelation, used to
 *  model three-point correlation functions of any kind
 *
 *  This file contains the implementation of the methods of the class
 *  Modelling_ThreePointCorrelation
 *
 *  @authors Federico Marulli, Michele Moresco
 *
 *  @authors federico.marulli3@unbo.it, michele.moresco@unibo.it
 */


#include "Modelling_ThreePointCorrelation.h"
#include "Modelling_ThreePointCorrelation_angular_connected.h"
#include "Modelling_ThreePointCorrelation_angular_reduced.h"
#include "Modelling_ThreePointCorrelation_comoving_connected.h"
#include "Modelling_ThreePointCorrelation_comoving_reduced.h"

using namespace cosmobl;
using namespace measure::threept;
using namespace modelling::threept;


// ============================================================================================


shared_ptr<modelling::threept::Modelling_ThreePointCorrelation> modelling::threept::Modelling_ThreePointCorrelation::Create (const shared_ptr<measure::threept::ThreePointCorrelation> threept)
{
  if (threept->threePType()==measure::threept::ThreePType::_angular_connected_)
    return move(unique_ptr<Modelling_ThreePointCorrelation_angular_connected> (new Modelling_ThreePointCorrelation_angular_connected(threept)));
  
  else if (threept->threePType()==measure::threept::ThreePType::_angular_reduced_)
    return move(unique_ptr<modelling::threept::Modelling_ThreePointCorrelation_angular_reduced> (new Modelling_ThreePointCorrelation_angular_reduced(threept)));

  else if (threept->threePType()==measure::threept::ThreePType::_comoving_connected_)
    return move(unique_ptr<Modelling_ThreePointCorrelation_comoving_connected> (new Modelling_ThreePointCorrelation_comoving_connected(threept)));

  else if (threept->threePType()==measure::threept::ThreePType::_comoving_reduced_)
    return move(unique_ptr<Modelling_ThreePointCorrelation_comoving_reduced> (new Modelling_ThreePointCorrelation_comoving_reduced(threept)));

  else ErrorCBL("Error in cosmobl::modelling::threept::Modelling_ThreePointCorrelation::Create of Modelling_ThreePointCorrelation.cpp: no such type of object, or error in the input parameters!");
 
  return NULL;
}


// ============================================================================================


shared_ptr<modelling::threept::Modelling_ThreePointCorrelation> modelling::threept::Modelling_ThreePointCorrelation::Create (const measure::threept::ThreePType threePType, const shared_ptr<data::Data> threept_dataset)
{
  if (threePType==measure::threept::ThreePType::_angular_connected_)
    return move(unique_ptr<Modelling_ThreePointCorrelation_angular_connected> (new Modelling_ThreePointCorrelation_angular_connected(threept_dataset)));

  else if (threePType==measure::threept::ThreePType::_angular_reduced_)
    return move(unique_ptr<Modelling_ThreePointCorrelation_angular_reduced> (new Modelling_ThreePointCorrelation_angular_reduced(threept_dataset)));

  else if (threePType==measure::threept::ThreePType::_comoving_connected_)
    return move(unique_ptr<Modelling_ThreePointCorrelation_comoving_connected> (new Modelling_ThreePointCorrelation_comoving_connected(threept_dataset)));

  else if (threePType==measure::threept::ThreePType::_comoving_reduced_)
    return move(unique_ptr<Modelling_ThreePointCorrelation_comoving_reduced> (new Modelling_ThreePointCorrelation_comoving_reduced(threept_dataset)));
  
  else ErrorCBL("Error in cosmobl::modelling::threept::Modelling_ThreePointCorrelation::Create of Modelling_ThreePointCorrelation.cpp: no such type of object, or error in the input parameters!");

  return NULL;
}

// ============================================================================================


void cosmobl::modelling::threept::Modelling_ThreePointCorrelation::set_data_model (const vector<double> Q_DM)
{
  m_data_model.Q_DM = Q_DM;
}


// ============================================================================================


void cosmobl::modelling::threept::Modelling_ThreePointCorrelation::set_data_Q_nonlocal (const cosmology::Cosmology cosmology, const double r1, const double r2, const vector<double> theta, const string model, const vector<double> kk, const vector<double> Pk_DM)
{
  m_data_model.cosmology = make_shared<cosmology::Cosmology>(cosmology);
  m_data_model.r1 = r1;
  m_data_model.r2 = r2;
  m_data_model.theta = theta;
  m_data_model.model = model;
  m_data_model.kk = kk;
  m_data_model.Pk_DM = Pk_DM;
}

