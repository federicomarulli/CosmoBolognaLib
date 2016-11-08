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
 *  @file Modelling/Modelling_TwoPointCorrelation.cpp
 *
 *  @brief Methods of the class Modelling_TwoPointCorrelation, used to
 *  model two-point correlation functions of any kind
 *
 *  This file contains the implementation of the methods of the class
 *  Modelling_TwoPointCorrelation
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */


#include "Modelling_TwoPointCorrelation.h"
#include "Modelling_TwoPointCorrelation_monopole.h"
#include "Modelling_TwoPointCorrelation_cartesian.h"
#include "Modelling_TwoPointCorrelation_projected.h"
#include "Modelling_TwoPointCorrelation_deprojected.h"

using namespace cosmobl;


// ============================================================================================


shared_ptr<modelling::Modelling_TwoPointCorrelation> modelling::Modelling_TwoPointCorrelation::Create (const shared_ptr<twopt::TwoPointCorrelation> twop)
{
  if (twop->twoPType()==twopt::TwoPType::_1D_monopole_)
    return move(unique_ptr<Modelling_TwoPointCorrelation_monopole> (new Modelling_TwoPointCorrelation_monopole(twop)));
  
  else if (twop->twoPType()==twopt::TwoPType::_2D_Cartesian_)
    return move(unique_ptr<Modelling_TwoPointCorrelation_cartesian> (new Modelling_TwoPointCorrelation_cartesian(twop)));

  else if (twop->twoPType()==twopt::TwoPType::_1D_projected_)
    return move(unique_ptr<Modelling_TwoPointCorrelation_projected> (new Modelling_TwoPointCorrelation_projected(twop)));

  else if (twop->twoPType()==twopt::TwoPType::_1D_deprojected_)
    return move(unique_ptr<Modelling_TwoPointCorrelation_deprojected> (new Modelling_TwoPointCorrelation_deprojected(twop)));

  else ErrorCBL("Error in cosmobl::modelling::Modelling_TwoPointCorrelation::Create of Modelling_TwoPointCorrelation.cpp: no such type of object, or error in the input parameters!");
 
  return NULL;
}


// ============================================================================================


shared_ptr<modelling::Modelling_TwoPointCorrelation> modelling::Modelling_TwoPointCorrelation::Create (const twopt::TwoPType twoPType, const shared_ptr<data::Data> twop_dataset)
{
  if (twoPType==twopt::TwoPType::_1D_monopole_)
    return move(unique_ptr<Modelling_TwoPointCorrelation_monopole> (new Modelling_TwoPointCorrelation_monopole(twop_dataset)));

  else if (twoPType==twopt::TwoPType::_2D_Cartesian_)
    return move(unique_ptr<Modelling_TwoPointCorrelation_cartesian> (new Modelling_TwoPointCorrelation_cartesian(twop_dataset)));

  else if (twoPType==twopt::TwoPType::_1D_projected_)
    return move(unique_ptr<Modelling_TwoPointCorrelation_projected> (new Modelling_TwoPointCorrelation_projected(twop_dataset)));

  else if (twoPType==twopt::TwoPType::_1D_deprojected_)
    return move(unique_ptr<Modelling_TwoPointCorrelation_deprojected> (new Modelling_TwoPointCorrelation_deprojected(twop_dataset)));
  
  else ErrorCBL("Error in cosmobl::modelling::Modelling_TwoPointCorrelation::Create of Modelling_TwoPointCorrelation.cpp: no such type of object, or error in the input parameters!");

  return NULL;
}


// ============================================================================================


void cosmobl::modelling::Modelling_TwoPointCorrelation::set_parameters_xiDM (const vector<double> fiducial_radDM, const cosmology::Cosmology cosmology, const double redshift, const string method_Pk, const double sigmaNL, const bool NL, const double pimax, const double r_min, const double r_max, const string output_root, const int norm, const double k_min, const double k_max, const double aa, const bool GSL, const double prec, const string file_par)
{
  m_twop_parameters.fiducial_radDM = fiducial_radDM;
  m_twop_parameters.cosmology = make_shared<cosmology::Cosmology>(cosmology);
  m_twop_parameters.redshift = redshift;
  m_twop_parameters.method_Pk = method_Pk;
  m_twop_parameters.sigmaNL = sigmaNL;
  m_twop_parameters.NL = NL;
  m_twop_parameters.pi_max = pimax;
  m_twop_parameters.r_min = r_min;
  m_twop_parameters.r_max = r_max;
  m_twop_parameters.output_root = output_root;
  m_twop_parameters.norm = norm;
  m_twop_parameters.k_min = k_min;
  m_twop_parameters.k_max = k_max;
  m_twop_parameters.aa = aa;
  m_twop_parameters.GSL = GSL;
  m_twop_parameters.prec = prec;
  m_twop_parameters.file_par = file_par;
}
	

// ============================================================================================


void cosmobl::modelling::Modelling_TwoPointCorrelation::set_parameters_xi2D_DM (const vector<double> fiducial_radDM, const cosmology::Cosmology cosmology, const double redshift, const string method_Pk, const double sigmaNL, const bool NL, const int FV, const vector<double> Xi, const vector<double> Xi_, const vector<double> Xi__, const string output_root, const bool bias_nl, const double bA, const bool xiType, const double k_star, const bool xiNL, const double v_min, const double v_max, const int step_v, const int norm, const double r_min, const double r_max, const double k_min, const double k_max, const double aa, const bool GSL, const double prec, const string file_par)
{
  m_twop_parameters.fiducial_radDM = fiducial_radDM;
  m_twop_parameters.cosmology = make_shared<cosmology::Cosmology>(cosmology);
  m_twop_parameters.redshift = redshift;
  m_twop_parameters.method_Pk = method_Pk;
  m_twop_parameters.sigmaNL = sigmaNL;
  m_twop_parameters.NL = NL;
  m_twop_parameters.FV = FV;
  m_twop_parameters.Xi = Xi;
  m_twop_parameters.Xi_ = Xi_;
  m_twop_parameters.Xi__ = Xi__;
  m_twop_parameters.output_root = output_root;
  m_twop_parameters.bias_nl = bias_nl;
  m_twop_parameters.bA = bA;
  m_twop_parameters.xiType = xiType;
  m_twop_parameters.k_star = k_star;
  m_twop_parameters.xiNL = xiNL;
  m_twop_parameters.v_min = v_min;
  m_twop_parameters.v_max = v_max;
  m_twop_parameters.step_v = step_v;
  m_twop_parameters.norm = norm;
  m_twop_parameters.r_min = r_min;
  m_twop_parameters.r_max = r_max;
  m_twop_parameters.k_min = k_min;
  m_twop_parameters.k_max = k_max;
  m_twop_parameters.aa = aa;
  m_twop_parameters.GSL = GSL;
  m_twop_parameters.prec = prec;
  m_twop_parameters.file_par = file_par;
}
