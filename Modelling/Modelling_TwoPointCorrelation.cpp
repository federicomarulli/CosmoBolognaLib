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


shared_ptr<modelling::Modelling_TwoPointCorrelation> modelling::Modelling_TwoPointCorrelation::Create (const shared_ptr<twopt::TwoPointCorrelation> twop)
{
  if (twop->twoPType()==twopt::TwoPType::_1D_monopole_) return move(unique_ptr<Modelling_TwoPointCorrelation_monopole> (new Modelling_TwoPointCorrelation_monopole(twop)));
  else if (twop->twoPType()==twopt::TwoPType::_1D_projected_) return move(unique_ptr<Modelling_TwoPointCorrelation_projected> (new Modelling_TwoPointCorrelation_projected(twop)));
  else if (twop->twoPType()==twopt::TwoPType::_1D_deprojected_) return move(unique_ptr<Modelling_TwoPointCorrelation_deprojected> (new Modelling_TwoPointCorrelation_deprojected(twop)));
  else ErrorMsg("Error in cosmobl::modelling::Modelling_TwoPointCorrelation::Create of Modelling_TwoPointCorrelation.cpp: no such type of object, or error in the input parameters!");

  return NULL;
}


// ============================================================================================


shared_ptr<modelling::Modelling_TwoPointCorrelation> modelling::Modelling_TwoPointCorrelation::Create (const shared_ptr<data::Data> twop_dataset, const twopt::TwoPType twoPType)
{
  if (twoPType==twopt::TwoPType::_1D_monopole_) 
  {
    auto mtwop = unique_ptr<Modelling_TwoPointCorrelation_monopole> (new Modelling_TwoPointCorrelation_monopole());
    mtwop->set_data(twop_dataset);
    return move(mtwop);
  }
  else if (twoPType==twopt::TwoPType::_1D_projected_) 
  {
    auto mtwop = unique_ptr<Modelling_TwoPointCorrelation_projected> (new Modelling_TwoPointCorrelation_projected());
    mtwop->set_data(twop_dataset);
    return move(mtwop);
  }
  else if (twoPType==twopt::TwoPType::_1D_deprojected_) 
  {
    auto mtwop = unique_ptr<Modelling_TwoPointCorrelation_deprojected> (new Modelling_TwoPointCorrelation_deprojected());
    mtwop->set_data(twop_dataset);
    return move(mtwop);
  }
  else ErrorMsg("Error in cosmobl::modelling::Modelling_TwoPointCorrelation::Create of Modelling_TwoPointCorrelation.cpp: no such type of object, or error in the input parameters!");

  return NULL;
}


// ============================================================================================


void modelling::Modelling_TwoPointCorrelation::set_parameters_twop_DM(const vector<double> model_scales, const cosmology::Cosmology cosmology, const double redshift, const string method, const double sigmaNL, const bool NL, const double pimax, const double r_min, const double r_max, const string output_root, const int norm, const double k_min, const double k_max, const double aa, const bool GSL, const double prec, const string file_par)
{
  m_twop_parameters.model_scales = model_scales;
  m_twop_parameters.cosmology = make_shared<cosmology::Cosmology>(cosmology);
  m_twop_parameters.redshift = redshift;
  m_twop_parameters.method = method;
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
