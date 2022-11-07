/********************************************************************
 *  Copyright (C) 2015 by Federico Marulli, Carlo Giocoli           *
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
 *  @file TwoPointCorrelation/TwoPointCorrelationCross.cpp
 *
 *  @brief Methods of the class TwoPointCorrelationCross
 *
 *  This file contains the implementation of the methods of the class
 *  TwoPointCorrelationCross, used to measure the cross two-point
 *  correlation function
 *
 *  @author Federico Marulli, Carlo Giocoli
 *
 *  @author federico.marulli3@unibo.it, carlo.giocoli@unibo.it
 */

#include "TwoPointCorrelationCross.h"

using namespace std;

using namespace cbl;
using namespace catalogue;
using namespace chainmesh;
using namespace data; 
using namespace pairs;
using namespace measure;
using namespace twopt;


// ============================================================================


void cbl::measure::twopt::TwoPointCorrelationCross::count_allPairs (const TwoPType type, const string dir_output_pairs, const vector<string> dir_input_pairs, const bool count_d1d2, const bool count_rr, const bool count_d1r, const bool count_d2r, const bool tcount, const Estimator estimator, const double fact)  
{
  // ----------- compute polar coordinates, if necessary ----------- 
 
  if (!m_data->isSetVar(Var::_RA_) || !m_data->isSetVar(Var::_Dec_) || !m_data->isSetVar(Var::_Dc_)) 
    m_data->computePolarCoordinates();
  
  if (!m_data2->isSetVar(Var::_RA_) || !m_data2->isSetVar(Var::_Dec_) || !m_data2->isSetVar(Var::_Dc_)) 
    m_data2->computePolarCoordinates();
  
  if (!m_random->isSetVar(Var::_RA_) || !m_random->isSetVar(Var::_Dec_) || !m_random->isSetVar(Var::_Dc_)) 
    m_random->computePolarCoordinates();
  
  if (type==TwoPType::_angular_) {
    m_data->normalizeComovingCoordinates();
    m_data2->normalizeComovingCoordinates();
    m_random->normalizeComovingCoordinates();
  }

  
  // ----------- dilute the random catalogue used to compute the RR pairs (to improve the performance) ----------- 
  
  if (estimator==Estimator::_natural_ && m_random_dilution_fraction!=1.) {
    m_random_dilution_fraction = 1.;
    WarningMsgCBL("m_random_dilution_fraction = 1, since the random catalogue is not diluted when using the natural estimator!", "count_allPairs", "TwoPointCorrelation.cpp");
  }

  auto random_dil = make_shared<catalogue::Catalogue>(catalogue::Catalogue(move(m_random->diluted_catalogue(m_random_dilution_fraction))));
  
  
  // ----------- create the chain-mesh ----------- 

  double rMAX;
  double rMIN;

  if (type==TwoPType::_monopole_ || type==TwoPType::_filtered_ || type==TwoPType::_multipoles_direct_) {
    rMAX = m_d1d2->sMax();
    rMIN = m_d1d2->sMin();
  }

  else if (type==TwoPType::_angular_) {
    double xx, yy, zz;
    cartesian_coord(radians(m_d1d2->sMax(), m_d1d2->angularUnits()), radians(m_d1d2->sMax(), m_d1d2->angularUnits()), 1., xx, yy, zz);
    rMAX = sqrt(xx*xx+zz*zz);
    cartesian_coord(radians(m_d1d2->sMin(), m_d1d2->angularUnits()), radians(m_d1d2->sMin(), m_d1d2->angularUnits()), 1., xx, yy, zz);
    rMIN = sqrt(xx*xx+zz*zz);
  }

  else if (type==TwoPType::_2D_polar_ || type==TwoPType::_multipoles_integrated_ || type ==TwoPType::_wedges_) {
    rMAX = m_d1d2->sMax_D1();
    rMIN = m_d1d2->sMin_D1();
  }
  
  else if (type==TwoPType::_2D_Cartesian_ || type==TwoPType::_projected_ || type==TwoPType::_deprojected_) {
    rMAX = max(m_d1d2->sMax_D1(), m_d1d2->sMax_D2())*sqrt(2.);
    rMIN = max(m_d1d2->sMin_D1(), m_d1d2->sMin_D2())*sqrt(2.);
  }
  
  else
    ErrorCBL("the chosen two-point correlation function type is uknown!", "count_allPairs", "TwoPointCorrelationCross.cpp");
  
  double cell_size = fact*rMAX; // to be optimized!!!
  
  ChainMesh_Catalogue ChM_data1data2, ChM_random_dil, ChM_data1random, ChM_data2random;
  
  if (count_d1d2)
    ChM_data1data2.set_par(cell_size, m_data2, rMAX, rMIN);

  if (count_rr)
    ChM_random_dil.set_par(cell_size, random_dil, rMAX, rMIN);    

  if (count_d1r)
    ChM_data1random.set_par(cell_size, m_random, rMAX, rMIN);
  
  if (count_d2r)
    ChM_data2random.set_par(cell_size, m_random, rMAX, rMIN);
  
  
  // ----------- count the number of pairs or read them from file -----------

  string file;
  
  cout << endl; coutCBL << par::col_green << "data1-data2" << par::col_default << endl;
  file = "d1d2.dat";
 
  if (count_d1d2) {
    count_pairs(m_data, ChM_data1data2, m_d1d2, true, tcount);
    if (dir_output_pairs!=par::defaultString) write_pairs(m_d1d2, dir_output_pairs, file);
  }
  else read_pairs(m_d1d2, dir_input_pairs, file);
  
  
  cout << endl; coutCBL << par::col_green << "random-random" << par::col_default << endl;
  file = "rr.dat";
 
  if (count_rr) {
    count_pairs(random_dil, ChM_random_dil, m_rr, false, tcount);
    if (dir_output_pairs!=par::defaultString) write_pairs(m_rr, dir_output_pairs, file);
  }
  else read_pairs(m_rr, dir_input_pairs, file);
  
    
  cout << endl; coutCBL << par::col_green << "data1-random" << par::col_default << endl; 
  file = "d1r.dat";
    
  if (count_d1r) {
    count_pairs(m_data, ChM_data1random, m_d1r, true, tcount);
    if (dir_output_pairs!=par::defaultString) write_pairs(m_d1r, dir_output_pairs, file);
  }
  else read_pairs(m_d1r, dir_input_pairs, file);

  
  cout << endl; coutCBL << par::col_green << "data2-random" << par::col_default << endl; 
  file = "d2r.dat";
    
  if (count_d2r) {
    count_pairs(m_data2, ChM_data2random, m_d2r, true, tcount);
    if (dir_output_pairs!=par::defaultString) write_pairs(m_d2r, dir_output_pairs, file);
  }
  else read_pairs(m_d2r, dir_input_pairs, file);

 
  if (count_d1d2) 
    m_data2->Order();
  
  if (count_rr || count_d1r || count_d2r)
    m_random->Order();

  if (type==TwoPType::_angular_) {
    m_data->restoreComovingCoordinates();
    m_data2->restoreComovingCoordinates();
    m_random->restoreComovingCoordinates();
  }
  
}


// ============================================================================


double TwoPointCorrelationCross::PoissonError (const Estimator estimator, const double d1d2, const double rr, const double d1r, const double d2r, const int nData1, const int nData2, const int nRandom) const
{
  (void)estimator; (void)d1d2; (void)rr;  (void)d1r; (void)d2r; (void)nData1; (void)nData2; (void)nRandom; 
  
  WarningMsgCBL("the PoissonError has still to be implemented for the cross-correlation function", "PoissonError", "TwoPointCorrelation.cpp");

  return -1000.;
}
