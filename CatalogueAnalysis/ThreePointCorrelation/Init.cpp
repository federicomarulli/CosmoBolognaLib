/********************************************************************
 *  Copyright (C) 2010 by Federico Marulli, Michele Moresco         *
 *  and Alfonso Veropalumbo                                         *
 *                                                                  *
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
 *  @file CatalogueAnalysis/ThreePointCorrelation/Init.cpp
 *
 *  @brief Methods of the class ThreePointCorrelation used to
 *  initialize the private members
 *
 *  This file contains the implementation of the methods of the class
 *  ThreePointCorrelation used to initialize the private members
 *
 *  @authors Federico Marulli, Michele Moresco, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, michele.moresco@unibo.it,
 *  alfonso.veropalumbo@unibo.it
 */

#include "ThreePointCorrelation.h"
using namespace cosmobl;


// ============================================================================


void cosmobl::ThreePointCorrelation::setParameters3p (double &side_s, double &side_u, double &perc_increase, string &type_binning, int &nbins)
{
  m_side_s = side_s;
  m_side_u = side_u;
  m_perc_increase = perc_increase;
  m_nbins = nbins;
  m_type_binning = type_binning;

  if (type_binning=="ang") 
    m_binsize = par::pi/nbins;
  else if (type_binning=="lin") 
    m_binsize = (side_s*(1+side_u)-side_s*(side_u-1))/nbins;
  else cosmobl::ErrorMsg("Error in binning type: choose ang or lin!");
}


// ============================================================================

 
void cosmobl::ThreePointCorrelation::allocate_vectors_zeta ()
{
  m_ggg.erase(m_ggg.begin(), m_ggg.end());
  m_rrr.erase(m_rrr.begin(), m_rrr.end());
  m_ggr.erase(m_ggr.begin(), m_ggr.end());
  m_grr.erase(m_grr.begin(), m_grr.end());
  m_zeta.erase(m_zeta.begin(), m_zeta.end());
  m_error_zeta.erase(m_error_zeta.begin(), m_error_zeta.end());
  m_zeta_red.erase(m_zeta_red.begin(), m_zeta_red.end());
  m_error_zeta_red.erase(m_error_zeta_red.begin(), m_error_zeta_red.end());

  for (int i=0; i<m_nbins+1; i++) {
    m_ggg.push_back(0.);
    m_rrr.push_back(0.);
    m_ggr.push_back(0.);
    m_grr.push_back(0.);
    m_zeta.push_back(-1.e30);
    m_error_zeta.push_back(-1.e30);
    m_zeta_red.push_back(-1.e30);
    m_error_zeta_red.push_back(-1.e30);
  }

}
