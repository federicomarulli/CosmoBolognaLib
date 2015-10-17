/*******************************************************************
 *  Copyright (C) 2010 by Federico Marulli                         *
 *  federico.marulli3@unibo.it                                     *
 *                                                                 *
 *  This program is free software; you can redistribute it and/or  *
 *  modify it under the terms of the GNU General Public License as *
 *  published by the Free Software Foundation; either version 2 of *
 *  the License, or (at your option) any later version.            *
 *                                                                 *
 *  This program is distributed in the hope that it will be useful,*
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the  *
 *  GNU General Public License for more details.                   *
 *                                                                 *
 *  You should have received a copy of the GNU General Public      *
 *  License along with this program; if not, write to the Free     *
 *  Software Foundation, Inc.,                                     *
 *  59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.      *
 *******************************************************************/

/**
 *  @file CatalogueAnalysis/ModelTwoPointCorrelation/HaloHost.cpp
 *
 *  @brief Methods of the class ModelTwoPointCorrelation used to
 *  estimate the mass of the host haloes
 *
 *  This file contains the implementation of the methods of the class
 *  ModelTwoPointCorrelation used to estimate the mass of the host
 *  dark matter haloes
 *
 *  @author Federico Marulli
 *
 *  @author federico.marulli3@unbo.it
 */

#include "ModelTwoPointCorrelation.h"
using namespace cosmobl;


// ============================================================================


double cosmobl::ModelTwoPointCorrelation::MhaloMin (Cosmology &cosm, double &lgM1_guess, double &lgM2_guess, double &mean_redshift, string &author_bias, string &author_MF, string &author_SS, string &Model, double &Delta, double &rMin, double &rMax, double &Mmax, bool proj, bool NL)
{
  // ------------------------------------------
  // ----- measure the bias, if necessary -----
  // ------------------------------------------
 
  bool measure = 0;
  if ((NL==0 && proj==0 && m_TwoP->sizeof_bias_lin_xi()==0) || 
      (NL==0 && proj==1 && m_TwoP->sizeof_bias_lin_wp()==0) || 
      (NL==1 && proj==0 && m_TwoP->sizeof_bias_nl_xi()==0) ||
      (NL==1 && proj==1 && m_TwoP->sizeof_bias_nl_wp()==0)) measure = 1;
  if (measure) m_TwoP->measure_bias(cosm, mean_redshift, author_SS, Model);

  double mean_bias = m_TwoP->mean_bias(rMin, rMax);

  double DeltaR = cosm.DeltaR(Delta, mean_redshift);

  cosmobl::classfunc::func_bias func(cosm, mean_bias, Mmax, mean_redshift, author_bias, author_MF, author_SS, Model, DeltaR);
  
  double prec = 0.0001;
  double lgM = zbrent(func, lgM1_guess, lgM2_guess, prec);

  return pow(10.,lgM);
}
