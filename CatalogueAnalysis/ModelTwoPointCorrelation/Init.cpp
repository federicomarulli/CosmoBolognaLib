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
 *  @file CatalogueAnalysis/ModelTwoPointCorrelation/Init.cpp
 *
 *  @brief Methods of the class ModelTwoPointCorrelation used to
 *  initialize the private members
 *
 *  This file contains the implementation of the methods of the class
 *  ModelTwoPointCorrelation used to initialize the private members
 *
 *  @author Federico Marulli
 *
 *  @author federico.marulli3@unbo.it
 */

#include "ModelTwoPointCorrelation.h"
using namespace cosmobl;


// ============================================================================


void cosmobl::ModelTwoPointCorrelation::setLimit (double lim1, double lim2)
{
  m_lim_index_fit.erase(m_lim_index_fit.begin(),m_lim_index_fit.end());
  m_lim_index_fit.push_back(max(int(lim1/m_TwoP->linbinsz()),0));
  m_lim_index_fit.push_back(int(lim2/m_TwoP->linbinsz()));
  
  m_lim_fit.erase(m_lim_fit.begin(),m_lim_fit.end());
  m_lim_fit.push_back(lim1);
  m_lim_fit.push_back(lim2);
}


// ============================================================================


void cosmobl::ModelTwoPointCorrelation::setLimitLog (double lim1, double lim2)
{
  m_lim_index_fit.erase(m_lim_index_fit.begin(),m_lim_index_fit.end());

  int indexMin = int(ceil((log10(lim1)-m_TwoP->shift_log()-log10(m_TwoP->rMIN()))/m_TwoP->logbinsz()));
  int indexMax = int(ceil((log10(lim2)-m_TwoP->shift_log()-log10(m_TwoP->rMIN()))/m_TwoP->logbinsz()));

  m_lim_index_fit.push_back(indexMin);
  m_lim_index_fit.push_back(indexMax);
  
  m_lim_fit.erase(m_lim_fit.begin(),m_lim_fit.end());
  m_lim_fit.push_back(lim1);
  m_lim_fit.push_back(lim2);
}


// ============================================================================


void cosmobl::ModelTwoPointCorrelation::setLimit (double lim1, double lim2, double lim3, double lim4)
{
  m_lim_index_fit.erase(m_lim_index_fit.begin(),m_lim_index_fit.end());
  m_lim_index_fit.push_back(max(int(lim1/m_TwoP->linbinsz()),0));
  m_lim_index_fit.push_back(int(lim2/m_TwoP->linbinsz()));
  m_lim_index_fit.push_back(max(int(lim3/m_TwoP->linbinsz()),0));
  m_lim_index_fit.push_back(int(lim4/m_TwoP->linbinsz()));

  m_lim_fit.erase(m_lim_fit.begin(),m_lim_fit.end());
  m_lim_fit.push_back(lim1);
  m_lim_fit.push_back(lim2);
  m_lim_fit.push_back(lim3);
  m_lim_fit.push_back(lim4);
}

