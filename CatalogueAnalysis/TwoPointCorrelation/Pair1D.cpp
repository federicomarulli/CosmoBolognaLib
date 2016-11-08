/*******************************************************************
 *  Copyright (C) 2015 by Federico Marulli                         *
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
 *  @file CatalogueAnalysis/TwoPointCorrelation/Pair1D.cpp
 *
 *  @brief Methods of the classes Pair1D*  
 *
 *  This file contains the implementation of all the methods of the
 *  classes Pair1D*, used to handle pairs of 1D objects of any kind
 *
 *  @author Federico Marulli 
 *
 *  @author federico.marulli3@unbo.it
 */

#include "Pair1D.h"

using namespace cosmobl;
using namespace catalogue;
using namespace pairs;


// ============================================================================================


void cosmobl::pairs::Pair1D_angular_lin::m_set_parameters_nbins ()
{
  double binSize = ((m_thetaMax-m_thetaMin)/m_nbins);
  m_binSize_inv = 1./binSize;
  
  m_scale.resize(m_nbins);
  for (int i=0; i<m_nbins; i++)
    m_scale[i] = binSize*(i+m_shift)+m_thetaMin;
}


// ============================================================================================


void cosmobl::pairs::Pair1D_angular_lin::m_set_parameters_binSize ()
{
  m_nbins = nint((m_thetaMax-m_thetaMin)*m_binSize_inv);
  m_thetaMax = m_nbins/m_binSize_inv+m_thetaMin;

  m_scale.resize(m_nbins);
  for (int i=0; i<m_nbins; i++)
    m_scale[i] = (i+m_shift)/m_binSize_inv+m_thetaMin;
}


// ============================================================================================


void cosmobl::pairs::Pair1D_angular_log::m_set_parameters_nbins ()
{
  if (m_thetaMin<1.e-30) ErrorCBL("Error in cosmobl::pairs::Pair1D_angular_log::m_set_parameters_nbins of Pair.cpp: Min must be >0!");
  
  double binSize = ((log10(m_thetaMax)-log10(m_thetaMin))/m_nbins);
  m_binSize_inv = 1./binSize;

  m_scale.resize(m_nbins);
  for (int i=0; i<m_nbins; i++)
    m_scale[i] = pow(10.,(i+m_shift)*binSize+log10(m_thetaMin));
}


// ============================================================================================


void cosmobl::pairs::Pair1D_angular_log::m_set_parameters_binSize ()
{
  if (m_thetaMin<1.e-30) ErrorCBL("Error in cosmobl::pairs::Pair1D_angular_log::m_set_parameters_binSize of Pair.cpp: Min must be >0!");

  m_nbins = nint((log10(m_thetaMax)-log10(m_thetaMin))*m_binSize_inv);
  m_thetaMax = pow(10.,m_nbins/m_binSize_inv+log10(m_thetaMin));

  m_scale.resize(m_nbins);
  for (int i=0; i<m_nbins; i++)
    m_scale[i] = pow(10.,(i+m_shift)/m_binSize_inv+log10(m_thetaMin));
}


// ============================================================================================


void cosmobl::pairs::Pair1D_comoving_lin::m_set_parameters_nbins ()
{
  double binSize = ((m_rMax-m_rMin)/m_nbins);
  m_binSize_inv = 1./binSize;

  m_scale.resize(m_nbins);
  for (int i=0; i<m_nbins; i++)
    m_scale[i] = (i+m_shift)*binSize+m_rMin;
}


// ============================================================================================


void cosmobl::pairs::Pair1D_comoving_lin::m_set_parameters_binSize ()
{
  m_nbins = nint((m_rMax-m_rMin)*m_binSize_inv);
  m_rMax = m_nbins/m_binSize_inv+m_rMin;
  
  m_scale.resize(m_nbins);
  for (int i=0; i<m_nbins; i++)
    m_scale[i] = (i+m_shift)/m_binSize_inv+m_rMin;
}


// ============================================================================================


void cosmobl::pairs::Pair1D_comoving_log::m_set_parameters_nbins ()
{
  if (m_rMin<1.e-30) ErrorCBL("Error in cosmobl::pairs::Pair1D_comoving_log::m_set_parameters_nbins of Pair.cpp: Min must be >0!");
  
  double binSize = ((log10(m_rMax)-log10(m_rMin))/m_nbins);
  m_binSize_inv = 1./binSize;

  m_scale.resize(m_nbins);
  for (int i=0; i<m_nbins; i++)
    m_scale[i] = pow(10.,(i+m_shift)*binSize+log10(m_rMin));
}


// ============================================================================================


void cosmobl::pairs::Pair1D_comoving_log::m_set_parameters_binSize ()
{
  if (m_rMin<1.e-30) ErrorCBL("Error in cosmobl::pairs::Pair1D_comoving_log::m_set_parameters_binSize of Pair.cpp: Min must be >0!");
  
  m_nbins = nint((log10(m_rMax)-log10(m_rMin))*m_binSize_inv);
  m_rMax = pow(10.,m_nbins/m_binSize_inv+log10(m_rMin));

  m_scale.resize(m_nbins);
  for (int i=0; i<m_nbins; i++)
    m_scale[i] = pow(10.,(i+m_shift)/m_binSize_inv+log10(m_rMin));
}


// ============================================================================================


void cosmobl::pairs::Pair1D_angular_lin::put (const shared_ptr<Object> obj1, const shared_ptr<Object> obj2)
{
  double dist = (m_angularUnits==_radians_) ? angular_distance(obj1->xx(), obj2->xx(), obj1->yy(), obj2->yy(), obj1->zz(), obj2->zz()) 
	     : converted_angle(angular_distance(obj1->xx(), obj2->xx(), obj1->yy(), obj2->yy(), obj1->zz(), obj2->zz()), _radians_, m_angularUnits);

  if (m_thetaMin < dist && dist < m_thetaMax) {

    int kk = max(0, min(int((dist-m_thetaMin)*m_binSize_inv), m_nbins));

    m_PP1D[kk] += obj1->weight()*obj2->weight();
    
  } 
}


// ============================================================================================


void cosmobl::pairs::Pair1D_angular_log::put (const shared_ptr<Object> obj1, const shared_ptr<Object> obj2)
{
  double dist = (m_angularUnits==_radians_) ? angular_distance(obj1->xx(), obj2->xx(), obj1->yy(), obj2->yy(), obj1->zz(), obj2->zz()) 
	     : converted_angle(angular_distance(obj1->xx(), obj2->xx(), obj1->yy(), obj2->yy(), obj1->zz(), obj2->zz()), _radians_, m_angularUnits);
  
  if (m_thetaMin < dist && dist < m_thetaMax) {

    int kk = max(0, min(int((log10(dist)-log10(m_thetaMin))*m_binSize_inv), m_nbins));

    m_PP1D[kk] += obj1->weight()*obj2->weight();
    
  }
}


// ============================================================================================


void cosmobl::pairs::Pair1D_comoving_lin::put (const shared_ptr<Object> obj1, const shared_ptr<Object> obj2)
{
  double dist = Euclidean_distance(obj1->xx(), obj2->xx(), obj1->yy(), obj2->yy(), obj1->zz(), obj2->zz()); 
  
  if (m_rMin < dist && dist < m_rMax) {

    int kk = max(0, min(int((dist-m_rMin)*m_binSize_inv), m_nbins));
    
    double angWeight = (m_angularWeight==nullptr) ? 1.
      : m_angularWeight(converted_angle(angular_distance(obj1->xx()/obj1->dc(), obj2->xx()/obj2->dc(), obj1->yy()/obj1->dc(), obj2->yy()/obj2->dc(), obj1->zz()/obj1->dc(), obj2->zz()/obj2->dc()), _radians_, m_angularUnits));
    
    m_PP1D[kk] += obj1->weight()*obj2->weight()*angWeight;
    
  }
}


// ============================================================================================


void cosmobl::pairs::Pair1D_comoving_log::put (const shared_ptr<Object> obj1, const shared_ptr<Object> obj2)
{
  double dist = Euclidean_distance(obj1->xx(), obj2->xx(), obj1->yy(), obj2->yy(), obj1->zz(), obj2->zz()); 
  
  if (m_rMin < dist && dist < m_rMax) {

    int kk = max(0, min(int((log10(dist)-log10(m_rMin))*m_binSize_inv), m_nbins));

    double angWeight = (m_angularWeight==nullptr) ? 1.
      : m_angularWeight(converted_angle(angular_distance(obj1->xx()/obj1->dc(), obj2->xx()/obj2->dc(), obj1->yy()/obj1->dc(), obj2->yy()/obj2->dc(), obj1->zz()/obj1->dc(), obj2->zz()/obj2->dc()), _radians_, m_angularUnits));

    m_PP1D[kk] += obj1->weight() * obj2->weight() * angWeight;
  }
}


// ============================================================================================


void cosmobl::pairs::Pair1D::add_data1D (const int i, const vector<double> data)
{
  checkDim(m_PP1D, i, "m_PP1D", false);
  checkDim(data, 1, "data");
  
  m_PP1D[i] += data[0];
}


// ============================================================================================


void cosmobl::pairs::Pair1D::add_data1D (const int i, const shared_ptr<pairs::Pair> pair, const double ww) 
{
  checkDim(m_PP1D, i, "m_PP1D", false);
  
  m_PP1D[i] += ww*pair->PP1D(i);
}


// ============================================================================================


void cosmobl::pairs::Pair1D::Sum (const shared_ptr<Pair> pp, const double ww)
{
  if (m_nbins != pp->nbins()) 
    ErrorCBL("Error in cosmobl::pairs::Pair1D::Sum of Pair.cpp: dimension problems!");
  
  for (int i=0; i<m_nbins; ++i) 
    m_PP1D[i] += ww*pp->PP1D(i);
}
