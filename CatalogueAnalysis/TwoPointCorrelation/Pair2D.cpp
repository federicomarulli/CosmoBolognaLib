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
 *  @file CatalogueAnalysis/TwoPointCorrelation/Pair2D.cpp
 *
 *  @brief Methods of the classes Pair2D*  
 *
 *  This file contains the implementation of all the methods of the
 *  classes Pair2D*, used to handle 2D pairs of objects of any kind
 *
 *  @author Federico Marulli 
 *
 *  @author federico.marulli3@unbo.it
 */

#include "Pair2D.h"

using namespace cosmobl;
using namespace catalogue;
using namespace pairs;


// ============================================================================================


void cosmobl::pairs::Pair2D_comovingCartesian_linlin::m_set_parameters_nbins ()
{
  const double binSize_D1 = ((m_rpMax-m_rpMin)/m_nbins_D1);
  m_binSize_inv_D1 = 1./binSize_D1;

  const double binSize_D2 = ((m_piMax-m_piMin)/m_nbins_D2);
  m_binSize_inv_D2 = 1./binSize_D2;

  m_scale_D1.resize(m_nbins_D1); m_scale_D2.resize(m_nbins_D2);
  for (int i=0; i<m_nbins_D1; i++) 
    m_scale_D1[i] = (i+m_shift_D1)*binSize_D1+m_rpMin;
  for (int i=0; i<m_nbins_D2; i++) 
    m_scale_D2[i] = (i+m_shift_D2)*binSize_D2+m_piMin;
}


// ============================================================================================


void cosmobl::pairs::Pair2D_comovingCartesian_linlin::m_set_parameters_binSize ()
{
  m_nbins_D1 = nint((m_rpMax-m_rpMin)*m_binSize_inv_D1);
  m_rpMax = m_nbins_D1/m_binSize_inv_D1+m_rpMin;

  m_nbins_D2 = nint((m_piMax-m_piMin)*m_binSize_inv_D2);
  m_piMax = m_nbins_D2/m_binSize_inv_D2+m_piMin;

  m_scale_D1.resize(m_nbins_D1); m_scale_D2.resize(m_nbins_D2);
  for (int i=0; i<m_nbins_D1; i++) 
    m_scale_D1[i] = (i+m_shift_D1)/m_binSize_inv_D1+m_rpMin;
  for (int i=0; i<m_nbins_D2; i++) 
    m_scale_D2[i] = (i+m_shift_D2)/m_binSize_inv_D2+m_piMin;
}


// ============================================================================================


void cosmobl::pairs::Pair2D_comovingCartesian_linlog::m_set_parameters_nbins ()
{
  if (m_piMin<1.e-30)
    ErrorCBL("Error in cosmobl::pairs::Pair2D_comovingCartesian_linlog::m_set_parameters_nbins of Pair.cpp: m_piMin must be >0!");
  
  const double binSize_D1 = ((m_rpMax-m_rpMin)/m_nbins_D1);
  m_binSize_inv_D1 = 1./binSize_D1;

  const double binSize_D2 = ((log10(m_piMax)-log10(m_piMin))/m_nbins_D2);
  m_binSize_inv_D2 = 1./binSize_D2;

  m_scale_D1.resize(m_nbins_D1); m_scale_D2.resize(m_nbins_D2);
  for (int i=0; i<m_nbins_D1; i++) 
    m_scale_D1[i] = (i+m_shift_D1)*binSize_D1+m_rpMin;
  for (int i=0; i<m_nbins_D2; i++) 
    m_scale_D2[i] = pow(10.,(i+m_shift_D2)*binSize_D2+log10(m_piMin));
}


// ============================================================================================


void cosmobl::pairs::Pair2D_comovingCartesian_linlog::m_set_parameters_binSize ()
{
  if (m_piMin<1.e-30)
    ErrorCBL("Error in cosmobl::pairs::Pair2D_comovingCartesian_linlog::m_set_parameters_binSize of Pair.cpp: m_piMin must be >0!");
  
  m_nbins_D1 = nint((m_rpMax-m_rpMin)*m_binSize_inv_D1);
  m_rpMax = m_nbins_D1/m_binSize_inv_D1+m_rpMin;

  m_nbins_D2 = nint((log10(m_piMax)-log10(m_piMin))*m_binSize_inv_D2);
  m_piMax = pow(10.,m_nbins_D2/m_binSize_inv_D2+log10(m_piMin));

  m_scale_D1.resize(m_nbins_D1); m_scale_D2.resize(m_nbins_D2);
  for (int i=0; i<m_nbins_D1; i++) 
    m_scale_D1[i] = (i+m_shift_D1)/m_binSize_inv_D1+m_rpMin;
  for (int i=0; i<m_nbins_D2; i++) 
    m_scale_D2[i] = pow(10.,(i+m_shift_D2)/m_binSize_inv_D2+log10(m_piMin));
}


// ============================================================================================


void cosmobl::pairs::Pair2D_comovingCartesian_loglin::m_set_parameters_nbins ()
{
  if (m_rpMin<1.e-30)
    ErrorCBL("Error in cosmobl::pairs::Pair2D_comovingCartesian_loglin::m_set_parameters_nbins of Pair.cpp: m_rpMin must be >0!");
  
  const double binSize_D1 = ((log10(m_rpMax)-log10(m_rpMin))/m_nbins_D1);
  m_binSize_inv_D1 = 1./binSize_D1;
  
  const double binSize_D2 = ((m_piMax-m_piMin)/m_nbins_D2);
  m_binSize_inv_D2 = 1./binSize_D2;

  m_scale_D1.resize(m_nbins_D1); m_scale_D2.resize(m_nbins_D2);
  for (int i=0; i<m_nbins_D1; i++) 
    m_scale_D1[i] = pow(10.,(i+m_shift_D1)*binSize_D1+log10(m_rpMin)); 
  for (int i=0; i<m_nbins_D2; i++) 
    m_scale_D2[i] = (i+m_shift_D2)*binSize_D2+m_piMin;
}


// ============================================================================================


void cosmobl::pairs::Pair2D_comovingCartesian_loglin::m_set_parameters_binSize ()
{
  if (m_rpMin<1.e-30)
    ErrorCBL("Error in cosmobl::pairs::Pair2D_comovingCartesian_linlin::m_set_parameters_binSize of Pair.cpp: m_rpMin must be >0!");
  
  m_nbins_D1 = nint((log10(m_rpMax)-log10(m_rpMin))*m_binSize_inv_D1);
  m_rpMax = pow(10.,m_nbins_D1/m_binSize_inv_D1+log10(m_rpMin));
  
  m_nbins_D2 = nint((m_piMax-m_piMin)*m_binSize_inv_D2);
  m_piMax = m_nbins_D2/m_binSize_inv_D2+m_piMin;

  m_scale_D1.resize(m_nbins_D1); m_scale_D2.resize(m_nbins_D2);
  for (int i=0; i<m_nbins_D1; i++) 
    m_scale_D1[i] = pow(10.,(i+m_shift_D1)/m_binSize_inv_D1+log10(m_rpMin));
  for (int i=0; i<m_nbins_D2; i++) 
    m_scale_D2[i] = (i+m_shift_D2)/m_binSize_inv_D2+m_piMin;
}


// ============================================================================================


void cosmobl::pairs::Pair2D_comovingCartesian_loglog::m_set_parameters_nbins ()
{
  if (m_rpMin<1.e-30 || m_piMin<1.e-30)
    ErrorCBL("Error in cosmobl::pairs::Pair2D_comovingCartesian_loglog::m_set_parameters_nbins of Pair.cpp: m_rpMin and m_piMin must be >0!");
  
  const double binSize_D1 = ((log10(m_rpMax)-log10(m_rpMin))/m_nbins_D1);
  m_binSize_inv_D1 = 1./binSize_D1;
  
  const double binSize_D2 = ((log10(m_piMax)-log10(m_piMin))/m_nbins_D2);
  m_binSize_inv_D2 = 1./binSize_D2;

  m_scale_D1.resize(m_nbins_D1); m_scale_D2.resize(m_nbins_D2);
  for (int i=0; i<m_nbins_D1; i++) 
    m_scale_D1[i] = pow(10.,(i+m_shift_D1)*binSize_D1+log10(m_rpMin)); 
  for (int i=0; i<m_nbins_D2; i++) 
    m_scale_D2[i] = pow(10.,(i+m_shift_D2)*binSize_D2+log10(m_piMin)); 
}


// ============================================================================================


void cosmobl::pairs::Pair2D_comovingCartesian_loglog::m_set_parameters_binSize ()
{
  if (m_rpMin<1.e-30 || m_piMin<1.e-30)
    ErrorCBL("Error in cosmobl::pairs::Pair2D_comovingCartesian_loglog::m_set_parameters_binSize of Pair.cpp: m_rpMin and m_piMin must be >0!");
  
  m_nbins_D1 = nint((log10(m_rpMax)-log10(m_rpMin))*m_binSize_inv_D1);
  m_rpMax = pow(10.,m_nbins_D1/m_binSize_inv_D1+log10(m_rpMin));
  
  m_nbins_D2 = nint((log10(m_piMax)-log10(m_piMin))*m_binSize_inv_D2);
  m_piMax = pow(10.,m_nbins_D2/m_binSize_inv_D2+log10(m_piMin));

  m_scale_D1.resize(m_nbins_D1); m_scale_D2.resize(m_nbins_D2);
  for (int i=0; i<m_nbins_D1; i++) 
    m_scale_D1[i] = pow(10.,(i+m_shift_D1)/m_binSize_inv_D1+log10(m_rpMin));
  for (int i=0; i<m_nbins_D2; i++) 
    m_scale_D2[i] = pow(10.,(i+m_shift_D2)/m_binSize_inv_D2+log10(m_piMin));
}


// ============================================================================================


void cosmobl::pairs::Pair2D_comovingPolar_linlin::m_set_parameters_nbins ()
{ 
  const double binSize_D1 = ((m_rMax-m_rMin)/m_nbins_D1);
  m_binSize_inv_D1 = 1./binSize_D1;

  const double binSize_D2 = ((m_muMax-m_muMin)/m_nbins_D2);
  m_binSize_inv_D2 = 1./binSize_D2;

  m_scale_D1.resize(m_nbins_D1); m_scale_D2.resize(m_nbins_D2);
  for (int i=0; i<m_nbins_D1; i++) 
    m_scale_D1[i] = (i+m_shift_D1)*binSize_D1+m_rMin;
  for (int i=0; i<m_nbins_D2; i++) 
    m_scale_D2[i] = (i+m_shift_D2)*binSize_D2+m_muMin;
}


// ============================================================================================


void cosmobl::pairs::Pair2D_comovingPolar_linlin::m_set_parameters_binSize ()
{
  m_nbins_D1 = nint((m_rMax-m_rMin)*m_binSize_inv_D1);
  m_rMax = m_nbins_D1/m_binSize_inv_D1+m_rMin;

  m_nbins_D2 = nint((m_muMax-m_muMin)*m_binSize_inv_D2);
  m_muMax = m_nbins_D2/m_binSize_inv_D2+m_muMin;

  m_scale_D1.resize(m_nbins_D1); m_scale_D2.resize(m_nbins_D2);
  for (int i=0; i<m_nbins_D1; i++) 
    m_scale_D1[i] = (i+m_shift_D1)/m_binSize_inv_D1+m_rMin;
  for (int i=0; i<m_nbins_D2; i++) 
    m_scale_D2[i] = (i+m_shift_D2)/m_binSize_inv_D2+m_muMin;
}


// ============================================================================================


void cosmobl::pairs::Pair2D_comovingPolar_linlog::m_set_parameters_nbins ()
{
  if (m_muMin<1.e-30)
    ErrorCBL("Error in cosmobl::pairs::Pair2D_comovingPolar_linlog::m_set_parameters_nbins of Pair.cpp: m_muMin must be >0!");
  
  const double binSize_D1 = ((m_rMax-m_rMin)/m_nbins_D1);
  m_binSize_inv_D1 = 1./binSize_D1;

  const double binSize_D2 = ((log10(m_muMax)-log10(m_muMin))/m_nbins_D2);
  m_binSize_inv_D2 = 1./binSize_D2;

  m_scale_D1.resize(m_nbins_D1); m_scale_D2.resize(m_nbins_D2);
  for (int i=0; i<m_nbins_D1; i++) 
    m_scale_D1[i] = (i+m_shift_D1)*binSize_D1+m_rMin;
  for (int i=0; i<m_nbins_D2; i++) 
    m_scale_D2[i] = pow(10.,(i+m_shift_D2)*binSize_D2+log10(m_muMin));
}


// ============================================================================================


void cosmobl::pairs::Pair2D_comovingPolar_linlog::m_set_parameters_binSize ()
{
  if (m_muMin<1.e-30)
    ErrorCBL("Error in cosmobl::pairs::Pair2D_comovingPolar_linlog::m_set_parameters_linlog of Pair.cpp: m_muMin must be >0!");
  
  m_nbins_D1 = nint((m_rMax-m_rMin)*m_binSize_inv_D1);
  m_rMax = m_nbins_D1/m_binSize_inv_D1+m_rMin;

  m_nbins_D2 = nint((log10(m_muMax)-log10(m_muMin))*m_binSize_inv_D2);
  m_muMax = pow(10.,(m_nbins_D2-m_shift_D2)/m_binSize_inv_D2+log10(m_muMin));

  m_scale_D1.resize(m_nbins_D1); m_scale_D2.resize(m_nbins_D2);
  for (int i=0; i<m_nbins_D1; i++) 
    m_scale_D1[i] = (i+m_shift_D1)/m_binSize_inv_D1+m_rMin;
  for (int i=0; i<m_nbins_D2; i++) 
    m_scale_D2[i] = pow(10.,(i+m_shift_D2)/m_binSize_inv_D2+log10(m_muMin));
}


// ============================================================================================


void cosmobl::pairs::Pair2D_comovingPolar_loglin::m_set_parameters_nbins ()
{
  if (m_rMin<1.e-30)
    ErrorCBL("Error in cosmobl::pairs::Pair2D_comovingPolar_loglin::m_set_parameters_nbins of Pair.cpp: m_rMin must be >0!");
  
  const double binSize_D1 = ((log10(m_rMax)-log10(m_rMin))/m_nbins_D1);
  m_binSize_inv_D1 = 1./binSize_D1;
  
  const double binSize_D2 = ((m_muMax-m_muMin)/m_nbins_D2);
  m_binSize_inv_D2 = 1./binSize_D2;

  m_scale_D1.resize(m_nbins_D1); m_scale_D2.resize(m_nbins_D2);
  for (int i=0; i<m_nbins_D1; i++) 
    m_scale_D1[i] = pow(10.,(i+m_shift_D1)*binSize_D1+log10(m_rMin)); 
  for (int i=0; i<m_nbins_D2; i++) 
    m_scale_D2[i] = (i+m_shift_D2)*binSize_D2+m_muMin;
}


// ============================================================================================


void cosmobl::pairs::Pair2D_comovingPolar_loglin::m_set_parameters_binSize ()
{
  if (m_rMin<1.e-30)
    ErrorCBL("Error in cosmobl::pairs::Pair2D_comovingPolar_loglin::m_set_parameters_binSize of Pair.cpp: m_rMin must be >0!");
  
  m_nbins_D1 = nint((log10(m_rMax)-log10(m_rMin))*m_binSize_inv_D1);
  m_rMax = pow(10.,m_nbins_D1/m_binSize_inv_D1+log10(m_rMin));
  
  m_nbins_D2 = nint((m_muMax-m_muMin)*m_binSize_inv_D2);
  m_muMax = m_nbins_D2/m_binSize_inv_D2+m_muMin;

  m_scale_D1.resize(m_nbins_D1); m_scale_D2.resize(m_nbins_D2);
  for (int i=0; i<m_nbins_D1; i++) 
    m_scale_D1[i] = pow(10.,(i+m_shift_D1)/m_binSize_inv_D1+log10(m_rMin));
  for (int i=0; i<m_nbins_D2; i++) 
    m_scale_D2[i] = (i+m_shift_D2)/m_binSize_inv_D2+m_muMin;
}


// ============================================================================================


void cosmobl::pairs::Pair2D_comovingPolar_loglog::m_set_parameters_nbins ()
{
  if (m_rMin<1.e-30 || m_muMin<1.e-30)
    ErrorCBL("Error in cosmobl::pairs::Pair2D_comovingPolar_loglog::m_set_parameters_nbins of Pair.cpp: m_rMin and m_muMin must be >0!");
  
  const double binSize_D1 = ((log10(m_rMax)-log10(m_rMin))/m_nbins_D1);
  m_binSize_inv_D1 = 1./binSize_D1;
  
  const double binSize_D2 = ((log10(m_muMax)-log10(m_muMin))/m_nbins_D2);
  m_binSize_inv_D2 = 1./binSize_D2;

  m_scale_D1.resize(m_nbins_D1); m_scale_D2.resize(m_nbins_D2);
  for (int i=0; i<m_nbins_D1; i++) 
    m_scale_D1[i] = pow(10.,(i+m_shift_D1)*binSize_D1+log10(m_rMin)); 
  for (int i=0; i<m_nbins_D2; i++) 
    m_scale_D2[i] = pow(10.,(i+m_shift_D2)*binSize_D2+log10(m_muMin)); 
}


// ============================================================================================


void cosmobl::pairs::Pair2D_comovingPolar_loglog::m_set_parameters_binSize ()
{
  if (m_rMin<1.e-30 || m_muMin<1.e-30)
    ErrorCBL("Error in cosmobl::pairs::Pair2D_comovingPolar_loglog::m_set_parameters_binSize of Pair.cpp: m_rMin and m_muMin must be >0!");
  
  m_nbins_D1 = nint((log10(m_rMax)-log10(m_rMin))*m_binSize_inv_D1);
  m_rMax = pow(10.,m_nbins_D1/m_binSize_inv_D1+log10(m_rMin));
  
  m_nbins_D2 = nint((log10(m_muMax)-log10(m_muMin))*m_binSize_inv_D2);
  m_muMax = pow(10.,m_nbins_D2/m_binSize_inv_D2+log10(m_muMin));

  m_scale_D1.resize(m_nbins_D1); m_scale_D2.resize(m_nbins_D2);
  for (int i=0; i<m_nbins_D1; i++) 
    m_scale_D1[i] = pow(10.,(i+m_shift_D1)/m_binSize_inv_D1+log10(m_rMin));
  for (int i=0; i<m_nbins_D2; i++) 
    m_scale_D2[i] = pow(10.,(i+m_shift_D2)/m_binSize_inv_D2+log10(m_muMin));
}


// ============================================================================================


void cosmobl::pairs::Pair2D_comovingCartesian_linlin::put (const shared_ptr<Object> obj1, const shared_ptr<Object> obj2) 
{
  const double rp = perpendicular_distance(obj1->ra(), obj2->ra(), obj1->dec(), obj2->dec(), obj1->dc(), obj2->dc());
  const double pi = fabs(obj1->dc()-obj2->dc());
  
  if (m_rpMin<rp && rp<m_rpMax && m_piMin<pi && pi<m_piMax) {

    const int ir = max(0, min(int((rp-m_rpMin)*m_binSize_inv_D1), m_nbins_D1));
    const int jr = max(0, min(int((pi-m_piMin)*m_binSize_inv_D2), m_nbins_D2));

    const double angWeight = (m_angularWeight==nullptr) ? 1.
      : max(0., m_angularWeight(converted_angle(angular_distance(obj1->xx()/obj1->dc(), obj2->xx()/obj2->dc(), obj1->yy()/obj1->dc(), obj2->yy()/obj2->dc(), obj1->zz()/obj1->dc(), obj2->zz()/obj2->dc()), _radians_, m_angularUnits)));
    
    m_PP2D[ir][jr] ++;
    m_PP2D_weighted[ir][jr] += obj1->weight()*obj2->weight()*angWeight;
    
  }
}


// ============================================================================================


void cosmobl::pairs::Pair2D_comovingCartesian_linlog::put (const shared_ptr<Object> obj1, const shared_ptr<Object> obj2) 
{
  const double rp = perpendicular_distance(obj1->ra(), obj2->ra(), obj1->dec(), obj2->dec(), obj1->dc(), obj2->dc());
  const double pi = fabs(obj1->dc()-obj2->dc());
  
  if (m_rpMin<rp && rp<m_rpMax && m_piMin<pi && pi<m_piMax) {
    
    const int ir = max(0, min(int((rp-m_rpMin)*m_binSize_inv_D1), m_nbins_D1));
    const int jr = max(0, min(int((log10(pi)-log10(m_piMin))*m_binSize_inv_D2), m_nbins_D2));

    const double angWeight = (m_angularWeight==nullptr) ? 1.
      : max(0., m_angularWeight(converted_angle(angular_distance(obj1->xx()/obj1->dc(), obj2->xx()/obj2->dc(), obj1->yy()/obj1->dc(), obj2->yy()/obj2->dc(), obj1->zz()/obj1->dc(), obj2->zz()/obj2->dc()), _radians_, m_angularUnits)));
    
    m_PP2D[ir][jr] ++;
    m_PP2D_weighted[ir][jr] += obj1->weight()*obj2->weight()*angWeight;
    
  }
}


// ============================================================================================


void cosmobl::pairs::Pair2D_comovingCartesian_loglin::put (const shared_ptr<Object> obj1, const shared_ptr<Object> obj2) 
{
  const double rp = perpendicular_distance(obj1->ra(), obj2->ra(), obj1->dec(), obj2->dec(), obj1->dc(), obj2->dc());
  const double pi = fabs(obj1->dc()-obj2->dc());

  if (m_rpMin<rp && rp<m_rpMax && m_piMin<pi && pi<m_piMax) {

    const int ir = max(0, min(int((log10(rp)-log10(m_rpMin))*m_binSize_inv_D1), m_nbins_D1)); 
    const int jr = max(0, min(int((pi-m_piMin)*m_binSize_inv_D2), m_nbins_D2));

    const double angWeight = (m_angularWeight==nullptr) ? 1.
      : max(0., m_angularWeight(converted_angle(angular_distance(obj1->xx()/obj1->dc(), obj2->xx()/obj2->dc(), obj1->yy()/obj1->dc(), obj2->yy()/obj2->dc(), obj1->zz()/obj1->dc(), obj2->zz()/obj2->dc()), _radians_, m_angularUnits)));
 
    m_PP2D[ir][jr] ++;
    m_PP2D_weighted[ir][jr] += obj1->weight()*obj2->weight()*angWeight;
    
  }
}


// ============================================================================================


void cosmobl::pairs::Pair2D_comovingCartesian_loglog::put (const shared_ptr<Object> obj1, const shared_ptr<Object> obj2) 
{
  const double rp = perpendicular_distance(obj1->ra(), obj2->ra(), obj1->dec(), obj2->dec(), obj1->dc(), obj2->dc());
  const double pi = fabs(obj1->dc()-obj2->dc());
  
  if (m_rpMin<rp && rp<m_rpMax && m_piMin<pi && pi<m_piMax) {

    const int ir = max(0, min(int((log10(rp)-log10(m_rpMin))*m_binSize_inv_D1), m_nbins_D1)); 
    const int jr = max(0, min(int((log10(pi)-log10(m_piMin))*m_binSize_inv_D2), m_nbins_D2)); 

    const double angWeight = (m_angularWeight==nullptr) ? 1.
      : max(0., m_angularWeight(converted_angle(angular_distance(obj1->xx()/obj1->dc(), obj2->xx()/obj2->dc(), obj1->yy()/obj1->dc(), obj2->yy()/obj2->dc(), obj1->zz()/obj1->dc(), obj2->zz()/obj2->dc()), _radians_, m_angularUnits)));
 
    m_PP2D[ir][jr] ++;
    m_PP2D_weighted[ir][jr] += obj1->weight()*obj2->weight()*angWeight;
    
  }
}


// ============================================================================================


void cosmobl::pairs::Pair2D_comovingPolar_linlin::put (const shared_ptr<Object> obj1, const shared_ptr<Object> obj2) 
{
  const double rr = Euclidean_distance(obj1->xx(), obj2->xx(), obj1->yy(), obj2->yy(), obj1->zz(), obj2->zz());
  const double mu = fabs(obj1->dc()-obj2->dc())/rr;
  
  if (m_rMin<rr && rr<m_rMax && m_muMin<mu && mu<m_muMax) {
    
    const int ir = max(0, min(int((rr-m_rMin)*m_binSize_inv_D1), m_nbins_D1));
    const int jr = max(0, min(int((mu-m_muMin)*m_binSize_inv_D2), m_nbins_D2));

    const double angWeight = (m_angularWeight==nullptr) ? 1.
      : max(0., m_angularWeight(converted_angle(angular_distance(obj1->xx()/obj1->dc(), obj2->xx()/obj2->dc(), obj1->yy()/obj1->dc(), obj2->yy()/obj2->dc(), obj1->zz()/obj1->dc(), obj2->zz()/obj2->dc()), _radians_, m_angularUnits)));
 
    m_PP2D[ir][jr] ++;
    m_PP2D_weighted[ir][jr] += obj1->weight()*obj2->weight()*angWeight;
    
  }
}


// ============================================================================================


void cosmobl::pairs::Pair2D_comovingPolar_linlog::put (const shared_ptr<Object> obj1, const shared_ptr<Object> obj2) 
{
  const double rr = Euclidean_distance(obj1->xx(), obj2->xx(), obj1->yy(), obj2->yy(), obj1->zz(), obj2->zz());
  const double mu = fabs(obj1->dc()-obj2->dc())/rr;
  
  if (m_rMin<rr && rr<m_rMax && m_muMin<mu && mu<m_muMax) {
    
    const int ir = max(0, min(int((rr-m_rMin)*m_binSize_inv_D1), m_nbins_D1));
    const int jr = max(0, min(int((log10(mu)-log10(m_muMin))*m_binSize_inv_D2), m_nbins_D2)); 

    const double angWeight = (m_angularWeight==nullptr) ? 1.
      : max(0., m_angularWeight(converted_angle(angular_distance(obj1->xx()/obj1->dc(), obj2->xx()/obj2->dc(), obj1->yy()/obj1->dc(), obj2->yy()/obj2->dc(), obj1->zz()/obj1->dc(), obj2->zz()/obj2->dc()), _radians_, m_angularUnits)));
    
    m_PP2D[ir][jr] ++;
    m_PP2D_weighted[ir][jr] += obj1->weight()*obj2->weight()*angWeight;
    
  }
}


// ============================================================================================


void cosmobl::pairs::Pair2D_comovingPolar_loglin::put (const shared_ptr<Object> obj1, const shared_ptr<Object> obj2) 
{
  const double rr = Euclidean_distance(obj1->xx(), obj2->xx(), obj1->yy(), obj2->yy(), obj1->zz(), obj2->zz());
  const double mu = fabs(obj1->dc()-obj2->dc())/rr;

  if (m_rMin<rr && rr<m_rMax && m_muMin<mu && mu<m_muMax) {

    const int ir = max(0, min(int((log10(rr)-log10(m_rMin))*m_binSize_inv_D1), m_nbins_D1)); 
    const int jr = max(0, min(int((mu-m_muMin)*m_binSize_inv_D2), m_nbins_D2));

    const double angWeight = (m_angularWeight==nullptr) ? 1.
      : max(0., m_angularWeight(converted_angle(angular_distance(obj1->xx()/obj1->dc(), obj2->xx()/obj2->dc(), obj1->yy()/obj1->dc(), obj2->yy()/obj2->dc(), obj1->zz()/obj1->dc(), obj2->zz()/obj2->dc()), _radians_, m_angularUnits)));
    
    m_PP2D[ir][jr] ++;
    m_PP2D_weighted[ir][jr] += obj1->weight()*obj2->weight()*angWeight;
    
  }
}


// ============================================================================================


void cosmobl::pairs::Pair2D_comovingPolar_loglog::put (const shared_ptr<Object> obj1, const shared_ptr<Object> obj2) 
{
  const double rr = Euclidean_distance(obj1->xx(), obj2->xx(), obj1->yy(), obj2->yy(), obj1->zz(), obj2->zz());
  const double mu = fabs(obj1->dc()-obj2->dc())/rr;
  
  if (m_rMin<rr && rr<m_rMax && m_muMin<mu && mu<m_muMax) {

    const int ir = max(0, min(int((log10(rr)-log10(m_rMin))*m_binSize_inv_D1), m_nbins_D1)); 
    const int jr = max(0, min(int((log10(mu)-log10(m_muMin))*m_binSize_inv_D2), m_nbins_D2)); 

    const double angWeight = (m_angularWeight==nullptr) ? 1.
      : max(0., m_angularWeight(converted_angle(angular_distance(obj1->xx()/obj1->dc(), obj2->xx()/obj2->dc(), obj1->yy()/obj1->dc(), obj2->yy()/obj2->dc(), obj1->zz()/obj1->dc(), obj2->zz()/obj2->dc()), _radians_, m_angularUnits)));
    
    m_PP2D[ir][jr] ++;
    m_PP2D_weighted[ir][jr] += obj1->weight()*obj2->weight()*angWeight;
    
  }
}


// ============================================================================================


void cosmobl::pairs::Pair2D::add_data2D (const int i, const int j, const vector<double> data)
{
  /*
    checkDim(m_PP2D, i, j, "m_PP2D", false);
    checkDim(m_PP2D_weighted, i, j, "m_PP2D_weighted", false);
    checkDim(data, 1, "data");
  */
  
  m_PP2D[i][j] += data[0];
  m_PP2D_weighted[i][j] += data[1];
}


// ============================================================================================


void cosmobl::pairs::Pair2D::add_data2D (const int i, const int j, const shared_ptr<pairs::Pair> pair, const double ww)
{
  add_data2D(i, j, {ww*pair->PP2D(i, j), ww*pair->PP2D_weighted(i, j)});
}


// ============================================================================================


void cosmobl::pairs::Pair2D::Sum (const shared_ptr<Pair> pair, const double ww)
{
  if (m_nbins_D1 != pair->nbins_D1() || m_nbins_D2 != pair->nbins_D2()) 
    ErrorCBL("Error in cosmobl::pairs::Pair2D::Sum of Pair.cpp: dimension problems!");
    
  for (int i=0; i<m_nbins_D1; i++)
    for (int j=0; j<m_nbins_D2; j++) 
      add_data2D(i, j, pair, ww);
}

