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
 *  @file Measure/TwoPointCorrelation/Pair1D.cpp
 *
 *  @brief Methods of the classes Pair1D*  
 *
 *  This file contains the implementation of all the methods of the
 *  classes Pair1D*, used to handle pairs of 1D objects of any kind
 *
 *  @author Federico Marulli 
 *
 *  @author federico.marulli3@unibo.it
 */

#include "Pair1D.h"

using namespace std;

using namespace cbl;
using namespace catalogue;
using namespace pairs;


// ============================================================================================


void cbl::pairs::Pair1D::reset ()
{
  for (int i=0; i<m_nbins; i++) {
    m_PP1D[i] = 0.;
    m_PP1D_weighted[i] = 0;
  }
}


// ============================================================================================


void cbl::pairs::Pair1D_angular_lin::m_set_parameters_nbins ()
{
  const double binSize = ((m_thetaMax-m_thetaMin)/m_nbins);
  m_binSize_inv = 1./binSize;
  
  m_scale.resize(m_nbins);
  for (int i=0; i<m_nbins; i++)
    m_scale[i] = binSize*(i+m_shift)+m_thetaMin;
}


// ============================================================================================


void cbl::pairs::Pair1D_angular_lin::m_set_parameters_binSize ()
{
  m_nbins = nint((m_thetaMax-m_thetaMin)*m_binSize_inv);
  m_thetaMax = m_nbins/m_binSize_inv+m_thetaMin;

  m_scale.resize(m_nbins);
  for (int i=0; i<m_nbins; i++)
    m_scale[i] = (i+m_shift)/m_binSize_inv+m_thetaMin;
}


// ============================================================================================


void cbl::pairs::Pair1D_angular_log::m_set_parameters_nbins ()
{
  if (m_thetaMin<1.e-30) ErrorCBL("m_thetaMin must be >0!", "m_set_parameters_nbins", "Pair1D.cpp");
  
  const double binSize = ((log10(m_thetaMax)-log10(m_thetaMin))/m_nbins);
  m_binSize_inv = 1./binSize;

  m_scale.resize(m_nbins);
  for (int i=0; i<m_nbins; i++)
    m_scale[i] = pow(10.,(i+m_shift)*binSize+log10(m_thetaMin));
}


// ============================================================================================


void cbl::pairs::Pair1D_angular_log::m_set_parameters_binSize ()
{
  if (m_thetaMin<1.e-30) ErrorCBL("m_thetaMin must be >0!", "m_set_parameters_binSize", "Pair1D.cpp");

  m_nbins = nint((log10(m_thetaMax)-log10(m_thetaMin))*m_binSize_inv);
  m_thetaMax = pow(10.,m_nbins/m_binSize_inv+log10(m_thetaMin));

  m_scale.resize(m_nbins);
  for (int i=0; i<m_nbins; i++)
    m_scale[i] = pow(10.,(i+m_shift)/m_binSize_inv+log10(m_thetaMin));
}


// ============================================================================================


void cbl::pairs::Pair1D_comoving_lin::m_set_parameters_nbins ()
{
  const double binSize = ((m_rMax-m_rMin)/m_nbins);
  m_binSize_inv = 1./binSize;

  m_scale.resize(m_nbins);
  for (int i=0; i<m_nbins; i++)
    m_scale[i] = (i+m_shift)*binSize+m_rMin;
}


// ============================================================================================


void cbl::pairs::Pair1D_comoving_lin::m_set_parameters_binSize ()
{
  m_nbins = nint((m_rMax-m_rMin)*m_binSize_inv);
  m_rMax = m_nbins/m_binSize_inv+m_rMin;
  
  m_scale.resize(m_nbins);
  for (int i=0; i<m_nbins; i++)
    m_scale[i] = (i+m_shift)/m_binSize_inv+m_rMin;
}


// ============================================================================================


void cbl::pairs::Pair1D_comoving_log::m_set_parameters_nbins ()
{
  if (m_rMin<1.e-30) ErrorCBL("m_rMin must be >0!", "m_set_parameters_nbins", "Pair1D.cpp");
  
  const double binSize = ((log10(m_rMax)-log10(m_rMin))/m_nbins);
  m_binSize_inv = 1./binSize;

  m_scale.resize(m_nbins);
  for (int i=0; i<m_nbins; i++)
    m_scale[i] = pow(10.,(i+m_shift)*binSize+log10(m_rMin));
}


// ============================================================================================


void cbl::pairs::Pair1D_comoving_log::m_set_parameters_binSize ()
{
  if (m_rMin<1.e-30) ErrorCBL("m_rMin must be >0!", "m_set_parameters_binSize", "Pair1D.cpp");
  
  m_nbins = nint((log10(m_rMax)-log10(m_rMin))*m_binSize_inv);
  m_rMax = pow(10.,m_nbins/m_binSize_inv+log10(m_rMin));

  m_scale.resize(m_nbins);
  for (int i=0; i<m_nbins; i++)
    m_scale[i] = pow(10.,(i+m_shift)/m_binSize_inv+log10(m_rMin));
}


// ============================================================================================


void cbl::pairs::Pair1D_comoving_multipoles_lin::m_set_parameters_nbins ()
{
  const double binSize = ((m_rMax-m_rMin)/m_nbins);
  m_binSize_inv = 1./binSize;

  m_scale.resize(3*m_nbins);
  for (int l=0; l<3; l++)
    for (int i=0; i<m_nbins; i++)
      m_scale[l*m_nbins+i] = (i+m_shift)*binSize+m_rMin;
}


// ============================================================================================


void cbl::pairs::Pair1D_comoving_multipoles_lin::m_set_parameters_binSize ()
{
  m_nbins = nint((m_rMax-m_rMin)*m_binSize_inv);
  m_rMax = m_nbins/m_binSize_inv+m_rMin;

  m_scale.resize(3*m_nbins);
  for (int l=0; l<3; l++)
    for (int i=0; i<m_nbins; i++)
      m_scale[l*m_nbins+i] = (i+m_shift)/m_binSize_inv+m_rMin;
}


// ============================================================================================


void cbl::pairs::Pair1D_comoving_multipoles_log::m_set_parameters_nbins ()
{
  if (m_rMin<1.e-30) ErrorCBL("m_rMin must be >0!", "m_set_parameters_nbins", "Pair1D.cpp");
  
  const double binSize = ((log10(m_rMax)-log10(m_rMin))/m_nbins);
  m_binSize_inv = 1./binSize;

  m_scale.resize(3*m_nbins);
  for (int l=0; l<3; l++)
    for (int i=0; i<m_nbins; i++)
      m_scale[l*m_nbins+i] = pow(10.,(i+m_shift)*binSize+log10(m_rMin));
}


// ============================================================================================


void cbl::pairs::Pair1D_comoving_multipoles_log::m_set_parameters_binSize ()
{
  if (m_rMin<1.e-30) ErrorCBL("m_rMin must be >0!", "m_set_parameters_binSize", "Pair1D.cpp");
  
  m_nbins = nint((log10(m_rMax)-log10(m_rMin))*m_binSize_inv);
  m_rMax = pow(10.,m_nbins/m_binSize_inv+log10(m_rMin));

  m_scale.resize(3*m_nbins);
  for (int l=0; l<3; l++)
    for (int i=0; i<m_nbins; i++)
      m_scale[l*m_nbins+i] = pow(10.,(i+m_shift)/m_binSize_inv+log10(m_rMin));
}


// ============================================================================================


void cbl::pairs::Pair1D_angular_lin::get (const shared_ptr<Object> obj1, const shared_ptr<Object> obj2, int &kk, double &wkk)
{
  kk = -1;
  wkk = 0;

  const double dist = (m_angularUnits==CoordinateUnits::_radians_) ? angular_distance(obj1->xx(), obj2->xx(), obj1->yy(), obj2->yy(), obj1->zz(), obj2->zz()) 
	     : converted_angle(angular_distance(obj1->xx(), obj2->xx(), obj1->yy(), obj2->yy(), obj1->zz(), obj2->zz()), CoordinateUnits::_radians_, m_angularUnits);

  if (m_thetaMin < dist && dist < m_thetaMax) {
    kk = max(0, min(int((dist-m_thetaMin)*m_binSize_inv), m_nbins));

    wkk += obj1->weight()*obj2->weight();
  }
}


// ============================================================================================


void cbl::pairs::Pair1D_angular_log::get (const shared_ptr<Object> obj1, const shared_ptr<Object> obj2, int &kk, double &wkk)
{
  kk = -1;
  wkk = 0;

  const double dist = (m_angularUnits==CoordinateUnits::_radians_) ? angular_distance(obj1->xx(), obj2->xx(), obj1->yy(), obj2->yy(), obj1->zz(), obj2->zz()) 
	     : converted_angle(angular_distance(obj1->xx(), obj2->xx(), obj1->yy(), obj2->yy(), obj1->zz(), obj2->zz()), CoordinateUnits::_radians_, m_angularUnits);
  
  if (m_thetaMin < dist && dist < m_thetaMax) {
    kk = max(0, min(int((log10(dist)-log10(m_thetaMin))*m_binSize_inv), m_nbins));

    wkk = obj1->weight()*obj2->weight();
  }
}


// ============================================================================================


void cbl::pairs::Pair1D_comoving_lin::get (const shared_ptr<Object> obj1, const shared_ptr<Object> obj2, int &kk, double &wkk)
{
  kk = -1;
  wkk = 0;

  const double dist = Euclidean_distance(obj1->xx(), obj2->xx(), obj1->yy(), obj2->yy(), obj1->zz(), obj2->zz()); 
  
  if (m_rMin < dist && dist < m_rMax) {
    kk = max(0, min(int((dist-m_rMin)*m_binSize_inv), m_nbins));
    
    const double angWeight = (m_angularWeight==nullptr) ? 1.
      : max(0., m_angularWeight(converted_angle(angular_distance(obj1->xx()/obj1->dc(), obj2->xx()/obj2->dc(), obj1->yy()/obj1->dc(), obj2->yy()/obj2->dc(), obj1->zz()/obj1->dc(), obj2->zz()/obj2->dc()), CoordinateUnits::_radians_, m_angularUnits)));
    
    wkk = obj1->weight()*obj2->weight()*angWeight;
  }
}


// ============================================================================================


void cbl::pairs::Pair1D_comoving_log::get (const shared_ptr<Object> obj1, const shared_ptr<Object> obj2, int &kk, double &wkk)
{
  kk = -1;
  wkk = 0;

  const double dist = Euclidean_distance(obj1->xx(), obj2->xx(), obj1->yy(), obj2->yy(), obj1->zz(), obj2->zz()); 
  
  if (m_rMin < dist && dist < m_rMax) {
    kk = max(0, min(int((log10(dist)-log10(m_rMin))*m_binSize_inv), m_nbins));

    const double angWeight = (m_angularWeight==nullptr) ? 1.
      : max(0., m_angularWeight(converted_angle(angular_distance(obj1->xx()/obj1->dc(), obj2->xx()/obj2->dc(), obj1->yy()/obj1->dc(), obj2->yy()/obj2->dc(), obj1->zz()/obj1->dc(), obj2->zz()/obj2->dc()), CoordinateUnits::_radians_, m_angularUnits)));

    wkk = obj1->weight()*obj2->weight()*angWeight;
  }
}


// ============================================================================================


void cbl::pairs::Pair1D_comoving_multipoles_lin::get (const shared_ptr<Object> obj1, const shared_ptr<Object> obj2, int &kk, double &cosmu, double &wkk)
{
  kk = -1;
  wkk = 0;

  const double dist = Euclidean_distance(obj1->xx(), obj2->xx(), obj1->yy(), obj2->yy(), obj1->zz(), obj2->zz()); 
  
  if (m_rMin < dist && dist < m_rMax) {
    kk = max(0, min(int((dist-m_rMin)*m_binSize_inv), m_nbins));
    
    const double angWeight = (m_angularWeight==nullptr) ? 1.
      : max(0., m_angularWeight(converted_angle(angular_distance(obj1->xx()/obj1->dc(), obj2->xx()/obj2->dc(), obj1->yy()/obj1->dc(), obj2->yy()/obj2->dc(), obj1->zz()/obj1->dc(), obj2->zz()/obj2->dc()), CoordinateUnits::_radians_, m_angularUnits)));
    
    wkk = obj1->weight()*obj2->weight()*angWeight;

    cosmu = (obj2->dc()-obj1->dc())/dist;
  }
}


// ============================================================================================


void cbl::pairs::Pair1D_comoving_multipoles_log::get (const shared_ptr<Object> obj1, const shared_ptr<Object> obj2, int &kk, double &cosmu, double &wkk)
{
  kk = -1;
  wkk = 0;

  const double dist = Euclidean_distance(obj1->xx(), obj2->xx(), obj1->yy(), obj2->yy(), obj1->zz(), obj2->zz()); 
  
  if (m_rMin < dist && dist < m_rMax) {
    kk = max(0, min(int((log10(dist)-log10(m_rMin))*m_binSize_inv), m_nbins));

    const double angWeight = (m_angularWeight==nullptr) ? 1.
      : max(0., m_angularWeight(converted_angle(angular_distance(obj1->xx()/obj1->dc(), obj2->xx()/obj2->dc(), obj1->yy()/obj1->dc(), obj2->yy()/obj2->dc(), obj1->zz()/obj1->dc(), obj2->zz()/obj2->dc()), CoordinateUnits::_radians_, m_angularUnits)));

    wkk = obj1->weight()*obj2->weight()*angWeight;

    cosmu = (obj2->dc()-obj1->dc())/dist;
  }
}

// ============================================================================================


void cbl::pairs::Pair1D_angular_lin::set (const int kk, const double wkk, const double weight)
{
  if (kk>-1) {
    m_PP1D[kk] += weight;
    m_PP1D_weighted[kk] += wkk*weight;
  }
}


// ============================================================================================


void cbl::pairs::Pair1D_angular_log::set (const int kk, const double wkk, const double weight)
{
  if (kk>-1) {
    m_PP1D[kk] += weight;
    m_PP1D_weighted[kk] += wkk*weight;
  }
}


// ============================================================================================


void cbl::pairs::Pair1D_comoving_lin::set (const int kk, const double wkk, const double weight)
{
  if (kk>-1) {
    m_PP1D[kk] += weight;
    m_PP1D_weighted[kk] += wkk*weight;
  }
}


// ============================================================================================


void cbl::pairs::Pair1D_comoving_log::set (const int kk, const double wkk, const double weight)
{
  if (kk>-1) {
    m_PP1D[kk] += weight;
    m_PP1D_weighted[kk] += wkk*weight;
  }
}


// ============================================================================================


void cbl::pairs::Pair1D_comoving_multipoles_lin::set (const double cosmu, const int kk, const double wkk, const double weight)
{
  if (kk>-1) {
    double cosmu2 = cosmu*cosmu;
    m_PP1D[kk] += weight;
    m_PP1D_weighted[kk] += wkk*weight;

    double leg_pol_2 = 0.5*(3.*cosmu2-1);
    m_PP1D[(m_nbins+1)+kk] += 5.*leg_pol_2*weight ;
    m_PP1D_weighted[(m_nbins+1)+kk] += 5.*wkk*leg_pol_2*weight;

    double leg_pol_4 = 0.125*(35.*cosmu2*cosmu2-30.*cosmu2+3.);
    m_PP1D[2.*(m_nbins+1)+kk] += 9.*leg_pol_4*weight;
    m_PP1D_weighted[2.*(m_nbins+1)+kk] += 9.*wkk*leg_pol_4*weight;
  }
}


// ============================================================================================


void cbl::pairs::Pair1D_comoving_multipoles_log::set (const double cosmu, const int kk, const double wkk, const double weight)
{
  if (kk>-1) {
    double cosmu2 = cosmu*cosmu;
    m_PP1D[kk] += weight;
    m_PP1D_weighted[kk] += wkk*weight;

    double leg_pol_2 = 0.5*(3.*cosmu2-1);
    m_PP1D[(m_nbins+1)+kk] += 5.*leg_pol_2*weight ;
    m_PP1D_weighted[(m_nbins+1)+kk] += 5.*wkk*leg_pol_2*weight;

    double leg_pol_4 = 0.125*(35.*cosmu2*cosmu2-30.*cosmu2+3.);
    m_PP1D[2.*(m_nbins+1)+kk] += 9.*leg_pol_4*weight;
    m_PP1D_weighted[2.*(m_nbins+1)+kk] += 9.*wkk*leg_pol_4*weight;
  }
}

// ============================================================================================


void cbl::pairs::Pair1D_angular_lin::put (const shared_ptr<Object> obj1, const shared_ptr<Object> obj2)
{
  const double dist = (m_angularUnits==CoordinateUnits::_radians_) ? angular_distance(obj1->xx(), obj2->xx(), obj1->yy(), obj2->yy(), obj1->zz(), obj2->zz()) 
	     : converted_angle(angular_distance(obj1->xx(), obj2->xx(), obj1->yy(), obj2->yy(), obj1->zz(), obj2->zz()), CoordinateUnits::_radians_, m_angularUnits);

  if (m_thetaMin < dist && dist < m_thetaMax) {

    const int kk = max(0, min(int((dist-m_thetaMin)*m_binSize_inv), m_nbins));

    m_PP1D[kk] ++;
    m_PP1D_weighted[kk] += obj1->weight()*obj2->weight();
    
  }
}


// ============================================================================================


void cbl::pairs::Pair1D_angular_log::put (const shared_ptr<Object> obj1, const shared_ptr<Object> obj2)
{
  const double dist = (m_angularUnits==CoordinateUnits::_radians_) ? angular_distance(obj1->xx(), obj2->xx(), obj1->yy(), obj2->yy(), obj1->zz(), obj2->zz()) 
	     : converted_angle(angular_distance(obj1->xx(), obj2->xx(), obj1->yy(), obj2->yy(), obj1->zz(), obj2->zz()), CoordinateUnits::_radians_, m_angularUnits);
  
  if (m_thetaMin < dist && dist < m_thetaMax) {

    const int kk = max(0, min(int((log10(dist)-log10(m_thetaMin))*m_binSize_inv), m_nbins));

    m_PP1D[kk] ++;
    m_PP1D_weighted[kk] += obj1->weight()*obj2->weight();

  }
}


// ============================================================================================


void cbl::pairs::Pair1D_comoving_lin::put (const shared_ptr<Object> obj1, const shared_ptr<Object> obj2)
{
  const double dist = Euclidean_distance(obj1->xx(), obj2->xx(), obj1->yy(), obj2->yy(), obj1->zz(), obj2->zz()); 
  
  if (m_rMin < dist && dist < m_rMax) {

    const int kk = max(0, min(int((dist-m_rMin)*m_binSize_inv), m_nbins));
    
    const double angWeight = (m_angularWeight==nullptr) ? 1.
      : max(0., m_angularWeight(converted_angle(angular_distance(obj1->xx()/obj1->dc(), obj2->xx()/obj2->dc(), obj1->yy()/obj1->dc(), obj2->yy()/obj2->dc(), obj1->zz()/obj1->dc(), obj2->zz()/obj2->dc()), CoordinateUnits::_radians_, m_angularUnits)));
    
    m_PP1D[kk] ++;
    m_PP1D_weighted[kk] += obj1->weight()*obj2->weight()*angWeight;
    
  }
}


// ============================================================================================


void cbl::pairs::Pair1D_comoving_log::put (const shared_ptr<Object> obj1, const shared_ptr<Object> obj2)
{
  const double dist = Euclidean_distance(obj1->xx(), obj2->xx(), obj1->yy(), obj2->yy(), obj1->zz(), obj2->zz()); 
  
  if (m_rMin < dist && dist < m_rMax) {

    const int kk = max(0, min(int((log10(dist)-log10(m_rMin))*m_binSize_inv), m_nbins));

    const double angWeight = (m_angularWeight==nullptr) ? 1.
      : max(0., m_angularWeight(converted_angle(angular_distance(obj1->xx()/obj1->dc(), obj2->xx()/obj2->dc(), obj1->yy()/obj1->dc(), obj2->yy()/obj2->dc(), obj1->zz()/obj1->dc(), obj2->zz()/obj2->dc()), CoordinateUnits::_radians_, m_angularUnits)));

    m_PP1D[kk] ++;
    m_PP1D_weighted[kk] += obj1->weight()*obj2->weight()*angWeight;

  }
}


// ============================================================================================


void cbl::pairs::Pair1D_comoving_multipoles_lin::put (const shared_ptr<Object> obj1, const shared_ptr<Object> obj2)
{
  const double dist = Euclidean_distance(obj1->xx(), obj2->xx(), obj1->yy(), obj2->yy(), obj1->zz(), obj2->zz()); 
  
  if (m_rMin < dist && dist < m_rMax) {

    const int kk = max(0, min(int((dist-m_rMin)*m_binSize_inv), m_nbins));
    
    const double angWeight = (m_angularWeight==nullptr) ? 1.
      : max(0., m_angularWeight(converted_angle(angular_distance(obj1->xx()/obj1->dc(), obj2->xx()/obj2->dc(), obj1->yy()/obj1->dc(), obj2->yy()/obj2->dc(), obj1->zz()/obj1->dc(), obj2->zz()/obj2->dc()), CoordinateUnits::_radians_, m_angularUnits)));

    const double ww = obj1->weight()*obj2->weight()*angWeight;
    double cosmu2 = pow((obj2->dc()-obj1->dc())/dist, 2);

    m_PP1D[kk] += 1.;
    m_PP1D_weighted[kk] += ww;

    double leg_pol_2 = 0.5*(3.*cosmu2-1.);
    m_PP1D[(m_nbins+1)+kk] += 5.*leg_pol_2 ;
    m_PP1D_weighted[(m_nbins+1)+kk] += 5.*ww*leg_pol_2;

    double leg_pol_4 = 0.125*(35.*cosmu2*cosmu2-30.*cosmu2+3.);
    m_PP1D[2.*(m_nbins+1)+kk] += 9.*leg_pol_4;
    m_PP1D_weighted[2.*(m_nbins+1)+kk] += 9.*ww*leg_pol_4;

    
  }
}


// ============================================================================================


void cbl::pairs::Pair1D_comoving_multipoles_log::put (const shared_ptr<Object> obj1, const shared_ptr<Object> obj2)
{
  const double dist = Euclidean_distance(obj1->xx(), obj2->xx(), obj1->yy(), obj2->yy(), obj1->zz(), obj2->zz()); 

  if (m_rMin < dist && dist < m_rMax) {

    const int kk = max(0, min(int((log10(dist)-log10(m_rMin))*m_binSize_inv), m_nbins));

    const double angWeight = (m_angularWeight==nullptr) ? 1.
      : max(0., m_angularWeight(converted_angle(angular_distance(obj1->xx()/obj1->dc(), obj2->xx()/obj2->dc(), obj1->yy()/obj1->dc(), obj2->yy()/obj2->dc(), obj1->zz()/obj1->dc(), obj2->zz()/obj2->dc()), CoordinateUnits::_radians_, m_angularUnits)));

    const double ww = obj1->weight()*obj2->weight()*angWeight;
    double cosmu2 = pow((obj2->dc()-obj1->dc())/dist, 2);

    m_PP1D[kk] += 1.;
    m_PP1D_weighted[kk] += ww;

    double leg_pol_2 = 0.5*(3.*cosmu2-1);
    m_PP1D[(m_nbins+1)+kk] += 5.*leg_pol_2 ;
    m_PP1D_weighted[(m_nbins+1)+kk] += 5.*ww*leg_pol_2;

    double leg_pol_4 = 0.125*(35.*cosmu2*cosmu2-30.*cosmu2+3.);
    m_PP1D[2.*(m_nbins+1)+kk] += 9.*leg_pol_4;
    m_PP1D_weighted[2.*(m_nbins+1)+kk] += 9.*ww*leg_pol_4;

  }
}


// ============================================================================================


void cbl::pairs::Pair1D::add_data1D (const int i, const std::vector<double> data)
{
  m_PP1D[i] += data[0];
  m_PP1D_weighted[i] += data[1];
}


// ============================================================================================


void cbl::pairs::Pair1D::add_data1D (const int i, const std::shared_ptr<pairs::Pair> pair, const double ww) 
{
  add_data1D(i, {ww*pair->PP1D(i), ww*pair->PP1D_weighted(i)});
}


// ============================================================================================


void cbl::pairs::Pair1D::Sum (const std::shared_ptr<Pair> pair, const double ww)
{
  if (m_nbins != pair->nbins()) 
    ErrorCBL("dimension problems!", "Sum", "Pair1D.cpp");
  
  for (int i=0; i<m_nbins; ++i)
    add_data1D(i, pair, ww);
}


// ============================================================================================


void cbl::pairs::Pair1D_comoving_multipoles_lin::Sum (const std::shared_ptr<Pair> pair, const double ww)
{
  if (m_nbins != pair->nbins()) 
    ErrorCBL("dimension problems!", "Sum", "Pair1D.cpp");
  
  for (int l=0; l<3; ++l)
    for (int i=0; i<m_nbins; ++i)
      add_data1D(l*(m_nbins+1)+i, pair, ww);
}


// ============================================================================================


void cbl::pairs::Pair1D_comoving_multipoles_log::Sum (const std::shared_ptr<Pair> pair, const double ww)
{
  if (m_nbins != pair->nbins()) 
    ErrorCBL("dimension problems!", "Sum", "Pair1D.cpp");
  
  for (int l=0; l<3; ++l)
    for (int i=0; i<m_nbins; ++i)
      add_data1D(l*(m_nbins+1)+i, pair, ww);
}
