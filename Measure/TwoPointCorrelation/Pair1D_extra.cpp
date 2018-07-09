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
 *  @file Measure/TwoPointCorrelation/Pair1D_extra.cpp
 *
 *  @brief Methods of the classes Pair1D_extra*  
 *
 *  This file contains the implementation of all the methods of the
 *  classes Pair1D_extra*, used to handle 1D pairs of objects of any
 *  kind, with extra information
 *
 *  @author Federico Marulli 
 *
 *  @author federico.marulli3@unbo.it
 */

#include "Pair1D_extra.h"

using namespace std;

using namespace cbl;
using namespace catalogue;
using namespace pairs;


// ============================================================================================


void cbl::pairs::Pair1D_angular_lin_extra::put (const shared_ptr<Object> obj1, const shared_ptr<Object> obj2)
{
  const double dist = (m_angularUnits==CoordinateUnits::_radians_) ? angular_distance(obj1->xx(), obj2->xx(), obj1->yy(), obj2->yy(), obj1->zz(), obj2->zz()) 
	     : converted_angle(angular_distance(obj1->xx(), obj2->xx(), obj1->yy(), obj2->yy(), obj1->zz(), obj2->zz()), CoordinateUnits::_radians_, m_angularUnits);

  if (m_thetaMin < dist && dist < m_thetaMax) {

    const int kk = max(0, min(int((dist-m_thetaMin)*m_binSize_inv), m_nbins));

    const double WeightTOT = obj1->weight()*obj2->weight();
    
    m_PP1D[kk] ++;
    m_PP1D_weighted[kk] += WeightTOT;
    
    if (m_PP1D_weighted[kk]>0) {

      const double scale_mean_p = m_scale_mean[kk];    
      m_scale_mean[kk] += WeightTOT/m_PP1D_weighted[kk]*(dist-scale_mean_p);
      m_scale_S[kk] += WeightTOT*(dist-scale_mean_p)*(dist-m_scale_mean[kk]);
      
      const double pair_redshift = (obj1->redshift()>0 && obj2->redshift()>0) ? (obj1->redshift()+obj2->redshift())*0.5 : -1.;
      const double z_mean_p = m_z_mean[kk];
      m_z_mean[kk] += WeightTOT/m_PP1D_weighted[kk]*(pair_redshift-z_mean_p); 
      m_z_S[kk] += WeightTOT*(pair_redshift-z_mean_p)*(pair_redshift-m_z_mean[kk]);
      
    }
    
  } 
}


// ============================================================================================


void cbl::pairs::Pair1D_angular_log_extra::put (const shared_ptr<Object> obj1, const shared_ptr<Object> obj2)
{
  const double dist = (m_angularUnits==CoordinateUnits::_radians_) ? angular_distance(obj1->xx(), obj2->xx(), obj1->yy(), obj2->yy(), obj1->zz(), obj2->zz()) 
	     : converted_angle(angular_distance(obj1->xx(), obj2->xx(), obj1->yy(), obj2->yy(), obj1->zz(), obj2->zz()), CoordinateUnits::_radians_, m_angularUnits);
  
  if (m_thetaMin < dist && dist < m_thetaMax) {

    const int kk = max(0, min(int((log10(dist)-log10(m_thetaMin))*m_binSize_inv), m_nbins));

    const double WeightTOT = obj1->weight()*obj2->weight();
    
    m_PP1D[kk] ++;
    m_PP1D_weighted[kk] += WeightTOT;

    if (m_PP1D_weighted[kk]>0) {

      const double scale_mean_p = m_scale_mean[kk];    
      m_scale_mean[kk] += WeightTOT/m_PP1D_weighted[kk]*(dist-scale_mean_p);
      m_scale_S[kk] += WeightTOT*(dist-scale_mean_p)*(dist-m_scale_mean[kk]);
      
      const double pair_redshift = (obj1->redshift()>0 && obj2->redshift()>0) ? (obj1->redshift()+obj2->redshift())*0.5 : -1.;
      const double z_mean_p = m_z_mean[kk];
      m_z_mean[kk] += WeightTOT/m_PP1D_weighted[kk]*(pair_redshift-z_mean_p); 
      m_z_S[kk] += WeightTOT*(pair_redshift-z_mean_p)*(pair_redshift-m_z_mean[kk]);
      
    }
    
  }
}


// ============================================================================================


void cbl::pairs::Pair1D_comoving_lin_extra::put (const shared_ptr<Object> obj1, const shared_ptr<Object> obj2)
{
  const double dist = Euclidean_distance(obj1->xx(), obj2->xx(), obj1->yy(), obj2->yy(), obj1->zz(), obj2->zz()); 
  
  if (m_rMin < dist && dist < m_rMax) {

    const int kk = max(0, min(int((dist-m_rMin)*m_binSize_inv), m_nbins));
    
    const double angWeight = (m_angularWeight==nullptr) ? 1.
      : max(0., m_angularWeight(converted_angle(angular_distance(obj1->xx()/obj1->dc(), obj2->xx()/obj2->dc(), obj1->yy()/obj1->dc(), obj2->yy()/obj2->dc(), obj1->zz()/obj1->dc(), obj2->zz()/obj2->dc()), CoordinateUnits::_radians_, m_angularUnits)));

    const double WeightTOT = obj1->weight()*obj2->weight()*angWeight;
 
    m_PP1D[kk] ++;
    m_PP1D_weighted[kk] += WeightTOT;

    if (m_PP1D_weighted[kk]>0) {
      
      const double scale_mean_p = m_scale_mean[kk];    
      m_scale_mean[kk] += WeightTOT/m_PP1D_weighted[kk]*(dist-scale_mean_p);
      m_scale_S[kk] += WeightTOT*(dist-scale_mean_p)*(dist-m_scale_mean[kk]);
      
      const double pair_redshift = (obj1->redshift()>0 && obj2->redshift()>0) ? (obj1->redshift()+obj2->redshift())*0.5 : -1.;
      const double z_mean_p = m_z_mean[kk];
      m_z_mean[kk] += WeightTOT/m_PP1D_weighted[kk]*(pair_redshift-z_mean_p); 
      m_z_S[kk] += WeightTOT*(pair_redshift-z_mean_p)*(pair_redshift-m_z_mean[kk]);
      
    }
    
  }
}


// ============================================================================================


void cbl::pairs::Pair1D_comoving_log_extra::put (const shared_ptr<Object> obj1, const shared_ptr<Object> obj2)
{ 
  const double dist = Euclidean_distance(obj1->xx(), obj2->xx(), obj1->yy(), obj2->yy(), obj1->zz(), obj2->zz()); 

  if (m_rMin < dist && dist < m_rMax) {

    const int kk = max(0, min(int((log10(dist)-log10(m_rMin))*m_binSize_inv), m_nbins));

    const double angWeight = (m_angularWeight==nullptr) ? 1.
      : max(0., m_angularWeight(converted_angle(angular_distance(obj1->xx()/obj1->dc(), obj2->xx()/obj2->dc(), obj1->yy()/obj1->dc(), obj2->yy()/obj2->dc(), obj1->zz()/obj1->dc(), obj2->zz()/obj2->dc()), CoordinateUnits::_radians_, m_angularUnits)));

    const double WeightTOT = obj1->weight()*obj2->weight()*angWeight;
    
    m_PP1D[kk] ++;
    m_PP1D_weighted[kk] += WeightTOT;

    if (m_PP1D_weighted[kk]>0) {

      const double scale_mean_p = m_scale_mean[kk];    
      m_scale_mean[kk] += WeightTOT/m_PP1D_weighted[kk]*(dist-scale_mean_p);
      m_scale_S[kk] += WeightTOT*(dist-scale_mean_p)*(dist-m_scale_mean[kk]);
      
      const double pair_redshift = (obj1->redshift()>0 && obj2->redshift()>0) ? (obj1->redshift()+obj2->redshift())*0.5 : -1.;
      const double z_mean_p = m_z_mean[kk];
      m_z_mean[kk] += WeightTOT/m_PP1D_weighted[kk]*(pair_redshift-z_mean_p); 
      m_z_S[kk] += WeightTOT*(pair_redshift-z_mean_p)*(pair_redshift-m_z_mean[kk]);
      
    }
    
  }
}


// ============================================================================================


void cbl::pairs::Pair1D_comoving_multipoles_lin_extra::put (const shared_ptr<Object> obj1, const shared_ptr<Object> obj2)
{
  const double dist = Euclidean_distance(obj1->xx(), obj2->xx(), obj1->yy(), obj2->yy(), obj1->zz(), obj2->zz()); 
  
  if (m_rMin < dist && dist < m_rMax) {

    const int kk = max(0, min(int((dist-m_rMin)*m_binSize_inv), m_nbins));
    
    const double angWeight = (m_angularWeight==nullptr) ? 1.
      : max(0., m_angularWeight(converted_angle(angular_distance(obj1->xx()/obj1->dc(), obj2->xx()/obj2->dc(), obj1->yy()/obj1->dc(), obj2->yy()/obj2->dc(), obj1->zz()/obj1->dc(), obj2->zz()/obj2->dc()), CoordinateUnits::_radians_, m_angularUnits)));

    const double WeightTOT = obj1->weight()*obj2->weight()*angWeight;
    const double cosmu2 = pow((obj2->dc()-obj1->dc())/dist, 2);
 
    m_PP1D[kk] ++;
    m_PP1D_weighted[kk] += WeightTOT;

    double leg_pol_2 = 0.5*(3.*cosmu2-1.);
    m_PP1D[(m_nbins+1)+kk] += 5.*leg_pol_2 ;
    m_PP1D_weighted[(m_nbins+1)+kk] += 5.*WeightTOT*leg_pol_2;

    double leg_pol_4 = 0.125*(35.*cosmu2*cosmu2-30.*cosmu2+3.);
    m_PP1D[2.*(m_nbins+1)+kk] += 9.*leg_pol_4;
    m_PP1D_weighted[2.*(m_nbins+1)+kk] += 9.*WeightTOT*leg_pol_4;
    
    if (m_PP1D_weighted[kk]>0) {
      
      const double scale_mean_p = m_scale_mean[kk];    
      m_scale_mean[kk] += WeightTOT/m_PP1D_weighted[kk]*(dist-scale_mean_p);
      m_scale_S[kk] += WeightTOT*(dist-scale_mean_p)*(dist-m_scale_mean[kk]);
      
      const double pair_redshift = (obj1->redshift()>0 && obj2->redshift()>0) ? (obj1->redshift()+obj2->redshift())*0.5 : -1.;
      const double z_mean_p = m_z_mean[kk];
      m_z_mean[kk] += WeightTOT/m_PP1D_weighted[kk]*(pair_redshift-z_mean_p); 
      m_z_S[kk] += WeightTOT*(pair_redshift-z_mean_p)*(pair_redshift-m_z_mean[kk]);
      
      m_scale_mean[kk+(m_nbins+1)] += WeightTOT/m_PP1D_weighted[kk]*(dist-scale_mean_p);
      m_scale_S[kk+(m_nbins+1)] += WeightTOT*(dist-scale_mean_p)*(dist-m_scale_mean[kk]);

      m_scale_mean[kk+2*(m_nbins+1)] += WeightTOT/m_PP1D_weighted[kk]*(dist-scale_mean_p);
      m_scale_S[kk+2*(m_nbins+1)] += WeightTOT*(dist-scale_mean_p)*(dist-m_scale_mean[kk]);
    }
  }
}


// ============================================================================================


void cbl::pairs::Pair1D_comoving_multipoles_log_extra::put (const shared_ptr<Object> obj1, const shared_ptr<Object> obj2)
{
  const double dist = Euclidean_distance(obj1->xx(), obj2->xx(), obj1->yy(), obj2->yy(), obj1->zz(), obj2->zz()); 
  
  if (m_rMin < dist && dist < m_rMax) {

    const int kk = max(0, min(int((log10(dist)-log10(m_rMin))*m_binSize_inv), m_nbins));
    
    const double angWeight = (m_angularWeight==nullptr) ? 1.
      : max(0., m_angularWeight(converted_angle(angular_distance(obj1->xx()/obj1->dc(), obj2->xx()/obj2->dc(), obj1->yy()/obj1->dc(), obj2->yy()/obj2->dc(), obj1->zz()/obj1->dc(), obj2->zz()/obj2->dc()), CoordinateUnits::_radians_, m_angularUnits)));

    const double WeightTOT = obj1->weight()*obj2->weight()*angWeight;
    const double cosmu2 = pow((obj2->dc()-obj1->dc())/dist, 2);
 
    m_PP1D[kk] ++;
    m_PP1D_weighted[kk] += WeightTOT;

    double leg_pol_2 = 0.5*(3.*cosmu2-1.);
    m_PP1D[(m_nbins+1)+kk] += 5.*leg_pol_2 ;
    m_PP1D_weighted[(m_nbins+1)+kk] += 5.*WeightTOT*leg_pol_2;

    double leg_pol_4 = 0.125*(35.*cosmu2*cosmu2-30.*cosmu2+3.);
    m_PP1D[2.*(m_nbins+1)+kk] += 9.*leg_pol_4;
    m_PP1D_weighted[2.*(m_nbins+1)+kk] += 9.*WeightTOT*leg_pol_4;
    
    if (m_PP1D_weighted[kk]>0) {
      
      const double scale_mean_p = m_scale_mean[kk];    
      m_scale_mean[kk] += WeightTOT/m_PP1D_weighted[kk]*(dist-scale_mean_p);
      m_scale_S[kk] += WeightTOT*(dist-scale_mean_p)*(dist-m_scale_mean[kk]);
      
      const double pair_redshift = (obj1->redshift()>0 && obj2->redshift()>0) ? (obj1->redshift()+obj2->redshift())*0.5 : -1.;
      const double z_mean_p = m_z_mean[kk];
      m_z_mean[kk] += WeightTOT/m_PP1D_weighted[kk]*(pair_redshift-z_mean_p); 
      m_z_S[kk] += WeightTOT*(pair_redshift-z_mean_p)*(pair_redshift-m_z_mean[kk]);
            
      m_scale_mean[kk+(m_nbins+1)] += WeightTOT/m_PP1D_weighted[kk]*(dist-scale_mean_p);
      m_scale_S[kk+(m_nbins+1)] += WeightTOT*(dist-scale_mean_p)*(dist-m_scale_mean[kk]);

      m_scale_mean[kk+2*(m_nbins+1)] += WeightTOT/m_PP1D_weighted[kk]*(dist-scale_mean_p);
      m_scale_S[kk+2*(m_nbins+1)] += WeightTOT*(dist-scale_mean_p)*(dist-m_scale_mean[kk]);
    }

  }
}


// ============================================================================================


void cbl::pairs::Pair1D_extra::add_data1D (const int i, const vector<double> data)
{
  /*
    checkDim(m_PP1D, i, "m_PP1D", false);
    checkDim(m_PP1D_weighted, i, "m_PP1D_weighted", false);
    checkDim(m_scale_mean, i, "m_scale_mean", false);
    checkDim(m_scale_sigma, i, "m_scale_sigma", false);
    checkDim(m_z_mean, i, "m_z_mean", false);
    checkDim(m_z_sigma, i, "m_z_sigma", false); 
    checkDim(data, 6, "data");
  */

  // mean scale up to this computation step
  const double scale_mean_temp_p = m_scale_mean[i];
    
  // mean redshift up to this computation step
  const double z_mean_temp_p = m_z_mean[i];

  // number of pairs
  m_PP1D[i] += data[0];

  // number of weighted pairs
  m_PP1D_weighted[i] += data[1];
  
 
  if (m_PP1D_weighted[i]>0) {
    
    // mean separation
    m_scale_mean[i] += data[1]/m_PP1D_weighted[i]*(data[2]-scale_mean_temp_p);
    
    // mean redshift
    m_z_mean[i] += data[1]/m_PP1D_weighted[i]*(data[4]-z_mean_temp_p);

    // compute the weighted standard deviation of the scale distribution
    m_scale_S[i] += data[3]+pow(data[2]-scale_mean_temp_p, 2)*data[1]*(m_PP1D_weighted[i]-data[1])/m_PP1D_weighted[i];
    m_scale_sigma[i] = sqrt(m_scale_S[i]/m_PP1D_weighted[i]);

    // compute the weighted standard deviation of the redshift distribution
    m_z_S[i] += data[5]+pow(data[4]-z_mean_temp_p, 2)*data[1]*(m_PP1D_weighted[i]-data[1])/m_PP1D_weighted[i];
    m_z_sigma[i] = sqrt(m_z_S[i]/m_PP1D_weighted[i]);
    
  }
  
}


// ============================================================================================
     

void cbl::pairs::Pair1D_extra::add_data1D (const int i, const shared_ptr<pairs::Pair> pair, const double ww)
{
  add_data1D(i, {ww*pair->PP1D(i), ww*pair->PP1D_weighted(i), pair->scale_mean(i), pair->scale_S(i), pair->z_mean(i), pair->z_S(i)});
}


// ============================================================================================


void cbl::pairs::Pair1D_extra::Sum (const shared_ptr<Pair> pair, const double ww)
{
  if (m_nbins != pair->nbins()) 
    ErrorCBL("Error in cbl::pairs::Pair1D::Sum of Pair.cpp: dimension problems!");

  for (int i=0; i<m_nbins; ++i) 
    add_data1D(i, pair, ww);   
}


// ============================================================================================


void cbl::pairs::Pair1D_comoving_multipoles_lin_extra::Sum (const shared_ptr<Pair> pair, const double ww)
{
  if (m_nbins != pair->nbins()) 
    ErrorCBL("Error in cbl::pairs::Pair1D_comoving_multipoles_lin_extra::Sum of Pair.cpp: dimension problems!");
  
  for (int l=0; l<3; ++l)
    for (int i=0; i<m_nbins; ++i)
      add_data1D(l*(m_nbins+1)+i, pair, ww);
}


// ============================================================================================


void cbl::pairs::Pair1D_comoving_multipoles_log_extra::Sum (const shared_ptr<Pair> pair, const double ww)
{
  if (m_nbins != pair->nbins()) 
    ErrorCBL("Error in cbl::pairs::Pair1D_comoving_multipoles_log_extra::Sum of Pair.cpp: dimension problems!");
  
  for (int l=0; l<3; ++l)
    for (int i=0; i<m_nbins; ++i)
      add_data1D(l*(m_nbins+1)+i, pair, ww);
}
