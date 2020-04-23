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
 *  @file Measure/TwoPointCorrelation/Pair2D_extra.cpp
 *
 *  @brief Methods of the classese Pair2D_extra*
 *
 *  This file contains the implementation of all the methods of the
 *  classes Pair2D_extra*, used to handle 2D pairs of objects of any
 *  kind, with extra information
 *
 *  @author Federico Marulli 
 *
 *  @author federico.marulli3@unibo.it
 */

#include "Pair2D_extra.h"

using namespace std;

using namespace cbl;
using namespace catalogue;
using namespace pairs;


// ============================================================================================


void cbl::pairs::Pair2D_comovingCartesian_linlin_extra::put (const std::shared_ptr<Object> obj1, const std::shared_ptr<Object> obj2) 
{
  const double rp = perpendicular_distance(obj1->ra(), obj2->ra(), obj1->dec(), obj2->dec(), obj1->dc(), obj2->dc());
  const double pi = fabs(obj1->dc()-obj2->dc());
  
  if (m_rpMin<rp && rp<m_rpMax && m_piMin<pi && pi<m_piMax) {

    const int ir = max(0, min(int((rp-m_rpMin)*m_binSize_inv_D1), m_nbins_D1));
    const int jr = max(0, min(int((pi-m_piMin)*m_binSize_inv_D2), m_nbins_D2));

    const double angWeight = (m_angularWeight==nullptr) ? 1.
      : max(0., m_angularWeight(converted_angle(angular_distance(obj1->xx()/obj1->dc(), obj2->xx()/obj2->dc(), obj1->yy()/obj1->dc(), obj2->yy()/obj2->dc(), obj1->zz()/obj1->dc(), obj2->zz()/obj2->dc()), CoordinateUnits::_radians_, m_angularUnits)));

    const double WeightTOT = obj1->weight()*obj2->weight()*angWeight;
    
    m_PP2D[ir][jr] ++;
    m_PP2D_weighted[ir][jr] += WeightTOT;
   
    if (m_PP2D_weighted[ir][jr]>0) {
  
      const double scale_D1_mean_p = m_scale_D1_mean[ir][jr];
      const double scale_D2_mean_p = m_scale_D2_mean[ir][jr];
      m_scale_D1_mean[ir][jr] += WeightTOT/m_PP2D_weighted[ir][jr]*(rp-scale_D1_mean_p);
      m_scale_D2_mean[ir][jr] += WeightTOT/m_PP2D_weighted[ir][jr]*(pi-scale_D2_mean_p);
      m_scale_D1_S[ir][jr] += WeightTOT*(rp-scale_D1_mean_p)*(rp-m_scale_D1_mean[ir][jr]);
      m_scale_D2_S[ir][jr] += WeightTOT*(pi-scale_D2_mean_p)*(pi-m_scale_D2_mean[ir][jr]);
      
      const double pair_redshift = (obj1->redshift()>0 && obj2->redshift()>0) ? (obj1->redshift()+obj2->redshift())*0.5 : -1.;
      const double z_mean_p = m_z_mean[ir][jr];
      m_z_mean[ir][jr] += WeightTOT/m_PP2D_weighted[ir][jr]*(pair_redshift-z_mean_p); 
      m_z_S[ir][jr] += WeightTOT*(pair_redshift-z_mean_p)*(pair_redshift-m_z_mean[ir][jr]);
      
    }

  }
}

// ============================================================================================


void cbl::pairs::Pair2D_comovingCartesian_linlog_extra::put (const std::shared_ptr<Object> obj1, const std::shared_ptr<Object> obj2) 
{
  const double rp = perpendicular_distance(obj1->ra(), obj2->ra(), obj1->dec(), obj2->dec(), obj1->dc(), obj2->dc());
  const double pi = fabs(obj1->dc()-obj2->dc());
  
  if (m_rpMin<rp && rp<m_rpMax && m_piMin<pi && pi<m_piMax) {
    
    const int ir = max(0, min(int((rp-m_rpMin)*m_binSize_inv_D1), m_nbins_D1));
    const int jr = max(0, min(int((log10(pi)-log10(m_piMin))*m_binSize_inv_D2), m_nbins_D2));

    const double angWeight = (m_angularWeight==nullptr) ? 1.
      : max(0., m_angularWeight(converted_angle(angular_distance(obj1->xx()/obj1->dc(), obj2->xx()/obj2->dc(), obj1->yy()/obj1->dc(), obj2->yy()/obj2->dc(), obj1->zz()/obj1->dc(), obj2->zz()/obj2->dc()), CoordinateUnits::_radians_, m_angularUnits)));

    const double WeightTOT = obj1->weight()*obj2->weight()*angWeight;
    
    m_PP2D[ir][jr] ++;
    m_PP2D_weighted[ir][jr] += WeightTOT;

    if (m_PP2D_weighted[ir][jr]>0) {

      const double scale_D1_mean_p = m_scale_D1_mean[ir][jr];
      const double scale_D2_mean_p = m_scale_D2_mean[ir][jr];
      m_scale_D1_mean[ir][jr] += WeightTOT/m_PP2D_weighted[ir][jr]*(rp-scale_D1_mean_p);
      m_scale_D2_mean[ir][jr] += WeightTOT/m_PP2D_weighted[ir][jr]*(pi-scale_D2_mean_p);
      m_scale_D1_S[ir][jr] += WeightTOT*(rp-scale_D1_mean_p)*(rp-m_scale_D1_mean[ir][jr]);
      m_scale_D2_S[ir][jr] += WeightTOT*(pi-scale_D2_mean_p)*(pi-m_scale_D2_mean[ir][jr]);

      const double pair_redshift = (obj1->redshift()>0 && obj2->redshift()>0) ? (obj1->redshift()+obj2->redshift())*0.5 : -1.;
      const double z_mean_p = m_z_mean[ir][jr];
      m_z_mean[ir][jr] += WeightTOT/m_PP2D_weighted[ir][jr]*(pair_redshift-z_mean_p); 
      m_z_S[ir][jr] += WeightTOT*(pair_redshift-z_mean_p)*(pair_redshift-m_z_mean[ir][jr]);

    }
    
  }
}


// ============================================================================================


void cbl::pairs::Pair2D_comovingCartesian_loglin_extra::put (const std::shared_ptr<Object> obj1, const std::shared_ptr<Object> obj2) 
{
  const double rp = perpendicular_distance(obj1->ra(), obj2->ra(), obj1->dec(), obj2->dec(), obj1->dc(), obj2->dc());
  const double pi = fabs(obj1->dc()-obj2->dc());

  if (m_rpMin<rp && rp<m_rpMax && m_piMin<pi && pi<m_piMax) {

    const int ir = max(0, min(int((log10(rp)-log10(m_rpMin))*m_binSize_inv_D1), m_nbins_D1)); 
    const int jr = max(0, min(int((pi-m_piMin)*m_binSize_inv_D2), m_nbins_D2));

    const double angWeight = (m_angularWeight==nullptr) ? 1.
      : max(0., m_angularWeight(converted_angle(angular_distance(obj1->xx()/obj1->dc(), obj2->xx()/obj2->dc(), obj1->yy()/obj1->dc(), obj2->yy()/obj2->dc(), obj1->zz()/obj1->dc(), obj2->zz()/obj2->dc()), CoordinateUnits::_radians_, m_angularUnits)));

    const double WeightTOT = obj1->weight()*obj2->weight()*angWeight;
 
    m_PP2D[ir][jr] ++;
    m_PP2D_weighted[ir][jr] += WeightTOT;

    if (m_PP2D_weighted[ir][jr]>0) {
    
      const double scale_D1_mean_p = m_scale_D1_mean[ir][jr];
      const double scale_D2_mean_p = m_scale_D2_mean[ir][jr];
      m_scale_D1_mean[ir][jr] += WeightTOT/m_PP2D_weighted[ir][jr]*(rp-scale_D1_mean_p);
      m_scale_D2_mean[ir][jr] += WeightTOT/m_PP2D_weighted[ir][jr]*(pi-scale_D2_mean_p);
      m_scale_D1_S[ir][jr] += WeightTOT*(rp-scale_D1_mean_p)*(rp-m_scale_D1_mean[ir][jr]);
      m_scale_D2_S[ir][jr] += WeightTOT*(pi-scale_D2_mean_p)*(pi-m_scale_D2_mean[ir][jr]);

      const double pair_redshift = (obj1->redshift()>0 && obj2->redshift()>0) ? (obj1->redshift()+obj2->redshift())*0.5 : -1.;
      const double z_mean_p = m_z_mean[ir][jr];
      m_z_mean[ir][jr] += WeightTOT/m_PP2D_weighted[ir][jr]*(pair_redshift-z_mean_p); 
      m_z_S[ir][jr] += WeightTOT*(pair_redshift-z_mean_p)*(pair_redshift-m_z_mean[ir][jr]);

    }
    
  }
}


// ============================================================================================


void cbl::pairs::Pair2D_comovingCartesian_loglog_extra::put (const std::shared_ptr<Object> obj1, const std::shared_ptr<Object> obj2) 
{
  const double rp = perpendicular_distance(obj1->ra(), obj2->ra(), obj1->dec(), obj2->dec(), obj1->dc(), obj2->dc());
  const double pi = fabs(obj1->dc()-obj2->dc());
  
  if (m_rpMin<rp && rp<m_rpMax && m_piMin<pi && pi<m_piMax) {

    const int ir = max(0, min(int((log10(rp)-log10(m_rpMin))*m_binSize_inv_D1), m_nbins_D1)); 
    const int jr = max(0, min(int((log10(pi)-log10(m_piMin))*m_binSize_inv_D2), m_nbins_D2)); 

    const double angWeight = (m_angularWeight==nullptr) ? 1.
      : max(0., m_angularWeight(converted_angle(angular_distance(obj1->xx()/obj1->dc(), obj2->xx()/obj2->dc(), obj1->yy()/obj1->dc(), obj2->yy()/obj2->dc(), obj1->zz()/obj1->dc(), obj2->zz()/obj2->dc()), CoordinateUnits::_radians_, m_angularUnits)));

    const double WeightTOT = obj1->weight()*obj2->weight()*angWeight;
 
    m_PP2D[ir][jr] ++;
    m_PP2D_weighted[ir][jr] += WeightTOT;

    if (m_PP2D_weighted[ir][jr]>0) {

      const double scale_D1_mean_p = m_scale_D1_mean[ir][jr];
      const double scale_D2_mean_p = m_scale_D2_mean[ir][jr];
      m_scale_D1_mean[ir][jr] += WeightTOT/m_PP2D_weighted[ir][jr]*(rp-scale_D1_mean_p);
      m_scale_D2_mean[ir][jr] += WeightTOT/m_PP2D_weighted[ir][jr]*(pi-scale_D2_mean_p);
      m_scale_D1_S[ir][jr] += WeightTOT*(rp-scale_D1_mean_p)*(rp-m_scale_D1_mean[ir][jr]);
      m_scale_D2_S[ir][jr] += WeightTOT*(pi-scale_D2_mean_p)*(pi-m_scale_D2_mean[ir][jr]);

      const double pair_redshift = (obj1->redshift()>0 && obj2->redshift()>0) ? (obj1->redshift()+obj2->redshift())*0.5 : -1.;
      const double z_mean_p = m_z_mean[ir][jr];
      m_z_mean[ir][jr] += WeightTOT/m_PP2D_weighted[ir][jr]*(pair_redshift-z_mean_p); 
      m_z_S[ir][jr] += WeightTOT*(pair_redshift-z_mean_p)*(pair_redshift-m_z_mean[ir][jr]);

    }
    
  }
}


// ============================================================================================


void cbl::pairs::Pair2D_comovingPolar_linlin_extra::put (const std::shared_ptr<Object> obj1, const std::shared_ptr<Object> obj2) 
{
  double rr = Euclidean_distance(obj1->xx(), obj2->xx(), obj1->yy(), obj2->yy(), obj1->zz(), obj2->zz());
  double mu = fabs(obj1->dc()-obj2->dc())/rr;
  
  if (m_rMin<rr && rr<m_rMax && m_muMin<mu && mu<m_muMax) {
    
    const int ir = max(0, min(int((rr-m_rMin)*m_binSize_inv_D1), m_nbins_D1));
    const int jr = max(0, min(int((mu-m_muMin)*m_binSize_inv_D2), m_nbins_D2));

    const double angWeight = (m_angularWeight==nullptr) ? 1.
      : max(0., m_angularWeight(converted_angle(angular_distance(obj1->xx()/obj1->dc(), obj2->xx()/obj2->dc(), obj1->yy()/obj1->dc(), obj2->yy()/obj2->dc(), obj1->zz()/obj1->dc(), obj2->zz()/obj2->dc()), CoordinateUnits::_radians_, m_angularUnits)));

    const double WeightTOT = obj1->weight()*obj2->weight()*angWeight;
 
    m_PP2D[ir][jr] ++;
    m_PP2D_weighted[ir][jr] += WeightTOT;

    if (m_PP2D_weighted[ir][jr]>0) {

      const double scale_D1_mean_p = m_scale_D1_mean[ir][jr];
      const double scale_D2_mean_p = m_scale_D2_mean[ir][jr];
      m_scale_D1_mean[ir][jr] += WeightTOT/m_PP2D_weighted[ir][jr]*(rr-scale_D1_mean_p);
      m_scale_D2_mean[ir][jr] += WeightTOT/m_PP2D_weighted[ir][jr]*(mu-scale_D2_mean_p);
      m_scale_D1_S[ir][jr] += WeightTOT*(rr-scale_D1_mean_p)*(rr-m_scale_D1_mean[ir][jr]);
      m_scale_D2_S[ir][jr] += WeightTOT*(mu-scale_D2_mean_p)*(mu-m_scale_D2_mean[ir][jr]);

      const double pair_redshift = (obj1->redshift()>0 && obj2->redshift()>0) ? (obj1->redshift()+obj2->redshift())*0.5 : -1.;
      const double z_mean_p = m_z_mean[ir][jr];
      m_z_mean[ir][jr] += WeightTOT/m_PP2D_weighted[ir][jr]*(pair_redshift-z_mean_p); 
      m_z_S[ir][jr] += WeightTOT*(pair_redshift-z_mean_p)*(pair_redshift-m_z_mean[ir][jr]);

    }
    
  }
}


// ============================================================================================


void cbl::pairs::Pair2D_comovingPolar_linlog_extra::put (const std::shared_ptr<Object> obj1, const std::shared_ptr<Object> obj2) 
{
  double rr = Euclidean_distance(obj1->xx(), obj2->xx(), obj1->yy(), obj2->yy(), obj1->zz(), obj2->zz());
  double mu = fabs(obj1->dc()-obj2->dc())/rr;
  
  if (m_rMin<rr && rr<m_rMax && m_muMin<mu && mu<m_muMax) {
    
    const int ir = max(0, min(int((rr-m_rMin)*m_binSize_inv_D1), m_nbins_D1));
    const int jr = max(0, min(int((log10(mu)-log10(m_muMin))*m_binSize_inv_D2), m_nbins_D2)); 

    const double angWeight = (m_angularWeight==nullptr) ? 1.
      : max(0., m_angularWeight(converted_angle(angular_distance(obj1->xx()/obj1->dc(), obj2->xx()/obj2->dc(), obj1->yy()/obj1->dc(), obj2->yy()/obj2->dc(), obj1->zz()/obj1->dc(), obj2->zz()/obj2->dc()), CoordinateUnits::_radians_, m_angularUnits)));

    const double WeightTOT = obj1->weight()*obj2->weight()*angWeight;
    
    m_PP2D[ir][jr] ++;
    m_PP2D_weighted[ir][jr] += WeightTOT;

    if (m_PP2D_weighted[ir][jr]>0) {

      const double scale_D1_mean_p = m_scale_D1_mean[ir][jr];
      const double scale_D2_mean_p = m_scale_D2_mean[ir][jr];
      m_scale_D1_mean[ir][jr] += WeightTOT/m_PP2D_weighted[ir][jr]*(rr-scale_D1_mean_p);
      m_scale_D2_mean[ir][jr] += WeightTOT/m_PP2D_weighted[ir][jr]*(mu-scale_D2_mean_p);
      m_scale_D1_S[ir][jr] += WeightTOT*(rr-scale_D1_mean_p)*(rr-m_scale_D1_mean[ir][jr]);
      m_scale_D2_S[ir][jr] += WeightTOT*(mu-scale_D2_mean_p)*(mu-m_scale_D2_mean[ir][jr]);

      const double pair_redshift = (obj1->redshift()>0 && obj2->redshift()>0) ? (obj1->redshift()+obj2->redshift())*0.5 : -1.;
      const double z_mean_p = m_z_mean[ir][jr];
      m_z_mean[ir][jr] += WeightTOT/m_PP2D_weighted[ir][jr]*(pair_redshift-z_mean_p); 
      m_z_S[ir][jr] += WeightTOT*(pair_redshift-z_mean_p)*(pair_redshift-m_z_mean[ir][jr]);

    }
    
  }
}


// ============================================================================================


void cbl::pairs::Pair2D_comovingPolar_loglin_extra::put (const std::shared_ptr<Object> obj1, const std::shared_ptr<Object> obj2) 
{
  double rr = Euclidean_distance(obj1->xx(), obj2->xx(), obj1->yy(), obj2->yy(), obj1->zz(), obj2->zz());
  double mu = fabs(obj1->dc()-obj2->dc())/rr;

  if (m_rMin<rr && rr<m_rMax && m_muMin<mu && mu<m_muMax) {

    const int ir = max(0, min(int((log10(rr)-log10(m_rMin))*m_binSize_inv_D1), m_nbins_D1)); 
    const int jr = max(0, min(int((mu-m_muMin)*m_binSize_inv_D2), m_nbins_D2));

    const double angWeight = (m_angularWeight==nullptr) ? 1.
      : max(0., m_angularWeight(converted_angle(angular_distance(obj1->xx()/obj1->dc(), obj2->xx()/obj2->dc(), obj1->yy()/obj1->dc(), obj2->yy()/obj2->dc(), obj1->zz()/obj1->dc(), obj2->zz()/obj2->dc()), CoordinateUnits::_radians_, m_angularUnits)));

    const double WeightTOT = obj1->weight()*obj2->weight()*angWeight;
    
    m_PP2D[ir][jr] ++;
    m_PP2D_weighted[ir][jr] += WeightTOT;

    if (m_PP2D_weighted[ir][jr]>0) {

      const double scale_D1_mean_p = m_scale_D1_mean[ir][jr];
      const double scale_D2_mean_p = m_scale_D2_mean[ir][jr];
      m_scale_D1_mean[ir][jr] += WeightTOT/m_PP2D_weighted[ir][jr]*(rr-scale_D1_mean_p);
      m_scale_D2_mean[ir][jr] += WeightTOT/m_PP2D_weighted[ir][jr]*(mu-scale_D2_mean_p);
      m_scale_D1_S[ir][jr] += WeightTOT*(rr-scale_D1_mean_p)*(rr-m_scale_D1_mean[ir][jr]);
      m_scale_D2_S[ir][jr] += WeightTOT*(mu-scale_D2_mean_p)*(mu-m_scale_D2_mean[ir][jr]);

      const double pair_redshift = (obj1->redshift()>0 && obj2->redshift()>0) ? (obj1->redshift()+obj2->redshift())*0.5 : -1.;
      const double z_mean_p = m_z_mean[ir][jr];
      m_z_mean[ir][jr] += WeightTOT/m_PP2D_weighted[ir][jr]*(pair_redshift-z_mean_p); 
      m_z_S[ir][jr] += WeightTOT*(pair_redshift-z_mean_p)*(pair_redshift-m_z_mean[ir][jr]);
      
    }
    
  }
}


// ============================================================================================


void cbl::pairs::Pair2D_comovingPolar_loglog_extra::put (const std::shared_ptr<Object> obj1, const std::shared_ptr<Object> obj2) 
{
  double rr = Euclidean_distance(obj1->xx(), obj2->xx(), obj1->yy(), obj2->yy(), obj1->zz(), obj2->zz());
  double mu = fabs(obj1->dc()-obj2->dc())/rr;
  
  if (m_rMin<rr && rr<m_rMax && m_muMin<mu && mu<m_muMax) {

    const int ir = max(0, min(int((log10(rr)-log10(m_rMin))*m_binSize_inv_D1), m_nbins_D1)); 
    const int jr = max(0, min(int((log10(mu)-log10(m_muMin))*m_binSize_inv_D2), m_nbins_D2)); 

    const double angWeight = (m_angularWeight==nullptr) ? 1.
      : max(0., m_angularWeight(converted_angle(angular_distance(obj1->xx()/obj1->dc(), obj2->xx()/obj2->dc(), obj1->yy()/obj1->dc(), obj2->yy()/obj2->dc(), obj1->zz()/obj1->dc(), obj2->zz()/obj2->dc()), CoordinateUnits::_radians_, m_angularUnits)));

    const double WeightTOT = obj1->weight()*obj2->weight()*angWeight;
    
    m_PP2D[ir][jr] ++;
    m_PP2D_weighted[ir][jr] += WeightTOT;

    if (m_PP2D_weighted[ir][jr]>0) {

      const double scale_D1_mean_p = m_scale_D1_mean[ir][jr];
      const double scale_D2_mean_p = m_scale_D2_mean[ir][jr];
      m_scale_D1_mean[ir][jr] += WeightTOT/m_PP2D_weighted[ir][jr]*(rr-scale_D1_mean_p);
      m_scale_D2_mean[ir][jr] += WeightTOT/m_PP2D_weighted[ir][jr]*(mu-scale_D2_mean_p);
      m_scale_D1_S[ir][jr] += WeightTOT*(rr-scale_D1_mean_p)*(rr-m_scale_D1_mean[ir][jr]);
      m_scale_D2_S[ir][jr] += WeightTOT*(mu-scale_D2_mean_p)*(mu-m_scale_D2_mean[ir][jr]);

      const double pair_redshift = (obj1->redshift()>0 && obj2->redshift()>0) ? (obj1->redshift()+obj2->redshift())*0.5 : -1.;
      const double z_mean_p = m_z_mean[ir][jr];
      m_z_mean[ir][jr] += WeightTOT/m_PP2D_weighted[ir][jr]*(pair_redshift-z_mean_p); 
      m_z_S[ir][jr] += WeightTOT*(pair_redshift-z_mean_p)*(pair_redshift-m_z_mean[ir][jr]);
    
    }
    
  }
}


// ============================================================================================


void cbl::pairs::Pair2D_extra::add_data2D (const int i, const int j, const std::vector<double> data)
{
  /*
    checkDim(m_PP2D, i, j, "m_PP2D", false);
    checkDim(m_PP2D_weighted, i, j, "m_PP2D_weighted", false);
    checkDim(m_scale_D1_mean, i, j, "m_scale_D1_mean", false);
    checkDim(m_scale_D2_mean, i, j, "m_scale_D2_mean", false);
    checkDim(m_scale_D1_sigma, i, j, "m_scale_D1_sigma", false);
    checkDim(m_scale_D2_sigma, i, j, "m_scale_D2_sigma", false);
    checkDim(m_z_mean, i, j, "m_z_mean", false);
    checkDim(m_z_sigma, i, j, "m_z_sigma", false);
    checkDim(data, 8, "data");
  */
      
  // mean scale in the first dimension up to this computation step
  const double scale_D1_mean_temp_p = m_scale_D1_mean[i][j];

  // mean scale in the second dimention up to this computation step
  const double scale_D2_mean_temp_p = m_scale_D2_mean[i][j];

  // mean redshift up to this computation step
  const double z_mean_temp_p = m_z_mean[i][j];
      
  // number of pairs 
  m_PP2D[i][j] += data[0];

  // number of weighted pairs
  m_PP2D_weighted[i][j] += data[1];


  if (m_PP2D_weighted[i][j]>0) {

    // mean separation in the first dimension
    m_scale_D1_mean[i][j] += data[1]/m_PP2D_weighted[i][j]*(data[2]-scale_D1_mean_temp_p);
    
    // mean separation in the first dimension
    m_scale_D2_mean[i][j] += data[1]/m_PP2D_weighted[i][j]*(data[4]-scale_D2_mean_temp_p);
    
    // mean redshift
    m_z_mean[i][j] += data[1]/m_PP2D_weighted[i][j]*(data[6]-z_mean_temp_p);
    
    // compute the weighted standard deviation of the scale distribution in the first dimension
    m_scale_D1_S[i][j] += data[3]+pow(data[2]-scale_D1_mean_temp_p, 2)*data[1]*(m_PP2D_weighted[i][j]-data[1])/m_PP2D_weighted[i][j];
    m_scale_D1_sigma[i][j] = sqrt(m_scale_D1_S[i][j]/m_PP2D_weighted[i][j]);

    // compute the weighted standard deviation of the scale distribution in the second dimension
    m_scale_D2_S[i][j] += data[5]+pow(data[4]-scale_D2_mean_temp_p, 2)*data[1]*(m_PP2D_weighted[i][j]-data[1])/m_PP2D_weighted[i][j];
    m_scale_D2_sigma[i][j] = sqrt(m_scale_D2_S[i][j]/m_PP2D_weighted[i][j]);

    // compute the weighted standard deviation of the redshift distribution
    m_z_S[i][j] += data[7]+pow(data[6]-z_mean_temp_p, 2)*data[1]*(m_PP2D_weighted[i][j]-data[1])/m_PP2D_weighted[i][j];
    m_z_sigma[i][j] = sqrt(m_z_S[i][j]/m_PP2D_weighted[i][j]);
    
  }

}


// ============================================================================================


void cbl::pairs::Pair2D_extra::add_data2D (const int i, const int j, const std::shared_ptr<pairs::Pair> pair, const double ww)
{
  add_data2D(i, j, {ww*pair->PP2D(i, j), ww*pair->PP2D_weighted(i, j), pair->scale_D1_mean(i, j), pair->scale_D1_S(i, j), pair->scale_D2_mean(i, j), pair->scale_D2_S(i, j), pair->z_mean(i, j), pair->z_S(i, j)});
}


// ============================================================================================


void cbl::pairs::Pair2D_extra::Sum (const std::shared_ptr<Pair> pair, const double ww)
{
  if (m_nbins_D1 != pair->nbins_D1() || m_nbins_D2 != pair->nbins_D2()) 
    ErrorCBL("dimension problems!", "Sum", "Pair2D_extra.cpp");
  
  for (int i=0; i<m_nbins_D1; ++i)
    for (int j=0; j<m_nbins_D2; ++j)
      add_data2D(i, j, pair, ww);
}

