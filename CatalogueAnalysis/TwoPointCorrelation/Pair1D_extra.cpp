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
 *  @file CatalogueAnalysis/TwoPointCorrelation/Pair1D_extra.cpp
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

using namespace cosmobl;
using namespace catalogue;
using namespace pairs;


// ============================================================================================


void cosmobl::pairs::Pair1D_angular_lin_extra::put (const shared_ptr<Object> obj1, const shared_ptr<Object> obj2)
{
  double dist = (m_angularUnits==_radians_) ? angular_distance(obj1->xx(), obj2->xx(), obj1->yy(), obj2->yy(), obj1->zz(), obj2->zz()) 
	     : converted_angle(angular_distance(obj1->xx(), obj2->xx(), obj1->yy(), obj2->yy(), obj1->zz(), obj2->zz()), _radians_, m_angularUnits);

  if (m_thetaMin < dist && dist < m_thetaMax) {

    int kk = max(0, min(int((dist-m_thetaMin)*m_binSize_inv), m_nbins));

    m_PP1D[kk] += obj1->weight()*obj2->weight();

    m_scale_mean[kk] += dist*obj1->weight()*obj2->weight();
    m_scale_sigma[kk] += pow(dist, 2)*obj1->weight()*obj2->weight();

    double pair_redshift = (obj1->redshift()>0 && obj2->redshift()>0) ? (obj1->redshift()+obj2->redshift())*0.5 : -1.;
    m_z_mean[kk] += pair_redshift*obj1->weight()*obj2->weight();
    m_z_sigma[kk] += pow(pair_redshift, 2)*obj1->weight()*obj2->weight();

  } 
}


// ============================================================================================


void cosmobl::pairs::Pair1D_angular_log_extra::put (const shared_ptr<Object> obj1, const shared_ptr<Object> obj2)
{
  double dist = (m_angularUnits==_radians_) ? angular_distance(obj1->xx(), obj2->xx(), obj1->yy(), obj2->yy(), obj1->zz(), obj2->zz()) 
	     : converted_angle(angular_distance(obj1->xx(), obj2->xx(), obj1->yy(), obj2->yy(), obj1->zz(), obj2->zz()), _radians_, m_angularUnits);
  
  if (m_thetaMin < dist && dist < m_thetaMax) {

    int kk = max(0, min(int((log10(dist)-log10(m_thetaMin))*m_binSize_inv), m_nbins));

    m_PP1D[kk] += obj1->weight()*obj2->weight();

    m_scale_mean[kk] += dist*obj1->weight()*obj2->weight();
    m_scale_sigma[kk] += pow(dist, 2)*obj1->weight()*obj2->weight();

    double pair_redshift = (obj1->redshift()>0 && obj2->redshift()>0) ? (obj1->redshift()+obj2->redshift())*0.5 : -1.;
    m_z_mean[kk] += pair_redshift*obj1->weight()*obj2->weight();
    m_z_sigma[kk] += pow(pair_redshift, 2)*obj1->weight()*obj2->weight();
    
  }
}


// ============================================================================================


void cosmobl::pairs::Pair1D_comoving_lin_extra::put (const shared_ptr<Object> obj1, const shared_ptr<Object> obj2)
{
  double dist = Euclidean_distance(obj1->xx(), obj2->xx(), obj1->yy(), obj2->yy(), obj1->zz(), obj2->zz()); 
  
  if (m_rMin < dist && dist < m_rMax) {

    int kk = max(0, min(int((dist-m_rMin)*m_binSize_inv), m_nbins));
    
    double angWeight = (m_angularWeight==nullptr) ? 1.
      : m_angularWeight(converted_angle(angular_distance(obj1->xx()/obj1->dc(), obj2->xx()/obj2->dc(), obj1->yy()/obj1->dc(), obj2->yy()/obj2->dc(), obj1->zz()/obj1->dc(), obj2->zz()/obj2->dc()), _radians_, m_angularUnits));
    
    m_PP1D[kk] += obj1->weight()*obj2->weight()*angWeight;

    m_scale_mean[kk] += dist*obj1->weight()*obj2->weight()*angWeight;
    m_scale_sigma[kk] += pow(dist, 2)*obj1->weight()*obj2->weight()*angWeight;

    double pair_redshift = (obj1->redshift()>0 && obj2->redshift()>0) ? (obj1->redshift()+obj2->redshift())*0.5 : -1.;
    m_z_mean[kk] += pair_redshift*obj1->weight()*obj2->weight();
    m_z_sigma[kk] += pow(pair_redshift, 2)*obj1->weight()*obj2->weight();
    
  }
}


// ============================================================================================


void cosmobl::pairs::Pair1D_comoving_log_extra::put (const shared_ptr<Object> obj1, const shared_ptr<Object> obj2)
{ 
  double dist = Euclidean_distance(obj1->xx(), obj2->xx(), obj1->yy(), obj2->yy(), obj1->zz(), obj2->zz()); 

  if (m_rMin < dist && dist < m_rMax) {

    int kk = max(0, min(int((log10(dist)-log10(m_rMin))*m_binSize_inv), m_nbins));

    double angWeight = (m_angularWeight==nullptr) ? 1.
      : m_angularWeight(converted_angle(angular_distance(obj1->xx()/obj1->dc(), obj2->xx()/obj2->dc(), obj1->yy()/obj1->dc(), obj2->yy()/obj2->dc(), obj1->zz()/obj1->dc(), obj2->zz()/obj2->dc()), _radians_, m_angularUnits));

    m_PP1D[kk] += obj1->weight()*obj2->weight()*angWeight;

    m_scale_mean[kk] += dist*obj1->weight()*obj2->weight()*angWeight;
    m_scale_sigma[kk] += pow(dist, 2)*obj1->weight()*obj2->weight()*angWeight;

    double pair_redshift = (obj1->redshift()>0 && obj2->redshift()>0) ? (obj1->redshift()+obj2->redshift())*0.5 : -1.;
    m_z_mean[kk] += pair_redshift*obj1->weight()*obj2->weight();
    m_z_sigma[kk] += pow(pair_redshift, 2)*obj1->weight()*obj2->weight();
  }
}


// ============================================================================================


void cosmobl::pairs::Pair1D_extra::add_data1D (const int i, const vector<double> data)
{
  /*
  checkDim(m_PP1D, i, "m_PP1D", false);
  checkDim(m_scale_mean, i, "m_scale_mean", false);
  checkDim(m_scale_sigma, i, "m_scale_sigma", false);
  checkDim(m_z_mean, i, "m_z_mean", false);
  checkDim(m_z_sigma, i, "m_z_sigma", false); 
  checkDim(data, 5, "data");
  */
  
  m_PP1D[i] += data[0];
  m_scale_mean[i] += data[1];
  m_scale_sigma[i] += data[2];
  m_z_mean[i] += data[3];
  m_z_sigma[i] += data[4];
}


// ============================================================================================
     

void cosmobl::pairs::Pair1D_extra::add_data1D (const int i, const shared_ptr<pairs::Pair> pair, const double ww)
{
  /*
  checkDim(m_PP1D, i, "m_PP1D", false); 
  checkDim(m_scale_mean, i, "m_scale_mean", false);
  checkDim(m_scale_sigma, i, "m_scale_sigma", false); 
  checkDim(m_z_mean, i, "m_z_mean", false);
  checkDim(m_z_sigma, i, "m_z_sigma", false); 
  */
  
  m_PP1D[i] += ww*pair->PP1D(i);
  m_scale_mean[i] += ww*pair->scale_mean(i);
  m_scale_sigma[i] += ww*pair->scale_sigma(i);
  m_z_mean[i] += ww*pair->z_mean(i);
  m_z_sigma[i] += ww*pair->z_sigma(i);
}


// ============================================================================================


void cosmobl::pairs::Pair1D_extra::Sum (const shared_ptr<Pair> pp, const double ww)
{
  if (m_nbins != pp->nbins()) 
    ErrorCBL("Error in cosmobl::pairs::Pair1D::Sum of Pair.cpp: dimension problems!");
  
  for (int i=0; i<m_nbins; ++i) {
    m_PP1D[i] += ww*pp->PP1D(i);
    m_scale_mean[i] += ww*pp->scale_mean(i);
    m_scale_sigma[i] += ww*pp->scale_sigma(i);
    m_z_mean[i] += ww*pp->z_mean(i);
    m_z_sigma[i] += ww*pp->z_sigma(i);
  }
}


// ============================================================================================


void cosmobl::pairs::Pair1D_extra::finalise ()
{
  for (int i=0; i<m_nbins; ++i) {
    m_scale_mean[i] = (m_PP1D[i]>0) ? m_scale_mean[i]/m_PP1D[i] : 0.;
    m_scale_sigma[i] = (m_PP1D[i]>0 && m_scale_sigma[i]/m_PP1D[i]-pow(m_scale_mean[i], 2)>0) ? sqrt(m_scale_sigma[i]/m_PP1D[i]-pow(m_scale_mean[i], 2)) : 0.;
    m_z_mean[i] = (m_PP1D[i]>0) ? m_z_mean[i]/m_PP1D[i] : 0.;
    m_z_sigma[i] = (m_PP1D[i]>0 && m_z_sigma[i]/m_PP1D[i]-pow(m_z_mean[i], 2)>0) ? sqrt(m_z_sigma[i]/m_PP1D[i]-pow(m_z_mean[i], 2)) : 0.;
  }
}
