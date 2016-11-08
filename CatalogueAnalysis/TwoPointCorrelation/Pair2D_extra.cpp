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
 *  @file CatalogueAnalysis/TwoPointCorrelation/Pair2D_extra.cpp
 *
 *  @brief Methods of the classese Pair2D_extra*
 *
 *  This file contains the implementation of all the methods of the
 *  classes Pair2D_extra*, used to handle 2D pairs of objects of any
 *  kind, with extra information
 *
 *  @author Federico Marulli 
 *
 *  @author federico.marulli3@unbo.it
 */

#include "Pair2D_extra.h"

using namespace cosmobl;
using namespace catalogue;
using namespace pairs;


// ============================================================================================


void cosmobl::pairs::Pair2D_comovingCartesian_linlin_extra::put (const shared_ptr<Object> obj1, const shared_ptr<Object> obj2) 
{
  double rp = perpendicular_distance(obj1->ra(), obj2->ra(), obj1->dec(), obj2->dec(), obj1->dc(), obj2->dc());
  double pi = fabs(obj1->dc()-obj2->dc());
  
  if (m_rpMin<rp && rp<m_rpMax && m_piMin<pi && pi<m_piMax) {

    int ir = max(0, min(int((rp-m_rpMin)*m_binSize_inv_D1), m_nbins_D1));
    int jr = max(0, min(int((pi-m_piMin)*m_binSize_inv_D2), m_nbins_D2));

    double angWeight = (m_angularWeight==nullptr) ? 1.
      : m_angularWeight(converted_angle(angular_distance(obj1->xx()/obj1->dc(), obj2->xx()/obj2->dc(), obj1->yy()/obj1->dc(), obj2->yy()/obj2->dc(), obj1->zz()/obj1->dc(), obj2->zz()/obj2->dc()), _radians_, m_angularUnits));
    
    m_PP2D[ir][jr] += obj1->weight()*obj2->weight()*angWeight;

    m_scale_D1_mean[ir][jr] += rp*obj1->weight()*obj2->weight()*angWeight;
    m_scale_D2_mean[ir][jr] += pi*obj1->weight()*obj2->weight()*angWeight;
    
    m_scale_D1_sigma[ir][jr] += pow(rp*obj1->weight()*obj2->weight(), 2);
    m_scale_D2_sigma[ir][jr] += pow(pi*obj1->weight()*obj2->weight(), 2);

    double pair_redshift = (obj1->redshift()>0 && obj2->redshift()>0) ? (obj1->redshift()+obj2->redshift())*0.5 : -1.;
    m_z_mean[ir][jr] += pair_redshift*obj1->weight()*obj2->weight()*angWeight;
    m_z_sigma[ir][jr] += pow(pair_redshift*obj1->weight()*obj2->weight(), 2);
    
  }
}

// ============================================================================================


void cosmobl::pairs::Pair2D_comovingCartesian_linlog_extra::put (const shared_ptr<Object> obj1, const shared_ptr<Object> obj2) 
{
  double rp = perpendicular_distance(obj1->ra(), obj2->ra(), obj1->dec(), obj2->dec(), obj1->dc(), obj2->dc());
  double pi = fabs(obj1->dc()-obj2->dc());
  
  if (m_rpMin<rp && rp<m_rpMax && m_piMin<pi && pi<m_piMax) {
    
    int ir = max(0, min(int((rp-m_rpMin)*m_binSize_inv_D1), m_nbins_D1));
    int jr = max(0, min(int((log10(pi)-log10(m_piMin))*m_binSize_inv_D2), m_nbins_D2));

    double angWeight = (m_angularWeight==nullptr) ? 1.
      : m_angularWeight(converted_angle(angular_distance(obj1->xx()/obj1->dc(), obj2->xx()/obj2->dc(), obj1->yy()/obj1->dc(), obj2->yy()/obj2->dc(), obj1->zz()/obj1->dc(), obj2->zz()/obj2->dc()), _radians_, m_angularUnits));
    
    m_PP2D[ir][jr] += obj1->weight()*obj2->weight()*angWeight;

    m_scale_D1_mean[ir][jr] += rp*obj1->weight()*obj2->weight()*angWeight;
    m_scale_D2_mean[ir][jr] += pi*obj1->weight()*obj2->weight()*angWeight;
    
    m_scale_D1_sigma[ir][jr] += pow(rp*obj1->weight()*obj2->weight(), 2);
    m_scale_D2_sigma[ir][jr] += pow(pi*obj1->weight()*obj2->weight(), 2);

    double pair_redshift = (obj1->redshift()>0 && obj2->redshift()>0) ? (obj1->redshift()+obj2->redshift())*0.5 : -1.;
    m_z_mean[ir][jr] += pair_redshift*obj1->weight()*obj2->weight()*angWeight;
    m_z_sigma[ir][jr] += pow(pair_redshift*obj1->weight()*obj2->weight(), 2);
    
  }
}


// ============================================================================================


void cosmobl::pairs::Pair2D_comovingCartesian_loglin_extra::put (const shared_ptr<Object> obj1, const shared_ptr<Object> obj2) 
{
  double rp = perpendicular_distance(obj1->ra(), obj2->ra(), obj1->dec(), obj2->dec(), obj1->dc(), obj2->dc());
  double pi = fabs(obj1->dc()-obj2->dc());

  if (m_rpMin<rp && rp<m_rpMax && m_piMin<pi && pi<m_piMax) {

    int ir = max(0, min(int((log10(rp)-log10(m_rpMin))*m_binSize_inv_D1), m_nbins_D1)); 
    int jr = max(0, min(int((pi-m_piMin)*m_binSize_inv_D2), m_nbins_D2));

    double angWeight = (m_angularWeight==nullptr) ? 1.
      : m_angularWeight(converted_angle(angular_distance(obj1->xx()/obj1->dc(), obj2->xx()/obj2->dc(), obj1->yy()/obj1->dc(), obj2->yy()/obj2->dc(), obj1->zz()/obj1->dc(), obj2->zz()/obj2->dc()), _radians_, m_angularUnits));
 
    m_PP2D[ir][jr] += obj1->weight()*obj2->weight()*angWeight;
    
    m_scale_D1_mean[ir][jr] += rp*obj1->weight()*obj2->weight()*angWeight;
    m_scale_D2_mean[ir][jr] += pi*obj1->weight()*obj2->weight()*angWeight;
    
    m_scale_D1_sigma[ir][jr] += pow(rp*obj1->weight()*obj2->weight(), 2);
    m_scale_D2_sigma[ir][jr] += pow(pi*obj1->weight()*obj2->weight(), 2);

    double pair_redshift = (obj1->redshift()>0 && obj2->redshift()>0) ? (obj1->redshift()+obj2->redshift())*0.5 : -1.;
    m_z_mean[ir][jr] += pair_redshift*obj1->weight()*obj2->weight()*angWeight;
    m_z_sigma[ir][jr] += pow(pair_redshift*obj1->weight()*obj2->weight(), 2);
    
  }
}


// ============================================================================================


void cosmobl::pairs::Pair2D_comovingCartesian_loglog_extra::put (const shared_ptr<Object> obj1, const shared_ptr<Object> obj2) 
{
  double rp = perpendicular_distance(obj1->ra(), obj2->ra(), obj1->dec(), obj2->dec(), obj1->dc(), obj2->dc());
  double pi = fabs(obj1->dc()-obj2->dc());
  
  if (m_rpMin<rp && rp<m_rpMax && m_piMin<pi && pi<m_piMax) {

    int ir = max(0, min(int((log10(rp)-log10(m_rpMin))*m_binSize_inv_D1), m_nbins_D1)); 
    int jr = max(0, min(int((log10(pi)-log10(m_piMin))*m_binSize_inv_D2), m_nbins_D2)); 

    double angWeight = (m_angularWeight==nullptr) ? 1.
      : m_angularWeight(converted_angle(angular_distance(obj1->xx()/obj1->dc(), obj2->xx()/obj2->dc(), obj1->yy()/obj1->dc(), obj2->yy()/obj2->dc(), obj1->zz()/obj1->dc(), obj2->zz()/obj2->dc()), _radians_, m_angularUnits));
 
    m_PP2D[ir][jr] += obj1->weight()*obj2->weight()*angWeight;

    m_scale_D1_mean[ir][jr] += rp*obj1->weight()*obj2->weight()*angWeight;
    m_scale_D2_mean[ir][jr] += pi*obj1->weight()*obj2->weight()*angWeight;
    
    m_scale_D1_sigma[ir][jr] += pow(rp*obj1->weight()*obj2->weight(), 2);
    m_scale_D2_sigma[ir][jr] += pow(pi*obj1->weight()*obj2->weight(), 2);

    double pair_redshift = (obj1->redshift()>0 && obj2->redshift()>0) ? (obj1->redshift()+obj2->redshift())*0.5 : -1.;
    m_z_mean[ir][jr] += pair_redshift*obj1->weight()*obj2->weight()*angWeight;
    m_z_sigma[ir][jr] += pow(pair_redshift*obj1->weight()*obj2->weight(), 2);
    
  }
}


// ============================================================================================


void cosmobl::pairs::Pair2D_comovingPolar_linlin_extra::put (const shared_ptr<Object> obj1, const shared_ptr<Object> obj2) 
{
  double rr = Euclidean_distance(obj1->xx(), obj2->xx(), obj1->yy(), obj2->yy(), obj1->zz(), obj2->zz());
  double mu = fabs(obj1->dc()-obj2->dc())/rr;
  
  if (m_rMin<rr && rr<m_rMax && m_muMin<mu && mu<m_muMax) {
    
    int ir = max(0, min(int((rr-m_rMin)*m_binSize_inv_D1), m_nbins_D1));
    int jr = max(0, min(int((mu-m_muMin)*m_binSize_inv_D2), m_nbins_D2));

    double angWeight = (m_angularWeight==nullptr) ? 1.
      : m_angularWeight(converted_angle(angular_distance(obj1->xx()/obj1->dc(), obj2->xx()/obj2->dc(), obj1->yy()/obj1->dc(), obj2->yy()/obj2->dc(), obj1->zz()/obj1->dc(), obj2->zz()/obj2->dc()), _radians_, m_angularUnits));
 
    m_PP2D[ir][jr] += obj1->weight()*obj2->weight()*angWeight;

    m_scale_D1_mean[ir][jr] += rr*obj1->weight()*obj2->weight()*angWeight;
    m_scale_D2_mean[ir][jr] += mu*obj1->weight()*obj2->weight()*angWeight;
    
    m_scale_D1_sigma[ir][jr] += pow(rr*obj1->weight()*obj2->weight(), 2);
    m_scale_D2_sigma[ir][jr] += pow(mu*obj1->weight()*obj2->weight(), 2);

    double pair_redshift = (obj1->redshift()>0 && obj2->redshift()>0) ? (obj1->redshift()+obj2->redshift())*0.5 : -1.;
    m_z_mean[ir][jr] += pair_redshift*obj1->weight()*obj2->weight()*angWeight;
    m_z_sigma[ir][jr] += pow(pair_redshift*obj1->weight()*obj2->weight(), 2);
    
  }
}


// ============================================================================================


void cosmobl::pairs::Pair2D_comovingPolar_linlog_extra::put (const shared_ptr<Object> obj1, const shared_ptr<Object> obj2) 
{
  double rr = Euclidean_distance(obj1->xx(), obj2->xx(), obj1->yy(), obj2->yy(), obj1->zz(), obj2->zz());
  double mu = fabs(obj1->dc()-obj2->dc())/rr;
  
  if (m_rMin<rr && rr<m_rMax && m_muMin<mu && mu<m_muMax) {
    
    int ir = max(0, min(int((rr-m_rMin)*m_binSize_inv_D1), m_nbins_D1));
    int jr = max(0, min(int((log10(mu)-log10(m_muMin))*m_binSize_inv_D2), m_nbins_D2)); 

    double angWeight = (m_angularWeight==nullptr) ? 1.
      : m_angularWeight(converted_angle(angular_distance(obj1->xx()/obj1->dc(), obj2->xx()/obj2->dc(), obj1->yy()/obj1->dc(), obj2->yy()/obj2->dc(), obj1->zz()/obj1->dc(), obj2->zz()/obj2->dc()), _radians_, m_angularUnits));
    
    m_PP2D[ir][jr] += obj1->weight()*obj2->weight()*angWeight;

    m_scale_D1_mean[ir][jr] += rr*obj1->weight()*obj2->weight()*angWeight;
    m_scale_D2_mean[ir][jr] += mu*obj1->weight()*obj2->weight()*angWeight;
    
    m_scale_D1_sigma[ir][jr] += pow(rr*obj1->weight()*obj2->weight(), 2);
    m_scale_D2_sigma[ir][jr] += pow(mu*obj1->weight()*obj2->weight(), 2);

    double pair_redshift = (obj1->redshift()>0 && obj2->redshift()>0) ? (obj1->redshift()+obj2->redshift())*0.5 : -1.;
    m_z_mean[ir][jr] += pair_redshift*obj1->weight()*obj2->weight()*angWeight;
    m_z_sigma[ir][jr] += pow(pair_redshift*obj1->weight()*obj2->weight(), 2);
    
  }
}


// ============================================================================================


void cosmobl::pairs::Pair2D_comovingPolar_loglin_extra::put (const shared_ptr<Object> obj1, const shared_ptr<Object> obj2) 
{
  double rr = Euclidean_distance(obj1->xx(), obj2->xx(), obj1->yy(), obj2->yy(), obj1->zz(), obj2->zz());
  double mu = fabs(obj1->dc()-obj2->dc())/rr;

  if (m_rMin<rr && rr<m_rMax && m_muMin<mu && mu<m_muMax) {

    int ir = max(0, min(int((log10(rr)-log10(m_rMin))*m_binSize_inv_D1), m_nbins_D1)); 
    int jr = max(0, min(int((mu-m_muMin)*m_binSize_inv_D2), m_nbins_D2));

    double angWeight = (m_angularWeight==nullptr) ? 1.
      : m_angularWeight(converted_angle(angular_distance(obj1->xx()/obj1->dc(), obj2->xx()/obj2->dc(), obj1->yy()/obj1->dc(), obj2->yy()/obj2->dc(), obj1->zz()/obj1->dc(), obj2->zz()/obj2->dc()), _radians_, m_angularUnits));
    
    m_PP2D[ir][jr] += obj1->weight()*obj2->weight()*angWeight;

    m_scale_D1_mean[ir][jr] += rr*obj1->weight()*obj2->weight()*angWeight;
    m_scale_D2_mean[ir][jr] += mu*obj1->weight()*obj2->weight()*angWeight;
    
    m_scale_D1_sigma[ir][jr] += pow(rr*obj1->weight()*obj2->weight(), 2);
    m_scale_D2_sigma[ir][jr] += pow(mu*obj1->weight()*obj2->weight(), 2);

    double pair_redshift = (obj1->redshift()>0 && obj2->redshift()>0) ? (obj1->redshift()+obj2->redshift())*0.5 : -1.;
    m_z_mean[ir][jr] += pair_redshift*obj1->weight()*obj2->weight()*angWeight;
    m_z_sigma[ir][jr] += pow(pair_redshift*obj1->weight()*obj2->weight(), 2);
    
  }
}


// ============================================================================================


void cosmobl::pairs::Pair2D_comovingPolar_loglog_extra::put (const shared_ptr<Object> obj1, const shared_ptr<Object> obj2) 
{
  double rr = Euclidean_distance(obj1->xx(), obj2->xx(), obj1->yy(), obj2->yy(), obj1->zz(), obj2->zz());
  double mu = fabs(obj1->dc()-obj2->dc())/rr;
  
  if (m_rMin<rr && rr<m_rMax && m_muMin<mu && mu<m_muMax) {

    int ir = max(0, min(int((log10(rr)-log10(m_rMin))*m_binSize_inv_D1), m_nbins_D1)); 
    int jr = max(0, min(int((log10(mu)-log10(m_muMin))*m_binSize_inv_D2), m_nbins_D2)); 

    double angWeight = (m_angularWeight==nullptr) ? 1.
      : m_angularWeight(converted_angle(angular_distance(obj1->xx()/obj1->dc(), obj2->xx()/obj2->dc(), obj1->yy()/obj1->dc(), obj2->yy()/obj2->dc(), obj1->zz()/obj1->dc(), obj2->zz()/obj2->dc()), _radians_, m_angularUnits));
    
    m_PP2D[ir][jr] += obj1->weight()*obj2->weight()*angWeight;

    m_scale_D1_mean[ir][jr] += rr*obj1->weight()*obj2->weight()*angWeight;
    m_scale_D2_mean[ir][jr] += mu*obj1->weight()*obj2->weight()*angWeight;
    
    m_scale_D1_sigma[ir][jr] += pow(rr*obj1->weight()*obj2->weight(), 2);
    m_scale_D2_sigma[ir][jr] += pow(mu*obj1->weight()*obj2->weight(), 2);

    double pair_redshift = (obj1->redshift()>0 && obj2->redshift()>0) ? (obj1->redshift()+obj2->redshift())*0.5 : -1.;
    m_z_mean[ir][jr] += pair_redshift*obj1->weight()*obj2->weight()*angWeight;
    m_z_sigma[ir][jr] += pow(pair_redshift*obj1->weight()*obj2->weight(), 2);
    
  }
}


// ============================================================================================


void cosmobl::pairs::Pair2D_extra::add_data2D (const int i, const int j, const vector<double> data)
{
  checkDim(m_PP2D, i, j, "m_PP2D", false);
  checkDim(m_scale_D1_mean, i, j, "m_scale_D1_mean", false);
  checkDim(m_scale_D2_mean, i, j, "m_scale_D2_mean", false);
  checkDim(m_scale_D1_sigma, i, j, "m_scale_D1_sigma", false);
  checkDim(m_scale_D2_sigma, i, j, "m_scale_D2_sigma", false);
  checkDim(m_z_mean, i, j, "m_z_mean", false);
  checkDim(m_z_sigma, i, j, "m_z_sigma", false);
  checkDim(data, 7, "data");
  
  m_PP2D[i][j] += data[0];
  m_scale_D1_mean[i][j] += data[1];
  m_scale_D1_sigma[i][j] += data[2];
  m_scale_D2_mean[i][j] += data[3];
  m_scale_D2_sigma[i][j] += data[4];
  m_z_mean[i][j] += data[5];
  m_z_sigma[i][j] += data[6];
}


// ============================================================================================


void cosmobl::pairs::Pair2D_extra::add_data2D (const int i, const int j, const shared_ptr<pairs::Pair> pair, const double ww)
{
  checkDim(m_PP2D, i, j, "m_PP2D", false);
  checkDim(m_scale_D1_mean, i, j, "m_scale_D1_mean", false);
  checkDim(m_scale_D2_mean, i, j, "m_scale_D2_mean", false);
  checkDim(m_scale_D1_sigma, i, j, "m_scale_D1_sigma", false);
  checkDim(m_scale_D2_sigma, i, j, "m_scale_D2_sigma", false); 
  checkDim(m_z_mean, i, j, "m_z_mean", false);
  checkDim(m_z_sigma, i, j, "m_z_sigma", false);
  
  m_PP2D[i][j] += ww*pair->PP2D(i, j);
  m_scale_D1_mean[i][j] += ww*pair->scale_D1_mean(i, j);
  m_scale_D1_sigma[i][j] += ww*pair->scale_D1_sigma(i, j);
  m_scale_D2_mean[i][j] += ww*pair->scale_D2_mean(i, j);
  m_scale_D2_sigma[i][j] += ww*pair->scale_D2_sigma(i, j);
  m_z_mean[i][j] += ww*pair->z_mean(i, j);
  m_z_sigma[i][j] += ww*pair->z_sigma(i, j);
}


// ============================================================================================


void cosmobl::pairs::Pair2D_extra::Sum (const shared_ptr<Pair> pp, const double ww)
{
  if (m_nbins_D1 != pp->nbins_D1() || m_nbins_D2 != pp->nbins_D2()) 
    ErrorCBL("Error in cosmobl::pairs::Pair2D::Sum of Pair.cpp: dimension problems!");
  
  for (int i=0; i<m_nbins_D1; i++) 
    for (int j=0; j<m_nbins_D2; j++) {
      m_PP2D[i][j] += ww*pp->PP2D(i, j);
      m_scale_D1_mean[i][j] += ww*pp->scale_D1_mean(i, j);
      m_scale_D1_sigma[i][j] += ww*pp->scale_D1_sigma(i, j);
      m_scale_D2_mean[i][j] += ww*pp->scale_D2_mean(i, j);
      m_scale_D2_sigma[i][j] += ww*pp->scale_D2_sigma(i, j);
      m_z_mean[i][j] += ww*pp->z_mean(i, j);
      m_z_sigma[i][j] += ww*pp->z_sigma(i, j);
    }
}

// ============================================================================================


void cosmobl::pairs::Pair2D_extra::finalise ()
{
  for (int i=0; i<m_nbins_D1; i++) 
    for (int j=0; j<m_nbins_D2; j++) {
      m_scale_D1_mean[i][j] = (m_PP2D[i][j]>0) ? m_scale_D1_mean[i][j]/m_PP2D[i][j] : 0.;
      m_scale_D1_sigma[i][j] = (m_PP2D[i][j]>0) ? sqrt(m_scale_D1_sigma[i][j]/m_PP2D[i][j]-pow(m_scale_D1_mean[i][j], 2)) : 0.;
      m_scale_D2_mean[i][j] = (m_PP2D[i][j]>0) ? m_scale_D2_mean[i][j]/m_PP2D[i][j] : 0.;
      m_scale_D2_sigma[i][j] = (m_PP2D[i][j]>0) ? sqrt(m_scale_D2_sigma[i][j]/m_PP2D[i][j]-pow(m_scale_D2_mean[i][j], 2)) : 0.;
      m_z_mean[i][j] = (m_PP2D[i][j]>0) ? m_z_mean[i][j]/m_PP2D[i][j] : 0.;
      m_z_sigma[i][j] = (m_PP2D[i][j]>0) ? sqrt(m_z_sigma[i][j]/m_PP2D[i][j]-pow(m_z_mean[i][j], 2)) : 0.;
    }
}

