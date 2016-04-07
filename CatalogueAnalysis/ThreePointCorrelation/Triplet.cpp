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
 *  @file CatalogueAnalysis/ThreePointCorrelation/Triplet.cpp
 *
 *  @brief Methods of the class Triplet  
 *
 *  This file contains the implementation of methods of the class Triplet,
 *  used to handle pairs of objects to compute the two-point correlation function
 *
 *  @authors Federico Marulli, Michele Moresco, Alfonso Veropalumbo 
 *
 *  @authors federico.marulli3@unbo.it, michele.moresco@unibo.it,
 *  alfonso.veropalumbo@unibo.it
 */

#include "Triplet1D.h"
#include "Triplet2D.h"

using namespace cosmobl;
using namespace triplets;


// ============================================================================================


shared_ptr<Triplet> cosmobl::triplets::Triplet::Create (const TripletType type, const double side_s, const int side_u, const double perc_increase, const int nbins)
{
  if (type==_comoving_theta_) return move(unique_ptr<Triplet1D_comoving_theta>{new Triplet1D_comoving_theta(side_s, side_u, perc_increase, nbins)});
  else if (type==_comoving_side_)  return move(unique_ptr<Triplet1D_comoving_side>{new Triplet1D_comoving_side(side_s, side_u, perc_increase, nbins)});
  
  else ErrorMsg("Error in cosmobl::triplets::Create of Triplet.cpp: no such type of object!");
  
  return NULL;
}


// ============================================================================================
// ============================================================================================


void cosmobl::triplets::Triplet1D::Sum (const shared_ptr<Triplet> tt, const double ww)  
{
  for (size_t i=0; i<m_TT1D.size(); i++) 
    m_TT1D[i] += ww*tt->TT1D(i);	 
}


// ============================================================================================
// ============================================================================================


void cosmobl::triplets::Triplet1D_comoving_side::set_parameters () 
{
  m_binSize = (m_side_s*(1+m_side_u)-m_side_s*(m_side_u-1))/m_nbins;

  m_scale.resize(m_nbins);
  for (int i=0; i<m_nbins; i++)
    m_scale[i] = m_side_s+(i+0.5)*m_binSize;
}


// ============================================================================================


void cosmobl::triplets::Triplet1D_comoving_side::put (const double r12, const double r13, const double r23, const double ww) 
{
  int klin = int((r23-m_side_s)/m_binSize);

  m_TT1D[klin] += ww;
}


// ============================================================================================


void cosmobl::triplets::Triplet1D_comoving_side::put (const shared_ptr<catalogue::Object> obj1, const shared_ptr<catalogue::Object> obj2, const shared_ptr<catalogue::Object> obj3)
{
  double x2 = obj2->xx(), y2 = obj2->yy(), z2 = obj2->zz(), w2 = obj2->weight();
  double x3 = obj3->xx(), y3 = obj3->yy(), z3 = obj3->zz(), w3 = obj3->weight();

  double r23 = cosmobl::Euclidean_distance(x2, x3, y2, y3, z2, z3); 
  
  double ww = w2*w3;

  int klin = int((r23-m_side_s)/m_binSize);

  m_TT1D[klin] += ww;
}


// ============================================================================================
// ============================================================================================


void cosmobl::triplets::Triplet1D_comoving_theta::set_parameters () 
{
  m_binSize = par::pi/m_nbins;

  m_scale.resize(m_nbins);
  for (int i=0; i<m_nbins; i++)
    m_scale[i] = (i+0.5)*m_binSize/par::pi;
}


// ============================================================================================


void cosmobl::triplets::Triplet1D_comoving_theta::put (const double r12, const double r13, const double r23, const double ww) 
{
  double shift_angle = ((((r12*r12)+(r13*r13)-(r23*r23))/(2*r12*r13))<0.) ? 0.00000001 : -0.00000001;
  double angle = (abs((((r12*r12)+(r13*r13)-(r23*r23))/(2*r12*r13)))>0.99999999) ? acos((((r12*r12)+(r13*r13)-(r23*r23))/(2*r12*r13))+shift_angle) : acos(((r12*r12)+(r13*r13)-(r23*r23))/(2*r12*r13));
  //angle = acos(((r12*r12)+(r13*r13)-(r23*r23))/(2*r12*r13));
  int klin = int(angle/m_binSize);

  m_TT1D[klin] += ww;
}


// ============================================================================================


void cosmobl::triplets::Triplet1D_comoving_theta::put (const shared_ptr<catalogue::Object> obj1, const shared_ptr<catalogue::Object> obj2, const shared_ptr<catalogue::Object> obj3) 
{
  double x1 = obj1->xx(), y1 = obj1->yy(), z1 = obj1->zz(), w1 = obj1->weight();
  double x2 = obj2->xx(), y2 = obj2->yy(), z2 = obj2->zz(), w2 = obj2->weight();
  double x3 = obj3->xx(), y3 = obj3->yy(), z3 = obj3->zz(), w3 = obj3->weight();

  double r12 = cosmobl::Euclidean_distance(x1, x2, y1, y2, z1, z2); 
  double r13 = cosmobl::Euclidean_distance(x1, x3, y1, y3, z1, z3); 
  double r23 = cosmobl::Euclidean_distance(x2, x3, y2, y3, z2, z3); 
  
  double ww = w1*w2*w3;

  double shift_angle = ((((r12*r12)+(r13*r13)-(r23*r23))/(2*r12*r13))<0.) ? 0.00000001 : -0.00000001;
  double angle = (abs((((r12*r12)+(r13*r13)-(r23*r23))/(2*r12*r13)))>0.99999999) ? acos((((r12*r12)+(r13*r13)-(r23*r23))/(2*r12*r13))+shift_angle) : acos(((r12*r12)+(r13*r13)-(r23*r23))/(2*r12*r13));
  //angle = acos(((r12*r12)+(r13*r13)-(r23*r23))/(2*r12*r13));
  int klin = int(angle/m_binSize);

  m_TT1D[klin] += ww;
}


