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
 *  @file Measure/ThreePointCorrelation/Triplet.cpp
 *
 *  @brief Methods of the class Triplet  
 *
 *  This file contains the implementation of methods of the class
 *  Triplet, used to handle pairs of objects to compute the two-point
 *  correlation function
 *
 *  @authors Federico Marulli, Michele Moresco, Alfonso Veropalumbo 
 *
 *  @authors federico.marulli3@unbo.it, michele.moresco@unibo.it,
 *  alfonso.veropalumbo@unibo.it
 */

#include "Triplet.h"
#include "Triplet1D.h"
#include "Triplet2D.h"

using namespace std;

using namespace cbl;
using namespace triplets;


// ============================================================================================


shared_ptr<Triplet> cbl::triplets::Triplet::Create (const cbl::triplets::TripletType type, const double r12, const double r12_binSize, const double r13, const double r13_binSize, const int nbins)
{
  if (type==TripletType::_comoving_theta_) return move(unique_ptr<Triplet1D_comoving_theta>{new Triplet1D_comoving_theta(r12, r12_binSize, r13, r13_binSize, nbins)});
  else if (type==TripletType::_comoving_side_)  return move(unique_ptr<Triplet1D_comoving_side>{new Triplet1D_comoving_side(r12, r12_binSize, r13, r13_binSize, nbins)});
  else if (type==TripletType::_comoving_costheta_)  return move(unique_ptr<Triplet1D_comoving_costheta>{new Triplet1D_comoving_costheta(r12, r12_binSize, r13, r13_binSize, nbins)});
  else if (type==TripletType::_multipoles_direct_)  return move(unique_ptr<Triplet1D_multipoles_direct>{new Triplet1D_multipoles_direct(r12, r12_binSize, r13, r13_binSize, nbins)});
  
  else ErrorCBL("Error in cbl::triplets::Create of Triplet.cpp: no such type of object!");
  
  return NULL;
}


// ============================================================================================
// ============================================================================================


void cbl::triplets::Triplet1D::Sum (const shared_ptr<Triplet> tt, const double ww)  
{
  for (size_t i=0; i<m_TT1D.size(); i++) 
    m_TT1D[i] += ww*tt->TT1D(i);	 
}


// ============================================================================================
// ============================================================================================


void cbl::triplets::Triplet1D_comoving_side::set_parameters () 
{
  double Lmin = (m_r13-0.5*m_r13_binSize)-(m_r12-0.5*m_r12_binSize);
  double Lmax = (m_r13+0.5*m_r13_binSize)+(m_r12+0.5*m_r12_binSize);
  m_binSize = (Lmax-Lmin)/m_nbins;

  m_scale.resize(m_nbins);
  for (int i=0; i<m_nbins; i++)
    m_scale[i] = (m_r12)*(m_r13-1)+((i+0.5)*m_binSize);
}


// ============================================================================================


void cbl::triplets::Triplet1D_comoving_side::get_triplet (const double r12, const double r13, const double r23, int &klin) 
{
  (void)r12; (void)r13;
  
  klin = int((r23-m_r12)/m_binSize);
}

// ============================================================================================


void cbl::triplets::Triplet1D_comoving_side::set_triplet (const int klin, const double ww) 
{
  m_TT1D[klin] += ww;
}

// ============================================================================================


void cbl::triplets::Triplet1D_comoving_side::put (const double r12, const double r13, const double r23, const double ww) 
{
  (void)r12; (void)r13;
  
  int klin = int((r23-m_r12)/m_binSize);

  m_TT1D[klin] += ww;
}


// ============================================================================================


void cbl::triplets::Triplet1D_comoving_side::put (const shared_ptr<catalogue::Object> obj1, const shared_ptr<catalogue::Object> obj2, const shared_ptr<catalogue::Object> obj3)
{
  (void)obj1;
  
  double x2 = obj2->xx(), y2 = obj2->yy(), z2 = obj2->zz(), w2 = obj2->weight();
  double x3 = obj3->xx(), y3 = obj3->yy(), z3 = obj3->zz(), w3 = obj3->weight();

  double r23 = cbl::Euclidean_distance(x2, x3, y2, y3, z2, z3); 
  
  double ww = w2*w3;

  int klin = int((r23-m_r12)/m_binSize);

  m_TT1D[klin] += ww;
}


// ============================================================================================
// ============================================================================================


void cbl::triplets::Triplet1D_comoving_theta::set_parameters () 
{
  m_binSize = par::pi/m_nbins;

  m_scale.resize(m_nbins);
  for (int i=0; i<m_nbins; i++)
    m_scale[i] = (i+0.5)*m_binSize/par::pi;
}


// ============================================================================================


void cbl::triplets::Triplet1D_comoving_theta::get_triplet (const double r12, const double r13, const double r23, int &klin) 
{
  double shift_angle = ((((r12*r12)+(r13*r13)-(r23*r23))/(2*r12*r13))<0.) ? 0.00000001 : -0.00000001;
  double angle = (abs((((r12*r12)+(r13*r13)-(r23*r23))/(2*r12*r13)))>0.99999999) ? acos((((r12*r12)+(r13*r13)-(r23*r23))/(2*r12*r13))+shift_angle) : acos(((r12*r12)+(r13*r13)-(r23*r23))/(2*r12*r13));
  //angle = acos(((r12*r12)+(r13*r13)-(r23*r23))/(2*r12*r13));
  klin = int(angle/m_binSize);

}


// ============================================================================================


void cbl::triplets::Triplet1D_comoving_theta::set_triplet (const int klin, const double ww) 
{
  m_TT1D[klin] += ww;
}


// ============================================================================================


void cbl::triplets::Triplet1D_comoving_theta::put (const double r12, const double r13, const double r23, const double ww) 
{
  double shift_angle = ((((r12*r12)+(r13*r13)-(r23*r23))/(2*r12*r13))<0.) ? 0.00000001 : -0.00000001;
  double angle = (abs((((r12*r12)+(r13*r13)-(r23*r23))/(2*r12*r13)))>0.99999999) ? acos((((r12*r12)+(r13*r13)-(r23*r23))/(2*r12*r13))+shift_angle) : acos(((r12*r12)+(r13*r13)-(r23*r23))/(2*r12*r13));
  //angle = acos(((r12*r12)+(r13*r13)-(r23*r23))/(2*r12*r13));
  int klin = int(angle/m_binSize);

  m_TT1D[klin] += ww;
}


// ============================================================================================


void cbl::triplets::Triplet1D_comoving_theta::put (const shared_ptr<catalogue::Object> obj1, const shared_ptr<catalogue::Object> obj2, const shared_ptr<catalogue::Object> obj3) 
{
  double x1 = obj1->xx(), y1 = obj1->yy(), z1 = obj1->zz(), w1 = obj1->weight();
  double x2 = obj2->xx(), y2 = obj2->yy(), z2 = obj2->zz(), w2 = obj2->weight();
  double x3 = obj3->xx(), y3 = obj3->yy(), z3 = obj3->zz(), w3 = obj3->weight();

  double r12 = cbl::Euclidean_distance(x1, x2, y1, y2, z1, z2); 
  double r13 = cbl::Euclidean_distance(x1, x3, y1, y3, z1, z3); 
  double r23 = cbl::Euclidean_distance(x2, x3, y2, y3, z2, z3); 
  
  double ww = w1*w2*w3;

  double shift_angle = ((((r12*r12)+(r13*r13)-(r23*r23))/(2*r12*r13))<0.) ? 0.00000001 : -0.00000001;
  double angle = (abs((((r12*r12)+(r13*r13)-(r23*r23))/(2*r12*r13)))>0.99999999) ? acos((((r12*r12)+(r13*r13)-(r23*r23))/(2*r12*r13))+shift_angle) : acos(((r12*r12)+(r13*r13)-(r23*r23))/(2*r12*r13));
  //angle = acos(((r12*r12)+(r13*r13)-(r23*r23))/(2*r12*r13));
  int klin = int(angle/m_binSize);

  m_TT1D[klin] += ww;
}


// ============================================================================================
// ============================================================================================


void cbl::triplets::Triplet1D_comoving_costheta::set_parameters () 
{
  m_binSize = 2./m_nbins;

  m_scale.resize(m_nbins);
  for (int i=0; i<m_nbins; i++)
    m_scale[i] = (i+0.5)*m_binSize-1;
}


// ============================================================================================


void cbl::triplets::Triplet1D_comoving_costheta::get_triplet (const double r12, const double r13, const double r23, int &klin) 
{
  double shift_angle = ((((r12*r12)+(r13*r13)-(r23*r23))/(2*r12*r13))<0.) ? 0.00000001 : -0.00000001;
  double cos_angle = (abs((((r12*r12)+(r13*r13)-(r23*r23))/(2*r12*r13)))>0.99999999) ? (((r12*r12)+(r13*r13)-(r23*r23))/(2*r12*r13))+shift_angle : ((r12*r12)+(r13*r13)-(r23*r23))/(2*r12*r13);
  //double cos_angle = ((r12*r12)+(r13*r13)-(r23*r23))/(2*r12*r13);
  klin = int((cos_angle+1)/m_binSize);

}


// ============================================================================================


void cbl::triplets::Triplet1D_comoving_costheta::set_triplet (const int klin, const double ww) 
{
  m_TT1D[klin] += ww;
}


// ============================================================================================


void cbl::triplets::Triplet1D_comoving_costheta::put (const double r12, const double r13, const double r23, const double ww) 
{
  double shift_angle = ((((r12*r12)+(r13*r13)-(r23*r23))/(2*r12*r13))<0.) ? 0.00000001 : -0.00000001;
  double cos_angle = (abs((((r12*r12)+(r13*r13)-(r23*r23))/(2*r12*r13)))>0.99999999) ? (((r12*r12)+(r13*r13)-(r23*r23))/(2*r12*r13))+shift_angle : ((r12*r12)+(r13*r13)-(r23*r23))/(2*r12*r13);
  //double cos_angle = ((r12*r12)+(r13*r13)-(r23*r23))/(2*r12*r13);
  int klin = int((cos_angle+1)/m_binSize);

  m_TT1D[klin] += ww;
}


// ============================================================================================


void cbl::triplets::Triplet1D_comoving_costheta::put (const shared_ptr<catalogue::Object> obj1, const shared_ptr<catalogue::Object> obj2, const shared_ptr<catalogue::Object> obj3) 
{
  double x1 = obj1->xx(), y1 = obj1->yy(), z1 = obj1->zz(), w1 = obj1->weight();
  double x2 = obj2->xx(), y2 = obj2->yy(), z2 = obj2->zz(), w2 = obj2->weight();
  double x3 = obj3->xx(), y3 = obj3->yy(), z3 = obj3->zz(), w3 = obj3->weight();

  double r12 = cbl::Euclidean_distance(x1, x2, y1, y2, z1, z2); 
  double r13 = cbl::Euclidean_distance(x1, x3, y1, y3, z1, z3); 
  double r23 = cbl::Euclidean_distance(x2, x3, y2, y3, z2, z3); 
  
  double ww = w1*w2*w3;

  double shift_angle = ((((r12*r12)+(r13*r13)-(r23*r23))/(2*r12*r13))<0.) ? 0.00000001 : -0.00000001;
  double cos_angle = (abs((((r12*r12)+(r13*r13)-(r23*r23))/(2*r12*r13)))>0.99999999) ? (((r12*r12)+(r13*r13)-(r23*r23))/(2*r12*r13))+shift_angle : ((r12*r12)+(r13*r13)-(r23*r23))/(2*r12*r13);
  //double cos_angle = ((r12*r12)+(r13*r13)-(r23*r23))/(2*r12*r13);
  int klin = int((cos_angle+1)/m_binSize);

  m_TT1D[klin] += ww;
}


// ============================================================================================


void cbl::triplets::Triplet1D_multipoles_direct::set_parameters () 
{
  m_scale.resize(m_nbins);
  for (int i=0; i<m_nbins; i++)
    m_scale[i] = i;
}


// ============================================================================================


void cbl::triplets::Triplet1D_multipoles_direct::put (const double r12, const double r13, const double r23, const double ww) 
{
  //double shift_angle = ((((r12*r12)+(r13*r13)-(r23*r23))/(2*r12*r13))<0.) ? 0.00000001 : -0.00000001;
  //double cos_angle = (abs((((r12*r12)+(r13*r13)-(r23*r23))/(2*r12*r13)))>0.99999999) ? (((r12*r12)+(r13*r13)-(r23*r23))/(2*r12*r13))+shift_angle : ((r12*r12)+(r13*r13)-(r23*r23))/(2*r12*r13);
  double cos_angle = ((r12*r12)+(r13*r13)-(r23*r23))/(2*r12*r13);

  for(int i=0; i<m_nbins; i++)
    m_TT1D[i] += ww*legendre_polynomial(cos_angle, i);

}


// ============================================================================================


void cbl::triplets::Triplet1D_multipoles_direct::put (const shared_ptr<catalogue::Object> obj1, const shared_ptr<catalogue::Object> obj2, const shared_ptr<catalogue::Object> obj3) 
{
  double x1 = obj1->xx(), y1 = obj1->yy(), z1 = obj1->zz(), w1 = obj1->weight();
  double x2 = obj2->xx(), y2 = obj2->yy(), z2 = obj2->zz(), w2 = obj2->weight();
  double x3 = obj3->xx(), y3 = obj3->yy(), z3 = obj3->zz(), w3 = obj3->weight();

  double r12 = cbl::Euclidean_distance(x1, x2, y1, y2, z1, z2); 
  double r13 = cbl::Euclidean_distance(x1, x3, y1, y3, z1, z3); 
  double r23 = cbl::Euclidean_distance(x2, x3, y2, y3, z2, z3); 
  
  double ww = w1*w2*w3;

  //double shift_angle = ((((r12*r12)+(r13*r13)-(r23*r23))/(2*r12*r13))<0.) ? 0.00000001 : -0.00000001;
  //double cos_angle = (abs((((r12*r12)+(r13*r13)-(r23*r23))/(2*r12*r13)))>0.99999999) ? (((r12*r12)+(r13*r13)-(r23*r23))/(2*r12*r13))+shift_angle : ((r12*r12)+(r13*r13)-(r23*r23))/(2*r12*r13);

  double cos_angle = ((r12*r12)+(r13*r13)-(r23*r23))/(2*r12*r13);

  for(int i=0; i<m_nbins; i++) 
    m_TT1D[i] += ww*legendre_polynomial(cos_angle, i);
  
}
