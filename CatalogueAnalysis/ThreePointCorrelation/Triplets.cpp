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
 *  @file CatalogueAnalysis/ThreePointCorrelation/Triplets.cpp
 *
 *  @brief Methods of the class Triplets  
 *
 *  This file contains the implementation of methods of the class Triplets,
 *  used to handle pairs of objects to compute the two-point correlation function
 *
 *  @authors Federico Marulli, Michele Moresco, Alfonso Veropalumbo 
 *
 *  @authors federico.marulli3@unbo.it, michele.moresco@unibo.it,
 *  alfonso.veropalumbo@unibo.it
 */

#include "Triplets.h"
using namespace cosmobl;


// ============================================================================================


void cosmobl::Triplets2D::sum (shared_ptr<Triplets> tt, double ww) 
{
  if (fabs(m_binsize-tt->binsize())>1.e-30)
    cosmobl::ErrorMsg("Error in cosmobl::Triplets2D::sum of Triplets.cpp: dimension problems!");
    
  for (size_t i=0; i<m_TT.size(); i++) 
    m_TT[i] += ww*tt->TT(i);	 
}


// ============================================================================================


void cosmobl::Triplets2D::put (double &r12, double &r13, double &r23, double ww) 
{
  int klin = -1;
  
  if (m_type_binning=="ang") {
    double shift_angle = ((((r12*r12)+(r13*r13)-(r23*r23))/(2*r12*r13))<0.) ? 0.00000001 : -0.00000001;
    double angle = (abs((((r12*r12)+(r13*r13)-(r23*r23))/(2*r12*r13)))>0.99999999) ? acos((((r12*r12)+(r13*r13)-(r23*r23))/(2*r12*r13))+shift_angle) : acos(((r12*r12)+(r13*r13)-(r23*r23))/(2*r12*r13));
    //angle = acos(((r12*r12)+(r13*r13)-(r23*r23))/(2*r12*r13));
    klin = int(angle/m_binsize);
  }

  else if (m_type_binning=="lin") 
    klin = int((r23-m_side_s)/m_binsize);

#pragma omp critical
  {
    m_TT[klin] += ww;
  }
}


// ============================================================================================


void cosmobl::Triplets2D::put (shared_ptr<Object> obj1, shared_ptr<Object> obj2, shared_ptr<Object> obj3) 
{
  double x1 = obj1->xx(), y1 = obj1->yy(), z1 = obj1->zz(), w1 = obj1->weight();
  double x2 = obj2->xx(), y2 = obj2->yy(), z2 = obj2->zz(), w2 = obj2->weight();
  double x3 = obj3->xx(), y3 = obj3->yy(), z3 = obj3->zz(), w3 = obj3->weight();

  double r12 = cosmobl::angsep_xyz(x1-x2, y1-y2, z1-z2); 
  double r13 = cosmobl::angsep_xyz(x1-x3, y1-y3, z1-z3); 
  double r23 = cosmobl::angsep_xyz(x2-x3, y2-y3, z2-z3); 
  
  double ww = w1*w2*w3;

  int klin = -1;
  
  if (m_type_binning=="ang") {
    double shift_angle = ((((r12*r12)+(r13*r13)-(r23*r23))/(2*r12*r13))<0.) ? 0.00000001 : -0.00000001;
    double angle = (abs((((r12*r12)+(r13*r13)-(r23*r23))/(2*r12*r13)))>0.99999999) ? acos((((r12*r12)+(r13*r13)-(r23*r23))/(2*r12*r13))+shift_angle) : acos(((r12*r12)+(r13*r13)-(r23*r23))/(2*r12*r13));
    //angle = acos(((r12*r12)+(r13*r13)-(r23*r23))/(2*r12*r13));
    klin = int(angle/m_binsize);
  }

  else if (m_type_binning=="lin") 
    klin = int((r23-m_side_s)/m_binsize);

#pragma omp critical
  {
    m_TT[klin] += ww;
  }
}


// ============================================================================================


void cosmobl::Triplets3D::sum (shared_ptr<Triplets> tt, double ww) 
{
  if (fabs(m_binsize-tt->binsize())>1.e-30) 
    cosmobl::ErrorMsg("Error in cosmobl::Triplets3D::sum of Triplets.cpp: dimension problems!");
    
  for (size_t i=0; i<m_TT.size(); i++) 
    m_TT[i] += ww*tt->TT(i);	 
}



// ============================================================================================


void cosmobl::Triplets3D::put (double &r12, double &r13, double &r23, double ww) 
{
  int klin = -1;

  if (m_type_binning=="ang") {
    double shift_angle = ((((r12*r12)+(r13*r13)-(r23*r23))/(2*r12*r13))<0.) ? 0.00000001 : -0.00000001;
    double angle = (abs((((r12*r12)+(r13*r13)-(r23*r23))/(2*r12*r13)))>0.99999999) ? acos((((r12*r12)+(r13*r13)-(r23*r23))/(2*r12*r13))+shift_angle) : acos(((r12*r12)+(r13*r13)-(r23*r23))/(2*r12*r13));
    //angle = acos(((r12*r12)+(r13*r13)-(r23*r23))/(2*r12*r13));
    klin = int(angle/m_binsize);
  }

  else if (m_type_binning=="lin") 
    klin = int((r23-m_side_s)/m_binsize);

#pragma omp critical
  {
    m_TT[klin] += ww;
  }
}


// ============================================================================================


void cosmobl::Triplets3D::put (shared_ptr<Object> obj1, shared_ptr<Object> obj2, shared_ptr<Object> obj3) 
{
  double x1 = obj1->xx(), y1 = obj1->yy(), z1 = obj1->zz(), w1 = obj1->weight();
  double x2 = obj2->xx(), y2 = obj2->yy(), z2 = obj2->zz(), w2 = obj2->weight();
  double x3 = obj3->xx(), y3 = obj3->yy(), z3 = obj3->zz(), w3 = obj3->weight();

  double r12 = cosmobl::distance(x1-x2, y1-y2, z1-z2); 
  double r13 = cosmobl::distance(x1-x3, y1-y3, z1-z3); 
  double r23 = cosmobl::distance(x2-x3, y2-y3, z2-z3); 
  
  double ww = w1*w2*w3;

  int klin = -1;
  
  if (m_type_binning=="ang") {
    double shift_angle = ((((r12*r12)+(r13*r13)-(r23*r23))/(2*r12*r13))<0.) ? 0.00000001 : -0.00000001;
    double angle = (abs((((r12*r12)+(r13*r13)-(r23*r23))/(2*r12*r13)))>0.99999999) ? acos((((r12*r12)+(r13*r13)-(r23*r23))/(2*r12*r13))+shift_angle) : acos(((r12*r12)+(r13*r13)-(r23*r23))/(2*r12*r13));
    //angle = acos(((r12*r12)+(r13*r13)-(r23*r23))/(2*r12*r13));
    klin = int(angle/m_binsize);
  }
  else if (m_type_binning=="lin") 
    klin = int((r23-m_side_s)/m_binsize);

#pragma omp critical
  {
    m_TT[klin] += ww;
  }
}
