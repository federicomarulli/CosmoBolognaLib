/*******************************************************************
 *  Copyright (C) 2015 by Federico Marulli and Alfonso Veropalumbo *
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
 *  @file CatalogueAnalysis/TwoPointCorrelation/Pairs.cpp
 *
 *  @brief Methods of the class Pairs  
 *
 *  This file contains the implementation of methods of the class Pairs,
 *  used to handle pairs of objects to compute the two-point correlation function
 *
 *  @author Federico Marulli 
 *
 *  @author federico.marulli3@unbo.it
 */

#include "Pairs.h"
using namespace cosmobl;


// ============================================================================================


void cosmobl::Pairs2D::sum (shared_ptr<Pairs> pp, double ww) 
{
  if (m_nlog != pp->nlog() || m_nlin != pp->nlin()) 
    cosmobl::ErrorMsg("Error in cosmobl::Pairs2D::sum of Pairs.cpp: dimension problems!");

  for (int i=0; i<m_nlog; i++)
    m_PPlog[i] += ww*pp->PPlog(i);

  for (int i=0; i<m_nlin; i++) 
    m_PPlin[i] += ww*pp->PPlin(i);	 
}


// ============================================================================================


void cosmobl::Pairs2D::put (shared_ptr<Object> obj1, shared_ptr<Object> obj2) 
{
  double x1 = obj1->xx(), y1 = obj1->yy(), z1 = obj1->zz(), w1 = obj1->weight();
  double x2 = obj2->xx(), y2 = obj2->yy(), z2 = obj2->zz(), w2 = obj2->weight();

  double tt = cosmobl::angsep_xyz(x1-x2, y1-y2, z1-z2); 

  if (m_thetaMin <= tt && tt <= m_thetaMax) {

    double ww = w1*w2;

    // linear binning
    int klin = max(0,min(int((tt-m_thetaMin)*m_linbinsz_inv),m_nlin));

    // log binning
    int klog = min(int((Log(tt)-Log(m_thetaMin))*m_logbinsz_inv),m_nlog);

#pragma omp critical
    {
      if (klog>=0) m_PPlog[klog] += ww;

      m_PPlin[klin] += ww;
    }
  }
}


// ============================================================================================


void cosmobl::Pairs3D::sum (shared_ptr<Pairs> pp, double ww) 
{    
  if (m_nlog != pp->nlog() || m_nlin != pp->nlin() || m_ncos != pp->ncos()) 
    cosmobl::ErrorMsg("Error in cosmobl::Pairs3D::sum of Pairs.cpp: dimension problems!");

  
  // 1D
	 
  for (int i=0; i<m_nlog; i++)
    m_PPlog[i] += ww*pp->PPlog(i);
    
  for (int i=0; i<m_nlin; i++) 
    m_PPlin[i] += ww*pp->PPlin(i);
  

  // 2D
    
  for (int i=0; i<m_nlin; i++)
    for (int j=0; j<m_nlin; j++)
      m_PP2d[i][j] += ww*pp->PP2d(i,j);

  for (int i=0; i<m_nlog; i++)
    for (int j=0; j<m_nlin; j++)
      m_PPslog[i][j] += ww*pp->PPslog(i,j);

  for (int i=0; i<m_nlog; i++)
    for (int j=0; j<m_ncos; j++)
      m_PPcoslog[i][j] += ww*pp->PPcoslog(i,j);

  for (int i=0; i<m_nlin; i++)
    for (int j=0; j<m_ncos; j++)
      m_PPcoslin[i][j] += ww*pp->PPcoslin(i,j);
}
  

// ============================================================================================


void cosmobl::Pairs3D::put_log (shared_ptr<Object> obj1, shared_ptr<Object> obj2) 
{
  double x1 = obj1->xx(), y1 = obj1->yy(), z1 = obj1->zz(), w1 = obj1->weight(); 
  double x2 = obj2->xx(), y2 = obj2->yy(), z2 = obj2->zz(), w2 = obj2->weight(); 

  double rr = cosmobl::distance(x1-x2, y1-y2, z1-z2); 

  if (m_rMin<rr && rr<m_rMax) {

    double ww = w1*w2;
      
    // log binning
    int klog = max(0, min(int((Log(rr)-Log(m_rMin))*m_logbinsz_inv), int(m_nlog-1))); 

#pragma omp critical
    { 
      m_PPlog[klog] += ww;
    }
  }
}


// ============================================================================================


void cosmobl::Pairs3D::put (shared_ptr<Object> obj1, shared_ptr<Object> obj2) 
{
  double x1 = obj1->xx(), y1 = obj1->yy(), z1 = obj1->zz(), w1 = obj1->weight(), ra1 = obj1->ra(), dec1 = obj1->dec(), d1 = obj1->dc(); 
  double x2 = obj2->xx(), y2 = obj2->yy(), z2 = obj2->zz(), w2 = obj2->weight(), ra2 = obj2->ra(), dec2 = obj2->dec(), d2 = obj2->dc(); 

  double rr = cosmobl::distance(x1-x2, y1-y2, z1-z2); 

  if (m_rMin<rr && rr<m_rMax) {

    double ww = w1*w2;
    

    // ----- 1D ----- 

    // lin binning
    int klin = max(0, min(int((rr-m_rMin)*m_linbinsz_inv), int(m_nlin-1))); 
      
    // log binning
    int klog = max(0, min(int((Log(rr)-Log(m_rMin))*m_logbinsz_inv), int(m_nlog-1))); 
 

    // ----- 2D ----- 
	     
    double costheta = sin(dec1)*sin(dec2)+cos(dec1)*cos(dec2)*(cos(ra1)*cos(ra2)+sin(ra1)*sin(ra2));
    double theta = 0.;
    if (fabs(costheta)<1.-1.e-30) theta = acos(costheta);
    else if (costheta>=1.-1.e-30) theta = 0.;
    else theta = par::pi;
             
    double rp = (d1+d2)*tan(theta*0.5);
    rp *= 4.*d1*d2/((d1+d2)*(d1+d2));
    double pi = fabs(d1-d2);

    // 2D lin-lin binning
    int ir = max(0, min(int((rp-m_rMin)*m_linbinsz_inv), m_nlin));
    int jr = max(0, min(int((pi-m_rMin)*m_linbinsz_inv), m_nlin));
    int irlog = max(0, min(int((Log(rp)-Log(m_rMin))*m_logbinsz_inv), m_nlog)); 

    // polar 2D r,cos binning (log-lin)
    int clog = max(0, min(int((Log(sqrt(rp*rp+pi*pi))-Log(m_rMin))*m_logbinsz_inv), m_nlog));
    int ccos = max(0, min(int(pi/sqrt(rp*rp+pi*pi)*m_cosbinsz_inv), m_ncos));

    // polar 2D r,cos binning (lin-lin)
    int clin = max(0, min(int((sqrt(rp*rp+pi*pi)-m_rMin)*m_linbinsz_inv), m_nlin));

#pragma omp critical
    {
      m_PPlog[klog] += ww;
      m_PPlin[klin] += ww;

      m_PP2d[ir][jr] += ww;
      m_PPslog[irlog][jr] += ww;
      m_PPcoslog[clog][ccos] += ww;
      m_PPcoslin[clin][ccos] += ww;
    }
  }
}
