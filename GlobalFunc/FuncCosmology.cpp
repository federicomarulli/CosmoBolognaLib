/*******************************************************************
 *  Copyright (C) 2010 by Federico Marulli                         *
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
 *  @file GlobalFunc/FuncCosmology.cpp
 *
 *  @brief Generic functions that use the class Cosmology
 *
 *  This file contains the implementation of a set of generic
 *  functions that use the class Cosmology
 *
 *  @author Federico Marulli 
 *
 *  @author federico.marulli3@unibo.it
 */

#include "GlobalFunc.h"

using namespace std;

using namespace cbl;


// ============================================================================


void cbl::Vmax_DC_distribution (std::vector<double> &dc, std::vector<double> &nObj, const std::vector<double> D_C, const std::vector<double> zobj_min, const std::vector<double> zobj_max, const double z_min, const double z_max, const double zbin_min, const double zbin_max, cosmology::Cosmology &cosm, const double Area, const int nObjRan, const bool norm, const std::string file_Vmax, const double delta_dc_Vmax, const int seed)
{
  if (dc.size()>0 || nObj.size()>0) ErrorCBL("", "Vmax_DC_distribution", "GlobalFunc/FuncCosmology.cpp");

  random::UniformRandomNumbers ran(0., 1., seed);
  
  vector<double> err, dc_Vmax, ww;
  double Volume, zz;

  for (unsigned int i=0; i<D_C.size(); i++) { 
    for (int j=0; j<nObjRan; j++) {
      Volume = ran()*cosm.Volume(zobj_min[i], zobj_max[i], Area);
      zz = cosm.max_redshift(Volume, Area, zobj_min[i]);
      if (z_min<zz && zz<z_max) { 
	dc_Vmax.push_back(cosm.D_C(zz));
	ww.push_back(1.);
      }
    }
  }

  double fact = (norm) ? double(dc_Vmax.size())/double(D_C.size()) : dc_Vmax.size();
  
  double dc1 = cosm.D_C(zbin_min), dc2 = cosm.D_C(zbin_max);

  int nbin = nint((dc2-dc1))/delta_dc_Vmax;

  bool linear = 1;
  
  distribution(dc, nObj, err, dc_Vmax, ww, nbin, linear, file_Vmax, fact, dc1, dc2);
}


// ============================================================================================


double cbl::AP_shift_r (const double redshift, const cosmology::Cosmology &cosm1, const cosmology::Cosmology &cosm2)
{
  return cosm2.D_V(redshift)/cosm1.D_V(redshift);
}

double cbl::AP_shift_rp (const double redshift, const cosmology::Cosmology &cosm1, const cosmology::Cosmology &cosm2)
{
  return cosm1.D_A(redshift)/cosm2.D_A(redshift);
}

double cbl::AP_shift_pi (const double redshift, const cosmology::Cosmology &cosm1, const cosmology::Cosmology &cosm2)
{
  return cosm2.HH(redshift)/cosm1.HH(redshift);
}


// ============================================================================================


void cbl::max_separations_AP (const double Rp_max, const double Pi_max, const double redshift, const cosmology::Cosmology &cosm1, const std::vector<cosmology::Cosmology> &cosm2, double &rpM_AP, double &piM_AP, double &rM_AP) 
{
  vector<double> rp(cosm2.size()), pi(cosm2.size());

  for (unsigned int i=0; i<cosm2.size(); i++) {
    rp[i] = Rp_max*AP_shift_rp(redshift,cosm1,cosm2[i]);
    pi[i] = Pi_max*AP_shift_pi(redshift,cosm1,cosm2[i]);
  }

  rpM_AP = Max(rp);
  piM_AP = Max(pi);
  
  rM_AP = sqrt(pow(rpM_AP,2)+pow(piM_AP,2));
}


// ============================================================================================


double cbl::converted_xi (const double RR, const double redshift, const std::vector<double> rr, const std::vector<double> Xi, const cosmology::Cosmology &cosm1, const cosmology::Cosmology &cosm2, const bool direction) 
{
  if (RR==0) ErrorCBL("RR must be >0!", "converted_xi", "GlobalFunc/FuncCosmology.cpp");

  double gamma = AP_shift_r(redshift, cosm1, cosm2);

  // direction: 0:cosm2->cosm1, 1:cosm1->cosm2
  double lgRR = (direction) ? log10(RR/gamma) : log10(RR*gamma);

  vector<double> lgrr, lgXi;
  for (unsigned int i=0; i<rr.size(); i++) 
    if (rr[i]>0 && Xi[i]>0) {
      lgrr.push_back(log10(rr[i]));
      lgXi.push_back(log10(Xi[i]));
    }

  return pow(10., interpolated(lgRR, lgrr, lgXi, "Poly"));
}


// ============================================================================


double cbl::converted_xi (const double RP, const double PI, const double redshift, const std::vector<double> rp, const std::vector<double> pi, const std::vector<std::vector<double> > Xi, const cosmology::Cosmology &cosm1, const cosmology::Cosmology &cosm2, const bool direction) 
{
  double fDA = AP_shift_rp(redshift, cosm1, cosm2);
  double fH = AP_shift_pi(redshift, cosm1, cosm2);
  
  // direction: 0:cosm2->cosm1, 1:cosm1->cosm2
  double _RP = (direction) ? RP*fDA : RP/fDA;
  double _PI = (direction) ? PI*fH : PI/fH;

  return interpolated_2D(_RP, _PI, rp, pi, Xi, "Poly");
}
