/*******************************************************************
 *  Copyright (C) 2022 by Federico Marulli and Giorgio Lesci       *
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
 *  @file CAMBwrapper/CAMB.cpp
 *
 *  @brief functions that wrap CAMB routines
 *
 *  This file contains the implementation of 
 *  wrappers of CAMB routines
 *
 *  @author Giorgio Lesci
 *
 *  @author giorgio.lesci2@unibo.it
 */

#include "CAMB.h"

using namespace std;


// ============================================================================


std::vector<double> cbl::wrapper::camb::Pk_CAMB (const bool nonlinear, const double redshift, const double kmin, const double kmax, const int npoints, const double ombh2, const double omch2, const double omnuh2, const double massless_nu, const int massive_nu, const double omk, const double H0, const double ns, const double As, const double pivot_scalar, const double w, const double wa, const double tau, const bool accurate_massive_nu, const int neutrino_hierarchy, const int dark_energy_model, const double cs2_lam, const double T_cmb, const double helium_fraction)
{

  if (kmin < 5.e-5)
    ErrorCBL("the minimum k must be > 5.e-5.", "Pk_CAMB", "CAMB.cpp");
  
  void* params = GetCAMBparams();
  void* data = GetCAMBdata();

  // Set the cosmological parameters
  SetCAMBparams(params, ombh2, omch2, omnuh2, massless_nu, massive_nu, neutrino_hierarchy, omk, H0, dark_energy_model, w, wa, tau, cs2_lam, T_cmb, helium_fraction);
  SetCAMBPk(params, redshift, ns, As, pivot_scalar, accurate_massive_nu, kmax, nonlinear);

  // Get CAMB's P(k)
  GetCAMBresults(params, data);

  std::vector<double> Pk (npoints, 0.);
  const double h = (H0/100.);
  const double minkh = kmin/h;
  const double dlnkh = (log(kmax/h)-log(minkh)) / (npoints-1);

  GetCAMBPk(data, Pk.data(), minkh, dlnkh, npoints);

  // De-allocate
  ReleaseCAMBparams(params);
  ReleaseCAMBdata(data);

  return Pk;
  
}

