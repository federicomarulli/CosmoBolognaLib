/********************************************************************
 *  Copyright (C) 2014 by Federico Marulli and Alfonso Veropalumbo  *
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
 *  @file Prior/Prior.cpp
 *
 *  @brief Methods of the class Prior
 *
 *  This file contains the implementation of the methods of the class
 *  Prior
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#include "Prior.h"
using namespace cosmobl;


// ======================================================================================


cosmobl::Prior::Prior (vector<double> _Limits) 
{
  xmin = _Limits[0]; 
  xmax = _Limits[1];
  func = &identity;  
}


// =====================================================================================


cosmobl::Prior::Prior (vector<double> _Limits, double &mean, double &sigma) 
{
  xmin = _Limits[0]; 
  xmax = _Limits[1];
  prior_func_pars.push_back(mean);
  prior_func_pars.push_back(sigma);
  func = &gaussian;
}


// =====================================================================================


cosmobl::Prior::Prior (vector<double> _Limits, prior_func _func, vector<double> _prior_func_pars) 
: func(_func), prior_func_pars(_prior_func_pars) 
{
   xmin = _Limits[0]; 
   xmax = _Limits[1];
}


// =====================================================================================


void cosmobl::Prior::put_limits (vector<double> _Limits) {
  xmin = _Limits[0]; 
  xmax = _Limits[1];
}


// =====================================================================================


void cosmobl::Prior::put_gaussian_parameters (double &mean, double &sigma) 
{
  func = &gaussian;
  prior_func_pars.erase(prior_func_pars.begin(), prior_func_pars.end());
  prior_func_pars.push_back(mean);
  prior_func_pars.push_back(sigma);
}


// =====================================================================================


void cosmobl::Prior::put_func_parameters (prior_func _func, vector<double> pars) 
{
   func = _func;
   prior_func_pars.erase(prior_func_pars.begin(),prior_func_pars.end());

   for (unsigned int i=0; i<pars.size(); i++)
      prior_func_pars.push_back(pars[i]);
}
