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
 *  @file Wrappers/CUBAwrapper.cpp
 *
 *  @brief functions that wrap CUBA routines for 
 *  multidimensional integration,
 *  root finding and minimization
 *
 *  This file contains the implementation of 
 *  wrappers of  CUBA routines for multidimensional
 *  integration
 *
 *  @author Alfonso Veropalumbo
 *
 *  @author alfonso.veropalumbo@unbo.it
 */

#include "CUBAwrapper.h"

using namespace std;


// ============================================================================


int cbl::wrapper::cuba::CUBAIntegrand (const int *ndim, const cubareal xx[], const int *ncomp, cubareal ff[], void *userdata)
{
  cbl::wrapper::cuba::STR_CUBA_integrand *pp = static_cast<cbl::wrapper::cuba::STR_CUBA_integrand *>(userdata);

  vector<double> var;
  double fact = 1.;
  for (int i=0; i<*ndim; i++) {
    var.push_back(xx[i]*(pp->integration_limits[i][1]-pp->integration_limits[i][0])+pp->integration_limits[i][0]);
    fact *= (pp->integration_limits[i][1]-pp->integration_limits[i][0]);
  }

  ff[*ncomp-1] = pp->func(var)*fact;

  return 0;
}


// ============================================================================


cbl::wrapper::cuba::CUBAwrapper::CUBAwrapper (FunctionDoubleVectorPtrVectorRef func, const std::shared_ptr<void> function_parameters, std::vector<double> &parameters, const int ndim)
{
  set_integrand(func, function_parameters, parameters, ndim);
}


// ============================================================================


cbl::wrapper::cuba::CUBAwrapper::CUBAwrapper (FunctionDoubleVector func, const int ndim)
{
  set_integrand(func, ndim);
}


// ============================================================================


void cbl::wrapper::cuba::CUBAwrapper::set_integrand (FunctionDoubleVectorPtrVectorRef func, const std::shared_ptr<void> function_parameters, std::vector<double> &parameters, const int ndim)
{
  m_integrand = std::bind(func, std::placeholders::_1, function_parameters, parameters);
  m_ndim = ndim;
}


// ============================================================================


void cbl::wrapper::cuba::CUBAwrapper::set_integrand (FunctionDoubleVector func, const int ndim)
{
  m_integrand = func;
  m_ndim = ndim;
}


// ============================================================================


double cbl::wrapper::cuba::CUBAwrapper::IntegrateVegas (vector<vector<double>> integration_limits)
{
  int comp, neval, fail;
  cubareal integral[NCOMP], error[NCOMP], prob[NCOMP];

  cbl::wrapper::cuba::STR_CUBA_integrand *userdata = new cbl::wrapper::cuba::STR_CUBA_integrand;
  userdata->func = m_integrand;
  userdata->integration_limits = integration_limits;

  Vegas(m_ndim, NCOMP, cbl::wrapper::cuba::CUBAIntegrand, userdata, NVEC,
    EPSREL, EPSABS, VERBOSE, SEED,
    MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH,
    GRIDNO, STATEFILE, SPIN,
    &neval, &fail, integral, error, prob);

  if(VERBOSE>0){
    printf("VEGAS RESULT:\tneval %d\tfail %d\n", neval, fail);
    for( comp = 0; comp < NCOMP; ++comp )
      printf("VEGAS RESULT:\t%.8f +- %.8f\tp = %.3f\n", (double)integral[comp], (double)error[comp], (double)prob[comp]);
  }

  return integral[0];
}


// ============================================================================


double cbl::wrapper::cuba::CUBAwrapper::IntegrateSuave (vector<vector<double>> integration_limits)
{
  int comp, neval, fail, nregions;
  cubareal integral[NCOMP], error[NCOMP], prob[NCOMP];

  cbl::wrapper::cuba::STR_CUBA_integrand *userdata = new cbl::wrapper::cuba::STR_CUBA_integrand;
  userdata->func = m_integrand;
  userdata->integration_limits = integration_limits;
  Suave(m_ndim, NCOMP, cbl::wrapper::cuba::CUBAIntegrand, userdata, NVEC,
    EPSREL, EPSABS, VERBOSE | LAST, SEED,
    MINEVAL, MAXEVAL, NNEW, NMIN, FLATNESS,
    STATEFILE, SPIN,
    &nregions, &neval, &fail, integral, error, prob);

  if(VERBOSE>0){
    printf("SUAVE RESULT:\tnregions %d\tneval %d\tfail %d\n", nregions, neval, fail);
    for( comp = 0; comp < NCOMP; ++comp )
      printf("SUAVE RESULT:\t%.8f +- %.8f\tp = %.3f\n", (double)integral[comp], (double)error[comp], (double)prob[comp]);
  }

  return integral[0];
}


// ============================================================================


double cbl::wrapper::cuba::CUBAwrapper::IntegrateDivonne (vector<vector<double>> integration_limits)
{
  int comp, neval, fail, nregions;
  cubareal integral[NCOMP], error[NCOMP], prob[NCOMP];

  cbl::wrapper::cuba::STR_CUBA_integrand *userdata = new cbl::wrapper::cuba::STR_CUBA_integrand;
  userdata->func = m_integrand;
  userdata->integration_limits = integration_limits;

  Divonne(m_ndim, NCOMP, cbl::wrapper::cuba::CUBAIntegrand, userdata, NVEC,
    EPSREL, EPSABS, VERBOSE, SEED,
    MINEVAL, MAXEVAL, KEY1, KEY2, KEY3, MAXPASS,
    BORDER, MAXCHISQ, MINDEVIATION,
    NGIVEN, LDXGIVEN, NULL, NEXTRA, NULL,
    STATEFILE, SPIN,
    &nregions, &neval, &fail, integral, error, prob);

  if(VERBOSE>0){
    printf("DIVONNE RESULT:\tnregions %d\tneval %d\tfail %d\n", nregions, neval, fail);
    for( comp = 0; comp < NCOMP; ++comp ) 
      printf("DIVONNE RESULT:\t%.8f +- %.8f\tp = %.3f\n", (double)integral[comp], (double)error[comp], (double)prob[comp]);
  }

  return integral[0];
}


// ============================================================================


double cbl::wrapper::cuba::CUBAwrapper::IntegrateCuhre (vector<vector<double>> integration_limits)
{
  int comp, neval, fail, nregions;
  cubareal integral[NCOMP], error[NCOMP], prob[NCOMP];

  cbl::wrapper::cuba::STR_CUBA_integrand *userdata = new cbl::wrapper::cuba::STR_CUBA_integrand;
  userdata->func = m_integrand;
  userdata->integration_limits = integration_limits;

  Cuhre(m_ndim, NCOMP, cbl::wrapper::cuba::CUBAIntegrand, userdata, NVEC,
    EPSREL, EPSABS, VERBOSE | LAST,
    MINEVAL, MAXEVAL, KEY,
    STATEFILE, SPIN,
    &nregions, &neval, &fail, integral, error, prob);

  if(VERBOSE>0){
    printf("CUHRE RESULT:\tnregions %d\tneval %d\tfail %d\n", nregions, neval, fail);
    for( comp = 0; comp < NCOMP; ++comp )
      printf("CUHRE RESULT:\t%.8f +- %.8f\tp = %.3f\n", (double)integral[comp], (double)error[comp], (double)prob[comp]);
  }

  return integral[0];
}
