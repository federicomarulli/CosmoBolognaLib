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
 *  @author alfonso.veropalumbo@unibo.it
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


double cbl::wrapper::cuba::CUBAwrapper::IntegrateVegas (vector<vector<double>> integration_limits, const bool parallelize)
{
  int comp, neval, fail;
  vector<double> integral(m_inputs.NCOMP, 0), error(m_inputs.NCOMP, 0), prob(m_inputs.NCOMP, 0);

  cbl::wrapper::cuba::STR_CUBA_integrand *userdata = new cbl::wrapper::cuba::STR_CUBA_integrand;
  userdata->func = m_integrand;
  userdata->integration_limits = integration_limits;

  if (!parallelize) {
    int zero = 0;
    cubacores(&zero, &zero);
  }

  Vegas(m_ndim, m_inputs.NCOMP, cbl::wrapper::cuba::CUBAIntegrand, userdata, m_inputs.NVEC,
    m_inputs.EPSREL, m_inputs.EPSABS, m_inputs.VERBOSE, m_inputs.SEED,
    m_inputs.MINEVAL, m_inputs.MAXEVAL, m_inputs.NSTART, m_inputs.NINCREASE, m_inputs.NBATCH,
    m_inputs.GRIDNO, m_inputs.STATEFILE, m_inputs.SPIN.get(),
    &neval, &fail, integral.data(), error.data(), prob.data());

  if(m_inputs.VERBOSE>0){
    printf("VEGAS RESULT:\tneval %d\tfail %d\n", neval, fail);
    for( comp = 0; comp < m_inputs.NCOMP; ++comp )
      printf("VEGAS RESULT:\t%.8f +- %.8f\tp = %.3f\n", (double)integral[comp], (double)error[comp], (double)prob[comp]);
  }

  return integral[0];
}


// ============================================================================


double cbl::wrapper::cuba::CUBAwrapper::IntegrateSuave (vector<vector<double>> integration_limits, const bool parallelize)
{
  int comp, neval, fail, nregions;
  vector<double> integral(m_inputs.NCOMP, 0), error(m_inputs.NCOMP, 0), prob(m_inputs.NCOMP, 0);

  cbl::wrapper::cuba::STR_CUBA_integrand *userdata = new cbl::wrapper::cuba::STR_CUBA_integrand;
  userdata->func = m_integrand;
  userdata->integration_limits = integration_limits;

  if (!parallelize) {
    int zero = 0;
    cubacores(&zero, &zero);
  }
  
  Suave(m_ndim, m_inputs.NCOMP, cbl::wrapper::cuba::CUBAIntegrand, userdata, m_inputs.NVEC,
    m_inputs.EPSREL, m_inputs.EPSABS, m_inputs.VERBOSE | m_inputs.LAST, m_inputs.SEED,
    m_inputs.MINEVAL, m_inputs.MAXEVAL, m_inputs.NNEW, m_inputs.NMIN, m_inputs.FLATNESS,
    m_inputs.STATEFILE, m_inputs.SPIN.get(),
    &nregions, &neval, &fail, integral.data(), error.data(), prob.data());

  if(m_inputs.VERBOSE>0){
    printf("SUAVE RESULT:\tnregions %d\tneval %d\tfail %d\n", nregions, neval, fail);
    for( comp = 0; comp < m_inputs.NCOMP; ++comp )
      printf("SUAVE RESULT:\t%.8f +- %.8f\tp = %.3f\n", (double)integral[comp], (double)error[comp], (double)prob[comp]);
  }

  return integral[0];
}


// ============================================================================


double cbl::wrapper::cuba::CUBAwrapper::IntegrateDivonne (vector<vector<double>> integration_limits, const bool parallelize)
{
  int comp, neval, fail, nregions;
  vector<double> integral(m_inputs.NCOMP, 0), error(m_inputs.NCOMP, 0), prob(m_inputs.NCOMP, 0);

  cbl::wrapper::cuba::STR_CUBA_integrand *userdata = new cbl::wrapper::cuba::STR_CUBA_integrand;
  userdata->func = m_integrand;
  userdata->integration_limits = integration_limits;

  if (!parallelize) {
    int zero = 0;
    cubacores(&zero, &zero);
  }

  Divonne(m_ndim, m_inputs.NCOMP, cbl::wrapper::cuba::CUBAIntegrand, userdata, m_inputs.NVEC,
    m_inputs.EPSREL, m_inputs.EPSABS, m_inputs.VERBOSE, m_inputs.SEED,
    m_inputs.MINEVAL, m_inputs.MAXEVAL, m_inputs.KEY1, m_inputs.KEY2, m_inputs.KEY3, m_inputs.MAXPASS,
    m_inputs.BORDER, m_inputs.MAXCHISQ, m_inputs.MINDEVIATION,
    m_inputs.NGIVEN, m_inputs.LDXGIVEN, NULL, m_inputs.NEXTRA, NULL,
    m_inputs.STATEFILE, m_inputs.SPIN.get(),
    &nregions, &neval, &fail, integral.data(), error.data(), prob.data());

  if(m_inputs.VERBOSE>0){
    printf("DIVONNE RESULT:\tnregions %d\tneval %d\tfail %d\n", nregions, neval, fail);
    for( comp = 0; comp < m_inputs.NCOMP; ++comp ) 
      printf("DIVONNE RESULT:\t%.8f +- %.8f\tp = %.3f\n", (double)integral[comp], (double)error[comp], (double)prob[comp]);
  }

  return integral[0];
}


// ============================================================================


double cbl::wrapper::cuba::CUBAwrapper::IntegrateCuhre (vector<vector<double>> integration_limits, const bool parallelize)
{
  int comp, neval, fail, nregions;
  vector<double> integral(m_inputs.NCOMP, 0), error(m_inputs.NCOMP, 0), prob(m_inputs.NCOMP, 0);

  cbl::wrapper::cuba::STR_CUBA_integrand *userdata = new cbl::wrapper::cuba::STR_CUBA_integrand;
  userdata->func = m_integrand;
  userdata->integration_limits = integration_limits;

  if (!parallelize) {
    int zero = 0;
    cubacores(&zero, &zero);
  }

  Cuhre(m_ndim, m_inputs.NCOMP, cbl::wrapper::cuba::CUBAIntegrand, userdata, m_inputs.NVEC,
    m_inputs.EPSREL, m_inputs.EPSABS, m_inputs.VERBOSE | m_inputs.LAST,
    m_inputs.MINEVAL, m_inputs.MAXEVAL, m_inputs.KEY,
    m_inputs.STATEFILE, m_inputs.SPIN.get(),
    &nregions, &neval, &fail, integral.data(), error.data(), prob.data());

  if(m_inputs.VERBOSE>0){
    printf("CUHRE RESULT:\tnregions %d\tneval %d\tfail %d\n", nregions, neval, fail);
    for( comp = 0; comp < m_inputs.NCOMP; ++comp )
      printf("CUHRE RESULT:\t%.8f +- %.8f\tp = %.3f\n", (double)integral[comp], (double)error[comp], (double)prob[comp]);
  }

  return integral[0];
}
