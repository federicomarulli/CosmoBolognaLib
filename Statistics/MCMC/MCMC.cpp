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
 *  @file MCMC/MCMC.cpp
 *
 *  @brief Methods of the class MCMC 
 *
 *  This file contains the implementation of the methods of the class
 *  MCMC, used for the Monte Carlo Markov Chain (MCMC) method
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo 
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#include "MCMC.h"
using namespace cosmobl;


// ======================================================================================


double cosmobl::MCMC::prior_eval (vector<double> &pr)
{
   double psum = 0.;
   for (int i=0; i<m_npar; i++)
      psum += -2.*log(m_prior[i](pr[i]));
   
   return psum;
}

// ======================================================================================


void cosmobl::MCMC::start_process (int idum, double logstep, int niter, int nsubiter, Chi2 &chi, vector<double> &start)
{
  if (int(start.size())!=m_npar) {
    string Err = "Error in start_process of cosmobl::MCMC::MCMC.cpp! start.size() != " + conv(m_npar,par::fINT) + "!";
    ErrorMsg(Err);
  }

  vector<double> temp(niter,0);

  m_par.erase(m_par.begin(), m_par.end()); m_par.resize(m_npar, temp);
  m_chi2.erase(m_chi2.begin(), m_chi2.end()); m_chi2.resize(niter, 0);
  m_accept.erase(m_accept.begin(), m_accept.end()); m_accept.resize(niter, 0);

  cosmobl::glob::State state;
  state.m_par = start;

  cosmobl::glob::Plog plog(chi, m_prior);

  cosmobl::glob::Proposal propose(idum, logstep);
   
  double acc = 1.;


  // write the best-fit
   
  for (int j=0; j<m_npar; j++) m_par[j][0] = state.m_par[j]; 
  m_chi2[0] = -2*plog(state);
  m_accept[0] = acc;


  for (int i=1; i<niter; i++) {
    acc = mcmc_step(nsubiter, state, plog, propose);
    for (int j=0; j<m_npar; j++) m_par[j][i] = state.m_par[j]; 
    m_chi2[i] = -2*state.m_plog;
    m_accept[i] = acc;
  }
  
  cout << "Accepted "<< 100*Average(m_accept) << "%" << endl;
  for (int i=0; i<m_npar; i++) {
    double ss = 0., err1 = 0., err2 = 0.;
    Quartile(m_par[i], err1, ss, err2);
    start[i] = ss;
    cout << start[i] << " - " << fabs(err1-start[i]) << " + " << fabs(err2-start[i]) << endl;
  }
   
}


// ======================================================================================


void cosmobl::MCMC::write_chain (string chain)
{
  ofstream fout (chain.c_str()); checkIO (chain,0);

  for (unsigned int i=0; i<m_chi2.size(); i++) {
    fout << i << " ";
    for (int j=0; j<m_npar; j++) 
      fout << m_par[j][i] << " " ;
    fout << m_chi2[i] << " " << m_accept[i] << endl;
  }

  fout.clear(); fout.close(); 
}


// ======================================================================================

/// @cond glob

double cosmobl::glob::mcmc_step (int mm, cosmobl::glob::State &state, cosmobl::glob::Plog &plog, cosmobl::glob::Proposal &propose)
{
  cosmobl::glob::State sprop;
  double qratio, alpha, ran;
  int accept = 0;

  plog(state);
  
  for (int i=0; i<mm; i++) {
    propose(state,sprop,qratio);
    alpha = min(1.,exp(plog(sprop)-state.m_plog));
    ran = propose.m_gau.doub();
    if (ran<alpha) {
      state = sprop;
      plog(state);
      accept ++;
    }
  }

  return accept/double(mm);
}

/// @endcond
