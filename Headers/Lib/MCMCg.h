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
 *  @file Headers/Lib/MCMCg.h
 *
 *  @brief Structures of the class MCMC 
 *
 *  This file contains some structures used for the Monte Carlo Markov
 *  Chain (MCMC) method
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#ifndef __MCMCg__
#define __MCMCg__

#include "Chi2.h"

// ======================================================================================


namespace cosmobl {

  namespace glob {
    
    class Prior {
      
    protected:
      
      typedef double(*prior_func) (double, void *, vector<double>);

      prior_func func;
      vector<double> prior_func_pars;
      double xmin, xmax;

    public:

      Prior () {
	xmin = -1.e30;
	xmax = -1.e30;
	func = &identity;
      }

      Prior (vector<double> _Limits) {
	xmin = _Limits[0]; 
	xmax = _Limits[1];
	func = &identity;  
      }

      Prior (vector<double> _Limits, double &mean, double &sigma) {
	xmin = _Limits[0]; xmax = _Limits[1];
	prior_func_pars.push_back(mean);
	prior_func_pars.push_back(sigma);
	func = &gaussian;
      }

      Prior (vector<double> _Limits, prior_func _func, vector<double> _prior_func_pars) 
	: func(_func), prior_func_pars(_prior_func_pars) 
      {
	xmin = _Limits[0]; 
	xmax = _Limits[1];
      }

      double operator() (double _xx) {
	void *pp = NULL;
	if (_xx <xmin || _xx >xmax) return 1.e30;
	else return func(_xx, pp, prior_func_pars);
      }

      void put_limits(vector<double> &Limits) {
	xmin = Limits[0]; 
	xmax = Limits[1];
      }

      void put_gaussian_parameter(double &mean, double &sigma) {
	func = &gaussian;
	prior_func_pars.erase(prior_func_pars.begin(), prior_func_pars.end());
	prior_func_pars.push_back(mean);
	prior_func_pars.push_back(sigma);
      }

      void put_func_parameter (prior_func _func, vector<double> pars) {
	func = _func;
	prior_func_pars.erase(prior_func_pars.begin(), prior_func_pars.end());
	for (unsigned int i=0; i<pars.size(); i++) 
	  prior_func_pars.push_back(pars[i]);
      }
    
    };


    // ======================================================================================

    
    struct State {
  
      vector<double> m_par;
      double m_plog, m_prob;

      State (vector<double> par) 
      : m_par(par), m_prob(1.) 
      {  
	for (unsigned int i=0; i<par.size(); i++) m_prob *= abs(par[i]);
      }

      State() {};
    };


    // ======================================================================================


    struct Plog
    {
      cosmobl::Chi2 &m_chi2;
      vector<Prior> m_prior;

      Plog (cosmobl::Chi2 &chi2) : m_chi2(chi2) {}
      Plog (Chi2 &_chi2, vector<Prior> _prior) : m_chi2(_chi2), m_prior(_prior) {}

      double operator() (State &state) {
	vector<double> spar = state.m_par;
	state.m_prob = 0;

	for (unsigned int i=0; i<m_prior.size(); i++) state.m_prob +=log(m_prior[i](spar[i]));

	return (state.m_plog=-0.5*m_chi2.get_chi2(spar)+state.m_prob);
      }
    };


    // ======================================================================================


    struct Proposal
    {
      Normaldev m_gau;
      double m_lstep;

      Proposal (int randseed, double logstep) 
      : m_gau(0.,1.,randseed), m_lstep(logstep) {}
  
      void operator() (const State &s1, State &s2, double &qratio) {
	vector<double> tPar = s1.m_par;
    
	for (unsigned int i=0; i<s1.m_par.size(); i++) tPar[i] = s1.m_par[i]*exp(m_lstep*m_gau.dev());

	State ts(tPar);
	qratio = 1.;
	s2 = ts;
      }
    };


    // ======================================================================================

    /// @cond glob
    double mcmc_step (int, State &, Plog &, Proposal &);
    /// @endcond

  }
}

#endif
