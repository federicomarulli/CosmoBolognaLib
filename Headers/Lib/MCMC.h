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
 *  @file Headers/Lib/MCMC.h
 *
 *  @brief The class MCMC 
 *
 *  This file defines the interface of the class MCMC, used for the
 *  Monte Carlo Markov Chain (MCMC) method
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#ifndef __MCMC__
#define __MCMC__

#include "MCMCg.h"

// ======================================================================================


namespace cosmobl {

  /**
   *  @class MCMC MCMC.h "Headers/Lib/MCMC.h"
   *
   *  @brief The class MCMC
   *
   *  This class is used to handle objects of type <EM> MCMC
   *  </EM>. It is used for Monte Carlo Markov Chain analyses.
   */
  class MCMC {
  
  protected:
 
    int m_npar;
  
    vector<vector<double> > m_par;
    vector<double> m_plog, m_chi2, m_accept;
    vector<cosmobl::glob::Prior> m_prior;

  public:

    MCMC (int npar) : m_npar(npar) 
    { m_par.erase(m_par.begin(),m_par.end()); m_par.resize(npar); } 

    MCMC (vector<cosmobl::glob::Prior> _prior) : m_npar(_prior.size()), m_prior(_prior)
     { 
	m_par.erase(m_par.begin(),m_par.end()); m_par.resize(m_npar); 
     } 

    void start_process (int, double , int , int , cosmobl::Chi2 &, vector<double> &);
    double prior_eval(vector<double> &);
    void write_chain (string);

  };
}

#endif
