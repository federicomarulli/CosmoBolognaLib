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
 *  @file Headers/Lib/Prior.h
 *
 *  @brief The class Prior 
 *
 *  This file defines the interface of the class Prior
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#ifndef __PRIOR__
#define __PRIOR__

#include "Func.h"


// ============================================================================================


namespace cosmobl {

  /**
   *  @class Chi2 Chi2.h "Headers/Lib/Prior.h"
   *
   *  @brief The class Prior
   *
   *  This class is used to define the prior
   */
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

    Prior (vector<double>);

    Prior (vector<double>, double &, double &);

    Prior (vector<double>, prior_func, vector<double>);

    double operator() (double _xx) {
      void *pp = NULL;
      if (_xx<xmin || _xx>xmax) return 1.e30;
      else return func(_xx, pp, prior_func_pars);
    }

    void put_limits (vector<double>);

    void put_gaussian_parameters (double &, double &);
        
    void put_func_parameters (prior_func, vector<double>);

  };
}

#endif
