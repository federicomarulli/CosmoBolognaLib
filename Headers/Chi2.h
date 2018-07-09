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
 *  @file Headers/Chi2.h
 *
 *  @brief The class Chi2 
 *
 *  This file defines the interface of the class Chi2, used for
 *  statistical analyses and Bayesian inference
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#ifndef __CHI2__
#define __CHI2__

#include "Likelihood.h"


// ============================================================================================


namespace cbl {

  namespace statistics {
    
    /**
     *  @class Chi2 Chi2.h "Headers/Chi2.h"
     *
     *  @brief The class Chi2
     *
     *  This class is used to handle objects of type chi2. It is
     *  used for all kind of chi2 minimization
     */
    class Chi2 : public Likelihood
    {

    public:

      /**
       *  @name Constructors/destructors
       */
      ///@{

      /**
       *  @brief default constructor
       *  @return object of class Chi2
       */
      Chi2 () : Likelihood() {}

      /**
       *  @brief constructor
       *  
       *  @param data pointers to the data container
       *  
       *  @param model pointers to the model 
       *
       *  @param x_index index(s) of the extra info std::vector containing the point(s) where to evaluate the model
       *
       *  @param w_index index of the extra info std::vector containing the data point weight 
       *
       *  @return object of class Chi2
       */
      Chi2 (const std::shared_ptr<data::Data> data, const std::shared_ptr<Model> model, const std::vector<size_t> x_index={0,2}, const int w_index=-1);

      /**
       *  @brief default destructor
       *
       *  @return none
       */
      ~Chi2 () = default;

      ///@}

      /**
       * @brief evaluate the \f$\chi^2\f$
       *
       * @param pp the parameters
       *
       * @return the \f$\chi^2\f$
       */
      double operator () (std::vector<double> &pp) const;

      /**
       *  @brief function that minimizes the \f$\chi^2\f$, finds the
       *  best-fit parameters and stores them in model parameters
       *
       *  @param start std::vector containing initial values for
       *  the likelihood maximization
       *
       *  @param parameter_limits limits for the parameters
       *
       *  @param use_covariance if true use the full data covariance matrix
       *  else use the diagonal
       *
       *  @param max_iter the maximum number of iterations
       *
       *  @param tol the tolerance in finding convergence 
       *
       *  @param epsilon the relative fraction of the interval size
       *
       *  @return none
       */
      void minimize (const std::vector<double> start, const std::vector<std::vector<double>> parameter_limits, const bool use_covariance=false, const unsigned int max_iter=10000, const double tol=1.e-6, const double epsilon=1.e-3); 

    };
  }
}

#endif
