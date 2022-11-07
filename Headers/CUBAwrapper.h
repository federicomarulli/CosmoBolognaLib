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
 *  @file Headers/CUBAwrapper.h
 *
 *  @brief class CUBAwrapper that wrap CUBA routines for 
 *  multidimensional integration
 *
 *  This file contains the wrappers of CUBA routines for 
 *  multidimensional integration
 *
 *  @author Alfonso Veropalumbo 
 *
 *  @author alfonso.veropalumbo@unibo.it
 */

#ifndef __CUBAwrap__
#define __CUBAwrap__

#include "cuba.h"
#include "Kernel.h"

namespace cbl {

  namespace wrapper {

    /**
     *  @brief The namespace of the <B> CUBA wrappers </B>
     *  
     *  The \e gsl namespace contains all the wrapper functions of the
     *  CUBA routines
     */
    namespace cuba {

      /**
       * @brief struct to manage cuba parameters
       */
      struct STR_CUBA_inputs {

	/// CUBA NDIM parameter
	int NDIM = 1;

	/// CUBA NCOMP parameter
	int NCOMP = 1;

	/// CUBA USERDATA parameter
	std::shared_ptr<void> USERDATA = NULL;

	/// CUBA NVEC parameter
	int NVEC = 1;

	/// CUBA EPSREL parameter
	double EPSREL = 1e-4;

	/// CUBA EPSABS parameter
	double EPSABS = 1e-12;

	/// CUBA VERBOSE parameter
	int VERBOSE = 0;

	/// CUBA LAST parameter
	int LAST = 4;

	/// CUBA SEED parameter
	int SEED = 0;

	/// CUBA MINEVAL parameter
	int MINEVAL = 0;

	/// CUBA MAXEVAL parameter
	int MAXEVAL = 50000;

	/// CUBA NSTART parameter
	int NSTART = 1000;

	/// CUBA NINCREASE parameter
	int NINCREASE = 500;

	/// CUBA NBATCH parameter
	int NBATCH = 1000;

	/// CUBA GRIDNO parameter
	int GRIDNO = 0;

	/// CUBA STATEFILE parameter
	char *STATEFILE = NULL;

	/// CUBA SPIN parameter
	std::shared_ptr<void> SPIN = NULL;

	/// CUBA NNEW parameter
	int NNEW = 1000;

	/// CUBA NMIN parameter
	int NMIN = 2;

	/// CUBA FLATNESS parameter
	double FLATNESS = 25.;

	/// CUBA KEY1 parameter
	int KEY1 = 47;

	/// CUBA KEY2 parameter
	int KEY2 = 1;

	/// CUBA KEY3 parameter
	int KEY3 = 1;

	/// CUBA MAXPASS parameter
	int MAXPASS = 5;

	/// CUBA parameter
	double BORDER = 0.;

	/// CUBA MAXCHISQ parameter
	double MAXCHISQ = 10.;

	/// CUBA MINDEVIATION parameter
	double MINDEVIATION = 0.25;

	/// CUBA NGIVEN parameter
	int NGIVEN = 0;

	/// CUBA LDXGIVEN parameter
	int LDXGIVEN = NDIM;

	/// CUBA NEXTRA parameter
	int NEXTRA = 0;

	/// CUBA KEY parameter
	int KEY = 0;
      };

    
      /**
       *  @brief support object for cuba integrand
       */
      struct STR_CUBA_integrand
      {
	/// function to be integrated
	FunctionDoubleVector func;

	/// limits of the integration
	std::vector<std::vector<double>> integration_limits;
      };

      /**
       *  @brief generic CUBA integrand
       *
       *  @param ndim integral dimension
       *  @param xx integral limits
       *  @param ncomp dimension of the function
       *  @param ff integrand
       *  @param userdata data for integrand evaluation
       *  @return 0
       */
      int CUBAIntegrand (const int *ndim, const cubareal xx[], const int *ncomp, cubareal ff[], void *userdata);

      /**
       *  @class CUBAwrapper CUBAwrapper.h "Headers/CUBAwrapper.h"
       *
       *  @brief The class CUBAwrapper
       *
       *  This class is used to handle objects of type <EM> CUBAwrapper
       *  </EM>. It can be used to estimate multidimensional integrals
       *  using montecarlo methods implemented in the CUBA library.  See
       *  http://www.feynarts.de/cuba/ for complete documentation
       */
      class CUBAwrapper {

      protected:

	/// integrand
	FunctionDoubleVector m_integrand;

	/// integral dimension
	int m_ndim;

	/// CUBA integration inputs
	STR_CUBA_inputs m_inputs;

      public:

	/**
	 *  @brief default constructor
	 */
	CUBAwrapper () = default;

	/**
	 *  @brief default constructor
	 *
	 *  @param func the integrand
	 *  @param function_parameters suppoer function parameters
	 *  @param parameters parameters
	 *  @param ndim the integral dimension
	 */
	CUBAwrapper (FunctionDoubleVectorPtrVectorRef func, const std::shared_ptr<void> function_parameters, std::vector<double> &parameters, const int ndim);

	/**
	 *  @brief default constructor
	 *  @param func the integrand
	 *  @param ndim the integral dimension
	 */
	CUBAwrapper (FunctionDoubleVector func, const int ndim);

	/**
	 *  @brief default destructor
	 */
	~CUBAwrapper () = default;

	/**
	 *  @brief set the integrand
	 *
	 *  @param func the integrand
	 *  @param function_parameters suppoer function parameters
	 *  @param parameters parameters
	 *  @param ndim the integral dimension
	 */
	void set_integrand (FunctionDoubleVectorPtrVectorRef func, const std::shared_ptr<void> function_parameters, std::vector<double> &parameters, const int ndim);

	/**
	 *  @brief set the integrand
	 *  @param func the integrand
	 *  @param ndim the integral dimension
	 */
	void set_integrand (FunctionDoubleVector func, const int ndim);
	
	/**
	 *  @brief set integration limits
	 *
	 *  @param integration_limits vector containing integration
	 *  limits
	 */
	void set_limits (std::vector<std::vector<double>> integration_limits)
	{ (void)integration_limits; ErrorCBL("", "set_limits", "CUBAwrapper.h", glob::ExitCode::_workInProgress_); };

	
	/**
	 * @brief return reference to integration parameters
	 * @return reference to integration parameters
	 */
	STR_CUBA_inputs& inputs () {return m_inputs;}

	/**
	 *  @brief integrate using the Vegas routine
	 *
	 *  @param integration_limits vector containing integration
	 *  limits
	 *  @param parallelize parallelize the integration
	 *
	 *  @return the integral
	 */
	double IntegrateVegas (std::vector<std::vector<double>> integration_limits, const bool parallelize=true);

	/**
	 *  @brief integrate using the Suave routine
	 *
	 *  @param integration_limits vector containing integration
	 *  limits
	 *  @param parallelize parallelize the integration
	 *
	 *  @return the integral
	 */
	double IntegrateSuave (std::vector<std::vector<double>> integration_limits, const bool parallelize=true);

	/**
	 *  @brief integrate using the Divonne routine
	 *
	 *  @param integration_limits vector containing integration
	 *  limits
	 *  @param parallelize parallelize the integration
	 *
	 *  @return the integral
	 */
	double IntegrateDivonne (std::vector<std::vector<double>> integration_limits, const bool parallelize=true);

	/**
	 *  @brief integrate using the Cuhre routine
	 *
	 *  @param integration_limits vector containing integration
	 *  limits
	 *  @param parallelize parallelize the integration
	 *
	 *  @return the integral
	 */
	double IntegrateCuhre (std::vector<std::vector<double>> integration_limits, const bool parallelize=true);

      };
    }
  }
}
 
#endif
