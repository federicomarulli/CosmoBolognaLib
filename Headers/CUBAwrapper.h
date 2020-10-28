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

      /// CUBA NDIM parameter
#define NDIM 2

      /// CUBA NCOMP parameter
#define NCOMP 1

      /// CUBA USERDATA parameter
#define USERDATA NULL

      /// CUBA NVEC parameter
#define NVEC 1

      /// CUBA EPSREL parameter
#define EPSREL 1e-4

      /// CUBA EPSABS parameter
#define EPSABS 1e-12

      /// CUBA VERBOSE parameter
#define VERBOSE 0

      /// CUBA LAST parameter
#define LAST 4

      /// CUBA SEED parameter
#define SEED 0

      /// CUBA MINEVAL parameter
#define MINEVAL 0

      /// CUBA MAXEVAL parameter
#define MAXEVAL 50000

      /// CUBA NSTART parameter
#define NSTART 1000

      /// CUBA NINCREASE parameter
#define NINCREASE 500

      /// CUBA NBATCH parameter
#define NBATCH 1000

      /// CUBA GRIDNO parameter
#define GRIDNO 0

      /// CUBA STATEFILE parameter
#define STATEFILE NULL

      /// CUBA SPIN parameter
#define SPIN NULL

      /// CUBA NNEW parameter
#define NNEW 1000

      /// CUBA NMIN parameter
#define NMIN 2

      /// CUBA FLATNESSv parameter
#define FLATNESS 25.

      /// CUBA KEY1 parameter
#define KEY1 47

      /// CUBA KEY2 parameter
#define KEY2 1

      /// CUBA KEY3 parameter
#define KEY3 1

      /// CUBA MAXPASS parameter
#define MAXPASS 5

      /// CUBA parameter
#define BORDER 0.

      /// CUBA MAXCHISQ parameter
#define MAXCHISQ 10.

      /// CUBA MINDEVIATION parameter
#define MINDEVIATION .25

      /// CUBA NGIVEN parameter
#define NGIVEN 0

      /// CUBA LDXGIVEN parameter
#define LDXGIVEN NDIM

      /// CUBA NEXTRA parameter
#define NEXTRA 0

      /// CUBA KEY parameter
#define KEY 0

    
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
	 *  default destructor
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
	 *  @brief integrate using the Vegas routine
	 *
	 *  @param integration_limits vector containing integration
	 *  limits
	 *
	 *  @return the integral
	 */
	double IntegrateVegas (std::vector<std::vector<double>> integration_limits);

	/**
	 *  @brief integrate using the Suave routine
	 *
	 *  @param integration_limits vector containing integration
	 *  limits
	 *
	 *  @return the integral
	 */
	double IntegrateSuave (std::vector<std::vector<double>> integration_limits);

	/**
	 *  @brief integrate using the Divonne routine
	 *
	 *  @param integration_limits vector containing integration
	 *  limits
	 *
	 *  @return the integral
	 */
	double IntegrateDivonne (std::vector<std::vector<double>> integration_limits);

	/**
	 *  @brief integrate using the Cuhre routine
	 *
	 *  @param integration_limits vector containing integration
	 *  limits
	 *
	 *  @return the integral
	 */
	double IntegrateCuhre (std::vector<std::vector<double>> integration_limits);

      };
    }
  }
}
 
#endif
