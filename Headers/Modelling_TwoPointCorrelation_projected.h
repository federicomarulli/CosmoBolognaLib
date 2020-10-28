/********************************************************************
 *  Copyright (C) 2016 by Federico Marulli and Alfonso Veropalumbo  *
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
 *  @file Headers/Modelling_TwoPointCorrelation_projected.h
 *
 *  @brief The class Modelling_TwoPointCorrelation_projected
 *
 *  This file defines the interface of the class
 *  Modelling_TwoPointCorrelation_projected, used to model the
 *  projected two-point correlation function
 *
 *  @author Federico Marulli, Alfonso Veropalumbo
 *
 *  @author federico.marulli3@unibo.it, alfonso.veropalumbo@unibo.it
 */

#ifndef __MODELLINGTWOPOINTPROJ__
#define __MODELLINGTWOPOINTPROJ__


#include "Modelling_TwoPointCorrelation1D_monopole.h"
#include "ModelFunction_TwoPointCorrelation_projected.h"


// ===================================================================================================


namespace cbl {

  namespace modelling {

    namespace twopt {

      /**
       *  @class Modelling_TwoPointCorrelation_projected
       *  Modelling_TwoPointCorrelation_projected.h
       *  "Headers/Modelling_TwoPointCorrelation_projected.h"
       *
       *  @brief The class Modelling_TwoPointCorrelation_projected
       *
       *  This file defines the interface of the base class
       *  Modelling_TwoPointCorrelation_projected, used to model the
       *  projected of the two-point correlation function
       *
       */
      class Modelling_TwoPointCorrelation_projected : public Modelling_TwoPointCorrelation1D_monopole {
	
	public:

	  /**
	   *  @name Constructors/destructors
	   */
	  ///@{

	  /**
	   *  @brief default constuctor
	   */
	  Modelling_TwoPointCorrelation_projected () = default;

	  /**
	   *  @brief constructor
	   *  
	   *  @param twop the two-point correlation function to model
	   */
	  Modelling_TwoPointCorrelation_projected (const std::shared_ptr<cbl::measure::twopt::TwoPointCorrelation> twop)
	    : Modelling_TwoPointCorrelation1D_monopole(twop) {}

	  /**
	   *  @brief constructor
	   *  
	   *  @param twop_dataset the dataset containing the two-point
	   *  correlation function to model
	   */
	  Modelling_TwoPointCorrelation_projected (const std::shared_ptr<data::Data> twop_dataset)
	    : Modelling_TwoPointCorrelation1D_monopole() { set_data(twop_dataset); }

	  /**
	   *  @brief default destructor
	   *  
	   */
	  ~Modelling_TwoPointCorrelation_projected () = default;

	  ///@}


	  /**
	   *  @name Member functions used to set the model parameters
	   */
	  ///@{

	  /**
	   *  @brief set the fiducial model for the dark matter
	   *  projected correlation function
	   *
	   *  
	   */
	  void set_fiducial_wpDM ();
	
	  /**
	   *  @brief set the model to fit the projected two-point
	   *  correlation function assuming a linear bias
	   *
	   *  the model is the following:
	   *
	   *  \f[w_p(r_p) = (b\sigma_8)^2 \cdot w_p^{\rm DM}(r_p)\f]
	   *
	   *  @param bsigma8_prior prior for the parameter
	   *  \f$b(z)\sigma_8(z)\f$
	   *
	   *  
	   */
	  void set_model_linearBias (const statistics::PriorDistribution bsigma8_prior);

	  /**
	   *  @brief set the HOD parameters used to model the full
	   *  shape of the projected two-point correlation function
	   *
	   *  the HOD model for \f$w_p(r_p)\f$ is the one implemented
	   *  by cbl::modelling::twopt::wp_HOD
	   *
	   *  the model has 5 free parameters:
	   *
	   *  - \f$M_{min}\f$: the mass scale at which 50% of haloes
	   *   host a satellite galaxy
	   *
	   *  - \f$\sigma_{\log M_h}\f$: transition width reflecting
	   *   the scatter in the luminosity-halo mass relation
	   *
	   *  - \f$M_0\f$: the cutoff mass
	   *
	   *  - \f$M_1\f$: the amplitude of the power law
	   *
	   *  - \f$\alpha\f$: the slope of the power law
	   *
	   *  @param Mmin_prior \f$M_{min}\f$ prior
	   *
	   *  @param sigmalgM_prior \f$\sigma_{\log M_h}\f$ prior
	   *
	   *  @param M0_prior \f$M_0\f$ prior
	   *
	   *  @param M1_prior \f$\alpha\f$ prior
	   *
	   *  @param alpha_prior \f$\alpha\f$ prior
	   *
	   *  
	   */
	  void set_model_HOD (const statistics::PriorDistribution Mmin_prior={}, const statistics::PriorDistribution sigmalgM_prior={}, const statistics::PriorDistribution M0_prior={}, const statistics::PriorDistribution M1_prior={}, const statistics::PriorDistribution alpha_prior={});
	  
	  ///@}

      };
    }
  }
}

#endif
