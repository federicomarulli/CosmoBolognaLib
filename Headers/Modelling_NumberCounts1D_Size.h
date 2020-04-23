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
 *  @file Headers/Modelling_NumberCounts1D_Size.h
 *
 *  @brief The class Modelling_NumberCounts1D_Size
 *
 *  This file defines the interface of the class Modelling_NumberCounts1D_Size, 
 *  used to  model number counts of any kind
 *
 *  @author Federico Marulli, Alfonso Veropalumbo
 *
 *  @author federico.marulli3@unibo.it, alfonso.veropalumbo@unibo.it
 */

#ifndef __MODELLINGNCS__
#define __MODELLINGNCS_


#include "Modelling_NumberCounts1D.h"


// ===================================================================================================


namespace cbl {

  namespace modelling {

    /**
     *  @brief The namespace of the <B> number counts
     *  modelling </B>
     *  
     *  The \e modelling::numbercounts namespace contains all the functions
     *  and classes to model number counts
     */
    namespace numbercounts {
    
      /**
       *  @class Modelling_NumberCounts1D_Size
       *  Modelling_NumberCounts1D_Size.h
       *  "Headers/Modelling_NumberCounts1D_Size.h"
       *
       *  @brief The class Modelling_NumberCounts1D_Size
       *
       *  This file defines the interface of the base class
       *  Modelling_NumberCounts1D_Size, used for modelling 
       *  redshift number counts measurements
       *
       */
      class Modelling_NumberCounts1D_Size : public Modelling_NumberCounts1D
      {

      public:

	/**
	 *  @name Constructors/destructors
	 */
	///@{

	/**
	 *  @brief default constuctor
	 *  @return object of class Modelling_NumberCounts1D_Size
	 */
	Modelling_NumberCounts1D_Size () = default;
	
	/**
	 *  @brief constuctor
	 *  @param nc the number counts to model
	 *  @return object of class Modelling_NumberCounts1D_Size
	 */
	Modelling_NumberCounts1D_Size (const std::shared_ptr<cbl::measure::numbercounts::NumberCounts> nc) 
	  : Modelling_NumberCounts1D (nc) {}
	
	/**
	 *  @brief constuctor
	 *  @param dataset the number counts dataset
	 *  @param hist_type the histogram type
	 *  @param fact the normalization factor
	 *  @return object of class Modelling_NumberCounts1D_Size
	 */
	Modelling_NumberCounts1D_Size (const std::shared_ptr<cbl::data::Data> dataset, glob::HistogramType hist_type, double fact)
	  : Modelling_NumberCounts1D (dataset, hist_type, fact) {}
	
	/**
	 *  @brief default destructor
	 *  @return none
	 */
	virtual ~Modelling_NumberCounts1D_Size () = default;
	
	///@}

	/**
	 *  @name Member functions used to set the model parameters
	 */
	///@{

	/**
	 *  @brief set the cosmological parameters used to model the 
	 *  size function
	 *
	 *  the model has N cosmological parameters
	 *
	 *  @param cosmo_params vector of enums containing cosmological
	 *  parameters
	 *
\	 *  @param cosmo_param_priors vector containing the priors for
	 *  the cosmological parameters
	 *
	 *  @return none
	 */
	void set_model_NumberCounts_cosmology (const std::vector<cbl::cosmology::CosmologicalParameter> cosmo_params={}, const std::vector<statistics::PriorDistribution> cosmo_param_priors={});

		///@{

	/**
	 *  @brief set the cosmological parameters and parameters
	 *  necessary to re-parametrise the void size function
	 *
	 *  the model has N cosmological parameters, the tracer bias
	 *  and the coefficients (b_slope and b_offset) of the linear
	 *  relation necessary to convert the tracer bias into the
	 *  value of the bias computed inside cosmic voids (see
	 *  Contarini et al. 2019)
	 *
	 *  @param cosmo_params vector of enums containing cosmological
	 *  parameters
	 *
	 *  @param cosmo_param_priors vector containing the priors for
	 *  the cosmological parameters
	 *
	 *  @param bias_priors the priors of the effective bias,
	 *  b_slope and b_offset
	 *
	 *  @return none
	 */
	void set_model_NumberCounts_cosmology_and_bias (const std::vector<cbl::cosmology::CosmologicalParameter> cosmo_params={}, const std::vector<statistics::PriorDistribution> cosmo_param_priors={}, const std::vector<statistics::PriorDistribution> bias_priors={});

	///@}

      };
    }
  }
}

#endif
