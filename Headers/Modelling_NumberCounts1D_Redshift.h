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
 *  @file Headers/Modelling_NumberCounts1D_Redshift.h
 *
 *  @brief The class Modelling_NumberCounts1D_Redshift
 *
 *  This file defines the interface of the class Modelling_NumberCounts1D_Redshift, 
 *  used to  model number counts of any kind
 *
 *  @author Federico Marulli, Alfonso Veropalumbo
 *
 *  @author federico.marulli3@unibo.it, alfonso.veropalumbo@unibo.it
 */

#ifndef __MODELLINGNCR__
#define __MODELLINGNCR__


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
       *  @class Modelling_NumberCounts1D_Redshift
       *  Modelling_NumberCounts1D_Redshift.h
       *  "Headers/Modelling_NumberCounts1D_Redshift.h"
       *
       *  @brief The class Modelling_NumberCounts1D_Redshift
       *
       *  This file defines the interface of the base class
       *  Modelling_NumberCounts1D_Redshift, used for modelling 
       *  redshift number counts measurements
       *
       */
      class Modelling_NumberCounts1D_Redshift : public Modelling_NumberCounts1D 
      {

      public:

	/**
	 *  @name Constructors/destructors
	 */
	///@{

	/**
	 *  @brief default constuctor
	 *  @return object of class Modelling_NumberCounts1D_Redshift
	 */
	Modelling_NumberCounts1D_Redshift () = default;
	
	/**
	 *  @brief constuctor
	 *  @param nc the number counts to model
	 *  @return object of class Modelling_NumberCounts1D_Redshift
	 */
	Modelling_NumberCounts1D_Redshift (const std::shared_ptr<cbl::measure::numbercounts::NumberCounts> nc) 
	  : Modelling_NumberCounts1D (nc) {}
	
	/**
	 *  @brief constuctor
	 *  @param dataset the number counts dataset
	 *  @param hist_type the histogram type
	 *  @param fact the normalization factor
	 *  @return object of class Modelling_NumberCounts1D_Redshift
	 */
	Modelling_NumberCounts1D_Redshift (const std::shared_ptr<cbl::data::Data> dataset, glob::HistogramType hist_type, double fact)
	  : Modelling_NumberCounts1D (dataset, hist_type, fact) {}
	
	/**
	 *  @brief default destructor
	 *  @return none
	 */
	virtual ~Modelling_NumberCounts1D_Redshift () = default;

	///@}

	/**
	 *  @name Member functions used to set the model parameters
	 */
	///@{

	/**
	 *  @brief set the cosmological parameters used to model the 
	 *  mass function
	 *
	 *  the model has N cosmological parameters
	 *
	 *  @param cosmo_param vector of enums containing cosmological
	 *  parameters
	 *
	 *  @param cosmo_param_prior vector containing the priors for
	 *  the cosmological parameters
	 *
	 *  @return none
	 */
	void set_model_NumberCounts_cosmology (const std::vector<cbl::cosmology::CosmologicalParameter> cosmo_param={}, const std::vector<statistics::PriorDistribution> cosmo_param_prior={});

	///@}
      };
    }
  }
}

#endif
