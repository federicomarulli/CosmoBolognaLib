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
 *  @file Headers/Modelling_NumberCounts1D_Mass.h
 *
 *  @brief The class Modelling_NumberCounts1D_Mass
 *
 *  This file defines the interface of the class Modelling_NumberCounts1D_Mass, 
 *  used to  model number counts of any kind
 *
 *  @author Federico Marulli, Alfonso Veropalumbo
 *
 *  @author federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#ifndef __MODELLINGNCM__
#define __MODELLINGNCM__


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
       *  @class Modelling_NumberCounts1D_Mass
       *  Modelling_NumberCounts1D_Mass.h
       *  "Headers/Modelling_NumberCounts1D_Mass.h"
       *
       *  @brief The class Modelling_NumberCounts1D_Mass
       *
       *  This file defines the interface of the base class
       *  Modelling_NumberCounts1D_Mass, used for modelling 
       *  redshift number counts measurements
       *
       */
      class Modelling_NumberCounts1D_Mass : public Modelling_NumberCounts1D
      {

      public:

	/**
	 *  @name Constructors/destructors
	 */
	///@{

	/**
	 *  @brief default constuctor
	 *  @return object of class Modelling_NumberCounts1D_Mass
	 */
	Modelling_NumberCounts1D_Mass () = default;
	
	/**
	 *  @brief constuctor
	 *  @param nc the number counts to model
	 *  @return object of class Modelling_NumberCounts1D_Mass
	 */
	Modelling_NumberCounts1D_Mass (const std::shared_ptr<cbl::measure::numbercounts::NumberCounts> nc) 
	  : Modelling_NumberCounts1D (nc) {}
	
	/**
	 *  @brief constuctor
	 *  @param dataset the number counts dataset
	 *  @param hist_type the histogram type
	 *  @param fact the normalization factor
	 *  @return object of class Modelling_NumberCounts1D_Mass
	 */
	Modelling_NumberCounts1D_Mass (const std::shared_ptr<cbl::data::Data> dataset, glob::HistogramType hist_type, double fact)
	  : Modelling_NumberCounts1D (dataset, hist_type, fact) {}
	
	/**
	 *  @brief default destructor
	 *  @return none
	 */
	virtual ~Modelling_NumberCounts1D_Mass () = default;

	///@}


		
	/**
	 *  @name Member functions used to set the model parameters
	 */
	///@{

	/**
	 *  @brief set the data used to construct mass
	 *  number counts of simulation snapshots
         *  
	 *  @param cosmology the cosmological model used to compute
	 *  &xi;<SUB>DM</SUB>
	 *
	 *  @param redshift redshift
	 *
	 *  @param method_Pk method used to compute the power
	 *  spectrum; valid choices for method_Pk are: CAMB
	 *  [http://camb.info/], CLASS [http://class-code.net/],
	 *  MPTbreeze-v1 [http://arxiv.org/abs/1207.1465],
	 *  EisensteinHu
	 *  [http://background.uchicago.edu/~whu/transfer/transferpage.html]
	 *    
	 *  @param k_min minimum wave vector module up to which the
	 *  binned dark matter power spectrum is computed
	 *  
	 *  @param k_max maximum wave vector module up to which the
	 *  binned dark matter power spectrum is computed
	 *
	 *  @param step number of steps used to compute the binned
	 *  power spectrum
	 *
	 *  @param output_dir the output_dir directory
	 *  where the output of external codes are written
	 *  
	 *  @param norm 0 &rarr; don't normalize the power spectrum; 1
	 *  &rarr; normalize the power spectrum
	 *
	 *  @param Delta \f$\Delta\f$: the overdensity, defined as the
	 *  mean interior density relative to the background
	 *
	 *  @param isDelta_vir \f$\rightarrow\f$ \f$\Delta\f$ is the
	 *  virial overdensity
	 *
	 *  @param model_MF author(s) who proposed the mass function
	 *
	 *  @param Volume the snapshot volume
	 *
	 *  @param Mass_min minimum halo mass
	 *
	 *  @param Mass_max maximum halo mass
	 *
	 *  @param Mass_step the number of bins for mass vector
       	 *
	 *  @param prec the precision
	 *
	 *  @return none
	 */
	void set_data_model_snapshot (const cbl::cosmology::Cosmology cosmology={}, const double redshift=0., const std::string method_Pk="CAMB", const double k_min=1.e-4, const double k_max=100., const int step=500,  const std::string output_dir=par::defaultString, const int norm=-1, const double Delta=200., const bool isDelta_vir=true, const std::string model_MF="Tinker", const double Volume=par::defaultDouble, const double Mass_min=par::defaultDouble, const double Mass_max=par::defaultDouble, const int Mass_step=100, const double prec=1.e-4); 
	
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
