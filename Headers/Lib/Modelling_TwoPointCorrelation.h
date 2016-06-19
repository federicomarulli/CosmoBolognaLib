/********************************************************************
 *  Copyright (C) 2010 by Federico Marulli                          *
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
 *  @file Headers/Lib/Modelling_TwoPointCorrelation.h
 *
 *  @brief The class Modelling_TwoPointCorrelation
 *
 *  This file defines the interface of the class
 *  Modelling, used for modelling any kind of 
 *  2pcf measurements
 *
 *  @author Federico Marulli
 *
 *  @author federico.marulli3@unbo.it
 */

#ifndef __MODELLING2P__
#define __MODELLING2P__


#include "Cosmology.h"
#include "TwoPointCorrelation.h"
#include "Modelling.h"
#include "ModelBias.h"



// ===================================================================================================


namespace cosmobl {

  /**
   *  @brief The namespace of functions and classes used for modelling
   *  
   * The \e modelling namespace contains all the functions and classes
   * used to model any kind of measurements
   */
  namespace modelling {

    /**
     *  @class Modelling_TwoPointCorrelation Modelling_TwoPointCorrelation.h
     *  "Headers/Lib/Modelling_TwoPointCorrelation.h"
     *
     *  @brief The class Modelling_TwoPointCorrelation
     *
     *  This file defines the interface of the base class Modelling_TwoPointCorrelation,
     *  used for modelling any kind of measurements
     *
     */
    class Modelling_TwoPointCorrelation : public Modelling 
    {
      protected:
	
	/// two-point correlation function type
	twopt::TwoPType m_twoPType;

	/// The redshift of the two point correlation
	double m_redshift;

	/// The fiducial cosmology
	shared_ptr<Cosmology> m_cosmology;

      public:

	/**
	 *  @name Constructors/destructors
	 */
	///@{

	/**
	 *  @brief default constuctor
	 *  @return object of class ModellingTwoPointCorrelation
	 */
	Modelling_TwoPointCorrelation () {}

	/**
	 *  @brief default destructor
	 *  @return none
	 */
	virtual ~Modelling_TwoPointCorrelation () {}

	/**
	 *  @brief static factory used to construct modelling of 
	 *  two-point correlation functions of any type
	 *
	 *  @param twop the two-point correlation function to model
	 *
	 *  @param redshift the redhisft of the two-point correlation
	 *  signal
	 *
	 *  @param cosmology the fiducial cosmology
	 *
	 *  @return a pointer to an object of class 
	 *  Modelling_TwoPointCorrelation of a given type
	 */
	static shared_ptr<Modelling_TwoPointCorrelation> Create (const shared_ptr<cosmobl::twopt::TwoPointCorrelation> twop, const double redshift, const Cosmology cosmology);

	///@}

	/**
	 * @brief return the type of correlation function
	 * @return the type of correlation function
	 */
	twopt::TwoPType twoPType () {return m_twoPType;}

	/**
	 * @brief
	 *
	 * @param likelihoodType the likelihood type
	 * @param xlimits range for the fit
	 * @param bias_value prior for the bias 
	 * @param bias_prior prior for the bias
	 * @param nChains number of chains
	 * @param chain_size size of the chain
	 * @param dir_output output directory for the chains
	 * @param start starting position in the chains 
	 * @param stop final position in the chains
	 * @param thin interval of parameter in output
	 *
	 * @return
	 */
	virtual void fit_bias(const statistics::LikelihoodType likelihoodType, const vector<double> xlimits, const double bias_value, const statistics::Prior bias_prior, const int nChains, const int chain_size, const string dir_output, const double start=0.5, const double stop=1, const int thin=1)
	{ErrorMsg("Error in fit_bias of Modelling_TwoPointCorrelation.h");}

	/**
	 * @brief
	 *
	 * @param likelihoodType the likelihood type
	 * @param xlimits range for the fit
	 * @param bias_value value of the bias
	 * @param bias_prior prior for the bias
	 * @param CosmoPars vector containing the cosmological parameters
	 * @param prior_CosmoPars vector containing the priors on cosmological parameters
	 * @param nChains number of chains
	 * @param chain_size size of the chain
	 * @param dir_output output directory for the chains
	 * @param start starting position in the chains
	 * @param stop final position in the chains
	 * @param thin interval of parameter in output
	 *
	 * @return
	 */
	virtual void fit_bias_cosmology(const statistics::LikelihoodType likelihoodType, const vector<double> xlimits, const double bias_value, const statistics::Prior bias_prior, const vector<CosmoPar> CosmoPars, const vector<statistics::Prior> prior_CosmoPars, const int nChains, const int chain_size, const string dir_output, const double start=0.5, const double stop=1, const int thin=1)
	{ErrorMsg("Error in fit_bias_cosmology of Modelling_TwoPointCorrelation.h");}

//	virtual void fit_BAO(const string LikelihoodType, const vector<double> xlimits, const double bias_value, const statistics::Prior bias_prior, const int nChains, const int chain_size, const string dir_output, const double start=0.5, const double stop=1, const int thin=1, const int nbins_model, )
//	{ErrorMsg("Error in fit_BAO_cosmology of Modelling_TwoPointCorrelation.h");}
    };
  }
}

#endif
