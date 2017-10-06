/********************************************************************
 *  Copyright (C) 2017 by Federico Marulli                          *
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
 *  @file
 *  Headers/Lib/Modelling_ThreePointCorrelation_comoving_reduced.h
 *
 *  @brief The class Modelling_ThreePointCorrelation_comoving_reduced
 *
 *  This file defines the interface of the class
 *  Modelling_ThreePointCorrelation_comoving_reduced, that contains
 *  all the methods to model the connected three-point correlation
 *  function in angular coordinates
 *
 *  @author Federico Marulli
 *
 *  @author federico.marulli3@unbo.it
 */

#ifndef __MODELLINGTHREEPOINTCOMRED__
#define __MODELLINGTHREEPOINTCOMRED__


#include "Modelling_ThreePointCorrelation_comoving_connected.h"
#include "ModelFunction_ThreePointCorrelation_comoving_reduced.h"


// ===================================================================================================


namespace cosmobl {

  namespace modelling {
    
    namespace threept {
      
      /**
       *  @class Modelling_ThreePointCorrelation_comoving_reduced
       *  Modelling_ThreePointCorrelation_comoving_reduced.h
       *  "Headers/Lib/Modelling_ThreePointCorrelation_comoving_reduced.h"
       *
       *  @brief The class
       *  Modelling_ThreePointCorrelation_comoving_reduced
       *
       *  This file defines the interface of the class
       *  Modelling_ThreePointCorrelation_comoving_reduced, that
       *  contains all the methods to model the connected three-point
       *  correlation function in angular coordinates
       *
       */
      class Modelling_ThreePointCorrelation_comoving_reduced : public Modelling_ThreePointCorrelation_comoving_connected {

      public:

	/**
	 *  @name Constructors/destructors
	 */
	///@{

	/**
	 *  @brief default constuctor
	 *  @return object of class
	 *  Modelling_ThreePointCorrelation_comoving_reduced
	 */
	Modelling_ThreePointCorrelation_comoving_reduced () = default;

	/**
	 *  @brief constructor
	 *  
	 *  @param threep the three-point correlation function to
	 *  model
	 *
	 *  @return object of type
	 *  Modelling_ThreePointCorrelation_comoving_reduced
	 */
	Modelling_ThreePointCorrelation_comoving_reduced (const shared_ptr<cosmobl::measure::threept::ThreePointCorrelation> threep)
	  : Modelling_ThreePointCorrelation_comoving_connected(threep) {}

	/**
	 *  @brief constructor
	 *  
	 *  @param threep_dataset the dataset containing the
	 *  three-point correlation function to model
	 *
	 *  @return object of type
	 *  Modelling_ThreePointCorrelation_comoving_reduced
	 */
	Modelling_ThreePointCorrelation_comoving_reduced (const shared_ptr<data::Data> threep_dataset)
	  : Modelling_ThreePointCorrelation_comoving_connected() { set_data(threep_dataset); }
	
	/**
	 *  @brief default destructor
	 *  @return none
	 */
	virtual ~Modelling_ThreePointCorrelation_comoving_reduced () = default;
	
	///@}

	/**
	 *  @brief set the parameters used to model the full shape
	 *  of the monopole of the two-point correlation function
	 *
	 *  the model is the following:
	 *
	 *  \f[\xi_0(s) = \left[ (b\sigma_8)^2 + \frac{2}{3}
	 *  f\sigma_8 \cdot b\sigma_8 + \frac{1}{5}(f\sigma_8)^2
	 *  \right] \cdot \xi_{\rm DM}(\alpha\cdot s)/\sigma_8^2 +
	 *  \sum_{i=0}^N \frac{A_i}{r^{i}}\f]
	 *
	 *  the model has 3+N parameters: 
	 *    - \f$\alpha\f$
	 *    - \f$f(z)\sigma_8(z)\f$
	 *    - \f$b(z)\sigma_8(z)\f$
	 *    - \f$A_i\f$
	 *
	 *  the dark matter two-point correlation function is computed
	 *  using the input cosmological parameters
	 *
	 *  @param alpha_value if alpha_value>par::defaultDouble
	 *  then the \f$\alpha\f$ value is fixed at alpha_value
	 *
	 *  @param alpha_prior prior for the parameter \f$\alpha\f$
	 *
	 *  @param fsigma8_value if fsigma8_value>par::defaultDouble
	 *  then the \f$f\sigma_8\f$ value is fixed at fsigma8_value
	 *
	 *  @param fsigma8_prior prior for the parameter
	 *  \f$f(z)\sigma_8(z)\f$
	 *
	 *  @param bsigma8_value if bsigma8_value>par::defaultDouble
	 *  then the \f$b\sigma_8\f$ value is fixed at bsigma8_value
	 *
	 *  @param bsigma8_prior prior for the parameter
	 *  \f$b(z)\sigma_8(z)\f$
	 *
	 *  @param polynomial_value if
	 *  polynomial_value_value[i]>par::defaultDouble then the
	 *  the polynomial values are fixed at polynomial_value[i]
	 *
	 *  @param polynomial_prior vector containing the priors
	 *  for the polynomial part: the order of the polynomial is
	 *  the size of the vector \f$-1\f$
	 *
	 *  @param compute_xiDM true \f$rightarrow\f$ compute the
	 *  fiducial model of the dark matter correlation function
	 *
	 *  @return none
	 */
	void set_model_nonlinear_localbias (const statistics::Prior bias1_prior={}, const statistics::Prior bias2_prior={});

      
      };
    }
  }
}

#endif
