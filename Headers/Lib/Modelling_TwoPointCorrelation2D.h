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
 *  @file Headers/Lib/Modelling_TwoPointCorrelation2D.h
 *
 *  @brief The class Modelling_TwoPointCorrelation2D
 *
 *  This file defines the interface of the class
 *  Modelling_TwoPointCorrelation2D, that contains all the common
 *  methods to model 2D two-point correlation functions
 *
 *  @author Federico Marulli, Alfonso Veropalumbo
 *
 *  @author federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#ifndef __MODELLING2D__
#define __MODELLING2D__


#include "Modelling_TwoPointCorrelation.h"


// ===================================================================================================


namespace cosmobl {
  
  namespace modelling {
    
    /**
     *  @class Modelling_TwoPointCorrelation2D
     *  Modelling_TwoPointCorrelation2D.h
     *  "Headers/Lib/Modelling_TwoPointCorrelation2D.h"
     *
     *  @brief The class Modelling_TwoPointCorrelation2D
     *
     *  This file defines the interface of the base class
     *  Modelling_TwoPointCorrelation2D, used for modelling
     *  the 2D two-point correlation function in cartesian coordinates
     *
     */
    class Modelling_TwoPointCorrelation2D : public Modelling_TwoPointCorrelation {

    public:

      /**
       *  @name Constructors/destructors
       */
      ///@{

      /**
       *  @brief default constuctor
       *  @return object of class Modelling_TwoPointCorrelation2D
       */
      Modelling_TwoPointCorrelation2D () = default;

      /**
       *  @brief constructor
       *  
       *  @param twop the two-point correlation function to model
       *
       *  @return object of type Modelling_TwoPointCorrelation2D
       */
      Modelling_TwoPointCorrelation2D (const shared_ptr<cosmobl::twopt::TwoPointCorrelation> twop)
	: Modelling_TwoPointCorrelation(twop) {}
      
      /**
       *  @brief default destructor
       *  @return none
       */
      virtual ~Modelling_TwoPointCorrelation2D () = default;
	
      ///@}


      /**
       *  @brief set the fiducial model for dark matter two-point
       *  correlation function
       *
       *  @return none
       */
      virtual void set_fiducial_xiDM () = 0;
      
      /**
       *  @brief set the parameters used to model the 2D two-point
       *  correlation function, in Cartesian coordinates
       *
       *  @param alpha_perp_guess guess value for the parameter
       *  \f$\alpha_\perp=\frac{D_{\rm A,1}(z)}{D_{\rm A,2}(z)}\f$
       *
       *  @param alpha_perp_prior prior for the parameter
       *  \f$\alpha_\perp=\frac{D_{\rm A,1}(z)}{D_{\rm A,2}(z)}\f$
       *
       *  @param alpha_perp_type type of the parameter
       *  \f$\alpha_\perp=\frac{D_{\rm A,1}(z)}{D_{\rm A,2}(z)}\f$:
       *  it can be either _free_ or _fixed_
       *
       *  @param alpha_par_guess guess value for the parameter
       *  \f$\alpha_\parallel=\frac{H_2(z)}{H_1(z)}\f$
       *
       *  @param alpha_par_prior prior for the parameter
       *  \f$\alpha_\parallel=\frac{H_2(z)}{H_1(z)}\f$
       *
       *  @param alpha_par_type type of the parameter
       *  \f$\alpha_\parallel=\frac{H_2(z)}{H_1(z)}\f$: it can be
       *  either _free_ or _fixed_
       *
       *  @param fsigma8_guess guess value for the parameter
       *  \f$f(z)\sigma_8(z)\f$
       *
       *  @param fsigma8_prior prior for the parameter
       *  \f$f(z)\sigma_8(z)\f$
       *
       *  @param fsigma8_type type of the parameter
       *  \f$f(z)\sigma_8(z)\f$: it can be either _free_ or _fixed_
       *
       *  @param bsigma8_guess guess value for the parameter
       *  \f$b(z)\sigma_8(z)\f$
       *
       *  @param bsigma8_prior prior for the parameter
       *  \f$b(z)\sigma_8(z)\f$
       *
       *  @param bsigma8_type type of the parameter
       *  \f$b(z)\sigma_8(z)\f$: it can be either _free_ or _fixed_
       *
       *  @param sigma12_guess guess value for the parameter
       *  \f$\sigma_{12}(z)\f$
       *
       *  @param sigma12_prior prior for the parameter \f$\sigma_{12}(z)\f$
       *
       *  @param sigma12_type type of the parameter
       *  \f$\sigma_{12}(z)\f$: it can be either _free_ or _fixed_
       *
       *  @return none
       */
      void set_model_dispersionModel (const double alpha_perp_guess=par::defaultDouble, const statistics::Prior alpha_perp_prior={}, const statistics::ParameterType alpha_perp_type=statistics::_free_, const double alpha_par_guess=par::defaultDouble, const statistics::Prior alpha_par_prior={}, const statistics::ParameterType alpha_par_type=statistics::_free_, const double fsigma8_guess=par::defaultDouble, const statistics::Prior fsigma8_prior={}, const statistics::ParameterType fsigma8_type=statistics::_free_, const double bsigma8_guess=par::defaultDouble, const statistics::Prior bsigma8_prior={}, const statistics::ParameterType bsigma8_type=statistics::_free_, const double sigma12_guess=par::defaultDouble, const statistics::Prior sigma12_prior={}, const statistics::ParameterType sigma12_type=statistics::_free_);


      /**
       *  @brief overloading of the function used to set the
       *  parameters used to model the 2D two-point correlation
       *  function, in Cartesian coordinates
       *
       *  @param fsigma8_prior prior for the parameter
       *  \f$f(z)\sigma_8(z)\f$
       *
       *  @param bsigma8_prior prior for the parameter
       *  \f$b(z)\sigma_8(z)\f$
       *
       *  @param sigma12_prior prior for the parameter
       *  \f$\sigma_{12}(z)\f$
       *
       *  @return none
       */
      void set_model_dispersionModel (const statistics::Prior fsigma8_prior, const statistics::Prior bsigma8_prior, const statistics::Prior sigma12_prior)
      {
	set_model_dispersionModel(1., {}, statistics::_fixed_, 1., {}, statistics::_fixed_, par::defaultDouble, fsigma8_prior, statistics::_free_, par::defaultDouble, bsigma8_prior, statistics::_free_, par::defaultDouble, sigma12_prior, statistics::_free_);
      }
      
      
      /**
       *  @brief compute and write the model 
       *
       *  @param xx vector of points at which the model is computed,
       *  first axis
       *  @param yy vector of points at which the model is computed,
       *  second axis
       *  @param output_dir the output directory
       *  @param output_file the output file
       *
       *  @param parameter vector containing input parameters used to
       *  compute the model; if it is not provided, internal
       *  parameters will be used
       *
       *  @return none
       */
      void write_model (const vector<double> xx, const vector<double> yy, const string output_dir, const string output_file, const vector<double> parameter={}) const override
      { m_model->write_model(xx, yy, output_dir, output_file, parameter); }
      
    };
  }
}

#endif
